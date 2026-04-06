//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::should_implement_trait)]
use super::functions::*;
use crate::{Error, Result};
use oxiphysics_core::math::Vec3;
use std::io::{BufRead, BufReader, Read, Write};

/// Describes a custom column set for a LAMMPS dump.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DumpColumnSpec {
    /// Column names in order.
    pub columns: Vec<String>,
}
impl DumpColumnSpec {
    /// Standard full column set: `id type x y z vx vy vz`.
    pub fn standard() -> Self {
        Self {
            columns: vec![
                "id".into(),
                "type".into(),
                "x".into(),
                "y".into(),
                "z".into(),
                "vx".into(),
                "vy".into(),
                "vz".into(),
            ],
        }
    }
    /// Positions only: `id type x y z`.
    pub fn positions_only() -> Self {
        Self {
            columns: vec![
                "id".into(),
                "type".into(),
                "x".into(),
                "y".into(),
                "z".into(),
            ],
        }
    }
    /// Return the `ITEM: ATOMS` header line.
    pub fn header_line(&self) -> String {
        format!("ITEM: ATOMS {}", self.columns.join(" "))
    }
    /// Number of columns.
    pub fn len(&self) -> usize {
        self.columns.len()
    }
    /// True if no columns are defined.
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }
    /// Check if a column name is present.
    pub fn contains(&self, name: &str) -> bool {
        self.columns.iter().any(|c| c == name)
    }
    /// Index of a column by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|c| c == name)
    }
}
/// LAMMPS-style restart file: stores atom positions and velocities in a
/// compact binary format for simulation restart.
pub struct LammpsRestartWriter;
impl LammpsRestartWriter {
    /// Write restart data to `writer`.
    ///
    /// Layout:
    /// ```text
    /// [magic: b"LRST"] 4 bytes
    /// [version: i32 LE = 1]
    /// [timestep: i64 LE]
    /// [n_atoms: i64 LE]
    /// [n_types: i32 LE]
    /// [masses: n_types × f64 LE]
    /// for each atom: [id i32][type i32][x f64][y f64][z f64][vx f64][vy f64][vz f64]
    /// ```
    pub fn write<W: std::io::Write>(
        writer: &mut W,
        timestep: i64,
        masses: &[f64],
        atoms: &[LammpsAtom],
    ) -> Result<()> {
        writer.write_all(b"LRST")?;
        writer.write_all(&1_i32.to_le_bytes())?;
        writer.write_all(&timestep.to_le_bytes())?;
        writer.write_all(&(atoms.len() as i64).to_le_bytes())?;
        writer.write_all(&(masses.len() as i32).to_le_bytes())?;
        for &m in masses {
            writer.write_all(&m.to_le_bytes())?;
        }
        for a in atoms {
            writer.write_all(&(a.id as i32).to_le_bytes())?;
            writer.write_all(&(a.type_id as i32).to_le_bytes())?;
            writer.write_all(&a.position.x.to_le_bytes())?;
            writer.write_all(&a.position.y.to_le_bytes())?;
            writer.write_all(&a.position.z.to_le_bytes())?;
            writer.write_all(&a.velocity.x.to_le_bytes())?;
            writer.write_all(&a.velocity.y.to_le_bytes())?;
            writer.write_all(&a.velocity.z.to_le_bytes())?;
        }
        Ok(())
    }
}
/// Simulation box geometry for a LAMMPS run.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LammpsBox {
    /// Low corner `[xlo, ylo, zlo]`.
    pub lo: [f64; 3],
    /// High corner `[xhi, yhi, zhi]`.
    pub hi: [f64; 3],
    /// Tilt factors `[xy, xz, yz]` for triclinic cells (0 for orthogonal).
    pub tilt: [f64; 3],
}
#[allow(dead_code)]
impl LammpsBox {
    /// Create an orthogonal box with given extents.
    pub fn orthogonal(lo: [f64; 3], hi: [f64; 3]) -> Self {
        Self {
            lo,
            hi,
            tilt: [0.0; 3],
        }
    }
    /// Box lengths in each dimension.
    pub fn lengths(&self) -> [f64; 3] {
        [
            self.hi[0] - self.lo[0],
            self.hi[1] - self.lo[1],
            self.hi[2] - self.lo[2],
        ]
    }
    /// Volume of an orthogonal box.
    pub fn volume(&self) -> f64 {
        let l = self.lengths();
        l[0] * l[1] * l[2]
    }
    /// Number density given an atom count.
    pub fn number_density(&self, n_atoms: usize) -> f64 {
        let vol = self.volume();
        if vol < 1e-30 {
            return 0.0;
        }
        n_atoms as f64 / vol
    }
    /// Generate LAMMPS `region` command for a box region.
    pub fn to_region_command(&self, region_id: &str) -> String {
        format!(
            "region {} block {} {} {} {} {} {}\n",
            region_id, self.lo[0], self.hi[0], self.lo[1], self.hi[1], self.lo[2], self.hi[2],
        )
    }
    /// Check whether a position lies inside the box (for orthogonal boxes).
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        (0..3).all(|d| pos[d] >= self.lo[d] && pos[d] < self.hi[d])
    }
    /// Wrap a position into the primary box using periodic boundary conditions.
    pub fn wrap_pbc(&self, pos: [f64; 3]) -> [f64; 3] {
        let mut out = pos;
        for d in 0..3 {
            let l = self.hi[d] - self.lo[d];
            if l > 1e-30 {
                out[d] = out[d] - (((out[d] - self.lo[d]) / l).floor()) * l + self.lo[d];
            }
        }
        out
    }
}
/// Generator for common LAMMPS input-script sections.
#[allow(dead_code)]
pub struct LammpsInputScript;
#[allow(dead_code)]
impl LammpsInputScript {
    /// Generate a minimisation script section.
    ///
    /// * `force_tol`  — force convergence tolerance (eV/Å or kcal/mol/Å).
    /// * `energy_tol` — energy convergence tolerance.
    pub fn minimize_script(force_tol: f64, energy_tol: f64) -> String {
        format!("# Energy minimisation\nminimize {energy_tol} {force_tol} 1000 10000\n")
    }
    /// Generate an NPT run script section.
    ///
    /// * `temp`    — target temperature (K).
    /// * `press`   — target pressure (bar).
    /// * `t_damp`  — temperature damping parameter (time units).
    /// * `p_damp`  — pressure damping parameter (time units).
    /// * `n_steps` — number of MD steps.
    pub fn npt_script(temp: f64, press: f64, t_damp: f64, p_damp: f64, n_steps: u64) -> String {
        format!(
            "# NPT ensemble\nfix 1 all npt temp {temp} {temp} {t_damp} iso {press} {press} {p_damp}\nrun {n_steps}\nunfix 1\n"
        )
    }
}
#[allow(dead_code)]
impl LammpsInputScript {
    /// Generate a complete NVE run section.
    pub fn nve_script(n_steps: u64, dt: f64) -> String {
        let mut s = String::new();
        s.push_str(&format!("timestep {dt}\n"));
        s.push_str("fix 1 all nve\n");
        s.push_str(&format!("run {n_steps}\n"));
        s.push_str("unfix 1\n");
        s
    }
    /// Generate a thermalisation script with Langevin dynamics.
    pub fn thermalise_script(temp: f64, damp: f64, n_steps: u64, seed: u64) -> String {
        let mut s = String::new();
        s.push_str(&format!(
            "fix therm all langevin {temp} {temp} {damp} {seed}\n"
        ));
        s.push_str("fix nve_int all nve\n");
        s.push_str(&format!("run {n_steps}\n"));
        s.push_str("unfix therm\n");
        s.push_str("unfix nve_int\n");
        s
    }
    /// Generate a dump command.
    pub fn dump_command(
        dump_id: &str,
        group: &str,
        style: &str,
        every: u64,
        filename: &str,
        fields: &str,
    ) -> String {
        format!("dump {dump_id} {group} {style} {every} {filename} {fields}\n")
    }
    /// Generate a thermo output command.
    pub fn thermo_output(every: u64, keywords: &[&str]) -> String {
        let mut s = format!("thermo {every}\n");
        if !keywords.is_empty() {
            s.push_str(&format!("thermo_style custom {}\n", keywords.join(" ")));
        }
        s
    }
    /// Generate a units command.
    pub fn units(style: &str) -> String {
        format!("units {style}\n")
    }
    /// Generate a boundary command.
    pub fn boundary(style: &str) -> String {
        format!("boundary {style}\n")
    }
    /// Generate a read_data command.
    pub fn read_data(filename: &str) -> String {
        format!("read_data {filename}\n")
    }
}
/// Writer for LAMMPS binary dump format (little-endian).
///
/// Produces a compact binary representation of a single frame with the
/// standard `id type x y z vx vy vz` per-atom layout.
pub struct LammpsBinaryDumpWriter;
impl LammpsBinaryDumpWriter {
    /// Write one binary frame to `writer`.
    ///
    /// Layout per frame:
    /// ```text
    /// [timestep: i64 LE]
    /// [n_atoms:  i64 LE]
    /// [6 × f64 LE: xlo, xhi, ylo, yhi, zlo, zhi]
    /// for each atom:
    ///   [id:    i32 LE]
    ///   [type:  i32 LE]
    ///   [x y z vx vy vz: 6 × f64 LE]
    /// ```
    pub fn write_frame<W: std::io::Write>(
        writer: &mut W,
        timestep: i64,
        box_bounds: [[f64; 2]; 3],
        atoms: &[LammpsAtom],
    ) -> Result<()> {
        writer.write_all(&timestep.to_le_bytes())?;
        writer.write_all(&(atoms.len() as i64).to_le_bytes())?;
        for b in &box_bounds {
            writer.write_all(&b[0].to_le_bytes())?;
            writer.write_all(&b[1].to_le_bytes())?;
        }
        for a in atoms {
            writer.write_all(&(a.id as i32).to_le_bytes())?;
            writer.write_all(&(a.type_id as i32).to_le_bytes())?;
            writer.write_all(&a.position.x.to_le_bytes())?;
            writer.write_all(&a.position.y.to_le_bytes())?;
            writer.write_all(&a.position.z.to_le_bytes())?;
            writer.write_all(&a.velocity.x.to_le_bytes())?;
            writer.write_all(&a.velocity.y.to_le_bytes())?;
            writer.write_all(&a.velocity.z.to_le_bytes())?;
        }
        Ok(())
    }
    /// Byte size of one binary frame (header + atoms).
    pub fn frame_byte_size(n_atoms: usize) -> usize {
        8 + 8 + 6 * 8 + n_atoms * (4 + 4 + 6 * 8)
    }
}
/// Writer for the LAMMPS data-file format (used with `read_data`).
#[allow(dead_code)]
pub struct LammpsDataWriter;
#[allow(dead_code)]
impl LammpsDataWriter {
    /// Generate the header section of a LAMMPS data file.
    ///
    /// `box_lo` / `box_hi` are the lower/upper corners of the simulation box.
    pub fn write_header(
        n_atoms: usize,
        n_bonds: usize,
        box_lo: [f64; 3],
        box_hi: [f64; 3],
    ) -> String {
        let mut s = String::new();
        s.push_str("LAMMPS data file\n\n");
        s.push_str(&format!("{} atoms\n", n_atoms));
        if n_bonds > 0 {
            s.push_str(&format!("{} bonds\n", n_bonds));
        }
        s.push('\n');
        s.push_str(&format!("{} {} xlo xhi\n", box_lo[0], box_hi[0]));
        s.push_str(&format!("{} {} ylo yhi\n", box_lo[1], box_hi[1]));
        s.push_str(&format!("{} {} zlo zhi\n", box_lo[2], box_hi[2]));
        s
    }
    /// Generate the `Atoms` section for `atom_style atomic`.
    ///
    /// `types` is 1-indexed atom type per atom.
    pub fn write_atoms_atomic(positions: &[[f64; 3]], masses: &[f64], types: &[u32]) -> String {
        let _ = masses;
        let mut s = String::from("\nAtoms  # atomic\n\n");
        for (i, (pos, &t)) in positions.iter().zip(types.iter()).enumerate() {
            s.push_str(&format!(
                "{} {} {} {} {}\n",
                i + 1,
                t,
                pos[0],
                pos[1],
                pos[2]
            ));
        }
        s
    }
    /// Generate the `Atoms` section for `atom_style charge`.
    pub fn write_atoms_charge(
        positions: &[[f64; 3]],
        masses: &[f64],
        charges: &[f64],
        types: &[u32],
    ) -> String {
        let _ = masses;
        let mut s = String::from("\nAtoms  # charge\n\n");
        for (i, ((pos, &q), &t)) in positions
            .iter()
            .zip(charges.iter())
            .zip(types.iter())
            .enumerate()
        {
            s.push_str(&format!(
                "{} {} {} {} {} {}\n",
                i + 1,
                t,
                q,
                pos[0],
                pos[1],
                pos[2]
            ));
        }
        s
    }
    /// Generate the `Bonds` section.
    ///
    /// Each entry is `(bond_type, atom_i, atom_j)` with 1-indexed atoms.
    pub fn write_bonds(bonds: &[(u32, u32, u32)]) -> String {
        let mut s = String::from("\nBonds\n\n");
        for (idx, &(btype, i, j)) in bonds.iter().enumerate() {
            s.push_str(&format!("{} {} {} {}\n", idx + 1, btype, i, j));
        }
        s
    }
    /// Generate a complete LAMMPS data file string (atomic style, no bonds).
    pub fn write_complete(
        positions: &[[f64; 3]],
        masses: &[f64],
        types: &[u32],
        box_lo: [f64; 3],
        box_hi: [f64; 3],
    ) -> String {
        let mut s = Self::write_header(positions.len(), 0, box_lo, box_hi);
        let max_type = *types.iter().max().unwrap_or(&1) as usize;
        s.push_str(&format!("{} atom types\n", max_type));
        s.push_str("\nMasses\n\n");
        let mut type_mass = vec![1.0_f64; max_type + 1];
        for (&t, &m) in types.iter().zip(masses.iter()) {
            type_mass[t as usize] = m;
        }
        for ti in 1..=max_type {
            s.push_str(&format!("{} {}\n", ti, type_mass[ti]));
        }
        s.push_str(&Self::write_atoms_atomic(positions, masses, types));
        s
    }
}
#[allow(dead_code)]
impl LammpsDataWriter {
    /// Generate the `Atoms` section for `atom_style full`.
    ///
    /// Format: `id mol-id type charge x y z`
    pub fn write_atoms_full(
        positions: &[[f64; 3]],
        charges: &[f64],
        types: &[u32],
        mol_ids: &[u32],
    ) -> String {
        let mut s = String::from("\nAtoms  # full\n\n");
        for (i, (((&pos, &q), &t), &mol)) in positions
            .iter()
            .zip(charges.iter())
            .zip(types.iter())
            .zip(mol_ids.iter())
            .enumerate()
        {
            s.push_str(&format!(
                "{} {} {} {} {} {} {}\n",
                i + 1,
                mol,
                t,
                q,
                pos[0],
                pos[1],
                pos[2]
            ));
        }
        s
    }
    /// Generate the `Atoms` section for `atom_style molecular`.
    ///
    /// Format: `id mol-id type x y z`
    pub fn write_atoms_molecular(positions: &[[f64; 3]], types: &[u32], mol_ids: &[u32]) -> String {
        let mut s = String::from("\nAtoms  # molecular\n\n");
        for (i, ((&pos, &t), &mol)) in positions
            .iter()
            .zip(types.iter())
            .zip(mol_ids.iter())
            .enumerate()
        {
            s.push_str(&format!(
                "{} {} {} {} {} {}\n",
                i + 1,
                mol,
                t,
                pos[0],
                pos[1],
                pos[2]
            ));
        }
        s
    }
}
/// Reader for the LAMMPS data-file format.
#[allow(dead_code)]
pub struct LammpsDataReader {
    pub(super) positions: Vec<[f64; 3]>,
    pub(super) masses: Vec<f64>,
    pub(super) types: Vec<u32>,
    pub(super) box_lo: [f64; 3],
    pub(super) box_hi: [f64; 3],
}
#[allow(dead_code)]
impl LammpsDataReader {
    /// Parse a LAMMPS data file from a string.
    pub fn from_str(data: &str) -> Result<Self> {
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut masses_map: Vec<(usize, f64)> = Vec::new();
        let mut types: Vec<u32> = Vec::new();
        let mut box_lo = [0.0_f64; 3];
        let mut box_hi = [10.0_f64; 3];
        let mut in_atoms = false;
        let mut in_masses = false;
        for raw in data.lines() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line.starts_with("Atoms") {
                in_atoms = true;
                in_masses = false;
                continue;
            }
            if line.starts_with("Masses") {
                in_masses = true;
                in_atoms = false;
                continue;
            }
            if line.starts_with("Bonds") || line.starts_with("Velocities") {
                in_atoms = false;
                in_masses = false;
                continue;
            }
            if line.contains("xlo xhi") {
                let p = line.split_whitespace().collect::<Vec<_>>();
                box_lo[0] = p[0]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                box_hi[0] = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                continue;
            }
            if line.contains("ylo yhi") {
                let p = line.split_whitespace().collect::<Vec<_>>();
                box_lo[1] = p[0]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                box_hi[1] = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                continue;
            }
            if line.contains("zlo zhi") {
                let p = line.split_whitespace().collect::<Vec<_>>();
                box_lo[2] = p[0]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                box_hi[2] = p[1]
                    .parse::<f64>()
                    .map_err(|e| Error::Parse(e.to_string()))?;
                continue;
            }
            if in_masses {
                let p = line.split_whitespace().collect::<Vec<_>>();
                if p.len() >= 2 {
                    let tid = p[0]
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let m = p[1]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    masses_map.push((tid, m));
                }
                continue;
            }
            if in_atoms {
                let p = line.split_whitespace().collect::<Vec<_>>();
                if p.len() >= 5 {
                    let t = p[1]
                        .parse::<u32>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let x = p[2]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let y = p[3]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let z = p[4]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    positions.push([x, y, z]);
                    types.push(t);
                }
                continue;
            }
        }
        let max_type = types.iter().copied().max().unwrap_or(1) as usize;
        let mut type_mass = vec![1.0_f64; max_type + 1];
        for (tid, m) in masses_map {
            if tid <= max_type {
                type_mass[tid] = m;
            }
        }
        let masses: Vec<f64> = types.iter().map(|&t| type_mass[t as usize]).collect();
        Ok(Self {
            positions,
            masses,
            types,
            box_lo,
            box_hi,
        })
    }
    /// Return parsed positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }
    /// Return per-atom masses.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }
    /// Return per-atom type IDs (1-indexed).
    pub fn types(&self) -> &[u32] {
        &self.types
    }
    /// Return the simulation box bounds as `(lo, hi)`.
    pub fn box_bounds(&self) -> ([f64; 3], [f64; 3]) {
        (self.box_lo, self.box_hi)
    }
}
/// Represents a single atom in a LAMMPS dump frame.
#[derive(Debug, Clone)]
pub struct LammpsAtom {
    /// Atom ID.
    pub id: usize,
    /// Atom type ID.
    pub type_id: usize,
    /// 3D position.
    pub position: Vec3,
    /// 3D velocity.
    pub velocity: Vec3,
}
/// Generator for LAMMPS fix commands.
#[allow(dead_code)]
pub struct LammpsFix;
#[allow(dead_code)]
impl LammpsFix {
    /// Generate an NVE fix.
    pub fn nve(fix_id: &str, group: &str) -> String {
        format!("fix {fix_id} {group} nve\n")
    }
    /// Generate an NVT (Nose-Hoover thermostat) fix.
    pub fn nvt(fix_id: &str, group: &str, t_start: f64, t_end: f64, t_damp: f64) -> String {
        format!("fix {fix_id} {group} nvt temp {t_start} {t_end} {t_damp}\n")
    }
    /// Generate a Langevin thermostat fix.
    pub fn langevin(
        fix_id: &str,
        group: &str,
        t_start: f64,
        t_end: f64,
        damp: f64,
        seed: u64,
    ) -> String {
        format!("fix {fix_id} {group} langevin {t_start} {t_end} {damp} {seed}\n")
    }
    /// Generate a wall/reflect fix.
    pub fn wall_reflect(fix_id: &str, group: &str, face: &str, coord: f64) -> String {
        format!("fix {fix_id} {group} wall/reflect {face} EDGE\n")
            .replace("EDGE", &format!("{coord}"))
    }
    /// Generate a `fix setforce` command.
    pub fn setforce(fix_id: &str, group: &str, fx: f64, fy: f64, fz: f64) -> String {
        format!("fix {fix_id} {group} setforce {fx} {fy} {fz}\n")
    }
    /// Generate a `fix rigid` command.
    pub fn rigid(fix_id: &str, group: &str, style: &str) -> String {
        format!("fix {fix_id} {group} rigid {style}\n")
    }
}
/// Lennard-Jones pair interaction parameters for one type pair.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsLJPair {
    /// Type index I (1-based).
    pub type_i: u32,
    /// Type index J (1-based).
    pub type_j: u32,
    /// Well depth epsilon (energy units).
    pub epsilon: f64,
    /// Diameter sigma (distance units).
    pub sigma: f64,
    /// Cutoff distance. If `None`, uses the global cutoff.
    pub cutoff: Option<f64>,
}
#[allow(dead_code)]
impl LammpsLJPair {
    /// Compute the 12-6 LJ potential energy at separation `r`.
    ///
    /// `U(r) = 4ε[(σ/r)^12 − (σ/r)^6]`
    pub fn potential_energy(&self, r: f64) -> f64 {
        let sr = self.sigma / r;
        let sr6 = sr * sr * sr * sr * sr * sr;
        let sr12 = sr6 * sr6;
        4.0 * self.epsilon * (sr12 - sr6)
    }
    /// Compute the magnitude of the 12-6 LJ force at separation `r`.
    ///
    /// `F(r) = 24ε/r [2(σ/r)^12 − (σ/r)^6]`
    pub fn force_magnitude(&self, r: f64) -> f64 {
        let sr = self.sigma / r;
        let sr6 = sr * sr * sr * sr * sr * sr;
        let sr12 = sr6 * sr6;
        24.0 * self.epsilon / r * (2.0 * sr12 - sr6)
    }
    /// Return the distance where the LJ potential is at its minimum (r_min = 2^(1/6) σ).
    pub fn r_min(&self) -> f64 {
        2.0_f64.powf(1.0 / 6.0) * self.sigma
    }
    /// Write LAMMPS `pair_coeff` line for this pair.
    pub fn to_pair_coeff_line(&self) -> String {
        if let Some(rc) = self.cutoff {
            format!(
                "pair_coeff {} {} {} {} {}\n",
                self.type_i, self.type_j, self.epsilon, self.sigma, rc
            )
        } else {
            format!(
                "pair_coeff {} {} {} {}\n",
                self.type_i, self.type_j, self.epsilon, self.sigma
            )
        }
    }
}
/// A single thermo output record from a LAMMPS log file.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LammpsThermoRecord {
    /// Step number.
    pub step: u64,
    /// Temperature (K), if present.
    pub temp: Option<f64>,
    /// Total energy (energy units), if present.
    pub etotal: Option<f64>,
    /// Potential energy, if present.
    pub pe: Option<f64>,
    /// Kinetic energy, if present.
    pub ke: Option<f64>,
    /// Pressure (pressure units), if present.
    pub press: Option<f64>,
    /// Volume, if present.
    pub vol: Option<f64>,
}
/// Builder for the coefficient sections of a LAMMPS data file.
///
/// Generates the `Masses`, `Pair Coeffs`, and `Bond Coeffs` sections.
#[allow(dead_code)]
pub struct LammpsDataSectionBuilder {
    /// (type_id, mass)
    pub masses: Vec<(usize, f64)>,
    /// (type_i, type_j, epsilon, sigma)
    pub pair_coeffs: Vec<(usize, usize, f64, f64)>,
    /// (bond_type, k, r0)
    pub bond_coeffs: Vec<(usize, f64, f64)>,
    /// (angle_type, k, theta0_deg)
    pub angle_coeffs: Vec<(usize, f64, f64)>,
}
impl LammpsDataSectionBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self {
            masses: Vec::new(),
            pair_coeffs: Vec::new(),
            bond_coeffs: Vec::new(),
            angle_coeffs: Vec::new(),
        }
    }
    /// Add a mass entry.
    pub fn add_mass(&mut self, type_id: usize, mass: f64) {
        self.masses.push((type_id, mass));
    }
    /// Add a Lennard-Jones pair coefficient.
    pub fn add_pair_coeff(&mut self, type_i: usize, type_j: usize, epsilon: f64, sigma: f64) {
        self.pair_coeffs.push((type_i, type_j, epsilon, sigma));
    }
    /// Add a harmonic bond coefficient.
    pub fn add_bond_coeff(&mut self, bond_type: usize, k: f64, r0: f64) {
        self.bond_coeffs.push((bond_type, k, r0));
    }
    /// Add a harmonic angle coefficient.
    pub fn add_angle_coeff(&mut self, angle_type: usize, k: f64, theta0_deg: f64) {
        self.angle_coeffs.push((angle_type, k, theta0_deg));
    }
    /// Generate the `Masses` section text.
    pub fn masses_section(&self) -> String {
        let mut s = "Masses\n\n".to_string();
        for &(id, m) in &self.masses {
            s.push_str(&format!("{} {:.6}\n", id, m));
        }
        s
    }
    /// Generate the `Pair Coeffs` section text.
    pub fn pair_coeffs_section(&self) -> String {
        let mut s = "Pair Coeffs\n\n".to_string();
        for &(i, j, eps, sig) in &self.pair_coeffs {
            s.push_str(&format!("{} {} {:.6} {:.6}\n", i, j, eps, sig));
        }
        s
    }
    /// Generate the `Bond Coeffs` section text.
    pub fn bond_coeffs_section(&self) -> String {
        let mut s = "Bond Coeffs\n\n".to_string();
        for &(bt, k, r0) in &self.bond_coeffs {
            s.push_str(&format!("{} {:.6} {:.6}\n", bt, k, r0));
        }
        s
    }
    /// Generate the `Angle Coeffs` section text.
    pub fn angle_coeffs_section(&self) -> String {
        let mut s = "Angle Coeffs\n\n".to_string();
        for &(at, k, th) in &self.angle_coeffs {
            s.push_str(&format!("{} {:.6} {:.6}\n", at, k, th));
        }
        s
    }
    /// Combine all sections into one string.
    pub fn build(&self) -> String {
        let mut out = String::new();
        if !self.masses.is_empty() {
            out.push_str(&self.masses_section());
            out.push('\n');
        }
        if !self.pair_coeffs.is_empty() {
            out.push_str(&self.pair_coeffs_section());
            out.push('\n');
        }
        if !self.bond_coeffs.is_empty() {
            out.push_str(&self.bond_coeffs_section());
            out.push('\n');
        }
        if !self.angle_coeffs.is_empty() {
            out.push_str(&self.angle_coeffs_section());
            out.push('\n');
        }
        out
    }
}
/// Harmonic bond potential parameters.
///
/// `U(r) = K (r − r₀)²`
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsHarmonicBond {
    /// Bond type index (1-based).
    pub bond_type: u32,
    /// Spring constant K (energy/distance²).
    pub k: f64,
    /// Equilibrium bond length r₀.
    pub r0: f64,
}
#[allow(dead_code)]
impl LammpsHarmonicBond {
    /// Compute bond potential energy.
    pub fn potential_energy(&self, r: f64) -> f64 {
        self.k * (r - self.r0) * (r - self.r0)
    }
    /// Compute bond force magnitude (restoring).
    pub fn force_magnitude(&self, r: f64) -> f64 {
        2.0 * self.k * (r - self.r0)
    }
    /// Generate LAMMPS `bond_coeff` line.
    pub fn to_bond_coeff_line(&self) -> String {
        format!("bond_coeff {} {} {}\n", self.bond_type, self.k, self.r0)
    }
}
/// Builder for a LAMMPS input script `run` / `timestep` / `fix` block.
#[allow(dead_code)]
pub struct LammpsRunBlock {
    /// Integration timestep in fs.
    pub timestep: f64,
    /// Number of MD steps.
    pub run_steps: u64,
    /// Thermo output frequency.
    pub thermo_freq: u64,
    /// Dump output frequency.
    pub dump_freq: u64,
    /// Ensemble: "nve", "nvt", "npt".
    pub ensemble: String,
    /// Temperature target (for NVT/NPT).
    pub temp_target: f64,
    /// Pressure target in bar (for NPT).
    pub pressure_target: Option<f64>,
}
impl LammpsRunBlock {
    /// Create a new NVE run block.
    pub fn nve(timestep: f64, run_steps: u64) -> Self {
        Self {
            timestep,
            run_steps,
            thermo_freq: 100,
            dump_freq: 1000,
            ensemble: "nve".to_string(),
            temp_target: 300.0,
            pressure_target: None,
        }
    }
    /// Create a new NVT run block.
    pub fn nvt(timestep: f64, run_steps: u64, temp: f64) -> Self {
        Self {
            timestep,
            run_steps,
            thermo_freq: 100,
            dump_freq: 1000,
            ensemble: "nvt".to_string(),
            temp_target: temp,
            pressure_target: None,
        }
    }
    /// Create a new NPT run block.
    pub fn npt(timestep: f64, run_steps: u64, temp: f64, pressure: f64) -> Self {
        Self {
            timestep,
            run_steps,
            thermo_freq: 100,
            dump_freq: 1000,
            ensemble: "npt".to_string(),
            temp_target: temp,
            pressure_target: Some(pressure),
        }
    }
    /// Set thermo output frequency.
    pub fn thermo_every(mut self, freq: u64) -> Self {
        self.thermo_freq = freq;
        self
    }
    /// Set dump frequency.
    pub fn dump_every(mut self, freq: u64) -> Self {
        self.dump_freq = freq;
        self
    }
    /// Total simulated time in fs.
    pub fn total_time_fs(&self) -> f64 {
        self.timestep * self.run_steps as f64
    }
    /// Generate the LAMMPS input script block.
    pub fn to_script(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("timestep {:.6}\n", self.timestep));
        s.push_str(&format!("thermo {}\n", self.thermo_freq));
        s.push_str("thermo_style custom step temp pe ke etotal press vol\n");
        s.push_str(&format!(
            "dump 1 all custom {} dump.lammpstrj id type x y z vx vy vz\n",
            self.dump_freq
        ));
        match self.ensemble.as_str() {
            "nve" => {
                s.push_str("fix 1 all nve\n");
            }
            "nvt" => {
                s.push_str(&format!(
                    "fix 1 all nvt temp {t:.1} {t:.1} $(100.0*dt)\n",
                    t = self.temp_target
                ));
            }
            "npt" => {
                let p = self.pressure_target.unwrap_or(1.0);
                s.push_str(&format!(
                    "fix 1 all npt temp {t:.1} {t:.1} $(100.0*dt) iso {p:.1} {p:.1} $(1000.0*dt)\n",
                    t = self.temp_target,
                    p = p
                ));
            }
            _ => {}
        }
        s.push_str(&format!("run {}\n", self.run_steps));
        s
    }
}
/// Reader for LAMMPS custom dump format.
pub struct LammpsDumpReader;
impl LammpsDumpReader {
    /// Read a single frame from a LAMMPS dump stream.
    ///
    /// Returns `(timestep, atoms)`. The stream must start at the beginning of a frame.
    pub fn read_frame<R: Read>(reader: R) -> Result<(u64, Vec<LammpsAtom>)> {
        let buf = BufReader::new(reader);
        let mut lines = buf.lines();
        let mut timestep: u64 = 0;
        let mut n_atoms: usize = 0;
        let mut atoms: Vec<LammpsAtom> = Vec::new();
        let mut reading_atoms = false;
        let mut atom_count = 0;
        while let Some(line) = lines.next() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed == "ITEM: TIMESTEP" {
                if let Some(ts_line) = lines.next() {
                    timestep = ts_line?
                        .trim()
                        .parse::<u64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                }
            } else if trimmed == "ITEM: NUMBER OF ATOMS" {
                if let Some(n_line) = lines.next() {
                    n_atoms = n_line?
                        .trim()
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    atoms.reserve(n_atoms);
                }
            } else if trimmed.starts_with("ITEM: BOX BOUNDS") {
                for _ in 0..3 {
                    lines.next();
                }
            } else if trimmed.starts_with("ITEM: ATOMS") {
                reading_atoms = true;
            } else if reading_atoms && atom_count < n_atoms {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 8 {
                    let id = parts[0]
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let type_id = parts[1]
                        .parse::<usize>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let x = parts[2]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let y = parts[3]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let z = parts[4]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let vx = parts[5]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let vy = parts[6]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    let vz = parts[7]
                        .parse::<f64>()
                        .map_err(|e| Error::Parse(e.to_string()))?;
                    atoms.push(LammpsAtom {
                        id,
                        type_id,
                        position: Vec3::new(x, y, z),
                        velocity: Vec3::new(vx, vy, vz),
                    });
                    atom_count += 1;
                }
            }
        }
        Ok((timestep, atoms))
    }
}
/// Generator for LAMMPS pair style commands.
#[allow(dead_code)]
pub struct LammpsPairStyle;
#[allow(dead_code)]
impl LammpsPairStyle {
    /// Generate a `pair_style lj/cut` command.
    pub fn lj_cut(cutoff: f64) -> String {
        format!("pair_style lj/cut {cutoff}\n")
    }
    /// Generate a `pair_coeff` command for LJ parameters.
    pub fn pair_coeff_lj(type_i: u32, type_j: u32, epsilon: f64, sigma: f64) -> String {
        format!("pair_coeff {type_i} {type_j} {epsilon} {sigma}\n")
    }
    /// Generate a `pair_style lj/cut/coul/long` command (for charged systems).
    pub fn lj_cut_coul_long(lj_cutoff: f64, coul_cutoff: f64) -> String {
        format!("pair_style lj/cut/coul/long {lj_cutoff} {coul_cutoff}\n")
    }
    /// Generate a `pair_style hybrid` command.
    pub fn hybrid(styles: &[&str]) -> String {
        format!("pair_style hybrid {}\n", styles.join(" "))
    }
    /// Generate a `pair_modify` command.
    pub fn pair_modify(options: &str) -> String {
        format!("pair_modify {options}\n")
    }
}
/// Atom style for a LAMMPS data file.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LammpsAtomStyle {
    /// `atomic` — id type x y z
    Atomic,
    /// `full` — id mol-id type charge x y z
    Full,
    /// `molecular` — id mol-id type x y z
    Molecular,
    /// `charge` — id type charge x y z
    Charge,
}
/// Parse LAMMPS thermo output from a log file string.
///
/// Looks for lines starting with a header of whitespace-separated column names,
/// then reads subsequent numeric rows.
#[allow(dead_code)]
pub struct LammpsLogParser;
#[allow(dead_code)]
impl LammpsLogParser {
    /// Parse thermo records from a LAMMPS log file string.
    ///
    /// Returns a vector of thermo records. Columns are matched by name (case-insensitive).
    /// Unknown columns are ignored.
    pub fn parse_thermo(log_text: &str) -> Vec<LammpsThermoRecord> {
        let mut records = Vec::new();
        let mut headers: Vec<String> = Vec::new();
        let mut in_thermo = false;
        for line in log_text.lines() {
            let line = line.trim();
            if line.starts_with("Step") || line.starts_with("step") {
                headers = line.split_whitespace().map(|s| s.to_lowercase()).collect();
                in_thermo = true;
                continue;
            }
            if line.starts_with("Loop") || line.starts_with("WARNING") || line.is_empty() {
                if in_thermo && !records.is_empty() {
                    in_thermo = false;
                }
                continue;
            }
            if in_thermo && !headers.is_empty() {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() < headers.len() {
                    continue;
                }
                let vals: Vec<Option<f64>> = parts.iter().map(|p| p.parse::<f64>().ok()).collect();
                if vals[0].is_none() {
                    in_thermo = false;
                    continue;
                }
                let get = |name: &str| -> Option<f64> {
                    headers
                        .iter()
                        .position(|h| h == name)
                        .and_then(|i| vals.get(i).and_then(|v| *v))
                };
                let step = match get("step").map(|v| v as u64) {
                    Some(s) => s,
                    None => continue,
                };
                records.push(LammpsThermoRecord {
                    step,
                    temp: get("temp"),
                    etotal: get("etotal"),
                    pe: get("pe"),
                    ke: get("ke"),
                    press: get("press"),
                    vol: get("vol"),
                });
            }
        }
        records
    }
    /// Extract the total simulation time from a log file.
    ///
    /// Looks for a line of the form `Total wall time: HH:MM:SS`.
    /// Returns `None` if not found.
    pub fn parse_wall_time(log_text: &str) -> Option<String> {
        for line in log_text.lines() {
            let t = line.trim();
            if t.starts_with("Total wall time:") {
                return Some(t.to_string());
            }
        }
        None
    }
    /// Count the number of MD steps in a LAMMPS log file.
    ///
    /// Looks for lines matching `run N` and sums up N.
    pub fn count_total_steps(log_text: &str) -> u64 {
        let mut total = 0u64;
        for line in log_text.lines() {
            let t = line.trim().to_lowercase();
            if t.starts_with("run ") {
                let parts: Vec<&str> = t.split_whitespace().collect();
                if parts.len() >= 2
                    && let Ok(n) = parts[1].parse::<u64>()
                {
                    total += n;
                }
            }
        }
        total
    }
}
/// Reader for the binary dump format written by [`LammpsBinaryDumpWriter`].
pub struct LammpsBinaryDumpReader;
impl LammpsBinaryDumpReader {
    /// Read one binary frame from the byte slice starting at `offset`.
    ///
    /// Returns `(timestep, box_bounds, atoms, bytes_consumed)` on success.
    #[allow(clippy::type_complexity)]
    pub fn read_frame(
        data: &[u8],
        offset: usize,
    ) -> Result<(i64, [[f64; 2]; 3], Vec<LammpsAtom>, usize)> {
        let mut pos = offset;
        let timestep = Self::read_i64(data, &mut pos)?;
        let n_atoms = Self::read_i64(data, &mut pos)? as usize;
        let mut box_bounds = [[0.0_f64; 2]; 3];
        for b in &mut box_bounds {
            b[0] = Self::read_f64(data, &mut pos)?;
            b[1] = Self::read_f64(data, &mut pos)?;
        }
        let mut atoms = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let id = Self::read_i32(data, &mut pos)? as usize;
            let type_id = Self::read_i32(data, &mut pos)? as usize;
            let x = Self::read_f64(data, &mut pos)?;
            let y = Self::read_f64(data, &mut pos)?;
            let z = Self::read_f64(data, &mut pos)?;
            let vx = Self::read_f64(data, &mut pos)?;
            let vy = Self::read_f64(data, &mut pos)?;
            let vz = Self::read_f64(data, &mut pos)?;
            atoms.push(LammpsAtom {
                id,
                type_id,
                position: Vec3::new(x, y, z),
                velocity: Vec3::new(vx, vy, vz),
            });
        }
        Ok((timestep, box_bounds, atoms, pos - offset))
    }
    fn read_i32(data: &[u8], pos: &mut usize) -> Result<i32> {
        if *pos + 4 > data.len() {
            return Err(Error::Parse("binary dump: unexpected EOF (i32)".into()));
        }
        let v = i32::from_le_bytes([data[*pos], data[*pos + 1], data[*pos + 2], data[*pos + 3]]);
        *pos += 4;
        Ok(v)
    }
    fn read_i64(data: &[u8], pos: &mut usize) -> Result<i64> {
        if *pos + 8 > data.len() {
            return Err(Error::Parse("binary dump: unexpected EOF (i64)".into()));
        }
        let b = &data[*pos..*pos + 8];
        let v = i64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
        *pos += 8;
        Ok(v)
    }
    fn read_f64(data: &[u8], pos: &mut usize) -> Result<f64> {
        if *pos + 8 > data.len() {
            return Err(Error::Parse("binary dump: unexpected EOF (f64)".into()));
        }
        let b = &data[*pos..*pos + 8];
        let v = f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
        *pos += 8;
        Ok(v)
    }
}
/// A table of LJ pair parameters for a simulation.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct LammpsLJTable {
    /// All pair interactions.
    pub pairs: Vec<LammpsLJPair>,
    /// Global LJ cutoff (distance units).
    pub global_cutoff: f64,
}
#[allow(dead_code)]
impl LammpsLJTable {
    /// Create an empty table with a given global cutoff.
    pub fn new(cutoff: f64) -> Self {
        Self {
            pairs: Vec::new(),
            global_cutoff: cutoff,
        }
    }
    /// Add a pair interaction.
    pub fn add_pair(&mut self, pair: LammpsLJPair) {
        self.pairs.push(pair);
    }
    /// Look up parameters for types (i, j). Symmetric lookup.
    pub fn get(&self, type_i: u32, type_j: u32) -> Option<&LammpsLJPair> {
        self.pairs.iter().find(|p| {
            (p.type_i == type_i && p.type_j == type_j) || (p.type_i == type_j && p.type_j == type_i)
        })
    }
    /// Generate `pair_style lj/cut` + all `pair_coeff` lines.
    pub fn to_lammps_input(&self) -> String {
        let mut s = LammpsPairStyle::lj_cut(self.global_cutoff);
        for pair in &self.pairs {
            s.push_str(&pair.to_pair_coeff_line());
        }
        s
    }
    /// Apply Lorentz-Berthelot mixing rules to fill in cross-interactions.
    ///
    /// For types `(i, j)` where `i != j`, if no cross-term exists,
    /// create one from the like-pairs using:
    /// - `ε_ij = sqrt(ε_ii * ε_jj)`
    /// - `σ_ij = (σ_ii + σ_jj) / 2`
    pub fn apply_lorentz_berthelot_mixing(&mut self) {
        let like: Vec<(u32, f64, f64)> = self
            .pairs
            .iter()
            .filter(|p| p.type_i == p.type_j)
            .map(|p| (p.type_i, p.epsilon, p.sigma))
            .collect();
        let mut new_pairs = Vec::new();
        for i in 0..like.len() {
            for j in (i + 1)..like.len() {
                let (ti, eps_i, sig_i) = like[i];
                let (tj, eps_j, sig_j) = like[j];
                if self.get(ti, tj).is_none() {
                    new_pairs.push(LammpsLJPair {
                        type_i: ti,
                        type_j: tj,
                        epsilon: (eps_i * eps_j).sqrt(),
                        sigma: (sig_i + sig_j) / 2.0,
                        cutoff: None,
                    });
                }
            }
        }
        self.pairs.extend(new_pairs);
    }
}
/// Harmonic angle potential parameters.
///
/// `U(θ) = K (θ − θ₀)²`
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct LammpsHarmonicAngle {
    /// Angle type index (1-based).
    pub angle_type: u32,
    /// Spring constant K (energy/rad²).
    pub k: f64,
    /// Equilibrium angle θ₀ (degrees).
    pub theta0_deg: f64,
}
#[allow(dead_code)]
impl LammpsHarmonicAngle {
    /// Compute angle potential energy given `theta` in degrees.
    pub fn potential_energy_deg(&self, theta_deg: f64) -> f64 {
        let dt = (theta_deg - self.theta0_deg).to_radians();
        self.k * dt * dt
    }
    /// Generate LAMMPS `angle_coeff` line.
    pub fn to_angle_coeff_line(&self) -> String {
        format!(
            "angle_coeff {} {} {}\n",
            self.angle_type, self.k, self.theta0_deg
        )
    }
}
/// Generator for LAMMPS compute commands.
#[allow(dead_code)]
pub struct LammpsCompute;
#[allow(dead_code)]
impl LammpsCompute {
    /// Generate a `compute` for temperature.
    pub fn temp(compute_id: &str, group: &str) -> String {
        format!("compute {compute_id} {group} temp\n")
    }
    /// Generate a `compute` for potential energy.
    pub fn pe(compute_id: &str, group: &str) -> String {
        format!("compute {compute_id} {group} pe\n")
    }
    /// Generate a `compute` for kinetic energy.
    pub fn ke(compute_id: &str, group: &str) -> String {
        format!("compute {compute_id} {group} ke\n")
    }
    /// Generate a `compute` for pressure.
    pub fn pressure(compute_id: &str, group: &str, temp_compute: &str) -> String {
        format!("compute {compute_id} {group} pressure {temp_compute}\n")
    }
    /// Generate a `compute` for RDF.
    pub fn rdf(compute_id: &str, group: &str, n_bins: usize, cutoff: f64) -> String {
        format!("compute {compute_id} {group} rdf {n_bins} cutoff {cutoff}\n")
    }
}
/// Writer for LAMMPS custom dump format.
pub struct LammpsDumpWriter;
impl LammpsDumpWriter {
    /// Write a single frame in LAMMPS custom dump format.
    ///
    /// `box_bounds` is `[[xlo, xhi\], [ylo, yhi], [zlo, zhi]]`.
    pub fn write_frame<W: Write>(
        writer: &mut W,
        timestep: u64,
        box_bounds: [[f64; 2]; 3],
        atoms: &[LammpsAtom],
    ) -> Result<()> {
        writeln!(writer, "ITEM: TIMESTEP")?;
        writeln!(writer, "{}", timestep)?;
        writeln!(writer, "ITEM: NUMBER OF ATOMS")?;
        writeln!(writer, "{}", atoms.len())?;
        writeln!(writer, "ITEM: BOX BOUNDS pp pp pp")?;
        for b in &box_bounds {
            writeln!(writer, "{} {}", b[0], b[1])?;
        }
        writeln!(writer, "ITEM: ATOMS id type x y z vx vy vz")?;
        for a in atoms {
            writeln!(
                writer,
                "{} {} {} {} {} {} {} {}",
                a.id,
                a.type_id,
                a.position.x,
                a.position.y,
                a.position.z,
                a.velocity.x,
                a.velocity.y,
                a.velocity.z,
            )?;
        }
        Ok(())
    }
}
/// Reader for LAMMPS restart files written by [`LammpsRestartWriter`].
pub struct LammpsRestartReader;
impl LammpsRestartReader {
    /// Parse a restart file from bytes.
    ///
    /// Returns `(timestep, masses, atoms)` on success.
    pub fn read(data: &[u8]) -> Result<(i64, Vec<f64>, Vec<LammpsAtom>)> {
        if data.len() < 4 || &data[..4] != b"LRST" {
            return Err(Error::Parse("restart: missing LRST magic".into()));
        }
        let mut pos = 4usize;
        let _version = read_i32_le(data, &mut pos)?;
        let timestep = read_i64_le(data, &mut pos)?;
        let n_atoms = read_i64_le(data, &mut pos)? as usize;
        let n_types = read_i32_le(data, &mut pos)? as usize;
        let mut masses = Vec::with_capacity(n_types);
        for _ in 0..n_types {
            masses.push(read_f64_le(data, &mut pos)?);
        }
        let mut atoms = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let id = read_i32_le(data, &mut pos)? as usize;
            let type_id = read_i32_le(data, &mut pos)? as usize;
            let x = read_f64_le(data, &mut pos)?;
            let y = read_f64_le(data, &mut pos)?;
            let z = read_f64_le(data, &mut pos)?;
            let vx = read_f64_le(data, &mut pos)?;
            let vy = read_f64_le(data, &mut pos)?;
            let vz = read_f64_le(data, &mut pos)?;
            atoms.push(LammpsAtom {
                id,
                type_id,
                position: Vec3::new(x, y, z),
                velocity: Vec3::new(vx, vy, vz),
            });
        }
        Ok((timestep, masses, atoms))
    }
}
