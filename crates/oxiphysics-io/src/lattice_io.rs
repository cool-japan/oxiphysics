// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lattice structure I/O for LBM and crystallography.
//!
//! Provides `LatticeGrid` for LBM fields and `CrystalLattice` for crystal
//! structures, together with binary, VTK, XYZ, POSCAR, and CIF output routines.

use std::fs::File;
use std::io::{Read, Write};

// ── Lattice type ──────────────────────────────────────────────────────────────

/// Lattice connectivity / velocity-set type.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum LatticeType {
    /// 2-D, 9-velocity lattice (LBM standard).
    D2Q9,
    /// 2-D, 5-velocity lattice.
    D2Q5,
    /// 3-D, 15-velocity lattice.
    D3Q15,
    /// 3-D, 19-velocity lattice (most common 3-D LBM).
    D3Q19,
    /// 3-D, 27-velocity lattice.
    D3Q27,
    /// Simple cubic crystal.
    SC,
    /// Face-centred cubic crystal.
    FCC,
    /// Body-centred cubic crystal.
    BCC,
    /// Hexagonal close-packed crystal.
    HCP,
    /// User-defined velocity set.
    Custom(Vec<[i32; 3]>),
}

impl LatticeType {
    /// Number of discrete velocities in this lattice.
    pub fn velocity_count(&self) -> usize {
        match self {
            LatticeType::D2Q9 => 9,
            LatticeType::D2Q5 => 5,
            LatticeType::D3Q15 => 15,
            LatticeType::D3Q19 => 19,
            LatticeType::D3Q27 => 27,
            LatticeType::SC => 6,
            LatticeType::FCC => 12,
            LatticeType::BCC => 8,
            LatticeType::HCP => 12,
            LatticeType::Custom(v) => v.len(),
        }
    }

    /// Numeric type identifier for binary serialisation.
    fn type_id(&self) -> u8 {
        match self {
            LatticeType::D2Q9 => 0,
            LatticeType::D2Q5 => 1,
            LatticeType::D3Q15 => 2,
            LatticeType::D3Q19 => 3,
            LatticeType::D3Q27 => 4,
            LatticeType::SC => 5,
            LatticeType::FCC => 6,
            LatticeType::BCC => 7,
            LatticeType::HCP => 8,
            LatticeType::Custom(_) => 9,
        }
    }

    /// Reconstruct from type id (Custom returns D3Q19 as fallback).
    fn from_type_id(id: u8) -> Self {
        match id {
            0 => LatticeType::D2Q9,
            1 => LatticeType::D2Q5,
            2 => LatticeType::D3Q15,
            3 => LatticeType::D3Q19,
            4 => LatticeType::D3Q27,
            5 => LatticeType::SC,
            6 => LatticeType::FCC,
            7 => LatticeType::BCC,
            8 => LatticeType::HCP,
            _ => LatticeType::D3Q19,
        }
    }
}

// ── LatticeGrid ───────────────────────────────────────────────────────────────

/// A 3-D labelled lattice field.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct LatticeGrid {
    /// Connectivity / velocity-set type.
    pub lattice_type: LatticeType,
    /// Grid size in x.
    pub nx: usize,
    /// Grid size in y.
    pub ny: usize,
    /// Grid size in z.
    pub nz: usize,
    /// Physical lattice spacing (m or dimensionless).
    pub dx: f64,
    /// Per-node data: `node_data[node_idx]` is a slice of length `n_components`.
    pub node_data: Vec<Vec<f64>>,
    /// Number of scalar components stored at each node.
    pub n_components: usize,
}

impl LatticeGrid {
    /// Create a new lattice grid initialised to zero.
    pub fn new(lt: LatticeType, nx: usize, ny: usize, nz: usize, dx: f64, n_comp: usize) -> Self {
        let nc = nx * ny * nz;
        Self {
            lattice_type: lt,
            nx,
            ny,
            nz,
            dx,
            node_data: vec![vec![0.0; n_comp]; nc],
            n_components: n_comp,
        }
    }

    /// Flat node index for node (i, j, k).
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.nx * self.ny + j * self.nx + i
    }

    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Set the data at node (i, j, k).
    pub fn set_node(&mut self, i: usize, j: usize, k: usize, data: &[f64]) {
        let idx = self.index(i, j, k);
        let nc = self.n_components.min(data.len());
        self.node_data[idx][..nc].copy_from_slice(&data[..nc]);
    }

    /// Get the data slice at node (i, j, k).
    pub fn get_node(&self, i: usize, j: usize, k: usize) -> &[f64] {
        let idx = self.index(i, j, k);
        &self.node_data[idx]
    }
}

// ── Binary I/O ─────────────────────────────────────────────────────────────────

/// Write a `LatticeGrid` to a custom binary file.
///
/// Format: `[magic: 4 bytes][type_id: 1 byte][nx,ny,nz: 3×u64 LE][dx: f64 LE]
/// [n_comp: u64 LE][data: n_nodes × n_comp × f64 LE]`
pub fn write_lattice_binary(grid: &LatticeGrid, path: &str) -> Result<(), String> {
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    // Magic
    f.write_all(b"LATT").map_err(|e| e.to_string())?;
    // Type id
    f.write_all(&[grid.lattice_type.type_id()])
        .map_err(|e| e.to_string())?;
    // Dimensions
    for dim in [grid.nx as u64, grid.ny as u64, grid.nz as u64] {
        f.write_all(&dim.to_le_bytes()).map_err(|e| e.to_string())?;
    }
    // dx
    f.write_all(&grid.dx.to_le_bytes())
        .map_err(|e| e.to_string())?;
    // n_components
    f.write_all(&(grid.n_components as u64).to_le_bytes())
        .map_err(|e| e.to_string())?;
    // Data
    for node in &grid.node_data {
        for &val in node {
            f.write_all(&val.to_le_bytes()).map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

/// Read a `LatticeGrid` from a custom binary file.
pub fn read_lattice_binary(path: &str) -> Result<LatticeGrid, String> {
    let mut f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let mut buf4 = [0u8; 4];
    f.read_exact(&mut buf4).map_err(|e| e.to_string())?;
    if &buf4 != b"LATT" {
        return Err("invalid magic bytes".into());
    }
    let mut type_id_buf = [0u8; 1];
    f.read_exact(&mut type_id_buf).map_err(|e| e.to_string())?;
    let lt = LatticeType::from_type_id(type_id_buf[0]);

    let read_u64 = |f: &mut File| -> Result<u64, String> {
        let mut b = [0u8; 8];
        f.read_exact(&mut b).map_err(|e| e.to_string())?;
        Ok(u64::from_le_bytes(b))
    };
    let nx = read_u64(&mut f)? as usize;
    let ny = read_u64(&mut f)? as usize;
    let nz = read_u64(&mut f)? as usize;

    let mut dx_buf = [0u8; 8];
    f.read_exact(&mut dx_buf).map_err(|e| e.to_string())?;
    let dx = f64::from_le_bytes(dx_buf);

    let n_comp = read_u64(&mut f)? as usize;
    let nc = nx * ny * nz;
    let mut node_data: Vec<Vec<f64>> = vec![vec![0.0; n_comp]; nc];
    let mut val_buf = [0u8; 8];
    for node in node_data.iter_mut() {
        for comp in node.iter_mut() {
            f.read_exact(&mut val_buf).map_err(|e| e.to_string())?;
            *comp = f64::from_le_bytes(val_buf);
        }
    }
    Ok(LatticeGrid {
        lattice_type: lt,
        nx,
        ny,
        nz,
        dx,
        node_data,
        n_components: n_comp,
    })
}

// ── VTK export ────────────────────────────────────────────────────────────────

/// Write the first component of a `LatticeGrid` field as a VTK structured grid.
pub fn write_lattice_vtk(grid: &LatticeGrid, path: &str, field_name: &str) -> Result<(), String> {
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    writeln!(f, "# vtk DataFile Version 3.0").map_err(|e| e.to_string())?;
    writeln!(f, "LatticeGrid").map_err(|e| e.to_string())?;
    writeln!(f, "ASCII").map_err(|e| e.to_string())?;
    writeln!(f, "DATASET STRUCTURED_POINTS").map_err(|e| e.to_string())?;
    writeln!(f, "DIMENSIONS {} {} {}", grid.nx, grid.ny, grid.nz).map_err(|e| e.to_string())?;
    writeln!(f, "ORIGIN 0 0 0").map_err(|e| e.to_string())?;
    writeln!(f, "SPACING {} {} {}", grid.dx, grid.dx, grid.dx).map_err(|e| e.to_string())?;
    writeln!(f, "POINT_DATA {}", grid.node_count()).map_err(|e| e.to_string())?;
    writeln!(f, "SCALARS {} double 1", field_name).map_err(|e| e.to_string())?;
    writeln!(f, "LOOKUP_TABLE default").map_err(|e| e.to_string())?;
    for node in &grid.node_data {
        let val = node.first().copied().unwrap_or(0.0);
        writeln!(f, "{val}").map_err(|e| e.to_string())?;
    }
    Ok(())
}

// ── XYZ export ────────────────────────────────────────────────────────────────

/// Write lattice node positions to an XYZ file for visualisation.
pub fn lattice_to_xyz(grid: &LatticeGrid, path: &str) -> Result<(), String> {
    let nc = grid.node_count();
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    writeln!(f, "{nc}").map_err(|e| e.to_string())?;
    writeln!(f, "LatticeGrid XYZ export").map_err(|e| e.to_string())?;
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let x = i as f64 * grid.dx;
                let y = j as f64 * grid.dx;
                let z = k as f64 * grid.dx;
                writeln!(f, "Latt {x:.6} {y:.6} {z:.6}").map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(())
}

// ── Crystal lattice ───────────────────────────────────────────────────────────

/// A crystal structure defined by basis atoms and lattice vectors.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CrystalLattice {
    /// Fractional or Cartesian coordinates of basis atoms.
    pub basis_atoms: Vec<[f64; 3]>,
    /// Lattice vectors as rows: `[[a1x, a1y, a1z\], [a2x, …], [a3x, …]]`.
    pub lattice_vectors: [[f64; 3]; 3],
    /// Element symbol for each basis atom.
    pub atom_types: Vec<String>,
}

impl CrystalLattice {
    /// Create an FCC lattice with lattice constant `a`.
    ///
    /// FCC has 4 basis atoms in the conventional cell.
    pub fn fcc(a: f64) -> Self {
        Self {
            basis_atoms: vec![
                [0.0, 0.0, 0.0],
                [0.5 * a, 0.5 * a, 0.0],
                [0.5 * a, 0.0, 0.5 * a],
                [0.0, 0.5 * a, 0.5 * a],
            ],
            lattice_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
            atom_types: vec!["X".into(); 4],
        }
    }

    /// Create a BCC lattice with lattice constant `a`.
    ///
    /// BCC has 2 basis atoms in the conventional cell.
    pub fn bcc(a: f64) -> Self {
        Self {
            basis_atoms: vec![[0.0, 0.0, 0.0], [0.5 * a, 0.5 * a, 0.5 * a]],
            lattice_vectors: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
            atom_types: vec!["X".into(); 2],
        }
    }

    /// Create an HCP lattice with lattice constants `a` and `c`.
    ///
    /// HCP has 2 atoms per primitive cell.
    pub fn hcp(a: f64, c: f64) -> Self {
        Self {
            basis_atoms: vec![
                [0.0, 0.0, 0.0],
                [0.5 * a, a * (3.0_f64).sqrt() / 6.0, 0.5 * c],
            ],
            lattice_vectors: [
                [a, 0.0, 0.0],
                [0.5 * a, a * (3.0_f64).sqrt() / 2.0, 0.0],
                [0.0, 0.0, c],
            ],
            atom_types: vec!["X".into(); 2],
        }
    }

    /// Number of basis atoms.
    pub fn atom_count(&self) -> usize {
        self.basis_atoms.len()
    }
}

/// Build a supercell with `n1 × n2 × n3` repetitions.
pub fn supercell(crystal: &CrystalLattice, n1: usize, n2: usize, n3: usize) -> CrystalLattice {
    let lv = crystal.lattice_vectors;
    let mut new_atoms: Vec<[f64; 3]> = Vec::new();
    let mut new_types: Vec<String> = Vec::new();
    for k in 0..n3 {
        for j in 0..n2 {
            for i in 0..n1 {
                let shift = [
                    i as f64 * lv[0][0] + j as f64 * lv[1][0] + k as f64 * lv[2][0],
                    i as f64 * lv[0][1] + j as f64 * lv[1][1] + k as f64 * lv[2][1],
                    i as f64 * lv[0][2] + j as f64 * lv[1][2] + k as f64 * lv[2][2],
                ];
                for (atom, atype) in crystal.basis_atoms.iter().zip(crystal.atom_types.iter()) {
                    new_atoms.push([atom[0] + shift[0], atom[1] + shift[1], atom[2] + shift[2]]);
                    new_types.push(atype.clone());
                }
            }
        }
    }
    let new_lv = [
        [
            lv[0][0] * n1 as f64,
            lv[0][1] * n1 as f64,
            lv[0][2] * n1 as f64,
        ],
        [
            lv[1][0] * n2 as f64,
            lv[1][1] * n2 as f64,
            lv[1][2] * n2 as f64,
        ],
        [
            lv[2][0] * n3 as f64,
            lv[2][1] * n3 as f64,
            lv[2][2] * n3 as f64,
        ],
    ];
    CrystalLattice {
        basis_atoms: new_atoms,
        lattice_vectors: new_lv,
        atom_types: new_types,
    }
}

// ── POSCAR / CIF output ───────────────────────────────────────────────────────

/// Write crystal structure to a VASP POSCAR file with `n × n × n` supercell.
pub fn write_poscar(crystal: &CrystalLattice, n: usize, path: &str) -> Result<(), String> {
    let sc = supercell(crystal, n, n, n);
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    writeln!(f, "Generated by OxiPhysics lattice_io").map_err(|e| e.to_string())?;
    writeln!(f, "1.0").map_err(|e| e.to_string())?;
    for row in &sc.lattice_vectors {
        writeln!(f, "  {:.8}  {:.8}  {:.8}", row[0], row[1], row[2]).map_err(|e| e.to_string())?;
    }
    // Unique element types
    let types: Vec<&str> = {
        let mut seen: Vec<&str> = Vec::new();
        for t in &sc.atom_types {
            if !seen.contains(&t.as_str()) {
                seen.push(t.as_str());
            }
        }
        seen
    };
    let type_counts: Vec<usize> = types
        .iter()
        .map(|t| sc.atom_types.iter().filter(|a| a.as_str() == *t).count())
        .collect();
    writeln!(f, "{}", types.join(" ")).map_err(|e| e.to_string())?;
    writeln!(
        f,
        "{}",
        type_counts
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    )
    .map_err(|e| e.to_string())?;
    writeln!(f, "Cartesian").map_err(|e| e.to_string())?;
    for atom in &sc.basis_atoms {
        writeln!(f, "  {:.8}  {:.8}  {:.8}", atom[0], atom[1], atom[2])
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Write a minimal CIF file for the crystal structure.
pub fn write_cif_minimal(crystal: &CrystalLattice, path: &str) -> Result<(), String> {
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    writeln!(f, "data_structure").map_err(|e| e.to_string())?;
    writeln!(f, "_cell_length_a {:.6}", crystal.lattice_vectors[0][0])
        .map_err(|e| e.to_string())?;
    writeln!(f, "_cell_length_b {:.6}", crystal.lattice_vectors[1][1])
        .map_err(|e| e.to_string())?;
    writeln!(f, "_cell_length_c {:.6}", crystal.lattice_vectors[2][2])
        .map_err(|e| e.to_string())?;
    writeln!(f, "_cell_angle_alpha 90").map_err(|e| e.to_string())?;
    writeln!(f, "_cell_angle_beta 90").map_err(|e| e.to_string())?;
    writeln!(f, "_cell_angle_gamma 90").map_err(|e| e.to_string())?;
    writeln!(f, "loop_").map_err(|e| e.to_string())?;
    writeln!(f, "_atom_site_label").map_err(|e| e.to_string())?;
    writeln!(f, "_atom_site_Cartn_x").map_err(|e| e.to_string())?;
    writeln!(f, "_atom_site_Cartn_y").map_err(|e| e.to_string())?;
    writeln!(f, "_atom_site_Cartn_z").map_err(|e| e.to_string())?;
    for (atom, atype) in crystal.basis_atoms.iter().zip(crystal.atom_types.iter()) {
        writeln!(f, "{} {:.6} {:.6} {:.6}", atype, atom[0], atom[1], atom[2])
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── LatticeType tests ─────────────────────────────────────────────────

    #[test]
    fn test_d2q9_velocity_count() {
        assert_eq!(LatticeType::D2Q9.velocity_count(), 9);
    }

    #[test]
    fn test_d3q19_velocity_count() {
        assert_eq!(LatticeType::D3Q19.velocity_count(), 19);
    }

    #[test]
    fn test_d3q27_velocity_count() {
        assert_eq!(LatticeType::D3Q27.velocity_count(), 27);
    }

    #[test]
    fn test_custom_velocity_count() {
        let vels: Vec<[i32; 3]> = (0..7).map(|i| [i, 0, 0]).collect();
        assert_eq!(LatticeType::Custom(vels).velocity_count(), 7);
    }

    #[test]
    fn test_fcc_velocity_count() {
        assert_eq!(LatticeType::FCC.velocity_count(), 12);
    }

    // ── LatticeGrid tests ──────────────────────────────────────────────────

    #[test]
    fn test_node_count() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 4, 5, 6, 1.0, 1);
        assert_eq!(g.node_count(), 120);
    }

    #[test]
    fn test_node_count_equal_dims() {
        let g = LatticeGrid::new(LatticeType::D2Q9, 3, 3, 3, 0.1, 2);
        assert_eq!(g.node_count(), 27);
    }

    #[test]
    fn test_index_origin() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 4, 4, 4, 1.0, 1);
        assert_eq!(g.index(0, 0, 0), 0);
    }

    #[test]
    fn test_index_last() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 4, 4, 4, 1.0, 1);
        assert_eq!(g.index(3, 3, 3), 63);
    }

    #[test]
    fn test_set_get_node_roundtrip() {
        let mut g = LatticeGrid::new(LatticeType::D3Q19, 3, 3, 3, 1.0, 3);
        g.set_node(1, 2, 0, &[1.0, 2.0, 3.0]);
        let data = g.get_node(1, 2, 0);
        assert!((data[0] - 1.0).abs() < 1e-12);
        assert!((data[1] - 2.0).abs() < 1e-12);
        assert!((data[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_set_node_initial_zero() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 2, 2, 2, 1.0, 2);
        let data = g.get_node(0, 0, 0);
        assert!(data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_node_data_length() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 2, 3, 4, 0.5, 5);
        assert_eq!(g.node_data.len(), 24);
        assert!(g.node_data.iter().all(|d| d.len() == 5));
    }

    // ── Binary I/O tests ──────────────────────────────────────────────────

    #[test]
    fn test_write_read_binary_roundtrip() {
        let mut g = LatticeGrid::new(LatticeType::D3Q19, 2, 2, 2, 0.5, 2);
        g.set_node(0, 1, 1, &[3.125, 2.72]);
        let path = "/tmp/lattice_test_roundtrip.bin";
        write_lattice_binary(&g, path).expect("write failed");
        let g2 = read_lattice_binary(path).expect("read failed");
        assert_eq!(g2.nx, 2);
        assert_eq!(g2.ny, 2);
        assert_eq!(g2.nz, 2);
        assert_eq!(g2.n_components, 2);
        let data = g2.get_node(0, 1, 1);
        assert!((data[0] - 3.125).abs() < 1e-12);
        assert!((data[1] - 2.72).abs() < 1e-12);
    }

    #[test]
    fn test_write_binary_creates_file() {
        let g = LatticeGrid::new(LatticeType::D2Q9, 2, 2, 1, 1.0, 1);
        write_lattice_binary(&g, "/tmp/lattice_create_test.bin").expect("write failed");
        assert!(std::path::Path::new("/tmp/lattice_create_test.bin").exists());
    }

    #[test]
    fn test_read_binary_invalid_magic() {
        use std::io::Write;
        let path = "/tmp/lattice_bad_magic.bin";
        let mut f = File::create(path).unwrap();
        f.write_all(b"XXXX").unwrap();
        let res = read_lattice_binary(path);
        assert!(res.is_err());
    }

    #[test]
    fn test_binary_roundtrip_node_count() {
        let g = LatticeGrid::new(LatticeType::D2Q9, 3, 4, 5, 1.0, 1);
        let path = "/tmp/lattice_nc_roundtrip.bin";
        write_lattice_binary(&g, path).expect("write failed");
        let g2 = read_lattice_binary(path).expect("read failed");
        assert_eq!(g2.node_count(), 60);
    }

    // ── VTK export tests ──────────────────────────────────────────────────

    #[test]
    fn test_write_lattice_vtk_creates_file() {
        let g = LatticeGrid::new(LatticeType::D3Q19, 3, 3, 3, 1.0, 1);
        write_lattice_vtk(&g, "/tmp/lattice_test.vtk", "density").expect("vtk failed");
        assert!(std::path::Path::new("/tmp/lattice_test.vtk").exists());
    }

    // ── XYZ export tests ──────────────────────────────────────────────────

    #[test]
    fn test_lattice_to_xyz_creates_file() {
        let g = LatticeGrid::new(LatticeType::SC, 2, 2, 2, 1.0, 1);
        lattice_to_xyz(&g, "/tmp/lattice_test.xyz").expect("xyz failed");
        assert!(std::path::Path::new("/tmp/lattice_test.xyz").exists());
    }

    // ── CrystalLattice tests ──────────────────────────────────────────────

    #[test]
    fn test_fcc_has_4_basis_atoms() {
        let c = CrystalLattice::fcc(4.05);
        assert_eq!(c.atom_count(), 4);
    }

    #[test]
    fn test_bcc_has_2_basis_atoms() {
        let c = CrystalLattice::bcc(3.3);
        assert_eq!(c.atom_count(), 2);
    }

    #[test]
    fn test_hcp_has_2_basis_atoms() {
        let c = CrystalLattice::hcp(3.2, 5.2);
        assert_eq!(c.atom_count(), 2);
    }

    #[test]
    fn test_fcc_lattice_constant() {
        let a = 4.05;
        let c = CrystalLattice::fcc(a);
        assert!((c.lattice_vectors[0][0] - a).abs() < 1e-12);
    }

    // ── Supercell tests ───────────────────────────────────────────────────

    #[test]
    fn test_supercell_atom_count_fcc_2x2x2() {
        let c = CrystalLattice::fcc(4.05);
        let sc = supercell(&c, 2, 2, 2);
        assert_eq!(sc.atom_count(), 4 * 8);
    }

    #[test]
    fn test_supercell_atom_count_bcc_3x1x1() {
        let c = CrystalLattice::bcc(3.3);
        let sc = supercell(&c, 3, 1, 1);
        assert_eq!(sc.atom_count(), 2 * 3);
    }

    #[test]
    fn test_supercell_1x1x1_equals_unit_cell() {
        let c = CrystalLattice::fcc(4.05);
        let sc = supercell(&c, 1, 1, 1);
        assert_eq!(sc.atom_count(), c.atom_count());
    }

    #[test]
    fn test_supercell_lattice_vectors_scale() {
        let c = CrystalLattice::bcc(3.3);
        let sc = supercell(&c, 2, 2, 2);
        assert!((sc.lattice_vectors[0][0] - 2.0 * c.lattice_vectors[0][0]).abs() < 1e-12);
    }

    // ── POSCAR / CIF tests ────────────────────────────────────────────────

    #[test]
    fn test_write_poscar_creates_file() {
        let c = CrystalLattice::fcc(4.05);
        write_poscar(&c, 1, "/tmp/lattice_test.poscar").expect("poscar failed");
        assert!(std::path::Path::new("/tmp/lattice_test.poscar").exists());
    }

    #[test]
    fn test_write_cif_creates_file() {
        let c = CrystalLattice::bcc(3.3);
        write_cif_minimal(&c, "/tmp/lattice_test.cif").expect("cif failed");
        assert!(std::path::Path::new("/tmp/lattice_test.cif").exists());
    }

    #[test]
    fn test_type_id_roundtrip() {
        for id in 0..=9u8 {
            let lt = LatticeType::from_type_id(id);
            if id < 9 {
                assert_eq!(lt.type_id(), id);
            }
        }
    }

    #[test]
    fn test_d3q19_type_id() {
        assert_eq!(LatticeType::D3Q19.type_id(), 3);
    }
}
