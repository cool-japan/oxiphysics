//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::types::{
    Atom, CifAtomSite, CifBlock, MolFormat, MolVizError, MoldenMo, PsfAtom, XsfAtom,
};

/// Convenience result alias.
pub type Result<T> = std::result::Result<T, MolVizError>;
/// A 3-D position in angstroms or Bohr radii (context-dependent).
pub type Vec3 = [f64; 3];
/// A 3×3 lattice matrix stored row-major: rows are lattice vectors a, b, c.
pub type Lattice = [[f64; 3]; 3];
pub(super) fn parse_psf_atom_line(line: &str, _xplor: bool) -> Result<PsfAtom> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < 8 {
        return Err(MolVizError::ParseError {
            message: format!("PSF atom line too short: '{line}'"),
        });
    }
    let serial: u32 = toks[0].parse().map_err(|_| MolVizError::ParseError {
        message: format!("bad serial: '{}'", toks[0]),
    })?;
    let segname = toks[1].to_string();
    let resid: i32 = toks[2].parse().unwrap_or(0);
    let resname = toks[3].to_string();
    let atomname = toks[4].to_string();
    let atom_type = toks[5].to_string();
    let charge: f64 = toks[6].parse().unwrap_or(0.0);
    let mass: f64 = toks[7].parse().unwrap_or(1.0);
    Ok(PsfAtom {
        serial,
        segname,
        resid,
        resname,
        atomname,
        atom_type,
        charge,
        mass,
    })
}
pub(super) fn voxel_volume(dx: Vec3, dy: Vec3, dz: Vec3) -> f64 {
    let cross = [
        dy[1] * dz[2] - dy[2] * dz[1],
        dy[2] * dz[0] - dy[0] * dz[2],
        dy[0] * dz[1] - dy[1] * dz[0],
    ];
    (dx[0] * cross[0] + dx[1] * cross[1] + dx[2] * cross[2]).abs()
}
pub(super) fn read_3x3_matrix(lines: &[&str], idx: &mut usize) -> Result<Lattice> {
    let mut mat = [[0.0f64; 3]; 3];
    for row in mat.iter_mut() {
        if *idx >= lines.len() {
            return Err(MolVizError::ParseError {
                message: "unexpected EOF reading 3×3 matrix".into(),
            });
        }
        let line = lines[*idx].trim();
        *idx += 1;
        let toks: Vec<f64> = line
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        if toks.len() >= 3 {
            row[0] = toks[0];
            row[1] = toks[1];
            row[2] = toks[2];
        }
    }
    Ok(mat)
}
pub(super) fn parse_xsf_atom(line: &str) -> Option<XsfAtom> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if toks.len() < 4 {
        return None;
    }
    let atomic_number: i32 = toks[0].parse().ok()?;
    let pos = [
        toks[1].parse().ok()?,
        toks[2].parse().ok()?,
        toks[3].parse().ok()?,
    ];
    let force = if toks.len() >= 7 {
        Some([
            toks[4].parse().unwrap_or(0.0),
            toks[5].parse().unwrap_or(0.0),
            toks[6].parse().unwrap_or(0.0),
        ])
    } else {
        None
    };
    Some(XsfAtom {
        atomic_number,
        position: pos,
        force,
    })
}
/// Parse a simplified CIF file into a list of data blocks.
pub fn parse_cif(content: &str) -> Result<Vec<CifBlock>> {
    let mut blocks: Vec<CifBlock> = Vec::new();
    let mut current: Option<CifBlock> = None;
    let mut in_loop = false;
    let mut loop_headers: Vec<String> = Vec::new();
    let mut loop_rows: Vec<Vec<String>> = Vec::new();
    for raw in content.lines() {
        let line = raw.trim();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if line.starts_with("data_") {
            if let Some(blk) = current.take() {
                blocks.push(blk);
            }
            let name = line.trim_start_matches("data_").to_string();
            current = Some(CifBlock {
                name,
                tags: HashMap::new(),
                atom_sites: Vec::new(),
            });
            in_loop = false;
            loop_headers.clear();
            loop_rows.clear();
            continue;
        }
        let Some(ref mut blk) = current else {
            continue;
        };
        if line == "loop_" {
            if !loop_headers.is_empty() {
                process_cif_loop(&loop_headers, &loop_rows, blk);
            }
            in_loop = true;
            loop_headers.clear();
            loop_rows.clear();
            continue;
        }
        if in_loop {
            if line.starts_with('_') {
                loop_headers.push(line.to_string());
            } else if !line.is_empty() {
                let row: Vec<String> = line.split_whitespace().map(|s| s.to_string()).collect();
                if !row.is_empty() {
                    loop_rows.push(row);
                }
            }
        } else if line.starts_with('_') {
            let mut parts = line.splitn(2, char::is_whitespace);
            let key = parts.next().unwrap_or("").to_string();
            let value = parts.next().unwrap_or("").trim().to_string();
            blk.tags.insert(key, value);
        }
    }
    if let Some(ref mut blk) = current
        && !loop_headers.is_empty()
    {
        process_cif_loop(&loop_headers, &loop_rows, blk);
    }
    if let Some(blk) = current {
        blocks.push(blk);
    }
    Ok(blocks)
}
pub(super) fn process_cif_loop(headers: &[String], rows: &[Vec<String>], blk: &mut CifBlock) {
    let is_atom_site = headers.iter().any(|h| h.contains("_atom_site"));
    if !is_atom_site {
        return;
    }
    let idx_label = headers.iter().position(|h| h.contains("label"));
    let idx_type = headers
        .iter()
        .position(|h| h.contains("type_symbol") || h.contains("element"));
    let idx_x = headers.iter().position(|h| h.ends_with("_x"));
    let idx_y = headers.iter().position(|h| h.ends_with("_y"));
    let idx_z = headers.iter().position(|h| h.ends_with("_z"));
    let idx_occ = headers.iter().position(|h| h.contains("occupancy"));
    for row in rows {
        let label = idx_label
            .and_then(|i| row.get(i))
            .cloned()
            .unwrap_or_else(|| "X".to_string());
        let element = idx_type
            .and_then(|i| row.get(i))
            .cloned()
            .unwrap_or_else(|| label.chars().take(2).collect());
        let fx = idx_x
            .and_then(|i| row.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let fy = idx_y
            .and_then(|i| row.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let fz = idx_z
            .and_then(|i| row.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let occupancy = idx_occ
            .and_then(|i| row.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0);
        blk.atom_sites.push(CifAtomSite {
            label,
            element,
            fx,
            fy,
            fz,
            occupancy,
        });
    }
}
pub(super) fn parse_molden_kv(line: &str, mo: &mut MoldenMo) {
    if let Some(rest) = line.strip_prefix("Sym=") {
        mo.symmetry = rest.trim().to_string();
    } else if let Some(rest) = line.strip_prefix("Ene=") {
        mo.energy = rest.trim().parse().unwrap_or(0.0);
    } else if let Some(rest) = line.strip_prefix("Spin=") {
        mo.spin = rest.trim().to_string();
    } else if let Some(rest) = line.strip_prefix("Occup=") {
        mo.occupation = rest.trim().parse().unwrap_or(0.0);
    }
}
pub(super) fn lattice_volume(lat: Lattice) -> f64 {
    let a = lat[0];
    let b = lat[1];
    let c = lat[2];
    let cross = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]).abs()
}
/// Attempt to detect the molecular file format from filename and/or content.
///
/// Detection priority: filename extension > content magic.
pub fn detect_format(filename: &str, content: &str) -> MolFormat {
    let lower = filename.to_lowercase();
    if lower.ends_with(".psf") {
        return MolFormat::Psf;
    }
    if lower.ends_with(".dx") {
        return MolFormat::Dx;
    }
    if lower.ends_with(".cube") {
        return MolFormat::Cube;
    }
    if lower.ends_with(".xsf") {
        return MolFormat::Xsf;
    }
    if lower.ends_with(".cif") {
        return MolFormat::Cif;
    }
    if lower.ends_with(".molden") || lower.ends_with(".mold") {
        return MolFormat::Molden;
    }
    let basename = lower
        .split('/')
        .next_back()
        .unwrap_or(&lower)
        .split('\\')
        .next_back()
        .unwrap_or(&lower);
    if basename == "poscar" || basename == "contcar" {
        return MolFormat::Poscar;
    }
    if basename == "chgcar" || basename == "parchg" || basename == "elfcar" {
        return MolFormat::Chgcar;
    }
    let first = content.lines().next().unwrap_or("").trim();
    if first.starts_with("PSF") {
        return MolFormat::Psf;
    }
    if content.contains("object 1 class gridpositions") {
        return MolFormat::Dx;
    }
    if content.contains("[Molden Format]") || content.contains("[MOLDEN FORMAT]") {
        return MolFormat::Molden;
    }
    if content.contains("CRYSTAL") || content.contains("PRIMVEC") {
        return MolFormat::Xsf;
    }
    if content.contains("data_") && content.contains("_cell_length") {
        return MolFormat::Cif;
    }
    MolFormat::Unknown
}
/// Write a list of atoms to a minimal XSF molecule file.
pub fn write_xsf_molecule(atoms: &[Atom]) -> String {
    let mut out = String::from("MOLECULE\nATOMS\n");
    for atom in atoms {
        let z = element_to_z(&atom.element);
        out.push_str(&format!(
            "{} {:.6} {:.6} {:.6}\n",
            z, atom.position[0], atom.position[1], atom.position[2]
        ));
    }
    out
}
/// Write atoms and lattice to a minimal XSF crystal file.
pub fn write_xsf_crystal(atoms: &[Atom], lattice: Lattice) -> String {
    let mut out = String::from("CRYSTAL\nPRIMVEC\n");
    for row in &lattice {
        out.push_str(&format!(" {:.6} {:.6} {:.6}\n", row[0], row[1], row[2]));
    }
    out.push_str(&format!("PRIMCOORD\n{} 1\n", atoms.len()));
    for atom in atoms {
        let z = element_to_z(&atom.element);
        out.push_str(&format!(
            "{} {:.6} {:.6} {:.6}\n",
            z, atom.position[0], atom.position[1], atom.position[2]
        ));
    }
    out
}
/// Write a simple POSCAR file from atoms and lattice.
pub fn write_poscar(comment: &str, lattice: Lattice, atoms: &[Atom]) -> String {
    let mut out = format!("{comment}\n1.0\n");
    for row in &lattice {
        out.push_str(&format!("  {:.10} {:.10} {:.10}\n", row[0], row[1], row[2]));
    }
    let mut element_order: Vec<String> = Vec::new();
    for a in atoms {
        if !element_order.contains(&a.element) {
            element_order.push(a.element.clone());
        }
    }
    out.push_str(&element_order.join("  "));
    out.push('\n');
    let counts: Vec<usize> = element_order
        .iter()
        .map(|e| atoms.iter().filter(|a| &a.element == e).count())
        .collect();
    out.push_str(
        &counts
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join("  "),
    );
    out.push_str("\nDirect\n");
    out.push_str("Cartesian\n");
    for e in &element_order {
        for a in atoms.iter().filter(|a| &a.element == e) {
            out.push_str(&format!(
                "  {:.10} {:.10} {:.10}\n",
                a.position[0], a.position[1], a.position[2]
            ));
        }
    }
    out
}
/// Map common element symbols to atomic numbers (returns 0 if unknown).
pub fn element_to_z(symbol: &str) -> i32 {
    match symbol.trim() {
        "H" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Fe" => 26,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Au" => 79,
        "Pt" => 78,
        _ => 0,
    }
}
/// Map an atomic number to an element symbol.
pub fn z_to_element(z: i32) -> &'static str {
    match z {
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        10 => "Ne",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        26 => "Fe",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        78 => "Pt",
        79 => "Au",
        _ => "X",
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecular_visualization_io::types::*;
    fn sample_psf() -> &'static str {
        "PSF\n\n       1 !NTITLE\n REMARKS test PSF file\n\n       3 !NATOM\n       1 SEG1     1 ALA  CA  CT1   -0.270000     12.0110\n       2 SEG1     1 ALA  N   NH1   -0.470000     14.0070\n       3 SEG1     1 ALA  C   C      0.510000     12.0110\n\n       2 !NBOND\n       1       2       2       3\n\n       1 !NTHETA\n       1       2       3\n"
    }
    #[test]
    fn test_psf_parse_atoms() {
        let psf = PsfFile::parse(sample_psf()).unwrap();
        assert_eq!(psf.atoms.len(), 3);
        assert_eq!(psf.atoms[0].atomname, "CA");
        assert_eq!(psf.atoms[1].atomname, "N");
    }
    #[test]
    fn test_psf_parse_bonds() {
        let psf = PsfFile::parse(sample_psf()).unwrap();
        assert_eq!(psf.bonds.len(), 2);
        assert_eq!(psf.bonds[0], PsfBond { i: 1, j: 2 });
        assert_eq!(psf.bonds[1], PsfBond { i: 2, j: 3 });
    }
    #[test]
    fn test_psf_parse_angles() {
        let psf = PsfFile::parse(sample_psf()).unwrap();
        assert_eq!(psf.angles.len(), 1);
        assert_eq!(psf.angles[0].i, 1);
        assert_eq!(psf.angles[0].j, 2);
        assert_eq!(psf.angles[0].k, 3);
    }
    #[test]
    fn test_psf_total_charge() {
        let psf = PsfFile::parse(sample_psf()).unwrap();
        let total = psf.total_charge();
        assert!((total - (-0.23)).abs() < 1e-6, "total charge = {total}");
    }
    #[test]
    fn test_psf_segments() {
        let psf = PsfFile::parse(sample_psf()).unwrap();
        let segs = psf.segments();
        assert!(segs.contains_key("SEG1"));
        assert_eq!(segs["SEG1"].len(), 3);
    }
    #[test]
    fn test_psf_parse_non_psf_returns_error() {
        let result = PsfFile::parse("NOT A PSF FILE\nsome data\n");
        assert!(result.is_err());
    }
    fn sample_dx() -> &'static str {
        "# Data from APBS\nobject 1 class gridpositions counts 2 2 2\norigin 0.0 0.0 0.0\ndelta 1.0 0.0 0.0\ndelta 0.0 1.0 0.0\ndelta 0.0 0.0 1.0\nobject 2 class array type double rank 0 items 8 data follows\n1.0 2.0 3.0\n4.0 5.0 6.0\n7.0 8.0\nattribute \"dep\" string \"positions\"\n"
    }
    #[test]
    fn test_dx_parse_dims() {
        let dx = DxGrid::parse(sample_dx()).unwrap();
        assert_eq!(dx.nx, 2);
        assert_eq!(dx.ny, 2);
        assert_eq!(dx.nz, 2);
    }
    #[test]
    fn test_dx_parse_data_count() {
        let dx = DxGrid::parse(sample_dx()).unwrap();
        assert_eq!(dx.data.len(), 8);
    }
    #[test]
    fn test_dx_get_set() {
        let mut dx = DxGrid::empty(2, 2, 2);
        dx.set(1, 0, 1, 42.0);
        assert!((dx.get(1, 0, 1) - 42.0).abs() < 1e-12);
    }
    #[test]
    fn test_dx_range() {
        let dx = DxGrid::parse(sample_dx()).unwrap();
        let (mn, mx) = dx.range();
        assert!((mn - 1.0).abs() < 1e-12);
        assert!((mx - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_dx_coord() {
        let dx = DxGrid::empty(3, 3, 3);
        let c = dx.coord(1, 1, 1);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_dx_roundtrip() {
        let dx = DxGrid::parse(sample_dx()).unwrap();
        let s = dx.to_dx_string();
        let dx2 = DxGrid::parse(&s).unwrap();
        assert_eq!(dx2.nx, dx.nx);
        assert_eq!(dx2.data.len(), dx.data.len());
    }
    fn sample_cube() -> &'static str {
        "Electron density\nGenerated by OxiPhysics\n2 0.0 0.0 0.0\n2 0.5 0.0 0.0\n2 0.0 0.5 0.0\n2 0.0 0.0 0.5\n6 6.0 0.0 0.0 0.0\n8 8.0 1.2 0.0 0.0\n1.0 2.0 3.0 4.0\n5.0 6.0 7.0 8.0\n"
    }
    #[test]
    fn test_cube_parse_atoms() {
        let cube = CubeFile::parse(sample_cube()).unwrap();
        assert_eq!(cube.n_atoms, 2);
        assert_eq!(cube.atoms[0].0, 6);
        assert_eq!(cube.atoms[1].0, 8);
    }
    #[test]
    fn test_cube_parse_grid() {
        let cube = CubeFile::parse(sample_cube()).unwrap();
        assert_eq!(cube.nx, 2);
        assert_eq!(cube.ny, 2);
        assert_eq!(cube.nz, 2);
    }
    #[test]
    fn test_cube_data_len() {
        let cube = CubeFile::parse(sample_cube()).unwrap();
        assert_eq!(cube.data.len(), 8);
        assert!((cube.data[0] - 1.0).abs() < 1e-12);
        assert!((cube.data[7] - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_cube_integrate_positive() {
        let cube = CubeFile::parse(sample_cube()).unwrap();
        let integral = cube.integrate();
        assert!(integral > 0.0);
    }
    #[test]
    fn test_cube_position() {
        let cube = CubeFile::parse(sample_cube()).unwrap();
        let pos = cube.position(0, 0, 0);
        assert!((pos[0] - 0.0).abs() < 1e-12);
        assert!((pos[1] - 0.0).abs() < 1e-12);
        assert!((pos[2] - 0.0).abs() < 1e-12);
    }
    fn sample_xsf_crystal() -> &'static str {
        "CRYSTAL\nPRIMVEC\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nPRIMCOORD\n2 1\n26 0.0 0.0 0.0\n26 2.0 2.0 2.0\n"
    }
    #[test]
    fn test_xsf_parse_crystal() {
        let xsf = XsfFile::parse(sample_xsf_crystal()).unwrap();
        assert_eq!(xsf.periodicity, XsfPeriodicity::Crystal);
        assert!(xsf.lattice.is_some());
    }
    #[test]
    fn test_xsf_parse_atoms() {
        let xsf = XsfFile::parse(sample_xsf_crystal()).unwrap();
        assert_eq!(xsf.atoms.len(), 2);
        assert_eq!(xsf.atoms[0].atomic_number, 26);
    }
    #[test]
    fn test_xsf_write_and_detect() {
        let atoms = vec![Atom::new("Fe", [0.0, 0.0, 0.0])];
        let lat = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let s = write_xsf_crystal(&atoms, lat);
        assert!(s.contains("CRYSTAL"));
        assert!(s.contains("PRIMVEC"));
    }
    #[test]
    fn test_xsf_molecule() {
        let atoms = vec![
            Atom::new("C", [0.0, 0.0, 0.0]),
            Atom::new("H", [1.0, 0.0, 0.0]),
        ];
        let s = write_xsf_molecule(&atoms);
        assert!(s.contains("MOLECULE"));
        assert!(s.contains("ATOMS"));
    }
    fn sample_poscar() -> &'static str {
        "BCC Fe\n1.0\n2.87 0.00 0.00\n0.00 2.87 0.00\n0.00 0.00 2.87\nFe\n2\nDirect\n0.00 0.00 0.00\n0.50 0.50 0.50\n"
    }
    #[test]
    fn test_poscar_parse_lattice() {
        let p = PoscarFile::parse(sample_poscar()).unwrap();
        assert!((p.lattice[0][0] - 2.87).abs() < 1e-10);
    }
    #[test]
    fn test_poscar_parse_natoms() {
        let p = PoscarFile::parse(sample_poscar()).unwrap();
        assert_eq!(p.n_atoms(), 2);
    }
    #[test]
    fn test_poscar_cartesian_positions() {
        let p = PoscarFile::parse(sample_poscar()).unwrap();
        let cart = p.cartesian_positions();
        assert_eq!(cart.len(), 2);
        assert!(cart[0][0].abs() < 1e-10);
        assert!((cart[1][0] - 1.435).abs() < 1e-3);
    }
    #[test]
    fn test_poscar_write_contains_lattice() {
        let atoms = vec![Atom::new("Fe", [0.0, 0.0, 0.0])];
        let lat = [[2.87, 0.0, 0.0], [0.0, 2.87, 0.0], [0.0, 0.0, 2.87]];
        let s = write_poscar("test", lat, &atoms);
        assert!(s.contains("2.87"));
        assert!(s.contains("Fe"));
    }
    fn sample_cif() -> &'static str {
        "data_test\n_cell_length_a 4.0\n_cell_length_b 4.0\n_cell_length_c 4.0\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\nFe1 Fe 0.0 0.0 0.0 1.0\nFe2 Fe 0.5 0.5 0.5 1.0\n"
    }
    #[test]
    fn test_cif_parse_block() {
        let blocks = parse_cif(sample_cif()).unwrap();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].name, "test");
    }
    #[test]
    fn test_cif_parse_atom_sites() {
        let blocks = parse_cif(sample_cif()).unwrap();
        assert_eq!(blocks[0].atom_sites.len(), 2);
    }
    #[test]
    fn test_cif_atom_site_values() {
        let blocks = parse_cif(sample_cif()).unwrap();
        let site0 = &blocks[0].atom_sites[0];
        assert_eq!(site0.element, "Fe");
        assert!((site0.fx - 0.0).abs() < 1e-12);
        assert!((site0.occupancy - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_cif_parse_tags() {
        let blocks = parse_cif(sample_cif()).unwrap();
        assert!(blocks[0].tags.contains_key("_cell_length_a"));
    }
    fn make_frame(idx: usize, t: f64, n: usize) -> TrajectoryFrame {
        let positions = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        TrajectoryFrame {
            index: idx,
            time: t,
            positions,
            velocities: None,
            cell: None,
        }
    }
    #[test]
    fn test_trajectory_push_and_len() {
        let mut traj = Trajectory::new(3);
        traj.push_frame(make_frame(0, 0.0, 3)).unwrap();
        traj.push_frame(make_frame(1, 1.0, 3)).unwrap();
        assert_eq!(traj.len(), 2);
    }
    #[test]
    fn test_trajectory_push_wrong_natoms() {
        let mut traj = Trajectory::new(3);
        let result = traj.push_frame(make_frame(0, 0.0, 2));
        assert!(result.is_err());
    }
    #[test]
    fn test_trajectory_slice() {
        let mut traj = Trajectory::new(2);
        for i in 0..5 {
            traj.push_frame(make_frame(i, i as f64, 2)).unwrap();
        }
        let sliced = traj.slice(1, 4).unwrap();
        assert_eq!(sliced.len(), 3);
    }
    #[test]
    fn test_trajectory_concat() {
        let mut t1 = Trajectory::new(2);
        let mut t2 = Trajectory::new(2);
        for i in 0..3 {
            t1.push_frame(make_frame(i, i as f64, 2)).unwrap();
        }
        for i in 0..2 {
            t2.push_frame(make_frame(i, i as f64, 2)).unwrap();
        }
        let combined = t1.concat(&t2).unwrap();
        assert_eq!(combined.len(), 5);
    }
    #[test]
    fn test_trajectory_stride() {
        let mut traj = Trajectory::new(1);
        for i in 0..10 {
            traj.push_frame(make_frame(i, i as f64, 1)).unwrap();
        }
        let strided = traj.stride(3);
        assert_eq!(strided.len(), 4);
    }
    #[test]
    fn test_trajectory_msd_stationary() {
        let mut traj = Trajectory::new(1);
        for i in 0..5 {
            let frame = TrajectoryFrame {
                index: i,
                time: i as f64,
                positions: vec![[0.0, 0.0, 0.0]],
                velocities: None,
                cell: None,
            };
            traj.push_frame(frame).unwrap();
        }
        let msd = traj.msd_atom(0).unwrap();
        for v in &msd {
            assert!(v.abs() < 1e-12);
        }
    }
    #[test]
    fn test_trajectory_msd_moving() {
        let mut traj = Trajectory::new(1);
        for i in 0..4 {
            let x = i as f64;
            let frame = TrajectoryFrame {
                index: i,
                time: i as f64,
                positions: vec![[x, 0.0, 0.0]],
                velocities: None,
                cell: None,
            };
            traj.push_frame(frame).unwrap();
        }
        let msd = traj.msd_atom(0).unwrap();
        assert!((msd[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_detect_psf_by_extension() {
        assert_eq!(detect_format("protein.psf", ""), MolFormat::Psf);
    }
    #[test]
    fn test_detect_dx_by_extension() {
        assert_eq!(detect_format("grid.dx", ""), MolFormat::Dx);
    }
    #[test]
    fn test_detect_cube_by_extension() {
        assert_eq!(detect_format("density.cube", ""), MolFormat::Cube);
    }
    #[test]
    fn test_detect_xsf_by_extension() {
        assert_eq!(detect_format("crystal.xsf", ""), MolFormat::Xsf);
    }
    #[test]
    fn test_detect_poscar_by_name() {
        assert_eq!(detect_format("POSCAR", ""), MolFormat::Poscar);
        assert_eq!(detect_format("CONTCAR", ""), MolFormat::Poscar);
    }
    #[test]
    fn test_detect_cif_by_extension() {
        assert_eq!(detect_format("struct.cif", ""), MolFormat::Cif);
    }
    #[test]
    fn test_detect_molden_by_extension() {
        assert_eq!(detect_format("orbs.molden", ""), MolFormat::Molden);
    }
    #[test]
    fn test_detect_chgcar_by_name() {
        assert_eq!(detect_format("CHGCAR", ""), MolFormat::Chgcar);
    }
    #[test]
    fn test_detect_by_content_psf() {
        assert_eq!(detect_format("unknown", "PSF\n\n"), MolFormat::Psf);
    }
    #[test]
    fn test_detect_unknown() {
        assert_eq!(
            detect_format("file.xyz", "some random content"),
            MolFormat::Unknown
        );
    }
    #[test]
    fn test_element_to_z() {
        assert_eq!(element_to_z("H"), 1);
        assert_eq!(element_to_z("C"), 6);
        assert_eq!(element_to_z("Fe"), 26);
        assert_eq!(element_to_z("Au"), 79);
        assert_eq!(element_to_z("Xx"), 0);
    }
    #[test]
    fn test_z_to_element() {
        assert_eq!(z_to_element(1), "H");
        assert_eq!(z_to_element(6), "C");
        assert_eq!(z_to_element(26), "Fe");
        assert_eq!(z_to_element(999), "X");
    }
    #[test]
    fn test_atom_distance() {
        let a = Atom::new("C", [0.0, 0.0, 0.0]);
        let b = Atom::new("C", [3.0, 4.0, 0.0]);
        let d = a.distance_to(&b);
        assert!((d - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_voxel_volume_unit_cell() {
        let dx = [1.0, 0.0, 0.0];
        let dy = [0.0, 1.0, 0.0];
        let dz = [0.0, 0.0, 1.0];
        let vol = voxel_volume(dx, dy, dz);
        assert!((vol - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_lattice_volume_cube() {
        let lat = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        let vol = lattice_volume(lat);
        assert!((vol - 64.0).abs() < 1e-10);
    }
    #[test]
    fn test_trajectory_slice_out_of_bounds() {
        let traj = Trajectory::new(1);
        let result = traj.slice(0, 5);
        assert!(result.is_err());
    }
    #[test]
    fn test_trajectory_concat_natoms_mismatch() {
        let t1 = Trajectory::new(2);
        let t2 = Trajectory::new(3);
        assert!(t1.concat(&t2).is_err());
    }
    #[test]
    fn test_dx_empty_grid() {
        let dx = DxGrid::empty(4, 4, 4);
        assert_eq!(dx.data.len(), 64);
        let (mn, mx) = dx.range();
        assert!((mn - 0.0).abs() < 1e-12);
        assert!((mx - 0.0).abs() < 1e-12);
    }
}
