//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::types::{AmberCoordinates, AmberRestraint, AmberTopology};

#[cfg(test)]
mod tests_amber_ext {
    use super::*;
    use crate::amber::AmberAtom;
    use crate::amber::AmberBox;
    use crate::amber::AmberDcd;
    use crate::amber::AmberMask;
    use crate::amber::AmberRst7;
    use crate::amber::FrcmodFile;
    use crate::amber::*;
    use std::str::FromStr;
    fn sample_frcmod() -> &'static str {
        "Modified force field for custom ligand\n\
         MASS\n\
         CX  12.011\n\
         \n\
         BOND\n\
         CX-HC  340.0  1.090\n\
         CX-CX  310.0  1.526\n\
         \n\
         ANGL\n\
         HC-CX-HC  35.0  109.5\n\
         CX-CX-HC  50.0  109.5\n\
         \n\
         DIHE\n\
         X-CX-CX-X   9  1.4  0.0  3.0\n\
         \n\
         NONB\n\
         CX    1.9080  0.1094\n\
         \n\
         END\n"
    }
    #[test]
    fn test_frcmod_parse_title() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert!(frc.title.contains("Modified"), "title={}", frc.title);
    }
    #[test]
    fn test_frcmod_bonds_count() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert_eq!(frc.n_bonds(), 2, "bonds={}", frc.n_bonds());
    }
    #[test]
    fn test_frcmod_bond_values() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        let b = frc.get_bond("CX-HC").unwrap();
        assert!((b.k - 340.0).abs() < 1.0, "k={}", b.k);
        assert!((b.r0 - 1.090).abs() < 0.001, "r0={}", b.r0);
    }
    #[test]
    fn test_frcmod_bond_reverse_lookup() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert!(frc.get_bond("HC-CX").is_some());
    }
    #[test]
    fn test_frcmod_angles_count() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert_eq!(frc.n_angles(), 2, "angles={}", frc.n_angles());
    }
    #[test]
    fn test_frcmod_angle_values() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        let a = &frc.angles[0];
        assert!((a.k - 35.0).abs() < 1.0, "k={}", a.k);
        assert!(
            (a.theta0_deg - 109.5).abs() < 0.1,
            "theta0={}",
            a.theta0_deg
        );
    }
    #[test]
    fn test_frcmod_dihedrals_count() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert_eq!(frc.n_dihedrals(), 1, "dihedrals={}", frc.n_dihedrals());
    }
    #[test]
    fn test_frcmod_nonbonded_count() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        assert_eq!(frc.n_nonbonded(), 1, "nonbonded={}", frc.n_nonbonded());
    }
    #[test]
    fn test_frcmod_nonbonded_values() {
        let frc = FrcmodFile::from_str(sample_frcmod()).unwrap();
        let nb = frc.get_nonbonded("CX").unwrap();
        assert!((nb.r_star - 1.9080).abs() < 0.001, "r_star={}", nb.r_star);
        assert!(
            (nb.epsilon - 0.1094).abs() < 0.001,
            "epsilon={}",
            nb.epsilon
        );
    }
    #[test]
    fn test_frcmod_empty_file() {
        let frc = FrcmodFile::from_str("").unwrap();
        assert_eq!(frc.n_bonds(), 0);
        assert_eq!(frc.n_angles(), 0);
    }
    #[test]
    fn test_mask_single_residue() {
        let m = AmberMask::parse(":5").unwrap();
        assert_eq!(m.residues, vec![5]);
        assert!(m.atom_names.is_empty());
    }
    #[test]
    fn test_mask_residue_range() {
        let m = AmberMask::parse(":1-5").unwrap();
        assert_eq!(m.residues, vec![1, 2, 3, 4, 5]);
    }
    #[test]
    fn test_mask_residue_list() {
        let m = AmberMask::parse(":1,3,5").unwrap();
        assert_eq!(m.residues, vec![1, 3, 5]);
    }
    #[test]
    fn test_mask_atom_name() {
        let m = AmberMask::parse("@CA").unwrap();
        assert!(m.residues.is_empty());
        assert_eq!(m.atom_names, vec!["CA"]);
    }
    #[test]
    fn test_mask_multiple_atom_names() {
        let m = AmberMask::parse("@CA,N,C").unwrap();
        assert_eq!(m.atom_names.len(), 3);
        assert!(m.atom_names.contains(&"CA".to_string()));
        assert!(m.atom_names.contains(&"N".to_string()));
    }
    #[test]
    fn test_mask_residue_and_atom() {
        let m = AmberMask::parse(":1-10@CA").unwrap();
        assert_eq!(m.residues.len(), 10);
        assert_eq!(m.atom_names, vec!["CA"]);
    }
    #[test]
    fn test_mask_matches() {
        let m = AmberMask::parse(":1-3@CA").unwrap();
        assert!(m.matches(1, "CA"));
        assert!(m.matches(3, "CA"));
        assert!(!m.matches(4, "CA"));
        assert!(!m.matches(1, "N"));
    }
    #[test]
    fn test_mask_empty() {
        let m = AmberMask::parse("").unwrap();
        assert!(m.matches_residue(1));
        assert!(m.matches_atom("CA"));
    }
    #[test]
    fn test_mask_invalid_prefix() {
        assert!(AmberMask::parse("1-5@CA").is_err());
    }
    #[test]
    fn test_amber_box_cubic() {
        let b = AmberBox::cubic(50.0);
        assert!(b.is_cubic());
        assert!(b.is_orthorhombic());
        assert!((b.volume() - 50.0_f64.powi(3)).abs() < 1e-5);
    }
    #[test]
    fn test_amber_box_orthorhombic() {
        let b = AmberBox::orthorhombic(10.0, 20.0, 30.0);
        assert!(!b.is_cubic());
        assert!(b.is_orthorhombic());
        assert!((b.volume() - 6000.0).abs() < 1e-6);
    }
    #[test]
    fn test_amber_box_volume_cubic() {
        let b = AmberBox::cubic(10.0);
        assert!((b.volume() - 1000.0).abs() < 1e-8);
    }
    #[test]
    fn test_amber_box_in_nm() {
        let b = AmberBox::cubic(50.0);
        let nm = b.in_nm();
        assert!((nm[0] - 5.0).abs() < 1e-10);
    }
    fn sample_mdcrd() -> &'static str {
        "test trajectory\n\
           1.0000000  2.0000000  3.0000000  4.0000000  5.0000000  6.0000000\n"
    }
    #[test]
    fn test_mdcrd_parse_one_frame() {
        let frames = parse_mdcrd(sample_mdcrd(), 2, false);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].n_atoms(), 2);
    }
    #[test]
    fn test_mdcrd_position() {
        let frames = parse_mdcrd(sample_mdcrd(), 2, false);
        let p0 = frames[0].position(0).unwrap();
        assert!((p0[0] - 1.0).abs() < 1e-6);
        assert!((p0[1] - 2.0).abs() < 1e-6);
        assert!((p0[2] - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_mdcrd_position_atom1() {
        let frames = parse_mdcrd(sample_mdcrd(), 2, false);
        let p1 = frames[0].position(1).unwrap();
        assert!((p1[0] - 4.0).abs() < 1e-6);
    }
    #[test]
    fn test_mdcrd_no_frames_on_empty() {
        let frames = parse_mdcrd("title\n", 3, false);
        assert!(frames.is_empty());
    }
    #[test]
    fn test_amber99sb_bond_ct_ct() {
        let bonds = amber99sb_bonds();
        let ct_ct = bonds.iter().find(|b| b.type_a == "CT" && b.type_b == "CT");
        assert!(ct_ct.is_some());
        let b = ct_ct.unwrap();
        assert!((b.r0 - 1.526).abs() < 0.001);
    }
    #[test]
    fn test_amber99sb_angle_hc_ct_hc() {
        let angles = amber99sb_angles();
        let a = angles
            .iter()
            .find(|a| a.type_i == "HC" && a.type_j == "CT" && a.type_k == "HC");
        assert!(a.is_some());
        assert!((a.unwrap().theta0_deg - 109.5).abs() < 0.1);
    }
    #[test]
    fn test_write_prmtop_contains_title() {
        let topo = AmberTopology {
            title: "TestMol".into(),
            atoms: vec![],
            bonds: vec![],
            angles: vec![],
            n_atoms: 0,
            n_bonds: 0,
        };
        let text = write_prmtop(&topo);
        assert!(text.contains("%FLAG TITLE"), "should contain FLAG TITLE");
        assert!(text.contains("TestMol"), "should contain the title string");
    }
    #[test]
    fn test_write_prmtop_atom_names() {
        let topo = AmberTopology {
            title: "mol".into(),
            atoms: vec![
                AmberAtom {
                    name: "CA".into(),
                    residue_name: "ALA".into(),
                    charge: 0.1,
                    mass: 12.0,
                    atom_type: "CT".into(),
                },
                AmberAtom {
                    name: "N".into(),
                    residue_name: "ALA".into(),
                    charge: -0.4,
                    mass: 14.0,
                    atom_type: "N".into(),
                },
            ],
            bonds: vec![],
            angles: vec![],
            n_atoms: 2,
            n_bonds: 0,
        };
        let text = write_prmtop(&topo);
        assert!(text.contains("%FLAG ATOM_NAME"));
        assert!(text.contains("CA  ") || text.contains("CA"));
    }
    #[test]
    fn test_write_prmtop_pointers_section() {
        let topo = AmberTopology {
            title: "t".into(),
            atoms: vec![AmberAtom {
                name: "C".into(),
                residue_name: "R".into(),
                charge: 0.0,
                mass: 12.0,
                atom_type: "C".into(),
            }],
            bonds: vec![],
            angles: vec![],
            n_atoms: 1,
            n_bonds: 0,
        };
        let text = write_prmtop(&topo);
        assert!(text.contains("%FLAG POINTERS"));
    }
    #[test]
    fn test_amber_rst7_write_velocity_ok() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let mut rst7 = AmberRst7::new("test", positions);
        let velocities = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        assert!(rst7.write_velocity(velocities).is_ok());
        assert!(rst7.has_velocities());
    }
    #[test]
    fn test_amber_rst7_write_velocity_mismatch_error() {
        let positions = vec![[1.0, 2.0, 3.0]];
        let mut rst7 = AmberRst7::new("test", positions);
        let result = rst7.write_velocity(vec![[0.0; 3], [0.0; 3]]);
        assert!(result.is_err());
    }
    #[test]
    fn test_amber_rst7_to_string_with_velocity() {
        let pos = vec![[1.0, 2.0, 3.0]];
        let mut rst7 = AmberRst7::new("mol", pos);
        rst7.write_velocity(vec![[0.5, 0.6, 0.7]]).unwrap();
        let text = rst7.to_string_repr();
        assert!(text.contains("mol"));
        assert!(text.contains("1.0000000") || text.contains("1."));
    }
    #[test]
    fn test_amber_dcd_write_frame_roundtrip() {
        let mut dcd = AmberDcd::new(2);
        let frame0 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        dcd.write_frame(&frame0).unwrap();
        assert_eq!(dcd.n_frames(), 1);
        let pos = dcd.get_position(0, 0).unwrap();
        assert!((pos[0] - 1.0f32).abs() < 1e-5);
    }
    #[test]
    fn test_amber_dcd_write_frame_wrong_atom_count() {
        let mut dcd = AmberDcd::new(3);
        let result = dcd.write_frame(&[[0.0; 3], [0.0; 3]]);
        assert!(result.is_err());
    }
    #[test]
    fn test_amber_dcd_bytes_roundtrip() {
        let mut dcd = AmberDcd::new(2);
        dcd.write_frame(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
            .unwrap();
        dcd.write_frame(&[[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
            .unwrap();
        let bytes = dcd.to_bytes();
        let loaded = AmberDcd::from_bytes(&bytes).unwrap();
        assert_eq!(loaded.n_frames(), 2);
        assert_eq!(loaded.n_atoms, 2);
        let p = loaded.get_position(1, 0).unwrap();
        assert!((p[0] - 2.0f32).abs() < 1e-5);
    }
}
/// Validate that an `AmberTopology` is self-consistent.
///
/// Checks include: atom indices in bonds are in range, bond force constants
/// are non-negative, and angle equilibrium values are in [0, π].
pub fn validate_topology(topo: &AmberTopology) -> Vec<String> {
    let mut issues = Vec::new();
    let n = topo.atoms.len();
    for (idx, bond) in topo.bonds.iter().enumerate() {
        if bond.i >= n || bond.j >= n {
            issues.push(format!(
                "Bond {}: atom index out of range (i={}, j={}, n_atoms={})",
                idx, bond.i, bond.j, n
            ));
        }
        if bond.k < 0.0 {
            issues.push(format!("Bond {}: negative force constant {}", idx, bond.k));
        }
        if bond.r0 <= 0.0 {
            issues.push(format!(
                "Bond {}: non-positive equilibrium length {}",
                idx, bond.r0
            ));
        }
    }
    for (idx, angle) in topo.angles.iter().enumerate() {
        if angle.i >= n || angle.j >= n {
            issues.push(format!("Angle {}: atom index out of range", idx));
        }
        if angle.theta0 < 0.0 || angle.theta0 > std::f64::consts::PI {
            issues.push(format!(
                "Angle {}: equilibrium angle {} rad out of range [0, π]",
                idx, angle.theta0
            ));
        }
    }
    for atom in &topo.atoms {
        if atom.mass < 0.0 {
            issues.push(format!(
                "Atom '{}' has negative mass {}",
                atom.name, atom.mass
            ));
        }
    }
    issues
}
/// Count the number of bonds for each atom index.
pub fn atom_bond_counts(topo: &AmberTopology) -> Vec<usize> {
    let mut counts = vec![0usize; topo.atoms.len()];
    for bond in &topo.bonds {
        if bond.i < counts.len() {
            counts[bond.i] += 1;
        }
        if bond.j < counts.len() {
            counts[bond.j] += 1;
        }
    }
    counts
}
/// Return sorted list of (atom_idx, bond_count) pairs for atoms with ≥ min_bonds bonds.
pub fn highly_connected_atoms(topo: &AmberTopology, min_bonds: usize) -> Vec<(usize, usize)> {
    atom_bond_counts(topo)
        .into_iter()
        .enumerate()
        .filter(|(_, count)| *count >= min_bonds)
        .collect()
}
/// Compute the Coulomb energy (kcal mol⁻¹) between two point charges.
///
/// Uses the AMBER convention: `E = 332.0636 * q_i * q_j / r` where `r` is
/// in Ångströms and charges are in elementary charge units.
pub fn coulomb_energy(q_i: f64, q_j: f64, r: f64) -> f64 {
    pub(super) const COULOMB_FACTOR: f64 = 332.0636;
    if r < 1e-10 {
        return 0.0;
    }
    COULOMB_FACTOR * q_i * q_j / r
}
/// Compute the Lennard-Jones 12-6 energy between two atoms.
///
/// `e_ij = A / r^12 - B / r^6`  where A and B come from combining rules.
pub fn lj_energy(a_coeff: f64, b_coeff: f64, r: f64) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    let r6 = r.powi(6);
    let r12 = r6 * r6;
    a_coeff / r12 - b_coeff / r6
}
/// Lorentz-Berthelot combining rules: compute A and B for a pair (i, j) from
/// individual sigma_i, sigma_j, eps_i, eps_j values.
///
/// Returns `(a_coeff, b_coeff)`.
pub fn lorentz_berthelot(sigma_i: f64, sigma_j: f64, eps_i: f64, eps_j: f64) -> (f64, f64) {
    let sigma_ij = (sigma_i + sigma_j) * 0.5;
    let eps_ij = (eps_i * eps_j).sqrt();
    let s6 = sigma_ij.powi(6);
    let s12 = s6 * s6;
    let b_coeff = 4.0 * eps_ij * s6;
    let a_coeff = 4.0 * eps_ij * s12;
    (a_coeff, b_coeff)
}
/// Compute the harmonic bond potential energy for a single bond.
///
/// `E = k * (r - r0)^2`
pub fn bond_energy(k: f64, r0: f64, r: f64) -> f64 {
    k * (r - r0).powi(2)
}
/// Compute the harmonic angle potential energy.
///
/// `E = k * (theta - theta0)^2`
pub fn angle_energy(k: f64, theta0: f64, theta: f64) -> f64 {
    k * (theta - theta0).powi(2)
}
/// Compute the AMBER dihedral (torsion) energy.
///
/// `E = V_n/2 * (1 + cos(n*phi - gamma))`
pub fn dihedral_energy(v_n: f64, n: f64, phi: f64, gamma: f64) -> f64 {
    0.5 * v_n * (1.0 + (n * phi - gamma).cos())
}
/// Compute Euclidean distance between two 3D points (Å or any consistent unit).
pub fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Compute minimum-image distance between two points in an orthorhombic box.
pub fn min_image_distance(a: [f64; 3], b: [f64; 3], box_a: f64, box_b: f64, box_c: f64) -> f64 {
    let dx = wrap(a[0] - b[0], box_a);
    let dy = wrap(a[1] - b[1], box_b);
    let dz = wrap(a[2] - b[2], box_c);
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Wrap a displacement into `\[-L/2, L/2)`.
pub(super) fn wrap(d: f64, l: f64) -> f64 {
    if l < 1e-30 {
        return d;
    }
    let mut v = d;
    while v > l * 0.5 {
        v -= l;
    }
    while v < -l * 0.5 {
        v += l;
    }
    v
}
/// Compute the angle (radians) A–B–C.
pub fn compute_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
    let mag_ba = (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]).sqrt();
    let mag_bc = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]).sqrt();
    if mag_ba < 1e-10 || mag_bc < 1e-10 {
        return 0.0;
    }
    let cos_theta = (dot / (mag_ba * mag_bc)).clamp(-1.0, 1.0);
    cos_theta.acos()
}
/// Compute the dihedral angle (radians) A–B–C–D.
pub fn compute_dihedral(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b1 = sub3(b, a);
    let b2 = sub3(c, b);
    let b3 = sub3(d, c);
    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let mag1 = mag3(n1);
    let mag2 = mag3(n2);
    if mag1 < 1e-10 || mag2 < 1e-10 {
        return 0.0;
    }
    let cos_phi = dot3(n1, n2) / (mag1 * mag2);
    let cos_phi = cos_phi.clamp(-1.0, 1.0);
    let sign = dot3(cross3(n1, n2), b2);
    let phi = cos_phi.acos();
    if sign < 0.0 { -phi } else { phi }
}
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn mag3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
/// Generate a geometric temperature ladder for AMBER REMD simulations.
///
/// Returns `n_replicas` temperatures between `t_low` and `t_high` (K),
/// spaced geometrically (equal exchange acceptance rate approximation).
pub fn remd_temperature_ladder(t_low: f64, t_high: f64, n_replicas: usize) -> Vec<f64> {
    if n_replicas == 0 {
        return Vec::new();
    }
    if n_replicas == 1 {
        return vec![t_low];
    }
    let ratio = (t_high / t_low).powf(1.0 / (n_replicas - 1) as f64);
    (0..n_replicas)
        .map(|i| t_low * ratio.powi(i as i32))
        .collect()
}
/// Compute the (approximate) REMD exchange probability between two replicas.
///
/// `P = min(1, exp(-delta))` where
/// `delta = (1/kT_i - 1/kT_j) * (E_j - E_i)`.
///
/// `k_b` = Boltzmann constant (kcal mol⁻¹ K⁻¹) ≈ 0.001987207.
pub fn remd_exchange_probability(e_i: f64, e_j: f64, t_i: f64, t_j: f64) -> f64 {
    pub(super) const K_B: f64 = 0.001987207;
    if t_i < 1e-10 || t_j < 1e-10 {
        return 0.0;
    }
    let beta_i = 1.0 / (K_B * t_i);
    let beta_j = 1.0 / (K_B * t_j);
    let delta = (beta_j - beta_i) * (e_i - e_j);
    let exponent = -delta;
    if exponent >= 0.0 { 1.0 } else { exponent.exp() }
}
/// Parse a minimal AMBER NMR distance restraint file (RST format).
///
/// Each restraint block looks like:
/// ```text
///  &rst  iat=1,2, r1=1.5,r2=2.0,r3=3.0,r4=4.0, rk2=2.0,rk3=2.0 /
/// ```
pub fn parse_restraints(s: &str) -> Vec<AmberRestraint> {
    let mut restraints = Vec::new();
    let mut in_rst = false;
    let mut block = String::new();
    for line in s.lines() {
        let t = line.trim();
        if t.starts_with("&rst") || t.starts_with("&RST") {
            in_rst = true;
            block.clear();
            block.push_str(t);
            block.push(' ');
            if t.contains('/') {
                in_rst = false;
                if let Some(rst) = parse_rst_block(&block) {
                    restraints.push(rst);
                }
                block.clear();
            }
            continue;
        }
        if in_rst {
            block.push_str(t);
            block.push(' ');
            if t.ends_with('/') || t.contains('/') {
                in_rst = false;
                if let Some(rst) = parse_rst_block(&block) {
                    restraints.push(rst);
                }
                block.clear();
            }
        }
    }
    restraints
}
pub(super) fn parse_rst_block(block: &str) -> Option<AmberRestraint> {
    fn get_val(block: &str, key: &str) -> Option<f64> {
        let klen = key.len();
        let pos = block.find(key)?;
        let rest = block[pos + klen..].trim_start_matches(|c: char| c == '=' || c.is_whitespace());
        let end = rest
            .find(|c: char| c == ',' || c == '/' || c.is_whitespace())
            .unwrap_or(rest.len());
        rest[..end].parse().ok()
    }
    fn get_int_pair(block: &str, key: &str) -> Option<(usize, usize)> {
        let klen = key.len();
        let pos = block.find(key)?;
        let rest = block[pos + klen..].trim_start_matches(|c: char| c == '=' || c.is_whitespace());
        let pair_str: String = rest
            .chars()
            .take_while(|c| c.is_ascii_digit() || *c == ',')
            .collect();
        let parts: Vec<&str> = pair_str.split(',').filter(|s| !s.is_empty()).collect();
        if parts.len() >= 2 {
            let a: usize = parts[0].trim().parse().ok()?;
            let b: usize = parts[1].trim().parse().ok()?;
            Some((a.saturating_sub(1), b.saturating_sub(1)))
        } else if parts.len() == 1 {
            let a: usize = parts[0].trim().parse().ok()?;
            Some((a.saturating_sub(1), 0))
        } else {
            None
        }
    }
    let (iat1, iat2) = get_int_pair(block, "iat")?;
    let r1 = get_val(block, "r1").unwrap_or(0.0);
    let r2 = get_val(block, "r2").unwrap_or(0.0);
    let r3 = get_val(block, "r3").unwrap_or(0.0);
    let r4 = get_val(block, "r4").unwrap_or(0.0);
    let rk2 = get_val(block, "rk2").unwrap_or(1.0);
    let rk3 = get_val(block, "rk3").unwrap_or(1.0);
    Some(AmberRestraint {
        iat1,
        iat2,
        r1,
        r2,
        r3,
        r4,
        rk2,
        rk3,
    })
}
/// List all `%FLAG` section names found in a prmtop string.
pub fn list_prmtop_flags(prmtop: &str) -> Vec<String> {
    prmtop
        .lines()
        .filter_map(|line| {
            line.trim()
                .strip_prefix("%FLAG ")
                .map(|rest| rest.trim().to_string())
        })
        .collect()
}
/// Check whether a given `%FLAG` section is present in a prmtop string.
pub fn has_prmtop_flag(prmtop: &str, flag: &str) -> bool {
    let target = format!("%FLAG {}", flag.to_uppercase());
    prmtop.lines().any(|l| l.trim() == target.as_str())
}
/// Compute the centre of mass (CoM) of all atoms in a topology given their
/// positions from an `AmberCoordinates` record.
pub fn centre_of_mass(topo: &AmberTopology, coords: &AmberCoordinates) -> [f64; 3] {
    let mut total_mass = 0.0_f64;
    let mut com = [0.0_f64; 3];
    for (i, atom) in topo.atoms.iter().enumerate() {
        if i >= coords.n_atoms {
            break;
        }
        let pos = coords.position(i);
        let m = atom.mass;
        com[0] += m * pos[0];
        com[1] += m * pos[1];
        com[2] += m * pos[2];
        total_mass += m;
    }
    if total_mass > 1e-30 {
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
    }
    com
}
/// Compute the radius of gyration (Å) of atoms in a topology.
pub fn radius_of_gyration(topo: &AmberTopology, coords: &AmberCoordinates) -> f64 {
    let com = centre_of_mass(topo, coords);
    let mut total_mass = 0.0_f64;
    let mut sum = 0.0_f64;
    for (i, atom) in topo.atoms.iter().enumerate() {
        if i >= coords.n_atoms {
            break;
        }
        let pos = coords.position(i);
        let m = atom.mass;
        let d2 = (pos[0] - com[0]).powi(2) + (pos[1] - com[1]).powi(2) + (pos[2] - com[2]).powi(2);
        sum += m * d2;
        total_mass += m;
    }
    if total_mass > 1e-30 {
        (sum / total_mass).sqrt()
    } else {
        0.0
    }
}
#[cfg(test)]
mod tests_amber_geo {
    use super::*;
    use crate::amber::AmberAtom;
    use crate::amber::types::*;
    use std::f64::consts::PI;
    use std::str::FromStr;
    #[test]
    fn test_validate_topology_clean() {
        let topo = AmberTopology {
            title: "valid".into(),
            atoms: vec![
                AmberAtom {
                    name: "C".into(),
                    residue_name: "ALA".into(),
                    charge: 0.0,
                    mass: 12.0,
                    atom_type: "CT".into(),
                },
                AmberAtom {
                    name: "N".into(),
                    residue_name: "ALA".into(),
                    charge: 0.0,
                    mass: 14.0,
                    atom_type: "N".into(),
                },
            ],
            bonds: vec![AmberBond {
                i: 0,
                j: 1,
                k: 300.0,
                r0: 1.5,
            }],
            angles: vec![],
            n_atoms: 2,
            n_bonds: 1,
        };
        let issues = validate_topology(&topo);
        assert!(issues.is_empty(), "unexpected issues: {:?}", issues);
    }
    #[test]
    fn test_validate_topology_bad_bond_index() {
        let topo = AmberTopology {
            title: "bad".into(),
            atoms: vec![AmberAtom {
                name: "C".into(),
                residue_name: "R".into(),
                charge: 0.0,
                mass: 12.0,
                atom_type: "C".into(),
            }],
            bonds: vec![AmberBond {
                i: 0,
                j: 99,
                k: 300.0,
                r0: 1.5,
            }],
            angles: vec![],
            n_atoms: 1,
            n_bonds: 1,
        };
        let issues = validate_topology(&topo);
        assert!(!issues.is_empty(), "should flag out-of-range atom index");
    }
    #[test]
    fn test_validate_topology_negative_mass() {
        let topo = AmberTopology {
            title: "neg".into(),
            atoms: vec![AmberAtom {
                name: "X".into(),
                residue_name: "R".into(),
                charge: 0.0,
                mass: -1.0,
                atom_type: "X".into(),
            }],
            bonds: vec![],
            angles: vec![],
            n_atoms: 1,
            n_bonds: 0,
        };
        let issues = validate_topology(&topo);
        assert!(!issues.is_empty(), "should flag negative mass");
    }
    #[test]
    fn test_atom_bond_counts_basic() {
        let topo = AmberTopology {
            title: "t".into(),
            atoms: vec![
                AmberAtom {
                    name: "A".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "A".into(),
                },
                AmberAtom {
                    name: "B".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "B".into(),
                },
                AmberAtom {
                    name: "C".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "C".into(),
                },
            ],
            bonds: vec![
                AmberBond {
                    i: 0,
                    j: 1,
                    k: 1.0,
                    r0: 1.0,
                },
                AmberBond {
                    i: 1,
                    j: 2,
                    k: 1.0,
                    r0: 1.0,
                },
            ],
            angles: vec![],
            n_atoms: 3,
            n_bonds: 2,
        };
        let counts = atom_bond_counts(&topo);
        assert_eq!(counts[0], 1);
        assert_eq!(counts[1], 2);
        assert_eq!(counts[2], 1);
    }
    #[test]
    fn test_coulomb_energy_unit_charges() {
        let e = coulomb_energy(1.0, 1.0, 1.0);
        assert!((e - 332.0636).abs() < 1e-3, "e={}", e);
    }
    #[test]
    fn test_coulomb_energy_opposite_charges() {
        let e = coulomb_energy(1.0, -1.0, 1.0);
        assert!(e < 0.0, "opposite charges should give negative energy");
    }
    #[test]
    fn test_lj_energy_at_r_min() {
        let a = 4.0_f64;
        let b = 4.0_f64;
        let r_min = (2.0 * a / b).powf(1.0 / 6.0);
        let e = lj_energy(a, b, r_min);
        let expected = -b * b / (4.0 * a);
        assert!(
            (e - expected).abs() < 1e-8,
            "e={}, expected={}",
            e,
            expected
        );
    }
    #[test]
    fn test_lorentz_berthelot_same_atoms() {
        let sigma = 3.0;
        let eps = 0.1;
        let (a, b) = lorentz_berthelot(sigma, sigma, eps, eps);
        let s6 = sigma.powi(6);
        let s12 = s6 * s6;
        let expected_b = 4.0 * eps * s6;
        let expected_a = 4.0 * eps * s12;
        assert!((a - expected_a).abs() < 1e-8, "a={}", a);
        assert!((b - expected_b).abs() < 1e-8, "b={}", b);
    }
    #[test]
    fn test_bond_energy_at_equilibrium() {
        assert!((bond_energy(300.0, 1.5, 1.5)).abs() < 1e-12);
    }
    #[test]
    fn test_bond_energy_displaced() {
        let e = bond_energy(300.0, 1.5, 1.6);
        assert!((e - 3.0).abs() < 1e-8, "e={}", e);
    }
    #[test]
    fn test_angle_energy_at_equilibrium() {
        let theta0 = 109.5 * PI / 180.0;
        assert!((angle_energy(50.0, theta0, theta0)).abs() < 1e-12);
    }
    #[test]
    fn test_dihedral_energy_at_zero_phase() {
        let e = dihedral_energy(2.0, 1.0, 0.0, 0.0);
        assert!((e - 2.0).abs() < 1e-8, "e={}", e);
    }
    #[test]
    fn test_distance_basic() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        assert!((distance(a, b) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_distance_3d() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 6.0, 3.0];
        let d = distance(a, b);
        assert!((d - 5.0).abs() < 1e-10, "d={}", d);
    }
    #[test]
    fn test_min_image_distance_no_wrap() {
        let a = [1.0, 0.0, 0.0];
        let b = [3.0, 0.0, 0.0];
        let d = min_image_distance(a, b, 100.0, 100.0, 100.0);
        assert!((d - 2.0).abs() < 1e-10, "d={}", d);
    }
    #[test]
    fn test_min_image_distance_wrap() {
        let a = [0.5, 0.0, 0.0];
        let b = [9.5, 0.0, 0.0];
        let d = min_image_distance(a, b, 10.0, 10.0, 10.0);
        assert!((d - 1.0).abs() < 1e-8, "d={}", d);
    }
    #[test]
    fn test_compute_angle_linear() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0];
        let theta = compute_angle(a, b, c);
        assert!((theta - PI).abs() < 1e-8, "theta={}", theta);
    }
    #[test]
    fn test_compute_angle_right() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let theta = compute_angle(a, b, c);
        assert!((theta - PI / 2.0).abs() < 1e-8, "theta={}", theta);
    }
    #[test]
    fn test_remd_temperature_ladder_two_replicas() {
        let temps = remd_temperature_ladder(300.0, 600.0, 2);
        assert_eq!(temps.len(), 2);
        assert!((temps[0] - 300.0).abs() < 1e-6);
        assert!((temps[1] - 600.0).abs() < 1e-6);
    }
    #[test]
    fn test_remd_temperature_ladder_geometric() {
        let temps = remd_temperature_ladder(300.0, 1200.0, 3);
        assert_eq!(temps.len(), 3);
        assert!((temps[1] / temps[0] - temps[2] / temps[1]).abs() < 1e-6);
    }
    #[test]
    fn test_remd_exchange_probability_identical() {
        let p = remd_exchange_probability(-100.0, -100.0, 300.0, 300.0);
        assert!((p - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_remd_exchange_probability_unfavorable() {
        let p = remd_exchange_probability(0.0, -1000.0, 300.0, 310.0);
        assert!((0.0..=1.0).contains(&p), "p={}", p);
    }
    #[test]
    fn test_parse_restraints_basic() {
        let rst_text = " &rst  iat=1,2, r1=1.5,r2=2.0,r3=3.0,r4=4.0, rk2=2.0,rk3=2.0 /\n";
        let rests = parse_restraints(rst_text);
        assert_eq!(rests.len(), 1);
        assert_eq!(rests[0].iat1, 0);
        assert_eq!(rests[0].iat2, 1);
        assert!((rests[0].r2 - 2.0).abs() < 1e-6);
        assert!((rests[0].rk2 - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_restraint_energy_flat_bottom() {
        let r = AmberRestraint {
            iat1: 0,
            iat2: 1,
            r1: 1.5,
            r2: 2.0,
            r3: 3.0,
            r4: 4.0,
            rk2: 1.0,
            rk3: 1.0,
        };
        assert!((r.energy(2.5)).abs() < 1e-10);
        let e_out = r.energy(3.5);
        assert!((e_out - 0.25).abs() < 1e-8, "e={}", e_out);
    }
    #[test]
    fn test_list_prmtop_flags() {
        let prmtop = "%FLAG TITLE\n%FORMAT(20a4)\nhello\n%FLAG POINTERS\n%FORMAT(10I8)\n1\n";
        let flags = list_prmtop_flags(prmtop);
        assert!(flags.contains(&"TITLE".to_string()));
        assert!(flags.contains(&"POINTERS".to_string()));
    }
    #[test]
    fn test_has_prmtop_flag_true() {
        let prmtop = "%FLAG CHARGE\n%FORMAT(5E16.8)\n1.0\n";
        assert!(has_prmtop_flag(prmtop, "CHARGE"));
    }
    #[test]
    fn test_has_prmtop_flag_false() {
        let prmtop = "%FLAG CHARGE\n%FORMAT(5E16.8)\n1.0\n";
        assert!(!has_prmtop_flag(prmtop, "MISSING_SECTION"));
    }
    #[test]
    fn test_centre_of_mass_equal_masses() {
        let topo = AmberTopology {
            title: "t".into(),
            atoms: vec![
                AmberAtom {
                    name: "A".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "A".into(),
                },
                AmberAtom {
                    name: "B".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "B".into(),
                },
            ],
            bonds: vec![],
            angles: vec![],
            n_atoms: 2,
            n_bonds: 0,
        };
        let coords_str = "test\n    2\n   0.0000000   0.0000000   0.0000000   2.0000000   0.0000000   0.0000000\n";
        let coords = AmberCoordinates::from_str(coords_str).unwrap();
        let com = centre_of_mass(&topo, &coords);
        assert!((com[0] - 1.0).abs() < 1e-8, "com_x={}", com[0]);
        assert!((com[1]).abs() < 1e-8);
        assert!((com[2]).abs() < 1e-8);
    }
    #[test]
    fn test_radius_of_gyration_dumbbell() {
        let topo = AmberTopology {
            title: "t".into(),
            atoms: vec![
                AmberAtom {
                    name: "A".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "A".into(),
                },
                AmberAtom {
                    name: "B".into(),
                    residue_name: "R".into(),
                    charge: 0.0,
                    mass: 1.0,
                    atom_type: "B".into(),
                },
            ],
            bonds: vec![],
            angles: vec![],
            n_atoms: 2,
            n_bonds: 0,
        };
        let coords_str = "test\n    2\n  -1.0000000   0.0000000   0.0000000   1.0000000   0.0000000   0.0000000\n";
        let coords = AmberCoordinates::from_str(coords_str).unwrap();
        let rg = radius_of_gyration(&topo, &coords);
        assert!((rg - 1.0).abs() < 1e-8, "rg={}", rg);
    }
}
