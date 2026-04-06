//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{AmberLjParams, DsspLabel, GyrRadius, ProteinTopology, SecondaryStructure};

/// Compute phi backbone dihedral angle (C_prev, N, CA, C) in degrees.
pub fn phi_angle(c_prev: [f64; 3], n: [f64; 3], ca: [f64; 3], c: [f64; 3]) -> f64 {
    dihedral_angle(c_prev, n, ca, c)
}
/// Compute psi backbone dihedral angle (N, CA, C, N_next) in degrees.
pub fn psi_angle(n: [f64; 3], ca: [f64; 3], c: [f64; 3], n_next: [f64; 3]) -> f64 {
    dihedral_angle(n, ca, c, n_next)
}
/// Compute omega backbone dihedral angle (CA, C, N_next, CA_next) in degrees.
pub fn omega_angle(ca: [f64; 3], c: [f64; 3], n_next: [f64; 3], ca_next: [f64; 3]) -> f64 {
    dihedral_angle(ca, c, n_next, ca_next)
}
/// Compute H-bond score from distance and angle.
pub fn hbond_score(dist: f64, angle_deg: f64) -> f64 {
    let r_factor = (2.5 / dist.max(1e-5)).powi(6);
    let a_factor = ((angle_deg - 150.0) / 30.0).powf(2.0);
    r_factor * (-a_factor).exp()
}
/// Return DSSP label name as a char.
pub fn dssp_label(label: DsspLabel) -> char {
    match label {
        DsspLabel::AlphaHelix => 'H',
        DsspLabel::Helix310 => 'G',
        DsspLabel::PiHelix => 'I',
        DsspLabel::BetaSheet => 'E',
        DsspLabel::BetaBridge => 'B',
        DsspLabel::Bend => 'S',
        DsspLabel::Turn => 'T',
        DsspLabel::Coil => 'C',
    }
}
/// Compute dihedral angle between four points in degrees.
pub fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let b1 = sub3(p2, p1);
    let b2 = sub3(p3, p2);
    let b3 = sub3(p4, p3);
    let n1 = cross3(b1, b2);
    let n2 = cross3(b2, b3);
    let m1 = cross3(n1, b2);
    let y = dot3(m1, n2) / (mag3(m1).max(1e-30) * mag3(n2).max(1e-30));
    let x = dot3(n1, n2) / (mag3(n1).max(1e-30) * mag3(n2).max(1e-30));
    y.atan2(x).to_degrees()
}
/// H-bond angle D-H...A in degrees.
pub fn hbond_angle(d: [f64; 3], h: [f64; 3], a: [f64; 3]) -> f64 {
    let dh = sub3(h, d);
    let ha = sub3(a, h);
    let cos_angle = dot3(dh, ha) / (mag3(dh).max(1e-30) * mag3(ha).max(1e-30));
    cos_angle.clamp(-1.0, 1.0).acos().to_degrees()
}
/// 3D distance.
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}
/// 3D subtraction.
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// 3D dot product.
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// 3D cross product.
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// 3D magnitude.
pub(super) fn mag3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
/// Compute AMBER 12-6 Lennard-Jones energy between two atoms.
///
/// Uses the Lorentz-Berthelot combining rules:
/// - epsilon_ij = sqrt(eps_i * eps_j)
/// - Rmin_ij = Rmin2_i + Rmin2_j
///
/// E_LJ = eps_ij * \[(Rmin_ij/r)^12 - 2*(Rmin_ij/r)^6\]
pub fn amber_lj_energy(params_i: &AmberLjParams, params_j: &AmberLjParams, r: f64) -> f64 {
    if r < 1e-6 {
        return 1e30;
    }
    let eps = (params_i.epsilon * params_j.epsilon).sqrt();
    let rmin = params_i.rmin2 + params_j.rmin2;
    let ratio = rmin / r;
    let r6 = ratio.powi(6);
    eps * (r6 * r6 - 2.0 * r6)
}
/// Compute the 1-4 scaled AMBER LJ energy (scale factor 0.5 for torsion pairs).
///
/// In AMBER, 1-4 van der Waals interactions are scaled by 0.5 to avoid
/// double-counting with the dihedral term.
pub fn amber_lj_14_energy(params_i: &AmberLjParams, params_j: &AmberLjParams, r: f64) -> f64 {
    0.5 * amber_lj_energy(params_i, params_j, r)
}
/// Compute phi and psi backbone dihedral angles for all residues in a topology.
///
/// Returns a vector of (phi, psi) pairs in degrees. Residues at the termini
/// that cannot have both angles yield `None`.
pub fn compute_ramachandran(topo: &ProteinTopology) -> Vec<Option<(f64, f64)>> {
    let n_res = topo.n_residues();
    if n_res < 2 {
        return vec![None; n_res];
    }
    let mut result = Vec::with_capacity(n_res);
    for res in 0..n_res {
        let base = res * 4;
        if base + 3 >= topo.n_atoms {
            result.push(None);
            continue;
        }
        let phi_opt = if res > 0 {
            let prev_c = topo.positions[(res - 1) * 4 + 2];
            let n_pos = topo.positions[base];
            let ca_pos = topo.positions[base + 1];
            let c_pos = topo.positions[base + 2];
            Some(dihedral_angle(prev_c, n_pos, ca_pos, c_pos))
        } else {
            None
        };
        let psi_opt = if res + 1 < n_res && (res + 1) * 4 < topo.n_atoms {
            let n_pos = topo.positions[base];
            let ca_pos = topo.positions[base + 1];
            let c_pos = topo.positions[base + 2];
            let next_n = topo.positions[(res + 1) * 4];
            Some(dihedral_angle(n_pos, ca_pos, c_pos, next_n))
        } else {
            None
        };
        match (phi_opt, psi_opt) {
            (Some(phi), Some(psi)) => result.push(Some((phi, psi))),
            _ => result.push(None),
        }
    }
    result
}
/// Detect alpha-helical residues from backbone H-bond geometry.
///
/// A residue `i` is considered helical if the O(i) ... N(i+4) distance < 3.5 Å
/// and the N(i+4)-H ... O(i) angle > 120° (simplified).
pub fn detect_helix_residues(
    n_pos: &[[f64; 3]],
    _ca_pos: &[[f64; 3]],
    _c_pos: &[[f64; 3]],
    o_pos: &[[f64; 3]],
) -> Vec<bool> {
    let n = n_pos.len().min(o_pos.len());
    let mut helix = vec![false; n];
    for i in 0..n {
        if i + 4 < n {
            let d = dist3(o_pos[i], n_pos[i + 4]);
            if d < 3.5 {
                helix[i] = true;
            }
        }
    }
    helix
}
/// Detect beta-sheet residues from backbone H-bond geometry.
///
/// A residue `i` is flagged as beta if it has an O ... N distance < 3.5 Å
/// to any non-adjacent residue |j - i| >= 3 (excluding helix).
pub fn detect_sheet_residues(
    n_pos: &[[f64; 3]],
    _ca_pos: &[[f64; 3]],
    o_pos: &[[f64; 3]],
    min_sep: usize,
) -> Vec<bool> {
    let n = n_pos.len().min(o_pos.len());
    let mut sheet = vec![false; n];
    for i in 0..n {
        for j in 0..n {
            if (i as isize - j as isize).unsigned_abs() < min_sep {
                continue;
            }
            let d = dist3(o_pos[i], n_pos[j]);
            if d < 3.5 {
                sheet[i] = true;
                break;
            }
        }
    }
    sheet
}
/// Compute the radius of gyration Rg of a protein in Angstroms.
///
/// Rg = sqrt( sum_i m_i |r_i - r_com|^2 / sum_i m_i )
pub fn radius_of_gyration(positions: &[[f64; 3]], masses: &[f64]) -> f64 {
    let gr = GyrRadius::new(positions.to_vec(), masses.to_vec());
    gr.rg_squared().sqrt()
}
/// Compute solvent-accessible surface area (SASA) using a simplified
/// Lee-Richards sphere-overlap approach.
///
/// Each atom contributes (4π r_probe²) minus the fraction occluded by neighbours.
/// `probe_radius` is typically 1.4 Å (water probe).
///
/// Returns per-atom SASA in Å².
pub fn compute_sasa(positions: &[[f64; 3]], radii: &[f64], probe_radius: f64) -> Vec<f64> {
    let n = positions.len().min(radii.len());
    let mut sasa = Vec::with_capacity(n);
    for i in 0..n {
        let ri = radii[i] + probe_radius;
        let area_i = 4.0 * std::f64::consts::PI * ri * ri;
        let mut overlap_frac = 0.0_f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let rj = radii[j] + probe_radius;
            let d = dist3(positions[i], positions[j]);
            if d < ri + rj && d > (ri - rj).abs() {
                let cos_alpha = (ri * ri + d * d - rj * rj) / (2.0 * ri * d).max(1e-20);
                let frac = 0.5 * (1.0 - cos_alpha.clamp(-1.0, 1.0));
                overlap_frac += frac;
            }
        }
        sasa.push(area_i * (1.0 - overlap_frac.min(1.0)));
    }
    sasa
}
/// Total SASA (sum over all atoms) in Å².
pub fn total_sasa(positions: &[[f64; 3]], radii: &[f64], probe_radius: f64) -> f64 {
    compute_sasa(positions, radii, probe_radius).iter().sum()
}
/// Compute hydrogen bond energy using the Kabsch-Sander electrostatic model.
///
/// E_HB = f * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332 kcal/mol
///
/// where f = 0.084, and the distances are in Å.
/// Returns energy in kcal/mol. Hydrogen bonds have E < -0.5 kcal/mol.
pub fn kabsch_sander_hbond(
    o_pos: [f64; 3],
    c_pos: [f64; 3],
    n_pos: [f64; 3],
    h_pos: [f64; 3],
) -> f64 {
    let r_on = dist3(o_pos, n_pos).max(0.01);
    let r_ch = dist3(c_pos, h_pos).max(0.01);
    let r_oh = dist3(o_pos, h_pos).max(0.01);
    let r_cn = dist3(c_pos, n_pos).max(0.01);
    0.084 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn) * 332.0
}
/// Detect all Kabsch-Sander hydrogen bonds in a secondary structure object.
///
/// Returns list of (donor_idx, acceptor_idx, energy) for bonds with E < cutoff.
pub fn find_kabsch_sander_hbonds(
    ss: &SecondaryStructure,
    energy_cutoff: f64,
) -> Vec<(usize, usize, f64)> {
    let n = ss.n_pos.len();
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let h_approx = ss.n_pos[j];
            let e = kabsch_sander_hbond(ss.o_pos[i], ss.c_pos[i], ss.n_pos[j], h_approx);
            if e < energy_cutoff {
                bonds.push((i, j, e));
            }
        }
    }
    bonds
}
/// Compute the folding order parameter Q based on fraction of native contacts.
///
/// Q = (1/Q0) * sum_{i<j, |i-j|>1} exp(-((r_ij - r_ij^0) / (2*sigma))^2)
///
/// where r_ij are current distances, r_ij^0 native distances, sigma = 1 Å.
pub fn q_fold_order_parameter(
    positions: &[[f64; 3]],
    native_positions: &[[f64; 3]],
    sigma: f64,
    min_sep: usize,
) -> f64 {
    let n = positions.len().min(native_positions.len());
    if n < 2 {
        return 1.0;
    }
    let mut sum = 0.0;
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + min_sep)..n {
            let r_curr = dist3(positions[i], positions[j]);
            let r_nat = dist3(native_positions[i], native_positions[j]);
            let x = (r_curr - r_nat) / (2.0 * sigma);
            sum += (-x * x).exp();
            count += 1;
        }
    }
    if count == 0 {
        return 1.0;
    }
    sum / count as f64
}
/// Detect disulfide bonds from sulfur atom positions.
///
/// A disulfide bond S-S is detected when d(S_i, S_j) < cutoff_angstrom (typically 2.2 Å).
/// Returns list of (atom_i, atom_j, distance) pairs.
pub fn detect_disulfide_bonds(
    s_positions: &[[f64; 3]],
    s_atom_indices: &[usize],
    cutoff_angstrom: f64,
) -> Vec<(usize, usize, f64)> {
    let n = s_positions.len();
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist3(s_positions[i], s_positions[j]);
            if d < cutoff_angstrom {
                let idx_i = if i < s_atom_indices.len() {
                    s_atom_indices[i]
                } else {
                    i
                };
                let idx_j = if j < s_atom_indices.len() {
                    s_atom_indices[j]
                } else {
                    j
                };
                bonds.push((idx_i, idx_j, d));
            }
        }
    }
    bonds
}
/// Disulfide bond detector operating directly on ProteinTopology.
///
/// Scans all atoms named "SG" (cysteine sulfur) and checks for S-S contacts.
pub fn find_disulfide_in_topology(
    topo: &ProteinTopology,
    cutoff_angstrom: f64,
) -> Vec<(usize, usize, f64)> {
    let sg_indices: Vec<usize> = topo
        .atom_names
        .iter()
        .enumerate()
        .filter(|(_, name)| name.as_str() == "SG")
        .map(|(i, _)| i)
        .collect();
    let sg_positions: Vec<[f64; 3]> = sg_indices.iter().map(|&i| topo.positions[i]).collect();
    detect_disulfide_bonds(&sg_positions, &sg_indices, cutoff_angstrom)
}
/// Identify salt bridges from positions of charged groups.
///
/// A salt bridge is a close contact between a positive and negative
/// charged group with centre-to-centre distance < cutoff (typically 4.0 Å).
///
/// Returns list of (pos_idx, neg_idx, distance) pairs.
pub fn identify_salt_bridges(
    pos_centers: &[[f64; 3]],
    neg_centers: &[[f64; 3]],
    cutoff_angstrom: f64,
) -> Vec<(usize, usize, f64)> {
    let mut bridges = Vec::new();
    for (i, &p) in pos_centers.iter().enumerate() {
        for (j, &n) in neg_centers.iter().enumerate() {
            let d = dist3(p, n);
            if d < cutoff_angstrom {
                bridges.push((i, j, d));
            }
        }
    }
    bridges
}
/// Salt bridge strength estimated from Coulomb energy (kcal/mol).
///
/// E_salt = 332 * q_pos * q_neg / (epsilon_r * r)
/// with q_pos = +1, q_neg = -1, epsilon_r = 4 (protein interior).
pub fn salt_bridge_energy(r: f64, epsilon_r: f64) -> f64 {
    if r < 0.1 {
        return -1e10;
    }
    -332.0 / (epsilon_r * r)
}
/// 3×3 matrix multiply A * B^T.
pub fn mat3_mul_t(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[j][k];
            }
        }
    }
    c
}
/// 3×3 matrix determinant.
pub fn mat3_det(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Jacobi eigendecomposition for symmetric 3×3 matrix.
///
/// Returns (eigenvalues, eigenvectors) after up to 50 sweeps.
pub fn jacobi3(a: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut d = a;
    let mut v = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let mut max_val = 0.0_f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if d[i][j].abs() > max_val {
                    max_val = d[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }
        let theta = 0.5 * (d[q][q] - d[p][p]) / d[p][q];
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            1.0 / (theta - (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let dpq = d[p][q];
        d[p][p] -= t * dpq;
        d[q][q] += t * dpq;
        d[p][q] = 0.0;
        d[q][p] = 0.0;
        for r in 0..3 {
            if r != p && r != q {
                let drp = d[r][p];
                let drq = d[r][q];
                d[r][p] = c * drp - s * drq;
                d[p][r] = d[r][p];
                d[r][q] = c * drq + s * drp;
                d[q][r] = d[r][q];
            }
        }
        for r in 0..3 {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = c * vrp - s * vrq;
            v[r][q] = c * vrq + s * vrp;
        }
    }
    ([d[0][0], d[1][1], d[2][2]], v)
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::biomolecular::ProteinElasticity;
    use crate::biomolecular::RgProtein;
    use crate::biomolecular::RmsdCalculation;
    use crate::biomolecular::SecStructLabel;
    use crate::biomolecular::SecondaryStructureStats;
    #[test]
    fn test_protein_topology_add_residue() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0, 0.0, 0.0],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        assert_eq!(topo.n_residues(), 1);
        assert_eq!(topo.n_atoms, 4);
    }
    #[test]
    fn test_protein_topology_two_residues() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0; 3],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        topo.add_residue(
            "GLY",
            [3.8, 0.0, 0.0],
            [5.2, 0.5, 0.0],
            [6.0, 1.5, 0.0],
            [6.5, 2.5, 0.0],
        );
        assert_eq!(topo.n_residues(), 2);
        assert!(topo.bonds.len() > 3);
    }
    #[test]
    fn test_ff_amino_lj_energy() {
        let ff = ForceFieldAmino::new_charmm36();
        let e = ff.lj_energy("N", "CA", 4.0);
        assert!(e.is_finite());
    }
    #[test]
    fn test_ff_amino_bond_energy_at_eq() {
        let ff = ForceFieldAmino::new_charmm36();
        let r0 = *ff.bond_lengths.get("N-CA").unwrap();
        let e = ff.bond_energy("N-CA", r0);
        assert!(e.abs() < 1e-10);
    }
    #[test]
    fn test_water_tip3p_oh_length() {
        let w = WaterModel::tip3p();
        assert!((w.oh_length - 0.9572).abs() < 1e-4);
    }
    #[test]
    fn test_water_tip4p_virtual() {
        let w = WaterModel::tip4p();
        assert!(w.d_virtual > 0.0);
    }
    #[test]
    fn test_water_spce_lj() {
        let w = WaterModel::spce();
        let e = w.lj_oo(3.5);
        assert!(e.is_finite());
    }
    #[test]
    fn test_water_dipole_moment() {
        let w = WaterModel::tip3p();
        let mu = w.dipole_moment();
        assert!(mu > 0.0);
    }
    #[test]
    fn test_water_h_positions() {
        let w = WaterModel::tip3p();
        let [h1, h2] = w.h_positions([0.0, 0.0, 0.0]);
        assert!((dist3(h1, [0.0, 0.0, 0.0]) - w.oh_length).abs() < 0.001);
        assert!((dist3(h2, [0.0, 0.0, 0.0]) - w.oh_length).abs() < 0.001);
    }
    #[test]
    fn test_protein_folding_q_factor() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0; 3],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        let mut fold = ProteinFolding::new(topo, 300.0, 8.0, 0.002);
        fold.step();
        let q = fold.contact_fraction_q();
        assert!((0.0..=1.0).contains(&q));
    }
    #[test]
    fn test_protein_folding_rmsd() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0; 3],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        let fold = ProteinFolding::new(topo, 300.0, 8.0, 0.002);
        let rmsd = fold.rmsd();
        assert!(rmsd >= 0.0);
    }
    #[test]
    fn test_ramachandran_helical() {
        let plot = RamachandranPlot::new(10.0);
        assert!(plot.is_helical(-60.0, -40.0));
        assert!(!plot.is_helical(60.0, 60.0));
    }
    #[test]
    fn test_ramachandran_add_point() {
        let mut plot = RamachandranPlot::new(10.0);
        plot.add_point(-60.0, -40.0);
        assert_eq!(plot.phi_angles.len(), 1);
    }
    #[test]
    fn test_ramachandran_helix_fraction() {
        let mut plot = RamachandranPlot::new(10.0);
        plot.add_point(-60.0, -40.0);
        plot.add_point(60.0, 60.0);
        let frac = plot.helix_fraction();
        assert!((frac - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_solvation_shell_accumulate() {
        let mut ss = SolvationShell::new(10.0, 0.1, 3.5);
        ss.solute_positions = vec![[0.0; 3]];
        ss.solvent_positions = vec![[2.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        ss.accumulate();
        assert!(ss.n_snapshots == 1);
    }
    #[test]
    fn test_gyr_radius_com() {
        let gr = GyrRadius::new(vec![[0.0; 3], [2.0, 0.0, 0.0]], vec![1.0, 1.0]);
        let com = gr.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_gyr_radius_rg_sq() {
        let gr = GyrRadius::new(vec![[0.0; 3], [2.0, 0.0, 0.0]], vec![1.0, 1.0]);
        let rg2 = gr.rg_squared();
        assert!((rg2 - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_hbond_analysis_detect() {
        let mut hb = HBondAnalysis::new(3.5, 150.0);
        hb.add_donor([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        hb.add_acceptor([2.5, 0.1, 0.0]);
        hb.detect();
        let _ = hb.n_hbonds();
    }
    #[test]
    fn test_dssp_creation() {
        let n = vec![[0.0; 3], [4.0, 0.0, 0.0]];
        let ca = vec![[1.46, 0.0, 0.0], [5.46, 0.0, 0.0]];
        let c = vec![[2.3, 1.2, 0.0], [6.3, 1.2, 0.0]];
        let o = vec![[2.3, 2.4, 0.0], [6.3, 2.4, 0.0]];
        let ss = SecondaryStructure::new(n, ca, c, o);
        assert_eq!(ss.labels.len(), 2);
    }
    #[test]
    fn test_binding_free_energy_mmpbsa() {
        let mut bfe = BindingFreeEnergy::new(300.0, 0.0054);
        bfe.add_frame(-100.0, -60.0, -30.0, -5.0);
        bfe.add_frame(-102.0, -61.0, -31.0, -5.5);
        let dg = bfe.mmpbsa();
        assert!(dg.is_finite());
    }
    #[test]
    fn test_phi_angle_planar() {
        let c_prev = [0.0, 0.0, 0.0];
        let n = [1.46, 0.0, 0.0];
        let ca = [2.0, 1.2, 0.0];
        let c = [3.0, 1.2, 0.0];
        let phi = phi_angle(c_prev, n, ca, c);
        assert!(phi.is_finite());
    }
    #[test]
    fn test_hbond_score_decreases_with_distance() {
        let s1 = hbond_score(2.5, 170.0);
        let s2 = hbond_score(3.0, 170.0);
        assert!(s1 > s2);
    }
    #[test]
    fn test_dssp_label_helix() {
        assert_eq!(dssp_label(DsspLabel::AlphaHelix), 'H');
        assert_eq!(dssp_label(DsspLabel::Coil), 'C');
    }
    #[test]
    fn test_ff_amino_angle_energy() {
        let ff = ForceFieldAmino::new_charmm36();
        let theta0 = *ff.angle_eq.get("N-CA-C").unwrap();
        let e = ff.angle_energy("N-CA-C", theta0);
        assert!(e.abs() < 1e-10);
    }
    #[test]
    fn test_amber_lj_minimum_at_rmin() {
        let pi = AmberLjParams::new("C", 0.086, 1.908);
        let pj = AmberLjParams::new("C", 0.086, 1.908);
        let r_min = pi.rmin2 + pj.rmin2;
        let e_min = amber_lj_energy(&pi, &pj, r_min);
        let eps = (pi.epsilon * pj.epsilon).sqrt();
        assert!(
            (e_min - (-eps)).abs() < 1e-8,
            "e_min={e_min} expected {}",
            -eps
        );
    }
    #[test]
    fn test_amber_lj_repulsive_at_short_range() {
        let pi = AmberLjParams::new("C", 0.086, 1.908);
        let pj = AmberLjParams::new("C", 0.086, 1.908);
        let e = amber_lj_energy(&pi, &pj, 1.0);
        assert!(e > 0.0, "energy should be repulsive at short range");
    }
    #[test]
    fn test_amber_lj_14_scale_factor() {
        let pi = AmberLjParams::new("N", 0.17, 1.824);
        let pj = AmberLjParams::new("O", 0.21, 1.661);
        let e_full = amber_lj_energy(&pi, &pj, 3.0);
        let e_14 = amber_lj_14_energy(&pi, &pj, 3.0);
        assert!((e_14 - 0.5 * e_full).abs() < 1e-12);
    }
    #[test]
    fn test_amber_lj_defaults_count() {
        let params = AmberLjParams::amber99sb_defaults();
        assert!(params.len() >= 7, "should have at least 7 atom types");
    }
    #[test]
    fn test_charmm_bond_energy_at_eq() {
        let bond = CharmmBond::new("C", "N", 370.0, 1.335);
        assert!(bond.energy(1.335).abs() < 1e-12);
    }
    #[test]
    fn test_charmm_bond_energy_displaced() {
        let bond = CharmmBond::new("CA", "C", 317.0, 1.522);
        let e = bond.energy(1.622);
        assert!((e - 317.0 * 0.01).abs() < 1e-8);
    }
    #[test]
    fn test_charmm_angle_energy_at_eq() {
        let ang = CharmmAngle::new("N", "CA", "C", 63.0, 111.2);
        assert!(ang.energy(111.2).abs() < 1e-10);
    }
    #[test]
    fn test_charmm_angle_energy_urey_bradley() {
        let mut ang = CharmmAngle::new("N", "CA", "C", 63.0, 111.2);
        ang.k_ub = 40.0;
        ang.s0 = 2.24;
        let e = ang.energy_with_ub(111.2, 2.34);
        assert!(e > 0.0, "UB term should add energy when displaced");
    }
    #[test]
    fn test_charmm_dihedral_zero_at_delta() {
        let dih = CharmmDihedral::new("C", "CA", "N", "C", 1.0, 1, 0.0);
        let e = dih.energy(180.0);
        assert!(
            e.abs() < 1e-8,
            "energy at phi=180 with delta=0, n=1 should be 0"
        );
    }
    #[test]
    fn test_charmm_dihedral_max_at_zero() {
        let dih = CharmmDihedral::new("C", "CA", "N", "C", 1.5, 1, 0.0);
        let e = dih.energy(0.0);
        assert!((e - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_charmm_improper_energy_at_eq() {
        let imp = CharmmImproper::new("C", "CA", "N", "O", 80.0, 0.0);
        assert!(imp.energy(0.0).abs() < 1e-12);
    }
    #[test]
    fn test_compute_ramachandran_single_residue() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0; 3],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        let angles = compute_ramachandran(&topo);
        assert_eq!(angles.len(), 1);
        assert!(angles[0].is_none());
    }
    #[test]
    fn test_compute_ramachandran_two_residues() {
        let mut topo = ProteinTopology::new();
        topo.add_residue(
            "ALA",
            [0.0; 3],
            [1.46, 0.0, 0.0],
            [2.3, 1.2, 0.0],
            [2.3, 2.4, 0.0],
        );
        topo.add_residue(
            "GLY",
            [3.8, 0.0, 0.0],
            [5.2, 0.5, 0.0],
            [6.0, 1.5, 0.0],
            [6.5, 2.5, 0.0],
        );
        let angles = compute_ramachandran(&topo);
        assert_eq!(angles.len(), 2);
    }
    #[test]
    fn test_detect_helix_residues_none() {
        let n = vec![[0.0; 3], [10.0, 0.0, 0.0]];
        let ca = vec![[1.0; 3], [11.0; 3]];
        let c = vec![[2.0; 3], [12.0; 3]];
        let o = vec![[2.5; 3], [12.5; 3]];
        let helix = detect_helix_residues(&n, &ca, &c, &o);
        assert_eq!(helix.len(), 2);
        assert!(!helix[0]);
    }
    #[test]
    fn test_detect_sheet_residues_empty() {
        let n: Vec<[f64; 3]> = vec![];
        let ca: Vec<[f64; 3]> = vec![];
        let o: Vec<[f64; 3]> = vec![];
        let sheet = detect_sheet_residues(&n, &ca, &o, 3);
        assert!(sheet.is_empty());
    }
    #[test]
    fn test_radius_of_gyration_two_atoms() {
        let positions = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let rg = radius_of_gyration(&positions, &masses);
        assert!((rg - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sasa_single_atom() {
        let positions = vec![[0.0; 3]];
        let radii = vec![1.7_f64];
        let probe = 1.4;
        let sasa = compute_sasa(&positions, &radii, probe);
        let r_eff = radii[0] + probe;
        let expected = 4.0 * std::f64::consts::PI * r_eff * r_eff;
        assert!(
            (sasa[0] - expected).abs() < 1e-8,
            "single atom SASA mismatch"
        );
    }
    #[test]
    fn test_total_sasa_positive() {
        let positions = vec![[0.0; 3], [5.0, 0.0, 0.0]];
        let radii = vec![1.7, 1.7];
        let s = total_sasa(&positions, &radii, 1.4);
        assert!(s > 0.0);
    }
    #[test]
    fn test_hydrophobic_potential_positive() {
        let hp = HydrophobicPotential::default_carbon(3);
        let positions = vec![[0.0; 3], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0]];
        let e = hp.energy(&positions);
        assert!(e >= 0.0, "hydrophobic energy should be non-negative");
    }
    #[test]
    fn test_kabsch_sander_hbond_negative_for_hbond() {
        let o = [0.0, 0.0, 0.0];
        let c = [0.5, 0.0, 0.0];
        let n = [0.0, 2.0, 0.0];
        let h = [0.0, 1.0, 0.0];
        let e = kabsch_sander_hbond(o, c, n, h);
        assert!(e.is_finite());
    }
    #[test]
    fn test_find_ks_hbonds_empty_ss() {
        let ss = SecondaryStructure::new(vec![], vec![], vec![], vec![]);
        let bonds = find_kabsch_sander_hbonds(&ss, -0.5);
        assert!(bonds.is_empty());
    }
    #[test]
    fn test_q_fold_identical_conformations() {
        let positions = vec![[0.0; 3], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let q = q_fold_order_parameter(&positions, &positions, 1.0, 2);
        assert!(
            (q - 1.0).abs() < 1e-8,
            "identical conformation should give Q=1"
        );
    }
    #[test]
    fn test_q_fold_different_conformations() {
        let native = vec![[0.0; 3], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let denatured = vec![[0.0; 3], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let q = q_fold_order_parameter(&denatured, &native, 1.0, 2);
        assert!(q < 1.0, "Q should be less than 1 for denatured state");
    }
    #[test]
    fn test_detect_disulfide_bonds_found() {
        let positions = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        let indices = vec![5usize, 12usize];
        let bonds = detect_disulfide_bonds(&positions, &indices, 2.2);
        assert_eq!(bonds.len(), 1);
        assert!((bonds[0].2 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_detect_disulfide_bonds_none() {
        let positions = vec![[0.0; 3], [5.0, 0.0, 0.0]];
        let indices = vec![3usize, 9usize];
        let bonds = detect_disulfide_bonds(&positions, &indices, 2.2);
        assert!(bonds.is_empty());
    }
    #[test]
    fn test_identify_salt_bridges_found() {
        let pos = vec![[0.0; 3]];
        let neg = vec![[3.5, 0.0, 0.0]];
        let bridges = identify_salt_bridges(&pos, &neg, 4.0);
        assert_eq!(bridges.len(), 1);
    }
    #[test]
    fn test_identify_salt_bridges_none() {
        let pos = vec![[0.0; 3]];
        let neg = vec![[10.0, 0.0, 0.0]];
        let bridges = identify_salt_bridges(&pos, &neg, 4.0);
        assert!(bridges.is_empty());
    }
    #[test]
    fn test_salt_bridge_energy_attractive() {
        let e = salt_bridge_energy(3.0, 4.0);
        assert!(
            e < 0.0,
            "salt bridge energy should be negative (attractive)"
        );
    }
    #[test]
    fn test_protein_secondary_structure_classify_helix() {
        let mut pss = ProteinSecondaryStructure::new(5);
        pss.set_angles(1, -57.0, -47.0);
        let label = pss.classify(1);
        assert_eq!(label, SecStructLabel::Helix);
    }
    #[test]
    fn test_protein_secondary_structure_classify_sheet() {
        let mut pss = ProteinSecondaryStructure::new(5);
        pss.set_angles(1, -119.0, 113.0);
        let label = pss.classify(1);
        assert_eq!(label, SecStructLabel::Sheet);
    }
    #[test]
    fn test_protein_secondary_structure_classify_coil() {
        let mut pss = ProteinSecondaryStructure::new(5);
        pss.set_angles(1, 60.0, 60.0);
        let label = pss.classify(1);
        assert_eq!(label, SecStructLabel::Coil);
    }
    #[test]
    fn test_protein_secondary_structure_all_residues() {
        let pss = ProteinSecondaryStructure::new(10);
        let all = pss.classify_all();
        assert_eq!(all.len(), 10);
    }
    #[test]
    fn test_ramachandran_analysis_helix_allowed() {
        let ra = RamachandranAnalysis::new();
        assert!(ra.is_helix_region(-57.0, -47.0));
    }
    #[test]
    fn test_ramachandran_analysis_sheet_allowed() {
        let ra = RamachandranAnalysis::new();
        assert!(ra.is_sheet_region(-119.0, 113.0));
    }
    #[test]
    fn test_ramachandran_disallowed_region() {
        let ra = RamachandranAnalysis::new();
        assert!(!ra.is_allowed(60.0, -120.0));
    }
    #[test]
    fn test_ramachandran_analysis_left_helix() {
        let ra = RamachandranAnalysis::new();
        assert!(ra.is_allowed(57.0, 47.0));
    }
    #[test]
    fn test_hydrogen_bond_detected_short_distance() {
        let hb = HydrogenBond::new([0.0; 3], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]);
        assert!(hb.is_hbond());
    }
    #[test]
    fn test_hydrogen_bond_not_detected_long_distance() {
        let hb = HydrogenBond::new([0.0; 3], [0.0, 1.0, 0.0], [0.0, 5.0, 0.0]);
        assert!(!hb.is_hbond());
    }
    #[test]
    fn test_hydrogen_bond_energy_negative() {
        let hb = HydrogenBond::new([0.0; 3], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]);
        assert!(hb.energy() < 0.0);
    }
    #[test]
    fn test_hydrogen_bond_distance() {
        let hb = HydrogenBond::new([0.0; 3], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]);
        assert!((hb.donor_acceptor_distance() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_sas_single_atom() {
        let sas = SolventAccessibleSurface::new(vec![[0.0; 3]], vec![1.7], 1.4);
        let area = sas.compute_atom(0);
        let expected = 4.0 * std::f64::consts::PI * (1.7_f64 + 1.4_f64).powi(2);
        assert!((area - expected).abs() < 1.0);
    }
    #[test]
    fn test_sas_two_distant_atoms() {
        let sas =
            SolventAccessibleSurface::new(vec![[0.0; 3], [100.0, 0.0, 0.0]], vec![1.7, 1.7], 1.4);
        let total = sas.total();
        let single = 4.0 * std::f64::consts::PI * (1.7_f64 + 1.4_f64).powi(2);
        assert!((total - 2.0 * single).abs() < 5.0);
    }
    #[test]
    fn test_sas_total_positive() {
        let sas =
            SolventAccessibleSurface::new(vec![[0.0; 3], [3.0, 0.0, 0.0]], vec![1.7, 1.7], 1.4);
        assert!(sas.total() > 0.0);
    }
    #[test]
    fn test_debye_huckel_screened_potential() {
        let dh = DebyeHuckelProtein::new(0.15, 80.0, 300.0);
        let phi = dh.screened_potential(1.0, 5.0);
        assert!(phi.is_finite());
    }
    #[test]
    fn test_debye_huckel_decays_with_distance() {
        let dh = DebyeHuckelProtein::new(0.15, 80.0, 300.0);
        let phi1 = dh.screened_potential(1.0, 3.0);
        let phi2 = dh.screened_potential(1.0, 10.0);
        assert!(
            phi1.abs() > phi2.abs(),
            "Screened potential should decay with distance"
        );
    }
    #[test]
    fn test_debye_huckel_length() {
        let dh = DebyeHuckelProtein::new(0.15, 80.0, 300.0);
        let kappa = dh.inverse_debye_length();
        assert!(kappa > 0.0);
    }
    #[test]
    fn test_rg_protein_single_residue() {
        let rg = RgProtein::new(vec![[0.0; 3]]);
        assert_eq!(rg.rg(), 0.0);
    }
    #[test]
    fn test_rg_protein_symmetric() {
        let rg = RgProtein::new(vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert!((rg.rg() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rg_protein_center_of_mass() {
        let rg = RgProtein::new(vec![[0.0; 3], [2.0, 0.0, 0.0]]);
        let com = rg.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rmsd_identical_structures() {
        let pos = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let rmsd = RmsdCalculation::rmsd(&pos, &pos);
        assert!(rmsd.abs() < 1e-12);
    }
    #[test]
    fn test_rmsd_shifted_structure() {
        let pos1 = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        let pos2 = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let rmsd = RmsdCalculation::rmsd(&pos1, &pos2);
        assert!((rmsd - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rmsd_kabsch_superposition() {
        let pos1 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let pos2 = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        let rmsd = RmsdCalculation::kabsch_rmsd(&pos1, &pos2);
        assert!(
            rmsd < 1e-8,
            "Kabsch RMSD after optimal superposition should be ~0"
        );
    }
    #[test]
    fn test_secondary_structure_stats_fractions_sum() {
        let labels = vec![
            SecStructLabel::Helix,
            SecStructLabel::Sheet,
            SecStructLabel::Coil,
            SecStructLabel::Helix,
            SecStructLabel::Coil,
        ];
        let stats = SecondaryStructureStats::from_labels(&labels);
        let total = stats.helix_fraction + stats.sheet_fraction + stats.coil_fraction;
        assert!((total - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_secondary_structure_stats_all_helix() {
        let labels = vec![SecStructLabel::Helix; 5];
        let stats = SecondaryStructureStats::from_labels(&labels);
        assert!((stats.helix_fraction - 1.0).abs() < 1e-12);
        assert_eq!(stats.sheet_fraction, 0.0);
    }
    #[test]
    fn test_secondary_structure_stats_empty() {
        let stats = SecondaryStructureStats::from_labels(&[]);
        assert_eq!(stats.helix_fraction, 0.0);
    }
    #[test]
    fn test_protein_elasticity_spring_constant_cutoff() {
        let pe = ProteinElasticity::new(7.0, 1.0);
        let k = pe.spring_constant(5.0);
        assert!(k > 0.0);
    }
    #[test]
    fn test_protein_elasticity_no_spring_beyond_cutoff() {
        let pe = ProteinElasticity::new(7.0, 1.0);
        let k = pe.spring_constant(10.0);
        assert_eq!(k, 0.0);
    }
    #[test]
    fn test_protein_elasticity_network_energy() {
        let pe = ProteinElasticity::new(7.0, 1.0);
        let positions = vec![[0.0; 3], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let e = pe.elastic_energy(&positions);
        assert!(e >= 0.0);
    }
    #[test]
    fn test_protein_elasticity_pair_count() {
        let pe = ProteinElasticity::new(7.0, 1.0);
        let positions = vec![[0.0; 3], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]];
        let pairs = pe.contact_pairs(&positions);
        assert_eq!(pairs.len(), 2);
    }
}
