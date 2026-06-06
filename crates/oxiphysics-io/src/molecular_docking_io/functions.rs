//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{DockingPose, PdbqtAtom, PoseCluster, SmilesAtom, VinaResult};

/// Parse a single ATOM/HETATM line from a PDBQT file.
///
/// Returns `None` if the line is too short or unparseable.
pub fn parse_pdbqt_atom_line(line: &str) -> Option<PdbqtAtom> {
    if line.len() < 60 {
        return None;
    }
    let serial = line.get(6..11)?.trim().parse::<u32>().ok()?;
    let name = line.get(12..16)?.trim().to_string();
    let res_name = line.get(17..20)?.trim().to_string();
    let chain_id = line.chars().nth(21).unwrap_or('A');
    let res_seq = line.get(22..26)?.trim().parse::<i32>().ok()?;
    let x = line.get(30..38)?.trim().parse::<f64>().ok()?;
    let y = line.get(38..46)?.trim().parse::<f64>().ok()?;
    let z = line.get(46..54)?.trim().parse::<f64>().ok()?;
    let occupancy = line
        .get(54..60)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(1.0);
    let temp_factor = line
        .get(60..66)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);
    let charge = line
        .get(70..76)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);
    let atom_type = line.get(77..).unwrap_or("C").trim().to_string();
    Some(PdbqtAtom::new(
        serial,
        &name,
        &res_name,
        chain_id,
        res_seq,
        [x, y, z],
        occupancy,
        temp_factor,
        charge,
        &atom_type,
    ))
}
/// Tokenize element symbols from a SMILES string (organic subset).
pub(super) fn tokenize_smiles_atoms(smiles: &str) -> Vec<SmilesAtom> {
    let organic_subset = ["Cl", "Br", "C", "N", "O", "S", "P", "F", "I", "B"];
    let mut atoms = Vec::new();
    let chars: Vec<char> = smiles.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let ch = chars[i];
        if i + 1 < chars.len() {
            let two: String = [ch, chars[i + 1]].iter().collect();
            if organic_subset.contains(&two.as_str()) {
                atoms.push(SmilesAtom::new(&two, 0, false));
                i += 2;
                continue;
            }
        }
        if "cnospb".contains(ch) {
            let el = ch.to_ascii_uppercase().to_string();
            atoms.push(SmilesAtom::new(&el, 0, true));
            i += 1;
            continue;
        }
        if ch.is_ascii_uppercase() {
            atoms.push(SmilesAtom::new(&ch.to_string(), 0, false));
            i += 1;
            continue;
        }
        if ch == '['
            && let Some(end) = smiles[i..].find(']')
        {
            let bracket = &smiles[i + 1..i + end];
            let el = bracket
                .chars()
                .filter(|c| c.is_ascii_alphabetic())
                .take(2)
                .collect::<String>();
            let charge = if bracket.contains('+') {
                1
            } else if bracket.contains('-') {
                -1
            } else {
                0
            };
            atoms.push(SmilesAtom::new(&el, charge, false));
            i += end + 1;
            continue;
        }
        i += 1;
    }
    atoms
}
/// Parse Vina log output into a list of `VinaResult` entries.
///
/// Looks for lines of the form: `   1      -9.0      0.000      0.000`
pub fn parse_vina_log(log: &str) -> Vec<VinaResult> {
    let mut results = Vec::new();
    let mut in_results = false;
    for line in log.lines() {
        if line.contains("mode") && line.contains("affinity") {
            in_results = true;
            continue;
        }
        if in_results {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4
                && let (Ok(mode), Ok(aff), Ok(lb), Ok(ub)) = (
                    parts[0].parse::<u32>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                )
            {
                results.push(VinaResult::new(mode, aff, lb, ub));
            }
        }
    }
    results
}
/// Perform greedy RMSD-based clustering of docking poses.
///
/// * `poses`       – list of docking poses.
/// * `rmsd_cutoff` – maximum RMSD (Angstrom) to be in the same cluster.
///
/// Returns a list of `PoseCluster` instances.
pub fn cluster_poses_by_rmsd(poses: &[DockingPose], rmsd_cutoff: f64) -> Vec<PoseCluster> {
    let mut clusters: Vec<PoseCluster> = Vec::new();
    let mut assigned = vec![false; poses.len()];
    for i in 0..poses.len() {
        if assigned[i] {
            continue;
        }
        let mut cluster = PoseCluster::new(clusters.len(), rmsd_cutoff);
        cluster.pose_indices.push(i);
        cluster.best_pose_index = i;
        assigned[i] = true;
        for j in (i + 1)..poses.len() {
            if assigned[j] {
                continue;
            }
            if let Some(rmsd) = poses[i].rmsd_to(&poses[j])
                && rmsd <= rmsd_cutoff
            {
                cluster.pose_indices.push(j);
                assigned[j] = true;
                if poses[j].binding_affinity < poses[cluster.best_pose_index].binding_affinity {
                    cluster.best_pose_index = j;
                }
            }
        }
        clusters.push(cluster);
    }
    clusters
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecular_docking_io::BindingAffinityScore;
    use crate::molecular_docking_io::DockingLog;
    use crate::molecular_docking_io::DockingRestraint;
    use crate::molecular_docking_io::GridBox;
    use crate::molecular_docking_io::Mol2Atom;
    use crate::molecular_docking_io::Mol2Bond;
    use crate::molecular_docking_io::Mol2Molecule;
    use crate::molecular_docking_io::PdbqtMolecule;
    use crate::molecular_docking_io::PharmacophoreFeature;
    use crate::molecular_docking_io::PharmacophoreFeatureType;
    use crate::molecular_docking_io::PharmacophoreHypothesis;
    use crate::molecular_docking_io::ReceptorLigandComplex;
    use crate::molecular_docking_io::RestraintType;
    use crate::molecular_docking_io::SdfAtom;
    use crate::molecular_docking_io::SdfEntry;
    use crate::molecular_docking_io::SmilesMolecule;
    use crate::molecular_docking_io::VirtualScreeningResult;
    use crate::molecular_docking_io::VirtualScreeningResults;
    #[test]
    fn test_pdbqt_atom_line_format() {
        let atom = PdbqtAtom::new(
            1,
            "CA",
            "ALA",
            'A',
            10,
            [1.0, 2.0, 3.0],
            1.0,
            0.0,
            -0.15,
            "C",
        );
        let line = atom.to_pdbqt_line();
        assert!(line.starts_with("ATOM"), "line={line}");
        assert!(line.contains("ALA"));
        assert!(line.contains("CA"));
    }
    #[test]
    fn test_parse_pdbqt_atom_line_short() {
        let result = parse_pdbqt_atom_line("ATOM    1  CA  ALA");
        assert!(result.is_none());
    }
    #[test]
    fn test_pdbqt_molecule_centre_of_mass_empty() {
        let mol = PdbqtMolecule::new("empty");
        let com = mol.centre_of_mass();
        for c in com {
            assert!(c.abs() < 1e-10);
        }
    }
    #[test]
    fn test_pdbqt_molecule_centre_of_mass() {
        let mut mol = PdbqtMolecule::new("test");
        mol.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [2.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        mol.add_atom(PdbqtAtom::new(
            2,
            "C2",
            "LIG",
            'A',
            1,
            [4.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let com = mol.centre_of_mass();
        assert!((com[0] - 3.0).abs() < 1e-10, "com={:?}", com);
    }
    #[test]
    fn test_pdbqt_molecule_to_pdbqt_string() {
        let mol = PdbqtMolecule::new("LIG");
        let s = mol.to_pdbqt_string();
        assert!(s.contains("REMARK"));
        assert!(s.contains("END"));
    }
    #[test]
    fn test_docking_pose_rmsd_same() {
        let mut p1 = DockingPose::new(1, -8.0, 0.0, 0.0);
        let mut p2 = DockingPose::new(2, -7.5, 0.5, 1.0);
        p1.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        p2.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let rmsd = p1.rmsd_to(&p2).unwrap();
        assert!(rmsd.abs() < 1e-10, "rmsd={rmsd}");
    }
    #[test]
    fn test_docking_pose_rmsd_different() {
        let mut p1 = DockingPose::new(1, -8.0, 0.0, 0.0);
        let mut p2 = DockingPose::new(2, -7.5, 0.5, 1.0);
        p1.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [0.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        p2.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [3.0, 4.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let rmsd = p1.rmsd_to(&p2).unwrap();
        assert!((rmsd - 5.0).abs() < 1e-10, "rmsd={rmsd}");
    }
    #[test]
    fn test_docking_pose_rmsd_mismatched_atoms() {
        let mut p1 = DockingPose::new(1, -8.0, 0.0, 0.0);
        let p2 = DockingPose::new(2, -7.5, 0.5, 1.0);
        p1.add_atom(PdbqtAtom::new(
            1, "C1", "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, 0.0, "C",
        ));
        assert!(p1.rmsd_to(&p2).is_none());
    }
    #[test]
    fn test_binding_affinity_delta_g() {
        let s = BindingAffinityScore::new(-10.0, -8.0, -1.5, -0.5, -1.0);
        assert!((s.estimated_delta_g() - (-8.5)).abs() < 1e-10);
    }
    #[test]
    fn test_binding_affinity_remark_block() {
        let s = BindingAffinityScore::new(-9.0, -7.5, -1.0, -0.3, -1.2);
        let r = s.to_remark_block();
        assert!(r.contains("VINA RESULT"));
    }
    #[test]
    fn test_smiles_molecule_from_simple() {
        let mol = SmilesMolecule::from_smiles("CC", "ethane");
        assert_eq!(mol.atoms.len(), 2);
    }
    #[test]
    fn test_smiles_molecule_formula() {
        let mol = SmilesMolecule::from_smiles("CCO", "ethanol");
        let formula = mol.molecular_formula();
        assert!(formula.contains('C'), "formula={formula}");
    }
    #[test]
    fn test_smiles_roundtrip() {
        let smiles = "c1ccccc1";
        let mol = SmilesMolecule::from_smiles(smiles, "benzene");
        assert_eq!(mol.to_smiles(), smiles);
    }
    #[test]
    fn test_mol2_atom_line() {
        let a = Mol2Atom::new(1, "C1", [1.0, 2.0, 3.0], "C.3", 1, "LIG", -0.1);
        let line = a.to_mol2_line();
        assert!(line.contains("C.3"));
    }
    #[test]
    fn test_mol2_bond_line() {
        let b = Mol2Bond::new(1, 1, 2, "1");
        let line = b.to_mol2_line();
        assert!(line.contains('1'));
    }
    #[test]
    fn test_mol2_molecule_to_string() {
        let mut mol = Mol2Molecule::new("LIG", "SMALL", "GASTEIGER");
        mol.atoms
            .push(Mol2Atom::new(1, "C1", [0.0; 3], "C.3", 1, "LIG", 0.0));
        let s = mol.to_mol2_string();
        assert!(s.contains("@<TRIPOS>MOLECULE"));
        assert!(s.contains("@<TRIPOS>ATOM"));
    }
    #[test]
    fn test_sdf_entry_to_string() {
        let mut entry = SdfEntry::new("Aspirin");
        entry.atoms.push(SdfAtom::new([1.0, 0.0, 0.0], "C", 0));
        entry.add_tag("BINDING_SCORE", "-8.5");
        let s = entry.to_sdf_string();
        assert!(s.contains("Aspirin"));
        assert!(s.contains("BINDING_SCORE"));
        assert!(s.contains("$$$$"));
    }
    #[test]
    fn test_sdf_atom_v2000_line() {
        let a = SdfAtom::new([1.234, 5.678, 9.012], "O", 0);
        let line = a.to_v2000_line();
        assert!(line.contains('O'));
    }
    #[test]
    fn test_pharmacophore_match() {
        let f =
            PharmacophoreFeature::new(PharmacophoreFeatureType::HbDonor, [0.0, 0.0, 0.0], 2.0, 1.0);
        assert!(f.is_matched([1.0, 0.0, 0.0]));
        assert!(!f.is_matched([5.0, 0.0, 0.0]));
    }
    #[test]
    fn test_pharmacophore_feature_type_str() {
        let f = PharmacophoreFeature::new(PharmacophoreFeatureType::Aromatic, [0.0; 3], 1.5, 1.0);
        assert_eq!(f.feature_type_str(), "AR");
    }
    #[test]
    fn test_pharmacophore_hypothesis_count_matches() {
        let mut hyp = PharmacophoreHypothesis::new("Hyp1");
        hyp.add_feature(PharmacophoreFeature::new(
            PharmacophoreFeatureType::HbAcceptor,
            [0.0, 0.0, 0.0],
            1.0,
            1.0,
        ));
        hyp.add_feature(PharmacophoreFeature::new(
            PharmacophoreFeatureType::Hydrophobic,
            [10.0, 0.0, 0.0],
            1.0,
            1.0,
        ));
        let positions = vec![[0.5, 0.0, 0.0], [9.8, 0.0, 0.0]];
        let count = hyp.count_matches(&positions);
        assert_eq!(count, 2);
    }
    #[test]
    fn test_grid_box_contains() {
        let gb = GridBox::new([0.0, 0.0, 0.0], [20.0, 20.0, 20.0], 0.375);
        assert!(gb.contains([5.0, 5.0, 5.0]));
        assert!(!gb.contains([15.0, 0.0, 0.0]));
    }
    #[test]
    fn test_grid_box_total_points() {
        let gb = GridBox::new([0.0; 3], [10.0, 10.0, 10.0], 1.0);
        let gp = gb.grid_points();
        let total = gb.total_grid_points();
        assert_eq!(total, gp[0] as u64 * gp[1] as u64 * gp[2] as u64);
    }
    #[test]
    fn test_parse_vina_log() {
        let log = "\
Vina 1.2.0
   mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1      -9.0      0.000      0.000
   2      -8.5      1.234      2.456
";
        let results = parse_vina_log(log);
        assert_eq!(results.len(), 2);
        assert!((results[0].affinity - (-9.0)).abs() < 1e-6);
        assert!((results[1].rmsd_lb - 1.234).abs() < 1e-6);
    }
    #[test]
    fn test_cluster_poses_same_position() {
        let mut p1 = DockingPose::new(1, -9.0, 0.0, 0.0);
        let mut p2 = DockingPose::new(2, -8.5, 0.0, 0.0);
        p1.add_atom(PdbqtAtom::new(
            1, "C1", "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, 0.0, "C",
        ));
        p2.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [0.5, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let poses = vec![p1, p2];
        let clusters = cluster_poses_by_rmsd(&poses, 2.0);
        assert_eq!(
            clusters.len(),
            1,
            "Expected 1 cluster, got {}",
            clusters.len()
        );
        assert_eq!(clusters[0].size(), 2);
    }
    #[test]
    fn test_cluster_poses_separate() {
        let mut p1 = DockingPose::new(1, -9.0, 0.0, 0.0);
        let mut p2 = DockingPose::new(2, -8.5, 5.0, 6.0);
        p1.add_atom(PdbqtAtom::new(
            1, "C1", "LIG", 'A', 1, [0.0; 3], 1.0, 0.0, 0.0, "C",
        ));
        p2.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'A',
            1,
            [20.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let poses = vec![p1, p2];
        let clusters = cluster_poses_by_rmsd(&poses, 2.0);
        assert_eq!(
            clusters.len(),
            2,
            "Expected 2 clusters, got {}",
            clusters.len()
        );
    }
    #[test]
    fn test_complex_contact_count() {
        let mut receptor = PdbqtMolecule::new("REC");
        receptor.add_atom(PdbqtAtom::new(
            1, "CA", "ALA", 'A', 1, [0.0; 3], 1.0, 0.0, 0.0, "C",
        ));
        let mut ligand = PdbqtMolecule::new("LIG");
        ligand.add_atom(PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'L',
            1,
            [1.5, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        ));
        let complex = ReceptorLigandComplex::new(receptor, ligand);
        assert_eq!(complex.contact_count(2.0), 1);
        assert_eq!(complex.contact_count(1.0), 0);
    }
    #[test]
    fn test_restraint_distance_satisfied() {
        let r = DockingRestraint::new(
            1,
            RestraintType::Distance {
                receptor_atom: 0,
                ligand_atom: 0,
                target_distance: 3.0,
                tolerance: 1.0,
            },
            10.0,
        );
        let rec = vec![PdbqtAtom::new(
            1, "CA", "ALA", 'A', 1, [0.0; 3], 1.0, 0.0, 0.0, "C",
        )];
        let lig = vec![PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'L',
            1,
            [3.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        )];
        let penalty = r.evaluate_penalty(&rec, &lig);
        assert!(penalty.abs() < 1e-10, "penalty={penalty}");
    }
    #[test]
    fn test_restraint_excluded_sphere_violated() {
        let r = DockingRestraint::new(
            2,
            RestraintType::ExcludedSphere {
                ligand_atom: 0,
                centre: [0.0; 3],
                radius: 5.0,
            },
            100.0,
        );
        let rec: Vec<PdbqtAtom> = Vec::new();
        let lig = vec![PdbqtAtom::new(
            1,
            "C1",
            "LIG",
            'L',
            1,
            [1.0, 0.0, 0.0],
            1.0,
            0.0,
            0.0,
            "C",
        )];
        let penalty = r.evaluate_penalty(&rec, &lig);
        assert!(penalty > 0.0, "penalty={penalty}");
    }
    #[test]
    fn test_virtual_screening_sort() {
        let mut vsr = VirtualScreeningResults::new("receptor");
        vsr.add_result(VirtualScreeningResult::new(
            "CPDA", "CC", -7.0, 3, true, 300.0,
        ));
        vsr.add_result(VirtualScreeningResult::new(
            "CPDB", "CCO", -9.5, 5, true, 400.0,
        ));
        vsr.sort_by_score();
        assert!((vsr.results[0].best_score - (-9.5)).abs() < 1e-10);
    }
    #[test]
    fn test_virtual_screening_top_n() {
        let mut vsr = VirtualScreeningResults::new("receptor");
        for i in 0..5 {
            vsr.add_result(VirtualScreeningResult::new(
                &format!("CPD{i}"),
                "CC",
                -(i as f64),
                1,
                true,
                200.0,
            ));
        }
        let top3 = vsr.top_n(3);
        assert_eq!(top3.len(), 3);
    }
    #[test]
    fn test_virtual_screening_csv() {
        let mut vsr = VirtualScreeningResults::new("receptor");
        vsr.add_result(VirtualScreeningResult::new(
            "CPD1", "CCO", -8.0, 3, true, 320.0,
        ));
        let csv = vsr.to_csv();
        assert!(csv.contains("CPD1"));
        assert!(csv.contains("compound_id"));
    }
    #[test]
    fn test_docking_log_from_vina_log() {
        let log = "\
Vina 1.2.0
   mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1      -9.0      0.000      0.000
";
        let dl = DockingLog::from_vina_log(log);
        assert_eq!(dl.vina_results.len(), 1);
        assert!((dl.best_affinity() - (-9.0)).abs() < 1e-6);
    }
    #[test]
    fn test_docking_log_default_best_affinity_empty() {
        let dl = DockingLog::new();
        assert!(dl.best_affinity().is_infinite());
    }
}
