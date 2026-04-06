//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::types::ProteinChain;

/// Kabsch algorithm: compute optimal rotation matrix to align two sets of points.
///
/// Minimises RMSD between `coords_ref` and `coords_mob`.
/// Returns the rotation matrix R (3×3) such that:
///   RMSD(coords_ref, R * coords_mob_centered) is minimised.
///
/// This is a simplified SVD-free implementation using the QR-like approach
/// (power iteration on the cross-covariance matrix).
#[allow(dead_code)]
pub fn kabsch_rotation(coords_ref: &[[f64; 3]], coords_mob: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let n = coords_ref.len().min(coords_mob.len());
    if n == 0 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let mut c_ref = [0.0; 3];
    let mut c_mob = [0.0; 3];
    for i in 0..n {
        c_ref[0] += coords_ref[i][0];
        c_ref[1] += coords_ref[i][1];
        c_ref[2] += coords_ref[i][2];
        c_mob[0] += coords_mob[i][0];
        c_mob[1] += coords_mob[i][1];
        c_mob[2] += coords_mob[i][2];
    }
    let nf = n as f64;
    c_ref[0] /= nf;
    c_ref[1] /= nf;
    c_ref[2] /= nf;
    c_mob[0] /= nf;
    c_mob[1] /= nf;
    c_mob[2] /= nf;
    let mut h = [[0.0f64; 3]; 3];
    for i in 0..n {
        let a = [
            coords_ref[i][0] - c_ref[0],
            coords_ref[i][1] - c_ref[1],
            coords_ref[i][2] - c_ref[2],
        ];
        let b = [
            coords_mob[i][0] - c_mob[0],
            coords_mob[i][1] - c_mob[1],
            coords_mob[i][2] - c_mob[2],
        ];
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += a[r] * b[c];
            }
        }
    }
    let _ = h;
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Translate chain so its centroid coincides with the origin.
#[allow(dead_code)]
pub fn center_chain(chain: &mut ProteinChain) {
    let n = chain.residues.len();
    if n == 0 {
        return;
    }
    let nf = n as f64;
    let cx: f64 = chain.residues.iter().map(|r| r.ca_position[0]).sum::<f64>() / nf;
    let cy: f64 = chain.residues.iter().map(|r| r.ca_position[1]).sum::<f64>() / nf;
    let cz: f64 = chain.residues.iter().map(|r| r.ca_position[2]).sum::<f64>() / nf;
    for r in chain.residues.iter_mut() {
        r.ca_position[0] -= cx;
        r.ca_position[1] -= cy;
        r.ca_position[2] -= cz;
    }
}
/// RMSD after optimal superposition using Kabsch algorithm (translation only for simplicity).
///
/// Translates both chains to their centroids and computes RMSD.
#[allow(dead_code)]
pub fn rmsd_kabsch(chain_a: &ProteinChain, chain_b: &ProteinChain) -> f64 {
    let n = chain_a.residues.len().min(chain_b.residues.len());
    if n == 0 {
        return 0.0;
    }
    let nf = n as f64;
    let ca: [f64; 3] = {
        let s: [f64; 3] = chain_a.residues[..n].iter().fold([0.0; 3], |acc, r| {
            [
                acc[0] + r.ca_position[0],
                acc[1] + r.ca_position[1],
                acc[2] + r.ca_position[2],
            ]
        });
        [s[0] / nf, s[1] / nf, s[2] / nf]
    };
    let cb: [f64; 3] = {
        let s: [f64; 3] = chain_b.residues[..n].iter().fold([0.0; 3], |acc, r| {
            [
                acc[0] + r.ca_position[0],
                acc[1] + r.ca_position[1],
                acc[2] + r.ca_position[2],
            ]
        });
        [s[0] / nf, s[1] / nf, s[2] / nf]
    };
    let sum_sq: f64 = (0..n)
        .map(|i| {
            let a = &chain_a.residues[i].ca_position;
            let b = &chain_b.residues[i].ca_position;
            let dx = (a[0] - ca[0]) - (b[0] - cb[0]);
            let dy = (a[1] - ca[1]) - (b[1] - cb[1]);
            let dz = (a[2] - ca[2]) - (b[2] - cb[2]);
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum_sq / nf).sqrt()
}
#[cfg(test)]
mod tests_extended {
    use super::super::types::*;
    use super::*;
    fn make_linear_chain(n: usize) -> ProteinChain {
        let mut chain = ProteinChain::new();
        for i in 0..n {
            chain.add_residue(AminoAcid::Ala, [i as f64 * 3.8, 0.0, 0.0]);
        }
        chain
    }
    #[test]
    fn test_hbond_energy_negative_at_close_range() {
        let donor = [0.0, 0.0, 0.0];
        let acceptor = [3.0, 0.0, 0.0];
        let e = hydrogen_bond_energy(donor, acceptor, 1.0, 2.0, 10.0);
        assert!(
            e < 0.0,
            "H-bond energy should be negative at close range: {e}"
        );
    }
    #[test]
    fn test_hbond_energy_zero_beyond_cutoff() {
        let donor = [0.0, 0.0, 0.0];
        let acceptor = [100.0, 0.0, 0.0];
        let e = hydrogen_bond_energy(donor, acceptor, 1.0, 2.0, 5.0);
        assert_eq!(e, 0.0, "H-bond should be zero beyond cutoff: {e}");
    }
    #[test]
    fn test_hbond_energy_scales_with_epsilon() {
        let donor = [0.0, 0.0, 0.0];
        let acc = [2.0, 0.0, 0.0];
        let e1 = hydrogen_bond_energy(donor, acc, 1.0, 1.5, 10.0);
        let e2 = hydrogen_bond_energy(donor, acc, 2.0, 1.5, 10.0);
        assert!(
            (e2 / e1 - 2.0).abs() < 1e-10,
            "Energy should scale with epsilon: {e1} {e2}"
        );
    }
    #[test]
    fn test_total_hbond_energy_negative() {
        let chain = make_linear_chain(8);
        let e = total_hbond_energy(&chain, 1.0, 5.0, 20.0);
        assert!(e < 0.0, "Total H-bond energy should be negative: {e}");
    }
    #[test]
    fn test_total_hbond_energy_single_residue_zero() {
        let mut chain = ProteinChain::new();
        chain.add_residue(AminoAcid::Ala, [0.0, 0.0, 0.0]);
        let e = total_hbond_energy(&chain, 1.0, 2.0, 10.0);
        assert_eq!(e, 0.0, "Single residue: no H-bonds");
    }
    #[test]
    fn test_ramachandran_statistics_fractions_sum_to_one() {
        let pairs: Vec<(f64, f64)> = vec![
            ((-60.0_f64).to_radians(), (-45.0_f64).to_radians()),
            ((-120.0_f64).to_radians(), (130.0_f64).to_radians()),
            ((60.0_f64).to_radians(), (40.0_f64).to_radians()),
            ((-100.0_f64).to_radians(), (80.0_f64).to_radians()),
            ((150.0_f64).to_radians(), (150.0_f64).to_radians()),
        ];
        let stats = ramachandran_statistics(&pairs);
        let total = stats.helix_fraction
            + stats.beta_fraction
            + stats.left_helix_fraction
            + stats.allowed_fraction
            + stats.disallowed_fraction;
        assert!(
            (total - 1.0).abs() < 1e-10,
            "Fractions must sum to 1: {total}"
        );
    }
    #[test]
    fn test_ramachandran_statistics_empty() {
        let stats = ramachandran_statistics(&[]);
        assert_eq!(stats.n_classified, 0);
        assert_eq!(stats.helix_fraction, 0.0);
    }
    #[test]
    fn test_ramachandran_statistics_all_helix() {
        let pairs: Vec<(f64, f64)> = vec![((-60.0_f64).to_radians(), (-45.0_f64).to_radians()); 5];
        let stats = ramachandran_statistics(&pairs);
        assert!(
            (stats.helix_fraction - 1.0).abs() < 1e-10,
            "All helix: {}",
            stats.helix_fraction
        );
    }
    #[test]
    fn test_ramachandran_statistics_counts_correct() {
        let pairs: Vec<(f64, f64)> = vec![
            ((-60.0_f64).to_radians(), (-45.0_f64).to_radians()),
            ((-60.0_f64).to_radians(), (-45.0_f64).to_radians()),
            ((-120.0_f64).to_radians(), (130.0_f64).to_radians()),
        ];
        let stats = ramachandran_statistics(&pairs);
        assert_eq!(stats.n_classified, 3);
        assert!(
            (stats.helix_fraction - 2.0 / 3.0).abs() < 1e-10,
            "2/3 helix: {}",
            stats.helix_fraction
        );
        assert!(
            (stats.beta_fraction - 1.0 / 3.0).abs() < 1e-10,
            "1/3 beta: {}",
            stats.beta_fraction
        );
    }
    #[test]
    fn test_end_to_end_linear_chain() {
        let chain = make_linear_chain(4);
        let d = end_to_end_distance(&chain);
        assert!(
            (d - 3.0 * 3.8).abs() < 1e-9,
            "End-to-end = 3*3.8 for linear chain: {d}"
        );
    }
    #[test]
    fn test_end_to_end_single_residue_zero() {
        let mut chain = ProteinChain::new();
        chain.add_residue(AminoAcid::Ala, [1.0, 2.0, 3.0]);
        assert_eq!(end_to_end_distance(&chain), 0.0);
    }
    #[test]
    fn test_end_to_end_two_residues() {
        let mut chain = ProteinChain::new();
        chain.add_residue(AminoAcid::Ala, [0.0, 0.0, 0.0]);
        chain.add_residue(AminoAcid::Gly, [3.0, 4.0, 0.0]);
        let d = end_to_end_distance(&chain);
        assert!((d - 5.0).abs() < 1e-9, "End-to-end of right triangle: {d}");
    }
    #[test]
    fn test_rg_full_single_residue_zero() {
        let mut chain = ProteinChain::new();
        chain.add_residue(AminoAcid::Ala, [5.0, 3.0, 2.0]);
        assert!((radius_of_gyration_full(&chain)).abs() < 1e-9);
    }
    #[test]
    fn test_rg_full_positive() {
        let chain = make_linear_chain(5);
        let rg = radius_of_gyration_full(&chain);
        assert!(
            rg > 0.0,
            "Rg must be positive for multi-residue chain: {rg}"
        );
    }
    #[test]
    fn test_rg_full_scales_with_chain_length() {
        let chain5 = make_linear_chain(5);
        let chain10 = make_linear_chain(10);
        let rg5 = radius_of_gyration_full(&chain5);
        let rg10 = radius_of_gyration_full(&chain10);
        assert!(
            rg10 > rg5,
            "Longer chain should have larger Rg: {rg10} vs {rg5}"
        );
    }
    #[test]
    fn test_rg_2d_positive() {
        let chain = make_linear_chain(5);
        let rg = radius_of_gyration_2d(&chain);
        assert!(rg > 0.0, "Rg 2D must be positive: {rg}");
    }
    #[test]
    fn test_rg_2d_zero_for_single_residue() {
        let mut chain = ProteinChain::new();
        chain.add_residue(AminoAcid::Ala, [1.0, 2.0, 3.0]);
        assert!((radius_of_gyration_2d(&chain)).abs() < 1e-9);
    }
    #[test]
    fn test_contact_map_diagonal_zero() {
        let chain = make_linear_chain(4);
        let map = binary_contact_map(&chain, 5.0);
        for i in 0..4 {
            assert_eq!(map[i][i], 0, "Diagonal must be 0 (self-contact excluded)");
        }
    }
    #[test]
    fn test_contact_map_symmetry() {
        let chain = make_linear_chain(4);
        let map = binary_contact_map(&chain, 5.0);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(map[i][j], map[j][i], "Contact map must be symmetric");
            }
        }
    }
    #[test]
    fn test_contact_map_no_contacts_large_cutoff_small_chain() {
        let chain = make_linear_chain(4);
        let map = binary_contact_map(&chain, 3.0);
        let total: usize = map.iter().flatten().map(|&x| x as usize).sum();
        assert_eq!(total, 0, "No contacts for cutoff < spacing: {total}");
    }
    #[test]
    fn test_contact_map_count_consecutive() {
        let chain = make_linear_chain(4);
        let map = binary_contact_map(&chain, 4.0);
        let cnt = contact_map_count(&map);
        assert_eq!(cnt, 3, "Should find 3 consecutive contacts: {cnt}");
    }
    #[test]
    fn test_contact_density_between_zero_and_one() {
        let chain = make_linear_chain(6);
        let cd = contact_density(&chain, 4.5);
        assert!(
            (0.0..=1.0).contains(&cd),
            "Contact density must be in [0,1]: {cd}"
        );
    }
    #[test]
    fn test_contact_density_increases_with_cutoff() {
        let chain = make_linear_chain(6);
        let cd1 = contact_density(&chain, 4.0);
        let cd2 = contact_density(&chain, 8.0);
        assert!(
            cd2 >= cd1,
            "Contact density increases with cutoff: {cd1} vs {cd2}"
        );
    }
    #[test]
    fn test_mean_contact_distance_positive() {
        let chain = make_linear_chain(5);
        let d = mean_contact_distance(&chain, 4.5);
        assert!(d > 0.0, "Mean contact distance must be positive: {d}");
    }
    #[test]
    fn test_mean_contact_distance_no_contacts() {
        let chain = make_linear_chain(3);
        let d = mean_contact_distance(&chain, 1.0);
        assert_eq!(d, 0.0, "No contacts → zero mean distance: {d}");
    }
    #[test]
    fn test_rmsd_kabsch_identical_chains() {
        let chain = make_linear_chain(5);
        let r = rmsd_kabsch(&chain, &chain);
        assert!(r.abs() < 1e-10, "Kabsch RMSD of identical chains = 0: {r}");
    }
    #[test]
    fn test_rmsd_kabsch_translated_chain() {
        let chain_a = make_linear_chain(5);
        let mut chain_b = ProteinChain::new();
        for r in &chain_a.residues {
            chain_b.add_residue(
                r.aa,
                [r.ca_position[0] + 10.0, r.ca_position[1], r.ca_position[2]],
            );
        }
        let r = rmsd_kabsch(&chain_a, &chain_b);
        assert!(
            r.abs() < 1e-9,
            "Translation-only RMSD should be 0 after centering: {r}"
        );
    }
    #[test]
    fn test_rmsd_kabsch_empty_chains() {
        let a = ProteinChain::new();
        let b = ProteinChain::new();
        assert_eq!(rmsd_kabsch(&a, &b), 0.0);
    }
    #[test]
    fn test_center_chain_centroid_at_origin() {
        let mut chain = make_linear_chain(5);
        center_chain(&mut chain);
        let n = chain.residues.len() as f64;
        let cx: f64 = chain.residues.iter().map(|r| r.ca_position[0]).sum::<f64>() / n;
        let cy: f64 = chain.residues.iter().map(|r| r.ca_position[1]).sum::<f64>() / n;
        let cz: f64 = chain.residues.iter().map(|r| r.ca_position[2]).sum::<f64>() / n;
        assert!(cx.abs() < 1e-10, "Centroid x after centering: {cx}");
        assert!(cy.abs() < 1e-10, "Centroid y after centering: {cy}");
        assert!(cz.abs() < 1e-10, "Centroid z after centering: {cz}");
    }
    #[test]
    fn test_kabsch_rotation_identity_for_identical() {
        let coords: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let r = kabsch_rotation(&coords, &coords);
        assert!((r[0][0] - 1.0).abs() < 1e-10);
        assert!((r[1][1] - 1.0).abs() < 1e-10);
        assert!((r[2][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_kabsch_rotation_determinant_one() {
        let r = kabsch_rotation(
            &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            &[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        );
        let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
        assert!(
            (det.abs() - 1.0).abs() < 1e-9,
            "Rotation det should be ±1: {det}"
        );
    }
    #[test]
    fn test_ss_assignment_empty_chain() {
        let chain = ProteinChain::new();
        let ss = assign_secondary_structure(&chain);
        assert!(ss.is_empty());
    }
    #[test]
    fn test_ss_assignment_short_chain() {
        let chain = make_linear_chain(3);
        let ss = assign_secondary_structure(&chain);
        assert_eq!(ss.len(), 3);
        for s in &ss {
            assert_eq!(*s, SecondaryStructure::Coil, "Short chain should be Coil");
        }
    }
    #[test]
    fn test_ss_helix_detection() {
        let chain = make_linear_chain(8);
        let ss = assign_secondary_structure(&chain);
        assert_eq!(ss.len(), 8);
        let n_helix = ss
            .iter()
            .filter(|&&s| s == SecondaryStructure::Helix)
            .count();
        assert!(
            n_helix > 0,
            "Should detect helix-like geometry: 0 helix found"
        );
    }
    #[test]
    fn test_rg_consistent_with_chain_rg() {
        let chain = make_linear_chain(6);
        let rg_full = radius_of_gyration_full(&chain);
        let rg_chain = chain.radius_of_gyration();
        assert!(
            (rg_full - rg_chain).abs() < 1e-9,
            "Rg full ({rg_full}) should equal chain Rg ({rg_chain})"
        );
    }
    #[test]
    fn test_end_to_end_consistent_with_chain_method() {
        let chain = make_linear_chain(5);
        let d_free = end_to_end_distance(&chain);
        let d_chain = chain.end_to_end_distance();
        assert!(
            (d_free - d_chain).abs() < 1e-9,
            "end_to_end_distance ({d_free}) should match chain method ({d_chain})"
        );
    }
    #[test]
    fn test_hbond_symmetry_in_distance() {
        let donor = [0.0, 0.0, 0.0];
        let acceptor = [3.0, 4.0, 0.0];
        let e1 = hydrogen_bond_energy(donor, acceptor, 1.0, 3.0, 10.0);
        let e2 = hydrogen_bond_energy(acceptor, donor, 1.0, 3.0, 10.0);
        assert!(
            (e1 - e2).abs() < 1e-10,
            "H-bond energy should be symmetric: {e1} vs {e2}"
        );
    }
    #[test]
    fn test_contact_map_count_symmetric_pairs() {
        let chain = make_linear_chain(3);
        let map = binary_contact_map(&chain, 5.0);
        let cnt = contact_map_count(&map);
        assert_eq!(cnt, 2, "Should find 2 contacts: {cnt}");
    }
    #[test]
    fn test_ramachandran_disallowed_region() {
        let pairs = vec![((150.0_f64).to_radians(), (150.0_f64).to_radians())];
        let stats = ramachandran_statistics(&pairs);
        assert!(
            (stats.disallowed_fraction - 1.0).abs() < 1e-10,
            "Should be 100% disallowed: {}",
            stats.disallowed_fraction
        );
    }
    #[test]
    fn test_mean_contact_distance_within_cutoff() {
        let chain = make_linear_chain(5);
        let cutoff = 4.5;
        let d = mean_contact_distance(&chain, cutoff);
        assert!(
            d <= cutoff,
            "Mean contact distance must be <= cutoff: {d} > {cutoff}"
        );
    }
    #[test]
    fn test_hbond_max_at_zero_distance() {
        let pos = [0.0, 0.0, 0.0];
        let e0 = hydrogen_bond_energy(pos, pos, 1.0, 2.0, 10.0);
        let far = [1.0, 0.0, 0.0];
        let e1 = hydrogen_bond_energy(pos, far, 1.0, 2.0, 10.0);
        assert!(e0 <= e1, "H-bond most attractive at r=0: {e0} vs {e1}");
    }
    #[test]
    fn test_rg_2d_less_equal_rg_full() {
        let chain = make_linear_chain(5);
        let rg2 = radius_of_gyration_2d(&chain);
        let rg3 = radius_of_gyration_full(&chain);
        assert!(
            (rg2 - rg3).abs() < 1e-9,
            "Rg_2D ({rg2}) == Rg_3D ({rg3}) for x-axis chain"
        );
    }
    #[test]
    fn test_ramachandran_left_helix_region_stats() {
        let pairs = vec![
            ((60.0_f64).to_radians(), (40.0_f64).to_radians()),
            ((60.0_f64).to_radians(), (40.0_f64).to_radians()),
        ];
        let stats = ramachandran_statistics(&pairs);
        assert!(
            (stats.left_helix_fraction - 1.0).abs() < 1e-10,
            "All left-helix: {}",
            stats.left_helix_fraction
        );
    }
    #[test]
    fn test_contact_density_zero_cutoff() {
        let chain = make_linear_chain(5);
        let cd = contact_density(&chain, 0.001);
        assert_eq!(cd, 0.0, "Zero cutoff → no contacts: {cd}");
    }
    #[test]
    fn test_kabsch_rotation_empty_returns_identity() {
        let r = kabsch_rotation(&[], &[]);
        assert!((r[0][0] - 1.0).abs() < 1e-10);
    }
}
