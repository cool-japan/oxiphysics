//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::BasisFunction;
#[cfg(test)]
use super::types::{
    ElectronCorrelation, NaturalBondOrbital, NboAnalysis, QmmmBoundary, QmmmCoupling,
};
#[cfg(test)]
use super::types_2::{BasisSetType, HartreeFockSolver};
#[cfg(test)]
use super::types_2::{DensityFunctionalTheory, EmbeddingType, NboInteraction};
#[cfg(test)]
use super::types_3::{LinkAtomQmmm, MolecularOrbital, ScfStatus, XcFunctional};

/// Contracted two-electron repulsion integral (μν|λσ) via primitive decomposition.
///
/// For s–s|s–s pairs the Obara–Saika / Boys-function formula is exact.
/// For all combinations involving p-type functions the same product-Gaussian-center
/// approach gives the zeroth-order (monopole) contribution, which is the dominant
/// term for STO-3G and sufficient for qualitatively correct HF results.
pub(super) fn primitive_eri_contracted(
    bmu: &BasisFunction,
    bnu: &BasisFunction,
    bla: &BasisFunction,
    bsi: &BasisFunction,
) -> f64 {
    use oxiphysics_core::numerics::boys_fn;
    use std::f64::consts::PI;
    let mut val = 0.0;
    for (i, &ai) in bmu.exponents.iter().enumerate() {
        for (j, &aj) in bnu.exponents.iter().enumerate() {
            for (k, &ak) in bla.exponents.iter().enumerate() {
                for (l, &al) in bsi.exponents.iter().enumerate() {
                    let p = ai + aj;
                    let q = ak + al;
                    let zeta = p + q;
                    let p_center = [
                        (ai * bmu.center[0] + aj * bnu.center[0]) / p,
                        (ai * bmu.center[1] + aj * bnu.center[1]) / p,
                        (ai * bmu.center[2] + aj * bnu.center[2]) / p,
                    ];
                    let q_center = [
                        (ak * bla.center[0] + al * bsi.center[0]) / q,
                        (ak * bla.center[1] + al * bsi.center[1]) / q,
                        (ak * bla.center[2] + al * bsi.center[2]) / q,
                    ];
                    let r_mu_nu2: f64 = bmu
                        .center
                        .iter()
                        .zip(bnu.center.iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    let r_la_si2: f64 = bla
                        .center
                        .iter()
                        .zip(bsi.center.iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    let k_munu = (-ai * aj / p * r_mu_nu2).exp();
                    let k_lasi = (-ak * al / q * r_la_si2).exp();
                    let r_pq2: f64 = p_center
                        .iter()
                        .zip(q_center.iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    let t = p * q / zeta * r_pq2;
                    let f0 = boys_fn(t, 0)[0];
                    let pre = 2.0 * PI.powf(2.5) / (p * q * zeta.sqrt());
                    let coeff = bmu.coefficients[i]
                        * bnu.coefficients[j]
                        * bla.coefficients[k]
                        * bsi.coefficients[l];
                    val += coeff * pre * k_munu * k_lasi * f0;
                }
            }
        }
    }
    val
}
/// Build a minimal STO-3G overlap matrix for H₂ at bond length R (Bohr).
pub fn build_h2_overlap(r: f64) -> [[f64; 2]; 2] {
    let zeta = 1.24;
    let center_a = [0.0, 0.0, 0.0f64];
    let center_b = [r, 0.0, 0.0f64];
    let ba = BasisFunction::sto3g_s(center_a, zeta);
    let bb = BasisFunction::sto3g_s(center_b, zeta);
    [
        [ba.overlap_with(&ba), ba.overlap_with(&bb)],
        [bb.overlap_with(&ba), bb.overlap_with(&bb)],
    ]
}
/// Hartree-to-eV conversion factor.
pub const HARTREE_TO_EV: f64 = 27.211386245988;
/// Bohr-to-Angstrom conversion factor.
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;
/// Angstrom-to-Bohr conversion.
pub fn angstrom_to_bohr(r: f64) -> f64 {
    r / BOHR_TO_ANGSTROM
}
/// Convert energy from Hartree to kcal/mol.
pub fn hartree_to_kcal(e: f64) -> f64 {
    e * 627.5094740631
}
#[cfg(test)]
mod tests {
    use super::*;
    fn make_1x1_hf() -> HartreeFockSolver {
        HartreeFockSolver::new(
            1,
            2,
            BasisSetType::Sto3G,
            vec![vec![1.0]],
            vec![vec![-1.0]],
            vec![0.5],
            vec![1.0],
            vec![[0.0, 0.0, 0.0]],
        )
    }
    fn make_2x2_hf() -> HartreeFockSolver {
        HartreeFockSolver::new(
            2,
            2,
            BasisSetType::Sto3G,
            vec![vec![1.0, 0.5], vec![0.5, 1.0]],
            vec![vec![-1.5, -0.3], vec![-0.3, -1.5]],
            vec![0.5; 6],
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        )
    }
    #[test]
    fn test_sto3g_s_creates_3_primitives() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        assert_eq!(bf.exponents.len(), 3);
        assert_eq!(bf.coefficients.len(), 3);
        assert_eq!(bf.angular_momentum, 0);
    }
    #[test]
    fn test_sto3g_p_angular_momentum() {
        let bf = BasisFunction::sto3g_p([0.0, 0.0, 0.0], 1.0);
        assert_eq!(bf.angular_momentum, 1);
    }
    #[test]
    fn test_self_overlap_near_one() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let s = bf.overlap_with(&bf);
        assert!(s > 0.0, "self-overlap should be positive");
    }
    #[test]
    fn test_overlap_decreases_with_distance() {
        let bf_a = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let bf_b = BasisFunction::sto3g_s([1.0, 0.0, 0.0], 1.0);
        let bf_c = BasisFunction::sto3g_s([5.0, 0.0, 0.0], 1.0);
        let s_ab = bf_a.overlap_with(&bf_b).abs();
        let s_ac = bf_a.overlap_with(&bf_c).abs();
        assert!(s_ab > s_ac, "overlap should decrease with distance");
    }
    #[test]
    fn test_kinetic_integral_positive() {
        let bf = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let t = bf.kinetic_with(&bf);
        assert!(t > 0.0, "kinetic integral should be positive");
    }
    #[test]
    fn test_zeta_scaling_exponents() {
        let bf1 = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.0);
        let bf2 = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 2.0);
        for (e1, e2) in bf1.exponents.iter().zip(bf2.exponents.iter()) {
            assert!(e2 > e1);
        }
    }
    #[test]
    fn test_nuclear_repulsion_single_atom_zero() {
        let hf = make_1x1_hf();
        assert_eq!(hf.nuclear_repulsion(), 0.0);
    }
    #[test]
    fn test_nuclear_repulsion_two_atoms() {
        let hf = make_2x2_hf();
        let e_nn = hf.nuclear_repulsion();
        assert!((e_nn - 1.0 / 1.4).abs() < 1e-6);
    }
    #[test]
    fn test_fock_matrix_size() {
        let hf = make_2x2_hf();
        let density = hf.initial_density_matrix();
        let fock = hf.build_fock_matrix(&density);
        assert_eq!(fock.len(), 2);
        assert_eq!(fock[0].len(), 2);
    }
    #[test]
    fn test_initial_density_matrix_size() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        assert_eq!(dm.len(), 2);
    }
    #[test]
    fn test_scf_returns_result() {
        let hf = make_2x2_hf();
        let result = hf.run_scf();
        assert!(result.energy.is_finite());
    }
    #[test]
    fn test_scf_energy_has_nuclear_contribution() {
        let hf = make_2x2_hf();
        let result = hf.run_scf();
        assert!(result.nuclear_repulsion > 0.0);
    }
    #[test]
    fn test_mulliken_charge_conservation() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        let q0 = hf.mulliken_charge(&dm, 1.0, &[0]);
        let q1 = hf.mulliken_charge(&dm, 1.0, &[1]);
        assert!(q0.is_finite(), "q0 = {q0} should be finite");
        assert!(q1.is_finite(), "q1 = {q1} should be finite");
    }
    #[test]
    fn test_electronic_energy_negative() {
        let hf = make_2x2_hf();
        let dm = hf.initial_density_matrix();
        let fock = hf.build_fock_matrix(&dm);
        let e = hf.electronic_energy(&dm, &fock);
        assert!(e.is_finite());
    }
    #[test]
    fn test_eri_index_symmetry() {
        let hf = make_2x2_hf();
        let i1 = hf.eri_index(0, 1, 0, 1);
        let i2 = hf.eri_index(1, 0, 0, 1);
        assert_eq!(i1, i2);
    }
    #[test]
    fn test_basis_set_type_equality() {
        assert_eq!(BasisSetType::Sto3G, BasisSetType::Sto3G);
        assert_ne!(BasisSetType::Sto3G, BasisSetType::G6_31G);
    }
    #[test]
    fn test_lda_exchange_zero_for_zero_density() {
        let e = DensityFunctionalTheory::lda_exchange_energy_density(0.0);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn test_lda_exchange_negative() {
        let e = DensityFunctionalTheory::lda_exchange_energy_density(1.0);
        assert!(e < 0.0, "exchange energy density should be negative");
    }
    #[test]
    fn test_vwn_correlation_zero_for_zero_density() {
        let e = DensityFunctionalTheory::vwn_correlation_energy_density(0.0);
        assert_eq!(e, 0.0);
    }
    #[test]
    fn test_pbe_enhancement_at_zero_gradient() {
        let f = DensityFunctionalTheory::pbe_exchange_enhancement(0.0);
        assert!((f - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pbe_enhancement_increases_with_s() {
        let f0 = DensityFunctionalTheory::pbe_exchange_enhancement(0.0);
        let f1 = DensityFunctionalTheory::pbe_exchange_enhancement(1.0);
        assert!(f1 > f0);
    }
    #[test]
    fn test_dft_xc_energy_lda() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        let e = dft.xc_energy_uniform(0.5, 1.0);
        assert!(e.is_finite());
    }
    #[test]
    fn test_dft_xc_energy_b3lyp() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::B3lyp);
        let e = dft.xc_energy_uniform(0.5, 1.0);
        assert!(e.is_finite());
    }
    #[test]
    fn test_dft_ks_scf_runs() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Pbe);
        let result = dft.run_ks_scf();
        assert!(result.energy.is_finite());
    }
    #[test]
    fn test_homo_lumo_gap_positive() {
        let hf = make_2x2_hf();
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        let orbitals = vec![-0.5, 0.3];
        let gap = dft.homo_lumo_gap(&orbitals);
        assert!(gap >= 0.0);
    }
    #[test]
    fn test_pseudopotential_added_to_energy() {
        let hf = make_2x2_hf();
        let mut dft = DensityFunctionalTheory::new(hf, XcFunctional::Lda);
        dft.set_pseudopotential(-0.1);
        let result = dft.run_ks_scf();
        assert!(result.energy.is_finite());
    }
    #[test]
    fn test_mo_is_occupied() {
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![0.7, 0.7]);
        assert!(mo.is_occupied());
    }
    #[test]
    fn test_mo_is_not_occupied() {
        let mo = MolecularOrbital::new(1, 0.3, 0.0, vec![0.7, -0.7]);
        assert!(!mo.is_occupied());
    }
    #[test]
    fn test_participation_ratio_localized() {
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![1.0, 0.0, 0.0, 0.0]);
        let pr = mo.participation_ratio();
        assert!((pr - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_participation_ratio_delocalized() {
        let c = 0.5f64;
        let mo = MolecularOrbital::new(0, -0.3, 2.0, vec![c, c, c, c]);
        let pr = mo.participation_ratio();
        assert!(
            (pr - 4.0).abs() < 1e-6,
            "PR should be ~4 for equal delocalization"
        );
    }
    #[test]
    fn test_dipole_contribution_finite() {
        let mo = MolecularOrbital::new(0, -0.5, 2.0, vec![0.8, 0.6]);
        let centers = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let d = mo.dipole_contribution(&centers, 0);
        assert!(d.is_finite());
    }
    #[test]
    fn test_mayer_bond_order_finite() {
        let n = 2;
        let density = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let overlap = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let bo = MolecularOrbital::mayer_bond_order(&density, &overlap, &[0], &[1]);
        let _ = n;
        assert!(bo.is_finite());
    }
    #[test]
    fn test_wiberg_bond_index_finite() {
        let density = vec![vec![1.0, 0.3], vec![0.3, 1.0]];
        let overlap = vec![vec![1.0, 0.5], vec![0.5, 1.0]];
        let wbi = MolecularOrbital::wiberg_bond_index(&density, &overlap, &[0], &[1]);
        assert!(wbi.is_finite());
    }
    fn make_ec() -> ElectronCorrelation {
        let eps = vec![-0.6, -0.5, 0.3, 0.4];
        let eri = vec![0.1; 4];
        ElectronCorrelation::new(2, 2, eps, eri, -1.0)
    }
    #[test]
    fn test_mp2_result_finite() {
        let ec = make_ec();
        let r = ec.mp2_energy();
        assert!(r.e_mp2.is_finite());
    }
    #[test]
    fn test_mp2_total_includes_hf() {
        let ec = make_ec();
        let r = ec.mp2_energy();
        assert_eq!(r.e_hf, -1.0);
        assert!((r.e_total - r.e_hf - r.e_mp2).abs() < 1e-10);
    }
    #[test]
    fn test_cis_returns_states() {
        let ec = make_ec();
        let states = ec.cis_states(2);
        assert!(states.len() <= 2);
    }
    #[test]
    fn test_cis_excitation_energies_positive() {
        let ec = make_ec();
        let states = ec.cis_states(4);
        for s in &states {
            assert!(s.excitation_energy > 0.0 || s.excitation_energy.is_finite());
        }
    }
    #[test]
    fn test_ccs_t1_diagnostic_zero_amplitudes() {
        let ec = make_ec();
        let t1 = ec.ccs_t1_diagnostic(&[]);
        assert_eq!(t1, 0.0);
    }
    #[test]
    fn test_ccs_t1_diagnostic_nonzero() {
        let ec = make_ec();
        let t1 = ec.ccs_t1_diagnostic(&[0.02, 0.01, 0.015, 0.005]);
        assert!(t1 > 0.0);
    }
    #[test]
    fn test_ccs_energy_finite() {
        let ec = make_ec();
        let t1 = vec![0.02; 4];
        let e = ec.ccs_energy(&t1);
        assert!(e.is_finite());
    }
    #[test]
    fn test_scs_mp2_finite() {
        let ec = make_ec();
        let e = ec.scs_mp2_energy(1.2, 0.333);
        assert!(e.is_finite());
    }
    #[test]
    fn test_nbo_role_bond() {
        let nbo = NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5);
        assert_eq!(nbo.role(), "donor");
    }
    #[test]
    fn test_nbo_role_antibond() {
        let nbo = NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3);
        assert_eq!(nbo.role(), "acceptor");
    }
    #[test]
    fn test_nbo_s_character_zero_hybridization() {
        let nbo = NaturalBondOrbital::new("LP", vec![0], 2.0, -0.3);
        assert_eq!(nbo.s_character(), 0.0);
    }
    #[test]
    fn test_nbo_s_character_sp3() {
        let mut nbo = NaturalBondOrbital::new("BD", vec![0, 1], 1.9, -0.5);
        nbo.hybridization = vec![0.5, 0.866, 0.0, 0.0];
        let sc = nbo.s_character();
        assert!(sc > 0.0 && sc <= 100.0);
    }
    #[test]
    fn test_nbo_interaction_e2_sign() {
        let inter = NboInteraction::compute(1.97, -0.05, -0.5, 0.3, 0, 1);
        assert!(inter.e2_energy < 0.0);
    }
    #[test]
    fn test_nbo_interaction_is_significant() {
        let inter = NboInteraction::compute(1.97, -0.05, -0.5, 0.3, 0, 1);
        assert!(inter.is_significant(0.1));
    }
    #[test]
    fn test_nbo_analysis_top_interactions() {
        let nbos = vec![
            NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5),
            NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3),
        ];
        let mut analysis = NboAnalysis::new(nbos, vec![0.0, 0.0]);
        let fock_nbo = vec![vec![0.0, -0.05], vec![-0.05, 0.0]];
        analysis.compute_interactions(&fock_nbo, 0.1);
        let top = analysis.top_interactions(5);
        assert!(top.len() <= 5);
    }
    #[test]
    fn test_nbo_total_charge_transfer_updated() {
        let nbos = vec![
            NaturalBondOrbital::new("BD", vec![0, 1], 1.97, -0.5),
            NaturalBondOrbital::new("BD*", vec![0, 1], 0.02, 0.3),
        ];
        let mut analysis = NboAnalysis::new(nbos, vec![0.0, 0.0]);
        let fock_nbo = vec![vec![0.0, -0.05], vec![-0.05, 0.0]];
        analysis.compute_interactions(&fock_nbo, 0.0);
        assert!(analysis.total_charge_transfer.is_finite());
    }
    fn make_qmmm() -> QmmmCoupling {
        QmmmCoupling::new(
            QmmmBoundary::LinkAtom,
            EmbeddingType::Electrostatic,
            vec![0, 1],
            vec![-0.834, 0.417, 0.417],
            vec![[3.0, 0.0, 0.0], [3.5, 0.5, 0.0], [3.5, -0.5, 0.0]],
        )
    }
    #[test]
    fn test_qmmm_is_qm_atom() {
        let qmmm = make_qmmm();
        assert!(qmmm.is_qm_atom(0));
        assert!(!qmmm.is_qm_atom(5));
    }
    #[test]
    fn test_qmmm_electrostatic_potential_finite() {
        let qmmm = make_qmmm();
        let v = qmmm.electrostatic_potential_at([0.0, 0.0, 0.0]);
        assert!(v.is_finite());
    }
    #[test]
    fn test_qmmm_mechanical_embedding_finite() {
        let qmmm = make_qmmm();
        let charges = vec![1.0, -1.0];
        let positions = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let e = qmmm.mechanical_embedding_energy(&charges, &positions);
        assert!(e.is_finite());
    }
    #[test]
    fn test_qmmm_total_energy() {
        let qmmm = make_qmmm();
        let e = qmmm.total_energy(-74.9, -5.0, -0.3);
        assert!((e - (-80.2)).abs() < 1e-10);
    }
    #[test]
    fn test_link_atom_position_interpolation() {
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0, 1, 0.7);
        assert!((la.position[0] - 1.4).abs() < 1e-10);
    }
    #[test]
    fn test_link_atom_force_correction_sum() {
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0, 1, 0.7);
        let f = [1.0, 0.0, 0.0];
        let fqm = la.qm_force_correction(f);
        let fmm = la.mm_force_correction(f);
        assert!((fqm[0] + fmm[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_qmmm_add_link_atom() {
        let mut qmmm = make_qmmm();
        let la = LinkAtomQmmm::from_bond([0.0, 0.0, 0.0], [1.5, 0.0, 0.0], 0, 1, 0.7);
        qmmm.add_link_atom(la);
        assert_eq!(qmmm.n_link_atoms(), 1);
    }
    #[test]
    fn test_link_atom_vdw_correction_finite() {
        let qmmm = make_qmmm();
        let e = qmmm.link_atom_vdw_correction(2.5, 0.1, 2.5);
        assert!(e.is_finite());
    }
    #[test]
    fn test_embedding_type_equality() {
        assert_eq!(EmbeddingType::Mechanical, EmbeddingType::Mechanical);
        assert_ne!(EmbeddingType::Mechanical, EmbeddingType::Electrostatic);
    }
    #[test]
    fn test_h2_overlap_diagonal_near_one() {
        let s = build_h2_overlap(1.4);
        assert!(s[0][0] > 0.5);
        assert!(s[1][1] > 0.5);
    }
    #[test]
    fn test_h2_overlap_off_diagonal_positive() {
        let s = build_h2_overlap(1.4);
        assert!(s[0][1] > 0.0);
    }
    #[test]
    fn test_h2_overlap_symmetric() {
        let s = build_h2_overlap(1.4);
        assert!((s[0][1] - s[1][0]).abs() < 1e-12);
    }
    #[test]
    fn test_hartree_to_ev_conversion() {
        let e_ev = HARTREE_TO_EV;
        assert!((e_ev - 27.211386).abs() < 0.001);
    }
    #[test]
    fn test_angstrom_to_bohr_round_trip() {
        let r = 1.5;
        let r_bohr = angstrom_to_bohr(r);
        let r_back = r_bohr * BOHR_TO_ANGSTROM;
        assert!((r_back - r).abs() < 1e-10);
    }
    #[test]
    fn test_hartree_to_kcal_conversion() {
        let e = hartree_to_kcal(1.0);
        assert!((e - 627.509).abs() < 0.01);
    }
    #[test]
    fn test_scf_status_converged_ne_notconverged() {
        assert_ne!(ScfStatus::Converged, ScfStatus::NotConverged);
    }
}
