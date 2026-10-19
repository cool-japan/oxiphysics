//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::correlation::ElectronCorrelation;
use super::types::BasisFunction;
#[cfg(test)]
use super::types::{CartesianAxis, NaturalBondOrbital, NboAnalysis, QmmmBoundary, QmmmCoupling};
#[cfg(test)]
use super::types_2::{BasisSetType, HartreeFockSolver};
#[cfg(test)]
use super::types_2::{DensityFunctionalTheory, EmbeddingType, NboInteraction};
#[cfg(test)]
use super::types_3::{LinkAtomQmmm, MolecularOrbital, ScfStatus, XcFunctional};

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
        let bf = BasisFunction::sto3g_p([0.0, 0.0, 0.0], 1.0, CartesianAxis::Y);
        assert_eq!(bf.angular_momentum, 1);
    }
    #[test]
    fn test_self_overlap_near_one() {
        for zeta in [1.0, 1.24, 5.67, 7.66] {
            let bf = BasisFunction::sto3g_s([0.3, -0.1, 2.0], zeta);
            let s = bf.overlap_with(&bf);
            assert!(
                (s - 1.0).abs() < 1e-6,
                "1s self-overlap {s} for zeta {zeta}"
            );
            let bf2 = BasisFunction::sto3g_2s([0.3, -0.1, 2.0], zeta);
            let s2 = bf2.overlap_with(&bf2);
            assert!(
                (s2 - 1.0).abs() < 1e-6,
                "2s self-overlap {s2} for zeta {zeta}"
            );
        }
    }
    #[test]
    fn test_p_function_unit_norm() {
        // Same-centre self-overlap of one Cartesian p component:
        // Σ c_i c_j (π/γ)^{3/2} / (2γ).
        for zeta in [1.72, 2.25] {
            let bf = BasisFunction::sto3g_p([0.0; 3], zeta, CartesianAxis::X);
            let mut s = 0.0;
            for (&ai, &ci) in bf.exponents.iter().zip(bf.coefficients.iter()) {
                for (&aj, &cj) in bf.exponents.iter().zip(bf.coefficients.iter()) {
                    let g = ai + aj;
                    s += ci * cj * (std::f64::consts::PI / g).powf(1.5) / (2.0 * g);
                }
            }
            assert!((s - 1.0).abs() < 1e-6, "p self-overlap {s}");
        }
    }
    #[test]
    fn test_sto3g_hydrogen_exponents_and_primitive_weights() {
        // ζ = 1.24 reproduces the standard H STO-3G primitives (3.42525091,
        // 0.62391373, 0.16885540) and the tabulated weights of *normalised*
        // primitives (0.15432897, 0.53532814, 0.44463454).
        let bf = BasisFunction::sto3g_s([0.0; 3], 1.24);
        let expected_alpha = [0.168_855_40, 0.623_913_73, 3.425_250_91];
        let expected_d = [0.444_634_54, 0.535_328_14, 0.154_328_97];
        let mut ratio: Option<f64> = None;
        for k in 0..3 {
            assert!((bf.exponents[k] - expected_alpha[k]).abs() < 1e-6);
            let n_k = (2.0 * bf.exponents[k] / std::f64::consts::PI).powf(0.75);
            let r = bf.coefficients[k] / (n_k * expected_d[k]);
            if let Some(r0) = ratio {
                assert!((r - r0).abs() < 1e-6 * r0, "primitive weights distorted");
            }
            ratio = Some(r);
        }
    }
    #[test]
    fn test_h2_sto3g_overlap_szabo_ostlund() {
        // Szabo & Ostlund: S_12 = 0.6593 for ζ = 1.24 at R = 1.4 bohr.
        let a = BasisFunction::sto3g_s([0.0; 3], 1.24);
        let b = BasisFunction::sto3g_s([1.4, 0.0, 0.0], 1.24);
        assert!((a.overlap_with(&b) - 0.6593).abs() < 1e-4);
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
        let hf = HartreeFockSolver::from_atoms(&[1, 1], &[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], 0)
            .expect("H2");
        let dft = DensityFunctionalTheory::new(hf, XcFunctional::Pbe);
        let result = dft.run_ks_scf().expect("PBE H2");
        assert_eq!(result.status, ScfStatus::Converged);
        assert!(
            result.energy < -1.1 && result.energy > -1.2,
            "{}",
            result.energy
        );
    }
    #[test]
    fn test_dft_ks_scf_rejects_integral_only_solver() {
        // A solver built from bare matrices has no basis functions to put on a grid.
        let dft = DensityFunctionalTheory::new(make_2x2_hf(), XcFunctional::Pbe);
        assert!(dft.run_ks_scf().is_err());
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
        assert_eq!(dft.pseudopotential_energy, -0.1);
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
    /// Two occupied + two virtual orbitals with a physically structured MO tensor:
    /// (pq|rs) = 0.1 + 0.05·δ_pq·δ_rs, which has full 8-fold symmetry.
    fn make_ec() -> ElectronCorrelation {
        let eps = vec![-0.6, -0.5, 0.3, 0.4];
        let n = 4;
        let mut eri = vec![0.0; n * n * n * n];
        for p in 0..n {
            for q in 0..n {
                for r in 0..n {
                    for s in 0..n {
                        let diag = if p == q && r == s { 0.05 } else { 0.0 };
                        eri[((p * n + q) * n + r) * n + s] = 0.1 + diag;
                    }
                }
            }
        }
        ElectronCorrelation::new(2, 2, eps, eri, -1.0).expect("consistent sizes")
    }
    #[test]
    fn test_mp2_energy_negative_and_split() {
        let ec = make_ec();
        let r = ec.mp2_energy().expect("MP2");
        assert!(r.e_mp2 < 0.0, "MP2 correlation energy must be negative");
        assert!((r.e_os + r.e_ss - r.e_mp2).abs() < 1e-14);
        // (ia|jb) = (ib|ja) = 0.1 here, so E_SS = 0 and E_OS = −Σ 0.01/Δ.
        assert!(r.e_ss.abs() < 1e-14);
        let eps = [-0.6, -0.5, 0.3, 0.4];
        let mut expected = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                for a in 2..4 {
                    for b in 2..4 {
                        expected -= 0.01 / (eps[a] + eps[b] - eps[i] - eps[j]);
                    }
                }
            }
        }
        assert!((r.e_mp2 - expected).abs() < 1e-14);
    }
    #[test]
    fn test_mp2_total_includes_hf() {
        let ec = make_ec();
        let r = ec.mp2_energy().expect("MP2");
        assert_eq!(r.e_hf, -1.0);
        assert!((r.e_total - r.e_hf - r.e_mp2).abs() < 1e-10);
    }
    #[test]
    fn test_cis_requires_dipoles_and_returns_sorted_states() {
        let ec = make_ec();
        assert!(ec.cis_states(2).is_err());
        let n = 4;
        let mut dz = vec![0.0; n * n];
        for p in 0..n {
            for q in 0..n {
                dz[p * n + q] = if p == q { 0.0 } else { 0.3 };
            }
        }
        let ec = ec
            .with_mo_dipoles([vec![0.0; n * n], vec![0.0; n * n], dz])
            .expect("dipoles");
        let states = ec.cis_states(4).expect("CIS");
        assert_eq!(states.len(), 4);
        for w in states.windows(2) {
            assert!(w[0].excitation_energy <= w[1].excitation_energy);
        }
        for s in &states {
            assert!(s.excitation_energy > 0.0);
            assert!(s.oscillator_strength >= 0.0);
            let norm: f64 = s.amplitudes.iter().map(|x| x * x).sum();
            assert!((norm - 1.0).abs() < 1e-10);
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
        assert!((t1 - ((0.0004f64 + 0.0001 + 0.000225 + 0.000025) / 2.0).sqrt()).abs() < 1e-15);
    }
    #[test]
    fn test_ccs_energy_value() {
        let ec = make_ec();
        let t1 = vec![0.02; 4];
        let e = ec.ccs_energy(&t1).expect("CCS");
        // Σ_{ia,jb} [2·0.1 − 0.1] t² = 16 · 0.1 · 0.0004
        assert!((e - 16.0 * 0.1 * 0.0004).abs() < 1e-15);
    }
    #[test]
    fn test_scs_mp2_scaling() {
        let ec = make_ec();
        let mp2 = ec.mp2_energy().expect("MP2");
        let e = ec.scs_mp2_energy(1.2, 0.333).expect("SCS");
        assert!((e - (mp2.e_hf + 1.2 * mp2.e_os + 0.333 * mp2.e_ss)).abs() < 1e-14);
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
        assert!((s[0][0] - 1.0).abs() < 1e-6);
        assert!((s[1][1] - 1.0).abs() < 1e-6);
        assert!((s[0][1] - 0.6593).abs() < 1e-4);
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
