// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! All-electron STO-3G basis set for hydrogen through argon.
//!
//! Every shell is a least-squares fit of a Slater orbital by three Gaussians
//! (Hehre, Stewart & Pople, J. Chem. Phys. 51, 2657 (1969); Hehre, Ditchfield,
//! Stewart & Pople, J. Chem. Phys. 52, 2769 (1970)).  The primitive exponents
//! below are the standard tabulated values (as distributed by the EMSL/BSE basis
//! set exchange); the contraction coefficients are universal per shell type.
//!
//! A first-row atom gets `1s, 2s, 2px, 2py, 2pz` (5 functions), a second-row
//! atom additionally `3s, 3px, 3py, 3pz` (9 functions), H and He a single `1s`.

use super::types::BasisFunction;

/// 1s contraction coefficients (for normalised primitives), exponents descending.
const COEF_1S: [f64; 3] = [0.154_328_967_3, 0.535_328_142_3, 0.444_634_542_2];
/// 2s contraction coefficients of the 2sp shell.
const COEF_2S: [f64; 3] = [-0.099_967_229_19, 0.399_512_826_1, 0.700_115_468_9];
/// 2p contraction coefficients of the 2sp shell.
const COEF_2P: [f64; 3] = [0.155_916_275, 0.607_683_718_6, 0.391_957_393_1];
/// 3s contraction coefficients of the 3sp shell.
const COEF_3S: [f64; 3] = [-0.219_620_369, 0.225_595_433_6, 0.900_398_426];
/// 3p contraction coefficients of the 3sp shell.
const COEF_3P: [f64; 3] = [0.010_587_604_29, 0.595_167_005_3, 0.462_001_012];

/// Primitive exponents of one atom: `(1s, 2sp, 3sp)`.
type AtomExponents = ([f64; 3], Option<[f64; 3]>, Option<[f64; 3]>);

/// Tabulated STO-3G primitive exponents (bohr⁻²) for Z = 1 … 18.
fn sto3g_exponents(z: u32) -> Option<AtomExponents> {
    let table: AtomExponents = match z {
        1 => ([3.425_250_91, 0.623_913_73, 0.168_855_40], None, None),
        2 => ([6.362_421_39, 1.158_923, 0.313_649_79], None, None),
        3 => (
            [16.119_575, 2.936_200_7, 0.794_650_5],
            Some([0.636_289_7, 0.147_860_1, 0.048_088_7]),
            None,
        ),
        4 => (
            [30.167_871, 5.495_115_3, 1.487_192_7],
            Some([1.314_833_1, 0.305_538_9, 0.099_370_7]),
            None,
        ),
        5 => (
            [48.791_113, 8.887_362_2, 2.405_267],
            Some([2.236_956_1, 0.519_820_5, 0.169_061_8]),
            None,
        ),
        6 => (
            [71.616_837, 13.045_096, 3.530_512_2],
            Some([2.941_249_4, 0.683_483_1, 0.222_289_9]),
            None,
        ),
        7 => (
            [99.106_169, 18.052_312, 4.885_660_2],
            Some([3.780_455_9, 0.878_496_6, 0.285_714_4]),
            None,
        ),
        8 => (
            [130.709_32, 23.808_861, 6.443_608_3],
            Some([5.033_151_3, 1.169_596_1, 0.380_389]),
            None,
        ),
        9 => (
            [166.679_13, 30.360_812, 8.216_820_7],
            Some([6.464_803_2, 1.502_281_2, 0.488_588_5]),
            None,
        ),
        10 => (
            [207.015_61, 37.708_151, 10.205_297],
            Some([8.246_315_1, 1.916_266_2, 0.623_229_3]),
            None,
        ),
        11 => (
            [250.772_43, 45.678_511, 12.362_388],
            Some([12.040_193, 2.797_881_9, 0.909_958]),
            Some([1.478_740_6, 0.412_564_9, 0.161_475_1]),
        ),
        12 => (
            [299.237_4, 54.506_47, 14.751_58],
            Some([15.121_82, 3.513_987, 1.142_857]),
            Some([1.395_448, 0.389_326, 0.152_38]),
        ),
        13 => (
            [351.421_476_7, 64.011_860_67, 17.324_107_61],
            Some([18.899_396_21, 4.391_813_233, 1.428_353_97]),
            Some([1.395_448_293, 0.389_326_531_8, 0.152_379_765_9]),
        ),
        14 => (
            [407.797_551_4, 74.280_833_05, 20.103_292_29],
            Some([23.193_656_06, 5.389_706_871, 1.752_899_952]),
            Some([1.478_740_622, 0.412_564_880_1, 0.161_475_097_9]),
        ),
        15 => (
            [468.365_637_8, 85.313_385_59, 23.089_131_56],
            Some([28.032_639_58, 6.514_182_577, 2.118_614_352]),
            Some([1.743_103_231, 0.486_321_377_1, 0.190_342_890_9]),
        ),
        16 => (
            [533.125_735_9, 97.109_518_3, 26.281_625_42],
            Some([33.329_751_73, 7.745_117_521, 2.518_952_599]),
            Some([2.029_194_274, 0.566_140_051_8, 0.221_583_379_2]),
        ),
        17 => (
            [601.345_613_6, 109.535_854_2, 29.644_676_86],
            Some([38.960_418_89, 9.053_563_477, 2.944_499_834]),
            Some([2.129_386_495, 0.594_093_427_4, 0.232_524_141]),
        ),
        18 => (
            [674.446_518_4, 122.851_275_3, 33.248_349_45],
            Some([45.164_243_92, 10.495_199, 3.413_364_448]),
            Some([2.621_366_518, 0.731_354_605, 0.286_247_235_6]),
        ),
        _ => return None,
    };
    Some(table)
}

/// Whether an STO-3G basis is tabulated for atomic number `z`.
pub fn sto3g_supports_element(z: u32) -> bool {
    sto3g_exponents(z).is_some()
}

/// Append an sp shell (one s and three Cartesian p functions sharing exponents).
fn push_sp_shell(
    out: &mut Vec<BasisFunction>,
    center: [f64; 3],
    exponents: &[f64; 3],
    coef_s: &[f64; 3],
    coef_p: &[f64; 3],
) {
    out.push(BasisFunction::contracted(
        center,
        [0, 0, 0],
        exponents,
        coef_s,
    ));
    for cart in [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
        out.push(BasisFunction::contracted(center, cart, exponents, coef_p));
    }
}

/// All-electron STO-3G basis functions of an atom with atomic number `z`
/// centred at `center` \[Bohr\], ordered `1s, 2s, 2px, 2py, 2pz, 3s, 3px, 3py, 3pz`.
///
/// Returns `None` for elements beyond argon.
pub fn sto3g_atom_basis(z: u32, center: [f64; 3]) -> Option<Vec<BasisFunction>> {
    let (core, shell2, shell3) = sto3g_exponents(z)?;
    let mut out = Vec::with_capacity(9);
    out.push(BasisFunction::contracted(
        center,
        [0, 0, 0],
        &core,
        &COEF_1S,
    ));
    if let Some(exps) = shell2 {
        push_sp_shell(&mut out, center, &exps, &COEF_2S, &COEF_2P);
    }
    if let Some(exps) = shell3 {
        push_sp_shell(&mut out, center, &exps, &COEF_3S, &COEF_3P);
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shell_counts_follow_the_periodic_table() {
        assert_eq!(sto3g_atom_basis(1, [0.0; 3]).map(|b| b.len()), Some(1));
        assert_eq!(sto3g_atom_basis(2, [0.0; 3]).map(|b| b.len()), Some(1));
        assert_eq!(sto3g_atom_basis(8, [0.0; 3]).map(|b| b.len()), Some(5));
        assert_eq!(sto3g_atom_basis(16, [0.0; 3]).map(|b| b.len()), Some(9));
        assert!(sto3g_atom_basis(19, [0.0; 3]).is_none());
        assert!(sto3g_supports_element(18));
        assert!(!sto3g_supports_element(0));
    }

    #[test]
    fn every_function_is_normalised_and_p_components_are_orthogonal() {
        let basis = sto3g_atom_basis(16, [0.3, -0.2, 0.1]).expect("S basis");
        for (i, bf) in basis.iter().enumerate() {
            assert!((bf.overlap_with(bf) - 1.0).abs() < 1e-10, "function {i}");
        }
        // 2px ⟂ 2py ⟂ 2pz, 1s ⟂ 2px on the same centre.
        assert!(basis[2].overlap_with(&basis[3]).abs() < 1e-14);
        assert!(basis[3].overlap_with(&basis[4]).abs() < 1e-14);
        assert!(basis[0].overlap_with(&basis[2]).abs() < 1e-14);
        // 1s and 2s overlap (STO-3G 2s is not orthogonalised against the core).
        assert!(basis[0].overlap_with(&basis[1]).abs() > 0.1);
    }

    #[test]
    fn first_row_exponents_equal_zeta_scaled_expansions() {
        // Standard molecular ζ (Hehre et al. 1969): O ζ_1s = 7.66, ζ_2sp = 2.25.
        let (core, shell2, _) = sto3g_exponents(8).expect("O");
        assert!((core[0] - 2.227_660_584 * 7.66 * 7.66).abs() < 1e-3);
        let sp = shell2.expect("O 2sp");
        assert!((sp[0] - 0.994_202_729 * 2.25 * 2.25).abs() < 1e-4);
    }
}
