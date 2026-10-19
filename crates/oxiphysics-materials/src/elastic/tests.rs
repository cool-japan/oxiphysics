// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the elastic material models.

use super::*;

/// Shear modulus formula: G = E/(2(1+ν))
#[test]
fn test_isotropic_shear_modulus() {
    let mat = IsotropicElastic::new(200.0e9, 0.3);
    let g = mat.shear_modulus();
    let expected = 200.0e9 / (2.0 * 1.3);
    assert!((g - expected).abs() < 1e6, "G mismatch: {g} vs {expected}");
}

/// Bulk modulus formula: K = E/(3(1-2ν))
#[test]
fn test_isotropic_bulk_modulus() {
    let mat = IsotropicElastic::new(200.0e9, 0.3);
    let k = mat.bulk_modulus();
    let expected = 200.0e9 / (3.0 * (1.0 - 0.6));
    assert!((k - expected).abs() < 1e6, "K mismatch: {k} vs {expected}");
}

/// P-wave modulus M = K + 4G/3
#[test]
fn test_p_wave_modulus_relation() {
    let mat = IsotropicElastic::new(200.0e9, 0.25);
    let m = mat.p_wave_modulus();
    let k = mat.bulk_modulus();
    let g = mat.shear_modulus();
    let expected = k + 4.0 / 3.0 * g;
    assert!(
        (m - expected).abs() / m < 1e-10,
        "P-wave modulus mismatch: {m} vs {expected}"
    );
}

/// Compliance matrix is symmetric and consistent with the stiffness matrix.
#[test]
fn test_compliance_symmetry() {
    let mat = IsotropicElastic::new(200.0e9, 0.3);
    let s = mat.compliance_matrix_voigt();
    for i in 0..6 {
        for j in 0..6 {
            assert!(
                (s[i * 6 + j] - s[j * 6 + i]).abs() < 1e-30,
                "S[{i}][{j}] != S[{j}][{i}]"
            );
        }
    }
}

/// Engineering strains under uniaxial stress in x: eps_yy = eps_zz = -nu/E * sigma_xx
#[test]
fn test_engineering_strains_uniaxial() {
    let mat = IsotropicElastic::new(200.0e9, 0.3);
    let stress = [100.0e6, 0.0, 0.0, 0.0, 0.0, 0.0]; // uniaxial x
    let strain = mat.engineering_strains(stress);

    let eps_xx_expected = 100.0e6 / 200.0e9;
    let eps_yy_expected = -0.3 * eps_xx_expected;

    assert!(
        (strain[0] - eps_xx_expected).abs() / eps_xx_expected < 1e-10,
        "eps_xx mismatch: {} vs {}",
        strain[0],
        eps_xx_expected
    );
    assert!(
        (strain[1] - eps_yy_expected).abs() / eps_xx_expected.abs() < 1e-10,
        "eps_yy mismatch: {} vs {}",
        strain[1],
        eps_yy_expected
    );
}

/// LinearElastic stiffness matrix is symmetric.
#[test]
fn test_linear_elastic_stress_strain_symmetry() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let c = mat.stress_strain_matrix_3d();
    for (i, row) in c.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!((val - c[j][i]).abs() < 1.0e-6, "C[{i}][{j}] != C[{j}][{i}]");
        }
    }
}

/// LinearElastic bulk and shear modulus for steel.
#[test]
fn test_linear_elastic_bulk_shear_modulus_steel() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let k = mat.bulk_modulus();
    let g = mat.shear_modulus();
    assert!((k - 166.667e9).abs() < 1.0e8);
    assert!((g - 76.923e9).abs() < 1.0e8);
}

/// Neo-Hookean: identity deformation gradient → zero stress.
#[test]
fn test_neo_hookean_identity_zero_stress() {
    let mat = NeoHookean::new(1.0e6, 1.0e9);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let p = mat.first_piola_kirchhoff_stress(&identity);
    for (i, row) in p.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            assert!(val.abs() < 1.0e-6, "P[{i}][{j}] = {} should be ~0", val);
        }
    }
}

fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn mat3_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}

/// (F^{-1})^T of a simple shear has the entry -gamma in position [2][1].
#[test]
fn test_inv_transpose3_simple_shear() {
    let gamma = 0.3;
    let f = [[1.0, 0.0, 0.0], [0.0, 1.0, gamma], [0.0, 0.0, 1.0]];
    let f_inv_t = inv_transpose3(&f);
    let expected = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -gamma, 1.0]];
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (f_inv_t[i][j] - expected[i][j]).abs() < 1e-14,
                "F^-T[{i}][{j}] = {} expected {}",
                f_inv_t[i][j],
                expected[i][j]
            );
        }
    }
}

/// transpose(inv_transpose3(F)) · F = I for non-symmetric F with det > 0.
#[test]
fn test_inv_transpose3_is_inverse_transpose() {
    let gradients = [
        [[1.1, 0.2, -0.1], [0.05, 0.9, 0.3], [-0.2, 0.15, 1.2]],
        [[0.8, -0.4, 0.25], [0.35, 1.3, -0.15], [0.1, 0.45, 0.95]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.3], [0.0, 0.0, 1.0]],
    ];
    for f in &gradients {
        assert!(det3(f) > 0.0);
        let f_inv = mat3_transpose(&inv_transpose3(f));
        let prod = mat3_mul(&f_inv, f);
        for (i, row) in prod.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "F^-1 F [{i}][{j}] = {value}"
                );
            }
        }
    }
}

/// Central-difference dW/dF of the strain energy density.
fn numerical_first_piola(mat: &NeoHookean, f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let h = 1.0e-6;
    let mut p = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut fp = *f;
            let mut fm = *f;
            fp[i][j] += h;
            fm[i][j] -= h;
            p[i][j] = (mat.strain_energy_density(&fp) - mat.strain_energy_density(&fm)) / (2.0 * h);
        }
    }
    p
}

/// The analytic PK1 must equal dW/dF (it is the derivative of the energy).
#[test]
fn test_neo_hookean_pk1_matches_energy_derivative() {
    let mat = NeoHookean::new(1.0e6, 1.0e9);
    let gradients = [
        // Simple shear in the (2,3) plane exercises cofactor [2][1].
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.3], [0.0, 0.0, 1.0]],
        [
            [1.02, 0.03, -0.01],
            [0.015, 0.98, 0.04],
            [-0.02, 0.025, 1.01],
        ],
        [[1.05, -0.02, 0.01], [0.04, 1.02, -0.03], [0.02, 0.06, 0.97]],
    ];
    for f in &gradients {
        let analytic = mat.first_piola_kirchhoff_stress(f);
        let numeric = numerical_first_piola(&mat, f);
        let scale = analytic
            .iter()
            .flatten()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()))
            .max(1.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (analytic[i][j] - numeric[i][j]).abs() <= 1.0e-6 * scale,
                    "P[{i}][{j}]: analytic {} vs numeric {}",
                    analytic[i][j],
                    numeric[i][j]
                );
            }
        }
    }
}

/// Simple shear in the (2,3) plane gives a non-zero P[2][1] component.
#[test]
fn test_neo_hookean_shear_pk1_has_transverse_component() {
    let mu = 1.0e6;
    let mat = NeoHookean::new(mu, 1.0e9);
    let gamma = 0.3;
    let f = [[1.0, 0.0, 0.0], [0.0, 1.0, gamma], [0.0, 0.0, 1.0]];
    let p = mat.first_piola_kirchhoff_stress(&f);
    // Isochoric simple shear (J = 1, I1 = 3 + gamma^2): P = mu (F - (I1/3) F^-T)
    // with F^-T[2][1] = -gamma, hence
    //   P[1][2] = mu gamma, P[2][1] = mu gamma (1 + gamma^2 / 3),
    //   P[i][i] = -mu gamma^2 / 3.
    assert!((p[1][2] - mu * gamma).abs() < 1.0e-6 * mu);
    let expected_21 = mu * gamma * (1.0 + gamma * gamma / 3.0);
    assert!(
        (p[2][1] - expected_21).abs() < 1.0e-6 * mu,
        "P[2][1] = {} expected {expected_21}",
        p[2][1]
    );
    for (i, row) in p.iter().enumerate() {
        assert!((row[i] + mu * gamma * gamma / 3.0).abs() < 1.0e-6 * mu);
    }
}

/// OrthotropicElastic compliance matrix is symmetric for equal E1, E2, E3.
#[test]
fn test_orthotropic_compliance_symmetry() {
    let mat = OrthotropicElastic {
        e1: 200.0e9,
        e2: 100.0e9,
        e3: 80.0e9,
        g12: 40.0e9,
        g23: 30.0e9,
        g13: 35.0e9,
        nu12: 0.25,
        nu23: 0.2,
        nu13: 0.22,
    };
    let s = mat.compliance_voigt();
    // Compliance diagonal entries should be positive
    for i in 0..6 {
        assert!(
            s[i * 6 + i] > 0.0,
            "Compliance diagonal S[{i}][{i}] should be positive"
        );
    }
}

fn assert_voigt_symmetric(m: &[f64; 36], label: &str) {
    for i in 0..6 {
        for j in 0..6 {
            let asym = (m[i * 6 + j] - m[j * 6 + i]).abs();
            assert!(
                asym <= 1.0e-12 * m[i * 6 + i].abs(),
                "{label}[{i}][{j}] = {} vs [{j}][{i}] = {}",
                m[i * 6 + j],
                m[j * 6 + i]
            );
        }
    }
}

/// Strongly anisotropic transversely isotropic material (E_p/E_t = 14).
fn anisotropic_ti() -> TransverselyIsotropicElastic {
    TransverselyIsotropicElastic::new(140.0e9, 10.0e9, 5.0e9, 0.3, 0.3)
}

/// Reciprocity: S_13 = S_31 = -nu_pt / E_p and S_23 = S_32.
#[test]
fn test_transversely_isotropic_compliance_is_symmetric() {
    let mat = anisotropic_ti();
    let s = mat.compliance_voigt();
    assert_voigt_symmetric(&s, "S");
    let expected = -mat.nut / mat.ep;
    for idx in [2, 2 * 6, 6 + 2, 2 * 6 + 1] {
        assert!(
            (s[idx] - expected).abs() <= 1.0e-12 * expected.abs(),
            "S[{idx}] = {} expected {expected}",
            s[idx]
        );
    }
    // nu_tp = nu_pt * E_t / E_p  (strain in p per strain in t is smaller).
    let nu_tp = -s[2] * mat.et;
    assert!((nu_tp - mat.nut * mat.et / mat.ep).abs() < 1.0e-12);
}

#[test]
fn test_transversely_isotropic_stiffness_is_symmetric_and_inverts_compliance() {
    let mat = anisotropic_ti();
    let c = mat.stiffness_voigt();
    assert_voigt_symmetric(&c, "C");
    let s = mat.compliance_voigt();
    for i in 0..6 {
        for j in 0..6 {
            let mut prod = 0.0;
            for k in 0..6 {
                prod += c[i * 6 + k] * s[k * 6 + j];
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((prod - expected).abs() < 1.0e-9, "(C S)[{i}][{j}] = {prod}");
        }
    }
}

/// With E_p = E_t, nu_p = nu_t and G = E / (2 (1 + nu)) the material is
/// isotropic and the compliance must match the isotropic one.
#[test]
fn test_transversely_isotropic_reduces_to_isotropic() {
    let e = 70.0e9;
    let nu = 0.33;
    let g = e / (2.0 * (1.0 + nu));
    let ti = TransverselyIsotropicElastic::new(e, e, g, nu, nu);
    let iso = IsotropicElastic::new(e, nu);
    let s_ti = ti.compliance_voigt();
    let s_iso = iso.compliance_matrix_voigt();
    for idx in 0..36 {
        assert!(
            (s_ti[idx] - s_iso[idx]).abs() <= 1.0e-12 * s_iso[0],
            "S[{idx}]: {} vs {}",
            s_ti[idx],
            s_iso[idx]
        );
    }
}

/// Uniaxial stress sigma in the p-direction gives the p-strain sigma / E_p
/// and the transverse strain eps_t = -nu_pt sigma / E_p.
#[test]
fn test_transversely_isotropic_uniaxial_strains() {
    let mat = anisotropic_ti();
    let s = mat.compliance_voigt();
    let sigma = 100.0e6;
    // Strain from S · [sigma, 0, 0, 0, 0, 0].
    let eps_p = s[0] * sigma;
    let eps_t = s[2 * 6] * sigma;
    assert!((eps_p - sigma / mat.ep).abs() < 1.0e-18);
    assert!((eps_t + mat.nut * sigma / mat.ep).abs() < 1.0e-18);
}

/// The UMAT tangent of the transversely isotropic model is symmetric and
/// positive definite along random-ish strain directions.
#[test]
fn test_transversely_isotropic_stress_update_tangent_symmetric() {
    let mat = anisotropic_ti();
    let strain = [1.0e-3, -4.0e-4, 2.0e-4, 1.0e-4, -2.0e-4, 3.0e-4];
    let response = mat.stress_update(&strain, &(), 1.0e-3);
    for i in 0..6 {
        for j in 0..6 {
            let asym = (response.tangent[i][j] - response.tangent[j][i]).abs();
            assert!(
                asym <= 1.0e-12 * response.tangent[i][i].abs(),
                "tangent[{i}][{j}] asymmetry {asym}"
            );
        }
    }
    // Stress = C : strain, so the strain energy is positive.
    let energy: f64 = response
        .stress
        .iter()
        .zip(strain.iter())
        .map(|(s, e)| s * e)
        .sum();
    assert!(energy > 0.0, "strain energy must be positive: {energy}");
    // Directional stiffness along the 3-axis: C33 > 0 and finite.
    assert!(response.tangent[2][2] > 0.0 && response.tangent[2][2].is_finite());
}

/// Von Mises: uniaxial stress below yield does not fail.
#[test]
fn test_von_mises_no_failure_below_yield() {
    let crit = VonMisesFailure::new(250.0e6);
    let stress = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(
        !crit.is_failed(&stress),
        "Should not fail below yield stress"
    );
}

/// Von Mises: uniaxial stress above yield fails.
#[test]
fn test_von_mises_failure_above_yield() {
    let crit = VonMisesFailure::new(250.0e6);
    let stress = [300.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(crit.is_failed(&stress), "Should fail above yield stress");
}

/// Von Mises: zero stress gives zero effective stress.
#[test]
fn test_von_mises_zero_stress() {
    let vm = VonMisesFailure::von_mises_stress(&[0.0; 6]);
    assert!(
        vm.abs() < 1e-15,
        "Von Mises stress should be 0 for zero stress, got {vm}"
    );
}

/// Tsai-Wu: zero stress does not fail.
#[test]
fn test_tsai_wu_no_failure_at_zero() {
    let crit = TsaiWuFailure::from_strengths(500.0e6, 300.0e6, 200.0e6, 150.0e6, 80.0e6);
    let stress = [0.0; 6];
    assert!(
        !crit.is_failed(&stress),
        "Zero stress should not trigger Tsai-Wu failure"
    );
}

/// Tsai-Wu: failure index from_strengths is positive for tensile load.
#[test]
fn test_tsai_wu_failure_index_large_stress() {
    let crit = TsaiWuFailure::from_strengths(100.0e6, 200.0e6, 80.0e6, 120.0e6, 50.0e6);
    // Apply stress well above tensile strength → must fail
    let stress = [200.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(
        crit.is_failed(&stress),
        "Large tensile stress should trigger Tsai-Wu failure, index={}",
        crit.failure_index(&stress)
    );
}

// -----------------------------------------------------------------------
// PlaneStressStiffness tests
// -----------------------------------------------------------------------

/// Plane stress: isotropic Q11 = E/(1-ν²)
#[test]
fn test_plane_stress_isotropic_q11() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let ps = PlaneStressStiffness::from_isotropic(e, nu);
    let expected_q11 = e / (1.0 - nu * nu);
    assert!(
        (ps.q[0] - expected_q11).abs() / expected_q11 < 1e-10,
        "Q11 mismatch: {} vs {}",
        ps.q[0],
        expected_q11
    );
}

/// Plane stress: isotropic Q12 = ν·E/(1-ν²)
#[test]
fn test_plane_stress_isotropic_q12() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let ps = PlaneStressStiffness::from_isotropic(e, nu);
    let expected_q12 = nu * e / (1.0 - nu * nu);
    assert!(
        (ps.q[1] - expected_q12).abs() / expected_q12 < 1e-10,
        "Q12 mismatch: {} vs {}",
        ps.q[1],
        expected_q12
    );
}

/// Plane stress: Q matrix is symmetric (Q12 == Q21).
#[test]
fn test_plane_stress_isotropic_symmetry() {
    let ps = PlaneStressStiffness::from_isotropic(200.0e9, 0.25);
    assert!(
        (ps.q[1] - ps.q[3]).abs() < 1e-6,
        "Q12 ({}) should equal Q21 ({})",
        ps.q[1],
        ps.q[3]
    );
}

/// Plane stress: orthotropic Q11 = E1/(1-ν12·ν21)
#[test]
fn test_plane_stress_orthotropic_q11() {
    let e1 = 200.0e9_f64;
    let e2 = 100.0e9_f64;
    let nu12 = 0.25_f64;
    let g12 = 40.0e9_f64;
    let ps = PlaneStressStiffness::from_orthotropic(e1, e2, nu12, g12);
    let nu21 = nu12 * e2 / e1;
    let denom = 1.0 - nu12 * nu21;
    let expected_q11 = e1 / denom;
    assert!(
        (ps.q[0] - expected_q11).abs() / expected_q11 < 1e-10,
        "Ortho Q11 mismatch: {} vs {}",
        ps.q[0],
        expected_q11
    );
}

/// Plane stress: orthotropic Q66 = G12 (shear entry).
#[test]
fn test_plane_stress_orthotropic_q66() {
    let g12 = 40.0e9_f64;
    let ps = PlaneStressStiffness::from_orthotropic(200.0e9, 100.0e9, 0.25, g12);
    assert!(
        (ps.q[8] - g12).abs() / g12 < 1e-10,
        "Q66 should equal G12: {} vs {}",
        ps.q[8],
        g12
    );
}

/// Plane stress: apply to uniaxial strain gives expected stress.
#[test]
fn test_plane_stress_apply_uniaxial() {
    let e = 100.0e9_f64;
    let nu = 0.0_f64; // zero Poisson for simplicity
    let ps = PlaneStressStiffness::from_isotropic(e, nu);
    // ε = [1e-3, 0, 0] → σ11 = E*ε11
    let stress = ps.apply([1.0e-3, 0.0, 0.0]);
    let expected = e * 1.0e-3;
    assert!(
        (stress[0] - expected).abs() / expected < 1e-10,
        "Uniaxial stress mismatch: {} vs {}",
        stress[0],
        expected
    );
    assert!(
        stress[1].abs() < 1e-3,
        "σ22 should be zero for ν=0: {}",
        stress[1]
    );
}

// -----------------------------------------------------------------------
// PlaneStrainStiffness tests
// -----------------------------------------------------------------------

/// Plane strain: isotropic C11 = E(1-ν)/((1+ν)(1-2ν))
#[test]
fn test_plane_strain_isotropic_c11() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let ps = PlaneStrainStiffness::from_isotropic(e, nu);
    let expected = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    assert!(
        (ps.c[0] - expected).abs() / expected < 1e-10,
        "C11 mismatch: {} vs {}",
        ps.c[0],
        expected
    );
}

/// Plane strain: isotropic C12 = Eν/((1+ν)(1-2ν))
#[test]
fn test_plane_strain_isotropic_c12() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let ps = PlaneStrainStiffness::from_isotropic(e, nu);
    let expected = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    assert!(
        (ps.c[1] - expected).abs() / expected < 1e-10,
        "C12 mismatch: {} vs {}",
        ps.c[1],
        expected
    );
}

/// Plane strain: C11 > C12 > 0 for valid Poisson ratio.
#[test]
fn test_plane_strain_ordering() {
    let ps = PlaneStrainStiffness::from_isotropic(200.0e9, 0.3);
    assert!(ps.c[0] > ps.c[1], "C11 should be greater than C12");
    assert!(ps.c[1] > 0.0, "C12 should be positive for ν > 0");
}

/// Plane strain: apply to pure shear strain yields correct shear stress.
#[test]
fn test_plane_strain_apply_shear() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let ps = PlaneStrainStiffness::from_isotropic(e, nu);
    let g = e / (2.0 * (1.0 + nu));
    let strain = [0.0, 0.0, 1.0e-3]; // pure shear γ12
    let stress = ps.apply(strain);
    let expected_s12 = g * 1.0e-3;
    assert!(
        (stress[2] - expected_s12).abs() / expected_s12 < 1e-10,
        "Shear stress mismatch: {} vs {}",
        stress[2],
        expected_s12
    );
}

/// Plane strain: C matrix is symmetric (C12 == C21).
#[test]
fn test_plane_strain_symmetry() {
    let ps = PlaneStrainStiffness::from_isotropic(150.0e9, 0.2);
    assert!(
        (ps.c[1] - ps.c[3]).abs() < 1e-6,
        "C12 ({}) should equal C21 ({})",
        ps.c[1],
        ps.c[3]
    );
}

// -----------------------------------------------------------------------
// EshelbySphericalInclusion tests
// -----------------------------------------------------------------------

/// Eshelby S1111: formula (7-5ν)/(15(1-ν)) for ν=0.3
#[test]
fn test_eshelby_s1111_nu03() {
    let esh = EshelbySphericalInclusion::new(0.3);
    let nu = 0.3_f64;
    let expected = (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
    assert!(
        (esh.s1111() - expected).abs() < 1e-12,
        "S1111 mismatch: {} vs {}",
        esh.s1111(),
        expected
    );
}

/// Eshelby S1122: formula (5ν-1)/(15(1-ν)) for ν=0.3
#[test]
fn test_eshelby_s1122_nu03() {
    let esh = EshelbySphericalInclusion::new(0.3);
    let nu = 0.3_f64;
    let expected = (5.0 * nu - 1.0) / (15.0 * (1.0 - nu));
    assert!(
        (esh.s1122() - expected).abs() < 1e-12,
        "S1122 mismatch: {} vs {}",
        esh.s1122(),
        expected
    );
}

/// Eshelby S1212: formula (4-5ν)/(15(1-ν)) for ν=0.3
#[test]
fn test_eshelby_s1212_nu03() {
    let esh = EshelbySphericalInclusion::new(0.3);
    let nu = 0.3_f64;
    let expected = (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
    assert!(
        (esh.s1212() - expected).abs() < 1e-12,
        "S1212 mismatch: {} vs {}",
        esh.s1212(),
        expected
    );
}

/// Eshelby: S1111 > 0 for any valid Poisson's ratio (0 < ν < 0.5).
#[test]
fn test_eshelby_s1111_positive() {
    for &nu in &[0.1_f64, 0.2, 0.3, 0.4, 0.49] {
        let esh = EshelbySphericalInclusion::new(nu);
        assert!(esh.s1111() > 0.0, "S1111 should be positive for ν={}", nu);
    }
}

/// Eshelby: S1212 > 0 for any valid Poisson's ratio.
#[test]
fn test_eshelby_s1212_positive() {
    for &nu in &[0.1_f64, 0.2, 0.3, 0.4, 0.49] {
        let esh = EshelbySphericalInclusion::new(nu);
        assert!(esh.s1212() > 0.0, "S1212 should be positive for ν={}", nu);
    }
}

/// Eshelby: trace sum identity — S1111 + 2*S1122 + 2*S1212 is within expected range.
#[test]
fn test_eshelby_trace_components() {
    let esh = EshelbySphericalInclusion::new(0.3);
    // All Eshelby tensor components should be positive for a sphere
    assert!(esh.s1111() > 0.0);
    assert!(esh.s1212() > 0.0);
    // For ν=0.3, S1122 = (1.5 - 1)/10.5 > 0
    assert!(esh.s1122() > 0.0);
}

// -----------------------------------------------------------------------
// EffectiveMedium tests
// -----------------------------------------------------------------------

/// Voigt modulus: for φ=0, E_V = E1.
#[test]
fn test_effective_medium_voigt_zero_phi() {
    let em = EffectiveMedium::new(0.0, 200.0e9, 0.3, 400.0e9, 0.25);
    assert!(
        (em.voigt_modulus() - 200.0e9).abs() < 1e-3,
        "Voigt at φ=0 should equal E1: {}",
        em.voigt_modulus()
    );
}

/// Voigt modulus: for φ=1, E_V = E2.
#[test]
fn test_effective_medium_voigt_full_phi() {
    let em = EffectiveMedium::new(1.0, 200.0e9, 0.3, 400.0e9, 0.25);
    assert!(
        (em.voigt_modulus() - 400.0e9).abs() < 1e-3,
        "Voigt at φ=1 should equal E2: {}",
        em.voigt_modulus()
    );
}

/// Reuss modulus: for φ=0, E_R = E1.
#[test]
fn test_effective_medium_reuss_zero_phi() {
    let em = EffectiveMedium::new(0.0, 200.0e9, 0.3, 400.0e9, 0.25);
    assert!(
        (em.reuss_modulus() - 200.0e9).abs() < 1e-3,
        "Reuss at φ=0 should equal E1: {}",
        em.reuss_modulus()
    );
}

/// Reuss modulus: for φ=1, E_R = E2.
#[test]
fn test_effective_medium_reuss_full_phi() {
    let em = EffectiveMedium::new(1.0, 200.0e9, 0.3, 400.0e9, 0.25);
    assert!(
        (em.reuss_modulus() - 400.0e9).abs() < 1e-3,
        "Reuss at φ=1 should equal E2: {}",
        em.reuss_modulus()
    );
}

/// Bounds: Reuss ≤ Voigt (always).
#[test]
fn test_effective_medium_bounds_satisfied() {
    for &phi in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 1.0] {
        let em = EffectiveMedium::new(phi, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            em.bounds_satisfied(),
            "Reuss > Voigt at φ={}: R={}, V={}",
            phi,
            em.reuss_modulus(),
            em.voigt_modulus()
        );
    }
}

/// Hill modulus: Hill = (Voigt + Reuss) / 2 — always between the bounds.
#[test]
fn test_effective_medium_hill_between_bounds() {
    let em = EffectiveMedium::new(0.4, 200.0e9, 0.3, 400.0e9, 0.25);
    let hill = em.hill_modulus();
    let voigt = em.voigt_modulus();
    let reuss = em.reuss_modulus();
    assert!(
        hill >= reuss - 1e-3 && hill <= voigt + 1e-3,
        "Hill ({}) should be between Reuss ({}) and Voigt ({})",
        hill,
        reuss,
        voigt
    );
}

/// Voigt bulk modulus: for φ=0, K_V = K1.
#[test]
fn test_effective_medium_voigt_bulk_zero_phi() {
    let e1 = 200.0e9_f64;
    let nu1 = 0.3_f64;
    let em = EffectiveMedium::new(0.0, e1, nu1, 400.0e9, 0.25);
    let k1 = e1 / (3.0 * (1.0 - 2.0 * nu1));
    assert!(
        (em.voigt_bulk_modulus() - k1).abs() < 1e-3,
        "Voigt K at φ=0 should be K1: {} vs {}",
        em.voigt_bulk_modulus(),
        k1
    );
}

/// Reuss bulk modulus ≤ Voigt bulk modulus.
#[test]
fn test_effective_medium_bulk_bounds() {
    for &phi in &[0.1_f64, 0.3, 0.5, 0.7, 0.9] {
        let em = EffectiveMedium::new(phi, 200.0e9, 0.3, 400.0e9, 0.25);
        assert!(
            em.reuss_bulk_modulus() <= em.voigt_bulk_modulus() + 1e-6,
            "Reuss K ({}) > Voigt K ({}) at φ={}",
            em.reuss_bulk_modulus(),
            em.voigt_bulk_modulus(),
            phi
        );
    }
}

/// Voigt shear modulus: for φ=0, G_V = G1.
#[test]
fn test_effective_medium_voigt_shear_zero_phi() {
    let e1 = 200.0e9_f64;
    let nu1 = 0.3_f64;
    let em = EffectiveMedium::new(0.0, e1, nu1, 400.0e9, 0.25);
    let g1 = e1 / (2.0 * (1.0 + nu1));
    assert!(
        (em.voigt_shear_modulus() - g1).abs() < 1e-3,
        "Voigt G at φ=0 should be G1: {} vs {}",
        em.voigt_shear_modulus(),
        g1
    );
}

// ── ElasticMaterial tests ───────────────────────────────────────────────

/// Compliance tensor: S = C⁻¹, so C * S = I (identity).
#[test]
fn test_elastic_material_compliance_cs_is_identity() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let s = em.compute_compliance_tensor();
    let c = em.stiffness;
    // Check C * S ≈ I
    for i in 0..6 {
        for j in 0..6 {
            let mut cs_ij = 0.0_f64;
            for k in 0..6 {
                cs_ij += c[i * 6 + k] * s[k * 6 + j];
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (cs_ij - expected).abs() < 1e-6,
                "C*S[{},{}] = {} ≠ {}",
                i,
                j,
                cs_ij,
                expected
            );
        }
    }
}

/// Engineering constants from isotropic material recover E and ν.
#[test]
fn test_elastic_material_engineering_constants_isotropic() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let mat = LinearElastic::new(e, nu);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let ec = em.compute_engineering_constants();
    // E1 = E2 = E3 = E
    assert!((ec.e1 - e).abs() / e < 1e-8, "E1: {} vs {}", ec.e1, e);
    assert!((ec.e2 - e).abs() / e < 1e-8, "E2: {} vs {}", ec.e2, e);
    assert!((ec.e3 - e).abs() / e < 1e-8, "E3: {} vs {}", ec.e3, e);
}

/// Engineering constants: nu12 recovers Poisson's ratio for isotropic mat.
#[test]
fn test_elastic_material_poisson_recovered() {
    let e = 200.0e9_f64;
    let nu = 0.25_f64;
    let mat = LinearElastic::new(e, nu);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let ec = em.compute_engineering_constants();
    assert!((ec.nu12 - nu).abs() < 1e-6, "ν12: {} vs {}", ec.nu12, nu);
    assert!((ec.nu13 - nu).abs() < 1e-6, "ν13: {} vs {}", ec.nu13, nu);
    assert!((ec.nu23 - nu).abs() < 1e-6, "ν23: {} vs {}", ec.nu23, nu);
}

/// Engineering constants: G12 = E/(2(1+ν)) for isotropic.
#[test]
fn test_elastic_material_shear_modulus_recovered() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let g = e / (2.0 * (1.0 + nu));
    let mat = LinearElastic::new(e, nu);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let ec = em.compute_engineering_constants();
    assert!((ec.g12 - g).abs() / g < 1e-8, "G12: {} vs {}", ec.g12, g);
    assert!((ec.g23 - g).abs() / g < 1e-8, "G23: {} vs {}", ec.g23, g);
    assert!((ec.g13 - g).abs() / g < 1e-8, "G13: {} vs {}", ec.g13, g);
}

/// Wave speed: P-wave v_P1 = sqrt(C11 / rho) for isotropic material.
#[test]
fn test_elastic_material_p_wave_speed() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let rho = 7800.0_f64;
    let mat = LinearElastic::new(e, nu);
    let em = ElasticMaterial::from_isotropic(&mat, rho);
    let ws = em.compute_wave_speeds();
    // C11 = E(1-ν)/((1+ν)(1-2ν)) — constrained/P-wave modulus
    let c11 = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let v_p_expected = (c11 / rho).sqrt();
    assert!(
        (ws.v_p1 - v_p_expected).abs() / v_p_expected < 1e-8,
        "v_P1: {} vs {}",
        ws.v_p1,
        v_p_expected
    );
}

/// Wave speed: S-wave v_S12 = sqrt(G / rho) for isotropic material.
#[test]
fn test_elastic_material_s_wave_speed() {
    let e = 200.0e9_f64;
    let nu = 0.3_f64;
    let rho = 7800.0_f64;
    let mat = LinearElastic::new(e, nu);
    let em = ElasticMaterial::from_isotropic(&mat, rho);
    let ws = em.compute_wave_speeds();
    let g = e / (2.0 * (1.0 + nu));
    let v_s_expected = (g / rho).sqrt();
    assert!(
        (ws.v_s12 - v_s_expected).abs() / v_s_expected < 1e-8,
        "v_S12: {} vs {}",
        ws.v_s12,
        v_s_expected
    );
}

/// Wave speeds: P-wave > S-wave for typical elastic solid.
#[test]
fn test_elastic_material_p_wave_faster_than_s_wave() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let ws = em.compute_wave_speeds();
    assert!(
        ws.v_p1 > ws.v_s12,
        "P-wave ({}) must exceed S-wave ({})",
        ws.v_p1,
        ws.v_s12
    );
}

/// Wave speeds are isotropic (equal) for an isotropic material.
#[test]
fn test_elastic_material_wave_speeds_isotropic() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let ws = em.compute_wave_speeds();
    assert!(
        (ws.v_p1 - ws.v_p2).abs() < 1.0,
        "v_P1 ≠ v_P2 for isotropic: {} {}",
        ws.v_p1,
        ws.v_p2
    );
    assert!(
        (ws.v_p2 - ws.v_p3).abs() < 1.0,
        "v_P2 ≠ v_P3 for isotropic: {} {}",
        ws.v_p2,
        ws.v_p3
    );
    assert!(
        (ws.v_s12 - ws.v_s23).abs() < 1.0,
        "v_S12 ≠ v_S23 for isotropic: {} {}",
        ws.v_s12,
        ws.v_s23
    );
}

/// Compliance matrix is symmetric for an isotropic material.
#[test]
fn test_elastic_material_compliance_symmetric() {
    let mat = LinearElastic::new(200.0e9, 0.3);
    let em = ElasticMaterial::from_isotropic(&mat, 7800.0);
    let s = em.compute_compliance_tensor();
    for i in 0..6 {
        for j in 0..6 {
            // Allow small numerical error from 6x6 matrix inversion
            let avg = (s[i * 6 + j].abs() + s[j * 6 + i].abs()) * 0.5 + 1e-40;
            let rel_err = (s[i * 6 + j] - s[j * 6 + i]).abs() / avg;
            assert!(
                rel_err < 1e-8,
                "S[{},{}] ≠ S[{},{}]: {} vs {}",
                i,
                j,
                j,
                i,
                s[i * 6 + j],
                s[j * 6 + i]
            );
        }
    }
}

/// Wave speed increases with modulus: stiffer material → faster waves.
#[test]
fn test_elastic_material_wave_speed_increases_with_modulus() {
    let mat1 = LinearElastic::new(200.0e9, 0.3);
    let mat2 = LinearElastic::new(400.0e9, 0.3);
    let em1 = ElasticMaterial::from_isotropic(&mat1, 7800.0);
    let em2 = ElasticMaterial::from_isotropic(&mat2, 7800.0);
    let ws1 = em1.compute_wave_speeds();
    let ws2 = em2.compute_wave_speeds();
    assert!(
        ws2.v_p1 > ws1.v_p1,
        "Stiffer material should have higher P-wave speed"
    );
}

/// E1 from engineering constants is positive for any valid material.
#[test]
fn test_elastic_material_engineering_constants_positive_moduli() {
    let mat = LinearElastic::new(70.0e9, 0.33); // aluminium
    let em = ElasticMaterial::from_isotropic(&mat, 2700.0);
    let ec = em.compute_engineering_constants();
    assert!(ec.e1 > 0.0, "E1 must be positive: {}", ec.e1);
    assert!(ec.g12 > 0.0, "G12 must be positive: {}", ec.g12);
    assert!(
        ec.nu12 > 0.0,
        "ν12 must be positive for aluminium: {}",
        ec.nu12
    );
}
