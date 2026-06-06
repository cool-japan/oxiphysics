//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Compute the trace of a 3x3 matrix.
#[inline]
pub(super) fn trace3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}
/// Determinant of a 3x3 matrix.
#[inline]
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Transpose a 3x3 matrix.
#[inline]
pub(super) fn transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Multiply two 3x3 matrices.
pub(super) fn mat_mul3(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for (r_row, a_row) in r.iter_mut().zip(a.iter()) {
        for (j, r_ij) in r_row.iter_mut().enumerate() {
            *r_ij = a_row
                .iter()
                .zip(b.iter())
                .map(|(&a_ik, b_col)| a_ik * b_col[j])
                .sum();
        }
    }
    r
}
/// Identity 3x3.
pub(super) fn identity3() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Right Cauchy-Green tensor C = F^T F.
pub(super) fn right_cauchy_green(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ft = transpose3(f);
    mat_mul3(&ft, f)
}
/// Green-Lagrange strain E = 0.5*(C - I).
pub(super) fn green_lagrange(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = right_cauchy_green(f);
    let mut e = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            e[i][j] = 0.5 * c[i][j];
            if i == j {
                e[i][j] -= 0.5;
            }
        }
    }
    e
}
/// Invariant I1 = tr(C).
pub(super) fn invariant_i1(c: &[[f64; 3]; 3]) -> f64 {
    trace3(c)
}
/// Invariant I2 = 0.5*(tr(C)^2 - tr(C^2)).
pub(super) fn invariant_i2(c: &[[f64; 3]; 3]) -> f64 {
    let tr_c = trace3(c);
    let c2 = mat_mul3(c, c);
    let tr_c2 = trace3(&c2);
    0.5 * (tr_c * tr_c - tr_c2)
}
/// Invariant I3 = det(C).
pub(super) fn invariant_i3(c: &[[f64; 3]; 3]) -> f64 {
    det3(c)
}
/// Jacobian J = det(F).
pub(super) fn jacobian(f: &[[f64; 3]; 3]) -> f64 {
    det3(f)
}
/// 3D vector dot product.
#[inline]
pub(super) fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// 3D vector norm.
#[inline]
pub(super) fn norm3(v: &[f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Scale a 3D vector.
#[inline]
pub(super) fn scale3(v: &[f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Outer product of two 3-vectors.
pub(super) fn outer3(a: &[f64; 3], b: &[f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}
/// Contract two symmetric 3x3 tensors: A:B = sum A_ij * B_ij.
pub(super) fn double_contract(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += a[i][j] * b[i][j];
        }
    }
    s
}
/// Clamp a value to \[lo, hi\].
#[inline]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}
/// Inverse of a 3x3 matrix.
pub(super) fn inv3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    if d.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let inv_d = 1.0 / d;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
        ],
    ]
}
/// Solve the cubic equation for eigenvalues of a 3x3 symmetric matrix.
///
/// x^3 - I1*x^2 + I2*x - I3 = 0
pub(super) fn solve_cubic_symmetric(i1: f64, i2: f64, i3: f64) -> [f64; 3] {
    let p = i1 * i1 - 3.0 * i2;
    let q = 2.0 * i1 * i1 * i1 - 9.0 * i1 * i2 + 27.0 * i3;
    if p.abs() < 1e-30 {
        let val = i1 / 3.0;
        return [val, val, val];
    }
    let discriminant = (q * q - 4.0 * p * p * p).max(0.0);
    if discriminant < 1e-30 {
        let r = q / (2.0 * p * p.sqrt());
        let r_clamped = clamp(r, -1.0, 1.0);
        let theta = r_clamped.acos();
        let sqrt_p = (p / 3.0).sqrt();
        let x1 = i1 / 3.0 + 2.0 * sqrt_p * (theta / 3.0).cos();
        let x2 = i1 / 3.0 + 2.0 * sqrt_p * ((theta + 2.0 * PI) / 3.0).cos();
        let x3 = i1 / 3.0 + 2.0 * sqrt_p * ((theta + 4.0 * PI) / 3.0).cos();
        let mut vals = [x1, x2, x3];
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        vals
    } else {
        let r = q / (2.0 * p * p.sqrt());
        let r_clamped = clamp(r, -1.0, 1.0);
        let theta = r_clamped.acos();
        let sqrt_p = (p / 3.0).sqrt();
        let x1 = i1 / 3.0 + 2.0 * sqrt_p * (theta / 3.0).cos();
        let x2 = i1 / 3.0 + 2.0 * sqrt_p * ((theta + 2.0 * PI) / 3.0).cos();
        let x3 = i1 / 3.0 + 2.0 * sqrt_p * ((theta + 4.0 * PI) / 3.0).cos();
        let mut vals = [x1, x2, x3];
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        vals
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::biomechanical_materials::BoneRemodeling;
    use crate::biomechanical_materials::CartilageBiphasic;
    use crate::biomechanical_materials::FungHyperelastic;
    use crate::biomechanical_materials::HillMuscleActivation;
    use crate::biomechanical_materials::HolzapfelGasserOgden;
    use crate::biomechanical_materials::MooneyRivlinBiological;
    use crate::biomechanical_materials::MultiLayerTissue;
    use crate::biomechanical_materials::QuasiLinearViscoelastic;
    use crate::biomechanical_materials::SkinMechanicsOgden;
    use crate::biomechanical_materials::TissueDamage;
    use crate::biomechanical_materials::TissueGrowth;
    pub(super) const EPS: f64 = 1e-8;
    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    fn identity_f() -> [[f64; 3]; 3] {
        identity3()
    }
    fn uniaxial_stretch(lambda: f64) -> [[f64; 3]; 3] {
        let lat = 1.0 / lambda.sqrt();
        [[lambda, 0.0, 0.0], [0.0, lat, 0.0], [0.0, 0.0, lat]]
    }
    #[test]
    fn test_fung_zero_strain_energy_at_identity() {
        let fung = FungHyperelastic::isotropic(100.0, 10.0);
        let w = fung.strain_energy(&identity_f());
        assert!(
            w.abs() < EPS,
            "Fung energy at identity should be 0, got {w}"
        );
    }
    #[test]
    fn test_fung_positive_strain_energy_under_stretch() {
        let fung = FungHyperelastic::isotropic(100.0, 10.0);
        let f = uniaxial_stretch(1.2);
        let w = fung.strain_energy(&f);
        assert!(
            w > 0.0,
            "Fung energy should be positive under stretch, got {w}"
        );
    }
    #[test]
    fn test_fung_exponential_stiffening() {
        let fung = FungHyperelastic::isotropic(100.0, 10.0);
        let w1 = fung.strain_energy(&uniaxial_stretch(1.1));
        let w2 = fung.strain_energy(&uniaxial_stretch(1.2));
        let w3 = fung.strain_energy(&uniaxial_stretch(1.3));
        let dw12 = w2 - w1;
        let dw23 = w3 - w2;
        assert!(
            dw23 > dw12,
            "Fung should stiffen exponentially: dW12={dw12:.6}, dW23={dw23:.6}"
        );
    }
    #[test]
    fn test_fung_stress_zero_at_identity() {
        let fung = FungHyperelastic::isotropic(100.0, 10.0);
        let s = fung.second_piola_kirchhoff(&identity_f());
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < EPS,
                    "Fung stress at identity should be 0, got S[{i}][{j}]={val:.6}"
                );
            }
        }
    }
    #[test]
    fn test_hgo_zero_energy_at_identity() {
        let hgo = HolzapfelGasserOgden::arterial_wall();
        let w = hgo.strain_energy(&identity_f());
        assert!(
            w.abs() < 1e-6,
            "HGO energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_hgo_positive_energy_under_stretch() {
        let hgo = HolzapfelGasserOgden::arterial_wall();
        let f = uniaxial_stretch(1.1);
        let w = hgo.strain_energy(&f);
        assert!(
            w > 0.0,
            "HGO energy should be positive under stretch, got {w}"
        );
    }
    #[test]
    fn test_hgo_anisotropy() {
        let hgo = HolzapfelGasserOgden::new(
            10.0e3,
            5.0e3,
            1.0,
            0.0,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            100.0e3,
        );
        let f_along = uniaxial_stretch(1.15);
        let w_along = hgo.strain_energy(&f_along);
        let lat = 1.0 / 1.15_f64.sqrt();
        let f_perp = [[lat, 0.0, 0.0], [0.0, 1.15, 0.0], [0.0, 0.0, lat]];
        let w_perp = hgo.strain_energy(&f_perp);
        assert!(
            w_along > w_perp,
            "HGO energy along fibers ({w_along:.6}) should exceed perpendicular ({w_perp:.6})"
        );
    }
    #[test]
    fn test_hgo_volumetric_penalty() {
        let hgo = HolzapfelGasserOgden::arterial_wall();
        let f = [[1.1, 0.0, 0.0], [0.0, 1.1, 0.0], [0.0, 0.0, 1.1]];
        let vs = hgo.volumetric_stress(&f);
        assert!(
            vs > 0.0,
            "Volumetric stress should be positive under dilation, got {vs}"
        );
    }
    #[test]
    fn test_mr_bio_zero_energy_at_identity() {
        let mr = MooneyRivlinBiological::soft_tissue();
        let w = mr.strain_energy(&identity_f());
        assert!(
            w.abs() < 1e-6,
            "MR bio energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_mr_bio_shear_modulus() {
        let mr = MooneyRivlinBiological::new(1000.0, 500.0, 0.0, 50000.0);
        let mu = mr.shear_modulus();
        assert!(
            approx_eq(mu, 3000.0, 1e-6),
            "Shear modulus should be 2*(c10+c01) = 3000, got {mu}"
        );
    }
    #[test]
    fn test_mr_bio_brain_tissue() {
        let brain = MooneyRivlinBiological::brain_tissue();
        let soft = MooneyRivlinBiological::soft_tissue();
        assert!(
            brain.shear_modulus() < soft.shear_modulus(),
            "Brain tissue should be softer than generic soft tissue"
        );
    }
    #[test]
    fn test_qlv_instantaneous_modulus() {
        let qlv = QuasiLinearViscoelastic::tendon_default();
        let ratio = qlv.instantaneous_modulus_ratio();
        assert!(
            approx_eq(ratio, 1.0, 1e-6),
            "Instantaneous modulus ratio should be 1.0, got {ratio}"
        );
    }
    #[test]
    fn test_qlv_relaxation_at_zero() {
        let qlv = QuasiLinearViscoelastic::tendon_default();
        let g0 = qlv.relaxation_function(0.0);
        assert!(
            approx_eq(g0, 1.0, 1e-6),
            "Relaxation at t=0 should be 1.0, got {g0}"
        );
    }
    #[test]
    fn test_qlv_relaxation_decays() {
        let qlv = QuasiLinearViscoelastic::tendon_default();
        let g1 = qlv.relaxation_function(1.0);
        let g10 = qlv.relaxation_function(10.0);
        let g100 = qlv.relaxation_function(100.0);
        assert!(g1 < 1.0, "G(1) should be less than 1.0");
        assert!(g10 < g1, "G(10) should be less than G(1)");
        assert!(g100 < g10, "G(100) should be less than G(10)");
    }
    #[test]
    fn test_qlv_long_term_limit() {
        let qlv = QuasiLinearViscoelastic::tendon_default();
        let g_long = qlv.relaxation_function(1e10);
        assert!(
            approx_eq(g_long, qlv.g_inf, 1e-3),
            "G(inf) should approach g_inf={:.6}, got {g_long:.6}",
            qlv.g_inf
        );
    }
    #[test]
    fn test_hill_force_length_peak() {
        let hill = HillMuscleActivation::skeletal_muscle(1000.0, 0.1);
        let fl = hill.force_length(0.1);
        assert!(
            approx_eq(fl, 1.0, 1e-6),
            "Force-length at optimal should be 1.0, got {fl}"
        );
    }
    #[test]
    fn test_hill_force_length_decreases_away() {
        let hill = HillMuscleActivation::skeletal_muscle(1000.0, 0.1);
        let fl_opt = hill.force_length(0.1);
        let fl_short = hill.force_length(0.07);
        let fl_long = hill.force_length(0.13);
        assert!(
            fl_short < fl_opt,
            "Force-length should decrease when shorter"
        );
        assert!(fl_long < fl_opt, "Force-length should decrease when longer");
    }
    #[test]
    fn test_hill_activation_dynamics() {
        let mut hill = HillMuscleActivation::skeletal_muscle(1000.0, 0.1);
        assert!(
            hill.activation().abs() < EPS,
            "Initial activation should be 0"
        );
        for _ in 0..100 {
            hill.update_activation(1.0, 0.001);
        }
        assert!(
            hill.activation() > 0.5,
            "Activation should increase, got {}",
            hill.activation()
        );
    }
    #[test]
    fn test_hill_zero_force_at_zero_activation() {
        let hill = HillMuscleActivation::skeletal_muscle(1000.0, 0.1);
        let force = hill.compute_force(0.1, 0.0);
        assert!(
            force.abs() < EPS,
            "Force at zero activation and optimal length should be ~0, got {force}"
        );
    }
    #[test]
    fn test_bone_cortical_defaults() {
        let bone = BoneRemodeling::cortical_bone();
        assert!(
            approx_eq(bone.density, 1900.0, EPS),
            "Cortical bone density should be 1900"
        );
        assert!(!bone.is_osteoporotic());
    }
    #[test]
    fn test_bone_elastic_modulus_scaling() {
        let bone = BoneRemodeling::cortical_bone();
        let e = bone.elastic_modulus();
        assert!(
            approx_eq(e, bone.e_ref, 1e-3),
            "Elastic modulus at ref density should be E_ref={:.6}, got {e:.6}",
            bone.e_ref
        );
    }
    #[test]
    fn test_bone_remodeling_increases_density() {
        let mut bone = BoneRemodeling::cortical_bone();
        let initial = bone.density;
        for _ in 0..1000 {
            bone.update_density(0.01, 0.1);
        }
        assert!(
            bone.density > initial,
            "High stimulus should increase density: initial={initial}, final={}",
            bone.density
        );
    }
    #[test]
    fn test_bone_remodeling_decreases_density() {
        let mut bone = BoneRemodeling::cortical_bone();
        let initial = bone.density;
        for _ in 0..1000 {
            bone.update_density(0.001, 0.1);
        }
        assert!(
            bone.density < initial,
            "Low stimulus should decrease density: initial={initial}, final={}",
            bone.density
        );
    }
    #[test]
    fn test_cartilage_zones_stiffness_order() {
        let sup = CartilageBiphasic::superficial_zone();
        let mid = CartilageBiphasic::middle_zone();
        let deep = CartilageBiphasic::deep_zone();
        assert!(
            sup.aggregate_modulus < mid.aggregate_modulus,
            "Superficial should be softer than middle"
        );
        assert!(
            mid.aggregate_modulus < deep.aggregate_modulus,
            "Middle should be softer than deep"
        );
    }
    #[test]
    fn test_cartilage_consolidation_time() {
        let cart = CartilageBiphasic::middle_zone();
        let tau = cart.consolidation_time(0.001);
        assert!(
            tau > 0.0,
            "Consolidation time should be positive, got {tau}"
        );
        assert!(
            tau < 1e6,
            "Consolidation time should be reasonable, got {tau}"
        );
    }
    #[test]
    fn test_cartilage_porosity() {
        let cart = CartilageBiphasic::middle_zone();
        let phi = cart.porosity();
        assert!(
            phi > 0.0 && phi < 1.0,
            "Porosity should be in (0,1), got {phi}"
        );
        assert!(
            approx_eq(phi, 1.0 - cart.solid_fraction, EPS),
            "Porosity should be 1 - solid_fraction"
        );
    }
    #[test]
    fn test_skin_zero_energy_at_identity() {
        let skin = SkinMechanicsOgden::human_skin();
        let w = skin.strain_energy_from_stretches(1.0, 1.0, 1.0);
        assert!(
            w.abs() < 1e-6,
            "Skin energy at identity should be ~0, got {w}"
        );
    }
    #[test]
    fn test_skin_positive_energy_under_stretch() {
        let skin = SkinMechanicsOgden::human_skin();
        let w = skin.strain_energy_from_stretches(1.2, 1.0 / 1.2_f64.sqrt(), 1.0 / 1.2_f64.sqrt());
        assert!(
            w > 0.0,
            "Skin energy should be positive under stretch, got {w}"
        );
    }
    #[test]
    fn test_skin_aged_softer() {
        let young = SkinMechanicsOgden::human_skin();
        let old = SkinMechanicsOgden::aged_skin();
        let mu_young = young.initial_shear_modulus();
        let mu_old = old.initial_shear_modulus();
        assert!(
            mu_old.abs() < mu_young.abs(),
            "Aged skin should be softer: young={mu_young:.6}, old={mu_old:.6}"
        );
    }
    #[test]
    fn test_multilayer_total_thickness() {
        let ml = MultiLayerTissue::skin_fat_muscle();
        let t = ml.total_thickness();
        assert!(
            approx_eq(t, 0.042, 1e-6),
            "Total thickness should be 0.042, got {t}"
        );
    }
    #[test]
    fn test_multilayer_voigt_reuss_bounds() {
        let ml = MultiLayerTissue::skin_fat_muscle();
        let e_v = ml.effective_youngs_modulus_voigt();
        let e_r = ml.effective_youngs_modulus_reuss();
        assert!(
            e_v >= e_r,
            "Voigt bound ({e_v}) should be >= Reuss bound ({e_r})"
        );
    }
    #[test]
    fn test_damage_initially_zero() {
        let td = TissueDamage::soft_tissue_default();
        assert!(td.damage.abs() < EPS, "Initial damage should be 0");
        assert!(!td.is_ruptured());
    }
    #[test]
    fn test_damage_evolves_above_threshold() {
        let mut td = TissueDamage::soft_tissue_default();
        td.update(0.5);
        assert!(
            td.damage > 0.0,
            "Damage should increase above threshold, got {}",
            td.damage
        );
    }
    #[test]
    fn test_damage_no_evolution_below_threshold() {
        let mut td = TissueDamage::soft_tissue_default();
        td.update(0.1);
        assert!(
            td.damage.abs() < EPS,
            "Damage should not evolve below threshold, got {}",
            td.damage
        );
    }
    #[test]
    fn test_growth_initial_volume_ratio() {
        let tg = TissueGrowth::arterial_growth();
        let vr = tg.growth_volume_ratio();
        assert!(
            approx_eq(vr, 1.0, EPS),
            "Initial growth volume ratio should be 1.0, got {vr}"
        );
    }
    #[test]
    fn test_growth_increases_under_high_stress() {
        let mut tg = TissueGrowth::arterial_growth();
        let high_stress = 200.0e3;
        for _ in 0..1000 {
            tg.update_isotropic(high_stress, 0.01);
        }
        let vr = tg.growth_volume_ratio();
        assert!(
            vr > 1.0,
            "High stress should cause growth: volume ratio = {vr}"
        );
    }
    #[test]
    fn test_growth_elastic_deformation() {
        let tg = TissueGrowth::arterial_growth();
        let f = identity_f();
        let fe = tg.elastic_deformation(&f);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    approx_eq(fe[i][j], f[i][j], EPS),
                    "Fe should equal F when Fg=I"
                );
            }
        }
    }
}
