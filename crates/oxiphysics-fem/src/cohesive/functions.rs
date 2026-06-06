//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::BilinearCzmParams;

/// Compute the 12×12 cohesive element stiffness matrix for a 2D (4-node)
/// cohesive interface element.
///
/// The cohesive element has 4 nodes (2 per face), each with 3 DOFs (x,y,z).
/// The element resists relative opening/sliding through the cohesive law.
///
/// The stiffness is: K_coh = ∫ Bᵀ D_coh B dA
/// where B extracts the displacement jump and D_coh is the cohesive tangent.
///
/// # Arguments
/// * `k_n` – normal cohesive stiffness (N/m³)
/// * `k_t` – tangential cohesive stiffness (N/m³)
/// * `area` – interface element area (m²)
///
/// Returns a 12×12 symmetric matrix.
pub fn cohesive_element_stiffness(k_n: f64, k_t: f64, area: f64) -> [[f64; 12]; 12] {
    let factor = area / 4.0;
    let kn = k_n * factor;
    let kt = k_t * factor;
    let d = [kn, kt, kt];
    let mut ke = [[0.0f64; 12]; 12];
    let pairs = [(0usize, 2usize), (1usize, 3usize)];
    for &(nu, nl) in &pairs {
        for (dof, &k) in d.iter().enumerate() {
            let row_u = nu * 3 + dof;
            let row_l = nl * 3 + dof;
            ke[row_u][row_u] += k;
            ke[row_l][row_l] += k;
            ke[row_u][row_l] -= k;
            ke[row_l][row_u] -= k;
        }
    }
    ke
}
/// Compute the mode-I energy release rate G_I from the traction-separation law.
///
/// G_I = area under the T_n – Δ_n curve up to complete failure.
/// For the bilinear law: G_I = ½ T_max * δ_f.
pub fn energy_release_rate_mode1(params: &BilinearCzmParams) -> f64 {
    params.fracture_toughness()
}
/// Compute the mode-II energy release rate G_II from a bilinear shear law.
///
/// `t_shear_max` – peak shear traction (Pa).
/// `delta_shear_f` – shear failure separation (m).
pub fn energy_release_rate_mode2(t_shear_max: f64, delta_shear_f: f64) -> f64 {
    0.5 * t_shear_max * delta_shear_f
}
/// Compute the mode-III energy release rate G_III from a bilinear tearing law.
pub fn energy_release_rate_mode3(t_tear_max: f64, delta_tear_f: f64) -> f64 {
    0.5 * t_tear_max * delta_tear_f
}
/// Total energy release rate for mixed-mode fracture using power law interaction.
///
/// The Benzeggagh-Kenane (B-K) criterion:
/// G_c(β) = G_Ic + (G_IIc - G_Ic) * β^η
///
/// where β = G_shear / G_total (mixed-mode ratio) and η is a material exponent.
pub fn bk_critical_energy_release_rate(g_ic: f64, g_iic: f64, beta: f64, eta: f64) -> f64 {
    g_ic + (g_iic - g_ic) * beta.clamp(0.0, 1.0).powf(eta)
}
/// Power-law mixed-mode failure criterion:
/// (G_I/G_Ic)^α + (G_II/G_IIc)^α + (G_III/G_IIIc)^α = 1
///
/// Returns the total index; failure when ≥ 1.
pub fn power_law_failure_index(
    g1: f64,
    g_ic: f64,
    g2: f64,
    g_iic: f64,
    g3: f64,
    g_iiic: f64,
    alpha: f64,
) -> f64 {
    let term1 = if g_ic > f64::EPSILON {
        (g1 / g_ic).powf(alpha)
    } else {
        0.0
    };
    let term2 = if g_iic > f64::EPSILON {
        (g2 / g_iic).powf(alpha)
    } else {
        0.0
    };
    let term3 = if g_iiic > f64::EPSILON {
        (g3 / g_iiic).powf(alpha)
    } else {
        0.0
    };
    (term1 + term2 + term3).powf(1.0 / alpha)
}
/// Estimate the cohesive zone length (Irwin-type) ahead of the crack tip.
///
/// L_cz ≈ E * G_c / (π * σ_c²)
///
/// where E is Young's modulus, G_c is the fracture toughness, and σ_c is
/// the cohesive strength.
pub fn cohesive_zone_length(young_modulus: f64, g_c: f64, sigma_c: f64) -> f64 {
    if sigma_c < f64::EPSILON {
        return f64::INFINITY;
    }
    young_modulus * g_c / (std::f64::consts::PI * sigma_c * sigma_c)
}
/// Estimate the number of elements needed to resolve the cohesive zone.
///
/// At least 3–5 elements per cohesive zone are recommended.
pub fn elements_per_cohesive_zone(cohesive_length: f64, element_size: f64) -> f64 {
    if element_size < f64::EPSILON {
        return f64::INFINITY;
    }
    cohesive_length / element_size
}
/// Compute the crack mouth opening displacement (CMOD) for a centre-cracked panel.
///
/// Analytical estimate for plane-stress / plane-strain:
/// CMOD ≈ 4 * σ * a / E * sqrt(1 - (a/W)^2)  (Tada handbook)
///
/// # Arguments
/// * `sigma` – remote stress (Pa)
/// * `a` – half crack length (m)
/// * `w` – half specimen width (m)
/// * `young_modulus` – E (Pa)
pub fn cmod_centre_crack(sigma: f64, a: f64, w: f64, young_modulus: f64) -> f64 {
    let aw = a / w.max(f64::EPSILON);
    let factor = (1.0 - aw * aw).sqrt();
    4.0 * sigma * a / young_modulus.max(f64::EPSILON) * factor
}
/// Dugdale (strip yield) cohesive zone length estimate.
///
/// For a central crack in plane stress under remote stress σ:
/// L_Dugdale = π/8 * (K_I / σ_yield)²
///
/// where K_I = σ * sqrt(π * a) and σ_yield is the yield stress.
pub fn dugdale_zone_length(sigma_remote: f64, a: f64, sigma_yield: f64) -> f64 {
    if sigma_yield < f64::EPSILON {
        return f64::INFINITY;
    }
    let k1 = sigma_remote * (std::f64::consts::PI * a).sqrt();
    std::f64::consts::PI / 8.0 * (k1 / sigma_yield).powi(2)
}
/// Barenblatt cohesive zone length (linear softening).
///
/// L_B = π * E * G_c / (2 * σ_c²)  (plane stress)
///
/// This is proportional to the Irwin cohesive zone length.
pub fn barenblatt_zone_length(young_modulus: f64, g_c: f64, sigma_c: f64) -> f64 {
    if sigma_c < f64::EPSILON {
        return f64::INFINITY;
    }
    std::f64::consts::PI * young_modulus * g_c / (2.0 * sigma_c * sigma_c)
}
/// Ratio of Dugdale to Barenblatt cohesive zone lengths.
///
/// Returns L_Dugdale / L_Barenblatt — independent of E and G_c.
pub fn dugdale_to_barenblatt_ratio(
    sigma_remote: f64,
    a: f64,
    sigma_yield: f64,
    sigma_c: f64,
) -> f64 {
    let ld = dugdale_zone_length(sigma_remote, a, sigma_yield);
    let lb = barenblatt_zone_length(1.0, sigma_remote * sigma_remote * a / (2.0), sigma_c);
    if lb < f64::EPSILON {
        return f64::INFINITY;
    }
    ld / lb
}
/// Compute the relative sliding / opening displacements for a 2-node pair
/// (1D cohesive element) from nodal displacements.
///
/// # Arguments
/// * `u_upper` – displacement of the upper face node \[ux, uy, uz\]
/// * `u_lower` – displacement of the lower face node \[ux, uy, uz\]
/// * `normal` – outward unit normal of the interface \[nx, ny, nz\]
///
/// Returns `(delta_n, delta_t1, delta_t2)`:
/// - `delta_n`  = normal opening = (u_upper - u_lower) · n
/// - `delta_t1` = in-plane tangential sliding component 1
/// - `delta_t2` = in-plane tangential sliding component 2
pub fn interface_displacement_jump(
    u_upper: [f64; 3],
    u_lower: [f64; 3],
    normal: [f64; 3],
) -> (f64, f64, f64) {
    let du = [
        u_upper[0] - u_lower[0],
        u_upper[1] - u_lower[1],
        u_upper[2] - u_lower[2],
    ];
    let delta_n = du[0] * normal[0] + du[1] * normal[1] + du[2] * normal[2];
    let (t1, t2) = compute_interface_tangents(normal);
    let delta_t1 = du[0] * t1[0] + du[1] * t1[1] + du[2] * t1[2];
    let delta_t2 = du[0] * t2[0] + du[1] * t2[1] + du[2] * t2[2];
    (delta_n, delta_t1, delta_t2)
}
/// Compute two unit tangent vectors forming an orthonormal frame with `normal`.
pub(super) fn compute_interface_tangents(normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let aux = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let t1_raw = [
        aux[1] * normal[2] - aux[2] * normal[1],
        aux[2] * normal[0] - aux[0] * normal[2],
        aux[0] * normal[1] - aux[1] * normal[0],
    ];
    let len1 = (t1_raw[0] * t1_raw[0] + t1_raw[1] * t1_raw[1] + t1_raw[2] * t1_raw[2])
        .sqrt()
        .max(f64::EPSILON);
    let t1 = [t1_raw[0] / len1, t1_raw[1] / len1, t1_raw[2] / len1];
    let t2 = [
        normal[1] * t1[2] - normal[2] * t1[1],
        normal[2] * t1[0] - normal[0] * t1[2],
        normal[0] * t1[1] - normal[1] * t1[0],
    ];
    (t1, t2)
}
/// Effective mixed-mode separation Δ_eff for use in scalar damage tracking.
///
/// Δ_eff = sqrt(max(0, Δ_n)² + β² * Δ_t²)
///
/// where β is a mode-mixity weight (β = 1: equal weighting; β = 0: mode-I only).
pub fn effective_separation(delta_n: f64, delta_t: f64, beta_mix: f64) -> f64 {
    let dn = delta_n.max(0.0);
    (dn * dn + beta_mix * beta_mix * delta_t * delta_t).sqrt()
}
/// Compute mode mixity angle ψ in degrees.
///
/// ψ = atan2(K_II, K_I) * 180/π, where K_I, K_II are stress intensity factors.
pub fn mode_mixity_angle_degrees(k1: f64, k2: f64) -> f64 {
    k2.atan2(k1) * 180.0 / std::f64::consts::PI
}
/// Polynomial stiffness degradation model.
///
/// The secant modulus degrades as:
/// E_eff = E_0 * (1 - D)^n_deg
///
/// where n_deg is the degradation exponent (n_deg = 1: linear, n_deg = 2: quadratic).
pub fn polynomial_stiffness_degradation(e0: f64, damage: f64, n_deg: f64) -> f64 {
    let d = damage.clamp(0.0, 1.0);
    e0 * (1.0 - d).powf(n_deg)
}
/// Compute the tangent modulus for a cohesive element in the softening branch.
///
/// For the bilinear law:
/// K_tangent = -T_max / (δ_f - δ_0)  (slope of the softening branch)
pub fn bilinear_tangent_modulus_softening(params: &BilinearCzmParams) -> f64 {
    let d0 = params.delta_onset();
    let df = params.delta_failure;
    let denom = (df - d0).max(f64::EPSILON);
    -params.traction_max / denom
}
/// Secant unloading stiffness for the bilinear law at separation `delta`.
///
/// K_sec = T(delta) / delta
/// In the softening branch: K_sec = T_max * (δ_f - δ) / (δ * (δ_f - δ_0))
pub fn bilinear_secant_unloading(params: &BilinearCzmParams, delta: f64) -> f64 {
    if delta < f64::EPSILON {
        return params.penalty_stiffness;
    }
    let d0 = params.delta_onset();
    let df = params.delta_failure;
    if delta <= d0 {
        params.penalty_stiffness
    } else if delta < df {
        let t = params.traction(delta);
        t / delta
    } else {
        0.0
    }
}
/// Weibull weakest-link model for interface strength.
///
/// The survival probability for an interface element of area A under traction T:
/// P_s = exp(-(A/A_0) * (T/T_0)^m_w)
///
/// where m_w is the Weibull modulus and T_0, A_0 are reference values.
pub fn weibull_survival_probability(
    traction: f64,
    area: f64,
    t_ref: f64,
    a_ref: f64,
    m_weibull: f64,
) -> f64 {
    if traction <= 0.0 {
        return 1.0;
    }
    let ratio_t = traction / t_ref.max(f64::EPSILON);
    let ratio_a = area / a_ref.max(f64::EPSILON);
    (-ratio_a * ratio_t.powf(m_weibull)).exp()
}
/// Weibull characteristic strength for a given element area and target survival probability.
///
/// T_char = T_0 * (-ln(P_s) * A_0 / A)^(1/m_w)
pub fn weibull_characteristic_strength(
    target_ps: f64,
    area: f64,
    t_ref: f64,
    a_ref: f64,
    m_weibull: f64,
) -> f64 {
    let ps = target_ps.clamp(f64::EPSILON, 1.0 - f64::EPSILON);
    let ratio_a = area / a_ref.max(f64::EPSILON);
    t_ref * (-ps.ln() * (1.0 / ratio_a)).powf(1.0 / m_weibull.max(f64::EPSILON))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohesive::*;
    #[test]
    fn test_bilinear_onset_separation() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let d0 = p.delta_onset();
        assert!((d0 - 50e6 / 1e9).abs() < 1e-12, "d0 = {d0}");
    }
    #[test]
    fn test_bilinear_traction_linear_branch() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let t = p.traction(d0 * 0.5);
        let expected = 1e9 * d0 * 0.5;
        assert!(
            (t - expected).abs() / expected < 1e-10,
            "linear branch t = {t}"
        );
    }
    #[test]
    fn test_bilinear_traction_peak() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let t = p.traction(d0 * (1.0 + 1e-10));
        assert!(
            (t - p.traction_max).abs() / p.traction_max < 1e-3,
            "peak t = {t}"
        );
    }
    #[test]
    fn test_bilinear_traction_zero_after_failure() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let t = p.traction(1e-3);
        assert!(t.abs() < 1e-12, "should be zero after failure, got {t}");
    }
    #[test]
    fn test_bilinear_fracture_toughness() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let gc = p.fracture_toughness();
        let expected = 0.5 * 50e6 * 1e-4;
        assert!((gc - expected).abs() / expected < 1e-10, "G_c = {gc}");
    }
    #[test]
    fn test_bilinear_damage_undamaged() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let d = p.damage(d0 * 0.5);
        assert!(d.abs() < 1e-12, "below onset: D = {d}");
    }
    #[test]
    fn test_bilinear_damage_fully_failed() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d = p.damage(p.delta_failure * 2.0);
        assert!((d - 1.0).abs() < 1e-12, "after failure: D = {d}");
    }
    #[test]
    fn test_bilinear_damage_mid() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let df = p.delta_failure;
        let delta_mid = 0.5 * (d0 + df);
        let d = p.damage(delta_mid);
        assert!(d > 0.0 && d < 1.0, "mid-range damage D = {d}");
    }
    #[test]
    fn test_xu_needleman_sigma_max() {
        let p = XuNeedlemanParams::new(500.0, 250.0, 1e-4, 1e-4);
        let s = p.sigma_max();
        let expected = 500.0 / (std::f64::consts::E * 1e-4);
        assert!((s - expected).abs() / expected < 1e-10, "sigma_max = {s}");
    }
    #[test]
    fn test_xu_needleman_normal_traction_at_zero() {
        let p = XuNeedlemanParams::new(500.0, 250.0, 1e-4, 1e-4);
        let t = p.normal_traction(0.0, 0.0);
        assert!(t.is_finite(), "should be finite, got {t}");
    }
    #[test]
    fn test_xu_needleman_tangential_zero_at_zero_tangential() {
        let p = XuNeedlemanParams::new(500.0, 250.0, 1e-4, 1e-4);
        let t = p.tangential_traction(0.0, 0.0);
        assert!(t.abs() < 1e-20, "tangential at zero separation = {t}");
    }
    #[test]
    fn test_xu_needleman_potential_zero_at_zero() {
        let p = XuNeedlemanParams::new(500.0, 250.0, 1e-4, 1e-4);
        let phi = p.potential(0.0, 0.0);
        assert!(phi.abs() < 1e-20, "potential at zero = {phi}");
    }
    #[test]
    fn test_mixed_mode_not_initiated_below_critical() {
        let c = MixedModeCriterion::new(50e6, 30e6, 20e6);
        assert!(!c.is_initiated(10e6, 10e6, 5e6), "below critical");
    }
    #[test]
    fn test_mixed_mode_initiated_at_critical_normal() {
        let c = MixedModeCriterion::new(50e6, 1e20, 1e20);
        assert!(
            c.is_initiated(50e6, 0.0, 0.0),
            "at critical normal traction"
        );
    }
    #[test]
    fn test_mixed_mode_compressive_not_initiated() {
        let c = MixedModeCriterion::new(50e6, 30e6, 20e6);
        assert!(
            !c.is_initiated(-100e6, 0.0, 0.0),
            "compressive normal should not initiate"
        );
    }
    #[test]
    fn test_mixed_mode_ratio_pure_modes() {
        let beta = MixedModeCriterion::mixed_mode_ratio(1.0, 0.0);
        assert!(beta.abs() < 1e-12, "pure mode I β = {beta}");
        let beta2 = MixedModeCriterion::mixed_mode_ratio(0.0, 1.0);
        assert!((beta2 - 1.0).abs() < 1e-12, "pure mode II β = {beta2}");
    }
    #[test]
    fn test_cohesive_stiffness_symmetry() {
        let ke = cohesive_element_stiffness(1e9, 5e8, 0.01);
        for (i, row) in ke.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!(
                    (v - ke[j][i]).abs() < 1e-8,
                    "not symmetric at ({i},{j}): {v} vs {}",
                    ke[j][i]
                );
            }
        }
    }
    #[test]
    fn test_cohesive_stiffness_positive_diagonal() {
        let ke = cohesive_element_stiffness(1e9, 5e8, 0.01);
        for (i, row) in ke.iter().enumerate() {
            assert!(row[i] > 0.0, "diagonal [{i}] = {}", row[i]);
        }
    }
    #[test]
    fn test_damage_state_initial() {
        let state = CohesiveDamageState::new();
        assert!((state.damage).abs() < 1e-12);
        assert!(!state.failed);
    }
    #[test]
    fn test_damage_state_update_undamaged() {
        let mut state = CohesiveDamageState::new();
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let d = state.update(p.delta_onset() * 0.5, &p);
        assert!(d.abs() < 1e-12, "undamaged D = {d}");
    }
    #[test]
    fn test_damage_state_update_damaged() {
        let mut state = CohesiveDamageState::new();
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d = state.update(p.delta_failure * 0.5, &p);
        assert!(d > 0.0 && d < 1.0, "partial damage D = {d}");
    }
    #[test]
    fn test_damage_state_irreversibility() {
        let mut state = CohesiveDamageState::new();
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d1 = state.update(p.delta_failure * 0.7, &p);
        let d2 = state.update(p.delta_failure * 0.3, &p);
        assert!(d2 >= d1, "damage must not decrease: d1={d1} d2={d2}");
    }
    #[test]
    fn test_damage_state_full_failure() {
        let mut state = CohesiveDamageState::new();
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        state.update(p.delta_failure * 2.0, &p);
        assert!(state.failed, "should be failed");
        assert!((state.damage - 1.0).abs() < 1e-12, "D = {}", state.damage);
    }
    #[test]
    fn test_damage_state_traction_after_failure() {
        let mut state = CohesiveDamageState::new();
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        state.update(p.delta_failure * 2.0, &p);
        let t = state.traction(1e-3, &p);
        assert!(t.abs() < 1e-12, "traction after failure = {t}");
    }
    #[test]
    fn test_cohesive_interface_creation() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let iface = CohesiveInterface::new(vec![0, 1, 2], vec![3, 4, 5], p);
        assert_eq!(iface.n_nodes(), 3);
        assert_eq!(iface.damage_states.len(), 3);
        assert!(!iface.is_fully_failed());
    }
    #[test]
    fn test_cohesive_interface_max_damage_zero() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let iface = CohesiveInterface::new(vec![0, 1], vec![2, 3], p);
        assert!(iface.max_damage().abs() < 1e-12);
    }
    #[test]
    fn test_energy_release_rate_mode1() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let g = energy_release_rate_mode1(&p);
        assert!(
            (g - 0.5 * 50e6 * 1e-4).abs() / (0.5 * 50e6 * 1e-4) < 1e-10,
            "G_I = {g}"
        );
    }
    #[test]
    fn test_energy_release_rate_mode2() {
        let g = energy_release_rate_mode2(30e6, 2e-4);
        assert!(
            (g - 0.5 * 30e6 * 2e-4).abs() / (0.5 * 30e6 * 2e-4) < 1e-10,
            "G_II = {g}"
        );
    }
    #[test]
    fn test_bk_criterion_pure_mode_i() {
        let gc = bk_critical_energy_release_rate(100.0, 200.0, 0.0, 2.0);
        assert!((gc - 100.0).abs() < 1e-10, "pure mode I G_c = {gc}");
    }
    #[test]
    fn test_bk_criterion_pure_mode_ii() {
        let gc = bk_critical_energy_release_rate(100.0, 200.0, 1.0, 2.0);
        assert!((gc - 200.0).abs() < 1e-10, "pure mode II G_c = {gc}");
    }
    #[test]
    fn test_power_law_failure_index() {
        let idx = power_law_failure_index(0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 2.0);
        let expected = (3.0 * 0.25_f64).sqrt();
        assert!((idx - expected).abs() < 1e-10, "idx = {idx}");
    }
    #[test]
    fn test_r_curve_initial_toughness() {
        let rc = RCurveModel::new(10.0, 30.0, 1e-3);
        let kr = rc.k_resistance(0.0);
        assert!((kr - 10.0).abs() < 1e-10, "K_R(0) = {kr}");
    }
    #[test]
    fn test_r_curve_steady_state_limit() {
        let rc = RCurveModel::new(10.0, 30.0, 1e-3);
        let kr = rc.k_resistance(1.0);
        assert!((kr - 30.0).abs() < 0.01, "K_R(∞) ≈ {kr}");
    }
    #[test]
    fn test_r_curve_propagation() {
        let rc = RCurveModel::new(10.0, 30.0, 1e-3);
        assert!(rc.propagates(50.0, 0.0), "K > K_R → propagates");
        assert!(!rc.propagates(5.0, 0.0), "K < K_R → does not propagate");
    }
    #[test]
    fn test_r_curve_stable_extension() {
        let rc = RCurveModel::new(10.0, 30.0, 1e-3);
        let result = rc.stable_crack_extension(20.0, 0.1, 1e-8);
        assert!(result.is_some(), "should find stable extension");
        let da = result.unwrap();
        let kr = rc.k_resistance(da);
        assert!((kr - 20.0).abs() < 1e-4, "K_R({da:.4e}) = {kr:.4} ≈ 20");
    }
    #[test]
    fn test_r_curve_unstable_returns_none() {
        let rc = RCurveModel::new(10.0, 30.0, 1e-3);
        let result = rc.stable_crack_extension(35.0, 1.0, 1e-8);
        assert!(result.is_none(), "K > K_ss → unstable (None)");
    }
    #[test]
    fn test_cohesive_zone_length_positive() {
        let lc = cohesive_zone_length(200e9, 500.0, 50e6);
        assert!(lc > 0.0 && lc.is_finite(), "L_cz = {lc}");
    }
    #[test]
    fn test_elements_per_cz() {
        let lc = cohesive_zone_length(200e9, 500.0, 50e6);
        let n = elements_per_cohesive_zone(lc, lc / 5.0);
        assert!((n - 5.0).abs() < 1e-10, "n_el = {n}");
    }
    #[test]
    fn test_cmod_centre_crack_positive() {
        let cmod = cmod_centre_crack(100e6, 0.05, 0.1, 200e9);
        assert!(cmod > 0.0 && cmod.is_finite(), "CMOD = {cmod}");
    }
    #[test]
    fn test_cmod_proportional_to_stress() {
        let c1 = cmod_centre_crack(100e6, 0.05, 0.1, 200e9);
        let c2 = cmod_centre_crack(200e6, 0.05, 0.1, 200e9);
        assert!(
            (c2 / c1 - 2.0).abs() < 1e-10,
            "CMOD should scale linearly with stress"
        );
    }
    #[test]
    fn test_ppr_delta_n_peak() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.1);
        let expected = 0.1 * 1e-4;
        assert!(
            (p.delta_n_peak() - expected).abs() < 1e-16,
            "δ_n_peak = {}",
            p.delta_n_peak()
        );
    }
    #[test]
    fn test_ppr_delta_t_peak() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 5e-5, 0.2, 0.15);
        let expected = 0.15 * 5e-5;
        assert!(
            (p.delta_t_peak() - expected).abs() < 1e-18,
            "δ_t_peak = {}",
            p.delta_t_peak()
        );
    }
    #[test]
    fn test_ppr_sigma_max_positive() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.1);
        let s = p.sigma_max();
        assert!(s > 0.0 && s.is_finite(), "σ_max = {s}");
    }
    #[test]
    fn test_ppr_normal_traction_zero_at_ends() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.1);
        assert!(p.normal_traction(0.0).abs() < 1e-30, "T_n(0) should be 0");
        assert!(
            p.normal_traction(p.delta_n).abs() < 1e-30,
            "T_n(δ_n) should be 0"
        );
    }
    #[test]
    fn test_ppr_normal_traction_peak_at_lambda() {
        let p = PprParams::new(500.0, 250.0, 5.0, 5.0, 1e-3, 1e-3, 0.2, 0.2);
        let peak = p.delta_n_peak();
        let t_peak = p.normal_traction(peak);
        assert!(t_peak > 0.0, "T at peak = {t_peak}");
    }
    #[test]
    fn test_ppr_mode_i_energy_positive() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.1);
        let energy = p.mode_i_energy_numerical();
        assert!(energy > 0.0 && energy.is_finite(), "G_I_num = {energy}");
    }
    #[test]
    fn test_ppr_normal_traction_always_nonneg() {
        let p = PprParams::new(500.0, 250.0, 4.0, 4.0, 1e-4, 1e-4, 0.1, 0.1);
        for i in 0..20 {
            let d = p.delta_n * (i as f64) / 20.0;
            assert!(p.normal_traction(d) >= 0.0, "T_n({d:.2e}) < 0");
        }
    }
    #[test]
    fn test_viscous_traction_zero_rate() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let vp = ViscousCzmParams::new(base, 1e4);
        let d0 = base.delta_onset();
        let t_static = base.traction(d0 * 0.5);
        let t_visc = vp.traction(d0 * 0.5, 0.0);
        assert!(
            (t_static - t_visc).abs() < 1e-10,
            "zero-rate traction mismatch"
        );
    }
    #[test]
    fn test_viscous_traction_nonzero_rate() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let vp = ViscousCzmParams::new(base, 1e4);
        let t1 = vp.traction(1e-5, 0.0);
        let t2 = vp.traction(1e-5, 1.0);
        assert!(
            t2 > t1,
            "viscous traction should be higher with positive rate"
        );
    }
    #[test]
    fn test_viscous_traction_viscous_term() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let vp = ViscousCzmParams::new(base, 500.0);
        let t_v = vp.viscous_traction(0.01);
        assert!(
            (t_v - 500.0 * 0.01).abs() < 1e-12,
            "viscous traction = {t_v}"
        );
    }
    #[test]
    fn test_viscous_dissipation_positive() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let vp = ViscousCzmParams::new(base, 500.0);
        let diss = vp.viscous_dissipation(0.001, 1e-3);
        assert!(diss > 0.0, "dissipation = {diss}");
    }
    #[test]
    fn test_viscous_consistent_tangent_dt_effect() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let vp = ViscousCzmParams::new(base, 1e3);
        let kt1 = vp.consistent_tangent(1e-5, 1.0);
        let kt2 = vp.consistent_tangent(1e-5, 0.1);
        assert!(
            kt2 > kt1,
            "smaller dt should give larger tangent: {kt1} vs {kt2}"
        );
    }
    #[test]
    fn test_thermal_factor_at_ref_temp() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 20.0, 1.0, 1000.0);
        assert!(
            (tp.thermal_factor(20.0) - 1.0).abs() < 1e-12,
            "f(T_ref) = 1"
        );
    }
    #[test]
    fn test_thermal_factor_below_ref() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 20.0, 1.0, 1000.0);
        assert!(
            (tp.thermal_factor(0.0) - 1.0).abs() < 1e-12,
            "f(T < T_ref) = 1"
        );
    }
    #[test]
    fn test_thermal_factor_at_melt() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 0.0, 1.0, 100.0);
        assert!(tp.thermal_factor(100.0).abs() < 1e-12, "f(T_melt) = 0");
    }
    #[test]
    fn test_thermal_traction_reduced_at_high_temp() {
        let base = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 20.0, 1.0, 200.0);
        let d0 = base.delta_onset();
        let t_cold = tp.traction(d0 * 0.5, 20.0);
        let t_hot = tp.traction(d0 * 0.5, 120.0);
        assert!(
            t_hot < t_cold,
            "traction should reduce at higher temperature"
        );
    }
    #[test]
    fn test_thermal_fracture_toughness_reduces() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 0.0, 1.0, 200.0);
        let gc_cold = tp.fracture_toughness(0.0);
        let gc_hot = tp.fracture_toughness(100.0);
        assert!(gc_hot < gc_cold, "G_c should decrease at high temperature");
    }
    #[test]
    fn test_thermal_heat_flux() {
        let base = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let tp = ThermalCzmParams::new(base, 20.0, 1.0, 1000.0);
        let q = tp.heat_flux(1000.0, 50.0);
        assert!((q - 50000.0).abs() < 1e-8, "heat flux = {q}");
    }
    #[test]
    fn test_fatigue_state_initial() {
        let st = FatigueCzmState::new();
        assert!((st.total_damage()).abs() < 1e-12, "initial damage = 0");
        assert!(!st.is_failed(), "initial: not failed");
        assert_eq!(st.n_cycles, 0);
    }
    #[test]
    fn test_fatigue_record_separation() {
        let mut st = FatigueCzmState::new();
        st.record_separation(1e-5);
        st.record_separation(5e-5);
        st.record_separation(3e-5);
        assert!(
            (st.delta_max_cycle - 5e-5).abs() < 1e-20,
            "max = {}",
            st.delta_max_cycle
        );
        assert!(
            (st.delta_min_cycle - 1e-5).abs() < 1e-20,
            "min = {}",
            st.delta_min_cycle
        );
    }
    #[test]
    fn test_fatigue_advance_cycle_positive() {
        let mut st = FatigueCzmState::new();
        st.record_separation(1e-6);
        st.record_separation(1e-5);
        let delta_f = 1e-4;
        st.advance_cycle(1e3, 2.0, delta_f);
        assert!(st.d_fatigue > 0.0, "fatigue damage should accumulate");
        assert_eq!(st.n_cycles, 1);
    }
    #[test]
    fn test_fatigue_monotonic_update() {
        let mut st = FatigueCzmState::new();
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        st.update_monotonic(p.delta_failure * 0.5, &p);
        assert!(st.d_mono > 0.0, "monotonic damage should be > 0");
    }
    #[test]
    fn test_fatigue_total_damage_capped_at_one() {
        let mut st = FatigueCzmState::new();
        st.d_mono = 0.8;
        st.d_fatigue = 0.5;
        assert!(
            (st.total_damage() - 1.0).abs() < 1e-12,
            "damage capped at 1"
        );
    }
    #[test]
    fn test_fatigue_cycles_reset_after_advance() {
        let mut st = FatigueCzmState::new();
        st.record_separation(5e-5);
        st.advance_cycle(1e3, 2.0, 1e-4);
        assert!(
            st.delta_max_cycle < 1e-20,
            "delta_max_cycle should reset: {}",
            st.delta_max_cycle
        );
    }
    #[test]
    fn test_fatigue_params_cycles_to_failure_positive() {
        let base = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let fp = FatigueCzmParams::new(base, 1e3, 2.0, 1e-7);
        let nf = fp.cycles_to_failure(1e-5, 0.0);
        assert!(nf.is_some(), "should compute cycles to failure");
        assert!(nf.unwrap() > 0.0, "N_f must be positive");
    }
    #[test]
    fn test_fatigue_params_cycles_to_failure_below_threshold() {
        let base = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let fp = FatigueCzmParams::new(base, 1e3, 2.0, 1e-4);
        let nf = fp.cycles_to_failure(1e-5, 0.0);
        assert!(nf.is_none(), "below threshold → None");
    }
    #[test]
    fn test_fatigue_params_fewer_cycles_at_higher_delta() {
        let base = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let fp = FatigueCzmParams::new(base, 1e3, 2.0, 1e-7);
        let nf_low = fp.cycles_to_failure(1e-5, 0.0).unwrap();
        let nf_high = fp.cycles_to_failure(5e-5, 0.0).unwrap();
        assert!(nf_high < nf_low, "larger Δδ → fewer cycles to failure");
    }
    #[test]
    fn test_frictional_normal_traction_open() {
        let normal = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let shear = BilinearCzmParams::new(5e8, 30e6, 2e-4);
        let fp = FrictionalCzmParams::new(normal, shear, 0.3);
        let d0 = normal.delta_onset();
        let t = fp.normal_traction(d0 * 0.5, 1e9);
        assert!(t > 0.0, "normal opening traction = {t}");
    }
    #[test]
    fn test_frictional_normal_traction_contact() {
        let normal = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let shear = BilinearCzmParams::new(5e8, 30e6, 2e-4);
        let fp = FrictionalCzmParams::new(normal, shear, 0.3);
        let t = fp.normal_traction(-1e-4, 1e10);
        assert!(t > 0.0, "contact compressive traction = {t}");
    }
    #[test]
    fn test_frictional_shear_increases_in_contact() {
        let normal = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let shear = BilinearCzmParams::new(3e12, 30e6, 2e-4);
        let fp = FrictionalCzmParams::new(normal, shear, 0.4);
        let t_open = fp.shear_traction(1e-5, 1e-5, 1e9);
        let t_contact = fp.shear_traction(1e-5, -1e-5, 1e9);
        assert!(
            t_contact > t_open,
            "friction increases shear: {t_open} vs {t_contact}"
        );
    }
    #[test]
    fn test_frictional_max_shear_capacity() {
        let normal = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let shear = BilinearCzmParams::new(5e8, 30e6, 2e-4);
        let fp = FrictionalCzmParams::new(normal, shear, 0.3);
        let cap = fp.max_shear_capacity(100e6);
        let expected = 30e6 + 0.3 * 100e6;
        assert!((cap - expected).abs() < 1e-6, "capacity = {cap}");
    }
    #[test]
    fn test_dugdale_zone_length_positive() {
        let ld = dugdale_zone_length(100e6, 0.01, 200e6);
        assert!(ld > 0.0 && ld.is_finite(), "Dugdale L = {ld}");
    }
    #[test]
    fn test_dugdale_zone_scales_with_stress_squared() {
        let ld1 = dugdale_zone_length(100e6, 0.01, 200e6);
        let ld2 = dugdale_zone_length(200e6, 0.01, 200e6);
        assert!((ld2 / ld1 - 4.0).abs() < 1e-10, "L_D scale: {}", ld2 / ld1);
    }
    #[test]
    fn test_barenblatt_zone_length_positive() {
        let lb = barenblatt_zone_length(200e9, 500.0, 50e6);
        assert!(lb > 0.0 && lb.is_finite(), "Barenblatt L = {lb}");
    }
    #[test]
    fn test_barenblatt_zero_sigma_c_gives_infinity() {
        let lb = barenblatt_zone_length(200e9, 500.0, 0.0);
        assert!(lb.is_infinite(), "σ_c = 0 → L_B = ∞");
    }
    #[test]
    fn test_interface_displacement_jump_pure_normal() {
        let u_up = [0.0, 0.0, 1e-4];
        let u_lo = [0.0, 0.0, 0.0];
        let n = [0.0, 0.0, 1.0];
        let (dn, dt1, dt2) = interface_displacement_jump(u_up, u_lo, n);
        assert!((dn - 1e-4).abs() < 1e-15, "δ_n = {dn}");
        assert!(dt1.abs() < 1e-15, "δ_t1 = {dt1}");
        assert!(dt2.abs() < 1e-15, "δ_t2 = {dt2}");
    }
    #[test]
    fn test_interface_displacement_jump_pure_tangential() {
        let u_up = [1e-4, 0.0, 0.0];
        let u_lo = [0.0, 0.0, 0.0];
        let n = [0.0, 0.0, 1.0];
        let (dn, dt1, dt2) = interface_displacement_jump(u_up, u_lo, n);
        assert!(dn.abs() < 1e-15, "δ_n = {dn}");
        let dt_mag = (dt1 * dt1 + dt2 * dt2).sqrt();
        assert!((dt_mag - 1e-4).abs() < 1e-15, "‖δ_t‖ = {dt_mag}");
    }
    #[test]
    fn test_interface_displacement_jump_zero() {
        let n = [0.0, 1.0, 0.0];
        let (dn, dt1, dt2) = interface_displacement_jump([0.0; 3], [0.0; 3], n);
        assert!(dn.abs() < 1e-30, "δ_n = {dn}");
        assert!(
            dt1.abs() < 1e-30 && dt2.abs() < 1e-30,
            "δ_t = ({dt1}, {dt2})"
        );
    }
    #[test]
    fn test_effective_separation_mode_i_only() {
        let deff = effective_separation(1e-4, 0.0, 1.0);
        assert!((deff - 1e-4).abs() < 1e-20, "mode-I: δ_eff = {deff}");
    }
    #[test]
    fn test_effective_separation_mode_ii_only() {
        let deff = effective_separation(0.0, 1e-4, 1.0);
        assert!((deff - 1e-4).abs() < 1e-20, "mode-II: δ_eff = {deff}");
    }
    #[test]
    fn test_effective_separation_compressive_ignored() {
        let deff = effective_separation(-1e-4, 0.0, 1.0);
        assert!(deff.abs() < 1e-20, "compressive: δ_eff = {deff}");
    }
    #[test]
    fn test_mode_mixity_angle_pure_mode_i() {
        let psi = mode_mixity_angle_degrees(1.0, 0.0);
        assert!(psi.abs() < 1e-12, "pure mode I ψ = {psi}°");
    }
    #[test]
    fn test_mode_mixity_angle_pure_mode_ii() {
        let psi = mode_mixity_angle_degrees(0.0, 1.0);
        assert!((psi - 90.0).abs() < 1e-10, "pure mode II ψ = {psi}°");
    }
    #[test]
    fn test_polynomial_degradation_no_damage() {
        let e = polynomial_stiffness_degradation(200e9, 0.0, 2.0);
        assert!((e - 200e9).abs() < 1e-3, "E(D=0) = {e}");
    }
    #[test]
    fn test_polynomial_degradation_full_damage() {
        let e = polynomial_stiffness_degradation(200e9, 1.0, 2.0);
        assert!(e.abs() < 1e-6, "E(D=1) = {e}");
    }
    #[test]
    fn test_polynomial_degradation_monotone() {
        let e1 = polynomial_stiffness_degradation(200e9, 0.3, 2.0);
        let e2 = polynomial_stiffness_degradation(200e9, 0.7, 2.0);
        assert!(e1 > e2, "E should decrease with damage");
    }
    #[test]
    fn test_bilinear_tangent_modulus_softening_negative() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let kt = bilinear_tangent_modulus_softening(&p);
        assert!(
            kt < 0.0,
            "tangent modulus in softening should be negative: {kt}"
        );
    }
    #[test]
    fn test_bilinear_secant_unloading_in_linear_branch() {
        let p = BilinearCzmParams::new(1e9, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let ks = bilinear_secant_unloading(&p, d0 * 0.5);
        assert!(
            (ks - p.penalty_stiffness).abs() / p.penalty_stiffness < 1e-10,
            "secant in linear branch = {ks}"
        );
    }
    #[test]
    fn test_bilinear_secant_unloading_softening_less_than_penalty() {
        let p = BilinearCzmParams::new(5e12, 50e6, 1e-4);
        let d0 = p.delta_onset();
        let df = p.delta_failure;
        let ks = bilinear_secant_unloading(&p, 0.5 * (d0 + df));
        assert!(ks < p.penalty_stiffness, "secant in softening < K_0: {ks}");
    }
    #[test]
    fn test_weibull_survival_prob_no_traction() {
        let ps = weibull_survival_probability(0.0, 1.0, 50e6, 1e-4, 5.0);
        assert!((ps - 1.0).abs() < 1e-12, "P_s(T=0) = 1");
    }
    #[test]
    fn test_weibull_survival_prob_at_ref() {
        let ps = weibull_survival_probability(50e6, 1e-4, 50e6, 1e-4, 1.0);
        assert!(
            (ps - (-1.0_f64).exp()).abs() < 1e-12,
            "P_s(T_ref,A_ref) = {ps}"
        );
    }
    #[test]
    fn test_weibull_survival_decreases_with_traction() {
        let ps1 = weibull_survival_probability(10e6, 1e-4, 50e6, 1e-4, 5.0);
        let ps2 = weibull_survival_probability(40e6, 1e-4, 50e6, 1e-4, 5.0);
        assert!(ps2 < ps1, "higher traction → lower survival probability");
    }
    #[test]
    fn test_weibull_characteristic_strength_positive() {
        let ts = weibull_characteristic_strength(0.368, 1e-4, 50e6, 1e-4, 5.0);
        assert!(ts > 0.0 && ts.is_finite(), "characteristic strength = {ts}");
    }
    #[test]
    fn test_weibull_characteristic_strength_at_ref() {
        let ts = weibull_characteristic_strength((-1.0_f64).exp(), 1e-4, 50e6, 1e-4, 1.0);
        assert!((ts - 50e6).abs() / 50e6 < 1e-10, "T_char at ref = {ts}");
    }
    #[test]
    fn test_weibull_characteristic_strength_smaller_area_higher_strength() {
        let ts_large = weibull_characteristic_strength(0.5, 1e-3, 50e6, 1e-4, 5.0);
        let ts_small = weibull_characteristic_strength(0.5, 1e-5, 50e6, 1e-4, 5.0);
        assert!(
            ts_small > ts_large,
            "smaller area → higher characteristic strength"
        );
    }
}
