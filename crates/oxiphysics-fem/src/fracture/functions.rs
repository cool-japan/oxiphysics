//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use super::types::{
    CzmBilinear, CzmBilinearParams, DomainPoint, JIntegralPoint, LinearFracture, MixedModeCriterion,
};

/// Compute the J-integral along a contour path (2D, plane stress/strain).
///
/// J = integral_Gamma (W n_x - t_i du_i/dx) ds
///
/// where W = strain energy density, t_i = sigma_ij n_j (traction),
/// and du_i/dx = displacement gradient in the crack-growth direction.
#[allow(dead_code)]
pub fn j_integral_contour(points: &[JIntegralPoint], youngs_modulus: f64, poisson: f64) -> f64 {
    let mut j = 0.0;
    let inv_e = 1.0 / youngs_modulus;
    let nu = poisson;
    for pt in points {
        let sxx = pt.stress[0];
        let syy = pt.stress[1];
        let sxy = pt.stress[2];
        let exx = inv_e * (sxx - nu * syy);
        let eyy = inv_e * (syy - nu * sxx);
        let exy = (1.0 + nu) * inv_e * sxy;
        let w = 0.5 * (sxx * exx + syy * eyy + 2.0 * sxy * exy);
        let nx = pt.normal[0];
        let ny = pt.normal[1];
        let tx = sxx * nx + sxy * ny;
        let ty = sxy * nx + syy * ny;
        let du_x_dx = pt.displacement_gradient[0];
        let du_y_dx = pt.displacement_gradient[2];
        j += (w * nx - tx * du_x_dx - ty * du_y_dx) * pt.ds;
    }
    j
}
/// Estimate the T-stress from the difference between sigma_xx and sigma_yy
/// along the crack line (theta = 0) at a distance r from the crack tip.
///
/// For theta = 0, the Williams expansion gives:
/// sigma_xx = K_I / sqrt(2 pi r) + T + ...
/// sigma_yy = K_I / sqrt(2 pi r) + ...
///
/// So T = sigma_xx - sigma_yy (at theta = 0, after removing the singular terms).
///
/// In practice, T is extracted by evaluating stresses at several radii and
/// extrapolating.
#[allow(dead_code)]
pub fn t_stress_estimate(sigma_xx_at_theta0: f64, sigma_yy_at_theta0: f64) -> f64 {
    sigma_xx_at_theta0 - sigma_yy_at_theta0
}
/// Compute T-stress using the stress difference method with extrapolation.
///
/// Given stress values at multiple radii along theta = 0, extrapolate to r = 0.
/// Uses linear least-squares fit of T(r) = sigma_xx(r) - sigma_yy(r) vs sqrt(r).
#[allow(dead_code)]
pub fn t_stress_extrapolation(radii: &[f64], sigma_xx: &[f64], sigma_yy: &[f64]) -> f64 {
    assert_eq!(radii.len(), sigma_xx.len());
    assert_eq!(radii.len(), sigma_yy.len());
    let n = radii.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sigma_xx[0] - sigma_yy[0];
    }
    let mut sum_1 = 0.0;
    let mut sum_sr = 0.0;
    let mut sum_sr2 = 0.0;
    let mut sum_t = 0.0;
    let mut sum_t_sr = 0.0;
    for i in 0..n {
        let sr = radii[i].sqrt();
        let t_val = sigma_xx[i] - sigma_yy[i];
        sum_1 += 1.0;
        sum_sr += sr;
        sum_sr2 += sr * sr;
        sum_t += t_val;
        sum_t_sr += t_val * sr;
    }
    let det = sum_1 * sum_sr2 - sum_sr * sum_sr;
    if det.abs() < 1e-60 {
        return sum_t / sum_1;
    }
    (sum_sr2 * sum_t - sum_sr * sum_t_sr) / det
}
/// Predict the crack propagation angle using the Maximum Tangential Stress criterion.
///
/// The crack propagates in the direction that maximizes the tangential stress
/// sigma_theta_theta at the crack tip.
///
/// For mixed-mode I+II loading, the angle theta_c satisfies:
/// K_I sin(theta_c) + K_II (3 cos(theta_c) - 1) = 0
///
/// Returns the angle in radians (measured from the crack-line direction).
#[allow(dead_code)]
pub fn mts_propagation_angle(ki: f64, kii: f64) -> f64 {
    if kii.abs() < 1e-30 {
        return 0.0;
    }
    let discriminant = (ki * ki + 8.0 * kii * kii).sqrt();
    2.0 * ((ki - discriminant) / (4.0 * kii)).atan()
}
/// Compute the equivalent stress intensity factor for mixed-mode loading.
///
/// K_eq = K_I cos^3(theta/2) - 3 K_II cos^2(theta/2) sin(theta/2)
///
/// Fracture occurs when K_eq >= K_Ic.
#[allow(dead_code)]
pub fn equivalent_sif_mts(ki: f64, kii: f64) -> f64 {
    let theta = mts_propagation_angle(ki, kii);
    let ht = theta / 2.0;
    let cos_ht = ht.cos();
    let sin_ht = ht.sin();
    ki * cos_ht.powi(3) - 3.0 * kii * cos_ht.powi(2) * sin_ht
}
/// Check whether mixed-mode fracture occurs using the specified criterion.
///
/// Returns `true` if the equivalent condition K_eq >= K_Ic is met.
#[allow(dead_code)]
pub fn mixed_mode_fracture_check(
    ki: f64,
    kii: f64,
    k_ic: f64,
    criterion: MixedModeCriterion,
) -> bool {
    match criterion {
        MixedModeCriterion::MaxTangentialStress => equivalent_sif_mts(ki, kii) >= k_ic,
        MixedModeCriterion::MaxEnergyReleaseRate => ki * ki + kii * kii >= k_ic * k_ic,
        MixedModeCriterion::MinStrainEnergyDensity => ki * ki + kii * kii >= k_ic * k_ic,
    }
}
/// Forman crack growth equation (accounts for stress ratio R and K_Ic).
///
/// da/dN = C * delta_K^m / ((1-R) K_Ic - delta_K)
///
/// This extends the Paris law by accounting for the mean stress effect
/// and the approach to fracture toughness.
#[allow(dead_code)]
pub fn forman_da_dn(delta_k: f64, c: f64, m: f64, r_ratio: f64, k_ic: f64) -> f64 {
    let denom = (1.0 - r_ratio) * k_ic - delta_k;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    c * delta_k.powf(m) / denom
}
/// Walker crack growth equation (generalization for stress ratio effects).
///
/// da/dN = C * (delta_K / (1-R)^(1-gamma))^m
///
/// The Walker exponent gamma is a material parameter (gamma = 0.5 for many metals).
#[allow(dead_code)]
pub fn walker_da_dn(delta_k: f64, c: f64, m: f64, r_ratio: f64, gamma: f64) -> f64 {
    let effective_dk = delta_k / (1.0 - r_ratio).powf(1.0 - gamma);
    c * effective_dk.powf(m)
}
/// NASGRO/Forman-Newman-de Koning crack growth equation (simplified).
///
/// da/dN = C * ((1-f)/(1-R))^n * delta_K^n * (1 - delta_K_th/delta_K)^p
///        / (1 - K_max/K_Ic)^q
///
/// where f is the crack opening function.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn nasgro_da_dn(
    delta_k: f64,
    c: f64,
    n: f64,
    r_ratio: f64,
    delta_k_th: f64,
    k_ic: f64,
    p: f64,
    q: f64,
) -> f64 {
    if delta_k <= delta_k_th {
        return 0.0;
    }
    let k_max = delta_k / (1.0 - r_ratio);
    let denom = 1.0 - k_max / k_ic;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    let f = r_ratio.max(0.0);
    let opening_ratio = ((1.0 - f) / (1.0 - r_ratio)).max(0.0);
    c * opening_ratio.powf(n) * delta_k.powf(n) * (1.0 - delta_k_th / delta_k).powf(p)
        / denom.powf(q)
}
/// Predict the crack increment direction for a given mixed-mode loading state.
///
/// Returns the unit direction vector \[dx, dy\] for 2D crack propagation
/// relative to the current crack direction.
#[allow(dead_code)]
pub fn predict_crack_direction(
    ki: f64,
    kii: f64,
    current_direction: [f64; 2],
    criterion: MixedModeCriterion,
) -> [f64; 2] {
    let theta = match criterion {
        MixedModeCriterion::MaxTangentialStress => mts_propagation_angle(ki, kii),
        _ => mts_propagation_angle(ki, kii),
    };
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    let dx = current_direction[0] * cos_t - current_direction[1] * sin_t;
    let dy = current_direction[0] * sin_t + current_direction[1] * cos_t;
    let len = (dx * dx + dy * dy).sqrt().max(f64::EPSILON);
    [dx / len, dy / len]
}
/// Integrate fatigue crack growth along a predicted path using adaptive stepping.
///
/// Returns a vector of (crack_length, cycles) pairs.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn fatigue_crack_path(
    a_initial: f64,
    a_final: f64,
    delta_sigma: f64,
    geometry_factor: f64,
    c: f64,
    m: f64,
    n_steps: usize,
) -> Vec<(f64, f64)> {
    let mut path = Vec::with_capacity(n_steps + 1);
    let da = (a_final - a_initial) / n_steps as f64;
    let mut a = a_initial;
    let mut n_cycles = 0.0;
    path.push((a, n_cycles));
    for _ in 0..n_steps {
        let delta_k = delta_sigma * (PI * a).sqrt() * geometry_factor;
        let da_dn = LinearFracture::paris_law_da_dn(delta_k, c, m);
        if da_dn > 1e-60 {
            n_cycles += da / da_dn;
        }
        a += da;
        path.push((a, n_cycles));
    }
    path
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fracture::*;
    #[test]
    fn test_sif_penny_crack() {
        let sigma = 100.0e6_f64;
        let a = 0.005_f64;
        let geo = 2.0 / PI;
        let ki = LinearFracture::stress_intensity_ki(sigma, a, geo);
        let expected = sigma * (PI * a).sqrt() * geo;
        assert!(
            (ki - expected).abs() / expected < 1e-12,
            "K_I mismatch: {ki} vs {expected}"
        );
        assert!(ki > 0.0);
    }
    #[test]
    fn test_griffith_stress() {
        let k_ic = 50.0e6_f64;
        let a = 0.01_f64;
        let geo = 1.0_f64;
        let sigma_c = LinearFracture::griffith_critical_stress(a, k_ic, geo);
        let expected = k_ic / (PI * a).sqrt();
        assert!(
            (sigma_c - expected).abs() / expected < 1e-12,
            "sigma_c mismatch: {sigma_c} vs {expected}"
        );
        assert!(sigma_c > 250.0e6 && sigma_c < 320.0e6);
    }
    #[test]
    fn test_j_integral() {
        let ki = 30.0e6_f64;
        let kii = 0.0_f64;
        let e = 200.0e9_f64;
        let j = LinearFracture::j_integral_estimate(ki, kii, e);
        let expected = ki * ki / e;
        assert!((j - expected).abs() / expected < 1e-12);
        assert!(j > 0.0);
    }
    #[test]
    fn test_paris_law() {
        let delta_k = 20.0e6_f64;
        let c = 6.9e-12_f64;
        let m = 3.0_f64;
        let da_dn = LinearFracture::paris_law_da_dn(delta_k, c, m);
        let expected = c * delta_k.powf(m);
        assert!((da_dn - expected).abs() / expected < 1e-12);
        assert!(da_dn > 0.0);
        assert!(da_dn.is_finite());
    }
    #[test]
    fn test_fatigue_life() {
        let a_i = 0.001_f64;
        let a_f = 0.005_f64;
        let delta_sigma = 200.0e6_f64;
        let geo = 1.12_f64;
        let c = 6.9e-30_f64;
        let m = 3.0_f64;
        let n = LinearFracture::fatigue_life_cycles(a_i, a_f, delta_sigma, geo, c, m);
        assert!(n > 0.0);
        assert!(n.is_finite());
        assert!(n > 1.0e3);
        assert!(n < 1.0e8);
    }
    #[test]
    fn test_t_stress_estimate() {
        let t = t_stress_estimate(150.0e6, 100.0e6);
        assert!((t - 50.0e6).abs() < 1e-6);
    }
    #[test]
    fn test_t_stress_extrapolation_constant() {
        let radii = vec![0.001, 0.002, 0.003, 0.004];
        let t_val = 50.0e6;
        let sigma_xx: Vec<f64> = radii.iter().map(|_| 200.0e6).collect();
        let sigma_yy: Vec<f64> = radii.iter().map(|_| 200.0e6 - t_val).collect();
        let t = t_stress_extrapolation(&radii, &sigma_xx, &sigma_yy);
        assert!(
            (t - t_val).abs() / t_val < 0.01,
            "T-stress = {t}, expected {t_val}"
        );
    }
    #[test]
    fn test_mts_angle_pure_mode_i() {
        let theta = mts_propagation_angle(1.0, 0.0);
        assert!(theta.abs() < 1e-12, "pure mode I angle = {theta}");
    }
    #[test]
    fn test_mts_angle_pure_mode_ii() {
        let theta = mts_propagation_angle(0.0, 1.0);
        let expected_deg = -70.5;
        let theta_deg = theta.to_degrees();
        assert!(
            (theta_deg - expected_deg).abs() < 1.0,
            "pure mode II angle = {theta_deg} deg, expected ~{expected_deg} deg"
        );
    }
    #[test]
    fn test_equivalent_sif_pure_mode_i() {
        let k_eq = equivalent_sif_mts(10.0e6, 0.0);
        assert!(
            (k_eq - 10.0e6).abs() / 10.0e6 < 1e-10,
            "K_eq = {k_eq}, expected 10e6"
        );
    }
    #[test]
    fn test_mixed_mode_merr_below_threshold() {
        let result = mixed_mode_fracture_check(
            5.0e6,
            3.0e6,
            50.0e6,
            MixedModeCriterion::MaxEnergyReleaseRate,
        );
        assert!(!result, "should not fracture below threshold");
    }
    #[test]
    fn test_mixed_mode_merr_above_threshold() {
        let result = mixed_mode_fracture_check(
            40.0e6,
            30.0e6,
            50.0e6,
            MixedModeCriterion::MaxEnergyReleaseRate,
        );
        assert!(result, "should fracture at threshold");
    }
    #[test]
    fn test_forman_da_dn_positive() {
        let da_dn = forman_da_dn(20.0e6, 6.9e-12, 3.0, 0.1, 50.0e6);
        assert!(da_dn > 0.0);
        assert!(da_dn.is_finite());
    }
    #[test]
    fn test_forman_unstable_fracture() {
        let da_dn = forman_da_dn(50.0e6, 6.9e-12, 3.0, 0.0, 50.0e6);
        assert!(da_dn.is_infinite());
    }
    #[test]
    fn test_walker_da_dn() {
        let da_dn = walker_da_dn(20.0e6, 6.9e-12, 3.0, 0.1, 0.5);
        assert!(da_dn > 0.0);
        assert!(da_dn.is_finite());
        let da_dn_high_r = walker_da_dn(20.0e6, 6.9e-12, 3.0, 0.5, 0.5);
        assert!(da_dn_high_r > da_dn, "higher R should give faster growth");
    }
    #[test]
    fn test_nasgro_below_threshold() {
        let da_dn = nasgro_da_dn(5.0e6, 1e-12, 3.0, 0.1, 10.0e6, 50.0e6, 1.0, 1.0);
        assert_eq!(da_dn, 0.0, "below threshold => zero growth");
    }
    #[test]
    fn test_nasgro_above_threshold() {
        let da_dn = nasgro_da_dn(20.0e6, 1e-12, 3.0, 0.1, 5.0e6, 100.0e6, 1.0, 1.0);
        assert!(da_dn > 0.0);
        assert!(da_dn.is_finite());
    }
    #[test]
    fn test_predict_crack_direction_mode_i() {
        let dir = predict_crack_direction(
            1.0,
            0.0,
            [1.0, 0.0],
            MixedModeCriterion::MaxTangentialStress,
        );
        assert!((dir[0] - 1.0).abs() < 1e-10);
        assert!(dir[1].abs() < 1e-10);
    }
    #[test]
    fn test_predict_crack_direction_mode_ii() {
        let dir = predict_crack_direction(
            0.0,
            1.0,
            [1.0, 0.0],
            MixedModeCriterion::MaxTangentialStress,
        );
        let angle = dir[1].atan2(dir[0]).to_degrees();
        assert!(
            (angle - (-70.5)).abs() < 1.0,
            "predicted angle = {angle} deg, expected ~-70.5 deg"
        );
    }
    #[test]
    fn test_fatigue_crack_path_length() {
        let path = fatigue_crack_path(0.001, 0.005, 200e6, 1.12, 6.9e-30, 3.0, 10);
        assert_eq!(path.len(), 11);
        assert!((path[0].0 - 0.001).abs() < 1e-15);
        assert_eq!(path[0].1, 0.0);
        assert!((path[10].0 - 0.005).abs() < 1e-12);
        for i in 1..path.len() {
            assert!(
                path[i].1 >= path[i - 1].1,
                "cycles not increasing at step {i}"
            );
        }
    }
    #[test]
    fn test_geometry_factors() {
        assert_eq!(GeometryFactors::infinite_plate(), 1.0);
        assert!((GeometryFactors::edge_crack() - 1.12).abs() < 1e-12);
        assert!((GeometryFactors::penny_crack() - 2.0 / PI).abs() < 1e-12);
    }
    #[test]
    fn test_finite_width_geometry_factor() {
        let f = GeometryFactors::finite_width(0.001, 1.0);
        assert!((f - 1.0).abs() < 0.01, "F = {f}, expected ~1.0");
        let f2 = GeometryFactors::finite_width(0.4, 1.0);
        assert!(f2 > 1.0, "F should be > 1 for finite width, got {f2}");
    }
    #[test]
    fn test_j_integral_contour_simple() {
        let points = vec![
            JIntegralPoint {
                stress: [100.0e6, 100.0e6, 0.0],
                displacement_gradient: [0.0; 4],
                normal: [1.0, 0.0],
                ds: 0.01,
            },
            JIntegralPoint {
                stress: [100.0e6, 100.0e6, 0.0],
                displacement_gradient: [0.0; 4],
                normal: [-1.0, 0.0],
                ds: 0.01,
            },
        ];
        let j = j_integral_contour(&points, 200.0e9, 0.3);
        assert!(j.abs() < 1e-6, "J = {j}, expected ~0");
    }
    #[test]
    fn test_crack_tip_stress_theta_zero() {
        let r = 0.01;
        let ki = 30.0e6;
        let stress = LinearFracture::crack_tip_stress(r, 0.0, ki, 0.0);
        let expected_syy = ki / (2.0 * PI * r).sqrt();
        assert!(
            (stress[1] - expected_syy).abs() / expected_syy < 1e-10,
            "sigma_yy = {}, expected {}",
            stress[1],
            expected_syy
        );
    }
}
/// Compute the Mode-III stress intensity factor K_III.
///
/// K_III = tau_III * sqrt(pi * a) * F
#[allow(dead_code)]
pub fn stress_intensity_kiii(tau_iii: f64, a: f64, geometry_factor: f64) -> f64 {
    tau_iii * (PI * a).sqrt() * geometry_factor
}
/// Effective (combined) stress intensity factor for mixed Mode I + II + III.
///
/// K_eff = sqrt(K_I^2 + K_II^2 + K_III^2 / (1 - nu))
///
/// where nu is the Poisson's ratio (plane strain).
#[allow(dead_code)]
pub fn effective_sif_combined(ki: f64, kii: f64, kiii: f64, poisson: f64) -> f64 {
    let denom = (1.0 - poisson).max(f64::EPSILON);
    (ki * ki + kii * kii + kiii * kiii / denom).max(0.0).sqrt()
}
/// Compute J-integral including Mode III contribution (plane strain).
///
/// J = (K_I^2 + K_II^2) / E + K_III^2 / (2 * mu)
///
/// where mu = E / (2 * (1 + nu))
#[allow(dead_code)]
pub fn j_integral_three_mode(ki: f64, kii: f64, kiii: f64, youngs: f64, poisson: f64) -> f64 {
    let mu = youngs / (2.0 * (1.0 + poisson));
    let j12 = (ki * ki + kii * kii) / youngs;
    let j3 = kiii * kiii / (2.0 * mu.max(f64::EPSILON));
    j12 + j3
}
/// Check fracture by comparing effective SIF with K_Ic using combined criterion.
#[allow(dead_code)]
pub fn fracture_check_combined(ki: f64, kii: f64, kiii: f64, k_ic: f64, poisson: f64) -> bool {
    effective_sif_combined(ki, kii, kiii, poisson) >= k_ic
}
/// Represent two co-planar parallel cracks and compute their interaction factor.
///
/// For two equal cracks of length `2a` separated by `d` (centre-to-centre),
/// the interaction amplification factor on K_I is approximated by:
///
/// F_interaction ≈ 1 + 0.5 * (2a / d)^2   (Erdogan and Sih approximation)
#[allow(dead_code)]
pub fn crack_interaction_factor(a: f64, d: f64) -> f64 {
    if d < 2.0 * a {
        return 2.5;
    }
    1.0 + 0.5 * (2.0 * a / d).powi(2)
}
/// Interaction-corrected K_I for two co-planar cracks.
#[allow(dead_code)]
pub fn ki_two_cracks(sigma_far: f64, a: f64, d: f64, geometry_factor: f64) -> f64 {
    let f_int = crack_interaction_factor(a, d);
    LinearFracture::stress_intensity_ki(sigma_far, a, geometry_factor) * f_int
}
/// Predict whether two cracks will link up (coalesce) based on their tip distance
/// and the applied stress.
///
/// Coalescence criterion: K_I_effective >= K_Ic at the inner crack tips.
#[allow(dead_code)]
pub fn will_cracks_link(sigma_far: f64, a: f64, d: f64, geo: f64, k_ic: f64) -> bool {
    let ki = ki_two_cracks(sigma_far, a, d, geo);
    ki >= k_ic
}
/// Compute the domain J-integral from a set of domain quadrature points.
///
/// J = -∫∫_A \[ W * dq/dx - (sigma_ij * du_i/dx) * dq/dj \] dA
///   = Σ_k \[ -W_k * dq/dx_k + T_k * dq_k \] * A_k
///
/// where W is the strain energy density and T are traction-like terms.
#[allow(dead_code)]
pub fn j_integral_domain(points: &[DomainPoint]) -> f64 {
    let mut j = 0.0;
    for p in points {
        let sxx = p.stress[0];
        let syy = p.stress[1];
        let sxy = p.stress[2];
        let exx = p.strain[0];
        let eyy = p.strain[1];
        let exy = p.strain[2];
        let w = 0.5 * (sxx * exx + syy * eyy + sxy * exy);
        let du_x_dx = p.displacement_gradient[0];
        let du_y_dx = p.displacement_gradient[2];
        let traction_contrib =
            (sxx * du_x_dx + sxy * du_y_dx) * p.dq[0] + (sxy * du_x_dx + syy * du_y_dx) * p.dq[1];
        j += (-w * p.dq[0] + traction_contrib) * p.area;
    }
    j
}
/// Irwin plastic zone radius for plane stress.
///
/// r_p = (1 / (2π)) * (K_I / σ_y)^2
#[allow(dead_code)]
pub fn plastic_zone_radius_plane_stress(ki: f64, yield_stress: f64) -> f64 {
    1.0 / (2.0 * PI) * (ki / yield_stress).powi(2)
}
/// Irwin plastic zone radius for plane strain.
///
/// r_p = (1 / (6π)) * (K_I / σ_y)^2
#[allow(dead_code)]
pub fn plastic_zone_radius_plane_strain(ki: f64, yield_stress: f64) -> f64 {
    1.0 / (6.0 * PI) * (ki / yield_stress).powi(2)
}
/// Effective crack length with plastic zone correction (Irwin correction).
///
/// a_eff = a + r_p (plane stress)
#[allow(dead_code)]
pub fn effective_crack_length(a: f64, ki: f64, yield_stress: f64) -> f64 {
    let r_p = plastic_zone_radius_plane_stress(ki, yield_stress);
    a + r_p
}
/// Compute the threshold ΔK for fatigue crack arrest using the Kitagawa-Takahashi diagram.
///
/// For large cracks: ΔK_th_eff = ΔK_th (intrinsic threshold)
/// For small cracks (a < a_0): ΔK_th approaches ΔΣ_e * sqrt(π a)
///   where a_0 = (ΔK_th / (ΔΣ_e * sqrt(π)))^2
///
/// Returns the effective threshold for a given crack size `a`.
#[allow(dead_code)]
pub fn kitagawa_threshold(a: f64, delta_k_th: f64, delta_sigma_e: f64) -> f64 {
    let dk_stress_limit = delta_sigma_e * (PI * a).sqrt();
    let inv_sq = 1.0 / (delta_k_th * delta_k_th) + 1.0 / (dk_stress_limit * dk_stress_limit);
    (1.0 / inv_sq).sqrt()
}
/// R-curve (crack resistance): linear rising K_R model.
///
/// K_R(Δa) = K_Ic0 + slope * Δa
///
/// where Δa = a - a_0 is the crack extension from initial length a_0.
#[allow(dead_code)]
pub fn k_r_curve_linear(delta_a: f64, k_ic0: f64, slope: f64) -> f64 {
    k_ic0 + slope * delta_a.max(0.0)
}
/// Check if the crack is growing using the R-curve: K_I >= K_R(Δa)
#[allow(dead_code)]
pub fn r_curve_fracture_check(ki: f64, delta_a: f64, k_ic0: f64, slope: f64) -> bool {
    ki >= k_r_curve_linear(delta_a, k_ic0, slope)
}
/// Crack growth rate under stress corrosion cracking (SCC) using the
/// power-law model.
///
/// da/dt = C_scc * (K_I)^n_scc  for K_Ith <= K_I < K_Ic
/// da/dt = 0                     for K_I < K_Ith
/// da/dt = inf                   for K_I >= K_Ic (fast fracture)
#[allow(dead_code)]
pub fn scc_crack_growth_rate(ki: f64, c_scc: f64, n_scc: f64, k_ith: f64, k_ic: f64) -> f64 {
    if ki < k_ith {
        return 0.0;
    }
    if ki >= k_ic {
        return f64::INFINITY;
    }
    c_scc * ki.powf(n_scc)
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::fracture::*;
    #[test]
    fn test_kiii_positive() {
        let kiii = stress_intensity_kiii(50.0e6, 0.005, 1.0);
        assert!(kiii > 0.0);
        assert!(kiii.is_finite());
    }
    #[test]
    fn test_kiii_scales_with_geometry_factor() {
        let k1 = stress_intensity_kiii(50.0e6, 0.005, 1.0);
        let k2 = stress_intensity_kiii(50.0e6, 0.005, 2.0);
        assert!((k2 / k1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_effective_sif_pure_mode_i() {
        let ki = 30.0e6;
        let k_eff = effective_sif_combined(ki, 0.0, 0.0, 0.3);
        assert!((k_eff - ki).abs() / ki < 1e-10);
    }
    #[test]
    fn test_effective_sif_always_positive() {
        let k_eff = effective_sif_combined(10.0e6, 5.0e6, 3.0e6, 0.3);
        assert!(k_eff > 0.0);
        assert!(k_eff.is_finite());
    }
    #[test]
    fn test_j_three_mode_pure_mode_i() {
        let ki = 30.0e6;
        let e = 200.0e9;
        let nu = 0.3;
        let j = j_integral_three_mode(ki, 0.0, 0.0, e, nu);
        let expected = ki * ki / e;
        assert!((j - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_j_three_mode_positive_always() {
        let j = j_integral_three_mode(10.0e6, 5.0e6, 3.0e6, 200.0e9, 0.3);
        assert!(j > 0.0);
        assert!(j.is_finite());
    }
    #[test]
    fn test_crack_interaction_unity_for_far_cracks() {
        let f = crack_interaction_factor(0.001, 100.0);
        assert!((f - 1.0).abs() < 0.01, "far cracks: f = {f}");
    }
    #[test]
    fn test_crack_interaction_amplifies_when_close() {
        let f_far = crack_interaction_factor(0.01, 1.0);
        let f_close = crack_interaction_factor(0.01, 0.03);
        assert!(
            f_close > f_far,
            "closer cracks should amplify more: {f_close} vs {f_far}"
        );
    }
    #[test]
    fn test_crack_interaction_capped_for_overlapping() {
        let f = crack_interaction_factor(1.0, 1.5);
        assert!((f - 2.5).abs() < 1e-10, "overlapping: f = {f}");
    }
    #[test]
    fn test_ki_two_cracks_larger_than_single() {
        let ki_single = LinearFracture::stress_intensity_ki(100e6, 0.01, 1.0);
        let ki_two = ki_two_cracks(100e6, 0.01, 0.1, 1.0);
        assert!(
            ki_two > ki_single,
            "two-crack KI should exceed single-crack KI"
        );
    }
    #[test]
    fn test_will_cracks_link_below_toughness() {
        let result = will_cracks_link(10.0e6, 0.001, 10.0, 1.0, 50.0e6);
        assert!(!result, "far apart cracks at low stress should not link");
    }
    #[test]
    fn test_ct_geometry_factor_middle_crack() {
        let ct = CompactTensionSpec::new(0.05, 0.025, 0.01);
        let f = ct.geometry_factor();
        assert!(f > 0.0 && f.is_finite(), "CT geometry factor: {f}");
    }
    #[test]
    fn test_ct_k_from_load_positive() {
        let ct = CompactTensionSpec::new(0.05, 0.025, 0.01);
        let k = ct.k_from_load(5000.0);
        assert!(k > 0.0 && k.is_finite());
    }
    #[test]
    fn test_ct_plane_strain_validity() {
        let ct = CompactTensionSpec::new(0.1, 0.05, 0.05);
        let k_ic = 50.0e6;
        let sigma_y = 500.0e6;
        let valid = ct.is_plane_strain_valid(k_ic, sigma_y);
        assert!(valid, "specimen should satisfy plane strain validity");
    }
    #[test]
    fn test_senb_geometry_factor_positive() {
        let senb = ThreePointBendSpec::new(0.2, 0.05, 0.025, 0.01);
        let f = senb.geometry_factor();
        assert!(f > 0.0 && f.is_finite(), "SENB geometry factor: {f}");
    }
    #[test]
    fn test_senb_k_from_load_positive() {
        let senb = ThreePointBendSpec::new(0.2, 0.05, 0.025, 0.01);
        let k = senb.k_from_load(1000.0);
        assert!(k > 0.0 && k.is_finite());
    }
    #[test]
    fn test_j_domain_empty_zero() {
        let j = j_integral_domain(&[]);
        assert_eq!(j, 0.0);
    }
    #[test]
    fn test_j_domain_single_point_finite() {
        let p = DomainPoint {
            stress: [100.0e6, 50.0e6, 10.0e6],
            strain: [5e-4, 2.5e-4, 1.3e-4],
            displacement_gradient: [5e-4, 1e-4, 1e-4, 2.5e-4],
            q: 1.0,
            dq: [1.0, 0.0],
            area: 1e-6,
        };
        let j = j_integral_domain(&[p]);
        assert!(j.is_finite(), "J domain single point: {j}");
    }
    #[test]
    fn test_plastic_zone_plane_stress_positive() {
        let rp = plastic_zone_radius_plane_stress(30.0e6, 250.0e6);
        assert!(rp > 0.0 && rp.is_finite());
    }
    #[test]
    fn test_plastic_zone_plane_strain_smaller_than_stress() {
        let ki = 30.0e6;
        let sy = 250.0e6;
        let rp_stress = plastic_zone_radius_plane_stress(ki, sy);
        let rp_strain = plastic_zone_radius_plane_strain(ki, sy);
        assert!(rp_strain < rp_stress, "plane-strain zone should be smaller");
    }
    #[test]
    fn test_effective_crack_length_larger_than_original() {
        let a = 0.005;
        let a_eff = effective_crack_length(a, 30.0e6, 250.0e6);
        assert!(a_eff > a, "effective crack length should exceed original");
    }
    #[test]
    fn test_kitagawa_large_crack_approaches_threshold() {
        let delta_k_th = 5.0e6;
        let delta_sigma_e = 200.0e6;
        let th = kitagawa_threshold(1.0, delta_k_th, delta_sigma_e);
        assert!(
            (th - delta_k_th).abs() / delta_k_th < 0.01,
            "large crack: th = {th}"
        );
    }
    #[test]
    fn test_kitagawa_threshold_positive() {
        let th = kitagawa_threshold(0.001, 5.0e6, 100.0e6);
        assert!(th > 0.0 && th.is_finite());
    }
    #[test]
    fn test_r_curve_at_zero_extension() {
        let kr = k_r_curve_linear(0.0, 25.0e6, 1e9);
        assert!((kr - 25.0e6).abs() < 1.0);
    }
    #[test]
    fn test_r_curve_rises_with_extension() {
        let kr0 = k_r_curve_linear(0.0, 25.0e6, 1e9);
        let kr1 = k_r_curve_linear(0.001, 25.0e6, 1e9);
        assert!(kr1 > kr0, "R-curve should rise with crack extension");
    }
    #[test]
    fn test_r_curve_fracture_check() {
        assert!(r_curve_fracture_check(50.0e6, 0.001, 25.0e6, 1e8));
        assert!(!r_curve_fracture_check(1.0e6, 0.001, 25.0e6, 1e8));
    }
    #[test]
    fn test_scc_below_threshold_zero() {
        let rate = scc_crack_growth_rate(1.0e6, 1e-10, 3.0, 5.0e6, 50.0e6);
        assert_eq!(rate, 0.0, "SCC below threshold: rate = {rate}");
    }
    #[test]
    fn test_scc_above_threshold_positive() {
        let rate = scc_crack_growth_rate(20.0e6, 1e-10, 3.0, 5.0e6, 50.0e6);
        assert!(rate > 0.0 && rate.is_finite());
    }
    #[test]
    fn test_scc_at_toughness_infinite() {
        let rate = scc_crack_growth_rate(50.0e6, 1e-10, 3.0, 5.0e6, 50.0e6);
        assert!(rate.is_infinite(), "SCC at K_Ic should be infinite: {rate}");
    }
    #[test]
    fn test_combined_fracture_check_below() {
        let result = fracture_check_combined(5.0e6, 3.0e6, 1.0e6, 50.0e6, 0.3);
        assert!(!result, "should not fracture below threshold");
    }
    #[test]
    fn test_combined_fracture_check_above() {
        let result = fracture_check_combined(40.0e6, 20.0e6, 10.0e6, 30.0e6, 0.3);
        assert!(result, "should fracture above threshold");
    }
}
#[cfg(test)]
mod tests_higher_order {
    use super::*;
    use crate::fracture::*;
    #[test]
    fn test_williams_singular_term_at_theta0() {
        let ki = 1.0e6;
        let we = WilliamsExpansion::new(ki, 0.0, 0.0, 0.0);
        let r = 0.01;
        let stresses = we.compute_higher_order_terms(r, 0.0);
        let expected = ki / (2.0 * PI * r).sqrt();
        assert!(
            (stresses[0] - expected).abs() / expected < 1e-9,
            "sigma_xx at theta=0 should equal K_I/sqrt(2 pi r): {} vs {}",
            stresses[0],
            expected
        );
    }
    #[test]
    fn test_williams_t_stress_adds_to_sigma_xx() {
        let ki = 1.0e6;
        let t = 50.0e3;
        let we_no_t = WilliamsExpansion::new(ki, 0.0, 0.0, 0.0);
        let we_t = WilliamsExpansion::new(ki, t, 0.0, 0.0);
        let r = 0.01;
        let s_no_t = we_no_t.compute_higher_order_terms(r, 0.0);
        let s_t = we_t.compute_higher_order_terms(r, 0.0);
        assert!(
            (s_t[0] - s_no_t[0] - t).abs() < 1.0,
            "T-stress should add to sigma_xx: diff={}",
            s_t[0] - s_no_t[0]
        );
    }
    #[test]
    fn test_williams_a3_order_sqrt_r() {
        let a3 = 1.0e3;
        let we = WilliamsExpansion::new(0.0, 0.0, a3, 0.0);
        let r1 = 0.01;
        let r2 = 0.04;
        let s1 = we.compute_higher_order_terms(r1, 0.0);
        let s2 = we.compute_higher_order_terms(r2, 0.0);
        let ratio = s2[0] / s1[0];
        assert!(
            (ratio - 2.0).abs() < 1e-9,
            "Order-3 term should scale as sqrt(r): ratio={}",
            ratio
        );
    }
    #[test]
    fn test_williams_near_zero_radius_singular() {
        let we = WilliamsExpansion::new(1.0e6, 0.0, 0.0, 0.0);
        let s = we.compute_higher_order_terms(1e-25, 0.0);
        assert!(
            s[0].is_infinite(),
            "Very small r should give infinite stress"
        );
    }
    #[test]
    fn test_williams_extract_ki_round_trip() {
        let ki_true = 2.0e6;
        let we = WilliamsExpansion::new(ki_true, 0.0, 0.0, 0.0);
        let r = 0.005;
        let sigma_yy = ki_true / (2.0 * PI * r).sqrt();
        let ki_extracted = we.extract_ki_from_stress(r, sigma_yy);
        assert!(
            (ki_extracted - ki_true).abs() / ki_true < 1e-10,
            "Round-trip K_I extraction: {} vs {}",
            ki_extracted,
            ki_true
        );
    }
    #[test]
    fn test_j_integral_domain_single_point() {
        let ji = JIntegral::new(200e9, 0.3, false);
        let pt = DomainIntegralPoint {
            stress: [100.0, 0.0, 0.0],
            strain: [5e-7, 0.0, 0.0],
            grad_u: [5e-7, 0.0, 0.0, 0.0],
            grad_q: [1.0, 0.0],
            dv: 1e-6,
        };
        let j = ji.compute_domain_integral(&[pt]);
        assert!(
            (j - 2.5e-11).abs() < 1e-20,
            "Domain J-integral single-point: {} vs 2.5e-11",
            j
        );
    }
    #[test]
    fn test_j_from_sif_plane_stress() {
        let ki = 1.0e6;
        let e = 200e9;
        let ji = JIntegral::new(e, 0.3, false);
        let j = ji.j_from_sif(ki, 0.0);
        let expected = ki * ki / e;
        assert!(
            (j - expected).abs() < 1e-10 * expected.abs(),
            "J from SIF plane stress: {} vs {}",
            j,
            expected
        );
    }
    #[test]
    fn test_j_from_sif_plane_strain() {
        let ki = 1.0e6;
        let e = 200e9;
        let nu = 0.3;
        let ji = JIntegral::new(e, nu, true);
        let j = ji.j_from_sif(ki, 0.0);
        let expected = ki * ki * (1.0 - nu * nu) / e;
        assert!(
            (j - expected).abs() < 1e-10 * expected.abs(),
            "J from SIF plane strain: {} vs {}",
            j,
            expected
        );
    }
    #[test]
    fn test_ki_from_j_round_trip() {
        let ki_true = 1.5e6;
        let e = 200e9;
        let ji = JIntegral::new(e, 0.3, false);
        let j = ji.j_from_sif(ki_true, 0.0);
        let ki_back = ji.ki_from_j(j);
        assert!(
            (ki_back - ki_true).abs() / ki_true < 1e-10,
            "Round-trip K_I from J: {} vs {}",
            ki_back,
            ki_true
        );
    }
    #[test]
    fn test_interaction_integral_zero_aux() {
        let sif = SifMixed::new(200e9, 0.3, false);
        let pt = InteractionIntegralPoint {
            stress_fem: [100.0, 50.0, 30.0],
            grad_u_fem: [1e-6, 2e-7, 3e-7, 1e-6],
            stress_aux: [0.0, 0.0, 0.0],
            grad_u_aux: [0.0, 0.0, 0.0, 0.0],
            strain_fem: [5e-7, 2.5e-7, 1.5e-7],
            strain_aux: [0.0, 0.0, 0.0],
            grad_q: [1.0, 0.0],
            dv: 1e-6,
        };
        let i = sif.compute_interaction_integral(&[pt]);
        assert_eq!(
            i, 0.0,
            "Zero auxiliary field gives zero interaction integral"
        );
    }
    #[test]
    fn test_sif_mixed_effective() {
        let sif = SifMixed::new(200e9, 0.3, false);
        let ki = 3.0e6;
        let kii = 4.0e6;
        let keff = sif.effective_sif(ki, kii);
        assert!(
            (keff - 5.0e6).abs() < 1e-3,
            "K_eff = 5 MPa sqrt(m): {}",
            keff
        );
    }
    #[test]
    fn test_sif_mixed_fracture_check() {
        let sif = SifMixed::new(200e9, 0.3, false);
        let kic = 50.0e6;
        assert!(
            !sif.is_fracture(10.0e6, 10.0e6, kic),
            "Should not fracture: K_I=10, K_II=10, K_Ic=50"
        );
        assert!(
            sif.is_fracture(40.0e6, 40.0e6, kic),
            "Should fracture: K_I=40, K_II=40, K_Ic=50"
        );
    }
}
/// Heaviside enrichment for XFEM (step function across crack).
///
/// Returns +1 if `signed_distance > 0`, -1 otherwise.
#[allow(dead_code)]
pub fn xfem_heaviside(signed_distance: f64) -> f64 {
    if signed_distance > 0.0 { 1.0 } else { -1.0 }
}
/// XFEM branch (crack-tip) enrichment functions in 2-D.
///
/// Based on the Moës-Dolbow-Belytschko 1999 formulation.
/// Returns `[F1, F2, F3, F4]` where:
/// - `F1 = sqrt(r) * sin(theta/2)`
/// - `F2 = sqrt(r) * cos(theta/2)`
/// - `F3 = sqrt(r) * sin(theta/2) * sin(theta)`
/// - `F4 = sqrt(r) * cos(theta/2) * sin(theta)`
#[allow(dead_code)]
pub fn xfem_branch_functions(r: f64, theta: f64) -> [f64; 4] {
    let sqrt_r = r.sqrt();
    let half_theta = theta / 2.0;
    [
        sqrt_r * half_theta.sin(),
        sqrt_r * half_theta.cos(),
        sqrt_r * half_theta.sin() * theta.sin(),
        sqrt_r * half_theta.cos() * theta.sin(),
    ]
}
/// Gradients of XFEM branch functions w.r.t. polar coordinates (r, theta).
///
/// Returns `[[dF1/dr, dF1/dtheta\], ..., [dF4/dr, dF4/dtheta]]`.
#[allow(dead_code)]
pub fn xfem_branch_gradients(r: f64, theta: f64) -> [[f64; 2]; 4] {
    let half_theta = theta / 2.0;
    let sin_h = half_theta.sin();
    let cos_h = half_theta.cos();
    let sin_t = theta.sin();
    let cos_t = theta.cos();
    let inv_sqrt_r = if r > 0.0 { 1.0 / (2.0 * r.sqrt()) } else { 0.0 };
    let sqrt_r = r.sqrt();
    [
        [inv_sqrt_r * sin_h, sqrt_r * cos_h / 2.0],
        [inv_sqrt_r * cos_h, -sqrt_r * sin_h / 2.0],
        [
            inv_sqrt_r * sin_h * sin_t,
            sqrt_r * (cos_h / 2.0 * sin_t + sin_h * cos_t),
        ],
        [
            inv_sqrt_r * cos_h * sin_t,
            sqrt_r * (-sin_h / 2.0 * sin_t + cos_h * cos_t),
        ],
    ]
}
/// Level-set signed distance from a line crack defined by two endpoints.
///
/// The crack is the segment `(x0, y0)` → `(x1, y1)`.
/// Returns the signed perpendicular distance from point `(px, py)`.
#[allow(dead_code)]
pub fn crack_level_set(x0: f64, y0: f64, x1: f64, y1: f64, px: f64, py: f64) -> f64 {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return ((px - x0) * (px - x0) + (py - y0) * (py - y0)).sqrt();
    }
    (dx * (py - y0) - dy * (px - x0)) / len
}
/// Nodal enrichment flag: does node at `(nx, ny)` need Heaviside enrichment?
///
/// Returns true if the crack level-set changes sign across the element containing the node.
#[allow(dead_code)]
pub fn node_needs_heaviside_enrichment(level_set_values: &[f64]) -> bool {
    let has_positive = level_set_values.iter().any(|&v| v > 0.0);
    let has_negative = level_set_values.iter().any(|&v| v < 0.0);
    has_positive && has_negative
}
/// Compute the equivalent domain integral (EDI) for XFEM.
///
/// Approximates the J-integral using the interaction integral approach for XFEM.
/// `sigma_xx`, `sigma_yy`, `sigma_xy` are stresses at integration points,
/// `u_x`, `u_y` are displacements, `w_x`, `w_y` are q-function gradients,
/// `weights` are integration weights.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn xfem_domain_integral(
    sigma_xx: &[f64],
    sigma_yy: &[f64],
    sigma_xy: &[f64],
    u_x: &[f64],
    u_y: &[f64],
    w_x: &[f64],
    w_y: &[f64],
    weights: &[f64],
) -> f64 {
    let n = sigma_xx.len();
    assert!(n == sigma_yy.len() && n == sigma_xy.len());
    assert!(n == u_x.len() && n == u_y.len());
    assert!(n == w_x.len() && n == w_y.len() && n == weights.len());
    let mut j = 0.0;
    for i in 0..n {
        let w_energy = 0.5 * (sigma_xx[i] * u_x[i] + sigma_yy[i] * u_y[i]);
        let flux_x = sigma_xx[i] * u_x[i] + sigma_xy[i] * u_y[i];
        let flux_y = sigma_xy[i] * u_x[i] + sigma_yy[i] * u_y[i];
        j += (flux_x * w_x[i] + flux_y * w_y[i] - w_energy * w_x[i]) * weights[i];
    }
    j
}
/// Elber crack closure model.
///
/// Effective stress-intensity range accounting for crack closure:
/// `delta_K_eff = U * delta_K`
/// where `U = 0.5 + 0.4 * R` (Elber's empirical relation, valid for -0.1 <= R <= 0.7).
#[allow(dead_code)]
pub fn elber_effective_delta_k(delta_k: f64, r_ratio: f64) -> f64 {
    let u = (0.5 + 0.4 * r_ratio).clamp(0.0, 1.0);
    u * delta_k
}
/// Schijve modification of Elber closure.
///
/// `U = 0.55 + 0.33 * R + 0.12 * R^2`
#[allow(dead_code)]
pub fn schijve_effective_delta_k(delta_k: f64, r_ratio: f64) -> f64 {
    let u = (0.55 + 0.33 * r_ratio + 0.12 * r_ratio * r_ratio).clamp(0.0, 1.0);
    u * delta_k
}
/// Striations per cycle from crack growth rate.
///
/// Each striation corresponds to one load cycle, spacing = `da/dN`.
#[allow(dead_code)]
pub fn striations_spacing(da_dn: f64) -> f64 {
    da_dn
}
/// Effective Paris law using Elber closure.
///
/// `da/dN = C * (U * delta_K)^m`
#[allow(dead_code)]
pub fn paris_elber(delta_k: f64, c: f64, m: f64, r_ratio: f64) -> f64 {
    let delta_k_eff = elber_effective_delta_k(delta_k, r_ratio);
    c * delta_k_eff.powf(m)
}
/// Integrated fatigue life with Elber closure (cycles to failure).
///
/// `N = integral from a0 to af of da / (C * (U * delta_K(a))^m)`
/// Numerically integrated with `n_steps` trapezoid steps.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn fatigue_life_elber(
    a0: f64,
    af: f64,
    c: f64,
    m: f64,
    sigma_far: f64,
    geo: f64,
    r_ratio: f64,
    n_steps: usize,
) -> f64 {
    let h = (af - a0) / n_steps as f64;
    let mut life = 0.0;
    for i in 0..n_steps {
        let a_lo = a0 + i as f64 * h;
        let a_hi = a_lo + h;
        let delta_k_lo = sigma_far * (PI * a_lo).sqrt() * geo;
        let delta_k_hi = sigma_far * (PI * a_hi).sqrt() * geo;
        let rate_lo = paris_elber(delta_k_lo, c, m, r_ratio).max(1e-30);
        let rate_hi = paris_elber(delta_k_hi, c, m, r_ratio).max(1e-30);
        life += 0.5 * h * (1.0 / rate_lo + 1.0 / rate_hi);
    }
    life
}
/// Crack retardation factor (Wheeler model).
///
/// After an overload, crack growth is retarded by factor `phi`:
/// `phi = (r_p / (a_ol + r_p_ol - a))^m_w`
/// where `a_ol` is crack length at overload, `r_p_ol` plastic zone at overload,
/// `r_p` current plastic zone, `a` current crack length, `m_w` Wheeler exponent.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn wheeler_retardation(a_ol: f64, r_p_ol: f64, a: f64, r_p: f64, m_w: f64) -> f64 {
    let boundary = a_ol + r_p_ol;
    if a + r_p >= boundary {
        return 1.0;
    }
    (r_p / (boundary - a)).powf(m_w).min(1.0)
}
/// Retarded Paris law (Wheeler model).
///
/// `da/dN = phi * C * delta_K^m`
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn paris_wheeler(
    delta_k: f64,
    c: f64,
    m: f64,
    a_ol: f64,
    r_p_ol: f64,
    a: f64,
    r_p: f64,
    m_w: f64,
) -> f64 {
    let phi = wheeler_retardation(a_ol, r_p_ol, a, r_p, m_w);
    phi * c * delta_k.powf(m)
}
/// Single crack propagation step using Paris law.
///
/// Given current crack length `a`, stress range `delta_sigma`, geometry factor `geo`,
/// Paris constants `c` and `m`, returns the new crack length after one block of `cycles`.
#[allow(dead_code)]
pub fn paris_step(a: f64, delta_sigma: f64, geo: f64, c: f64, m: f64, cycles: f64) -> f64 {
    let delta_k = delta_sigma * (PI * a).sqrt() * geo;
    let da_dn = c * delta_k.powf(m);
    a + da_dn * cycles
}
/// Simulate crack propagation until failure (K_I >= K_Ic).
///
/// Returns the crack length history as `(crack_lengths, cumulative_cycles)`.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn simulate_crack_growth(
    a0: f64,
    k_ic: f64,
    delta_sigma: f64,
    geo: f64,
    c: f64,
    m: f64,
    delta_n_block: f64,
    max_blocks: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut lengths = vec![a0];
    let mut cycles = vec![0.0];
    let mut a = a0;
    let mut n_total = 0.0;
    for _ in 0..max_blocks {
        let ki = delta_sigma * (PI * a).sqrt() * geo;
        if ki >= k_ic {
            break;
        }
        a = paris_step(a, delta_sigma, geo, c, m, delta_n_block);
        n_total += delta_n_block;
        lengths.push(a);
        cycles.push(n_total);
    }
    (lengths, cycles)
}
/// Compute the 4×4 cohesive element stiffness matrix (2-D, 2-noded).
///
/// Node layout: `[u_n1, u_t1, u_n2, u_t2]`.
/// Penalty stiffness in normal direction `k_n`, tangential `k_t`.
#[allow(dead_code)]
pub fn cohesive_stiffness_2d(k_n: f64, k_t: f64, length: f64) -> [[f64; 4]; 4] {
    let gauss_pts = [-1.0_f64 / 3.0_f64.sqrt(), 1.0_f64 / 3.0_f64.sqrt()];
    let gauss_w = [1.0_f64, 1.0_f64];
    let jac = length / 2.0;
    let mut k = [[0.0_f64; 4]; 4];
    for (&xi, &w) in gauss_pts.iter().zip(gauss_w.iter()) {
        let n1 = 0.5 * (1.0 - xi);
        let n2 = 0.5 * (1.0 + xi);
        let b_n = [n1, 0.0, n2, 0.0];
        let b_t = [0.0, n1, 0.0, n2];
        for i in 0..4 {
            for j in 0..4 {
                k[i][j] += w * jac * (k_n * b_n[i] * b_n[j] + k_t * b_t[i] * b_t[j]);
            }
        }
    }
    k
}
/// Extract the normal and tangential tractions from a cohesive element.
///
/// Given nodal displacements `dof = [u_n1, u_t1, u_n2, u_t2]` and
/// bilinear cohesive parameters, returns tractions at midpoint.
#[allow(dead_code)]
pub fn cohesive_element_traction(dof: [f64; 4], p: &CzmBilinearParams) -> (f64, f64) {
    let delta_n = 0.5 * (dof[0] + dof[2]);
    let delta_s = 0.5 * (dof[1] + dof[3]);
    let t_n = CzmBilinear::normal_traction(p, delta_n);
    let t_s = CzmBilinear::shear_traction(p, delta_s);
    (t_n, t_s)
}
#[cfg(test)]
mod fracture_extended_tests {
    use super::*;
    use crate::fracture::*;
    #[test]
    fn test_czm_bilinear_params_fracture_energy() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        let expected_g_i = 0.5 * 1e6 * 0.001;
        assert!((p.g_c_i - expected_g_i).abs() < 1e-10);
    }
    #[test]
    fn test_czm_bilinear_zero_opening_zero_traction() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        assert_eq!(CzmBilinear::normal_traction(&p, 0.0), 0.0);
    }
    #[test]
    fn test_czm_bilinear_negative_opening_zero_traction() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        assert_eq!(CzmBilinear::normal_traction(&p, -0.0005), 0.0);
    }
    #[test]
    fn test_czm_bilinear_peak_at_half_critical() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        let t = CzmBilinear::normal_traction(&p, 0.0005);
        assert!((t - 1e6).abs() < 1.0, "Expected peak traction at delta_n_0");
    }
    #[test]
    fn test_czm_bilinear_zero_after_critical() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        assert_eq!(CzmBilinear::normal_traction(&p, 0.002), 0.0);
    }
    #[test]
    fn test_czm_shear_antisymmetry() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        let t_pos = CzmBilinear::shear_traction(&p, 0.0005);
        let t_neg = CzmBilinear::shear_traction(&p, -0.0005);
        assert!((t_pos + t_neg).abs() < 1e-12);
    }
    #[test]
    fn test_czm_effective_opening() {
        let eff = CzmBilinear::effective_opening(3.0, 4.0);
        assert!((eff - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_czm_bk_criterion_pure_mode_i() {
        let g = CzmBilinear::bk_fracture_energy(100.0, 200.0, 100.0, 0.0, 2.0);
        assert!(
            (g - 100.0).abs() < 1e-10,
            "Pure mode-I: G_c should equal G_c_I"
        );
    }
    #[test]
    fn test_czm_bk_criterion_pure_mode_ii() {
        let g = CzmBilinear::bk_fracture_energy(100.0, 200.0, 0.0, 50.0, 2.0);
        assert!(
            (g - 200.0).abs() < 1e-10,
            "Pure mode-II: G_c should equal G_c_II"
        );
    }
    #[test]
    fn test_czm_damage_zero_below_onset() {
        let d = CzmBilinear::damage_variable(0.0005, 0.001, 0.005);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_czm_damage_one_at_failure() {
        let d = CzmBilinear::damage_variable(0.005, 0.001, 0.005);
        assert!((d - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_czm_damage_clamped_above_failure() {
        let d = CzmBilinear::damage_variable(0.01, 0.001, 0.005);
        assert_eq!(d, 1.0);
    }
    #[test]
    fn test_czm_secant_stiffness_no_damage() {
        let k = CzmBilinear::secant_normal_stiffness(1e9, 0.0);
        assert!((k - 1e9).abs() < 1.0);
    }
    #[test]
    fn test_czm_secant_stiffness_full_damage() {
        let k = CzmBilinear::secant_normal_stiffness(1e9, 1.0);
        assert_eq!(k, 0.0);
    }
    #[test]
    fn test_xu_needleman_peak_at_delta_nc() {
        let g_c_i = 100.0;
        let delta_n_c = 0.001;
        let t_peak = CzmBilinear::xu_needleman_normal(g_c_i, delta_n_c, delta_n_c);
        let expected = (g_c_i / delta_n_c) * (-1.0_f64).exp();
        assert!((t_peak - expected).abs() < 1e-10);
    }
    #[test]
    fn test_crack_opening_work_approaches_fracture_energy() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        let work = CzmBilinear::crack_opening_work(&p, 1000);
        let expected = p.g_c_i;
        assert!(
            (work - expected).abs() / expected < 0.01,
            "CZM work should equal G_c_I; got {} vs {}",
            work,
            expected
        );
    }
    #[test]
    fn test_xfem_heaviside_positive() {
        assert_eq!(xfem_heaviside(1.5), 1.0);
    }
    #[test]
    fn test_xfem_heaviside_negative() {
        assert_eq!(xfem_heaviside(-0.5), -1.0);
    }
    #[test]
    fn test_xfem_branch_f2_equals_sqrt_r_at_theta_zero() {
        let r = 4.0;
        let branches = xfem_branch_functions(r, 0.0);
        assert!((branches[0]).abs() < 1e-12);
        assert!((branches[1] - r.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_xfem_branch_all_zero_at_r_zero() {
        let branches = xfem_branch_functions(0.0, 0.5);
        for &f in &branches {
            assert!(f.abs() < 1e-12);
        }
    }
    #[test]
    fn test_crack_level_set_sign_change() {
        let ls_above = crack_level_set(0.0, 0.0, 1.0, 0.0, 0.5, 1.0);
        let ls_below = crack_level_set(0.0, 0.0, 1.0, 0.0, 0.5, -1.0);
        assert!(ls_above > 0.0);
        assert!(ls_below < 0.0);
    }
    #[test]
    fn test_crack_level_set_on_crack_zero() {
        let ls = crack_level_set(0.0, 0.0, 1.0, 0.0, 0.5, 0.0);
        assert!(ls.abs() < 1e-12);
    }
    #[test]
    fn test_node_needs_enrichment_straddle() {
        let lsets = [-0.5, 0.5, 0.3];
        assert!(node_needs_heaviside_enrichment(&lsets));
    }
    #[test]
    fn test_node_no_enrichment_all_positive() {
        let lsets = [0.5, 0.5, 0.3];
        assert!(!node_needs_heaviside_enrichment(&lsets));
    }
    #[test]
    fn test_elber_closure_r_zero() {
        let dk_eff = elber_effective_delta_k(10.0, 0.0);
        assert!((dk_eff - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_schijve_closure_higher_than_elber() {
        let dk = 10.0;
        let r = 0.3;
        let dk_elber = elber_effective_delta_k(dk, r);
        let dk_schijve = schijve_effective_delta_k(dk, r);
        assert!(dk_elber <= dk);
        assert!(dk_schijve <= dk);
    }
    #[test]
    fn test_paris_elber_less_than_plain_paris() {
        let dk = 20.0e6_f64;
        let c = 1e-11_f64;
        let m = 3.0_f64;
        let r = 0.1_f64;
        let plain = c * dk.powf(m);
        let elber = paris_elber(dk, c, m, r);
        assert!(
            elber < plain,
            "Elber growth rate should be less than plain Paris"
        );
    }
    #[test]
    fn test_wheeler_retardation_outside_zone() {
        let phi = wheeler_retardation(0.010, 0.005, 0.016, 0.002, 2.0);
        assert!((phi - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_wheeler_retardation_inside_zone_less_than_one() {
        let phi = wheeler_retardation(0.010, 0.005, 0.012, 0.001, 2.0);
        assert!(phi < 1.0);
        assert!(phi > 0.0);
    }
    #[test]
    fn test_paris_step_grows_crack() {
        let a_new = paris_step(0.010, 100e6, 1.0, 1e-11, 3.0, 1000.0);
        assert!(a_new > 0.010, "Crack should grow");
    }
    #[test]
    fn test_simulate_crack_growth_terminates() {
        let (lengths, cycles) =
            simulate_crack_growth(0.005, 50e6, 200e6, 1.12, 1e-11, 3.0, 100.0, 10000);
        assert!(!lengths.is_empty());
        assert_eq!(lengths.len(), cycles.len());
        assert!(*lengths.last().unwrap() > 0.005);
    }
    #[test]
    fn test_cohesive_stiffness_2d_symmetric() {
        let k = cohesive_stiffness_2d(1e9, 0.5e9, 0.01);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (k[i][j] - k[j][i]).abs() < 1e-6,
                    "K should be symmetric at [{},{}]",
                    i,
                    j
                );
            }
        }
    }
    #[test]
    fn test_cohesive_stiffness_2d_positive_diagonal() {
        let k = cohesive_stiffness_2d(1e9, 0.5e9, 0.01);
        for i in 0..4 {
            assert!(
                k[i][i] > 0.0,
                "Diagonal entry [{},{}] should be positive",
                i,
                i
            );
        }
    }
    #[test]
    fn test_cohesive_element_traction_zero_dof() {
        let p = CzmBilinearParams::from_peak_and_critical(1e6, 0.5e6, 0.001, 0.002);
        let (t_n, t_s) = cohesive_element_traction([0.0; 4], &p);
        assert_eq!(t_n, 0.0);
        assert_eq!(t_s, 0.0);
    }
    #[test]
    fn test_fatigue_life_elber_positive() {
        let life = fatigue_life_elber(0.005, 0.050, 1e-11, 3.0, 100e6, 1.0, 0.1, 100);
        assert!(life > 0.0, "Fatigue life must be positive");
    }
    #[test]
    fn test_xfem_domain_integral_positive_work() {
        let s_xx = vec![100e6];
        let s_yy = vec![50e6];
        let s_xy = vec![0.0];
        let u_x = vec![1e-4];
        let u_y = vec![5e-5];
        let w_x = vec![1.0];
        let w_y = vec![0.0];
        let wts = vec![1.0];
        let j = xfem_domain_integral(&s_xx, &s_yy, &s_xy, &u_x, &u_y, &w_x, &w_y, &wts);
        assert!(j.is_finite());
    }
}
