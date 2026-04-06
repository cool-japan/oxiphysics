//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ContactPair, FrictionStatus};

/// Estimate the contact stiffness from material properties using Hertz theory.
///
/// k_contact = 4/3 * E* * sqrt(R*)
///
/// This is the linearized contact stiffness at approach delta.
#[allow(dead_code)]
pub fn hertz_contact_stiffness(e_star: f64, r_star: f64, approach: f64) -> f64 {
    if approach <= 0.0 {
        return 0.0;
    }
    2.0 * e_star * (r_star * approach).sqrt()
}
/// Compute the contact energy for Hertz contact.
///
/// U = 8/15 * E* * sqrt(R*) * delta^(5/2)
#[allow(dead_code)]
pub fn hertz_contact_energy(e_star: f64, r_star: f64, approach: f64) -> f64 {
    if approach <= 0.0 {
        return 0.0;
    }
    (8.0 / 15.0) * e_star * r_star.sqrt() * approach.powf(2.5)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::contact::*;
    use std::f64::consts::PI;
    #[test]
    fn test_hertz_sphere_sphere_steel() {
        let r = 0.01_f64;
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let approach = 1.0e-6_f64;
        let res = HertzContact::sphere_sphere(r, r, e, nu, e, nu, approach);
        let e_star = e / (2.0 * (1.0 - nu * nu));
        let r_star = r / 2.0;
        let expected_a = (r_star * approach).sqrt();
        let expected_f = (4.0 / 3.0) * e_star * r_star.sqrt() * approach.powf(1.5);
        let expected_p0 = 3.0 * expected_f / (2.0 * PI * expected_a * expected_a);
        assert!(
            (res.contact_radius - expected_a).abs() / expected_a < 1e-12,
            "contact_radius mismatch"
        );
        assert!(
            (res.force - expected_f).abs() / expected_f < 1e-12,
            "force mismatch"
        );
        assert!(
            (res.peak_pressure - expected_p0).abs() / expected_p0 < 1e-12,
            "peak_pressure mismatch"
        );
        assert!(res.force > 0.0);
        assert!(res.peak_pressure > 0.0);
        assert!(res.contact_radius > 0.0);
    }
    #[test]
    fn test_hertz_force_displacement() {
        let r = 0.01_f64;
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let d1 = 1.0e-6_f64;
        let d2 = 2.0e-6_f64;
        let f1 = HertzContact::sphere_sphere(r, r, e, nu, e, nu, d1).force;
        let f2 = HertzContact::sphere_sphere(r, r, e, nu, e, nu, d2).force;
        let ratio = f2 / f1;
        let expected_ratio = (d2 / d1).powf(1.5);
        assert!((ratio - expected_ratio).abs() / expected_ratio < 1e-12);
    }
    #[test]
    fn test_hertz_inverse() {
        let r1 = 0.01_f64;
        let r2 = 0.02_f64;
        let e1 = 200.0e9_f64;
        let nu1 = 0.3_f64;
        let e2 = 70.0e9_f64;
        let nu2 = 0.33_f64;
        let approach_original = 2.5e-6_f64;
        let res_fwd = HertzContact::sphere_sphere(r1, r2, e1, nu1, e2, nu2, approach_original);
        let approach_recovered = ContactPenetrationDepth::from_force_sphere_sphere(
            res_fwd.force,
            r1,
            r2,
            e1,
            nu1,
            e2,
            nu2,
        );
        assert!((approach_recovered - approach_original).abs() / approach_original < 1e-12,);
    }
    #[test]
    fn test_hertz_sphere_flat_special_case() {
        let r = 0.01_f64;
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let approach = 1.0e-6_f64;
        let r2_large = 1.0e12_f64;
        let res_flat = HertzContact::sphere_flat(r, e, nu, e, nu, approach);
        let res_ss = HertzContact::sphere_sphere(r, r2_large, e, nu, e, nu, approach);
        assert!(
            (res_flat.contact_radius - res_ss.contact_radius).abs() / res_flat.contact_radius
                < 1e-6,
        );
        assert!((res_flat.force - res_ss.force).abs() / res_flat.force < 1e-6,);
    }
    #[test]
    fn test_penalty_parameters_from_material() {
        let params = PenaltyParameters::from_material(200.0e9, 0.01, 100.0);
        assert!(params.normal_stiffness > 0.0);
        assert!(params.tangential_stiffness > 0.0);
        assert!(
            params.tangential_stiffness < params.normal_stiffness,
            "kt should be less than kn"
        );
    }
    #[test]
    fn test_penalty_adaptive_stiffness() {
        let params = PenaltyParameters::from_material(200.0e9, 0.01, 100.0);
        let k_normal = params.adaptive_stiffness(0.0005);
        let k_high = params.adaptive_stiffness(0.002);
        assert!(
            k_high > k_normal,
            "should increase stiffness for excessive penetration"
        );
        assert!((k_high / k_normal - params.scaling_factor).abs() < 1e-10);
    }
    #[test]
    fn test_penalty_normal_force() {
        let params = PenaltyParameters::from_material(200.0e9, 0.01, 100.0);
        assert_eq!(params.normal_force(0.001), 0.0);
        let f = params.normal_force(-0.0005);
        assert!(f > 0.0);
    }
    #[test]
    fn test_penalty_friction_force() {
        let params = PenaltyParameters::from_material(200.0e9, 0.01, 100.0);
        let f_n = 100.0;
        let mu = 0.3;
        let f_t = params.friction_force(f_n, 0.0, mu);
        assert!(f_t.abs() < 1e-10);
        let f_t_full = params.friction_force(f_n, 1.0, mu);
        assert!((f_t_full - mu * f_n).abs() / (mu * f_n) < 0.01);
    }
    #[test]
    fn test_augmented_lagrangian_uzawa_convergence() {
        let mut solver = AugmentedLagrangianSolver::new(1000.0, 2, 1e-6, 100);
        let gaps = vec![-0.01, 0.05];
        let slidings = vec![0.001, 0.0];
        for _ in 0..20 {
            solver.uzawa_update(&gaps, &slidings, 0.3);
        }
        assert!(
            solver.lambda_n[0] < 0.0,
            "should have compressive multiplier"
        );
        assert_eq!(
            solver.lambda_n[1], 0.0,
            "separated pair should have zero multiplier"
        );
    }
    #[test]
    fn test_augmented_lagrangian_normal_force() {
        let mut solver = AugmentedLagrangianSolver::new(1000.0, 1, 1e-6, 100);
        solver.lambda_n[0] = -5.0;
        let f = solver.normal_force(0, -0.01);
        assert!((f - (-15.0)).abs() < 1e-10);
    }
    #[test]
    fn test_segment_distance_parallel() {
        let (dist, _t1, _t2) = SegmentToSegmentContact::segment_distance(
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        );
        assert!(
            (dist - 1.0).abs() < 1e-10,
            "distance = {dist}, expected 1.0"
        );
    }
    #[test]
    fn test_segment_distance_touching() {
        let (dist, _t1, _t2) = SegmentToSegmentContact::segment_distance(
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
        );
        assert!(
            dist < 1e-10,
            "touching segments should have zero distance, got {dist}"
        );
    }
    #[test]
    fn test_mortar_d_integral() {
        let d = MortarContact::mortar_d_integral(0.5, 2.0);
        assert!((d[0] - 1.0).abs() < 1e-12);
        assert!((d[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mortar_m_integral() {
        let m = MortarContact::mortar_m_integral(6.0);
        assert!((m[0][0] - 2.0).abs() < 1e-12);
        assert!((m[0][1] - 1.0).abs() < 1e-12);
        assert!((m[1][0] - 1.0).abs() < 1e-12);
        assert!((m[1][1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_mortar_project_slave() {
        let (xi, gap) = MortarContact::project_slave_to_master([0.5, 1.0], [0.0, 0.0], [1.0, 0.0]);
        assert!((xi - 0.5).abs() < 1e-12, "xi = {xi}");
        assert!((gap - 1.0).abs() < 1e-12, "gap = {gap}");
    }
    #[test]
    fn test_mortar_project_clamped() {
        let (xi, _gap) = MortarContact::project_slave_to_master([2.0, 0.0], [0.0, 0.0], [1.0, 0.0]);
        assert!(
            (xi - 1.0).abs() < 1e-12,
            "xi should be clamped to 1.0, got {xi}"
        );
    }
    #[test]
    fn test_hertz_contact_stiffness() {
        let e = 200.0e9;
        let nu = 0.3;
        let e_star = e / (2.0 * (1.0 - nu * nu));
        let r_star = 0.005;
        let approach = 1e-6;
        let k = hertz_contact_stiffness(e_star, r_star, approach);
        assert!(k > 0.0);
        let k2 = hertz_contact_stiffness(e_star, r_star, 2.0 * approach);
        assert!(k2 > k, "stiffness should increase with approach");
    }
    #[test]
    fn test_hertz_contact_energy() {
        let e_star = 100.0e9;
        let r_star = 0.005;
        let approach = 1e-6;
        let u = hertz_contact_energy(e_star, r_star, approach);
        assert!(u > 0.0, "energy should be positive");
        assert_eq!(hertz_contact_energy(e_star, r_star, 0.0), 0.0);
        assert_eq!(hertz_contact_energy(e_star, r_star, -1e-6), 0.0);
    }
    #[test]
    fn test_hertz_contact_energy_scales() {
        let e_star = 100.0e9;
        let r_star = 0.005;
        let d1 = 1e-6;
        let d2 = 2e-6;
        let u1 = hertz_contact_energy(e_star, r_star, d1);
        let u2 = hertz_contact_energy(e_star, r_star, d2);
        let ratio = u2 / u1;
        let expected = (d2 / d1).powf(2.5);
        assert!((ratio - expected).abs() / expected < 1e-10);
    }
}
#[cfg(test)]
mod tests_extended {

    use crate::contact::*;
    #[test]
    fn test_nts_project_midpoint() {
        let slave = [0.5, 1.0, 0.0];
        let (xi, gap) = NodeToSegmentContact::project_node_to_segment_3d(
            slave,
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        );
        assert!((xi - 0.5).abs() < 1e-12, "xi = {xi}");
        assert!((gap[1] - 1.0).abs() < 1e-12, "gap_y = {}", gap[1]);
    }
    #[test]
    fn test_nts_project_beyond_endpoint() {
        let slave = [2.0, 0.5, 0.0];
        let (xi, _gap) = NodeToSegmentContact::project_node_to_segment_3d(
            slave,
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        );
        assert!((xi - 1.0).abs() < 1e-12, "xi should clamp to 1.0, got {xi}");
    }
    #[test]
    fn test_nts_signed_gap_separated() {
        let gap = NodeToSegmentContact::signed_gap(
            [0.5, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!(gap > 0.0, "separated node should have positive gap: {gap}");
    }
    #[test]
    fn test_nts_signed_gap_penetrating() {
        let gap = NodeToSegmentContact::signed_gap(
            [0.5, -0.1, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!(
            gap < 0.0,
            "penetrating node should have negative gap: {gap}"
        );
    }
    #[test]
    fn test_nts_penalty_force_no_penetration() {
        let f = NodeToSegmentContact::penalty_force(
            [0.5, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            1e6,
        );
        assert!(
            f[0].abs() < 1e-12 && f[1].abs() < 1e-12,
            "no force when separated"
        );
    }
    #[test]
    fn test_nts_penalty_force_penetration() {
        let f = NodeToSegmentContact::penalty_force(
            [0.5, -0.1, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            1e6,
        );
        assert!(
            f[1] > 0.0,
            "penetrating node should get restoring force: {}",
            f[1]
        );
    }
    #[test]
    fn test_aabb_overlaps() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        let c = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(a.overlaps(&b), "overlapping AABBs should return true");
        assert!(!a.overlaps(&c), "disjoint AABBs should return false");
    }
    #[test]
    fn test_aabb_contains() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(a.contains([0.5, 0.5, 0.5]));
        assert!(!a.contains([1.5, 0.5, 0.5]));
    }
    #[test]
    fn test_aabb_from_points() {
        let pts = vec![[1.0, 2.0, 3.0], [-1.0, 0.0, 5.0], [0.0, 4.0, -2.0]];
        let aabb = Aabb::from_points(&pts).unwrap();
        assert!((aabb.min[0] - (-1.0)).abs() < 1e-12);
        assert!((aabb.max[0] - 1.0).abs() < 1e-12);
        assert!((aabb.min[2] - (-2.0)).abs() < 1e-12);
        assert!((aabb.max[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_fem_contact_detector_candidates() {
        let detector = FemContactDetector::new(0.5);
        let slaves = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let masters = vec![[0.3, 0.0, 0.0], [5.0, 0.5, 0.0]];
        let pairs = detector.find_candidates(&slaves, &masters);
        assert!(!pairs.is_empty(), "should find at least one candidate pair");
    }
    #[test]
    fn test_fem_contact_detector_no_candidates() {
        let detector = FemContactDetector::new(0.1);
        let slaves = vec![[0.0, 0.0, 0.0]];
        let masters = vec![[5.0, 0.0, 0.0]];
        let pairs = detector.find_candidates(&slaves, &masters);
        assert!(pairs.is_empty(), "far-away nodes should have no candidates");
    }
    #[test]
    fn test_hertz_validator_force() {
        let e = 200.0e9;
        let nu = 0.3;
        let r = 0.01;
        let approach = 1e-6;
        let e_star = e / (2.0 * (1.0 - nu * nu));
        let r_star = r / 2.0;
        let result = HertzContact::sphere_sphere(r, r, e, nu, e, nu, approach);
        let err = HertzValidator::validate_force(&result, e_star, r_star);
        assert!(err < 1e-10, "force validation error = {err}");
    }
    #[test]
    fn test_hertz_validator_contact_radius() {
        let e = 200.0e9;
        let nu = 0.3;
        let r = 0.01;
        let approach = 1e-6;
        let r_star = r / 2.0;
        let result = HertzContact::sphere_sphere(r, r, e, nu, e, nu, approach);
        let err = HertzValidator::validate_contact_radius(&result, r_star);
        assert!(err < 1e-10, "contact radius validation error = {err}");
    }
    #[test]
    fn test_hertz_validator_pressure_ratio() {
        let e = 200.0e9;
        let nu = 0.3;
        let r = 0.01;
        let approach = 1e-6;
        let result = HertzContact::sphere_sphere(r, r, e, nu, e, nu, approach);
        let err = HertzValidator::validate_pressure_ratio(&result);
        assert!(err < 1e-10, "pressure ratio error = {err}");
    }
    #[test]
    fn test_hertz_validator_peak_pressure() {
        let e = 200.0e9;
        let nu = 0.3;
        let r = 0.01;
        let approach = 1e-6;
        let result = HertzContact::sphere_sphere(r, r, e, nu, e, nu, approach);
        let err = HertzValidator::validate_peak_pressure(&result);
        assert!(err < 1e-10, "peak pressure validation error = {err}");
    }
}
/// Detect node-to-face contact by checking each node against every triangular
/// face.  Returns all pairs where the signed distance is negative (penetrating)
/// or within a small tolerance.
#[allow(dead_code)]
pub fn detect_node_face_contact(nodes: &[[f64; 3]], faces: &[[usize; 3]]) -> Vec<ContactPair> {
    let mut pairs = Vec::new();
    let tol = 1e-10_f64;
    for (node_idx, &node) in nodes.iter().enumerate() {
        for face in faces.iter() {
            let v0 = nodes[face[0]];
            let v1 = nodes[face[1]];
            let v2 = nodes[face[2]];
            let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let nx = e1[1] * e2[2] - e1[2] * e2[1];
            let ny = e1[2] * e2[0] - e1[0] * e2[2];
            let nz = e1[0] * e2[1] - e1[1] * e2[0];
            let n_len = (nx * nx + ny * ny + nz * nz).sqrt();
            if n_len < 1e-30 {
                continue;
            }
            let normal = [nx / n_len, ny / n_len, nz / n_len];
            let d = [node[0] - v0[0], node[1] - v0[1], node[2] - v0[2]];
            let gap = d[0] * normal[0] + d[1] * normal[1] + d[2] * normal[2];
            if gap <= tol {
                let proj = [
                    node[0] - gap * normal[0],
                    node[1] - gap * normal[1],
                    node[2] - gap * normal[2],
                ];
                let w0 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
                let w1 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
                let w2 = [proj[0] - v0[0], proj[1] - v0[1], proj[2] - v0[2]];
                let d00 = w0[0] * w0[0] + w0[1] * w0[1] + w0[2] * w0[2];
                let d01 = w0[0] * w1[0] + w0[1] * w1[1] + w0[2] * w1[2];
                let d11 = w1[0] * w1[0] + w1[1] * w1[1] + w1[2] * w1[2];
                let d20 = w2[0] * w0[0] + w2[1] * w0[1] + w2[2] * w0[2];
                let d21 = w2[0] * w1[0] + w2[1] * w1[1] + w2[2] * w1[2];
                let denom = d00 * d11 - d01 * d01;
                if denom.abs() < 1e-30 {
                    continue;
                }
                let v = (d11 * d20 - d01 * d21) / denom;
                let w = (d00 * d21 - d01 * d20) / denom;
                let u = 1.0 - v - w;
                if u >= -tol && v >= -tol && w >= -tol {
                    pairs.push(ContactPair {
                        node_a: node_idx,
                        node_b: face[0],
                        gap,
                        normal,
                        lambda: 0.0,
                    });
                }
            }
        }
    }
    pairs
}
/// Hertz contact force for two elastic spheres given effective modulus,
/// effective radius, and indentation depth `delta`.
///
/// F = (4/3) * E_eff * sqrt(R_eff) * delta^(3/2)
#[allow(dead_code)]
pub fn hertz_contact_force(
    e1: f64,
    nu1: f64,
    e2: f64,
    nu2: f64,
    r1: f64,
    r2: f64,
    delta: f64,
) -> f64 {
    if delta <= 0.0 {
        return 0.0;
    }
    let inv_e = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
    let e_eff = 1.0 / inv_e;
    let r_eff = 1.0 / (1.0 / r1 + 1.0 / r2);
    (4.0 / 3.0) * e_eff * r_eff.sqrt() * delta.powf(1.5)
}
/// Contact radius from Hertz theory given effective modulus, effective radius,
/// and applied force F.
///
/// a = ( (3 * F) / (4 * E_eff * sqrt(R_eff)) )^(2/3)
#[allow(dead_code)]
pub fn hertz_contact_radius(e_eff: f64, r_eff: f64, f: f64) -> f64 {
    if f <= 0.0 || e_eff <= 0.0 || r_eff <= 0.0 {
        return 0.0;
    }
    let delta = ((3.0 * f) / (4.0 * e_eff * r_eff.sqrt())).powf(2.0 / 3.0);
    (r_eff * delta).sqrt()
}
#[cfg(test)]
mod tests_contact_new {
    use super::*;
    use crate::contact::*;
    #[test]
    fn test_penalty_contact_zero_at_positive_gap() {
        let pc = PenaltyContact {
            penalty: 1e6,
            friction_coeff: 0.3,
        };
        assert_eq!(
            pc.normal_force(0.001),
            0.0,
            "positive gap should give zero force"
        );
        assert_eq!(pc.normal_force(0.0), 0.0, "zero gap should give zero force");
    }
    #[test]
    fn test_penalty_contact_nonzero_at_negative_gap() {
        let pc = PenaltyContact {
            penalty: 1e6,
            friction_coeff: 0.3,
        };
        let f = pc.normal_force(-0.001);
        assert!(f > 0.0, "negative gap should give positive force, got {f}");
        assert!((f - 1e3).abs() < 1e-6, "expected 1000 N, got {f}");
    }
    #[test]
    fn test_al_contact_multiplier_update_violated() {
        let mut alc = AugmentedLagrangianContact::new(1000.0, 2);
        let gaps = [-0.01, 0.05];
        alc.update_multipliers(&gaps);
        assert!(alc.lambda[0] > 0.0);
    }
    #[test]
    fn test_hertz_contact_force_positive() {
        let e = 200.0e9;
        let nu = 0.3;
        let r = 0.01;
        let delta = 1e-6;
        let f = hertz_contact_force(e, nu, e, nu, r, r, delta);
        assert!(
            f > 0.0,
            "Hertz force should be positive for delta>0, got {f}"
        );
        assert_eq!(hertz_contact_force(e, nu, e, nu, r, r, 0.0), 0.0);
    }
    #[test]
    fn test_contact_stiffness_triplets_count() {
        let pairs = vec![ContactPair {
            node_a: 0,
            node_b: 1,
            gap: -0.001,
            normal: [0.0, 1.0, 0.0],
            lambda: 0.0,
        }];
        let triplets = ContactStiffness::assemble(&pairs, 1e6, 6);
        assert_eq!(
            triplets.len(),
            12,
            "expected 12 triplets for one pair, got {}",
            triplets.len()
        );
    }
    #[test]
    fn test_contact_stiffness_no_triplets_for_open_pair() {
        let pairs = vec![ContactPair {
            node_a: 0,
            node_b: 1,
            gap: 0.001,
            normal: [0.0, 1.0, 0.0],
            lambda: 0.0,
        }];
        let triplets = ContactStiffness::assemble(&pairs, 1e6, 6);
        assert!(
            triplets.is_empty(),
            "open contact pair should produce no triplets"
        );
    }
}
/// Evaluate Coulomb friction status.
///
/// Returns [`FrictionStatus`] based on gap, tangential trial force, and friction limit.
#[allow(dead_code)]
pub fn coulomb_status(
    gap: f64,
    tangential_trial: f64,
    normal_force: f64,
    mu: f64,
) -> FrictionStatus {
    if gap >= 0.0 {
        return FrictionStatus::Open;
    }
    let limit = mu * normal_force.abs();
    if tangential_trial.abs() <= limit {
        FrictionStatus::Stick
    } else {
        FrictionStatus::Slip
    }
}
/// Regularized Coulomb friction using a hyperbolic tangent smoothing.
///
/// `t_f = mu * |f_n| * tanh(slip_vel / epsilon)`
#[allow(dead_code)]
pub fn regularized_friction_force(mu: f64, normal_force: f64, slip_vel: f64, epsilon: f64) -> f64 {
    mu * normal_force.abs() * (slip_vel / epsilon).tanh()
}
/// Penalty-based tangential contact force (stick-slip).
///
/// If the penalty tangential force exceeds `mu * |f_n|`, it is capped (slip).
#[allow(dead_code)]
pub fn penalty_tangential_force(slip: f64, k_t: f64, normal_force: f64, mu: f64) -> f64 {
    let trial = k_t * slip;
    let limit = mu * normal_force.abs();
    if trial.abs() <= limit {
        trial
    } else {
        limit * trial.signum()
    }
}
/// Return mapping for Coulomb friction (2-D slip direction).
///
/// Updates trial tangential force by projecting onto Coulomb cone.
/// Returns `(corrected_slip, corrected_force, slipped)`.
#[allow(dead_code)]
pub fn return_mapping_coulomb_2d(
    trial_t: [f64; 2],
    normal_force: f64,
    mu: f64,
) -> ([f64; 2], [f64; 2], bool) {
    let limit = mu * normal_force.abs();
    let norm = (trial_t[0] * trial_t[0] + trial_t[1] * trial_t[1]).sqrt();
    if norm <= limit {
        (trial_t, trial_t, false)
    } else {
        let scale = limit / norm;
        let corrected = [trial_t[0] * scale, trial_t[1] * scale];
        (trial_t, corrected, true)
    }
}
/// 2-D gap function between a node and a line segment.
///
/// Returns the signed perpendicular distance from point `p` to segment `(a, b)`.
/// Negative means penetration (node on wrong side).
#[allow(dead_code)]
pub fn gap_node_to_segment_2d(
    p: [f64; 2],
    seg_a: [f64; 2],
    seg_b: [f64; 2],
    inward_normal: [f64; 2],
) -> f64 {
    let dx = seg_b[0] - seg_a[0];
    let dy = seg_b[1] - seg_a[1];
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return 0.0;
    }
    let nx = -dy / len;
    let ny = dx / len;
    let sign = if nx * inward_normal[0] + ny * inward_normal[1] < 0.0 {
        -1.0
    } else {
        1.0
    };
    sign * ((p[0] - seg_a[0]) * nx + (p[1] - seg_a[1]) * ny)
}
/// Maximum shear stress at depth `z` beneath a circular Hertz contact.
///
/// Based on the Hertz elastic half-space solution:
/// `tau_max(z) = p0 * f(z/a)` where `f` is a tabulated function.
/// Here we use the approximate formula by Johnson (1985):
/// `tau_max ≈ 0.31 * p0` at `z ≈ 0.47 * a`.
#[allow(dead_code)]
pub fn hertz_max_shear_stress(p0: f64) -> f64 {
    0.31 * p0
}
/// Depth of maximum shear stress beneath circular Hertz contact.
#[allow(dead_code)]
pub fn hertz_max_shear_depth(a: f64) -> f64 {
    0.47 * a
}
/// Hertz subsurface stress field (approximate) at depth `z`, radius `r` from axis.
///
/// Returns `[sigma_r, sigma_z, tau_rz]` using the Sneddon approximation.
#[allow(dead_code)]
pub fn hertz_subsurface_stress(p0: f64, a: f64, z: f64, r: f64) -> [f64; 3] {
    let zeta = z / a;
    let rho = r / a;
    let denom = (rho * rho + (1.0 + zeta) * (1.0 + zeta)).sqrt();
    let sigma_r = -p0 * (1.0 + 2.0 * zeta * zeta / (denom * denom)).recip();
    let sigma_z = -p0 / (1.0 + zeta * zeta / (denom * denom * denom));
    let tau_rz = p0 * rho * zeta / (denom * denom * denom);
    [sigma_r, sigma_z, tau_rz]
}
/// Mindlin tangential compliance for a sphere-sphere contact.
///
/// Returns the tangential displacement `delta_t` under tangential force `Q`.
/// `Q < mu * P` is assumed (no slip).
#[allow(dead_code)]
pub fn mindlin_tangential_compliance(q: f64, _p: f64, a: f64, e_star: f64, _nu: f64) -> f64 {
    3.0 * q / (16.0 * e_star * a)
}
/// Rolling resistance torque (simple model).
///
/// `M = mu_r * P * R`
#[allow(dead_code)]
pub fn rolling_resistance_torque(mu_r: f64, normal_force: f64, r: f64) -> f64 {
    mu_r * normal_force * r
}
#[cfg(test)]
mod contact_extended_tests {
    use super::*;
    use crate::contact::*;
    #[test]
    fn test_coulomb_status_open_gap() {
        let status = coulomb_status(0.001, 500.0, 1000.0, 0.3);
        assert_eq!(status, FrictionStatus::Open);
    }
    #[test]
    fn test_coulomb_status_stick() {
        let status = coulomb_status(-0.001, 100.0, 1000.0, 0.3);
        assert_eq!(status, FrictionStatus::Stick);
    }
    #[test]
    fn test_coulomb_status_slip() {
        let status = coulomb_status(-0.001, 400.0, 1000.0, 0.3);
        assert_eq!(status, FrictionStatus::Slip);
    }
    #[test]
    fn test_regularized_friction_zero_velocity() {
        let f = regularized_friction_force(0.3, 1000.0, 0.0, 0.001);
        assert_eq!(f, 0.0);
    }
    #[test]
    fn test_regularized_friction_saturates() {
        let f = regularized_friction_force(0.3, 1000.0, 1e6, 0.001);
        assert!((f - 0.3 * 1000.0).abs() < 1.0);
    }
    #[test]
    fn test_penalty_tangential_stick() {
        let f = penalty_tangential_force(0.00001, 1e8, 1000.0, 0.3);
        assert!((f.abs() - 300.0).abs() < 1e-6);
    }
    #[test]
    fn test_penalty_tangential_slip_capped() {
        let f = penalty_tangential_force(0.001, 1e8, 500.0, 0.3);
        let limit = 0.3 * 500.0;
        assert!((f.abs() - limit).abs() < 1e-6);
    }
    #[test]
    fn test_return_mapping_no_slip() {
        let (_trial, corrected, slipped) = return_mapping_coulomb_2d([50.0, 0.0], 500.0, 0.3);
        assert!(!slipped);
        assert!((corrected[0] - 50.0).abs() < 1e-12);
    }
    #[test]
    fn test_return_mapping_slip_reduces_force() {
        let (_trial, corrected, slipped) = return_mapping_coulomb_2d([200.0, 0.0], 500.0, 0.3);
        assert!(slipped);
        assert!((corrected[0] - 150.0).abs() < 1e-10);
    }
    #[test]
    fn test_return_mapping_2d_direction_preserved() {
        let trial = [120.0_f64, 160.0_f64];
        let (_, corrected, slipped) = return_mapping_coulomb_2d(trial, 500.0, 0.3);
        assert!(slipped);
        let norm = (corrected[0] * corrected[0] + corrected[1] * corrected[1]).sqrt();
        assert!((norm - 150.0).abs() < 1e-8);
        let trial_norm = (trial[0] * trial[0] + trial[1] * trial[1]).sqrt();
        assert!((corrected[0] / norm - trial[0] / trial_norm).abs() < 1e-10);
    }
    #[test]
    fn test_elliptical_hertz_semi_axes_positive() {
        let e_star = 100e9_f64 / (2.0 * (1.0 - 0.09));
        let result = EllipticalHertz::contact(0.02, 0.02, e_star, 1000.0);
        assert!(result.semi_axis_a > 0.0);
        assert!(result.semi_axis_b > 0.0);
    }
    #[test]
    fn test_elliptical_hertz_pressure_positive() {
        let result = EllipticalHertz::contact(0.02, 0.02, 50e9, 1000.0);
        assert!(result.peak_pressure > 0.0);
    }
    #[test]
    fn test_elliptical_pressure_at_center() {
        let p = EllipticalHertz::pressure_distribution(0.0, 0.0, 100e6, 0.005, 0.003);
        assert!((p - 100e6).abs() < 1.0);
    }
    #[test]
    fn test_elliptical_pressure_zero_outside() {
        let p = EllipticalHertz::pressure_distribution(0.01, 0.01, 100e6, 0.005, 0.003);
        assert_eq!(p, 0.0);
    }
    #[test]
    fn test_aabb2d_overlap() {
        let a = Aabb2d::new([0.0, 0.0], [1.0, 1.0]);
        let b = Aabb2d::new([0.5, 0.5], [1.5, 1.5]);
        assert!(a.overlaps(&b));
    }
    #[test]
    fn test_aabb2d_no_overlap() {
        let a = Aabb2d::new([0.0, 0.0], [1.0, 1.0]);
        let b = Aabb2d::new([2.0, 2.0], [3.0, 3.0]);
        assert!(!a.overlaps(&b));
    }
    #[test]
    fn test_aabb2d_from_points() {
        let pts = [[0.0, 1.0], [2.0, -1.0], [1.0, 3.0]];
        let bb = Aabb2d::from_points(&pts).unwrap();
        assert!((bb.min[0] - 0.0).abs() < 1e-12);
        assert!((bb.min[1] - (-1.0)).abs() < 1e-12);
        assert!((bb.max[0] - 2.0).abs() < 1e-12);
        assert!((bb.max[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb2d_center() {
        let a = Aabb2d::new([0.0, 0.0], [2.0, 4.0]);
        let c = a.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb2d_expanded() {
        let a = Aabb2d::new([1.0, 1.0], [3.0, 3.0]);
        let b = a.expanded(0.5);
        assert!((b.min[0] - 0.5).abs() < 1e-12);
        assert!((b.max[0] - 3.5).abs() < 1e-12);
    }
    #[test]
    fn test_gap_node_to_segment_above() {
        let p = [0.5, 1.0];
        let seg_a = [0.0, 0.0];
        let seg_b = [1.0, 0.0];
        let inward = [0.0, 1.0];
        let g = gap_node_to_segment_2d(p, seg_a, seg_b, inward);
        assert!(
            g > 0.0,
            "Node above segment should have positive gap: {}",
            g
        );
    }
    #[test]
    fn test_gap_node_to_segment_on_segment_zero() {
        let p = [0.5, 0.0];
        let seg_a = [0.0, 0.0];
        let seg_b = [1.0, 0.0];
        let inward = [0.0, 1.0];
        let g = gap_node_to_segment_2d(p, seg_a, seg_b, inward);
        assert!(
            g.abs() < 1e-12,
            "Node on segment should have zero gap: {}",
            g
        );
    }
    #[test]
    fn test_dual_state_uzawa_negative_gap_activates() {
        let mut state = DualContactState::new(2, 1e6);
        state.uzawa_step(&[-0.001, -0.002]);
        assert!(state.lambda[0] < 0.0);
        assert!(state.lambda[1] < 0.0);
    }
    #[test]
    fn test_dual_state_uzawa_positive_gap_stays_zero() {
        let mut state = DualContactState::new(2, 1e6);
        state.uzawa_step(&[0.001, 0.002]);
        assert_eq!(state.lambda[0], 0.0);
        assert_eq!(state.lambda[1], 0.0);
    }
    #[test]
    fn test_dual_state_active_set() {
        let mut state = DualContactState::new(3, 1e6);
        state.lambda[0] = -100.0;
        state.lambda[2] = -50.0;
        let active = state.active_set();
        assert_eq!(active, vec![0, 2]);
    }
    #[test]
    fn test_dual_state_contact_forces_magnitude() {
        let mut state = DualContactState::new(1, 1e6);
        state.lambda[0] = -1000.0;
        let normals = [[0.0, 1.0, 0.0]];
        let forces = state.contact_forces(&normals);
        assert!((forces[0][1] - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_hertz_max_shear_stress_positive() {
        let tau = hertz_max_shear_stress(1e9);
        assert!(tau > 0.0);
        assert!((tau - 0.31e9).abs() < 1.0);
    }
    #[test]
    fn test_hertz_max_shear_depth_scales_with_a() {
        let d = hertz_max_shear_depth(0.005);
        assert!((d - 0.47 * 0.005).abs() < 1e-12);
    }
    #[test]
    fn test_hertz_subsurface_stress_finite() {
        let s = hertz_subsurface_stress(1e9, 0.005, 0.002, 0.001);
        for &v in &s {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn test_mindlin_compliance_positive() {
        let delta_t = mindlin_tangential_compliance(100.0, 1000.0, 0.005, 50e9, 0.3);
        assert!(delta_t > 0.0);
    }
    #[test]
    fn test_rolling_resistance_torque_proportional() {
        let t1 = rolling_resistance_torque(0.01, 1000.0, 0.05);
        let t2 = rolling_resistance_torque(0.01, 2000.0, 0.05);
        assert!((t2 - 2.0 * t1).abs() < 1e-12);
    }
    #[test]
    fn test_rolling_resistance_scales_with_radius() {
        let t1 = rolling_resistance_torque(0.01, 1000.0, 0.05);
        let t2 = rolling_resistance_torque(0.01, 1000.0, 0.10);
        assert!((t2 - 2.0 * t1).abs() < 1e-12);
    }
}
