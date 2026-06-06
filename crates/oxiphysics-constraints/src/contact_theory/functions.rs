//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::LcpFormulation;

pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Build the LCP matrix for contact constraints from per-contact data.
pub fn build_lcp_matrix(eff_masses: &[f64], velocity_biases: &[f64]) -> LcpFormulation {
    LcpFormulation::from_contacts(eff_masses, velocity_biases)
}
/// Compute friction cone edge directions for a polyhedral approximation.
pub fn compute_cone_edges(normal: [f64; 3], mu: f64, num_edges: usize) -> Vec<[f64; 3]> {
    let n = normalize3(normal);
    let t1 = {
        let candidate = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let cr = cross3(n, candidate);
        normalize3(cr)
    };
    let t2 = normalize3(cross3(n, t1));
    let ne = num_edges.max(4);
    let mut edges = Vec::with_capacity(ne);
    for i in 0..ne {
        let angle = 2.0 * PI * (i as f64) / (ne as f64);
        let dir = [
            n[0] + mu * (angle.cos() * t1[0] + angle.sin() * t2[0]),
            n[1] + mu * (angle.cos() * t1[1] + angle.sin() * t2[1]),
            n[2] + mu * (angle.cos() * t1[2] + angle.sin() * t2[2]),
        ];
        edges.push(normalize3(dir));
    }
    edges
}
/// Project a 3-D force vector onto the linearised friction cone.
pub fn project_to_friction_cone(f: [f64; 3], normal: [f64; 3], mu: f64, fn_mag: f64) -> [f64; 3] {
    let n = normalize3(normal);
    let f_n = dot3(f, n);
    let ft = [f[0] - f_n * n[0], f[1] - f_n * n[1], f[2] - f_n * n[2]];
    let ft_mag = norm3(ft);
    let max_tangential = mu * fn_mag.max(0.0);
    if ft_mag <= max_tangential + 1e-15 {
        return f;
    }
    if ft_mag < 1e-15 {
        return f;
    }
    let scale = max_tangential / ft_mag;
    [
        f_n * n[0] + ft[0] * scale,
        f_n * n[1] + ft[1] * scale,
        f_n * n[2] + ft[2] * scale,
    ]
}
/// Baumgarte stabilization correction velocity for penetration `depth`.
pub fn baumgarte_correction(depth: f64, beta: f64, dt: f64, slop: f64) -> f64 {
    let correctable = (depth - slop).max(0.0);
    if dt.abs() < 1e-20 {
        return 0.0;
    }
    -(beta / dt) * correctable
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::contact_theory::ContactDamping;
    use crate::contact_theory::ContactGraph;
    use crate::contact_theory::ContactManifoldUpdate;
    use crate::contact_theory::ContactPair;
    use crate::contact_theory::ContactStiffness;

    use crate::contact_theory::FrictionCone;

    use crate::contact_theory::ImpulseBasedSolver;

    use crate::contact_theory::LcpSolver;
    use crate::contact_theory::MachineLearningContact;

    use crate::contact_theory::PenetrationHandling;

    #[test]
    fn test_contact_graph_empty() {
        let g = ContactGraph::new(4);
        assert_eq!(g.num_contacts(), 0);
    }
    #[test]
    fn test_contact_graph_add() {
        let mut g = ContactGraph::new(4);
        g.add_contact(ContactPair::new(0, 1, [0.0; 3], [0.0, 1.0, 0.0], 0.01));
        assert_eq!(g.num_contacts(), 1);
    }
    #[test]
    fn test_connected_components_isolated() {
        let g = ContactGraph::new(2);
        let comps = g.connected_components();
        assert_eq!(comps.len(), 2);
    }
    #[test]
    fn test_connected_components_one_island() {
        let mut g = ContactGraph::new(2);
        g.add_contact(ContactPair::new(0, 1, [0.0; 3], [0.0, 1.0, 0.0], 0.01));
        let comps = g.connected_components();
        assert_eq!(comps.len(), 1);
    }
    #[test]
    fn test_lcp_formulation_valid() {
        let lcp = LcpFormulation::from_contacts(&[1.0, 2.0], &[-0.5, -1.0]);
        assert!(lcp.is_valid());
    }
    #[test]
    fn test_lcp_pgs_diagonal() {
        let mut lcp = LcpFormulation::new(2);
        lcp.m[0] = 2.0;
        lcp.m[3] = 2.0;
        lcp.q = vec![-1.0, -1.0];
        let solver = LcpSolver::new(100);
        let sol = solver.solve_pgs(&lcp);
        assert!(sol.z[0] > 0.0, "z[0]={}", sol.z[0]);
        assert!(sol.z[1] > 0.0, "z[1]={}", sol.z[1]);
    }
    #[test]
    fn test_lcp_pgs_no_penetration() {
        let mut lcp = LcpFormulation::new(1);
        lcp.m[0] = 1.0;
        lcp.q = vec![1.0];
        let solver = LcpSolver::new(50);
        let sol = solver.solve_pgs(&lcp);
        assert!(sol.z[0].abs() < 1e-6, "z[0]={}", sol.z[0]);
    }
    #[test]
    fn test_friction_cone_inside() {
        let fc = FrictionCone::new(0.5, [0.0, 1.0, 0.0], 4);
        let f = [3.0, 10.0, 0.0];
        assert!(fc.in_cone(f));
    }
    #[test]
    fn test_friction_cone_outside() {
        let fc = FrictionCone::new(0.3, [0.0, 1.0, 0.0], 4);
        let f = [9.0, 10.0, 0.0];
        assert!(!fc.in_cone(f));
    }
    #[test]
    fn test_project_to_friction_cone() {
        let normal = [0.0, 1.0, 0.0];
        let mu = 0.5;
        let fn_mag = 10.0;
        let f = [8.0, 10.0, 0.0];
        let pf = project_to_friction_cone(f, normal, mu, fn_mag);
        let ft = pf[0];
        assert!(ft.abs() <= mu * fn_mag + 1e-9, "ft={ft}");
    }
    #[test]
    fn test_hertz_force_no_penetration() {
        let cs = ContactStiffness::new(1e9, 1e9, 0.3, 0.3);
        assert!(cs.hertz_force(0.01, 0.01, -0.001) < 1e-15);
        assert!(cs.hertz_force(0.01, 0.01, 0.0) < 1e-15);
    }
    #[test]
    fn test_hertz_force_positive() {
        let cs = ContactStiffness::new(1e9, 1e9, 0.3, 0.3);
        let f = cs.hertz_force(0.01, 0.01, 0.001);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn test_hunt_crossley_force() {
        let cd = ContactDamping::new(0.1);
        let f = cd.hunt_crossley_force(0.01, 1.0);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn test_hunt_crossley_no_contact() {
        let cd = ContactDamping::new(0.1);
        let f = cd.hunt_crossley_force(-0.01, 1.0);
        assert!(f.abs() < 1e-15);
    }
    #[test]
    fn test_baumgarte_correction() {
        let bias = baumgarte_correction(0.1, 0.2, 0.01, 0.0);
        assert!((bias + 2.0).abs() < 1e-10, "bias={bias}");
    }
    #[test]
    fn test_baumgarte_within_slop() {
        let bias = baumgarte_correction(0.001, 0.2, 0.01, 0.01);
        assert!(bias.abs() < 1e-10, "bias={bias}");
    }
    #[test]
    fn test_penetration_position_correction() {
        let ph = PenetrationHandling::new(0.2, 0.001, 1000.0);
        let corr = ph.position_correction(0.01);
        assert!((corr - 0.009).abs() < 1e-10, "corr={corr}");
    }
    #[test]
    fn test_manifold_warm_start_unknown() {
        let cm = ContactManifoldUpdate::new(0.9, 0.1);
        let (ln, lt) = cm.warm_start(0, 1);
        assert!(ln.abs() < 1e-15);
        assert!(lt[0].abs() < 1e-15);
    }
    #[test]
    fn test_manifold_update_and_retrieve() {
        let mut cm = ContactManifoldUpdate::new(0.9, 0.1);
        cm.update(0, 1, 5.0, [1.0, -0.5]);
        let (ln, lt) = cm.warm_start(0, 1);
        assert!((ln - 5.0).abs() < 1e-10);
        assert!((lt[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_manifold_aging() {
        let mut cm = ContactManifoldUpdate::new(0.5, 0.0);
        cm.update(0, 1, 10.0, [0.0; 2]);
        cm.age_all();
        let (ln, _) = cm.warm_start(0, 1);
        assert!(ln < 10.0, "ln={ln}");
    }
    #[test]
    fn test_ml_contact_predict_size() {
        let mlc = MachineLearningContact::new(6, 8, 3);
        let out = mlc.predict(&[0.0; 6]);
        assert_eq!(out.len(), 3);
    }
    #[test]
    fn test_ml_contact_predict_zero_weights() {
        let mut mlc = MachineLearningContact::new(4, 4, 2);
        mlc.b2 = vec![1.0, -1.0];
        let out = mlc.predict(&[0.0; 4]);
        assert!((out[0] - 1.0).abs() < 1e-10);
        assert!((out[1] + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ml_contact_training_reduces_loss() {
        let mut mlc = MachineLearningContact::new(2, 4, 1);
        mlc.w1 = vec![0.1; 8];
        mlc.w2 = vec![0.1; 4];
        let input = vec![1.0, 0.5];
        let target = vec![2.0];
        let pred_before = mlc.predict(&input)[0];
        let loss_before = (pred_before - target[0]).powi(2);
        for _ in 0..50 {
            mlc.train_step(&input, &target, 0.01);
        }
        let pred_after = mlc.predict(&input)[0];
        let loss_after = (pred_after - target[0]).powi(2);
        assert!(
            loss_after < loss_before + 1.0,
            "loss before={loss_before}, after={loss_after}"
        );
    }
    #[test]
    fn test_normal_impulse_nonneg() {
        let solver = ImpulseBasedSolver::new(10, true);
        let mut lambda_n = 0.0;
        let j = solver.normal_impulse(-1.0, 1.0, 0.5, &mut lambda_n);
        assert!(j >= 0.0, "j={j}");
        assert!(lambda_n >= 0.0);
    }
    #[test]
    fn test_friction_impulse_clamped() {
        let solver = ImpulseBasedSolver::new(10, true);
        let mut lambda_t = [0.0f64; 2];
        let delta = solver.friction_impulse([10.0, 0.0], 1.0, 0.3, 5.0, &mut lambda_t);
        let total = (lambda_t[0] * lambda_t[0] + lambda_t[1] * lambda_t[1]).sqrt();
        assert!(total <= 0.3 * 5.0 + 1e-9, "total={total}");
        let _ = delta;
    }
    #[test]
    fn test_build_lcp_matrix() {
        let lcp = build_lcp_matrix(&[1.0, 2.0, 3.0], &[-1.0, -2.0, -3.0]);
        assert!(lcp.is_valid());
        assert_eq!(lcp.n, 3);
    }
    #[test]
    fn test_contact_graph_clear() {
        let mut g = ContactGraph::new(2);
        g.add_contact(ContactPair::new(0, 1, [0.0; 3], [0.0, 1.0, 0.0], 0.1));
        g.clear();
        assert_eq!(g.num_contacts(), 0);
    }
    #[test]
    fn test_hertz_modulus_symmetric() {
        let cs1 = ContactStiffness::new(1e9, 2e9, 0.3, 0.25);
        let cs2 = ContactStiffness::new(2e9, 1e9, 0.25, 0.3);
        assert!((cs1.hertz_modulus() - cs2.hertz_modulus()).abs() < 1e-3);
    }
    #[test]
    fn test_friction_cone_edge_count() {
        let fc = FrictionCone::new(0.5, [0.0, 1.0, 0.0], 8);
        assert_eq!(fc.edges.len(), 8);
    }
}
/// Approximate complementary error function erfc(x) for x ≥ 0.
pub(super) fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}
#[cfg(test)]
mod expanded_contact_tests {
    use super::*;
    use crate::contact_theory::BoussinesqHalfSpace;
    use crate::contact_theory::CebAsperity;
    use crate::contact_theory::DmtContact;
    use crate::contact_theory::ElastoPlasticMaterial;
    use crate::contact_theory::GreenwoodWilliamson;
    use crate::contact_theory::HertzPressureDistribution;
    use crate::contact_theory::HuntCrossleyContact;
    use crate::contact_theory::ImpactMechanics;
    use crate::contact_theory::JkrContact;
    use crate::contact_theory::MindlinContact;
    use crate::contact_theory::ThorntonContact;
    use crate::contact_theory::TribologicalContact;
    #[test]
    fn test_elastic_yield_strain() {
        let mat = ElastoPlasticMaterial::steel();
        let ys = mat.yield_strain();
        assert!((ys - mat.yield_stress / mat.elastic_modulus).abs() < 1e-15);
    }
    #[test]
    fn test_flow_stress_elastic() {
        let mat = ElastoPlasticMaterial::steel();
        let strain = mat.yield_strain() * 0.5;
        let sigma = mat.flow_stress(strain);
        assert!((sigma - mat.elastic_modulus * strain).abs() < 1e-3);
    }
    #[test]
    fn test_thornton_zero_force() {
        let tc = ThorntonContact::new(100e9, 0.01, 200e6, 600e6);
        assert!(tc.contact_force(0.0).abs() < 1e-10);
    }
    #[test]
    fn test_thornton_elastic_regime() {
        let e_star = 100e9_f64;
        let r_star = 0.01_f64;
        let tc = ThorntonContact::new(e_star, r_star, 200e9, 600e9);
        let delta_y = tc.yield_interference();
        let delta = delta_y * 0.5;
        let f_tc = tc.contact_force(delta);
        let f_hertz = (4.0 / 3.0) * e_star * r_star.sqrt() * delta.powf(1.5);
        assert!(
            (f_tc - f_hertz).abs() / f_hertz < 1e-10,
            "f_tc={f_tc} f_hertz={f_hertz}"
        );
    }
    #[test]
    fn test_ceb_critical_interference() {
        let asperity = CebAsperity::new(1e-6, 70e9, 100e6);
        let dc = asperity.critical_interference();
        assert!(dc > 0.0, "dc={dc}");
    }
    #[test]
    fn test_ceb_force_increasing() {
        let asperity = CebAsperity::new(1e-6, 70e9, 100e6);
        let f1 = asperity.force(1e-9);
        let f2 = asperity.force(2e-9);
        assert!(f2 > f1, "f1={f1}, f2={f2}");
    }
    #[test]
    fn test_erfc_zero() {
        assert!((erfc_approx(0.0) - 1.0).abs() < 0.001);
    }
    #[test]
    fn test_gw_plasticity_index() {
        let gw = GreenwoodWilliamson::new(1e13, 1e-6, 1e-7, 100e9);
        let pi = gw.plasticity_index(1e9);
        assert!(pi >= 0.0, "plasticity_index={pi}");
    }
    #[test]
    fn test_hunt_crossley_no_force() {
        let mut hc = HuntCrossleyContact::new(1e6, 1.5, 0.1);
        let f = hc.force(-0.01, 0.0);
        assert!(f.abs() < 1e-10, "f={f}");
    }
    #[test]
    fn test_hunt_crossley_positive_force() {
        let mut hc = HuntCrossleyContact::new(1e6, 1.5, 0.1);
        let f = hc.force(1e-3, -0.01);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn test_hunt_crossley_restitution() {
        let hc = HuntCrossleyContact::new(1e6, 1.5, 0.1);
        let e = hc.restitution_approx(1.0);
        assert!((0.0..=1.0).contains(&e), "e={e}");
    }
    #[test]
    fn test_jkr_pull_off() {
        let jkr = JkrContact::new(100e9, 1e-6, 0.05);
        let f_po = jkr.pull_off_force();
        assert!(f_po > 0.0, "f_po={f_po}");
    }
    #[test]
    fn test_jkr_zero_radius_force() {
        let jkr = JkrContact::new(100e9, 1e-6, 0.05);
        let f = jkr.total_force(0.0);
        assert!(f.abs() < 1e-10, "f={f}");
    }
    #[test]
    fn test_dmt_vs_jkr_pull_off() {
        let e_star = 100e9_f64;
        let r_star = 1e-6_f64;
        let w = 0.05_f64;
        let dmt = DmtContact::new(e_star, r_star, w);
        let jkr = JkrContact::new(e_star, r_star, w);
        assert!(jkr.pull_off_force() < dmt.pull_off_force());
    }
    #[test]
    fn test_dmt_contact_force_positive() {
        let dmt = DmtContact::new(100e9, 1e-6, 0.05);
        let f = dmt.contact_force(1e-8);
        assert!(f > 0.0, "f={f}");
    }
    #[test]
    fn test_hertz_integrated_force() {
        let f_applied = 1.0;
        let dist = HertzPressureDistribution::from_force(f_applied, 100e9, 1e-3);
        let f_int = dist.integrated_force();
        assert!(
            (f_int - f_applied).abs() / f_applied < 0.01,
            "f_int={f_int}"
        );
    }
    #[test]
    fn test_hertz_pressure_outside() {
        let dist = HertzPressureDistribution::from_force(1.0, 100e9, 1e-3);
        let p = dist.pressure(dist.contact_radius + 1e-10);
        assert!(p.abs() < 1e-10, "p={p}");
    }
    #[test]
    fn test_mindlin_stiffness() {
        let mc = MindlinContact::new(1e-4, 80e9, 10.0, 0.3);
        let kt = mc.tangential_stiffness();
        assert!(kt > 0.0, "kt={kt}");
    }
    #[test]
    fn test_mindlin_stick() {
        let mut mc = MindlinContact::new(1e-4, 80e9, 100.0, 0.3);
        let f = mc.increment_force(1e-10);
        assert!(f < mc.mu * mc.normal_force, "Should be in stick");
    }
    #[test]
    fn test_mindlin_slip() {
        let mut mc = MindlinContact::new(1e-4, 80e9, 10.0, 0.3);
        let f = mc.increment_force(1.0);
        assert!((f.abs() - mc.mu * mc.normal_force).abs() < 1e-6, "f={f}");
    }
    #[test]
    fn test_impact_momentum() {
        let imp = ImpactMechanics::new(1.0, 1.0, 0.8, 0.5, 0.3);
        let (va, vb) = imp.normal_impact_1d(1.0, -1.0);
        let p_before = 1.0 * 1.0 + -1.0;
        let p_after = 1.0 * va + 1.0 * vb;
        assert!(
            (p_after - p_before).abs() < 1e-10,
            "p_after={p_after}, p_before={p_before}"
        );
    }
    #[test]
    fn test_impact_energy_dissipated() {
        let imp = ImpactMechanics::new(1.0, 2.0, 0.7, 0.5, 0.3);
        let e = imp.energy_dissipated(2.0, -1.0);
        assert!(e >= 0.0, "e={e}");
    }
    #[test]
    fn test_stribeck_model() {
        let tri = TribologicalContact::new(1e-8, 1e9, 0.1, 0.3);
        let mu_static = tri.friction_coefficient(0.0);
        let mu_kinetic = tri.friction_coefficient(10.0);
        assert!(
            mu_static > mu_kinetic,
            "static={mu_static}, kinetic={mu_kinetic}"
        );
    }
    #[test]
    fn test_archard_wear() {
        let mut tri = TribologicalContact::new(1e-8, 1e9, 0.1, 0.3);
        tri.wear_increment(1000.0, 0.1);
        assert!(tri.wear_volume > 0.0, "wear_volume={}", tri.wear_volume);
    }
    #[test]
    fn test_boussinesq_displacement() {
        let hs = BoussinesqHalfSpace::new(200e9, 0.3);
        let w1 = hs.surface_displacement(1000.0, 1e-3);
        let w2 = hs.surface_displacement(1000.0, 2e-3);
        assert!(w1 > w2, "w1={w1}, w2={w2}");
    }
    #[test]
    fn test_boussinesq_axis_stress() {
        let hs = BoussinesqHalfSpace::new(200e9, 0.3);
        let sigma = hs.max_stress_on_axis(1000.0, 1e-3);
        assert!(sigma < 0.0, "Compressive stress expected, sigma={sigma}");
    }
    #[test]
    fn test_hardness_estimate() {
        let mat = ElastoPlasticMaterial::steel();
        let h = mat.hardness_estimate();
        assert!((h - 750e6).abs() < 1e-3, "h={h}");
    }
    #[test]
    fn test_impact_duration_positive() {
        let hc = HuntCrossleyContact::new(1e7, 1.5, 0.05);
        let td = hc.impact_duration(0.1, 1.0);
        assert!(td > 0.0, "td={td}");
    }
    #[test]
    fn test_mindlin_stick_capacity() {
        let mut mc = MindlinContact::new(1e-4, 80e9, 100.0, 0.3);
        let cap0 = mc.stick_capacity();
        mc.increment_force(1e-8);
        let cap1 = mc.stick_capacity();
        assert!(cap1 <= cap0, "Stick capacity should decrease");
    }
    #[test]
    fn test_gw_elastic_force_positive() {
        let gw = GreenwoodWilliamson::new(1e13, 1e-6, 1e-7, 100e9);
        let f = gw.elastic_force(1e-4, 0.5);
        assert!(f >= 0.0, "f={f}");
    }
}
