//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Multiply two 3×3 matrices (row-major).
pub(super) fn mat3_mul(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}
/// Compute the inverse of a 3×3 row-major matrix.
///
/// Returns the zero matrix if singular.
pub(super) fn mat3_inv(m: [f64; 9]) -> [f64; 9] {
    let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    if det.abs() < 1e-15 {
        return [0.0; 9];
    }
    let inv_det = 1.0 / det;
    [
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ]
}
/// Frobenius norm squared of a 3×3 matrix.
pub(super) fn mat3_frobenius_sq(m: [f64; 9]) -> f64 {
    m.iter().map(|x| x * x).sum()
}
/// Determinant of a 3×3 row-major matrix.
pub(super) fn mat3_det(m: [f64; 9]) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}
/// Symmetric product F^T F (returns the symmetric 3×3 right Cauchy-Green tensor).
pub(super) fn mat3_sym_product_fft(f: [f64; 9]) -> [f64; 9] {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = 0.0_f64;
            for k in 0..3 {
                s += f[k * 3 + i] * f[k * 3 + j];
            }
            c[i * 3 + j] = s;
        }
    }
    c
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::surgery_simulation::BleedingModel;
    use crate::surgery_simulation::BloodVessel;
    use crate::surgery_simulation::CollisionResponseSurgical;
    use crate::surgery_simulation::ContactParams;
    use crate::surgery_simulation::CuttingModel;
    use crate::surgery_simulation::CuttingSimulation;
    use crate::surgery_simulation::EdgeStatus;
    use crate::surgery_simulation::HapticFeedback;
    use crate::surgery_simulation::HealingState;
    use crate::surgery_simulation::MeshEdge;
    use crate::surgery_simulation::NeoHookeanParams;
    use crate::surgery_simulation::SimulationStep;
    use crate::surgery_simulation::SoftTissueModel;
    use crate::surgery_simulation::SurgicalTool;
    use crate::surgery_simulation::SuturModel;
    use crate::surgery_simulation::SutureMaterial;
    use crate::surgery_simulation::SutureThread;
    use crate::surgery_simulation::TissueEdge;
    use crate::surgery_simulation::TissueLayer;
    use crate::surgery_simulation::TissueModel;
    use crate::surgery_simulation::TissueNode;
    use crate::surgery_simulation::ToolGeometry;
    fn skin_layer() -> TissueLayer {
        TissueLayer::new("skin", 200_000.0, 0.45, 0.002, 1100.0, 500.0)
    }
    fn fat_layer() -> TissueLayer {
        TissueLayer::new("fat", 3_000.0, 0.49, 0.005, 920.0, 50.0)
    }
    fn simple_tissue() -> TissueModel {
        let layers = vec![skin_layer(), fat_layer()];
        let mut model = TissueModel::new(layers, 0.0);
        model.add_node(TissueNode::new_static([0.0, 0.0, 0.0], 0));
        model.add_node(TissueNode::new([1.0, 0.0, 0.0], 0.01, 0));
        model.add_node(TissueNode::new([0.0, 1.0, 0.0], 0.01, 0));
        model.add_node(TissueNode::new([0.0, 0.0, 1.0], 0.01, 0));
        model.add_element([0, 1, 2, 3]);
        model
    }
    #[test]
    fn test_layer_shear_modulus_positive() {
        let layer = skin_layer();
        assert!(layer.shear_modulus() > 0.0);
    }
    #[test]
    fn test_layer_bulk_modulus_positive() {
        let layer = skin_layer();
        assert!(layer.bulk_modulus() > 0.0);
    }
    #[test]
    fn test_layer_lame_lambda_positive() {
        let layer = skin_layer();
        assert!(layer.lame_lambda() > 0.0);
    }
    #[test]
    fn test_static_node_inv_mass_zero() {
        let node = TissueNode::new_static([0.0, 0.0, 0.0], 0);
        assert_eq!(node.inv_mass(), 0.0);
    }
    #[test]
    fn test_dynamic_node_inv_mass() {
        let node = TissueNode::new([1.0, 0.0, 0.0], 2.0, 0);
        assert!((node.inv_mass() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_static_node_does_not_move() {
        let mut model = simple_tissue();
        let init_pos = model.nodes[0].position;
        model.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        assert_eq!(model.nodes[0].position, init_pos);
    }
    #[test]
    fn test_dynamic_node_falls() {
        let mut model = simple_tissue();
        let init_y = model.nodes[1].position[1];
        for _ in 0..60 {
            model.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        }
        assert!(model.nodes[1].position[1] < init_y);
    }
    #[test]
    fn test_deformation_gradient_identity_at_rest() {
        let model = simple_tissue();
        let rest: Vec<[f64; 3]> = model.nodes.iter().map(|n| n.position).collect();
        let f = model.deformation_gradient(0, &rest);
        assert!((f[0] - 1.0).abs() < 1e-9, "F[0,0] = {}", f[0]);
        assert!((f[4] - 1.0).abs() < 1e-9, "F[1,1] = {}", f[4]);
        assert!((f[8] - 1.0).abs() < 1e-9, "F[2,2] = {}", f[8]);
    }
    #[test]
    fn test_elastic_energy_zero_at_rest() {
        let model = simple_tissue();
        let rest: Vec<[f64; 3]> = model.nodes.iter().map(|n| n.position).collect();
        let energy = model.elastic_energy_density(0, &rest);
        assert!(
            energy.abs() < 1e-8,
            "Energy at rest should be ≈ 0, got {energy}"
        );
    }
    #[test]
    fn test_tool_contact_inside() {
        let tool = SurgicalTool::new("scalpel", ToolGeometry::Blade, 0.003, [0.0; 3]);
        let node_inside = [0.001, 0.0, 0.0];
        let depth = tool.contact_depth(node_inside, 0.001);
        assert!(
            depth > 0.0,
            "Node inside tool should have positive depth, got {depth}"
        );
    }
    #[test]
    fn test_tool_contact_outside() {
        let tool = SurgicalTool::new("scalpel", ToolGeometry::Blade, 0.003, [0.0; 3]);
        let node_outside = [0.1, 0.0, 0.0];
        let depth = tool.contact_depth(node_outside, 0.001);
        assert!(
            depth <= 0.0,
            "Node outside tool should have non-positive depth"
        );
    }
    #[test]
    fn test_tool_contact_normal_direction() {
        let tool = SurgicalTool::new("needle", ToolGeometry::Needle, 0.001, [0.0; 3]);
        let node = [1.0, 0.0, 0.0];
        let normal = tool.contact_normal(node);
        assert!(
            normal[0] > 0.9,
            "Normal should point towards +X, got {:?}",
            normal
        );
    }
    #[test]
    fn test_tool_translate() {
        let mut tool = SurgicalTool::new("tool", ToolGeometry::Sphere, 0.005, [0.0; 3]);
        tool.translate([0.01, 0.02, 0.0]);
        assert!((tool.position[0] - 0.01).abs() < 1e-12);
        assert!((tool.position[1] - 0.02).abs() < 1e-12);
    }
    #[test]
    fn test_edge_damage_accumulates() {
        let mut edge = TissueEdge::new(0, 1, 0.01);
        edge.apply_damage(0.4);
        assert!((edge.damage - 0.4).abs() < 1e-10);
        assert!(matches!(edge.status, EdgeStatus::Partial(_)));
    }
    #[test]
    fn test_edge_severed_at_full_damage() {
        let mut edge = TissueEdge::new(0, 1, 0.01);
        edge.apply_damage(1.0);
        assert!(edge.is_severed());
    }
    #[test]
    fn test_edge_damage_clamped() {
        let mut edge = TissueEdge::new(0, 1, 0.01);
        edge.apply_damage(5.0);
        assert!((edge.damage - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cutting_no_initial_cuts() {
        let mut cm = CuttingModel::new(500.0, 0.9);
        cm.add_edge(TissueEdge::new(0, 1, 0.01));
        assert_eq!(cm.num_severed(), 0);
    }
    #[test]
    fn test_cutting_force_severs_edge() {
        let mut cm = CuttingModel::new(1.0, 1.0);
        let e = TissueEdge::new(0, 1, 0.01);
        cm.add_edge(e);
        let nodes = vec![[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]];
        for _ in 0..100 {
            cm.apply_cut([0.0, 0.0, 0.0], 1e6, 0.1, &nodes);
        }
        assert_eq!(
            cm.num_severed(),
            1,
            "Edge should be severed after repeated cutting"
        );
    }
    #[test]
    fn test_suture_zero_force_coincident() {
        let suture = SuturModel::new(0, 1, 0.0, 1000.0, 50.0, 1.0);
        let nodes = vec![
            TissueNode::new([0.0, 0.0, 0.0], 0.01, 0),
            TissueNode::new([0.0, 0.0, 0.0], 0.01, 0),
        ];
        let (_f, tension) = suture.tension_force(&nodes);
        assert!(tension.abs() < 1e-10);
    }
    #[test]
    fn test_suture_tension_positive() {
        let suture = SuturModel::new(0, 1, 0.0, 1000.0, 50.0, 1.0);
        let nodes = vec![
            TissueNode::new([0.0, 0.0, 0.0], 0.01, 0),
            TissueNode::new([0.01, 0.0, 0.0], 0.01, 0),
        ];
        let (_f, tension) = suture.tension_force(&nodes);
        assert!(
            tension > 0.0,
            "Separated nodes should have positive tension, got {tension}"
        );
    }
    #[test]
    fn test_suture_initial_healing_zero() {
        let suture = SuturModel::new(0, 1, 0.0, 1000.0, 50.0, 100.0);
        assert!((suture.healing_fraction()).abs() < 1e-10);
    }
    #[test]
    fn test_suture_healing_increases() {
        let mut suture = SuturModel::new(0, 1, 0.0, 1000.0, 50.0, 10.0);
        suture.advance_healing(5.0);
        let frac = suture.healing_fraction();
        assert!(
            frac > 0.0 && frac < 1.0,
            "Healing fraction should be between 0 and 1, got {frac}"
        );
    }
    #[test]
    fn test_suture_heals_fully() {
        let mut suture = SuturModel::new(0, 1, 0.0, 1000.0, 50.0, 1.0);
        suture.advance_healing(1000.0);
        assert_eq!(suture.healing, HealingState::Healed);
    }
    #[test]
    fn test_haptic_reset() {
        let mut hf = HapticFeedback::zero();
        hf.force = [1.0, 2.0, 3.0];
        hf.reset();
        assert!(hf.force_magnitude() < 1e-12);
    }
    #[test]
    fn test_haptic_force_magnitude() {
        let mut hf = HapticFeedback::zero();
        hf.force = [3.0, 4.0, 0.0];
        assert!((hf.force_magnitude() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_haptic_contact_energy() {
        let mut hf = HapticFeedback::zero();
        hf.add_contact([1.0, 0.0, 0.0], 0.002, 10000.0, [0.0; 3]);
        assert!(hf.contact_energy > 0.0, "Contact energy should be positive");
    }
    #[test]
    fn test_simulation_time_advances() {
        let tissue = simple_tissue();
        let tool = SurgicalTool::new("scalpel", ToolGeometry::Blade, 0.003, [5.0, 5.0, 5.0]);
        let cutting = CuttingModel::new(500.0, 0.9);
        let mut sim = SimulationStep::new(tissue, tool, cutting, 1e5, 0.001);
        sim.advance(1.0 / 60.0, [0.0, -9.81, 0.0]);
        assert!((sim.time - 1.0 / 60.0).abs() < 1e-12);
    }
    #[test]
    fn test_simulation_contact_haptic_force() {
        let tissue = simple_tissue();
        let tool = SurgicalTool::new("probe", ToolGeometry::Sphere, 0.01, [1.005, 0.0, 0.0]);
        let cutting = CuttingModel::new(500.0, 0.0);
        let mut sim = SimulationStep::new(tissue, tool, cutting, 1e5, 0.001);
        sim.advance(1.0 / 60.0, [0.0, 0.0, 0.0]);
        let fmag = sim.haptic.force_magnitude();
        assert!(
            fmag > 0.0,
            "Overlapping tool should produce haptic force, got {fmag}"
        );
    }
    #[test]
    fn test_mat3_mul_identity() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0_f64];
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let result = mat3_mul(a, identity);
        for (r, ai) in result.iter().zip(a.iter()) {
            assert!((r - ai).abs() < 1e-10, "A*I should equal A");
        }
    }
    #[test]
    fn test_mat3_inv_product_is_identity() {
        let a = [2.0, 1.0, 0.0, 0.0, 3.0, 1.0, 1.0, 0.0, 2.0_f64];
        let inv_a = mat3_inv(a);
        let product = mat3_mul(a, inv_a);
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        for (p, e) in product.iter().zip(identity.iter()) {
            assert!(
                (p - e).abs() < 1e-10,
                "A*inv(A) should be identity, got {p} vs {e}"
            );
        }
    }
    fn make_soft_tissue() -> SoftTissueModel {
        let mut model = SoftTissueModel::new(0.02);
        let mat = NeoHookeanParams::from_young(10_000.0, 0.45, 1.0);
        model.add_node([0.0, 0.0, 0.0], 0.0, true);
        model.add_node([1.0, 0.0, 0.0], 0.01, false);
        model.add_node([0.0, 1.0, 0.0], 0.01, false);
        model.add_node([0.0, 0.0, 1.0], 0.01, false);
        model.add_element([0, 1, 2, 3], mat);
        model
    }
    #[test]
    fn test_neo_hookean_params_positive() {
        let p = NeoHookeanParams::from_young(10_000.0, 0.45, 1.0);
        assert!(p.mu > 0.0, "mu should be positive");
        assert!(p.kappa > 0.0, "kappa should be positive");
    }
    #[test]
    fn test_neo_hookean_energy_zero_at_identity() {
        let p = NeoHookeanParams::from_young(10_000.0, 0.45, 1.0);
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let w = p.strain_energy_density(&identity);
        assert!(w.abs() < 1e-8, "Energy at identity should be ≈ 0, got {w}");
    }
    #[test]
    fn test_neo_hookean_energy_positive_under_stretch() {
        let p = NeoHookeanParams::from_young(10_000.0, 0.45, 1.0);
        let f_stretch = [1.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let w = p.strain_energy_density(&f_stretch);
        assert!(w > 0.0, "Energy under stretch should be positive, got {w}");
    }
    #[test]
    fn test_soft_tissue_fixed_node_stationary() {
        let mut model = make_soft_tissue();
        let pos0 = model.nodes[0];
        model.step(0.01, [0.0, -9.81, 0.0]);
        assert_eq!(model.nodes[0], pos0, "fixed node should not move");
    }
    #[test]
    fn test_soft_tissue_free_node_moves() {
        let mut model = make_soft_tissue();
        let y0 = model.nodes[1][1];
        model.step(0.1, [0.0, -9.81, 0.0]);
        assert!(
            model.nodes[1][1] < y0,
            "free node should fall under gravity"
        );
    }
    #[test]
    fn test_soft_tissue_element_volume_positive() {
        let model = make_soft_tissue();
        let vol = model.element_volume(0);
        assert!(vol > 0.0, "element volume should be positive, got {vol}");
    }
    #[test]
    fn test_soft_tissue_total_energy_finite() {
        let model = make_soft_tissue();
        let energy = model.total_energy();
        assert!(
            energy.is_finite(),
            "total energy should be finite, got {energy}"
        );
    }
    #[test]
    fn test_cutting_sim_no_initial_severed() {
        let mut cs = CuttingSimulation::new(500.0, 0.9);
        cs.add_edge(MeshEdge::new(0, 1, 1.0));
        assert_eq!(cs.n_severed(), 0);
    }
    #[test]
    fn test_cutting_sim_force_severs_edge() {
        let mut cs = CuttingSimulation::new(1.0, 1.0);
        cs.add_edge(MeshEdge::new(0, 1, 1.0));
        let nodes = vec![[0.0_f64; 3], [0.01, 0.0, 0.0]];
        for _ in 0..200 {
            cs.apply_cut([0.005, 0.0, 0.0], 1e5_f64, 0.05, &nodes);
        }
        assert!(cs.n_severed() >= 1, "edge should be severed after cutting");
    }
    #[test]
    fn test_cutting_sim_intact_count_decreases() {
        let mut cs = CuttingSimulation::new(1.0, 1.0);
        cs.add_edge(MeshEdge::new(0, 1, 1.0));
        cs.add_edge(MeshEdge::new(1, 2, 1.0));
        let nodes = vec![[0.0_f64; 3], [0.005, 0.0, 0.0], [0.01, 0.0, 0.0]];
        for _ in 0..500 {
            cs.apply_cut([0.005, 0.0, 0.0], 1e6_f64, 0.05, &nodes);
        }
        assert!(cs.n_intact() < 2, "intact count should decrease");
    }
    #[test]
    fn test_cutting_sim_crack_path_nonempty() {
        let mut cs = CuttingSimulation::new(1.0, 1.0);
        cs.add_edge(MeshEdge::new(0, 1, 0.5));
        let nodes = vec![[0.0_f64; 3], [0.01, 0.0, 0.0]];
        cs.apply_cut([0.005, 0.0, 0.0], 1e6_f64, 0.05, &nodes);
        cs.apply_cut([0.005, 0.0, 0.0], 1e6_f64, 0.05, &nodes);
        if cs.n_severed() > 0 {
            let path = cs.crack_path(&nodes);
            assert!(!path.is_empty(), "crack path should be non-empty");
        }
    }
    #[test]
    fn test_suture_material_axial_stiffness() {
        let mat = SutureMaterial::new(1e8_f64, 1e-6_f64, 5e7_f64, 0.9, 0.1);
        let expected = 1e8_f64 * 1e-6_f64;
        assert!(
            (mat.axial_stiffness() - expected).abs() < 1e-4,
            "axial stiffness mismatch"
        );
    }
    #[test]
    fn test_suture_material_max_force_reduced_by_knot() {
        let mat_perfect = SutureMaterial::new(1e8_f64, 1e-6_f64, 5e7_f64, 1.0, 0.0);
        let mat_knotted = SutureMaterial::new(1e8_f64, 1e-6_f64, 5e7_f64, 0.7, 0.0);
        assert!(
            mat_knotted.max_force() < mat_perfect.max_force(),
            "knotted suture should be weaker"
        );
    }
    #[test]
    fn test_suture_material_pretension_included() {
        let mat = SutureMaterial::new(1e8_f64, 1e-6_f64, 5e7_f64, 1.0, 0.5);
        let f = mat.elastic_force(0.0);
        assert!(
            (f - 0.5).abs() < 1e-10,
            "force at zero elongation should equal pretension, got {f}"
        );
    }
    #[test]
    fn test_suture_thread_cut_gives_zero_force() {
        let mat = SutureMaterial::new(1e6_f64, 1e-6_f64, 1e6_f64, 1.0, 0.0);
        let mut thread = SutureThread::new(0, 1, 0.0, mat);
        thread.cut();
        let nodes = vec![[0.0_f64; 3], [0.01, 0.0, 0.0]];
        let (_f, tension) = thread.force(&nodes);
        assert!(tension.abs() < 1e-12, "cut thread should have zero tension");
    }
    #[test]
    fn test_suture_thread_heals() {
        let mat = SutureMaterial::new(1e6_f64, 1e-6_f64, 1e6_f64, 1.0, 0.0);
        let mut thread = SutureThread::new(0, 1, 0.0, mat);
        thread.heal(100.0, 10.0);
        assert!(
            thread.is_healed(),
            "thread should be healed after sufficient time"
        );
    }
    #[test]
    fn test_vessel_no_bleed_intact() {
        let vessel = BloodVessel::new([0.0; 3], 0.001, 1e6_f64, 1200.0);
        assert!(
            (vessel.bleeding_rate(3e-3_f64)).abs() < 1e-20,
            "intact vessel should not bleed"
        );
    }
    #[test]
    fn test_vessel_bleeds_after_rupture() {
        let mut vessel = BloodVessel::new([0.0; 3], 0.001, 1e4_f64, 1200.0);
        vessel.check_rupture(2e4_f64);
        let rate = vessel.bleeding_rate(3e-3_f64);
        assert!(rate > 0.0, "ruptured vessel should bleed, rate={rate}");
    }
    #[test]
    fn test_bleeding_model_trauma_ruptures_vessel() {
        let vessel = BloodVessel::new([0.0; 3], 0.001, 1.0, 1200.0);
        let mut bm = BleedingModel::new(vec![vessel], 3e-3_f64);
        let new_rupt = bm.apply_trauma([0.0; 3], 1e8_f64, 0.05);
        assert!(
            new_rupt > 0 || bm.n_ruptured() > 0,
            "trauma should rupture vessel"
        );
    }
    #[test]
    fn test_collision_response_normal_force_positive() {
        let params = ContactParams::new(1e5_f64, 10.0, 0.3);
        let mut cr = CollisionResponseSurgical::new(params);
        let (force_node, _force_tool) = cr.compute_contact(0.002, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(force_node[1] > 0.0, "contact force should be along normal");
    }
}
