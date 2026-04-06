//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

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
pub(super) fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::cloth_advanced::ClothEdge;
    use crate::cloth_advanced::ClothGarment;
    use crate::cloth_advanced::ClothInteraction;
    use crate::cloth_advanced::ClothParticle;
    use crate::cloth_advanced::ClothRenderMesh;
    use crate::cloth_advanced::ClothSewing;
    use crate::cloth_advanced::ClothTearing;
    use crate::cloth_advanced::ClothWetting;

    use crate::cloth_advanced::ElasticFabric;

    use crate::cloth_advanced::FabricWeave;

    use crate::cloth_advanced::InteractionMode;

    use crate::cloth_advanced::LayeredCloth;
    use crate::cloth_advanced::SeamLine;

    use crate::cloth_advanced::WindDraping;

    #[test]
    fn test_layered_cloth_creation() {
        let lc = LayeredCloth::new(3, 10, 0.01, 0.1, 0.002);
        assert_eq!(lc.num_layers(), 3);
        assert_eq!(lc.layers[0].len(), 10);
    }
    #[test]
    fn test_layered_cloth_initial_ke_zero() {
        let lc = LayeredCloth::new(2, 4, 0.01, 0.1, 0.001);
        assert!(lc.kinetic_energy().abs() < 1e-10);
    }
    #[test]
    fn test_layered_cloth_friction_no_crash() {
        let mut lc = LayeredCloth::new(2, 4, 0.01, 0.1, 0.001);
        lc.layers[0][0].velocity = [1.0, 0.0, 0.0];
        lc.apply_inter_layer_friction(0.01);
    }
    #[test]
    fn test_fabric_weave_warp_stiffness() {
        let fw = FabricWeave::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1000.0, 500.0, 200.0);
        let k = fw.stiffness_in_direction([1.0, 0.0, 0.0]);
        assert!((k - 1000.0).abs() < 1.0, "k={k}");
    }
    #[test]
    fn test_fabric_weave_weft_stiffness() {
        let fw = FabricWeave::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1000.0, 500.0, 200.0);
        let k = fw.stiffness_in_direction([0.0, 1.0, 0.0]);
        assert!((k - 500.0).abs() < 1.0, "k={k}");
    }
    #[test]
    fn test_fabric_weave_edge_force_linear() {
        let fw = FabricWeave::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1000.0, 500.0, 200.0);
        let f1 = fw.edge_force([1.0, 0.0, 0.0], 0.01);
        let f2 = fw.edge_force([1.0, 0.0, 0.0], 0.02);
        assert!((f2 - 2.0 * f1).abs() < 1e-6, "f1={f1}, f2={f2}");
    }
    #[test]
    fn test_cloth_sewing_no_crash() {
        let mut cs = ClothSewing::new(2, 3, 0.01);
        let seam = SeamLine::new(vec![0, 1, 2], vec![3, 4, 5], 10.0);
        cs.add_seam(seam);
        cs.apply_seam_constraints(0.01);
    }
    #[test]
    fn test_garment_collision_resolution() {
        let mut garment = ClothGarment::new("test", 1, 1, 0.01);
        garment.sewing.particles[0].position = [0.2, 0.0, 0.0];
        garment.add_avatar_sphere([0.0, 0.0, 0.0], 0.5);
        garment.resolve_avatar_collisions();
        let p = garment.sewing.particles[0].position;
        let dist = norm3(p);
        assert!(dist >= 0.5 - 1e-6, "dist={dist}");
    }
    #[test]
    fn test_wind_draping_zero_wind() {
        let wd = WindDraping::new([0.0; 3], 1.225, 1.0);
        let f = wd.triangle_force([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        for fi in &f {
            for c in fi {
                assert!(c.abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_wind_draping_normal_force() {
        let wd = WindDraping::new([0.0, 10.0, 0.0], 1.225, 1.0);
        let f = wd.triangle_force([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        let fy: f64 = f.iter().map(|fi| fi[1]).sum();
        assert!(fy.abs() > 1e-10, "fy={fy}");
    }
    #[test]
    fn test_cloth_tearing_initially_intact() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([0.1, 0.0, 0.0], 0.01),
        ];
        let edges = vec![ClothEdge::new(0, 1, 0.1, 0.5)];
        let ct = ClothTearing::new(particles, edges);
        assert_eq!(ct.intact_edge_count(), 1);
        assert_eq!(ct.torn_edge_count(), 0);
    }
    #[test]
    fn test_cloth_tearing_rupture() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([0.2, 0.0, 0.0], 0.01),
        ];
        let edges = vec![ClothEdge::new(0, 1, 0.1, 0.5)];
        let mut ct = ClothTearing::new(particles, edges);
        ct.check_all_ruptures();
        assert_eq!(ct.torn_edge_count(), 1);
    }
    #[test]
    fn test_cloth_wetting_initially_dry() {
        let cw = ClothWetting::new(4, 0.3, 0.5);
        assert!(cw.total_water_mass() < 1e-10);
    }
    #[test]
    fn test_cloth_wetting_accumulates() {
        let mut cw = ClothWetting::new(4, 0.3, 0.5);
        cw.wet_particle(0, 0.1);
        assert!((cw.water_content[0] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_cloth_wetting_evaporation() {
        let mut cw = ClothWetting::new(2, 0.3, 0.5);
        cw.wet_particle(0, 0.3);
        let before = cw.water_content[0];
        cw.update_evaporation(1.0);
        assert!(cw.water_content[0] < before);
    }
    #[test]
    fn test_elastic_fabric_zero_strain() {
        let ef = ElasticFabric::new(1000.0, 500.0, 0.3, 200.0, 0.01);
        let f = ef.in_plane_force([0.0, 0.0]);
        assert!(f[0].abs() < 1e-10 && f[1].abs() < 1e-10);
    }
    #[test]
    fn test_elastic_fabric_tensor_symmetry() {
        let ef = ElasticFabric::new(1000.0, 500.0, 0.3, 200.0, 0.01);
        let c = ef.stiffness_tensor();
        assert!((c[0][1] - c[1][0]).abs() < 1e-10);
    }
    #[test]
    fn test_cloth_render_mesh_vertex_count() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([1.0, 0.0, 0.0], 0.01),
            ClothParticle::new([0.0, 1.0, 0.0], 0.01),
        ];
        let tris = vec![[0usize, 1, 2]];
        let rm = ClothRenderMesh::from_particles(&particles, &tris);
        assert_eq!(rm.num_vertices(), 3);
        assert_eq!(rm.num_triangles(), 1);
    }
    #[test]
    fn test_cloth_render_mesh_normals_unit() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([1.0, 0.0, 0.0], 0.01),
            ClothParticle::new([0.0, 1.0, 0.0], 0.01),
        ];
        let tris = vec![[0usize, 1, 2]];
        let rm = ClothRenderMesh::from_particles(&particles, &tris);
        for v in &rm.vertices {
            let len =
                (v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2])
                    .sqrt();
            assert!((len - 1.0).abs() < 1e-6, "len={len}");
        }
    }
    #[test]
    fn test_cloth_interaction_find_closest() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([1.0, 0.0, 0.0], 0.01),
            ClothParticle::new([2.0, 0.0, 0.0], 0.01),
        ];
        let ci = ClothInteraction::new(100.0);
        let idx = ci.find_closest_particle(&particles, [0.9, 0.0, 0.0]);
        assert_eq!(idx, Some(1));
    }
    #[test]
    fn test_cloth_interaction_grab_moves() {
        let mut particles = vec![ClothParticle::new([0.0, 0.0, 0.0], 0.01)];
        let mut ci = ClothInteraction::new(100.0);
        ci.target_position = [1.0, 0.0, 0.0];
        ci.select(0, InteractionMode::Pull);
        ci.apply(&mut particles, 0.01);
        assert!(
            particles[0].velocity[0] > 0.0,
            "vx={}",
            particles[0].velocity[0]
        );
    }
    #[test]
    fn test_cloth_interaction_cut() {
        let mut particles = vec![ClothParticle::new([0.0, 0.0, 0.0], 0.01)];
        let mut ci = ClothInteraction::new(100.0);
        ci.select(0, InteractionMode::Cut);
        ci.apply(&mut particles, 0.01);
        assert!(particles[0].is_pinned());
    }
    #[test]
    fn test_cloth_wetting_max_capacity() {
        let mut cw = ClothWetting::new(1, 0.3, 0.5);
        cw.wet_particle(0, 10.0);
        assert!(cw.water_content[0] <= 0.5 + 1e-10);
    }
    #[test]
    fn test_elastic_fabric_sed_positive() {
        let ef = ElasticFabric::new(1000.0, 500.0, 0.3, 200.0, 0.01);
        let sed = ef.strain_energy_density([0.01, 0.005]);
        assert!(sed > 0.0, "sed={sed}");
    }
    #[test]
    fn test_seam_line_num_pairs() {
        let seam = SeamLine::new(vec![0, 1, 2], vec![3, 4], 10.0);
        assert_eq!(seam.num_pairs(), 2);
    }
    #[test]
    fn test_cloth_particle_pinned() {
        let p = ClothParticle::pinned([0.0, 1.0, 2.0]);
        assert!(p.is_pinned());
        assert!(p.inv_mass.abs() < 1e-15);
    }
    #[test]
    fn test_wind_draping_apply_empty() {
        let wd = WindDraping::new([1.0, 0.0, 0.0], 1.225, 1.0);
        let mut particles = vec![ClothParticle::new([0.0; 3], 0.01)];
        wd.apply(&mut particles, &[]);
    }
    #[test]
    fn test_cloth_edge_strain_at_rest() {
        let particles = vec![
            ClothParticle::new([0.0, 0.0, 0.0], 0.01),
            ClothParticle::new([0.1, 0.0, 0.0], 0.01),
        ];
        let edge = ClothEdge::new(0, 1, 0.1, 0.5);
        assert!(edge.strain(&particles).abs() < 1e-10);
    }
}
/// Approximate erfc for fabric thermal calculations.
#[allow(dead_code)]
pub(super) fn erfc_fabric(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_fabric(-x);
    }
    if x > 5.0 {
        return 0.0;
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}
#[cfg(test)]
mod expanded_cloth_tests {

    use crate::cloth_advanced::ClothBvhLeaf;
    use crate::cloth_advanced::ClothCollisionBroadphase;

    use crate::cloth_advanced::DrapeSimParams;

    use crate::cloth_advanced::FabricFailureCriterion;
    use crate::cloth_advanced::FabricThermal;

    use crate::cloth_advanced::FitMetrics;
    use crate::cloth_advanced::GarmentPanel;

    use crate::cloth_advanced::KnittedFabric;

    use crate::cloth_advanced::WeavePattern;

    use crate::cloth_advanced::WovenFabric;
    use crate::cloth_advanced::YarnSection;
    #[test]
    fn test_yarn_axial_stiffness() {
        let yarn = YarnSection::circular(0.15e-3, 8e9, 3e9);
        let ea = yarn.axial_stiffness();
        assert!(ea > 0.0, "ea={ea}");
    }
    #[test]
    fn test_weave_float_length() {
        assert_eq!(WeavePattern::Plain.float_length(), 1);
        assert_eq!(WeavePattern::Twill(3, 1).float_length(), 3);
    }
    #[test]
    fn test_weave_crimp() {
        let c = WeavePattern::Plain.crimp_factor();
        assert!(c > 0.0 && c < 1.0, "c={c}");
    }
    #[test]
    fn test_woven_cotton_stiffness() {
        let fabric = WovenFabric::cotton_plain();
        assert!(fabric.warp_stiffness() > 0.0);
        assert!(fabric.weft_stiffness() > 0.0);
    }
    #[test]
    fn test_woven_thickness() {
        let fabric = WovenFabric::cotton_plain();
        assert!(fabric.thickness() > 0.0);
    }
    #[test]
    fn test_knit_stitch_density() {
        let knit = KnittedFabric::jersey_cotton();
        assert!(knit.stitch_density() > 0.0);
    }
    #[test]
    fn test_knit_extensibility() {
        let jersey = KnittedFabric::jersey_cotton();
        let mut rib = KnittedFabric::jersey_cotton();
        rib.is_jersey = false;
        assert!(jersey.extensibility() > rib.extensibility());
    }
    #[test]
    fn test_panel_area() {
        let panel = GarmentPanel::rectangle("front", 0.5, 0.8);
        assert!((panel.area() - 0.4).abs() < 1e-10, "area={}", panel.area());
    }
    #[test]
    fn test_panel_perimeter() {
        let panel = GarmentPanel::rectangle("back", 1.0, 1.0);
        assert!((panel.perimeter() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_panel_centroid() {
        let panel = GarmentPanel::rectangle("sleeve", 2.0, 2.0);
        let c = panel.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fit_score_range() {
        let mut fm = FitMetrics {
            ease: 0.03,
            avg_strain: 0.05,
            max_strain: 0.1,
            pressure_areas: 0,
            fit_score: 0.0,
        };
        fm.compute_score();
        assert!(
            fm.fit_score >= 0.0 && fm.fit_score <= 1.0,
            "score={}",
            fm.fit_score
        );
    }
    #[test]
    fn test_drape_coefficient() {
        let drape = DrapeSimParams::muslin();
        let dc = drape.drape_coefficient();
        assert!((0.0..=1.0).contains(&dc), "dc={dc}");
    }
    #[test]
    fn test_cotton_clo() {
        let cotton = FabricThermal::cotton();
        let clo = cotton.clo_value();
        assert!(clo > 0.0, "clo={clo}");
    }
    #[test]
    fn test_thermal_conductivity_comparison() {
        let cotton = FabricThermal::cotton();
        let poly = FabricThermal::polyester();
        assert!(poly.conductivity < cotton.conductivity);
    }
    #[test]
    fn test_thermal_depth_profile() {
        let cotton = FabricThermal::cotton();
        let t1 = cotton.temperature_at_depth(100.0, 20.0, 0.0, 1.0);
        let t2 = cotton.temperature_at_depth(100.0, 20.0, cotton.thickness, 1.0);
        assert!(t1 >= t2, "surface should be hotter: t1={t1}, t2={t2}");
    }
    #[test]
    fn test_failure_safe_at_zero() {
        let fc = FabricFailureCriterion::cotton_plain();
        assert!(!fc.is_failed(0.0, 0.0, 0.0));
    }
    #[test]
    fn test_failure_exceeds_strength() {
        let fc = FabricFailureCriterion::cotton_plain();
        assert!(fc.is_failed(fc.sigma_warp_t * 1.1, 0.0, 0.0));
    }
    #[test]
    fn test_bvh_leaf_overlap() {
        let l1 =
            ClothBvhLeaf::from_triangle(0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], 0.0);
        let l2 =
            ClothBvhLeaf::from_triangle(1, [0.4, 0.0, 0.0], [1.4, 0.0, 0.0], [0.9, 1.0, 0.0], 0.0);
        assert!(l1.overlaps(&l2), "Should overlap in XY region");
    }
    #[test]
    fn test_bvh_leaf_no_overlap() {
        let l1 =
            ClothBvhLeaf::from_triangle(0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], 0.0);
        let l2 =
            ClothBvhLeaf::from_triangle(1, [5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.5, 6.0, 5.0], 0.0);
        assert!(!l1.overlaps(&l2), "Should not overlap");
    }
    #[test]
    fn test_broadphase_pairs() {
        let mut bp = ClothCollisionBroadphase::new();
        bp.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], 0.01);
        bp.add_triangle([0.4, 0.0, 0.0], [1.4, 0.0, 0.0], [0.9, 1.0, 0.0], 0.01);
        bp.find_pairs();
        assert!(bp.pair_count() > 0, "Should find overlap pair");
    }
    #[test]
    fn test_woven_shear_stiffness() {
        let fabric = WovenFabric::cotton_plain();
        let ks = fabric.shear_stiffness();
        assert!(ks >= 0.0, "ks={ks}");
    }
    #[test]
    fn test_knit_bending_rigidity() {
        let knit = KnittedFabric::jersey_cotton();
        assert!(knit.bending_rigidity() > 0.0);
    }
    #[test]
    fn test_drape_fold_count() {
        let drape = DrapeSimParams::muslin();
        assert_eq!(drape.fold_count(), drape.petal_count);
    }
    #[test]
    fn test_panel_net_perimeter() {
        let panel = GarmentPanel::rectangle("test", 1.0, 1.0);
        assert!(panel.net_perimeter() < panel.perimeter());
    }
    #[test]
    fn test_failure_shear() {
        let fc = FabricFailureCriterion::cotton_plain();
        assert!(fc.is_failed(0.0, 0.0, fc.tau_s * 1.05));
    }
}
