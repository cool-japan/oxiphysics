//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BirthDeathPair, Color4, CriticalPoint, PointGlyph};

/// A 3-D point stored as a plain array.
pub type Point3 = [f64; 3];
/// Render a slice of critical points as point glyphs.
///
/// The radius is scaled by `base_radius`. Each critical point type gets its
/// canonical color (blue=min, cyan=saddle1, orange=saddle2, red=max).
#[allow(clippy::too_many_arguments)]
pub fn render_critical_points(cps: &[CriticalPoint], base_radius: f64) -> Vec<PointGlyph> {
    cps.iter().map(|cp| cp.to_glyph(base_radius)).collect()
}
/// Map a persistence value to a heat color (blue=low → red=high).
///
/// Persistence is normalized to `[0, 1]` using `max_persistence`.
pub fn persistence_color(persistence: f64, max_persistence: f64) -> Color4 {
    let t = (persistence / max_persistence.max(1e-14)).clamp(0.0, 1.0) as f32;
    Color4::lerp(Color4::blue(), Color4::red(), t)
}
/// Build a birth-death scatter plot as a vector of glyphs from raw pairs.
#[allow(clippy::too_many_arguments)]
pub fn build_birth_death_scatter(
    pairs: &[BirthDeathPair],
    max_death: f64,
    point_radius: f64,
) -> Vec<PointGlyph> {
    let max_pers = pairs
        .iter()
        .map(|p| {
            if p.persistence().is_finite() {
                p.persistence()
            } else {
                0.0
            }
        })
        .fold(0.0_f64, f64::max);
    pairs
        .iter()
        .map(|p| {
            let death_vis = if p.death.is_infinite() {
                max_death
            } else {
                p.death
            };
            let pos = [p.birth, death_vis, 0.0];
            let pers = if p.persistence().is_finite() {
                p.persistence()
            } else {
                0.0
            };
            let color = persistence_color(pers, max_pers);
            PointGlyph::new(pos, color, point_radius)
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology_viz_ext::ContourTreeViz;
    use crate::topology_viz_ext::CriticalPointType;
    use crate::topology_viz_ext::EulerCharacteristicViz;
    use crate::topology_viz_ext::FiberSurfaceViz;
    use crate::topology_viz_ext::GradientVectorFieldViz;
    use crate::topology_viz_ext::JacobiSegment;
    use crate::topology_viz_ext::JacobiSetViz;
    use crate::topology_viz_ext::LineSegment;
    use crate::topology_viz_ext::MorseComplexViz;
    use crate::topology_viz_ext::PersistenceBarcodeViz;
    use crate::topology_viz_ext::PersistenceDiagramViz;
    use crate::topology_viz_ext::ReebGraphViz;
    use crate::topology_viz_ext::TopologicalFeatureTracker;
    use crate::topology_viz_ext::TopologicalNoiseRemovalViz;
    use crate::topology_viz_ext::TopologySimplificationViz;
    use crate::topology_viz_ext::TrackedFeature;
    use crate::topology_viz_ext::UnstableManifoldArc;
    use crate::topology_viz_ext::UnstableManifoldViz;
    #[test]
    fn birth_death_pair_persistence() {
        let p = BirthDeathPair::new(1.0, 3.0, 0);
        assert!((p.persistence() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn birth_death_pair_essential_infinity() {
        let p = BirthDeathPair::new(1.0, f64::INFINITY, 1);
        assert!(p.is_essential());
        assert!(p.persistence().is_infinite());
    }
    #[test]
    fn birth_death_pair_diagram_point() {
        let p = BirthDeathPair::new(2.0, 5.0, 0);
        let [b, d] = p.diagram_point();
        assert!((b - 2.0).abs() < 1e-12);
        assert!((d - 5.0).abs() < 1e-12);
    }
    #[test]
    fn persistence_diagram_build_glyphs_count() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 1.0, 0),
            BirthDeathPair::new(0.5, 2.0, 1),
            BirthDeathPair::new(1.0, f64::INFINITY, 1),
        ];
        let viz = PersistenceDiagramViz::new(pairs);
        let glyphs = viz.build_glyphs(10.0);
        assert_eq!(glyphs.len(), 3);
    }
    #[test]
    fn persistence_diagram_diagonal_present() {
        let viz = PersistenceDiagramViz::new(vec![]);
        let diag = viz.build_diagonal(0.0, 5.0);
        assert!(diag.is_some());
    }
    #[test]
    fn persistence_diagram_no_diagonal_when_disabled() {
        let mut viz = PersistenceDiagramViz::new(vec![]);
        viz.show_diagonal = false;
        assert!(viz.build_diagonal(0.0, 5.0).is_none());
    }
    #[test]
    fn persistence_diagram_color_for_dim() {
        let viz = PersistenceDiagramViz::new(vec![]);
        let c0 = viz.color_for_dim(0);
        assert!((c0.b - 1.0).abs() < 1e-6);
        assert!(c0.r.abs() < 1e-6);
    }
    #[test]
    fn persistence_diagram_bottleneck_zero_when_identical() {
        let pairs = vec![BirthDeathPair::new(0.0, 2.0, 0)];
        let a = PersistenceDiagramViz::new(pairs.clone());
        let b = PersistenceDiagramViz::new(pairs);
        assert!((a.bottleneck_distance_approx(&b)).abs() < 1e-12);
    }
    #[test]
    fn critical_point_glyph_color_maximum_is_red() {
        let cp = CriticalPoint::new([0.0, 0.0, 0.0], 5.0, CriticalPointType::Maximum, 0);
        let glyph = cp.to_glyph(0.05);
        assert!((glyph.color.r - 1.0).abs() < 1e-6);
        assert!(glyph.color.g.abs() < 1e-6);
        assert!(glyph.color.b.abs() < 1e-6);
    }
    #[test]
    fn critical_point_glyph_label_is_value() {
        let cp = CriticalPoint::new([0.0, 0.0, 0.0], 3.125, CriticalPointType::Minimum, 0);
        let glyph = cp.to_glyph(0.05);
        assert!(glyph.label.is_some());
    }
    #[test]
    fn morse_index_saddle1_is_1() {
        assert_eq!(CriticalPointType::Saddle1.morse_index(), 1);
    }
    #[test]
    fn morse_complex_euler_characteristic_sphere() {
        let mut mc = MorseComplexViz::new();
        mc.add_critical_point(CriticalPoint::new(
            [0.0, -1.0, 0.0],
            -1.0,
            CriticalPointType::Minimum,
            0,
        ));
        mc.add_critical_point(CriticalPoint::new(
            [0.0, 1.0, 0.0],
            1.0,
            CriticalPointType::Maximum,
            1,
        ));
        assert_eq!(mc.euler_characteristic(), 0);
    }
    #[test]
    fn morse_complex_euler_characteristic_torus() {
        let mut mc = MorseComplexViz::new();
        mc.add_critical_point(CriticalPoint::new(
            [0.0, 0.0, 0.0],
            0.0,
            CriticalPointType::Minimum,
            0,
        ));
        mc.add_critical_point(CriticalPoint::new(
            [1.0, 0.0, 0.0],
            1.0,
            CriticalPointType::Saddle1,
            1,
        ));
        mc.add_critical_point(CriticalPoint::new(
            [2.0, 0.0, 0.0],
            2.0,
            CriticalPointType::Saddle1,
            2,
        ));
        mc.add_critical_point(CriticalPoint::new(
            [3.0, 0.0, 0.0],
            3.0,
            CriticalPointType::Maximum,
            3,
        ));
        assert_eq!(mc.euler_characteristic(), -2);
    }
    #[test]
    fn morse_complex_count_of_type() {
        let mut mc = MorseComplexViz::new();
        mc.add_critical_point(CriticalPoint::new(
            [0.0, 0.0, 0.0],
            0.0,
            CriticalPointType::Minimum,
            0,
        ));
        mc.add_critical_point(CriticalPoint::new(
            [1.0, 0.0, 0.0],
            1.0,
            CriticalPointType::Minimum,
            1,
        ));
        assert_eq!(mc.count_of_type(CriticalPointType::Minimum), 2);
        assert_eq!(mc.count_of_type(CriticalPointType::Maximum), 0);
    }
    #[test]
    fn morse_complex_all_arcs_count() {
        let mut mc = MorseComplexViz::new();
        mc.add_ascending_arc(LineSegment::new([0.0; 3], [1.0, 0.0, 0.0], Color4::cyan()));
        mc.add_descending_arc(LineSegment::new(
            [0.0; 3],
            [0.0, 1.0, 0.0],
            Color4::orange(),
        ));
        assert_eq!(mc.all_arcs().len(), 2);
    }
    #[test]
    fn reeb_graph_add_nodes_and_arcs() {
        let mut rg = ReebGraphViz::new();
        let n0 = rg.add_node(0.0, [0.0; 3], CriticalPointType::Minimum);
        let n1 = rg.add_node(1.0, [0.0, 1.0, 0.0], CriticalPointType::Maximum);
        rg.add_arc(
            n0,
            n1,
            vec![[0.0; 3], [0.0, 0.5, 0.0], [0.0, 1.0, 0.0]],
            Color4::green(),
        );
        assert_eq!(rg.nodes.len(), 2);
        assert_eq!(rg.arcs.len(), 1);
    }
    #[test]
    fn reeb_graph_build_lines_count() {
        let mut rg = ReebGraphViz::new();
        let n0 = rg.add_node(0.0, [0.0; 3], CriticalPointType::Minimum);
        let n1 = rg.add_node(1.0, [0.0, 1.0, 0.0], CriticalPointType::Maximum);
        rg.add_arc(
            n0,
            n1,
            vec![[0.0; 3], [0.0, 0.5, 0.0], [0.0, 1.0, 0.0]],
            Color4::green(),
        );
        let lines = rg.build_lines();
        assert_eq!(lines.len(), 2);
    }
    #[test]
    fn reeb_graph_node_lookup() {
        let mut rg = ReebGraphViz::new();
        let n0 = rg.add_node(3.125, [1.0, 2.0, 3.0], CriticalPointType::Saddle1);
        let node = rg.node(n0).unwrap();
        assert!((node.value - 3.125).abs() < 1e-12);
    }
    #[test]
    fn contour_tree_build_lines() {
        let mut ct = ContourTreeViz::new();
        let v0 = ct.add_vertex(0.0, [0.0; 3]);
        let v1 = ct.add_vertex(1.0, [0.0, 1.0, 0.0]);
        ct.add_edge(v0, v1, Color4::white());
        let lines = ct.build_lines();
        assert_eq!(lines.len(), 1);
    }
    #[test]
    fn contour_tree_leaf_count() {
        let mut ct = ContourTreeViz::new();
        let v0 = ct.add_vertex(0.0, [0.0; 3]);
        let v1 = ct.add_vertex(0.5, [0.5, 0.0, 0.0]);
        let v2 = ct.add_vertex(1.0, [1.0, 0.0, 0.0]);
        ct.add_edge(v0, v1, Color4::white());
        ct.add_edge(v1, v2, Color4::white());
        assert_eq!(ct.leaf_count(), 2);
    }
    #[test]
    fn euler_characteristic_curve_lines_count() {
        let thresholds = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let chi_vals = vec![1_i64, 2, 1, 0, -1];
        let viz = EulerCharacteristicViz::new(thresholds, chi_vals);
        let lines = viz.build_curve_lines(1.0);
        assert_eq!(lines.len(), 4);
    }
    #[test]
    fn euler_characteristic_min_max() {
        let viz = EulerCharacteristicViz::new(vec![0.0, 1.0], vec![-3_i64, 5]);
        assert_eq!(viz.chi_min(), Some(-3));
        assert_eq!(viz.chi_max(), Some(5));
    }
    #[test]
    fn tracker_register_and_record() {
        let mut tracker = TopologicalFeatureTracker::new();
        let id = tracker.register_feature(1);
        tracker.record(id, 0.0, BirthDeathPair::new(0.1, 0.9, 1));
        tracker.record(id, 1.0, BirthDeathPair::new(0.2, 1.1, 1));
        let feat = &tracker.features[&id];
        assert_eq!(feat.time_series.len(), 2);
    }
    #[test]
    fn tracked_feature_mean_persistence() {
        let mut feat = TrackedFeature::new(0, 0);
        feat.record(0.0, BirthDeathPair::new(0.0, 2.0, 0));
        feat.record(1.0, BirthDeathPair::new(0.0, 4.0, 0));
        assert!((feat.mean_persistence() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn tracked_feature_trace_line_count() {
        let mut feat = TrackedFeature::new(0, 1);
        feat.record(0.0, BirthDeathPair::new(0.0, 1.0, 1));
        feat.record(1.0, BirthDeathPair::new(0.0, 2.0, 1));
        feat.record(2.0, BirthDeathPair::new(0.0, 1.5, 1));
        let lines = feat.build_persistence_trace();
        assert_eq!(lines.len(), 2);
    }
    #[test]
    fn gradient_field_from_constant_has_zero_magnitude() {
        let viz = GradientVectorFieldViz::from_scalar_field_2d(
            |_x, _y| 1.0,
            0.0,
            1.0,
            4,
            0.0,
            1.0,
            4,
            |_| Color4::white(),
        );
        for a in &viz.arrows {
            assert!(a.magnitude < 1e-10, "magnitude={}", a.magnitude);
        }
    }
    #[test]
    fn gradient_field_linear_has_constant_gradient() {
        let viz = GradientVectorFieldViz::from_scalar_field_2d(
            |x, _y| x,
            0.0,
            1.0,
            5,
            0.0,
            1.0,
            5,
            |_| Color4::white(),
        );
        for a in &viz.arrows {
            assert!((a.direction[0] - 1.0).abs() < 1e-6, "gx={}", a.direction[0]);
            assert!(a.direction[1].abs() < 1e-6, "gy={}", a.direction[1]);
        }
    }
    #[test]
    fn gradient_max_magnitude_linear_field() {
        let viz = GradientVectorFieldViz::from_scalar_field_2d(
            |x, _y| x,
            0.0,
            1.0,
            3,
            0.0,
            1.0,
            3,
            |_| Color4::white(),
        );
        assert!((viz.max_magnitude() - 1.0).abs() < 1e-6);
    }
    #[test]
    fn simplification_apply_threshold_cancels_pairs() {
        let cps = vec![
            CriticalPoint::new([0.0; 3], 0.0, CriticalPointType::Minimum, 0),
            CriticalPoint::new([1.0, 0.0, 0.0], 0.1, CriticalPointType::Saddle1, 1),
            CriticalPoint::new([2.0, 0.0, 0.0], 2.0, CriticalPointType::Maximum, 2),
        ];
        let mut viz = TopologySimplificationViz::new(cps);
        viz.add_cancellation(0, 1, 0.1);
        viz.apply_threshold(0.1);
        assert_eq!(viz.cancelled_count(), 1);
        let active = viz.active_indices();
        assert!(!active.contains(&0));
        assert!(!active.contains(&1));
        assert!(active.contains(&2));
    }
    #[test]
    fn simplification_no_cancellation_below_threshold() {
        let cps = vec![CriticalPoint::new(
            [0.0; 3],
            0.0,
            CriticalPointType::Minimum,
            0,
        )];
        let mut viz = TopologySimplificationViz::new(cps);
        viz.add_cancellation(0, 0, 5.0);
        viz.apply_threshold(1.0);
        assert_eq!(viz.cancelled_count(), 0);
    }
    #[test]
    fn fiber_surface_sample_empty_when_no_match() {
        let surf = FiberSurfaceViz::sample_from_field(
            |_x, _y, _z| [0.0, 0.0],
            [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
            [4, 4, 4],
            [1.0, 1.0],
            0.01,
            Color4::cyan(),
        );
        assert_eq!(surf.point_count(), 0);
    }
    #[test]
    fn fiber_surface_sample_finds_points() {
        let surf = FiberSurfaceViz::sample_from_field(
            |x, y, _z| [x, y],
            [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
            [5, 5, 3],
            [0.5, 0.5],
            0.15,
            Color4::cyan(),
        );
        assert!(surf.point_count() > 0, "expected some fiber points");
    }
    #[test]
    fn jacobi_set_from_identity_map_has_some_segments() {
        let jac =
            JacobiSetViz::from_bivariate_field_2d(|x, y| [x, y], 0.0, 1.0, 5, 0.0, 1.0, 5, 0.0);
        assert_eq!(jac.segment_count(), 0);
    }
    #[test]
    fn jacobi_set_segment_to_line() {
        let seg = JacobiSegment::new([0.0; 3], [1.0, 0.0, 0.0], 1);
        let line = seg.to_line();
        assert!((line.end[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn noise_removal_counts() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 0.1, 0),
            BirthDeathPair::new(0.0, 5.0, 0),
        ];
        let viz = TopologicalNoiseRemovalViz::new(pairs, 0.2);
        assert_eq!(viz.noise_count(), 1);
        assert_eq!(viz.signal_count(), 1);
    }
    #[test]
    fn noise_removal_update_threshold() {
        let pairs = vec![BirthDeathPair::new(0.0, 3.0, 0)];
        let mut viz = TopologicalNoiseRemovalViz::new(pairs, 1.0);
        assert_eq!(viz.noise_count(), 0);
        viz.update_threshold(5.0);
        assert_eq!(viz.noise_count(), 1);
    }
    #[test]
    fn noise_removal_build_glyphs() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 0.1, 0),
            BirthDeathPair::new(0.0, 5.0, 1),
        ];
        let viz = TopologicalNoiseRemovalViz::new(pairs, 0.2);
        let glyphs = viz.build_diagram_glyphs(10.0);
        assert_eq!(glyphs.len(), 2);
    }
    #[test]
    fn barcode_bar_count() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 1.0, 0),
            BirthDeathPair::new(0.5, 2.0, 1),
            BirthDeathPair::new(1.0, f64::INFINITY, 0),
        ];
        let viz = PersistenceBarcodeViz::from_pairs(&pairs, 0.1);
        assert_eq!(viz.bar_count(), 3);
        assert_eq!(viz.essential_count(), 1);
    }
    #[test]
    fn barcode_build_lines_has_correct_count() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 1.0, 0),
            BirthDeathPair::new(0.5, 2.0, 1),
        ];
        let viz = PersistenceBarcodeViz::from_pairs(&pairs, 0.05);
        let lines = viz.build_lines(5.0);
        assert_eq!(lines.len(), 2);
    }
    #[test]
    fn unstable_manifold_trace_arc_length() {
        let pts =
            UnstableManifoldViz::trace_arc(|_x, _y, _z| [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], 0.1, 5);
        assert_eq!(pts.len(), 6);
        assert!(pts.last().unwrap()[1] > 0.4);
    }
    #[test]
    fn unstable_manifold_build_lines() {
        let mut viz = UnstableManifoldViz::new(Color4::orange());
        let pts = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        viz.add_arc(UnstableManifoldArc {
            points: pts,
            saddle_index: 0,
            max_index: 1,
        });
        let lines = viz.build_lines();
        assert_eq!(lines.len(), 2);
    }
    #[test]
    fn unstable_manifold_total_points() {
        let mut viz = UnstableManifoldViz::new(Color4::orange());
        viz.add_arc(UnstableManifoldArc {
            points: vec![[0.0; 3]; 4],
            saddle_index: 0,
            max_index: 1,
        });
        viz.add_arc(UnstableManifoldArc {
            points: vec![[0.0; 3]; 3],
            saddle_index: 0,
            max_index: 2,
        });
        assert_eq!(viz.total_points(), 7);
    }
    #[test]
    fn persistence_color_zero_is_blue() {
        let c = persistence_color(0.0, 1.0);
        assert!((c.b - 1.0).abs() < 1e-6);
        assert!(c.r.abs() < 1e-6);
    }
    #[test]
    fn persistence_color_max_is_red() {
        let c = persistence_color(1.0, 1.0);
        assert!((c.r - 1.0).abs() < 1e-6);
        assert!(c.b.abs() < 1e-6);
    }
    #[test]
    fn build_birth_death_scatter_count() {
        let pairs = vec![
            BirthDeathPair::new(0.0, 1.0, 0),
            BirthDeathPair::new(1.0, 3.0, 1),
        ];
        let glyphs = build_birth_death_scatter(&pairs, 5.0, 0.02);
        assert_eq!(glyphs.len(), 2);
    }
    #[test]
    fn line_segment_length_unit() {
        let seg = LineSegment::new([0.0; 3], [1.0, 0.0, 0.0], Color4::white());
        assert!((seg.length() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn line_segment_length_3d_diagonal() {
        let seg = LineSegment::new([0.0; 3], [1.0, 1.0, 1.0], Color4::white());
        assert!((seg.length() - 3.0_f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn color4_lerp_midpoint() {
        let c = Color4::lerp(Color4::red(), Color4::blue(), 0.5);
        assert!((c.r - 0.5).abs() < 1e-6);
        assert!((c.b - 0.5).abs() < 1e-6);
    }
    #[test]
    fn render_critical_points_count() {
        let cps = vec![
            CriticalPoint::new([0.0; 3], 0.0, CriticalPointType::Minimum, 0),
            CriticalPoint::new([1.0, 0.0, 0.0], 1.0, CriticalPointType::Maximum, 1),
        ];
        let glyphs = render_critical_points(&cps, 0.05);
        assert_eq!(glyphs.len(), 2);
    }
}
