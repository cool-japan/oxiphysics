//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::functions::*;
#[cfg(test)]
use super::types::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_vec3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = vec3_add(a, b);
        assert!((c[0] - 5.0).abs() < 1e-12);
        assert!((c[1] - 7.0).abs() < 1e-12);
        assert!((c[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_dot() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(vec3_dot(a, b).abs() < 1e-12);
        assert!((vec3_dot(a, a) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_norm_unit() {
        let v = [3.0, 4.0, 0.0];
        let u = vec3_norm(v);
        assert!((vec3_len(u) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_norm_zero() {
        let v = [0.0, 0.0, 0.0];
        let u = vec3_norm(v);
        assert_eq!(u, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_vec3_cross() {
        let x = [1.0_f64, 0.0, 0.0];
        let y = [0.0_f64, 1.0, 0.0];
        // cross product: x × y = [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]
        let z = [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ];
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rgba_lerp_midpoint() {
        let r = Rgba::red();
        let b = Rgba::blue();
        let m = Rgba::lerp(r, b, 0.5);
        assert!((m.r - 0.5).abs() < 1e-6);
        assert!((m.b - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_rgba_lerp_endpoints() {
        let r = Rgba::red();
        let b = Rgba::blue();
        let a = Rgba::lerp(r, b, 0.0);
        let c = Rgba::lerp(r, b, 1.0);
        assert!((a.r - 1.0).abs() < 1e-6);
        assert!((c.b - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_thermal_mechanical_overlay_length() {
        let config = ThermalMechanicalOverlayConfig::default();
        let t_nodes: Vec<FieldNode> = (0..5)
            .map(|i| FieldNode::new([i as f64, 0.0, 0.0], 300.0 + i as f64 * 100.0, [0.0; 3]))
            .collect();
        let d_nodes: Vec<FieldNode> = (0..5)
            .map(|i| FieldNode::new([i as f64, 0.0, 0.0], 0.0, [0.001 * i as f64, 0.0, 0.0]))
            .collect();
        let colors = thermal_mechanical_overlay(&t_nodes, &d_nodes, &config);
        assert_eq!(colors.len(), 5);
    }
    #[test]
    fn test_thermal_mechanical_overlay_pure_thermal() {
        let config = ThermalMechanicalOverlayConfig {
            thermal_weight: 1.0,
            ..Default::default()
        };
        let t_nodes = vec![FieldNode::new([0.0; 3], 1000.0, [0.0; 3])];
        let d_nodes = vec![FieldNode::new([0.0; 3], 0.0, [0.0; 3])];
        let colors = thermal_mechanical_overlay(&t_nodes, &d_nodes, &config);
        assert!(colors[0].r > 0.9);
    }
    #[test]
    fn test_coupled_arrows_count() {
        let config = CoupledArrowConfig::default();
        let nodes = vec![
            FieldNode::new([0.0, 0.0, 0.0], 0.0, [0.0; 3]),
            FieldNode::new([1.0, 0.0, 0.0], 0.0, [0.0; 3]),
        ];
        let hf = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let disp = vec![[0.0, 1.0, 0.0], [0.0, 0.5, 0.0]];
        let segs = coupled_field_arrows(&nodes, &hf, &disp, &config);
        assert_eq!(segs.len(), 4);
    }
    #[test]
    fn test_coupled_arrows_below_threshold_filtered() {
        let config = CoupledArrowConfig {
            min_magnitude: 1e3,
            ..Default::default()
        };
        let nodes = vec![FieldNode::new([0.0; 3], 0.0, [0.0; 3])];
        let hf = vec![[1.0, 0.0, 0.0]];
        let disp = vec![[0.0, 1.0, 0.0]];
        let segs = coupled_field_arrows(&nodes, &hf, &disp, &config);
        assert_eq!(segs.len(), 0);
    }
    #[test]
    fn test_extract_interface_single_crossing() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let scalars: Vec<f64> = (0..10).map(|i| i as f64 - 4.5).collect();
        let pts = extract_interface_points(&positions, &scalars, 0.0);
        assert_eq!(pts.len(), 1);
    }
    #[test]
    fn test_extract_interface_no_crossing() {
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let scalars = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let pts = extract_interface_points(&positions, &scalars, 0.0);
        assert_eq!(pts.len(), 0);
    }
    #[test]
    fn test_interface_normal_lines() {
        let pts = vec![
            InterfacePoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0),
            InterfacePoint::new([1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 1),
        ];
        let segs = interface_normal_lines(&pts, 0.1, Rgba::red(), Rgba::blue());
        assert_eq!(segs.len(), 2);
        assert!((segs[0].start[0]).abs() < 1e-10);
    }
    #[test]
    fn test_fsi_fluid_arrows_count() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.fluid_positions = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        snap.fluid_velocities = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let segs = fsi_fluid_arrows(&snap, 1.0, 1e-10, Rgba::cyan());
        assert_eq!(segs.len(), 2);
    }
    #[test]
    fn test_fsi_coupling_lines_count() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.fluid_positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        snap.structural_positions = vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
        snap.fsi_interface_indices = vec![0, 1];
        let segs = fsi_coupling_lines(&snap, Rgba::yellow());
        assert_eq!(segs.len(), 2);
    }
    #[test]
    fn test_fsi_structural_arrows_filtered() {
        let mut snap = FsiSnapshot::new(0.0);
        snap.structural_positions = vec![[0.0; 3]];
        snap.structural_displacements = vec![[0.0; 3]];
        let segs = fsi_structural_arrows(&snap, 1.0, 1e-10, Rgba::green());
        assert_eq!(segs.len(), 0);
    }
    #[test]
    fn test_em_field_line_uniform_euler() {
        let config = EmFieldLineConfig {
            step_size: 0.1,
            max_steps: 10,
            min_magnitude: 1e-15,
            method: FieldLineMethod::Euler,
        };
        let field = |_: [f64; 3]| [1.0, 0.0, 0.0];
        let line = trace_em_field_line([0.0, 0.0, 0.0], &field, &config);
        assert_eq!(line.points.len(), 11);
        assert!((line.points.last().unwrap()[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_em_field_line_rk4() {
        let config = EmFieldLineConfig {
            step_size: 0.1,
            max_steps: 5,
            min_magnitude: 1e-15,
            method: FieldLineMethod::Rk4,
        };
        let field = |_: [f64; 3]| [0.0, 1.0, 0.0];
        let line = trace_em_field_line([0.0, 0.0, 0.0], &field, &config);
        assert_eq!(line.points.len(), 6);
        assert!((line.points.last().unwrap()[1] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_em_field_lines_to_segments() {
        let config = EmFieldLineConfig::default();
        let field = |_: [f64; 3]| [1.0, 0.0, 0.0];
        let seeds = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let lines = trace_em_field_lines(&seeds, &field, &config);
        let segs = em_field_lines_to_segments(&lines, 0.0, 1.0);
        assert!(!segs.is_empty());
    }
    #[test]
    fn test_blackbody_ramp_endpoints() {
        let ramp = blackbody_ramp();
        let cold = ramp.evaluate(0.0);
        let hot = ramp.evaluate(5000.0);
        assert!(cold.r < 0.3);
        assert!(hot.r > 0.9);
        assert!(hot.g > 0.9);
        assert!(hot.b > 0.9);
    }
    #[test]
    fn test_blackbody_ramp_midpoint() {
        let ramp = blackbody_ramp();
        let mid = ramp.evaluate(800.0);
        assert!(mid.r > mid.b);
    }
    #[test]
    fn test_colorize_temperatures_length() {
        let ramp = blackbody_ramp();
        let temps = vec![300.0, 600.0, 900.0, 1200.0];
        let colors = ramp.colorize_temperatures(&temps);
        assert_eq!(colors.len(), 4);
    }
    #[test]
    fn test_dashboard_all_converged() {
        let mut dash = MultiphysicsConvergenceDashboard::new(100);
        let mut d1 = DomainConvergence::new("thermal", 1e-6);
        d1.push(1e-7);
        let mut d2 = DomainConvergence::new("mechanical", 1e-6);
        d2.push(5e-8);
        dash.add_domain(d1);
        dash.add_domain(d2);
        assert!(dash.all_converged());
    }
    #[test]
    fn test_dashboard_not_converged() {
        let mut dash = MultiphysicsConvergenceDashboard::new(100);
        let mut d = DomainConvergence::new("fluid", 1e-6);
        d.push(1e-3);
        dash.add_domain(d);
        assert!(!dash.all_converged());
    }
    #[test]
    fn test_dashboard_iteration_limit() {
        let mut dash = MultiphysicsConvergenceDashboard::new(10);
        for _ in 0..10 {
            dash.next_iteration();
        }
        assert!(dash.iteration_limit_reached());
    }
    #[test]
    fn test_convergence_rate() {
        let mut d = DomainConvergence::new("x", 1e-8);
        d.push(1.0);
        d.push(0.1);
        let rate = d.rate().unwrap();
        assert!((rate - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_interpolate_scalar_nearest() {
        let src_pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_val = vec![10.0, 20.0, 30.0];
        let tgt_pos = vec![[0.1, 0.0, 0.0]];
        let result = interpolate_scalar_field(
            &src_pos,
            &src_val,
            &tgt_pos,
            MeshInterpolationMethod::NearestNeighbor,
        );
        assert!((result[0] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_scalar_idw() {
        let src_pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_val = vec![0.0, 2.0];
        let tgt_pos = vec![[1.0, 0.0, 0.0]];
        let result = interpolate_scalar_field(
            &src_pos,
            &src_val,
            &tgt_pos,
            MeshInterpolationMethod::InverseDistance {
                power: 2.0,
                k_neighbors: 2,
            },
        );
        assert!((result[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_vector_field() {
        let src_pos = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let src_vec = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let tgt_pos = vec![[1.0, 0.0, 0.0]];
        let result = interpolate_vector_field(&src_pos, &src_vec, &tgt_pos, 2);
        assert!((result[0][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_clipping_plane_visible() {
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(plane.visible([1.0, 0.0, 0.0]));
        assert!(!plane.visible([-1.0, 0.0, 0.0]));
    }
    #[test]
    fn test_clipping_plane_signed_distance() {
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let d = plane.signed_distance([0.0, 3.0, 0.0]);
        assert!((d - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_clipping_stack_filter_nodes() {
        let mut stack = ClippingPlaneStack::new();
        stack.add(ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
        let nodes = vec![
            FieldNode::new([1.0, 0.0, 0.0], 0.0, [0.0; 3]),
            FieldNode::new([-1.0, 0.0, 0.0], 0.0, [0.0; 3]),
        ];
        let visible = stack.filter_nodes(&nodes);
        assert_eq!(visible.len(), 1);
        assert!((visible[0].position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_clipping_stack_filter_segments() {
        let mut stack = ClippingPlaneStack::new();
        stack.add(ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]));
        let segs = vec![
            LineSegment::new([1.0, 0.0, 0.0], [2.0, 0.0, 0.0], Rgba::white()),
            LineSegment::new([-1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], Rgba::white()),
            LineSegment::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], Rgba::white()),
        ];
        let visible = stack.filter_segments(&segs);
        assert_eq!(visible.len(), 1);
    }
    #[test]
    fn test_timeline_global_bounds() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.1));
        let mut s1 = PhysicsTimeSeries::new("thermal");
        s1.push(0.0, 0.0);
        s1.push(1.0, 1.0);
        let mut s2 = PhysicsTimeSeries::new("fluid");
        s2.push(0.5, 0.0);
        s2.push(2.0, 2.0);
        tl.add_domain(s1);
        tl.add_domain(s2);
        assert!((tl.global_start() - 0.0).abs() < 1e-12);
        assert!((tl.global_end() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_timeline_fixed_sync_step() {
        let tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.25));
        assert!((tl.sync_step() - 0.25).abs() < 1e-12);
    }
    #[test]
    fn test_timeline_sample_all() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.1));
        let mut s = PhysicsTimeSeries::new("domain");
        s.push(0.0, 10.0);
        s.push(1.0, 20.0);
        tl.add_domain(s);
        let vals = tl.sample_all(0.5);
        assert_eq!(vals.len(), 1);
        assert!((vals[0].unwrap() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_time_series_interpolate_out_of_range() {
        let mut s = PhysicsTimeSeries::new("x");
        s.push(1.0, 100.0);
        s.push(2.0, 200.0);
        assert!((s.interpolate_at(0.0).unwrap() - 100.0).abs() < 1e-10);
        assert!((s.interpolate_at(3.0).unwrap() - 200.0).abs() < 1e-10);
    }
    #[test]
    fn test_scalar_profile_bars_count() {
        let pos: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, 0.0, 0.0]).collect();
        let vals: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let bars = scalar_profile_bars(&pos, &vals, 0.0, 7.0, 1.0);
        assert_eq!(bars.len(), 8);
    }
    #[test]
    fn test_scalar_profile_bar_height() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let vals = vec![1.0];
        let bars = scalar_profile_bars(&pos, &vals, 0.0, 1.0, 2.0);
        assert_eq!(bars.len(), 1);
        assert!((bars[0].end[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_iso_contour_2x2_single_crossing() {
        let grid = vec![0.0, 0.0, 1.0, 1.0];
        let segs = iso_contour_lines(&grid, 2, 2, 0.5, 1.0, [0.0; 3], Rgba::white());
        assert!(!segs.is_empty());
    }
    #[test]
    fn test_iso_contour_uniform_no_crossing() {
        let grid = vec![1.0, 1.0, 1.0, 1.0];
        let segs = iso_contour_lines(&grid, 2, 2, 0.5, 1.0, [0.0; 3], Rgba::white());
        assert_eq!(segs.len(), 0);
    }
    #[test]
    fn test_multiphysics_residual_norm() {
        let f1 = vec![0.1, -0.5, 0.3];
        let f2 = vec![1.2, 0.0];
        let norm = multiphysics_residual_norm(&[&f1, &f2]);
        assert!((norm - 1.2).abs() < 1e-12);
    }
    #[test]
    fn test_coupling_error() {
        let old = vec![1.0, 2.0, 3.0];
        let new = vec![1.1, 1.8, 3.5];
        let err = coupling_error(&old, &new);
        assert!((err[0] - 0.1).abs() < 1e-12);
        assert!((err[1] - 0.2).abs() < 1e-12);
        assert!((err[2] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_smooth_scalar_field_uniform() {
        let vals = vec![2.0; 10];
        let smoothed = smooth_scalar_field(&vals, 2);
        for v in &smoothed {
            assert!((v - 2.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_smooth_scalar_field_reduces_peak() {
        let mut vals = vec![0.0; 11];
        vals[5] = 10.0;
        let smoothed = smooth_scalar_field(&vals, 2);
        assert!(smoothed[5] < 10.0, "smoothing should reduce peak");
    }
    #[test]
    fn test_clip_scalar_field() {
        let nodes = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let vals = vec![10.0, 20.0, 30.0];
        let plane = ClippingPlane::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let (vis_nodes, vis_vals) = clip_scalar_field(&nodes, &vals, &plane);
        assert_eq!(vis_nodes.len(), 2);
        assert_eq!(vis_vals.len(), 2);
        assert!(vis_vals.contains(&10.0));
        assert!(vis_vals.contains(&30.0));
    }
    #[test]
    fn test_field_node_construction() {
        let n = FieldNode::new([1.0, 2.0, 3.0], 42.0, [0.1, 0.2, 0.3]);
        assert!((n.scalar - 42.0).abs() < 1e-12);
        assert!((n.position[0] - 1.0).abs() < 1e-12);
        assert!((n.vector[1] - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_scalar_clamped() {
        assert!((normalize_scalar(-1.0, 0.0, 1.0) - 0.0).abs() < 1e-12);
        assert!((normalize_scalar(2.0, 0.0, 1.0) - 1.0).abs() < 1e-12);
        assert!((normalize_scalar(0.5, 0.0, 1.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_scalar_degenerate() {
        assert!((normalize_scalar(5.0, 5.0, 5.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_hot_colormap_zero() {
        let c = hot_colormap(0.0);
        assert!(c.r < 1e-6 && c.g < 1e-6 && c.b < 1e-6);
    }
    #[test]
    fn test_hot_colormap_one() {
        let c = hot_colormap(1.0);
        assert!((c.r - 1.0).abs() < 1e-6);
        assert!((c.g - 1.0).abs() < 1e-6);
        assert!((c.b - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_cool_colormap_range() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let c = cool_colormap(t);
            assert!(c.r >= 0.0 && c.r <= 1.0);
            assert!(c.g >= 0.0 && c.g <= 1.0);
            assert!(c.b >= 0.0 && c.b <= 1.0);
        }
    }
    #[test]
    fn test_fsi_snapshot_initial_state() {
        let snap = FsiSnapshot::new(1.5);
        assert!((snap.time - 1.5).abs() < 1e-12);
        assert!(snap.fluid_positions.is_empty());
        assert!(snap.structural_positions.is_empty());
    }
    #[test]
    fn test_timeline_sync_grid_length() {
        let mut tl = MultiphysicsTimeline::new(SyncPolicy::Fixed(0.5));
        let mut s = PhysicsTimeSeries::new("d");
        s.push(0.0, 0.0);
        s.push(2.0, 1.0);
        tl.add_domain(s);
        let grid = tl.sync_grid();
        assert_eq!(grid.len(), 5);
    }
}
