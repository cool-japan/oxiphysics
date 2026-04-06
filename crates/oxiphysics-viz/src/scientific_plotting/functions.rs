//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{PlotColor, PlotColormap};

/// Generate evenly-spaced values in `[start, end]` (inclusive).
pub fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * i as f64 / (n - 1) as f64)
        .collect()
}
/// Apply a colormap to an array of scalar values, returning colors.
pub fn colorize(values: &[f64], colormap: PlotColormap) -> Vec<PlotColor> {
    let mn = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let mx = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let range = mx - mn;
    values
        .iter()
        .map(|&v| {
            let t = if range.abs() < 1e-14 {
                0.5
            } else {
                ((v - mn) / range) as f32
            };
            colormap.map(t.clamp(0.0, 1.0))
        })
        .collect()
}
/// Compute a 2D histogram from `(x, y)` data.
///
/// Returns a flat row-major array of counts of size `bins_x × bins_y`.
pub fn histogram2d(
    x: &[f64],
    y: &[f64],
    bins_x: usize,
    bins_y: usize,
    x_range: [f64; 2],
    y_range: [f64; 2],
) -> Vec<u64> {
    let mut counts = vec![0u64; bins_x * bins_y];
    let dx = (x_range[1] - x_range[0]) / bins_x as f64;
    let dy = (y_range[1] - y_range[0]) / bins_y as f64;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        if xi < x_range[0] || xi >= x_range[1] {
            continue;
        }
        if yi < y_range[0] || yi >= y_range[1] {
            continue;
        }
        let bx = ((xi - x_range[0]) / dx) as usize;
        let by = ((yi - y_range[0]) / dy) as usize;
        let bx = bx.min(bins_x - 1);
        let by = by.min(bins_y - 1);
        counts[by * bins_x + bx] += 1;
    }
    counts
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::scientific_plotting::AxisConfig;
    use crate::scientific_plotting::ContourPlot;
    use crate::scientific_plotting::ExportFormat;
    use crate::scientific_plotting::ExportOptions;
    use crate::scientific_plotting::LinePlot;
    use crate::scientific_plotting::LineSeries;
    use crate::scientific_plotting::LineStyle2D;
    use crate::scientific_plotting::PhasePortrait;
    use crate::scientific_plotting::PhaseTrajectory;
    use crate::scientific_plotting::PlotLayout;
    use crate::scientific_plotting::PlotType;
    use crate::scientific_plotting::ScatterPlot;
    use crate::scientific_plotting::ScatterPoint;
    use crate::scientific_plotting::SurfaceLighting;
    use crate::scientific_plotting::SurfacePlot3D;
    use crate::scientific_plotting::SvgExporter;
    #[test]
    fn test_color_rgb_components() {
        let c = PlotColor::rgb(0.1, 0.2, 0.3);
        assert!((c.r - 0.1).abs() < 1e-6);
        assert!((c.g - 0.2).abs() < 1e-6);
        assert!((c.b - 0.3).abs() < 1e-6);
        assert!((c.a - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_color_rgba_alpha() {
        let c = PlotColor::rgba(1.0, 0.0, 0.0, 0.5);
        assert!((c.a - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_color_lerp_midpoint() {
        let c = PlotColor::lerp(PlotColor::black(), PlotColor::white(), 0.5);
        assert!((c.r - 0.5).abs() < 1e-6);
        assert!((c.g - 0.5).abs() < 1e-6);
        assert!((c.b - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_color_lerp_t0_is_a() {
        let a = PlotColor::red();
        let b = PlotColor::blue();
        let c = PlotColor::lerp(a, b, 0.0);
        assert!((c.r - 1.0).abs() < 1e-6);
        assert!(c.b.abs() < 1e-6);
    }
    #[test]
    fn test_color_lerp_t1_is_b() {
        let a = PlotColor::red();
        let b = PlotColor::blue();
        let c = PlotColor::lerp(a, b, 1.0);
        assert!(c.r.abs() < 1e-6);
        assert!((c.b - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_colormap_jet_t0_is_blue() {
        let c = PlotColormap::Jet.map(0.0);
        assert!(c.b > c.r, "jet at t=0 should be bluish");
    }
    #[test]
    fn test_colormap_jet_t1_is_red() {
        let c = PlotColormap::Jet.map(1.0);
        assert!(c.r > c.b, "jet at t=1 should be reddish");
    }
    #[test]
    fn test_colormap_greys_monotone() {
        let c0 = PlotColormap::Greys.map(0.0);
        let c1 = PlotColormap::Greys.map(1.0);
        assert!(c0.r < c1.r, "greys should go dark to light");
    }
    #[test]
    fn test_colormap_bwr_midpoint_is_white() {
        let c = PlotColormap::BlueWhiteRed.map(0.5);
        assert!(c.r > 0.9, "BWR midpoint r≈1");
        assert!(c.g > 0.9, "BWR midpoint g≈1");
        assert!(c.b > 0.9, "BWR midpoint b≈1");
    }
    #[test]
    fn test_colormap_clamp_below_zero() {
        let c = PlotColormap::Jet.map(-1.0);
        assert!(c.r >= 0.0 && c.r <= 1.0);
    }
    #[test]
    fn test_colormap_clamp_above_one() {
        let c = PlotColormap::Jet.map(2.0);
        assert!(c.r >= 0.0 && c.r <= 1.0);
    }
    #[test]
    fn test_colorize_array() {
        let vals = vec![0.0, 0.5, 1.0];
        let colors = colorize(&vals, PlotColormap::Greys);
        assert_eq!(colors.len(), 3);
        assert!(colors[0].r < colors[2].r);
    }
    #[test]
    fn test_axis_effective_range_auto() {
        let ax = AxisConfig::default();
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let (lo, hi) = ax.effective_range(&data);
        assert!((lo - 1.0).abs() < 1e-10);
        assert!((hi - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_axis_effective_range_fixed() {
        let ax = AxisConfig::default().with_range(-5.0, 5.0);
        let (lo, hi) = ax.effective_range(&[1.0, 2.0]);
        assert!((lo + 5.0).abs() < 1e-10);
        assert!((hi - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_axis_effective_range_single_value_expanded() {
        let ax = AxisConfig::default();
        let (lo, hi) = ax.effective_range(&[3.0]);
        assert!(hi > lo, "single value should be expanded");
    }
    #[test]
    fn test_line_series_len() {
        let s = LineSeries::new(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 4.0],
            LineStyle2D::default(),
        );
        assert_eq!(s.len(), 3);
    }
    #[test]
    fn test_line_series_y_range() {
        let s = LineSeries::new(vec![0.0, 1.0], vec![-3.0, 7.0], LineStyle2D::default());
        let (lo, hi) = s.y_range();
        assert!((lo + 3.0).abs() < 1e-10);
        assert!((hi - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_series_interpolate_at() {
        let s = LineSeries::new(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 4.0],
            LineStyle2D::default(),
        );
        let y = s.interpolate_at(0.5);
        assert!((y - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_line_series_interpolate_clamp_low() {
        let s = LineSeries::new(vec![1.0, 2.0], vec![5.0, 10.0], LineStyle2D::default());
        assert!((s.interpolate_at(0.0) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_series_interpolate_clamp_high() {
        let s = LineSeries::new(vec![1.0, 2.0], vec![5.0, 10.0], LineStyle2D::default());
        assert!((s.interpolate_at(3.0) - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_series_resample() {
        let s = LineSeries::new(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 0.0],
            LineStyle2D::default(),
        );
        let r = s.resample(5);
        assert_eq!(r.len(), 5);
    }
    #[test]
    fn test_line_plot_empty() {
        let p = LinePlot::new("test");
        assert_eq!(p.num_series(), 0);
    }
    #[test]
    fn test_line_plot_add_series() {
        let mut p = LinePlot::new("test");
        let s = LineSeries::new(vec![0.0, 1.0], vec![0.0, 1.0], LineStyle2D::default());
        p.add_series(s);
        assert_eq!(p.num_series(), 1);
    }
    #[test]
    fn test_line_plot_auto_range() {
        let mut p = LinePlot::new("test");
        p.add_xy(vec![0.0, 1.0], vec![-1.0, 1.0], PlotColor::blue(), "a");
        let (xlo, xhi) = p.auto_x_range();
        assert!((xlo - 0.0).abs() < 1e-10);
        assert!((xhi - 1.0).abs() < 1e-10);
        let (ylo, yhi) = p.auto_y_range();
        assert!((ylo + 1.0).abs() < 1e-10);
        assert!((yhi - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_plot_set_labels() {
        let mut p = LinePlot::new("test");
        p.set_labels("Time (s)", "Displacement (m)");
        assert_eq!(p.x_axis.label, "Time (s)");
        assert_eq!(p.y_axis.label, "Displacement (m)");
    }
    #[test]
    fn test_line_plot_multiple_series() {
        let mut p = LinePlot::new("multi");
        for i in 0..4 {
            p.add_xy(
                vec![0.0, 1.0],
                vec![i as f64, i as f64 + 1.0],
                PlotColor::blue(),
                "s",
            );
        }
        assert_eq!(p.num_series(), 4);
    }
    #[test]
    fn test_scatter_empty() {
        let s = ScatterPlot::new("test");
        assert!(s.is_empty());
    }
    #[test]
    fn test_scatter_add_point() {
        let mut s = ScatterPlot::new("test");
        s.add_point(ScatterPoint::new(1.0, 2.0));
        assert_eq!(s.len(), 1);
    }
    #[test]
    fn test_scatter_add_xy() {
        let mut s = ScatterPlot::new("test");
        s.add_xy(&[0.0, 1.0, 2.0], &[3.0, 4.0, 5.0]);
        assert_eq!(s.len(), 3);
    }
    #[test]
    fn test_scatter_color_for_no_colormap() {
        let mut s = ScatterPlot::new("test");
        s.add_point(ScatterPoint::new(0.0, 0.0));
        s.default_color = PlotColor::red();
        let c = s.color_for(0, 0.0, 1.0);
        assert!((c.r - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_scatter_color_for_with_colormap() {
        let mut s = ScatterPlot::new("test");
        s.use_colormap = true;
        s.colormap = PlotColormap::Greys;
        s.add_point(ScatterPoint::new(0.0, 0.0).with_color_value(1.0));
        s.add_point(ScatterPoint::new(1.0, 0.0).with_color_value(0.0));
        let (mn, mx) = s.auto_color_range();
        let c_bright = s.color_for(0, mn, mx);
        let c_dark = s.color_for(1, mn, mx);
        assert!(c_bright.r > c_dark.r);
    }
    #[test]
    fn test_scatter_x_range() {
        let mut s = ScatterPlot::new("test");
        s.add_xy(&[-2.0, 0.0, 3.0], &[0.0, 0.0, 0.0]);
        let (lo, hi) = s.x_range();
        assert!((lo + 2.0).abs() < 1e-10);
        assert!((hi - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_scatter_auto_color_range() {
        let mut s = ScatterPlot::new("test");
        s.add_point(ScatterPoint::new(0.0, 0.0).with_color_value(10.0));
        s.add_point(ScatterPoint::new(1.0, 0.0).with_color_value(20.0));
        let (mn, mx) = s.auto_color_range();
        assert!((mn - 10.0).abs() < 1e-10);
        assert!((mx - 20.0).abs() < 1e-10);
    }
    #[test]
    fn test_contour_from_function() {
        let cp = ContourPlot::from_function(
            |x, y| x * x + y * y,
            10,
            10,
            [-1.0, 1.0],
            [-1.0, 1.0],
            "circles",
        );
        assert_eq!(cp.values.len(), 100);
    }
    #[test]
    fn test_contour_data_range() {
        let cp = ContourPlot::from_function(|x, y| x + y, 5, 5, [0.0, 1.0], [0.0, 1.0], "linear");
        let (mn, mx) = cp.data_range();
        assert!(mn >= 0.0 - 1e-10);
        assert!(mx <= 2.0 + 1e-10);
    }
    #[test]
    fn test_contour_auto_levels() {
        let mut cp = ContourPlot::from_function(|x, _y| x, 5, 5, [0.0, 1.0], [0.0, 1.0], "test");
        cp.auto_levels(4);
        assert_eq!(cp.levels.len(), 4);
    }
    #[test]
    fn test_contour_extract_level_no_panic() {
        let mut cp = ContourPlot::from_function(
            |x, y| x * x + y * y,
            10,
            10,
            [-1.0, 1.0],
            [-1.0, 1.0],
            "test",
        );
        cp.auto_levels(3);
        for &l in &cp.levels.clone() {
            let level = cp.extract_level(l);
            assert_eq!(level.value, l);
        }
    }
    #[test]
    fn test_contour_value_at() {
        let cp = ContourPlot::from_function(|x, _y| x, 5, 5, [0.0, 4.0], [0.0, 4.0], "test");
        assert!(cp.value_at(0, 0).abs() < 1e-10);
    }
    #[test]
    fn test_contour_extract_all_levels() {
        let mut cp = ContourPlot::from_function(
            |x, y| x * x + y * y,
            8,
            8,
            [-1.0, 1.0],
            [-1.0, 1.0],
            "circles",
        );
        cp.auto_levels(5);
        let levels = cp.extract_all_levels();
        assert_eq!(levels.len(), 5);
    }
    #[test]
    fn test_surface_from_function_vertex_count() {
        let s = SurfacePlot3D::from_function(
            |x, y| x * x + y * y,
            5,
            5,
            [-1.0, 1.0],
            [-1.0, 1.0],
            "paraboloid",
        );
        assert_eq!(s.num_vertices(), 25);
    }
    #[test]
    fn test_surface_triangle_count() {
        let s = SurfacePlot3D::from_function(|x, y| x + y, 4, 4, [0.0, 1.0], [0.0, 1.0], "plane");
        assert_eq!(s.num_triangles(), 2 * 3 * 3);
    }
    #[test]
    fn test_surface_z_range() {
        let s = SurfacePlot3D::from_function(|x, _y| x, 5, 5, [0.0, 4.0], [0.0, 4.0], "ramp");
        let (mn, mx) = s.z_range();
        assert!(mn >= -1e-10);
        assert!((mx - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_surface_color_for_vertex() {
        let s = SurfacePlot3D::from_function(|x, _y| x, 5, 5, [0.0, 1.0], [0.0, 1.0], "test");
        let (mn, mx) = s.z_range();
        let c0 = s.color_for_vertex(0, mn, mx);
        assert!(c0.r >= 0.0 && c0.r <= 1.0);
    }
    #[test]
    fn test_surface_lighting_apply() {
        let light = SurfaceLighting::default();
        let color = PlotColor::rgb(1.0, 1.0, 1.0);
        let n = [0.0_f32, 0.0, 1.0];
        let lit = light.apply(color, n);
        assert!(lit.r >= 0.0 && lit.r <= 1.0);
    }
    #[test]
    fn test_phase_portrait_empty() {
        let pp = PhasePortrait::new("test");
        assert_eq!(pp.num_trajectories(), 0);
    }
    #[test]
    fn test_phase_trajectory_from_xy() {
        let t =
            PhaseTrajectory::from_xy(vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 0.0], PlotColor::blue());
        assert_eq!(t.len(), 3);
        assert_eq!(t.dim(), 2);
    }
    #[test]
    fn test_phase_trajectory_from_xyz() {
        let t = PhaseTrajectory::from_xyz(
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            PlotColor::red(),
        );
        assert_eq!(t.len(), 2);
        assert_eq!(t.dim(), 3);
    }
    #[test]
    fn test_phase_trajectory_arc_length() {
        let t =
            PhaseTrajectory::from_xy(vec![0.0, 1.0, 2.0], vec![0.0, 0.0, 0.0], PlotColor::blue());
        assert!((t.arc_length() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_phase_portrait_integrate_2d_uniform() {
        let t =
            PhasePortrait::integrate_2d(|_x, _y| (1.0, 0.0), 0.0, 0.0, 0.1, 5, PlotColor::blue());
        assert_eq!(t.len(), 6);
        assert!((t.points.last().unwrap()[0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_phase_portrait_integrate_3d() {
        let t = PhasePortrait::integrate_3d(
            |_x, _y, _z| (1.0, 0.0, 0.0),
            0.0,
            0.0,
            0.0,
            0.1,
            3,
            PlotColor::blue(),
        );
        assert_eq!(t.len(), 4);
        assert_eq!(t.dim(), 3);
    }
    #[test]
    fn test_phase_portrait_add_trajectory() {
        let mut pp = PhasePortrait::new("test");
        let t = PhaseTrajectory::from_xy(vec![0.0, 1.0], vec![0.0, 1.0], PlotColor::blue());
        pp.add_trajectory(t);
        assert_eq!(pp.num_trajectories(), 1);
    }
    #[test]
    fn test_phase_compute_nullcline() {
        let nc =
            PhasePortrait::compute_nullcline_x(|x, _y| (x - 1.0, 0.0), [0.0, 2.0], [0.0, 1.0], 20);
        assert!(!nc.points.is_empty(), "should find nullcline points");
    }
    #[test]
    fn test_layout_empty() {
        let layout = PlotLayout::new(2, 2, "test");
        assert_eq!(layout.num_panels(), 0);
    }
    #[test]
    fn test_layout_add_panels() {
        let mut layout = PlotLayout::new(2, 2, "test");
        layout.add_panel(0, 0, PlotType::Line);
        layout.add_panel(0, 1, PlotType::Scatter);
        assert_eq!(layout.num_panels(), 2);
    }
    #[test]
    fn test_layout_fill_empty() {
        let mut layout = PlotLayout::new(2, 2, "test");
        layout.add_panel(0, 0, PlotType::Line);
        layout.fill_empty();
        assert_eq!(layout.num_panels(), 4);
    }
    #[test]
    fn test_layout_panel_bounds() {
        let layout = PlotLayout::new(1, 2, "test");
        let b0 = layout.panel_bounds(0, 0);
        let b1 = layout.panel_bounds(0, 1);
        assert!(b1[0] > b0[0], "second panel should be to the right");
    }
    #[test]
    fn test_layout_axis_range() {
        let mut layout = PlotLayout::new(2, 2, "test");
        layout.set_axis_range(0, 0, (-1.0, 1.0), (-2.0, 2.0));
        let ((x0, x1), (y0, y1)) = layout.axis_ranges[&(0, 0)];
        assert!((x0 + 1.0).abs() < 1e-10);
        assert!((x1 - 1.0).abs() < 1e-10);
        assert!((y0 + 2.0).abs() < 1e-10);
        assert!((y1 - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_layout_describe_contains_title() {
        let layout = PlotLayout::new(2, 3, "MyFigure");
        let desc = layout.describe();
        assert!(desc.contains("MyFigure"));
    }
    #[test]
    fn test_layout_share_axes() {
        let mut layout = PlotLayout::new(2, 2, "test");
        layout.add_panel(1, 0, PlotType::Line);
        layout.share_x_axis(1, 0);
        assert!(layout.panels[0].shared_x);
    }
    #[test]
    fn test_svg_exporter_contains_svg_tag() {
        let mut p = LinePlot::new("SVG Test");
        p.add_xy(
            vec![0.0, 1.0, 2.0],
            vec![0.0, 1.0, 0.0],
            PlotColor::blue(),
            "wave",
        );
        let svg = SvgExporter::export_line_plot(&p);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }
    #[test]
    fn test_svg_exporter_contains_polyline() {
        let mut p = LinePlot::new("test");
        p.add_xy(vec![0.0, 1.0], vec![0.0, 1.0], PlotColor::red(), "line");
        let svg = SvgExporter::export_line_plot(&p);
        assert!(svg.contains("polyline"));
    }
    #[test]
    fn test_svg_exporter_title() {
        let p = LinePlot::new("My Plot");
        let svg = SvgExporter::export_line_plot(&p);
        assert!(svg.contains("My Plot"));
    }
    #[test]
    fn test_svg_exporter_empty_plot() {
        let p = LinePlot::new("empty");
        let svg = SvgExporter::export_line_plot(&p);
        assert!(svg.contains("<svg"));
        assert!(!svg.contains("polyline"));
    }
    #[test]
    fn test_linspace_count() {
        let v = linspace(0.0, 1.0, 11);
        assert_eq!(v.len(), 11);
    }
    #[test]
    fn test_linspace_endpoints() {
        let v = linspace(2.0, 5.0, 4);
        assert!((v[0] - 2.0).abs() < 1e-10);
        assert!((v[3] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_linspace_zero() {
        let v = linspace(0.0, 1.0, 0);
        assert!(v.is_empty());
    }
    #[test]
    fn test_linspace_one() {
        let v = linspace(3.0, 7.0, 1);
        assert_eq!(v.len(), 1);
        assert!((v[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_histogram2d_single_point() {
        let counts = histogram2d(&[0.5], &[0.5], 4, 4, [0.0, 1.0], [0.0, 1.0]);
        assert_eq!(counts.len(), 16);
        assert_eq!(counts.iter().sum::<u64>(), 1);
    }
    #[test]
    fn test_histogram2d_out_of_range_ignored() {
        let counts = histogram2d(&[2.0], &[2.0], 4, 4, [0.0, 1.0], [0.0, 1.0]);
        assert_eq!(counts.iter().sum::<u64>(), 0);
    }
    #[test]
    fn test_histogram2d_multiple_points() {
        let x = vec![0.1, 0.6, 0.1];
        let y = vec![0.1, 0.6, 0.9];
        let counts = histogram2d(&x, &y, 2, 2, [0.0, 1.0], [0.0, 1.0]);
        assert_eq!(counts.iter().sum::<u64>(), 3);
    }
    #[test]
    fn test_export_options_default() {
        let opts = ExportOptions::default();
        assert_eq!(opts.dpi, 300);
        assert_eq!(opts.format, ExportFormat::Svg);
        assert!(opts.tight_layout);
    }
    #[test]
    fn test_export_format_equality() {
        assert_eq!(ExportFormat::Png, ExportFormat::Png);
        assert_ne!(ExportFormat::Svg, ExportFormat::Pdf);
    }
}
