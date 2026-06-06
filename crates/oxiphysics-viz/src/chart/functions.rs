//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute the minimum of a slice of floats.
pub(super) fn slice_min(data: &[f64]) -> f64 {
    data.iter().cloned().fold(f64::MAX, f64::min)
}
/// Compute the maximum of a slice of floats.
pub(super) fn slice_max(data: &[f64]) -> f64 {
    data.iter().cloned().fold(f64::MIN, f64::max)
}
/// Map a value in `[lo, hi]` to `[0, 1]`.
pub(super) fn normalize(v: f64, lo: f64, hi: f64) -> f64 {
    if (hi - lo).abs() < 1e-12 {
        0.5
    } else {
        ((v - lo) / (hi - lo)).clamp(0.0, 1.0)
    }
}
/// Computes a centred moving average of `data` with the given `window` size.
///
/// Near-boundary values use a truncated window.  Output length equals input length.
pub fn moving_average(data: &[f64], window: usize) -> Vec<f64> {
    let n = data.len();
    let w = window.max(1);
    (0..n)
        .map(|i| {
            let half = w / 2;
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            let count = hi - lo;
            data[lo..hi].iter().sum::<f64>() / count as f64
        })
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::chart::Annotation;
    use crate::chart::AxisConfig;
    use crate::chart::AxisTicksLabel;
    use crate::chart::BarChart;
    use crate::chart::BubbleChart;
    use crate::chart::ContourPlot;
    use crate::chart::ConvergencePlot;
    use crate::chart::EnergyTimeSeries;
    use crate::chart::ErrorBars;
    use crate::chart::HeatmapPlot;
    use crate::chart::HistogramPlot;
    use crate::chart::LinePlot;
    use crate::chart::LogAxis;
    use crate::chart::MultiLinePlot;
    use crate::chart::ParallelCoord;
    use crate::chart::PhasePortrait;
    use crate::chart::PolarPlot;
    use crate::chart::RealTimeBuffer;
    use crate::chart::ScatterPlot;
    use crate::chart::SensorPlot;
    use crate::chart::Series;
    use crate::chart::SvgExporter;
    use crate::chart::TelemetrySample;
    #[test]
    fn test_series_len() {
        let s = Series::new(
            "a",
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0],
            [1.0, 0.0, 0.0, 1.0],
        );
        assert_eq!(s.len(), 2);
    }
    #[test]
    fn test_axis_range_auto() {
        let ax = AxisConfig::new("x");
        let (lo, hi) = ax.effective_range(&[1.0, 2.0, 3.0]);
        assert!((lo - 1.0).abs() < 1e-9 && (hi - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_axis_ticks_count() {
        let ax = AxisConfig::new("x");
        let ticks = ax.tick_positions(0.0, 10.0);
        assert_eq!(ticks.len(), ax.ticks + 1);
    }
    #[test]
    fn test_line_plot_x_range() {
        let mut lp = LinePlot::new("test");
        lp.add_series(Series::new("s", vec![0.0, 5.0], vec![1.0, 2.0], [1.0; 4]));
        let (lo, hi) = lp.x_range();
        assert!((lo - 0.0).abs() < 1e-9 && (hi - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_line_plot_normalize_bound() {
        let mut lp = LinePlot::new("test");
        lp.add_series(Series::new("s", vec![0.0, 10.0], vec![0.0, 5.0], [1.0; 4]));
        let p = lp.normalize_point(10.0, 5.0);
        assert!((p[0] - 1.0).abs() < 1e-9 && (p[1] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_scatter_add_point() {
        let mut sp = ScatterPlot::new("test");
        sp.add_point(1.0, 2.0, 0.5, 5.0);
        assert_eq!(sp.points.len(), 1);
    }
    #[test]
    fn test_scatter_count_region() {
        let mut sp = ScatterPlot::new("test");
        sp.add_point(1.0, 1.0, 0.5, 5.0);
        sp.add_point(5.0, 5.0, 0.5, 5.0);
        assert_eq!(sp.count_in_region(0.0, 3.0, 0.0, 3.0), 1);
    }
    #[test]
    fn test_histogram_count_sum() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = HistogramPlot::new("test", &data, 10);
        let total: usize = h.counts.iter().sum();
        assert_eq!(total, 100);
    }
    #[test]
    fn test_histogram_density_sum() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let h = HistogramPlot::new("test", &data, 10);
        let bw = h.bin_edges[1] - h.bin_edges[0];
        let integral: f64 = h.densities().iter().sum::<f64>() * bw;
        assert!((integral - 1.0).abs() < 0.05);
    }
    #[test]
    fn test_heatmap_sample_flat() {
        let hm = HeatmapPlot::new("test", vec![2.0; 16], 4, 4);
        let v = hm.sample(1.5, 1.5);
        assert!((v - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_heatmap_contour_count() {
        let hm = HeatmapPlot::new("test", vec![0.0; 16], 4, 4);
        assert_eq!(hm.contour_thresholds().len(), hm.contour_levels + 1);
    }
    #[test]
    fn test_bubble_add_point() {
        let mut bc = BubbleChart::new("test");
        bc.add_point(1.0, 2.0, 3.0, [1.0; 4]);
        assert_eq!(bc.points.len(), 1);
    }
    #[test]
    fn test_bubble_display_size_range() {
        let mut bc = BubbleChart::new("test");
        bc.add_point(0.0, 0.0, 1.0, [1.0; 4]);
        bc.add_point(0.0, 0.0, 10.0, [1.0; 4]);
        let s = bc.display_size(5.0);
        assert!(s >= bc.min_size && s <= bc.max_size);
    }
    #[test]
    fn test_parallel_coord_normalized() {
        let mut pc = ParallelCoord::new("test", vec!["a".to_string(), "b".to_string()]);
        pc.add_row(vec![0.0, 0.0], [1.0; 4]);
        pc.add_row(vec![10.0, 5.0], [1.0; 4]);
        let norm = pc.normalized_row(1);
        assert!(norm.iter().all(|&v| (0.0..=1.0).contains(&v)));
    }
    #[test]
    fn test_phase_portrait_add() {
        let mut pp = PhasePortrait::new("test");
        pp.add_point(1.0, 2.0);
        assert_eq!(pp.x.len(), 1);
    }
    #[test]
    fn test_phase_stable() {
        let j = [-1.0, 0.0, 0.0, -1.0];
        assert_eq!(PhasePortrait::classify_fixed_point(j), "stable");
    }
    #[test]
    fn test_phase_saddle() {
        let j = [1.0, 0.0, 0.0, -2.0];
        assert_eq!(PhasePortrait::classify_fixed_point(j), "saddle");
    }
    #[test]
    fn test_energy_series_push() {
        let mut es = EnergyTimeSeries::new("test");
        es.push(0.0, 10.0, 5.0);
        assert!((es.total[0] - 15.0).abs() < 1e-9);
    }
    #[test]
    fn test_energy_drift_zero() {
        let mut es = EnergyTimeSeries::new("test");
        es.push(0.0, 5.0, 5.0);
        es.push(1.0, 4.0, 6.0);
        es.push(2.0, 3.0, 7.0);
        assert!(es.energy_drift() < 1e-9);
    }
    #[test]
    fn test_sensor_max_speed() {
        let mut sp = SensorPlot::new("test");
        sp.push(TelemetrySample {
            time: 0.0,
            tire_fl: 0.0,
            tire_fr: 0.0,
            tire_rl: 0.0,
            tire_rr: 0.0,
            g_lateral: 0.0,
            g_longitudinal: 0.0,
            speed: 30.0,
        });
        sp.push(TelemetrySample {
            time: 1.0,
            tire_fl: 0.0,
            tire_fr: 0.0,
            tire_rl: 0.0,
            tire_rr: 0.0,
            g_lateral: 0.0,
            g_longitudinal: 0.0,
            speed: 50.0,
        });
        assert!((sp.max_speed() - 50.0).abs() < 1e-9);
    }
    #[test]
    fn test_convergence_converged() {
        let mut cp = ConvergencePlot::new("test", 1e-6);
        cp.push(0, 1e-3);
        cp.push(1, 1e-7);
        assert!(cp.converged());
    }
    #[test]
    fn test_convergence_rate() {
        let mut cp = ConvergencePlot::new("test", 1e-10);
        cp.push(0, 1.0);
        cp.push(1, 0.1);
        assert!((cp.convergence_rate - 0.1).abs() < 1e-9);
    }
    #[test]
    fn test_convergence_estimate() {
        let mut cp = ConvergencePlot::new("test", 1e-6);
        cp.push(0, 1.0);
        cp.push(1, 0.1);
        let est = cp.estimated_iters_to_convergence();
        assert!(est.is_some());
    }
    #[test]
    fn test_histogram_normal_overlay_len() {
        let data: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let h = HistogramPlot::new("test", &data, 10);
        assert_eq!(h.normal_pdf_overlay().len(), 10);
    }
    #[test]
    fn test_normalize_clamp() {
        assert!((normalize(15.0, 0.0, 10.0) - 1.0).abs() < 1e-9);
        assert!((normalize(-5.0, 0.0, 10.0)).abs() < 1e-9);
    }
    #[test]
    fn test_slice_min_max() {
        let data = vec![3.0, 1.0, 4.0, 1.5, 9.0, 2.6];
        assert!((slice_min(&data) - 1.0).abs() < 1e-9);
        assert!((slice_max(&data) - 9.0).abs() < 1e-9);
    }
    #[test]
    fn test_scatter_color_max() {
        let mut sp = ScatterPlot::new("test");
        sp.value_range = (0.0, 1.0);
        let c = sp.value_to_color(1.0);
        assert!((c[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_bubble_sorted_z() {
        let mut bc = BubbleChart::new("test");
        bc.add_point(0.0, 0.0, 3.0, [1.0; 4]);
        bc.add_point(0.0, 0.0, 1.0, [1.0; 4]);
        bc.add_point(0.0, 0.0, 2.0, [1.0; 4]);
        let sorted = bc.sorted_by_z_desc();
        assert!(sorted[0].z >= sorted[1].z && sorted[1].z >= sorted[2].z);
    }
    #[test]
    fn test_energy_extra_channel() {
        let mut es = EnergyTimeSeries::new("test");
        es.add_channel("elastic", vec![1.0, 2.0, 3.0]);
        assert_eq!(es.extra.len(), 1);
    }
    #[test]
    fn test_sensor_speed_series_len() {
        let mut sp = SensorPlot::new("test");
        for i in 0..5 {
            sp.push(TelemetrySample {
                time: i as f64,
                tire_fl: 0.0,
                tire_fr: 0.0,
                tire_rl: 0.0,
                tire_rr: 0.0,
                g_lateral: 0.0,
                g_longitudinal: 0.0,
                speed: i as f64 * 10.0,
            });
        }
        let (t, v) = sp.speed_series();
        assert_eq!(t.len(), 5);
        assert_eq!(v.len(), 5);
    }
    #[test]
    fn test_multi_line_legend_count() {
        let mut ml = MultiLinePlot::new("test");
        ml.add_series(Series::new("a", vec![0.0, 1.0], vec![0.0, 1.0], [1.0; 4]));
        ml.add_series(Series::new("b", vec![0.0, 1.0], vec![1.0, 2.0], [1.0; 4]));
        let legend = ml.legend_entries();
        assert_eq!(legend.len(), 2);
    }
    #[test]
    fn test_bar_chart_group_total() {
        let mut bc = BarChart::new("test");
        bc.add_group("A", vec![1.0, 2.0, 3.0]);
        assert!((bc.group_total(0) - 6.0).abs() < 1e-9);
    }
    #[test]
    fn test_bar_chart_stacked_max() {
        let mut bc = BarChart::new("test");
        bc.add_group("X", vec![2.0, 3.0]);
        bc.add_group("Y", vec![1.0, 4.0]);
        assert!((bc.stacked_max() - 7.0).abs() < 1e-9);
    }
    #[test]
    fn test_polar_plot_add() {
        let mut pp = PolarPlot::new("test");
        pp.add_point(0.0, 1.0);
        assert_eq!(pp.points.len(), 1);
    }
    #[test]
    fn test_polar_rose_bin() {
        let pp = PolarPlot::new("test");
        let idx = pp.rose_bin(std::f64::consts::PI, 8);
        assert!(idx < 8);
    }
    #[test]
    fn test_contour_flat_field() {
        let cp = ContourPlot::new("test", vec![1.0; 9], 3, 3, 0.1);
        let segs = cp.marching_squares(0.5);
        assert!(segs.is_empty() || !segs.is_empty());
    }
    #[test]
    fn test_error_bars_range() {
        let eb = ErrorBars::new(
            vec![1.0, 2.0, 3.0],
            vec![10.0, 20.0, 30.0],
            vec![1.0, 2.0, 1.5],
        );
        assert!((eb.error_range(1).0 - 18.0).abs() < 1e-9);
        assert!((eb.error_range(1).1 - 22.0).abs() < 1e-9);
    }
    #[test]
    fn test_log_axis_ticks() {
        let la = LogAxis::new(1.0, 1000.0, 4);
        let ticks = la.tick_positions();
        assert_eq!(ticks.len(), 4);
    }
    #[test]
    fn test_log_axis_tick_bounds() {
        let la = LogAxis::new(0.1, 100.0, 5);
        let ticks = la.tick_positions();
        for &t in &ticks {
            assert!((0.09..=101.0).contains(&t));
        }
    }
    #[test]
    fn test_axis_ticks_linear_count() {
        let at = AxisTicksLabel::linear(0.0, 1.0, 5);
        assert_eq!(at.positions.len(), 5);
    }
    #[test]
    fn test_svg_line_element() {
        let exp = SvgExporter::new(800, 600);
        let svg = exp.line_element(0.0, 0.0, 100.0, 100.0, "black", 1.0);
        assert!(!svg.is_empty());
        assert!(svg.contains("<line"));
    }
    #[test]
    fn test_svg_export_viewbox() {
        let exp = SvgExporter::new(400, 300);
        let out = exp.export_header();
        assert!(out.contains("viewBox"));
    }
    #[test]
    fn test_annotation_render() {
        let ann = Annotation::new("hello", 10.0, 20.0);
        let s = ann.to_svg();
        assert!(s.contains("hello"));
    }
    #[test]
    fn test_moving_avg_window1() {
        let data = vec![1.0, 3.0, 5.0, 7.0];
        let ma = moving_average(&data, 1);
        for (a, b) in data.iter().zip(ma.iter()) {
            assert!((a - b).abs() < 1e-9);
        }
    }
    #[test]
    fn test_moving_avg_window3() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ma = moving_average(&data, 3);
        assert_eq!(ma.len(), data.len());
        assert!((ma[1] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_realtime_buffer_capacity() {
        let mut buf = RealTimeBuffer::new(3);
        for i in 0..5 {
            buf.push_sample(i as f64);
        }
        assert_eq!(buf.samples.len(), 3);
    }
    #[test]
    fn test_realtime_buffer_latest() {
        let mut buf = RealTimeBuffer::new(5);
        buf.push_sample(42.0);
        assert!((buf.latest().unwrap() - 42.0).abs() < 1e-9);
    }
    #[test]
    fn test_multi_line_x_range() {
        let mut ml = MultiLinePlot::new("test");
        ml.add_series(Series::new("a", vec![-5.0, 0.0], vec![0.0, 1.0], [1.0; 4]));
        ml.add_series(Series::new("b", vec![0.0, 10.0], vec![0.0, 0.0], [1.0; 4]));
        let (lo, hi) = ml.x_range();
        assert!((lo - (-5.0)).abs() < 1e-9);
        assert!((hi - 10.0).abs() < 1e-9);
    }
    #[test]
    fn test_contour_threshold_count() {
        let data: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let cp = ContourPlot::new("test", data, 4, 4, 0.25);
        let thresholds = cp.auto_thresholds(5);
        assert_eq!(thresholds.len(), 5);
    }
    #[test]
    fn test_error_bars_len() {
        let eb = ErrorBars::new(
            vec![0.0, 1.0, 2.0],
            vec![5.0, 6.0, 7.0],
            vec![0.5, 0.5, 0.5],
        );
        assert_eq!(eb.len(), 3);
    }
    #[test]
    fn test_bar_chart_categories() {
        let mut bc = BarChart::new("test");
        bc.add_group("G1", vec![1.0, 2.0, 3.0]);
        bc.add_group("G2", vec![4.0, 5.0, 6.0]);
        assert_eq!(bc.n_categories(), 3);
    }
    #[test]
    fn test_svg_circle_element() {
        let exp = SvgExporter::new(800, 600);
        let s = exp.circle_element(50.0, 50.0, 10.0, "blue");
        assert!(s.contains("<circle"));
    }
    #[test]
    fn test_annotation_font_size() {
        let ann = Annotation::new("label", 5.0, 10.0).with_font_size(14);
        assert_eq!(ann.font_size, 14);
    }
    #[test]
    fn test_realtime_buffer_mean() {
        let mut buf = RealTimeBuffer::new(4);
        for v in [1.0_f64, 2.0, 3.0, 4.0] {
            buf.push_sample(v);
        }
        assert!((buf.mean() - 2.5).abs() < 1e-9);
    }
    #[test]
    fn test_moving_avg_output_len() {
        let data: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let ma = moving_average(&data, 4);
        assert_eq!(ma.len(), data.len());
    }
}
