//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::types::*;
#[allow(unused_imports)]
use super::types_core::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sim_dashboard_total_energy() {
        let mut d = SimDashboard::new();
        d.total_ke = 10.0;
        d.total_pe = 5.0;
        assert!((d.total_energy() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_sim_dashboard_step() {
        let mut d = SimDashboard::new();
        d.step(0.01);
        assert_eq!(d.step_count, 1);
        assert!((d.sim_time - 0.01).abs() < 1e-10);
    }
    #[test]
    fn test_sim_dashboard_paused() {
        let mut d = SimDashboard::new();
        d.paused = true;
        d.step(0.01);
        assert!(d.sim_time.abs() < 1e-10);
        assert_eq!(d.step_count, 0);
    }
    #[test]
    fn test_sim_dashboard_status_string() {
        let mut d = SimDashboard::new();
        d.fps = 60.0;
        let s = d.status_string();
        assert!(s.contains("fps=60.0"), "s={s}");
    }
    #[test]
    fn test_energy_panel_mean_ke() {
        let mut ep = EnergyPanel::new(100, 0.5);
        ep.add_sample(EnergySample {
            time: 0.0,
            ke: 10.0,
            pe: 5.0,
        });
        ep.add_sample(EnergySample {
            time: 0.01,
            ke: 20.0,
            pe: 5.0,
        });
        assert!((ep.mean_ke() - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_energy_panel_anomaly() {
        let mut ep = EnergyPanel::new(100, 0.1);
        ep.add_sample(EnergySample {
            time: 0.0,
            ke: 10.0,
            pe: 0.0,
        });
        ep.add_sample(EnergySample {
            time: 0.01,
            ke: 100.0,
            pe: 0.0,
        });
        assert!(ep.anomaly_count > 0);
    }
    #[test]
    fn test_contact_panel_mean_count() {
        let mut cp = ContactPanel::new(100, 10, 1000.0);
        cp.add_sample(ContactSample {
            time: 0.0,
            count: 4,
            mean_force: 10.0,
            max_force: 20.0,
        });
        cp.add_sample(ContactSample {
            time: 0.01,
            count: 6,
            mean_force: 15.0,
            max_force: 30.0,
        });
        assert!((cp.mean_contact_count() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_panel_histogram() {
        let mut cp = ContactPanel::new(100, 10, 100.0);
        cp.add_sample(ContactSample {
            time: 0.0,
            count: 2,
            mean_force: 50.0,
            max_force: 80.0,
        });
        let total: u32 = cp.histogram.iter().sum();
        assert!(total > 0);
    }
    #[test]
    fn test_performance_panel_avg_ms() {
        let mut pp = PerformancePanel::new(100);
        pp.record_frame(FrameTiming {
            broadphase_ms: 1.0,
            narrowphase_ms: 2.0,
            solver_ms: 5.0,
            integration_ms: 1.0,
            memory_bytes: 1024,
        });
        assert!((pp.avg_total_ms() - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_performance_panel_avg_fps() {
        let mut pp = PerformancePanel::new(100);
        pp.record_frame(FrameTiming {
            broadphase_ms: 0.0,
            narrowphase_ms: 0.0,
            solver_ms: 16.667,
            integration_ms: 0.0,
            memory_bytes: 0,
        });
        let fps = pp.avg_fps();
        assert!(fps > 50.0 && fps < 70.0, "fps={fps}");
    }
    #[test]
    fn test_performance_panel_phase_fractions_sum() {
        let mut pp = PerformancePanel::new(10);
        pp.record_frame(FrameTiming {
            broadphase_ms: 1.0,
            narrowphase_ms: 2.0,
            solver_ms: 4.0,
            integration_ms: 1.0,
            memory_bytes: 0,
        });
        let fracs = pp.phase_fractions();
        let sum: f64 = fracs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "sum={sum}");
    }
    #[test]
    fn test_body_inspector_ke() {
        let bs = BodyState {
            index: 0,
            mass: 2.0,
            inertia: [1.0; 3],
            velocity: [3.0, 4.0, 0.0],
            angular_velocity: [0.0; 3],
            net_force: [0.0; 3],
            net_torque: [0.0; 3],
            is_sleeping: false,
        };
        assert!((bs.kinetic_energy() - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_body_inspector_selected_ke() {
        let mut bi = BodyInspector::new();
        assert!(bi.selected_ke().abs() < 1e-10);
        let bs = BodyState {
            index: 0,
            mass: 1.0,
            inertia: [0.0; 3],
            velocity: [2.0, 0.0, 0.0],
            angular_velocity: [0.0; 3],
            net_force: [0.0; 3],
            net_torque: [0.0; 3],
            is_sleeping: false,
        };
        bi.select(bs);
        assert!((bi.selected_ke() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_constraint_state_normalized_position() {
        let cs = ConstraintState {
            index: 0,
            joint_type: JointType::Revolute,
            current_value: 0.5,
            current_force: 10.0,
            limit_lo: 0.0,
            limit_hi: 1.0,
            motor_speed: 0.0,
            motor_active: false,
        };
        assert!((cs.normalized_position() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_constraint_state_at_limits() {
        let mut cs = ConstraintState {
            index: 0,
            joint_type: JointType::Revolute,
            current_value: 0.0,
            current_force: 0.0,
            limit_lo: 0.0,
            limit_hi: 1.0,
            motor_speed: 0.0,
            motor_active: false,
        };
        assert!(cs.at_lower_limit());
        cs.current_value = 1.0;
        assert!(cs.at_upper_limit());
    }
    #[test]
    fn test_material_state_stress_at_strain() {
        let mut ms = MaterialState::new("steel", 200e9, 0.3, 7800.0);
        ms.generate_linear_curve(0.01, 11);
        let s = ms.stress_at_strain(0.005);
        assert!((s - 1e9).abs() < 1e6, "s={s}");
    }
    #[test]
    fn test_material_state_linear_curve() {
        let mut ms = MaterialState::new("rubber", 1e6, 0.49, 1100.0);
        ms.generate_linear_curve(0.5, 6);
        assert_eq!(ms.stress_strain_curve.len(), 6);
        assert!(ms.stress_strain_curve[0].stress.abs() < 1e-6);
    }
    #[test]
    fn test_parameter_sweep_scatter_pair() {
        let mut ps = ParameterSweep::new(vec!["a".into(), "b".into()], vec!["m".into()]);
        ps.add_result(SweepResult {
            parameters: vec![1.0, 2.0],
            metrics: vec![3.0],
        });
        ps.add_result(SweepResult {
            parameters: vec![4.0, 5.0],
            metrics: vec![6.0],
        });
        let pairs = ps.scatter_pair(0, 1);
        assert_eq!(pairs.len(), 2);
        assert!((pairs[0].0 - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sim_trace_sample_at() {
        let mut t = SimTrace::new("energy");
        t.push(0.0, 0.0);
        t.push(1.0, 10.0);
        let v = t.sample_at(0.5).unwrap();
        assert!((v - 5.0).abs() < 1e-10, "v={v}");
    }
    #[test]
    fn test_sim_trace_std_dev() {
        let mut t = SimTrace::new("x");
        for i in 0..5 {
            t.push(i as f64, i as f64 * 2.0);
        }
        let sd = t.std_dev();
        assert!(sd > 0.0, "sd={sd}");
    }
    #[test]
    fn test_sim_compare_mean_abs_diff() {
        let mut a = SimTrace::new("a");
        let mut b = SimTrace::new("b");
        a.push(0.0, 1.0);
        a.push(1.0, 2.0);
        b.push(0.0, 2.0);
        b.push(1.0, 3.0);
        let sc = SimCompare::new(a, b);
        assert!((sc.mean_abs_diff() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_alert_system_energy_spike() {
        let mut al = AlertSystem::new(0.1, 0.001);
        al.check_energy(100.0, 0.0);
        al.check_energy(300.0, 0.01);
        assert!(!al.alerts.is_empty());
        assert!(al.unacknowledged_count(AlertSeverity::Warning) > 0);
    }
    #[test]
    fn test_alert_system_penetration() {
        let mut al = AlertSystem::new(0.5, 0.001);
        al.check_penetration(0.01, 0.0);
        assert!(al.alerts.iter().any(|a| a.severity == AlertSeverity::Error));
    }
    #[test]
    fn test_alert_system_acknowledge_all() {
        let mut al = AlertSystem::new(0.1, 0.001);
        al.check_energy(100.0, 0.0);
        al.check_energy(500.0, 0.01);
        al.acknowledge_all();
        assert_eq!(al.unacknowledged_count(AlertSeverity::Info), 0);
    }
    #[test]
    fn test_alert_system_clear_acknowledged() {
        let mut al = AlertSystem::new(0.1, 0.001);
        al.check_energy(100.0, 0.0);
        al.check_energy(500.0, 0.01);
        al.acknowledge_all();
        al.clear_acknowledged();
        assert!(al.alerts.is_empty());
    }
    #[test]
    fn test_time_range_normalize() {
        let tr = TimeRange::new(0.0, 10.0);
        assert!((tr.normalize(5.0) - 0.5).abs() < 1e-10);
        assert!(tr.normalize(-1.0).abs() < 1e-10);
        assert!((tr.normalize(11.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_energy_sample_total() {
        let s = EnergySample {
            time: 0.0,
            ke: 3.0,
            pe: 7.0,
        };
        assert!((s.total() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_performance_panel_peak_memory() {
        let mut pp = PerformancePanel::new(10);
        pp.record_frame(FrameTiming {
            broadphase_ms: 0.0,
            narrowphase_ms: 0.0,
            solver_ms: 0.0,
            integration_ms: 0.0,
            memory_bytes: 1000,
        });
        pp.record_frame(FrameTiming {
            broadphase_ms: 0.0,
            narrowphase_ms: 0.0,
            solver_ms: 0.0,
            integration_ms: 0.0,
            memory_bytes: 5000,
        });
        assert_eq!(pp.peak_memory, 5000);
    }
}
#[cfg(test)]
mod expanded_dashboard_tests {
    use super::*;
    #[test]
    fn test_snapshot_buffer_ring() {
        let mut buf = SnapshotBuffer::new(3);
        for i in 0..5 {
            buf.push(SimSnapshot {
                time: i as f64,
                ke: 1.0,
                pe: 1.0,
                temperature: 300.0,
                pressure: 1e5,
                n_active: 1,
                max_force: 0.0,
                rms_velocity: 0.0,
            });
        }
        assert_eq!(buf.len(), 3);
        assert!((buf.latest().unwrap().time - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_energy_drift() {
        let mut buf = SnapshotBuffer::new(10);
        buf.push(SimSnapshot {
            time: 0.0,
            ke: 5.0,
            pe: 5.0,
            temperature: 0.0,
            pressure: 0.0,
            n_active: 0,
            max_force: 0.0,
            rms_velocity: 0.0,
        });
        buf.push(SimSnapshot {
            time: 1.0,
            ke: 5.5,
            pe: 4.6,
            temperature: 0.0,
            pressure: 0.0,
            n_active: 0,
            max_force: 0.0,
            rms_velocity: 0.0,
        });
        let drift = buf.energy_drift().unwrap();
        assert!((drift - 0.01).abs() < 0.01, "drift={drift}");
    }
    #[test]
    fn test_convergence_absolute() {
        let crit = ConvergenceCriterion::Absolute(1e-6);
        assert!(crit.satisfied(1e-7, 1.0));
        assert!(!crit.satisfied(1e-5, 1.0));
    }
    #[test]
    fn test_convergence_relative() {
        let crit = ConvergenceCriterion::Relative(0.01);
        assert!(crit.satisfied(0.001, 1.0));
        assert!(!crit.satisfied(0.1, 1.0));
    }
    #[test]
    fn test_convergence_monitor() {
        let mut mon = ConvergenceMonitor::new(ConvergenceCriterion::Absolute(1e-4), 100);
        let mut converged = false;
        for i in 0..20 {
            let r = (0.5f64).powi(i + 1);
            if mon.record(r) {
                converged = true;
                break;
            }
        }
        assert!(converged, "Monitor should have converged");
        assert!(mon.converge_iter.is_some());
    }
    #[test]
    fn test_convergence_sparkline() {
        let mut mon = ConvergenceMonitor::new(ConvergenceCriterion::Absolute(1e-12), 50);
        for i in 0..10 {
            mon.record(1.0 - i as f64 * 0.05);
        }
        let spark = mon.sparkline();
        assert!(!spark.is_empty(), "Sparkline should not be empty");
    }
    #[test]
    fn test_velocity_histogram_binning() {
        let mut hist = VelocityHistogram::new(10, 0.0, 10.0);
        hist.add(5.0);
        hist.add(5.5);
        hist.add(0.1);
        assert_eq!(hist.total, 3);
    }
    #[test]
    fn test_velocity_histogram_pdf() {
        let mut hist = VelocityHistogram::new(10, 0.0, 10.0);
        for v in (0..=100).map(|i| i as f64 * 0.1) {
            hist.add(v);
        }
        let pdf = hist.pdf();
        let sum: f64 = pdf
            .iter()
            .zip((0..10).map(|i| {
                let h = &hist;
                h.edges[i + 1] - h.edges[i]
            }))
            .map(|(p, dv)| p * dv)
            .sum();
        assert!((sum - 1.0).abs() < 0.01, "sum={sum}");
    }
    #[test]
    fn test_velocity_histogram_mean() {
        let mut hist = VelocityHistogram::new(20, 0.0, 20.0);
        hist.add_batch(&[5.0, 5.0, 5.0, 5.0]);
        let mean = hist.mean();
        assert!((mean - 5.0).abs() < 1.0, "mean={mean}");
    }
    #[test]
    fn test_param_sweep_optimal() {
        let mut sweep = ParamSweepViz::new("dt", "error");
        sweep.add(0.01, 0.001, None, true);
        sweep.add(0.1, 0.1, None, true);
        sweep.add(0.001, 0.0001, None, true);
        assert!((sweep.optimal_min().unwrap() - 0.001).abs() < 1e-10);
    }
    #[test]
    fn test_param_sweep_convergence() {
        let mut sweep = ParamSweepViz::new("stiffness", "force");
        sweep.add(1e3, 10.0, None, true);
        sweep.add(1e6, 100.0, None, true);
        sweep.add(1e9, 0.0, None, false);
        assert!((sweep.convergence_fraction() - 2.0 / 3.0).abs() < 0.01);
    }
    #[test]
    fn test_dashboard_standard_layout() {
        let layout = DashboardLayout::standard_4panel("Test Sim");
        assert_eq!(layout.panels.len(), 4);
        assert_eq!(layout.visible_count(), 4);
    }
    #[test]
    fn test_dashboard_out_of_bounds() {
        let mut layout = DashboardLayout::new(2, 2, "test");
        let panel = DashboardPanel::new("bad", PanelPosition::single(5, 5), PanelContent::Empty);
        assert!(!layout.add_panel(panel));
    }
    #[test]
    fn test_svg_builder_render() {
        let mut svg = SvgPlotBuilder::new(400, 300);
        svg.add_background("#1a1a2e");
        svg.add_line(&[0.0, 1.0, 2.0], &[0.0, 1.0, 0.5], "red", 1.5);
        let output = svg.render();
        assert!(output.contains("<svg"), "Should contain SVG tag");
        assert!(output.contains("polyline"), "Should contain polyline");
    }
    #[test]
    fn test_svg_coord_mapping() {
        let svg = SvgPlotBuilder::new(400, 300);
        let x_px = svg.map_x(0.5, 0.0, 1.0);
        let expected = svg.margins[0] as f64 + 0.5 * svg.plot_width() as f64;
        assert!(
            (x_px - expected).abs() < 0.1,
            "x_px={x_px}, expected={expected}"
        );
    }
    #[test]
    fn test_colormap_legend_ticks() {
        let legend = ColormapLegend::linear(0.0, 100.0, "Temperature (K)");
        let ticks = legend.tick_values();
        assert_eq!(ticks.len(), 5);
        assert!(ticks[0].abs() < 1e-10);
        assert!((ticks[4] - 100.0).abs() < 1e-10);
    }
    #[test]
    fn test_colormap_legend_normalize() {
        let legend = ColormapLegend::linear(0.0, 100.0, "Pressure");
        assert!((legend.normalize(0.0)).abs() < 1e-10);
        assert!((legend.normalize(100.0) - 1.0).abs() < 1e-10);
        assert!((legend.normalize(50.0) - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_colormap_legend_fit() {
        let mut legend = ColormapLegend::linear(0.0, 1.0, "Stress");
        legend.fit_range(&[5.0, 10.0, 3.0, 7.0]);
        assert!((legend.vmin - 3.0).abs() < 1e-10);
        assert!((legend.vmax - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_convergence_rate() {
        let mut mon = ConvergenceMonitor::new(ConvergenceCriterion::Absolute(1e-12), 100);
        for i in 0..5 {
            mon.record(10.0_f64.powi(-i));
        }
        let rate = mon.convergence_rate().unwrap();
        assert!(
            rate < 0.0,
            "Convergence rate should be negative for decreasing residuals"
        );
    }
    #[test]
    fn test_maxwell_boltzmann() {
        let mb = VelocityHistogram::maxwell_boltzmann(500.0, 1.381e-23 * 300.0, 3.34e-27);
        assert!(mb >= 0.0, "MB PDF should be non-negative");
    }
    #[test]
    fn test_real_time_metrics_fps() {
        let mut m = RealTimeMetrics::new(100);
        for _ in 0..10 {
            m.record_frame(16.667);
        }
        let fps = m.fps();
        assert!((fps - 60.0).abs() < 1.0, "fps={fps}");
    }
    #[test]
    fn test_real_time_metrics_step_time() {
        let mut m = RealTimeMetrics::new(50);
        m.record_frame(5.0);
        m.record_frame(15.0);
        m.record_frame(10.0);
        assert!((m.mean_step_ms() - 10.0).abs() < 1e-9);
        assert!((m.max_step_ms() - 15.0).abs() < 1e-9);
        assert!((m.min_step_ms() - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_time_series_buffer_capacity() {
        let mut buf = TimeSeriesBuffer::new(3);
        buf.push(1.0, 10.0);
        buf.push(2.0, 20.0);
        buf.push(3.0, 30.0);
        buf.push(4.0, 40.0);
        assert_eq!(buf.len(), 3);
        assert!((buf.times()[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_time_series_buffer_range() {
        let mut buf = TimeSeriesBuffer::new(20);
        for i in 0..10 {
            buf.push(i as f64, (i * i) as f64);
        }
        let slice = buf.range(2.0, 5.0);
        assert!(!slice.is_empty());
        for &(t, _) in &slice {
            assert!((2.0..=5.0).contains(&t));
        }
    }
    #[test]
    fn test_phase_space_diagram() {
        let mut ps = PhaseSpaceDiagram::new(1000);
        ps.add_point(1.0, 0.5);
        ps.add_point(-1.0, -0.3);
        ps.add_point(2.0, 1.0);
        let (xmin, xmax, ymin, ymax) = ps.bounds();
        assert!((xmin - (-1.0)).abs() < 1e-9);
        assert!((xmax - 2.0).abs() < 1e-9);
        assert!((ymin - (-0.3)).abs() < 1e-9);
        assert!((ymax - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_histogram_renderer_bins() {
        let data = vec![0.5, 1.5, 2.5, 1.2, 0.8];
        let hr = HistogramRenderer::from_data(&data, 3, 0.0, 3.0);
        assert_eq!(hr.bins.len(), 3);
        assert_eq!(hr.bins.iter().sum::<u32>(), 5);
    }
    #[test]
    fn test_histogram_renderer_normalize() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let hr = HistogramRenderer::from_data(&data, 4, 0.5, 4.5);
        let norm = hr.normalized();
        let total: f64 = norm.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "sum={total}");
    }
    #[test]
    fn test_heatmap_renderer_dimensions() {
        let hm = HeatmapRenderer::from_grid(vec![vec![1.0_f64, 2.0], vec![3.0, 4.0]], "pressure");
        assert_eq!(hm.rows(), 2);
        assert_eq!(hm.cols(), 2);
        assert!((hm.global_max() - 4.0).abs() < 1e-9);
        assert!((hm.global_min() - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_heatmap_cell_normalize() {
        let hm = HeatmapRenderer::from_grid(vec![vec![0.0_f64, 5.0], vec![10.0, 5.0]], "stress");
        assert!((hm.cell_normalized(0, 0) - 0.0).abs() < 1e-9);
        assert!((hm.cell_normalized(1, 0) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_dashboard_config_visibility() {
        let mut cfg = DashboardConfig::default_4panel();
        cfg.set_panel_visible("energy", false);
        assert!(!cfg.is_panel_visible("energy"));
        assert!(cfg.is_panel_visible("performance"));
    }
    #[test]
    fn test_data_stream_buffer_overflow() {
        let mut dsb = DataStreamBuffer::new(4, 2);
        for i in 0..10u64 {
            dsb.push_sample(i as f64, i);
        }
        assert!(dsb.len() <= 4);
    }
    #[test]
    fn test_data_stream_latest_n() {
        let mut dsb = DataStreamBuffer::new(20, 3);
        for i in 0..10 {
            dsb.push_sample(i as f64, i as u64);
        }
        let last3 = dsb.latest_n(3);
        assert_eq!(last3.len(), 3);
        assert!((last3[2].0 - 9.0).abs() < 1e-9);
    }
    #[test]
    fn test_chart_zoom() {
        let mut zp = ChartZoomPan::new(0.0, 10.0, -5.0, 5.0);
        zp.zoom(0.5);
        let (xmin, xmax, _, _) = zp.view();
        assert!(xmax - xmin < 10.0, "zoom should narrow view");
    }
    #[test]
    fn test_chart_pan() {
        let mut zp = ChartZoomPan::new(0.0, 20.0, 0.0, 5.0);
        zp.zoom(0.4);
        let (xmin_before, _, _, _) = zp.view();
        zp.pan(2.0, 0.0);
        let (xmin_after, _, _, _) = zp.view();
        assert!(
            (xmin_after - xmin_before - 2.0).abs() < 1e-9,
            "before={xmin_before} after={xmin_after}"
        );
    }
    #[test]
    fn test_multi_panel_layout_grid() {
        let layout = MultiPanelLayout::grid(2, 3, 800.0, 600.0);
        assert_eq!(layout.panels.len(), 6);
        let (w, h) = layout.panels[0].cell_size();
        assert!((w - 800.0 / 3.0).abs() < 1e-6);
        assert!((h - 600.0 / 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_performance_profiler_sum() {
        let mut pp = PerformanceProfiler::new(50);
        pp.record(StepBreakdown {
            broadphase_us: 100,
            narrowphase_us: 200,
            solver_us: 300,
            integration_us: 50,
            render_us: 50,
        });
        let avg = pp.average_breakdown();
        assert_eq!(avg.total_us(), 700);
    }
    #[test]
    fn test_performance_profiler_fps() {
        let mut pp = PerformanceProfiler::new(50);
        pp.record(StepBreakdown {
            broadphase_us: 0,
            narrowphase_us: 0,
            solver_us: 0,
            integration_us: 0,
            render_us: 10_000,
        });
        let fps = pp.estimated_fps();
        assert!((fps - 100.0).abs() < 1.0, "fps={fps}");
    }
    #[test]
    fn test_energy_time_series_window() {
        let mut ets = EnergyTimeSeries::new(10.0);
        for i in 0..20 {
            ets.record(i as f64, i as f64 * 2.0, i as f64 * 0.5);
        }
        let recent = ets.recent_ke(5);
        assert_eq!(recent.len(), 5);
    }
    #[test]
    fn test_momentum_time_series_conservation() {
        let mut mts = MomentumTimeSeries::new(100);
        mts.record(0.0, [1.0, 0.0, 0.0]);
        mts.record(1.0, [1.0, 0.0, 0.0]);
        mts.record(2.0, [1.001, 0.0, 0.0]);
        assert!(mts.drift() < 0.01);
    }
    #[test]
    fn test_phase_space_density() {
        let mut ps = PhaseSpaceDiagram::new(500);
        for _ in 0..100 {
            ps.add_point(0.0, 0.0);
        }
        ps.add_point(10.0, 10.0);
        assert!(!ps.is_empty());
        assert!(ps.len() <= 500);
    }
    #[test]
    fn test_heatmap_resample() {
        let hm = HeatmapRenderer::from_grid(
            vec![vec![1.0, 2.0, 3.0, 4.0], vec![5.0, 6.0, 7.0, 8.0]],
            "field",
        );
        let small = hm.resample_to(1, 2);
        assert_eq!(small.rows(), 1);
        assert_eq!(small.cols(), 2);
    }
    #[test]
    fn test_data_stream_statistics() {
        let mut dsb = DataStreamBuffer::new(100, 1);
        for i in 0..5 {
            dsb.push_sample(i as f64, i as u64);
        }
        let mean = dsb.mean_value();
        assert!((mean - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_dashboard_config_panel_count() {
        let cfg = DashboardConfig::default_4panel();
        assert_eq!(cfg.panels.len(), 4);
    }
    #[test]
    fn test_time_series_interpolate() {
        let mut buf = TimeSeriesBuffer::new(100);
        buf.push(0.0, 0.0);
        buf.push(1.0, 10.0);
        buf.push(2.0, 20.0);
        let v = buf.interpolate_at(0.5);
        assert!((v - 5.0).abs() < 1e-9, "v={v}");
    }
    #[test]
    fn test_multi_panel_by_name() {
        let layout = MultiPanelLayout::grid(2, 2, 400.0, 400.0);
        assert_eq!(layout.panels.len(), 4);
        for p in &layout.panels {
            let (w, h) = p.cell_size();
            assert!(w > 0.0 && h > 0.0);
        }
    }
}
