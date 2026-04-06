//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[inline]
pub(super) fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn vec3_dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = vec3_sub(a, b);
    vec3_dot(d, d)
}
#[cfg(test)]
mod tests {

    use crate::telemetry::BrakeTraceAnalyzer;
    use crate::telemetry::CornerAnalyzer;
    use crate::telemetry::DriverInputAnalyzer;
    use crate::telemetry::EnergyConsumption;
    use crate::telemetry::ErsEnergyLog;
    use crate::telemetry::FuelModel;
    use crate::telemetry::GForceRecorder;
    use crate::telemetry::LapComparison;
    use crate::telemetry::LapTimer;
    use crate::telemetry::PitEvent;
    use crate::telemetry::PitlaneDetector;
    use crate::telemetry::SectorComparison;
    use crate::telemetry::SessionSummary;
    use crate::telemetry::SpeedTrace;
    use crate::telemetry::TelemetryChannel;
    use crate::telemetry::TelemetryFrame;
    use crate::telemetry::TelemetryRecorder;
    use crate::telemetry::TelemetryReplay;
    use crate::telemetry::TelemetryStatistics;
    use crate::telemetry::TireSlipRecorder;
    use crate::telemetry::TireTempTrend;
    fn make_frame(time: f64, speed: f64) -> TelemetryFrame {
        TelemetryFrame {
            time,
            position: [speed * time, 0.0, 0.0],
            velocity: [speed, 0.0, 0.0],
            acceleration: [0.0; 3],
            steering_angle: 0.0,
            throttle: 0.5,
            brake: 0.0,
            gear: 3,
            rpm: 5000.0,
            wheel_speeds: [10.0; 4],
            lateral_g: 0.3,
            longitudinal_g: 0.5,
        }
    }
    #[test]
    fn lap_timer_sector_sum_equals_lap_time() {
        let mut timer = LapTimer::new(3);
        timer.start_lap(0.0);
        timer.complete_sector(30.0);
        timer.complete_sector(55.0);
        timer.complete_sector(80.0);
        let lap = timer.complete_lap(80.0);
        assert!((lap - 80.0).abs() < 1e-9, "lap time should be 80 s");
        let sectors = &timer.sector_times[0];
        assert_eq!(sectors.len(), 3);
        let sector_sum: f64 = sectors.iter().sum();
        assert!(
            (sector_sum - lap).abs() < 1e-9,
            "sector sum {sector_sum} != lap {lap}"
        );
        assert_eq!(timer.last_lap(), Some(80.0));
        assert_eq!(timer.best_lap(), Some(80.0));
    }
    #[test]
    fn telemetry_csv_line_count() {
        let mut rec = TelemetryRecorder::new(100);
        for i in 0..10 {
            rec.record(make_frame(i as f64 * 0.1, 20.0));
        }
        let csv = rec.export_csv();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 11, "expected 11 lines, got {}", lines.len());
    }
    #[test]
    fn speed_and_distance_consistency() {
        let mut rec = TelemetryRecorder::new(100);
        let speed = 20.0_f64;
        for i in 0..10 {
            rec.record(make_frame(i as f64 * 0.1, speed));
        }
        let avg = rec.average_speed();
        assert!(
            (avg - speed).abs() < 1e-6,
            "expected avg speed {speed}, got {avg}"
        );
        let dist = rec.distance_traveled();
        let expected_dist = 9.0 * 0.1 * speed;
        assert!(
            (dist - expected_dist).abs() < 1e-6,
            "expected distance {expected_dist}, got {dist}"
        );
    }
    #[test]
    fn fuel_model_formula_car_initial_full() {
        let m = FuelModel::formula_car();
        assert!((m.fuel_fraction() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn fuel_model_step_reduces_fuel() {
        let mut m = FuelModel::formula_car();
        let initial = m.fuel_kg;
        let consumed = m.step(0.5, 8000.0, 1.0);
        assert!(consumed > 0.0, "should consume fuel");
        assert!(m.fuel_kg < initial, "fuel level should decrease");
    }
    #[test]
    fn fuel_model_zero_throttle_base_rate_only() {
        let m = FuelModel::formula_car();
        let rate = m.consumption_rate(0.0, 0.0);
        assert!(
            (rate - m.base_rate).abs() < 1e-12,
            "zero throttle, zero RPM → base rate only"
        );
    }
    #[test]
    fn fuel_model_remaining_time_decreases_with_higher_throttle() {
        let m = FuelModel::formula_car();
        let t_low = m.remaining_time(0.1, 5000.0);
        let t_high = m.remaining_time(0.9, 5000.0);
        assert!(t_high < t_low, "higher throttle → less remaining time");
    }
    #[test]
    fn fuel_model_empty_tank_returns_zero_rate() {
        let mut m = FuelModel::formula_car();
        m.fuel_kg = 0.0;
        assert_eq!(m.consumption_rate(1.0, 10000.0), 0.0);
    }
    #[test]
    fn fuel_model_fraction_clamped_to_one() {
        let mut m = FuelModel::formula_car();
        m.fuel_kg = m.capacity_kg * 1.5;
        assert!((m.fuel_fraction() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn tire_temp_trend_empty_current_is_none() {
        let t = TireTempTrend::new(10);
        assert!(t.current(0).is_none());
    }
    #[test]
    fn tire_temp_trend_current_returns_latest() {
        let mut t = TireTempTrend::new(10);
        t.record(0.0, [80.0, 81.0, 82.0, 83.0]);
        t.record(1.0, [90.0, 91.0, 92.0, 93.0]);
        assert!((t.current(0).unwrap() - 90.0).abs() < 1e-9);
    }
    #[test]
    fn tire_temp_trend_rate_of_change_positive_when_warming() {
        let mut t = TireTempTrend::new(10);
        t.record(0.0, [60.0; 4]);
        t.record(10.0, [80.0; 4]);
        let rate = t.rate_of_change(0);
        assert!((rate - 2.0).abs() < 1e-9, "expected 2 °C/s, got {rate}");
        assert!(t.is_warming(0));
    }
    #[test]
    fn tire_temp_trend_peak_and_trough() {
        let mut t = TireTempTrend::new(10);
        t.record(0.0, [60.0; 4]);
        t.record(1.0, [100.0; 4]);
        t.record(2.0, [50.0; 4]);
        assert!((t.peak(0) - 100.0).abs() < 1e-9);
        assert!((t.trough(0) - 50.0).abs() < 1e-9);
    }
    #[test]
    fn tire_temp_trend_mean_correct() {
        let mut t = TireTempTrend::new(10);
        t.record(0.0, [100.0; 4]);
        t.record(1.0, [200.0; 4]);
        let mean = t.mean(0);
        assert!((mean - 150.0).abs() < 1e-9, "expected mean=150, got {mean}");
    }
    #[test]
    fn tire_temp_trend_rolling_window_evicts_oldest() {
        let mut t = TireTempTrend::new(3);
        t.record(0.0, [10.0; 4]);
        t.record(1.0, [20.0; 4]);
        t.record(2.0, [30.0; 4]);
        t.record(3.0, [40.0; 4]);
        assert!(
            (t.trough(0) - 20.0).abs() < 1e-9,
            "oldest sample should be 20, trough={}",
            t.trough(0)
        );
    }
    #[test]
    fn channel_linear_interp_exact_sample() {
        let mut ch = TelemetryChannel::new("speed");
        ch.record(0.0, 0.0);
        ch.record(1.0, 10.0);
        ch.record(2.0, 20.0);
        assert!((ch.interpolate_linear(1.0) - 10.0).abs() < 1e-12);
    }
    #[test]
    fn channel_linear_interp_between_samples() {
        let mut ch = TelemetryChannel::new("speed");
        ch.record(0.0, 0.0);
        ch.record(2.0, 20.0);
        assert!((ch.interpolate_linear(1.0) - 10.0).abs() < 1e-12);
    }
    #[test]
    fn channel_linear_interp_clamped_before() {
        let mut ch = TelemetryChannel::new("x");
        ch.record(1.0, 5.0);
        ch.record(2.0, 10.0);
        assert!((ch.interpolate_linear(0.0) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn channel_linear_interp_clamped_after() {
        let mut ch = TelemetryChannel::new("x");
        ch.record(1.0, 5.0);
        ch.record(2.0, 10.0);
        assert!((ch.interpolate_linear(99.0) - 10.0).abs() < 1e-12);
    }
    #[test]
    fn channel_cubic_interp_monotone_sequence() {
        let mut ch = TelemetryChannel::new("v");
        for i in 0..10 {
            ch.record(i as f64, (i * i) as f64);
        }
        let v = ch.interpolate_cubic(4.5);
        assert!(
            (16.0..=25.0).contains(&v),
            "cubic interp at 4.5 out of range: {v}"
        );
    }
    #[test]
    fn channel_is_empty_initially() {
        let ch = TelemetryChannel::new("empty");
        assert!(ch.is_empty());
        assert_eq!(ch.len(), 0);
    }
    #[test]
    fn gforce_recorder_peak_g() {
        let mut rec = GForceRecorder::new(100);
        rec.record(0.0, [1.0, 0.0, 0.0]);
        rec.record(1.0, [3.0, 4.0, 0.0]);
        rec.record(2.0, [0.0, 0.0, 2.0]);
        assert!((rec.peak_resultant_g() - 5.0).abs() < 1e-9);
    }
    #[test]
    fn gforce_recorder_mean_g() {
        let mut rec = GForceRecorder::new(100);
        rec.record(0.0, [3.0, 4.0, 0.0]);
        rec.record(1.0, [0.0, 0.0, 5.0]);
        assert!((rec.mean_resultant_g() - 5.0).abs() < 1e-9);
    }
    #[test]
    fn gforce_recorder_time_above_threshold() {
        let mut rec = GForceRecorder::new(100);
        rec.record(0.0, [0.0; 3]);
        rec.record(1.0, [5.0, 0.0, 0.0]);
        rec.record(2.0, [5.0, 0.0, 0.0]);
        rec.record(3.0, [0.0; 3]);
        let t = rec.time_above_g(3.0);
        assert!(
            (t - 2.0).abs() < 1e-9,
            "time above 3G should be 2s, got {t}"
        );
    }
    #[test]
    fn energy_consumption_net_energy() {
        let mut e = EnergyConsumption::new();
        e.record_step(100.0, 20.0, 10.0);
        e.record_step(50.0, 10.0, 5.0);
        assert!((e.net_energy() - 120.0).abs() < 1e-9);
    }
    #[test]
    fn energy_consumption_recovery_efficiency() {
        let mut e = EnergyConsumption::new();
        e.record_step(100.0, 25.0, 10.0);
        assert!((e.recovery_efficiency() - 0.25).abs() < 1e-9);
    }
    #[test]
    fn energy_consumption_zero_consumed() {
        let e = EnergyConsumption::new();
        assert_eq!(e.recovery_efficiency(), 0.0);
    }
    #[test]
    fn energy_consumption_total_consumed_accumulates() {
        let mut e = EnergyConsumption::new();
        e.record_step(50.0, 0.0, 0.0);
        e.record_step(30.0, 0.0, 0.0);
        assert!((e.total_consumed - 80.0).abs() < 1e-9);
    }
    #[test]
    fn replay_returns_first_frame_immediately() {
        let frames: Vec<TelemetryFrame> = (0..5).map(|i| make_frame(i as f64, 10.0)).collect();
        let mut replay = TelemetryReplay::new(frames, 1.0);
        let f = replay.advance(0.0).unwrap();
        assert!((f.time).abs() < 1e-12);
    }
    #[test]
    fn replay_advances_through_frames() {
        let frames: Vec<TelemetryFrame> = (0..5).map(|i| make_frame(i as f64, 10.0)).collect();
        let mut replay = TelemetryReplay::new(frames, 1.0);
        replay.advance(2.5);
        let f = replay.advance(0.0).unwrap();
        assert!(
            (f.time - 2.0).abs() < 1e-12,
            "expected frame t=2, got {}",
            f.time
        );
    }
    #[test]
    fn replay_finished_after_last_frame() {
        let frames: Vec<TelemetryFrame> = (0..3).map(|i| make_frame(i as f64, 10.0)).collect();
        let mut replay = TelemetryReplay::new(frames, 1.0);
        replay.advance(100.0);
        assert!(replay.is_finished());
    }
    #[test]
    fn replay_reset_goes_back_to_start() {
        let frames: Vec<TelemetryFrame> = (0..5).map(|i| make_frame(i as f64, 10.0)).collect();
        let mut replay = TelemetryReplay::new(frames, 1.0);
        replay.advance(4.0);
        replay.reset();
        let f = replay.advance(0.0).unwrap();
        assert!((f.time).abs() < 1e-12, "after reset, should be at t=0");
    }
    #[test]
    fn replay_duration_correct() {
        let frames: Vec<TelemetryFrame> =
            (0..11).map(|i| make_frame(i as f64 * 0.5, 10.0)).collect();
        let replay = TelemetryReplay::new(frames, 1.0);
        assert!((replay.duration() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn driver_input_overlap_zero_no_overlap() {
        let mut a = DriverInputAnalyzer::new(100);
        for i in 0..5 {
            a.record(i as f64, 0.8, 0.0, 0.1);
        }
        assert!(a.overlap_duration() < 1e-12, "no brake → no overlap");
    }
    #[test]
    fn driver_input_overlap_positive_when_both_active() {
        let mut a = DriverInputAnalyzer::new(100);
        a.record(0.0, 0.8, 0.5, 0.0);
        a.record(1.0, 0.9, 0.6, 0.0);
        a.record(2.0, 0.7, 0.4, 0.0);
        let ov = a.overlap_duration();
        assert!(ov > 0.0, "both > 0 → overlap detected");
    }
    #[test]
    fn driver_input_peak_throttle_correct() {
        let mut a = DriverInputAnalyzer::new(20);
        a.record(0.0, 0.5, 0.0, 0.0);
        a.record(1.0, 1.0, 0.0, 0.0);
        a.record(2.0, 0.3, 0.0, 0.0);
        assert!((a.peak_throttle() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn driver_input_full_throttle_fraction() {
        let mut a = DriverInputAnalyzer::new(20);
        a.record(0.0, 0.5, 0.0, 0.0);
        a.record(1.0, 1.0, 0.0, 0.0);
        a.record(2.0, 1.0, 0.0, 0.0);
        a.record(3.0, 0.2, 0.0, 0.0);
        let frac = a.full_throttle_fraction();
        assert!(frac > 0.0 && frac <= 1.0, "fraction should be in (0,1]");
    }
    #[test]
    fn sector_comparison_total_delta_zero_identical() {
        let ref_s = vec![30.0, 25.0, 28.0];
        let cur_s = vec![30.0, 25.0, 28.0];
        let cmp = SectorComparison::new(ref_s, cur_s);
        assert!(
            cmp.total_delta().abs() < 1e-12,
            "identical laps → zero delta"
        );
    }
    #[test]
    fn sector_comparison_total_delta_positive_slower() {
        let cmp = SectorComparison::new(vec![30.0, 25.0, 28.0], vec![31.0, 26.0, 29.0]);
        let delta = cmp.total_delta();
        assert!((delta - 3.0).abs() < 1e-9, "should be 3 s slower");
    }
    #[test]
    fn sector_comparison_worst_sector_index() {
        let cmp = SectorComparison::new(vec![30.0, 25.0, 28.0], vec![30.5, 27.0, 28.1]);
        let worst = cmp.worst_sector().unwrap();
        assert_eq!(worst, 1, "sector 1 has the largest time loss");
    }
    #[test]
    fn sector_comparison_best_sector_negative_delta() {
        let cmp = SectorComparison::new(vec![30.0, 25.0, 28.0], vec![31.0, 23.0, 29.0]);
        let best = cmp.best_sector().unwrap();
        assert_eq!(best, 1, "sector 1 should be the best improvement");
    }
    #[test]
    fn sector_comparison_cumulative_deltas_count() {
        let cmp = SectorComparison::new(vec![30.0, 25.0, 28.0], vec![31.0, 25.0, 27.0]);
        let cum = cmp.cumulative_deltas();
        assert_eq!(cum.len(), 3, "three sectors → three cumulative entries");
        assert!((cum[0] - 1.0).abs() < 1e-9);
        assert!((cum[1] - 1.0).abs() < 1e-9);
        assert!((cum[2] - 0.0).abs() < 1e-9);
    }
    #[test]
    fn speed_trace_min_max() {
        let mut tr = SpeedTrace::new();
        tr.record(0.0, 10.0);
        tr.record(100.0, 50.0);
        tr.record(200.0, 20.0);
        assert!((tr.min_speed() - 10.0).abs() < 1e-9);
        assert!((tr.max_speed() - 50.0).abs() < 1e-9);
    }
    #[test]
    fn speed_trace_max_speed_distance() {
        let mut tr = SpeedTrace::new();
        tr.record(0.0, 10.0);
        tr.record(150.0, 60.0);
        tr.record(300.0, 30.0);
        let d = tr.max_speed_distance().unwrap();
        assert!((d - 150.0).abs() < 1e-9, "max speed should be at 150 m");
    }
    #[test]
    fn speed_trace_speed_at_known_point() {
        let mut tr = SpeedTrace::new();
        tr.record(0.0, 0.0);
        tr.record(100.0, 50.0);
        let v = tr.speed_at(50.0);
        assert!((v - 25.0).abs() < 1e-9);
    }
    #[test]
    fn speed_trace_delta_versus_same_trace_is_zero() {
        let mut tr = SpeedTrace::new();
        for i in 0..5 {
            tr.record(i as f64 * 100.0, 30.0);
        }
        let tr2 = SpeedTrace {
            distances: tr.distances.clone(),
            speeds: tr.speeds.clone(),
        };
        let deltas = tr.delta_versus(&tr2);
        for (_, d) in &deltas {
            assert!(d.abs() < 1e-9, "delta should be zero");
        }
    }
    #[test]
    fn brake_trace_event_count() {
        let mut bt = BrakeTraceAnalyzer::new(0.1);
        for i in 0..=10 {
            let t = i as f64 * 0.5;
            let p = if (2..=6).contains(&i) { 0.8 } else { 0.0 };
            bt.record(t, p);
        }
        let events = bt.detect_events();
        assert_eq!(
            events.len(),
            1,
            "one braking event expected, got {}",
            events.len()
        );
    }
    #[test]
    fn brake_trace_mean_peak_pressure() {
        let mut bt = BrakeTraceAnalyzer::new(0.05);
        bt.record(0.0, 0.0);
        bt.record(1.0, 0.9);
        bt.record(2.0, 0.9);
        bt.record(3.0, 0.0);
        let mean_peak = bt.mean_peak_pressure();
        assert!((mean_peak - 0.9).abs() < 1e-9);
    }
    #[test]
    fn brake_trace_total_braking_time_positive() {
        let mut bt = BrakeTraceAnalyzer::new(0.1);
        bt.record(0.0, 0.0);
        bt.record(1.0, 0.8);
        bt.record(2.0, 0.9);
        bt.record(3.0, 0.7);
        bt.record(4.0, 0.0);
        let t = bt.total_braking_time();
        assert!(t > 0.0, "should have positive braking time");
    }
    #[test]
    fn corner_analyzer_detects_corners() {
        let mut ca = CornerAnalyzer::new(0.5);
        for i in 0..=20 {
            let t = i as f64 * 0.5;
            let g = if (4..=12).contains(&i) { 1.2 } else { 0.1 };
            ca.record(t, g);
        }
        assert_eq!(ca.corner_count(), 1, "one corner expected");
    }
    #[test]
    fn corner_analyzer_mean_corner_g_positive() {
        let mut ca = CornerAnalyzer::new(0.5);
        ca.record(0.0, 0.0);
        ca.record(1.0, 1.0);
        ca.record(2.0, 1.5);
        ca.record(3.0, 0.8);
        ca.record(4.0, 0.0);
        let mean_g = ca.mean_corner_g();
        assert!(mean_g > 0.0, "corners detected → mean G > 0");
    }
    #[test]
    fn corner_analyzer_total_corner_time_positive() {
        let mut ca = CornerAnalyzer::new(0.5);
        ca.record(0.0, 0.0);
        ca.record(0.5, 1.0);
        ca.record(1.0, 1.0);
        ca.record(1.5, 1.0);
        ca.record(2.0, 0.0);
        let t = ca.total_corner_time();
        assert!(t > 0.0, "should have positive corner time");
    }
    #[test]
    fn telemetry_statistics_from_empty_recorder() {
        let rec = TelemetryRecorder::new(100);
        let stats = TelemetryStatistics::from_recorder(&rec);
        assert_eq!(stats.total_distance_m, 0.0);
        assert_eq!(stats.max_speed_ms, 0.0);
    }
    #[test]
    fn telemetry_statistics_max_speed() {
        let mut rec = TelemetryRecorder::new(100);
        for i in 0..5 {
            rec.record(make_frame(i as f64, (i as f64 + 1.0) * 10.0));
        }
        let stats = TelemetryStatistics::from_recorder(&rec);
        assert!((stats.max_speed_ms - 50.0).abs() < 1e-6);
    }
    #[test]
    fn telemetry_statistics_summary_line_non_empty() {
        let mut rec = TelemetryRecorder::new(100);
        for i in 0..5 {
            rec.record(make_frame(i as f64, 20.0));
        }
        let stats = TelemetryStatistics::from_recorder(&rec);
        let line = stats.summary_line();
        assert!(!line.is_empty(), "summary line should not be empty");
        assert!(line.contains("dist="), "should contain distance");
    }
    #[test]
    fn telemetry_statistics_duration_matches_frames() {
        let mut rec = TelemetryRecorder::new(100);
        rec.record(make_frame(0.0, 10.0));
        rec.record(make_frame(5.0, 10.0));
        let stats = TelemetryStatistics::from_recorder(&rec);
        assert!((stats.duration_s - 5.0).abs() < 1e-9);
    }
    #[test]
    fn tire_slip_recorder_peak_slip_angle() {
        let mut rec = TireSlipRecorder::new(100);
        rec.record(0.0, [0.05, 0.04, 0.06, 0.03], [0.0; 4]);
        rec.record(1.0, [0.12, 0.10, 0.08, 0.09], [0.0; 4]);
        assert!((rec.peak_slip_angle() - 0.12).abs() < 1e-9);
    }
    #[test]
    fn tire_slip_recorder_mean_slip_angle() {
        let mut rec = TireSlipRecorder::new(100);
        rec.record(0.0, [0.0, 0.0, 0.0, 0.0], [0.0; 4]);
        rec.record(1.0, [0.2, 0.0, 0.0, 0.0], [0.0; 4]);
        let mean = rec.mean_slip_angle(0);
        assert!((mean - 0.1).abs() < 1e-9, "mean should be 0.1, got {mean}");
    }
    #[test]
    fn tire_slip_recorder_peak_long_slip() {
        let mut rec = TireSlipRecorder::new(100);
        rec.record(0.0, [0.0; 4], [0.1, -0.2, 0.05, 0.15]);
        assert!((rec.peak_long_slip() - 0.2).abs() < 1e-9);
    }
    #[test]
    fn tire_slip_recorder_time_above_threshold() {
        let mut rec = TireSlipRecorder::new(100);
        rec.record(0.0, [0.0; 4], [0.0; 4]);
        rec.record(1.0, [0.15, 0.0, 0.0, 0.0], [0.0; 4]);
        rec.record(2.0, [0.15, 0.0, 0.0, 0.0], [0.0; 4]);
        rec.record(3.0, [0.05, 0.0, 0.0, 0.0], [0.0; 4]);
        let t = rec.time_above_slip_threshold(0, 0.1);
        assert!(t > 0.0, "should have time above threshold");
    }
    #[test]
    fn ers_log_initial_soc_correct() {
        let log = ErsEnergyLog::new(4_000_000.0, 0.8);
        assert!((log.current_soc() - 0.8).abs() < 1e-9);
    }
    #[test]
    fn ers_log_deploy_reduces_soc() {
        let mut log = ErsEnergyLog::new(4_000_000.0, 1.0);
        log.step(0.0, 120_000.0, 1.0);
        assert!(log.current_soc() < 1.0, "SoC should have decreased");
    }
    #[test]
    fn ers_log_harvest_increases_soc() {
        let mut log = ErsEnergyLog::new(4_000_000.0, 0.5);
        let soc_before = log.current_soc();
        log.step(0.0, -80_000.0, 1.0);
        assert!(
            log.current_soc() > soc_before,
            "harvesting should increase SoC"
        );
    }
    #[test]
    fn ers_log_total_deployed_positive() {
        let mut log = ErsEnergyLog::new(4_000_000.0, 1.0);
        log.step(0.0, 100_000.0, 0.5);
        log.step(0.5, 100_000.0, 0.5);
        log.step(1.0, 100_000.0, 0.5);
        let deployed = log.total_deployed_j();
        assert!((deployed - 100_000.0).abs() < 1.0, "deployed={deployed}");
    }
    #[test]
    fn ers_log_soc_clamped_to_zero() {
        let mut log = ErsEnergyLog::new(1_000_000.0, 0.1);
        log.step(0.0, 10_000_000.0, 1.0);
        assert!(log.current_soc() < 1e-9, "SoC should be ~0 after over-draw");
    }
    #[test]
    fn lap_comparison_time_deltas_zero_same_laps() {
        let frames: Vec<TelemetryFrame> = (0..5).map(|i| make_frame(i as f64, 20.0)).collect();
        let cmp = LapComparison::new(frames.clone(), frames);
        let deltas = cmp.time_deltas();
        for d in &deltas {
            assert!(d.abs() < 1e-12, "identical laps → zero time delta");
        }
    }
    #[test]
    fn lap_comparison_total_time_delta_positive_when_slower() {
        let ref_frames: Vec<TelemetryFrame> =
            (0..5).map(|i| make_frame(i as f64 * 0.9, 20.0)).collect();
        let tgt_frames: Vec<TelemetryFrame> =
            (0..5).map(|i| make_frame(i as f64 * 1.0, 20.0)).collect();
        let cmp = LapComparison::new(ref_frames, tgt_frames);
        let delta = cmp.total_time_delta();
        assert!(delta > 0.0, "target is slower → positive total delta");
    }
    #[test]
    fn lap_comparison_speed_deltas_zero_same_speed() {
        let frames: Vec<TelemetryFrame> = (0..4).map(|i| make_frame(i as f64, 15.0)).collect();
        let cmp = LapComparison::new(frames.clone(), frames);
        for d in cmp.speed_deltas() {
            assert!(d.abs() < 1e-9, "same speed → zero speed delta");
        }
    }
    #[test]
    fn lap_comparison_frames_faster_count() {
        let ref_frames: Vec<TelemetryFrame> = (0..5).map(|i| make_frame(i as f64, 20.0)).collect();
        let tgt_frames: Vec<TelemetryFrame> =
            (0..5).map(|i| make_frame(i as f64 * 0.8, 20.0)).collect();
        let cmp = LapComparison::new(ref_frames, tgt_frames);
        let faster = cmp.frames_faster_count();
        assert!(faster > 0, "target is faster on some frames");
    }
    #[test]
    fn pitlane_detector_single_stop() {
        let mut det = PitlaneDetector::new(15.0);
        det.update(0.0, 30.0);
        let e = det.update(1.0, 10.0);
        assert_eq!(e, Some(PitEvent::Entry));
        det.update(10.0, 10.0);
        let x = det.update(20.0, 20.0);
        assert_eq!(x, Some(PitEvent::Exit));
        assert_eq!(det.num_pit_stops(), 1);
        assert!(
            (det.total_pit_time() - 19.0).abs() < 1e-9,
            "pit time should be 19 s"
        );
    }
    #[test]
    fn pitlane_detector_no_stop_all_fast() {
        let mut det = PitlaneDetector::new(15.0);
        for i in 0..10 {
            det.update(i as f64, 30.0);
        }
        assert_eq!(det.num_pit_stops(), 0);
    }
    #[test]
    fn pitlane_detector_mean_stop_duration() {
        let mut det = PitlaneDetector::new(15.0);
        det.update(0.0, 30.0);
        det.update(5.0, 10.0);
        det.update(15.0, 20.0);
        det.update(30.0, 10.0);
        det.update(35.0, 20.0);
        assert_eq!(det.num_pit_stops(), 2);
        let mean = det.mean_stop_duration();
        assert!(
            (mean - 7.5).abs() < 1e-9,
            "mean stop duration should be 7.5 s, got {mean}"
        );
    }
    #[test]
    fn session_summary_best_lap() {
        let mut s = SessionSummary::new();
        s.record_lap(95.0, vec![30.0, 32.0, 33.0]);
        s.record_lap(92.0, vec![29.0, 31.0, 32.0]);
        s.record_lap(94.0, vec![30.0, 31.5, 32.5]);
        assert!((s.best_lap().unwrap() - 92.0).abs() < 1e-9);
    }
    #[test]
    fn session_summary_theoretical_best() {
        let mut s = SessionSummary::new();
        s.record_lap(95.0, vec![30.0, 32.0, 33.0]);
        s.record_lap(92.5, vec![29.0, 31.5, 32.0]);
        let tbl = s.theoretical_best_lap();
        assert_eq!(tbl.len(), 3, "three sectors");
        assert!((tbl[0] - 29.0).abs() < 1e-9);
        assert!((tbl[1] - 31.5).abs() < 1e-9);
    }
    #[test]
    fn session_summary_report_non_empty() {
        let mut s = SessionSummary::new();
        s.record_lap(90.0, vec![28.0, 30.0, 32.0]);
        s.set_pit_data(25.0, 1);
        let report = s.report();
        assert!(!report.is_empty());
        assert!(
            report.contains("Pit stops:"),
            "report should mention pit stops"
        );
    }
    #[test]
    fn session_summary_lap_std_zero_single_lap() {
        let mut s = SessionSummary::new();
        s.record_lap(90.0, vec![30.0, 30.0, 30.0]);
        assert!(s.lap_time_std() < 1e-9, "single lap → zero std");
    }
    #[test]
    fn session_summary_mean_lap_multiple() {
        let mut s = SessionSummary::new();
        s.record_lap(90.0, vec![]);
        s.record_lap(100.0, vec![]);
        assert!(
            (s.mean_lap() - 95.0).abs() < 1e-9,
            "mean of 90 and 100 = 95"
        );
    }
}
