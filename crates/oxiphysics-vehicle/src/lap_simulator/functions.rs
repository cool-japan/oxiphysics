//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    FastestLapFinder, LapResult, LapSimConfig, OptimalLap, RaceSimResult, TireCompound,
    TireDegradationState, TrackLayout, TrackSegment,
};

/// Simulate a full lap given a track layout and configuration.
///
/// Returns a [`LapResult`] with estimated lap time and sector breakdown.
#[allow(clippy::too_many_arguments)]
pub fn simulate_lap(
    track: &TrackLayout,
    config: &mut LapSimConfig,
    lap_number: usize,
) -> LapResult {
    let mut total_time = 0.0_f64;
    let mut sector_times = vec![0.0_f64; 3];
    let mut max_speed = 0.0_f64;
    let mut total_distance = 0.0_f64;
    let fuel_kg = config.fuel.fuel_at_lap(lap_number);
    let fuel_penalty = config.fuel.laptime_penalty(fuel_kg);
    let opt = OptimalLap {
        perf: config.perf.clone(),
        driver: config.driver.clone(),
    };
    let sector_len = track.total_length / 3.0;
    let mut current_sector = 0usize;
    let mut v_current = 30.0_f64;
    for seg in &track.segments {
        let seg_len = seg.length();
        let seg_time = match seg {
            TrackSegment::Straight {
                length,
                drs_available,
            } => {
                let v_lim = if *drs_available {
                    config.perf.v_max
                } else {
                    config.perf.v_max * 0.92
                };
                let t = opt.straight_time(*length, v_current, v_lim);
                v_current = v_lim.min(config.perf.v_max);
                max_speed = max_speed.max(v_current);
                t
            }
            TrackSegment::Corner { radius, angle, .. } => {
                let v_c = opt.v_corner(*radius).max(config.min_speed);
                let d_brake = opt.braking_distance(v_current, v_c);
                let t_brake = if d_brake > 1e-3 && v_current > v_c {
                    (v_current - v_c) / opt.perf.brake_g_at(v_current).max(1e-3)
                } else {
                    0.0
                };
                let arc = radius * angle.abs();
                let t_corner = arc / v_c.max(1.0);
                v_current = v_c;
                t_brake + t_corner
            }
            TrackSegment::Chicane { radii, angles } => {
                let mut t_chicane = 0.0;
                for (r, a) in radii.iter().zip(angles.iter()) {
                    let v_c = opt.v_corner(*r).max(config.min_speed);
                    let arc = r * a.abs();
                    let t_brake = if v_current > v_c {
                        (v_current - v_c) / opt.perf.brake_g_at(v_current).max(1e-3)
                    } else {
                        0.0
                    };
                    t_chicane += t_brake + arc / v_c.max(1.0);
                    v_current = v_c;
                }
                t_chicane
            }
            TrackSegment::Elevation { length, gradient } => {
                let g = 9.81;
                let a_adjust = g * gradient;
                let a_eff = (opt.perf.traction_g_at(v_current) - a_adjust).max(1.0);
                let v_new = (v_current * v_current + 2.0 * a_eff * length).sqrt();
                let v_new = v_new.min(config.perf.v_max);
                let t = 2.0 * length / (v_current + v_new).max(1.0);
                v_current = v_new;
                t
            }
        };
        total_time += seg_time;
        total_distance += seg_len;
        max_speed = max_speed.max(v_current);
        let sector_for_this = if total_distance < sector_len {
            0
        } else if total_distance < 2.0 * sector_len {
            1
        } else {
            2
        };
        if sector_for_this < 3 {
            sector_times[sector_for_this] += seg_time;
            current_sector = sector_for_this;
        }
    }
    let _ = current_sector;
    total_time += fuel_penalty;
    let avg_speed = if total_time > 1e-3 {
        total_distance / total_time
    } else {
        0.0
    };
    let ers_deployed = config.ers.deploy(config.ers.max_deploy_power, 1.0);
    config.ers.reset_lap();
    LapResult {
        lap_time: total_time,
        sector_times,
        max_speed,
        avg_speed,
        tire_degradation: 0.02 * lap_number as f64,
        fuel_consumption: config.fuel.fuel_per_lap,
        ers_deployed,
    }
}
/// Simulate a full race of `total_laps` using the given track and config.
///
/// Returns a [`RaceSimResult`] with per-lap breakdowns.
pub fn simulate_race(
    track: &TrackLayout,
    config: &LapSimConfig,
    total_laps: usize,
    compound: TireCompound,
    track_temp: f64,
) -> RaceSimResult {
    let mut lap_times = Vec::with_capacity(total_laps);
    let mut tire_degradation_vec = Vec::with_capacity(total_laps);
    let mut fuel_masses = Vec::with_capacity(total_laps);
    let mut ers_deployed_vec = Vec::with_capacity(total_laps);
    let mut tire = TireDegradationState::new(compound, track_temp);
    let mut finder = FastestLapFinder::new();
    for lap_num in 1..=total_laps {
        let mut lap_config = config.clone();
        let fuel = lap_config.fuel.fuel_at_lap(lap_num);
        fuel_masses.push(fuel);
        let result = simulate_lap(track, &mut lap_config, lap_num);
        let tire_delta = tire.laptime_delta_s();
        let adjusted_time = result.lap_time + tire_delta;
        lap_times.push(adjusted_time);
        finder.add_lap(adjusted_time);
        tire.advance_lap();
        tire_degradation_vec.push(tire.degradation);
        ers_deployed_vec.push(result.ers_deployed);
    }
    let total_race_time = lap_times.iter().sum();
    let (fastest_idx, fastest_time) = finder.fastest().unwrap_or((0, 0.0));
    RaceSimResult {
        lap_times,
        tire_degradation: tire_degradation_vec,
        fuel_masses,
        ers_deployed: ers_deployed_vec,
        total_race_time,
        fastest_lap_number: fastest_idx + 1,
        fastest_lap_time: fastest_time,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::lap_simulator::CompoundOptimiser;
    use crate::lap_simulator::CornerApex;
    use crate::lap_simulator::DriverModel;
    use crate::lap_simulator::EnergyRecovery;
    use crate::lap_simulator::ErsZoneMap;
    use crate::lap_simulator::FuelEffect;
    use crate::lap_simulator::GGDiagram;
    use crate::lap_simulator::LapCompare;
    use crate::lap_simulator::LapDelta;
    use crate::lap_simulator::PitStopStrategy;
    use crate::lap_simulator::RaceConditions;
    use crate::lap_simulator::RacingLine;
    use crate::lap_simulator::SectorTiming;
    use crate::lap_simulator::VehiclePerformanceMap;
    #[test]
    fn straight_length() {
        let s = TrackSegment::Straight {
            length: 500.0,
            drs_available: false,
        };
        assert!((s.length() - 500.0).abs() < 1e-10);
    }
    #[test]
    fn corner_length_arc() {
        let c = TrackSegment::Corner {
            radius: 100.0,
            angle: std::f64::consts::FRAC_PI_2,
            direction: 1.0,
        };
        let expected = 100.0 * std::f64::consts::FRAC_PI_2;
        assert!((c.length() - expected).abs() < 1e-6);
    }
    #[test]
    fn chicane_min_radius() {
        let c = TrackSegment::Chicane {
            radii: vec![30.0, 50.0],
            angles: vec![0.5, 0.5],
        };
        assert_eq!(c.min_radius(), Some(30.0));
    }
    #[test]
    fn elevation_length() {
        let e = TrackSegment::Elevation {
            length: 200.0,
            gradient: 0.05,
        };
        assert!((e.length() - 200.0).abs() < 1e-10);
    }
    #[test]
    fn drs_zone_flag() {
        let s = TrackSegment::Straight {
            length: 600.0,
            drs_available: true,
        };
        assert!(s.is_drs_zone());
        let s2 = TrackSegment::Straight {
            length: 300.0,
            drs_available: false,
        };
        assert!(!s2.is_drs_zone());
    }
    #[test]
    fn straight_has_no_radius() {
        let s = TrackSegment::Straight {
            length: 600.0,
            drs_available: false,
        };
        assert!(s.min_radius().is_none());
    }
    #[test]
    fn track_layout_total_length_positive() {
        let t = TrackLayout::f1_example();
        assert!(t.total_length > 0.0, "total length: {}", t.total_length);
    }
    #[test]
    fn track_layout_drs_zones_detected() {
        let t = TrackLayout::f1_example();
        assert!(t.drs_zone_count() > 0, "should have DRS zones");
    }
    #[test]
    fn track_layout_corner_count() {
        let t = TrackLayout::f1_example();
        assert!(t.corner_count() > 0);
    }
    #[test]
    fn track_layout_sector_boundaries() {
        let t = TrackLayout::f1_example();
        assert_eq!(t.sector_boundaries.len(), 2);
    }
    #[test]
    fn perf_map_lat_g_positive() {
        let p = VehiclePerformanceMap::default_f1();
        assert!(p.lat_g_at(50.0) > 0.0);
    }
    #[test]
    fn perf_map_brake_g_positive() {
        let p = VehiclePerformanceMap::default_f1();
        assert!(p.brake_g_at(80.0) > 0.0);
    }
    #[test]
    fn perf_map_drag_increases_with_speed() {
        let p = VehiclePerformanceMap::default_f1();
        assert!(p.drag_at(80.0) > p.drag_at(40.0));
    }
    #[test]
    fn perf_map_traction_at_low_speed() {
        let p = VehiclePerformanceMap::default_f1();
        let trac_low = p.traction_g_at(5.0);
        assert!(trac_low > 0.0);
    }
    #[test]
    fn driver_perfect_grip_fraction() {
        let d = DriverModel::perfect();
        assert!((d.effective_grip_fraction() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn driver_experienced_less_than_perfect() {
        let d = DriverModel::experienced();
        assert!(d.effective_grip_fraction() < 1.0);
    }
    #[test]
    fn optimal_lap_v_corner_positive() {
        let opt = OptimalLap::new_f1();
        let v = opt.v_corner(100.0);
        assert!(v > 0.0, "corner speed should be positive: {v}");
    }
    #[test]
    fn optimal_lap_v_corner_increases_with_radius() {
        let opt = OptimalLap::new_f1();
        let v1 = opt.v_corner(50.0);
        let v2 = opt.v_corner(150.0);
        assert!(v2 > v1, "larger radius → higher corner speed: {v1} vs {v2}");
    }
    #[test]
    fn optimal_lap_braking_distance_positive() {
        let opt = OptimalLap::new_f1();
        let d = opt.braking_distance(80.0, 30.0);
        assert!(d > 0.0, "braking distance: {d}");
    }
    #[test]
    fn optimal_lap_braking_distance_zero_if_slower() {
        let opt = OptimalLap::new_f1();
        let d = opt.braking_distance(20.0, 30.0);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn optimal_lap_corner_time_positive() {
        let opt = OptimalLap::new_f1();
        let t = opt.corner_time(80.0, std::f64::consts::FRAC_PI_2);
        assert!(t > 0.0, "corner time: {t}");
    }
    #[test]
    fn fuel_at_lap_decreases() {
        let f = FuelEffect::default_f1_one_stop();
        assert!(f.fuel_at_lap(1) > f.fuel_at_lap(10));
    }
    #[test]
    fn fuel_lap_time_penalty_proportional() {
        let f = FuelEffect::default_f1_one_stop();
        let p1 = f.laptime_penalty(10.0);
        let p2 = f.laptime_penalty(20.0);
        assert!((p2 - 2.0 * p1).abs() < 1e-10);
    }
    #[test]
    fn fuel_total_race_time_positive() {
        let f = FuelEffect::default_f1_one_stop();
        let t = f.total_race_time(90.0, 50);
        assert!(t > 0.0 && t < 10_000.0, "race time out of range: {t}");
    }
    #[test]
    fn fuel_at_lap_floor_zero() {
        let f = FuelEffect::default_f1_one_stop();
        let fuel = f.fuel_at_lap(1_000);
        assert_eq!(fuel, 0.0);
    }
    #[test]
    fn ers_deploy_reduces_soc() {
        let mut ers = EnergyRecovery::default_f1_mguk();
        let initial = ers.soc_j;
        ers.deploy(120_000.0, 1.0);
        assert!(ers.soc_j < initial);
    }
    #[test]
    fn ers_harvest_increases_soc() {
        let mut ers = EnergyRecovery::default_f1_mguk();
        ers.soc_j = 0.0;
        ers.harvest(50_000.0, 1.0);
        assert!(ers.soc_j > 0.0);
    }
    #[test]
    fn ers_soc_never_exceeds_capacity() {
        let mut ers = EnergyRecovery::default_f1_mguk();
        ers.harvest(1_000_000.0, 100.0);
        assert!(ers.soc_j <= ers.capacity_j);
    }
    #[test]
    fn ers_soc_fraction_range() {
        let ers = EnergyRecovery::default_f1_mguk();
        let f = ers.soc_fraction();
        assert!((0.0..=1.0).contains(&f), "soc fraction: {f}");
    }
    #[test]
    fn ers_reset_lap_clears_counters() {
        let mut ers = EnergyRecovery::default_f1_mguk();
        ers.deploy(100_000.0, 0.5);
        ers.reset_lap();
        assert_eq!(ers.deployed_this_lap, 0.0);
        assert_eq!(ers.harvested_this_lap, 0.0);
    }
    #[test]
    fn lap_result_time_str_format() {
        let r = LapResult {
            lap_time: 90.123,
            sector_times: vec![30.0, 30.0, 30.123],
            max_speed: 95.0,
            avg_speed: 50.0,
            tire_degradation: 0.05,
            fuel_consumption: 2.2,
            ers_deployed: 100_000.0,
        };
        let s = r.lap_time_str();
        assert!(s.starts_with("1:"), "expected 1:... format, got: {s}");
    }
    #[test]
    fn lap_result_sector_time_some() {
        let r = LapResult {
            lap_time: 90.0,
            sector_times: vec![30.0, 30.0, 30.0],
            max_speed: 90.0,
            avg_speed: 45.0,
            tire_degradation: 0.0,
            fuel_consumption: 2.2,
            ers_deployed: 0.0,
        };
        assert_eq!(r.sector_time(1), Some(30.0));
        assert_eq!(r.sector_time(5), None);
    }
    #[test]
    fn lap_compare_delta_sign() {
        let r = LapResult {
            lap_time: 90.0,
            sector_times: vec![30.0, 30.0, 30.0],
            max_speed: 90.0,
            avg_speed: 45.0,
            tire_degradation: 0.0,
            fuel_consumption: 2.2,
            ers_deployed: 0.0,
        };
        let c = LapResult {
            lap_time: 89.5,
            sector_times: vec![29.8, 30.0, 29.7],
            max_speed: 91.0,
            avg_speed: 46.0,
            tire_degradation: 0.01,
            fuel_consumption: 2.2,
            ers_deployed: 100_000.0,
        };
        let cmp = LapCompare::new(r, c);
        assert!(cmp.lap_time_delta() < 0.0, "comparison should be faster");
        assert!(cmp.is_comparison_faster());
    }
    #[test]
    fn lap_compare_sector_deltas_count() {
        let r = LapResult {
            lap_time: 90.0,
            sector_times: vec![30.0, 30.0, 30.0],
            max_speed: 90.0,
            avg_speed: 45.0,
            tire_degradation: 0.0,
            fuel_consumption: 2.2,
            ers_deployed: 0.0,
        };
        let c = r.clone();
        let cmp = LapCompare::new(r, c);
        assert_eq!(cmp.sector_deltas().len(), 3);
    }
    #[test]
    fn simulate_lap_returns_positive_time() {
        let track = TrackLayout::f1_example();
        let mut config = LapSimConfig::default_f1();
        let result = simulate_lap(&track, &mut config, 1);
        assert!(result.lap_time > 0.0, "lap time: {}", result.lap_time);
    }
    #[test]
    fn simulate_lap_max_speed_positive() {
        let track = TrackLayout::f1_example();
        let mut config = LapSimConfig::default_f1();
        let result = simulate_lap(&track, &mut config, 1);
        assert!(result.max_speed > 0.0);
    }
    #[test]
    fn simulate_lap_tire_degradation_increases() {
        let track = TrackLayout::f1_example();
        let mut config = LapSimConfig::default_f1();
        let r1 = simulate_lap(&track, &mut config, 1);
        let r5 = simulate_lap(&track, &mut config, 5);
        assert!(r5.tire_degradation > r1.tire_degradation);
    }
    #[test]
    fn simulate_lap_fuel_consumption_positive() {
        let track = TrackLayout::f1_example();
        let mut config = LapSimConfig::default_f1();
        let result = simulate_lap(&track, &mut config, 1);
        assert!(result.fuel_consumption > 0.0);
    }
    #[test]
    fn corner_apex_detect_finds_corners() {
        let track = TrackLayout::f1_example();
        let opt = OptimalLap::new_f1();
        let apexes = CornerApex::detect_all(&track, &opt);
        assert!(!apexes.is_empty(), "should detect at least one apex");
    }
    #[test]
    fn corner_apex_speed_positive() {
        let track = TrackLayout::f1_example();
        let opt = OptimalLap::new_f1();
        let apexes = CornerApex::detect_all(&track, &opt);
        for apex in &apexes {
            assert!(apex.apex_speed > 0.0, "apex speed must be positive");
        }
    }
    #[test]
    fn corner_apex_in_braking_zone() {
        let track = TrackLayout::f1_example();
        let opt = OptimalLap::new_f1();
        let apexes = CornerApex::detect_all(&track, &opt);
        if let Some(apex) = apexes.first() {
            let pos = apex.braking_zone_start + 1.0;
            assert!(apex.in_braking_zone(pos));
        }
    }
    #[test]
    fn corner_apex_in_accel_zone() {
        let track = TrackLayout::f1_example();
        let opt = OptimalLap::new_f1();
        let apexes = CornerApex::detect_all(&track, &opt);
        if let Some(apex) = apexes.first() {
            let pos = apex.apex_distance + 1.0;
            if pos <= apex.accel_zone_end {
                assert!(apex.in_accel_zone(pos));
            }
        }
    }
    #[test]
    fn racing_line_total_curvature_nonnegative() {
        let track = TrackLayout::f1_example();
        let rl = RacingLine::new(&track, 10.0);
        assert!(rl.total_curvature_cost() >= 0.0);
    }
    #[test]
    fn racing_line_optimise_reduces_curvature() {
        let track = TrackLayout::f1_example();
        let mut rl = RacingLine::new(&track, 10.0);
        let before = rl.total_curvature_cost();
        rl.optimise(5, 0.1);
        let after = rl.total_curvature_cost();
        assert!(
            after <= before + 1e-6,
            "optimisation should not increase curvature"
        );
    }
    #[test]
    fn racing_line_effective_corner_speed_increases() {
        let track = TrackLayout::f1_example();
        let opt = OptimalLap::new_f1();
        let mut rl = RacingLine::new(&track, 10.0);
        let corner_idx = track
            .segments
            .iter()
            .position(|s| matches!(s, TrackSegment::Corner { .. }))
            .unwrap_or(0);
        let speed_before = rl.effective_corner_speed(corner_idx, &opt);
        rl.optimise(3, 0.5);
        let speed_after = rl.effective_corner_speed(corner_idx, &opt);
        assert!(speed_after >= speed_before - 1e-6);
    }
    #[test]
    fn gg_diagram_len_matches_samples() {
        let perf = VehiclePerformanceMap::default_f1();
        let gg = GGDiagram::build(&perf, 50.0, 10);
        assert_eq!(gg.len(), 21);
    }
    #[test]
    fn gg_diagram_feasible_zero_accel() {
        let perf = VehiclePerformanceMap::default_f1();
        let gg = GGDiagram::build(&perf, 50.0, 20);
        assert!(gg.is_feasible(0.0, 0.0));
    }
    #[test]
    fn gg_diagram_infeasible_huge_accel() {
        let perf = VehiclePerformanceMap::default_f1();
        let gg = GGDiagram::build(&perf, 50.0, 20);
        assert!(!gg.is_feasible(0.0, 1000.0));
    }
    #[test]
    fn soft_grip_higher_than_hard() {
        assert!(TireCompound::Soft.peak_grip() > TireCompound::Hard.peak_grip());
    }
    #[test]
    fn soft_degrades_faster_than_hard() {
        assert!(TireCompound::Soft.degradation_rate() > TireCompound::Hard.degradation_rate());
    }
    #[test]
    fn tire_degradation_advances_lap() {
        let mut state = TireDegradationState::new(TireCompound::Soft, 100.0);
        let initial_deg = state.degradation;
        state.advance_lap();
        assert!(state.degradation > initial_deg);
    }
    #[test]
    fn tire_at_cliff_after_many_laps() {
        let mut state = TireDegradationState::new(TireCompound::Soft, 110.0);
        for _ in 0..50 {
            state.advance_lap();
        }
        assert!(state.at_cliff());
    }
    #[test]
    fn tire_grip_fraction_one_when_new() {
        let state = TireDegradationState::new(TireCompound::Medium, 95.0);
        let grip = state.grip_fraction();
        assert!(grip > 0.9 && grip <= 1.1, "grip={grip}");
    }
    #[test]
    fn sector_timing_personal_best_updates() {
        let mut st = SectorTiming::new(0, 30.0);
        st.update(28.5);
        assert!((st.personal_best - 28.5).abs() < 1e-10);
    }
    #[test]
    fn sector_timing_delta_to_pb_nonnegative_after_slower() {
        let mut st = SectorTiming::new(0, 30.0);
        st.update(31.0);
        assert!(st.delta_to_pb() >= 0.0);
    }
    #[test]
    fn sector_timing_is_personal_best_initial() {
        let st = SectorTiming::new(0, 30.0);
        assert!(st.is_personal_best());
    }
    #[test]
    fn lap_delta_final_delta_correct() {
        let pos = vec![0.0, 100.0, 200.0];
        let ref_t = vec![10.0, 20.0, 30.0];
        let cmp_t = vec![9.5, 19.5, 29.5];
        let ld = LapDelta::from_segments(pos, &ref_t, &cmp_t);
        assert!((ld.final_delta - (-0.5)).abs() < 1e-10);
    }
    #[test]
    fn lap_delta_max_gain_negative() {
        let pos = vec![0.0, 100.0];
        let ref_t = vec![10.0, 20.0];
        let cmp_t = vec![9.0, 18.0];
        let ld = LapDelta::from_segments(pos, &ref_t, &cmp_t);
        assert!(ld.max_gain() < 0.0);
    }
    #[test]
    fn one_stop_strategy_has_two_stints() {
        let strategy = PitStopStrategy::one_stop_soft_hard(50, 22.0, 25);
        assert_eq!(strategy.stints.len(), 2);
        assert_eq!(strategy.pit_count(), 1);
    }
    #[test]
    fn two_stop_strategy_has_three_stints() {
        let strategy = PitStopStrategy::two_stop_soft_medium_soft(50, 22.0, 15, 35);
        assert_eq!(strategy.stints.len(), 3);
    }
    #[test]
    fn optimise_one_stop_returns_valid_pit_lap() {
        let (pit_lap, race_time) = PitStopStrategy::optimise_one_stop(50, 90.0, 22.0, 10, 40);
        assert!((10..=40).contains(&pit_lap), "pit_lap={pit_lap}");
        assert!(race_time > 0.0);
    }
    #[test]
    fn compound_at_lap_returns_correct() {
        let strategy = PitStopStrategy::one_stop_soft_hard(50, 22.0, 20);
        assert_eq!(strategy.compound_at_lap(10), Some(TireCompound::Soft));
        assert_eq!(strategy.compound_at_lap(30), Some(TireCompound::Hard));
    }
    #[test]
    fn fastest_lap_finder_finds_minimum() {
        let mut f = FastestLapFinder::new();
        f.add_lap(91.0);
        f.add_lap(89.5);
        f.add_lap(90.2);
        let (idx, t) = f.fastest().unwrap();
        assert_eq!(idx, 1);
        assert!((t - 89.5).abs() < 1e-10);
    }
    #[test]
    fn fastest_lap_finder_std_dev_positive() {
        let mut f = FastestLapFinder::new();
        for t in [90.0, 91.0, 92.0, 89.0] {
            f.add_lap(t);
        }
        assert!(f.lap_time_std_dev() > 0.0);
    }
    #[test]
    fn compound_optimiser_returns_wet_when_wet() {
        let cond = RaceConditions {
            track_temp_c: 10.0,
            remaining_laps: 30,
            wet: true,
            laps_on_current: 0,
        };
        let rec = CompoundOptimiser::select(&cond);
        assert_eq!(rec.compound, TireCompound::Wet);
    }
    #[test]
    fn compound_optimiser_dry_returns_valid_compound() {
        let cond = RaceConditions {
            track_temp_c: 40.0,
            remaining_laps: 20,
            wet: false,
            laps_on_current: 5,
        };
        let rec = CompoundOptimiser::select(&cond);
        assert!(matches!(
            rec.compound,
            TireCompound::Soft | TireCompound::Medium | TireCompound::Hard
        ));
    }
    #[test]
    fn ers_zone_map_has_deploy_zones() {
        let track = TrackLayout::f1_example();
        let map = ErsZoneMap::from_track(&track);
        assert!(map.deploy_zone_count() > 0);
    }
    #[test]
    fn ers_zone_map_has_harvest_zones() {
        let track = TrackLayout::f1_example();
        let map = ErsZoneMap::from_track(&track);
        assert!(map.harvest_zone_count() > 0);
    }
    #[test]
    fn simulate_race_total_time_positive() {
        let track = TrackLayout::f1_example();
        let config = LapSimConfig::default_f1();
        let result = simulate_race(&track, &config, 5, TireCompound::Medium, 40.0);
        assert!(result.total_race_time > 0.0);
    }
    #[test]
    fn simulate_race_fastest_lap_within_bounds() {
        let track = TrackLayout::f1_example();
        let config = LapSimConfig::default_f1();
        let result = simulate_race(&track, &config, 5, TireCompound::Soft, 35.0);
        assert!(result.fastest_lap_number >= 1 && result.fastest_lap_number <= 5);
    }
    #[test]
    fn simulate_race_lap_count_matches() {
        let track = TrackLayout::f1_example();
        let config = LapSimConfig::default_f1();
        let result = simulate_race(&track, &config, 10, TireCompound::Hard, 45.0);
        assert_eq!(result.lap_times.len(), 10);
    }
}
