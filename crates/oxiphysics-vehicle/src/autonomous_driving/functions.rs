//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Simple LCG random number generator (deterministic, no external deps).
pub(super) fn lcg_rand(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    (*state >> 33) as f64 / (1u64 << 31) as f64
}

/// Row-major 4x4 matrix product `a * b`.
pub(super) fn mat4_mul(a: &[f64; 16], b: &[f64; 16]) -> [f64; 16] {
    let mut c = [0.0f64; 16];
    for (c_row, a_row) in c.chunks_exact_mut(4).zip(a.chunks_exact(4)) {
        for (j, c_ij) in c_row.iter_mut().enumerate() {
            *c_ij = a_row
                .iter()
                .zip(b.chunks_exact(4))
                .map(|(a_ik, b_row)| a_ik * b_row[j])
                .sum();
        }
    }
    c
}
/// Row-major 4x4 matrix transpose.
pub(super) fn mat4_transpose(a: &[f64; 16]) -> [f64; 16] {
    let mut t = [0.0f64; 16];
    for (i, a_row) in a.chunks_exact(4).enumerate() {
        for (j, a_ij) in a_row.iter().enumerate() {
            t[j * 4 + i] = *a_ij;
        }
    }
    t
}
#[cfg(test)]
mod tests {

    use super::lcg_rand;
    use crate::autonomous_driving::AebSystem;
    use crate::autonomous_driving::AutonomousVehicle;
    use crate::autonomous_driving::BicycleModel;
    use crate::autonomous_driving::CameraDetector;
    use crate::autonomous_driving::ComfortEvent;
    use crate::autonomous_driving::ComfortMonitor;
    use crate::autonomous_driving::DubinsPath;
    use crate::autonomous_driving::DubinsSegment;
    use crate::autonomous_driving::GapAcceptance;
    use crate::autonomous_driving::KalmanTracker;
    use crate::autonomous_driving::LidarSensor;
    use crate::autonomous_driving::MpcPathTracker;
    use crate::autonomous_driving::ObjectClass;
    use crate::autonomous_driving::OccupancyGrid;
    use crate::autonomous_driving::ParkingPlanner;
    use crate::autonomous_driving::ParkingSlot;
    use crate::autonomous_driving::ParkingType;
    use crate::autonomous_driving::Pedestrian;
    use crate::autonomous_driving::PedestrianCrowd;
    use crate::autonomous_driving::PlanState;
    use crate::autonomous_driving::RadarSensor;
    use crate::autonomous_driving::RiskAssessor;
    use crate::autonomous_driving::RrtPlanner;
    use crate::autonomous_driving::SensorFusion;
    use crate::autonomous_driving::Vec2;
    use std::f64::consts::{PI, TAU};
    #[test]
    fn test_vec2_basic() {
        let a = Vec2::new(3.0, 4.0);
        assert!((a.norm() - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_vec2_dot() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        assert_eq!(a.dot(&b), 0.0);
    }
    #[test]
    fn test_vec2_cross() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        assert_eq!(a.cross(&b), 1.0);
    }
    #[test]
    fn test_vec2_normalize() {
        let a = Vec2::new(3.0, 4.0);
        let n = a.normalize();
        assert!((n.norm() - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_vec2_rotate() {
        let a = Vec2::new(1.0, 0.0);
        let r = a.rotate(std::f64::consts::FRAC_PI_2);
        assert!(r.x.abs() < 1e-9);
        assert!((r.y - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_lidar_sensor_creation() {
        let lidar = LidarSensor::automotive();
        assert_eq!(lidar.num_beams, 1800);
        assert!(lidar.max_range > 100.0);
    }
    #[test]
    fn test_lidar_scan_no_obstacles() {
        let lidar = LidarSensor::automotive();
        let scan = lidar.scan(&[]);
        assert_eq!(scan.len(), 1800);
        for (_, range) in &scan {
            assert!((*range - lidar.max_range).abs() < 1e-9);
        }
    }
    #[test]
    fn test_lidar_scan_with_obstacle() {
        let lidar = LidarSensor::automotive();
        let obs = (Vec2::new(10.0, 0.0), 1.0);
        let scan = lidar.scan(&[obs]);
        let front_beam = scan
            .iter()
            .find(|(a, _)| a.abs() < lidar.angular_resolution);
        if let Some((_, r)) = front_beam {
            assert!(*r < lidar.max_range);
        }
    }
    #[test]
    fn test_radar_detect_in_range() {
        let radar = RadarSensor::long_range();
        let pos = Vec2::new(50.0, 0.0);
        assert!(radar.detect(pos, -10.0).is_some());
    }
    #[test]
    fn test_radar_detect_out_of_range() {
        let radar = RadarSensor::long_range();
        let pos = Vec2::new(300.0, 0.0);
        assert!(radar.detect(pos, 0.0).is_none());
    }
    #[test]
    fn test_camera_project_in_front() {
        let cam = CameraDetector::forward_camera();
        let result = cam.project([0.0, 0.0, 10.0]);
        assert!(result.is_some());
    }
    #[test]
    fn test_camera_project_behind() {
        let cam = CameraDetector::forward_camera();
        let result = cam.project([0.0, 0.0, -1.0]);
        assert!(result.is_none());
    }
    #[test]
    fn test_kalman_tracker_init() {
        let t = KalmanTracker::new(0, 5.0, 3.0);
        assert!((t.x[0] - 5.0).abs() < 1e-9);
        assert!((t.x[1] - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_kalman_tracker_predict() {
        let mut t = KalmanTracker::new(0, 0.0, 0.0);
        t.x[2] = 1.0;
        t.predict(1.0);
        assert!((t.x[0] - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_kalman_tracker_update() {
        let mut t = KalmanTracker::new(0, 0.0, 0.0);
        t.update([10.0, 0.0]);
        assert!(t.x[0] > 0.0);
    }
    #[test]
    fn test_kalman_predict_builds_position_velocity_cross_covariance() {
        // P0 = 10 I, Q = diag(0.01, 0.01, 0.1, 0.1):  P' = F P F^T + Q dt.
        let mut t = KalmanTracker::new(0, 0.0, 0.0);
        let dt = 0.5;
        t.predict(dt);
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        assert!(
            close(t.p[0], 10.0 + dt * dt * 10.0 + 0.01 * dt),
            "{}",
            t.p[0]
        );
        assert!(
            close(t.p[5], 10.0 + dt * dt * 10.0 + 0.01 * dt),
            "{}",
            t.p[5]
        );
        assert!(close(t.p[10], 10.0 + 0.1 * dt), "{}", t.p[10]);
        assert!(close(t.p[15], 10.0 + 0.1 * dt), "{}", t.p[15]);
        // cov(x, vx) = dt * var(vx): previously exactly zero.
        assert!(close(t.p[2], dt * 10.0), "{}", t.p[2]);
        assert!(close(t.p[8], dt * 10.0), "{}", t.p[8]);
        assert!(close(t.p[7], dt * 10.0), "{}", t.p[7]);
        assert!(close(t.p[13], dt * 10.0), "{}", t.p[13]);
        // x and y stay uncorrelated.
        assert!(close(t.p[1], 0.0) && close(t.p[4], 0.0));
    }
    #[test]
    fn test_kalman_velocity_is_estimated_from_position_measurements() {
        let dt = 0.1;
        let mut t = KalmanTracker::new(7, 0.0, 0.0);
        for k in 1..=30 {
            t.predict(dt);
            t.update([2.0 * dt * f64::from(k), 0.0]);
        }
        let v = t.velocity();
        assert!((v.x - 2.0).abs() < 0.2, "vx = {}, expected 2.0", v.x);
        assert!(v.y.abs() < 0.2, "vy = {}, expected 0.0", v.y);
        // The estimated position follows the target.
        assert!((t.position().x - 6.0).abs() < 0.3, "x = {}", t.position().x);
    }
    #[test]
    fn test_kalman_velocity_estimated_for_diagonal_motion() {
        let dt = 0.1;
        let (vx, vy) = (-1.5, 3.0);
        let mut t = KalmanTracker::new(0, 5.0, -2.0);
        for k in 1..=40 {
            t.predict(dt);
            let time = dt * f64::from(k);
            t.update([5.0 + vx * time, -2.0 + vy * time]);
        }
        let v = t.velocity();
        assert!((v.x - vx).abs() < 0.2, "vx = {}", v.x);
        assert!((v.y - vy).abs() < 0.2, "vy = {}", v.y);
    }
    #[test]
    fn test_kalman_covariance_stays_symmetric_and_shrinks() {
        let dt = 0.1;
        let mut t = KalmanTracker::new(0, 0.0, 0.0);
        let initial_vel_var = t.p[10];
        for k in 1..=30 {
            t.predict(dt);
            t.update([2.0 * dt * f64::from(k), 0.0]);
        }
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (t.p[i * 4 + j] - t.p[j * 4 + i]).abs() < 1e-12,
                    "P not symmetric at ({i}, {j})"
                );
            }
            assert!(t.p[i * 4 + i] > 0.0, "variance {i} must stay positive");
        }
        assert!(
            t.p[10] < 0.1 * initial_vel_var,
            "velocity variance should collapse once velocity is observed: {}",
            t.p[10]
        );
    }
    #[test]
    fn test_sensor_fusion_track_reports_nonzero_velocity() {
        let dt = 0.1;
        let mut fusion = SensorFusion::new();
        let mut last = Vec::new();
        for k in 0..30 {
            let x = 5.0 + 2.0 * dt * f64::from(k);
            last = fusion.update(&[Vec2::new(x, -3.0)], dt);
        }
        assert_eq!(last.len(), 1, "one confirmed track expected");
        let (_, pos, vel) = last[0];
        assert!((vel.x - 2.0).abs() < 0.3, "vx = {}", vel.x);
        assert!(vel.y.abs() < 0.3, "vy = {}", vel.y);
        assert!((pos.y + 3.0).abs() < 0.3, "y = {}", pos.y);
    }
    #[test]
    fn test_sensor_fusion_empty() {
        let mut fusion = SensorFusion::new();
        let tracks = fusion.update(&[], 0.1);
        assert!(tracks.is_empty());
    }
    #[test]
    fn test_sensor_fusion_creates_tracks() {
        let mut fusion = SensorFusion::new();
        let dets = vec![Vec2::new(10.0, 0.0), Vec2::new(-5.0, 3.0)];
        fusion.update(&dets, 0.1);
        fusion.update(&dets, 0.1);
        let tracks = fusion.update(&dets, 0.1);
        assert!(!tracks.is_empty());
    }
    #[test]
    fn test_occupancy_grid_creation() {
        let grid = OccupancyGrid::new(100, 100, 0.5);
        assert_eq!(grid.log_odds.len(), 10000);
    }
    #[test]
    fn test_occupancy_grid_update() {
        let mut grid = OccupancyGrid::new(100, 100, 0.5);
        grid.update_cell(50, 50, 0.9);
        assert!(grid.probability(50, 50) > 0.5);
    }
    #[test]
    fn test_occupancy_grid_world_to_cell() {
        let grid = OccupancyGrid::new(100, 100, 0.5);
        let cell = grid.world_to_cell(Vec2::new(10.0, 10.0));
        assert!(cell.is_some());
    }
    #[test]
    fn test_bicycle_model_straight() {
        let mut bike = BicycleModel::new(2.7);
        bike.state[3] = 10.0;
        bike.step(0.0, 0.0, 1.0);
        assert!((bike.state[0] - 10.0).abs() < 1e-6);
    }
    #[test]
    fn test_bicycle_model_turn() {
        let mut bike = BicycleModel::new(2.7);
        bike.state[3] = 10.0;
        bike.step(0.3, 0.0, 1.0);
        assert!(bike.state[2].abs() > 0.0);
    }
    #[test]
    fn test_mpc_tracker_no_panic() {
        let mut tracker = MpcPathTracker::new(2.7, 10, 0.05);
        tracker.model.state[3] = 10.0;
        let path = vec![
            Vec2::new(5.0, 0.0),
            Vec2::new(10.0, 1.0),
            Vec2::new(20.0, 2.0),
        ];
        tracker.step(&path, 10.0);
    }
    #[test]
    fn test_pure_pursuit_returns_bounded_steer() {
        let tracker = MpcPathTracker::new(2.7, 10, 0.05);
        let path = vec![Vec2::new(5.0, 1.0), Vec2::new(10.0, 2.0)];
        let steer = tracker.pure_pursuit_steer(&path, 3.0);
        assert!(steer.abs() <= tracker.model.max_steer + 1e-9);
    }
    #[test]
    fn test_risk_assessor_ttc_infinite() {
        let ra = RiskAssessor::new();
        let ttc = ra.time_to_collision(10.0, 10.0, 50.0);
        assert_eq!(ttc, f64::INFINITY);
    }
    #[test]
    fn test_risk_assessor_ttc_finite() {
        let ra = RiskAssessor::new();
        let ttc = ra.time_to_collision(20.0, 10.0, 100.0);
        assert!((ttc - 10.0).abs() < 1e-6);
    }
    #[test]
    fn test_aeb_full_brake() {
        let mut aeb = AebSystem::new();
        let decel = aeb.update(0.5, ObjectClass::Vehicle);
        assert!(decel > 0.0);
        assert!(aeb.active);
    }
    #[test]
    fn test_aeb_no_brake() {
        let mut aeb = AebSystem::new();
        let decel = aeb.update(10.0, ObjectClass::Vehicle);
        assert_eq!(decel, 0.0);
        assert!(!aeb.active);
    }
    #[test]
    fn test_comfort_monitor_jerk() {
        let mut cm = ComfortMonitor::new(0.1);
        cm.update(0.0, 0.0);
        let event = cm.update(100.0, 0.0);
        assert_ne!(event, ComfortEvent::Ok);
    }
    #[test]
    fn test_gap_acceptance_safe() {
        let ga = GapAcceptance::new(0.5);
        assert!(ga.is_acceptable(50.0, 20.0, 30.0, 5.0, 15.0));
    }
    #[test]
    fn test_gap_acceptance_unsafe() {
        let ga = GapAcceptance::new(0.5);
        assert!(!ga.is_acceptable(0.5, 20.0, 0.5, 25.0, 15.0));
    }
    #[test]
    fn test_pedestrian_step() {
        let mut ped = Pedestrian::new(0, Vec2::new(0.0, 0.0), Vec2::new(10.0, 0.0));
        let others = vec![];
        let obstacles = vec![];
        ped.step(&others, &obstacles, 0.1);
        assert!(ped.pos.x > 0.0);
    }
    #[test]
    fn test_pedestrian_crowd_step() {
        let mut crowd = PedestrianCrowd::new();
        crowd.agents.push(Pedestrian::new(
            0,
            Vec2::new(0.0, 0.0),
            Vec2::new(10.0, 0.0),
        ));
        crowd.agents.push(Pedestrian::new(
            1,
            Vec2::new(1.0, 0.0),
            Vec2::new(10.0, 0.0),
        ));
        crowd.step(0.1);
    }
    #[test]
    fn test_parking_planner_can_park() {
        let planner = ParkingPlanner::new(2.7, 2.0);
        let slot = ParkingSlot {
            center: Vec2::new(10.0, 0.0),
            heading: 0.0,
            length: 8.0,
            width: 3.0,
            slot_type: ParkingType::Parallel,
        };
        assert!(planner.can_parallel_park(&slot));
    }
    #[test]
    fn test_parking_planner_too_small() {
        let planner = ParkingPlanner::new(2.7, 2.0);
        let slot = ParkingSlot {
            center: Vec2::new(10.0, 0.0),
            heading: 0.0,
            length: 4.0,
            width: 2.1,
            slot_type: ParkingType::Parallel,
        };
        assert!(!planner.can_parallel_park(&slot));
    }
    #[test]
    fn test_dubins_path_sample() {
        let path = DubinsPath {
            segments: [
                DubinsSegment::Straight,
                DubinsSegment::Left,
                DubinsSegment::Right,
            ],
            lengths: [10.0, 1.0, 0.5],
            radius: 5.0,
            start: [0.0, 0.0, 0.0],
            end: [10.0, 2.0, 0.3],
        };
        let pt = path.sample(5.0);
        assert!((pt[0] - 5.0).abs() < 1e-6);
    }
    /// Smallest absolute angle difference.
    fn angle_diff(a: f64, b: f64) -> f64 {
        ((a - b + PI).rem_euclid(TAU) - PI).abs()
    }
    /// Assert that `path` ends exactly at `end` (position and heading).
    fn assert_reaches_end(path: &DubinsPath, end: [f64; 3], what: &str) {
        let pose = path.sample(path.total_length());
        assert!(
            (pose[0] - end[0]).abs() < 1e-6 && (pose[1] - end[1]).abs() < 1e-6,
            "{what}: position {:?} != {:?}",
            &pose[..2],
            &end[..2]
        );
        assert!(
            angle_diff(pose[2], end[2]) < 1e-6,
            "{what}: heading {} != {}",
            pose[2],
            end[2]
        );
    }
    #[test]
    fn test_dubins_every_word_reaches_the_end_pose() {
        use std::collections::HashSet;
        let mut rng = 20_240_611_u64;
        let mut seen_words: HashSet<Vec<u8>> = HashSet::new();
        for trial in 0..3000 {
            let radius = 0.5 + 9.5 * lcg_rand(&mut rng);
            // Mix near and far goals so that the CCC words also become optimal.
            let reach = if trial % 2 == 0 { 4.0 * radius } else { 60.0 };
            let start = [
                reach * (2.0 * lcg_rand(&mut rng) - 1.0),
                reach * (2.0 * lcg_rand(&mut rng) - 1.0),
                TAU * lcg_rand(&mut rng) - PI,
            ];
            let end = [
                start[0] + reach * (2.0 * lcg_rand(&mut rng) - 1.0),
                start[1] + reach * (2.0 * lcg_rand(&mut rng) - 1.0),
                TAU * lcg_rand(&mut rng) - PI,
            ];
            let candidates = DubinsPath::candidates(start, end, radius);
            assert!(
                !candidates.is_empty(),
                "trial {trial}: Dubins' theorem guarantees a feasible word"
            );
            for path in &candidates {
                assert_reaches_end(path, end, &format!("trial {trial} {:?}", path.segments));
                assert!(
                    path.lengths.iter().all(|&l| l >= 0.0),
                    "negative segment length {:?}",
                    path.lengths
                );
                seen_words.insert(
                    path.segments
                        .iter()
                        .map(|s| match s {
                            DubinsSegment::Left => b'L',
                            DubinsSegment::Straight => b'S',
                            DubinsSegment::Right => b'R',
                        })
                        .collect(),
                );
            }
            let best = DubinsPath::shortest(start, end, radius).expect("path exists");
            let min_len = candidates
                .iter()
                .map(DubinsPath::total_length)
                .fold(f64::INFINITY, f64::min);
            assert!((best.total_length() - min_len).abs() < 1e-12);
            let chord = ((end[0] - start[0]).powi(2) + (end[1] - start[1]).powi(2)).sqrt();
            assert!(
                best.total_length() >= chord - 1e-9,
                "a path cannot be shorter than the straight line"
            );
        }
        // All six words were produced (and verified) at least once.
        for word in ["LSL", "RSR", "LSR", "RSL", "RLR", "LRL"] {
            assert!(
                seen_words.contains(word.as_bytes()),
                "word {word} was never feasible in the sample"
            );
        }
    }
    #[test]
    fn test_dubins_aligned_configurations_give_the_straight_line() {
        let path = DubinsPath::shortest([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 2.0).expect("path");
        assert!((path.total_length() - 10.0).abs() < 1e-9);
        assert_reaches_end(&path, [10.0, 0.0, 0.0], "straight");
        // Heading along the direction of travel, 45 degrees.
        let end = [3.0, 3.0, PI / 4.0];
        let path = DubinsPath::shortest([0.0, 0.0, PI / 4.0], end, 1.5).expect("path");
        assert!((path.total_length() - 18.0_f64.sqrt()).abs() < 1e-9);
        assert_reaches_end(&path, end, "diagonal straight");
    }
    #[test]
    fn test_dubins_u_turn_is_a_half_circle() {
        let radius = 3.0;
        let end = [0.0, 2.0 * radius, PI];
        let path = DubinsPath::shortest([0.0, 0.0, 0.0], end, radius).expect("path");
        assert!(
            (path.total_length() - PI * radius).abs() < 1e-9,
            "half circle expected, got {}",
            path.total_length()
        );
        assert_reaches_end(&path, end, "u-turn");
    }
    #[test]
    fn test_dubins_end_heading_is_honoured() {
        // The old implementation always returned a straight line towards the goal
        // and therefore ended with the start heading.
        let end = [10.0, 5.0, PI];
        let path = DubinsPath::shortest([0.0, 0.0, 0.0], end, 2.0).expect("path");
        assert_reaches_end(&path, end, "reverse heading");
        assert!(path.total_length() > (125.0_f64).sqrt());
    }
    #[test]
    fn test_dubins_same_position_different_heading_gives_a_loop() {
        let start = [1.0, 2.0, 0.0];
        for heading in [0.3, PI / 2.0, PI, -2.0, 6.0] {
            let end = [1.0, 2.0, heading];
            let path = DubinsPath::shortest(start, end, 2.0).expect("a loop path must exist");
            assert!(path.total_length() > 0.0);
            assert_reaches_end(&path, end, &format!("heading {heading}"));
        }
    }
    #[test]
    fn test_dubins_identical_configurations_and_bad_radius_give_none() {
        let pose = [3.0, -4.0, 1.0];
        assert!(DubinsPath::shortest(pose, pose, 2.0).is_none());
        // Equal heading modulo a full turn is the same configuration.
        assert!(DubinsPath::shortest(pose, [3.0, -4.0, 1.0 + TAU], 2.0).is_none());
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(DubinsPath::shortest([0.0; 3], [5.0, 5.0, 0.0], radius).is_none());
        }
        assert!(DubinsPath::shortest([0.0; 3], [f64::NAN, 5.0, 0.0], 1.0).is_none());
    }
    #[test]
    fn test_dubins_lengths_scale_with_radius() {
        let a = DubinsPath::shortest([0.0, 0.0, 0.3], [8.0, -3.0, 2.0], 1.0).expect("path");
        let b = DubinsPath::shortest([0.0, 0.0, 0.3], [16.0, -6.0, 2.0], 2.0).expect("path");
        assert!((b.total_length() - 2.0 * a.total_length()).abs() < 1e-9);
        assert_eq!(a.segments, b.segments);
    }
    #[test]
    fn test_autonomous_vehicle_step() {
        let mut av = AutonomousVehicle::new();
        av.tracker.model.state[3] = 10.0;
        av.planned_path = vec![
            Vec2::new(5.0, 0.0),
            Vec2::new(10.0, 0.0),
            Vec2::new(20.0, 0.0),
        ];
        let scan = vec![(0.0_f64, 100.0_f64)];
        let action = av.step(&[], &scan, 0.05);
        assert!(action.brake_demand >= 0.0);
    }
    #[test]
    fn test_rrt_planner_finds_path() {
        let planner = RrtPlanner::new(1.0, [-20.0, -20.0, 20.0, 20.0]);
        let start = PlanState::new(0.0, 0.0, 0.0, 0.0);
        let goal = PlanState::new(5.0, 5.0, 0.0, 0.0);
        let path = planner.plan(start, goal, &|_| false);
        assert!(path.is_some());
    }
}
