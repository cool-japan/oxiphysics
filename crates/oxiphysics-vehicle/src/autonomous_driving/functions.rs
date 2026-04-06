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
#[cfg(test)]
mod tests {

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
