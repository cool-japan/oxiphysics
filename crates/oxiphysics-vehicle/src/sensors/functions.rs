//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types_core::*;
use std::f64::consts::PI;

/// Radar range equation (Friis / monostatic radar).
///
/// `Pr = Pt * Gt * Gr * lambda² * sigma / ((4π)³ * R⁴)`
pub fn radar_range_equation(pt: f64, gt: f64, gr: f64, lambda: f64, sigma: f64, range: f64) -> f64 {
    let four_pi_cubed = (4.0 * PI).powi(3);
    pt * gt * gr * lambda * lambda * sigma / (four_pi_cubed * range.powi(4))
}
/// Compute the centroid of a point cloud.
pub fn point_cloud_centroid(points: &[LidarPoint]) -> [f64; 3] {
    if points.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let n = points.len() as f64;
    let sum = points.iter().fold([0.0_f64; 3], |acc, p| {
        [acc[0] + p.x, acc[1] + p.y, acc[2] + p.z]
    });
    [sum[0] / n, sum[1] / n, sum[2] / n]
}
/// Compute the axis-aligned bounding box of a point cloud.
///
/// Returns `(min_corner, max_corner)`.
pub fn point_cloud_bounding_box(points: &[LidarPoint]) -> ([f64; 3], [f64; 3]) {
    if points.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in points {
        let coords = [p.x, p.y, p.z];
        for i in 0..3 {
            if coords[i] < min[i] {
                min[i] = coords[i];
            }
            if coords[i] > max[i] {
                max[i] = coords[i];
            }
        }
    }
    (min, max)
}
/// Compute the number of cells per side for a square bird's-eye-view grid.
///
/// `range_m` is the half-extent of the grid; `resolution_m` is metres per cell.
pub fn bev_grid_size(range_m: f64, resolution_m: f64) -> usize {
    ((2.0 * range_m) / resolution_m).ceil() as usize
}
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[cfg(test)]
mod tests {
    use super::super::types_advanced::*;
    use super::*;
    use rand::SeedableRng;
    #[test]
    fn test_sensor_noise_zero_noise() {
        let n = SensorNoise::new(0.0, 0.0);
        assert_eq!(n.apply(5.0, 1.0), 5.0);
    }
    #[test]
    fn test_sensor_noise_bias_only() {
        let n = SensorNoise::new(0.0, 0.5);
        assert!((n.apply(3.0, 0.0) - 3.5).abs() < 1e-10);
    }
    #[test]
    fn test_sensor_noise_std_dev() {
        let n = SensorNoise::new(2.0, 0.0);
        assert!((n.apply(0.0, 1.0) - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_sensor_noise_combined() {
        let n = SensorNoise::new(1.0, 0.3);
        assert!((n.apply(10.0, 2.0) - 12.3).abs() < 1e-10);
    }
    #[test]
    fn test_radar_automotive_config() {
        let cfg = RadarConfig::automotive();
        assert_eq!(cfg.max_range, 200.0);
        assert_eq!(cfg.beam_width_deg, 6.0);
        assert_eq!(cfg.update_rate_hz, 20.0);
    }
    #[test]
    fn test_radar_detect_in_range() {
        let cfg = RadarConfig::automotive();
        let result = RadarSensor::detect(&cfg, 0.0, [50.0, 0.0], [0.0, 0.0], [0.0, 0.0], 0.0);
        assert!(result.is_some());
        let t = result.unwrap();
        assert!((t.range - 50.0).abs() < 0.5);
    }
    #[test]
    fn test_radar_detect_out_of_range() {
        let cfg = RadarConfig::automotive();
        let result = RadarSensor::detect(&cfg, 0.0, [300.0, 0.0], [0.0, 0.0], [0.0, 0.0], 0.0);
        assert!(result.is_none());
    }
    #[test]
    fn test_radar_detect_outside_beam() {
        let cfg = RadarConfig::automotive();
        let result = RadarSensor::detect(&cfg, 0.0, [50.0, 50.0], [0.0, 0.0], [0.0, 0.0], 0.0);
        assert!(result.is_none());
    }
    #[test]
    fn test_radar_detect_azimuth_value() {
        let cfg = RadarConfig::automotive();
        let result = RadarSensor::detect(&cfg, 0.0, [100.0, 0.0], [0.0, 0.0], [0.0, 0.0], 0.0);
        assert!(result.is_some());
        assert!(result.unwrap().azimuth.abs() < 1.0);
    }
    #[test]
    fn test_radar_doppler_approaching() {
        let los = [1.0, 0.0];
        let v = RadarSensor::doppler_velocity(0.0, [-20.0, 0.0], los);
        assert!((v - (-20.0)).abs() < 1e-10);
    }
    #[test]
    fn test_radar_doppler_receding() {
        let los = [1.0, 0.0];
        let v = RadarSensor::doppler_velocity(0.0, [30.0, 0.0], los);
        assert!((v - 30.0).abs() < 1e-10);
    }
    #[test]
    fn test_radar_range_equation_positive() {
        let pr = radar_range_equation(1000.0, 30.0, 30.0, 0.039, 1.0, 100.0);
        assert!(pr > 0.0);
    }
    #[test]
    fn test_radar_range_equation_r4_dependence() {
        let pr1 = radar_range_equation(1.0, 1.0, 1.0, 1.0, 1.0, 10.0);
        let pr2 = radar_range_equation(1.0, 1.0, 1.0, 1.0, 1.0, 20.0);
        let ratio = pr1 / pr2;
        assert!((ratio - 16.0).abs() < 1e-6);
    }
    #[test]
    fn test_lidar_vlp16_config() {
        let cfg = LidarConfig::velodyne_vlp16();
        assert_eq!(cfg.n_beams, 16);
        assert_eq!(cfg.max_range, 100.0);
        assert_eq!(cfg.vertical_fov_deg, 30.0);
    }
    #[test]
    fn test_ray_sphere_hit() {
        let origin = [0.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let center = [5.0, 0.0, 0.0];
        let t = LidarScan::ray_sphere_intersect(origin, dir, center, 1.0);
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-9);
    }
    #[test]
    fn test_ray_sphere_miss() {
        let origin = [0.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let center = [0.0, 10.0, 0.0];
        let t = LidarScan::ray_sphere_intersect(origin, dir, center, 1.0);
        assert!(t.is_none());
    }
    #[test]
    fn test_ray_sphere_behind_origin() {
        let origin = [0.0, 0.0, 0.0];
        let dir = [1.0, 0.0, 0.0];
        let center = [-5.0, 0.0, 0.0];
        let t = LidarScan::ray_sphere_intersect(origin, dir, center, 1.0);
        assert!(t.is_none());
    }
    #[test]
    fn test_generate_scan_hits_sphere() {
        let cfg = LidarConfig::velodyne_vlp16();
        let origin = [0.0, 0.0, 0.0];
        let obstacles = [([10.0, 0.0, 0.0], 1.0)];
        let pts = LidarScan::generate_scan(&cfg, origin, &obstacles);
        assert!(!pts.is_empty());
    }
    #[test]
    fn test_generate_scan_no_obstacles() {
        let cfg = LidarConfig::velodyne_vlp16();
        let origin = [0.0, 0.0, 0.0];
        let pts = LidarScan::generate_scan(&cfg, origin, &[]);
        assert!(pts.is_empty());
    }
    #[test]
    fn test_centroid_empty() {
        let c = point_cloud_centroid(&[]);
        assert_eq!(c, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_centroid_single_point() {
        let pts = vec![LidarPoint {
            x: 1.0,
            y: 2.0,
            z: 3.0,
            intensity: 1.0,
            ring: 0,
        }];
        let c = point_cloud_centroid(&pts);
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_bounding_box_empty() {
        let (mn, mx) = point_cloud_bounding_box(&[]);
        assert_eq!(mn, [0.0; 3]);
        assert_eq!(mx, [0.0; 3]);
    }
    #[test]
    fn test_bounding_box_two_points() {
        let pts = vec![
            LidarPoint {
                x: -1.0,
                y: 0.0,
                z: 0.0,
                intensity: 1.0,
                ring: 0,
            },
            LidarPoint {
                x: 1.0,
                y: 2.0,
                z: -1.0,
                intensity: 1.0,
                ring: 1,
            },
        ];
        let (mn, mx) = point_cloud_bounding_box(&pts);
        assert!((mn[0] - (-1.0)).abs() < 1e-10);
        assert!((mx[0] - 1.0).abs() < 1e-10);
        assert!((mx[1] - 2.0).abs() < 1e-10);
        assert!((mn[2] - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_pinhole_new_focal_length() {
        let cam = PinholeCamera::new(90.0, 640, 480);
        assert!((cam.fx - 320.0).abs() < 1.0);
    }
    #[test]
    fn test_pinhole_project_centre() {
        let cam = PinholeCamera::new(90.0, 640, 480);
        let px = cam.project([0.0, 0.0, 5.0]);
        assert!(px.is_some());
        let px = px.unwrap();
        assert!((px[0] - cam.cx).abs() < 1e-6);
        assert!((px[1] - cam.cy).abs() < 1e-6);
    }
    #[test]
    fn test_pinhole_project_behind_camera() {
        let cam = PinholeCamera::new(90.0, 640, 480);
        assert!(cam.project([0.0, 0.0, -1.0]).is_none());
    }
    #[test]
    fn test_pinhole_unproject_roundtrip() {
        let cam = PinholeCamera::new(60.0, 1280, 720);
        let pt3d = [1.5, -0.5, 4.0];
        if let Some(px) = cam.project(pt3d) {
            let back = cam.unproject(px, pt3d[2]);
            assert!((back[0] - pt3d[0]).abs() < 1e-9);
            assert!((back[1] - pt3d[1]).abs() < 1e-9);
            assert!((back[2] - pt3d[2]).abs() < 1e-9);
        }
    }
    #[test]
    fn test_bev_grid_size() {
        assert_eq!(bev_grid_size(50.0, 0.1), 1000);
    }
    #[test]
    fn test_bev_grid_size_ceil() {
        assert_eq!(bev_grid_size(5.0, 0.3), 34);
    }
    #[test]
    fn test_gps_automotive_config() {
        let cfg = GpsConfig::automotive();
        assert!(cfg.horizontal_accuracy_m > 0.0);
        assert_eq!(cfg.update_rate_hz, 10.0);
    }
    #[test]
    fn test_gps_measure_zero_noise_matches_true_pos() {
        let cfg = GpsConfig::rtk();
        let true_pos = [100.0, 200.0, 5.0];
        let true_vel = [10.0, 0.0, 0.0];
        let m = GpsSensor::measure(&cfg, true_pos, true_vel, 0.0, 0.0, 12);
        assert!((m.position[0] - 100.0).abs() < 1e-9);
        assert!((m.position[2] - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_gps_hdop_increases_with_fewer_satellites() {
        let cfg = GpsConfig::automotive();
        let m_good = GpsSensor::measure(&cfg, [0.0; 3], [0.0; 3], 0.0, 0.0, 12);
        let m_bad = GpsSensor::measure(&cfg, [0.0; 3], [0.0; 3], 0.0, 0.0, 4);
        assert!(m_bad.hdop > m_good.hdop, "fewer satellites → higher HDOP");
    }
    #[test]
    fn test_imu_measure_zero_noise_matches_true() {
        let cfg = ImuConfig::mems_automotive();
        let a = [0.0, 0.0, 9.81];
        let g = [0.0, 0.0, 0.1];
        let m = ImuSensor::measure(&cfg, a, g, 0.0);
        assert!((m.accel[2] - 9.81).abs() < 1e-9);
        assert!((m.gyro[2] - 0.1).abs() < 1e-9);
    }
    #[test]
    fn test_imu_accel_attitude_level() {
        let (roll, pitch) = ImuSensor::accel_attitude([0.0, 0.0, 9.81]);
        assert!(roll.abs() < 1e-9, "roll should be ~0, got {roll}");
        assert!(pitch.abs() < 1e-9, "pitch should be ~0, got {pitch}");
    }
    #[test]
    fn test_imu_accel_attitude_roll_90_deg() {
        let (roll, _pitch) = ImuSensor::accel_attitude([0.0, 9.81, 0.0]);
        assert!(
            (roll.abs() - std::f64::consts::FRAC_PI_2).abs() < 1e-6,
            "expected roll=π/2, got {roll}"
        );
    }
    #[test]
    fn test_wheel_speed_zero_noise() {
        let cfg = WheelSpeedConfig::abs_hall_effect();
        let omega = [10.0, 10.0, 10.0, 10.0];
        let m = WheelSpeedSensor::measure(&cfg, omega, 0.0);
        for i in 0..4 {
            assert!((m.angular_velocity[i] - 10.0).abs() < 1e-9);
            assert!((m.linear_speed[i] - 10.0 * cfg.wheel_radius).abs() < 1e-9);
        }
    }
    #[test]
    fn test_wheel_speed_linear_speed_proportional_to_radius() {
        let cfg = WheelSpeedConfig {
            wheel_radius: 0.5,
            ..WheelSpeedConfig::abs_hall_effect()
        };
        let m = WheelSpeedSensor::measure(&cfg, [4.0; 4], 0.0);
        for i in 0..4 {
            assert!((m.linear_speed[i] - 2.0).abs() < 1e-9);
        }
    }
    #[test]
    fn test_steering_angle_sensor_zero_noise() {
        let s = SteeringAngleSensor::optical_encoder();
        let angle = s.measure(0.1, 0.0);
        assert!((angle - 0.1).abs() < 1e-9);
    }
    #[test]
    fn test_steering_angle_sensor_bias() {
        let s = SteeringAngleSensor {
            noise_rad: 0.0,
            bias_rad: 0.01,
        };
        let angle = s.measure(0.0, 0.0);
        assert!((angle - 0.01).abs() < 1e-12);
    }
    #[test]
    fn test_pedal_sensor_throttle_clamped() {
        let s = PedalSensor::hall_effect();
        assert_eq!(s.measure_throttle(1.1, 0.0), 1.0);
        assert_eq!(s.measure_throttle(-0.1, 0.0), 0.0);
    }
    #[test]
    fn test_pedal_sensor_brake_clamped() {
        let s = PedalSensor::hall_effect();
        assert_eq!(s.measure_brake(1.2, 0.0), 1.0);
        assert_eq!(s.measure_brake(-0.5, 0.0), 0.0);
    }
    #[test]
    fn test_sensor_fusion_init_zero() {
        let sf = SensorFusion::new(0.98);
        assert_eq!(sf.roll, 0.0);
        assert_eq!(sf.pitch, 0.0);
    }
    #[test]
    fn test_sensor_fusion_update_level_vehicle() {
        let mut sf = SensorFusion::new(0.98);
        let g = 9.81;
        for _ in 0..100 {
            sf.update([0.0, 0.0, g], [0.0, 0.0, 0.0], 0.01);
        }
        assert!(
            sf.roll.abs() < 1e-6,
            "roll should converge to 0, got {}",
            sf.roll
        );
        assert!(
            sf.pitch.abs() < 0.01,
            "pitch should converge to 0, got {}",
            sf.pitch
        );
    }
    #[test]
    fn test_sensor_fusion_reset() {
        let mut sf = SensorFusion::new(0.98);
        sf.roll = 0.5;
        sf.pitch = 0.3;
        sf.reset();
        assert_eq!(sf.roll, 0.0);
        assert_eq!(sf.pitch, 0.0);
    }
    #[test]
    fn test_speed_sensor_adds_noise() {
        let sensor = SpeedSensor::new(0.31, 0.5);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(42);
        let omega = 50.0;
        let true_speed = omega * 0.31;
        let measured = sensor.measure(omega, &mut rng);
        assert!(
            (measured - true_speed).abs() < 5.0,
            "measured speed should be close to true, got {measured}"
        );
        let measured2 = sensor.measure(omega, &mut rng);
        let measured3 = sensor.measure(omega, &mut rng);
        let all_exact = (measured - true_speed).abs() < 1e-12
            && (measured2 - true_speed).abs() < 1e-12
            && (measured3 - true_speed).abs() < 1e-12;
        assert!(!all_exact, "noise should perturb the measurement");
    }
    #[test]
    fn test_gps_noise_within_3sigma() {
        let gps = GpsUnit::new(2.0);
        let mut rng = rand::rngs::SmallRng::seed_from_u64(99);
        let true_pos = [100.0, 200.0, 5.0];
        let mut max_err = 0.0f64;
        for _ in 0..100 {
            let m = gps.measure_position(true_pos, &mut rng);
            for i in 0..3 {
                let err = (m[i] - true_pos[i]).abs();
                if err > max_err {
                    max_err = err;
                }
            }
        }
        assert!(
            max_err < 6.0 * 3.0,
            "GPS noise too large: max_err={max_err}"
        );
    }
    #[test]
    fn test_lidar_ray_hits_plane() {
        let ray = LidarRay::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0], 100.0);
        let hit = ray.cast_against_plane([0.0, 0.0, 1.0], 0.0);
        assert!(hit.is_some(), "ray should hit z=0 plane");
        assert!(
            (hit.unwrap() - 5.0).abs() < 1e-9,
            "hit distance should be 5.0"
        );
    }
    #[test]
    fn test_lidar_ray_misses_behind_plane() {
        let ray = LidarRay::new([0.0, 0.0, -1.0], [0.0, 0.0, -1.0], 100.0);
        let hit = ray.cast_against_plane([0.0, 0.0, 1.0], 0.0);
        assert!(hit.is_none(), "ray pointing away should not hit plane");
    }
    #[test]
    fn test_lidar_ray_exceeds_max_range() {
        let ray = LidarRay::new([0.0, 0.0, 200.0], [0.0, 0.0, -1.0], 10.0);
        let hit = ray.cast_against_plane([0.0, 0.0, 1.0], 0.0);
        assert!(hit.is_none(), "distance 200 exceeds max_range 10");
    }
    #[test]
    fn test_telemetry_ring_buffer_wraps() {
        let mut buf = TelemetryBuffer::new(4);
        for i in 0..10u32 {
            buf.push(i as f64, i as f64 * 2.0);
        }
        assert_eq!(buf.len(), 4, "buffer should hold exactly 4 entries");
        let last2 = buf.last_n(2);
        assert_eq!(last2.len(), 2);
    }
    #[test]
    fn test_telemetry_ring_buffer_push_and_retrieve() {
        let mut buf = TelemetryBuffer::new(8);
        assert!(buf.is_empty());
        buf.push(1.0, 10.0);
        buf.push(2.0, 20.0);
        buf.push(3.0, 30.0);
        assert_eq!(buf.len(), 3);
        let last = buf.last_n(3);
        assert_eq!(last.len(), 3);
    }
    #[test]
    fn test_wheel_pulse_count() {
        let wp = WheelPulse::new(48, 0.0);
        let pulses = wp.pulse_count(2.0 * PI, 1.0);
        assert_eq!(pulses, 48);
    }
    #[test]
    fn test_wheel_pulse_velocity_from_pulses() {
        let wp = WheelPulse::new(48, 0.0);
        let circumference = 2.0 * PI * 0.31;
        let v = wp.velocity_from_pulses(48, 1.0, circumference);
        assert!((v - circumference).abs() < 1e-9, "got {v}");
    }
    #[test]
    fn test_imu_unit_bias_applied() {
        let mut rng = rand::rngs::SmallRng::seed_from_u64(0);
        let imu = ImuUnit::new(0.0, 0.0, [1.0, 2.0, 3.0]);
        let meas = imu.measure_accel([0.0, 0.0, 0.0], &mut rng);
        assert!((meas[0] - 1.0).abs() < 1e-9);
        assert!((meas[1] - 2.0).abs() < 1e-9);
        assert!((meas[2] - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_position_fusion_gps_weight_1_follows_gps() {
        let mut fuse = PositionFusion::new(1.0);
        let gps = [10.0, 20.0, 5.0];
        let result = fuse.update(gps, [0.0, 0.0, 0.0], 0.01);
        assert!((result[0] - 10.0).abs() < 1e-9);
        assert!((result[1] - 20.0).abs() < 1e-9);
    }
    #[test]
    fn test_barometric_altitude_zero_height() {
        let s = BarometricSensor::standard();
        let p = s.measure_pressure(0.0, 0.0);
        assert!((p - s.p0).abs() < 1e-6, "at h=0 pressure should equal p0");
    }
    #[test]
    fn test_barometric_altitude_roundtrip() {
        let s = BarometricSensor::standard();
        let h_true = 500.0;
        let p = s.measure_pressure(h_true, 0.0);
        let h_est = s.altitude_from_pressure(p);
        assert!(
            (h_est - h_true).abs() < 0.1,
            "roundtrip altitude error: {}",
            (h_est - h_true).abs()
        );
    }
    #[test]
    fn test_barometric_higher_altitude_lower_pressure() {
        let s = BarometricSensor::standard();
        let p_low = s.measure_pressure(100.0, 0.0);
        let p_high = s.measure_pressure(500.0, 0.0);
        assert!(p_high < p_low, "higher altitude → lower pressure");
    }
    #[test]
    fn test_temperature_sensor_zero_noise() {
        let s = TemperatureSensor {
            noise_c: 0.0,
            bias_c: 0.0,
            update_rate_hz: 10.0,
        };
        assert!((s.measure(100.0, 0.0) - 100.0).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_sensor_bias_applied() {
        let s = TemperatureSensor {
            noise_c: 0.0,
            bias_c: 2.5,
            update_rate_hz: 10.0,
        };
        assert!((s.measure(50.0, 0.0) - 52.5).abs() < 1e-12);
    }
    #[test]
    fn test_temperature_sensor_rtd_preset() {
        let s = TemperatureSensor::rtd_tyre();
        assert!((s.bias_c).abs() < 1e-9, "RTD should have zero bias");
        assert!(s.noise_c < 0.5, "RTD should be low-noise");
    }
    #[test]
    fn test_kalman1d_predict_advances_position() {
        let mut kf = KalmanFilter1D::new(0.01, 0.01, 1.0);
        kf.position = 0.0;
        kf.velocity = 5.0;
        kf.predict(1.0);
        assert!(
            (kf.position - 5.0).abs() < 1e-9,
            "position should advance by v*dt"
        );
    }
    #[test]
    fn test_kalman1d_update_converges_to_measurement() {
        let mut kf = KalmanFilter1D::new(0.01, 0.01, 0.001);
        kf.position = 0.0;
        for _ in 0..50 {
            kf.update(10.0);
        }
        assert!(
            (kf.position - 10.0).abs() < 0.5,
            "should converge to measured position"
        );
    }
    #[test]
    fn test_kalman1d_std_positive() {
        let kf = KalmanFilter1D::new(0.1, 0.1, 1.0);
        assert!(
            kf.position_std() >= 0.0,
            "position std must be non-negative"
        );
        assert!(
            kf.velocity_std() >= 0.0,
            "velocity std must be non-negative"
        );
    }
    #[test]
    fn test_kalman3d_tracks_constant_position() {
        let mut kf = Kalman3D::new(0.001, 0.001, 0.1);
        let true_pos = [5.0, -3.0, 2.0];
        for _ in 0..20 {
            kf.predict(0.05);
            kf.update(true_pos);
        }
        let est = kf.position();
        for i in 0..3 {
            assert!(
                (est[i] - true_pos[i]).abs() < 1.0,
                "axis {i}: estimated {}, true {}",
                est[i],
                true_pos[i]
            );
        }
    }
    #[test]
    fn test_kalman3d_position_std_non_negative() {
        let kf = Kalman3D::new(0.01, 0.01, 1.0);
        assert!(kf.position_std_rms() >= 0.0);
    }
    #[test]
    fn test_imu_dead_reckoning_straight_line() {
        let mut nav = ImuDeadReckoning::new();
        for _ in 0..10 {
            nav.step(1.0, 0.0, 0.0, 0.1);
        }
        assert!(nav.velocity[0] > 0.8, "forward velocity should have grown");
        assert!(nav.position[0] > 0.0, "should have moved forward");
    }
    #[test]
    fn test_imu_dead_reckoning_yaw_rate() {
        let mut nav = ImuDeadReckoning::new();
        let dt = 0.01;
        for _ in 0..100 {
            nav.step(0.0, 0.0, std::f64::consts::FRAC_PI_2, dt);
        }
        let expected_yaw = std::f64::consts::FRAC_PI_2;
        assert!(
            (nav.heading - expected_yaw).abs() < 0.05,
            "heading should be ~π/2, got {}",
            nav.heading
        );
    }
    #[test]
    fn test_imu_dead_reckoning_reset_position() {
        let mut nav = ImuDeadReckoning::new();
        nav.step(5.0, 0.0, 0.0, 1.0);
        nav.reset_position([100.0, 200.0]);
        assert!((nav.position[0] - 100.0).abs() < 1e-9);
        assert!((nav.position[1] - 200.0).abs() < 1e-9);
    }
    #[test]
    fn test_slip_angle_straight_no_lateral() {
        let est = SlipAngleEstimator::new(1.2, 1.5);
        let (af, ar) = est.estimate(20.0, 0.0, 0.0, 0.0);
        assert!(af.abs() < 1e-9, "front slip should be zero on straight");
        assert!(ar.abs() < 1e-9, "rear slip should be zero on straight");
    }
    #[test]
    fn test_slip_angle_steering_produces_front_slip() {
        let est = SlipAngleEstimator::new(1.2, 1.5);
        let delta = 0.1;
        let (af, _) = est.estimate(15.0, 0.0, 0.0, delta);
        assert!(
            (af - delta).abs() < 0.01,
            "front slip ≈ delta for slow corner"
        );
    }
    #[test]
    fn test_slip_angle_lateral_force() {
        let est = SlipAngleEstimator::new(1.2, 1.5);
        let c_alpha = 30_000.0;
        let fy = est.lateral_force(c_alpha, 0.05);
        assert!(
            (fy - 1500.0).abs() < 1.0,
            "lateral force should be C_alpha * alpha"
        );
    }
    #[test]
    fn test_optical_flow_velocity_roundtrip() {
        let sensor = OpticalFlowSensor::ground_vehicle_camera();
        let vx_true = 3.0;
        let vy_true = 0.5;
        let dt = 0.033;
        let (fx, fy) = sensor.simulate_flow(vx_true, vy_true, dt, 0.0, 0.0);
        let v_est = sensor.velocity_from_flow(fx, fy, dt);
        assert!(
            (v_est[0] - vx_true).abs() < 1e-9,
            "Vx roundtrip: {}",
            v_est[0]
        );
        assert!(
            (v_est[1] - vy_true).abs() < 1e-9,
            "Vy roundtrip: {}",
            v_est[1]
        );
    }
    #[test]
    fn test_optical_flow_zero_velocity_produces_zero_flow() {
        let sensor = OpticalFlowSensor::ground_vehicle_camera();
        let (fx, fy) = sensor.simulate_flow(0.0, 0.0, 0.033, 0.0, 0.0);
        assert!(fx.abs() < 1e-12);
        assert!(fy.abs() < 1e-12);
    }
    #[test]
    fn test_camera_frustum_point_directly_ahead() {
        let cam = CameraFrustum::new(60.0, 45.0, 0.5, 100.0);
        assert!(cam.is_visible([0.0, 0.0, 20.0]));
    }
    #[test]
    fn test_camera_frustum_behind_near_plane() {
        let cam = CameraFrustum::new(60.0, 45.0, 0.5, 100.0);
        assert!(!cam.is_visible([0.0, 0.0, 0.1]));
    }
    #[test]
    fn test_camera_frustum_beyond_far_plane() {
        let cam = CameraFrustum::new(60.0, 45.0, 0.5, 100.0);
        assert!(!cam.is_visible([0.0, 0.0, 150.0]));
    }
    #[test]
    fn test_camera_frustum_side_culled() {
        let cam = CameraFrustum::new(60.0, 45.0, 0.5, 100.0);
        assert!(!cam.is_visible([500.0, 0.0, 10.0]));
    }
    #[test]
    fn test_gps_bias_noise_zero_drift_stays_near_true() {
        let gps = GpsBiasNoise::new(0.5, 0.0);
        let true_pos = [10.0, 20.0, 5.0];
        let meas = gps.measure(true_pos, [0.0, 0.0, 0.0]);
        for i in 0..3 {
            assert!(
                (meas[i] - true_pos[i]).abs() < 1e-9,
                "axis {i}: zero bias+noise should be exact"
            );
        }
    }
    #[test]
    fn test_gps_bias_noise_drift_grows() {
        let mut gps = GpsBiasNoise::new(0.0, 0.01);
        let pos = [0.0, 0.0, 0.0];
        gps.measure(pos, [0.0, 0.0, 0.0]);
        gps.step_bias(1.0);
        gps.step_bias(1.0);
        let any_bias = gps.bias.iter().any(|&b| b.abs() > 1e-12);
        assert!(any_bias, "bias should have grown after steps");
    }
    #[test]
    fn test_allen_deviation_noise_floor_positive() {
        let imu = AllenDeviationImu::new(0.05, 0.001, 0.002, 0.0001);
        let sigma = imu.accel_noise_at_tau(1.0);
        assert!(sigma > 0.0, "noise floor must be positive");
    }
    #[test]
    fn test_allen_deviation_noise_increases_with_bias_instability() {
        let imu_no_bias = AllenDeviationImu::new(0.05, 0.001, 0.0, 0.0);
        let imu_bias = AllenDeviationImu::new(0.05, 0.001, 0.05, 0.001);
        let sigma_no = imu_no_bias.accel_noise_at_tau(100.0);
        let sigma_bi = imu_bias.accel_noise_at_tau(100.0);
        assert!(
            sigma_bi >= sigma_no,
            "more bias instability → higher long-tau noise"
        );
    }
    #[test]
    fn test_wheel_encoder_quantised_speed_exact_multiple() {
        let enc = WheelEncoder::new(100, 0.3);
        let omega = 2.0 * PI * 5.0;
        let speed = enc.quantised_speed(omega, 1.0);
        let expected = 2.0 * PI * 5.0 * 0.3;
        assert!(
            (speed - expected).abs() < 0.1,
            "got {speed}, expected {expected}"
        );
    }
    #[test]
    fn test_wheel_encoder_zero_omega_zero_speed() {
        let enc = WheelEncoder::new(48, 0.31);
        let speed = enc.quantised_speed(0.0, 0.1);
        assert!(speed.abs() < 1e-9);
    }
    #[test]
    fn test_complementary_filter_converges_to_level() {
        let mut cf = ComplementaryFilter::new(0.98);
        let g = 9.81;
        for _ in 0..200 {
            cf.update([0.0, 0.0, g], [0.0, 0.0, 0.0], 0.01);
        }
        assert!(
            cf.roll.abs() < 1e-4,
            "should converge to roll=0, got {}",
            cf.roll
        );
        assert!(
            cf.pitch.abs() < 1e-4,
            "should converge to pitch=0, got {}",
            cf.pitch
        );
    }
    #[test]
    fn test_complementary_filter_roll_from_lateral_g() {
        let mut cf = ComplementaryFilter::new(0.5);
        let g = 9.81;
        for _ in 0..50 {
            cf.update([0.0, g * 0.5, g * 0.866], [0.0, 0.0, 0.0], 0.01);
        }
        assert!(
            cf.roll > 0.0,
            "lateral gravity should produce positive roll"
        );
    }
    #[test]
    fn test_complementary_filter_reset() {
        let mut cf = ComplementaryFilter::new(0.98);
        cf.roll = 1.0;
        cf.pitch = 0.5;
        cf.reset();
        assert_eq!(cf.roll, 0.0);
        assert_eq!(cf.pitch, 0.0);
    }
}
/// Compute the full radar range-rate measurement vector `[range, range_rate]`
/// for a target at 3-D position and velocity.
///
/// # Arguments
/// * `radar_pos`    – radar 3-D position (m)
/// * `radar_vel`    – radar 3-D velocity (m/s)
/// * `target_pos`   – target 3-D position (m)
/// * `target_vel`   – target 3-D velocity (m/s)
/// * `noise_range`  – additive noise for range (m)
/// * `noise_rate`   – additive noise for range-rate (m/s)
///
/// Returns `[range_m, range_rate_m_s]`.
pub fn radar_range_rate_measurement(
    radar_pos: [f64; 3],
    radar_vel: [f64; 3],
    target_pos: [f64; 3],
    target_vel: [f64; 3],
    noise_range: f64,
    noise_rate: f64,
) -> [f64; 2] {
    let dx = target_pos[0] - radar_pos[0];
    let dy = target_pos[1] - radar_pos[1];
    let dz = target_pos[2] - radar_pos[2];
    let range = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9);
    let los = [dx / range, dy / range, dz / range];
    let rel_vx = target_vel[0] - radar_vel[0];
    let rel_vy = target_vel[1] - radar_vel[1];
    let rel_vz = target_vel[2] - radar_vel[2];
    let range_rate = rel_vx * los[0] + rel_vy * los[1] + rel_vz * los[2];
    [range + noise_range, range_rate + noise_rate]
}
#[cfg(test)]
mod tests_sensors_extended {
    use super::super::types_advanced::*;
    use super::*;
    #[test]
    fn test_dop_from_many_sats_is_low() {
        let dop = DopValues::from_satellite_count(12);
        assert!(dop.hdop < 3.0, "many satellites → low HDOP: {}", dop.hdop);
    }
    #[test]
    fn test_dop_from_few_sats_is_high() {
        let dop_good = DopValues::from_satellite_count(12);
        let dop_bad = DopValues::from_satellite_count(4);
        assert!(dop_bad.hdop > dop_good.hdop, "fewer sats → worse DOP");
    }
    #[test]
    fn test_dop_pdop_from_hdop_vdop() {
        let dop = DopValues::from_satellite_count(8);
        let expected_pdop = (dop.hdop * dop.hdop + dop.vdop * dop.vdop).sqrt();
        assert!((dop.pdop - expected_pdop).abs() < 1e-9);
    }
    #[test]
    fn test_gps_extended_zero_noise_exact() {
        let gps = GpsModelExtended::automotive();
        let dop = DopValues::ideal();
        let true_pos = [100.0, 200.0, 10.0];
        let m = gps.measure_position(true_pos, &dop, [0.0, 0.0, 0.0]);
        for i in 0..3 {
            assert!((m[i] - true_pos[i]).abs() < 1e-9, "axis {i}: got {}", m[i]);
        }
    }
    #[test]
    fn test_gps_extended_bias_grows_with_drift() {
        let mut gps = GpsModelExtended::automotive();
        gps.bias_drift_rate = 0.1;
        gps.bias_time_const = 1000.0;
        gps.step_bias(1.0, [1.0, 1.0, 1.0]);
        let any_bias = gps.bias.iter().any(|&b| b.abs() > 1e-12);
        assert!(any_bias, "bias should have grown after step");
    }
    #[test]
    fn test_gps_extended_bias_decays_to_zero() {
        let mut gps = GpsModelExtended::automotive();
        gps.bias = [5.0, -3.0, 1.0];
        gps.bias_time_const = 10.0;
        for _ in 0..200 {
            gps.step_bias(0.5, [0.0, 0.0, 0.0]);
        }
        let mag: f64 = gps.bias.iter().map(|b| b * b).sum::<f64>().sqrt();
        assert!(mag < 0.01, "bias should have decayed, magnitude={mag}");
    }
    #[test]
    fn test_gps_extended_velocity_estimation() {
        let mut gps = GpsModelExtended::automotive();
        let _first = gps.estimate_velocity([0.0, 0.0, 0.0], 0.0);
        assert!(_first.is_none(), "no velocity on first call");
        let vel = gps.estimate_velocity([10.0, 0.0, 0.0], 1.0);
        assert!(vel.is_some());
        let v = vel.unwrap();
        assert!(
            (v[0] - 10.0).abs() < 1e-9,
            "velocity x should be 10 m/s: {}",
            v[0]
        );
    }
    #[test]
    fn test_gps_extended_noise_scales_with_hdop() {
        let gps = GpsModelExtended::new(1.0, 2.0, 0.0);
        let dop_good = DopValues {
            hdop: 1.0,
            vdop: 1.5,
            pdop: 1.8,
        };
        let dop_bad = DopValues {
            hdop: 4.0,
            vdop: 6.0,
            pdop: 7.2,
        };
        let sigma_good = gps.horizontal_1sigma(&dop_good);
        let sigma_bad = gps.horizontal_1sigma(&dop_bad);
        assert!(sigma_bad > sigma_good, "higher HDOP → larger sigma");
    }
    #[test]
    fn test_imu_bias_initial_zero() {
        let bm = ImuBiasModel::automotive_mems();
        assert_eq!(bm.accel_bias, [0.0; 3]);
        assert_eq!(bm.gyro_bias, [0.0; 3]);
    }
    #[test]
    fn test_imu_bias_step_changes_bias() {
        let mut bm = ImuBiasModel::automotive_mems();
        bm.step(0.01, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let accel_mag = bm.accel_bias_rms();
        let gyro_mag = bm.gyro_bias_rms();
        assert!(accel_mag >= 0.0);
        assert!(gyro_mag >= 0.0);
    }
    #[test]
    fn test_imu_bias_decays_with_zero_noise() {
        let mut bm = ImuBiasModel::automotive_mems();
        bm.accel_bias = [0.1, 0.1, 0.1];
        for _ in 0..1000 {
            bm.step(0.1, [0.0; 3], [0.0; 3]);
        }
        let mag = bm.accel_bias_rms();
        assert!(mag < 0.05, "bias should decay with zero noise: {mag}");
    }
    #[test]
    fn test_imu_bias_corrupt_accel_adds_bias() {
        let mut bm = ImuBiasModel::automotive_mems();
        bm.accel_bias = [1.0, 0.0, 0.0];
        let raw = bm.corrupt_accel([0.0, 0.0, 9.81], [0.0; 3], 0.0);
        assert!((raw[0] - 1.0).abs() < 1e-9);
        assert!((raw[2] - 9.81).abs() < 1e-9);
    }
    #[test]
    fn test_imu_bias_reset() {
        let mut bm = ImuBiasModel::automotive_mems();
        bm.accel_bias = [1.0, 2.0, 3.0];
        bm.reset();
        assert_eq!(bm.accel_bias, [0.0; 3]);
        assert_eq!(bm.gyro_bias, [0.0; 3]);
    }
    #[test]
    fn test_radar_extended_cfar_threshold_positive() {
        let radar = RadarExtended::long_range_77ghz();
        let t = radar.cfar.threshold_factor();
        assert!(t > 0.0, "threshold factor must be positive: {t}");
    }
    #[test]
    fn test_radar_extended_quantise_range() {
        let radar = RadarExtended::long_range_77ghz();
        let q = radar.quantise_range(50.15);
        let expected = (50.15_f64 / 0.2).round() * 0.2;
        assert!((q - expected).abs() < 1e-9);
    }
    #[test]
    fn test_radar_extended_quantise_velocity() {
        let radar = RadarExtended::long_range_77ghz();
        let q = radar.quantise_velocity(-5.35);
        let expected = (-5.35_f64 / 0.1).round() * 0.1;
        assert!((q - expected).abs() < 1e-9);
    }
    #[test]
    fn test_radar_extended_detects_close_target() {
        let radar = RadarExtended::long_range_77ghz();
        let det = radar.detect_cfar(0.0, [10.0, 0.0], [0.0, 0.0], [0.0, 0.0], 0.0);
        assert!(det.is_some(), "close target should be detected");
    }
    #[test]
    fn test_radar_extended_unambiguous_range_positive() {
        let radar = RadarExtended::long_range_77ghz();
        let r = radar.unambiguous_range();
        assert!(r > 0.0, "unambiguous range must be positive: {r}");
    }
    #[test]
    fn test_radar_extended_unambiguous_velocity_positive() {
        let radar = RadarExtended::long_range_77ghz();
        let v = radar.unambiguous_velocity();
        assert!(v > 0.0, "unambiguous velocity must be positive: {v}");
    }
    #[test]
    fn test_kalman2d_predict_advances_position() {
        let mut kf = KalmanFusion2D::new(0.01, 0.01, 1.0);
        kf.state[2] = 10.0;
        kf.predict(0.0, 0.0, 1.0);
        assert!(
            (kf.state[0] - 10.0).abs() < 1e-9,
            "x should advance: {}",
            kf.state[0]
        );
    }
    #[test]
    fn test_kalman2d_update_gps_moves_position_estimate() {
        let mut kf = KalmanFusion2D::new(0.1, 0.1, 1.0);
        kf.update_gps([100.0, 0.0]);
        assert!(
            kf.state[0] > 0.0,
            "state should move toward GPS: {}",
            kf.state[0]
        );
    }
    #[test]
    fn test_kalman2d_repeated_gps_converges() {
        let mut kf = KalmanFusion2D::new(0.01, 0.01, 0.5);
        let true_pos = [50.0, -30.0];
        for _ in 0..50 {
            kf.predict(0.0, 0.0, 0.1);
            kf.update_gps(true_pos);
        }
        assert!(
            (kf.state[0] - true_pos[0]).abs() < 5.0,
            "should converge to GPS x"
        );
        assert!(
            (kf.state[1] - true_pos[1]).abs() < 5.0,
            "should converge to GPS y"
        );
    }
    #[test]
    fn test_kalman2d_position_std_non_negative() {
        let kf = KalmanFusion2D::new(0.1, 0.1, 1.0);
        let std = kf.position_std();
        assert!(std[0] >= 0.0 && std[1] >= 0.0);
    }
    #[test]
    fn test_kalman2d_velocity_std_non_negative() {
        let kf = KalmanFusion2D::new(0.1, 0.1, 1.0);
        let std = kf.velocity_std();
        assert!(std[0] >= 0.0 && std[1] >= 0.0);
    }
    #[test]
    fn test_kalman2d_predict_with_accel() {
        let mut kf = KalmanFusion2D::new(0.01, 0.01, 1.0);
        kf.state = [0.0, 0.0, 0.0, 0.0];
        kf.predict(1.0, 0.0, 2.0);
        assert!(kf.state[2] > 0.0, "vx should increase: {}", kf.state[2]);
        assert!(kf.state[0] > 0.0, "x should increase: {}", kf.state[0]);
    }
    #[test]
    fn test_wheel_speed_fusion_all_equal() {
        let wf = WheelSpeedFusion::typical();
        let speed = wf.fuse([10.0, 10.0, 10.0, 10.0]);
        assert!(
            (speed - 10.0 * 0.31).abs() < 1e-9,
            "all-equal should average: {speed}"
        );
    }
    #[test]
    fn test_wheel_speed_fusion_rejects_spinning_wheel() {
        let wf = WheelSpeedFusion::new(0.31, 20.0);
        let omega = [10.0, 10.0, 10.0, 100.0];
        let fused = wf.fuse(omega);
        assert!(
            (fused - 3.1).abs() < 1.0,
            "spinning wheel should be rejected: fused={fused}"
        );
    }
    #[test]
    fn test_wheel_speed_fusion_detect_slip() {
        let wf = WheelSpeedFusion::new(0.31, 3.0);
        let omega = [10.0, 10.0, 10.0, 80.0];
        let slip = wf.detect_slip(omega);
        assert!(slip[3], "last wheel should be detected as slipping");
        assert!(!slip[0], "first three should not be slip");
    }
    #[test]
    fn test_wheel_speed_fusion_zero_speed() {
        let wf = WheelSpeedFusion::typical();
        let speed = wf.fuse([0.0, 0.0, 0.0, 0.0]);
        assert!(speed.abs() < 1e-9);
    }
}
