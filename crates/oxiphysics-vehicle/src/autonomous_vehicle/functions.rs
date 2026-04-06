//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::*;
use rand::RngExt;

/// Synthesise a synthetic LiDAR scan by ray-casting against axis-aligned boxes.
///
/// Each obstacle is an AABB `(min, max)`.  For demonstration purposes the
/// ray is cast along the XY plane only (azimuth sweep, elevation ignored in
/// this simplified version).
pub fn synthesise_lidar_scan(
    config: &LidarConfig,
    sensor_pos: Vec3,
    obstacles: &[(Vec3, Vec3)],
) -> Vec<LidarPoint> {
    let mut rng = rand::rng();
    let mut points = Vec::new();
    let n_az = config.rays_per_revolution;
    let n_el = config.num_channels;
    for el_i in 0..n_el {
        let el_frac = if n_el > 1 {
            el_i as f64 / (n_el - 1) as f64
        } else {
            0.5
        };
        let elevation = -config.vfov_half + 2.0 * config.vfov_half * el_frac;
        let noise_el: f64 = rng.random_range(-config.angle_noise_std..config.angle_noise_std);
        let el_noisy = elevation + noise_el;
        for az_i in 0..n_az {
            let azimuth = 2.0 * std::f64::consts::PI * az_i as f64 / n_az as f64;
            let noise_az: f64 = rng.random_range(-config.angle_noise_std..config.angle_noise_std);
            let az_noisy = azimuth + noise_az;
            let dir = Vec3::new(
                az_noisy.cos() * el_noisy.cos(),
                az_noisy.sin() * el_noisy.cos(),
                el_noisy.sin(),
            );
            let mut hit_range = config.range_max;
            let mut hit = false;
            for &(amin, amax) in obstacles {
                if let Some(t) = ray_aabb_intersect(sensor_pos, dir, amin, amax)
                    && t >= config.range_min
                    && t < hit_range
                {
                    hit_range = t;
                    hit = true;
                }
            }
            if hit {
                let range_noise: f64 =
                    rng.random_range(-config.range_noise_std..config.range_noise_std);
                let r = (hit_range + range_noise).clamp(config.range_min, config.range_max);
                let pos = sensor_pos.add(dir.scale(r));
                let intensity: f64 = rng.random_range(0.3_f64..1.0_f64);
                points.push(LidarPoint {
                    position: pos,
                    intensity,
                    channel: el_i,
                });
            }
        }
    }
    points
}
/// Ray–AABB intersection (slab method).  Returns `Some(t)` if hit, else `None`.
pub fn ray_aabb_intersect(origin: Vec3, dir: Vec3, amin: Vec3, amax: Vec3) -> Option<f64> {
    let inv_dx = if dir.x.abs() > 1e-15 {
        1.0 / dir.x
    } else {
        f64::INFINITY
    };
    let inv_dy = if dir.y.abs() > 1e-15 {
        1.0 / dir.y
    } else {
        f64::INFINITY
    };
    let inv_dz = if dir.z.abs() > 1e-15 {
        1.0 / dir.z
    } else {
        f64::INFINITY
    };
    let t1x = (amin.x - origin.x) * inv_dx;
    let t2x = (amax.x - origin.x) * inv_dx;
    let t1y = (amin.y - origin.y) * inv_dy;
    let t2y = (amax.y - origin.y) * inv_dy;
    let t1z = (amin.z - origin.z) * inv_dz;
    let t2z = (amax.z - origin.z) * inv_dz;
    let tmin = t1x.min(t2x).max(t1y.min(t2y)).max(t1z.min(t2z));
    let tmax = t1x.max(t2x).min(t1y.max(t2y)).min(t1z.max(t2z));
    if tmax >= tmin && tmax >= 0.0 {
        Some(tmin.max(0.0))
    } else {
        None
    }
}
/// Simulate radar detections from a list of moving objects.
pub fn simulate_radar(
    config: &RadarConfig,
    ego_velocity: Vec2,
    objects: &[(Vec2, Vec2, f64)],
) -> Vec<RadarTarget> {
    let mut rng = rand::rng();
    let mut targets = Vec::new();
    for &(pos, vel, rcs) in objects {
        if rcs < config.min_rcs {
            continue;
        }
        let range = pos.norm();
        if range > config.range_max || range < 1e-3 {
            continue;
        }
        let angle = pos.y.atan2(pos.x);
        if angle.abs() > config.hfov_half {
            continue;
        }
        let r_hat = pos.normalise();
        let rel_vel = vel.sub(ego_velocity);
        let v_r = rel_vel.dot(r_hat);
        let rn: f64 = rng.random_range(-config.range_noise_std..config.range_noise_std);
        let an: f64 = rng.random_range(-config.angle_noise_std..config.angle_noise_std);
        let r_noisy = (range + rn).max(0.0);
        let theta_noisy = angle + an;
        let pos_noisy = Vec2::new(r_noisy * theta_noisy.cos(), r_noisy * theta_noisy.sin());
        let snr_db = 20.0 * (rcs / config.min_rcs.max(1e-12)).log10()
            - 40.0 * (range / 1.0_f64.max(1e-12)).log10();
        targets.push(RadarTarget {
            position: pos_noisy,
            radial_velocity: v_r,
            rcs,
            snr_db,
        });
    }
    targets
}
/// Bresenham line algorithm returning integer cell coordinates.
pub fn bresenham_line(x0: isize, y0: isize, x1: isize, y1: isize) -> Vec<(isize, isize)> {
    let mut cells = Vec::new();
    let mut x = x0;
    let mut y = y0;
    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();
    let sx: isize = if x0 < x1 { 1 } else { -1 };
    let sy: isize = if y0 < y1 { 1 } else { -1 };
    let mut err = dx - dy;
    loop {
        cells.push((x, y));
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 > -dy {
            err -= dy;
            x += sx;
        }
        if e2 < dx {
            err += dx;
            y += sy;
        }
    }
    cells
}
/// A* path planner on an occupancy grid.
///
/// Cells with probability > `occ_threshold` are treated as obstacles.
pub fn astar_plan(
    grid: &OccupancyGrid,
    start: (usize, usize),
    goal: (usize, usize),
    occ_threshold: f64,
) -> Option<AStarResult> {
    use std::collections::BinaryHeap;
    use std::collections::HashMap;
    #[derive(PartialEq)]
    pub(super) struct State {
        f: ordered_float::NotNanF64,
        cell: (usize, usize),
    }
    impl Eq for State {}
    impl PartialOrd for State {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }
    impl Ord for State {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            other.f.cmp(&self.f)
        }
    }
    let heuristic = |(c, r): (usize, usize)| {
        let dc = (c as isize - goal.0 as isize).abs() as f64;
        let dr = (r as isize - goal.1 as isize).abs() as f64;
        (dc * dc + dr * dr).sqrt()
    };
    let mut open: BinaryHeap<(i64, (usize, usize))> = BinaryHeap::new();
    let mut g_cost: HashMap<(usize, usize), f64> = HashMap::new();
    let mut came_from: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
    g_cost.insert(start, 0.0);
    let h0 = heuristic(start);
    open.push((-(h0.to_bits() as i64), start));
    pub(super) const DIRS: [(isize, isize); 8] = [
        (1, 0),
        (-1, 0),
        (0, 1),
        (0, -1),
        (1, 1),
        (1, -1),
        (-1, 1),
        (-1, -1),
    ];
    while let Some((_, current)) = open.pop() {
        if current == goal {
            let mut path = vec![goal];
            let mut cur = goal;
            while let Some(&prev) = came_from.get(&cur) {
                path.push(prev);
                cur = prev;
            }
            path.reverse();
            let cost = *g_cost.get(&goal).unwrap_or(&0.0);
            return Some(AStarResult { path, cost });
        }
        let g_cur = *g_cost.get(&current).unwrap_or(&f64::INFINITY);
        for &(dc, dr) in &DIRS {
            let nc = current.0 as isize + dc;
            let nr = current.1 as isize + dr;
            if nc < 0 || nr < 0 {
                continue;
            }
            let (nc, nr) = (nc as usize, nr as usize);
            if nc >= grid.width || nr >= grid.height {
                continue;
            }
            if grid.probability(nc, nr) > occ_threshold {
                continue;
            }
            let step = if dc != 0 && dr != 0 {
                2.0_f64.sqrt()
            } else {
                1.0
            };
            let g_new = g_cur + step;
            let old_g = *g_cost.get(&(nc, nr)).unwrap_or(&f64::INFINITY);
            if g_new < old_g {
                g_cost.insert((nc, nr), g_new);
                came_from.insert((nc, nr), current);
                let f = g_new + heuristic((nc, nr));
                open.push((-(f.to_bits() as i64), (nc, nr)));
            }
        }
    }
    None
}
mod ordered_float {
    /// A newtype for `f64` that implements `Ord` via total ordering.
    #[derive(PartialEq, Clone, Copy, Debug)]
    pub struct NotNanF64(pub f64);
    impl Eq for NotNanF64 {}
    impl PartialOrd for NotNanF64 {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }
    impl Ord for NotNanF64 {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            self.0
                .partial_cmp(&other.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        }
    }
}
pub(super) fn mat5x5_mul(a: &[f64; 25], b: &[f64; 25]) -> [f64; 25] {
    let mut c = [0.0f64; 25];
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                c[i * 5 + j] += a[i * 5 + k] * b[k * 5 + j];
            }
        }
    }
    c
}
pub(super) fn mat5x5_transpose(a: &[f64; 25]) -> [f64; 25] {
    let mut t = [0.0f64; 25];
    for i in 0..5 {
        for j in 0..5 {
            t[j * 5 + i] = a[i * 5 + j];
        }
    }
    t
}
pub(super) fn mat2x5_mul_5x5_mul_5x2(h: &[f64; 10], p: &[f64; 25]) -> [f64; 4] {
    let mut hp = [0.0f64; 10];
    for i in 0..2 {
        for j in 0..5 {
            for k in 0..5 {
                hp[i * 5 + j] += h[i * 5 + k] * p[k * 5 + j];
            }
        }
    }
    let mut res = [0.0f64; 4];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..5 {
                res[i * 2 + j] += hp[i * 5 + k] * h[j * 5 + k];
            }
        }
    }
    res
}
pub(super) fn mat2x5_transpose(h: &[f64; 10]) -> [f64; 10] {
    let mut t = [0.0f64; 10];
    for i in 0..2 {
        for j in 0..5 {
            t[j * 2 + i] = h[i * 5 + j];
        }
    }
    t
}
pub(super) fn mat5x5_mul_5x2(p: &[f64; 25], ht: &[f64; 10]) -> [f64; 10] {
    let mut res = [0.0f64; 10];
    for i in 0..5 {
        for j in 0..2 {
            for k in 0..5 {
                res[i * 2 + j] += p[i * 5 + k] * ht[k * 2 + j];
            }
        }
    }
    res
}
pub(super) fn mat2x2_inv(m: &[f64; 4]) -> [f64; 4] {
    let det = m[0] * m[3] - m[1] * m[2];
    let inv_det = if det.abs() < 1e-30 { 0.0 } else { 1.0 / det };
    [
        m[3] * inv_det,
        -m[1] * inv_det,
        -m[2] * inv_det,
        m[0] * inv_det,
    ]
}
pub(super) fn mat5x2_mul_2x2(k: &[f64; 10], s: &[f64; 4]) -> [f64; 10] {
    let mut res = [0.0f64; 10];
    for i in 0..5 {
        for j in 0..2 {
            for k2 in 0..2 {
                res[i * 2 + j] += k[i * 2 + k2] * s[k2 * 2 + j];
            }
        }
    }
    res
}
pub(super) fn mat5x2_mul_2x5(a: &[f64; 10], b: &[f64; 10]) -> [f64; 25] {
    let mut res = [0.0f64; 25];
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..2 {
                res[i * 5 + j] += a[i * 2 + k] * b[k * 5 + j];
            }
        }
    }
    res
}
/// Solve a 3×3 linear system by Cramer's rule.
pub(super) fn solve3x3(a: &[f64; 9], b: &[f64; 3]) -> Option<[f64; 3]> {
    let det = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6]);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv = 1.0 / det;
    let x0 = inv
        * (b[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (b[1] * a[8] - a[5] * b[2])
            + a[2] * (b[1] * a[7] - a[4] * b[2]));
    let x1 = inv
        * (a[0] * (b[1] * a[8] - a[5] * b[2]) - b[0] * (a[3] * a[8] - a[5] * a[6])
            + a[2] * (a[3] * b[2] - b[1] * a[6]));
    let x2 = inv
        * (a[0] * (a[4] * b[2] - b[1] * a[7]) - a[1] * (a[3] * b[2] - b[1] * a[6])
            + b[0] * (a[3] * a[7] - a[4] * a[6]));
    Some([x0, x1, x2])
}
/// Simplified traffic sign recognition via template matching on 2-D pixel descriptors.
///
/// `descriptors`: each entry is a feature vector of length 8 (HOG-like, normalised).
/// `templates`: (sign_type, feature_vector) pairs.
pub fn recognise_traffic_signs(
    descriptors: &[(Vec2, [f64; 8])],
    templates: &[(TrafficSignType, [f64; 8])],
    ego_pos: Vec2,
    confidence_threshold: f64,
) -> Vec<TrafficSign> {
    let mut results = Vec::new();
    for &(det_pos, ref det_feat) in descriptors {
        let mut best_score = f64::NEG_INFINITY;
        let mut best_type = TrafficSignType::Unknown;
        for (tmpl_type, tmpl_feat) in templates {
            let dot: f64 = det_feat
                .iter()
                .zip(tmpl_feat.iter())
                .map(|(a, b)| a * b)
                .sum();
            let norm_a: f64 = det_feat.iter().map(|a| a * a).sum::<f64>().sqrt();
            let norm_b: f64 = tmpl_feat.iter().map(|b| b * b).sum::<f64>().sqrt();
            let sim = if norm_a > 1e-12 && norm_b > 1e-12 {
                dot / (norm_a * norm_b)
            } else {
                0.0
            };
            if sim > best_score {
                best_score = sim;
                best_type = tmpl_type.clone();
            }
        }
        if best_score >= confidence_threshold {
            let dist = ego_pos.dist(det_pos);
            results.push(TrafficSign {
                position: det_pos,
                sign_type: best_type,
                confidence: best_score,
                distance: dist,
            });
        }
    }
    results
}
#[cfg(test)]
mod tests {
    use super::*;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_vec2_norm() {
        let v = Vec2::new(3.0, 4.0);
        assert!((v.norm() - 5.0).abs() < EPS);
    }
    #[test]
    fn test_vec2_dot() {
        let a = Vec2::new(1.0, 2.0);
        let b = Vec2::new(3.0, 4.0);
        assert!((a.dot(b) - 11.0).abs() < EPS);
    }
    #[test]
    fn test_vec2_cross2() {
        let a = Vec2::new(1.0, 0.0);
        let b = Vec2::new(0.0, 1.0);
        assert!((a.cross2(b) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_ray_aabb_miss() {
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let amin = Vec3::new(5.0, 0.0, 0.0);
        let amax = Vec3::new(6.0, 1.0, 1.0);
        let hit = ray_aabb_intersect(origin, dir, amin, amax);
        assert!(
            hit.is_none(),
            "Ray shooting upward should miss box at x=5..6"
        );
    }
    #[test]
    fn test_ray_aabb_hit() {
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let amin = Vec3::new(5.0, -1.0, -1.0);
        let amax = Vec3::new(6.0, 1.0, 1.0);
        let hit = ray_aabb_intersect(origin, dir, amin, amax);
        assert!(hit.is_some(), "Ray shooting +x should hit box at x=5..6");
        let t = hit.unwrap();
        assert!((t - 5.0).abs() < EPS, "Hit at t=5, got {t}");
    }
    #[test]
    fn test_camera_projection_inside() {
        let cam = CameraIntrinsics {
            fx: 500.0,
            fy: 500.0,
            cx: 320.0,
            cy: 240.0,
            width: 640,
            height: 480,
            near: 0.1,
            far: 100.0,
        };
        let p = Vec3::new(0.0, 0.0, 10.0);
        let proj = cam.project(p);
        assert!(proj.is_some());
        let (u, v) = proj.unwrap();
        assert!((u - 320.0).abs() < EPS && (v - 240.0).abs() < EPS);
    }
    #[test]
    fn test_camera_projection_behind() {
        let cam = CameraIntrinsics {
            fx: 500.0,
            fy: 500.0,
            cx: 320.0,
            cy: 240.0,
            width: 640,
            height: 480,
            near: 0.1,
            far: 100.0,
        };
        let p = Vec3::new(0.0, 0.0, -5.0);
        assert!(cam.project(p).is_none());
    }
    #[test]
    fn test_occupancy_grid_initial() {
        let grid = OccupancyGrid::new(10, 10, 1.0, Vec2::zero());
        assert!((grid.probability(5, 5) - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_occupancy_grid_occupied() {
        let mut grid = OccupancyGrid::new(10, 10, 1.0, Vec2::zero());
        for _ in 0..5 {
            grid.update_occupied(3, 3);
        }
        assert!(grid.probability(3, 3) > 0.5);
    }
    #[test]
    fn test_occupancy_grid_free() {
        let mut grid = OccupancyGrid::new(10, 10, 1.0, Vec2::zero());
        for _ in 0..5 {
            grid.update_free(3, 3);
        }
        assert!(grid.probability(3, 3) < 0.5);
    }
    #[test]
    fn test_bresenham_single_point() {
        let pts = bresenham_line(3, 4, 3, 4);
        assert_eq!(pts.len(), 1);
        assert_eq!(pts[0], (3, 4));
    }
    #[test]
    fn test_bresenham_horizontal() {
        let pts = bresenham_line(0, 0, 4, 0);
        assert_eq!(pts.len(), 5);
        for (i, &(x, y)) in pts.iter().enumerate() {
            assert_eq!(x, i as isize);
            assert_eq!(y, 0);
        }
    }
    #[test]
    fn test_astar_no_obstacles() {
        let grid = OccupancyGrid::new(5, 5, 1.0, Vec2::zero());
        let result = astar_plan(&grid, (0, 0), (4, 4), 0.9);
        assert!(result.is_some(), "A* should find path with no obstacles");
        let path = result.unwrap().path;
        assert_eq!(*path.first().unwrap(), (0, 0));
        assert_eq!(*path.last().unwrap(), (4, 4));
    }
    #[test]
    fn test_astar_blocked_goal() {
        let mut grid = OccupancyGrid::new(5, 5, 1.0, Vec2::zero());
        for i in 3..=4 {
            for j in 3..=4 {
                for _ in 0..20 {
                    grid.update_occupied(i, j);
                }
            }
        }
        let result = astar_plan(&grid, (0, 0), (4, 4), 0.5);
        assert!(result.is_none(), "A* should fail when goal is blocked");
    }
    #[test]
    fn test_rrt_open_space() {
        let start = Vec2::new(0.0, 0.0);
        let goal = Vec2::new(5.0, 5.0);
        let mut planner = RrtPlanner::new(
            start,
            0.5,
            0.6,
            5000,
            Vec2::new(-1.0, -1.0),
            Vec2::new(6.0, 6.0),
        );
        let result = planner.plan(goal, &|_, _| true);
        assert!(result.is_some(), "RRT should find a path in open space");
    }
    #[test]
    fn test_potential_field_start() {
        let pf = PotentialFieldPlanner {
            k_att: 1.0,
            k_rep: 10.0,
            d_safe: 1.0,
            step_size: 0.1,
            max_iter: 50,
        };
        let start = Vec2::new(0.0, 0.0);
        let path = pf.plan(start, Vec2::new(5.0, 5.0), &[]);
        assert!(!path.is_empty());
        assert!((path[0].x - start.x).abs() < EPS && (path[0].y - start.y).abs() < EPS);
    }
    #[test]
    fn test_bicycle_zero_velocity() {
        let s = BicycleState {
            x: 1.0,
            y: 2.0,
            yaw: 0.0,
            v: 0.0,
        };
        let s2 = s.step(0.0, 0.0, 2.8, 0.1);
        assert!((s2.x - 1.0).abs() < EPS && (s2.y - 2.0).abs() < EPS);
    }
    #[test]
    fn test_bicycle_forward() {
        let s = BicycleState {
            x: 0.0,
            y: 0.0,
            yaw: 0.0,
            v: 10.0,
        };
        let s2 = s.step(0.0, 0.0, 2.8, 0.1);
        assert!(
            s2.x > 0.0,
            "Vehicle moving forward should have x > 0: {}",
            s2.x
        );
    }
    #[test]
    fn test_mpc_horizon_length() {
        let mpc = Mpc::default_params();
        let init = BicycleState {
            x: 0.0,
            y: 0.0,
            yaw: 0.0,
            v: 5.0,
        };
        let reference: Vec<BicycleState> = (0..mpc.horizon)
            .map(|k| BicycleState {
                x: k as f64 * 0.5,
                y: 0.0,
                yaw: 0.0,
                v: 5.0,
            })
            .collect();
        let result = mpc.solve(init, &reference, 50);
        assert_eq!(result.controls.len(), mpc.horizon);
        assert_eq!(result.states.len(), mpc.horizon + 1);
    }
    #[test]
    fn test_ekf_predict_covariance_grows() {
        let mut ekf = EkfState::new();
        let p0 = ekf.p[0];
        let q = [0.1f64; 25];
        ekf.predict(0.1, &q);
        assert!(
            ekf.p[0] > p0,
            "Covariance should grow after predict: {} > {}",
            ekf.p[0],
            p0
        );
    }
    #[test]
    fn test_ekf_update_covariance_shrinks() {
        let mut ekf = EkfState::new();
        for i in 0..5 {
            ekf.p[i * 5 + i] = 10.0;
        }
        let p0_diag = ekf.p[0] + ekf.p[6] + ekf.p[12] + ekf.p[18] + ekf.p[24];
        let r = [0.01, 0.0, 0.0, 0.01];
        ekf.update_position([0.0, 0.0], &r);
        let p1_diag = ekf.p[0] + ekf.p[6] + ekf.p[12] + ekf.p[18] + ekf.p[24];
        assert!(
            p1_diag < p0_diag,
            "Trace should shrink after update: {p1_diag} < {p0_diag}"
        );
    }
    #[test]
    fn test_lane_fit_straight() {
        let points: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, 0.5)).collect();
        let lane = LaneModel::fit(&points);
        assert!(lane.is_some(), "Should fit a straight lane");
        let l = lane.unwrap();
        assert!(
            (l.lateral_offset(0.0) - 0.5).abs() < 0.1,
            "Lateral offset at x=0 should be ~0.5"
        );
    }
    #[test]
    fn test_lane_curvature_straight() {
        let points: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, 0.0)).collect();
        let lane = LaneModel::fit(&points).unwrap();
        let curv = lane.curvature(0.0).abs();
        assert!(
            curv < 0.01,
            "Curvature of straight lane should be ~0: {curv}"
        );
    }
    #[test]
    fn test_pedestrian_moves_toward_goal() {
        let desired = Vec2::new(1.0, 0.0);
        let mut ped = Pedestrian::new(Vec2::zero(), desired, 70.0, 0.3);
        for _ in 0..20 {
            ped.step(&[], &[], 0.1, 0.5, 2000.0, 0.08);
        }
        assert!(
            ped.position.x > 0.0,
            "Pedestrian should move in +x: {}",
            ped.position.x
        );
    }
    #[test]
    fn test_traffic_sign_exact_match() {
        let feat = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let descriptors = vec![(Vec2::new(10.0, 0.0), feat)];
        let templates = vec![(TrafficSignType::Stop, feat)];
        let signs = recognise_traffic_signs(&descriptors, &templates, Vec2::zero(), 0.9);
        assert_eq!(signs.len(), 1);
        assert_eq!(signs[0].sign_type, TrafficSignType::Stop);
        assert!((signs[0].confidence - 1.0).abs() < EPS);
    }
    #[test]
    fn test_traffic_sign_below_threshold() {
        let feat_det = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let feat_tmpl = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let descriptors = vec![(Vec2::new(10.0, 0.0), feat_det)];
        let templates = vec![(TrafficSignType::Stop, feat_tmpl)];
        let signs = recognise_traffic_signs(&descriptors, &templates, Vec2::zero(), 0.5);
        assert!(signs.is_empty(), "Orthogonal features should not match");
    }
    #[test]
    fn test_slam_initial_pose() {
        let slam = ParticleFilterSlam::new(100, 0.01, 0.01, 0.1);
        let (x, y, _yaw) = slam.mean_pose();
        assert!(x.abs() < EPS && y.abs() < EPS);
    }
    #[test]
    fn test_slam_predict_forward() {
        let mut slam = ParticleFilterSlam::new(100, 0.001, 0.001, 0.1);
        slam.predict(1.0, 0.0);
        let (x, _y, _yaw) = slam.mean_pose();
        assert!(x > 0.5, "Mean x should be ~1 after forward motion: {x}");
    }
    #[test]
    fn test_slam_resample_count() {
        let mut slam = ParticleFilterSlam::new(50, 0.01, 0.01, 0.1);
        slam.resample();
        assert_eq!(slam.particles.len(), 50);
    }
    #[test]
    fn test_safety_envelope_initial() {
        let state = BicycleState {
            x: 0.0,
            y: 0.0,
            yaw: 0.0,
            v: 10.0,
        };
        let env = SafetyEnvelope::compute(state, 1.0, 0.1, 2.0, 5.0);
        assert!(!env.intervals.is_empty());
        let (_xmin, xmax, ymin, ymax) = env.intervals[0];
        // At v=10 m/s, after dt=0.1s the vehicle travels ~1m forward
        assert!(
            xmax >= 0.5,
            "x_max should include forward motion: xmax={xmax}"
        );
        assert!(
            ymin <= 0.0 && ymax >= 0.0,
            "y range brackets origin: [{ymin},{ymax}]"
        );
    }
    #[test]
    fn test_safety_envelope_outside() {
        let state = BicycleState {
            x: 0.0,
            y: 0.0,
            yaw: 0.0,
            v: 5.0,
        };
        let env = SafetyEnvelope::compute(state, 1.0, 0.1, 0.5, 3.0);
        assert!(!env.contains(9, 1000.0, 0.0));
    }
    #[test]
    fn test_ukf_sigma_mean() {
        let ukf = UkfState::new();
        let pts = ukf.sigma_points();
        assert!((pts[0][0] - ukf.mean[0]).abs() < EPS);
        assert!((pts[0][1] - ukf.mean[1]).abs() < EPS);
        assert!((pts[0][2] - ukf.mean[2]).abs() < EPS);
    }
    #[test]
    fn test_ukf_predict_covariance() {
        let mut ukf = UkfState::new();
        let c0 = ukf.cov[0];
        ukf.predict(0.1, 0.01);
        assert!(
            ukf.cov[0] > c0,
            "Covariance should grow after predict: {} > {}",
            ukf.cov[0],
            c0
        );
    }
    #[test]
    fn test_radar_out_of_fov() {
        let cfg = RadarConfig::automotive_77ghz();
        let ego_vel = Vec2::zero();
        let target_pos = Vec2::new(0.0, 50.0);
        let objects = vec![(target_pos, Vec2::zero(), 10.0)];
        let targets = simulate_radar(&cfg, ego_vel, &objects);
        assert!(
            targets.is_empty(),
            "Target at 90° should be outside 10° FOV"
        );
    }
    #[test]
    fn test_radar_in_fov() {
        let cfg = RadarConfig::automotive_77ghz();
        let ego_vel = Vec2::zero();
        let target_pos = Vec2::new(50.0, 0.0);
        let objects = vec![(target_pos, Vec2::zero(), 5.0)];
        let targets = simulate_radar(&cfg, ego_vel, &objects);
        assert!(
            !targets.is_empty(),
            "Target straight ahead should be detected"
        );
    }
}
