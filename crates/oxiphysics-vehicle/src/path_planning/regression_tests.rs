// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the RRT goal edge, the path smoother and the Dubins
//! planner.

use super::*;
use std::f64::consts::PI;

/// Deterministic xorshift generator so the random poses are reproducible.
struct XorShift(u64);

impl XorShift {
    fn next_unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn pose(&mut self, extent: f64) -> PathPlannerState {
        PathPlannerState {
            x: self.uniform(-extent, extent),
            y: self.uniform(-extent, extent),
            theta: self.uniform(-PI, PI),
            ..Default::default()
        }
    }
}

fn pose(x: f64, y: f64, theta: f64) -> PathPlannerState {
    PathPlannerState {
        x,
        y,
        theta,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// RRT: the edge to the goal must be collision free
// ---------------------------------------------------------------------------

/// 40 x 40 m map with a one-cell-thick vertical wall at x = 30..31 m.  The
/// wall is complete except for a 3 m door at y = 2..5 m when `door` is true.
fn walled_map(door: bool) -> GridMap {
    let mut map = GridMap::new(40, 40, 1.0, [0.0, 0.0]);
    for y in 0..40 {
        if door && (2..5).contains(&y) {
            continue;
        }
        map.set_obstacle(30, y);
    }
    map
}

fn assert_path_is_collision_free(map: &GridMap, path: &[[f64; 2]]) {
    for pair in path.windows(2) {
        assert!(
            map.line_of_sight(pair[0], pair[1]),
            "edge {:?} -> {:?} crosses an obstacle",
            pair[0],
            pair[1]
        );
    }
}

#[test]
fn rrt_does_not_connect_through_a_closed_wall_near_the_goal() {
    let map = walled_map(false);
    let planner = RrtPlanner::new(3000, 2.0, 0.2);
    // The goal sits right behind the wall, within three wall thicknesses of
    // tree nodes on the near side; no collision-free route exists.
    let path = planner.plan(&map, [5.5, 20.5], [31.5, 20.5], 3.0);
    assert!(
        path.is_none(),
        "a path through the closed wall was returned: {path:?}"
    );
}

#[test]
fn rrt_path_around_the_wall_is_collision_free_including_the_goal_edge() {
    let map = walled_map(true);
    let planner = RrtPlanner::new(20_000, 2.0, 0.1);
    let goal = [31.5, 20.5];
    let path = planner
        .plan(&map, [5.5, 20.5], goal, 3.0)
        .expect("a route through the door exists");
    assert_eq!(path.first().copied(), Some([5.5, 20.5]));
    assert_eq!(path.last().copied(), Some(goal));
    assert_path_is_collision_free(&map, &path);
    // The only way through the wall is the door at y = 2..5 m.
    assert!(
        path.iter().any(|p| p[0] > 29.0 && p[1] < 6.0),
        "path does not pass the door: {path:?}"
    );
}

#[test]
fn rrt_open_field_path_still_reaches_the_goal() {
    let map = GridMap::new(40, 40, 1.0, [0.0, 0.0]);
    let planner = RrtPlanner::new(5000, 2.0, 0.2);
    let path = planner
        .plan(&map, [2.5, 2.5], [35.5, 35.5], 2.0)
        .expect("open field");
    assert_path_is_collision_free(&map, &path);
    assert_eq!(path.first().copied(), Some([2.5, 2.5]));
    assert!(dist(*path.last().unwrap(), [35.5, 35.5]) < 1e-9);
}

// ---------------------------------------------------------------------------
// PathSmoother::smooth_path terminates
// ---------------------------------------------------------------------------

fn zigzag() -> Vec<[f64; 2]> {
    (0..10)
        .map(|i| [i as f64, if i % 2 == 0 { 0.0 } else { 1.0 }])
        .collect()
}

fn roughness(path: &[[f64; 2]]) -> f64 {
    path.windows(3)
        .map(|w| {
            let ddx = w[0][0] - 2.0 * w[1][0] + w[2][0];
            let ddy = w[0][1] - 2.0 * w[1][1] + w[2][1];
            ddx * ddx + ddy * ddy
        })
        .sum()
}

#[test]
fn smooth_path_with_zero_tolerance_terminates_and_smooths() {
    let path = zigzag();
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, 0.0);
    assert_eq!(smoothed.len(), path.len());
    assert_eq!(smoothed[0], path[0]);
    assert_eq!(smoothed[9], path[9]);
    assert!(
        smoothed
            .iter()
            .all(|p| p[0].is_finite() && p[1].is_finite())
    );
    assert!(roughness(&smoothed) < 0.5 * roughness(&path));
}

#[test]
fn smooth_path_with_negative_tolerance_terminates() {
    let path = zigzag();
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, -1.0);
    assert_eq!(smoothed.len(), path.len());
    assert!(roughness(&smoothed) < roughness(&path));
}

#[test]
fn smooth_path_with_divergent_weights_returns_the_input() {
    let path = zigzag();
    // weight_smooth = 1.0 makes the in-place update unstable.
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 1.0, 1e-6);
    assert_eq!(smoothed, path);
}

#[test]
fn smooth_path_with_over_relaxed_weights_returns_the_input() {
    // The in-place sweep is an SOR iteration with the factor w_data + 2 w_smooth.
    // From 2 on it oscillates or blows up: 2.02 reaches about 1e86 within the
    // sweep limit without ever overflowing, so only the factor tells the
    // difference between a usable and a ruined path.
    let path = zigzag();
    for (weight_data, weight_smooth) in [
        (0.5, 0.76),
        (0.0, 1.01),
        (0.1, 0.95),
        (1.0, 0.5),
        (2.5, 0.0),
        (0.2, 0.95),
    ] {
        assert_eq!(
            PathSmoother::smooth_path(&path, weight_data, weight_smooth, 1e-6),
            path,
            "weights ({weight_data}, {weight_smooth})"
        );
    }
}

#[test]
fn smooth_path_just_inside_the_stability_bound_converges() {
    let path = zigzag();
    for (weight_data, weight_smooth) in [(0.5, 0.74), (0.0, 0.99), (0.9, 0.5)] {
        let smoothed = PathSmoother::smooth_path(&path, weight_data, weight_smooth, 1e-12);
        assert_eq!(smoothed.len(), path.len());
        for i in 1..path.len() - 1 {
            for k in 0..2 {
                // The minimiser is a weighted mean of the input, so it stays in
                // the bounding box of the input points.
                assert!(
                    (0.0..=9.0).contains(&smoothed[i][0]) && (0.0..=1.0).contains(&smoothed[i][1]),
                    "weights ({weight_data}, {weight_smooth}): point {i} left the input box: {:?}",
                    smoothed[i]
                );
                let residual = weight_data * (path[i][k] - smoothed[i][k])
                    + weight_smooth
                        * (smoothed[i - 1][k] - 2.0 * smoothed[i][k] + smoothed[i + 1][k]);
                assert!(
                    residual.abs() < 1e-8,
                    "weights ({weight_data}, {weight_smooth}): point {i}, axis {k}: {residual}"
                );
            }
        }
    }
}

#[test]
fn smooth_path_with_non_finite_point_returns_in_bounded_time() {
    let mut path = zigzag();
    path[4][1] = f64::NAN;
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, 1e-6);
    assert_eq!(smoothed.len(), path.len());
    assert!(smoothed[4][1].is_nan());
    assert_eq!(smoothed[3], path[3]);
}

#[test]
fn smooth_path_with_invalid_parameters_returns_the_input() {
    let path = zigzag();
    for (wd, ws, tol) in [
        (-0.5, 0.1, 1e-6),
        (0.5, -0.1, 1e-6),
        (f64::NAN, 0.1, 1e-6),
        (0.5, f64::INFINITY, 1e-6),
        (0.5, 0.1, f64::NAN),
        (0.5, 0.1, f64::INFINITY),
    ] {
        assert_eq!(
            PathSmoother::smooth_path(&path, wd, ws, tol),
            path,
            "weights ({wd}, {ws}), tolerance {tol}"
        );
    }
}

#[test]
fn smooth_path_converges_to_a_stationary_point() {
    let path = zigzag();
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, 1e-12);
    for i in 1..path.len() - 1 {
        for k in 0..2 {
            let residual = 0.5 * (path[i][k] - smoothed[i][k])
                + 0.1 * (smoothed[i - 1][k] - 2.0 * smoothed[i][k] + smoothed[i + 1][k]);
            assert!(residual.abs() < 1e-9, "point {i}, axis {k}: {residual}");
        }
    }
}

#[test]
fn smooth_path_stops_after_the_first_sweep_for_a_huge_tolerance() {
    let path = zigzag();
    let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, 1e9);
    // One Gauss-Seidel sweep, computed by hand.
    let mut expected = path.clone();
    for i in 1..path.len() - 1 {
        for k in 0..2 {
            let delta = 0.5 * (path[i][k] - expected[i][k])
                + 0.1 * (expected[i - 1][k] - 2.0 * expected[i][k] + expected[i + 1][k]);
            expected[i][k] += delta;
        }
    }
    assert_eq!(smoothed, expected);
}

// ---------------------------------------------------------------------------
// Dubins planner
// ---------------------------------------------------------------------------

/// End pose of a word, from the closed-form integral of a constant-curvature
/// arc (independent of the incremental sampling in `sample_dubins`).
fn word_end_pose(
    start: PathPlannerState,
    word: DubinsWord,
    segs: [f64; 3],
    radius: f64,
) -> (f64, f64, f64) {
    let turns: [f64; 3] = match word {
        DubinsWord::Lsl => [1.0, 0.0, 1.0],
        DubinsWord::Rsr => [-1.0, 0.0, -1.0],
        DubinsWord::Lsr => [1.0, 0.0, -1.0],
        DubinsWord::Rsl => [-1.0, 0.0, 1.0],
        DubinsWord::Rlr => [-1.0, 1.0, -1.0],
        DubinsWord::Lrl => [1.0, -1.0, 1.0],
    };
    let (mut x, mut y, mut theta) = (start.x, start.y, start.theta);
    for (seg, turn) in segs.iter().zip(turns) {
        let length = seg * radius;
        if turn == 0.0 {
            x += length * theta.cos();
            y += length * theta.sin();
        } else {
            let curvature = turn / radius;
            let theta_new = theta + curvature * length;
            x += (theta_new.sin() - theta.sin()) / curvature;
            y -= (theta_new.cos() - theta.cos()) / curvature;
            theta = theta_new;
        }
    }
    (x, y, theta)
}

fn heading_error(a: f64, b: f64) -> f64 {
    (mod2pi(a - b + PI) - PI).abs()
}

#[test]
fn dubins_words_reach_the_goal_pose_for_random_poses() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let mut words_seen = std::collections::HashSet::new();
    for case in 0..300 {
        let radius = rng.uniform(0.5, 6.0);
        let start = rng.pose(30.0);
        let goal = rng.pose(30.0);
        let (word, segs) =
            shortest_dubins_word(start, goal, radius).expect("a Dubins word always exists");
        let (x, y, theta) = word_end_pose(start, word, segs, radius);
        assert!(
            (x - goal.x).hypot(y - goal.y) < 1e-8,
            "case {case}: {word:?} misses the goal by {} m",
            (x - goal.x).hypot(y - goal.y)
        );
        assert!(
            heading_error(theta, goal.theta) < 1e-8,
            "case {case}: {word:?} heading {theta} vs {}",
            goal.theta
        );
        // The path cannot be shorter than the straight line.
        let length = (segs[0] + segs[1] + segs[2]) * radius;
        assert!(length >= (goal.x - start.x).hypot(goal.y - start.y) - 1e-9);
        words_seen.insert(format!("{word:?}"));
    }
    // Both CSC families are exercised by random poses.
    assert!(words_seen.len() >= 4, "words seen: {words_seen:?}");
}

#[test]
fn every_dubins_word_formula_reaches_the_goal_when_valid() {
    let mut rng = XorShift(0x1234_5678_9abc_def1);
    for case in 0..200 {
        let radius = rng.uniform(0.5, 4.0);
        let start = rng.pose(20.0);
        let goal = rng.pose(20.0);
        let dx = goal.x - start.x;
        let dy = goal.y - start.y;
        let d = dx.hypot(dy) / radius;
        let phi = dy.atan2(dx);
        let alpha = mod2pi(start.theta - phi);
        let beta = mod2pi(goal.theta - phi);
        for (word, segs) in [
            (DubinsWord::Lsl, dubins_lsl(d, alpha, beta)),
            (DubinsWord::Rsr, dubins_rsr(d, alpha, beta)),
            (DubinsWord::Lsr, dubins_lsr(d, alpha, beta)),
            (DubinsWord::Rsl, dubins_rsl(d, alpha, beta)),
            (DubinsWord::Rlr, dubins_rlr(d, alpha, beta)),
            (DubinsWord::Lrl, dubins_lrl(d, alpha, beta)),
        ] {
            let Some(segs) = segs else { continue };
            let (x, y, theta) = word_end_pose(start, word, segs, radius);
            assert!(
                (x - goal.x).hypot(y - goal.y) < 1e-8,
                "case {case}: {word:?} misses the goal"
            );
            assert!(heading_error(theta, goal.theta) < 1e-8);
        }
    }
}

#[test]
fn dubins_plan_samples_end_at_the_goal() {
    let mut rng = XorShift(0x0bad_cafe_dead_beef);
    for case in 0..20 {
        let radius = rng.uniform(1.0, 5.0);
        let start = rng.pose(25.0);
        let goal = rng.pose(25.0);
        let path = DubinsPath::plan(start, goal, radius);
        let first = path[0];
        let last = *path.last().unwrap();
        assert!(dist(first, [start.x, start.y]) < 1e-12);
        assert!(
            dist(last, [goal.x, goal.y]) < 0.06 * radius,
            "case {case}: ends {} m from the goal",
            dist(last, [goal.x, goal.y])
        );
        // Waypoints are at most one sampling step apart.
        for pair in path.windows(2) {
            assert!(dist(pair[0], pair[1]) <= DUBINS_SAMPLE_SPACING + 1e-9);
        }
        // Final heading, from the last chord: the chord of a 0.05 m step
        // deviates from the tangent by less than 0.05 / (2 r) rad.
        let prev = path[path.len() - 2];
        let chord_heading = (last[1] - prev[1]).atan2(last[0] - prev[0]);
        assert!(
            heading_error(chord_heading, goal.theta) < 1e-2 + DUBINS_SAMPLE_SPACING / radius,
            "case {case}: chord heading {chord_heading} vs goal {}",
            goal.theta
        );
    }
}

#[test]
fn dubins_word_for_aligned_poses_is_the_straight_line() {
    // alpha = beta = 0 and d = 10: LSL degenerates to a straight run of 10 radii.
    let segs = dubins_lsl(10.0, 0.0, 0.0).expect("valid");
    assert!(segs[0].abs() < 1e-12);
    assert!((segs[1] - 10.0).abs() < 1e-12);
    assert!(segs[2].abs() < 1e-12);

    let path = DubinsPath::plan(pose(0.0, 0.0, 0.0), pose(10.0, 0.0, 0.0), 1.0);
    assert!(path.iter().all(|p| p[1].abs() < 1e-9));
    assert!((PathSmoother::path_length(&path) - 10.0).abs() < 1e-9);
}

#[test]
fn dubins_frame_follows_the_start_goal_chord() {
    // The chord points along +y: heading +y at both ends is a straight line.
    let up = DubinsPath::plan(pose(0.0, 0.0, PI / 2.0), pose(0.0, 10.0, PI / 2.0), 1.0);
    assert!((PathSmoother::path_length(&up) - 10.0).abs() < 1e-9);
    assert!(up.iter().all(|p| p[0].abs() < 1e-9));

    // The same along a 45 degree chord.
    let diag = DubinsPath::plan(
        pose(1.0, 2.0, PI / 4.0),
        pose(1.0 + 7.0, 2.0 + 7.0, PI / 4.0),
        2.0,
    );
    let expected = 7.0 * 2.0_f64.sqrt();
    assert!((PathSmoother::path_length(&diag) - expected).abs() < 1e-9);
}

#[test]
fn dubins_u_turn_has_the_textbook_length() {
    // Start at the origin heading +x, goal 2 r to the left heading -x:
    // a half circle of radius r (length pi r) joins them.
    let radius = 3.0;
    let path = DubinsPath::plan(pose(0.0, 0.0, 0.0), pose(0.0, 2.0 * radius, PI), radius);
    let length = PathSmoother::path_length(&path);
    assert!((length - PI * radius).abs() < 1e-2, "length {length}");
}

#[test]
fn dubins_invalid_radius_falls_back_to_the_chord() {
    let start = pose(1.0, 2.0, 0.3);
    let goal = pose(4.0, 6.0, -0.2);
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let path = DubinsPath::plan(start, goal, radius);
        assert_eq!(path, vec![[1.0, 2.0], [4.0, 6.0]], "radius {radius}");
    }
}

#[test]
fn mod2pi_wraps_into_one_turn() {
    assert!((mod2pi(-0.5) - (2.0 * PI - 0.5)).abs() < 1e-12);
    assert!((mod2pi(7.0) - (7.0 - 2.0 * PI)).abs() < 1e-12);
    assert_eq!(mod2pi(0.0), 0.0);
    // A tiny negative angle must not become a full spurious turn.
    assert_eq!(mod2pi(-1e-17), 0.0);
    for x in [-20.0, -6.3, -1e-9, 0.0, 3.0, 6.2, 25.0] {
        let v = mod2pi(x);
        assert!((0.0..2.0 * PI).contains(&v), "mod2pi({x}) = {v}");
    }
}

// ---------------------------------------------------------------------------
// Nearest-neighbour index of the RRT
// ---------------------------------------------------------------------------

#[test]
fn node_grid_agrees_with_a_linear_scan() {
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    for (count, cell, extent) in [
        (1_usize, 2.0, 40.0),
        (5, 2.0, 40.0),
        (200, 2.0, 40.0),
        (2000, 2.0, 40.0),
        (2000, 0.05, 40.0),
        (400, 3.0, 1.0e5),
    ] {
        let mut grid = NodeGrid::new(cell);
        let mut nodes: Vec<[f64; 2]> = Vec::new();
        for k in 0..count {
            let p = [
                rng.uniform(-0.2 * extent, 0.8 * extent),
                rng.uniform(-0.5 * extent, 0.5 * extent),
            ];
            nodes.push(p);
            grid.insert(k, p);
            if k.is_multiple_of(5) {
                for _ in 0..3 {
                    let q = [
                        rng.uniform(-0.75 * extent, 0.75 * extent),
                        rng.uniform(-0.75 * extent, 0.75 * extent),
                    ];
                    assert_eq!(
                        grid.nearest(&nodes, q),
                        nearest_node_linear(&nodes, q),
                        "{k} nodes, cell {cell}, query {q:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn node_grid_breaks_ties_towards_the_lowest_index() {
    // identical nodes
    let mut grid = NodeGrid::new(2.0);
    let mut nodes = Vec::new();
    for k in 0..10 {
        nodes.push([1.0, 1.0]);
        grid.insert(k, [1.0, 1.0]);
    }
    assert_eq!(grid.nearest(&nodes, [1.5, 1.0]), Some(0));

    // three nodes in three different cells at exactly the same distance
    let nodes = vec![[3.0, 0.0], [-3.0, 0.0], [0.0, 3.0]];
    let mut grid = NodeGrid::new(2.0);
    for (k, p) in nodes.iter().enumerate() {
        grid.insert(k, *p);
    }
    assert_eq!(grid.nearest(&nodes, [0.0, 0.0]), Some(0));
}

#[test]
fn node_grid_survives_degenerate_cell_sizes_and_coordinates() {
    let mut rng = XorShift(0x0123_4567_89ab_cdef);
    let nodes: Vec<[f64; 2]> = (0..300)
        .map(|_| [rng.uniform(0.0, 40.0), rng.uniform(0.0, 40.0)])
        .collect();
    for cell in [f64::NAN, -1.0, 0.0, f64::INFINITY, 1.0e-12] {
        let mut grid = NodeGrid::new(cell);
        for (k, p) in nodes.iter().enumerate() {
            grid.insert(k, *p);
        }
        for q in [[20.0, 20.0], [-5.0, 50.0], [1.0e9, -1.0e9]] {
            assert_eq!(
                grid.nearest(&nodes, q),
                nearest_node_linear(&nodes, q),
                "cell size {cell}, query {q:?}"
            );
        }
    }

    let mut grid = NodeGrid::new(2.0);
    for (k, p) in nodes.iter().enumerate() {
        grid.insert(k, *p);
    }
    for q in [
        [f64::INFINITY, 0.0],
        [f64::NEG_INFINITY, f64::INFINITY],
        [f64::NAN, 1.0],
        [1.0e300, -1.0e300],
        [2.0e16, 3.0],
        [-1.0e14, 5.0],
    ] {
        assert_eq!(
            grid.nearest(&nodes, q),
            nearest_node_linear(&nodes, q),
            "query {q:?}"
        );
    }

    // a node that cannot be bucketed is still found
    let mut far = nodes.clone();
    far.push([1.0e300, 0.0]);
    let mut grid = NodeGrid::new(2.0);
    for (k, p) in far.iter().enumerate() {
        grid.insert(k, *p);
    }
    assert_eq!(grid.nearest(&far, [1.0e300, 1.0]), Some(far.len() - 1));
    assert_eq!(
        grid.nearest(&far, [20.0, 20.0]),
        nearest_node_linear(&far, [20.0, 20.0])
    );
}

#[test]
fn node_grid_of_nothing_has_no_nearest_node() {
    let grid = NodeGrid::new(2.0);
    assert_eq!(grid.nearest(&[], [0.0, 0.0]), None);
    assert_eq!(nearest_node_linear(&[], [0.0, 0.0]), None);
}

#[test]
fn rrt_with_a_degenerate_step_size_terminates() {
    let map = GridMap::new(20, 20, 1.0, [0.0, 0.0]);
    let (start, goal) = ([1.5, 1.5], [18.5, 18.5]);
    for step in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let planner = RrtPlanner::new(200, step, 0.1);
        // Whatever such a step does to the tree, planning has to return and any
        // path it reports has to connect the two ends.
        if let Some(path) = planner.plan(&map, start, goal, 1.0) {
            assert_eq!(path.first().copied(), Some(start), "step {step}");
            assert_eq!(path.last().copied(), Some(goal), "step {step}");
        }
    }
}
