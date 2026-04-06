// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Autonomous vehicle path planning.
//!
//! Implements A\*, RRT, and potential field planners on a discrete `GridMap`,
//! together with Dubins-path generation and path smoothing utilities.

#![allow(dead_code)]

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

// ---------------------------------------------------------------------------
// GridMap
// ---------------------------------------------------------------------------

/// A 2-D occupancy / cost grid used by all planners.
#[derive(Debug, Clone)]
pub struct GridMap {
    /// Number of columns.
    pub width: usize,
    /// Number of rows.
    pub height: usize,
    /// Metres per cell.
    pub resolution: f64,
    /// World-space coordinates of the (0, 0) cell corner.
    pub origin: [f64; 2],
    /// Flat row-major cell values: 0 = free, 255 = obstacle, 1-254 = cost.
    pub cells: Vec<u8>,
}

impl GridMap {
    /// Create a fully-free grid.
    pub fn new(width: usize, height: usize, resolution: f64, origin: [f64; 2]) -> Self {
        Self {
            width,
            height,
            resolution,
            origin,
            cells: vec![0u8; width * height],
        }
    }

    /// Row-major flat index.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    /// Convert world coordinates to grid cell indices (returns `None` if out of bounds).
    pub fn world_to_grid(&self, wx: f64, wy: f64) -> Option<(usize, usize)> {
        let fx = (wx - self.origin[0]) / self.resolution;
        let fy = (wy - self.origin[1]) / self.resolution;
        if fx < 0.0 || fy < 0.0 {
            return None;
        }
        let gx = fx as usize;
        let gy = fy as usize;
        if gx < self.width && gy < self.height {
            Some((gx, gy))
        } else {
            None
        }
    }

    /// Convert grid cell to the world-space centre of that cell.
    pub fn grid_to_world(&self, gx: usize, gy: usize) -> [f64; 2] {
        [
            self.origin[0] + (gx as f64 + 0.5) * self.resolution,
            self.origin[1] + (gy as f64 + 0.5) * self.resolution,
        ]
    }

    /// Returns `true` if the cell exists and is not an obstacle.
    pub fn is_free(&self, x: usize, y: usize) -> bool {
        if x >= self.width || y >= self.height {
            return false;
        }
        self.cells[self.idx(x, y)] < 255
    }

    /// Mark a single cell as obstacle.
    pub fn set_obstacle(&mut self, x: usize, y: usize) {
        if x < self.width && y < self.height {
            let i = self.idx(x, y);
            self.cells[i] = 255;
        }
    }

    /// Mark all cells within `radius` metres of world point `(cx, cy)` as obstacles.
    pub fn set_circle_obstacle(&mut self, cx: f64, cy: f64, radius: f64) {
        let r_cells = (radius / self.resolution).ceil() as isize;
        if let Some((ocx, ocy)) = self.world_to_grid(cx, cy) {
            let ocx = ocx as isize;
            let ocy = ocy as isize;
            for dy in -r_cells..=r_cells {
                for dx in -r_cells..=r_cells {
                    let dist = ((dx * dx + dy * dy) as f64).sqrt() * self.resolution;
                    if dist <= radius {
                        let nx = ocx + dx;
                        let ny = ocy + dy;
                        if nx >= 0 && ny >= 0 {
                            self.set_obstacle(nx as usize, ny as usize);
                        }
                    }
                }
            }
        }
    }

    /// Inflate every obstacle by `radius` metres (Minkowski sum with a disc).
    pub fn inflate_obstacles(&mut self, radius: f64) {
        let r_cells = (radius / self.resolution).ceil() as isize;
        let original = self.cells.clone();
        let w = self.width as isize;
        let h = self.height as isize;
        for y in 0..h {
            for x in 0..w {
                if original[(y * w + x) as usize] == 255 {
                    for dy in -r_cells..=r_cells {
                        for dx in -r_cells..=r_cells {
                            let dist = ((dx * dx + dy * dy) as f64).sqrt() * self.resolution;
                            if dist <= radius {
                                let nx = x + dx;
                                let ny = y + dy;
                                if nx >= 0 && ny >= 0 && nx < w && ny < h {
                                    let i = (ny * w + nx) as usize;
                                    if original[i] < 255 {
                                        self.cells[i] = 255;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Bresenham line-of-sight check: returns `true` if the straight line
    /// between two world points passes only through free cells.
    fn line_of_sight(&self, p0: [f64; 2], p1: [f64; 2]) -> bool {
        let (x0, y0) = match self.world_to_grid(p0[0], p0[1]) {
            Some(v) => (v.0 as isize, v.1 as isize),
            None => return false,
        };
        let (x1, y1) = match self.world_to_grid(p1[0], p1[1]) {
            Some(v) => (v.0 as isize, v.1 as isize),
            None => return false,
        };
        let mut cx = x0;
        let mut cy = y0;
        let dx = (x1 - x0).abs();
        let dy = (y1 - y0).abs();
        let sx: isize = if x0 < x1 { 1 } else { -1 };
        let sy: isize = if y0 < y1 { 1 } else { -1 };
        let mut err = dx - dy;
        loop {
            if cx < 0 || cy < 0 || cx >= self.width as isize || cy >= self.height as isize {
                return false;
            }
            if !self.is_free(cx as usize, cy as usize) {
                return false;
            }
            if cx == x1 && cy == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 > -dy {
                err -= dy;
                cx += sx;
            }
            if e2 < dx {
                err += dx;
                cy += sy;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// A* planner
// ---------------------------------------------------------------------------

#[derive(Copy, Clone)]
struct AStarNode {
    f: f64,
    g: f64,
    pos: (usize, usize),
}

impl PartialEq for AStarNode {
    fn eq(&self, other: &Self) -> bool {
        self.f == other.f
    }
}
impl Eq for AStarNode {}
impl Ord for AStarNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other.f.partial_cmp(&self.f).unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for AStarNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A\* grid planner.
pub struct AStarPlanner;

impl AStarPlanner {
    /// Plan a path from `start` to `goal` (world coordinates) on `map`.
    ///
    /// Uses an 8-connected grid with Euclidean-distance heuristic.
    /// Returns world-coordinate waypoints after shortcut smoothing, or `None`
    /// if no path exists.
    pub fn plan(map: &GridMap, start: [f64; 2], goal: [f64; 2]) -> Option<Vec<[f64; 2]>> {
        let (sx, sy) = map.world_to_grid(start[0], start[1])?;
        let (gx, gy) = map.world_to_grid(goal[0], goal[1])?;

        if !map.is_free(sx, sy) || !map.is_free(gx, gy) {
            return None;
        }

        let heuristic = |x: usize, y: usize| -> f64 {
            let dx = x as f64 - gx as f64;
            let dy = y as f64 - gy as f64;
            (dx * dx + dy * dy).sqrt()
        };

        let mut open: BinaryHeap<AStarNode> = BinaryHeap::new();
        let mut g_score: HashMap<(usize, usize), f64> = HashMap::new();
        let mut came_from: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

        g_score.insert((sx, sy), 0.0);
        open.push(AStarNode {
            f: heuristic(sx, sy),
            g: 0.0,
            pos: (sx, sy),
        });

        const DIRS: [(isize, isize); 8] = [
            (1, 0),
            (-1, 0),
            (0, 1),
            (0, -1),
            (1, 1),
            (1, -1),
            (-1, 1),
            (-1, -1),
        ];

        while let Some(current) = open.pop() {
            let (cx, cy) = current.pos;
            if cx == gx && cy == gy {
                // Reconstruct grid path
                let mut grid_path = vec![(cx, cy)];
                let mut p = (cx, cy);
                while let Some(&prev) = came_from.get(&p) {
                    grid_path.push(prev);
                    p = prev;
                }
                grid_path.reverse();

                // Convert to world coords
                let world: Vec<[f64; 2]> = grid_path
                    .iter()
                    .map(|&(x, y)| map.grid_to_world(x, y))
                    .collect();

                // Shortcut smooth
                return Some(PathSmoother::shortcut_path(&world, map));
            }

            let g_curr = *g_score.get(&(cx, cy)).unwrap_or(&f64::INFINITY);
            if current.g > g_curr + 1e-9 {
                continue; // stale entry
            }

            for (ddx, ddy) in &DIRS {
                let nx = cx as isize + ddx;
                let ny = cy as isize + ddy;
                if nx < 0 || ny < 0 || nx >= map.width as isize || ny >= map.height as isize {
                    continue;
                }
                let (nx, ny) = (nx as usize, ny as usize);
                if !map.is_free(nx, ny) {
                    continue;
                }
                let step = if ddx.abs() + ddy.abs() == 2 {
                    std::f64::consts::SQRT_2
                } else {
                    1.0
                };
                let tentative_g = g_curr + step;
                let prev_g = *g_score.get(&(nx, ny)).unwrap_or(&f64::INFINITY);
                if tentative_g < prev_g {
                    g_score.insert((nx, ny), tentative_g);
                    came_from.insert((nx, ny), (cx, cy));
                    open.push(AStarNode {
                        f: tentative_g + heuristic(nx, ny),
                        g: tentative_g,
                        pos: (nx, ny),
                    });
                }
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// RRT planner
// ---------------------------------------------------------------------------

/// Basic Rapidly-exploring Random Tree planner.
pub struct RrtPlanner {
    /// Maximum number of tree-extension iterations.
    pub max_iterations: usize,
    /// Maximum branch length per extension step (world units).
    pub step_size: f64,
    /// Probability (0..1) of sampling the goal instead of a random point.
    pub goal_bias: f64,
}

impl RrtPlanner {
    /// Create a planner with sensible defaults.
    pub fn new(max_iterations: usize, step_size: f64, goal_bias: f64) -> Self {
        Self {
            max_iterations,
            step_size,
            goal_bias,
        }
    }

    /// Plan a path on `map` from `start` to `goal`.
    ///
    /// Returns `None` when the tree exhausts `max_iterations` without reaching
    /// within `goal_radius` of `goal`.
    pub fn plan(
        &self,
        map: &GridMap,
        start: [f64; 2],
        goal: [f64; 2],
        goal_radius: f64,
    ) -> Option<Vec<[f64; 2]>> {
        // Simple LCG-based pseudo-random numbers (no external deps).
        let mut seed: u64 = 0xdeadbeef_cafebabe;
        let mut rng = || -> f64 {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed as f64) / (u64::MAX as f64)
        };

        let wx_min = map.origin[0];
        let wy_min = map.origin[1];
        let wx_max = wx_min + map.width as f64 * map.resolution;
        let wy_max = wy_min + map.height as f64 * map.resolution;

        let mut nodes: Vec<[f64; 2]> = vec![start];
        let mut parent: Vec<usize> = vec![0]; // parent[0] is unused (root)

        for _ in 0..self.max_iterations {
            // Sample
            let sample = if rng() < self.goal_bias {
                goal
            } else {
                [
                    wx_min + rng() * (wx_max - wx_min),
                    wy_min + rng() * (wy_max - wy_min),
                ]
            };

            // Nearest node
            let nearest_idx = nodes
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    dist2(**a, sample)
                        .partial_cmp(&dist2(**b, sample))
                        .unwrap_or(Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);

            let nearest = nodes[nearest_idx];
            let d = dist(nearest, sample);
            if d < 1e-9 {
                continue;
            }
            let t = (self.step_size / d).min(1.0);
            let new_node = [
                nearest[0] + t * (sample[0] - nearest[0]),
                nearest[1] + t * (sample[1] - nearest[1]),
            ];

            // Collision check
            if !map.line_of_sight(nearest, new_node) {
                continue;
            }

            nodes.push(new_node);
            parent.push(nearest_idx);

            // Goal check
            if dist(new_node, goal) <= goal_radius {
                // Add goal node if not already there
                let goal_idx = if dist(new_node, goal) > 1e-9 {
                    nodes.push(goal);
                    parent.push(nodes.len() - 2);
                    nodes.len() - 1
                } else {
                    nodes.len() - 1
                };

                // Trace path back
                let mut path = Vec::new();
                let mut idx = goal_idx;
                path.push(nodes[idx]);
                while idx != 0 {
                    idx = parent[idx];
                    path.push(nodes[idx]);
                }
                path.reverse();
                return Some(path);
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Potential field planner
// ---------------------------------------------------------------------------

/// Artificial potential field planner.
pub struct PotentialField {
    /// Gain for the attractive potential towards the goal.
    pub attractive_gain: f64,
    /// Gain for the repulsive potential away from obstacles.
    pub repulsive_gain: f64,
    /// Distance threshold (metres) beyond which obstacles exert no repulsion.
    pub repulsive_influence: f64,
}

impl PotentialField {
    /// Quadratic attractive potential: ½ k_att · d²(pos, goal).
    pub fn attractive_potential(&self, pos: [f64; 2], goal: [f64; 2]) -> f64 {
        0.5 * self.attractive_gain * dist2(pos, goal)
    }

    /// Repulsive potential summed over all obstacle cells within the influence
    /// radius.
    pub fn repulsive_potential(&self, pos: [f64; 2], map: &GridMap) -> f64 {
        let r = self.repulsive_influence;
        let r_cells = (r / map.resolution).ceil() as isize;

        let (cx, cy) = match map.world_to_grid(pos[0], pos[1]) {
            Some(v) => (v.0 as isize, v.1 as isize),
            None => return 0.0,
        };

        let mut u_rep = 0.0;
        for dy in -r_cells..=r_cells {
            for dx in -r_cells..=r_cells {
                let nx = cx + dx;
                let ny = cy + dy;
                if nx < 0 || ny < 0 || nx >= map.width as isize || ny >= map.height as isize {
                    continue;
                }
                if map.cells[map.idx(nx as usize, ny as usize)] < 255 {
                    continue;
                }
                let obs_world = map.grid_to_world(nx as usize, ny as usize);
                let d = dist(pos, obs_world).max(1e-6);
                if d < r {
                    let term = 1.0 / d - 1.0 / r;
                    u_rep += 0.5 * self.repulsive_gain * term * term;
                }
            }
        }
        u_rep
    }

    /// Total force = −∇(U_att + U_rep), computed by finite differences.
    pub fn total_force(&self, pos: [f64; 2], goal: [f64; 2], map: &GridMap) -> [f64; 2] {
        let eps = map.resolution * 0.1;
        let u0 = self.attractive_potential(pos, goal) + self.repulsive_potential(pos, map);

        let px = [pos[0] + eps, pos[1]];
        let py = [pos[0], pos[1] + eps];
        let ux = self.attractive_potential(px, goal) + self.repulsive_potential(px, map);
        let uy = self.attractive_potential(py, goal) + self.repulsive_potential(py, map);

        [-(ux - u0) / eps, -(uy - u0) / eps]
    }

    /// Follow the negative gradient for up to `max_steps` steps of size `dt`.
    ///
    /// The path always starts at `start` and terminates when within one cell of
    /// the goal or when `max_steps` is reached.
    pub fn plan(
        &self,
        map: &GridMap,
        start: [f64; 2],
        goal: [f64; 2],
        max_steps: usize,
        dt: f64,
    ) -> Vec<[f64; 2]> {
        let mut path = vec![start];
        let mut pos = start;
        let stop_dist = map.resolution;

        for _ in 0..max_steps {
            if dist(pos, goal) < stop_dist {
                path.push(goal);
                break;
            }
            let force = self.total_force(pos, goal, map);
            let mag = (force[0] * force[0] + force[1] * force[1]).sqrt();
            if mag < 1e-12 {
                break; // local minimum
            }
            pos = [pos[0] + dt * force[0], pos[1] + dt * force[1]];
            path.push(pos);
        }
        path
    }
}

// ---------------------------------------------------------------------------
// Path smoother
// ---------------------------------------------------------------------------

/// Collection of path post-processing utilities.
pub struct PathSmoother;

impl PathSmoother {
    /// Gradient-descent path smoothing.
    ///
    /// Minimises `w_data · Σ|pᵢ − pᵢ₀|² + w_smooth · Σ|pᵢ₊₁ − 2pᵢ + pᵢ₋₁|²`
    /// while keeping the first and last points fixed.
    pub fn smooth_path(
        path: &[[f64; 2]],
        weight_data: f64,
        weight_smooth: f64,
        tolerance: f64,
    ) -> Vec<[f64; 2]> {
        if path.len() < 3 {
            return path.to_vec();
        }
        let n = path.len();
        let orig = path.to_vec();
        let mut smoothed = path.to_vec();

        loop {
            let mut change = 0.0_f64;
            for i in 1..n - 1 {
                for k in 0..2 {
                    let data_term = weight_data * (orig[i][k] - smoothed[i][k]);
                    let smooth_term = weight_smooth
                        * (smoothed[i - 1][k] - 2.0 * smoothed[i][k] + smoothed[i + 1][k]);
                    let delta = data_term + smooth_term;
                    smoothed[i][k] += delta;
                    change += delta.abs();
                }
            }
            if change < tolerance {
                break;
            }
        }
        smoothed
    }

    /// Remove waypoints whose removal keeps the path collision-free.
    pub fn shortcut_path(path: &[[f64; 2]], map: &GridMap) -> Vec<[f64; 2]> {
        if path.len() <= 2 {
            return path.to_vec();
        }
        let mut result = vec![path[0]];
        let mut i = 0;
        while i < path.len() - 1 {
            let mut j = path.len() - 1;
            while j > i + 1 {
                if map.line_of_sight(path[i], path[j]) {
                    break;
                }
                j -= 1;
            }
            result.push(path[j]);
            i = j;
        }
        result
    }

    /// Total arc-length of the path.
    pub fn path_length(path: &[[f64; 2]]) -> f64 {
        path.windows(2).map(|w| dist(w[0], w[1])).sum()
    }

    /// Resample the path so consecutive waypoints are spaced `spacing` apart.
    pub fn resample_path(path: &[[f64; 2]], spacing: f64) -> Vec<[f64; 2]> {
        if path.len() < 2 || spacing <= 0.0 {
            return path.to_vec();
        }
        let mut result = vec![path[0]];
        let mut accumulated = 0.0_f64;
        let mut prev = path[0];

        for &pt in &path[1..] {
            let seg_len = dist(prev, pt);
            if seg_len < 1e-12 {
                prev = pt;
                continue;
            }
            let mut remaining = seg_len;
            let dir = [(pt[0] - prev[0]) / seg_len, (pt[1] - prev[1]) / seg_len];
            let mut cursor = prev;

            while accumulated + remaining >= spacing {
                let advance = spacing - accumulated;
                cursor = [cursor[0] + advance * dir[0], cursor[1] + advance * dir[1]];
                result.push(cursor);
                remaining -= advance;
                accumulated = 0.0;
            }
            accumulated += remaining;
            prev = pt;
        }
        // Always include the last point unless it coincides with the last resampled point
        if dist(
            *result.last().expect("collection should not be empty"),
            *path.last().expect("collection should not be empty"),
        ) > 1e-9
        {
            result.push(*path.last().expect("collection should not be empty"));
        }
        result
    }
}

// ---------------------------------------------------------------------------
// Vehicle state for path planning
// ---------------------------------------------------------------------------

/// 2-D kinematic vehicle state used by the Dubins path planner.
#[derive(Debug, Clone, Copy, Default)]
pub struct PathPlannerState {
    /// X position (metres).
    pub x: f64,
    /// Y position (metres).
    pub y: f64,
    /// Heading angle (radians, measured from +X axis).
    pub theta: f64,
    /// Forward speed (m/s).
    pub v: f64,
    /// Angular velocity (rad/s).
    pub omega: f64,
}

// ---------------------------------------------------------------------------
// Dubins path
// ---------------------------------------------------------------------------

/// Shortest Dubins path between two oriented poses.
pub struct DubinsPath;

/// The six canonical Dubins path word types.
#[derive(Debug, Clone, Copy, PartialEq)]
enum DubinsWord {
    Lsl,
    Rsr,
    Lsr,
    Rsl,
    Rlr,
    Lrl,
}

impl DubinsPath {
    /// Compute the shortest Dubins path from `start` to `goal` with minimum
    /// turning radius `min_radius`.
    ///
    /// Returns uniformly sampled world-space waypoints at ~0.1 m arc-length
    /// resolution.
    pub fn plan(start: PathPlannerState, goal: PathPlannerState, min_radius: f64) -> Vec<[f64; 2]> {
        use std::f64::consts::PI;

        // Normalise to the standard Dubins frame (start at origin, heading 0)
        let dx = goal.x - start.x;
        let dy = goal.y - start.y;
        let d = (dx * dx + dy * dy).sqrt() / min_radius;
        let alpha = mod2pi(start.theta - dx.atan2(dy) + PI / 2.0); // local start heading
        let beta = mod2pi(goal.theta - dx.atan2(dy) + PI / 2.0); // local goal heading

        let candidates = [
            (DubinsWord::Lsl, dubins_lsl(d, alpha, beta)),
            (DubinsWord::Rsr, dubins_rsr(d, alpha, beta)),
            (DubinsWord::Lsr, dubins_lsr(d, alpha, beta)),
            (DubinsWord::Rsl, dubins_rsl(d, alpha, beta)),
            (DubinsWord::Rlr, dubins_rlr(d, alpha, beta)),
            (DubinsWord::Lrl, dubins_lrl(d, alpha, beta)),
        ];

        // Pick shortest valid path
        let best = candidates
            .iter()
            .filter_map(|(word, seg)| seg.map(|s| (*word, s)))
            .min_by(|(_, a), (_, b)| {
                let la: f64 = a[0] + a[1] + a[2];
                let lb: f64 = b[0] + b[1] + b[2];
                la.partial_cmp(&lb).unwrap_or(Ordering::Equal)
            });

        let (word, segs) = match best {
            Some(v) => v,
            None => return vec![[start.x, start.y], [goal.x, goal.y]],
        };

        sample_dubins(start, goal, word, segs, min_radius, 0.05)
    }
}

// ---------------------------------------------------------------------------
// Dubins helper functions
// ---------------------------------------------------------------------------

fn mod2pi(x: f64) -> f64 {
    let mut v = x % (2.0 * std::f64::consts::PI);
    if v < 0.0 {
        v += 2.0 * std::f64::consts::PI;
    }
    v
}

fn dubins_lsl(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp0 = d + sa - sb;
    let p_sq = 2.0 + d * d - 2.0 * (d * (sa - sb) - 2.0 * ca * cb + 2.0);
    // LSL: t = mod2pi(-alpha + atan2(cb-ca, d+sa-sb)), p = sqrt(p_sq), q = mod2pi(beta - ...)
    let tmp1 = (cb - ca).atan2(tmp0);
    let t = mod2pi(-alpha + tmp1);
    let p_sq2 = 2.0 + d * d - 2.0 * ca * cb * 2.0 + 2.0 * d * (sa - sb);
    if p_sq2 < 0.0 {
        return None;
    }
    let _ = p_sq; // suppress unused warning
    let p = p_sq2.max(0.0).sqrt();
    let q = mod2pi(beta - tmp1);
    Some([t, p, q])
}

fn dubins_rsr(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp0 = d - sa + sb;
    let tmp1 = (ca - cb).atan2(tmp0);
    let t = mod2pi(alpha - tmp1);
    let p_sq = 2.0 + d * d - 2.0 * ca * cb * 2.0 - 2.0 * d * (sa - sb);
    if p_sq < 0.0 {
        return None;
    }
    let p = p_sq.max(0.0).sqrt();
    let q = mod2pi(-beta + tmp1);
    Some([t, p, q])
}

fn dubins_lsr(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let p_sq = -2.0 + d * d + 2.0 * ca * cb * 2.0 + 2.0 * d * (sa + sb);
    if p_sq < 0.0 {
        return None;
    }
    let p = p_sq.max(0.0).sqrt();
    let tmp0 = (-ca - cb).atan2(d + sa + sb) - (-2.0_f64).atan2(p);
    let t = mod2pi(-alpha + tmp0);
    let q = mod2pi(-beta + tmp0);
    Some([t, p, q])
}

fn dubins_rsl(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let p_sq = d * d - 2.0 + 2.0 * ca * cb * 2.0 - 2.0 * d * (sa + sb);
    if p_sq < 0.0 {
        return None;
    }
    let p = p_sq.max(0.0).sqrt();
    let tmp0 = (ca + cb).atan2(d - sa - sb) - (2.0_f64).atan2(p);
    let t = mod2pi(alpha - tmp0);
    let q = mod2pi(beta - tmp0);
    Some([t, p, q])
}

fn dubins_rlr(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp0 = (6.0 - d * d + 2.0 * ca * cb * 2.0 + 2.0 * d * (sa - sb)) / 8.0;
    if tmp0.abs() > 1.0 {
        return None;
    }
    let p = mod2pi(2.0 * std::f64::consts::PI - tmp0.acos());
    let t = mod2pi(alpha - (ca - cb).atan2(d - sa + sb) + p / 2.0);
    let q = mod2pi(alpha - beta - t + p);
    Some([t, p, q])
}

fn dubins_lrl(d: f64, alpha: f64, beta: f64) -> Option<[f64; 3]> {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp0 = (6.0 - d * d + 2.0 * ca * cb * 2.0 - 2.0 * d * (sa - sb)) / 8.0;
    if tmp0.abs() > 1.0 {
        return None;
    }
    let p = mod2pi(2.0 * std::f64::consts::PI - tmp0.acos());
    let t = mod2pi(-alpha + (cb - ca).atan2(d + sa - sb) + p / 2.0);
    let q = mod2pi(beta - alpha - t + p);
    Some([t, p, q])
}

/// Sample world-space waypoints along a Dubins path.
fn sample_dubins(
    start: PathPlannerState,
    _goal: PathPlannerState,
    word: DubinsWord,
    segs: [f64; 3],
    r: f64,
    ds: f64,
) -> Vec<[f64; 2]> {
    // Each segment type: L = left turn (+), R = right turn (-), S = straight
    let types: [(i8, i8, i8); 6] = [
        (1, 0, 1),   // LSL
        (-1, 0, -1), // RSR
        (1, 0, -1),  // LSR
        (-1, 0, 1),  // RSL
        (-1, 1, -1), // RLR
        (1, -1, 1),  // LRL
    ];
    let idx = match word {
        DubinsWord::Lsl => 0,
        DubinsWord::Rsr => 1,
        DubinsWord::Lsr => 2,
        DubinsWord::Rsl => 3,
        DubinsWord::Rlr => 4,
        DubinsWord::Lrl => 5,
    };
    let (s0, s1, s2) = (types[idx].0, types[idx].1, types[idx].2);
    // Arc lengths in world units
    let lengths = [segs[0] * r, segs[1] * r, segs[2] * r];

    let mut pts: Vec<[f64; 2]> = Vec::new();
    let mut x = start.x;
    let mut y = start.y;
    let mut theta = start.theta;

    for (seg_idx, &seg_len) in lengths.iter().enumerate() {
        let turn = match seg_idx {
            0 => s0,
            1 => s1,
            _ => s2,
        };
        let steps = ((seg_len.abs() / ds).ceil() as usize).max(1);
        let step_len = seg_len / steps as f64;
        for _ in 0..steps {
            pts.push([x, y]);
            if turn == 0 {
                // Straight
                x += step_len * theta.cos();
                y += step_len * theta.sin();
            } else {
                // Arc: turn == 1 → left (CCW), turn == -1 → right (CW)
                let dtheta = step_len / r * turn as f64;
                // Centre of curvature
                let cx = x - r * (turn as f64) * theta.sin();
                let cy = y + r * (turn as f64) * theta.cos();
                theta += dtheta;
                x = cx + r * (turn as f64) * theta.sin();
                y = cy - r * (turn as f64) * theta.cos();
            }
        }
    }
    pts.push([x, y]);
    pts
}

// ---------------------------------------------------------------------------
// Small geometry helpers
// ---------------------------------------------------------------------------

#[inline]
fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

#[inline]
fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    dist2(a, b).sqrt()
}

// ---------------------------------------------------------------------------
// Road graph for Dijkstra
// ---------------------------------------------------------------------------

/// A road network represented as a weighted directed graph.
///
/// Nodes are identified by integer IDs.  Edges store the world-space distance
/// between adjacent nodes.
#[derive(Debug, Clone)]
pub struct RoadGraph {
    /// Node positions (world coordinates).
    pub nodes: Vec<[f64; 2]>,
    /// Adjacency list: `edges[i]` is a list of `(j, weight)` pairs for edges
    /// from node `i` to node `j`.
    pub edges: Vec<Vec<(usize, f64)>>,
}

impl Default for RoadGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl RoadGraph {
    /// Create an empty road graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Add a node at world position `pos`.  Returns the new node index.
    pub fn add_node(&mut self, pos: [f64; 2]) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(pos);
        self.edges.push(Vec::new());
        idx
    }

    /// Add a directed edge from `from` to `to` with the given `weight`.
    ///
    /// If `weight` is negative the Euclidean distance is used.
    pub fn add_edge(&mut self, from: usize, to: usize, weight: f64) {
        let w = if weight < 0.0 && from < self.nodes.len() && to < self.nodes.len() {
            dist(self.nodes[from], self.nodes[to])
        } else {
            weight
        };
        if from < self.edges.len() {
            self.edges[from].push((to, w));
        }
    }

    /// Add an undirected edge (two directed edges with the same weight).
    pub fn add_undirected_edge(&mut self, a: usize, b: usize, weight: f64) {
        self.add_edge(a, b, weight);
        self.add_edge(b, a, weight);
    }

    /// Build a road graph from a `GridMap` where free cells become nodes and
    /// 8-connected neighbours become edges (diagonal cost = √2).
    pub fn from_grid(map: &GridMap) -> Self {
        let n = map.width * map.height;
        let mut nodes = Vec::with_capacity(n);
        let mut edges = vec![Vec::new(); n];

        for y in 0..map.height {
            for x in 0..map.width {
                nodes.push(map.grid_to_world(x, y));
            }
        }

        for y in 0..map.height {
            for x in 0..map.width {
                if map.cells[map.idx(x, y)] == 255 {
                    continue;
                }
                let from = map.idx(x, y);
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        if dx == 0 && dy == 0 {
                            continue;
                        }
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;
                        if nx < 0 || ny < 0 || nx >= map.width as i32 || ny >= map.height as i32 {
                            continue;
                        }
                        let to = map.idx(nx as usize, ny as usize);
                        if map.cells[to] == 255 {
                            continue;
                        }
                        let w = if dx != 0 && dy != 0 {
                            map.resolution * 2_f64.sqrt()
                        } else {
                            map.resolution
                        };
                        edges[from].push((to, w));
                    }
                }
            }
        }

        Self { nodes, edges }
    }

    /// Number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

// ---------------------------------------------------------------------------
// Dijkstra planner on road graph
// ---------------------------------------------------------------------------

/// Dijkstra shortest-path planner on a `RoadGraph`.
pub struct DijkstraPlanner;

impl DijkstraPlanner {
    /// Find the shortest path from `start_node` to `goal_node`.
    ///
    /// Returns `Some(node_indices)` or `None` if unreachable.
    pub fn plan(graph: &RoadGraph, start: usize, goal: usize) -> Option<Vec<usize>> {
        let n = graph.nodes.len();
        if start >= n || goal >= n {
            return None;
        }

        let mut dist_arr = vec![f64::INFINITY; n];
        let mut prev = vec![usize::MAX; n];
        dist_arr[start] = 0.0;

        // Min-heap: (cost_as_ordered_float, node_index)
        let mut heap: BinaryHeap<(ordered_float::OrderedFloat, usize)> = BinaryHeap::new();
        heap.push((ordered_float::OrderedFloat(-0.0), start));

        while let Some((neg_d, u)) = heap.pop() {
            let d = -neg_d.0;
            if d > dist_arr[u] {
                continue;
            }
            if u == goal {
                break;
            }
            for &(v, w) in &graph.edges[u] {
                let nd = d + w;
                if nd < dist_arr[v] {
                    dist_arr[v] = nd;
                    prev[v] = u;
                    heap.push((ordered_float::OrderedFloat(-nd), v));
                }
            }
        }

        if dist_arr[goal].is_infinite() {
            return None;
        }

        let mut path = Vec::new();
        let mut cur = goal;
        while cur != usize::MAX {
            path.push(cur);
            if cur == start {
                break;
            }
            cur = prev[cur];
        }
        path.reverse();
        if path.first() == Some(&start) {
            Some(path)
        } else {
            None
        }
    }

    /// Plan and return world-space waypoints for the found path.
    pub fn plan_waypoints(graph: &RoadGraph, start: usize, goal: usize) -> Option<Vec<[f64; 2]>> {
        let nodes = Self::plan(graph, start, goal)?;
        Some(nodes.iter().map(|&i| graph.nodes[i]).collect())
    }
}

/// Minimal ordered-float wrapper so we can push negative f64s into a BinaryHeap.
mod ordered_float {
    #[derive(Clone, Copy, PartialEq)]
    pub struct OrderedFloat(pub f64);

    impl Eq for OrderedFloat {}

    impl PartialOrd for OrderedFloat {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for OrderedFloat {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            self.0
                .partial_cmp(&other.0)
                .unwrap_or(std::cmp::Ordering::Equal)
        }
    }
}

// ---------------------------------------------------------------------------
// Bézier path smoothing
// ---------------------------------------------------------------------------

/// Path smoothing utilities using cubic Bézier curves.
///
/// Given a set of waypoints, this module fits a C1-continuous piecewise
/// cubic Bézier spline and re-samples it at uniform arc-length intervals.
pub struct BezierSmoother;

impl BezierSmoother {
    /// Evaluate a cubic Bézier at parameter `t ∈ [0,1]`.
    ///
    /// Points: `p0, p1, p2, p3` (control points).
    pub fn cubic_bezier(
        p0: [f64; 2],
        p1: [f64; 2],
        p2: [f64; 2],
        p3: [f64; 2],
        t: f64,
    ) -> [f64; 2] {
        let u = 1.0 - t;
        let b0 = u * u * u;
        let b1 = 3.0 * u * u * t;
        let b2 = 3.0 * u * t * t;
        let b3 = t * t * t;
        [
            b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0],
            b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1],
        ]
    }

    /// Smooth a polyline by fitting Catmull-Rom splines converted to Bézier form.
    ///
    /// Returns a densely-sampled path.  `samples_per_segment` controls output
    /// resolution.
    pub fn smooth(path: &[[f64; 2]], samples_per_segment: usize) -> Vec<[f64; 2]> {
        if path.len() < 2 {
            return path.to_vec();
        }
        let n = path.len();
        let sps = samples_per_segment.max(2);
        let mut result = Vec::with_capacity((n - 1) * sps + 1);

        for i in 0..n - 1 {
            let p0 = if i == 0 { path[0] } else { path[i - 1] };
            let p1 = path[i];
            let p2 = path[i + 1];
            let p3 = if i + 2 < n { path[i + 2] } else { path[n - 1] };

            // Catmull-Rom to Bézier control points (tension = 0.5)
            let cp1 = [p1[0] + (p2[0] - p0[0]) / 6.0, p1[1] + (p2[1] - p0[1]) / 6.0];
            let cp2 = [p2[0] - (p3[0] - p1[0]) / 6.0, p2[1] - (p3[1] - p1[1]) / 6.0];

            let skip_last = i < n - 2;
            let count = if skip_last { sps } else { sps + 1 };
            for j in 0..count {
                let t = j as f64 / sps as f64;
                result.push(Self::cubic_bezier(p1, cp1, cp2, p2, t));
            }
        }
        result
    }

    /// Compute the curvature at parameter `t` of a cubic Bézier.
    pub fn curvature(p0: [f64; 2], p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], t: f64) -> f64 {
        // First and second derivatives
        let u = 1.0 - t;
        let d1 = [
            3.0 * (u * u * (p1[0] - p0[0])
                + 2.0 * u * t * (p2[0] - p1[0])
                + t * t * (p3[0] - p2[0])),
            3.0 * (u * u * (p1[1] - p0[1])
                + 2.0 * u * t * (p2[1] - p1[1])
                + t * t * (p3[1] - p2[1])),
        ];
        let d2 = [
            6.0 * (u * (p2[0] - 2.0 * p1[0] + p0[0]) + t * (p3[0] - 2.0 * p2[0] + p1[0])),
            6.0 * (u * (p2[1] - 2.0 * p1[1] + p0[1]) + t * (p3[1] - 2.0 * p2[1] + p1[1])),
        ];
        let cross = d1[0] * d2[1] - d1[1] * d2[0];
        let mag1 = (d1[0] * d1[0] + d1[1] * d1[1]).powf(1.5);
        if mag1 < 1e-30 { 0.0 } else { cross / mag1 }
    }
}

// ---------------------------------------------------------------------------
// Collision-aware path planning
// ---------------------------------------------------------------------------

/// A path planner that checks for moving obstacles at each step and replans
/// if a collision is imminent.
#[derive(Debug, Clone)]
pub struct CollisionAwarePlanner {
    /// Safety clearance around each obstacle (m).
    pub safety_radius: f64,
    /// Maximum lookahead horizon (m).
    pub lookahead: f64,
}

/// A circular obstacle with a position and radius.
#[derive(Debug, Clone, Copy)]
pub struct CircularObstacle {
    /// Centre position (world frame).
    pub center: [f64; 2],
    /// Radius of the obstacle (m).
    pub radius: f64,
}

impl CollisionAwarePlanner {
    /// Create a planner with given safety clearance and lookahead.
    pub fn new(safety_radius: f64, lookahead: f64) -> Self {
        Self {
            safety_radius,
            lookahead,
        }
    }

    /// Check whether the given `path` segment (from `start` to `end`) is
    /// clear of all `obstacles`.
    pub fn segment_is_clear(
        &self,
        start: [f64; 2],
        end: [f64; 2],
        obstacles: &[CircularObstacle],
    ) -> bool {
        for obs in obstacles {
            // Find closest point on segment to obstacle centre
            let dx = end[0] - start[0];
            let dy = end[1] - start[1];
            let len_sq = dx * dx + dy * dy;
            let t = if len_sq > 1e-20 {
                let dot = (obs.center[0] - start[0]) * dx + (obs.center[1] - start[1]) * dy;
                (dot / len_sq).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let closest = [start[0] + t * dx, start[1] + t * dy];
            let d = dist(closest, obs.center);
            if d < obs.radius + self.safety_radius {
                return false;
            }
        }
        true
    }

    /// Filter a path to remove waypoints inside obstacles.
    ///
    /// Keeps only waypoints that lie outside the inflated obstacle radii.
    pub fn filter_path(&self, path: &[[f64; 2]], obstacles: &[CircularObstacle]) -> Vec<[f64; 2]> {
        path.iter()
            .filter(|&&pt| {
                obstacles
                    .iter()
                    .all(|obs| dist(pt, obs.center) > obs.radius + self.safety_radius)
            })
            .copied()
            .collect()
    }

    /// Plan around obstacles: if the direct path is blocked, insert a detour
    /// waypoint that avoids the obstacle.
    ///
    /// Returns waypoints in world coordinates.
    pub fn plan_around(
        &self,
        start: [f64; 2],
        goal: [f64; 2],
        obstacles: &[CircularObstacle],
    ) -> Vec<[f64; 2]> {
        if self.segment_is_clear(start, goal, obstacles) {
            return vec![start, goal];
        }

        // Find the first blocking obstacle
        let blocking = obstacles
            .iter()
            .find(|obs| !self.segment_is_clear(start, goal, std::slice::from_ref(obs)));

        match blocking {
            None => vec![start, goal],
            Some(obs) => {
                // Create a tangential detour to the left of the obstacle
                let dx = goal[0] - start[0];
                let dy = goal[1] - start[1];
                let d = (dx * dx + dy * dy).sqrt().max(1e-12);
                let perp = [-dy / d, dx / d]; // left perpendicular
                let offset = obs.radius + self.safety_radius + 0.5;
                let detour = [
                    obs.center[0] + perp[0] * offset,
                    obs.center[1] + perp[1] * offset,
                ];
                vec![start, detour, goal]
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Dynamic replanning
// ---------------------------------------------------------------------------

/// Lightweight dynamic replanner that triggers A* replanning when the tracked
/// path deviates from the current vehicle position beyond a threshold.
pub struct DynamicReplanner {
    /// Distance threshold that triggers replanning (m).
    pub deviation_threshold: f64,
    /// Cached current path.
    pub current_path: Vec<[f64; 2]>,
    /// Current waypoint index the vehicle is heading toward.
    pub target_idx: usize,
}

impl DynamicReplanner {
    /// Create a new dynamic replanner with the given deviation threshold.
    pub fn new(deviation_threshold: f64) -> Self {
        Self {
            deviation_threshold,
            current_path: Vec::new(),
            target_idx: 0,
        }
    }

    /// Update the planner given the current vehicle position.
    ///
    /// Returns `true` if replanning was triggered.
    pub fn update(&mut self, pos: [f64; 2], goal: [f64; 2], map: &GridMap) -> bool {
        // Advance target index to nearest on-path waypoint
        if !self.current_path.is_empty() {
            let best = self
                .current_path
                .iter()
                .enumerate()
                .skip(self.target_idx)
                .min_by(|(_, a), (_, b)| {
                    dist2(pos, **a)
                        .partial_cmp(&dist2(pos, **b))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            if let Some((idx, _)) = best {
                self.target_idx = idx;
            }
        }

        // Check deviation
        let path_point = if self.current_path.is_empty() {
            pos
        } else {
            self.current_path[self.target_idx]
        };
        let deviation = dist(pos, path_point);

        if deviation > self.deviation_threshold || self.current_path.is_empty() {
            // Trigger replan
            if let Some(new_path) = AStarPlanner::plan(map, pos, goal) {
                self.current_path = new_path;
                self.target_idx = 0;
            }
            return true;
        }
        false
    }

    /// Returns the next target waypoint for the vehicle to steer toward.
    pub fn next_waypoint(&self) -> Option<[f64; 2]> {
        let idx = (self.target_idx + 1).min(self.current_path.len().saturating_sub(1));
        self.current_path.get(idx).copied()
    }

    /// Check whether the goal has been reached.
    pub fn goal_reached(&self, pos: [f64; 2], goal_radius: f64) -> bool {
        match self.current_path.last() {
            Some(&last) => dist(pos, last) <= goal_radius,
            None => false,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- GridMap ---

    #[test]
    fn test_grid_world_roundtrip() {
        let map = GridMap::new(20, 20, 0.5, [1.0, 2.0]);
        let wx = 3.75;
        let wy = 4.25;
        let (gx, gy) = map.world_to_grid(wx, wy).expect("in bounds");
        let world = map.grid_to_world(gx, gy);
        // grid_to_world returns cell centre, so check within half a cell
        assert!((world[0] - wx).abs() < map.resolution);
        assert!((world[1] - wy).abs() < map.resolution);
    }

    #[test]
    fn test_world_to_grid_oob_returns_none() {
        let map = GridMap::new(10, 10, 1.0, [0.0, 0.0]);
        assert!(map.world_to_grid(-1.0, 0.0).is_none());
        assert!(map.world_to_grid(0.0, 11.0).is_none());
    }

    // --- A* ---

    #[test]
    fn test_astar_finds_path_no_obstacles() {
        let map = GridMap::new(5, 5, 1.0, [0.0, 0.0]);
        let start = [0.5, 0.5];
        let goal = [4.5, 4.5];
        let path = AStarPlanner::plan(&map, start, goal);
        assert!(path.is_some(), "A* must find a path in an empty 5x5 grid");
        let path = path.unwrap();
        assert!(path.len() >= 2);
        // First waypoint close to start, last close to goal
        assert!(dist(path[0], start) < 1.5);
        assert!(dist(*path.last().unwrap(), goal) < 1.5);
    }

    #[test]
    fn test_astar_none_when_start_is_obstacle() {
        let mut map = GridMap::new(5, 5, 1.0, [0.0, 0.0]);
        map.set_obstacle(0, 0);
        let path = AStarPlanner::plan(&map, [0.5, 0.5], [4.5, 4.5]);
        assert!(path.is_none());
    }

    #[test]
    fn test_astar_none_when_goal_is_obstacle() {
        let mut map = GridMap::new(5, 5, 1.0, [0.0, 0.0]);
        map.set_obstacle(4, 4);
        let path = AStarPlanner::plan(&map, [0.5, 0.5], [4.5, 4.5]);
        assert!(path.is_none());
    }

    // --- PathSmoother ---

    #[test]
    fn test_path_length_positive() {
        let path: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]];
        let len = PathSmoother::path_length(&path);
        assert!(len > 0.0);
    }

    #[test]
    fn test_smooth_path_endpoints_fixed() {
        let path: Vec<[f64; 2]> = (0..10).map(|i| [i as f64, 0.0]).collect();
        let smoothed = PathSmoother::smooth_path(&path, 0.5, 0.1, 1e-6);
        // Endpoints must remain unchanged
        let eps = 1e-9;
        assert!((smoothed[0][0] - path[0][0]).abs() < eps);
        assert!((smoothed[0][1] - path[0][1]).abs() < eps);
        let n = path.len() - 1;
        assert!((smoothed[n][0] - path[n][0]).abs() < eps);
        assert!((smoothed[n][1] - path[n][1]).abs() < eps);
    }

    #[test]
    fn test_resample_path_uniform_spacing() {
        let path: Vec<[f64; 2]> = vec![[0.0, 0.0], [10.0, 0.0]];
        let spacing = 1.0;
        let resampled = PathSmoother::resample_path(&path, spacing);
        // All consecutive spacings should be approximately equal (within 1%)
        for w in resampled.windows(2) {
            let d = dist(w[0], w[1]);
            // Last segment may be shorter; check all but the last
            let _ = d;
        }
        // Internal spacings should be very close to `spacing`
        let n = resampled.len();
        assert!(n >= 2, "must have at least two points");
        for w in resampled[..n - 1].windows(2) {
            let d = dist(w[0], w[1]);
            assert!(
                (d - spacing).abs() < 1e-9,
                "spacing was {d}, expected {spacing}"
            );
        }
    }

    // --- PotentialField ---

    #[test]
    fn test_potential_field_attractive_force_toward_goal() {
        let pf = PotentialField {
            attractive_gain: 1.0,
            repulsive_gain: 0.0,
            repulsive_influence: 1.0,
        };
        let map = GridMap::new(10, 10, 1.0, [0.0, 0.0]);
        let pos = [2.0, 2.0];
        let goal = [5.0, 5.0];
        let force = pf.total_force(pos, goal, &map);
        // Force should point from pos towards goal
        assert!(force[0] > 0.0, "x-component of force must be positive");
        assert!(force[1] > 0.0, "y-component of force must be positive");
    }

    // --- RoadGraph ---

    #[test]
    fn test_road_graph_add_node_and_edge() {
        let mut g = RoadGraph::new();
        let a = g.add_node([0.0, 0.0]);
        let b = g.add_node([1.0, 0.0]);
        g.add_undirected_edge(a, b, 1.0);
        assert_eq!(g.node_count(), 2);
        assert!(!g.edges[a].is_empty());
        assert!(!g.edges[b].is_empty());
    }

    #[test]
    fn test_road_graph_from_grid_node_count() {
        let map = GridMap::new(3, 3, 1.0, [0.0, 0.0]);
        let g = RoadGraph::from_grid(&map);
        assert_eq!(g.node_count(), 9, "3x3 free grid → 9 nodes");
    }

    // --- Dijkstra ---

    #[test]
    fn test_dijkstra_finds_path_in_chain() {
        let mut g = RoadGraph::new();
        let a = g.add_node([0.0, 0.0]);
        let b = g.add_node([1.0, 0.0]);
        let c = g.add_node([2.0, 0.0]);
        g.add_undirected_edge(a, b, 1.0);
        g.add_undirected_edge(b, c, 1.0);
        let path = DijkstraPlanner::plan(&g, a, c);
        assert!(path.is_some(), "Dijkstra should find path in chain");
        assert_eq!(path.unwrap(), vec![a, b, c]);
    }

    #[test]
    fn test_dijkstra_no_path_disconnected() {
        let mut g = RoadGraph::new();
        let a = g.add_node([0.0, 0.0]);
        let b = g.add_node([10.0, 0.0]);
        let _ = (a, b);
        // No edges
        let path = DijkstraPlanner::plan(&g, a, b);
        assert!(path.is_none(), "disconnected nodes should return None");
    }

    #[test]
    fn test_dijkstra_shortest_path_two_routes() {
        let mut g = RoadGraph::new();
        let a = g.add_node([0.0, 0.0]);
        let b = g.add_node([1.0, 0.0]);
        let c = g.add_node([2.0, 0.0]);
        let d = g.add_node([1.0, 1.0]);
        g.add_undirected_edge(a, b, 1.0);
        g.add_undirected_edge(b, c, 1.0);
        g.add_undirected_edge(a, d, 5.0);
        g.add_undirected_edge(d, c, 5.0);
        let path = DijkstraPlanner::plan(&g, a, c).unwrap();
        // Shortest should go a→b→c (cost=2), not a→d→c (cost=10)
        assert_eq!(path, vec![a, b, c], "should prefer short route");
    }

    // --- BezierSmoother ---

    #[test]
    fn test_bezier_cubic_at_endpoints() {
        let p0 = [0.0_f64, 0.0];
        let p1 = [1.0, 0.0];
        let p2 = [2.0, 1.0];
        let p3 = [3.0, 0.0];
        let s = BezierSmoother::cubic_bezier(p0, p1, p2, p3, 0.0);
        let e = BezierSmoother::cubic_bezier(p0, p1, p2, p3, 1.0);
        assert!((s[0] - p0[0]).abs() < 1e-9 && (s[1] - p0[1]).abs() < 1e-9);
        assert!((e[0] - p3[0]).abs() < 1e-9 && (e[1] - p3[1]).abs() < 1e-9);
    }

    #[test]
    fn test_bezier_smooth_more_points_than_input() {
        let path: Vec<[f64; 2]> = vec![[0.0, 0.0], [1.0, 0.5], [2.0, 0.0], [3.0, 0.5]];
        let smoothed = BezierSmoother::smooth(&path, 10);
        assert!(
            smoothed.len() > path.len(),
            "smooth should produce more points"
        );
    }

    #[test]
    fn test_bezier_curvature_zero_on_straight() {
        // A straight line: all points collinear → curvature should be ~0
        let p0 = [0.0_f64, 0.0];
        let p1 = [1.0, 0.0];
        let p2 = [2.0, 0.0];
        let p3 = [3.0, 0.0];
        let k = BezierSmoother::curvature(p0, p1, p2, p3, 0.5).abs();
        assert!(k < 1e-6, "straight line curvature should be ~0, got {k}");
    }

    // --- CollisionAwarePlanner ---

    #[test]
    fn test_collision_planner_clear_segment() {
        let planner = CollisionAwarePlanner::new(0.5, 10.0);
        let obs = CircularObstacle {
            center: [5.0, 5.0],
            radius: 1.0,
        };
        // Path well away from obstacle
        let clear = planner.segment_is_clear([0.0, 0.0], [1.0, 0.0], &[obs]);
        assert!(clear, "segment far from obstacle should be clear");
    }

    #[test]
    fn test_collision_planner_blocked_segment() {
        let planner = CollisionAwarePlanner::new(0.5, 10.0);
        let obs = CircularObstacle {
            center: [1.0, 0.0],
            radius: 0.3,
        };
        // Path passes right through obstacle
        let blocked = !planner.segment_is_clear([0.0, 0.0], [2.0, 0.0], &[obs]);
        assert!(blocked, "segment through obstacle should be blocked");
    }

    #[test]
    fn test_collision_planner_plan_around_inserts_detour() {
        let planner = CollisionAwarePlanner::new(0.5, 10.0);
        let obs = CircularObstacle {
            center: [1.0, 0.0],
            radius: 0.4,
        };
        let path = planner.plan_around([0.0, 0.0], [2.0, 0.0], &[obs]);
        // Should have more than 2 points (i.e. a detour was inserted)
        assert!(path.len() > 2, "detour should add an intermediate waypoint");
    }

    #[test]
    fn test_collision_planner_filter_removes_blocked() {
        let planner = CollisionAwarePlanner::new(0.1, 5.0);
        let obs = CircularObstacle {
            center: [1.0, 0.0],
            radius: 0.2,
        };
        let path = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let filtered = planner.filter_path(&path, &[obs]);
        assert!(
            filtered
                .iter()
                .all(|&pt| dist(pt, obs.center) > obs.radius + 0.1),
            "all remaining points should be outside obstacle"
        );
    }

    // --- DynamicReplanner ---

    #[test]
    fn test_dynamic_replanner_triggers_on_empty_path() {
        let mut replanner = DynamicReplanner::new(1.0);
        let map = GridMap::new(10, 10, 1.0, [0.0, 0.0]);
        let triggered = replanner.update([0.5, 0.5], [9.5, 9.5], &map);
        assert!(triggered, "empty path should always trigger replan");
    }

    #[test]
    fn test_dynamic_replanner_has_path_after_update() {
        let mut replanner = DynamicReplanner::new(1.0);
        let map = GridMap::new(10, 10, 1.0, [0.0, 0.0]);
        replanner.update([0.5, 0.5], [9.5, 9.5], &map);
        assert!(
            !replanner.current_path.is_empty(),
            "should have a path after replan"
        );
    }

    #[test]
    fn test_dynamic_replanner_goal_not_reached_at_start() {
        let mut replanner = DynamicReplanner::new(1.0);
        let map = GridMap::new(10, 10, 1.0, [0.0, 0.0]);
        replanner.update([0.5, 0.5], [9.5, 9.5], &map);
        assert!(
            !replanner.goal_reached([0.5, 0.5], 0.5),
            "should not be at goal from start"
        );
    }
}
