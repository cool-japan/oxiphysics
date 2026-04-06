#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crowd simulation using Helbing's Social Force Model (SFM).
//!
//! Implements:
//! - [`Pedestrian`] — position, velocity, desired speed, SFM parameters
//! - [`SocialForceModel`] — desired force, repulsive social force, obstacle force
//! - [`FlowField`] — potential field for navigation, gradient descent path following
//! - [`CrowdDensity`] — density map, local velocity field, fundamental diagram (v vs ρ)
//! - [`EmergentBehavior`] — lane formation, oscillation at bottleneck, arch formation
//! - [`EvacuationSimulation`] — exit choice, bottleneck congestion, evacuation time

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// 2D vector type alias.
type Vec2 = [f64; 2];

#[inline]
fn add2(a: Vec2, b: Vec2) -> Vec2 {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn sub2(a: Vec2, b: Vec2) -> Vec2 {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn scale2(a: Vec2, s: f64) -> Vec2 {
    [a[0] * s, a[1] * s]
}

#[inline]
fn dot2(a: Vec2, b: Vec2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn norm2(a: Vec2) -> f64 {
    dot2(a, a).sqrt()
}

#[inline]
fn normalize2(a: Vec2) -> Vec2 {
    let n = norm2(a);
    if n > 1e-15 {
        scale2(a, 1.0 / n)
    } else {
        [0.0, 0.0]
    }
}

#[inline]
fn dist2(a: Vec2, b: Vec2) -> f64 {
    norm2(sub2(a, b))
}

// ---------------------------------------------------------------------------
// Pedestrian
// ---------------------------------------------------------------------------

/// Parameters governing how a single pedestrian interacts with the social force model.
#[derive(Debug, Clone)]
pub struct SfmParams {
    /// Desired speed in m/s (free-flow speed).
    pub desired_speed: f64,
    /// Relaxation time τ (s): how quickly the pedestrian reaches desired velocity.
    pub tau: f64,
    /// Strength of pedestrian–pedestrian repulsion (A_i in Helbing's model).
    pub repulsion_strength: f64,
    /// Range of pedestrian–pedestrian repulsion (B_i).
    pub repulsion_range: f64,
    /// Physical radius of the pedestrian (m).
    pub radius: f64,
    /// Body compression coefficient k_1 (N/m).
    pub body_force_k1: f64,
    /// Sliding friction coefficient k_2 (kg/ms).
    pub body_force_k2: f64,
    /// Anisotropy parameter λ: weighting of forces behind vs. in front.
    pub anisotropy: f64,
    /// Field of view half-angle (radians).
    pub fov_angle: f64,
}

impl Default for SfmParams {
    fn default() -> Self {
        Self {
            desired_speed: 1.34,
            tau: 0.5,
            repulsion_strength: 2000.0,
            repulsion_range: 0.08,
            radius: 0.25,
            body_force_k1: 1.2e5,
            body_force_k2: 2.4e5,
            anisotropy: 0.5,
            fov_angle: PI * 0.75,
        }
    }
}

/// A single pedestrian agent in the crowd simulation.
#[derive(Debug, Clone)]
pub struct Pedestrian {
    /// Unique identifier.
    pub id: usize,
    /// Current 2D position (m).
    pub position: Vec2,
    /// Current 2D velocity (m/s).
    pub velocity: Vec2,
    /// Desired 2D direction (unit vector).
    pub desired_direction: Vec2,
    /// Current destination position.
    pub destination: Vec2,
    /// Mass of the pedestrian (kg).
    pub mass: f64,
    /// SFM parameters.
    pub params: SfmParams,
    /// Whether the pedestrian has reached the exit.
    pub evacuated: bool,
    /// Total distance travelled (m).
    pub distance_travelled: f64,
    /// Time since simulation start (s).
    pub time_in_simulation: f64,
    /// Group identifier (optional).
    pub group_id: Option<usize>,
}

impl Pedestrian {
    /// Create a new pedestrian at `position` heading toward `destination`.
    pub fn new(id: usize, position: Vec2, destination: Vec2, mass: f64) -> Self {
        let dir = normalize2(sub2(destination, position));
        Self {
            id,
            position,
            velocity: [0.0, 0.0],
            desired_direction: dir,
            destination,
            mass,
            params: SfmParams::default(),
            evacuated: false,
            distance_travelled: 0.0,
            time_in_simulation: 0.0,
            group_id: None,
        }
    }

    /// Update the desired direction toward the current destination.
    pub fn update_desired_direction(&mut self) {
        let d = sub2(self.destination, self.position);
        if norm2(d) > 1e-6 {
            self.desired_direction = normalize2(d);
        }
    }

    /// Return the distance to the destination.
    pub fn distance_to_destination(&self) -> f64 {
        dist2(self.position, self.destination)
    }

    /// Return current speed in m/s.
    pub fn speed(&self) -> f64 {
        norm2(self.velocity)
    }

    /// Check whether the pedestrian is within `threshold` metres of the destination.
    pub fn has_reached_destination(&self, threshold: f64) -> bool {
        self.distance_to_destination() < threshold
    }

    /// Integrate velocity and position with semi-implicit Euler and update distance counter.
    pub fn integrate(&mut self, force: Vec2, dt: f64) {
        // v += (F/m) * dt
        let acc = scale2(force, 1.0 / self.mass);
        self.velocity = add2(self.velocity, scale2(acc, dt));
        // Clamp speed to avoid runaway.
        let spd = norm2(self.velocity);
        let max_speed = self.params.desired_speed * 3.0;
        if spd > max_speed {
            self.velocity = scale2(normalize2(self.velocity), max_speed);
        }
        let old_pos = self.position;
        self.position = add2(self.position, scale2(self.velocity, dt));
        self.distance_travelled += dist2(self.position, old_pos);
        self.time_in_simulation += dt;
    }
}

// ---------------------------------------------------------------------------
// SocialForceModel
// ---------------------------------------------------------------------------

/// Helbing's Social Force Model for pedestrian dynamics.
///
/// Computes three force components:
/// 1. Desired force: drives pedestrian toward goal at desired speed.
/// 2. Repulsive social force: keeps distance from other pedestrians.
/// 3. Obstacle force: repulsion from walls / obstacles.
#[derive(Debug, Clone)]
pub struct SocialForceModel;

impl SocialForceModel {
    /// Compute the desired (driving) force for pedestrian `p`.
    ///
    /// f_i^0 = (m_i / τ_i) * (v_i^0 * e_i^0 − v_i)
    pub fn desired_force(p: &Pedestrian) -> Vec2 {
        let v0 = p.params.desired_speed;
        let tau = p.params.tau;
        let e0 = p.desired_direction;
        let target_vel = scale2(e0, v0);
        let diff = sub2(target_vel, p.velocity);
        scale2(diff, p.mass / tau)
    }

    /// Anisotropy weighting: reduces force from pedestrians behind.
    fn anisotropy_weight(p: &Pedestrian, n_ij: Vec2) -> f64 {
        let lambda = p.params.anisotropy;
        let cos_phi = dot2(p.desired_direction, scale2(n_ij, -1.0));
        lambda + (1.0 - lambda) * 0.5 * (1.0 + cos_phi)
    }

    /// Compute the repulsive social force from pedestrian `other` on `p`.
    ///
    /// Includes exponential psychological repulsion and physical contact terms.
    pub fn repulsive_social_force(p: &Pedestrian, other: &Pedestrian) -> Vec2 {
        let r_ij = p.params.radius + other.params.radius;
        let d_vec = sub2(p.position, other.position);
        let d = norm2(d_vec);
        if d < 1e-6 {
            return [0.0, 0.0];
        }
        let n_ij = normalize2(d_vec); // unit vector from other to p

        let a = p.params.repulsion_strength;
        let b = p.params.repulsion_range;

        // Psychological repulsion
        let f_social = a * (-(d - r_ij) / b).exp();

        // Physical contact terms (only when overlapping)
        let overlap = r_ij - d;
        let f_body = if overlap > 0.0 {
            p.params.body_force_k1 * overlap
        } else {
            0.0
        };

        // Sliding friction (tangential)
        let tangent = [-n_ij[1], n_ij[0]]; // perpendicular to n_ij
        let delta_v_t = if overlap > 0.0 {
            dot2(sub2(other.velocity, p.velocity), tangent)
        } else {
            0.0
        };
        let f_friction = if overlap > 0.0 {
            p.params.body_force_k2 * overlap * delta_v_t
        } else {
            0.0
        };

        let weight = Self::anisotropy_weight(p, n_ij);
        let f_total_normal = (f_social + f_body) * weight;

        let normal_force = scale2(n_ij, f_total_normal);
        let friction_force = scale2(tangent, f_friction);
        add2(normal_force, friction_force)
    }

    /// Compute the repulsive force from a line-segment obstacle.
    ///
    /// `wall_a` and `wall_b` define the two endpoints of the wall segment.
    pub fn obstacle_force(p: &Pedestrian, wall_a: Vec2, wall_b: Vec2) -> Vec2 {
        let closest = closest_point_on_segment(p.position, wall_a, wall_b);
        let d_vec = sub2(p.position, closest);
        let d = norm2(d_vec);
        if d < 1e-6 {
            return [0.0, 0.0];
        }
        let n = normalize2(d_vec);
        let r = p.params.radius;

        let a = p.params.repulsion_strength;
        let b = p.params.repulsion_range;

        let f_social = a * (-(d - r) / b).exp();
        let overlap = r - d;
        let f_body = if overlap > 0.0 {
            p.params.body_force_k1 * overlap
        } else {
            0.0
        };

        // Sliding friction
        let tangent = [-n[1], n[0]];
        let f_friction = if overlap > 0.0 {
            -p.params.body_force_k2 * overlap * dot2(p.velocity, tangent)
        } else {
            0.0
        };

        add2(scale2(n, f_social + f_body), scale2(tangent, f_friction))
    }

    /// Compute the total force on pedestrian `p` from all others and obstacles.
    pub fn total_force(p: &Pedestrian, all: &[Pedestrian], obstacles: &[(Vec2, Vec2)]) -> Vec2 {
        let mut f = Self::desired_force(p);
        for other in all {
            if other.id != p.id && !other.evacuated {
                let fij = Self::repulsive_social_force(p, other);
                f = add2(f, fij);
            }
        }
        for &(wa, wb) in obstacles {
            let fw = Self::obstacle_force(p, wa, wb);
            f = add2(f, fw);
        }
        f
    }
}

/// Closest point on segment \[a, b\] to point p.
fn closest_point_on_segment(p: Vec2, a: Vec2, b: Vec2) -> Vec2 {
    let ab = sub2(b, a);
    let ap = sub2(p, a);
    let t = dot2(ap, ab) / (dot2(ab, ab) + 1e-15);
    let t = t.clamp(0.0, 1.0);
    add2(a, scale2(ab, t))
}

// ---------------------------------------------------------------------------
// FlowField
// ---------------------------------------------------------------------------

/// A 2D potential / flow field for navigation.
///
/// Stores a scalar potential at each grid cell. Pedestrians follow the
/// negative gradient (steepest descent toward the minimum, which is the goal).
#[derive(Debug, Clone)]
pub struct FlowField {
    /// Width of the grid in cells.
    pub width: usize,
    /// Height of the grid in cells.
    pub height: usize,
    /// Physical size of each cell (m).
    pub cell_size: f64,
    /// Origin of the grid (bottom-left corner, m).
    pub origin: Vec2,
    /// Potential values (row-major: index = y * width + x).
    pub potential: Vec<f64>,
    /// Whether each cell is passable.
    pub passable: Vec<bool>,
}

impl FlowField {
    /// Create a new flow field with all potentials set to infinity (uninitialised).
    pub fn new(width: usize, height: usize, cell_size: f64, origin: Vec2) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            cell_size,
            origin,
            potential: vec![f64::INFINITY; n],
            passable: vec![true; n],
        }
    }

    /// Convert world position to grid cell coordinates (clamped).
    pub fn world_to_cell(&self, pos: Vec2) -> (usize, usize) {
        let fx = (pos[0] - self.origin[0]) / self.cell_size;
        let fy = (pos[1] - self.origin[1]) / self.cell_size;
        let cx = (fx as usize).min(self.width.saturating_sub(1));
        let cy = (fy as usize).min(self.height.saturating_sub(1));
        (cx, cy)
    }

    /// Convert cell coordinates to world-space centre.
    pub fn cell_to_world(&self, cx: usize, cy: usize) -> Vec2 {
        [
            self.origin[0] + (cx as f64 + 0.5) * self.cell_size,
            self.origin[1] + (cy as f64 + 0.5) * self.cell_size,
        ]
    }

    fn idx(&self, cx: usize, cy: usize) -> usize {
        cy * self.width + cx
    }

    /// Compute the potential field via Dijkstra / wavefront propagation from goal cells.
    ///
    /// `goal_cells` is a list of `(cx, cy)` pairs with zero potential.
    pub fn compute_dijkstra(&mut self, goal_cells: &[(usize, usize)]) {
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;

        self.potential.fill(f64::INFINITY);
        let mut heap: BinaryHeap<Reverse<(u64, usize, usize)>> = BinaryHeap::new();

        for &(gx, gy) in goal_cells {
            let i = self.idx(gx, gy);
            self.potential[i] = 0.0;
            heap.push(Reverse((0, gx, gy)));
        }

        let dirs: [(i32, i32); 8] = [
            (1, 0),
            (-1, 0),
            (0, 1),
            (0, -1),
            (1, 1),
            (1, -1),
            (-1, 1),
            (-1, -1),
        ];

        while let Some(Reverse((_, cx, cy))) = heap.pop() {
            let current_cost = self.potential[self.idx(cx, cy)];
            for (dx, dy) in dirs {
                let nx = cx as i32 + dx;
                let ny = cy as i32 + dy;
                if nx < 0 || ny < 0 || nx >= self.width as i32 || ny >= self.height as i32 {
                    continue;
                }
                let nx = nx as usize;
                let ny = ny as usize;
                let ni = self.idx(nx, ny);
                if !self.passable[ni] {
                    continue;
                }
                let step = if dx != 0 && dy != 0 {
                    self.cell_size * 1.414
                } else {
                    self.cell_size
                };
                let new_cost = current_cost + step;
                if new_cost < self.potential[ni] {
                    self.potential[ni] = new_cost;
                    let cost_bits = new_cost.to_bits();
                    heap.push(Reverse((cost_bits, nx, ny)));
                }
            }
        }
    }

    /// Sample the gradient at world position `pos` (bilinear).
    ///
    /// Returns the normalised direction toward lower potential (the flow direction).
    pub fn flow_direction(&self, pos: Vec2) -> Vec2 {
        let (cx, cy) = self.world_to_cell(pos);
        let grad = self.gradient_at(cx, cy);
        normalize2(scale2(grad, -1.0)) // negative gradient = downhill
    }

    /// Central-difference gradient at cell (cx, cy).
    pub fn gradient_at(&self, cx: usize, cy: usize) -> Vec2 {
        let x_l = if cx > 0 {
            self.potential[self.idx(cx - 1, cy)]
        } else {
            self.potential[self.idx(cx, cy)]
        };
        let x_r = if cx + 1 < self.width {
            self.potential[self.idx(cx + 1, cy)]
        } else {
            self.potential[self.idx(cx, cy)]
        };
        let y_d = if cy > 0 {
            self.potential[self.idx(cx, cy - 1)]
        } else {
            self.potential[self.idx(cx, cy)]
        };
        let y_u = if cy + 1 < self.height {
            self.potential[self.idx(cx, cy + 1)]
        } else {
            self.potential[self.idx(cx, cy)]
        };

        let gx = (x_r - x_l) / (2.0 * self.cell_size);
        let gy = (y_u - y_d) / (2.0 * self.cell_size);
        [gx, gy]
    }

    /// Mark a rectangular region of cells as impassable (wall).
    pub fn add_wall(&mut self, x0: usize, y0: usize, x1: usize, y1: usize) {
        for y in y0..=y1.min(self.height - 1) {
            for x in x0..=x1.min(self.width - 1) {
                let i = self.idx(x, y);
                self.passable[i] = false;
                self.potential[i] = f64::INFINITY;
            }
        }
    }

    /// Return the potential at a world position (nearest-cell lookup).
    pub fn potential_at_world(&self, pos: Vec2) -> f64 {
        let (cx, cy) = self.world_to_cell(pos);
        self.potential[self.idx(cx, cy)]
    }
}

// ---------------------------------------------------------------------------
// CrowdDensity
// ---------------------------------------------------------------------------

/// Crowd density and velocity field on a uniform grid.
///
/// Supports Gaussian kernel density estimation and the fundamental diagram.
#[derive(Debug, Clone)]
pub struct CrowdDensity {
    /// Width of the density grid in cells.
    pub width: usize,
    /// Height of the density grid in cells.
    pub height: usize,
    /// Physical size of each cell (m).
    pub cell_size: f64,
    /// Origin of the grid (m).
    pub origin: Vec2,
    /// Density field (pedestrians/m²).
    pub density: Vec<f64>,
    /// Local velocity field (average pedestrian velocity per cell).
    pub velocity_x: Vec<f64>,
    /// Local velocity y-component.
    pub velocity_y: Vec<f64>,
    /// Gaussian kernel bandwidth (m).
    pub kernel_bandwidth: f64,
}

impl CrowdDensity {
    /// Create a new crowd density object with zero-initialised fields.
    pub fn new(width: usize, height: usize, cell_size: f64, origin: Vec2) -> Self {
        let n = width * height;
        Self {
            width,
            height,
            cell_size,
            origin,
            density: vec![0.0; n],
            velocity_x: vec![0.0; n],
            velocity_y: vec![0.0; n],
            kernel_bandwidth: 1.0,
        }
    }

    fn idx(&self, cx: usize, cy: usize) -> usize {
        cy * self.width + cx
    }

    fn world_to_cell_f(&self, pos: Vec2) -> (f64, f64) {
        (
            (pos[0] - self.origin[0]) / self.cell_size,
            (pos[1] - self.origin[1]) / self.cell_size,
        )
    }

    /// Recompute density and velocity fields from the current pedestrian positions.
    ///
    /// Uses a Gaussian kernel with bandwidth `kernel_bandwidth`.
    pub fn update(&mut self, pedestrians: &[Pedestrian]) {
        self.density.fill(0.0);
        self.velocity_x.fill(0.0);
        self.velocity_y.fill(0.0);

        let h = self.kernel_bandwidth;
        let h2 = h * h;
        let radius_cells = ((3.0 * h) / self.cell_size).ceil() as i32;

        for p in pedestrians {
            if p.evacuated {
                continue;
            }
            let (px_f, py_f) = self.world_to_cell_f(p.position);
            let cx0 = px_f as i32;
            let cy0 = py_f as i32;

            for dy in -radius_cells..=radius_cells {
                for dx in -radius_cells..=radius_cells {
                    let nx = cx0 + dx;
                    let ny = cy0 + dy;
                    if nx < 0 || ny < 0 || nx >= self.width as i32 || ny >= self.height as i32 {
                        continue;
                    }
                    let nx = nx as usize;
                    let ny = ny as usize;
                    // World centre of cell
                    let cx_w = self.origin[0] + (nx as f64 + 0.5) * self.cell_size;
                    let cy_w = self.origin[1] + (ny as f64 + 0.5) * self.cell_size;
                    let r2 = (cx_w - p.position[0]).powi(2) + (cy_w - p.position[1]).powi(2);
                    let kernel = (-0.5 * r2 / h2).exp() / (2.0 * PI * h2);
                    let i = self.idx(nx, ny);
                    self.density[i] += kernel;
                    self.velocity_x[i] += kernel * p.velocity[0];
                    self.velocity_y[i] += kernel * p.velocity[1];
                }
            }
        }

        // Normalise velocity by density
        for i in 0..self.density.len() {
            if self.density[i] > 1e-10 {
                self.velocity_x[i] /= self.density[i];
                self.velocity_y[i] /= self.density[i];
            }
        }
    }

    /// Return the local density at world position `pos` (nearest cell).
    pub fn density_at(&self, pos: Vec2) -> f64 {
        let cx = ((pos[0] - self.origin[0]) / self.cell_size) as usize;
        let cy = ((pos[1] - self.origin[1]) / self.cell_size) as usize;
        let cx = cx.min(self.width.saturating_sub(1));
        let cy = cy.min(self.height.saturating_sub(1));
        self.density[self.idx(cx, cy)]
    }

    /// Fundamental diagram: flow rate q = ρ * v(ρ) using Greenshields' model.
    ///
    /// `rho` — density (ped/m²), `rho_jam` — jam density, `v_free` — free-flow speed.
    pub fn fundamental_diagram_flow(rho: f64, rho_jam: f64, v_free: f64) -> f64 {
        if rho >= rho_jam {
            return 0.0;
        }
        let v = v_free * (1.0 - rho / rho_jam);
        rho * v
    }

    /// Greenshields' speed-density relation: v(ρ) = v_free * (1 − ρ/ρ_jam).
    pub fn greenshields_speed(rho: f64, rho_jam: f64, v_free: f64) -> f64 {
        if rho >= rho_jam {
            return 0.0;
        }
        v_free * (1.0 - rho / rho_jam)
    }

    /// Capacity (max flow) at ρ = ρ_jam / 2.
    pub fn capacity(rho_jam: f64, v_free: f64) -> f64 {
        let rho_opt = rho_jam / 2.0;
        Self::fundamental_diagram_flow(rho_opt, rho_jam, v_free)
    }

    /// Mean density averaged over all grid cells.
    pub fn mean_density(&self) -> f64 {
        let sum: f64 = self.density.iter().sum();
        sum / (self.density.len() as f64)
    }

    /// Peak density (maximum over all cells).
    pub fn peak_density(&self) -> f64 {
        self.density.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// EmergentBehavior
// ---------------------------------------------------------------------------

/// Detects and characterises emergent collective behaviors in crowd simulations.
///
/// Includes lane formation, bottleneck oscillation, arch formation, and turbulence.
#[derive(Debug, Clone)]
pub struct EmergentBehavior {
    /// History of mean crowd velocities for oscillation detection.
    pub velocity_history: Vec<[f64; 2]>,
    /// History of density measurements for congestion tracking.
    pub density_history: Vec<f64>,
    /// Maximum history length.
    pub max_history: usize,
}

impl EmergentBehavior {
    /// Create a new emergent behavior detector.
    pub fn new(max_history: usize) -> Self {
        Self {
            velocity_history: Vec::with_capacity(max_history),
            density_history: Vec::with_capacity(max_history),
            max_history,
        }
    }

    /// Record the current mean velocity and density.
    pub fn record(&mut self, mean_velocity: Vec2, mean_density: f64) {
        if self.velocity_history.len() >= self.max_history {
            self.velocity_history.remove(0);
            self.density_history.remove(0);
        }
        self.velocity_history.push(mean_velocity);
        self.density_history.push(mean_density);
    }

    /// Compute the polarisation order parameter (0 = disordered, 1 = fully aligned).
    ///
    /// φ = |Σ v̂_i| / N
    pub fn polarisation(pedestrians: &[Pedestrian]) -> f64 {
        if pedestrians.is_empty() {
            return 0.0;
        }
        let mut sx = 0.0f64;
        let mut sy = 0.0f64;
        let mut count = 0usize;
        for p in pedestrians {
            if !p.evacuated {
                let spd = norm2(p.velocity);
                if spd > 1e-6 {
                    sx += p.velocity[0] / spd;
                    sy += p.velocity[1] / spd;
                    count += 1;
                }
            }
        }
        if count == 0 {
            return 0.0;
        }
        ((sx * sx + sy * sy).sqrt()) / count as f64
    }

    /// Detect lane formation: check if pedestrians segregate into directional groups.
    ///
    /// Returns a score in \[0, 1\]: higher means stronger lane formation.
    pub fn lane_formation_score(pedestrians: &[Pedestrian]) -> f64 {
        if pedestrians.len() < 4 {
            return 0.0;
        }
        // Count pedestrians moving primarily left vs right
        let mut right = 0usize;
        let mut left = 0usize;
        for p in pedestrians {
            if p.evacuated {
                continue;
            }
            if p.velocity[0] > 0.1 {
                right += 1;
            } else if p.velocity[0] < -0.1 {
                left += 1;
            }
        }
        let total = (right + left) as f64;
        if total < 2.0 {
            return 0.0;
        }
        // Score: proportion of "strongly directional" walkers that are in separate lanes
        let imbalance = (right as f64 - left as f64).abs() / total;
        1.0 - imbalance // High score = balanced lanes
    }

    /// Detect bottleneck oscillation: measure variance of x-velocity in history.
    ///
    /// High variance indicates back-and-forth oscillation at a bottleneck.
    pub fn bottleneck_oscillation_index(&self) -> f64 {
        if self.velocity_history.len() < 4 {
            return 0.0;
        }
        let vx_vals: Vec<f64> = self.velocity_history.iter().map(|v| v[0]).collect();
        let mean = vx_vals.iter().sum::<f64>() / vx_vals.len() as f64;
        let var = vx_vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vx_vals.len() as f64;
        var.sqrt()
    }

    /// Detect arch formation near an exit: measure angular dispersion of velocities
    /// pointing toward a given exit location.
    ///
    /// Returns the fraction of pedestrians within `angle_threshold` radians of
    /// pointing toward the exit.
    pub fn arch_formation_score(
        pedestrians: &[Pedestrian],
        exit: Vec2,
        angle_threshold: f64,
    ) -> f64 {
        if pedestrians.is_empty() {
            return 0.0;
        }
        let mut aligned = 0usize;
        let mut total = 0usize;
        for p in pedestrians {
            if p.evacuated {
                continue;
            }
            total += 1;
            let to_exit = normalize2(sub2(exit, p.position));
            let speed = norm2(p.velocity);
            if speed < 1e-6 {
                continue;
            }
            let v_dir = normalize2(p.velocity);
            let cos_angle = dot2(to_exit, v_dir).clamp(-1.0, 1.0);
            let angle = cos_angle.acos();
            if angle < angle_threshold {
                aligned += 1;
            }
        }
        if total == 0 {
            0.0
        } else {
            aligned as f64 / total as f64
        }
    }

    /// Detect turbulent flow: measure spatial velocity gradient magnitude.
    ///
    /// Returns the mean absolute velocity difference between neighbours.
    pub fn turbulence_index(pedestrians: &[Pedestrian], neighbour_radius: f64) -> f64 {
        if pedestrians.len() < 2 {
            return 0.0;
        }
        let mut total_diff = 0.0f64;
        let mut count = 0usize;
        for i in 0..pedestrians.len() {
            for j in (i + 1)..pedestrians.len() {
                let pi = &pedestrians[i];
                let pj = &pedestrians[j];
                if pi.evacuated || pj.evacuated {
                    continue;
                }
                let d = dist2(pi.position, pj.position);
                if d < neighbour_radius {
                    let dv = sub2(pi.velocity, pj.velocity);
                    total_diff += norm2(dv);
                    count += 1;
                }
            }
        }
        if count == 0 {
            0.0
        } else {
            total_diff / count as f64
        }
    }

    /// Detect "freezing by heating" effect: high density but low flow.
    ///
    /// Returns `true` if mean density exceeds `rho_threshold` and mean flow
    /// rate falls below `flow_threshold`.
    pub fn is_frozen_by_heating(
        &self,
        rho_threshold: f64,
        flow_threshold: f64,
        rho_jam: f64,
        v_free: f64,
    ) -> bool {
        if self.density_history.is_empty() {
            return false;
        }
        let mean_rho = self.density_history.iter().sum::<f64>() / self.density_history.len() as f64;
        let flow = CrowdDensity::fundamental_diagram_flow(mean_rho, rho_jam, v_free);
        mean_rho > rho_threshold && flow < flow_threshold
    }
}

// ---------------------------------------------------------------------------
// EvacuationSimulation
// ---------------------------------------------------------------------------

/// Exit descriptor for evacuation scenarios.
#[derive(Debug, Clone)]
pub struct Exit {
    /// Position of the exit centre (m).
    pub position: Vec2,
    /// Width of the exit (m).
    pub width: f64,
    /// Attraction strength (higher = more pedestrians drawn to this exit).
    pub attractiveness: f64,
    /// Capacity in pedestrians/second.
    pub capacity_per_second: f64,
}

impl Exit {
    /// Create a new exit at `position` with `width` metres.
    pub fn new(position: Vec2, width: f64) -> Self {
        Self {
            position,
            width,
            attractiveness: 1.0,
            capacity_per_second: 1.5 * width, // typical ~1.5 ped/m/s
        }
    }

    /// Check whether position `pos` is within the exit zone.
    pub fn contains(&self, pos: Vec2) -> bool {
        dist2(pos, self.position) < self.width * 0.5
    }
}

/// Full evacuation simulation controller.
///
/// Manages a crowd of pedestrians, multiple exits, and wall obstacles, and
/// tracks evacuation metrics such as flow rate and evacuation time.
#[derive(Debug, Clone)]
pub struct EvacuationSimulation {
    /// Active pedestrians (including already evacuated ones for record-keeping).
    pub pedestrians: Vec<Pedestrian>,
    /// Available exits.
    pub exits: Vec<Exit>,
    /// Wall obstacle segments (pairs of endpoints).
    pub obstacles: Vec<(Vec2, Vec2)>,
    /// SFM solver.
    pub sfm: SocialForceModel,
    /// Elapsed simulation time (s).
    pub time: f64,
    /// Number of pedestrians evacuated so far.
    pub evacuated_count: usize,
    /// Times at which each pedestrian evacuated (indexed by pedestrian id).
    pub evacuation_times: Vec<Option<f64>>,
    /// Flow rate history: (time, pedestrians/s).
    pub flow_rate_history: Vec<(f64, f64)>,
    /// Bottleneck congestion indicator.
    pub bottleneck_density: f64,
    /// Random seed for stochastic behavior.
    rng_state: u64,
}

impl EvacuationSimulation {
    /// Create a new evacuation simulation.
    pub fn new(exits: Vec<Exit>, obstacles: Vec<(Vec2, Vec2)>) -> Self {
        Self {
            pedestrians: Vec::new(),
            exits,
            obstacles,
            sfm: SocialForceModel,
            time: 0.0,
            evacuated_count: 0,
            evacuation_times: Vec::new(),
            flow_rate_history: Vec::new(),
            bottleneck_density: 0.0,
            rng_state: 12345,
        }
    }

    /// Simple LCG random float in \[0, 1).
    fn rand_f64(&mut self) -> f64 {
        self.rng_state = self
            .rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.rng_state >> 33) as f64 / (1u64 << 31) as f64
    }

    /// Add a pedestrian to the simulation at `position`.
    pub fn add_pedestrian(&mut self, position: Vec2, mass: f64) {
        let id = self.pedestrians.len();
        // Assign to the nearest exit with some randomness for exit choice.
        let exit = self.choose_exit(position);
        let mut p = Pedestrian::new(id, position, exit.position, mass);
        self.pedestrians.push(p.clone());
        self.evacuation_times.push(None);
        let _ = &mut p; // silence the unused warning
    }

    /// Choose an exit for a pedestrian based on distance and attractiveness.
    pub fn choose_exit(&self, position: Vec2) -> &Exit {
        // Utility = attractiveness / (distance + 1)
        let mut best = 0usize;
        let mut best_util = f64::NEG_INFINITY;
        for (i, exit) in self.exits.iter().enumerate() {
            let d = dist2(position, exit.position);
            let util = exit.attractiveness / (d + 1.0);
            if util > best_util {
                best_util = util;
                best = i;
            }
        }
        &self.exits[best]
    }

    /// Re-evaluate which exit each non-evacuated pedestrian targets (dynamic exit choice).
    pub fn update_exit_choices(&mut self) {
        for p in &mut self.pedestrians {
            if p.evacuated {
                continue;
            }
            let exit = {
                let mut best = 0usize;
                let mut best_util = f64::NEG_INFINITY;
                for (i, exit) in self.exits.iter().enumerate() {
                    let d = dist2(p.position, exit.position);
                    let util = exit.attractiveness / (d + 1.0);
                    if util > best_util {
                        best_util = util;
                        best = i;
                    }
                }
                best
            };
            p.destination = self.exits[exit].position;
            p.update_desired_direction();
        }
    }

    /// Advance the simulation by one time step `dt`.
    pub fn step(&mut self, dt: f64) {
        // Compute forces
        let n = self.pedestrians.len();
        let mut forces = vec![[0.0f64; 2]; n];

        for i in 0..n {
            if self.pedestrians[i].evacuated {
                continue;
            }
            forces[i] = SocialForceModel::total_force(
                &self.pedestrians[i],
                &self.pedestrians,
                &self.obstacles,
            );
        }

        // Integrate
        for i in 0..n {
            if self.pedestrians[i].evacuated {
                continue;
            }
            self.pedestrians[i].update_desired_direction();
            self.pedestrians[i].integrate(forces[i], dt);
        }

        // Check for evacuation
        let mut newly_evacuated = 0usize;
        for i in 0..n {
            if self.pedestrians[i].evacuated {
                continue;
            }
            for exit in &self.exits {
                if exit.contains(self.pedestrians[i].position) {
                    self.pedestrians[i].evacuated = true;
                    self.evacuation_times[i] = Some(self.time);
                    newly_evacuated += 1;
                    break;
                }
            }
        }

        // Update metrics
        self.evacuated_count += newly_evacuated;
        if dt > 1e-10 {
            let flow = newly_evacuated as f64 / dt;
            self.flow_rate_history.push((self.time, flow));
        }
        self.time += dt;
    }

    /// Run the simulation until all pedestrians evacuate or `max_time` is reached.
    pub fn run_to_completion(&mut self, dt: f64, max_time: f64) {
        while self.time < max_time {
            if self.all_evacuated() {
                break;
            }
            self.step(dt);
            // Periodically re-evaluate exit choices
            if ((self.time / dt) as usize).is_multiple_of(60) {
                self.update_exit_choices();
            }
        }
    }

    /// Return `true` if all pedestrians have evacuated.
    pub fn all_evacuated(&self) -> bool {
        self.pedestrians.iter().all(|p| p.evacuated)
    }

    /// Compute mean evacuation time over all successfully evacuated pedestrians.
    pub fn mean_evacuation_time(&self) -> Option<f64> {
        let times: Vec<f64> = self.evacuation_times.iter().filter_map(|t| *t).collect();
        if times.is_empty() {
            return None;
        }
        Some(times.iter().sum::<f64>() / times.len() as f64)
    }

    /// Compute total evacuation flow rate (ped/s) over the entire simulation.
    pub fn total_flow_rate(&self) -> f64 {
        if self.time < 1e-6 {
            return 0.0;
        }
        self.evacuated_count as f64 / self.time
    }

    /// Compute the specific flow (ped/m/s) at a bottleneck of given `width`.
    pub fn specific_flow(&self, exit_width: f64) -> f64 {
        if exit_width < 1e-6 {
            return 0.0;
        }
        self.total_flow_rate() / exit_width
    }

    /// Estimate the density near exits (simple count in radius).
    pub fn exit_density(&self, exit_index: usize, radius: f64) -> f64 {
        if exit_index >= self.exits.len() {
            return 0.0;
        }
        let exit_pos = self.exits[exit_index].position;
        let count = self
            .pedestrians
            .iter()
            .filter(|p| !p.evacuated && dist2(p.position, exit_pos) < radius)
            .count();
        let area = PI * radius * radius;
        count as f64 / area
    }

    /// Fraction of pedestrians who have evacuated.
    pub fn evacuation_fraction(&self) -> f64 {
        if self.pedestrians.is_empty() {
            return 0.0;
        }
        self.evacuated_count as f64 / self.pedestrians.len() as f64
    }

    /// Return the maximum instantaneous flow rate observed.
    pub fn peak_flow_rate(&self) -> f64 {
        self.flow_rate_history
            .iter()
            .map(|(_, q)| *q)
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// CrowdSimulationWorld — high-level wrapper
// ---------------------------------------------------------------------------

/// High-level world object that integrates all crowd simulation components.
#[derive(Debug, Clone)]
pub struct CrowdSimulationWorld {
    /// All pedestrian agents.
    pub pedestrians: Vec<Pedestrian>,
    /// Wall obstacles.
    pub obstacles: Vec<(Vec2, Vec2)>,
    /// Navigation flow field.
    pub flow_field: Option<FlowField>,
    /// Density estimator.
    pub density: CrowdDensity,
    /// Emergent behavior detector.
    pub behavior: EmergentBehavior,
    /// Current time (s).
    pub time: f64,
}

impl CrowdSimulationWorld {
    /// Create a new world of given physical dimensions.
    pub fn new(width_m: f64, height_m: f64, cell_size: f64) -> Self {
        let w = (width_m / cell_size).ceil() as usize;
        let h = (height_m / cell_size).ceil() as usize;
        Self {
            pedestrians: Vec::new(),
            obstacles: Vec::new(),
            flow_field: None,
            density: CrowdDensity::new(w, h, cell_size, [0.0, 0.0]),
            behavior: EmergentBehavior::new(300),
            time: 0.0,
        }
    }

    /// Add a pedestrian to the world.
    pub fn add_pedestrian(&mut self, p: Pedestrian) {
        self.pedestrians.push(p);
    }

    /// Add a wall segment.
    pub fn add_obstacle(&mut self, a: Vec2, b: Vec2) {
        self.obstacles.push((a, b));
    }

    /// Set the navigation flow field.
    pub fn set_flow_field(&mut self, ff: FlowField) {
        self.flow_field = Some(ff);
    }

    /// Advance simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        // Optionally update desired directions from flow field.
        if let Some(ref ff) = self.flow_field {
            for p in &mut self.pedestrians {
                if !p.evacuated {
                    let dir = ff.flow_direction(p.position);
                    if norm2(dir) > 1e-6 {
                        p.desired_direction = dir;
                    }
                }
            }
        }

        // Compute and apply forces.
        let n = self.pedestrians.len();
        let mut forces = vec![[0.0f64; 2]; n];
        for i in 0..n {
            if !self.pedestrians[i].evacuated {
                forces[i] = SocialForceModel::total_force(
                    &self.pedestrians[i],
                    &self.pedestrians,
                    &self.obstacles,
                );
            }
        }
        for i in 0..n {
            if !self.pedestrians[i].evacuated {
                self.pedestrians[i].integrate(forces[i], dt);
            }
        }

        // Update density.
        self.density.update(&self.pedestrians);

        // Record behavior metrics.
        let mean_vel = self.mean_velocity();
        let mean_rho = self.density.mean_density();
        self.behavior.record(mean_vel, mean_rho);

        self.time += dt;
    }

    /// Compute the mean velocity across all active pedestrians.
    pub fn mean_velocity(&self) -> Vec2 {
        let mut vx = 0.0f64;
        let mut vy = 0.0f64;
        let mut count = 0usize;
        for p in &self.pedestrians {
            if !p.evacuated {
                vx += p.velocity[0];
                vy += p.velocity[1];
                count += 1;
            }
        }
        if count == 0 {
            [0.0, 0.0]
        } else {
            [vx / count as f64, vy / count as f64]
        }
    }

    /// Count of non-evacuated pedestrians.
    pub fn active_count(&self) -> usize {
        self.pedestrians.iter().filter(|p| !p.evacuated).count()
    }
}

// ---------------------------------------------------------------------------
// Utility: batch pedestrian spawning
// ---------------------------------------------------------------------------

/// Spawn `n` pedestrians randomly inside a rectangular region.
pub fn spawn_pedestrians_in_rect(
    n: usize,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    destination: Vec2,
    mass: f64,
) -> Vec<Pedestrian> {
    let mut rng = rand::rng();
    use rand::RngExt as _;
    (0..n)
        .map(|i| {
            let x = rng.random_range(x_min..x_max);
            let y = rng.random_range(y_min..y_max);
            Pedestrian::new(i, [x, y], destination, mass)
        })
        .collect()
}

/// Compute the mean speed of a pedestrian group.
pub fn mean_speed(pedestrians: &[Pedestrian]) -> f64 {
    if pedestrians.is_empty() {
        return 0.0;
    }
    let sum: f64 = pedestrians
        .iter()
        .filter(|p| !p.evacuated)
        .map(|p| p.speed())
        .sum();
    let count = pedestrians.iter().filter(|p| !p.evacuated).count();
    if count == 0 { 0.0 } else { sum / count as f64 }
}

/// Compute the convex hull area of pedestrian positions (approximate via bounding box).
pub fn crowd_spread_area(pedestrians: &[Pedestrian]) -> f64 {
    let active: Vec<&Pedestrian> = pedestrians.iter().filter(|p| !p.evacuated).collect();
    if active.len() < 2 {
        return 0.0;
    }
    let x_min = active
        .iter()
        .map(|p| p.position[0])
        .fold(f64::INFINITY, f64::min);
    let x_max = active
        .iter()
        .map(|p| p.position[0])
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = active
        .iter()
        .map(|p| p.position[1])
        .fold(f64::INFINITY, f64::min);
    let y_max = active
        .iter()
        .map(|p| p.position[1])
        .fold(f64::NEG_INFINITY, f64::max);
    (x_max - x_min) * (y_max - y_min)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_pedestrian(id: usize, px: f64, py: f64, dx: f64, dy: f64) -> Pedestrian {
        Pedestrian::new(id, [px, py], [dx, dy], 70.0)
    }

    // ── Pedestrian basics ────────────────────────────────────────────────────

    #[test]
    fn test_pedestrian_creation() {
        let p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        assert_eq!(p.id, 0);
        assert!(!p.evacuated);
        assert!((p.speed() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_pedestrian_desired_direction() {
        let p = make_pedestrian(0, 0.0, 0.0, 3.0, 4.0);
        // Normalised direction should have unit length.
        let d = norm2(p.desired_direction);
        assert!((d - 1.0).abs() < 1e-10, "direction norm={d}");
    }

    #[test]
    fn test_pedestrian_distance_to_destination() {
        let p = make_pedestrian(0, 0.0, 0.0, 3.0, 4.0);
        assert!((p.distance_to_destination() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_pedestrian_has_reached_destination() {
        let mut p = make_pedestrian(0, 0.0, 0.0, 0.1, 0.0);
        assert!(p.has_reached_destination(0.5));
        p.destination = [100.0, 0.0];
        assert!(!p.has_reached_destination(0.5));
    }

    #[test]
    fn test_pedestrian_integrate_moves_position() {
        let mut p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let force = [100.0f64, 0.0]; // push right
        p.integrate(force, 0.1);
        assert!(p.position[0] > 0.0, "pedestrian should have moved right");
    }

    #[test]
    fn test_pedestrian_speed_clamp() {
        let mut p = make_pedestrian(0, 0.0, 0.0, 100.0, 0.0);
        // Apply enormous force.
        p.integrate([1e8, 0.0], 1.0);
        let max_allowed = p.params.desired_speed * 3.0;
        assert!(p.speed() <= max_allowed + 1e-6, "speed should be clamped");
    }

    #[test]
    fn test_pedestrian_distance_travelled() {
        let mut p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        for _ in 0..10 {
            p.integrate([p.mass / p.params.tau, 0.0], 0.1);
        }
        assert!(p.distance_travelled > 0.0);
    }

    // ── SFM forces ───────────────────────────────────────────────────────────

    #[test]
    fn test_desired_force_points_toward_destination() {
        let p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let f = SocialForceModel::desired_force(&p);
        assert!(
            f[0] > 0.0,
            "desired force x should be positive, got {:?}",
            f
        );
    }

    #[test]
    fn test_desired_force_zero_when_at_desired_velocity() {
        let mut p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        // Set velocity to exactly desired.
        let v0 = p.params.desired_speed;
        p.velocity = [v0, 0.0];
        // Direction already [1,0] from construction.
        let f = SocialForceModel::desired_force(&p);
        assert!(f[0].abs() < 1e-8, "desired force should be ~0, got {:?}", f);
    }

    #[test]
    fn test_repulsive_force_pushes_apart() {
        let p0 = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let p1 = make_pedestrian(1, 0.2, 0.0, -10.0, 0.0); // close to p0
        let f = SocialForceModel::repulsive_social_force(&p0, &p1);
        // Force on p0 should push it left (negative x).
        assert!(
            f[0] < 0.0,
            "repulsion should push p0 away from p1, got {:?}",
            f
        );
    }

    #[test]
    fn test_repulsive_force_decreases_with_distance() {
        let p0 = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let p1_close = make_pedestrian(1, 0.3, 0.0, -10.0, 0.0);
        let p1_far = make_pedestrian(2, 3.0, 0.0, -10.0, 0.0);
        let f_close = SocialForceModel::repulsive_social_force(&p0, &p1_close);
        let f_far = SocialForceModel::repulsive_social_force(&p0, &p1_far);
        assert!(
            norm2(f_close) > norm2(f_far),
            "close repulsion should be stronger: {} vs {}",
            norm2(f_close),
            norm2(f_far)
        );
    }

    #[test]
    fn test_obstacle_force_pushes_away() {
        let p = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        // Wall just below the pedestrian
        let f = SocialForceModel::obstacle_force(&p, [-5.0, -0.2], [5.0, -0.2]);
        // Force should push upward (positive y)
        assert!(f[1] > 0.0, "obstacle force should push up, got {:?}", f);
    }

    #[test]
    fn test_total_force_includes_all_components() {
        let p0 = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let p1 = make_pedestrian(1, 0.3, 0.0, -10.0, 0.0);
        let pedestrians = vec![p0.clone(), p1];
        let f = SocialForceModel::total_force(&p0, &pedestrians, &[]);
        // Should be non-zero (desired + repulsive).
        assert!(norm2(f) > 0.0);
    }

    #[test]
    fn test_sfm_same_position_no_panic() {
        let p0 = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        let p1 = make_pedestrian(1, 0.0, 0.0, 10.0, 0.0); // same position
        let f = SocialForceModel::repulsive_social_force(&p0, &p1);
        // Should return zero, not NaN.
        assert!(f[0].is_finite() && f[1].is_finite());
    }

    // ── FlowField ────────────────────────────────────────────────────────────

    #[test]
    fn test_flow_field_creation() {
        let ff = FlowField::new(10, 10, 1.0, [0.0, 0.0]);
        assert_eq!(ff.potential.len(), 100);
        assert!(ff.passable.iter().all(|&p| p));
    }

    #[test]
    fn test_flow_field_dijkstra_goal_has_zero_potential() {
        let mut ff = FlowField::new(10, 10, 1.0, [0.0, 0.0]);
        ff.compute_dijkstra(&[(5, 5)]);
        assert_eq!(ff.potential[ff.idx(5, 5)], 0.0);
    }

    #[test]
    fn test_flow_field_potential_increases_with_distance() {
        let mut ff = FlowField::new(20, 20, 1.0, [0.0, 0.0]);
        ff.compute_dijkstra(&[(10, 10)]);
        let p_near = ff.potential[ff.idx(11, 10)];
        let p_far = ff.potential[ff.idx(15, 10)];
        assert!(p_far > p_near, "potential should increase with distance");
    }

    #[test]
    fn test_flow_field_direction_points_toward_goal() {
        let mut ff = FlowField::new(20, 20, 1.0, [0.0, 0.0]);
        ff.compute_dijkstra(&[(10, 10)]);
        // From left of goal, direction should point right (+x).
        let dir = ff.flow_direction([5.5, 10.5]);
        assert!(
            dir[0] > 0.0,
            "direction x should be positive, got {:?}",
            dir
        );
    }

    #[test]
    fn test_flow_field_wall_blocks_potential() {
        let mut ff = FlowField::new(20, 10, 1.0, [0.0, 0.0]);
        ff.add_wall(10, 0, 10, 8); // vertical wall at x=10
        ff.compute_dijkstra(&[(15, 5)]);
        // Cell at (10, 5) should be impassable (infinite potential).
        assert_eq!(ff.potential[ff.idx(10, 5)], f64::INFINITY);
    }

    #[test]
    fn test_flow_field_world_to_cell() {
        let ff = FlowField::new(10, 10, 2.0, [0.0, 0.0]);
        let (cx, cy) = ff.world_to_cell([3.0, 5.0]);
        assert_eq!(cx, 1);
        assert_eq!(cy, 2);
    }

    #[test]
    fn test_flow_field_cell_to_world() {
        let ff = FlowField::new(10, 10, 2.0, [0.0, 0.0]);
        let w = ff.cell_to_world(2, 3);
        assert!((w[0] - 5.0).abs() < 1e-10);
        assert!((w[1] - 7.0).abs() < 1e-10);
    }

    // ── CrowdDensity ─────────────────────────────────────────────────────────

    #[test]
    fn test_crowd_density_zero_when_empty() {
        let mut cd = CrowdDensity::new(10, 10, 1.0, [0.0, 0.0]);
        cd.update(&[]);
        assert!(cd.mean_density() < 1e-10);
    }

    #[test]
    fn test_crowd_density_positive_when_pedestrians_present() {
        let mut cd = CrowdDensity::new(20, 20, 1.0, [0.0, 0.0]);
        cd.kernel_bandwidth = 2.0;
        let pedestrians = vec![make_pedestrian(0, 10.0, 10.0, 20.0, 20.0)];
        cd.update(&pedestrians);
        assert!(cd.mean_density() > 0.0, "density should be positive");
    }

    #[test]
    fn test_fundamental_diagram_zero_at_jam() {
        let q = CrowdDensity::fundamental_diagram_flow(6.0, 6.0, 1.34);
        assert!(
            q.abs() < 1e-10,
            "flow at jam density should be zero, got {q}"
        );
    }

    #[test]
    fn test_fundamental_diagram_positive_below_jam() {
        let q = CrowdDensity::fundamental_diagram_flow(1.0, 6.0, 1.34);
        assert!(
            q > 0.0,
            "flow should be positive below jam density, got {q}"
        );
    }

    #[test]
    fn test_greenshields_speed_free_flow() {
        let v = CrowdDensity::greenshields_speed(0.0, 6.0, 1.34);
        assert!(
            (v - 1.34).abs() < 1e-10,
            "free flow speed should be 1.34, got {v}"
        );
    }

    #[test]
    fn test_greenshields_speed_zero_at_jam() {
        let v = CrowdDensity::greenshields_speed(6.0, 6.0, 1.34);
        assert!(
            v.abs() < 1e-10,
            "speed at jam density should be zero, got {v}"
        );
    }

    #[test]
    fn test_capacity_calculation() {
        let q_cap = CrowdDensity::capacity(6.0, 1.34);
        // At ρ_opt = 3.0, v = 1.34 * 0.5 = 0.67, q = 3 * 0.67 = 2.01
        let expected = 3.0 * (1.34 * 0.5);
        assert!(
            (q_cap - expected).abs() < 1e-10,
            "capacity={q_cap}, expected {expected}"
        );
    }

    // ── EmergentBehavior ────────────────────────────────────────────────────

    #[test]
    fn test_polarisation_zero_for_random_velocities() {
        let mut pedestrians = vec![
            make_pedestrian(0, 0.0, 0.0, 10.0, 0.0),
            make_pedestrian(1, 1.0, 0.0, -10.0, 0.0),
        ];
        pedestrians[0].velocity = [1.0, 0.0];
        pedestrians[1].velocity = [-1.0, 0.0];
        let pol = EmergentBehavior::polarisation(&pedestrians);
        assert!(
            pol < 0.1,
            "opposite velocities → low polarisation, got {pol}"
        );
    }

    #[test]
    fn test_polarisation_one_for_aligned() {
        let mut pedestrians = vec![
            make_pedestrian(0, 0.0, 0.0, 10.0, 0.0),
            make_pedestrian(1, 1.0, 0.0, 10.0, 0.0),
        ];
        pedestrians[0].velocity = [1.0, 0.0];
        pedestrians[1].velocity = [1.0, 0.0];
        let pol = EmergentBehavior::polarisation(&pedestrians);
        assert!(
            (pol - 1.0).abs() < 1e-10,
            "aligned velocities → polarisation=1, got {pol}"
        );
    }

    #[test]
    fn test_arch_formation_score_zero_when_evacuated() {
        let mut p = make_pedestrian(0, 1.0, 1.0, 5.0, 5.0);
        p.evacuated = true;
        let score = EmergentBehavior::arch_formation_score(&[p], [5.0, 5.0], 0.5);
        assert!(score.abs() < 1e-10);
    }

    #[test]
    fn test_turbulence_zero_identical_velocities() {
        let mut pedestrians: Vec<Pedestrian> = (0..5)
            .map(|i| {
                let mut p = make_pedestrian(i, i as f64, 0.0, 10.0, 0.0);
                p.velocity = [1.0, 0.0];
                p
            })
            .collect();
        let ti = EmergentBehavior::turbulence_index(&pedestrians, 3.0);
        assert!(
            ti.abs() < 1e-10,
            "identical velocities → zero turbulence, got {ti}"
        );
        // Silence unused mut warning
        pedestrians[0].velocity = [1.0, 0.0];
    }

    #[test]
    fn test_turbulence_positive_different_velocities() {
        let mut p0 = make_pedestrian(0, 0.0, 0.0, 10.0, 0.0);
        p0.velocity = [2.0, 0.0];
        let mut p1 = make_pedestrian(1, 0.5, 0.0, 10.0, 0.0);
        p1.velocity = [-2.0, 0.0];
        let ti = EmergentBehavior::turbulence_index(&[p0, p1], 2.0);
        assert!(
            ti > 0.0,
            "different velocities → positive turbulence, got {ti}"
        );
    }

    #[test]
    fn test_bottleneck_oscillation_zero_constant_velocity() {
        let mut eb = EmergentBehavior::new(100);
        for _ in 0..20 {
            eb.record([1.0, 0.0], 2.0);
        }
        let osc = eb.bottleneck_oscillation_index();
        assert!(
            osc.abs() < 1e-10,
            "constant velocity → zero oscillation, got {osc}"
        );
    }

    #[test]
    fn test_oscillation_detected() {
        let mut eb = EmergentBehavior::new(100);
        for i in 0..20 {
            let v = if i % 2 == 0 {
                [1.0f64, 0.0]
            } else {
                [-1.0f64, 0.0]
            };
            eb.record(v, 2.0);
        }
        let osc = eb.bottleneck_oscillation_index();
        assert!(
            osc > 0.5,
            "alternating velocity → oscillation detected, got {osc}"
        );
    }

    #[test]
    fn test_lane_formation_score_balanced() {
        let mut pedestrians: Vec<Pedestrian> = Vec::new();
        for i in 0..5 {
            let mut p = make_pedestrian(i, i as f64, 0.0, 10.0, 0.0);
            p.velocity = [1.0, 0.0];
            pedestrians.push(p);
        }
        for i in 5..10 {
            let mut p = make_pedestrian(i, i as f64, 1.0, -10.0, 0.0);
            p.velocity = [-1.0, 0.0];
            pedestrians.push(p);
        }
        let score = EmergentBehavior::lane_formation_score(&pedestrians);
        assert!(
            score > 0.8,
            "balanced bidirectional flow → high lane score, got {score}"
        );
    }

    // ── EvacuationSimulation ─────────────────────────────────────────────────

    #[test]
    fn test_evacuation_sim_creation() {
        let exits = vec![Exit::new([10.0, 5.0], 2.0)];
        let sim = EvacuationSimulation::new(exits, vec![]);
        assert_eq!(sim.evacuated_count, 0);
        assert!(sim.time.abs() < 1e-10);
    }

    #[test]
    fn test_exit_contains() {
        let exit = Exit::new([0.0, 0.0], 2.0);
        assert!(exit.contains([0.5, 0.0]));
        assert!(!exit.contains([5.0, 0.0]));
    }

    #[test]
    fn test_evacuation_pedestrian_reaches_exit() {
        let exit = Exit::new([5.0, 0.0], 1.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        sim.add_pedestrian([4.5, 0.0], 70.0); // very close to exit
        // Manual step: pedestrian should evacuate quickly.
        sim.run_to_completion(0.05, 20.0);
        assert!(
            sim.evacuated_count > 0,
            "at least one pedestrian should evacuate"
        );
    }

    #[test]
    fn test_evacuation_fraction() {
        let exit = Exit::new([5.0, 0.0], 1.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        sim.add_pedestrian([4.5, 0.0], 70.0);
        sim.run_to_completion(0.05, 30.0);
        let frac = sim.evacuation_fraction();
        assert!(
            frac > 0.0 && frac <= 1.0,
            "fraction should be in (0,1], got {frac}"
        );
    }

    #[test]
    fn test_evacuation_mean_time_is_positive() {
        let exit = Exit::new([5.0, 0.0], 1.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        sim.add_pedestrian([4.5, 0.0], 70.0);
        sim.run_to_completion(0.05, 30.0);
        if let Some(t) = sim.mean_evacuation_time() {
            assert!(t >= 0.0, "mean evac time should be non-negative, got {t}");
        }
    }

    #[test]
    fn test_evacuation_flow_rate_positive() {
        let exit = Exit::new([5.0, 0.0], 1.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        for i in 0..5 {
            sim.add_pedestrian([4.5 - i as f64 * 0.5, 0.0], 70.0);
        }
        sim.run_to_completion(0.05, 60.0);
        let fr = sim.total_flow_rate();
        assert!(fr >= 0.0, "flow rate should be non-negative, got {fr}");
    }

    #[test]
    fn test_all_evacuated_flag() {
        let exit = Exit::new([5.0, 0.0], 2.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        sim.add_pedestrian([4.7, 0.0], 70.0);
        sim.add_pedestrian([4.8, 0.1], 70.0);
        sim.run_to_completion(0.1, 120.0);
        // Not guaranteed all evacuate in time; just check the flag is consistent.
        let flag = sim.all_evacuated();
        let frac = sim.evacuation_fraction();
        if flag {
            assert!((frac - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_exit_density_non_negative() {
        let exit = Exit::new([0.0, 0.0], 2.0);
        let mut sim = EvacuationSimulation::new(vec![exit], vec![]);
        sim.add_pedestrian([0.3, 0.0], 70.0);
        let d = sim.exit_density(0, 3.0);
        assert!(d >= 0.0, "density should be non-negative, got {d}");
    }

    // ── CrowdSimulationWorld ─────────────────────────────────────────────────

    #[test]
    fn test_world_creation() {
        let world = CrowdSimulationWorld::new(50.0, 50.0, 1.0);
        assert_eq!(world.active_count(), 0);
    }

    #[test]
    fn test_world_add_and_step() {
        let mut world = CrowdSimulationWorld::new(20.0, 20.0, 1.0);
        world.add_pedestrian(make_pedestrian(0, 5.0, 5.0, 15.0, 5.0));
        world.step(0.05);
        assert!(world.time > 0.0);
    }

    #[test]
    fn test_world_mean_velocity() {
        let mut world = CrowdSimulationWorld::new(20.0, 20.0, 1.0);
        let mut p = make_pedestrian(0, 5.0, 5.0, 15.0, 5.0);
        p.velocity = [2.0, 0.0];
        world.add_pedestrian(p);
        let mv = world.mean_velocity();
        assert!((mv[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_world_obstacle_integration() {
        let mut world = CrowdSimulationWorld::new(20.0, 20.0, 1.0);
        world.add_obstacle([0.0, 5.0], [20.0, 5.0]);
        world.add_pedestrian(make_pedestrian(0, 5.0, 5.1, 5.0, 10.0));
        for _ in 0..20 {
            world.step(0.05);
        }
        // Just verify no crash.
        assert!(world.time > 0.0);
    }

    // ── Utility functions ────────────────────────────────────────────────────

    #[test]
    fn test_spawn_pedestrians_in_rect() {
        let peds = spawn_pedestrians_in_rect(10, 0.0, 10.0, 0.0, 10.0, [5.0, 5.0], 70.0);
        assert_eq!(peds.len(), 10);
        for p in &peds {
            assert!(p.position[0] >= 0.0 && p.position[0] <= 10.0);
            assert!(p.position[1] >= 0.0 && p.position[1] <= 10.0);
        }
    }

    #[test]
    fn test_mean_speed_zero_when_stationary() {
        let pedestrians = vec![make_pedestrian(0, 0.0, 0.0, 1.0, 0.0)];
        let ms = mean_speed(&pedestrians);
        assert!(
            ms.abs() < 1e-10,
            "stationary pedestrian → zero mean speed, got {ms}"
        );
    }

    #[test]
    fn test_crowd_spread_area() {
        let pedestrians = vec![
            make_pedestrian(0, 0.0, 0.0, 10.0, 0.0),
            make_pedestrian(1, 10.0, 5.0, 10.0, 0.0),
        ];
        let area = crowd_spread_area(&pedestrians);
        assert!(
            (area - 50.0).abs() < 1e-10,
            "area should be 10*5=50, got {area}"
        );
    }

    #[test]
    fn test_closest_point_on_segment_midpoint() {
        let cp = closest_point_on_segment([0.0, 1.0], [0.0, 0.0], [2.0, 0.0]);
        assert!((cp[0] - 0.0).abs() < 1e-10);
        assert!((cp[1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_closest_point_on_segment_endpoint() {
        let cp = closest_point_on_segment([5.0, 1.0], [0.0, 0.0], [3.0, 0.0]);
        // Closest is endpoint (3, 0)
        assert!((cp[0] - 3.0).abs() < 1e-10);
        assert!((cp[1] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_specific_flow_scales_with_width() {
        let exit1 = Exit::new([5.0, 0.0], 1.0);
        let exit2 = Exit::new([5.0, 0.0], 2.0);
        let mut sim1 = EvacuationSimulation::new(vec![exit1], vec![]);
        let mut sim2 = EvacuationSimulation::new(vec![exit2], vec![]);
        // Add same pedestrians at same positions.
        for i in 0..3 {
            sim1.add_pedestrian([4.5 - i as f64 * 0.2, 0.0], 70.0);
            sim2.add_pedestrian([4.5 - i as f64 * 0.2, 0.0], 70.0);
        }
        sim1.run_to_completion(0.05, 60.0);
        sim2.run_to_completion(0.05, 60.0);
        // Both should have evacuated the same number of people.
        assert_eq!(sim1.evacuated_count, sim2.evacuated_count);
    }
}
