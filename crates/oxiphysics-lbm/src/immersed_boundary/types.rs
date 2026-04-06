//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use std::f64::consts::PI;

/// A Lagrangian marker point on an immersed boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct IbMarker {
    /// Position \[x, y\]
    pub position: [f64; 2],
    /// Velocity \[ux, uy\]
    pub velocity: [f64; 2],
    /// Body force \[fx, fy\]
    pub force: [f64; 2],
}
impl IbMarker {
    /// Create a new IbMarker at the given position with zero velocity and force.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            position: [x, y],
            velocity: [0.0, 0.0],
            force: [0.0, 0.0],
        }
    }
}
/// Statistics about the immersed boundary configuration.
#[derive(Clone, Debug)]
pub struct IbmStats {
    /// Net force on all markers \[fx, fy\]
    pub total_force: [f64; 2],
    /// Maximum displacement of any marker from its rest position
    pub max_marker_displacement: f64,
}
impl IbmStats {
    /// Compute stats from current marker positions and forces.
    pub fn compute(markers: &[IbMarker], rest_positions: &[[f64; 2]]) -> Self {
        let mut total_force = [0.0_f64; 2];
        let mut max_disp = 0.0_f64;
        for (i, m) in markers.iter().enumerate() {
            total_force[0] += m.force[0];
            total_force[1] += m.force[1];
            if i < rest_positions.len() {
                let dx = m.position[0] - rest_positions[i][0];
                let dy = m.position[1] - rest_positions[i][1];
                let disp = (dx * dx + dy * dy).sqrt();
                if disp > max_disp {
                    max_disp = disp;
                }
            }
        }
        Self {
            total_force,
            max_marker_displacement: max_disp,
        }
    }
}
/// A collection of immersed boundaries for multi-body coupling.
#[derive(Clone, Debug)]
pub struct MultiIbSystem {
    /// List of immersed bodies, each represented by a set of markers.
    pub bodies: Vec<Vec<IbMarker>>,
    /// Rest positions for each body.
    pub rest_positions: Vec<Vec<[f64; 2]>>,
    /// Stiffness for each body.
    pub stiffnesses: Vec<f64>,
}
impl MultiIbSystem {
    /// Create an empty multi-IB system.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            rest_positions: Vec::new(),
            stiffnesses: Vec::new(),
        }
    }
}

impl Default for MultiIbSystem {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiIbSystem {
    /// Add a body (circle) to the system.
    pub fn add_circle(&mut self, circle: &IbCircle) {
        let markers = circle.to_markers();
        let rest: Vec<[f64; 2]> = markers.iter().map(|m| m.position).collect();
        self.bodies.push(markers);
        self.rest_positions.push(rest);
        self.stiffnesses.push(circle.stiffness);
    }
    /// Compute elastic forces for all bodies.
    pub fn compute_all_forces(&mut self) {
        for (i, body) in self.bodies.iter_mut().enumerate() {
            let stiffness = self.stiffnesses[i];
            let rest = &self.rest_positions[i];
            compute_elastic_forces(body, rest, stiffness);
        }
    }
    /// Spread forces from all bodies onto a single grid.
    pub fn spread_all_forces(&self, nx: usize, ny: usize, dx: f64) -> Vec<[f64; 2]> {
        let mut total_force = vec![[0.0_f64; 2]; nx * ny];
        for body in &self.bodies {
            let body_force = spread_force(body, nx, ny, dx);
            for k in 0..(nx * ny) {
                total_force[k][0] += body_force[k][0];
                total_force[k][1] += body_force[k][1];
            }
        }
        total_force
    }
    /// Total number of markers across all bodies.
    pub fn total_markers(&self) -> usize {
        self.bodies.iter().map(|b| b.len()).sum()
    }
}
/// Virtual boundary / feedback forcing for IBM.
///
/// The feedback forcing drives the marker velocity toward a target velocity:
///
/// F_k = α_f * (U_target - U_marker) + β_f * ∫(U_target - U_marker) dt
///
/// This struct accumulates the integral term.
#[allow(dead_code)]
pub struct FeedbackForcing {
    /// Proportional gain α_f.
    pub alpha_f: f64,
    /// Integral gain β_f.
    pub beta_f: f64,
    /// Accumulated velocity error integral per marker \[fx, fy\].
    pub integral: Vec<[f64; 2]>,
}
impl FeedbackForcing {
    /// Create a new FeedbackForcing with given gains and `n` markers.
    #[allow(dead_code)]
    pub fn new(alpha_f: f64, beta_f: f64, n: usize) -> Self {
        Self {
            alpha_f,
            beta_f,
            integral: vec![[0.0; 2]; n],
        }
    }
    /// Compute forcing for each marker given interpolated and target velocities.
    ///
    /// Also updates the integral accumulator.
    #[allow(dead_code)]
    pub fn compute(
        &mut self,
        u_interp: &[[f64; 2]],
        u_target: &[[f64; 2]],
        dt: f64,
    ) -> Vec<[f64; 2]> {
        let n = u_interp.len().min(u_target.len()).min(self.integral.len());
        let mut forces = vec![[0.0_f64; 2]; n];
        for k in 0..n {
            for d in 0..2 {
                let err = u_target[k][d] - u_interp[k][d];
                self.integral[k][d] += err * dt;
                forces[k][d] = self.alpha_f * err + self.beta_f * self.integral[k][d];
            }
        }
        forces
    }
    /// Reset integral accumulators to zero.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        for s in self.integral.iter_mut() {
            *s = [0.0; 2];
        }
    }
}
/// A 3D Lagrangian marker for immersed boundary method.
#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq)]
pub struct IbMarker3D {
    /// Position \[x, y, z\]
    pub position: [f64; 3],
    /// Velocity \[ux, uy, uz\]
    pub velocity: [f64; 3],
    /// Body force \[fx, fy, fz\]
    pub force: [f64; 3],
}
impl IbMarker3D {
    /// Create a new 3D IBM marker at position with zero velocity and force.
    #[allow(dead_code)]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            position: [x, y, z],
            velocity: [0.0; 3],
            force: [0.0; 3],
        }
    }
    /// Distance from this marker to another marker.
    #[allow(dead_code)]
    pub fn distance_to(&self, other: &IbMarker3D) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Configuration parameters for the IBM-LBM coupling.
#[derive(Clone, Debug)]
pub struct IbmConfig {
    /// Grid spacing (lattice unit)
    pub delta_x: f64,
    /// Time step
    pub delta_t: f64,
    /// Fluid density
    pub rho_f: f64,
}
impl IbmConfig {
    /// Create a new IbmConfig.
    pub fn new(delta_x: f64, delta_t: f64, rho_f: f64) -> Self {
        Self {
            delta_x,
            delta_t,
            rho_f,
        }
    }
}
/// A rigid immersed body with surface marker points.
///
/// Tracks body centroid position, velocity, orientation, angular velocity,
/// and accumulates forces and torques from the IBM coupling.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct IbmBody {
    /// Centroid position \[x, y\]
    pub centroid: [f64; 2],
    /// Centroid velocity \[vx, vy\]
    pub velocity: [f64; 2],
    /// Angular velocity (rad/s)
    pub angular_velocity: f64,
    /// Orientation angle (rad)
    pub angle: f64,
    /// Surface marker positions relative to centroid (body frame)
    pub surface_offsets: Vec<[f64; 2]>,
    /// Surface marker forces \[fx, fy\]
    pub surface_forces: Vec<[f64; 2]>,
    /// Mass of the rigid body
    pub mass: f64,
    /// Moment of inertia
    pub moment_of_inertia: f64,
}
impl IbmBody {
    /// Create a new IbmBody centred at `(cx, cy)` with given mass and inertia.
    #[allow(dead_code)]
    pub fn new(cx: f64, cy: f64, mass: f64, moment: f64) -> Self {
        Self {
            centroid: [cx, cy],
            velocity: [0.0, 0.0],
            angular_velocity: 0.0,
            angle: 0.0,
            surface_offsets: Vec::new(),
            surface_forces: Vec::new(),
            mass,
            moment_of_inertia: moment,
        }
    }
    /// Add a surface marker at offset `(dx, dy)` from centroid in body frame.
    #[allow(dead_code)]
    pub fn add_surface_point(&mut self, dx: f64, dy: f64) {
        self.surface_offsets.push([dx, dy]);
        self.surface_forces.push([0.0, 0.0]);
    }
    /// Get world-frame positions of all surface markers.
    #[allow(dead_code)]
    pub fn surface_positions(&self) -> Vec<[f64; 2]> {
        let cos_a = self.angle.cos();
        let sin_a = self.angle.sin();
        self.surface_offsets
            .iter()
            .map(|&[ox, oy]| {
                let rx = cos_a * ox - sin_a * oy;
                let ry = sin_a * ox + cos_a * oy;
                [self.centroid[0] + rx, self.centroid[1] + ry]
            })
            .collect()
    }
    /// Get surface marker velocities (body + rotational contribution).
    #[allow(dead_code)]
    pub fn surface_velocities(&self) -> Vec<[f64; 2]> {
        let cos_a = self.angle.cos();
        let sin_a = self.angle.sin();
        self.surface_offsets
            .iter()
            .map(|&[ox, oy]| {
                let rx = cos_a * ox - sin_a * oy;
                let ry = sin_a * ox + cos_a * oy;
                [
                    self.velocity[0] - self.angular_velocity * ry,
                    self.velocity[1] + self.angular_velocity * rx,
                ]
            })
            .collect()
    }
    /// Compute net force on the body from surface forces.
    #[allow(dead_code)]
    pub fn net_force(&self) -> [f64; 2] {
        self.surface_forces
            .iter()
            .fold([0.0, 0.0], |acc, &f| [acc[0] + f[0], acc[1] + f[1]])
    }
    /// Compute net torque about centroid from surface forces.
    #[allow(dead_code)]
    pub fn net_torque(&self) -> f64 {
        let cos_a = self.angle.cos();
        let sin_a = self.angle.sin();
        self.surface_offsets
            .iter()
            .zip(self.surface_forces.iter())
            .map(|(&[ox, oy], &[fx, fy])| {
                let rx = cos_a * ox - sin_a * oy;
                let ry = sin_a * ox + cos_a * oy;
                rx * fy - ry * fx
            })
            .sum()
    }
    /// Advance body dynamics using Newton's 2nd law (Euler step).
    ///
    /// F = m * a  →  v += (F/m) * dt
    /// τ = I * α  →  ω += (τ/I) * dt
    /// x += v * dt,  θ += ω * dt
    #[allow(dead_code)]
    pub fn advance(&mut self, dt: f64) {
        let [fx, fy] = self.net_force();
        self.velocity[0] += fx / self.mass * dt;
        self.velocity[1] += fy / self.mass * dt;
        let torque = self.net_torque();
        self.angular_velocity += torque / self.moment_of_inertia * dt;
        self.centroid[0] += self.velocity[0] * dt;
        self.centroid[1] += self.velocity[1] * dt;
        self.angle += self.angular_velocity * dt;
    }
    /// Create a circular IbmBody with `n` surface points.
    #[allow(dead_code)]
    pub fn circle(cx: f64, cy: f64, r: f64, n: usize, rho_body: f64) -> Self {
        let area = PI * r * r;
        let mass = rho_body * area;
        let moment = 0.5 * mass * r * r;
        let mut body = Self::new(cx, cy, mass, moment);
        for i in 0..n {
            let theta = 2.0 * PI * i as f64 / n as f64;
            body.add_surface_point(r * theta.cos(), r * theta.sin());
        }
        body
    }
}
/// A circular elastic immersed boundary defined by center, radius, and stiffness.
#[derive(Clone, Debug)]
pub struct IbCircle {
    /// Center position \[cx, cy\]
    pub center: [f64; 2],
    /// Radius
    pub radius: f64,
    /// Number of Lagrangian markers
    pub n_markers: usize,
    /// Spring stiffness
    pub stiffness: f64,
}
impl IbCircle {
    /// Create a new IbCircle.
    pub fn new(cx: f64, cy: f64, radius: f64, n_markers: usize, stiffness: f64) -> Self {
        Self {
            center: [cx, cy],
            radius,
            n_markers,
            stiffness,
        }
    }
    /// Distribute `n_markers` evenly on the circle and return them.
    pub fn to_markers(&self) -> Vec<IbMarker> {
        (0..self.n_markers)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / self.n_markers as f64;
                IbMarker::new(
                    self.center[0] + self.radius * theta.cos(),
                    self.center[1] + self.radius * theta.sin(),
                )
            })
            .collect()
    }
}
/// Goldstein feedback forcing structure.
///
/// Implements the virtual boundary method with proportional-integral control:
///   F_k = α * (U_B - U_f) + β * ∫(U_B - U_f) dt
#[allow(dead_code)]
pub struct GoldsteinForcing {
    /// Proportional gain α (large, ~10⁴ to 10⁶)
    pub alpha: f64,
    /// Integral gain β
    pub beta: f64,
    /// Accumulated velocity error (integral term)
    pub error_integral: Vec<[f64; 2]>,
}
impl GoldsteinForcing {
    /// Create a new Goldstein forcing controller for `n` markers.
    #[allow(dead_code)]
    pub fn new(alpha: f64, beta: f64, n: usize) -> Self {
        Self {
            alpha,
            beta,
            error_integral: vec![[0.0; 2]; n],
        }
    }
    /// Compute feedback forces given desired and interpolated velocities.
    ///
    /// Updates integral accumulators and returns forces per marker.
    #[allow(dead_code)]
    pub fn compute_forces(
        &mut self,
        u_desired: &[[f64; 2]],
        u_interp: &[[f64; 2]],
        dt: f64,
    ) -> Vec<[f64; 2]> {
        let n = u_desired
            .len()
            .min(u_interp.len())
            .min(self.error_integral.len());
        let mut forces = vec![[0.0_f64; 2]; n];
        for k in 0..n {
            for d in 0..2 {
                let err = u_desired[k][d] - u_interp[k][d];
                self.error_integral[k][d] += err * dt;
                forces[k][d] = self.alpha * err + self.beta * self.error_integral[k][d];
            }
        }
        forces
    }
    /// Reset error integral to zero.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        for s in self.error_integral.iter_mut() {
            *s = [0.0; 2];
        }
    }
}
