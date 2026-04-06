//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// A simple inlet buffer that maintains layers of particles at a fixed x-plane.
///
/// Particles are stored with the prescribed inflow velocity and density.
/// Layers are added at `inflow_pos` along the x-axis and shift each time step.
#[derive(Debug, Clone)]
pub struct InletBuffer {
    /// Stored buffer particles (positions only; all carry `velocity` and `density`).
    pub particles: Vec<[f64; 3]>,
    /// Prescribed inflow velocity (m/s).
    pub velocity: [f64; 3],
    /// Prescribed density (kg/m³).
    pub density: f64,
    /// Number of layers to maintain.
    pub n_layers: usize,
    /// Particle spacing Δx (m).
    pub dx: f64,
}
impl InletBuffer {
    /// Create a new (empty) inlet buffer.
    pub fn new(velocity: [f64; 3], density: f64, n_layers: usize, dx: f64) -> Self {
        Self {
            particles: Vec::new(),
            velocity,
            density,
            n_layers,
            dx,
        }
    }
    /// Add a new layer of particles at `inflow_pos` along the x-axis.
    ///
    /// A layer consists of `n_layers` particles spaced `dx` apart in the
    /// y-direction, centred at y = 0.
    pub fn generate_layer(&mut self, inflow_pos: f64) {
        for k in 0..self.n_layers {
            let y = (k as f64 - (self.n_layers as f64 - 1.0) * 0.5) * self.dx;
            self.particles.push([inflow_pos, y, 0.0]);
        }
    }
    /// Advance all buffer particles by `velocity * dt`.
    pub fn shift_buffer(&mut self, dt: f64) {
        let v = self.velocity;
        for p in &mut self.particles {
            p[0] += v[0] * dt;
            p[1] += v[1] * dt;
            p[2] += v[2] * dt;
        }
    }
}
/// Convective (advective) outflow boundary condition.
///
/// ∂u/∂t + c_conv ∂u/∂x = 0  where c_conv = mean advection speed.
///
/// Discretised as: u^{n+1}_out = u^n_out - c_conv Δt/Δx (u^n_out - u^n_{in})
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConvectiveOutflowBC {
    /// Convective speed (m/s).
    pub c_conv: f64,
    /// Spacing between outlet and interior layer (m).
    pub dx: f64,
    /// Time step Δt (s).
    pub dt: f64,
    /// Previous outlet velocities.
    pub vel_outlet_prev: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl ConvectiveOutflowBC {
    /// Create a new convective outflow BC.
    pub fn new(c_conv: f64, dx: f64, dt: f64, n_particles: usize) -> Self {
        Self {
            c_conv,
            dx,
            dt,
            vel_outlet_prev: vec![[0.0; 3]; n_particles],
        }
    }
    /// Apply the convective BC.
    ///
    /// - `vel_outlet` : current outlet velocities.
    /// - `vel_interior` : velocities one layer inside.
    ///   Returns the updated outlet velocities.
    pub fn apply(&mut self, vel_outlet: &[[f64; 3]], vel_interior: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = vel_outlet
            .len()
            .min(vel_interior.len())
            .min(self.vel_outlet_prev.len());
        let courant = (self.c_conv * self.dt / self.dx).min(1.0);
        let mut result = vel_outlet.to_vec();
        for i in 0..n {
            for c in 0..3 {
                result[i][c] = vel_outlet[i][c] - courant * (vel_outlet[i][c] - vel_interior[i][c]);
            }
        }
        self.vel_outlet_prev = vel_outlet.to_vec();
        result
    }
    /// Resize internal buffer to match a new particle count.
    pub fn resize(&mut self, n: usize) {
        self.vel_outlet_prev.resize(n, [0.0; 3]);
    }
}
/// PID controller for open boundary flow-rate regulation.
///
/// Adjusts the prescribed inlet velocity to match a target volumetric
/// flow rate Q_target (m³/s).  The error is e = Q_actual - Q_target.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FlowRateController {
    /// Target volumetric flow rate Q_target (m³/s).
    pub q_target: f64,
    /// Proportional gain K_p.
    pub kp: f64,
    /// Integral gain K_i.
    pub ki: f64,
    /// Derivative gain K_d.
    pub kd: f64,
    /// Accumulated integral of the error (m³).
    pub integral: f64,
    /// Previous error for derivative (m³/s).
    pub prev_error: f64,
    /// Maximum magnitude of velocity correction (m/s).
    pub max_correction: f64,
}
#[allow(dead_code)]
impl FlowRateController {
    /// Create a new flow-rate PID controller.
    pub fn new(q_target: f64, kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            q_target,
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            max_correction: 5.0,
        }
    }
    /// Compute the velocity correction given the measured flow rate.
    ///
    /// Returns the signed velocity increment Δu to add to the inlet velocity.
    pub fn update(&mut self, q_actual: f64, dt: f64) -> f64 {
        let error = q_actual - self.q_target;
        self.integral += error * dt;
        let derivative = (error - self.prev_error) / dt.max(1e-30);
        self.prev_error = error;
        let correction = -(self.kp * error + self.ki * self.integral + self.kd * derivative);
        correction.clamp(-self.max_correction, self.max_correction)
    }
    /// Reset the controller state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
    /// Estimate volumetric flow rate from particle velocities crossing a plane.
    ///
    /// Q = Σ_p A_p (v_p · n̂)  where A_p = m_p / (ρ_p · L) is the particle
    /// cross-sectional area and L is the layer thickness.
    pub fn estimate_flow_rate(
        velocities: &[[f64; 3]],
        normal: [f64; 3],
        masses: &[f64],
        densities: &[f64],
        layer_thickness: f64,
    ) -> f64 {
        let mut q = 0.0_f64;
        for i in 0..velocities.len() {
            let v_dot_n = velocities[i][0] * normal[0]
                + velocities[i][1] * normal[1]
                + velocities[i][2] * normal[2];
            let area_p = if densities[i] > 1e-30 {
                masses[i] / (densities[i] * layer_thickness.max(1e-30))
            } else {
                0.0
            };
            q += v_dot_n * area_p;
        }
        q
    }
}
/// Riemann solver state for a buffer zone particle at an open boundary.
///
/// Uses the Roe-averaged characteristics to prescribe non-reflecting
/// inflow conditions (Riemann invariant method).
#[derive(Debug, Clone)]
pub struct RiemannBufferParticle {
    /// Position.
    pub position: [f64; 3],
    /// Prescribed velocity (characteristic from exterior state).
    pub velocity: [f64; 3],
    /// Prescribed pressure (from exterior state).
    pub pressure: f64,
    /// Prescribed density.
    pub density: f64,
    /// Local sound speed at this particle.
    pub sound_speed: f64,
}
impl RiemannBufferParticle {
    /// Create a Riemann buffer particle with a given exterior state.
    #[allow(dead_code)]
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        pressure: f64,
        density: f64,
        sound_speed: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            pressure,
            density,
            sound_speed,
        }
    }
    /// Riemann invariant in the normal direction (incoming wave).
    ///
    /// `W+ = v_n + 2*c / (gamma-1)` for the incoming Riemann invariant.
    #[allow(dead_code)]
    pub fn riemann_invariant_plus(&self, normal: [f64; 3], gamma: f64) -> f64 {
        let v_n = self.velocity[0] * normal[0]
            + self.velocity[1] * normal[1]
            + self.velocity[2] * normal[2];
        v_n + 2.0 * self.sound_speed / (gamma - 1.0)
    }
    /// Riemann invariant for the outgoing wave.
    ///
    /// `W- = v_n - 2*c / (gamma-1)`
    #[allow(dead_code)]
    pub fn riemann_invariant_minus(&self, normal: [f64; 3], gamma: f64) -> f64 {
        let v_n = self.velocity[0] * normal[0]
            + self.velocity[1] * normal[1]
            + self.velocity[2] * normal[2];
        v_n - 2.0 * self.sound_speed / (gamma - 1.0)
    }
}
/// A ghost particle residing in the inflow reservoir.
#[derive(Debug, Clone)]
pub struct BufferParticle {
    /// Current position of the ghost particle.
    pub position: [f64; 3],
    /// Prescribed velocity of this buffer particle.
    pub velocity: [f64; 3],
}
/// Applies ramped damping forces near an outlet and flags particles for removal.
#[derive(Debug, Clone)]
pub struct OutflowDamper {
    /// The outlet boundary zone.
    pub zone: BoundaryZone,
    /// Length of the damping ramp (m).
    pub damping_length: f64,
    /// Maximum damping coefficient α_max (1/s).
    pub alpha_max: f64,
}
impl OutflowDamper {
    /// Create a new outflow damper.
    pub fn new(zone: BoundaryZone, damping_length: f64, alpha_max: f64) -> Self {
        Self {
            zone,
            damping_length,
            alpha_max,
        }
    }
    /// Compute the ramped damping force for a particle at `pos` with velocity `vel`.
    ///
    /// The damping coefficient ramps from 0 (at distance `damping_length` before
    /// the boundary) to `alpha_max` (at the boundary).
    /// `f_damp = -α(x) · vel`
    pub fn apply_damping_force(&self, pos: [f64; 3], vel: [f64; 3]) -> [f64; 3] {
        let sd = self.zone.signed_distance(pos);
        if sd < 0.0 || sd > self.damping_length {
            return [0.0, 0.0, 0.0];
        }
        let t = 1.0 - sd / self.damping_length;
        let alpha = self.alpha_max * t;
        [-alpha * vel[0], -alpha * vel[1], -alpha * vel[2]]
    }
    /// Compute the quadratic damping coefficient at a given position.
    ///
    /// Uses quadratic ramping: α = α_max * t² for smoother transition.
    pub fn quadratic_damping_coefficient(&self, pos: [f64; 3]) -> f64 {
        let sd = self.zone.signed_distance(pos);
        if sd < 0.0 || sd > self.damping_length {
            return 0.0;
        }
        let t = 1.0 - sd / self.damping_length;
        self.alpha_max * t * t
    }
    /// Check if a particle should be removed (past the boundary).
    pub fn should_remove(&self, pos: [f64; 3]) -> bool {
        self.zone.signed_distance(pos) < 0.0
    }
}
/// A sponge (absorbing) layer that damps waves near an outlet.
///
/// The sponge layer gradually relaxes solution variables toward
/// a target (reference) state to prevent wave reflection.
#[derive(Debug, Clone)]
pub struct SpongeLayer {
    /// Zone defining the sponge region.
    pub zone: BoundaryZone,
    /// Thickness of the sponge layer (m).
    pub thickness: f64,
    /// Maximum damping rate (1/s).
    pub sigma_max: f64,
    /// Polynomial order for ramping (1 = linear, 2 = quadratic, 3 = cubic).
    pub polynomial_order: u32,
}
impl SpongeLayer {
    /// Create a new sponge layer.
    #[allow(dead_code)]
    pub fn new(zone: BoundaryZone, thickness: f64, sigma_max: f64, polynomial_order: u32) -> Self {
        Self {
            zone,
            thickness,
            sigma_max,
            polynomial_order,
        }
    }
    /// Compute the local damping coefficient at `pos`.
    ///
    /// Returns 0 outside the sponge, rising to `sigma_max` at the outlet.
    #[allow(dead_code)]
    pub fn damping_at(&self, pos: [f64; 3]) -> f64 {
        let sd = self.zone.signed_distance(pos);
        if sd < 0.0 || sd > self.thickness {
            return 0.0;
        }
        let t = 1.0 - sd / self.thickness;
        self.sigma_max * t.powi(self.polynomial_order as i32)
    }
    /// Apply the sponge forcing to a velocity field.
    ///
    /// Returns the relaxation force: `f = -sigma * (v - v_ref)`
    #[allow(dead_code)]
    pub fn apply_sponge_forcing(&self, pos: [f64; 3], vel: [f64; 3], v_ref: [f64; 3]) -> [f64; 3] {
        let sigma = self.damping_at(pos);
        [
            -sigma * (vel[0] - v_ref[0]),
            -sigma * (vel[1] - v_ref[1]),
            -sigma * (vel[2] - v_ref[2]),
        ]
    }
}
/// Turbulence injection parameters for SPH inlet boundaries.
#[derive(Debug, Clone)]
pub struct TurbulenceInjectionConfig {
    /// Turbulence intensity (ratio of RMS fluctuation to mean velocity, 0–1).
    pub intensity: f64,
    /// Turbulence length scale (m).
    pub length_scale: f64,
    /// Random seed for reproducibility.
    pub seed: u64,
}
/// Manages a reservoir of ghost (buffer) particles at an inlet.
///
/// When a buffer particle crosses the boundary it is promoted to a real
/// fluid particle.
#[derive(Debug, Clone)]
pub struct InflowBuffer {
    /// Zone describing this inlet.
    pub zone: BoundaryZone,
    /// Pre-filled ghost particles in the inlet reservoir.
    pub reservoir: Vec<BufferParticle>,
    /// Prescribed inflow velocity (m/s).
    pub inflow_velocity: [f64; 3],
    /// Nominal particle spacing (m).
    pub particle_spacing: f64,
}
impl InflowBuffer {
    /// Create a new inflow buffer.
    pub fn new(zone: BoundaryZone, inflow_velocity: [f64; 3], particle_spacing: f64) -> Self {
        Self {
            zone,
            reservoir: Vec::new(),
            inflow_velocity,
            particle_spacing,
        }
    }
    /// Fill a rectangular slab of ghost particles.
    ///
    /// The slab extends from `z_start` to `z_end` along the boundary normal,
    /// with `nx × ny` particles laid out on a regular lattice in the plane.
    /// Returns the positions of all generated particles.
    pub fn generate_layer(&self, z_start: f64, z_end: f64, nx: usize, ny: usize) -> Vec<[f64; 3]> {
        let n = self.zone.normal;
        let up = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let t1_raw = [
            n[1] * up[2] - n[2] * up[1],
            n[2] * up[0] - n[0] * up[2],
            n[0] * up[1] - n[1] * up[0],
        ];
        let t1_len = (t1_raw[0].powi(2) + t1_raw[1].powi(2) + t1_raw[2].powi(2)).sqrt();
        let t1 = if t1_len > 1e-14 {
            [t1_raw[0] / t1_len, t1_raw[1] / t1_len, t1_raw[2] / t1_len]
        } else {
            [1.0, 0.0, 0.0]
        };
        let t2 = [
            n[1] * t1[2] - n[2] * t1[1],
            n[2] * t1[0] - n[0] * t1[2],
            n[0] * t1[1] - n[1] * t1[0],
        ];
        let dx = self.particle_spacing;
        let n_layers = (((z_end - z_start) / dx).ceil() as usize).max(1);
        let mut positions = Vec::with_capacity(nx * ny * n_layers);
        let cx = self.zone.position;
        for lz in 0..n_layers {
            let z = z_start + lz as f64 * dx;
            for ix in 0..nx {
                let xi = (ix as f64 - (nx as f64 - 1.0) * 0.5) * dx;
                for iy in 0..ny {
                    let yi = (iy as f64 - (ny as f64 - 1.0) * 0.5) * dx;
                    let pos = [
                        cx[0] + xi * t1[0] + yi * t2[0] + z * n[0],
                        cx[1] + xi * t1[1] + yi * t2[1] + z * n[1],
                        cx[2] + xi * t1[2] + yi * t2[2] + z * n[2],
                    ];
                    positions.push(pos);
                }
            }
        }
        positions
    }
    /// Advance buffer particles and inject those that have crossed the boundary.
    ///
    /// Buffer particles move with `inflow_velocity * dt`.  Any buffer particle
    /// whose signed distance from the boundary plane becomes ≥ 0 (i.e. it has
    /// entered the fluid domain) is promoted: its data is appended to
    /// `positions` / `velocities`, it is removed from the reservoir, and a new
    /// buffer particle is placed one `particle_spacing` behind the boundary to
    /// refill the reservoir.
    ///
    /// Returns the indices (into the updated `positions`/`velocities` slices)
    /// of newly injected particles.
    pub fn inject_particles(
        &mut self,
        dt: f64,
        positions: &mut Vec<[f64; 3]>,
        velocities: &mut Vec<[f64; 3]>,
    ) -> Vec<usize> {
        let v = self.inflow_velocity;
        let dx = self.particle_spacing;
        let n = self.zone.normal;
        for bp in &mut self.reservoir {
            bp.position[0] += v[0] * dt;
            bp.position[1] += v[1] * dt;
            bp.position[2] += v[2] * dt;
        }
        let mut injected = Vec::new();
        let mut new_reservoir = Vec::new();
        for bp in self.reservoir.drain(..) {
            let sd = self.zone.signed_distance(bp.position);
            if sd >= 0.0 {
                let idx = positions.len();
                positions.push(bp.position);
                velocities.push(bp.velocity);
                injected.push(idx);
                let new_pos = [
                    bp.position[0] - n[0] * dx,
                    bp.position[1] - n[1] * dx,
                    bp.position[2] - n[2] * dx,
                ];
                new_reservoir.push(BufferParticle {
                    position: new_pos,
                    velocity: bp.velocity,
                });
            } else {
                new_reservoir.push(bp);
            }
        }
        self.reservoir = new_reservoir;
        injected
    }
    /// Set inflow velocity with a parabolic profile.
    pub fn set_parabolic_profile(&mut self, v_max: f64, radius: f64) {
        let center = self.zone.position;
        let axis = self.zone.normal;
        for bp in &mut self.reservoir {
            bp.velocity = parabolic_velocity_profile(center, axis, radius, v_max, bp.position);
        }
    }
}
/// Characteristic boundary conditions for SPH inlets/outlets.
///
/// Uses acoustic Riemann invariants to set boundary values that minimise
/// spurious wave reflections.
#[derive(Debug, Clone)]
pub struct CharacteristicBc {
    /// Reference velocity magnitude (m/s).
    pub u_ref: f64,
    /// Reference density (kg/m³).
    pub rho_ref: f64,
    /// Reference speed of sound (m/s).
    pub c_ref: f64,
}
impl CharacteristicBc {
    /// Create a new characteristic BC.
    pub fn new(u_ref: f64, rho_ref: f64, c_ref: f64) -> Self {
        Self {
            u_ref,
            rho_ref,
            c_ref,
        }
    }
    /// Apply inlet characteristic BC.
    ///
    /// Blends the interior state `(rho, u)` with the prescribed free-stream
    /// `u_inf` to set non-reflecting inlet values.
    ///
    /// Returns `(rho_bc, u_bc)`.
    pub fn apply_inlet(&self, rho: f64, u: [f64; 3], u_inf: [f64; 3]) -> (f64, [f64; 3]) {
        let z = self.rho_ref * self.c_ref;
        let u_n = u[0];
        let u_inf_n = u_inf[0];
        let p_int = rho * self.c_ref * self.c_ref;
        let p_inf = self.rho_ref * self.c_ref * self.c_ref;
        let j_plus = u_inf_n + p_inf / z.max(1e-14);
        let j_minus = u_n - p_int / z.max(1e-14);
        let u_n_bc = 0.5 * (j_plus + j_minus);
        let p_bc = 0.5 * z * (j_plus - j_minus);
        let rho_bc = (p_bc / (self.c_ref * self.c_ref)).max(0.0) + self.rho_ref;
        let u_bc = [u_n_bc, u[1], u[2]];
        let _ = u_inf;
        (rho_bc, u_bc)
    }
}
/// Locally One-Dimensional Inviscid (LODI) relations for inflow/outflow.
///
/// Implements the LODI system of Poinsot & Lele (1992) for subsonic inlets
/// and outlets, providing wave amplitudes for the characteristic variables.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct LodiBC {
    /// Speed of sound c (m/s).
    pub c_sound: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Relaxation coefficient σ for soft boundary (prevents reflection).
    pub sigma: f64,
    /// Reference pressure p_ref (Pa).
    pub p_ref: f64,
    /// Domain length scale L (m) for non-dimensionalisation.
    pub length: f64,
}
#[allow(dead_code)]
impl LodiBC {
    /// Create a new LODI BC.
    pub fn new(c_sound: f64, density: f64, sigma: f64, p_ref: f64, length: f64) -> Self {
        Self {
            c_sound,
            density,
            sigma,
            p_ref,
            length,
        }
    }
    /// Wave amplitude for the outgoing acoustic wave at the inlet
    /// (subsonic inlet: u > 0):
    /// L₁ = σ · (p - p_ref) · (1 - Ma²) · c / L
    pub fn l1_inlet(&self, pressure: f64, mach: f64) -> f64 {
        let k = self.sigma * self.c_sound * (1.0 - mach * mach) / self.length.max(1e-30);
        k * (pressure - self.p_ref)
    }
    /// Entropy wave amplitude:
    /// L₂ = u · (dp/dx - ρ c du/dx) ≈ 0 for LODI (simplified).
    pub fn l2_entropy(&self) -> f64 {
        0.0
    }
    /// Pressure update from LODI wave amplitudes.
    ///
    /// dp/dt = -(L₁ + L₅) / 2
    pub fn pressure_update(&self, l1: f64, l5: f64, dt: f64) -> f64 {
        -0.5 * (l1 + l5) * dt
    }
    /// Velocity update from LODI wave amplitudes.
    ///
    /// du/dt = -(L₅ - L₁) / (2 ρ c)
    pub fn velocity_update(&self, l1: f64, l5: f64, dt: f64) -> f64 {
        let rc = self.density * self.c_sound;
        if rc < 1e-30 {
            return 0.0;
        }
        -(l5 - l1) / (2.0 * rc) * dt
    }
    /// Density update: dρ/dt = -(L₂ + (L₁+L₅)/c²) ≈ -(L₁+L₅)/c²
    pub fn density_update(&self, l1: f64, l5: f64, dt: f64) -> f64 {
        let c2 = self.c_sound * self.c_sound;
        if c2 < 1e-30 {
            return 0.0;
        }
        -(l1 + l5) / c2 * dt
    }
}
/// Non-reflecting outlet using Richardson extrapolation to estimate and
/// remove reflected waves.
///
/// Stores two previous time levels of velocity at outlet particles to
/// extrapolate the outgoing characteristic.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AbsorbingOutlet {
    /// Outlet normal (points outward from the domain).
    pub normal: [f64; 3],
    /// Speed of sound c (m/s) at the outlet.
    pub c_sound: f64,
    /// Velocity at the previous time step.
    pub vel_prev: Vec<[f64; 3]>,
    /// Velocity two steps back.
    pub vel_prev2: Vec<[f64; 3]>,
    /// Grid spacing at the outlet (for the finite-difference stencil).
    pub dx: f64,
    /// Time-step size Δt.
    pub dt: f64,
}
#[allow(dead_code)]
impl AbsorbingOutlet {
    /// Create a new absorbing outlet.
    pub fn new(normal: [f64; 3], c_sound: f64, n_particles: usize, dx: f64, dt: f64) -> Self {
        Self {
            normal,
            c_sound,
            vel_prev: vec![[0.0; 3]; n_particles],
            vel_prev2: vec![[0.0; 3]; n_particles],
            dx,
            dt,
        }
    }
    /// Apply 1st-order Sommerfeld radiation condition:
    ///   ∂u/∂t + c ∂u/∂n = 0  →  u^{n+1} = u^n - c Δt/Δx (u^n - u^{n-1})
    ///
    /// Returns the corrected outlet velocities.
    pub fn apply(&mut self, vel_current: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = vel_current.len().min(self.vel_prev.len());
        let courant = (self.c_sound * self.dt / self.dx).min(1.0);
        let mut vel_new = vel_current.to_vec();
        for i in 0..n {
            for c in 0..3 {
                vel_new[i][c] =
                    vel_current[i][c] - courant * (vel_current[i][c] - self.vel_prev[i][c]);
            }
        }
        self.vel_prev2 = self.vel_prev.clone();
        self.vel_prev = vel_current.to_vec();
        vel_new
    }
    /// Resize internal buffers to match a new particle count.
    pub fn resize(&mut self, n: usize) {
        self.vel_prev.resize(n, [0.0; 3]);
        self.vel_prev2.resize(n, [0.0; 3]);
    }
}
/// An inflow boundary with turbulence injection.
#[derive(Debug, Clone)]
pub struct TurbulentInflowBoundary {
    /// Underlying inflow buffer.
    pub inflow: InflowBuffer,
    /// Turbulence configuration.
    pub turb: TurbulenceInjectionConfig,
}
impl TurbulentInflowBoundary {
    /// Create a new turbulent inflow boundary.
    #[allow(dead_code)]
    pub fn new(
        zone: BoundaryZone,
        mean_velocity: [f64; 3],
        spacing: f64,
        turb: TurbulenceInjectionConfig,
    ) -> Self {
        Self {
            inflow: InflowBuffer::new(zone, mean_velocity, spacing),
            turb,
        }
    }
    /// Update buffer particle velocities with mean + turbulent fluctuation.
    #[allow(dead_code)]
    pub fn apply_turbulent_velocities(&mut self) {
        let v_mean = self.inflow.inflow_velocity;
        let config = self.turb.clone();
        for bp in &mut self.inflow.reservoir {
            let perturb = turbulent_velocity_perturbation(bp.position, v_mean, &config);
            bp.velocity = [
                v_mean[0] + perturb[0],
                v_mean[1] + perturb[1],
                v_mean[2] + perturb[2],
            ];
        }
    }
}
/// Recycling/rescaling turbulence inlet following Lund et al. (1998).
///
/// Recycles velocity fluctuations from a downstream plane and rescales
/// them for injection at the inlet.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RecyclingTurbulenceInlet {
    /// Mean inlet velocity magnitude (m/s).
    pub u_mean: f64,
    /// Boundary layer thickness δ at the inlet (m).
    pub delta_inlet: f64,
    /// Boundary layer thickness δ at the recycling station (m).
    pub delta_recycle: f64,
    /// Rescaling factor γ = (δ_in/δ_re)^((n+1)/(2n)) for power-law profile.
    pub gamma: f64,
    /// Stored recycled velocity fluctuations for current step.
    pub fluctuations: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl RecyclingTurbulenceInlet {
    /// Create a new recycling turbulence inlet.
    pub fn new(u_mean: f64, delta_inlet: f64, delta_recycle: f64) -> Self {
        let gamma = if delta_recycle > 1e-30 {
            (delta_inlet / delta_recycle).powf(4.0 / 7.0)
        } else {
            1.0
        };
        Self {
            u_mean,
            delta_inlet,
            delta_recycle,
            gamma,
            fluctuations: Vec::new(),
        }
    }
    /// Store recycled velocity fluctuations from the recycling plane.
    ///
    /// The fluctuations are `u'(x_re, y) = u(x_re, y) - U_mean(y)`.
    pub fn store_fluctuations(&mut self, u_recycle: &[[f64; 3]], u_mean_profile: &[[f64; 3]]) {
        self.fluctuations = u_recycle
            .iter()
            .zip(u_mean_profile.iter())
            .map(|(u, m)| [u[0] - m[0], u[1] - m[1], u[2] - m[2]])
            .collect();
    }
    /// Generate inlet velocity for particle `i` at wall-normal distance `y`.
    ///
    /// u_inlet = U_mean(y) + γ · u'(x_re, y · δ_re/δ_in)
    pub fn generate_inlet_velocity(&self, idx: usize, mean_vel: [f64; 3]) -> [f64; 3] {
        if idx >= self.fluctuations.len() {
            return mean_vel;
        }
        let u_prime = self.fluctuations[idx];
        [
            mean_vel[0] + self.gamma * u_prime[0],
            mean_vel[1] + self.gamma * u_prime[1],
            mean_vel[2] + self.gamma * u_prime[2],
        ]
    }
    /// Root-mean-square of stored fluctuations (turbulence intensity measure).
    pub fn rms_fluctuation(&self) -> f64 {
        if self.fluctuations.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = self
            .fluctuations
            .iter()
            .map(|u| u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
            .sum();
        (sum_sq / self.fluctuations.len() as f64 / 3.0).sqrt()
    }
}
/// Type of an open boundary zone.
#[derive(Debug, Clone)]
pub enum BoundaryType {
    /// Inlet: prescribes velocity and density for incoming particles.
    Inflow {
        /// Prescribed inflow velocity (m/s).
        velocity: [f64; 3],
        /// Prescribed inflow density (kg/m³).
        density: f64,
    },
    /// Outlet: particles are removed when they cross the boundary.
    Outflow,
    /// Solid wall (no-slip / free-slip).
    Wall,
}
/// A buffer zone manages a collection of ghost particles near an open boundary.
///
/// Particles in the buffer zone are used to maintain proper SPH support
/// near boundaries and are gradually promoted to real particles or retired.
#[derive(Debug, Clone)]
pub struct BufferZoneManager {
    /// The boundary zone this buffer serves.
    pub zone: BoundaryZone,
    /// Buffer particles (ghost particles extending beyond the boundary).
    pub particles: Vec<BufferParticle>,
    /// Number of buffer layers behind the boundary.
    pub n_layers: usize,
    /// Particle spacing.
    pub spacing: f64,
}
impl BufferZoneManager {
    /// Create a new buffer zone manager.
    pub fn new(zone: BoundaryZone, n_layers: usize, spacing: f64) -> Self {
        Self {
            zone,
            particles: Vec::new(),
            n_layers,
            spacing,
        }
    }
    /// Fill the buffer zone with particles arranged in layers behind the boundary.
    pub fn fill_buffer(&mut self, nx: usize, ny: usize) {
        let n = self.zone.normal;
        let up = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let t1_raw = [
            n[1] * up[2] - n[2] * up[1],
            n[2] * up[0] - n[0] * up[2],
            n[0] * up[1] - n[1] * up[0],
        ];
        let t1_len = (t1_raw[0].powi(2) + t1_raw[1].powi(2) + t1_raw[2].powi(2)).sqrt();
        let t1 = if t1_len > 1e-14 {
            [t1_raw[0] / t1_len, t1_raw[1] / t1_len, t1_raw[2] / t1_len]
        } else {
            [1.0, 0.0, 0.0]
        };
        let t2 = [
            n[1] * t1[2] - n[2] * t1[1],
            n[2] * t1[0] - n[0] * t1[2],
            n[0] * t1[1] - n[1] * t1[0],
        ];
        let dx = self.spacing;
        let cx = self.zone.position;
        self.particles.clear();
        for layer in 0..self.n_layers {
            let z = -(layer as f64 + 1.0) * dx;
            for ix in 0..nx {
                let xi = (ix as f64 - (nx as f64 - 1.0) * 0.5) * dx;
                for iy in 0..ny {
                    let yi = (iy as f64 - (ny as f64 - 1.0) * 0.5) * dx;
                    let pos = [
                        cx[0] + xi * t1[0] + yi * t2[0] + z * n[0],
                        cx[1] + xi * t1[1] + yi * t2[1] + z * n[1],
                        cx[2] + xi * t1[2] + yi * t2[2] + z * n[2],
                    ];
                    self.particles.push(BufferParticle {
                        position: pos,
                        velocity: [0.0; 3],
                    });
                }
            }
        }
    }
    /// Update buffer particle velocities based on a velocity profile function.
    pub fn update_velocities<F: Fn([f64; 3]) -> [f64; 3]>(&mut self, profile: F) {
        for bp in &mut self.particles {
            bp.velocity = profile(bp.position);
        }
    }
    /// Remove buffer particles that have moved too far from the boundary.
    pub fn prune_distant_particles(&mut self, max_distance: f64) {
        self.particles.retain(|bp| {
            let sd = self.zone.signed_distance(bp.position).abs();
            sd <= max_distance
        });
    }
    /// Count the number of buffer particles.
    pub fn count(&self) -> usize {
        self.particles.len()
    }
}
/// Riemann invariants for the 1-D Euler equations.
///
/// R⁺ = u + 2c/(γ-1)  (outgoing / right-running)
/// R⁻ = u - 2c/(γ-1)  (incoming / left-running)
///
/// At an inlet the outgoing wave R⁺ is determined internally; R⁻ is prescribed.
/// At an outlet the incoming wave R⁻ is determined internally; R⁺ is prescribed.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RiemannInvariantBc {
    /// Ratio of specific heats γ (default 1.4 for air).
    pub gamma: f64,
    /// Reference speed of sound c_ref (m/s).
    pub c_ref: f64,
    /// Reference density ρ_ref (kg/m³).
    pub rho_ref: f64,
}
#[allow(dead_code)]
impl RiemannInvariantBc {
    /// Create a new Riemann-invariant BC.
    pub fn new(gamma: f64, c_ref: f64, rho_ref: f64) -> Self {
        Self {
            gamma,
            c_ref,
            rho_ref,
        }
    }
    /// Local speed of sound from density via isentropic relation:
    /// c = c_ref * (ρ/ρ_ref)^((γ-1)/2)
    pub fn sound_speed(&self, rho: f64) -> f64 {
        if rho < 1e-30 {
            return self.c_ref;
        }
        self.c_ref * (rho / self.rho_ref).powf((self.gamma - 1.0) / 2.0)
    }
    /// Outgoing Riemann invariant: R⁺ = u + 2c/(γ-1).
    pub fn r_plus(&self, u: f64, c: f64) -> f64 {
        u + 2.0 * c / (self.gamma - 1.0)
    }
    /// Incoming Riemann invariant: R⁻ = u - 2c/(γ-1).
    pub fn r_minus(&self, u: f64, c: f64) -> f64 {
        u - 2.0 * c / (self.gamma - 1.0)
    }
    /// Inlet boundary velocity from Riemann invariants.
    ///
    /// Given R⁺ from interior and prescribed R⁻_bc:
    ///   u_bc = (R⁺ + R⁻_bc) / 2
    ///   c_bc = (R⁺ - R⁻_bc) * (γ-1) / 4
    pub fn inlet_velocity(&self, r_plus_interior: f64, r_minus_bc: f64) -> (f64, f64) {
        let u_bc = (r_plus_interior + r_minus_bc) * 0.5;
        let c_bc = (r_plus_interior - r_minus_bc) * (self.gamma - 1.0) * 0.25;
        (u_bc, c_bc.max(0.0))
    }
    /// Outlet boundary velocity from Riemann invariants.
    ///
    /// Given R⁻ from interior and prescribed R⁺_bc (from far field):
    ///   u_bc = (R⁺_bc + R⁻_interior) / 2
    ///   c_bc = (R⁺_bc - R⁻_interior) * (γ-1) / 4
    pub fn outlet_velocity(&self, r_minus_interior: f64, r_plus_bc: f64) -> (f64, f64) {
        let u_bc = (r_plus_bc + r_minus_interior) * 0.5;
        let c_bc = (r_plus_bc - r_minus_interior) * (self.gamma - 1.0) * 0.25;
        (u_bc, c_bc.max(0.0))
    }
    /// Density from local sound speed via isentropic relation:
    /// ρ = ρ_ref * (c/c_ref)^(2/(γ-1))
    pub fn density_from_sound_speed(&self, c: f64) -> f64 {
        if c < 1e-30 || self.c_ref < 1e-30 {
            return self.rho_ref;
        }
        self.rho_ref * (c / self.c_ref).powf(2.0 / (self.gamma - 1.0))
    }
}
/// Outlet zone defined by a maximum x-coordinate.
///
/// Particles beyond `x_max` are flagged for removal after applying a
/// linear damping factor to their velocity.
#[derive(Debug, Clone)]
pub struct OutletZone {
    /// Maximum x-coordinate of the domain (m).
    pub x_max: f64,
    /// Damping coefficient (0 = no damping, 1 = full stop).
    pub damping: f64,
}
impl OutletZone {
    /// Create a new outlet zone.
    pub fn new(x_max: f64, damping: f64) -> Self {
        Self { x_max, damping }
    }
    /// Apply outlet damping and check for removal.
    ///
    /// Multiplies the velocity by `(1 - damping)` when the particle is near
    /// or past the outlet.  Returns `true` if the particle should be removed
    /// (i.e. it has crossed `x_max`).
    pub fn apply_outlet_damping(&self, vel: &mut [f64; 3], pos: [f64; 3]) -> bool {
        if pos[0] >= self.x_max {
            let factor = (1.0 - self.damping).max(0.0);
            vel[0] *= factor;
            vel[1] *= factor;
            vel[2] *= factor;
            return true;
        }
        false
    }
}
/// A circular/spherical zone that marks an open boundary.
#[derive(Debug, Clone)]
pub struct BoundaryZone {
    /// Center of the boundary zone.
    pub position: [f64; 3],
    /// Inward-pointing unit normal (pointing into the fluid domain).
    pub normal: [f64; 3],
    /// Radius of the zone.
    pub radius: f64,
    /// What kind of boundary this zone represents.
    pub zone_type: BoundaryType,
}
impl BoundaryZone {
    /// Create a new boundary zone with an automatically normalised normal.
    pub fn new(position: [f64; 3], normal: [f64; 3], radius: f64, zone_type: BoundaryType) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let n = if len > 1e-14 {
            [normal[0] / len, normal[1] / len, normal[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            position,
            normal: n,
            radius,
            zone_type,
        }
    }
    /// Signed distance from `pos` to the boundary plane (positive = inside fluid).
    pub fn signed_distance(&self, pos: [f64; 3]) -> f64 {
        let dx = pos[0] - self.position[0];
        let dy = pos[1] - self.position[1];
        let dz = pos[2] - self.position[2];
        dx * self.normal[0] + dy * self.normal[1] + dz * self.normal[2]
    }
    /// Check if a position is within the zone's cylinder (within radius of
    /// the boundary normal axis).
    pub fn within_radius(&self, pos: [f64; 3]) -> bool {
        let dx = pos[0] - self.position[0];
        let dy = pos[1] - self.position[1];
        let dz = pos[2] - self.position[2];
        let proj = dx * self.normal[0] + dy * self.normal[1] + dz * self.normal[2];
        let perp_x = dx - proj * self.normal[0];
        let perp_y = dy - proj * self.normal[1];
        let perp_z = dz - proj * self.normal[2];
        let perp_sq = perp_x * perp_x + perp_y * perp_y + perp_z * perp_z;
        perp_sq <= self.radius * self.radius
    }
}
/// Prescribes pressure at an outlet boundary using Tait equation of state.
///
/// Weakly-compressible SPH uses: p = B * ((rho/rho0)^gamma - 1)
#[derive(Debug, Clone)]
pub struct PressureOutletBC {
    /// Reference density (kg/m³).
    pub rho0: f64,
    /// Tait stiffness coefficient B (Pa).
    pub b_coeff: f64,
    /// Tait exponent (γ ≈ 7 for water).
    pub gamma: f64,
    /// Target outlet pressure (Pa).
    pub p_outlet: f64,
}
impl PressureOutletBC {
    /// Create a new pressure outlet BC.
    #[allow(dead_code)]
    pub fn new(rho0: f64, b_coeff: f64, gamma: f64, p_outlet: f64) -> Self {
        Self {
            rho0,
            b_coeff,
            gamma,
            p_outlet,
        }
    }
    /// Compute the target density at the outlet to achieve `p_outlet`.
    ///
    /// Inverts the Tait equation: rho = rho0 * (p/B + 1)^(1/gamma)
    #[allow(dead_code)]
    pub fn target_density(&self) -> f64 {
        let ratio = self.p_outlet / self.b_coeff + 1.0;
        if ratio < 0.0 {
            return self.rho0;
        }
        self.rho0 * ratio.powf(1.0 / self.gamma)
    }
    /// Compute pressure from density using the Tait equation.
    #[allow(dead_code)]
    pub fn pressure_from_density(&self, rho: f64) -> f64 {
        self.b_coeff * ((rho / self.rho0).powf(self.gamma) - 1.0)
    }
    /// Apply pressure correction to particles near the outlet.
    ///
    /// Returns the corrected pressure for a particle with density `rho`
    /// at signed distance `sd` from the boundary (0 = at boundary, <0 = outside).
    #[allow(dead_code)]
    pub fn corrected_pressure(&self, rho: f64, sd: f64, damping_length: f64) -> f64 {
        let p_interior = self.pressure_from_density(rho);
        if sd >= 0.0 && sd < damping_length {
            let t = 1.0 - sd / damping_length;
            p_interior * (1.0 - t) + self.p_outlet * t
        } else if sd < 0.0 {
            self.p_outlet
        } else {
            p_interior
        }
    }
}
/// High-level manager that drives all open boundary operations each time step.
#[derive(Debug, Clone, Default)]
pub struct OpenBoundaryManager {
    /// Inflow buffers.
    pub inflows: Vec<InflowBuffer>,
    /// Outflow dampers.
    pub outflows: Vec<OutflowDamper>,
}
impl OpenBoundaryManager {
    /// Create an empty manager.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add an inflow buffer.
    pub fn add_inflow(&mut self, inflow: InflowBuffer) {
        self.inflows.push(inflow);
    }
    /// Add an outflow damper.
    pub fn add_outflow(&mut self, outflow: OutflowDamper) {
        self.outflows.push(outflow);
    }
    /// Perform one open-boundary step.
    ///
    /// 1. Inject new particles from inlet buffers.
    /// 2. Identify particles that have exited past any outlet boundary.
    ///
    /// Returns `(injected, removed)` index lists.  The caller is responsible
    /// for actually removing particles from their particle arrays (in
    /// descending index order to preserve validity).
    pub fn step(
        &mut self,
        dt: f64,
        positions: &mut Vec<[f64; 3]>,
        velocities: &mut Vec<[f64; 3]>,
        _densities: &[f64],
    ) -> (Vec<usize>, Vec<usize>) {
        let mut injected = Vec::new();
        for inflow in &mut self.inflows {
            let new_idx = inflow.inject_particles(dt, positions, velocities);
            injected.extend(new_idx);
        }
        let mut to_remove: Vec<usize> = Vec::new();
        for outflow in &self.outflows {
            for (i, &pos) in positions.iter().enumerate() {
                if outflow.zone.signed_distance(pos) < 0.0 {
                    to_remove.push(i);
                }
            }
        }
        to_remove.sort_unstable();
        to_remove.dedup();
        (injected, to_remove)
    }
    /// Apply damping forces from all outflow dampers.
    ///
    /// Returns force vectors for each particle.
    pub fn compute_damping_forces(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        for outflow in &self.outflows {
            for i in 0..n {
                let f = outflow.apply_damping_force(positions[i], velocities[i]);
                forces[i][0] += f[0];
                forces[i][1] += f[1];
                forces[i][2] += f[2];
            }
        }
        forces
    }
    /// Count total buffer particles across all inflows.
    pub fn total_buffer_count(&self) -> usize {
        self.inflows.iter().map(|inf| inf.reservoir.len()).sum()
    }
}
/// State vector for characteristic boundary conditions.
#[derive(Debug, Clone)]
pub struct CharacteristicState {
    /// Normal velocity (m/s).
    pub v_normal: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Sound speed (m/s).
    pub sound_speed: f64,
}
impl CharacteristicState {
    /// Compute the acoustic Riemann invariants.
    ///
    /// Returns `(J_plus, J_minus)` = `(v_n + p/(rho*c), v_n - p/(rho*c))`.
    #[allow(dead_code)]
    pub fn acoustic_invariants(&self) -> (f64, f64) {
        let denom = self.density * self.sound_speed;
        if denom < 1e-30 {
            return (self.v_normal, self.v_normal);
        }
        let j_plus = self.v_normal + self.pressure / denom;
        let j_minus = self.v_normal - self.pressure / denom;
        (j_plus, j_minus)
    }
    /// Apply a non-reflecting outlet BC by zeroing the incoming wave amplitude.
    ///
    /// For subsonic outflow: incoming wave `J_minus` comes from exterior (set to zero).
    /// Returns the corrected normal velocity and pressure.
    #[allow(dead_code)]
    pub fn non_reflecting_outlet(&self, reference: &CharacteristicState) -> (f64, f64) {
        let (j_plus, _) = self.acoustic_invariants();
        let denom = reference.density * reference.sound_speed;
        if denom < 1e-30 {
            return (self.v_normal, reference.pressure);
        }
        let j_minus_ref = reference.v_normal - reference.pressure / denom;
        let v_n_bc = 0.5 * (j_plus + j_minus_ref);
        let p_bc = 0.5 * denom * (j_plus - j_minus_ref);
        (v_n_bc, p_bc.max(0.0))
    }
}
/// Manages particle generation at inflow and deletion at outflow boundaries.
pub struct ParticleLifecycleManager {
    /// Maximum number of particles allowed in the simulation.
    pub max_particles: usize,
    /// Total particles injected over the simulation lifetime.
    pub total_injected: usize,
    /// Total particles removed over the simulation lifetime.
    pub total_removed: usize,
}
impl ParticleLifecycleManager {
    /// Create a new lifecycle manager.
    pub fn new(max_particles: usize) -> Self {
        Self {
            max_particles,
            total_injected: 0,
            total_removed: 0,
        }
    }
    /// Check if more particles can be injected.
    pub fn can_inject(&self, current_count: usize, batch_size: usize) -> bool {
        current_count + batch_size <= self.max_particles
    }
    /// Record injection of particles.
    pub fn record_injection(&mut self, count: usize) {
        self.total_injected += count;
    }
    /// Record removal of particles.
    pub fn record_removal(&mut self, count: usize) {
        self.total_removed += count;
    }
    /// Net particle change (injected - removed).
    pub fn net_change(&self) -> i64 {
        self.total_injected as i64 - self.total_removed as i64
    }
    /// Compute the flow rate (particles per second) given time elapsed.
    pub fn injection_rate(&self, elapsed_time: f64) -> f64 {
        if elapsed_time < 1e-30 {
            return 0.0;
        }
        self.total_injected as f64 / elapsed_time
    }
}
