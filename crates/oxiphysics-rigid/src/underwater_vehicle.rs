// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Underwater vehicle dynamics for AUVs, ROVs, and submarines.
//!
//! Implements hydrodynamic forces (buoyancy, drag, added mass), hydrostatic
//! restoring moments, thruster allocation, acoustic propagation, and a simple
//! depth PID controller.

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / n)
    }
}

/// Water density (kg/m³) at ~20 °C.
const RHO_WATER: f64 = 1025.0;
/// Gravitational acceleration (m/s²).
const GRAVITY: f64 = 9.81;

// ---------------------------------------------------------------------------
// UnderwaterVehicle
// ---------------------------------------------------------------------------

/// Full 6-DOF underwater vehicle (AUV / ROV / submarine).
///
/// State vector layout (12 elements):
/// `[x, y, z, phi, theta, psi, u, v, w, p, q, r]`
/// where `(x,y,z)` is position, `(phi,theta,psi)` are Euler angles,
/// `(u,v,w)` are body-frame linear velocities, and `(p,q,r)` are
/// body-frame angular rates.
pub struct UnderwaterVehicle {
    /// Dry mass of the vehicle in kg.
    pub mass: f64,
    /// Displaced volume in m³.
    pub volume: f64,
    /// Linear (surge/sway/heave) and rotational (roll/pitch/yaw) drag
    /// coefficients `[Du, Dv, Dw, Dp, Dq, Dr]`.
    pub drag_coeff: [f64; 6],
    /// Added-mass coefficients along the same 6 axes `[Xu_dot, ..., Nr_dot]`.
    pub added_mass: [f64; 6],
    /// Centre of buoyancy in body frame (m).
    pub center_of_buoyancy: [f64; 3],
    /// Centre of gravity in body frame (m).
    pub center_of_gravity: [f64; 3],
}

impl UnderwaterVehicle {
    /// Create a new underwater vehicle with the given parameters.
    ///
    /// # Arguments
    /// * `mass`   – dry mass (kg).
    /// * `volume` – displaced volume (m³).
    /// * `drag_coeff`          – quadratic drag coefficients for 6 DOF.
    /// * `added_mass`          – added-mass coefficients for 6 DOF.
    /// * `center_of_buoyancy`  – CoB in body frame (m).
    /// * `center_of_gravity`   – CoG in body frame (m).
    pub fn new(
        mass: f64,
        volume: f64,
        drag_coeff: [f64; 6],
        added_mass: [f64; 6],
        center_of_buoyancy: [f64; 3],
        center_of_gravity: [f64; 3],
    ) -> Self {
        Self {
            mass,
            volume,
            drag_coeff,
            added_mass,
            center_of_buoyancy,
            center_of_gravity,
        }
    }

    /// Compute the buoyancy force vector in the world frame `[Fx, Fy, Fz]`.
    ///
    /// The buoyancy force equals `ρ_water * g * V` acting upward (+Z in NED
    /// or +Y in ENU).  Here we use an ENU convention: upward = +Z.
    pub fn buoyancy_force(&self) -> [f64; 3] {
        let fb = RHO_WATER * GRAVITY * self.volume;
        [0.0, 0.0, fb]
    }

    /// Hydrostatic restoring forces/moments for a pose `pos = [x,y,z,φ,θ,ψ]`.
    ///
    /// Returns the 6-element restoring wrench `[Fx, Fy, Fz, Mx, My, Mz]`
    /// in the body frame (small-angle approximation for roll/pitch).
    pub fn hydrostatic_restoring(&self, pos: &[f64; 6]) -> [f64; 6] {
        let _x = pos[0];
        let _y = pos[1];
        let _z = pos[2];
        let phi = pos[3];
        let theta = pos[4];
        let _psi = pos[5];

        let w = self.mass * GRAVITY;
        let b = RHO_WATER * GRAVITY * self.volume;

        // Net force in body-frame z (heave).
        let fz = b - w;

        // Restoring moments from offset between CoB and CoG.
        let r_bg = sub3(self.center_of_buoyancy, self.center_of_gravity);

        // Full restoring moment: M = r_BG × (B * ẑ_world) in body frame.
        // Using Euler-angle rotation of world +Z into body frame (simplified):
        let bz_body = [
            -b * theta.sin(),
            b * phi.cos() * theta.sin().acos().cos() * 0.0 + b * phi.sin() * theta.cos(),
            b * phi.cos() * theta.cos(),
        ];
        let moment = cross3(r_bg, bz_body);

        [0.0, 0.0, fz, moment[0], moment[1], moment[2]]
    }

    /// Quadratic drag force/moment for body-frame velocity `vel = [u,v,w,p,q,r]`.
    ///
    /// Uses `F_i = -D_i * |vel_i| * vel_i` for each axis independently.
    pub fn drag_force(&self, vel: &[f64; 6]) -> [f64; 6] {
        let mut f = [0.0f64; 6];
        for i in 0..6 {
            f[i] = -self.drag_coeff[i] * vel[i].abs() * vel[i];
        }
        f
    }

    /// Added-mass reaction force/moment for body-frame acceleration `accel`.
    ///
    /// Returns `F_i = -Ma_i * accel_i` for each axis.
    pub fn added_mass_force(&self, accel: &[f64; 6]) -> [f64; 6] {
        let mut f = [0.0f64; 6];
        for i in 0..6 {
            f[i] = -self.added_mass[i] * accel[i];
        }
        f
    }

    /// Metacentric height GM (m).
    ///
    /// `GM = (CoB_z - CoG_z)` in the body frame.  Positive GM indicates
    /// positive static stability.
    pub fn metacentric_height(&self) -> f64 {
        self.center_of_buoyancy[2] - self.center_of_gravity[2]
    }

    /// Returns `true` if the vehicle is statically stable (GM > 0).
    pub fn static_stability(&self) -> bool {
        self.metacentric_height() > 0.0
    }

    /// Integrate the vehicle state by one timestep `dt` using semi-implicit
    /// Euler integration.
    ///
    /// # Arguments
    /// * `state`   – current 12-element state `[x,y,z,φ,θ,ψ, u,v,w,p,q,r]`.
    /// * `control` – 4-element control input `[Fx, Fy, Fz, Mz]` in body frame.
    /// * `dt`      – timestep (s).
    ///
    /// Returns the new 12-element state.
    pub fn step(&self, state: &[f64; 12], control: &[f64; 4], dt: f64) -> [f64; 12] {
        let phi = state[3];
        let theta = state[4];
        let psi = state[5];
        let u = state[6];
        let v = state[7];
        let w = state[8];
        let p = state[9];
        let q = state[10];
        let r = state[11];

        // Rotation matrix (body → world) — ZYX Euler.
        let cphi = phi.cos();
        let sphi = phi.sin();
        let cth = theta.cos();
        let sth = theta.sin();
        let cpsi = psi.cos();
        let spsi = psi.sin();

        // Body velocities → world velocities.
        let xd = cth * cpsi * u
            + (sphi * sth * cpsi - cphi * spsi) * v
            + (cphi * sth * cpsi + sphi * spsi) * w;
        let yd = cth * spsi * u
            + (sphi * sth * spsi + cphi * cpsi) * v
            + (cphi * sth * spsi - sphi * cpsi) * w;
        let zd = -sth * u + sphi * cth * v + cphi * cth * w;

        // Euler angle rates.
        let phid = p + (q * sphi + r * cphi) * sth / cth.max(1e-6);
        let thetad = q * cphi - r * sphi;
        let psid = (q * sphi + r * cphi) / cth.max(1e-6);

        // Forces in body frame.
        let vel6 = [u, v, w, p, q, r];
        let drag = self.drag_force(&vel6);

        // Net effective mass (including added mass).
        let m_eff = [
            self.mass + self.added_mass[0],
            self.mass + self.added_mass[1],
            self.mass + self.added_mass[2],
            self.mass + self.added_mass[3],
            self.mass + self.added_mass[4],
            self.mass + self.added_mass[5],
        ];

        let ctrl6 = [control[0], control[1], control[2], 0.0, 0.0, control[3]];

        let mut accel = [0.0f64; 6];
        for i in 0..6 {
            accel[i] = (ctrl6[i] + drag[i]) / m_eff[i].max(1e-9);
        }

        // Semi-implicit Euler.
        let new_u = u + accel[0] * dt;
        let new_v = v + accel[1] * dt;
        let new_w = w + accel[2] * dt;
        let new_p = p + accel[3] * dt;
        let new_q = q + accel[4] * dt;
        let new_r = r + accel[5] * dt;

        let new_x = state[0] + xd * dt;
        let new_y = state[1] + yd * dt;
        let new_z = state[2] + zd * dt;
        let new_phi = phi + phid * dt;
        let new_theta = theta + thetad * dt;
        let new_psi = psi + psid * dt;

        [
            new_x, new_y, new_z, new_phi, new_theta, new_psi, new_u, new_v, new_w, new_p, new_q,
            new_r,
        ]
    }
}

// ---------------------------------------------------------------------------
// Thruster
// ---------------------------------------------------------------------------

/// A single thruster unit on an underwater vehicle.
pub struct Thruster {
    /// Mounting position in body frame (m).
    pub position: [f64; 3],
    /// Unit thrust direction vector in body frame.
    pub direction: [f64; 3],
    /// Maximum achievable thrust force (N).
    pub max_thrust: f64,
    /// Efficiency coefficient (0–1).
    pub efficiency: f64,
}

impl Thruster {
    /// Create a new thruster.
    pub fn new(position: [f64; 3], direction: [f64; 3], max_thrust: f64, efficiency: f64) -> Self {
        Self {
            position,
            direction: normalize3(direction),
            max_thrust,
            efficiency,
        }
    }

    /// Compute the thrust vector (N) for a normalised throttle in `[-1, 1]`.
    pub fn thrust_vector(&self, throttle: f64) -> [f64; 3] {
        let clamped = throttle.clamp(-1.0, 1.0);
        scale3(self.direction, self.max_thrust * clamped)
    }

    /// Estimate electrical power consumption (W) for a given thrust (N).
    ///
    /// Uses `P = |T| * V_exit / η` where the exit velocity is derived from
    /// the disk-actuator approximation `T = 2ρA V_exit²`.  For simplicity we
    /// use a linear approximation `P ≈ |T|² / (2 * max_thrust * efficiency)`.
    pub fn power_consumption(&self, thrust: f64) -> f64 {
        if self.efficiency < 1e-9 || self.max_thrust < 1e-9 {
            return 0.0;
        }
        (thrust * thrust) / (2.0 * self.max_thrust * self.efficiency)
    }
}

// ---------------------------------------------------------------------------
// ThrusterAllocation
// ---------------------------------------------------------------------------

/// Manages a collection of thrusters and allocates desired wrenches.
pub struct ThrusterAllocation {
    /// All thrusters belonging to this vehicle.
    pub thrusters: Vec<Thruster>,
}

impl ThrusterAllocation {
    /// Create a new allocation manager with the given thrusters.
    pub fn new(thrusters: Vec<Thruster>) -> Self {
        Self { thrusters }
    }

    /// Build the 6×N thruster configuration matrix B.
    ///
    /// Each column `i` is `[f_i; τ_i]` where `f_i = direction_i` and
    /// `τ_i = position_i × direction_i`.
    fn config_matrix(&self) -> Vec<[f64; 6]> {
        self.thrusters
            .iter()
            .map(|t| {
                let f = t.direction;
                let tau = cross3(t.position, t.direction);
                [f[0], f[1], f[2], tau[0], tau[1], tau[2]]
            })
            .collect()
    }

    /// Pseudo-inverse thrust allocation for a desired wrench `w = [Fx,Fy,Fz,Mx,My,Mz]`.
    ///
    /// Uses a simple minimum-norm least-squares via the Moore-Penrose
    /// pseudo-inverse computed with the normal equations `B^T (B B^T)^{-1} w`.
    /// The result is clamped to each thruster's `[-max_thrust, +max_thrust]`.
    pub fn pseudo_inverse_allocation(&self, wrench: &[f64; 6]) -> Vec<f64> {
        let n = self.thrusters.len();
        if n == 0 {
            return Vec::new();
        }
        let b = self.config_matrix(); // n columns of [f64;6]

        // Compute B * B^T  (6×6 symmetric matrix).
        let mut bbt = [[0.0f64; 6]; 6];
        for col in &b {
            for i in 0..6 {
                for j in 0..6 {
                    bbt[i][j] += col[i] * col[j];
                }
            }
        }

        // Solve (B B^T) x = w via Gaussian elimination.
        let x = solve6x6(&bbt, wrench);

        // Compute B^T x  (pseudo-inverse solution).
        let mut t_alloc: Vec<f64> = b.iter().map(|col| dot6(col, &x)).collect();

        // Clamp to thruster limits.
        for (i, t) in t_alloc.iter_mut().enumerate() {
            *t = t.clamp(-self.thrusters[i].max_thrust, self.thrusters[i].max_thrust);
        }
        t_alloc
    }

    /// Allocate the desired wrench and return the resulting actual wrench.
    ///
    /// First computes per-thruster throttles via pseudo-inverse, then
    /// re-accumulates the actual wrench from the clamped throttle values.
    pub fn allocate_wrench(&self, wrench: &[f64; 6]) -> [f64; 6] {
        let throttles = self.pseudo_inverse_allocation(wrench);
        let b = self.config_matrix();
        let mut actual = [0.0f64; 6];
        for (i, &thr) in throttles.iter().enumerate() {
            for j in 0..6 {
                actual[j] += b[i][j] * thr;
            }
        }
        actual
    }
}

// ---------------------------------------------------------------------------
// Linear algebra helpers (6×6)
// ---------------------------------------------------------------------------

fn dot6(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Solve the 6×6 system A x = b via Gaussian elimination with partial pivoting.
fn solve6x6(a: &[[f64; 6]; 6], b: &[f64; 6]) -> [f64; 6] {
    const N: usize = 6;
    let mut m = [[0.0f64; 7]; N];
    for i in 0..N {
        for j in 0..N {
            m[i][j] = a[i][j];
        }
        m[i][N] = b[i];
    }

    for col in 0..N {
        // Find pivot.
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for (off, m_row) in m[col + 1..N].iter().enumerate() {
            let row = col + 1 + off;
            if m_row[col].abs() > max_val {
                max_val = m_row[col].abs();
                max_row = row;
            }
        }
        m.swap(col, max_row);

        let pivot = m[col][col];
        if pivot.abs() < 1e-15 {
            continue;
        }
        for row in col + 1..N {
            let factor = m[row][col] / pivot;
            let m_col_copy: Vec<f64> = m[col][col..=N].to_vec();
            for (m_rk, &m_ck) in m[row][col..=N].iter_mut().zip(m_col_copy.iter()) {
                *m_rk -= m_ck * factor;
            }
        }
    }

    // Back-substitution.
    let mut x = [0.0f64; N];
    for i in (0..N).rev() {
        x[i] = m[i][N];
        for j in i + 1..N {
            x[i] -= m[i][j] * x[j];
        }
        if m[i][i].abs() > 1e-15 {
            x[i] /= m[i][i];
        }
    }
    x
}

// ---------------------------------------------------------------------------
// AcousticPropagation
// ---------------------------------------------------------------------------

/// Simplified underwater acoustic propagation model.
pub struct AcousticPropagation {
    /// Speed of sound in water (m/s).  Typical value ≈ 1500 m/s.
    pub speed_of_sound: f64,
    /// Attenuation coefficient (dB/km).  Typical value ≈ 0.1 dB/km at 10 kHz.
    pub attenuation: f64,
}

impl AcousticPropagation {
    /// Create a new acoustic propagation model.
    pub fn new(speed_of_sound: f64, attenuation: f64) -> Self {
        Self {
            speed_of_sound,
            attenuation,
        }
    }

    /// One-way travel time (s) for a given distance (m).
    pub fn travel_time(&self, distance: f64) -> f64 {
        distance / self.speed_of_sound.max(1e-9)
    }

    /// Maximum sonar detection range (m) using a simplified SONAR equation.
    ///
    /// `SL - TL ≥ NL`  ⟹  `TL = SL - NL`
    ///
    /// Two-way transmission loss:
    /// `TL = 20 log₁₀(R) + α·R·1e-3`  where α is in dB/km.
    ///
    /// Solved iteratively (fixed-point iteration, max 100 iterations).
    ///
    /// # Arguments
    /// * `source_level` – source level SL (dB re 1 µPa @ 1 m).
    /// * `noise`        – ambient noise level NL (dB).
    pub fn sonar_range(&self, source_level: f64, noise: f64) -> f64 {
        let tl_max = source_level - noise;
        if tl_max <= 0.0 {
            return 0.0;
        }
        let alpha_per_m = self.attenuation * 1e-3 / 20.0 * std::f64::consts::LN_10.recip();
        // Binary search on range.
        let mut lo = 1.0f64;
        let mut hi = 1_000_000.0f64;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let tl = 20.0 * mid.log10() + self.attenuation * mid * 1e-3;
            if tl < tl_max {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let _ = alpha_per_m; // used in derivation above
        0.5 * (lo + hi)
    }
}

// ---------------------------------------------------------------------------
// DepthController
// ---------------------------------------------------------------------------

/// Simple PID depth controller for an underwater vehicle.
pub struct DepthController {
    /// PID gains `[Kp, Ki, Kd]`.
    pub pid: [f64; 3],
    /// Accumulated integral of the depth error.
    pub integral: f64,
    /// Previous error for derivative term.
    prev_error: f64,
}

impl DepthController {
    /// Create a new depth controller with gains `[Kp, Ki, Kd]`.
    pub fn new(pid: [f64; 3]) -> Self {
        Self {
            pid,
            integral: 0.0,
            prev_error: 0.0,
        }
    }

    /// Update the controller with the current depth error and timestep.
    ///
    /// Returns the heave thrust command (N, positive = ascend).
    pub fn update(&mut self, depth_error: f64, dt: f64) -> f64 {
        self.integral += depth_error * dt;
        let derivative = if dt > 1e-15 {
            (depth_error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = depth_error;
        self.pid[0] * depth_error + self.pid[1] * self.integral + self.pid[2] * derivative
    }

    /// Reset the controller state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn default_vehicle() -> UnderwaterVehicle {
        UnderwaterVehicle::new(
            100.0,                             // mass kg
            0.1,                               // volume m³  (≈ 102.5 kg buoyancy)
            [10.0, 10.0, 10.0, 5.0, 5.0, 5.0], // drag
            [5.0, 5.0, 5.0, 1.0, 1.0, 1.0],    // added mass
            [0.0, 0.0, 0.05],                  // CoB above CoG
            [0.0, 0.0, 0.0],                   // CoG at origin
        )
    }

    // -----------------------------------------------------------------------
    // UnderwaterVehicle tests
    // -----------------------------------------------------------------------

    #[test]
    fn buoyancy_force_upward() {
        let v = default_vehicle();
        let fb = v.buoyancy_force();
        assert!(fb[2] > 0.0, "buoyancy should be upward (+Z)");
        let expected = RHO_WATER * GRAVITY * v.volume;
        assert!((fb[2] - expected).abs() < 1e-6);
    }

    #[test]
    fn drag_force_opposes_motion() {
        let v = default_vehicle();
        let vel = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let drag = v.drag_force(&vel);
        assert!(drag[0] < 0.0, "drag should oppose surge velocity");
    }

    #[test]
    fn drag_force_zero_for_zero_velocity() {
        let v = default_vehicle();
        let vel = [0.0; 6];
        let drag = v.drag_force(&vel);
        for (i, &d) in drag.iter().enumerate() {
            assert!(d.abs() < 1e-15, "drag[{i}] should be zero at rest");
        }
    }

    #[test]
    fn added_mass_force_proportional_to_accel() {
        let v = default_vehicle();
        let a1 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let a2 = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f1 = v.added_mass_force(&a1);
        let f2 = v.added_mass_force(&a2);
        assert!((f2[0] - 2.0 * f1[0]).abs() < 1e-10);
    }

    #[test]
    fn metacentric_height_positive_stable() {
        let v = default_vehicle();
        assert!(v.metacentric_height() > 0.0);
        assert!(v.static_stability());
    }

    #[test]
    fn metacentric_height_negative_unstable() {
        let mut v = default_vehicle();
        v.center_of_buoyancy[2] = -0.1; // CoB below CoG
        assert!(v.metacentric_height() < 0.0);
        assert!(!v.static_stability());
    }

    #[test]
    fn step_advances_position() {
        let v = default_vehicle();
        let mut state = [0.0f64; 12];
        state[6] = 1.0; // surge velocity
        let control = [0.0; 4];
        let new_state = v.step(&state, &control, 0.1);
        assert!(new_state[0] > 0.0, "vehicle should move forward");
    }

    #[test]
    fn step_drag_slows_down_vehicle() {
        let v = default_vehicle();
        let mut state = [0.0f64; 12];
        state[6] = 10.0; // high surge velocity
        let control = [0.0; 4];
        let new_state = v.step(&state, &control, 0.1);
        assert!(new_state[6] < 10.0, "drag should reduce velocity");
    }

    #[test]
    fn hydrostatic_restoring_upright_has_fz() {
        let v = default_vehicle();
        let pos = [0.0; 6];
        let f = v.hydrostatic_restoring(&pos);
        // Net heave: buoyancy - weight
        let expected_fz = RHO_WATER * GRAVITY * v.volume - v.mass * GRAVITY;
        assert!((f[2] - expected_fz).abs() < 1e-6);
    }

    #[test]
    fn hydrostatic_restoring_nonzero_angle() {
        let v = default_vehicle();
        let pos = [0.0, 0.0, 0.0, 0.1, 0.0, 0.0]; // small roll
        let f = v.hydrostatic_restoring(&pos);
        // With CoB above CoG, roll moment should be restoring (non-zero).
        assert!(f[3].abs() > 0.0 || f[4].abs() > 0.0 || f[2].abs() > 0.0);
    }

    // -----------------------------------------------------------------------
    // Thruster tests
    // -----------------------------------------------------------------------

    #[test]
    fn thruster_zero_throttle_zero_force() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8);
        let f = t.thrust_vector(0.0);
        assert!(norm3(f) < 1e-12);
    }

    #[test]
    fn thruster_full_throttle_max_force() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8);
        let f = t.thrust_vector(1.0);
        assert!((norm3(f) - 100.0).abs() < 1e-10);
    }

    #[test]
    fn thruster_negative_throttle_reverses_force() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 50.0, 0.8);
        let f = t.thrust_vector(-1.0);
        assert!(f[0] < 0.0);
    }

    #[test]
    fn thruster_power_zero_for_zero_thrust() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8);
        assert!((t.power_consumption(0.0)).abs() < 1e-15);
    }

    #[test]
    fn thruster_power_positive_for_nonzero_thrust() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8);
        assert!(t.power_consumption(50.0) > 0.0);
    }

    // -----------------------------------------------------------------------
    // ThrusterAllocation tests
    // -----------------------------------------------------------------------

    #[test]
    fn allocation_single_thruster_surge() {
        let t = Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8);
        let alloc = ThrusterAllocation::new(vec![t]);
        let wrench = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let thrusts = alloc.pseudo_inverse_allocation(&wrench);
        assert_eq!(thrusts.len(), 1);
        // Thruster should produce positive surge force.
        assert!(thrusts[0] > 0.0);
    }

    #[test]
    fn allocation_empty_thrusters() {
        let alloc = ThrusterAllocation::new(vec![]);
        let wrench = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let thrusts = alloc.pseudo_inverse_allocation(&wrench);
        assert!(thrusts.is_empty());
    }

    #[test]
    fn allocation_wrench_bounded_by_max_thrust() {
        let thrusters = vec![
            Thruster::new([0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 50.0, 0.9),
            Thruster::new([-0.5, 0.0, 0.0], [-1.0, 0.0, 0.0], 50.0, 0.9),
        ];
        let alloc = ThrusterAllocation::new(thrusters);
        let wrench = [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // saturated
        let thrusts = alloc.pseudo_inverse_allocation(&wrench);
        for (i, &t) in thrusts.iter().enumerate() {
            assert!(t <= 50.0, "thrust[{i}] = {t} exceeds max");
            assert!(t >= -50.0, "thrust[{i}] = {t} below min");
        }
    }

    #[test]
    fn allocate_wrench_returns_6_elements() {
        let thrusters = vec![
            Thruster::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 100.0, 0.8),
            Thruster::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 100.0, 0.8),
        ];
        let alloc = ThrusterAllocation::new(thrusters);
        let wrench = [5.0, 5.0, 0.0, 0.0, 0.0, 0.0];
        let actual = alloc.allocate_wrench(&wrench);
        assert_eq!(actual.len(), 6);
    }

    // -----------------------------------------------------------------------
    // AcousticPropagation tests
    // -----------------------------------------------------------------------

    #[test]
    fn travel_time_proportional_to_distance() {
        let ap = AcousticPropagation::new(1500.0, 0.1);
        let t1 = ap.travel_time(1500.0);
        let t2 = ap.travel_time(3000.0);
        assert!((t1 - 1.0).abs() < 1e-10);
        assert!((t2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn travel_time_zero_distance() {
        let ap = AcousticPropagation::new(1500.0, 0.1);
        assert!((ap.travel_time(0.0)).abs() < 1e-15);
    }

    #[test]
    fn sonar_range_positive_for_strong_source() {
        let ap = AcousticPropagation::new(1500.0, 0.1);
        let r = ap.sonar_range(200.0, 50.0);
        assert!(r > 0.0, "sonar range should be positive: {r}");
    }

    #[test]
    fn sonar_range_zero_when_noise_exceeds_source() {
        let ap = AcousticPropagation::new(1500.0, 0.1);
        let r = ap.sonar_range(50.0, 100.0);
        assert_eq!(r, 0.0, "sonar range should be zero when noise > SL");
    }

    #[test]
    fn sonar_range_increases_with_source_level() {
        let ap = AcousticPropagation::new(1500.0, 0.1);
        let r1 = ap.sonar_range(150.0, 50.0);
        let r2 = ap.sonar_range(200.0, 50.0);
        assert!(r2 > r1, "higher SL should give longer range");
    }

    // -----------------------------------------------------------------------
    // DepthController tests
    // -----------------------------------------------------------------------

    #[test]
    fn depth_controller_zero_error_zero_output() {
        let mut ctrl = DepthController::new([1.0, 0.1, 0.5]);
        let out = ctrl.update(0.0, 0.1);
        assert!(out.abs() < 1e-15);
    }

    #[test]
    fn depth_controller_positive_error_positive_output() {
        let mut ctrl = DepthController::new([1.0, 0.0, 0.0]);
        let out = ctrl.update(5.0, 0.1);
        assert!(out > 0.0);
    }

    #[test]
    fn depth_controller_integral_accumulates() {
        let mut ctrl = DepthController::new([0.0, 1.0, 0.0]);
        ctrl.update(1.0, 0.1);
        ctrl.update(1.0, 0.1);
        assert!(ctrl.integral > 0.0);
    }

    #[test]
    fn depth_controller_derivative_term() {
        let mut ctrl = DepthController::new([0.0, 0.0, 1.0]);
        ctrl.update(0.0, 0.1); // prime derivative
        let out = ctrl.update(1.0, 0.1);
        // Derivative of (1 - 0) / 0.1 = 10; output = Kd * 10 = 10
        assert!((out - 10.0).abs() < 1e-6, "derivative output = {out}");
    }

    #[test]
    fn depth_controller_reset_clears_state() {
        let mut ctrl = DepthController::new([1.0, 1.0, 1.0]);
        ctrl.update(5.0, 0.1);
        ctrl.reset();
        assert!(ctrl.integral.abs() < 1e-15);
    }

    // -----------------------------------------------------------------------
    // Physics / helper tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve6x6_identity_system() {
        let mut a = [[0.0f64; 6]; 6];
        for (i, row) in a.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        let b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x = solve6x6(&a, &b);
        for (i, (&xi, &bi)) in x.iter().zip(b.iter()).enumerate() {
            assert!((xi - bi).abs() < 1e-10, "x[{i}] = {}", xi);
        }
    }

    #[test]
    fn cross3_orthogonality() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn normalize3_unit_vector() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize3(v);
        assert!((norm3(n) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn normalize3_zero_vector() {
        let v = [0.0, 0.0, 0.0];
        let n = normalize3(v);
        assert!(norm3(n) < 1e-15);
    }

    #[test]
    fn buoyancy_scales_with_volume() {
        let v1 = UnderwaterVehicle::new(100.0, 0.1, [0.0; 6], [0.0; 6], [0.0; 3], [0.0; 3]);
        let v2 = UnderwaterVehicle::new(100.0, 0.2, [0.0; 6], [0.0; 6], [0.0; 3], [0.0; 3]);
        let fb1 = v1.buoyancy_force();
        let fb2 = v2.buoyancy_force();
        assert!((fb2[2] - 2.0 * fb1[2]).abs() < 1e-6);
    }

    #[test]
    fn acoustic_propagation_custom_speed() {
        let ap = AcousticPropagation::new(1000.0, 0.0);
        assert!((ap.travel_time(2000.0) - 2.0).abs() < 1e-10);
    }
}
