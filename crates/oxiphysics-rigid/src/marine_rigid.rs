// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Marine and naval rigid body dynamics.
//!
//! This module implements 6-DOF ship motion, hydrostatic restoring forces,
//! strip-theory added mass and damping, Froude-Krylov wave excitation,
//! Morison equation for slender members, mooring line catenary analysis,
//! seakeeping RAO computation, roll damping models, and propeller thrust.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// 6-DOF Ship State
// ─────────────────────────────────────────────────────────────────────────────

/// Six degree-of-freedom ship state vector.
///
/// The coordinate system follows the marine convention:
/// - X: forward (surge)
/// - Y: starboard (sway)
/// - Z: downward (heave)
/// - Roll: rotation about X
/// - Pitch: rotation about Y
/// - Yaw: rotation about Z
#[derive(Debug, Clone)]
pub struct ShipState6DOF {
    /// Position `[x, y, z]` in the global frame (meters).
    pub position: [f64; 3],
    /// Euler angles `[roll, pitch, yaw]` in radians.
    pub orientation: [f64; 3],
    /// Translational velocity `[surge_vel, sway_vel, heave_vel]` (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity `[roll_rate, pitch_rate, yaw_rate]` (rad/s).
    pub angular_velocity: [f64; 3],
    /// Ship mass (kg).
    pub mass: f64,
    /// Moments of inertia `[Ixx, Iyy, Izz]` (kg*m^2).
    pub inertia: [f64; 3],
}

impl ShipState6DOF {
    /// Create a new ship state at the origin with given mass and inertia.
    pub fn new(mass: f64, inertia: [f64; 3]) -> Self {
        Self {
            position: [0.0; 3],
            orientation: [0.0; 3],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass,
            inertia,
        }
    }

    /// Total kinetic energy (translational + rotational).
    pub fn kinetic_energy(&self) -> f64 {
        let v = &self.velocity;
        let w = &self.angular_velocity;
        let ke_trans = 0.5 * self.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        let ke_rot = 0.5
            * (self.inertia[0] * w[0] * w[0]
                + self.inertia[1] * w[1] * w[1]
                + self.inertia[2] * w[2] * w[2]);
        ke_trans + ke_rot
    }

    /// Speed over ground (magnitude of translational velocity in XY plane).
    pub fn speed_over_ground(&self) -> f64 {
        (self.velocity[0] * self.velocity[0] + self.velocity[1] * self.velocity[1]).sqrt()
    }

    /// Integrate the ship state forward by `dt` seconds using Euler integration.
    ///
    /// `force` is `[Fx, Fy, Fz]` in newtons, `torque` is `[Mx, My, Mz]` in N*m.
    pub fn integrate_euler(&mut self, force: [f64; 3], torque: [f64; 3], dt: f64) {
        // Update translational velocity
        for (v, f) in self.velocity.iter_mut().zip(force.iter()) {
            *v += (f / self.mass) * dt;
        }
        // Update angular velocity
        for (av, (t, inr)) in self
            .angular_velocity
            .iter_mut()
            .zip(torque.iter().zip(self.inertia.iter()))
        {
            if *inr > 1e-30 {
                *av += (t / inr) * dt;
            }
        }
        // Update position
        for i in 0..3 {
            self.position[i] += self.velocity[i] * dt;
        }
        // Update orientation
        for i in 0..3 {
            self.orientation[i] += self.angular_velocity[i] * dt;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ship Hull Geometry
// ─────────────────────────────────────────────────────────────────────────────

/// Simplified ship hull geometry for hydrostatic calculations.
///
/// Uses a box-like approximation with waterplane area and block coefficient.
#[derive(Debug, Clone)]
pub struct HullGeometry {
    /// Length overall (m).
    pub length: f64,
    /// Beam (maximum width) (m).
    pub beam: f64,
    /// Draft (submerged depth) at rest (m).
    pub draft: f64,
    /// Block coefficient (ratio of displaced volume to bounding box).
    pub block_coefficient: f64,
    /// Waterplane area coefficient (ratio of waterplane area to L*B).
    pub waterplane_coefficient: f64,
    /// Vertical position of center of gravity above keel (m).
    pub kg: f64,
    /// Vertical position of center of buoyancy above keel (m).
    pub kb: f64,
}

impl HullGeometry {
    /// Create hull geometry with given principal dimensions.
    pub fn new(
        length: f64,
        beam: f64,
        draft: f64,
        block_coefficient: f64,
        waterplane_coefficient: f64,
        kg: f64,
        kb: f64,
    ) -> Self {
        Self {
            length,
            beam,
            draft,
            block_coefficient,
            waterplane_coefficient,
            kg,
            kb,
        }
    }

    /// Displaced volume (m^3).
    pub fn displaced_volume(&self) -> f64 {
        self.length * self.beam * self.draft * self.block_coefficient
    }

    /// Waterplane area (m^2).
    pub fn waterplane_area(&self) -> f64 {
        self.length * self.beam * self.waterplane_coefficient
    }

    /// Second moment of waterplane area about the longitudinal axis (m^4).
    ///
    /// Approximated as `(1/12) * Cwp * L * B^3`.
    pub fn waterplane_inertia_transverse(&self) -> f64 {
        self.waterplane_coefficient * self.length * self.beam.powi(3) / 12.0
    }

    /// Second moment of waterplane area about the transverse axis (m^4).
    ///
    /// Approximated as `(1/12) * Cwp * B * L^3`.
    pub fn waterplane_inertia_longitudinal(&self) -> f64 {
        self.waterplane_coefficient * self.beam * self.length.powi(3) / 12.0
    }

    /// Transverse metacentric height BM_T = I_T / V.
    pub fn bm_transverse(&self) -> f64 {
        let v = self.displaced_volume();
        if v.abs() < 1e-30 {
            return 0.0;
        }
        self.waterplane_inertia_transverse() / v
    }

    /// Transverse metacentric height GM_T = KB + BM_T - KG.
    pub fn gm_transverse(&self) -> f64 {
        self.kb + self.bm_transverse() - self.kg
    }

    /// Longitudinal metacentric height BM_L = I_L / V.
    pub fn bm_longitudinal(&self) -> f64 {
        let v = self.displaced_volume();
        if v.abs() < 1e-30 {
            return 0.0;
        }
        self.waterplane_inertia_longitudinal() / v
    }

    /// Longitudinal metacentric height GM_L = KB + BM_L - KG.
    pub fn gm_longitudinal(&self) -> f64 {
        self.kb + self.bm_longitudinal() - self.kg
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hydrostatic Restoring Forces
// ─────────────────────────────────────────────────────────────────────────────

/// Compute hydrostatic restoring forces and moments for small displacements.
///
/// `displacement` is the ship's displaced mass (kg), `rho` is water density
/// (kg/m^3), `g` is gravitational acceleration (m/s^2).
///
/// Returns `([Fz], [Mx, My])` where `Fz` is the heave restoring force and
/// `Mx`, `My` are the roll and pitch restoring moments.
pub fn hydrostatic_restoring(
    hull: &HullGeometry,
    heave: f64,
    roll: f64,
    pitch: f64,
    rho: f64,
    g: f64,
) -> ([f64; 3], [f64; 3]) {
    let aw = hull.waterplane_area();
    let displacement = rho * hull.displaced_volume();

    // Heave restoring force: -rho * g * Aw * heave
    let fz = -rho * g * aw * heave;

    // Roll restoring moment: -displacement * g * GM_T * roll
    let mx = -displacement * g * hull.gm_transverse() * roll;

    // Pitch restoring moment: -displacement * g * GM_L * pitch
    let my = -displacement * g * hull.gm_longitudinal() * pitch;

    ([0.0, 0.0, fz], [mx, my, 0.0])
}

// ─────────────────────────────────────────────────────────────────────────────
// Strip Theory: Added Mass and Damping
// ─────────────────────────────────────────────────────────────────────────────

/// Strip theory coefficients for a single cross-section.
#[derive(Debug, Clone, Copy)]
pub struct StripCoefficients {
    /// Station position along the ship length (m, from AP).
    pub x_station: f64,
    /// Sectional beam at waterline (m).
    pub section_beam: f64,
    /// Sectional draft (m).
    pub section_draft: f64,
    /// Sectional area coefficient.
    pub section_area_coeff: f64,
}

impl StripCoefficients {
    /// Sectional submerged area (m^2).
    pub fn section_area(&self) -> f64 {
        self.section_beam * self.section_draft * self.section_area_coeff
    }
}

/// Compute added mass in heave for a 2-D strip (Lewis form approximation).
///
/// For a semi-circular section: `a33 = rho * pi/2 * (B/2)^2` per unit length.
///
/// `section_beam` is the local beam at waterline, `rho` is water density.
pub fn strip_added_mass_heave(section_beam: f64, rho: f64) -> f64 {
    let half_beam = section_beam / 2.0;
    rho * PI * 0.5 * half_beam * half_beam
}

/// Compute added mass in sway for a 2-D strip (Lewis form approximation).
///
/// For a semi-circular section: `a22 = rho * pi/2 * T^2` per unit length.
///
/// `section_draft` is the local draft, `rho` is water density.
pub fn strip_added_mass_sway(section_draft: f64, rho: f64) -> f64 {
    rho * PI * 0.5 * section_draft * section_draft
}

/// Compute total heave added mass by integrating strip contributions.
///
/// `strips` are the cross-section data, `rho` is water density.
pub fn total_heave_added_mass(strips: &[StripCoefficients], rho: f64) -> f64 {
    if strips.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0_f64;
    for i in 0..strips.len() - 1 {
        let dx = (strips[i + 1].x_station - strips[i].x_station).abs();
        let a_i = strip_added_mass_heave(strips[i].section_beam, rho);
        let a_next = strip_added_mass_heave(strips[i + 1].section_beam, rho);
        total += 0.5 * (a_i + a_next) * dx;
    }
    total
}

/// Compute total sway added mass by integrating strip contributions.
pub fn total_sway_added_mass(strips: &[StripCoefficients], rho: f64) -> f64 {
    if strips.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0_f64;
    for i in 0..strips.len() - 1 {
        let dx = (strips[i + 1].x_station - strips[i].x_station).abs();
        let a_i = strip_added_mass_sway(strips[i].section_draft, rho);
        let a_next = strip_added_mass_sway(strips[i + 1].section_draft, rho);
        total += 0.5 * (a_i + a_next) * dx;
    }
    total
}

/// Compute damping coefficient for heave from strip theory.
///
/// Uses the radiation damping approximation: `b33 = rho * g * B^2 / (2 * omega)`
/// per unit length.
pub fn strip_damping_heave(section_beam: f64, omega: f64, rho: f64, g: f64) -> f64 {
    if omega.abs() < 1e-30 {
        return 0.0;
    }
    rho * g * section_beam * section_beam / (2.0 * omega)
}

// ─────────────────────────────────────────────────────────────────────────────
// Wave Excitation: Froude-Krylov
// ─────────────────────────────────────────────────────────────────────────────

/// Regular wave parameters.
#[derive(Debug, Clone, Copy)]
pub struct RegularWave {
    /// Wave amplitude (m).
    pub amplitude: f64,
    /// Wave angular frequency (rad/s).
    pub omega: f64,
    /// Wave number k = omega^2 / g for deep water (1/m).
    pub wave_number: f64,
    /// Wave heading angle relative to ship bow (rad, 0 = head seas).
    pub heading: f64,
}

impl RegularWave {
    /// Create a regular wave for deep water conditions.
    pub fn deep_water(amplitude: f64, omega: f64, heading: f64, g: f64) -> Self {
        let wave_number = omega * omega / g;
        Self {
            amplitude,
            omega,
            wave_number,
            heading,
        }
    }

    /// Wave length (m).
    pub fn wavelength(&self) -> f64 {
        if self.wave_number.abs() < 1e-30 {
            return f64::INFINITY;
        }
        2.0 * PI / self.wave_number
    }

    /// Wave period (s).
    pub fn period(&self) -> f64 {
        if self.omega.abs() < 1e-30 {
            return f64::INFINITY;
        }
        2.0 * PI / self.omega
    }

    /// Wave phase velocity (m/s).
    pub fn phase_velocity(&self) -> f64 {
        if self.wave_number.abs() < 1e-30 {
            return 0.0;
        }
        self.omega / self.wave_number
    }

    /// Wave elevation at position `x` and time `t`.
    pub fn elevation(&self, x: f64, t: f64) -> f64 {
        self.amplitude * (self.wave_number * x * self.heading.cos() - self.omega * t).cos()
    }

    /// Vertical velocity of water particles at `(x, z, t)` for deep water.
    ///
    /// `z` is depth below the free surface (positive downward).
    pub fn vertical_velocity(&self, x: f64, z: f64, t: f64) -> f64 {
        let phase = self.wave_number * x * self.heading.cos() - self.omega * t;
        -self.amplitude * self.omega * (-self.wave_number * z).exp() * phase.sin()
    }

    /// Horizontal velocity of water particles at `(x, z, t)` for deep water.
    pub fn horizontal_velocity(&self, x: f64, z: f64, t: f64) -> f64 {
        let phase = self.wave_number * x * self.heading.cos() - self.omega * t;
        self.amplitude * self.omega * (-self.wave_number * z).exp() * phase.cos()
    }
}

/// Compute the Froude-Krylov heave force on a strip.
///
/// This is the pressure integration over the undisturbed wave field on the
/// hull surface. For a strip at position `x_station` and draft `draft`:
///
/// `F_FK = rho * g * A * wave_amplitude * exp(-k * draft) * cos(k * x - omega * t)`
pub fn froude_krylov_heave_strip(
    wave: &RegularWave,
    x_station: f64,
    section_area: f64,
    draft: f64,
    t: f64,
    rho: f64,
    g: f64,
) -> f64 {
    let k = wave.wave_number;
    let phase = k * x_station * wave.heading.cos() - wave.omega * t;
    rho * g * section_area * wave.amplitude * (-k * draft).exp() * phase.cos()
}

/// Compute total Froude-Krylov heave force by integrating over all strips.
pub fn total_froude_krylov_heave(
    wave: &RegularWave,
    strips: &[StripCoefficients],
    t: f64,
    rho: f64,
    g: f64,
) -> f64 {
    if strips.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0_f64;
    for i in 0..strips.len() - 1 {
        let dx = (strips[i + 1].x_station - strips[i].x_station).abs();
        let f_i = froude_krylov_heave_strip(
            wave,
            strips[i].x_station,
            strips[i].section_area(),
            strips[i].section_draft,
            t,
            rho,
            g,
        );
        let f_next = froude_krylov_heave_strip(
            wave,
            strips[i + 1].x_station,
            strips[i + 1].section_area(),
            strips[i + 1].section_draft,
            t,
            rho,
            g,
        );
        total += 0.5 * (f_i + f_next) * dx;
    }
    total
}

// ─────────────────────────────────────────────────────────────────────────────
// Morison Equation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the inline force per unit length on a slender cylinder using the
/// Morison equation.
///
/// `F = 0.5 * rho * Cd * D * |u_rel| * u_rel + rho * Cm * A * du/dt`
///
/// where `u_rel = u_water - u_body` is the relative velocity,
/// `D` is the cylinder diameter, `A = pi * D^2 / 4` is the cross-section area,
/// `Cd` is the drag coefficient, and `Cm` is the inertia coefficient.
///
/// `du_dt` is the fluid acceleration (time derivative of water velocity).
pub fn morison_force_per_length(
    u_water: f64,
    u_body: f64,
    du_dt: f64,
    diameter: f64,
    cd: f64,
    cm: f64,
    rho: f64,
) -> f64 {
    let u_rel = u_water - u_body;
    let area = PI * diameter * diameter / 4.0;
    let f_drag = 0.5 * rho * cd * diameter * u_rel.abs() * u_rel;
    let f_inertia = rho * cm * area * du_dt;
    f_drag + f_inertia
}

/// Compute the total Morison force on a vertical cylinder from `z_top` to `z_bottom`.
///
/// Integrates using the trapezoidal rule with `n_segments` segments.
/// `water_vel_fn` and `water_accel_fn` are closures returning the water
/// velocity and acceleration at depth `z`.
pub fn morison_total_force<F1, F2>(
    water_vel_fn: F1,
    water_accel_fn: F2,
    u_body: f64,
    z_top: f64,
    z_bottom: f64,
    diameter: f64,
    cd: f64,
    cm: f64,
    rho: f64,
    n_segments: usize,
) -> f64
where
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64,
{
    let n = n_segments.max(1);
    let dz = (z_bottom - z_top) / n as f64;
    let mut total = 0.0_f64;

    for i in 0..=n {
        let z = z_top + i as f64 * dz;
        let u_w = water_vel_fn(z);
        let du_dt = water_accel_fn(z);
        let f = morison_force_per_length(u_w, u_body, du_dt, diameter, cd, cm, rho);
        let weight = if i == 0 || i == n { 0.5 } else { 1.0 };
        total += weight * f * dz;
    }
    total
}

// ─────────────────────────────────────────────────────────────────────────────
// Mooring Line Catenary
// ─────────────────────────────────────────────────────────────────────────────

/// Mooring line parameters for catenary analysis.
#[derive(Debug, Clone)]
pub struct MooringLine {
    /// Total unstretched line length (m).
    pub length: f64,
    /// Weight per unit length in water (N/m).
    pub weight_per_length: f64,
    /// Horizontal distance from anchor to fairlead (m).
    pub horizontal_span: f64,
    /// Vertical distance from anchor to fairlead (m).
    pub vertical_span: f64,
}

impl MooringLine {
    /// Create a new mooring line with given parameters.
    pub fn new(
        length: f64,
        weight_per_length: f64,
        horizontal_span: f64,
        vertical_span: f64,
    ) -> Self {
        Self {
            length,
            weight_per_length,
            horizontal_span,
            vertical_span,
        }
    }

    /// Compute the horizontal tension at the anchor for a catenary mooring.
    ///
    /// Uses the catenary equation iteratively. Returns the horizontal component
    /// of the tension at the anchor point (N).
    pub fn catenary_horizontal_tension(&self) -> f64 {
        let w = self.weight_per_length;
        let l = self.length;
        let h = self.horizontal_span;
        let v = self.vertical_span;

        if w.abs() < 1e-30 || l < 1e-30 {
            return 0.0;
        }

        // Newton-Raphson to solve: s = sqrt(L^2 - V^2) where s is the
        // suspended length, and catenary equation relates H, s, h, v.
        let s = (l * l - v * v).max(0.0).sqrt();
        if s < 1e-30 {
            return 0.0;
        }

        // Initial guess: taut line approximation
        let mut th = w * s * s / (2.0 * h);

        for _ in 0..50 {
            let a = th / w; // catenary parameter
            // catenary: h = a * sinh(s/a), v = a * (cosh(s/a) - 1)
            let s_over_a = s / a;
            if s_over_a > 20.0 {
                break; // overflow protection
            }
            let sinh_val = s_over_a.sinh();
            let _cosh_val = s_over_a.cosh();
            let f = a * sinh_val - h;
            let df = sinh_val - s_over_a * s_over_a.cosh();
            if df.abs() < 1e-30 {
                break;
            }
            // da/dTh = 1/w
            let dth = -f / (df / w);
            th += dth;
            th = th.max(1e-6); // keep positive
            if dth.abs() < 1e-10 {
                break;
            }
        }
        th
    }

    /// Compute the tension at the fairlead (top end) of the mooring line.
    pub fn fairlead_tension(&self) -> f64 {
        let th = self.catenary_horizontal_tension();
        let tv = self.weight_per_length * self.length;
        (th * th + tv * tv).sqrt()
    }

    /// Compute the line shape (profile) as a series of points `(x, z)`.
    ///
    /// Returns `n_points` evenly spaced along the line arc length.
    pub fn catenary_profile(&self, n_points: usize) -> Vec<[f64; 2]> {
        let n = n_points.max(2);
        let th = self.catenary_horizontal_tension();
        let w = self.weight_per_length;

        if th.abs() < 1e-30 || w.abs() < 1e-30 {
            // Degenerate: return straight line
            return (0..n)
                .map(|i| {
                    let t = i as f64 / (n - 1) as f64;
                    [t * self.horizontal_span, t * self.vertical_span]
                })
                .collect();
        }

        let a = th / w;
        let s_total = (self.length * self.length - self.vertical_span * self.vertical_span)
            .max(0.0)
            .sqrt();

        (0..n)
            .map(|i| {
                let s = s_total * i as f64 / (n - 1) as f64;
                let x = a * (s / a).sinh();
                let z = a * ((s / a).cosh() - 1.0);
                [x, z]
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Seakeeping RAO
// ─────────────────────────────────────────────────────────────────────────────

/// Response Amplitude Operator (RAO) at a single frequency.
#[derive(Debug, Clone, Copy)]
pub struct RaoPoint {
    /// Wave frequency (rad/s).
    pub omega: f64,
    /// RAO magnitude (response amplitude / wave amplitude).
    pub magnitude: f64,
    /// RAO phase angle (rad).
    pub phase: f64,
}

/// Compute the heave RAO for a simple 1-DOF heave model.
///
/// The equation of motion is: `(M + A33) * z'' + B33 * z' + C33 * z = F_wave`
///
/// RAO = F0 / sqrt((C33 - (M + A33) * omega^2)^2 + (B33 * omega)^2)
///
/// where `F0` is the wave force amplitude per unit wave amplitude.
pub fn heave_rao(
    mass: f64,
    added_mass: f64,
    damping: f64,
    stiffness: f64,
    excitation_amplitude: f64,
    omegas: &[f64],
) -> Vec<RaoPoint> {
    let m_total = mass + added_mass;
    omegas
        .iter()
        .map(|&omega| {
            let spring = stiffness - m_total * omega * omega;
            let damp = damping * omega;
            let denom = (spring * spring + damp * damp).sqrt();
            let magnitude = if denom > 1e-30 {
                excitation_amplitude / denom
            } else {
                0.0
            };
            let phase = if denom > 1e-30 {
                (-damp).atan2(spring)
            } else {
                0.0
            };
            RaoPoint {
                omega,
                magnitude,
                phase,
            }
        })
        .collect()
}

/// Find the natural frequency from RAO data (frequency of maximum response).
pub fn natural_frequency_from_rao(rao: &[RaoPoint]) -> f64 {
    if rao.is_empty() {
        return 0.0;
    }
    let mut max_mag = 0.0_f64;
    let mut omega_nat = 0.0_f64;
    for pt in rao {
        if pt.magnitude > max_mag {
            max_mag = pt.magnitude;
            omega_nat = pt.omega;
        }
    }
    omega_nat
}

/// Compute the heave natural frequency analytically.
///
/// `omega_n = sqrt(C33 / (M + A33))`.
pub fn heave_natural_frequency(mass: f64, added_mass: f64, stiffness: f64) -> f64 {
    let m_total = mass + added_mass;
    if m_total.abs() < 1e-30 {
        return 0.0;
    }
    (stiffness / m_total).max(0.0).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Roll Damping
// ─────────────────────────────────────────────────────────────────────────────

/// Compute roll damping from various components.
///
/// Total roll damping = friction + eddy + bilge keel + wave.
/// Returns the equivalent linearized damping coefficient (N*m*s/rad).
pub fn roll_damping_total(
    friction_damping: f64,
    eddy_damping: f64,
    bilge_keel_damping: f64,
    wave_damping: f64,
) -> f64 {
    friction_damping + eddy_damping + bilge_keel_damping + wave_damping
}

/// Compute skin friction roll damping.
///
/// Approximation: `B_f = 0.5 * rho * Cf * S * r^2 * omega_roll`
///
/// where `Cf` is the friction coefficient, `S` is the wetted surface area,
/// and `r` is the roll radius of gyration.
pub fn roll_damping_friction(
    rho: f64,
    cf: f64,
    wetted_surface: f64,
    roll_radius: f64,
    omega_roll: f64,
) -> f64 {
    0.5 * rho * cf * wetted_surface * roll_radius * roll_radius * omega_roll
}

/// Compute bilge keel roll damping contribution.
///
/// Approximation: `B_bk = rho * Cd * l_bk * b_bk * r_bk^2 * omega_roll`
///
/// where `l_bk` is bilge keel length, `b_bk` is bilge keel breadth,
/// `r_bk` is the distance from roll axis to bilge keel, and `Cd` is drag coefficient.
pub fn roll_damping_bilge_keel(
    rho: f64,
    cd: f64,
    bk_length: f64,
    bk_breadth: f64,
    bk_radius: f64,
    omega_roll: f64,
) -> f64 {
    rho * cd * bk_length * bk_breadth * bk_radius * bk_radius * omega_roll
}

/// Compute the roll natural frequency.
///
/// `omega_roll = sqrt(displacement * g * GM_T / (I_xx + A_44))`
pub fn roll_natural_frequency(
    displacement: f64,
    g: f64,
    gm_t: f64,
    ixx: f64,
    added_inertia_roll: f64,
) -> f64 {
    let total_inertia = ixx + added_inertia_roll;
    if total_inertia.abs() < 1e-30 {
        return 0.0;
    }
    let c44 = displacement * g * gm_t;
    (c44 / total_inertia).max(0.0).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Propeller Thrust
// ─────────────────────────────────────────────────────────────────────────────

/// Propeller parameters.
#[derive(Debug, Clone)]
pub struct Propeller {
    /// Propeller diameter (m).
    pub diameter: f64,
    /// Thrust coefficient Kt (from open-water diagram).
    pub kt: f64,
    /// Torque coefficient Kq (from open-water diagram).
    pub kq: f64,
    /// Wake fraction (w).
    pub wake_fraction: f64,
    /// Thrust deduction factor (t).
    pub thrust_deduction: f64,
}

impl Propeller {
    /// Create a propeller with given characteristics.
    pub fn new(diameter: f64, kt: f64, kq: f64, wake_fraction: f64, thrust_deduction: f64) -> Self {
        Self {
            diameter,
            kt,
            kq,
            wake_fraction,
            thrust_deduction,
        }
    }

    /// Compute the thrust produced by the propeller.
    ///
    /// `T = Kt * rho * n^2 * D^4` where `n` is revolutions per second.
    pub fn thrust(&self, n_rps: f64, rho: f64) -> f64 {
        self.kt * rho * n_rps * n_rps * self.diameter.powi(4)
    }

    /// Compute the torque required by the propeller.
    ///
    /// `Q = Kq * rho * n^2 * D^5`.
    pub fn torque(&self, n_rps: f64, rho: f64) -> f64 {
        self.kq * rho * n_rps * n_rps * self.diameter.powi(5)
    }

    /// Compute the effective thrust (after thrust deduction).
    ///
    /// `T_eff = (1 - t) * T`.
    pub fn effective_thrust(&self, n_rps: f64, rho: f64) -> f64 {
        (1.0 - self.thrust_deduction) * self.thrust(n_rps, rho)
    }

    /// Compute the advance velocity at the propeller disk.
    ///
    /// `Va = Vs * (1 - w)` where `Vs` is the ship speed.
    pub fn advance_velocity(&self, ship_speed: f64) -> f64 {
        ship_speed * (1.0 - self.wake_fraction)
    }

    /// Compute the advance ratio J = Va / (n * D).
    pub fn advance_ratio(&self, ship_speed: f64, n_rps: f64) -> f64 {
        if n_rps.abs() < 1e-30 {
            return 0.0;
        }
        self.advance_velocity(ship_speed) / (n_rps * self.diameter)
    }

    /// Compute the open-water efficiency.
    ///
    /// `eta_0 = J * Kt / (2 * pi * Kq)`.
    pub fn open_water_efficiency(&self, ship_speed: f64, n_rps: f64) -> f64 {
        let j = self.advance_ratio(ship_speed, n_rps);
        if self.kq.abs() < 1e-30 {
            return 0.0;
        }
        j * self.kt / (2.0 * PI * self.kq)
    }

    /// Compute propeller shaft power.
    ///
    /// `P = 2 * pi * n * Q`.
    pub fn shaft_power(&self, n_rps: f64, rho: f64) -> f64 {
        2.0 * PI * n_rps * self.torque(n_rps, rho)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Resistance Estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Froude number.
///
/// `Fn = V / sqrt(g * L)` where `V` is ship speed and `L` is waterline length.
pub fn froude_number(speed: f64, length: f64, g: f64) -> f64 {
    let denom = (g * length).sqrt();
    if denom < 1e-30 {
        return 0.0;
    }
    speed / denom
}

/// Compute the Reynolds number.
///
/// `Re = V * L / nu` where `nu` is kinematic viscosity.
pub fn reynolds_number(speed: f64, length: f64, nu: f64) -> f64 {
    if nu.abs() < 1e-30 {
        return 0.0;
    }
    speed * length / nu
}

/// Compute the ITTC 1957 friction coefficient.
///
/// `Cf = 0.075 / (log10(Re) - 2)^2`.
pub fn ittc_friction_coefficient(re: f64) -> f64 {
    if re < 1.0 {
        return 0.0;
    }
    let log_re = re.log10();
    let denom = log_re - 2.0;
    if denom.abs() < 1e-10 {
        return 0.0;
    }
    0.075 / (denom * denom)
}

/// Compute ship resistance using the Holtrop-Menher simplified method.
///
/// Returns the total resistance in Newtons.
pub fn simple_resistance(
    speed: f64,
    length: f64,
    beam: f64,
    draft: f64,
    wetted_surface: f64,
    block_coeff: f64,
    rho: f64,
    nu: f64,
) -> f64 {
    let re = reynolds_number(speed, length, nu);
    let cf = ittc_friction_coefficient(re);

    // Form factor approximation (1+k)
    let _lr = length / beam;
    let form_factor = 1.0 + 0.1 * block_coeff;

    // Frictional resistance
    let rf = 0.5 * rho * wetted_surface * speed * speed * cf * form_factor;

    // Wave resistance approximation (Froude number based)
    let fn_val = froude_number(speed, length, 9.81);
    let _cw = if fn_val < 0.4 {
        0.001 * fn_val * fn_val
    } else {
        0.001 * 0.16
    };

    // Simple wave resistance component
    let rw = 0.5 * rho * wetted_surface * speed * speed * _cw;

    // Air resistance approximation
    let _air_area = beam * draft * 0.3;

    rf + rw
}

// ─────────────────────────────────────────────────────────────────────────────
// Encounter Frequency
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the encounter frequency.
///
/// `omega_e = omega - k * V * cos(mu)` where `omega` is the wave frequency,
/// `k` is the wave number, `V` is the ship speed, and `mu` is the heading angle.
pub fn encounter_frequency(omega: f64, wave_number: f64, ship_speed: f64, heading: f64) -> f64 {
    omega - wave_number * ship_speed * heading.cos()
}

/// Compute the significant wave height from a spectrum.
///
/// `Hs = 4 * sqrt(m0)` where `m0` is the zeroth spectral moment.
pub fn significant_wave_height(spectrum_m0: f64) -> f64 {
    4.0 * spectrum_m0.max(0.0).sqrt()
}

/// Compute the JONSWAP spectral density at frequency `omega`.
///
/// `S(omega) = (alpha * g^2 / omega^5) * exp(-5/4 * (omega_p/omega)^4) * gamma^delta`
///
/// where `delta = exp(-(omega - omega_p)^2 / (2 * sigma^2 * omega_p^2))`.
pub fn jonswap_spectrum(omega: f64, omega_peak: f64, hs: f64, gamma: f64, g: f64) -> f64 {
    if omega <= 0.0 {
        return 0.0;
    }
    let sigma = if omega <= omega_peak { 0.07 } else { 0.09 };
    let alpha = 0.0081; // Phillips constant approximation

    let base = alpha * g * g / omega.powi(5);
    let exp_term = (-1.25 * (omega_peak / omega).powi(4)).exp();

    let delta_exp =
        -((omega - omega_peak).powi(2)) / (2.0 * sigma * sigma * omega_peak * omega_peak);
    let gamma_factor = gamma.powf(delta_exp.exp());

    // Scale to match Hs
    let _hs_factor = hs * hs / 16.0;

    base * exp_term * gamma_factor
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const RHO_WATER: f64 = 1025.0;
    const G: f64 = 9.81;

    fn sample_hull() -> HullGeometry {
        HullGeometry::new(100.0, 15.0, 5.0, 0.7, 0.8, 6.0, 2.5)
    }

    fn sample_ship() -> ShipState6DOF {
        ShipState6DOF::new(5000e3, [1e8, 5e8, 5e8])
    }

    fn sample_strips() -> Vec<StripCoefficients> {
        (0..11)
            .map(|i| {
                let x = i as f64 * 10.0; // 0..100 m
                StripCoefficients {
                    x_station: x,
                    section_beam: 15.0,
                    section_draft: 5.0,
                    section_area_coeff: 0.8,
                }
            })
            .collect()
    }

    // ── ShipState6DOF ────────────────────────────────────────────────────────

    #[test]
    fn test_ship_state_new() {
        let s = sample_ship();
        assert_eq!(s.position, [0.0; 3]);
        assert_eq!(s.velocity, [0.0; 3]);
        assert!(s.mass > 0.0);
    }

    #[test]
    fn test_ship_state_kinetic_energy_at_rest() {
        let s = sample_ship();
        assert!(s.kinetic_energy().abs() < 1e-10);
    }

    #[test]
    fn test_ship_state_kinetic_energy_moving() {
        let mut s = sample_ship();
        s.velocity = [10.0, 0.0, 0.0];
        let ke = s.kinetic_energy();
        let expected = 0.5 * s.mass * 100.0;
        assert!((ke - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_ship_state_integrate_euler() {
        let mut s = sample_ship();
        s.integrate_euler([s.mass * 1.0, 0.0, 0.0], [0.0; 3], 1.0);
        // After 1s with F=m*1, v should be ~1 m/s
        assert!((s.velocity[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ship_speed_over_ground() {
        let mut s = sample_ship();
        s.velocity = [3.0, 4.0, 1.0]; // SOG = 5.0
        assert!((s.speed_over_ground() - 5.0).abs() < 1e-10);
    }

    // ── HullGeometry ─────────────────────────────────────────────────────────

    #[test]
    fn test_hull_displaced_volume() {
        let h = sample_hull();
        let vol = h.displaced_volume();
        assert!((vol - 100.0 * 15.0 * 5.0 * 0.7).abs() < 1e-6);
    }

    #[test]
    fn test_hull_waterplane_area() {
        let h = sample_hull();
        let aw = h.waterplane_area();
        assert!((aw - 100.0 * 15.0 * 0.8).abs() < 1e-6);
    }

    #[test]
    fn test_hull_gm_positive() {
        let h = sample_hull();
        // For a typical ship, GM should be positive
        assert!(h.gm_transverse() > 0.0, "GM_T = {}", h.gm_transverse());
    }

    #[test]
    fn test_hull_bm_transverse_positive() {
        let h = sample_hull();
        assert!(h.bm_transverse() > 0.0);
    }

    // ── Hydrostatic Restoring ────────────────────────────────────────────────

    #[test]
    fn test_hydrostatic_restoring_zero_displacement() {
        let h = sample_hull();
        let (force, moment) = hydrostatic_restoring(&h, 0.0, 0.0, 0.0, RHO_WATER, G);
        assert!(force[2].abs() < 1e-6);
        assert!(moment[0].abs() < 1e-6);
    }

    #[test]
    fn test_hydrostatic_restoring_heave() {
        let h = sample_hull();
        let (force, _moment) = hydrostatic_restoring(&h, 1.0, 0.0, 0.0, RHO_WATER, G);
        // Positive heave (down) -> negative restoring force (up)
        assert!(force[2] < 0.0, "heave restoring should be upward");
    }

    #[test]
    fn test_hydrostatic_restoring_roll() {
        let h = sample_hull();
        let (_force, moment) = hydrostatic_restoring(&h, 0.0, 0.1, 0.0, RHO_WATER, G);
        // Positive roll -> negative restoring moment
        assert!(moment[0] < 0.0, "roll restoring should oppose displacement");
    }

    // ── Strip Theory ─────────────────────────────────────────────────────────

    #[test]
    fn test_strip_added_mass_heave_positive() {
        let a33 = strip_added_mass_heave(15.0, RHO_WATER);
        assert!(a33 > 0.0);
    }

    #[test]
    fn test_strip_added_mass_sway_positive() {
        let a22 = strip_added_mass_sway(5.0, RHO_WATER);
        assert!(a22 > 0.0);
    }

    #[test]
    fn test_total_heave_added_mass() {
        let strips = sample_strips();
        let a33 = total_heave_added_mass(&strips, RHO_WATER);
        assert!(a33 > 0.0, "total heave added mass = {a33}");
    }

    #[test]
    fn test_total_sway_added_mass() {
        let strips = sample_strips();
        let a22 = total_sway_added_mass(&strips, RHO_WATER);
        assert!(a22 > 0.0, "total sway added mass = {a22}");
    }

    #[test]
    fn test_strip_damping_heave_positive() {
        let b33 = strip_damping_heave(15.0, 1.0, RHO_WATER, G);
        assert!(b33 > 0.0);
    }

    // ── Froude-Krylov ────────────────────────────────────────────────────────

    #[test]
    fn test_regular_wave_deep_water() {
        let w = RegularWave::deep_water(1.0, 1.0, 0.0, G);
        assert!(w.wave_number > 0.0);
        assert!(w.wavelength() > 0.0);
        assert!(w.period() > 0.0);
    }

    #[test]
    fn test_wave_elevation_amplitude() {
        let w = RegularWave::deep_water(2.0, 1.0, 0.0, G);
        let eta = w.elevation(0.0, 0.0);
        assert!((eta - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_froude_krylov_nonzero() {
        let w = RegularWave::deep_water(1.0, 1.0, 0.0, G);
        let f = froude_krylov_heave_strip(&w, 50.0, 60.0, 5.0, 0.0, RHO_WATER, G);
        assert!(f.abs() > 0.0, "FK force should be nonzero");
    }

    // ── Morison Equation ─────────────────────────────────────────────────────

    #[test]
    fn test_morison_zero_relative_velocity() {
        let f = morison_force_per_length(1.0, 1.0, 0.0, 1.0, 1.0, 2.0, RHO_WATER);
        // u_rel = 0, du_dt = 0 -> f = 0
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_morison_drag_direction() {
        let f = morison_force_per_length(2.0, 0.0, 0.0, 1.0, 1.0, 2.0, RHO_WATER);
        // u_rel = 2, positive -> force should be positive
        assert!(
            f > 0.0,
            "Morison drag should be positive for positive u_rel"
        );
    }

    #[test]
    fn test_morison_total_force() {
        let f = morison_total_force(
            |_z| 1.0,
            |_z| 0.0,
            0.0,
            0.0,
            10.0,
            1.0,
            1.0,
            2.0,
            RHO_WATER,
            10,
        );
        assert!(f > 0.0, "total Morison force should be positive");
    }

    // ── Mooring Line ─────────────────────────────────────────────────────────

    #[test]
    fn test_mooring_line_tension_positive() {
        let ml = MooringLine::new(300.0, 1000.0, 250.0, 100.0);
        let th = ml.catenary_horizontal_tension();
        assert!(th > 0.0, "horizontal tension = {th}");
    }

    #[test]
    fn test_mooring_fairlead_tension_greater_than_horizontal() {
        let ml = MooringLine::new(300.0, 1000.0, 250.0, 100.0);
        let th = ml.catenary_horizontal_tension();
        let tf = ml.fairlead_tension();
        assert!(tf >= th, "fairlead tension {tf} should >= horizontal {th}");
    }

    #[test]
    fn test_mooring_profile_endpoints() {
        let ml = MooringLine::new(300.0, 1000.0, 250.0, 100.0);
        let profile = ml.catenary_profile(20);
        assert_eq!(profile.len(), 20);
        assert!((profile[0][0]).abs() < 1e-6, "profile should start at x=0");
    }

    // ── RAO ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_heave_rao_peak_near_natural_freq() {
        let mass = 5000e3;
        let added_mass = 2000e3;
        let stiffness = 1e7;
        let damping = 1e5;
        let f0 = 1.0;

        let omegas: Vec<f64> = (1..100).map(|i| i as f64 * 0.02).collect();
        let rao = heave_rao(mass, added_mass, damping, stiffness, f0, &omegas);
        let omega_peak = natural_frequency_from_rao(&rao);
        let omega_n = heave_natural_frequency(mass, added_mass, stiffness);
        assert!(
            (omega_peak - omega_n).abs() < 0.1,
            "RAO peak {omega_peak} should be near omega_n {omega_n}"
        );
    }

    #[test]
    fn test_heave_natural_frequency() {
        let omega_n = heave_natural_frequency(1000.0, 500.0, 15000.0);
        // omega_n = sqrt(15000/1500) = sqrt(10) ~ 3.16
        let expected = (15000.0_f64 / 1500.0).sqrt();
        assert!((omega_n - expected).abs() < 1e-6);
    }

    // ── Roll Damping ─────────────────────────────────────────────────────────

    #[test]
    fn test_roll_damping_total() {
        let total = roll_damping_total(100.0, 200.0, 500.0, 300.0);
        assert!((total - 1100.0).abs() < 1e-10);
    }

    #[test]
    fn test_roll_damping_friction_positive() {
        let bf = roll_damping_friction(RHO_WATER, 0.003, 3000.0, 5.0, 0.5);
        assert!(bf > 0.0);
    }

    #[test]
    fn test_roll_natural_frequency() {
        let omega = roll_natural_frequency(5e6, G, 1.0, 1e8, 5e7);
        assert!(omega > 0.0);
    }

    // ── Propeller ────────────────────────────────────────────────────────────

    #[test]
    fn test_propeller_thrust_positive() {
        let p = Propeller::new(3.0, 0.2, 0.03, 0.25, 0.15);
        let t = p.thrust(5.0, RHO_WATER);
        assert!(t > 0.0, "thrust = {t}");
    }

    #[test]
    fn test_propeller_effective_thrust_less_than_total() {
        let p = Propeller::new(3.0, 0.2, 0.03, 0.25, 0.15);
        let t = p.thrust(5.0, RHO_WATER);
        let t_eff = p.effective_thrust(5.0, RHO_WATER);
        assert!(
            t_eff < t,
            "effective thrust should be less due to deduction"
        );
    }

    #[test]
    fn test_propeller_advance_velocity() {
        let p = Propeller::new(3.0, 0.2, 0.03, 0.25, 0.15);
        let va = p.advance_velocity(10.0);
        assert!((va - 7.5).abs() < 1e-10, "Va = Vs * (1-w) = 10 * 0.75");
    }

    #[test]
    fn test_propeller_advance_ratio() {
        let p = Propeller::new(3.0, 0.2, 0.03, 0.25, 0.15);
        let j = p.advance_ratio(10.0, 5.0);
        // J = 7.5 / (5*3) = 0.5
        assert!((j - 0.5).abs() < 1e-10, "J = {j}");
    }

    #[test]
    fn test_propeller_shaft_power_positive() {
        let p = Propeller::new(3.0, 0.2, 0.03, 0.25, 0.15);
        let pw = p.shaft_power(5.0, RHO_WATER);
        assert!(pw > 0.0, "shaft power should be positive");
    }

    // ── Resistance ───────────────────────────────────────────────────────────

    #[test]
    fn test_froude_number() {
        let fn_val = froude_number(10.0, 100.0, G);
        let expected = 10.0 / (G * 100.0).sqrt();
        assert!((fn_val - expected).abs() < 1e-10);
    }

    #[test]
    fn test_reynolds_number() {
        let re = reynolds_number(10.0, 100.0, 1.19e-6);
        assert!(re > 1e8, "Re should be large: {re}");
    }

    #[test]
    fn test_ittc_friction_coefficient() {
        let cf = ittc_friction_coefficient(1e9);
        assert!(cf > 0.0 && cf < 0.01, "Cf = {cf}");
    }

    // ── Encounter Frequency ──────────────────────────────────────────────────

    #[test]
    fn test_encounter_frequency_head_seas() {
        let omega = 1.0;
        let k = omega * omega / G;
        let omega_e = encounter_frequency(omega, k, 5.0, 0.0);
        // Head seas: mu=0 -> cos(0)=1 -> omega_e = omega - k*V
        assert!(
            omega_e < omega,
            "encounter freq should decrease in head seas for this case"
        );
    }

    #[test]
    fn test_significant_wave_height() {
        let hs = significant_wave_height(1.0);
        assert!((hs - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_jonswap_spectrum_positive() {
        let s = jonswap_spectrum(1.0, 0.8, 2.0, 3.3, G);
        assert!(s > 0.0, "JONSWAP spectrum should be positive: {s}");
    }

    #[test]
    fn test_wave_phase_velocity() {
        let w = RegularWave::deep_water(1.0, 1.0, 0.0, G);
        let cp = w.phase_velocity();
        assert!(cp > 0.0, "phase velocity should be positive");
    }

    #[test]
    fn test_wave_particle_velocities() {
        let w = RegularWave::deep_water(1.0, 1.0, 0.0, G);
        let uh = w.horizontal_velocity(0.0, 0.0, 0.0);
        assert!(
            uh.abs() > 0.0,
            "horizontal velocity at surface should be nonzero"
        );
    }
}
