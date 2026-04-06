#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Explosion and blast wave simulation using SPH.
//!
//! Provides TNT equivalent energy models, Friedlander blast wave profiles,
//! reflected pressure computation, Mach stem detection, fragment tracking,
//! JWL equation of state, shock tube (Sod) validation, blast loading on
//! structures, cratering models, air blast propagation, particle splitting,
//! and damage assessment.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Euclidean distance between two 3-D points.
#[inline]
fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Norm of a 3-D vector.
#[inline]
fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Subtract two 3-D vectors: `a - b`.
#[inline]
fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-D vectors.
#[inline]
fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-D vector.
#[inline]
fn scale3(s: f64, v: &[f64; 3]) -> [f64; 3] {
    [s * v[0], s * v[1], s * v[2]]
}

/// Normalize a 3-D vector. Returns the length.
#[inline]
fn normalize3(v: &mut [f64; 3]) -> f64 {
    let len = norm3(v);
    if len > f64::EPSILON {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
    len
}

// ---------------------------------------------------------------------------
// 1. TNT equivalent energy model
// ---------------------------------------------------------------------------

/// TNT detonation heat (J/kg).
pub const TNT_HEAT_OF_DETONATION: f64 = 4.184e6;

/// Ambient atmospheric pressure (Pa).
pub const AMBIENT_PRESSURE: f64 = 101_325.0;

/// Ambient air density (kg/m^3).
pub const AMBIENT_DENSITY: f64 = 1.225;

/// Speed of sound in air at standard conditions (m/s).
pub const SOUND_SPEED_AIR: f64 = 343.0;

/// Ratio of specific heats for air.
pub const GAMMA_AIR: f64 = 1.4;

/// Compute the TNT equivalent mass for an arbitrary explosive.
///
/// `energy_per_kg` is the specific energy (J/kg) of the explosive.
/// `mass` is the charge mass (kg).
pub fn tnt_equivalent(mass: f64, energy_per_kg: f64) -> f64 {
    mass * energy_per_kg / TNT_HEAT_OF_DETONATION
}

/// Hopkinson–Cranz scaled distance: `Z = R / W^{1/3}`.
///
/// `range` is in metres, `tnt_mass` in kg.
pub fn scaled_distance(range: f64, tnt_mass: f64) -> f64 {
    let w_cbrt = tnt_mass.cbrt().max(f64::EPSILON);
    range / w_cbrt
}

/// Peak incident overpressure from the Kingery–Bulmash polynomial fit
/// (simplified, valid for `Z` in `[0.5, 40]` m/kg^{1/3}).
///
/// Returns overpressure in Pa.
pub fn kingery_bulmash_overpressure(scaled_dist: f64) -> f64 {
    let z = scaled_dist.max(0.5);
    // Simplified Kingery–Bulmash: P_s = 808 * (1 + (z/4.5)^2) / sqrt(1+(z/0.048)^2)
    //                                / sqrt(1+(z/0.32)^2) / sqrt(1+(z/1.35)^2)  [kPa]
    let num = 808.0 * (1.0 + (z / 4.5).powi(2));
    let d1 = (1.0 + (z / 0.048).powi(2)).sqrt();
    let d2 = (1.0 + (z / 0.32).powi(2)).sqrt();
    let d3 = (1.0 + (z / 1.35).powi(2)).sqrt();
    let ps_kpa = num / (d1 * d2 * d3);
    ps_kpa * 1000.0 // convert to Pa
}

/// Peak overpressure at a given range from a TNT charge.
pub fn peak_overpressure(range: f64, tnt_mass: f64) -> f64 {
    let z = scaled_distance(range, tnt_mass);
    kingery_bulmash_overpressure(z)
}

// ---------------------------------------------------------------------------
// 2. Friedlander blast wave profile
// ---------------------------------------------------------------------------

/// Friedlander waveform parameters.
#[derive(Debug, Clone, Copy)]
pub struct FriedlanderParams {
    /// Peak overpressure (Pa).
    pub peak_pressure: f64,
    /// Positive phase duration (s).
    pub positive_duration: f64,
    /// Waveform decay coefficient (dimensionless).
    pub decay_coeff: f64,
    /// Arrival time (s).
    pub arrival_time: f64,
}

impl FriedlanderParams {
    /// Create from TNT equivalent charge at given range.
    pub fn from_blast(tnt_mass: f64, range: f64) -> Self {
        let z = scaled_distance(range, tnt_mass);
        let ps = kingery_bulmash_overpressure(z);
        // Positive phase duration: t_d ≈ 1.0 * W^{1/3} * Z^{0.5} ms (simplified)
        let w_cbrt = tnt_mass.cbrt().max(f64::EPSILON);
        let td = 1.0e-3 * w_cbrt * z.sqrt().max(0.1);
        // Arrival time: R / c (approximate)
        let ta = range / SOUND_SPEED_AIR;
        FriedlanderParams {
            peak_pressure: ps,
            positive_duration: td,
            decay_coeff: 1.5, // typical value
            arrival_time: ta,
        }
    }

    /// Evaluate the Friedlander overpressure at time `t`.
    ///
    /// `P(t) = P_s * (1 - tau) * exp(-b * tau)`, where `tau = (t - t_a) / t_d`.
    pub fn pressure(&self, t: f64) -> f64 {
        let tau = (t - self.arrival_time) / self.positive_duration;
        if !(0.0..=1.0).contains(&tau) {
            return 0.0;
        }
        self.peak_pressure * (1.0 - tau) * (-self.decay_coeff * tau).exp()
    }

    /// Positive impulse: integral of pressure over the positive phase.
    ///
    /// Computed analytically:
    /// `I = P_s * t_d * [1/b - (1/b^2) * (1 - exp(-b))]`
    pub fn positive_impulse(&self) -> f64 {
        let b = self.decay_coeff;
        let ps = self.peak_pressure;
        let td = self.positive_duration;
        if b.abs() < 1e-12 {
            return 0.5 * ps * td;
        }
        ps * td * (1.0 / b - (1.0 - (-b).exp()) / (b * b))
    }
}

// ---------------------------------------------------------------------------
// 3. Reflected pressure computation
// ---------------------------------------------------------------------------

/// Compute the reflected overpressure for a normal reflection.
///
/// Uses the Rankine–Hugoniot relation:
/// `P_r = 2 P_s * (7 P0 + 4 P_s) / (7 P0 + P_s)`
pub fn reflected_pressure(incident_pressure: f64) -> f64 {
    let p0 = AMBIENT_PRESSURE;
    2.0 * incident_pressure * (7.0 * p0 + 4.0 * incident_pressure) / (7.0 * p0 + incident_pressure)
}

/// Reflected pressure at an angle of incidence `alpha` (radians).
///
/// Uses the simplified relation with reflection coefficient.
pub fn oblique_reflected_pressure(incident_pressure: f64, alpha: f64) -> f64 {
    let pr_normal = reflected_pressure(incident_pressure);
    let ps = incident_pressure;
    // Interpolation between incident and reflected based on angle
    let cos_a = alpha.cos();
    let cr = pr_normal / ps.max(f64::EPSILON);
    let cr_oblique = 1.0 + (cr - 1.0) * cos_a * cos_a;
    ps * cr_oblique
}

/// Dynamic pressure behind the shock.
///
/// `q = 0.5 * rho * u^2 = 5 P_s^2 / (2 * (7 P0 + P_s))`
pub fn dynamic_pressure(incident_pressure: f64) -> f64 {
    let p0 = AMBIENT_PRESSURE;
    5.0 * incident_pressure.powi(2) / (2.0 * (7.0 * p0 + incident_pressure))
}

// ---------------------------------------------------------------------------
// 4. Mach stem formation detection
// ---------------------------------------------------------------------------

/// Configuration for Mach stem detection.
#[derive(Debug, Clone)]
pub struct MachStemDetector {
    /// Ground level z-coordinate.
    pub ground_z: f64,
    /// Blast center position.
    pub center: [f64; 3],
    /// Height of burst above ground.
    pub hob: f64,
    /// TNT equivalent mass (kg).
    pub tnt_mass: f64,
}

impl MachStemDetector {
    /// Create a new detector for an airburst.
    pub fn new(center: [f64; 3], tnt_mass: f64) -> Self {
        let hob = center[2];
        Self {
            ground_z: 0.0,
            center,
            hob,
            tnt_mass,
        }
    }

    /// Critical angle for Mach reflection (simplified).
    ///
    /// For strong shocks, Mach reflection occurs when the incidence angle
    /// exceeds ~40 degrees.
    pub fn critical_angle(&self) -> f64 {
        // Approximate: angle depends on shock strength
        let z = scaled_distance(self.hob, self.tnt_mass);
        let base = 39.0_f64.to_radians();
        let correction = (z - 1.0).max(0.0) * 0.5_f64.to_radians();
        (base + correction).min(50.0_f64.to_radians())
    }

    /// Compute the ground range at which the Mach stem begins to form.
    pub fn mach_stem_onset_range(&self) -> f64 {
        let alpha_c = self.critical_angle();
        self.hob * alpha_c.tan()
    }

    /// Check whether a ground point at range `r` from ground zero is in
    /// the Mach stem region.
    pub fn is_mach_stem(&self, ground_range: f64) -> bool {
        ground_range >= self.mach_stem_onset_range()
    }

    /// Estimate the Mach stem height at ground range `r`.
    ///
    /// Simplified linear growth from onset.
    pub fn mach_stem_height(&self, ground_range: f64) -> f64 {
        let onset = self.mach_stem_onset_range();
        if ground_range <= onset {
            return 0.0;
        }
        let growth_rate = 0.1; // simplified
        (ground_range - onset) * growth_rate
    }
}

// ---------------------------------------------------------------------------
// 5. Fragment tracking with ballistic trajectories
// ---------------------------------------------------------------------------

/// A single fragment with mass, position, and velocity.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Mass (kg).
    pub mass: f64,
    /// Position in 3-D (m).
    pub position: [f64; 3],
    /// Velocity in 3-D (m/s).
    pub velocity: [f64; 3],
    /// Drag coefficient (dimensionless).
    pub drag_coeff: f64,
    /// Cross-sectional area (m^2).
    pub cross_area: f64,
    /// Is the fragment still in flight?
    pub active: bool,
}

impl Fragment {
    /// Create a new fragment.
    pub fn new(
        mass: f64,
        position: [f64; 3],
        velocity: [f64; 3],
        drag_coeff: f64,
        cross_area: f64,
    ) -> Self {
        Self {
            mass,
            position,
            velocity,
            drag_coeff,
            cross_area,
            active: true,
        }
    }

    /// Advance the fragment by `dt` seconds under gravity and drag.
    pub fn step(&mut self, dt: f64) {
        if !self.active {
            return;
        }
        let speed = norm3(&self.velocity).max(f64::EPSILON);
        // Drag force: F_d = 0.5 * rho * Cd * A * v^2
        let drag_mag = 0.5 * AMBIENT_DENSITY * self.drag_coeff * self.cross_area * speed * speed;
        let drag_accel = drag_mag / self.mass;
        // Drag opposes motion
        let drag_dir = scale3(-1.0 / speed, &self.velocity);
        let drag_a = scale3(drag_accel, &drag_dir);
        // Gravity
        let gravity = [0.0, 0.0, -9.81];
        let accel = add3(&drag_a, &gravity);

        // Update velocity and position (Euler)
        self.velocity = add3(&self.velocity, &scale3(dt, &accel));
        self.position = add3(&self.position, &scale3(dt, &self.velocity));

        // Deactivate if below ground
        if self.position[2] < 0.0 {
            self.position[2] = 0.0;
            self.active = false;
        }
    }

    /// Kinetic energy of the fragment (J).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = dot3(&self.velocity, &self.velocity);
        0.5 * self.mass * v2
    }

    /// Speed (m/s).
    pub fn speed(&self) -> f64 {
        norm3(&self.velocity)
    }
}

/// Initial fragment velocity from the Gurney equation (cylindrical case).
///
/// `V = sqrt(2 E_g) * (M/C + 0.5)^{-0.5}`
///
/// `gurney_energy` in J/kg, `metal_mass` and `charge_mass` in kg.
pub fn gurney_velocity_cylinder(gurney_energy: f64, metal_mass: f64, charge_mass: f64) -> f64 {
    let mc_ratio = metal_mass / charge_mass.max(f64::EPSILON);
    (2.0 * gurney_energy).sqrt() * (mc_ratio + 0.5).powf(-0.5)
}

/// Initial fragment velocity from the Gurney equation (spherical case).
///
/// `V = sqrt(2 E_g) * (M/C + 0.6)^{-0.5}`
pub fn gurney_velocity_sphere(gurney_energy: f64, metal_mass: f64, charge_mass: f64) -> f64 {
    let mc_ratio = metal_mass / charge_mass.max(f64::EPSILON);
    (2.0 * gurney_energy).sqrt() * (mc_ratio + 0.6).powf(-0.5)
}

/// Generate random fragment velocities around the origin.
pub fn generate_fragments(
    n_fragments: usize,
    mass_each: f64,
    speed: f64,
    drag_coeff: f64,
    cross_area: f64,
    center: [f64; 3],
) -> Vec<Fragment> {
    let mut frags = Vec::with_capacity(n_fragments);
    let mut rng = rand::rng();

    use rand::RngExt;
    for _ in 0..n_fragments {
        // Random direction on unit sphere
        let theta: f64 = rng.random_range(0.0..PI);
        let phi: f64 = rng.random_range(0.0..(2.0 * PI));
        let vx = speed * theta.sin() * phi.cos();
        let vy = speed * theta.sin() * phi.sin();
        let vz = speed * theta.cos();
        frags.push(Fragment::new(
            mass_each,
            center,
            [vx, vy, vz],
            drag_coeff,
            cross_area,
        ));
    }
    frags
}

// ---------------------------------------------------------------------------
// 6. Jones–Wilkins–Lee (JWL) equation of state
// ---------------------------------------------------------------------------

/// JWL equation of state parameters for detonation products.
///
/// `P = A * (1 - omega / (R1 * V)) * exp(-R1 * V)
///    + B * (1 - omega / (R2 * V)) * exp(-R2 * V)
///    + omega * e / V`
///
/// where `V = rho_0 / rho` is the relative volume and `e` is the
/// specific internal energy.
#[derive(Debug, Clone, Copy)]
pub struct JwlEos {
    /// First pressure coefficient (Pa).
    pub a: f64,
    /// Second pressure coefficient (Pa).
    pub b: f64,
    /// First exponential coefficient.
    pub r1: f64,
    /// Second exponential coefficient.
    pub r2: f64,
    /// Gruneisen coefficient.
    pub omega: f64,
    /// Reference density (kg/m^3).
    pub rho_0: f64,
}

impl JwlEos {
    /// Standard TNT JWL parameters.
    pub fn tnt() -> Self {
        Self {
            a: 3.738e11,
            b: 3.747e9,
            r1: 4.15,
            r2: 0.9,
            omega: 0.35,
            rho_0: 1630.0,
        }
    }

    /// Compute pressure from density and specific internal energy.
    pub fn pressure(&self, rho: f64, energy: f64) -> f64 {
        let v = self.rho_0 / rho.max(f64::EPSILON);
        let term1 = self.a * (1.0 - self.omega / (self.r1 * v)) * (-self.r1 * v).exp();
        let term2 = self.b * (1.0 - self.omega / (self.r2 * v)) * (-self.r2 * v).exp();
        let term3 = self.omega * energy * rho;
        term1 + term2 + term3
    }

    /// Compute the sound speed (approximate).
    ///
    /// `c^2 = dP/drho |_e  ≈ (P(rho+eps) - P(rho-eps)) / (2*eps)`
    pub fn sound_speed(&self, rho: f64, energy: f64) -> f64 {
        let eps = rho * 1e-6;
        let pp = self.pressure(rho + eps, energy);
        let pm = self.pressure(rho - eps, energy);
        let dp_drho = (pp - pm) / (2.0 * eps);
        dp_drho.abs().sqrt()
    }

    /// Detonation Chapman–Jouguet pressure.
    pub fn cj_pressure(&self) -> f64 {
        // Approximate: use rho = rho_0 and typical CJ energy
        let cj_energy = TNT_HEAT_OF_DETONATION * 0.3; // fraction of total
        self.pressure(self.rho_0 * 1.3, cj_energy)
    }
}

/// Ideal gas equation of state.
///
/// `P = (gamma - 1) * rho * e`
#[derive(Debug, Clone, Copy)]
pub struct IdealGasEos {
    /// Ratio of specific heats.
    pub gamma: f64,
}

impl IdealGasEos {
    /// Standard air EOS.
    pub fn air() -> Self {
        Self { gamma: GAMMA_AIR }
    }

    /// Compute pressure.
    pub fn pressure(&self, rho: f64, energy: f64) -> f64 {
        (self.gamma - 1.0) * rho * energy
    }

    /// Compute sound speed.
    pub fn sound_speed(&self, rho: f64, pressure: f64) -> f64 {
        (self.gamma * pressure / rho.max(f64::EPSILON)).sqrt()
    }

    /// Compute specific internal energy from pressure and density.
    pub fn energy(&self, rho: f64, pressure: f64) -> f64 {
        pressure / ((self.gamma - 1.0) * rho.max(f64::EPSILON))
    }
}

// ---------------------------------------------------------------------------
// 7. Shock tube problem (Sod) for validation
// ---------------------------------------------------------------------------

/// Initial conditions for one side of a shock tube.
#[derive(Debug, Clone, Copy)]
pub struct ShockTubeState {
    /// Density (kg/m^3).
    pub rho: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Velocity (m/s).
    pub velocity: f64,
}

/// Sod shock tube problem definition.
#[derive(Debug, Clone)]
pub struct SodShockTube {
    /// Left state.
    pub left: ShockTubeState,
    /// Right state.
    pub right: ShockTubeState,
    /// Ratio of specific heats.
    pub gamma: f64,
    /// Domain length (m).
    pub length: f64,
    /// Diaphragm position (m).
    pub diaphragm: f64,
}

impl SodShockTube {
    /// Standard Sod problem.
    pub fn standard() -> Self {
        Self {
            left: ShockTubeState {
                rho: 1.0,
                pressure: 1.0,
                velocity: 0.0,
            },
            right: ShockTubeState {
                rho: 0.125,
                pressure: 0.1,
                velocity: 0.0,
            },
            gamma: 1.4,
            length: 1.0,
            diaphragm: 0.5,
        }
    }

    /// Compute the exact solution at time `t` and position `x`.
    ///
    /// Returns `(density, velocity, pressure)`.
    pub fn exact_solution(&self, x: f64, t: f64) -> (f64, f64, f64) {
        if t < f64::EPSILON {
            if x < self.diaphragm {
                return (self.left.rho, self.left.velocity, self.left.pressure);
            } else {
                return (self.right.rho, self.right.velocity, self.right.pressure);
            }
        }

        let g = self.gamma;
        let gm1 = g - 1.0;
        let gp1 = g + 1.0;

        let pl = self.left.pressure;
        let pr = self.right.pressure;
        let rl = self.left.rho;
        let rr = self.right.rho;
        let _ul = self.left.velocity;
        let _ur = self.right.velocity;

        let cl = (g * pl / rl).sqrt();
        let cr = (g * pr / rr).sqrt();

        // Iterative solver for the star-region pressure
        let mut p_star = 0.5 * (pl + pr);
        for _ in 0..50 {
            let fl = if p_star > pl {
                (p_star - pl) * (2.0 / (gp1 * rl) / (p_star + gm1 / gp1 * pl)).sqrt()
            } else {
                2.0 * cl / gm1 * ((p_star / pl).powf(gm1 / (2.0 * g)) - 1.0)
            };
            let fr = if p_star > pr {
                (p_star - pr) * (2.0 / (gp1 * rr) / (p_star + gm1 / gp1 * pr)).sqrt()
            } else {
                2.0 * cr / gm1 * ((p_star / pr).powf(gm1 / (2.0 * g)) - 1.0)
            };
            let dfl = if p_star > pl {
                let a = 2.0 / (gp1 * rl);
                let b_val = gm1 / gp1 * pl;
                (a / (b_val + p_star)).sqrt() * (1.0 - (p_star - pl) / (2.0 * (b_val + p_star)))
            } else {
                1.0 / (rl * cl) * (p_star / pl).powf(-(gp1) / (2.0 * g))
            };
            let dfr = if p_star > pr {
                let a = 2.0 / (gp1 * rr);
                let b_val = gm1 / gp1 * pr;
                (a / (b_val + p_star)).sqrt() * (1.0 - (p_star - pr) / (2.0 * (b_val + p_star)))
            } else {
                1.0 / (rr * cr) * (p_star / pr).powf(-(gp1) / (2.0 * g))
            };

            let residual = fl + fr;
            let deriv = dfl + dfr;
            if deriv.abs() < 1e-15 {
                break;
            }
            let dp = -residual / deriv;
            p_star = (p_star + dp).max(1e-10);
            if dp.abs() / p_star < 1e-10 {
                break;
            }
        }

        // Star-region velocity
        let fl_star = if p_star > pl {
            (p_star - pl) * (2.0 / (gp1 * rl) / (p_star + gm1 / gp1 * pl)).sqrt()
        } else {
            2.0 * cl / gm1 * ((p_star / pl).powf(gm1 / (2.0 * g)) - 1.0)
        };
        let fr_star = if p_star > pr {
            (p_star - pr) * (2.0 / (gp1 * rr) / (p_star + gm1 / gp1 * pr)).sqrt()
        } else {
            2.0 * cr / gm1 * ((p_star / pr).powf(gm1 / (2.0 * g)) - 1.0)
        };
        let u_star = 0.5 * (fl_star - fr_star);

        // Sample solution at x/t
        let s = (x - self.diaphragm) / t;

        if s < u_star {
            // Left of contact
            if p_star > pl {
                // Left shock
                let sl = -cl * ((gp1 / (2.0 * g)) * (p_star / pl) + gm1 / (2.0 * g)).sqrt();
                if s < sl {
                    (rl, 0.0, pl)
                } else {
                    let rho_star_l =
                        rl * (p_star / pl + gm1 / gp1) / (gm1 / gp1 * p_star / pl + 1.0);
                    (rho_star_l, u_star, p_star)
                }
            } else {
                // Left rarefaction
                let sh = -cl;
                let st = u_star - cl * (p_star / pl).powf(gm1 / (2.0 * g));
                if s < sh {
                    (rl, 0.0, pl)
                } else if s < st {
                    let cs = 2.0 / gp1 * (cl + gm1 / 2.0 * s.abs().copysign(-s));
                    // Simplified rarefaction fan
                    let rho_fan = rl * (cs / cl).powf(2.0 / gm1);
                    let p_fan = pl * (cs / cl).powf(2.0 * g / gm1);
                    let u_fan = 2.0 / gp1 * (cl + s);
                    (rho_fan, u_fan, p_fan)
                } else {
                    let rho_star_l = rl * (p_star / pl).powf(1.0 / g);
                    (rho_star_l, u_star, p_star)
                }
            }
        } else {
            // Right of contact
            if p_star > pr {
                // Right shock
                let sr = cr * ((gp1 / (2.0 * g)) * (p_star / pr) + gm1 / (2.0 * g)).sqrt();
                if s > sr {
                    (rr, 0.0, pr)
                } else {
                    let rho_star_r =
                        rr * (p_star / pr + gm1 / gp1) / (gm1 / gp1 * p_star / pr + 1.0);
                    (rho_star_r, u_star, p_star)
                }
            } else {
                // Right rarefaction
                let sh = cr;
                let st = u_star + cr * (p_star / pr).powf(gm1 / (2.0 * g));
                if s > sh {
                    (rr, 0.0, pr)
                } else if s > st {
                    let cs = 2.0 / gp1 * (cr - gm1 / 2.0 * s);
                    let rho_fan = rr * (cs / cr).powf(2.0 / gm1);
                    let p_fan = pr * (cs / cr).powf(2.0 * g / gm1);
                    let u_fan = 2.0 / gp1 * (-cr + s);
                    (rho_fan, u_fan, p_fan)
                } else {
                    let rho_star_r = rr * (p_star / pr).powf(1.0 / g);
                    (rho_star_r, u_star, p_star)
                }
            }
        }
    }

    /// Generate initial SPH particles for the shock tube.
    ///
    /// Returns `(positions, densities, pressures, velocities)` as 1-D arrays.
    pub fn generate_particles(
        &self,
        n_particles: usize,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let dx = self.length / n_particles as f64;
        let mut positions = Vec::with_capacity(n_particles);
        let mut densities = Vec::with_capacity(n_particles);
        let mut pressures = Vec::with_capacity(n_particles);
        let mut velocities = Vec::with_capacity(n_particles);

        for i in 0..n_particles {
            let x = (i as f64 + 0.5) * dx;
            positions.push(x);
            if x < self.diaphragm {
                densities.push(self.left.rho);
                pressures.push(self.left.pressure);
                velocities.push(self.left.velocity);
            } else {
                densities.push(self.right.rho);
                pressures.push(self.right.pressure);
                velocities.push(self.right.velocity);
            }
        }

        (positions, densities, pressures, velocities)
    }
}

// ---------------------------------------------------------------------------
// 8. Blast loading on structures (pressure–impulse)
// ---------------------------------------------------------------------------

/// Pressure–impulse (P–I) diagram point.
#[derive(Debug, Clone, Copy)]
pub struct PressureImpulse {
    /// Peak pressure (Pa).
    pub pressure: f64,
    /// Impulse (Pa.s).
    pub impulse: f64,
}

/// Compute the impulse on a flat wall at range `r` from a TNT charge.
///
/// Uses the Friedlander profile integrated analytically.
pub fn blast_impulse_on_wall(tnt_mass: f64, range: f64) -> PressureImpulse {
    let ps = peak_overpressure(range, tnt_mass);
    let pr = reflected_pressure(ps);
    let fw = FriedlanderParams::from_blast(tnt_mass, range);
    // Impulse with reflected pressure
    let b = fw.decay_coeff;
    let td = fw.positive_duration;
    let impulse = if b.abs() < 1e-12 {
        0.5 * pr * td
    } else {
        pr * td * (1.0 / b - (1.0 - (-b).exp()) / (b * b))
    };
    PressureImpulse {
        pressure: pr,
        impulse,
    }
}

/// Structural response category from overpressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DamageLevel {
    /// No structural damage.
    None,
    /// Light damage (window breakage).
    Light,
    /// Moderate damage (structural cracking).
    Moderate,
    /// Heavy damage (partial collapse).
    Heavy,
    /// Total destruction.
    Destroyed,
}

/// Classify structural damage from peak overpressure.
pub fn classify_damage(overpressure_pa: f64) -> DamageLevel {
    let p_kpa = overpressure_pa / 1000.0;
    if p_kpa < 3.5 {
        DamageLevel::None
    } else if p_kpa < 14.0 {
        DamageLevel::Light
    } else if p_kpa < 35.0 {
        DamageLevel::Moderate
    } else if p_kpa < 82.0 {
        DamageLevel::Heavy
    } else {
        DamageLevel::Destroyed
    }
}

/// Compute the standoff distance for a given damage level from TNT charge.
///
/// Iterates to find range where overpressure matches threshold.
pub fn standoff_distance(tnt_mass: f64, target_level: DamageLevel) -> f64 {
    let threshold_kpa = match target_level {
        DamageLevel::None => 3.5,
        DamageLevel::Light => 14.0,
        DamageLevel::Moderate => 35.0,
        DamageLevel::Heavy => 82.0,
        DamageLevel::Destroyed => 200.0,
    };
    let threshold = threshold_kpa * 1000.0;
    // Binary search for range
    let mut lo = 0.1;
    let mut hi = 1000.0;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let p = peak_overpressure(mid, tnt_mass);
        if p > threshold {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ---------------------------------------------------------------------------
// 9. Cratering model
// ---------------------------------------------------------------------------

/// Crater dimensions from an explosion.
#[derive(Debug, Clone, Copy)]
pub struct CraterDimensions {
    /// Crater radius (m).
    pub radius: f64,
    /// Crater depth (m).
    pub depth: f64,
    /// Volume of ejected material (m^3).
    pub volume: f64,
}

/// Empirical cratering model (surface burst).
///
/// Radius ~ 0.46 * W^{1/3} for TNT on moderate soil.
/// Depth ~ 0.35 * Radius.
pub fn crater_dimensions(tnt_mass: f64) -> CraterDimensions {
    let w_cbrt = tnt_mass.cbrt();
    let radius = 0.46 * w_cbrt;
    let depth = 0.35 * radius;
    let volume = PI / 6.0 * (3.0 * radius * radius * depth + depth.powi(3));
    CraterDimensions {
        radius,
        depth,
        volume,
    }
}

/// Cratering model for buried charges.
///
/// `depth_of_burial` in metres.
pub fn buried_crater(tnt_mass: f64, depth_of_burial: f64) -> CraterDimensions {
    let w_cbrt = tnt_mass.cbrt();
    let optimal_dob = 0.3 * w_cbrt;
    let scale = if depth_of_burial < optimal_dob {
        1.0 + 0.5 * depth_of_burial / optimal_dob
    } else {
        1.5 * (-((depth_of_burial - optimal_dob) / optimal_dob).powi(2)).exp()
    };
    let base = crater_dimensions(tnt_mass);
    CraterDimensions {
        radius: base.radius * scale,
        depth: base.depth * scale,
        volume: base.volume * scale * scale * scale,
    }
}

// ---------------------------------------------------------------------------
// 10. Air blast propagation (SPH particles)
// ---------------------------------------------------------------------------

/// An SPH particle for blast simulation.
#[derive(Debug, Clone)]
pub struct BlastParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Density (kg/m^3).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Specific internal energy (J/kg).
    pub energy: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub smoothing_length: f64,
    /// Is this particle part of the detonation products?
    pub is_detonation: bool,
}

/// Cubic spline kernel in 3-D.
///
/// `W(r, h) = sigma / h^3 * f(q)`, where `q = r / h`.
pub fn cubic_spline_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Gradient of the cubic spline kernel (magnitude) in 3-D.
pub fn cubic_spline_grad_3d(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h.powi(4));
    if q < f64::EPSILON {
        0.0
    } else if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        sigma * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    }
}

/// Density summation for a particle.
pub fn sph_density(particle_idx: usize, particles: &[BlastParticle], neighbors: &[usize]) -> f64 {
    let pi = &particles[particle_idx];
    let mut rho = 0.0;
    for &j in neighbors {
        let pj = &particles[j];
        let r = dist3(&pi.position, &pj.position);
        let w = cubic_spline_3d(r, pi.smoothing_length);
        rho += pj.mass * w;
    }
    rho
}

/// SPH momentum equation acceleration (Euler equation with artificial viscosity).
#[allow(clippy::too_many_arguments)]
pub fn sph_acceleration(
    particle_idx: usize,
    particles: &[BlastParticle],
    neighbors: &[usize],
    alpha_visc: f64,
    beta_visc: f64,
) -> [f64; 3] {
    let pi = &particles[particle_idx];
    let mut acc = [0.0; 3];
    let rhoi = pi.density.max(f64::EPSILON);
    let pi_p = pi.pressure;

    for &j in neighbors {
        let pj = &particles[j];
        let rhoj = pj.density.max(f64::EPSILON);
        let rij = sub3(&pi.position, &pj.position);
        let r = norm3(&rij).max(f64::EPSILON);
        let vij = sub3(&pi.velocity, &pj.velocity);
        let vr = dot3(&vij, &rij);

        // Artificial viscosity (Monaghan)
        let mu = if vr < 0.0 {
            let h_avg = 0.5 * (pi.smoothing_length + pj.smoothing_length);
            let c_avg =
                0.5 * ((GAMMA_AIR * pi_p / rhoi).sqrt() + (GAMMA_AIR * pj.pressure / rhoj).sqrt());
            let mu_val = h_avg * vr / (r * r + 0.01 * h_avg * h_avg);
            (-alpha_visc * c_avg * mu_val + beta_visc * mu_val * mu_val) / (0.5 * (rhoi + rhoj))
        } else {
            0.0
        };

        let grad_w = cubic_spline_grad_3d(r, pi.smoothing_length);
        let factor = pj.mass * (pi_p / (rhoi * rhoi) + pj.pressure / (rhoj * rhoj) + mu);
        let grad_dir = scale3(1.0 / r, &rij);
        acc[0] -= factor * grad_w * grad_dir[0];
        acc[1] -= factor * grad_w * grad_dir[1];
        acc[2] -= factor * grad_w * grad_dir[2];
    }
    acc
}

/// SPH energy equation (rate of change of specific internal energy).
#[allow(clippy::too_many_arguments)]
pub fn sph_energy_rate(
    particle_idx: usize,
    particles: &[BlastParticle],
    neighbors: &[usize],
    alpha_visc: f64,
    _beta_visc: f64,
) -> f64 {
    let pi = &particles[particle_idx];
    let rhoi = pi.density.max(f64::EPSILON);
    let pi_p = pi.pressure;
    let mut de = 0.0;

    for &j in neighbors {
        let pj = &particles[j];
        let rhoj = pj.density.max(f64::EPSILON);
        let rij = sub3(&pi.position, &pj.position);
        let r = norm3(&rij).max(f64::EPSILON);
        let vij = sub3(&pi.velocity, &pj.velocity);
        let vr = dot3(&vij, &rij);

        let h_avg = 0.5 * (pi.smoothing_length + pj.smoothing_length);
        let c_avg =
            0.5 * ((GAMMA_AIR * pi_p / rhoi).sqrt() + (GAMMA_AIR * pj.pressure / rhoj).sqrt());
        let visc = if vr < 0.0 {
            let mu_val = h_avg * vr / (r * r + 0.01 * h_avg * h_avg);
            alpha_visc * c_avg * mu_val.abs() / (0.5 * (rhoi + rhoj))
        } else {
            0.0
        };

        let grad_w = cubic_spline_grad_3d(r, pi.smoothing_length);
        let factor = pj.mass * (pi_p / (rhoi * rhoi) + 0.5 * visc);
        let grad_dir = scale3(1.0 / r, &rij);
        de += factor * grad_w * dot3(&vij, &grad_dir);
    }
    de
}

// ---------------------------------------------------------------------------
// 11. Particle splitting for blast front resolution
// ---------------------------------------------------------------------------

/// Split a particle into `n_children` daughter particles arranged in
/// a ring around the parent position.
pub fn split_particle(parent: &BlastParticle, n_children: usize) -> Vec<BlastParticle> {
    let offset = parent.smoothing_length * 0.25;
    let child_mass = parent.mass / n_children as f64;
    let child_h = parent.smoothing_length * (1.0 / n_children as f64).cbrt();
    let mut children = Vec::with_capacity(n_children);

    for k in 0..n_children {
        let angle = 2.0 * PI * k as f64 / n_children as f64;
        let dx = offset * angle.cos();
        let dy = offset * angle.sin();
        let pos = add3(&parent.position, &[dx, dy, 0.0]);
        children.push(BlastParticle {
            position: pos,
            velocity: parent.velocity,
            density: parent.density,
            pressure: parent.pressure,
            energy: parent.energy,
            mass: child_mass,
            smoothing_length: child_h,
            is_detonation: parent.is_detonation,
        });
    }
    children
}

/// Determine which particles should be split based on pressure gradient.
pub fn particles_to_split(
    particles: &[BlastParticle],
    pressure_gradient_threshold: f64,
) -> Vec<usize> {
    let mut to_split = Vec::new();
    let n = particles.len();
    if n < 2 {
        return to_split;
    }

    for i in 0..n {
        // Estimate pressure gradient from nearest neighbors (simplified)
        let pi = &particles[i];
        let mut max_grad = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let pj = &particles[j];
            let r = dist3(&pi.position, &pj.position);
            if r < 2.0 * pi.smoothing_length && r > f64::EPSILON {
                let grad = (pi.pressure - pj.pressure).abs() / r;
                if grad > max_grad {
                    max_grad = grad;
                }
            }
        }
        if max_grad > pressure_gradient_threshold {
            to_split.push(i);
        }
    }
    to_split
}

// ---------------------------------------------------------------------------
// 12. Blast simulation driver
// ---------------------------------------------------------------------------

/// Configuration for a blast SPH simulation.
#[derive(Debug, Clone)]
pub struct BlastSimConfig {
    /// TNT equivalent mass (kg).
    pub tnt_mass: f64,
    /// Charge center position.
    pub center: [f64; 3],
    /// Domain half-size.
    pub domain_size: f64,
    /// Initial particle spacing.
    pub dx: f64,
    /// Artificial viscosity alpha.
    pub alpha_visc: f64,
    /// Artificial viscosity beta.
    pub beta_visc: f64,
    /// CFL number for time stepping.
    pub cfl: f64,
    /// End time.
    pub end_time: f64,
}

impl BlastSimConfig {
    /// Default configuration for a 1 kg TNT charge.
    pub fn default_1kg() -> Self {
        Self {
            tnt_mass: 1.0,
            center: [0.0, 0.0, 0.0],
            domain_size: 5.0,
            dx: 0.1,
            alpha_visc: 1.0,
            beta_visc: 2.0,
            cfl: 0.3,
            end_time: 0.01,
        }
    }

    /// Estimate the initial smoothing length.
    pub fn smoothing_length(&self) -> f64 {
        1.3 * self.dx
    }

    /// Estimate the initial time step.
    pub fn initial_dt(&self) -> f64 {
        self.cfl * self.dx / SOUND_SPEED_AIR
    }
}

/// State of the blast simulation at a given time.
#[derive(Debug, Clone)]
pub struct BlastSimState {
    /// Current simulation time.
    pub time: f64,
    /// Particles.
    pub particles: Vec<BlastParticle>,
    /// Time step counter.
    pub step: usize,
}

impl BlastSimState {
    /// Create the initial state from configuration.
    pub fn initialize(config: &BlastSimConfig) -> Self {
        let h = config.smoothing_length();
        let mut particles = Vec::new();

        // Simple 1-D line of particles along x-axis for testing
        let n = (2.0 * config.domain_size / config.dx).round() as usize;
        let mass_per_particle = AMBIENT_DENSITY * config.dx;
        let ambient_energy = AMBIENT_PRESSURE / ((GAMMA_AIR - 1.0) * AMBIENT_DENSITY);

        for i in 0..n {
            let x = -config.domain_size + (i as f64 + 0.5) * config.dx;
            let r = (x - config.center[0]).abs();
            // Detonation region: within charge radius
            let charge_radius = 0.05 * config.tnt_mass.cbrt();
            let (rho, pressure, energy, is_det) = if r < charge_radius {
                let jwl = JwlEos::tnt();
                let det_rho = jwl.rho_0;
                let det_energy = TNT_HEAT_OF_DETONATION;
                let det_p = jwl.pressure(det_rho, det_energy);
                (det_rho, det_p, det_energy, true)
            } else {
                (AMBIENT_DENSITY, AMBIENT_PRESSURE, ambient_energy, false)
            };

            particles.push(BlastParticle {
                position: [x, 0.0, 0.0],
                velocity: [0.0, 0.0, 0.0],
                density: rho,
                pressure,
                energy,
                mass: mass_per_particle,
                smoothing_length: h,
                is_detonation: is_det,
            });
        }

        Self {
            time: 0.0,
            particles,
            step: 0,
        }
    }

    /// Number of active particles.
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Maximum pressure in the domain.
    pub fn max_pressure(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.pressure)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Maximum density in the domain.
    pub fn max_density(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.density)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

// ---------------------------------------------------------------------------
// 13. Overpressure damage criteria
// ---------------------------------------------------------------------------

/// Human injury thresholds from blast overpressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InjuryLevel {
    /// No injury.
    None,
    /// Threshold for eardrum rupture (~35 kPa).
    EardrumRupture,
    /// Lung damage threshold (~100 kPa).
    LungDamage,
    /// Lethal (~200 kPa).
    Lethal,
}

/// Classify human injury from peak overpressure.
pub fn classify_injury(overpressure_pa: f64) -> InjuryLevel {
    let p_kpa = overpressure_pa / 1000.0;
    if p_kpa < 35.0 {
        InjuryLevel::None
    } else if p_kpa < 100.0 {
        InjuryLevel::EardrumRupture
    } else if p_kpa < 200.0 {
        InjuryLevel::LungDamage
    } else {
        InjuryLevel::Lethal
    }
}

/// Probit function for lethality from overpressure.
///
/// `Pr = a + b * ln(P_s)` where `P_s` is in Pa.
pub fn lethality_probit(overpressure_pa: f64) -> f64 {
    let a = -77.1;
    let b = 6.91;
    a + b * overpressure_pa.max(1.0).ln()
}

/// Convert probit to probability.
pub fn probit_to_probability(probit: f64) -> f64 {
    // Approximation of the normal CDF
    let x = probit - 5.0;
    0.5 * (1.0 + (x / 2.0_f64.sqrt()).tanh())
}

// ---------------------------------------------------------------------------
// 14. Blast wave arrival time and positive phase duration
// ---------------------------------------------------------------------------

/// Estimate blast wave arrival time at range `r` from a TNT charge.
///
/// Uses the simplified model: `t_a ≈ R / c * (1 + 0.1 * (R/R_0)^{-1})`
pub fn arrival_time(range: f64, tnt_mass: f64) -> f64 {
    let r0 = tnt_mass.cbrt();
    let factor = 1.0 + 0.1 * r0 / range.max(f64::EPSILON);
    range / SOUND_SPEED_AIR * factor
}

/// Estimate positive phase duration.
///
/// `t_d ≈ W^{1/3} * f(Z)`, simplified.
pub fn positive_phase_duration(range: f64, tnt_mass: f64) -> f64 {
    let z = scaled_distance(range, tnt_mass);
    let w_cbrt = tnt_mass.cbrt();
    // Simplified fit: t_d ≈ 1e-3 * W^{1/3} * (1 + z)^0.5
    1.0e-3 * w_cbrt * (1.0 + z).sqrt()
}

/// Compute the negative phase parameters (simplified).
///
/// Returns `(peak negative pressure, negative duration)`.
pub fn negative_phase(range: f64, tnt_mass: f64) -> (f64, f64) {
    let ps = peak_overpressure(range, tnt_mass);
    let td = positive_phase_duration(range, tnt_mass);
    // Negative phase peak is typically 10-20% of positive
    let pn = 0.15 * ps;
    let tn = 1.5 * td; // negative phase is longer
    (pn, tn)
}

// ---------------------------------------------------------------------------
// 15. Blast energy partitioning
// ---------------------------------------------------------------------------

/// Compute the fraction of TNT energy partitioned into blast wave,
/// ground shock, and cratering.
///
/// Returns `(blast_fraction, ground_shock_fraction, crater_fraction)`.
pub fn energy_partition(height_of_burst: f64, tnt_mass: f64) -> (f64, f64, f64) {
    let r0 = tnt_mass.cbrt();
    let h_ratio = height_of_burst / r0.max(f64::EPSILON);
    if h_ratio > 3.0 {
        // Free air burst: all energy to blast
        (0.85, 0.10, 0.05)
    } else if h_ratio > 0.5 {
        // Near-surface
        (0.60, 0.25, 0.15)
    } else {
        // Surface or buried
        (0.40, 0.35, 0.25)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    #[test]
    fn test_tnt_equivalent_identity() {
        let eq = tnt_equivalent(1.0, TNT_HEAT_OF_DETONATION);
        assert!((eq - 1.0).abs() < TOL);
    }

    #[test]
    fn test_tnt_equivalent_double_energy() {
        let eq = tnt_equivalent(1.0, 2.0 * TNT_HEAT_OF_DETONATION);
        assert!((eq - 2.0).abs() < TOL);
    }

    #[test]
    fn test_scaled_distance() {
        let z = scaled_distance(10.0, 1.0);
        assert!((z - 10.0).abs() < TOL);
        let z8 = scaled_distance(10.0, 8.0);
        assert!((z8 - 5.0).abs() < TOL);
    }

    #[test]
    fn test_kingery_bulmash_positive() {
        let ps = kingery_bulmash_overpressure(1.0);
        assert!(ps > 0.0);
    }

    #[test]
    fn test_kingery_bulmash_decreasing() {
        let ps1 = kingery_bulmash_overpressure(1.0);
        let ps2 = kingery_bulmash_overpressure(5.0);
        assert!(ps1 > ps2, "Overpressure should decrease with distance");
    }

    #[test]
    fn test_friedlander_peak() {
        let fw = FriedlanderParams {
            peak_pressure: 100_000.0,
            positive_duration: 0.005,
            decay_coeff: 1.5,
            arrival_time: 0.01,
        };
        let p0 = fw.pressure(0.01); // at arrival
        assert!((p0 - 100_000.0).abs() < TOL);
    }

    #[test]
    fn test_friedlander_zero_before_arrival() {
        let fw = FriedlanderParams {
            peak_pressure: 100_000.0,
            positive_duration: 0.005,
            decay_coeff: 1.5,
            arrival_time: 0.01,
        };
        assert!(fw.pressure(0.005).abs() < TOL);
    }

    #[test]
    fn test_friedlander_zero_after_positive() {
        let fw = FriedlanderParams {
            peak_pressure: 100_000.0,
            positive_duration: 0.005,
            decay_coeff: 1.5,
            arrival_time: 0.01,
        };
        assert!(fw.pressure(0.02).abs() < TOL);
    }

    #[test]
    fn test_friedlander_impulse_positive() {
        let fw = FriedlanderParams::from_blast(1.0, 5.0);
        let impulse = fw.positive_impulse();
        assert!(impulse > 0.0);
    }

    #[test]
    fn test_reflected_pressure_gt_incident() {
        let ps = 100_000.0;
        let pr = reflected_pressure(ps);
        assert!(pr > ps);
    }

    #[test]
    fn test_reflected_pressure_factor() {
        // For very small overpressure, reflection factor -> 2
        let ps = 1.0; // very small
        let pr = reflected_pressure(ps);
        let factor = pr / ps;
        assert!((factor - 2.0).abs() < 0.1);
    }

    #[test]
    fn test_oblique_reflection_normal() {
        let ps = 100_000.0;
        let pr_norm = reflected_pressure(ps);
        let pr_obl = oblique_reflected_pressure(ps, 0.0); // alpha = 0 => normal
        assert!((pr_norm - pr_obl).abs() < TOL);
    }

    #[test]
    fn test_dynamic_pressure_positive() {
        let q = dynamic_pressure(100_000.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_mach_stem_detector() {
        let det = MachStemDetector::new([0.0, 0.0, 10.0], 100.0);
        assert!(det.hob > 0.0);
        let onset = det.mach_stem_onset_range();
        assert!(onset > 0.0);
        assert!(!det.is_mach_stem(0.0));
        assert!(det.is_mach_stem(onset + 1.0));
    }

    #[test]
    fn test_mach_stem_height() {
        let det = MachStemDetector::new([0.0, 0.0, 10.0], 100.0);
        let onset = det.mach_stem_onset_range();
        assert!(det.mach_stem_height(onset - 1.0).abs() < TOL);
        assert!(det.mach_stem_height(onset + 10.0) > 0.0);
    }

    #[test]
    fn test_fragment_creation() {
        let frag = Fragment::new(0.01, [0.0, 0.0, 10.0], [100.0, 0.0, 50.0], 0.5, 1e-4);
        assert!(frag.active);
        assert!(frag.kinetic_energy() > 0.0);
    }

    #[test]
    fn test_fragment_step() {
        let mut frag = Fragment::new(0.01, [0.0, 0.0, 100.0], [100.0, 0.0, 0.0], 0.5, 1e-4);
        frag.step(0.01);
        assert!(frag.position[0] > 0.0);
        assert!(frag.active);
    }

    #[test]
    fn test_fragment_ground_impact() {
        let mut frag = Fragment::new(0.01, [0.0, 0.0, 1.0], [0.0, 0.0, -100.0], 0.5, 1e-4);
        for _ in 0..100 {
            frag.step(0.01);
            if !frag.active {
                break;
            }
        }
        assert!(!frag.active);
        assert!(frag.position[2] >= 0.0 - TOL);
    }

    #[test]
    fn test_gurney_velocity_cylinder() {
        let v = gurney_velocity_cylinder(2.7e6, 1.0, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_gurney_velocity_sphere_gt_cylinder() {
        let vc = gurney_velocity_cylinder(2.7e6, 1.0, 1.0);
        let vs = gurney_velocity_sphere(2.7e6, 1.0, 1.0);
        // Gurney formula: sphere uses 0.6 correction vs cylinder's 0.5,
        // so sphere velocity is lower for same M/C ratio.
        assert!(
            vc > vs,
            "Cylinder should give higher velocity than sphere for same M/C"
        );
    }

    #[test]
    fn test_generate_fragments() {
        let frags = generate_fragments(10, 0.01, 500.0, 0.5, 1e-4, [0.0, 0.0, 0.0]);
        assert_eq!(frags.len(), 10);
        for f in &frags {
            assert!(f.active);
            assert!((f.speed() - 500.0).abs() < 1.0);
        }
    }

    #[test]
    fn test_jwl_eos_tnt() {
        let jwl = JwlEos::tnt();
        let p = jwl.pressure(jwl.rho_0, TNT_HEAT_OF_DETONATION);
        assert!(p > 0.0, "JWL pressure should be positive");
    }

    #[test]
    fn test_jwl_sound_speed_positive() {
        let jwl = JwlEos::tnt();
        let c = jwl.sound_speed(jwl.rho_0, TNT_HEAT_OF_DETONATION);
        assert!(c > 0.0);
    }

    #[test]
    fn test_ideal_gas_eos() {
        let eos = IdealGasEos::air();
        let p = eos.pressure(1.225, 2.5e5);
        assert!(p > 0.0);
        let c = eos.sound_speed(1.225, AMBIENT_PRESSURE);
        assert!((c - SOUND_SPEED_AIR).abs() < 10.0);
    }

    #[test]
    fn test_sod_standard_initial() {
        let sod = SodShockTube::standard();
        let (r, u, p) = sod.exact_solution(0.25, 0.0);
        assert!((r - 1.0).abs() < TOL);
        assert!(u.abs() < TOL);
        assert!((p - 1.0).abs() < TOL);
    }

    #[test]
    fn test_sod_right_state() {
        let sod = SodShockTube::standard();
        let (r, _u, p) = sod.exact_solution(0.75, 0.0);
        assert!((r - 0.125).abs() < TOL);
        assert!((p - 0.1).abs() < TOL);
    }

    #[test]
    fn test_sod_exact_at_time() {
        let sod = SodShockTube::standard();
        let (r, _u, p) = sod.exact_solution(0.25, 0.2);
        // Left of diaphragm at time 0.2 should still be approximately left state
        assert!(r > 0.1);
        assert!(p > 0.01);
    }

    #[test]
    fn test_sod_generate_particles() {
        let sod = SodShockTube::standard();
        let (pos, rho, _pres, _vel) = sod.generate_particles(100);
        assert_eq!(pos.len(), 100);
        assert_eq!(rho.len(), 100);
        // Left side should have density 1.0
        assert!((rho[0] - 1.0).abs() < TOL);
        // Right side should have density 0.125
        assert!((rho[99] - 0.125).abs() < TOL);
    }

    #[test]
    fn test_blast_impulse_positive() {
        let pi = blast_impulse_on_wall(1.0, 5.0);
        assert!(pi.pressure > 0.0);
        assert!(pi.impulse > 0.0);
    }

    #[test]
    fn test_classify_damage() {
        assert_eq!(classify_damage(1000.0), DamageLevel::None);
        assert_eq!(classify_damage(10_000.0), DamageLevel::Light);
        assert_eq!(classify_damage(25_000.0), DamageLevel::Moderate);
        assert_eq!(classify_damage(200_000.0), DamageLevel::Destroyed);
    }

    #[test]
    fn test_standoff_distance_positive() {
        let d = standoff_distance(100.0, DamageLevel::Light);
        assert!(d > 0.0);
    }

    #[test]
    fn test_crater_dimensions_positive() {
        let c = crater_dimensions(10.0);
        assert!(c.radius > 0.0);
        assert!(c.depth > 0.0);
        assert!(c.volume > 0.0);
    }

    #[test]
    fn test_crater_scaling() {
        let c1 = crater_dimensions(1.0);
        let c10 = crater_dimensions(10.0);
        assert!(c10.radius > c1.radius);
    }

    #[test]
    fn test_buried_crater() {
        let c = buried_crater(10.0, 0.5);
        assert!(c.radius > 0.0);
    }

    #[test]
    fn test_cubic_spline_3d_normalization() {
        // Integrate numerically: should be close to 1 in 3D
        let h = 1.0;
        let n = 100;
        let dr = 2.0 * h / n as f64;
        let mut integral = 0.0;
        for i in 0..n {
            let r = (i as f64 + 0.5) * dr;
            integral += cubic_spline_3d(r, h) * 4.0 * PI * r * r * dr;
        }
        assert!(
            (integral - 1.0).abs() < 0.05,
            "Kernel integral = {:.6}",
            integral
        );
    }

    #[test]
    fn test_split_particle() {
        let parent = BlastParticle {
            position: [0.0, 0.0, 0.0],
            velocity: [100.0, 0.0, 0.0],
            density: 1.225,
            pressure: 101325.0,
            energy: 2.5e5,
            mass: 0.1,
            smoothing_length: 0.1,
            is_detonation: false,
        };
        let children = split_particle(&parent, 4);
        assert_eq!(children.len(), 4);
        let total_mass: f64 = children.iter().map(|c| c.mass).sum();
        assert!((total_mass - parent.mass).abs() < TOL);
    }

    #[test]
    fn test_blast_sim_init() {
        let config = BlastSimConfig::default_1kg();
        let state = BlastSimState::initialize(&config);
        assert!(state.num_particles() > 0);
        assert!(state.max_pressure() > AMBIENT_PRESSURE);
    }

    #[test]
    fn test_classify_injury() {
        assert_eq!(classify_injury(10_000.0), InjuryLevel::None);
        assert_eq!(classify_injury(50_000.0), InjuryLevel::EardrumRupture);
        assert_eq!(classify_injury(150_000.0), InjuryLevel::LungDamage);
        assert_eq!(classify_injury(300_000.0), InjuryLevel::Lethal);
    }

    #[test]
    fn test_lethality_probit() {
        let pr = lethality_probit(200_000.0);
        // Should be a finite number
        assert!(pr.is_finite());
    }

    #[test]
    fn test_probit_to_probability_bounds() {
        let p1 = probit_to_probability(1.0);
        let p9 = probit_to_probability(9.0);
        assert!(p1 < 0.5);
        assert!(p9 > 0.5);
    }

    #[test]
    fn test_arrival_time_positive() {
        let ta = arrival_time(10.0, 1.0);
        assert!(ta > 0.0);
    }

    #[test]
    fn test_positive_phase_duration() {
        let td = positive_phase_duration(10.0, 1.0);
        assert!(td > 0.0);
    }

    #[test]
    fn test_negative_phase() {
        let (pn, tn) = negative_phase(10.0, 1.0);
        assert!(pn > 0.0);
        assert!(tn > 0.0);
    }

    #[test]
    fn test_energy_partition_sums() {
        let (b, g, c) = energy_partition(10.0, 1.0);
        let sum = b + g + c;
        assert!((sum - 1.0).abs() < TOL);
    }

    #[test]
    fn test_energy_partition_high_burst() {
        let (b, _g, _c) = energy_partition(100.0, 1.0);
        assert!(b > 0.7, "High burst should have most energy in blast");
    }
}
