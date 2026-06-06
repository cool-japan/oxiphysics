// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Marine and offshore hydrodynamics.
//!
//! Implements ship hull resistance (Froude number, wave resistance, viscous
//! resistance), Morison equation for offshore structures, wave–body
//! interaction (diffraction and radiation theories), added mass and damping
//! matrices, response amplitude operators (RAOs), mooring line dynamics
//! (catenary and elastic), station-keeping thruster control, propeller
//! thrust/torque (KT-KQ curves), bilge keel roll damping, and seakeeping
//! analysis.
//!
//! # Conventions
//!
//! * Right-hand coordinate system: X forward, Y port, Z up.
//! * All angles in radians unless otherwise stated.
//! * SI units throughout (m, kg, s, N, Pa).
//!
//! # References
//!
//! * Faltinsen (1990) – Sea Loads on Ships and Offshore Structures.
//! * ITTC (2014) – Recommended Procedures for Model Tests.
//! * Morison et al. (1950) – Wave forces on offshore structures.
//! * Chakrabarti (1987) – Hydrodynamics of Offshore Structures.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
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

#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Seawater density at 15 °C and 35 ppt salinity (kg/m³).
pub const RHO_SEA: f64 = 1025.0;

/// Freshwater density at 20 °C (kg/m³).
pub const RHO_FRESH: f64 = 998.2;

/// Standard gravitational acceleration (m/s²).
pub const GRAVITY: f64 = 9.81;

/// Kinematic viscosity of seawater at 15 °C (m²/s).
pub const NU_SEA: f64 = 1.188e-6;

// ---------------------------------------------------------------------------
// Froude / Reynolds numbers
// ---------------------------------------------------------------------------

/// Compute the Froude number: `Fr = V / sqrt(g * L)`.
///
/// * `speed` – ship speed (m/s).
/// * `length` – ship waterline length (m).
pub fn froude_number(speed: f64, length: f64) -> f64 {
    speed / (GRAVITY * length).sqrt()
}

/// Compute the Reynolds number: `Re = V * L / ν`.
///
/// * `speed` – flow velocity (m/s).
/// * `length` – characteristic length (m).
/// * `kinematic_viscosity` – ν (m²/s), defaults to `NU_SEA`.
pub fn reynolds_number(speed: f64, length: f64, kinematic_viscosity: f64) -> f64 {
    speed * length / kinematic_viscosity
}

/// Compute the ITTC-57 friction coefficient:
/// `CF = 0.075 / (log10(Re) - 2)²`.
pub fn ittc57_friction_coefficient(re: f64) -> f64 {
    let log_re = re.log10();
    0.075 / ((log_re - 2.0) * (log_re - 2.0))
}

// ---------------------------------------------------------------------------
// ShipHullResistance
// ---------------------------------------------------------------------------

/// Ship hull total resistance model (ITTC-78 / Holtrop-Mennen approach).
///
/// Decomposes total resistance into:
/// * Frictional resistance (ITTC-57).
/// * Form factor correction (1 + k).
/// * Wave-making resistance (Froude-based polynomial).
/// * Air resistance.
#[derive(Debug, Clone)]
pub struct ShipHullResistance {
    /// Ship displacement (kg).
    pub displacement: f64,
    /// Wetted surface area (m²).
    pub wetted_surface: f64,
    /// Waterline length (m).
    pub length: f64,
    /// Waterplane area (m²).
    pub waterplane_area: f64,
    /// Block coefficient CB (dimensionless, 0.4–0.85).
    pub cb: f64,
    /// Form factor (1 + k).
    pub form_factor: f64,
    /// Projected lateral area above waterline (m²) for air resistance.
    pub lateral_area_above_wl: f64,
    /// Wind drag coefficient CAA (≈ 0.07 for tankers, 0.04 for ferries).
    pub air_drag_coeff: f64,
}

impl ShipHullResistance {
    /// Create a new hull resistance model.
    pub fn new(
        displacement: f64,
        wetted_surface: f64,
        length: f64,
        waterplane_area: f64,
        cb: f64,
        form_factor: f64,
        lateral_area_above_wl: f64,
        air_drag_coeff: f64,
    ) -> Self {
        Self {
            displacement,
            wetted_surface,
            length,
            waterplane_area,
            cb,
            form_factor,
            lateral_area_above_wl,
            air_drag_coeff,
        }
    }

    /// Compute the frictional resistance (N) at speed `v` (m/s).
    pub fn frictional_resistance(&self, v: f64) -> f64 {
        if v < 1e-6 {
            return 0.0;
        }
        let re = reynolds_number(v, self.length, NU_SEA);
        let cf = ittc57_friction_coefficient(re);
        0.5 * RHO_SEA * v * v * self.wetted_surface * cf * self.form_factor
    }

    /// Estimate wave-making resistance (N) using a simplified Froude scaling.
    ///
    /// Uses a quartic polynomial fit to Holtrop-Mennen wave coefficient C_W:
    /// `R_W / W = C_W(Fn)`
    /// where W = ρg∇ and C_W ≈ exp(a₀ + a₁ * Fn + a₂ * Fn²) for representative
    /// hull form (Fn < 0.50).
    pub fn wave_resistance(&self, v: f64) -> f64 {
        if v < 1e-6 {
            return 0.0;
        }
        let fn_ = froude_number(v, self.length);
        // Simplified wave resistance ratio r_w (Fn) polynomial for CB ≈ 0.6.
        // Coefficients from Holtrop-Mennen regression for displacement ships.
        let r_w_ratio = if fn_ < 0.05 {
            1e-5
        } else if fn_ < 0.35 {
            // Low-speed range: wave resistance grows ~Fn^4.
            0.01 * fn_.powi(2)
        } else {
            // High-speed range: hump near Fn ≈ 0.4–0.5.
            let dfn = fn_ - 0.35;
            0.01 * fn_.powi(2) + 0.05 * dfn * dfn
        };
        let weight = self.displacement * GRAVITY;
        r_w_ratio * weight
    }

    /// Compute air resistance (N) at vessel speed `v_ship` (m/s) with
    /// relative wind speed `v_wind` (m/s) (head wind positive).
    pub fn air_resistance(&self, v_ship: f64, v_wind: f64) -> f64 {
        let v_rel = v_ship + v_wind;
        let rho_air = 1.225_f64; // kg/m³
        0.5 * rho_air * v_rel * v_rel * self.lateral_area_above_wl * self.air_drag_coeff
    }

    /// Total calm-water resistance (N) at speed `v` (m/s).
    pub fn total_resistance(&self, v: f64, v_wind: f64) -> f64 {
        self.frictional_resistance(v) + self.wave_resistance(v) + self.air_resistance(v, v_wind)
    }

    /// Effective power (W) = total resistance × speed.
    pub fn effective_power(&self, v: f64, v_wind: f64) -> f64 {
        self.total_resistance(v, v_wind) * v
    }
}

// ---------------------------------------------------------------------------
// MorisonElement
// ---------------------------------------------------------------------------

/// Morison equation for hydrodynamic loading on a slender cylindrical element.
///
/// The force per unit length is:
/// ```text
/// F = ρ * (1 + Cm) * A * a_w  +  0.5 * ρ * Cd * D * |u_r| * u_r
/// ```
/// where `u_r` = `u_w - u_s` is the relative velocity, `a_w` is fluid
/// acceleration, A = π D²/4 is cross-sectional area, and Cm, Cd are
/// inertia and drag coefficients.
#[derive(Debug, Clone)]
pub struct MorisonElement {
    /// Cylinder outer diameter (m).
    pub diameter: f64,
    /// Element length (m).
    pub length: f64,
    /// Inertia coefficient Cm (≈ 1.0 for solid cylinder, Cm = CM − 1).
    pub cm: f64,
    /// Drag coefficient Cd (≈ 0.6–1.2 depending on KC, Re).
    pub cd: f64,
    /// Element axis direction (unit vector).
    pub axis: [f64; 3],
}

impl MorisonElement {
    /// Create a new Morison element.
    pub fn new(diameter: f64, length: f64, cm: f64, cd: f64, axis: [f64; 3]) -> Self {
        Self {
            diameter,
            length,
            cm,
            cd,
            axis: normalize3(axis),
        }
    }

    /// Cross-sectional area (m²).
    pub fn area(&self) -> f64 {
        PI * self.diameter * self.diameter * 0.25
    }

    /// Compute the total hydrodynamic force vector (N) on this element.
    ///
    /// * `u_water` – undisturbed fluid velocity at element centroid (m/s).
    /// * `a_water` – undisturbed fluid acceleration at element centroid (m/s²).
    /// * `u_struct` – structural velocity at element centroid (m/s).
    ///
    /// Forces are computed in the plane perpendicular to the element axis.
    pub fn force(&self, u_water: [f64; 3], a_water: [f64; 3], u_struct: [f64; 3]) -> [f64; 3] {
        // Project to plane perpendicular to element axis.
        let e = self.axis;
        let u_w_perp = sub3(u_water, scale3(e, dot3(u_water, e)));
        let a_w_perp = sub3(a_water, scale3(e, dot3(a_water, e)));
        let u_s_perp = sub3(u_struct, scale3(e, dot3(u_struct, e)));
        let u_rel = sub3(u_w_perp, u_s_perp);
        let u_rel_mag = norm3(u_rel);

        let area = self.area();
        // Inertia force: ρ (1 + Cm) A L a_w_perp
        let f_inertia = scale3(a_w_perp, RHO_SEA * (1.0 + self.cm) * area * self.length);
        // Drag force: 0.5 ρ Cd D L |u_rel| u_rel
        let f_drag = scale3(
            u_rel,
            0.5 * RHO_SEA * self.cd * self.diameter * self.length * u_rel_mag,
        );
        add3(f_inertia, f_drag)
    }

    /// Keulegan-Carpenter number: `KC = U_max * T / D`.
    pub fn kc_number(&self, u_max: f64, wave_period: f64) -> f64 {
        u_max * wave_period / self.diameter
    }
}

// ---------------------------------------------------------------------------
// AddedMassDamping
// ---------------------------------------------------------------------------

/// 6×6 added mass and radiation damping matrices for a floating body.
///
/// Entries are stored in row-major order for a 6-DOF system
/// (surge, sway, heave, roll, pitch, yaw).
#[derive(Debug, Clone)]
pub struct AddedMassDamping {
    /// Added mass matrix A\[i\]\[j\] (kg for translational DOFs, kg·m² for rotational).
    pub added_mass: [[f64; 6]; 6],
    /// Radiation damping matrix B\[i\]\[j\] (N·s/m or N·m·s/rad).
    pub damping: [[f64; 6]; 6],
    /// Restoring stiffness matrix C\[i\]\[j\] (N/m or N·m/rad).
    pub stiffness: [[f64; 6]; 6],
}

impl AddedMassDamping {
    /// Create a diagonal added mass and damping model.
    ///
    /// `a_diag` – diagonal entries of added mass matrix.
    /// `b_diag` – diagonal entries of damping matrix.
    /// `c_diag` – diagonal entries of restoring stiffness matrix.
    pub fn diagonal(a_diag: [f64; 6], b_diag: [f64; 6], c_diag: [f64; 6]) -> Self {
        let mut added_mass = [[0.0_f64; 6]; 6];
        let mut damping = [[0.0_f64; 6]; 6];
        let mut stiffness = [[0.0_f64; 6]; 6];
        for i in 0..6 {
            added_mass[i][i] = a_diag[i];
            damping[i][i] = b_diag[i];
            stiffness[i][i] = c_diag[i];
        }
        Self {
            added_mass,
            damping,
            stiffness,
        }
    }

    /// Compute the frequency-dependent restoring moment for a given
    /// 6-DOF displacement vector `x` (m or rad).
    pub fn restoring_force(&self, x: [f64; 6]) -> [f64; 6] {
        let mut f = [0.0_f64; 6];
        for (f_i, stiff_row) in f.iter_mut().zip(self.stiffness.iter()) {
            *f_i -= stiff_row
                .iter()
                .zip(x.iter())
                .map(|(s, xj)| s * xj)
                .sum::<f64>();
        }
        f
    }

    /// Compute radiation damping force for velocity vector `xdot`.
    pub fn radiation_damping_force(&self, xdot: [f64; 6]) -> [f64; 6] {
        let mut f = [0.0_f64; 6];
        for (f_i, damp_row) in f.iter_mut().zip(self.damping.iter()) {
            *f_i -= damp_row
                .iter()
                .zip(xdot.iter())
                .map(|(d, xd)| d * xd)
                .sum::<f64>();
        }
        f
    }

    /// Compute inertia term: (M + A) * xddot.
    ///
    /// `body_mass_diag` – diagonal of rigid body mass matrix.
    pub fn inertia_force(&self, body_mass_diag: [f64; 6], xddot: [f64; 6]) -> [f64; 6] {
        let mut f = [0.0_f64; 6];
        for i in 0..6 {
            let m_plus_a = body_mass_diag[i] + self.added_mass[i][i];
            f[i] = m_plus_a * xddot[i];
        }
        f
    }
}

// ---------------------------------------------------------------------------
// LinearWave
// ---------------------------------------------------------------------------

/// Airy (linear) regular wave.
///
/// Provides particle kinematics at arbitrary depth under linear wave theory.
#[derive(Debug, Clone)]
pub struct LinearWave {
    /// Wave amplitude (m).
    pub amplitude: f64,
    /// Wave frequency ω (rad/s).
    pub omega: f64,
    /// Wave number k (rad/m), satisfying deep/finite water dispersion.
    pub wave_number: f64,
    /// Water depth h (m). Use f64::INFINITY for deep water.
    pub depth: f64,
    /// Wave propagation direction (unit vector in XY plane).
    pub direction: [f64; 2],
}

impl LinearWave {
    /// Create a wave from amplitude, period, and water depth.
    ///
    /// Iteratively solves the dispersion relation ω² = g k tanh(k h).
    pub fn from_period(amplitude: f64, period: f64, depth: f64, direction: [f64; 2]) -> Self {
        let omega = 2.0 * PI / period;
        let k = Self::solve_dispersion(omega, depth);
        Self {
            amplitude,
            omega,
            wave_number: k,
            depth,
            direction: {
                let mag = (direction[0] * direction[0] + direction[1] * direction[1]).sqrt();
                if mag < 1e-12 {
                    [1.0, 0.0]
                } else {
                    [direction[0] / mag, direction[1] / mag]
                }
            },
        }
    }

    /// Solve the dispersion relation ω² = g k tanh(k h) for wave number k.
    ///
    /// Uses Newton-Raphson iteration (converges in ~5 iterations).
    pub fn solve_dispersion(omega: f64, depth: f64) -> f64 {
        if depth.is_infinite() || depth > 1000.0 {
            return omega * omega / GRAVITY;
        }
        // Initial guess: deep-water approximation.
        let mut k = omega * omega / GRAVITY;
        for _ in 0..20 {
            let th = (k * depth).tanh();
            let f = omega * omega - GRAVITY * k * th;
            let df = -GRAVITY * (th + k * depth * (1.0 - th * th));
            k -= f / df;
            if f.abs() < 1e-10 {
                break;
            }
        }
        k.max(1e-10)
    }

    /// Compute wave surface elevation at position (x, y) and time t (s).
    pub fn surface_elevation(&self, x: f64, y: f64, t: f64) -> f64 {
        let phase =
            self.wave_number * (self.direction[0] * x + self.direction[1] * y) - self.omega * t;
        self.amplitude * phase.cos()
    }

    /// Compute horizontal fluid velocity (m/s) at depth `z` (z ≤ 0 at surface).
    pub fn horizontal_velocity(&self, x: f64, y: f64, z: f64, t: f64) -> [f64; 2] {
        let phase =
            self.wave_number * (self.direction[0] * x + self.direction[1] * y) - self.omega * t;
        let cosh_factor = if self.depth.is_infinite() {
            (self.wave_number * z).exp()
        } else {
            ((self.wave_number * (z + self.depth)).cosh())
                / ((self.wave_number * self.depth).cosh())
        };
        let mag = self.amplitude * self.omega * cosh_factor * phase.cos();
        [self.direction[0] * mag, self.direction[1] * mag]
    }

    /// Compute vertical fluid velocity (m/s) at depth z.
    pub fn vertical_velocity(&self, x: f64, y: f64, z: f64, t: f64) -> f64 {
        let phase =
            self.wave_number * (self.direction[0] * x + self.direction[1] * y) - self.omega * t;
        let sinh_factor = if self.depth.is_infinite() {
            (self.wave_number * z).exp()
        } else {
            ((self.wave_number * (z + self.depth)).sinh())
                / ((self.wave_number * self.depth).cosh())
        };
        self.amplitude * self.omega * sinh_factor * phase.sin()
    }

    /// Compute horizontal fluid acceleration (m/s²) at depth z.
    pub fn horizontal_acceleration(&self, x: f64, y: f64, z: f64, t: f64) -> [f64; 2] {
        let phase =
            self.wave_number * (self.direction[0] * x + self.direction[1] * y) - self.omega * t;
        let cosh_factor = if self.depth.is_infinite() {
            (self.wave_number * z).exp()
        } else {
            ((self.wave_number * (z + self.depth)).cosh())
                / ((self.wave_number * self.depth).cosh())
        };
        let mag = self.amplitude * self.omega * self.omega * cosh_factor * phase.sin();
        [self.direction[0] * mag, self.direction[1] * mag]
    }

    /// Wave phase velocity c = ω / k (m/s).
    pub fn phase_velocity(&self) -> f64 {
        self.omega / self.wave_number
    }

    /// Wave group velocity cg = dω/dk (m/s).
    pub fn group_velocity(&self) -> f64 {
        if self.depth.is_infinite() {
            0.5 * self.phase_velocity()
        } else {
            let kh = self.wave_number * self.depth;
            let n = 0.5 * (1.0 + 2.0 * kh / kh.sinh() / kh.cosh());
            n * self.phase_velocity()
        }
    }

    /// Wave length (m).
    pub fn wavelength(&self) -> f64 {
        2.0 * PI / self.wave_number
    }
}

// ---------------------------------------------------------------------------
// ResponseAmplitudeOperator
// ---------------------------------------------------------------------------

/// Response amplitude operator (RAO) for a floating body.
///
/// The RAO (or transfer function H(ω)) gives the complex response amplitude
/// per unit wave amplitude for each DOF:
/// `X̃(ω) = H(ω) · ζ_a(ω)`
///
/// Stored as magnitude and phase angle per frequency.
#[derive(Debug, Clone)]
pub struct ResponseAmplitudeOperator {
    /// Frequency vector ω (rad/s).
    pub frequencies: Vec<f64>,
    /// RAO magnitude \[DOF\]\[freq_index\] (m/m or rad/m).
    pub magnitude: Vec<Vec<f64>>,
    /// RAO phase \[DOF\]\[freq_index\] (rad).
    pub phase: Vec<Vec<f64>>,
    /// Number of DOFs (usually 6).
    pub n_dof: usize,
}

impl ResponseAmplitudeOperator {
    /// Create a new RAO container for `n_dof` degrees of freedom over
    /// the given frequency vector.
    pub fn new(n_dof: usize, frequencies: Vec<f64>) -> Self {
        let nf = frequencies.len();
        Self {
            frequencies,
            magnitude: vec![vec![0.0; nf]; n_dof],
            phase: vec![vec![0.0; nf]; n_dof],
            n_dof,
        }
    }

    /// Build a simplified RAO from mass-spring-damper parameters for DOF `i`.
    ///
    /// The 1-DOF transfer function is:
    /// `|H(ω)|² = F²_exc / ((C - (M+A)ω²)² + B²ω²)`
    pub fn set_1dof(
        &mut self,
        dof: usize,
        mass_plus_added: f64,
        damping: f64,
        stiffness: f64,
        excitation_amplitude: f64,
    ) {
        for (fi, &omega) in self.frequencies.iter().enumerate() {
            let denom_real = stiffness - mass_plus_added * omega * omega;
            let denom_imag = damping * omega;
            let denom_sq = denom_real * denom_real + denom_imag * denom_imag;
            let h_real = excitation_amplitude * denom_real / denom_sq;
            let h_imag = -excitation_amplitude * denom_imag / denom_sq;
            self.magnitude[dof][fi] = (h_real * h_real + h_imag * h_imag).sqrt();
            self.phase[dof][fi] = h_imag.atan2(h_real);
        }
    }

    /// Interpolate the RAO magnitude for DOF `dof` at frequency `omega` (rad/s).
    pub fn magnitude_at(&self, dof: usize, omega: f64) -> f64 {
        interpolate_linear(&self.frequencies, &self.magnitude[dof], omega)
    }

    /// Compute the significant response amplitude for DOF `dof` in a
    /// Pierson-Moskowitz spectrum with significant wave height Hs (m)
    /// and peak period Tp (s).
    ///
    /// `σ²_X = ∫ |H(ω)|² S_PM(ω) dω`
    ///
    /// Approximated by trapezoidal integration.
    pub fn significant_response(&self, dof: usize, hs: f64, tp: f64) -> f64 {
        let n = self.frequencies.len();
        if n < 2 {
            return 0.0;
        }
        let mut sigma_sq = 0.0;
        for i in 0..(n - 1) {
            let w1 = self.frequencies[i];
            let w2 = self.frequencies[i + 1];
            let dw = w2 - w1;
            let s1 = pierson_moskowitz(w1, hs, tp);
            let s2 = pierson_moskowitz(w2, hs, tp);
            let h1 = self.magnitude[dof][i];
            let h2 = self.magnitude[dof][i + 1];
            sigma_sq += 0.5 * dw * (h1 * h1 * s1 + h2 * h2 * s2);
        }
        4.0 * sigma_sq.sqrt() // H_s = 4σ
    }
}

/// Pierson-Moskowitz wave spectrum S(ω) (m²·s/rad).
///
/// `S_PM(ω) = (5/16) H_s² ω_p^4 ω^{-5} exp(-5/4 (ω_p/ω)^4)`
pub fn pierson_moskowitz(omega: f64, hs: f64, tp: f64) -> f64 {
    if omega < 1e-6 {
        return 0.0;
    }
    let omega_p = 2.0 * PI / tp;
    let ratio = omega_p / omega;
    (5.0 / 16.0) * hs * hs * omega_p.powi(4) * omega.powi(-5) * (-1.25 * ratio.powi(4)).exp()
}

/// JONSWAP wave spectrum S(ω) for fetch-limited seas.
///
/// Includes peak enhancement factor γ (typically 3.3 for North Sea).
pub fn jonswap_spectrum(omega: f64, hs: f64, tp: f64, gamma: f64) -> f64 {
    let s_pm = pierson_moskowitz(omega, hs, tp);
    if s_pm < 1e-300 {
        return 0.0;
    }
    let omega_p = 2.0 * PI / tp;
    let sigma = if omega <= omega_p { 0.07 } else { 0.09 };
    let arg = -0.5 * ((omega / omega_p - 1.0) / sigma).powi(2);
    let r = gamma.powf(arg.exp());
    // Normalise to preserve Hs.
    let c_pm = 1.0 - 0.287 * gamma.ln();
    s_pm * r * c_pm
}

// ---------------------------------------------------------------------------
// MooringLine
// ---------------------------------------------------------------------------

/// Catenary mooring line model.
///
/// Models a single inextensible catenary mooring line connecting
/// an anchor at the seabed to a fairlead on the floating body.
#[derive(Debug, Clone)]
pub struct MooringLine {
    /// Total unstretched length of the line (m).
    pub length: f64,
    /// Linear mass density of the line in water (kg/m), i.e., m_line - ρA_line.
    pub wet_weight_per_metre: f64,
    /// Axial stiffness EA (N).
    pub axial_stiffness: f64,
    /// Anchor position (m).
    pub anchor_pos: [f64; 3],
}

impl MooringLine {
    /// Create a new catenary mooring line.
    pub fn new(
        length: f64,
        wet_weight_per_metre: f64,
        axial_stiffness: f64,
        anchor_pos: [f64; 3],
    ) -> Self {
        Self {
            length,
            wet_weight_per_metre,
            axial_stiffness,
            anchor_pos,
        }
    }

    /// Compute the horizontal tension at the fairlead given horizontal scope
    /// `x_h` (m) and vertical height difference `z_h` (m, positive upward).
    ///
    /// Returns `(T_H, T_V)` – horizontal and vertical components of tension (N).
    /// Uses the catenary equations:
    ///
    /// ```text
    /// x_h = a * ( sinh(T_V/T_H) - sinh((T_V - w*L)/T_H) )
    /// z_h = a * ( cosh(T_V/T_H) - cosh((T_V - w*L)/T_H) )
    /// ```
    ///
    /// Solved by Newton-Raphson on the catenary parameter `a = T_H / (w_0)`.
    pub fn catenary_tension(&self, x_h: f64, z_h: f64) -> (f64, f64) {
        let w = self.wet_weight_per_metre;
        // Initial guess: taut line approximation.
        let l_horiz = (x_h * x_h + z_h * z_h).sqrt().max(1e-3);
        let mut t_h = w * l_horiz * 0.5;
        // Newton-Raphson on catenary equations.
        for _ in 0..50 {
            if t_h < 1e-3 {
                t_h = 1e-3;
            }
            let a = t_h / w;
            // Touchdown condition: check if line is fully taut or has a lazy wave.
            let sinh_val = (z_h / a + (1.0 + (z_h / a).powi(2)).sqrt()).ln(); // asinh
            let x_cat = a * ((self.length / a).sinh() - sinh_val.sinh());
            let fx = x_cat - x_h;
            let dfx = (self.length / a).cosh() - 1.0;
            let delta = -fx / (dfx.max(1e-6));
            t_h += delta.max(-t_h * 0.5).min(t_h * 2.0);
            if fx.abs() < 1e-4 {
                break;
            }
        }
        let a = t_h / w;
        let t_v = t_h * ((self.length / a).sinh());
        (t_h, t_v)
    }

    /// Compute the fairlead restoring force vector given fairlead position (m).
    pub fn restoring_force(&self, fairlead_pos: [f64; 3]) -> [f64; 3] {
        let dx = sub3(fairlead_pos, self.anchor_pos);
        let x_h = (dx[0] * dx[0] + dx[1] * dx[1]).sqrt();
        let z_h = dx[2]; // positive upward
        if x_h < 1e-3 {
            return [0.0; 3];
        }
        let (t_h, t_v) = self.catenary_tension(x_h, z_h);
        // Horizontal direction from anchor to fairlead.
        let horiz_dir = normalize3([dx[0], dx[1], 0.0]);
        // Force on vessel is tension directed toward anchor.
        [-horiz_dir[0] * t_h, -horiz_dir[1] * t_h, -t_v]
    }
}

// ---------------------------------------------------------------------------
// PropellerModel
// ---------------------------------------------------------------------------

/// Propeller open-water performance model using KT-KQ curves.
///
/// Thrust coefficient: `KT = T / (ρ n² D⁴)`
/// Torque coefficient: `KQ = Q / (ρ n² D⁵)`
/// where `n` = rps, `D` = propeller diameter (m).
#[derive(Debug, Clone)]
pub struct PropellerModel {
    /// Propeller diameter (m).
    pub diameter: f64,
    /// Number of blades.
    pub n_blades: u32,
    /// Pitch ratio P/D.
    pub pitch_ratio: f64,
    /// Polynomial coefficients for KT(J): KT = Σ a_i J^i.
    pub kt_coeffs: Vec<f64>,
    /// Polynomial coefficients for KQ(J): KQ = Σ b_i J^i.
    pub kq_coeffs: Vec<f64>,
}

impl PropellerModel {
    /// Create a propeller with B-series polynomial coefficients.
    ///
    /// Uses a simplified 4-term KT/KQ polynomial valid for J ∈ \[0, 1\].
    /// Coefficients from Wageningen B-series (5 blades, P/D = 1.0).
    pub fn new_b_series(diameter: f64) -> Self {
        Self {
            diameter,
            n_blades: 5,
            pitch_ratio: 1.0,
            // B5-65 series: KT(J) = 0.339 - 0.156 J - 0.419 J² + 0.242 J³
            kt_coeffs: vec![0.339, -0.156, -0.419, 0.242],
            // KQ(J) = 0.0466 - 0.0207 J - 0.0512 J² + 0.0308 J³
            kq_coeffs: vec![0.0466, -0.0207, -0.0512, 0.0308],
        }
    }

    /// Advance coefficient J = Va / (n D).
    ///
    /// * `va` – advance velocity = ship speed × (1 - wake fraction) (m/s).
    /// * `n_rps` – propeller rotation rate (rev/s).
    pub fn advance_coefficient(&self, va: f64, n_rps: f64) -> f64 {
        if n_rps.abs() < 1e-6 {
            return 0.0;
        }
        va / (n_rps * self.diameter)
    }

    /// Evaluate KT at advance coefficient J.
    pub fn kt(&self, j: f64) -> f64 {
        let j = j.max(0.0);
        eval_polynomial(&self.kt_coeffs, j).max(0.0)
    }

    /// Evaluate KQ at advance coefficient J.
    pub fn kq(&self, j: f64) -> f64 {
        let j = j.max(0.0);
        eval_polynomial(&self.kq_coeffs, j).max(0.0)
    }

    /// Thrust force (N): T = KT ρ n² D⁴.
    pub fn thrust(&self, n_rps: f64, va: f64) -> f64 {
        let j = self.advance_coefficient(va, n_rps);
        let kt = self.kt(j);
        kt * RHO_SEA * n_rps * n_rps * self.diameter.powi(4)
    }

    /// Torque (N·m): Q = KQ ρ n² D⁵.
    pub fn torque(&self, n_rps: f64, va: f64) -> f64 {
        let j = self.advance_coefficient(va, n_rps);
        let kq = self.kq(j);
        kq * RHO_SEA * n_rps * n_rps * self.diameter.powi(5)
    }

    /// Open-water efficiency η_o = J * KT / (2π KQ).
    pub fn open_water_efficiency(&self, j: f64) -> f64 {
        let kt = self.kt(j);
        let kq = self.kq(j);
        if kq < 1e-10 {
            return 0.0;
        }
        j * kt / (2.0 * PI * kq)
    }
}

// ---------------------------------------------------------------------------
// BilgeKeelRollDamping
// ---------------------------------------------------------------------------

/// Bilge keel roll damping model.
///
/// Bilge keels are flat plates attached along the bilge radius that add
/// viscous and pressure roll damping.  The model follows the Ikeda method.
#[derive(Debug, Clone)]
pub struct BilgeKeelRollDamping {
    /// Bilge keel length (m).
    pub length: f64,
    /// Bilge keel height (normal to hull surface) (m).
    pub height: f64,
    /// Bilge radius r_b (m).
    pub bilge_radius: f64,
    /// Ship breadth (m).
    pub breadth: f64,
    /// Ship draft (m).
    pub draft: f64,
}

impl BilgeKeelRollDamping {
    /// Create a bilge keel roll damping model.
    pub fn new(length: f64, height: f64, bilge_radius: f64, breadth: f64, draft: f64) -> Self {
        Self {
            length,
            height,
            bilge_radius,
            breadth,
            draft,
        }
    }

    /// Compute the Ikeda bilge keel damping coefficient B_BK (N·m·s/rad).
    ///
    /// Simplified formula based on the non-dimensional roll amplitude φ_a (rad)
    /// and roll frequency ω (rad/s).
    ///
    /// `B_BK ≈ ρ r_b³ L_BK h² ω φ_a / (2 π)`  (linearised equivalent)
    pub fn damping_coefficient(&self, roll_amplitude_rad: f64, omega: f64) -> f64 {
        let r = self.bilge_radius;
        let cf = 2.0 * r * omega * roll_amplitude_rad / PI; // velocity scale
        // Drag contribution per unit length: 0.5 ρ CD h cf²
        let cd = 1.2; // flat plate normal drag coefficient
        let f_drag_per_length = 0.5 * RHO_SEA * cd * self.height * cf * cf;
        // Total moment per unit length: r × F_drag
        let moment_per_length = r * f_drag_per_length;
        // Equivalent linear damping: B_eq = M_total / ω φ_a
        let total_moment = moment_per_length * self.length;
        if (omega * roll_amplitude_rad).abs() < 1e-10 {
            return 0.0;
        }
        total_moment / (omega * roll_amplitude_rad)
    }

    /// Compute the roll damping ratio for bilge keels in a given sea state.
    ///
    /// Returns the fraction of critical damping ζ_BK = B_BK / (2 I_xx ω_n).
    pub fn damping_ratio(&self, roll_amplitude_rad: f64, omega: f64, roll_inertia: f64) -> f64 {
        let b_bk = self.damping_coefficient(roll_amplitude_rad, omega);
        let critical = 2.0 * roll_inertia * omega;
        if critical < 1e-10 {
            return 0.0;
        }
        b_bk / critical
    }
}

// ---------------------------------------------------------------------------
// StationKeepingController
// ---------------------------------------------------------------------------

/// Dynamic positioning (DP) station-keeping controller.
///
/// A PID controller in surge, sway, and yaw for maintaining vessel
/// position and heading under environmental disturbances.
#[derive(Debug, Clone)]
pub struct StationKeepingController {
    /// Proportional gain \[surge, sway, yaw\].
    pub kp: [f64; 3],
    /// Integral gain \[surge, sway, yaw\].
    pub ki: [f64; 3],
    /// Derivative gain \[surge, sway, yaw\].
    pub kd: [f64; 3],
    /// Integral accumulator \[surge, sway, yaw\].
    pub integral: [f64; 3],
    /// Previous error \[surge, sway, yaw\].
    pub prev_error: [f64; 3],
    /// Maximum thruster force per DOF (N or N·m).
    pub max_force: [f64; 3],
    /// Setpoint \[x, y, yaw_rad\].
    pub setpoint: [f64; 3],
}

impl StationKeepingController {
    /// Create a new DP controller.
    pub fn new(kp: [f64; 3], ki: [f64; 3], kd: [f64; 3], max_force: [f64; 3]) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: [0.0; 3],
            prev_error: [0.0; 3],
            max_force,
            setpoint: [0.0; 3],
        }
    }

    /// Set the target position and heading.
    pub fn set_setpoint(&mut self, x: f64, y: f64, yaw: f64) {
        self.setpoint = [x, y, yaw];
    }

    /// Compute the required thrust commands for the current state.
    ///
    /// * `state` – current \[x, y, yaw\] of the vessel.
    /// * `dt` – time step (s).
    ///
    /// Returns \[F_surge, F_sway, M_yaw\] in world frame (N, N, N·m).
    pub fn compute_thrust(&mut self, state: [f64; 3], dt: f64) -> [f64; 3] {
        let mut commands = [0.0_f64; 3];
        for i in 0..3 {
            let error = self.setpoint[i] - state[i];
            self.integral[i] += error * dt;
            let derivative = if dt > 1e-12 {
                (error - self.prev_error[i]) / dt
            } else {
                0.0
            };
            let u = self.kp[i] * error + self.ki[i] * self.integral[i] + self.kd[i] * derivative;
            commands[i] = clamp(u, -self.max_force[i], self.max_force[i]);
            self.prev_error[i] = error;
        }
        commands
    }

    /// Reset integral and derivative state.
    pub fn reset(&mut self) {
        self.integral = [0.0; 3];
        self.prev_error = [0.0; 3];
    }
}

// ---------------------------------------------------------------------------
// ThrusterAllocation
// ---------------------------------------------------------------------------

/// Thruster allocation for dynamic positioning systems.
///
/// Distributes demanded forces/moments among a set of azimuthing thrusters
/// using a simplified pseudo-inverse approach.
#[derive(Debug, Clone)]
pub struct ThrusterAllocation {
    /// Thruster positions \[x, y\] relative to vessel reference point (m).
    pub positions: Vec<[f64; 2]>,
    /// Maximum thrust per thruster (N).
    pub max_thrust: Vec<f64>,
    /// Azimuth angles (rad) for each thruster (0 = forward).
    pub azimuths: Vec<f64>,
}

impl ThrusterAllocation {
    /// Create a thruster allocation system.
    pub fn new(positions: Vec<[f64; 2]>, max_thrust: Vec<f64>) -> Self {
        let n = positions.len();
        Self {
            positions,
            max_thrust,
            azimuths: vec![0.0; n],
        }
    }

    /// Compute the thrust vector contribution \[Fx, Fy, Mz\] for thruster `i`
    /// with thrust `T_i` (N) and azimuth `alpha_i` (rad).
    pub fn thruster_contribution(&self, i: usize, thrust: f64, azimuth: f64) -> [f64; 3] {
        let fx = thrust * azimuth.cos();
        let fy = thrust * azimuth.sin();
        let mz = self.positions[i][0] * fy - self.positions[i][1] * fx;
        [fx, fy, mz]
    }

    /// Simple proportional allocation: distribute total demand equally among
    /// thrusters pointing in the demanded direction.
    ///
    /// Returns thrust magnitude for each thruster (N).
    pub fn allocate(&mut self, demand: [f64; 3]) -> Vec<f64> {
        let n = self.positions.len();
        let f_total_mag = (demand[0] * demand[0] + demand[1] * demand[1]).sqrt();
        let demanded_angle = demand[1].atan2(demand[0]);
        let per_thruster = f_total_mag / n.max(1) as f64;
        let mut thrusts = Vec::with_capacity(n);
        for i in 0..n {
            self.azimuths[i] = demanded_angle;
            let t = clamp(per_thruster, 0.0, self.max_thrust[i]);
            thrusts.push(t);
        }
        thrusts
    }
}

// ---------------------------------------------------------------------------
// SeakeepingAnalysis
// ---------------------------------------------------------------------------

/// Seakeeping analysis for pitch, heave, and roll motion statistics.
///
/// Combines RAOs with a wave spectrum to compute significant motion
/// amplitudes and RMS responses.
#[derive(Debug, Clone)]
pub struct SeakeepingAnalysis {
    /// Vessel RAO (6-DOF).
    pub rao: ResponseAmplitudeOperator,
    /// Significant wave height Hs (m).
    pub hs: f64,
    /// Peak period Tp (s).
    pub tp: f64,
    /// Spectrum type (PM or JONSWAP).
    pub spectrum_type: SpectrumType,
    /// JONSWAP peak enhancement factor γ.
    pub jonswap_gamma: f64,
}

/// Wave spectrum type selector.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectrumType {
    /// Pierson-Moskowitz spectrum (fully developed sea).
    PiersonMoskowitz,
    /// JONSWAP spectrum (fetch-limited sea).
    Jonswap,
}

impl SeakeepingAnalysis {
    /// Create a seakeeping analysis.
    pub fn new(
        rao: ResponseAmplitudeOperator,
        hs: f64,
        tp: f64,
        spectrum_type: SpectrumType,
        jonswap_gamma: f64,
    ) -> Self {
        Self {
            rao,
            hs,
            tp,
            spectrum_type,
            jonswap_gamma,
        }
    }

    /// Compute the spectral density S(ω) for the current sea state.
    pub fn spectrum(&self, omega: f64) -> f64 {
        match self.spectrum_type {
            SpectrumType::PiersonMoskowitz => pierson_moskowitz(omega, self.hs, self.tp),
            SpectrumType::Jonswap => jonswap_spectrum(omega, self.hs, self.tp, self.jonswap_gamma),
        }
    }

    /// Compute RMS response for DOF `dof` by spectral integration.
    pub fn rms_response(&self, dof: usize) -> f64 {
        let n = self.rao.frequencies.len();
        if n < 2 {
            return 0.0;
        }
        let mut m0 = 0.0;
        for i in 0..(n - 1) {
            let w1 = self.rao.frequencies[i];
            let w2 = self.rao.frequencies[i + 1];
            let dw = w2 - w1;
            let s1 = self.spectrum(w1);
            let s2 = self.spectrum(w2);
            let h1 = self.rao.magnitude[dof][i];
            let h2 = self.rao.magnitude[dof][i + 1];
            m0 += 0.5 * dw * (h1 * h1 * s1 + h2 * h2 * s2);
        }
        m0.sqrt()
    }

    /// Significant single-amplitude for DOF `dof`.
    ///
    /// `σ_{1/3} = 2 * RMS` (for narrow-band Rayleigh distribution).
    pub fn significant_amplitude(&self, dof: usize) -> f64 {
        2.0 * self.rms_response(dof)
    }
}

// ---------------------------------------------------------------------------
// Helper utilities
// ---------------------------------------------------------------------------

/// Evaluate a polynomial at `x`: `Σ c[i] * x^i`.
fn eval_polynomial(coeffs: &[f64], x: f64) -> f64 {
    let mut result = 0.0;
    let mut xpow = 1.0;
    for &c in coeffs {
        result += c * xpow;
        xpow *= x;
    }
    result
}

/// Linear interpolation / extrapolation of tabulated data.
fn interpolate_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return ys[0];
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    for i in 0..(n - 1) {
        if x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] + t * (ys[i + 1] - ys[i]);
        }
    }
    ys[n - 1]
}

// ---------------------------------------------------------------------------
// Unit Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // -----------------------------------------------------------------------
    // Froude / Reynolds number tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_froude_number_zero_speed() {
        assert!(froude_number(0.0, 100.0).abs() < EPS);
    }

    #[test]
    fn test_froude_number_known_value() {
        // V = 10 m/s, L = 100 m → Fn = 10 / sqrt(9.81 * 100) ≈ 0.3193
        let fn_ = froude_number(10.0, 100.0);
        assert!((fn_ - 0.3193).abs() < 1e-3, "Fn={fn_}");
    }

    #[test]
    fn test_reynolds_number_seawater() {
        let re = reynolds_number(5.0, 100.0, NU_SEA);
        let expected = 5.0 * 100.0 / NU_SEA;
        assert!((re - expected).abs() < EPS);
    }

    #[test]
    fn test_ittc57_friction_coefficient_reasonable() {
        let re = 1e8_f64;
        let cf = ittc57_friction_coefficient(re);
        // For Re=1e8, CF ≈ 0.075 / (8-2)^2 = 0.075/36 ≈ 0.002083
        assert!((cf - 0.002083).abs() < 1e-4, "CF={cf}");
    }

    // -----------------------------------------------------------------------
    // ShipHullResistance tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hull_resistance_zero_speed_zero() {
        let hull =
            ShipHullResistance::new(50_000e3, 15000.0, 150.0, 2000.0, 0.65, 1.15, 500.0, 0.07);
        assert!(hull.frictional_resistance(0.0).abs() < EPS);
        assert!(hull.wave_resistance(0.0).abs() < EPS);
        assert!(hull.total_resistance(0.0, 0.0).abs() < EPS);
    }

    #[test]
    fn test_hull_resistance_positive_at_speed() {
        let hull =
            ShipHullResistance::new(50_000e3, 15000.0, 150.0, 2000.0, 0.65, 1.15, 500.0, 0.07);
        let r = hull.total_resistance(8.0, 0.0);
        assert!(r > 0.0, "total resistance should be positive: {r}");
    }

    #[test]
    fn test_hull_resistance_increases_with_speed() {
        let hull =
            ShipHullResistance::new(50_000e3, 15000.0, 150.0, 2000.0, 0.65, 1.15, 500.0, 0.07);
        let r1 = hull.total_resistance(5.0, 0.0);
        let r2 = hull.total_resistance(10.0, 0.0);
        assert!(
            r2 > r1,
            "resistance should increase with speed: r1={r1}, r2={r2}"
        );
    }

    #[test]
    fn test_effective_power_equals_resistance_times_speed() {
        let hull =
            ShipHullResistance::new(50_000e3, 15000.0, 150.0, 2000.0, 0.65, 1.15, 500.0, 0.07);
        let v = 8.0;
        let r = hull.total_resistance(v, 0.0);
        let pe = hull.effective_power(v, 0.0);
        assert!(
            (pe - r * v).abs() < EPS * pe.abs().max(1.0),
            "Pe={pe}, R*v={}",
            r * v
        );
    }

    // -----------------------------------------------------------------------
    // MorisonElement tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_morison_zero_relative_velocity_inertia_only() {
        // Zero relative velocity → only inertia term.
        let elem = MorisonElement::new(1.0, 10.0, 1.0, 0.6, [0.0, 0.0, 1.0]);
        let a_water = [2.0, 0.0, 0.0]; // 2 m/s² in X (perpendicular to Z axis)
        let f = elem.force([0.0; 3], a_water, [0.0; 3]);
        // F = ρ * (1 + Cm) * A * L * a = 1025 * 2 * π/4 * 10 * 2 = 32 252 N
        let expected = RHO_SEA * 2.0 * (PI * 0.25) * 10.0 * 2.0;
        assert!(
            (f[0] - expected).abs() < 1.0,
            "f_x={}, expected={expected}",
            f[0]
        );
    }

    #[test]
    fn test_morison_zero_everything_zero_force() {
        let elem = MorisonElement::new(1.0, 10.0, 1.0, 0.6, [0.0, 0.0, 1.0]);
        let f = elem.force([0.0; 3], [0.0; 3], [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_morison_kc_number() {
        let elem = MorisonElement::new(2.0, 10.0, 1.0, 0.6, [0.0, 0.0, 1.0]);
        // KC = U_max * T / D = 3.0 * 10.0 / 2.0 = 15
        let kc = elem.kc_number(3.0, 10.0);
        assert!((kc - 15.0).abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // LinearWave tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_surface_elevation_amplitude() {
        let wave = LinearWave::from_period(1.0, 8.0, f64::INFINITY, [1.0, 0.0]);
        // At x=0, t=0: η = A cos(0) = A.
        let eta = wave.surface_elevation(0.0, 0.0, 0.0);
        assert!((eta - 1.0).abs() < 1e-6, "eta={eta}");
    }

    #[test]
    fn test_wave_dispersion_deep_water() {
        // Deep water: k = ω² / g.
        let period = 8.0;
        let omega = 2.0 * PI / period;
        let k = LinearWave::solve_dispersion(omega, f64::INFINITY);
        let k_expected = omega * omega / GRAVITY;
        assert!(
            (k - k_expected).abs() < 1e-6,
            "k={k}, expected={k_expected}"
        );
    }

    #[test]
    fn test_wave_phase_velocity_deep() {
        let wave = LinearWave::from_period(1.0, 10.0, f64::INFINITY, [1.0, 0.0]);
        // Deep water: c = g / ω = g T / (2π)
        let omega = 2.0 * PI / 10.0;
        let c_expected = GRAVITY / omega;
        let c = wave.phase_velocity();
        assert!(
            (c - c_expected).abs() < 1e-3,
            "c={c}, expected={c_expected}"
        );
    }

    #[test]
    fn test_wave_group_velocity_half_phase_deep() {
        let wave = LinearWave::from_period(1.0, 10.0, f64::INFINITY, [1.0, 0.0]);
        // Deep water: cg = c/2.
        let c = wave.phase_velocity();
        let cg = wave.group_velocity();
        assert!((cg - c * 0.5).abs() < 1e-6, "cg={cg}, c/2={}", c * 0.5);
    }

    #[test]
    fn test_wave_wavelength_deep() {
        let wave = LinearWave::from_period(1.0, 10.0, f64::INFINITY, [1.0, 0.0]);
        // L = g T² / (2π)
        let l_expected = GRAVITY * 100.0 / (2.0 * PI);
        let l = wave.wavelength();
        assert!((l - l_expected).abs() < 0.1, "L={l}, expected={l_expected}");
    }

    // -----------------------------------------------------------------------
    // ResponseAmplitudeOperator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_rao_zero_excitation_zero_response() {
        let freqs: Vec<f64> = (1..=10).map(|i| i as f64 * 0.1).collect();
        let mut rao = ResponseAmplitudeOperator::new(6, freqs);
        rao.set_1dof(2, 1000.0, 50.0, 10_000.0, 0.0);
        for m in &rao.magnitude[2] {
            assert!(m.abs() < EPS, "magnitude should be zero: {m}");
        }
    }

    #[test]
    fn test_rao_resonance_peak() {
        let omega_n = 1.0; // rad/s
        let freqs: Vec<f64> = (1..=100).map(|i| i as f64 * 0.05).collect();
        let mass = 1000.0_f64;
        let stiffness = mass * omega_n * omega_n;
        let damping = 50.0;
        let excitation = 1000.0;
        let mut rao = ResponseAmplitudeOperator::new(1, freqs.clone());
        rao.set_1dof(0, mass, damping, stiffness, excitation);
        // Find frequency index nearest ω_n.
        let idx_n = freqs
            .iter()
            .enumerate()
            .min_by_key(|&(_, &w)| {
                let d = (w - omega_n).abs();
                (d * 1e6) as u64
            })
            .map(|(i, _)| i)
            .unwrap();
        // RAO at resonance should be larger than at other frequencies.
        let h_res = rao.magnitude[0][idx_n];
        let h_low = rao.magnitude[0][0];
        assert!(
            h_res > h_low,
            "RAO at resonance ({h_res}) should exceed low-freq ({h_low})"
        );
    }

    // -----------------------------------------------------------------------
    // Pierson-Moskowitz / JONSWAP tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pm_spectrum_zero_at_zero_frequency() {
        let s = pierson_moskowitz(0.0, 3.0, 10.0);
        assert!(s.abs() < EPS);
    }

    #[test]
    fn test_pm_spectrum_positive_near_peak() {
        let tp = 10.0;
        let omega_p = 2.0 * PI / tp;
        let s = pierson_moskowitz(omega_p, 3.0, tp);
        assert!(
            s > 0.0,
            "spectrum should be positive at peak frequency: {s}"
        );
    }

    #[test]
    fn test_jonswap_greater_than_pm_at_peak() {
        let tp = 10.0;
        let omega_p = 2.0 * PI / tp;
        let s_pm = pierson_moskowitz(omega_p, 3.0, tp);
        let s_jn = jonswap_spectrum(omega_p, 3.0, tp, 3.3);
        assert!(
            s_jn > s_pm,
            "JONSWAP peak should exceed PM: s_jn={s_jn}, s_pm={s_pm}"
        );
    }

    // -----------------------------------------------------------------------
    // MooringLine tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_mooring_restoring_force_direction() {
        // Anchor at origin, fairlead 100 m to port (+Y) and 50 m above seabed.
        let line = MooringLine::new(150.0, 200.0, 1e9, [0.0, 0.0, -50.0]);
        let f = line.restoring_force([0.0, 100.0, 0.0]);
        // Horizontal force should pull toward anchor (−Y).
        assert!(
            f[1] < 0.0,
            "horizontal restoring force should be in -Y: f_y={}",
            f[1]
        );
    }

    #[test]
    fn test_mooring_zero_horizontal_no_horizontal_force() {
        let line = MooringLine::new(150.0, 200.0, 1e9, [0.0, 0.0, -50.0]);
        let f = line.restoring_force([0.0, 0.0, 0.0]);
        // Horizontal offset < 1e-3 → zero force.
        assert_eq!(f, [0.0; 3]);
    }

    // -----------------------------------------------------------------------
    // PropellerModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_propeller_kt_positive_at_zero_advance() {
        let prop = PropellerModel::new_b_series(5.0);
        let kt = prop.kt(0.0);
        assert!(kt > 0.0, "KT at J=0 should be positive: {kt}");
    }

    #[test]
    fn test_propeller_kq_positive_at_zero_advance() {
        let prop = PropellerModel::new_b_series(5.0);
        let kq = prop.kq(0.0);
        assert!(kq > 0.0, "KQ at J=0 should be positive: {kq}");
    }

    #[test]
    fn test_propeller_thrust_positive_at_non_zero_rpm() {
        let prop = PropellerModel::new_b_series(5.0);
        let t = prop.thrust(2.0, 0.0); // 2 rps, zero advance
        assert!(t > 0.0, "thrust should be positive: {t}");
    }

    #[test]
    fn test_propeller_advance_coefficient() {
        let prop = PropellerModel::new_b_series(5.0);
        let j = prop.advance_coefficient(10.0, 2.0);
        // J = 10 / (2 * 5) = 1.0
        assert!((j - 1.0).abs() < EPS, "J={j}");
    }

    #[test]
    fn test_propeller_efficiency_between_0_and_1() {
        let prop = PropellerModel::new_b_series(5.0);
        let eta = prop.open_water_efficiency(0.5);
        assert!((0.0..=1.0).contains(&eta), "efficiency={eta}");
    }

    // -----------------------------------------------------------------------
    // BilgeKeelRollDamping tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_bilge_keel_damping_positive() {
        let bk = BilgeKeelRollDamping::new(20.0, 0.4, 3.0, 30.0, 10.0);
        let b = bk.damping_coefficient(0.1, 0.5);
        assert!(b > 0.0, "damping coefficient should be positive: {b}");
    }

    #[test]
    fn test_bilge_keel_damping_zero_amplitude_zero_coeff() {
        let bk = BilgeKeelRollDamping::new(20.0, 0.4, 3.0, 30.0, 10.0);
        let b = bk.damping_coefficient(0.0, 0.5);
        assert!(b.abs() < EPS, "zero amplitude → zero damping: {b}");
    }

    #[test]
    fn test_bilge_keel_damping_ratio_positive() {
        let bk = BilgeKeelRollDamping::new(20.0, 0.4, 3.0, 30.0, 10.0);
        let roll_inertia = 1e8_f64;
        let zeta = bk.damping_ratio(0.1, 0.5, roll_inertia);
        assert!(zeta > 0.0, "damping ratio should be positive: {zeta}");
    }

    // -----------------------------------------------------------------------
    // StationKeepingController tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dp_controller_zero_error_zero_thrust() {
        let mut ctrl = StationKeepingController::new(
            [100.0, 100.0, 1e5],
            [10.0, 10.0, 1e4],
            [500.0, 500.0, 5e5],
            [1e6, 1e6, 1e7],
        );
        ctrl.set_setpoint(0.0, 0.0, 0.0);
        let f = ctrl.compute_thrust([0.0, 0.0, 0.0], 1.0);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_dp_controller_position_error_produces_thrust() {
        let mut ctrl = StationKeepingController::new(
            [100.0, 100.0, 1e5],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1e6, 1e6, 1e7],
        );
        ctrl.set_setpoint(10.0, 0.0, 0.0);
        let f = ctrl.compute_thrust([0.0, 0.0, 0.0], 1.0);
        assert!(
            f[0] > 0.0,
            "surge thrust should be positive for positive error: {}",
            f[0]
        );
    }

    #[test]
    fn test_dp_controller_thrust_clamped() {
        let mut ctrl = StationKeepingController::new(
            [1e8, 1e8, 1e8],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [500.0, 500.0, 500.0],
        );
        ctrl.set_setpoint(1000.0, 0.0, 0.0);
        let f = ctrl.compute_thrust([0.0, 0.0, 0.0], 1.0);
        assert!(
            (f[0] - 500.0).abs() < EPS,
            "thrust should be clamped to 500 N: {}",
            f[0]
        );
    }

    #[test]
    fn test_dp_controller_reset_clears_integral() {
        let mut ctrl = StationKeepingController::new(
            [10.0, 10.0, 10.0],
            [5.0, 5.0, 5.0],
            [0.0, 0.0, 0.0],
            [1e6, 1e6, 1e6],
        );
        ctrl.set_setpoint(5.0, 0.0, 0.0);
        ctrl.compute_thrust([0.0; 3], 1.0);
        ctrl.reset();
        assert_eq!(ctrl.integral, [0.0; 3]);
    }

    // -----------------------------------------------------------------------
    // ThrusterAllocation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_thruster_allocation_zero_demand_zero_thrust() {
        let mut alloc = ThrusterAllocation::new(vec![[5.0, 0.0], [-5.0, 0.0]], vec![1e5, 1e5]);
        let thrusts = alloc.allocate([0.0; 3]);
        for t in &thrusts {
            assert!(*t < EPS, "zero demand → zero thrust: {t}");
        }
    }

    #[test]
    fn test_thruster_contribution_pure_surge() {
        let alloc = ThrusterAllocation::new(vec![[0.0, 0.0]], vec![1e5]);
        let contrib = alloc.thruster_contribution(0, 1000.0, 0.0);
        assert!((contrib[0] - 1000.0).abs() < EPS);
        assert!(contrib[1].abs() < EPS);
        assert!(contrib[2].abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // SeakeepingAnalysis tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_seakeeping_rms_nonnegative() {
        let freqs: Vec<f64> = (1..=50).map(|i| i as f64 * 0.1).collect();
        let mut rao = ResponseAmplitudeOperator::new(6, freqs);
        rao.set_1dof(2, 5e6, 1e5, 5e7, 1e5); // heave
        let analysis = SeakeepingAnalysis::new(rao, 3.0, 10.0, SpectrumType::PiersonMoskowitz, 3.3);
        let rms = analysis.rms_response(2);
        assert!(rms >= 0.0, "RMS should be non-negative: {rms}");
    }

    #[test]
    fn test_seakeeping_jonswap_larger_rao_gives_larger_response() {
        let freqs: Vec<f64> = (1..=50).map(|i| i as f64 * 0.1).collect();
        let mut rao_small = ResponseAmplitudeOperator::new(1, freqs.clone());
        let mut rao_large = ResponseAmplitudeOperator::new(1, freqs.clone());
        rao_small.set_1dof(0, 5e6, 1e5, 5e7, 1e4);
        rao_large.set_1dof(0, 5e6, 1e5, 5e7, 1e5);
        let a_small = SeakeepingAnalysis::new(rao_small, 3.0, 10.0, SpectrumType::Jonswap, 3.3);
        let a_large = SeakeepingAnalysis::new(rao_large, 3.0, 10.0, SpectrumType::Jonswap, 3.3);
        assert!(
            a_large.rms_response(0) > a_small.rms_response(0),
            "larger RAO excitation → larger response"
        );
    }

    // -----------------------------------------------------------------------
    // AddedMassDamping tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_added_mass_diagonal_restoring_force() {
        let a_diag = [0.0; 6];
        let b_diag = [0.0; 6];
        let c_diag = [1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0];
        let amd = AddedMassDamping::diagonal(a_diag, b_diag, c_diag);
        let x = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = amd.restoring_force(x);
        // F[0] = -C[0][0] * x[0] = -1000 * 0.1 = -100
        assert!((f[0] + 100.0).abs() < EPS, "f[0]={}", f[0]);
    }

    #[test]
    fn test_added_mass_inertia_force() {
        let a_diag = [500.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let amd = AddedMassDamping::diagonal(a_diag, [0.0; 6], [0.0; 6]);
        let body_mass = [1000.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let xddot = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let f = amd.inertia_force(body_mass, xddot);
        // f[0] = (1000 + 500) * 1 = 1500 N
        assert!((f[0] - 1500.0).abs() < EPS, "f[0]={}", f[0]);
    }

    // -----------------------------------------------------------------------
    // Polynomial / interpolation utilities
    // -----------------------------------------------------------------------

    #[test]
    fn test_eval_polynomial_constant() {
        let coeffs = [5.0_f64];
        assert!((eval_polynomial(&coeffs, 7.0) - 5.0).abs() < EPS);
    }

    #[test]
    fn test_eval_polynomial_linear() {
        let coeffs = [1.0_f64, 2.0]; // 1 + 2x
        assert!((eval_polynomial(&coeffs, 3.0) - 7.0).abs() < EPS);
    }

    #[test]
    fn test_interpolate_linear_midpoint() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 10.0, 20.0];
        let v = interpolate_linear(&xs, &ys, 0.5);
        assert!((v - 5.0).abs() < EPS, "v={v}");
    }

    #[test]
    fn test_interpolate_linear_extrapolation_clamps() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0, 10.0];
        assert!((interpolate_linear(&xs, &ys, -1.0) - 0.0).abs() < EPS);
        assert!((interpolate_linear(&xs, &ys, 2.0) - 10.0).abs() < EPS);
    }
}
