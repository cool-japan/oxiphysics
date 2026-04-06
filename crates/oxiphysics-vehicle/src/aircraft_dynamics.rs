#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aircraft flight dynamics and aerodynamics.
//!
//! This module implements a full 6-DOF rigid-body flight dynamics model
//! including:
//!
//! - **6-DOF equations of motion** (body-axis Newton-Euler)
//! - **Aerodynamic coefficients** CL, CD, CY, Cl, Cm, Cn and polar
//! - **Stability derivatives** Clα, Cmα, Cnβ, Clβ, …
//! - **Longitudinal stability** (phugoid and short-period modes)
//! - **Lateral-directional stability** (Dutch roll, spiral, roll modes)
//! - **Control surface effectiveness** (elevator, aileron, rudder)
//! - **Propeller / turbofan thrust** model
//! - **Flight envelope** (Vmin, Vmax, service ceiling)
//! - **Trim conditions** (straight-and-level, climbing)

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// Physical constants
// ============================================================================

/// Standard sea-level air density (kg/m³).
pub const RHO_SL: f64 = 1.225;

/// Standard gravity (m/s²).
pub const G: f64 = 9.80665;

/// Standard sea-level speed of sound (m/s).
pub const A_SL: f64 = 340.294;

/// Gas constant for dry air (J/(kg·K)).
pub const R_AIR: f64 = 287.058;

/// Sea-level standard temperature (K).
pub const T_SL: f64 = 288.15;

/// ISA lapse rate (K/m) — troposphere.
pub const LAPSE_RATE: f64 = 0.0065;

// ============================================================================
// ISA atmosphere model
// ============================================================================

/// International Standard Atmosphere (ISA) properties at a given altitude.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct AtmosphereState {
    /// Altitude above MSL (m).
    pub altitude: f64,
    /// Static temperature (K).
    pub temperature: f64,
    /// Static pressure (Pa).
    pub pressure: f64,
    /// Air density (kg/m³).
    pub density: f64,
    /// Speed of sound (m/s).
    pub speed_of_sound: f64,
    /// Dynamic viscosity (Pa·s).
    pub dynamic_viscosity: f64,
}

impl AtmosphereState {
    /// Compute ISA atmosphere up to the tropopause (11 km).
    ///
    /// Above 11 km the isothermal stratosphere is used.
    pub fn from_altitude(h: f64) -> Self {
        const P_SL: f64 = 101_325.0;
        const H_TROP: f64 = 11_000.0;
        const T_TROP: f64 = T_SL - LAPSE_RATE * H_TROP;
        const GAMMA: f64 = 1.4;

        let (temp, pressure) = if h <= H_TROP {
            let t = T_SL - LAPSE_RATE * h;
            let p = P_SL * (t / T_SL).powf(G / (LAPSE_RATE * R_AIR));
            (t, p)
        } else {
            let p_trop = P_SL * (T_TROP / T_SL).powf(G / (LAPSE_RATE * R_AIR));
            let t = T_TROP;
            let p = p_trop * (-(G / (R_AIR * T_TROP)) * (h - H_TROP)).exp();
            (t, p)
        };

        let density = pressure / (R_AIR * temp);
        let speed_of_sound = (GAMMA * R_AIR * temp).sqrt();
        // Sutherland's formula for dynamic viscosity
        let mu_ref = 1.716e-5;
        let t_ref = 273.15;
        let s = 110.4;
        let dynamic_viscosity = mu_ref * (temp / t_ref).powf(1.5) * (t_ref + s) / (temp + s);

        AtmosphereState {
            altitude: h,
            temperature: temp,
            pressure,
            density,
            speed_of_sound,
            dynamic_viscosity,
        }
    }

    /// Compute Mach number for a given true airspeed (m/s).
    pub fn mach(&self, tas: f64) -> f64 {
        tas / self.speed_of_sound
    }

    /// Compute dynamic pressure q = 0.5 ρ V² (Pa).
    pub fn dynamic_pressure(&self, tas: f64) -> f64 {
        0.5 * self.density * tas * tas
    }
}

// ============================================================================
// Aerodynamic coefficients and polar
// ============================================================================

/// Complete aerodynamic coefficient set for an aircraft.
///
/// Coefficients are dimensionless and referenced to the wing reference area
/// `S_ref` and mean aerodynamic chord `c_mac`.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct AeroCoefficients {
    /// Lift coefficient CL.
    pub cl: f64,
    /// Drag coefficient CD.
    pub cd: f64,
    /// Side-force coefficient CY.
    pub cy: f64,
    /// Rolling-moment coefficient Cl (sometimes written ℓ).
    pub c_roll: f64,
    /// Pitching-moment coefficient Cm.
    pub cm: f64,
    /// Yawing-moment coefficient Cn.
    pub cn: f64,
}

impl AeroCoefficients {
    /// Create a zero coefficient set.
    pub fn zero() -> Self {
        AeroCoefficients {
            cl: 0.0,
            cd: 0.0,
            cy: 0.0,
            c_roll: 0.0,
            cm: 0.0,
            cn: 0.0,
        }
    }
}

/// Lift-drag polar model.
///
/// Uses the classical parabolic polar:
/// `CD = CD0 + CL² / (π e AR)`
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct LiftDragPolar {
    /// Zero-lift drag coefficient CD0.
    pub cd0: f64,
    /// Oswald span efficiency factor (0 < e ≤ 1).
    pub oswald_e: f64,
    /// Wing aspect ratio AR = b² / S.
    pub aspect_ratio: f64,
    /// CL at zero angle of attack.
    pub cl0: f64,
    /// CL-α slope (per radian).
    pub cl_alpha: f64,
    /// Angle of attack at zero lift (rad).
    pub alpha_zl: f64,
    /// Maximum lift coefficient (stall).
    pub cl_max: f64,
}

impl LiftDragPolar {
    /// Create a typical subsonic transport polar.
    pub fn typical_transport() -> Self {
        LiftDragPolar {
            cd0: 0.020,
            oswald_e: 0.80,
            aspect_ratio: 8.0,
            cl0: 0.15,
            cl_alpha: 5.5,
            alpha_zl: -0.027,
            cl_max: 1.4,
        }
    }

    /// Compute CL for a given angle of attack (rad).
    ///
    /// Returns `cl_max` if the computed value exceeds it (stall limit).
    pub fn cl_from_alpha(&self, alpha: f64) -> f64 {
        let cl = self.cl0 + self.cl_alpha * (alpha - self.alpha_zl);
        cl.min(self.cl_max)
    }

    /// Compute induced drag coefficient from CL.
    ///
    /// `CDi = CL² / (π e AR)`
    pub fn induced_drag(&self, cl: f64) -> f64 {
        cl * cl / (PI * self.oswald_e * self.aspect_ratio)
    }

    /// Compute total drag coefficient.
    pub fn cd_from_cl(&self, cl: f64) -> f64 {
        self.cd0 + self.induced_drag(cl)
    }

    /// Lift-to-drag ratio.
    pub fn ld_ratio(&self, alpha: f64) -> f64 {
        let cl = self.cl_from_alpha(alpha);
        let cd = self.cd_from_cl(cl);
        if cd.abs() < 1e-15 {
            return 0.0;
        }
        cl / cd
    }

    /// Angle of attack for maximum L/D (rad).
    ///
    /// Occurs where `dCL/dα * CD = CL * dCD/dα`, solved analytically.
    pub fn alpha_at_max_ld(&self) -> f64 {
        // At max L/D: CL = sqrt(CD0 * π e AR)
        let cl_opt = (self.cd0 * PI * self.oswald_e * self.aspect_ratio).sqrt();
        (cl_opt - self.cl0) / self.cl_alpha + self.alpha_zl
    }
}

// ============================================================================
// Stability derivatives
// ============================================================================

/// Longitudinal stability derivatives.
///
/// All derivatives are per radian unless noted.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct LongitudinalDerivatives {
    /// CL-α: lift-curve slope (per rad).
    pub cl_alpha: f64,
    /// CD-α: drag-α derivative (per rad).
    pub cd_alpha: f64,
    /// Cm-α: pitch-stiffness derivative (per rad). Negative for stable.
    pub cm_alpha: f64,
    /// Cm-q: pitch-damping derivative (per rad/s, non-dimensionalised by c/2V).
    pub cm_q: f64,
    /// Cm-α_dot: pitch-rate-change derivative (per rad/s).
    pub cm_alphadot: f64,
    /// CL-q: lift due to pitch rate (per rad/s).
    pub cl_q: f64,
    /// Cm-δe: elevator effectiveness (per rad deflection).
    pub cm_de: f64,
    /// CL-δe: lift due to elevator (per rad).
    pub cl_de: f64,
}

impl LongitudinalDerivatives {
    /// Typical subsonic transport jet derivatives.
    pub fn typical_transport() -> Self {
        LongitudinalDerivatives {
            cl_alpha: 5.5,
            cd_alpha: 0.3,
            cm_alpha: -1.2,
            cm_q: -12.0,
            cm_alphadot: -4.0,
            cl_q: 3.5,
            cm_de: -1.5,
            cl_de: 0.4,
        }
    }
}

/// Lateral-directional stability derivatives.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct LateralDerivatives {
    /// Cy-β: side-force due to sideslip (per rad). Usually negative.
    pub cy_beta: f64,
    /// Cl-β: dihedral effect — roll due to sideslip (per rad). Negative for stable.
    pub cl_beta: f64,
    /// Cn-β: directional stability (per rad). Positive for stable.
    pub cn_beta: f64,
    /// Cl-p: roll damping (per rad/s, non-dim by b/2V).
    pub cl_p: f64,
    /// Cn-p: adverse yaw due to roll (per rad/s).
    pub cn_p: f64,
    /// Cl-r: roll due to yaw rate (per rad/s).
    pub cl_r: f64,
    /// Cn-r: yaw damping (per rad/s). Negative for stable.
    pub cn_r: f64,
    /// Cl-δa: aileron effectiveness (per rad).
    pub cl_da: f64,
    /// Cn-δa: aileron adverse yaw (per rad).
    pub cn_da: f64,
    /// Cy-δr: side force due to rudder (per rad).
    pub cy_dr: f64,
    /// Cl-δr: roll due to rudder (per rad).
    pub cl_dr: f64,
    /// Cn-δr: directional control from rudder (per rad). Negative.
    pub cn_dr: f64,
}

impl LateralDerivatives {
    /// Typical subsonic transport jet lateral derivatives.
    pub fn typical_transport() -> Self {
        LateralDerivatives {
            cy_beta: -0.60,
            cl_beta: -0.10,
            cn_beta: 0.12,
            cl_p: -0.45,
            cn_p: -0.05,
            cl_r: 0.08,
            cn_r: -0.15,
            cl_da: 0.18,
            cn_da: -0.03,
            cy_dr: 0.17,
            cl_dr: 0.02,
            cn_dr: -0.10,
        }
    }
}

// ============================================================================
// Longitudinal stability modes
// ============================================================================

/// Parameters of a second-order dynamic mode (complex eigenvalue pair).
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct DynamicMode {
    /// Natural frequency (rad/s).
    pub omega_n: f64,
    /// Damping ratio ζ (dimensionless).
    pub zeta: f64,
    /// Mode name.
    pub name: &'static str,
}

impl DynamicMode {
    /// Period of oscillation (s). Returns `f64::INFINITY` for overdamped modes.
    pub fn period(&self) -> f64 {
        if self.zeta >= 1.0 {
            return f64::INFINITY;
        }
        let omega_d = self.omega_n * (1.0 - self.zeta * self.zeta).sqrt();
        2.0 * PI / omega_d
    }

    /// Time to half-amplitude (s). Positive means stable.
    pub fn t_half(&self) -> f64 {
        if self.zeta * self.omega_n < 1e-15 {
            return f64::INFINITY;
        }
        0.693 / (self.zeta * self.omega_n)
    }

    /// Returns true if this mode is stable (ζ > 0).
    pub fn is_stable(&self) -> bool {
        self.zeta > 0.0
    }
}

/// Longitudinal stability analysis.
///
/// Provides analytical approximations for phugoid and short-period modes.
#[allow(dead_code)]
pub struct LongitudinalStability {
    /// Flight speed (m/s TAS).
    pub v0: f64,
    /// Air density (kg/m³).
    pub rho: f64,
    /// Wing area (m²).
    pub s_ref: f64,
    /// Mean aerodynamic chord (m).
    pub c_mac: f64,
    /// Aircraft total mass (kg).
    pub mass: f64,
    /// Pitch moment of inertia Iyy (kg·m²).
    pub i_yy: f64,
    /// Longitudinal derivatives.
    pub derivs: LongitudinalDerivatives,
    /// Polar model.
    pub polar: LiftDragPolar,
}

impl LongitudinalStability {
    /// Create a longitudinal stability model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        v0: f64,
        rho: f64,
        s_ref: f64,
        c_mac: f64,
        mass: f64,
        i_yy: f64,
        derivs: LongitudinalDerivatives,
        polar: LiftDragPolar,
    ) -> Self {
        LongitudinalStability {
            v0,
            rho,
            s_ref,
            c_mac,
            mass,
            i_yy,
            derivs,
            polar,
        }
    }

    /// Phugoid mode approximation (Lanchester).
    ///
    /// `ω_ph ≈ √2 * g / V₀`
    /// `ζ_ph ≈ CD0 / (√2 * CL0)`
    pub fn phugoid_mode(&self, alpha0: f64) -> DynamicMode {
        let cl0 = self.polar.cl_from_alpha(alpha0);
        let cd0 = self.polar.cd_from_cl(cl0);
        let omega_n = (2.0_f64).sqrt() * G / self.v0;
        let zeta = if cl0.abs() > 1e-6 {
            cd0 / ((2.0_f64).sqrt() * cl0)
        } else {
            0.0
        };
        DynamicMode {
            omega_n,
            zeta,
            name: "phugoid",
        }
    }

    /// Short-period mode approximation.
    ///
    /// Uses the simplified expression from aircraft dynamics texts.
    pub fn short_period_mode(&self, alpha0: f64) -> DynamicMode {
        let q_bar = 0.5 * self.rho * self.v0 * self.v0;
        let cl = self.polar.cl_from_alpha(alpha0);
        let _mu = self.mass / (self.rho * self.s_ref * self.c_mac);

        // ω_sp² ≈ q̄·S·c / (2·mu·I_yy) * (Cmα + CL_alpha * …)
        let factor = q_bar * self.s_ref * self.c_mac / (self.i_yy);
        let omega_sq = -factor * self.derivs.cm_alpha; // cm_alpha < 0 for stable
        let omega_n = if omega_sq > 0.0 { omega_sq.sqrt() } else { 0.0 };

        // ζ_sp from pitch damping Cmq and Cm_alphadot
        let zeta_num =
            -factor * (self.derivs.cm_q + self.derivs.cm_alphadot) * self.c_mac / (2.0 * self.v0);
        let zeta = if omega_n > 1e-6 {
            zeta_num / (2.0 * omega_n)
        } else {
            0.0
        };
        let _ = cl; // used in more detailed formulae; keep for doc

        DynamicMode {
            omega_n,
            zeta,
            name: "short_period",
        }
    }
}

// ============================================================================
// Lateral-directional stability modes
// ============================================================================

/// Lateral-directional stability analysis.
///
/// Provides analytical approximations for Dutch roll, spiral, and roll modes.
#[allow(dead_code)]
pub struct LateralStability {
    /// Flight speed (m/s TAS).
    pub v0: f64,
    /// Air density (kg/m³).
    pub rho: f64,
    /// Wing area (m²).
    pub s_ref: f64,
    /// Wing span (m).
    pub b_ref: f64,
    /// Aircraft total mass (kg).
    pub mass: f64,
    /// Roll moment of inertia Ixx (kg·m²).
    pub i_xx: f64,
    /// Yaw moment of inertia Izz (kg·m²).
    pub i_zz: f64,
    /// Product of inertia Ixz (kg·m²).
    pub i_xz: f64,
    /// Lateral-directional derivatives.
    pub derivs: LateralDerivatives,
}

impl LateralStability {
    /// Create a lateral stability model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        v0: f64,
        rho: f64,
        s_ref: f64,
        b_ref: f64,
        mass: f64,
        i_xx: f64,
        i_zz: f64,
        i_xz: f64,
        derivs: LateralDerivatives,
    ) -> Self {
        LateralStability {
            v0,
            rho,
            s_ref,
            b_ref,
            mass,
            i_xx,
            i_zz,
            i_xz,
            derivs,
        }
    }

    /// Non-dimensionalisation factor for lateral rates.
    fn q_bar(&self) -> f64 {
        0.5 * self.rho * self.v0 * self.v0
    }

    /// Roll mode time constant (s).
    ///
    /// `T_roll = -I_xx / (q̄ S b² Cl_p / (2V))`
    pub fn roll_mode_time_constant(&self) -> f64 {
        let qsb = self.q_bar() * self.s_ref * self.b_ref;
        let denom = qsb * self.b_ref * self.derivs.cl_p / (2.0 * self.v0);
        if denom.abs() < 1e-15 {
            return f64::INFINITY;
        }
        -self.i_xx / denom
    }

    /// Dutch roll mode (approximate second-order system).
    pub fn dutch_roll_mode(&self) -> DynamicMode {
        let qsb = self.q_bar() * self.s_ref * self.b_ref;
        let beta_factor = qsb / (self.mass * self.v0);
        let n_beta = qsb * self.derivs.cn_beta / self.i_zz;

        // ω_dr² ≈ N_β (yaw stiffness)
        let omega_sq = n_beta.max(0.0);
        let omega_n = omega_sq.sqrt();

        // ζ_dr from yaw damping Nr and roll coupling
        let n_r = qsb * self.b_ref * self.derivs.cn_r / (2.0 * self.v0 * self.i_zz);
        let zeta = if omega_n > 1e-6 {
            -n_r / (2.0 * omega_n)
        } else {
            0.0
        };
        let _ = beta_factor; // reserved for complete 4-DOF analysis

        DynamicMode {
            omega_n,
            zeta,
            name: "dutch_roll",
        }
    }

    /// Spiral mode time constant (s).
    ///
    /// Positive → divergent (unstable spiral); negative → convergent.
    pub fn spiral_mode_time_constant(&self) -> f64 {
        let qsb = self.q_bar() * self.s_ref * self.b_ref;
        let l_beta = qsb * self.derivs.cl_beta / self.i_xx;
        let n_beta = qsb * self.derivs.cn_beta / self.i_zz;
        let l_r = qsb * self.b_ref * self.derivs.cl_r / (2.0 * self.v0 * self.i_xx);
        let n_r = qsb * self.b_ref * self.derivs.cn_r / (2.0 * self.v0 * self.i_zz);

        // 1/T_spiral ≈ (L_β N_r - N_β L_r) / (L_β - N_β * something)
        // Simplified first-order approximation:
        let denom = l_beta * n_r - n_beta * l_r;
        if denom.abs() < 1e-15 {
            return f64::INFINITY;
        }
        (l_beta - n_beta) / denom
    }
}

// ============================================================================
// Control surface effectiveness
// ============================================================================

/// Control surface (elevator, aileron, rudder) model.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct ControlSurface {
    /// Surface name.
    pub name: &'static str,
    /// Maximum positive deflection (rad).
    pub deflection_max: f64,
    /// Minimum (most negative) deflection (rad).
    pub deflection_min: f64,
    /// Hinge-moment coefficient Chδ (per rad deflection).
    pub ch_delta: f64,
    /// Control moment coefficient per unit deflection (per rad).
    pub cm_delta: f64,
}

impl ControlSurface {
    /// Create an elevator definition.
    pub fn elevator(cm_de: f64) -> Self {
        ControlSurface {
            name: "elevator",
            deflection_max: 25.0_f64.to_radians(),
            deflection_min: -25.0_f64.to_radians(),
            ch_delta: -0.35,
            cm_delta: cm_de,
        }
    }

    /// Create an aileron definition.
    pub fn aileron(cl_da: f64) -> Self {
        ControlSurface {
            name: "aileron",
            deflection_max: 20.0_f64.to_radians(),
            deflection_min: -20.0_f64.to_radians(),
            ch_delta: -0.28,
            cm_delta: cl_da,
        }
    }

    /// Create a rudder definition.
    pub fn rudder(cn_dr: f64) -> Self {
        ControlSurface {
            name: "rudder",
            deflection_max: 30.0_f64.to_radians(),
            deflection_min: -30.0_f64.to_radians(),
            ch_delta: -0.22,
            cm_delta: cn_dr,
        }
    }

    /// Clamp a requested deflection to the travel limits.
    pub fn clamp(&self, delta: f64) -> f64 {
        delta.clamp(self.deflection_min, self.deflection_max)
    }

    /// Moment coefficient contribution for a given deflection.
    pub fn moment_coefficient(&self, delta: f64) -> f64 {
        self.cm_delta * self.clamp(delta)
    }

    /// Hinge moment coefficient for a given deflection.
    pub fn hinge_moment(&self, delta: f64) -> f64 {
        self.ch_delta * self.clamp(delta)
    }
}

// ============================================================================
// Propeller / turbofan thrust model
// ============================================================================

/// Propeller thrust model (momentum theory + Mach correction).
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct PropellerThrustModel {
    /// Diameter of the propeller (m).
    pub diameter: f64,
    /// Static thrust coefficient CT0 at zero advance ratio.
    pub ct0: f64,
    /// Advance-ratio coefficient slope (dCT/dJ).
    pub ct_slope: f64,
    /// Maximum RPM.
    pub rpm_max: f64,
    /// Propulsive efficiency at cruise.
    pub eta_p: f64,
}

impl PropellerThrustModel {
    /// Advance ratio J = V / (n D).
    pub fn advance_ratio(&self, v: f64, rpm: f64) -> f64 {
        let n = rpm / 60.0;
        if n < 1e-6 {
            return 0.0;
        }
        v / (n * self.diameter)
    }

    /// Thrust coefficient CT(J).
    pub fn thrust_coefficient(&self, j: f64) -> f64 {
        (self.ct0 + self.ct_slope * j).max(0.0)
    }

    /// Thrust force (N).
    ///
    /// `T = CT · ρ · n² · D⁴`
    pub fn thrust(&self, v: f64, rpm: f64, rho: f64) -> f64 {
        let n = rpm / 60.0;
        let j = self.advance_ratio(v, rpm);
        let ct = self.thrust_coefficient(j);
        ct * rho * n * n * self.diameter.powi(4)
    }

    /// Power absorbed (W).
    pub fn power(&self, v: f64, rpm: f64, rho: f64) -> f64 {
        let t = self.thrust(v, rpm, rho);
        if self.eta_p < 1e-6 {
            return 0.0;
        }
        t * v / self.eta_p
    }
}

/// Turbofan thrust model.
///
/// Models net thrust as a function of Mach number and altitude.
#[allow(dead_code)]
pub struct TurbofanThrustModelImpl {
    /// Sea-level static thrust (N).
    pub t_sls: f64,
    /// Bypass ratio (for fan thrust lapse approximation).
    pub bypass_ratio: f64,
    /// Thrust lapse exponent with density ratio σ.
    pub lapse_exponent: f64,
    /// Mach number penalty coefficient.
    pub mach_penalty: f64,
}

impl TurbofanThrustModelImpl {
    /// Create a turbofan with typical civil engine parameters.
    pub fn new(t_sls: f64) -> Self {
        TurbofanThrustModelImpl {
            t_sls,
            bypass_ratio: 8.0,
            lapse_exponent: 0.9,
            mach_penalty: 0.14,
        }
    }

    /// Net thrust at altitude and Mach number.
    ///
    /// Simplified lapse model:
    /// `T = T_SLS · σ^n · (1 - k_M · M)`
    pub fn net_thrust(&self, atm: &AtmosphereState, mach: f64) -> f64 {
        let sigma = atm.density / RHO_SL;
        let t_alt = self.t_sls * sigma.powf(self.lapse_exponent);
        t_alt * (1.0 - self.mach_penalty * mach).max(0.0)
    }

    /// Specific fuel consumption (kg/(N·s)) — constant approximation.
    pub fn sfc(&self) -> f64 {
        // Typical bypass ratio 8 engine SFC ≈ 1.5e-5 kg/(N·s)
        let base = 1.5e-5;
        // Lower BPR = higher SFC (roughly)
        base * (1.0 + 0.05 * (8.0 - self.bypass_ratio).max(0.0))
    }

    /// Fuel flow rate (kg/s) at given thrust setting.
    pub fn fuel_flow(&self, thrust: f64) -> f64 {
        self.sfc() * thrust.abs()
    }
}

// ============================================================================
// Flight envelope
// ============================================================================

/// Flight envelope boundaries for an aircraft.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct FlightEnvelope {
    /// Wing area (m²).
    pub s_ref: f64,
    /// Maximum take-off weight (N).
    pub mtow: f64,
    /// Lift-drag polar.
    pub polar: LiftDragPolar,
    /// Maximum level-flight thrust (N) at sea level static.
    pub t_max_sls: f64,
    /// Structural dive speed (EAS, m/s).
    pub v_d: f64,
    /// Design load factor limits \[n_min, n_max\].
    pub load_factor_limits: [f64; 2],
}

impl FlightEnvelope {
    /// Minimum speed (stall speed) at a given altitude and weight (N).
    ///
    /// `V_min = sqrt(2W / (ρ S CL_max))`
    pub fn v_min(&self, atm: &AtmosphereState, weight: f64) -> f64 {
        (2.0 * weight / (atm.density * self.s_ref * self.polar.cl_max)).sqrt()
    }

    /// Maximum level-flight speed at a given altitude.
    ///
    /// Iterates to find V where T_available = D(V).
    pub fn v_max(&self, atm: &AtmosphereState, thrust: f64) -> f64 {
        // Solve: 0.5 ρ V² S CD(V) = T
        // CD ≈ CD0 + CL² / (π e AR), CL = W / (0.5 ρ V² S)
        // This is a quartic in V; we solve numerically via bisection.
        let weight = self.mtow;
        let rho = atm.density;
        let s = self.s_ref;
        let cd0 = self.polar.cd0;
        let k = 1.0 / (PI * self.polar.oswald_e * self.polar.aspect_ratio);

        let drag = |v: f64| -> f64 {
            let q = 0.5 * rho * v * v;
            let cl = if q > 1e-6 { weight / (q * s) } else { 10.0 };
            let cd = cd0 + k * cl * cl;
            q * s * cd
        };

        let mut lo = 10.0_f64;
        let mut hi = self.v_d;
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if drag(mid) < thrust {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Service ceiling (m) — altitude where rate-of-climb < 0.508 m/s (100 ft/min).
    ///
    /// Uses ISA atmosphere and assumes thrust lapse with density.
    pub fn service_ceiling(&self, t_sls: f64, weight: f64) -> f64 {
        let target_roc = 0.508;
        let s = self.s_ref;
        let cd0 = self.polar.cd0;
        let k = 1.0 / (PI * self.polar.oswald_e * self.polar.aspect_ratio);

        let roc_at_alt = |h: f64| -> f64 {
            let atm = AtmosphereState::from_altitude(h);
            let sigma = atm.density / RHO_SL;
            let thrust = t_sls * sigma.powf(0.75);
            // Best ROC speed: V = sqrt(2W/(ρS) * sqrt(k/3CD0))
            let v = (2.0 * weight / (atm.density * s) * (k / (3.0 * cd0)).sqrt()).sqrt();
            let q = 0.5 * atm.density * v * v;
            let cl = weight / (q * s);
            let cd = cd0 + k * cl * cl;
            let drag = q * s * cd;
            let excess = thrust - drag;
            let _ = excess;
            (thrust - drag) * v / weight
        };

        let mut lo = 0.0_f64;
        let mut hi = 15_000.0_f64;
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if roc_at_alt(mid) > target_roc {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ============================================================================
// 6-DOF equations of motion
// ============================================================================

/// Full 6-DOF aircraft state in body axes.
///
/// Positions are in Earth (NED) frame; velocities, angular rates, and
/// orientations are in body frame.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct AircraftState6Dof {
    /// Position in NED frame (m): \[North, East, Down\].
    pub position_ned: [f64; 3],
    /// Body-axis velocity \[u, v, w\] (m/s).
    pub velocity_body: [f64; 3],
    /// Euler angles \[φ, θ, ψ\] (roll, pitch, yaw) in rad.
    pub euler_angles: [f64; 3],
    /// Body-axis angular rates \[p, q, r\] (rad/s).
    pub angular_rates: [f64; 3],
}

impl AircraftState6Dof {
    /// Create a level-flight initial state at speed `u0` (m/s).
    pub fn level_flight(u0: f64) -> Self {
        AircraftState6Dof {
            position_ned: [0.0, 0.0, 0.0],
            velocity_body: [u0, 0.0, 0.0],
            euler_angles: [0.0, 0.0, 0.0],
            angular_rates: [0.0, 0.0, 0.0],
        }
    }

    /// True airspeed magnitude (m/s).
    pub fn tas(&self) -> f64 {
        let [u, v, w] = self.velocity_body;
        (u * u + v * v + w * w).sqrt()
    }

    /// Angle of attack α (rad).
    pub fn alpha(&self) -> f64 {
        let [u, _v, w] = self.velocity_body;
        w.atan2(u)
    }

    /// Sideslip angle β (rad).
    pub fn beta(&self) -> f64 {
        let v_tot = self.tas();
        if v_tot < 1e-6 {
            return 0.0;
        }
        let [_u, v, _w] = self.velocity_body;
        (v / v_tot).asin()
    }

    /// Dynamic pressure (Pa).
    pub fn dynamic_pressure(&self, rho: f64) -> f64 {
        let v = self.tas();
        0.5 * rho * v * v
    }

    /// Rotation matrix from body to NED (3×3, row-major).
    ///
    /// Uses the standard ZYX Euler sequence.
    pub fn body_to_ned(&self) -> [[f64; 3]; 3] {
        let [phi, theta, psi] = self.euler_angles;
        let (sp, cp) = (phi.sin(), phi.cos());
        let (st, ct) = (theta.sin(), theta.cos());
        let (ss, cs) = (psi.sin(), psi.cos());
        [
            [ct * cs, sp * st * cs - cp * ss, cp * st * cs + sp * ss],
            [ct * ss, sp * st * ss + cp * cs, cp * st * ss - sp * cs],
            [-st, sp * ct, cp * ct],
        ]
    }
}

/// 6-DOF aircraft equations of motion integrator.
///
/// Integrates Newton-Euler EOM in body axes using 4th-order Runge-Kutta.
#[allow(dead_code)]
pub struct SixDofIntegrator {
    /// Aircraft total mass (kg).
    pub mass: f64,
    /// Inertia tensor diagonal \[Ixx, Iyy, Izz\] (kg·m²).
    pub inertia: [f64; 3],
    /// Product of inertia Ixz (kg·m²).
    pub i_xz: f64,
    /// Wing reference area (m²).
    pub s_ref: f64,
    /// Wing span (m).
    pub b_ref: f64,
    /// Mean aerodynamic chord (m).
    pub c_mac: f64,
    /// Longitudinal derivatives.
    pub long_derivs: LongitudinalDerivatives,
    /// Lateral derivatives.
    pub lat_derivs: LateralDerivatives,
    /// Drag polar.
    pub polar: LiftDragPolar,
}

impl SixDofIntegrator {
    /// Compute aerodynamic force and moment coefficients at the current state.
    #[allow(clippy::too_many_arguments)]
    pub fn aero_coefficients(
        &self,
        alpha: f64,
        beta: f64,
        _p: f64,
        _q: f64,
        _r: f64,
        delta_e: f64,
        delta_a: f64,
        delta_r: f64,
    ) -> AeroCoefficients {
        let cl = self.polar.cl_from_alpha(alpha) + self.long_derivs.cl_de * delta_e;
        let cd = self.polar.cd_from_cl(cl);
        let cy = self.lat_derivs.cy_beta * beta + self.lat_derivs.cy_dr * delta_r;
        let c_roll = self.lat_derivs.cl_beta * beta
            + self.lat_derivs.cl_da * delta_a
            + self.lat_derivs.cl_dr * delta_r;
        let cm = self.long_derivs.cm_alpha * alpha + self.long_derivs.cm_de * delta_e;
        let cn = self.lat_derivs.cn_beta * beta
            + self.lat_derivs.cn_da * delta_a
            + self.lat_derivs.cn_dr * delta_r;

        AeroCoefficients {
            cl,
            cd,
            cy,
            c_roll,
            cm,
            cn,
        }
    }

    /// Integrate one RK4 step.
    #[allow(clippy::too_many_arguments)]
    pub fn rk4_step(
        &self,
        state: &AircraftState6Dof,
        thrust: f64,
        delta_e: f64,
        delta_a: f64,
        delta_r: f64,
        rho: f64,
        dt: f64,
    ) -> AircraftState6Dof {
        let f = |s: &AircraftState6Dof| self.derivatives(s, thrust, delta_e, delta_a, delta_r, rho);

        let k1 = f(state);
        let s2 = self.advance(state, &k1, 0.5 * dt);
        let k2 = f(&s2);
        let s3 = self.advance(state, &k2, 0.5 * dt);
        let k3 = f(&s3);
        let s4 = self.advance(state, &k3, dt);
        let k4 = f(&s4);

        // Combine
        let mut new = state.clone();
        for i in 0..3 {
            new.position_ned[i] += dt / 6.0
                * (k1.position_ned[i]
                    + 2.0 * k2.position_ned[i]
                    + 2.0 * k3.position_ned[i]
                    + k4.position_ned[i]);
            new.velocity_body[i] += dt / 6.0
                * (k1.velocity_body[i]
                    + 2.0 * k2.velocity_body[i]
                    + 2.0 * k3.velocity_body[i]
                    + k4.velocity_body[i]);
            new.euler_angles[i] += dt / 6.0
                * (k1.euler_angles[i]
                    + 2.0 * k2.euler_angles[i]
                    + 2.0 * k3.euler_angles[i]
                    + k4.euler_angles[i]);
            new.angular_rates[i] += dt / 6.0
                * (k1.angular_rates[i]
                    + 2.0 * k2.angular_rates[i]
                    + 2.0 * k3.angular_rates[i]
                    + k4.angular_rates[i]);
        }
        new
    }

    /// Compute time derivatives of the full state vector (body-axis EOM).
    #[allow(clippy::too_many_arguments)]
    fn derivatives(
        &self,
        s: &AircraftState6Dof,
        thrust: f64,
        delta_e: f64,
        delta_a: f64,
        delta_r: f64,
        rho: f64,
    ) -> AircraftState6Dof {
        let [u, v, w] = s.velocity_body;
        let [p, q, r] = s.angular_rates;
        let [phi, theta, _psi] = s.euler_angles;

        let alpha = s.alpha();
        let beta = s.beta();
        let q_bar = s.dynamic_pressure(rho);
        let v_ref = self.s_ref;
        let coeffs = self.aero_coefficients(alpha, beta, p, q, r, delta_e, delta_a, delta_r);

        // Forces in stability axes → body axes
        let f_aero_x = q_bar * v_ref * (-coeffs.cd * alpha.cos() + coeffs.cl * alpha.sin());
        let f_aero_y = q_bar * v_ref * coeffs.cy;
        let f_aero_z = q_bar * v_ref * (-coeffs.cd * alpha.sin() - coeffs.cl * alpha.cos());

        // Gravity in body axes
        let gx = -G * theta.sin();
        let gy = G * theta.cos() * phi.sin();
        let gz = G * theta.cos() * phi.cos();

        // Thrust along body x-axis
        let fx = f_aero_x + thrust + self.mass * gx;
        let fy = f_aero_y + self.mass * gy;
        let fz = f_aero_z + self.mass * gz;

        // Translational accelerations
        let u_dot = fx / self.mass + r * v - q * w;
        let v_dot = fy / self.mass - r * u + p * w;
        let w_dot = fz / self.mass + q * u - p * v;

        // Moments
        let qsc = q_bar * v_ref * self.c_mac;
        let qsb = q_bar * v_ref * self.b_ref;
        let l_m = qsb * coeffs.c_roll;
        let m_m = qsc * coeffs.cm;
        let n_m = qsb * coeffs.cn;

        // Angular acceleration (simplified: ignore Ixz cross terms)
        let [ixx, iyy, izz] = self.inertia;
        let p_dot = (l_m + (iyy - izz) * q * r) / ixx;
        let q_dot = (m_m + (izz - ixx) * p * r) / iyy;
        let r_dot = (n_m + (ixx - iyy) * p * q) / izz;

        // Kinematic equations (Euler angle rates)
        let phi_dot = p + (q * phi.sin() + r * phi.cos()) * theta.tan();
        let theta_dot = q * phi.cos() - r * phi.sin();
        let psi_dot = (q * phi.sin() + r * phi.cos()) / theta.cos().max(1e-6);

        // NED velocity (rotate body velocity to NED)
        let dcm = s.body_to_ned();
        let vel_body = [u, v, w];
        let pos_dot = mat3_mul_v3(dcm, vel_body);

        AircraftState6Dof {
            position_ned: pos_dot,
            velocity_body: [u_dot, v_dot, w_dot],
            euler_angles: [phi_dot, theta_dot, psi_dot],
            angular_rates: [p_dot, q_dot, r_dot],
        }
    }

    /// Advance state by step `h` given derivatives `d`.
    fn advance(&self, s: &AircraftState6Dof, d: &AircraftState6Dof, h: f64) -> AircraftState6Dof {
        let add = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
            [a[0] + h * b[0], a[1] + h * b[1], a[2] + h * b[2]]
        };
        AircraftState6Dof {
            position_ned: add(s.position_ned, d.position_ned),
            velocity_body: add(s.velocity_body, d.velocity_body),
            euler_angles: add(s.euler_angles, d.euler_angles),
            angular_rates: add(s.angular_rates, d.angular_rates),
        }
    }
}

fn mat3_mul_v3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ============================================================================
// Trim conditions
// ============================================================================

/// Trim solution for straight-and-level flight.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct TrimCondition {
    /// True airspeed (m/s).
    pub v_tas: f64,
    /// Altitude (m).
    pub altitude: f64,
    /// Trim angle of attack (rad).
    pub alpha_trim: f64,
    /// Trim elevator deflection (rad).
    pub delta_e_trim: f64,
    /// Thrust required (N).
    pub thrust_required: f64,
    /// Lift coefficient at trim.
    pub cl_trim: f64,
    /// Drag coefficient at trim.
    pub cd_trim: f64,
    /// L/D ratio at trim.
    pub ld_trim: f64,
}

impl TrimCondition {
    /// Compute straight-and-level trim.
    ///
    /// Finds α such that `L = W` and `Cm = 0` simultaneously.
    pub fn straight_level(
        v_tas: f64,
        altitude: f64,
        weight: f64,
        s_ref: f64,
        polar: &LiftDragPolar,
        long_derivs: &LongitudinalDerivatives,
    ) -> Self {
        let atm = AtmosphereState::from_altitude(altitude);
        let q_bar = atm.dynamic_pressure(v_tas);

        // CL required for level flight
        let cl_trim = weight / (q_bar * s_ref);
        let cl_trim = cl_trim.min(polar.cl_max);

        // Alpha from CL
        let alpha_trim = (cl_trim - polar.cl0) / polar.cl_alpha + polar.alpha_zl;

        // Elevator for Cm = 0: delta_e = -Cm_alpha * alpha / Cm_de
        let delta_e_trim = if long_derivs.cm_de.abs() > 1e-10 {
            -long_derivs.cm_alpha * alpha_trim / long_derivs.cm_de
        } else {
            0.0
        };

        let cd_trim = polar.cd_from_cl(cl_trim);
        let thrust_required = q_bar * s_ref * cd_trim;
        let ld_trim = if cd_trim > 1e-10 {
            cl_trim / cd_trim
        } else {
            0.0
        };

        TrimCondition {
            v_tas,
            altitude,
            alpha_trim,
            delta_e_trim,
            thrust_required,
            cl_trim,
            cd_trim,
            ld_trim,
        }
    }

    /// Compute climbing trim at a given flight-path angle γ (rad).
    ///
    /// Weight component along flight path must be overcome by thrust excess.
    pub fn climbing(
        v_tas: f64,
        altitude: f64,
        weight: f64,
        s_ref: f64,
        gamma: f64,
        polar: &LiftDragPolar,
        long_derivs: &LongitudinalDerivatives,
    ) -> Self {
        // In a climb, L = W cos γ
        let atm = AtmosphereState::from_altitude(altitude);
        let q_bar = atm.dynamic_pressure(v_tas);
        let cl_trim = (weight * gamma.cos()) / (q_bar * s_ref);
        let cl_trim = cl_trim.min(polar.cl_max);
        let alpha_trim = (cl_trim - polar.cl0) / polar.cl_alpha + polar.alpha_zl + gamma;
        let delta_e_trim = if long_derivs.cm_de.abs() > 1e-10 {
            -long_derivs.cm_alpha * alpha_trim / long_derivs.cm_de
        } else {
            0.0
        };
        let cd_trim = polar.cd_from_cl(cl_trim);
        let thrust_required = q_bar * s_ref * cd_trim + weight * gamma.sin();
        let ld_trim = if cd_trim > 1e-10 {
            cl_trim / cd_trim
        } else {
            0.0
        };

        TrimCondition {
            v_tas,
            altitude,
            alpha_trim,
            delta_e_trim,
            thrust_required,
            cl_trim,
            cd_trim,
            ld_trim,
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // -------------------------------------------------------------------------
    // 1. ISA atmosphere: sea-level temperature is standard
    // -------------------------------------------------------------------------
    #[test]
    fn test_isa_sea_level_temp() {
        let atm = AtmosphereState::from_altitude(0.0);
        assert!(
            (atm.temperature - T_SL).abs() < 0.01,
            "T(0) = {}",
            atm.temperature
        );
    }

    // -------------------------------------------------------------------------
    // 2. ISA atmosphere: sea-level pressure
    // -------------------------------------------------------------------------
    #[test]
    fn test_isa_sea_level_pressure() {
        let atm = AtmosphereState::from_altitude(0.0);
        assert!(
            (atm.pressure - 101_325.0).abs() < 1.0,
            "P(0) = {}",
            atm.pressure
        );
    }

    // -------------------------------------------------------------------------
    // 3. ISA atmosphere: density decreases with altitude
    // -------------------------------------------------------------------------
    #[test]
    fn test_isa_density_decreases() {
        let atm0 = AtmosphereState::from_altitude(0.0);
        let atm5 = AtmosphereState::from_altitude(5_000.0);
        let atm11 = AtmosphereState::from_altitude(11_000.0);
        assert!(
            atm5.density < atm0.density,
            "Density should decrease with altitude"
        );
        assert!(atm11.density < atm5.density);
    }

    // -------------------------------------------------------------------------
    // 4. ISA atmosphere: stratosphere isothermal
    // -------------------------------------------------------------------------
    #[test]
    fn test_isa_stratosphere_isothermal() {
        let atm12 = AtmosphereState::from_altitude(12_000.0);
        let atm15 = AtmosphereState::from_altitude(15_000.0);
        assert!(
            (atm12.temperature - atm15.temperature).abs() < 0.01,
            "Stratosphere should be isothermal: {} vs {}",
            atm12.temperature,
            atm15.temperature
        );
    }

    // -------------------------------------------------------------------------
    // 5. Dynamic pressure formula
    // -------------------------------------------------------------------------
    #[test]
    fn test_dynamic_pressure() {
        let atm = AtmosphereState::from_altitude(0.0);
        let q = atm.dynamic_pressure(100.0);
        let expected = 0.5 * RHO_SL * 100.0 * 100.0;
        assert!((q - expected).abs() < 1.0, "q = {q}, expected {expected}");
    }

    // -------------------------------------------------------------------------
    // 6. Lift-drag polar: CL = CL0 at alpha = alpha_zl
    // -------------------------------------------------------------------------
    #[test]
    fn test_polar_cl_at_zero_lift_alpha() {
        let polar = LiftDragPolar::typical_transport();
        let cl = polar.cl_from_alpha(polar.alpha_zl);
        assert!((cl - polar.cl0).abs() < EPS, "CL at alpha_zl should be CL0");
    }

    // -------------------------------------------------------------------------
    // 7. Lift-drag polar: CD >= CD0
    // -------------------------------------------------------------------------
    #[test]
    fn test_polar_cd_min_is_cd0() {
        let polar = LiftDragPolar::typical_transport();
        let cd = polar.cd_from_cl(0.0);
        assert!((cd - polar.cd0).abs() < EPS, "CD(CL=0) should equal CD0");
    }

    // -------------------------------------------------------------------------
    // 8. Lift-drag polar: max L/D location
    // -------------------------------------------------------------------------
    #[test]
    fn test_polar_max_ld_alpha() {
        let polar = LiftDragPolar::typical_transport();
        let alpha_opt = polar.alpha_at_max_ld();
        let ld_opt = polar.ld_ratio(alpha_opt);
        // Check that slightly perturbed alphas give lower L/D
        let ld_p = polar.ld_ratio(alpha_opt + 0.01);
        let ld_m = polar.ld_ratio(alpha_opt - 0.01);
        assert!(ld_opt >= ld_p, "L/D should be maximum at optimal alpha");
        assert!(ld_opt >= ld_m, "L/D should be maximum at optimal alpha");
    }

    // -------------------------------------------------------------------------
    // 9. Induced drag grows with CL²
    // -------------------------------------------------------------------------
    #[test]
    fn test_induced_drag_quadratic() {
        let polar = LiftDragPolar::typical_transport();
        let cdi1 = polar.induced_drag(1.0);
        let cdi2 = polar.induced_drag(2.0);
        assert!((cdi2 / cdi1 - 4.0).abs() < EPS, "CDi should scale as CL²");
    }

    // -------------------------------------------------------------------------
    // 10. Stability derivatives: cm_alpha negative for stable aircraft
    // -------------------------------------------------------------------------
    #[test]
    fn test_long_derivs_stable() {
        let d = LongitudinalDerivatives::typical_transport();
        assert!(d.cm_alpha < 0.0, "Stable aircraft needs Cm_alpha < 0");
    }

    // -------------------------------------------------------------------------
    // 11. Lateral derivatives: cn_beta positive for stable aircraft
    // -------------------------------------------------------------------------
    #[test]
    fn test_lat_derivs_directional_stable() {
        let d = LateralDerivatives::typical_transport();
        assert!(
            d.cn_beta > 0.0,
            "Stable aircraft needs Cn_beta > 0 (weathercock)"
        );
    }

    // -------------------------------------------------------------------------
    // 12. Phugoid mode: expected frequency range
    // -------------------------------------------------------------------------
    #[test]
    fn test_phugoid_frequency() {
        let ls = LongitudinalStability::new(
            150.0,
            RHO_SL,
            120.0,
            4.0,
            70_000.0,
            5_000_000.0,
            LongitudinalDerivatives::typical_transport(),
            LiftDragPolar::typical_transport(),
        );
        let mode = ls.phugoid_mode(0.05);
        // ω_ph ≈ √2 * 9.81 / 150 ≈ 0.0924 rad/s
        let expected = (2.0_f64).sqrt() * G / 150.0;
        assert!(
            (mode.omega_n - expected).abs() < 1e-4,
            "Phugoid ω = {}, expected {}",
            mode.omega_n,
            expected
        );
    }

    // -------------------------------------------------------------------------
    // 13. Phugoid period > short-period period (phugoid is slow)
    // -------------------------------------------------------------------------
    #[test]
    fn test_phugoid_slower_than_sp() {
        let ls = LongitudinalStability::new(
            150.0,
            RHO_SL,
            120.0,
            4.0,
            70_000.0,
            5_000_000.0,
            LongitudinalDerivatives::typical_transport(),
            LiftDragPolar::typical_transport(),
        );
        let phugoid = ls.phugoid_mode(0.05);
        let sp = ls.short_period_mode(0.05);
        assert!(
            phugoid.omega_n < sp.omega_n,
            "Phugoid should have lower frequency than short-period"
        );
    }

    // -------------------------------------------------------------------------
    // 14. Short-period stability (ζ > 0)
    // -------------------------------------------------------------------------
    #[test]
    fn test_short_period_stable() {
        let ls = LongitudinalStability::new(
            150.0,
            RHO_SL,
            120.0,
            4.0,
            70_000.0,
            5_000_000.0,
            LongitudinalDerivatives::typical_transport(),
            LiftDragPolar::typical_transport(),
        );
        let sp = ls.short_period_mode(0.05);
        assert!(
            sp.is_stable(),
            "Short-period mode should be stable for typical transport"
        );
    }

    // -------------------------------------------------------------------------
    // 15. Dutch roll mode: omega_n > 0
    // -------------------------------------------------------------------------
    #[test]
    fn test_dutch_roll_positive_frequency() {
        let lat = LateralStability::new(
            150.0,
            RHO_SL,
            120.0,
            28.0,
            70_000.0,
            2_000_000.0,
            8_000_000.0,
            200_000.0,
            LateralDerivatives::typical_transport(),
        );
        let dr = lat.dutch_roll_mode();
        assert!(
            dr.omega_n > 0.0,
            "Dutch roll should have positive frequency"
        );
    }

    // -------------------------------------------------------------------------
    // 16. Roll mode time constant positive (damped)
    // -------------------------------------------------------------------------
    #[test]
    fn test_roll_mode_time_constant() {
        let lat = LateralStability::new(
            150.0,
            RHO_SL,
            120.0,
            28.0,
            70_000.0,
            2_000_000.0,
            8_000_000.0,
            200_000.0,
            LateralDerivatives::typical_transport(),
        );
        let t_roll = lat.roll_mode_time_constant();
        assert!(
            t_roll > 0.0,
            "Roll mode time constant should be positive, got {t_roll}"
        );
    }

    // -------------------------------------------------------------------------
    // 17. Control surface clamp
    // -------------------------------------------------------------------------
    #[test]
    fn test_control_surface_clamp() {
        let elev = ControlSurface::elevator(-1.5);
        let delta_big = 1.0; // well beyond ±25°
        let clamped = elev.clamp(delta_big);
        assert!(
            clamped <= elev.deflection_max + EPS,
            "Deflection should be clamped, got {clamped}"
        );
    }

    // -------------------------------------------------------------------------
    // 18. Elevator moment coefficient correct sign
    // -------------------------------------------------------------------------
    #[test]
    fn test_elevator_moment() {
        let elev = ControlSurface::elevator(-1.5);
        let cm = elev.moment_coefficient(0.1); // small positive deflection
        // cm_de = -1.5 → negative nose-down moment
        assert!(
            cm < 0.0,
            "Positive elevator should give negative Cm, got {cm}"
        );
    }

    // -------------------------------------------------------------------------
    // 19. Propeller thrust at zero velocity (static)
    // -------------------------------------------------------------------------
    #[test]
    fn test_propeller_static_thrust() {
        let prop = PropellerThrustModel {
            diameter: 2.0,
            ct0: 0.08,
            ct_slope: -0.05,
            rpm_max: 2400.0,
            eta_p: 0.80,
        };
        let t = prop.thrust(0.0, 2400.0, RHO_SL);
        assert!(t > 0.0, "Static thrust should be positive, got {t}");
    }

    // -------------------------------------------------------------------------
    // 20. Propeller thrust decreases at higher advance ratio
    // -------------------------------------------------------------------------
    #[test]
    fn test_propeller_thrust_vs_advance() {
        let prop = PropellerThrustModel {
            diameter: 2.0,
            ct0: 0.08,
            ct_slope: -0.05,
            rpm_max: 2400.0,
            eta_p: 0.80,
        };
        let t_low = prop.thrust(30.0, 2400.0, RHO_SL);
        let t_high = prop.thrust(80.0, 2400.0, RHO_SL);
        assert!(t_high < t_low, "Thrust should decrease at higher airspeed");
    }

    // -------------------------------------------------------------------------
    // 21. Turbofan thrust lapse with altitude
    // -------------------------------------------------------------------------
    #[test]
    fn test_turbofan_thrust_lapse() {
        let engine = TurbofanThrustModelImpl::new(300_000.0);
        let atm0 = AtmosphereState::from_altitude(0.0);
        let atm10k = AtmosphereState::from_altitude(10_000.0);
        let t0 = engine.net_thrust(&atm0, 0.0);
        let t10 = engine.net_thrust(&atm10k, 0.0);
        assert!(t10 < t0, "Turbofan thrust should decrease with altitude");
    }

    // -------------------------------------------------------------------------
    // 22. Turbofan SFC is positive
    // -------------------------------------------------------------------------
    #[test]
    fn test_turbofan_sfc_positive() {
        let engine = TurbofanThrustModelImpl::new(300_000.0);
        assert!(engine.sfc() > 0.0, "SFC should be positive");
    }

    // -------------------------------------------------------------------------
    // 23. Trim condition: CL = weight / (q S)
    // -------------------------------------------------------------------------
    #[test]
    fn test_trim_cl_level_flight() {
        let polar = LiftDragPolar::typical_transport();
        let derivs = LongitudinalDerivatives::typical_transport();
        let weight = 600_000.0; // N (≈ 61 t)
        let s_ref = 120.0;
        let v_tas = 200.0;
        let atm = AtmosphereState::from_altitude(0.0);
        let q = atm.dynamic_pressure(v_tas);
        let trim = TrimCondition::straight_level(v_tas, 0.0, weight, s_ref, &polar, &derivs);
        let cl_expected = (weight / (q * s_ref)).min(polar.cl_max);
        assert!(
            (trim.cl_trim - cl_expected).abs() < 1e-6,
            "Trim CL should equal W/(qS), got {}",
            trim.cl_trim
        );
    }

    // -------------------------------------------------------------------------
    // 24. Trim thrust = drag in level flight
    // -------------------------------------------------------------------------
    #[test]
    fn test_trim_thrust_equals_drag() {
        let polar = LiftDragPolar::typical_transport();
        let derivs = LongitudinalDerivatives::typical_transport();
        let weight = 600_000.0;
        let s_ref = 120.0;
        let v_tas = 200.0;
        let trim = TrimCondition::straight_level(v_tas, 0.0, weight, s_ref, &polar, &derivs);
        let atm = AtmosphereState::from_altitude(0.0);
        let q = atm.dynamic_pressure(v_tas);
        let drag = q * s_ref * trim.cd_trim;
        assert!(
            (trim.thrust_required - drag).abs() < 1.0,
            "Trim thrust should equal drag: T={}, D={}",
            trim.thrust_required,
            drag
        );
    }

    // -------------------------------------------------------------------------
    // 25. Climbing trim thrust > level trim thrust (excess power)
    // -------------------------------------------------------------------------
    #[test]
    fn test_climbing_trim_higher_thrust() {
        let polar = LiftDragPolar::typical_transport();
        let derivs = LongitudinalDerivatives::typical_transport();
        let weight = 600_000.0;
        let s_ref = 120.0;
        let v_tas = 200.0;
        let level = TrimCondition::straight_level(v_tas, 0.0, weight, s_ref, &polar, &derivs);
        let climb = TrimCondition::climbing(v_tas, 0.0, weight, s_ref, 0.05, &polar, &derivs);
        assert!(
            climb.thrust_required > level.thrust_required,
            "Climbing requires more thrust: climb={}, level={}",
            climb.thrust_required,
            level.thrust_required
        );
    }

    // -------------------------------------------------------------------------
    // 26. DynamicMode: stable mode has positive t_half
    // -------------------------------------------------------------------------
    #[test]
    fn test_dynamic_mode_t_half() {
        let mode = DynamicMode {
            omega_n: 1.0,
            zeta: 0.5,
            name: "test",
        };
        assert!(mode.is_stable());
        assert!(mode.t_half() > 0.0);
    }

    // -------------------------------------------------------------------------
    // 27. DynamicMode: period formula
    // -------------------------------------------------------------------------
    #[test]
    fn test_dynamic_mode_period() {
        let mode = DynamicMode {
            omega_n: 1.0,
            zeta: 0.0,
            name: "test",
        };
        let expected = 2.0 * PI;
        assert!(
            (mode.period() - expected).abs() < EPS,
            "Period = 2π for ζ=0, ω=1"
        );
    }

    // -------------------------------------------------------------------------
    // 28. Flight envelope: v_min at sea level for typical aircraft
    // -------------------------------------------------------------------------
    #[test]
    fn test_flight_envelope_v_min() {
        let polar = LiftDragPolar::typical_transport();
        let env = FlightEnvelope {
            s_ref: 120.0,
            mtow: 600_000.0,
            polar: polar.clone(),
            t_max_sls: 500_000.0,
            v_d: 300.0,
            load_factor_limits: [-1.0, 2.5],
        };
        let atm = AtmosphereState::from_altitude(0.0);
        let v_min = env.v_min(&atm, 600_000.0);
        assert!(
            v_min > 30.0 && v_min < 120.0,
            "V_min should be reasonable: {v_min}"
        );
    }

    // -------------------------------------------------------------------------
    // 29. Flight envelope: v_max > v_min
    // -------------------------------------------------------------------------
    #[test]
    fn test_flight_envelope_v_max_gt_vmin() {
        let polar = LiftDragPolar::typical_transport();
        let env = FlightEnvelope {
            s_ref: 120.0,
            mtow: 600_000.0,
            polar: polar.clone(),
            t_max_sls: 500_000.0,
            v_d: 300.0,
            load_factor_limits: [-1.0, 2.5],
        };
        let atm = AtmosphereState::from_altitude(0.0);
        let v_min = env.v_min(&atm, 600_000.0);
        let v_max = env.v_max(&atm, 480_000.0);
        assert!(
            v_max > v_min,
            "V_max ({v_max}) should exceed V_min ({v_min})"
        );
    }

    // -------------------------------------------------------------------------
    // 30. 6-DOF state: TAS consistent with velocity
    // -------------------------------------------------------------------------
    #[test]
    fn test_sixdof_tas() {
        let s = AircraftState6Dof::level_flight(100.0);
        assert!(
            (s.tas() - 100.0).abs() < EPS,
            "TAS should equal u0 in level flight"
        );
    }

    // -------------------------------------------------------------------------
    // 31. 6-DOF state: alpha = 0 in level flight (no w component)
    // -------------------------------------------------------------------------
    #[test]
    fn test_sixdof_alpha_zero_level() {
        let s = AircraftState6Dof::level_flight(100.0);
        assert!(
            s.alpha().abs() < EPS,
            "Alpha should be zero in axial flight"
        );
    }

    // -------------------------------------------------------------------------
    // 32. Body-to-NED DCM: at zero Euler angles is identity
    // -------------------------------------------------------------------------
    #[test]
    fn test_body_to_ned_identity() {
        let s = AircraftState6Dof::level_flight(100.0);
        let dcm = s.body_to_ned();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dcm[i][j] - expected).abs() < EPS,
                    "DCM[{i}][{j}] = {}, expected {}",
                    dcm[i][j],
                    expected
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // 33. 6-DOF RK4 step produces finite values
    // -------------------------------------------------------------------------
    #[test]
    fn test_rk4_step_finite() {
        let integrator = SixDofIntegrator {
            mass: 70_000.0,
            inertia: [2_000_000.0, 5_000_000.0, 8_000_000.0],
            i_xz: 200_000.0,
            s_ref: 120.0,
            b_ref: 28.0,
            c_mac: 4.0,
            long_derivs: LongitudinalDerivatives::typical_transport(),
            lat_derivs: LateralDerivatives::typical_transport(),
            polar: LiftDragPolar::typical_transport(),
        };
        let state = AircraftState6Dof::level_flight(150.0);
        let new_state = integrator.rk4_step(&state, 200_000.0, 0.0, 0.0, 0.0, RHO_SL, 0.01);
        for v in new_state
            .velocity_body
            .iter()
            .chain(new_state.angular_rates.iter())
        {
            assert!(
                v.is_finite(),
                "State values should be finite after RK4 step"
            );
        }
    }

    // -------------------------------------------------------------------------
    // 34. Mach number formula
    // -------------------------------------------------------------------------
    #[test]
    fn test_mach_number() {
        let atm = AtmosphereState::from_altitude(0.0);
        let m = atm.mach(A_SL);
        assert!(
            (m - 1.0).abs() < 0.001,
            "Mach number at speed of sound should be 1.0, got {m}"
        );
    }

    // -------------------------------------------------------------------------
    // 35. Service ceiling is positive altitude
    // -------------------------------------------------------------------------
    #[test]
    fn test_service_ceiling_positive() {
        let polar = LiftDragPolar::typical_transport();
        let env = FlightEnvelope {
            s_ref: 120.0,
            mtow: 600_000.0,
            polar,
            t_max_sls: 500_000.0,
            v_d: 300.0,
            load_factor_limits: [-1.0, 2.5],
        };
        let ceil = env.service_ceiling(500_000.0, 600_000.0);
        assert!(
            ceil > 5_000.0,
            "Service ceiling should be above 5000 m, got {ceil}"
        );
    }
}
