// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wind turbine aerodynamics and structural dynamics.
//!
//! Implements the Blade Element Momentum (BEM) method for wind turbine
//! aerodynamics, including:
//!
//! - [`BetzLimit`] — theoretical maximum power coefficient
//! - [`BladeElement`] — lift and drag forces on a blade cross-section
//! - [`BemSolver`] — full BEM iteration with tip and hub loss corrections
//! - [`PrandtlTipLoss`] — Prandtl tip-loss correction factor
//! - [`YawMisalignment`] — yaw error corrections to thrust and power
//! - [`VariablePitchControl`] — collective pitch controller (PI)
//! - [`GeneratorTorqueControl`] — MPPT and rated power torque controller
//! - [`TowerShadow`] — potential-flow tower shadow deficit model
//! - [`BladeLoads`] — flapwise and edgewise bending moments
//! - [`RotorPerformance`] — integrated rotor thrust, torque, and power
//!
//! # Conventions
//!
//! * SI units: m, kg, s, N, Pa, W, rad.
//! * Wind direction along the +X axis; rotor axis along +X; Z is up.
//! * Blade azimuth measured from the 12-o'clock position (straight up).
//! * All angles in radians.
//!
//! # References
//!
//! * Hansen (2008) – Aerodynamics of Wind Turbines, 2nd ed.
//! * Burton et al. (2011) – Wind Energy Handbook, 2nd ed.
//! * IEC 61400-1 (2019) – Wind turbines – Part 1: Design requirements.
//! * Manwell, McGowan & Rogers (2009) – Wind Energy Explained, 2nd ed.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Clamp a value to `[lo, hi]`.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

// ---------------------------------------------------------------------------
// 1. Betz Limit
// ---------------------------------------------------------------------------

/// Betz limit analysis for an ideal actuator-disc wind turbine.
///
/// The Betz limit (`C_P,max = 16/27 ≈ 0.593`) is the theoretical maximum
/// power that can be extracted from a free-stream wind by an ideal rotor.
#[derive(Debug, Clone)]
pub struct BetzLimit {
    /// Air density ρ (kg m⁻³).
    pub rho: f64,
    /// Rotor disc area A (m²).
    pub disc_area: f64,
    /// Free-stream wind speed U∞ (m s⁻¹).
    pub wind_speed: f64,
}

impl BetzLimit {
    /// Create a new Betz-limit analysis.
    ///
    /// * `rho` — air density (kg m⁻³), typically 1.225 at sea level.
    /// * `rotor_radius` — rotor radius (m).
    /// * `wind_speed` — free-stream wind speed (m s⁻¹).
    pub fn new(rho: f64, rotor_radius: f64, wind_speed: f64) -> Self {
        Self {
            rho,
            disc_area: PI * rotor_radius * rotor_radius,
            wind_speed,
        }
    }

    /// Theoretical maximum power coefficient: `C_P,max = 16/27`.
    pub fn max_power_coefficient() -> f64 {
        16.0 / 27.0
    }

    /// Available wind power: `P_wind = ½ ρ A U∞³`.
    pub fn wind_power(&self) -> f64 {
        0.5 * self.rho * self.disc_area * self.wind_speed.powi(3)
    }

    /// Maximum extractable power: `P_max = C_P,max * P_wind`.
    pub fn max_power(&self) -> f64 {
        Self::max_power_coefficient() * self.wind_power()
    }

    /// Optimal axial induction factor: `a* = 1/3`.
    pub fn optimal_induction() -> f64 {
        1.0 / 3.0
    }

    /// Thrust coefficient at Betz optimum: `C_T = 8/9`.
    pub fn optimal_thrust_coefficient() -> f64 {
        8.0 / 9.0
    }

    /// Thrust force at Betz optimum: `T = C_T * ½ ρ A U∞²`.
    pub fn optimal_thrust(&self) -> f64 {
        Self::optimal_thrust_coefficient()
            * 0.5
            * self.rho
            * self.disc_area
            * self.wind_speed.powi(2)
    }

    /// Power for arbitrary axial induction `a` (0 ≤ a < 0.5):
    /// `C_P(a) = 4 a (1-a)²`.
    pub fn power_coefficient(a: f64) -> f64 {
        4.0 * a * (1.0 - a).powi(2)
    }

    /// Thrust coefficient for induction `a`: `C_T(a) = 4 a (1-a)`.
    pub fn thrust_coefficient(a: f64) -> f64 {
        4.0 * a * (1.0 - a)
    }
}

// ---------------------------------------------------------------------------
// 2. Lift and Drag Coefficient Profiles
// ---------------------------------------------------------------------------

/// A simple piecewise-linear aerofoil polar (CL, CD vs α).
///
/// Uses linear interpolation between tabulated (α, CL, CD) entries.
#[derive(Debug, Clone)]
pub struct AerofoilPolar {
    /// Tabulated angles of attack (radians).
    pub alpha_table: Vec<f64>,
    /// Lift coefficients corresponding to each angle.
    pub cl_table: Vec<f64>,
    /// Drag coefficients corresponding to each angle.
    pub cd_table: Vec<f64>,
}

impl AerofoilPolar {
    /// Create a polar from parallel tables.
    ///
    /// * `alpha_table` — angles of attack (rad), must be monotone increasing.
    /// * `cl_table` — lift coefficients.
    /// * `cd_table` — drag coefficients.
    pub fn new(alpha_table: Vec<f64>, cl_table: Vec<f64>, cd_table: Vec<f64>) -> Self {
        assert_eq!(alpha_table.len(), cl_table.len());
        assert_eq!(alpha_table.len(), cd_table.len());
        Self {
            alpha_table,
            cl_table,
            cd_table,
        }
    }

    /// Build a simple flat-plate polar valid for small angles.
    ///
    /// `CL = 2π α`, `CD = CD_min + k * α²`.
    pub fn flat_plate(n_points: usize, alpha_max_deg: f64, cd_min: f64) -> Self {
        let n = n_points.max(2);
        let alpha_max = alpha_max_deg.to_radians();
        let alphas: Vec<f64> = (0..n)
            .map(|i| -alpha_max + 2.0 * alpha_max * i as f64 / (n - 1) as f64)
            .collect();
        let cls: Vec<f64> = alphas.iter().map(|&a| 2.0 * PI * a).collect();
        let cds: Vec<f64> = alphas.iter().map(|&a| cd_min + 0.1 * a * a).collect();
        Self::new(alphas, cls, cds)
    }

    /// Linearly interpolate `(CL, CD)` for angle of attack `alpha` (rad).
    pub fn interpolate(&self, alpha: f64) -> (f64, f64) {
        let n = self.alpha_table.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        if alpha <= self.alpha_table[0] {
            return (self.cl_table[0], self.cd_table[0]);
        }
        if alpha >= self.alpha_table[n - 1] {
            return (self.cl_table[n - 1], self.cd_table[n - 1]);
        }
        for i in 0..n - 1 {
            let a0 = self.alpha_table[i];
            let a1 = self.alpha_table[i + 1];
            if alpha >= a0 && alpha <= a1 {
                let t = (alpha - a0) / (a1 - a0);
                let cl = self.cl_table[i] * (1.0 - t) + self.cl_table[i + 1] * t;
                let cd = self.cd_table[i] * (1.0 - t) + self.cd_table[i + 1] * t;
                return (cl, cd);
            }
        }
        (
            *self
                .cl_table
                .last()
                .expect("collection should not be empty"),
            *self
                .cd_table
                .last()
                .expect("collection should not be empty"),
        )
    }
}

// ---------------------------------------------------------------------------
// 3. Blade Element
// ---------------------------------------------------------------------------

/// A single blade element (annular ring) used in the BEM method.
#[derive(Debug, Clone)]
pub struct BladeElement {
    /// Radial position of the element centre (m).
    pub r: f64,
    /// Radial width of the element (m).
    pub dr: f64,
    /// Local chord length (m).
    pub chord: f64,
    /// Local twist angle (rad) — positive nose-up.
    pub twist: f64,
    /// Aerofoil polar for this element.
    pub polar: AerofoilPolar,
}

impl BladeElement {
    /// Construct a blade element.
    ///
    /// * `r` — radial position (m).
    /// * `dr` — element width (m).
    /// * `chord` — chord (m).
    /// * `twist` — local twist (rad).
    /// * `polar` — aerofoil polar.
    pub fn new(r: f64, dr: f64, chord: f64, twist: f64, polar: AerofoilPolar) -> Self {
        Self {
            r,
            dr,
            chord,
            twist,
            polar,
        }
    }

    /// Relative inflow velocity components.
    ///
    /// Returns `(V_rel, phi)` where `phi` is the inflow angle (rad).
    ///
    /// * `u_axial` — axial flow velocity at the rotor plane (m s⁻¹).
    /// * `omega_r` — tangential velocity at this element: `Ω r (1 + a')`.
    pub fn inflow(&self, u_axial: f64, omega_r: f64) -> (f64, f64) {
        let v_rel = (u_axial.powi(2) + omega_r.powi(2)).sqrt();
        let phi = u_axial.atan2(omega_r);
        (v_rel, phi)
    }

    /// Local angle of attack: `α = φ - (twist + pitch)` (rad).
    pub fn angle_of_attack(&self, phi: f64, pitch: f64) -> f64 {
        phi - (self.twist + pitch)
    }

    /// Normal (thrust) force per unit span (N m⁻¹).
    ///
    /// `dFn/dr = ½ ρ V_rel² c (CL cos φ + CD sin φ)`.
    pub fn normal_force_density(&self, rho: f64, v_rel: f64, phi: f64, cl: f64, cd: f64) -> f64 {
        0.5 * rho * v_rel.powi(2) * self.chord * (cl * phi.cos() + cd * phi.sin())
    }

    /// Tangential (torque) force per unit span (N m⁻¹).
    ///
    /// `dFt/dr = ½ ρ V_rel² c (CL sin φ - CD cos φ)`.
    pub fn tangential_force_density(
        &self,
        rho: f64,
        v_rel: f64,
        phi: f64,
        cl: f64,
        cd: f64,
    ) -> f64 {
        0.5 * rho * v_rel.powi(2) * self.chord * (cl * phi.sin() - cd * phi.cos())
    }

    /// Contribution of this element to rotor thrust (N).
    ///
    /// `dT = B * dFn * dr`.
    pub fn thrust_contribution(
        &self,
        rho: f64,
        v_rel: f64,
        phi: f64,
        cl: f64,
        cd: f64,
        n_blades: u32,
    ) -> f64 {
        let dfn = self.normal_force_density(rho, v_rel, phi, cl, cd);
        dfn * self.dr * n_blades as f64
    }

    /// Contribution of this element to rotor torque (N·m).
    ///
    /// `dQ = B * dFt * r * dr`.
    pub fn torque_contribution(
        &self,
        rho: f64,
        v_rel: f64,
        phi: f64,
        cl: f64,
        cd: f64,
        n_blades: u32,
    ) -> f64 {
        let dft = self.tangential_force_density(rho, v_rel, phi, cl, cd);
        dft * self.r * self.dr * n_blades as f64
    }
}

// ---------------------------------------------------------------------------
// 4. Prandtl Tip-Loss Correction
// ---------------------------------------------------------------------------

/// Prandtl tip-loss and hub-loss correction factors for the BEM method.
///
/// The tip-loss factor `F` accounts for the finite number of blades and reduces
/// the effective induction near the blade tip.
#[derive(Debug, Clone)]
pub struct PrandtlTipLoss {
    /// Number of blades.
    pub n_blades: u32,
    /// Rotor radius (m).
    pub r_tip: f64,
    /// Hub radius (m).
    pub r_hub: f64,
}

impl PrandtlTipLoss {
    /// Create a tip-loss corrector.
    ///
    /// * `n_blades` — number of blades.
    /// * `r_tip` — rotor tip radius (m).
    /// * `r_hub` — hub radius (m).
    pub fn new(n_blades: u32, r_tip: f64, r_hub: f64) -> Self {
        Self {
            n_blades,
            r_tip,
            r_hub,
        }
    }

    /// Prandtl tip-loss factor at radial position `r` and inflow angle `phi`.
    ///
    /// `F_tip = (2/π) arccos[ exp(-B(R-r)/(2r sinφ)) ]`.
    pub fn tip_factor(&self, r: f64, phi: f64) -> f64 {
        let sin_phi = phi.sin().abs().max(1e-6);
        let exp_arg = -(self.n_blades as f64) * (self.r_tip - r) / (2.0 * r * sin_phi);
        let f = (2.0 / PI) * exp_arg.exp().min(1.0).acos();
        clamp(f, 0.0, 1.0)
    }

    /// Prandtl hub-loss factor at radial position `r` and inflow angle `phi`.
    ///
    /// `F_hub = (2/π) arccos[ exp(-B(r-r_hub)/(2 r_hub sinφ)) ]`.
    pub fn hub_factor(&self, r: f64, phi: f64) -> f64 {
        let sin_phi = phi.sin().abs().max(1e-6);
        let exp_arg = -(self.n_blades as f64) * (r - self.r_hub) / (2.0 * self.r_hub * sin_phi);
        let f = (2.0 / PI) * exp_arg.exp().min(1.0).acos();
        clamp(f, 0.0, 1.0)
    }

    /// Combined tip and hub loss factor: `F = F_tip * F_hub`.
    pub fn combined_factor(&self, r: f64, phi: f64) -> f64 {
        self.tip_factor(r, phi) * self.hub_factor(r, phi)
    }
}

// ---------------------------------------------------------------------------
// 5. BEM Solver
// ---------------------------------------------------------------------------

/// Result for a single blade element BEM computation.
#[derive(Debug, Clone)]
pub struct BemElementResult {
    /// Axial induction factor `a`.
    pub a: f64,
    /// Tangential induction factor `a'`.
    pub a_prime: f64,
    /// Inflow angle `φ` (rad).
    pub phi: f64,
    /// Angle of attack `α` (rad).
    pub alpha: f64,
    /// Lift coefficient.
    pub cl: f64,
    /// Drag coefficient.
    pub cd: f64,
    /// Prandtl tip-loss factor.
    pub f_tip: f64,
    /// Local thrust contribution `dT` (N).
    pub d_thrust: f64,
    /// Local torque contribution `dQ` (N·m).
    pub d_torque: f64,
}

/// Full-rotor BEM solver.
///
/// Iterates the axial and tangential induction factors for each blade element
/// until convergence, accounting for tip and hub losses.
#[derive(Debug, Clone)]
pub struct BemSolver {
    /// Number of blades.
    pub n_blades: u32,
    /// Air density ρ (kg m⁻³).
    pub rho: f64,
    /// Free-stream wind speed U∞ (m s⁻¹).
    pub wind_speed: f64,
    /// Rotor angular velocity Ω (rad s⁻¹).
    pub omega: f64,
    /// Collective pitch angle (rad).
    pub pitch: f64,
    /// Blade elements.
    pub elements: Vec<BladeElement>,
    /// Tip-loss corrector.
    pub tip_loss: PrandtlTipLoss,
    /// Maximum BEM iterations.
    pub max_iter: usize,
    /// Convergence tolerance on induction factors.
    pub tolerance: f64,
}

impl BemSolver {
    /// Construct a BEM solver.
    ///
    /// * `n_blades` — number of blades.
    /// * `rho` — air density (kg m⁻³).
    /// * `wind_speed` — free-stream wind speed (m s⁻¹).
    /// * `omega` — rotor angular velocity (rad s⁻¹).
    /// * `pitch` — collective pitch (rad).
    /// * `elements` — blade element discretisation.
    /// * `r_hub` — hub radius (m).
    pub fn new(
        n_blades: u32,
        rho: f64,
        wind_speed: f64,
        omega: f64,
        pitch: f64,
        elements: Vec<BladeElement>,
        r_hub: f64,
    ) -> Self {
        let r_tip = elements.last().map(|e| e.r + e.dr * 0.5).unwrap_or(1.0);
        let tip_loss = PrandtlTipLoss::new(n_blades, r_tip, r_hub);
        Self {
            n_blades,
            rho,
            wind_speed,
            omega,
            pitch,
            elements,
            tip_loss,
            max_iter: 100,
            tolerance: 1e-6,
        }
    }

    /// Solve BEM for a single blade element.
    ///
    /// Returns a [`BemElementResult`] after iterating to convergence.
    pub fn solve_element(&self, elem: &BladeElement) -> BemElementResult {
        let b = self.n_blades as f64;
        let sigma = b * elem.chord / (2.0 * PI * elem.r);

        let mut a = 0.1_f64;
        let mut a_prime = 0.01_f64;

        for _ in 0..self.max_iter {
            let u_axial = self.wind_speed * (1.0 - a);
            let omega_r = self.omega * elem.r * (1.0 + a_prime);
            let (_v_rel, phi) = elem.inflow(u_axial, omega_r);
            let alpha = elem.angle_of_attack(phi, self.pitch);
            let (cl, cd) = elem.polar.interpolate(alpha);
            let f = self.tip_loss.combined_factor(elem.r, phi);

            // Normal and tangential force coefficients
            let cn = cl * phi.cos() + cd * phi.sin();
            let ct = cl * phi.sin() - cd * phi.cos();

            // Update induction factors
            let sin2 = phi.sin().powi(2);
            let new_a = if f > 1e-6 && sin2 > 1e-10 {
                let denom = 4.0 * f * sin2 / (sigma * cn) + 1.0;
                clamp(1.0 / denom, 0.0, 0.95)
            } else {
                0.0
            };
            let sin_cos = (phi.sin() * phi.cos()).abs().max(1e-10);
            let new_a_prime = if f > 1e-6 && sin_cos > 1e-10 {
                let denom = 4.0 * f * sin_cos / (sigma * ct) - 1.0;
                if denom.abs() > 1e-10 {
                    clamp(1.0 / denom, 0.0, 1.0)
                } else {
                    a_prime
                }
            } else {
                0.0
            };

            if (new_a - a).abs() < self.tolerance && (new_a_prime - a_prime).abs() < self.tolerance
            {
                a = new_a;
                a_prime = new_a_prime;
                break;
            }
            a = new_a;
            a_prime = new_a_prime;
        }

        let u_axial = self.wind_speed * (1.0 - a);
        let omega_r = self.omega * elem.r * (1.0 + a_prime);
        let (v_rel, phi) = elem.inflow(u_axial, omega_r);
        let alpha = elem.angle_of_attack(phi, self.pitch);
        let (cl, cd) = elem.polar.interpolate(alpha);
        let f_tip = self.tip_loss.tip_factor(elem.r, phi);
        let d_thrust = elem.thrust_contribution(self.rho, v_rel, phi, cl, cd, self.n_blades);
        let d_torque = elem.torque_contribution(self.rho, v_rel, phi, cl, cd, self.n_blades);

        BemElementResult {
            a,
            a_prime,
            phi,
            alpha,
            cl,
            cd,
            f_tip,
            d_thrust,
            d_torque,
        }
    }

    /// Solve BEM for all elements and return results.
    pub fn solve_all(&self) -> Vec<BemElementResult> {
        self.elements
            .iter()
            .map(|e| self.solve_element(e))
            .collect()
    }

    /// Compute total rotor thrust (N).
    pub fn rotor_thrust(&self) -> f64 {
        self.solve_all().iter().map(|r| r.d_thrust).sum()
    }

    /// Compute total rotor torque (N·m).
    pub fn rotor_torque(&self) -> f64 {
        self.solve_all().iter().map(|r| r.d_torque).sum()
    }

    /// Aerodynamic power (W): `P = Q * Ω`.
    pub fn rotor_power(&self) -> f64 {
        self.rotor_torque() * self.omega
    }

    /// Power coefficient: `C_P = P / (½ ρ A U∞³)`.
    pub fn power_coefficient(&self) -> f64 {
        let r_tip = self.tip_loss.r_tip;
        let p_wind = 0.5 * self.rho * PI * r_tip * r_tip * self.wind_speed.powi(3);
        if p_wind > 1e-10 {
            self.rotor_power() / p_wind
        } else {
            0.0
        }
    }

    /// Thrust coefficient: `C_T = T / (½ ρ A U∞²)`.
    pub fn thrust_coefficient(&self) -> f64 {
        let r_tip = self.tip_loss.r_tip;
        let q_wind = 0.5 * self.rho * PI * r_tip * r_tip * self.wind_speed.powi(2);
        if q_wind > 1e-10 {
            self.rotor_thrust() / q_wind
        } else {
            0.0
        }
    }

    /// Tip-speed ratio: `λ = Ω R / U∞`.
    pub fn tip_speed_ratio(&self) -> f64 {
        if self.wind_speed > 1e-10 {
            self.omega * self.tip_loss.r_tip / self.wind_speed
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// 6. Yaw Misalignment
// ---------------------------------------------------------------------------

/// Yaw misalignment model for a wind turbine operating at non-zero yaw angle.
///
/// Accounts for reduced energy capture and modified thrust at yaw angles
/// using the cosine-power law and skewed-wake corrections.
#[derive(Debug, Clone)]
pub struct YawMisalignment {
    /// Yaw misalignment angle γ (rad).
    pub yaw_angle: f64,
    /// Yaw power exponent `p` (typically 1.8–2.0).
    pub power_exponent: f64,
    /// Thrust reduction exponent `q` (typically 2.0).
    pub thrust_exponent: f64,
}

impl YawMisalignment {
    /// Create a new yaw misalignment model.
    ///
    /// * `yaw_angle` — yaw error (rad).
    /// * `power_exponent` — cosine power law exponent for power reduction.
    /// * `thrust_exponent` — cosine exponent for thrust reduction.
    pub fn new(yaw_angle: f64, power_exponent: f64, thrust_exponent: f64) -> Self {
        Self {
            yaw_angle,
            power_exponent,
            thrust_exponent,
        }
    }

    /// Power reduction factor: `η_P = cos^p(γ)`.
    pub fn power_factor(&self) -> f64 {
        self.yaw_angle.cos().powf(self.power_exponent)
    }

    /// Thrust reduction factor: `η_T = cos^q(γ)`.
    pub fn thrust_factor(&self) -> f64 {
        self.yaw_angle.cos().powf(self.thrust_exponent)
    }

    /// Effective wind speed at the rotor: `U_eff = U∞ cos(γ)`.
    pub fn effective_wind_speed(&self, free_stream: f64) -> f64 {
        free_stream * self.yaw_angle.cos()
    }

    /// Wake skew angle χ from Jiménez et al. (2009):
    /// `χ ≈ (0.6 a + 1) γ` where a is axial induction.
    pub fn wake_skew_angle(&self, axial_induction: f64) -> f64 {
        (0.6 * axial_induction + 1.0) * self.yaw_angle
    }

    /// Lateral wake deflection at distance `x` downstream (m).
    ///
    /// Simple linear deflection model: `δ = χ * x`.
    pub fn wake_deflection(&self, x: f64, axial_induction: f64) -> f64 {
        self.wake_skew_angle(axial_induction) * x
    }
}

// ---------------------------------------------------------------------------
// 7. Variable Pitch Control
// ---------------------------------------------------------------------------

/// Collective blade pitch controller (PI, Region 3 above rated).
///
/// Adjusts pitch to limit rotor speed to the rated value when wind exceeds
/// rated speed. Uses a proportional-integral (PI) controller.
#[derive(Debug, Clone)]
pub struct VariablePitchControl {
    /// Rated rotor speed Ω_rated (rad s⁻¹).
    pub omega_rated: f64,
    /// Proportional gain K_P (rad per rad s⁻¹ error).
    pub kp: f64,
    /// Integral gain K_I (rad per rad error).
    pub ki: f64,
    /// Minimum pitch angle (rad) — typically 0°.
    pub pitch_min: f64,
    /// Maximum pitch angle (rad) — typically 90° (fully feathered).
    pub pitch_max: f64,
    /// Current pitch demand (rad).
    pub pitch: f64,
    /// Integral error accumulator.
    pub integral: f64,
    /// Maximum pitch rate (rad s⁻¹).
    pub pitch_rate_max: f64,
}

impl VariablePitchControl {
    /// Construct a pitch controller.
    ///
    /// * `omega_rated` — rated speed (rad s⁻¹).
    /// * `kp` — proportional gain.
    /// * `ki` — integral gain.
    /// * `pitch_min` — minimum pitch (rad).
    /// * `pitch_max` — maximum pitch (rad).
    /// * `pitch_rate_max` — maximum pitch rate (rad s⁻¹).
    pub fn new(
        omega_rated: f64,
        kp: f64,
        ki: f64,
        pitch_min: f64,
        pitch_max: f64,
        pitch_rate_max: f64,
    ) -> Self {
        Self {
            omega_rated,
            kp,
            ki,
            pitch_min,
            pitch_max,
            pitch: pitch_min,
            integral: 0.0,
            pitch_rate_max,
        }
    }

    /// Update the pitch demand given the current rotor speed and time step.
    ///
    /// * `omega` — current rotor speed (rad s⁻¹).
    /// * `dt` — time step (s).
    ///
    /// Returns the updated pitch angle (rad).
    pub fn update(&mut self, omega: f64, dt: f64) -> f64 {
        let error = omega - self.omega_rated;
        self.integral += error * dt;
        let demand = self.kp * error + self.ki * self.integral;
        let new_pitch = clamp(self.pitch + demand * dt, self.pitch_min, self.pitch_max);
        // Rate limiting
        let delta = clamp(
            new_pitch - self.pitch,
            -self.pitch_rate_max * dt,
            self.pitch_rate_max * dt,
        );
        self.pitch += delta;
        self.pitch = clamp(self.pitch, self.pitch_min, self.pitch_max);
        self.pitch
    }

    /// Reset the integral accumulator (e.g. after mode switch).
    pub fn reset_integral(&mut self) {
        self.integral = 0.0;
    }
}

// ---------------------------------------------------------------------------
// 8. Generator Torque Control
// ---------------------------------------------------------------------------

/// Region of operation for a variable-speed wind turbine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperationRegion {
    /// Region 1: below cut-in — zero torque.
    BelowCutIn,
    /// Region 2: MPPT — maximum power tracking.
    Mppt,
    /// Region 2.5: transition between MPPT and rated.
    Transition,
    /// Region 3: rated power — constant power.
    Rated,
}

/// Variable-speed generator torque controller.
///
/// Implements MPPT (Region 2) and rated power (Region 3) control.
#[derive(Debug, Clone)]
pub struct GeneratorTorqueControl {
    /// Optimal torque coefficient: `Q = K_opt Ω²`.
    pub k_opt: f64,
    /// Rated generator torque (N·m).
    pub torque_rated: f64,
    /// Rated rotor speed (rad s⁻¹).
    pub omega_rated: f64,
    /// Cut-in rotor speed (rad s⁻¹).
    pub omega_cut_in: f64,
    /// Current operation region.
    pub region: OperationRegion,
    /// Current generator torque demand (N·m).
    pub torque: f64,
}

impl GeneratorTorqueControl {
    /// Construct a generator torque controller.
    ///
    /// * `k_opt` — MPPT gain `K = ½ ρ π R⁵ C_P,max / λ_opt³`.
    /// * `torque_rated` — rated generator torque (N·m).
    /// * `omega_rated` — rated rotor speed (rad s⁻¹).
    /// * `omega_cut_in` — cut-in rotor speed (rad s⁻¹).
    pub fn new(k_opt: f64, torque_rated: f64, omega_rated: f64, omega_cut_in: f64) -> Self {
        Self {
            k_opt,
            torque_rated,
            omega_rated,
            omega_cut_in,
            region: OperationRegion::BelowCutIn,
            torque: 0.0,
        }
    }

    /// Compute and update the generator torque demand for the current speed.
    ///
    /// * `omega` — current rotor speed (rad s⁻¹).
    ///
    /// Returns the demanded torque (N·m).
    pub fn update(&mut self, omega: f64) -> f64 {
        if omega < self.omega_cut_in {
            self.region = OperationRegion::BelowCutIn;
            self.torque = 0.0;
        } else if omega < self.omega_rated {
            self.region = OperationRegion::Mppt;
            self.torque = self.k_opt * omega.powi(2);
        } else {
            self.region = OperationRegion::Rated;
            self.torque = self.torque_rated;
        }
        self.torque
    }

    /// MPPT gain from rotor parameters.
    ///
    /// `K = ½ ρ π R⁵ C_Pmax / λ_opt³`.
    pub fn mppt_gain(rho: f64, r: f64, cp_max: f64, lambda_opt: f64) -> f64 {
        0.5 * rho * PI * r.powi(5) * cp_max / lambda_opt.powi(3)
    }
}

// ---------------------------------------------------------------------------
// 9. Tower Shadow Effect
// ---------------------------------------------------------------------------

/// Potential-flow tower shadow model for a cylindrical tower.
///
/// The velocity deficit behind a cylinder is modelled using the potential-flow
/// solution: `u = U∞ [1 - R_t² (x²-y²)/(x²+y²)²]` along the downwind axis.
#[derive(Debug, Clone)]
pub struct TowerShadow {
    /// Tower cylinder radius (m).
    pub tower_radius: f64,
    /// Free-stream wind speed (m s⁻¹).
    pub wind_speed: f64,
}

impl TowerShadow {
    /// Construct a tower shadow model.
    ///
    /// * `tower_radius` — tower outer radius at rotor height (m).
    /// * `wind_speed` — free-stream wind speed (m s⁻¹).
    pub fn new(tower_radius: f64, wind_speed: f64) -> Self {
        Self {
            tower_radius,
            wind_speed,
        }
    }

    /// Local axial velocity at point `(x, y)` relative to tower centre.
    ///
    /// Uses the potential-flow cylinder solution for the streamwise component.
    pub fn velocity_at(&self, x: f64, y: f64) -> f64 {
        let r2 = x * x + y * y;
        if r2 < self.tower_radius * self.tower_radius {
            return 0.0; // inside the tower
        }
        let rt2 = self.tower_radius * self.tower_radius;
        let factor = 1.0 - rt2 * (x * x - y * y) / (r2 * r2);
        self.wind_speed * factor
    }

    /// Velocity deficit fraction at `(x, y)`: `Δu / U∞`.
    pub fn deficit_fraction(&self, x: f64, y: f64) -> f64 {
        1.0 - self.velocity_at(x, y) / self.wind_speed.max(1e-10)
    }

    /// Blade passage through tower shadow: azimuthal deficit at angle `psi`
    /// for a blade at radius `r_blade` upwind distance `d_upwind`.
    ///
    /// `psi = 0` corresponds to the blade pointing straight down (toward tower).
    pub fn azimuthal_deficit(&self, psi: f64, r_blade: f64, d_upwind: f64) -> f64 {
        let x = d_upwind;
        let y = r_blade * psi.sin();
        self.deficit_fraction(x, y)
    }
}

// ---------------------------------------------------------------------------
// 10. Structural Loads
// ---------------------------------------------------------------------------

/// Flapwise and edgewise bending moment distributions on a wind turbine blade.
///
/// Flapwise bending (out-of-plane) is driven by aerodynamic thrust.
/// Edgewise bending (in-plane) is driven by aerodynamic torque and gravity.
#[derive(Debug, Clone)]
pub struct BladeLoads {
    /// Radial positions of load stations (m).
    pub r_stations: Vec<f64>,
    /// Flapwise distributed force at each station (N m⁻¹).
    pub flap_force: Vec<f64>,
    /// Edgewise distributed force at each station (N m⁻¹).
    pub edge_force: Vec<f64>,
    /// Blade mass distribution (kg m⁻¹).
    pub mass_per_length: Vec<f64>,
}

impl BladeLoads {
    /// Construct a blade loads object.
    ///
    /// * `r_stations` — radial positions of load stations (m).
    /// * `flap_force` — aerodynamic flapwise load per span (N m⁻¹).
    /// * `edge_force` — aerodynamic edgewise load per span (N m⁻¹).
    /// * `mass_per_length` — blade mass distribution (kg m⁻¹).
    pub fn new(
        r_stations: Vec<f64>,
        flap_force: Vec<f64>,
        edge_force: Vec<f64>,
        mass_per_length: Vec<f64>,
    ) -> Self {
        Self {
            r_stations,
            flap_force,
            edge_force,
            mass_per_length,
        }
    }

    /// Flapwise bending moment at station index `i` (N·m).
    ///
    /// Integrates the distributed load from tip to station `i`.
    pub fn flapwise_moment(&self, i: usize) -> f64 {
        let n = self.r_stations.len();
        if i >= n {
            return 0.0;
        }
        let mut m = 0.0;
        for j in (i + 1)..n {
            let r_mid = 0.5 * (self.r_stations[j] + self.r_stations[j - 1]);
            let dr = self.r_stations[j] - self.r_stations[j - 1];
            let f_avg = 0.5 * (self.flap_force[j] + self.flap_force[j - 1]);
            let moment_arm = r_mid - self.r_stations[i];
            m += f_avg * dr * moment_arm;
        }
        m
    }

    /// Edgewise bending moment at station index `i` (N·m).
    ///
    /// Includes aerodynamic edgewise loads and gravity (for azimuth `psi`).
    pub fn edgewise_moment(&self, i: usize, psi: f64, gravity: f64) -> f64 {
        let n = self.r_stations.len();
        if i >= n {
            return 0.0;
        }
        let mut m = 0.0;
        for j in (i + 1)..n {
            let r_mid = 0.5 * (self.r_stations[j] + self.r_stations[j - 1]);
            let dr = self.r_stations[j] - self.r_stations[j - 1];
            let f_avg = 0.5 * (self.edge_force[j] + self.edge_force[j - 1]);
            let m_avg = 0.5 * (self.mass_per_length[j] + self.mass_per_length[j - 1]);
            let moment_arm = r_mid - self.r_stations[i];
            // Aerodynamic edgewise
            m += f_avg * dr * moment_arm;
            // Gravity (azimuth-dependent): M_grav = m * g * sin(ψ) * arm
            m += m_avg * dr * gravity * psi.sin() * moment_arm;
        }
        m
    }

    /// Root flapwise moment (N·m) — at blade root (station 0).
    pub fn root_flapwise_moment(&self) -> f64 {
        self.flapwise_moment(0)
    }

    /// Root edgewise moment at azimuth `psi` (N·m).
    pub fn root_edgewise_moment(&self, psi: f64, gravity: f64) -> f64 {
        self.edgewise_moment(0, psi, gravity)
    }

    /// Normalised fatigue loading indicator: RMS of flapwise moment.
    pub fn flapwise_rms(&self) -> f64 {
        let n = self.r_stations.len();
        if n == 0 {
            return 0.0;
        }
        let sum_sq: f64 = (0..n).map(|i| self.flapwise_moment(i).powi(2)).sum();
        (sum_sq / n as f64).sqrt()
    }
}

// ---------------------------------------------------------------------------
// 11. Rotor Performance Summary
// ---------------------------------------------------------------------------

/// Summary of rotor aerodynamic performance.
#[derive(Debug, Clone)]
pub struct RotorPerformance {
    /// Total thrust T (N).
    pub thrust: f64,
    /// Total torque Q (N·m).
    pub torque: f64,
    /// Aerodynamic power P (W).
    pub power: f64,
    /// Power coefficient C_P.
    pub cp: f64,
    /// Thrust coefficient C_T.
    pub ct: f64,
    /// Tip-speed ratio λ.
    pub tsr: f64,
}

impl RotorPerformance {
    /// Compute rotor performance from a configured [`BemSolver`].
    pub fn from_bem(bem: &BemSolver) -> Self {
        let thrust = bem.rotor_thrust();
        let torque = bem.rotor_torque();
        let power = torque * bem.omega;
        let cp = bem.power_coefficient();
        let ct = bem.thrust_coefficient();
        let tsr = bem.tip_speed_ratio();
        Self {
            thrust,
            torque,
            power,
            cp,
            ct,
            tsr,
        }
    }

    /// Specific power: power per unit swept area (W m⁻²).
    pub fn specific_power(&self, rotor_radius: f64) -> f64 {
        let area = PI * rotor_radius * rotor_radius;
        if area > 1e-10 { self.power / area } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// 12. Hub Wake Model
// ---------------------------------------------------------------------------

/// Simplified hub wake model with vortex core correction.
///
/// The hub vortex creates an induction region at the rotor centre that
/// modifies the axial induction for small radii.
#[derive(Debug, Clone)]
pub struct HubWakeModel {
    /// Hub radius (m).
    pub r_hub: f64,
    /// Rotor radius (m).
    pub r_tip: f64,
    /// Hub vortex strength Γ_hub (m² s⁻¹).
    pub hub_vortex_strength: f64,
    /// Rankine vortex core radius (m).
    pub core_radius: f64,
}

impl HubWakeModel {
    /// Construct a hub wake model.
    ///
    /// * `r_hub` — hub radius (m).
    /// * `r_tip` — rotor tip radius (m).
    /// * `hub_vortex_strength` — hub vortex circulation (m² s⁻¹).
    pub fn new(r_hub: f64, r_tip: f64, hub_vortex_strength: f64) -> Self {
        Self {
            r_hub,
            r_tip,
            hub_vortex_strength,
            core_radius: r_hub * 0.5,
        }
    }

    /// Induced tangential velocity from the hub vortex at radius `r` (m s⁻¹).
    ///
    /// Rankine vortex: solid rotation inside core, 1/r outside.
    pub fn induced_tangential_velocity(&self, r: f64) -> f64 {
        if r < self.core_radius {
            self.hub_vortex_strength * r / (2.0 * PI * self.core_radius.powi(2))
        } else {
            self.hub_vortex_strength / (2.0 * PI * r)
        }
    }

    /// Axial induction correction at radius `r` due to hub vortex.
    ///
    /// Uses the Biot-Savart-based correction: `Δa ≈ Γ_hub / (4π U∞ r)`.
    pub fn axial_induction_correction(&self, r: f64, wind_speed: f64) -> f64 {
        if r < 1e-10 || wind_speed < 1e-10 {
            return 0.0;
        }
        self.hub_vortex_strength / (4.0 * PI * wind_speed * r)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // --- Betz Limit ---

    #[test]
    fn test_betz_limit_constant() {
        let cp_max = BetzLimit::max_power_coefficient();
        assert!((cp_max - 16.0 / 27.0).abs() < EPS, "Betz cp_max={cp_max}");
    }

    #[test]
    fn test_betz_wind_power_scaling() {
        let b1 = BetzLimit::new(1.225, 50.0, 10.0);
        let b2 = BetzLimit::new(1.225, 50.0, 20.0);
        // P scales as U^3 => factor 8
        let ratio = b2.wind_power() / b1.wind_power();
        assert!((ratio - 8.0).abs() < 1e-6, "Wind power ratio={ratio}");
    }

    #[test]
    fn test_betz_max_power_less_than_wind_power() {
        let b = BetzLimit::new(1.225, 50.0, 12.0);
        assert!(
            b.max_power() < b.wind_power(),
            "Betz power must be less than total wind power"
        );
    }

    #[test]
    fn test_betz_power_coefficient_function() {
        let cp = BetzLimit::power_coefficient(1.0 / 3.0);
        assert!(
            (cp - 16.0 / 27.0).abs() < 1e-10,
            "C_P at a=1/3 should be Betz limit: {cp}"
        );
    }

    #[test]
    fn test_betz_thrust_coefficient_optimal() {
        let ct = BetzLimit::optimal_thrust_coefficient();
        assert!((ct - 8.0 / 9.0).abs() < EPS, "C_T,opt={ct}");
    }

    // --- AerofoilPolar ---

    #[test]
    fn test_polar_flat_plate_cl_at_zero_is_zero() {
        let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
        let (cl, _cd) = polar.interpolate(0.0);
        assert!(
            cl.abs() < 1e-6,
            "CL at alpha=0 should be ~0 for flat plate: {cl}"
        );
    }

    #[test]
    fn test_polar_cl_positive_at_positive_alpha() {
        let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
        let alpha = 5_f64.to_radians();
        let (cl, _cd) = polar.interpolate(alpha);
        assert!(cl > 0.0, "CL should be positive at positive alpha: {cl}");
    }

    #[test]
    fn test_polar_cd_always_positive() {
        let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
        for &alpha_deg in &[-15.0_f64, -5.0, 0.0, 5.0, 15.0] {
            let (_cl, cd) = polar.interpolate(alpha_deg.to_radians());
            assert!(cd > 0.0, "CD must be positive at alpha={alpha_deg}°: {cd}");
        }
    }

    // --- PrandtlTipLoss ---

    #[test]
    fn test_tip_loss_at_tip_approaches_zero() {
        let tl = PrandtlTipLoss::new(3, 50.0, 2.0);
        let phi = 0.1_f64;
        let f = tl.tip_factor(50.0, phi);
        // At r = R, exponent → 0, arccos(1) = 0, so F → 0
        assert!(f < 0.1, "Tip loss at r=R should be small: {f}");
    }

    #[test]
    fn test_tip_loss_increases_towards_hub() {
        let tl = PrandtlTipLoss::new(3, 50.0, 2.0);
        let phi = 0.2_f64;
        let f_near_tip = tl.tip_factor(48.0, phi);
        let f_mid = tl.tip_factor(25.0, phi);
        assert!(
            f_mid > f_near_tip,
            "Tip loss should be larger away from tip: f_mid={f_mid}, f_near_tip={f_near_tip}"
        );
    }

    #[test]
    fn test_combined_loss_factor_bounded() {
        let tl = PrandtlTipLoss::new(3, 50.0, 2.0);
        for r in [5.0, 10.0, 25.0, 40.0, 49.0] {
            let f = tl.combined_factor(r, 0.15);
            assert!(
                (0.0..=1.0).contains(&f),
                "Combined loss factor must be in [0,1]: f={f} at r={r}"
            );
        }
    }

    // --- BladeElement ---

    #[test]
    fn test_blade_element_inflow_angle() {
        let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
        let elem = BladeElement::new(25.0, 1.0, 2.5, 5_f64.to_radians(), polar);
        let (v_rel, phi) = elem.inflow(10.0, 50.0);
        assert!(v_rel > 0.0, "V_rel must be positive: {v_rel}");
        assert!(phi > 0.0, "phi must be positive: {phi}");
    }

    #[test]
    fn test_blade_element_normal_force_positive() {
        let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
        let elem = BladeElement::new(25.0, 1.0, 2.5, 0.0, polar);
        let (v_rel, phi) = elem.inflow(8.0, 50.0);
        let alpha = elem.angle_of_attack(phi, 0.0);
        let (cl, cd) = elem.polar.interpolate(alpha);
        let fn_ = elem.normal_force_density(1.225, v_rel, phi, cl, cd);
        assert!(fn_ != 0.0, "Normal force density should be non-zero");
    }

    // --- BemSolver ---

    fn build_simple_rotor() -> BemSolver {
        let n = 10;
        let r_tip = 50.0_f64;
        let r_hub = 2.0_f64;
        let elements: Vec<BladeElement> = (0..n)
            .map(|i| {
                let r = r_hub + (r_tip - r_hub) * (i as f64 + 0.5) / n as f64;
                let dr = (r_tip - r_hub) / n as f64;
                let chord = 3.0 - 2.0 * r / r_tip; // linearly tapered
                let twist = 10_f64.to_radians() * (1.0 - r / r_tip);
                let polar = AerofoilPolar::flat_plate(21, 20.0, 0.01);
                BladeElement::new(r, dr, chord, twist, polar)
            })
            .collect();
        BemSolver::new(3, 1.225, 10.0, 1.26, 0.0, elements, r_hub)
    }

    #[test]
    fn test_bem_thrust_positive() {
        let bem = build_simple_rotor();
        let t = bem.rotor_thrust();
        assert!(t > 0.0, "Rotor thrust should be positive: {t}");
    }

    #[test]
    fn test_bem_power_positive() {
        let bem = build_simple_rotor();
        let p = bem.rotor_power();
        assert!(p > 0.0, "Rotor power should be positive: {p}");
    }

    #[test]
    fn test_bem_cp_bounded() {
        let bem = build_simple_rotor();
        let cp = bem.power_coefficient();
        assert!(cp > 0.0 && cp < 1.0, "C_P must be in (0,1): {cp}");
    }

    #[test]
    fn test_bem_ct_bounded() {
        let bem = build_simple_rotor();
        let ct = bem.thrust_coefficient();
        assert!(
            ct > 0.0 && ct < 2.0,
            "C_T should be physically reasonable: {ct}"
        );
    }

    #[test]
    fn test_bem_tsr_correct() {
        let bem = build_simple_rotor();
        let tsr = bem.tip_speed_ratio();
        // λ = 1.26 * 50 / 10 = 6.3
        assert!((tsr - 6.3).abs() < 0.01, "TSR={tsr}");
    }

    #[test]
    fn test_bem_induction_converged() {
        let bem = build_simple_rotor();
        let results = bem.solve_all();
        for r in &results {
            assert!(r.a >= 0.0 && r.a <= 0.95, "a={} out of range", r.a);
            assert!(r.a_prime >= 0.0, "a'={} must be non-negative", r.a_prime);
        }
    }

    // --- YawMisalignment ---

    #[test]
    fn test_yaw_zero_no_reduction() {
        let yaw = YawMisalignment::new(0.0, 1.88, 2.0);
        assert!(
            (yaw.power_factor() - 1.0).abs() < EPS,
            "Zero yaw should give unit power factor"
        );
        assert!(
            (yaw.thrust_factor() - 1.0).abs() < EPS,
            "Zero yaw should give unit thrust factor"
        );
    }

    #[test]
    fn test_yaw_30deg_power_reduced() {
        let yaw = YawMisalignment::new(30_f64.to_radians(), 1.88, 2.0);
        let pf = yaw.power_factor();
        assert!(
            pf < 1.0 && pf > 0.0,
            "Power factor at 30° yaw should be in (0,1): {pf}"
        );
    }

    #[test]
    fn test_yaw_effective_wind_speed() {
        let yaw = YawMisalignment::new(60_f64.to_radians(), 1.88, 2.0);
        let u_eff = yaw.effective_wind_speed(10.0);
        assert!(
            (u_eff - 5.0).abs() < 1e-6,
            "u_eff at 60° yaw should be 5 m/s: {u_eff}"
        );
    }

    // --- VariablePitchControl ---

    #[test]
    fn test_pitch_controller_increases_pitch_on_overspeed() {
        let mut ctrl = VariablePitchControl::new(1.26, 0.5, 0.1, 0.0, PI / 2.0, 0.1);
        let omega_high = 1.5_f64; // above rated
        let mut pitch = ctrl.pitch;
        for _ in 0..50 {
            pitch = ctrl.update(omega_high, 0.1);
        }
        assert!(pitch > 0.0, "Pitch should increase on overspeed: {pitch}");
    }

    #[test]
    fn test_pitch_controller_bounds_respected() {
        let mut ctrl = VariablePitchControl::new(1.26, 0.5, 0.1, 0.0, PI / 2.0, 10.0);
        for _ in 0..1000 {
            ctrl.update(5.0, 0.01); // very high speed
        }
        assert!(
            ctrl.pitch <= PI / 2.0,
            "Pitch must not exceed max: {}",
            ctrl.pitch
        );
        assert!(
            ctrl.pitch >= 0.0,
            "Pitch must not go below min: {}",
            ctrl.pitch
        );
    }

    // --- GeneratorTorqueControl ---

    #[test]
    fn test_generator_torque_below_cut_in_is_zero() {
        let mut ctrl = GeneratorTorqueControl::new(1e5, 5e5, 1.26, 0.5);
        let q = ctrl.update(0.1);
        assert!(q == 0.0, "Below cut-in torque should be zero: {q}");
        assert_eq!(ctrl.region, OperationRegion::BelowCutIn);
    }

    #[test]
    fn test_generator_torque_mppt_region() {
        let mut ctrl = GeneratorTorqueControl::new(1e5, 5e5, 1.26, 0.5);
        let q = ctrl.update(1.0);
        assert!(
            (q - 1e5).abs() < EPS,
            "MPPT torque = K_opt * 1^2 = K_opt: {q}"
        );
        assert_eq!(ctrl.region, OperationRegion::Mppt);
    }

    #[test]
    fn test_generator_torque_rated_region() {
        let mut ctrl = GeneratorTorqueControl::new(1e5, 5e5, 1.26, 0.5);
        ctrl.update(2.0); // above rated
        assert_eq!(ctrl.region, OperationRegion::Rated);
        assert!(
            (ctrl.torque - 5e5).abs() < EPS,
            "Rated torque mismatch: {}",
            ctrl.torque
        );
    }

    // --- TowerShadow ---

    #[test]
    fn test_tower_shadow_deficit_downwind() {
        let ts = TowerShadow::new(2.0, 10.0);
        // Directly downwind on rotor axis (y=0, x > R)
        let u = ts.velocity_at(5.0, 0.0);
        // At x=5, y=0: factor = 1 - R^2*(x^2-0)/(x^2)^2 = 1 - R^2/x^2
        let expected = 10.0 * (1.0 - 4.0 / 25.0);
        assert!(
            (u - expected).abs() < 1e-10,
            "Tower shadow velocity: {u} vs {expected}"
        );
    }

    #[test]
    fn test_tower_shadow_zero_inside_tower() {
        let ts = TowerShadow::new(2.0, 10.0);
        let u = ts.velocity_at(1.0, 0.0); // r < tower_radius
        assert!(u == 0.0, "Velocity inside tower must be zero: {u}");
    }

    #[test]
    fn test_tower_shadow_full_speed_far_away() {
        let ts = TowerShadow::new(2.0, 10.0);
        let u = ts.velocity_at(1000.0, 0.0);
        // At x=1000, deficit ≈ R^2/x^2 ≈ 4e-6 → u ≈ U_inf
        assert!(
            (u - 10.0).abs() < 0.01,
            "Far from tower speed should be ~U_inf: {u}"
        );
    }

    // --- BladeLoads ---

    #[test]
    fn test_blade_loads_root_moment_positive() {
        let rs = vec![0.0, 10.0, 20.0, 30.0, 40.0, 50.0];
        let ff = vec![500.0, 450.0, 380.0, 300.0, 200.0, 50.0];
        let ef = vec![100.0, 90.0, 80.0, 60.0, 40.0, 10.0];
        let ml = vec![150.0, 120.0, 90.0, 60.0, 30.0, 10.0];
        let loads = BladeLoads::new(rs, ff, ef, ml);
        let m = loads.root_flapwise_moment();
        assert!(m > 0.0, "Root flapwise moment should be positive: {m}");
    }

    #[test]
    fn test_blade_loads_moment_decreases_toward_tip() {
        let rs = vec![0.0, 10.0, 25.0, 40.0, 50.0];
        let ff = vec![400.0, 350.0, 250.0, 150.0, 0.0];
        let ef = vec![80.0, 70.0, 50.0, 30.0, 0.0];
        let ml = vec![100.0, 80.0, 60.0, 30.0, 5.0];
        let loads = BladeLoads::new(rs, ff, ef, ml);
        let m0 = loads.flapwise_moment(0);
        let m2 = loads.flapwise_moment(2);
        let m4 = loads.flapwise_moment(4);
        assert!(
            m0 > m2,
            "Moment at root should exceed mid-span: m0={m0}, m2={m2}"
        );
        assert!(m4 == 0.0, "Moment at tip should be zero: {m4}");
    }

    #[test]
    fn test_blade_loads_edgewise_gravity_dependence() {
        let rs = vec![0.0, 25.0, 50.0];
        let ff = vec![300.0, 200.0, 0.0];
        let ef = vec![50.0, 30.0, 0.0];
        let ml = vec![100.0, 50.0, 5.0];
        let loads = BladeLoads::new(rs, ff, ef, ml);
        let m_0 = loads.root_edgewise_moment(0.0, 9.81);
        let m_pi2 = loads.root_edgewise_moment(PI / 2.0, 9.81);
        // At psi=PI/2 gravity contributes positively; at psi=0 gravity term is zero
        assert!(
            m_pi2 != m_0,
            "Edgewise moment should differ at different azimuths"
        );
    }

    // --- RotorPerformance ---

    #[test]
    fn test_rotor_performance_from_bem() {
        let bem = build_simple_rotor();
        let perf = RotorPerformance::from_bem(&bem);
        assert!(perf.power > 0.0, "Power must be positive: {}", perf.power);
        assert!(
            perf.cp > 0.0 && perf.cp < 1.0,
            "C_P must be in (0,1): {}",
            perf.cp
        );
    }

    #[test]
    fn test_rotor_specific_power_positive() {
        let bem = build_simple_rotor();
        let perf = RotorPerformance::from_bem(&bem);
        let sp = perf.specific_power(50.0);
        assert!(sp > 0.0, "Specific power must be positive: {sp}");
    }

    // --- HubWakeModel ---

    #[test]
    fn test_hub_wake_induced_velocity_positive() {
        let hw = HubWakeModel::new(2.0, 50.0, 10.0);
        let v = hw.induced_tangential_velocity(5.0);
        assert!(v > 0.0, "Induced velocity should be positive: {v}");
    }

    #[test]
    fn test_hub_wake_rankine_core_solid_rotation() {
        let hw = HubWakeModel::new(2.0, 50.0, 100.0);
        let v1 = hw.induced_tangential_velocity(0.3); // inside core
        let v2 = hw.induced_tangential_velocity(0.6); // inside core
        // Inside core: v ∝ r
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-6,
            "Solid rotation: v2/v1 should be 2: {}",
            v2 / v1
        );
    }

    #[test]
    fn test_hub_wake_axial_correction_decreases_with_r() {
        let hw = HubWakeModel::new(2.0, 50.0, 10.0);
        let c1 = hw.axial_induction_correction(5.0, 10.0);
        let c2 = hw.axial_induction_correction(25.0, 10.0);
        assert!(
            c1 > c2,
            "Axial correction should decrease with radius: c1={c1}, c2={c2}"
        );
    }
}
