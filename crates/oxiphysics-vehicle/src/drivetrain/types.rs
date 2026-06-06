//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
pub use super::types_ext::{Gearbox, GearboxLegacy};
use oxiphysics_core::math::Real;

/// Complete drivetrain: engine → clutch → gearbox → differential → wheels.
#[derive(Debug, Clone)]
pub struct Drivetrain {
    /// Engine model.
    pub engine: Engine,
    /// Clutch between engine and gearbox.
    pub clutch: Clutch,
    /// Gearbox (gear selection and ratio).
    pub gearbox: Gearbox,
    /// Differential (torque split to left/right wheels).
    pub differential: Differential,
    /// Driven wheel radius (m), used for speed ↔ torque conversion.
    pub wheel_radius: f64,
}
impl Drivetrain {
    /// Create a new drivetrain.
    pub fn new(engine: Engine, gearbox: Gearbox, wheel_radius: f64) -> Self {
        let max_t = engine.curve.torque_at_rpm(engine.idle_rpm) * 3.0;
        Self {
            engine,
            clutch: Clutch::new(max_t.max(500.0)),
            gearbox,
            differential: Differential::open(),
            wheel_radius,
        }
    }
    /// Torque delivered to the left and right drive wheels (N·m each).
    pub fn wheel_torque(&self) -> (f64, f64) {
        let engine_t = self.engine.available_torque();
        let clutch_t = self.clutch.transmitted_torque(engine_t);
        let axle_t = self.gearbox.output_torque(clutch_t);
        self.differential.split_torque(axle_t, 0.0)
    }
    /// Convert wheel RPM to engine RPM via current total gear ratio.
    pub fn wheel_rpm_to_engine_rpm(&self, wheel_rpm: f64) -> f64 {
        wheel_rpm * self.gearbox.total_ratio().abs()
    }
    /// Advance the full drivetrain by `dt` seconds.
    ///
    /// * `throttle`     – throttle position \[0, 1\]
    /// * `brake_torque` – opposing torque from brakes (N·m, ≥ 0)
    /// * `wheel_speed`  – current wheel angular velocity (rad/s)
    /// * `dt`           – timestep (seconds)
    pub fn step(&mut self, throttle: f64, brake_torque: f64, wheel_speed: f64, dt: f64) {
        self.engine.throttle = throttle.clamp(0.0, 1.0);
        let ratio = self.gearbox.total_ratio().abs();
        let engine_omega_from_wheels = wheel_speed * ratio;
        let engine_omega = self.engine.rpm * 2.0 * std::f64::consts::PI / 60.0;
        self.clutch.slip_speed = engine_omega - engine_omega_from_wheels;
        let (wl, wr) = self.wheel_torque();
        let wheel_load = (wl + wr) + brake_torque;
        let engine_load = if ratio > 1e-10 {
            wheel_load / ratio
        } else {
            0.0
        };
        self.engine.update(engine_load, dt);
    }
    /// Net vehicle acceleration (m/s²) at the current drivetrain state.
    ///
    /// `vehicle_mass` in kg, `drag` in N (opposing force at current speed).
    pub fn vehicle_acceleration(&self, vehicle_mass: f64, drag: f64) -> f64 {
        let (wl, wr) = self.wheel_torque();
        let drive_force = (wl + wr) / self.wheel_radius.max(1e-6);
        (drive_force - drag) / vehicle_mass.max(1.0)
    }
}
impl Drivetrain {
    /// Compute the total effective powertrain rotating inertia referred to the
    /// engine crankshaft (kg·m²).
    ///
    /// Contributions:
    /// 1. Engine inertia: `I_engine`
    /// 2. Gearbox input shaft: `k_gearbox * I_engine` (typical fraction)
    /// 3. Drive shaft and differential: `k_driveline * I_engine`
    /// 4. Wheel/tire inertia referred back through gear ratio:
    ///    `I_wheel_eff = I_wheel / (gear_ratio^2)` per driven wheel.
    ///
    /// # Arguments
    /// * `wheel_inertia` – rotational inertia of each driven wheel+tire (kg·m²)
    ///
    /// Returns total inertia in kg·m² referred to engine shaft.
    pub fn compute_powertrain_inertia(&self, wheel_inertia: f64) -> f64 {
        let i_engine = self.engine.inertia;
        let k_gearbox = 0.15_f64;
        let k_driveline = 0.10_f64;
        let i_gearbox = k_gearbox * i_engine;
        let i_driveline = k_driveline * i_engine;
        let ratio = self.gearbox.total_ratio().abs();
        let i_wheels_referred = if ratio > 1e-9 {
            2.0 * wheel_inertia / (ratio * ratio)
        } else {
            0.0
        };
        i_engine + i_gearbox + i_driveline + i_wheels_referred
    }
}
/// Driveshaft torsional compliance model.
///
/// Models the elastic wind-up of the driveshaft between transmission output
/// and driven wheels. Uses a spring-damper (Voigt) model.
#[derive(Debug, Clone)]
pub struct DriveshaftCompliance {
    /// Torsional stiffness (N·m/rad).
    pub stiffness: f64,
    /// Torsional damping coefficient (N·m·s/rad).
    pub damping: f64,
    /// Current torsional deflection (rad).
    pub deflection: f64,
    /// Current angular velocity of deflection (rad/s).
    pub deflection_rate: f64,
}
impl DriveshaftCompliance {
    /// Create a driveshaft compliance model.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self {
            stiffness: stiffness.max(1.0),
            damping: damping.max(0.0),
            deflection: 0.0,
            deflection_rate: 0.0,
        }
    }
    /// Typical steel driveshaft for a compact passenger car.
    pub fn typical_steel() -> Self {
        Self::new(8000.0, 80.0)
    }
    /// Transmitted torque based on current deflection and rate (N·m).
    ///
    /// `T = k * θ + c * dθ/dt`
    pub fn transmitted_torque(&self, engine_speed_rad_s: f64, wheel_speed_rad_s: f64) -> f64 {
        let speed_diff = engine_speed_rad_s - wheel_speed_rad_s;
        self.stiffness * self.deflection + self.damping * speed_diff
    }
    /// Advance the torsional state by `dt` seconds using Euler integration.
    ///
    /// `engine_speed_rad_s` and `wheel_speed_rad_s` in rad/s.
    pub fn step(&mut self, engine_speed_rad_s: f64, wheel_speed_rad_s: f64, dt: f64) {
        let speed_diff = engine_speed_rad_s - wheel_speed_rad_s;
        self.deflection += speed_diff * dt;
        self.deflection_rate = speed_diff;
    }
    /// Driveline natural frequency (Hz) for a given effective inertia (kg·m²).
    ///
    /// `f_n = (1/2π) * sqrt(k / J)`
    pub fn natural_frequency(&self, inertia_kg_m2: f64) -> f64 {
        if inertia_kg_m2 < 1e-12 {
            return 0.0;
        }
        (1.0 / (2.0 * std::f64::consts::PI)) * (self.stiffness / inertia_kg_m2).sqrt()
    }
    /// Reset deflection to zero (simulate slip-release event).
    pub fn reset_deflection(&mut self) {
        self.deflection = 0.0;
        self.deflection_rate = 0.0;
    }
}
/// Engine torque curve represented as a degree-N polynomial fit.
///
/// `T(rpm) = a0 + a1*rpm + a2*rpm² + ...`
///
/// Coefficients are stored from lowest to highest degree.
#[derive(Debug, Clone)]
pub struct PolynomialTorqueCurve {
    /// Polynomial coefficients `[a0, a1, a2, ...]`.
    pub coefficients: Vec<f64>,
    /// Minimum RPM (idle).
    pub min_rpm: f64,
    /// Maximum RPM (redline).
    pub max_rpm: f64,
}
impl PolynomialTorqueCurve {
    /// Create a new polynomial torque curve.
    pub fn new(coefficients: Vec<f64>, min_rpm: f64, max_rpm: f64) -> Self {
        Self {
            coefficients,
            min_rpm: min_rpm.max(0.0),
            max_rpm: max_rpm.max(min_rpm + 1.0),
        }
    }
    /// Typical naturally-aspirated 2.0 L polynomial fit.
    ///
    /// Approximate fit: T ≈ -0.0014·rpm² + 8.5·rpm - 5000 (shifted/scaled)
    /// Peak near 3000 RPM.
    pub fn typical_na() -> Self {
        Self::new(vec![-4800.0, 7.0, -0.0012], 800.0, 7000.0)
    }
    /// Evaluate torque at given RPM (clamped to \[min_rpm, max_rpm\]).
    pub fn torque_at_rpm(&self, rpm: f64) -> f64 {
        let r = rpm.clamp(self.min_rpm, self.max_rpm);
        let mut torque = 0.0;
        let mut power = 1.0;
        for &coeff in &self.coefficients {
            torque += coeff * power;
            power *= r;
        }
        torque.max(0.0)
    }
    /// RPM at peak torque (found by golden-section search over \[min_rpm, max_rpm\]).
    pub fn peak_torque_rpm(&self) -> f64 {
        let mut lo = self.min_rpm;
        let mut hi = self.max_rpm;
        let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
        for _ in 0..60 {
            let m1 = hi - gr * (hi - lo);
            let m2 = lo + gr * (hi - lo);
            if self.torque_at_rpm(m1) < self.torque_at_rpm(m2) {
                lo = m1;
            } else {
                hi = m2;
            }
        }
        (lo + hi) * 0.5
    }
}
/// Open differential: torque split equally, wheel speeds free to differ.
#[derive(Debug, Clone)]
pub struct OpenDifferential {
    /// Ring-gear ratio (typically 1.0).
    pub ring_gear_ratio: f64,
    /// Current angular velocity of the left output shaft (rad/s).
    pub left_wheel_speed: f64,
    /// Current angular velocity of the right output shaft (rad/s).
    pub right_wheel_speed: f64,
}
impl OpenDifferential {
    /// Create a new open differential with both wheel speeds at zero.
    pub fn new() -> Self {
        Self {
            ring_gear_ratio: 1.0,
            left_wheel_speed: 0.0,
            right_wheel_speed: 0.0,
        }
    }
    /// Split input torque equally between the two output shafts.
    pub fn torque_split(&self, input_torque: f64) -> (f64, f64) {
        let half = input_torque * 0.5;
        (half, half)
    }
    /// Average of both wheel speeds.
    pub fn average_speed(&self) -> f64 {
        (self.left_wheel_speed + self.right_wheel_speed) * 0.5
    }
    /// Speed difference: left − right.
    pub fn speed_diff(&self) -> f64 {
        self.left_wheel_speed - self.right_wheel_speed
    }
}
/// Manual H-pattern or sequential gearbox (companion to [`Gearbox`]).
#[derive(Debug, Clone)]
pub struct GearboxManual {
    /// Forward gear ratios (index 0 = 1st gear).
    pub gear_ratios: Vec<f64>,
    /// Reverse gear ratio (positive; direction encoded by sign of current_gear).
    pub reverse_ratio: f64,
    /// Final drive (axle differential) ratio.
    pub final_drive: f64,
    /// Currently selected gear: 0 = neutral, 1..N = forward, -1 = reverse.
    pub current_gear: i32,
}
impl GearboxManual {
    /// Create a new manual gearbox.
    pub fn new(gear_ratios: Vec<f64>, final_drive: f64) -> Self {
        Self {
            reverse_ratio: 3.5,
            gear_ratios,
            final_drive,
            current_gear: 1,
        }
    }
    /// Typical 6-speed gearbox with common sport-car ratios.
    pub fn six_speed() -> Self {
        Self::new(vec![3.82, 2.20, 1.52, 1.15, 0.89, 0.73], 3.73)
    }
    /// Combined gear ratio (gear × final drive).  Returns 0.0 for neutral.
    pub fn total_ratio(&self) -> f64 {
        match self.current_gear {
            0 => 0.0,
            -1 => -self.reverse_ratio * self.final_drive,
            g if g >= 1 && (g as usize) <= self.gear_ratios.len() => {
                self.gear_ratios[(g - 1) as usize] * self.final_drive
            }
            _ => 0.0,
        }
    }
    /// Shift up one gear.  Returns `true` if the gear actually changed.
    pub fn shift_up(&mut self) -> bool {
        let max = self.max_gear();
        if self.current_gear < max {
            self.current_gear += 1;
            true
        } else {
            false
        }
    }
    /// Shift down one gear.  Returns `true` if the gear actually changed.
    pub fn shift_down(&mut self) -> bool {
        if self.current_gear > -1 {
            self.current_gear -= 1;
            true
        } else {
            false
        }
    }
    /// Select an arbitrary gear (clamped to valid range).
    pub fn shift_to(&mut self, gear: i32) {
        let max = self.max_gear();
        self.current_gear = gear.clamp(-1, max);
    }
    /// Returns `true` when in neutral.
    pub fn is_neutral(&self) -> bool {
        self.current_gear == 0
    }
    /// Highest available forward gear number.
    pub fn max_gear(&self) -> i32 {
        self.gear_ratios.len() as i32
    }
    /// Wheel RPM for a given engine RPM and current gear.
    pub fn wheel_rpm(&self, engine_rpm: f64) -> f64 {
        let ratio = self.total_ratio().abs();
        if ratio < 1e-10 {
            return 0.0;
        }
        engine_rpm / ratio
    }
    /// Engine RPM that corresponds to a given wheel RPM in the current gear.
    pub fn engine_rpm_from_wheel(&self, wheel_rpm: f64) -> f64 {
        wheel_rpm * self.total_ratio().abs()
    }
}
/// A single gear with a reduction ratio and mechanical efficiency.
#[derive(Debug, Clone)]
pub struct Gear {
    /// Speed reduction ratio (input ÷ output speed).
    pub ratio: f64,
    /// Mechanical efficiency (0..1), typically ~0.97.
    pub efficiency: f64,
}
impl Gear {
    /// Create a new gear.
    pub fn new(ratio: f64, efficiency: f64) -> Self {
        Self { ratio, efficiency }
    }
}
/// Detailed engine torque curve with power/peak helpers.
#[derive(Debug, Clone)]
pub struct EngineTorqueCurve {
    /// RPM sample points (sorted ascending).
    pub rpm_points: Vec<f64>,
    /// Torque (N·m) at each corresponding RPM point.
    pub torque_points: Vec<f64>,
    /// Maximum operating RPM (redline).
    pub max_rpm: f64,
    /// Minimum operating RPM (idle).
    pub idle_rpm: f64,
}
impl EngineTorqueCurve {
    /// Create a new torque curve from parallel RPM and torque vectors.
    pub fn new(rpm_points: Vec<f64>, torque_points: Vec<f64>) -> Self {
        assert_eq!(
            rpm_points.len(),
            torque_points.len(),
            "rpm_points and torque_points must have equal length"
        );
        let mut pairs: Vec<(f64, f64)> = rpm_points.into_iter().zip(torque_points).collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let idle_rpm = pairs.first().map(|p| p.0).unwrap_or(800.0);
        let max_rpm = pairs.last().map(|p| p.0).unwrap_or(7000.0);
        let rpm_points = pairs.iter().map(|p| p.0).collect();
        let torque_points = pairs.iter().map(|p| p.1).collect();
        Self {
            rpm_points,
            torque_points,
            max_rpm,
            idle_rpm,
        }
    }
    /// Typical naturally aspirated 2.0 L engine torque curve.
    pub fn typical_na() -> Self {
        Self::new(
            vec![
                800.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0,
            ],
            vec![120.0, 160.0, 185.0, 200.0, 195.0, 185.0, 165.0, 130.0],
        )
    }
    /// Torque at given RPM via linear interpolation.
    pub fn torque_at_rpm(&self, rpm: f64) -> f64 {
        if self.rpm_points.is_empty() {
            return 0.0;
        }
        if self.rpm_points.len() == 1 {
            return self.torque_points[0];
        }
        let rpm_clamped = rpm.clamp(self.idle_rpm, self.max_rpm);
        for i in 0..self.rpm_points.len() - 1 {
            let r0 = self.rpm_points[i];
            let r1 = self.rpm_points[i + 1];
            let t0 = self.torque_points[i];
            let t1 = self.torque_points[i + 1];
            if rpm_clamped >= r0 && rpm_clamped <= r1 {
                if (r1 - r0).abs() < 1e-10 {
                    return t0;
                }
                let frac = (rpm_clamped - r0) / (r1 - r0);
                return t0 + frac * (t1 - t0);
            }
        }
        *self.torque_points.last().unwrap_or(&0.0)
    }
    /// Power at given RPM (watts).
    pub fn power_at_rpm(&self, rpm: f64) -> f64 {
        let torque = self.torque_at_rpm(rpm);
        torque * rpm * 2.0 * std::f64::consts::PI / 60.0
    }
    /// RPM at which torque is highest.
    pub fn peak_torque_rpm(&self) -> f64 {
        self.rpm_points
            .iter()
            .zip(self.torque_points.iter())
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(rpm, _)| *rpm)
            .unwrap_or(self.idle_rpm)
    }
    /// RPM at which power output is highest.
    pub fn peak_power_rpm(&self) -> f64 {
        self.rpm_points
            .iter()
            .max_by(|&&a, &&b| {
                self.power_at_rpm(a)
                    .partial_cmp(&self.power_at_rpm(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied()
            .unwrap_or(self.max_rpm)
    }
}
/// A piecewise-linear engine torque curve defined by parallel RPM and torque vectors.
///
/// Supports linear interpolation between sample points and clamping
/// to the defined RPM range.
#[derive(Debug, Clone)]
pub struct EngineCurve {
    /// RPM sample points (sorted ascending).
    pub rpm_points: Vec<f64>,
    /// Torque (N·m) at each corresponding RPM point.
    pub torque_points: Vec<f64>,
}
impl EngineCurve {
    /// Create a new engine curve from parallel RPM and torque vectors.
    ///
    /// Returns `Err` if the vectors have different lengths or are empty.
    pub fn new(rpm: Vec<f64>, torque: Vec<f64>) -> Result<Self, String> {
        if rpm.len() != torque.len() {
            return Err(format!(
                "rpm_points ({}) and torque_points ({}) must have equal length",
                rpm.len(),
                torque.len()
            ));
        }
        if rpm.is_empty() {
            return Err("rpm_points must not be empty".to_string());
        }
        let mut pairs: Vec<(f64, f64)> = rpm.into_iter().zip(torque).collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let rpm_points = pairs.iter().map(|p| p.0).collect();
        let torque_points = pairs.iter().map(|p| p.1).collect();
        Ok(Self {
            rpm_points,
            torque_points,
        })
    }
    /// Minimum RPM in the curve (idle).
    fn min_rpm(&self) -> f64 {
        *self.rpm_points.first().unwrap_or(&800.0)
    }
    /// Maximum RPM in the curve (redline).
    fn max_rpm(&self) -> f64 {
        *self.rpm_points.last().unwrap_or(&7000.0)
    }
    /// Torque at given RPM via linear interpolation, clamped to the defined range.
    pub fn torque_at_rpm(&self, rpm: f64) -> f64 {
        if self.rpm_points.is_empty() {
            return 0.0;
        }
        if self.rpm_points.len() == 1 {
            return self.torque_points[0];
        }
        let rpm_c = rpm.clamp(self.min_rpm(), self.max_rpm());
        for i in 0..self.rpm_points.len() - 1 {
            let r0 = self.rpm_points[i];
            let r1 = self.rpm_points[i + 1];
            let t0 = self.torque_points[i];
            let t1 = self.torque_points[i + 1];
            if rpm_c >= r0 && rpm_c <= r1 {
                if (r1 - r0).abs() < 1e-10 {
                    return t0;
                }
                let frac = (rpm_c - r0) / (r1 - r0);
                return t0 + frac * (t1 - t0);
            }
        }
        *self.torque_points.last().unwrap_or(&0.0)
    }
    /// Power at given RPM: P = τ · ω = τ · rpm · 2π / 60 (watts).
    pub fn power_at_rpm(&self, rpm: f64) -> f64 {
        self.torque_at_rpm(rpm) * rpm * 2.0 * std::f64::consts::PI / 60.0
    }
    /// Returns (rpm, torque_Nm) at peak torque.
    pub fn peak_torque(&self) -> (f64, f64) {
        self.rpm_points
            .iter()
            .zip(self.torque_points.iter())
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(&r, &t)| (r, t))
            .unwrap_or((self.min_rpm(), 0.0))
    }
    /// Returns (rpm, power_W) at peak power.
    pub fn peak_power(&self) -> (f64, f64) {
        let (rpm, _) = self
            .rpm_points
            .iter()
            .map(|&r| (r, self.power_at_rpm(r)))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((self.max_rpm(), 0.0));
        (rpm, self.power_at_rpm(rpm))
    }
    /// Representative naturally aspirated 4-cylinder torque curve.
    pub fn typical_4cylinder() -> Self {
        Self::new(
            vec![
                800.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0,
            ],
            vec![120.0, 160.0, 185.0, 200.0, 195.0, 185.0, 165.0, 130.0],
        )
        .expect("static curve is valid")
    }
}
/// Legacy differential enum.  Prefer [`Differential`] for new code.
#[derive(Debug, Clone, Default)]
pub enum DifferentialMode {
    /// Open differential: torque split equally regardless of wheel speeds.
    #[default]
    Open,
    /// Locked differential: both wheels rotate at the same speed.
    Locked,
    /// Limited-slip differential with a torque bias ratio.
    LimitedSlip {
        /// Maximum torque bias ratio.
        torque_bias: Real,
    },
}
impl DifferentialMode {
    /// Split total torque between two driven wheels.
    pub fn split_torque(
        &self,
        total_torque: Real,
        speed_left: Real,
        speed_right: Real,
    ) -> (Real, Real) {
        match self {
            DifferentialMode::Open | DifferentialMode::Locked => {
                let half = total_torque * 0.5;
                (half, half)
            }
            DifferentialMode::LimitedSlip { torque_bias } => {
                let speed_diff = (speed_left - speed_right).abs();
                let avg_speed = (speed_left.abs() + speed_right.abs()) * 0.5;
                if avg_speed < 1e-10 {
                    let half = total_torque * 0.5;
                    return (half, half);
                }
                let ratio = (speed_diff / avg_speed).min(1.0);
                let bias = 1.0 + ratio * (torque_bias - 1.0);
                let (slow_factor, fast_factor) = if bias > 1e-10 {
                    let fast = 1.0 / (1.0 + bias);
                    let slow = bias / (1.0 + bias);
                    (slow, fast)
                } else {
                    (0.5, 0.5)
                };
                if speed_left.abs() <= speed_right.abs() {
                    (total_torque * slow_factor, total_torque * fast_factor)
                } else {
                    (total_torque * fast_factor, total_torque * slow_factor)
                }
            }
        }
    }
}
/// Engine braking model.
///
/// Computes the resistive torque due to engine pumping losses when
/// the throttle is closed, proportional to RPM above idle.
#[derive(Debug, Clone)]
pub struct EngineBraking {
    /// Braking torque coefficient (N·m / RPM above idle).
    pub coefficient: f64,
    /// Idle RPM (below which no engine braking occurs).
    pub idle_rpm: f64,
}
impl EngineBraking {
    /// Create an engine braking model.
    ///
    /// `coefficient`: N·m per RPM above idle (typical: 0.005 – 0.02 N·m/RPM).
    pub fn new(coefficient: f64, idle_rpm: f64) -> Self {
        Self {
            coefficient: coefficient.max(0.0),
            idle_rpm: idle_rpm.max(0.0),
        }
    }
    /// Default naturally-aspirated 4-cylinder engine braking model.
    pub fn default_na() -> Self {
        Self::new(0.012, 800.0)
    }
    /// Braking torque at given RPM (N·m, always ≥ 0).
    ///
    /// Returns 0 when RPM is at or below idle.
    pub fn braking_torque(&self, rpm: f64) -> f64 {
        let excess_rpm = (rpm - self.idle_rpm).max(0.0);
        self.coefficient * excess_rpm
    }
    /// Net applied torque when throttle is fully released (negative = decelerating).
    ///
    /// Use this to oppose vehicle acceleration when coasting.
    pub fn deceleration_torque(&self, rpm: f64) -> f64 {
        -self.braking_torque(rpm)
    }
}
/// Torque converter with a lock-up clutch state machine.
///
/// Extends [`TorqueConverter`] with hysteresis-controlled lock-up logic.
/// The clutch engages above `lockup_speed_ratio` and disengages below
/// `lockup_speed_ratio - hysteresis`.
#[derive(Debug, Clone)]
pub struct TorqueConverterLockup {
    /// Base torque converter.
    pub converter: TorqueConverter,
    /// Speed ratio above which the lock-up clutch begins to close.
    pub lockup_speed_ratio: f64,
    /// Hysteresis band below `lockup_speed_ratio` at which the clutch opens.
    pub hysteresis: f64,
    /// Engagement time constant (s).
    pub engage_time: f64,
    /// Current lock-up state.
    pub state: LockupState,
    /// Current clutch engagement fraction \[0..1\].
    pub clutch_engagement: f64,
}
impl TorqueConverterLockup {
    /// Create a new torque converter with lock-up.
    pub fn new(stall_ratio: f64, capacity_factor_ref: f64) -> Self {
        Self {
            converter: TorqueConverter::new(stall_ratio, capacity_factor_ref),
            lockup_speed_ratio: 0.85,
            hysteresis: 0.05,
            engage_time: 0.4,
            state: LockupState::Open,
            clutch_engagement: 0.0,
        }
    }
    /// Typical passenger-car torque converter with lock-up (stall ~2.0).
    pub fn typical_passenger_car() -> Self {
        Self::new(2.0, 200.0)
    }
    /// Update lock-up state given current pump and turbine RPM and time step.
    ///
    /// Returns the effective speed ratio after clutch modulation.
    pub fn update(&mut self, pump_rpm: f64, turbine_rpm: f64, dt: f64) -> f64 {
        let sr = if pump_rpm.abs() < 1.0 {
            0.0
        } else {
            (turbine_rpm / pump_rpm).clamp(0.0, 1.0)
        };
        match self.state {
            LockupState::Open => {
                if sr >= self.lockup_speed_ratio {
                    self.state = LockupState::Engaging;
                }
            }
            LockupState::Engaging => {
                self.clutch_engagement = (self.clutch_engagement + dt / self.engage_time).min(1.0);
                if self.clutch_engagement >= 1.0 {
                    self.state = LockupState::Locked;
                }
                if sr < self.lockup_speed_ratio - self.hysteresis {
                    self.state = LockupState::Open;
                    self.clutch_engagement = 0.0;
                }
            }
            LockupState::Locked => {
                if sr < self.lockup_speed_ratio - self.hysteresis {
                    self.state = LockupState::Open;
                    self.clutch_engagement = 0.0;
                }
            }
        }
        sr
    }
    /// Effective output torque considering clutch engagement.
    ///
    /// When locked: passes input torque 1:1.
    /// When open/engaging: blends converter output with direct coupling.
    pub fn effective_output_torque(&self, input_torque: f64, speed_ratio: f64) -> f64 {
        let converter_out = self.converter.output_torque(input_torque, speed_ratio);
        converter_out * (1.0 - self.clutch_engagement) + input_torque * self.clutch_engagement
    }
    /// Returns `true` if the lock-up clutch is fully engaged.
    pub fn is_locked(&self) -> bool {
        self.state == LockupState::Locked
    }
}
/// Friction clutch between engine and gearbox input shaft.
#[derive(Debug, Clone)]
pub struct Clutch {
    /// Maximum torque the clutch can transmit when fully engaged (N·m).
    pub max_torque: f64,
    /// Engagement factor: 0.0 = fully open, 1.0 = fully engaged.
    pub engagement: f64,
    /// Angular velocity difference across the clutch (rad/s).
    pub slip_speed: f64,
}
impl Clutch {
    /// Create a new clutch with given maximum torque capacity, starting disengaged.
    pub fn new(max_torque: f64) -> Self {
        Self {
            max_torque,
            engagement: 0.0,
            slip_speed: 0.0,
        }
    }
    /// Fully engage the clutch.
    pub fn engage(&mut self) {
        self.engagement = 1.0;
    }
    /// Fully disengage the clutch.
    pub fn disengage(&mut self) {
        self.engagement = 0.0;
    }
    /// Set partial engagement (clamped 0..1).
    pub fn partial(&mut self, alpha: f64) {
        self.engagement = alpha.clamp(0.0, 1.0);
    }
    /// Torque transmitted through the clutch.
    ///
    /// Returns `min(engagement × max_torque, |engine_torque|)`, preserving sign.
    pub fn transmitted_torque(&self, engine_torque: f64) -> f64 {
        let capacity = self.engagement * self.max_torque;
        engine_torque.clamp(-capacity, capacity)
    }
    /// Heat generated by clutch slip (W).
    ///
    /// Q = |slip_speed| × |transmitted_torque|
    pub fn heat_generation(&self, engine_torque: f64) -> f64 {
        self.slip_speed.abs() * self.transmitted_torque(engine_torque).abs()
    }
}
/// Automatic transmission with a planetary gear set, torque converter, and
/// closed-loop shift strategy.
///
/// Models a simplified Ravigneaux or Simpson-type 4–6 speed automatic.
#[derive(Debug, Clone)]
pub struct AutomaticTransmission {
    /// Planetary gear set providing the gear ratios.
    pub planetary: PlanetaryGearSet,
    /// Torque converter with lock-up.
    pub torque_converter: TorqueConverterLockup,
    /// Forward gear ratios (index 0 = 1st).
    pub gear_ratios: Vec<f64>,
    /// Final drive ratio.
    pub final_drive: f64,
    /// Currently selected gear (1-based, 0 = neutral).
    pub current_gear: usize,
    /// Shift mode.
    pub shift_mode: ShiftMode,
    /// Current shift phase.
    pub shift_phase: ShiftPhase,
    /// Time remaining in current shift (s).
    pub shift_time_remaining: f64,
    /// Shift duration (s).
    pub shift_duration: f64,
    /// Pending target gear (used during shifts).
    pub pending_gear: usize,
    /// Upshift RPM thresholds for Economy mode (index 0 = 1→2, etc.).
    pub economy_upshift_rpm: Vec<f64>,
    /// Downshift RPM thresholds for Economy mode.
    pub economy_downshift_rpm: Vec<f64>,
    /// Upshift RPM thresholds for Sport mode.
    pub sport_upshift_rpm: Vec<f64>,
    /// Downshift RPM thresholds for Sport mode.
    pub sport_downshift_rpm: Vec<f64>,
    /// Mechanical efficiency of the transmission.
    pub efficiency: f64,
}
impl AutomaticTransmission {
    /// Create a new automatic transmission.
    pub fn new(
        planetary: PlanetaryGearSet,
        gear_ratios: Vec<f64>,
        final_drive: f64,
        economy_upshift_rpm: Vec<f64>,
        economy_downshift_rpm: Vec<f64>,
        sport_upshift_rpm: Vec<f64>,
        sport_downshift_rpm: Vec<f64>,
    ) -> Self {
        let n = gear_ratios.len();
        Self {
            planetary,
            torque_converter: TorqueConverterLockup::typical_passenger_car(),
            gear_ratios,
            final_drive,
            current_gear: 1,
            shift_mode: ShiftMode::Economy,
            shift_phase: ShiftPhase::Steady,
            shift_time_remaining: 0.0,
            shift_duration: 0.35,
            pending_gear: 1,
            economy_upshift_rpm: if economy_upshift_rpm.len() == n - 1 {
                economy_upshift_rpm
            } else {
                vec![2200.0; n.saturating_sub(1)]
            },
            economy_downshift_rpm: if economy_downshift_rpm.len() == n - 1 {
                economy_downshift_rpm
            } else {
                vec![1200.0; n.saturating_sub(1)]
            },
            sport_upshift_rpm: if sport_upshift_rpm.len() == n - 1 {
                sport_upshift_rpm
            } else {
                vec![5500.0; n.saturating_sub(1)]
            },
            sport_downshift_rpm: if sport_downshift_rpm.len() == n - 1 {
                sport_downshift_rpm
            } else {
                vec![2800.0; n.saturating_sub(1)]
            },
            efficiency: 0.92,
        }
    }
    /// Typical 6-speed automatic (similar to ZF 6HP).
    pub fn zf_6hp() -> Self {
        let planetary = PlanetaryGearSet::ravigneaux();
        let gear_ratios = vec![4.17, 2.34, 1.52, 1.14, 0.87, 0.69];
        let final_drive = 3.46;
        let eco_up = vec![1800.0, 1900.0, 2000.0, 2100.0, 2200.0];
        let eco_dn = vec![1200.0, 1300.0, 1400.0, 1500.0, 1600.0];
        let spt_up = vec![4800.0, 5000.0, 5200.0, 5400.0, 5600.0];
        let spt_dn = vec![2500.0, 2700.0, 2900.0, 3000.0, 3200.0];
        Self::new(
            planetary,
            gear_ratios,
            final_drive,
            eco_up,
            eco_dn,
            spt_up,
            spt_dn,
        )
    }
    /// Total gear ratio for the current gear × final drive.
    pub fn total_ratio(&self) -> f64 {
        if self.current_gear == 0 || self.current_gear > self.gear_ratios.len() {
            return 0.0;
        }
        self.gear_ratios[self.current_gear - 1] * self.final_drive
    }
    /// Output torque delivered to the driveshaft (N·m).
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        let ratio = self.total_ratio();
        input_torque * ratio * self.efficiency
    }
    /// Determine whether an upshift is requested based on current shift mode.
    fn wants_upshift(&self, engine_rpm: f64, throttle: f64) -> bool {
        if self.current_gear >= self.gear_ratios.len() {
            return false;
        }
        if self.shift_phase != ShiftPhase::Steady {
            return false;
        }
        let idx = self.current_gear.saturating_sub(1);
        match self.shift_mode {
            ShiftMode::Economy => {
                let threshold = self
                    .economy_upshift_rpm
                    .get(idx)
                    .copied()
                    .unwrap_or(f64::INFINITY);
                let adjusted = threshold - (1.0 - throttle) * 200.0;
                engine_rpm > adjusted
            }
            ShiftMode::Sport => {
                let threshold = self
                    .sport_upshift_rpm
                    .get(idx)
                    .copied()
                    .unwrap_or(f64::INFINITY);
                engine_rpm > threshold
            }
            ShiftMode::Manual => false,
        }
    }
    /// Determine whether a downshift is requested.
    fn wants_downshift(&self, engine_rpm: f64, throttle: f64) -> bool {
        if self.current_gear <= 1 {
            return false;
        }
        if self.shift_phase != ShiftPhase::Steady {
            return false;
        }
        let idx = self.current_gear.saturating_sub(2);
        match self.shift_mode {
            ShiftMode::Economy => {
                let threshold = self.economy_downshift_rpm.get(idx).copied().unwrap_or(0.0);
                let kick_down = throttle > 0.9 && engine_rpm < threshold * 1.8;
                engine_rpm < threshold || kick_down
            }
            ShiftMode::Sport => {
                let threshold = self.sport_downshift_rpm.get(idx).copied().unwrap_or(0.0);
                engine_rpm < threshold
            }
            ShiftMode::Manual => false,
        }
    }
    /// Run the automatic shift strategy for one time step.
    ///
    /// Updates current gear and shift phase. Returns `true` if a gear change
    /// was initiated this step.
    ///
    /// # Arguments
    /// * `engine_rpm`  – current engine speed (RPM)
    /// * `throttle`    – throttle position \[0, 1\]
    /// * `dt`          – time step (s)
    pub fn step_shift_strategy(&mut self, engine_rpm: f64, throttle: f64, dt: f64) -> bool {
        if self.shift_phase != ShiftPhase::Steady {
            self.shift_time_remaining -= dt;
            if self.shift_time_remaining <= 0.0 {
                self.current_gear = self.pending_gear;
                self.shift_phase = ShiftPhase::Steady;
                self.shift_time_remaining = 0.0;
            }
            return false;
        }
        if self.wants_upshift(engine_rpm, throttle) {
            let target = (self.current_gear + 1).min(self.gear_ratios.len());
            self.pending_gear = target;
            self.shift_phase = ShiftPhase::Torque;
            self.shift_time_remaining = self.shift_duration;
            return true;
        }
        if self.wants_downshift(engine_rpm, throttle) {
            let target = self.current_gear.saturating_sub(1).max(1);
            self.pending_gear = target;
            self.shift_phase = ShiftPhase::Inertia;
            self.shift_time_remaining = self.shift_duration * 0.8;
            return true;
        }
        false
    }
    /// Manual upshift (for ShiftMode::Manual or paddle shifters).
    pub fn manual_upshift(&mut self) {
        if self.current_gear < self.gear_ratios.len() && self.shift_phase == ShiftPhase::Steady {
            self.pending_gear = self.current_gear + 1;
            self.shift_phase = ShiftPhase::Torque;
            self.shift_time_remaining = self.shift_duration;
        }
    }
    /// Manual downshift.
    pub fn manual_downshift(&mut self) {
        if self.current_gear > 1 && self.shift_phase == ShiftPhase::Steady {
            self.pending_gear = self.current_gear - 1;
            self.shift_phase = ShiftPhase::Inertia;
            self.shift_time_remaining = self.shift_duration * 0.8;
        }
    }
    /// Returns `true` if a shift is currently in progress.
    pub fn is_shifting(&self) -> bool {
        self.shift_phase != ShiftPhase::Steady
    }
    /// Torque reduction factor during a shift (simulates momentary power interruption).
    ///
    /// Returns 1.0 in steady state, 0.3 during the torque phase.
    pub fn torque_reduction_factor(&self) -> f64 {
        match self.shift_phase {
            ShiftPhase::Steady | ShiftPhase::PowerOn => 1.0,
            ShiftPhase::Torque => 0.3,
            ShiftPhase::Inertia => 0.6,
        }
    }
}
/// Brake-specific fuel consumption lookup map.
///
/// BSFC (g/kWh) as a function of engine speed (RPM) and throttle position.
/// Minimum BSFC (best efficiency) occurs at mid-RPM, moderate load.
#[derive(Debug, Clone)]
pub struct BsfcMap {
    /// Optimal BSFC value (g/kWh) at the best efficiency point.
    pub bsfc_optimal: f64,
    /// RPM at the best efficiency point.
    pub best_rpm: f64,
    /// Throttle at the best efficiency point \[0, 1\].
    pub best_throttle: f64,
    /// BSFC penalty per RPM deviation from optimal (g/kWh per RPM).
    pub rpm_penalty: f64,
    /// BSFC penalty per throttle deviation from optimal (g/kWh per unit).
    pub throttle_penalty: f64,
}
impl BsfcMap {
    /// Create a BSFC map with given parameters.
    pub fn new(
        bsfc_optimal: f64,
        best_rpm: f64,
        best_throttle: f64,
        rpm_penalty: f64,
        throttle_penalty: f64,
    ) -> Self {
        Self {
            bsfc_optimal,
            best_rpm,
            best_throttle,
            rpm_penalty,
            throttle_penalty,
        }
    }
    /// Typical naturally-aspirated gasoline engine BSFC map.
    ///
    /// Best efficiency around 2500 RPM / 70% throttle: ~245 g/kWh.
    pub fn default_gasoline() -> Self {
        Self::new(245.0, 2500.0, 0.7, 0.008, 80.0)
    }
    /// Lookup BSFC (g/kWh) at given RPM and throttle.
    pub fn lookup(&self, rpm: f64, throttle: f64) -> f64 {
        let throttle_c = throttle.clamp(0.0, 1.0);
        let rpm_dev = (rpm - self.best_rpm).abs() / 1000.0;
        let thr_dev = (throttle_c - self.best_throttle).abs();
        self.bsfc_optimal
            + self.rpm_penalty * rpm_dev.powi(2) * 1000.0
            + self.throttle_penalty * thr_dev.powi(2)
    }
    /// Instantaneous fuel flow rate (g/s) at given operating point.
    ///
    /// `fuel_flow = BSFC * power / 3_600_000`  (BSFC in g/kWh, power in W → g/s)
    pub fn fuel_flow_g_per_s(&self, rpm: f64, throttle: f64, torque_nm: f64) -> f64 {
        let omega = rpm * 2.0 * std::f64::consts::PI / 60.0;
        let power_w = (torque_nm * omega).max(0.0);
        let bsfc = self.lookup(rpm, throttle);
        bsfc * power_w / 3_600_000.0
    }
}
/// Non-linear driveshaft model with angular backlash (lash) and a torsional
/// torque limit (shear-pin / overload protection).
///
/// The backlash is a dead-band in the torsional deflection: no torque is
/// transmitted until the deflection exceeds `±lash_rad`.
#[derive(Debug, Clone)]
pub struct DriveshaftNonlinear {
    /// Torsional stiffness (N·m/rad).
    pub stiffness: f64,
    /// Torsional damping (N·m·s/rad).
    pub damping: f64,
    /// Backlash half-band (rad); no torque for |θ| < lash_rad.
    pub lash_rad: f64,
    /// Maximum transmissible torque (N·m); torque is clamped at this level.
    pub torque_limit: f64,
    /// Current torsional deflection (rad).
    pub deflection: f64,
}
impl DriveshaftNonlinear {
    /// Create a non-linear driveshaft.
    pub fn new(stiffness: f64, damping: f64, lash_rad: f64, torque_limit: f64) -> Self {
        Self {
            stiffness: stiffness.max(1.0),
            damping: damping.max(0.0),
            lash_rad: lash_rad.abs(),
            torque_limit: torque_limit.max(0.0),
            deflection: 0.0,
        }
    }
    /// Typical front-drive half-shaft (lash ≈ 1° = 0.0175 rad).
    pub fn typical_half_shaft() -> Self {
        Self::new(6000.0, 60.0, 0.0175, 800.0)
    }
    /// Transmitted torque considering backlash and limit.
    ///
    /// `T = clamp(spring_torque + damping_torque, -limit, +limit)`
    /// where spring torque is 0 inside the lash band.
    pub fn transmitted_torque(&self, engine_speed_rad_s: f64, wheel_speed_rad_s: f64) -> f64 {
        let speed_diff = engine_speed_rad_s - wheel_speed_rad_s;
        let damping_t = self.damping * speed_diff;
        let spring_t = if self.deflection.abs() <= self.lash_rad {
            0.0
        } else {
            let active = self.deflection - self.deflection.signum() * self.lash_rad;
            self.stiffness * active
        };
        let total = spring_t + damping_t;
        total.clamp(-self.torque_limit, self.torque_limit)
    }
    /// Advance the deflection state by `dt` seconds (Euler integration).
    pub fn step(&mut self, engine_speed_rad_s: f64, wheel_speed_rad_s: f64, dt: f64) {
        let speed_diff = engine_speed_rad_s - wheel_speed_rad_s;
        self.deflection += speed_diff * dt;
    }
    /// Reset deflection to zero (simulates torque reversal / lash take-up).
    pub fn reset_deflection(&mut self) {
        self.deflection = 0.0;
    }
}
/// Shift strategy for the automatic transmission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShiftMode {
    /// Standard efficiency-oriented shift schedule.
    Economy,
    /// Sport/performance shift schedule (holds gears longer, higher shift RPM).
    Sport,
    /// Manual sequential override.
    Manual,
}
/// Continuously variable transmission (CVT) model.
///
/// A CVT adjusts its ratio continuously between a minimum (overdrive) and
/// maximum (low gear) ratio to keep the engine at an optimal operating point.
///
/// The ratio is controlled by a hydraulic actuator modelled as a first-order
/// lag system.
#[derive(Debug, Clone)]
pub struct CvtTransmission {
    /// Minimum gear ratio (highest speed / overdrive).
    pub ratio_min: f64,
    /// Maximum gear ratio (lowest speed / low gear).
    pub ratio_max: f64,
    /// Final drive ratio.
    pub final_drive: f64,
    /// Current CVT ratio.
    pub current_ratio: f64,
    /// First-order time constant for ratio change (s).
    pub ratio_time_const: f64,
    /// Mechanical efficiency.
    pub efficiency: f64,
    /// Target RPM for the engine (CVT tries to keep engine here).
    pub target_engine_rpm: f64,
}
impl CvtTransmission {
    /// Create a new CVT.
    pub fn new(ratio_min: f64, ratio_max: f64, final_drive: f64) -> Self {
        let ratio_min = ratio_min.max(0.3);
        let ratio_max = ratio_max.max(ratio_min + 0.1);
        Self {
            ratio_min,
            ratio_max,
            final_drive,
            current_ratio: ratio_max,
            ratio_time_const: 0.5,
            efficiency: 0.88,
            target_engine_rpm: 2200.0,
        }
    }
    /// Typical Jatco JF011E CVT (common in Nissan/Mitsubishi vehicles).
    pub fn jatco_jf011e() -> Self {
        Self::new(0.43, 2.35, 4.04)
    }
    /// Total ratio (CVT ratio × final drive).
    pub fn total_ratio(&self) -> f64 {
        self.current_ratio * self.final_drive
    }
    /// Output torque to the driveshaft (N·m).
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        input_torque * self.total_ratio() * self.efficiency
    }
    /// Compute the ideal CVT ratio to keep the engine at `target_engine_rpm`
    /// given the current wheel speed (rad/s).
    ///
    /// `ideal_ratio = target_rpm / (wheel_speed_rad_s * 60 / 2π) / final_drive`
    pub fn ideal_ratio_for_wheel_speed(&self, wheel_speed_rad_s: f64) -> f64 {
        if wheel_speed_rad_s.abs() < 1e-6 {
            return self.ratio_max;
        }
        let wheel_rpm = wheel_speed_rad_s * 60.0 / (2.0 * std::f64::consts::PI);
        let ideal = self.target_engine_rpm / (wheel_rpm * self.final_drive);
        ideal.clamp(self.ratio_min, self.ratio_max)
    }
    /// Advance the CVT ratio toward the target using a first-order lag.
    ///
    /// # Arguments
    /// * `target_ratio` – desired ratio \[ratio_min, ratio_max\]
    /// * `dt`           – time step (s)
    ///
    /// Returns the new ratio.
    pub fn step_ratio(&mut self, target_ratio: f64, dt: f64) -> f64 {
        let target_clamped = target_ratio.clamp(self.ratio_min, self.ratio_max);
        let alpha = (dt / self.ratio_time_const).min(1.0);
        self.current_ratio += alpha * (target_clamped - self.current_ratio);
        self.current_ratio
    }
    /// Step the CVT to optimise engine RPM given current wheel speed.
    ///
    /// Combines `ideal_ratio_for_wheel_speed` and `step_ratio`.
    pub fn step_optimise(&mut self, wheel_speed_rad_s: f64, dt: f64) -> f64 {
        let ideal = self.ideal_ratio_for_wheel_speed(wheel_speed_rad_s);
        self.step_ratio(ideal, dt)
    }
    /// Clamped ratio (convenience accessor).
    pub fn ratio_clamped(&self) -> f64 {
        self.current_ratio.clamp(self.ratio_min, self.ratio_max)
    }
    /// Ratio coverage: ratio_max / ratio_min (dimensionless spread).
    pub fn ratio_coverage(&self) -> f64 {
        self.ratio_max / self.ratio_min
    }
}
/// State of a gear shift in progress.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShiftPhase {
    /// No shift in progress.
    Steady,
    /// Power-off overlap phase (torque is ramping down).
    Torque,
    /// Inertia phase (speed synchronisation).
    Inertia,
    /// Power-on phase (torque is restoring).
    PowerOn,
}
/// Torque-splitting differential with a continuously variable lock factor.
#[derive(Debug, Clone)]
pub struct Differential {
    /// Base torque split ratio (not currently used; reserved for future asymmetric diffs).
    pub ratio: f64,
    /// Lock factor: 0.0 = fully open (50/50), 1.0 = fully locked.
    pub lock_factor: f64,
}
impl Differential {
    /// Create a new differential.
    pub fn new(lock_factor: f64) -> Self {
        Self {
            ratio: 1.0,
            lock_factor: lock_factor.clamp(0.0, 1.0),
        }
    }
    /// Open differential (pure 50/50 split).
    pub fn open() -> Self {
        Self::new(0.0)
    }
    /// Locked differential (equal speeds, all torque transferred).
    pub fn locked() -> Self {
        Self::new(1.0)
    }
    /// Split input torque between left and right output shafts.
    ///
    /// Open diff: 50/50.  Locked diff: also 50/50 but with speed constraint
    /// (simplified model).  `lock_factor` interpolates between the two.
    ///
    /// Returns `(torque_left, torque_right)`.
    pub fn split_torque(&self, input_torque: f64, _wheel_speed_diff: f64) -> (f64, f64) {
        let half = input_torque * 0.5;
        (half, half)
    }
    /// Limited-slip split with a preload bias term.
    ///
    /// Adds `preload` (N·m) of locking torque to the slower wheel when slipping.
    /// Returns `(torque_left, torque_right)`.
    pub fn limited_slip(&self, input: f64, speed_diff: f64, preload: f64) -> (f64, f64) {
        let half = input * 0.5;
        let bias = (preload * self.lock_factor * speed_diff.signum()).min(half.abs());
        if speed_diff > 0.0 {
            (half - bias, half + bias)
        } else {
            (half + bias, half - bias)
        }
    }
}
impl Differential {
    /// Torsen-style (worm-gear) limited-slip differential torque split.
    ///
    /// The Torsen LSD biases torque to the wheel with more traction.  When a
    /// speed difference exists between the two output shafts the slower wheel
    /// (higher traction) receives a greater share governed by the TBR
    /// (Torque-Bias Ratio).
    ///
    /// `T_slow = T_in * TBR / (1 + TBR)`
    /// `T_fast = T_in * 1   / (1 + TBR)`
    ///
    /// where `TBR = 1 + lock_factor * (|Δω| / (|Δω| + omega_ref))`.
    ///
    /// # Arguments
    /// * `input_torque`  – total input torque to the diff (N·m)
    /// * `speed_left`    – left shaft angular velocity (rad/s)
    /// * `speed_right`   – right shaft angular velocity (rad/s)
    ///
    /// Returns `(torque_left, torque_right)`.
    pub fn compute_torque_split_limited_slip(
        &self,
        input_torque: f64,
        speed_left: f64,
        speed_right: f64,
    ) -> (f64, f64) {
        let omega_ref = 5.0_f64;
        let delta_omega = (speed_left - speed_right).abs();
        let tbr = 1.0 + self.lock_factor * delta_omega / (delta_omega + omega_ref);
        let t_slow = input_torque * tbr / (1.0 + tbr);
        let t_fast = input_torque / (1.0 + tbr);
        if speed_left <= speed_right {
            (t_slow, t_fast)
        } else {
            (t_fast, t_slow)
        }
    }
}
/// Centre differential distributing torque to front and rear axles.
#[derive(Debug, Clone)]
pub struct AwdCenterDifferential {
    /// Current operating mode.
    pub mode: AwdCenterMode,
    /// Software-commanded front torque fraction for Active mode \[0, 1\].
    pub commanded_front: f64,
    /// Front axle angular speed (rad/s).
    pub front_speed: f64,
    /// Rear axle angular speed (rad/s).
    pub rear_speed: f64,
}
impl AwdCenterDifferential {
    /// Create a new centre differential in Fixed mode with the given front bias.
    pub fn new_fixed(front_fraction: f64) -> Self {
        Self {
            mode: AwdCenterMode::Fixed(front_fraction.clamp(0.0, 1.0)),
            commanded_front: 0.4,
            front_speed: 0.0,
            rear_speed: 0.0,
        }
    }
    /// 40/60 front/rear fixed split (common road AWD default).
    pub fn fixed_40_60() -> Self {
        Self::new_fixed(0.4)
    }
    /// Torsen-style centre differential (viscous coupling, base 50/50).
    pub fn torsen_centre() -> Self {
        Self {
            mode: AwdCenterMode::Viscous {
                base_front: 0.5,
                viscosity: 10.0,
            },
            commanded_front: 0.5,
            front_speed: 0.0,
            rear_speed: 0.0,
        }
    }
    /// Compute torque split (front, rear) from total input torque.
    ///
    /// Returns `(front_torque, rear_torque)`.
    pub fn split_torque(&self, total_torque: f64) -> (f64, f64) {
        let front_frac = match self.mode {
            AwdCenterMode::Fixed(f) => f,
            AwdCenterMode::Viscous {
                base_front,
                viscosity,
            } => {
                let speed_diff = (self.front_speed - self.rear_speed).abs();
                let extra_rear = (viscosity * speed_diff / (speed_diff + 5.0)) * 0.3;
                if self.front_speed > self.rear_speed {
                    (base_front - extra_rear).clamp(0.1, 0.9)
                } else {
                    (base_front + extra_rear).clamp(0.1, 0.9)
                }
            }
            AwdCenterMode::Active => self.commanded_front.clamp(0.0, 1.0),
        };
        (total_torque * front_frac, total_torque * (1.0 - front_frac))
    }
}
/// Driveline NVH (Noise, Vibration, Harshness) model.
///
/// Characterises the torsional resonance frequency and bandwidth of the
/// driveline, and detects when an excitation frequency falls within the
/// resonant band.
#[derive(Debug, Clone)]
pub struct NvhDriveline {
    /// Driveline torsional natural frequency (Hz).
    pub natural_freq_hz: f64,
    /// Half-power bandwidth of the resonance peak (Hz).
    pub bandwidth_hz: f64,
}
impl NvhDriveline {
    /// Create a driveline NVH model.
    pub fn new(natural_freq_hz: f64, bandwidth_hz: f64) -> Self {
        Self {
            natural_freq_hz: natural_freq_hz.max(1.0),
            bandwidth_hz: bandwidth_hz.max(0.1),
        }
    }
    /// Default compact-car driveline (resonance at ~40 Hz, bandwidth 4 Hz).
    pub fn default_compact() -> Self {
        Self::new(40.0, 4.0)
    }
    /// Returns `true` when the excitation frequency falls within the resonant band.
    ///
    /// Band is `[f_n - bandwidth/2, f_n + bandwidth/2]`.
    pub fn is_resonant(&self, excitation_hz: f64, _damping_ratio: f64) -> bool {
        let half_bw = self.bandwidth_hz * 0.5;
        (excitation_hz - self.natural_freq_hz).abs() <= half_bw
    }
    /// Vibration amplitude (relative, dimensionless) at a given excitation frequency.
    ///
    /// Uses a simplified Lorentzian frequency response:
    /// `A(f) = (f_n²) / sqrt((f_n² - f²)² + (B * f)²)`
    /// where `B = bandwidth_hz`.
    pub fn vibration_amplitude(&self, excitation_hz: f64) -> f64 {
        let fn2 = self.natural_freq_hz * self.natural_freq_hz;
        let f2 = excitation_hz * excitation_hz;
        let b = self.bandwidth_hz;
        let denominator = ((fn2 - f2).powi(2) + (b * excitation_hz).powi(2)).sqrt();
        if denominator < 1e-12 {
            return f64::INFINITY;
        }
        fn2 / denominator
    }
    /// Driveline RPM corresponding to the resonant frequency for a given order.
    ///
    /// `rpm = 60 * f_n / order` — e.g. order 4 for a 4-cylinder engine.
    pub fn resonant_rpm(&self, order: f64) -> f64 {
        if order < 1e-6 {
            return f64::INFINITY;
        }
        60.0 * self.natural_freq_hz / order
    }
}
/// Hydraulic torque converter model.
///
/// Models the torque multiplication as a function of speed ratio (turbine/pump).
/// At stall (speed_ratio = 0) the torque ratio is at maximum.
/// At coupling point (speed_ratio → 1) the torque ratio approaches 1.
#[derive(Debug, Clone)]
pub struct TorqueConverter {
    /// Stall torque ratio (torque multiplication at speed_ratio = 0).
    pub stall_ratio: f64,
    /// Speed ratio at the coupling point (where torque ratio → 1.0).
    pub coupling_ratio: f64,
    /// Reference capacity factor at 1000 RPM (N·m·s²/rad²).
    pub capacity_factor_ref: f64,
}
impl TorqueConverter {
    /// Create a torque converter with given stall ratio and capacity factor.
    ///
    /// `stall_ratio` = torque multiplication at zero turbine speed (typically 1.8–2.5).
    /// `capacity_factor_ref` = K-factor at 1000 RPM (Pa·s² / rad²).
    pub fn new(stall_ratio: f64, capacity_factor_ref: f64) -> Self {
        Self {
            stall_ratio: stall_ratio.max(1.0),
            coupling_ratio: 0.85,
            capacity_factor_ref: capacity_factor_ref.max(1.0),
        }
    }
    /// Typical automotive torque converter (stall ratio ~2.0).
    pub fn typical_auto() -> Self {
        Self::new(2.0, 200.0)
    }
    /// Speed ratio: `sr = n_turbine / n_pump` in \[0, 1\].
    fn speed_ratio(pump_rpm: f64, turbine_rpm: f64) -> f64 {
        if pump_rpm.abs() < 1.0 {
            return 0.0;
        }
        (turbine_rpm / pump_rpm).clamp(0.0, 1.0)
    }
    /// Torque ratio at a given speed ratio (piecewise linear approximation).
    ///
    /// At `sr = 0.0` → `stall_ratio`.
    /// At `sr = coupling_ratio` → 1.0.
    /// At `sr > coupling_ratio` → 1.0 (lock-up).
    pub fn torque_ratio(&self, speed_ratio: f64) -> f64 {
        let sr = speed_ratio.clamp(0.0, 1.0);
        if sr >= self.coupling_ratio {
            return 1.0;
        }
        let frac = sr / self.coupling_ratio;
        self.stall_ratio + frac * (1.0 - self.stall_ratio)
    }
    /// Output torque delivered to the turbine shaft (N·m).
    ///
    /// `input_torque`: pump torque (N·m, from engine).
    /// `speed_ratio`: turbine/pump speed ratio.
    pub fn output_torque(&self, input_torque: f64, speed_ratio: f64) -> f64 {
        input_torque * self.torque_ratio(speed_ratio)
    }
    /// Capacity factor (K-factor) at a given pump RPM (N·m/(rpm)²).
    ///
    /// Represents the ability to absorb torque. Real K-factor varies with
    /// operating point; here approximated as constant for simplicity.
    pub fn capacity_factor(&self, pump_rpm: f64) -> f64 {
        if pump_rpm.abs() < 1.0 {
            return self.capacity_factor_ref;
        }
        self.capacity_factor_ref * (1000.0 / pump_rpm.abs()).sqrt()
    }
    /// Pump torque at given RPM using the capacity factor.
    ///
    /// `T_pump = K * n²`
    pub fn pump_torque(&self, pump_rpm: f64) -> f64 {
        let k = self.capacity_factor(pump_rpm);
        k * (pump_rpm / 1000.0).powi(2)
    }
    /// Whether the converter is in lock-up (speed ratio ≥ coupling point).
    pub fn is_locked_up(&self, pump_rpm: f64, turbine_rpm: f64) -> bool {
        Self::speed_ratio(pump_rpm, turbine_rpm) >= self.coupling_ratio
    }
}
/// Centre differential for an all-wheel-drive vehicle.
///
/// Distributes engine torque between front and rear axles. Supports:
/// - Open (fixed front/rear bias ratio)
/// - Viscous coupling (speed-sensitive)
/// - Active (direct torque demand from software)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AwdCenterMode {
    /// Fixed front-to-total torque split fraction \[0, 1\].
    Fixed(f64),
    /// Viscous coupling: bias increases with axle speed difference.
    Viscous {
        /// Base front fraction.
        base_front: f64,
        /// Coupling viscosity coefficient (N·m·s/rad).
        viscosity: f64,
    },
    /// Active: software sets front fraction directly.
    Active,
}
/// A simple planetary gear set (epicyclic gear train).
///
/// Three elements: sun (S), planet carrier (C), ring (R).
/// The fundamental relation (Willis equation):
///
/// `(ω_R - ω_C) / (ω_S - ω_C) = -Z_S / Z_R`
///
/// where Z_S and Z_R are the tooth counts of the sun and ring gears.
///
/// Torque relations (ideal, lossless):
/// `T_R = -T_S * Z_R / Z_S` (reaction torque on ring from sun input)
#[derive(Debug, Clone)]
pub struct PlanetaryGearSet {
    /// Number of teeth on the sun gear.
    pub z_sun: u32,
    /// Number of teeth on the ring gear.
    pub z_ring: u32,
    /// Mechanical efficiency (0..1).
    pub efficiency: f64,
}
impl PlanetaryGearSet {
    /// Create a new planetary gear set.
    ///
    /// Typical Z_sun = 30, Z_ring = 90 → basic gear ratio ≈ 4.0.
    pub fn new(z_sun: u32, z_ring: u32, efficiency: f64) -> Self {
        Self {
            z_sun: z_sun.max(1),
            z_ring: z_ring.max(z_sun + 1),
            efficiency: efficiency.clamp(0.0, 1.0),
        }
    }
    /// Typical Simpson compound planetary (used in 3-speed automatics).
    pub fn simpson_3speed() -> Self {
        Self::new(30, 90, 0.97)
    }
    /// Ravigneaux compound planetary (used in many modern automatics).
    pub fn ravigneaux() -> Self {
        Self::new(24, 72, 0.975)
    }
    /// Basic (standing) gear ratio when carrier is fixed: ω_R / ω_S.
    ///
    /// `i_0 = -Z_R / Z_S`  (negative because ring rotates opposite to sun
    /// when carrier is fixed)
    pub fn basic_ratio(&self) -> f64 {
        -(self.z_ring as f64) / (self.z_sun as f64)
    }
    /// Gear ratio when ring is the output, sun is the input, carrier fixed.
    ///
    /// `i = Z_R / Z_S`  (absolute value; sign conventions handled by caller)
    pub fn ratio_sun_to_ring_fixed_carrier(&self) -> f64 {
        (self.z_ring as f64) / (self.z_sun as f64)
    }
    /// Gear ratio when carrier is the output, sun is the input, ring fixed.
    ///
    /// `i = 1 + Z_S / Z_R`
    pub fn ratio_sun_to_carrier_fixed_ring(&self) -> f64 {
        1.0 + (self.z_sun as f64) / (self.z_ring as f64)
    }
    /// Ring speed given sun speed and carrier speed (Willis equation).
    ///
    /// `ω_R = (1 + i_0) * ω_C - i_0 * ω_S`
    /// where `i_0 = -Z_R / Z_S`.
    pub fn ring_speed(&self, sun_speed: f64, carrier_speed: f64) -> f64 {
        let i0 = self.basic_ratio();
        (1.0 + i0) * carrier_speed - i0 * sun_speed
    }
    /// Sun speed given ring speed and carrier speed (Willis equation).
    pub fn sun_speed(&self, ring_speed: f64, carrier_speed: f64) -> f64 {
        let i0 = self.basic_ratio();
        if i0.abs() < 1e-10 {
            return carrier_speed;
        }
        ((1.0 + i0) * carrier_speed - ring_speed) / i0
    }
    /// Output torque on the ring when sun receives input torque.
    /// Carrier is fixed (reaction member).
    ///
    /// `T_R = T_S * Z_R / Z_S * efficiency`
    pub fn ring_torque_from_sun(&self, sun_torque: f64) -> f64 {
        let ratio = (self.z_ring as f64) / (self.z_sun as f64);
        sun_torque * ratio * self.efficiency
    }
    /// Carrier output torque when sun is input, ring is fixed.
    ///
    /// `T_C = T_S * (1 + Z_R / Z_S) * efficiency`
    pub fn carrier_torque_from_sun_ring_fixed(&self, sun_torque: f64) -> f64 {
        let ratio = 1.0 + (self.z_ring as f64) / (self.z_sun as f64);
        sun_torque * ratio * self.efficiency
    }
    /// Number of planet gears needed to mesh correctly.
    ///
    /// Planets must satisfy: `(Z_sun + Z_ring) mod n_planets == 0`.
    /// Returns the minimum number of evenly-spaced planets (typically 3).
    pub fn min_planet_count(&self) -> u32 {
        for n in 2u32..=8 {
            if (self.z_sun + self.z_ring).is_multiple_of(n) {
                return n;
            }
        }
        3
    }
}
/// Lock-up state of a torque converter's bypass clutch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LockupState {
    /// Converter slip mode — torque multiplication active.
    Open,
    /// Clutch is engaging (transition).
    Engaging,
    /// Clutch is fully locked — direct mechanical drive.
    Locked,
}
/// Which wheels are driven by the engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DriveLayout {
    /// Front-wheel drive.
    FrontWheelDrive,
    /// Rear-wheel drive.
    #[default]
    RearWheelDrive,
    /// All-wheel drive with a center differential split.
    AllWheelDrive,
}
/// Engine simulation with RPM dynamics driven by throttle and load torque.
#[derive(Debug, Clone)]
pub struct Engine {
    /// Engine torque curve.
    pub curve: EngineCurve,
    /// Idle RPM (minimum operating speed).
    pub idle_rpm: f64,
    /// Maximum RPM (redline).
    pub max_rpm: f64,
    /// Rotational moment of inertia (kg·m²).
    pub inertia: f64,
    /// Throttle position \[0, 1\].
    pub throttle: f64,
    /// Current engine speed (RPM).
    pub rpm: f64,
}
impl Engine {
    /// Create a new engine from a torque curve and basic parameters.
    pub fn new(curve: EngineCurve, idle_rpm: f64, max_rpm: f64, inertia: f64) -> Self {
        Self {
            curve,
            idle_rpm,
            max_rpm,
            inertia,
            throttle: 0.0,
            rpm: idle_rpm,
        }
    }
    /// Typical compact car engine (4-cylinder, ~2.0 L).
    pub fn typical_4cylinder() -> Self {
        Self::new(EngineCurve::typical_4cylinder(), 800.0, 7000.0, 0.25)
    }
    /// Torque available at current RPM after applying throttle.
    pub fn available_torque(&self) -> f64 {
        self.throttle.clamp(0.0, 1.0) * self.curve.torque_at_rpm(self.rpm)
    }
    /// Advance engine state by `dt` seconds under `load_torque` (N·m, opposing).
    ///
    /// Uses Euler integration: dω/dt = (T_engine − T_load) / I
    pub fn update(&mut self, load_torque: f64, dt: f64) {
        let omega = self.rpm * 2.0 * std::f64::consts::PI / 60.0;
        let t_engine = self.available_torque();
        let alpha = (t_engine - load_torque) / self.inertia.max(1e-6);
        let new_omega = (omega + alpha * dt).max(0.0);
        let new_rpm = new_omega * 60.0 / (2.0 * std::f64::consts::PI);
        self.rpm = new_rpm.clamp(self.idle_rpm, self.max_rpm);
    }
    /// Returns `true` when the engine is at or above its rev limit.
    pub fn is_rev_limiting(&self) -> bool {
        self.rpm >= self.max_rpm
    }
    /// Approximate brake-specific fuel consumption (g/kWh) using a simple BSFC map.
    ///
    /// A value around 250–300 g/kWh is typical at peak-efficiency conditions.
    pub fn fuel_consumption_rate(&self) -> f64 {
        let rpm_norm = (self.rpm - self.idle_rpm) / (self.max_rpm - self.idle_rpm).max(1.0);
        let throttle = self.throttle.clamp(0.0, 1.0);
        250.0 + 80.0 * (1.0 - throttle).powi(2) + 40.0 * (rpm_norm - 0.5).powi(2)
    }
}
impl Engine {
    /// Estimate exhaust gas temperature (K) based on operating conditions.
    ///
    /// Uses a simplified energy balance:
    /// - At full throttle and high RPM the EGT approaches ~1100 K.
    /// - At idle it drops toward ~700 K.
    /// - Rich mixture (high throttle) increases temperature.
    /// - High RPM increases temperature due to higher flow.
    ///
    /// `EGT = T_base + k_throttle * throttle + k_rpm * rpm_norm`
    ///
    /// where `rpm_norm = (rpm - idle_rpm) / (max_rpm - idle_rpm)`.
    ///
    /// Returns temperature in Kelvin.
    pub fn compute_exhaust_temperature(&self) -> f64 {
        let t_base = 700.0_f64;
        let k_throttle = 300.0_f64;
        let k_rpm = 100.0_f64;
        let throttle = self.throttle.clamp(0.0, 1.0);
        let rpm_range = (self.max_rpm - self.idle_rpm).max(1.0);
        let rpm_norm = ((self.rpm - self.idle_rpm) / rpm_range).clamp(0.0, 1.0);
        t_base + k_throttle * throttle + k_rpm * rpm_norm
    }
}
/// Launch control system for controlled standing starts.
///
/// Limits engine RPM / throttle during the launch phase to prevent excessive
/// wheel spin while maximising traction force at the drive wheels.
#[derive(Debug, Clone)]
pub struct LaunchControl {
    /// Target launch RPM (engine speed to hold during slip buildup).
    pub target_rpm: f64,
    /// Maximum allowed slip ratio at the driven wheels.
    pub max_slip: f64,
    /// Throttle ramp rate once launch is initiated (fraction per second).
    pub throttle_ramp_rate: f64,
    /// Proportional gain for RPM regulation.
    pub kp_rpm: f64,
    /// Current internal throttle demand (after LC correction).
    pub throttle_output: f64,
}
impl LaunchControl {
    /// Create a launch control system with default racing car parameters.
    pub fn default_racing() -> Self {
        Self {
            target_rpm: 4500.0,
            max_slip: 0.12,
            throttle_ramp_rate: 2.0,
            kp_rpm: 0.0005,
            throttle_output: 0.0,
        }
    }
    /// Compute the launch control throttle output.
    ///
    /// In the launch phase the system holds engine RPM near `target_rpm`.
    /// Once the vehicle begins rolling (above `launch_speed_threshold`) the
    /// throttle ramps up at `throttle_ramp_rate` per second until it reaches
    /// the driver-requested level.
    ///
    /// # Arguments
    /// * `driver_throttle`  – driver's throttle pedal position \[0, 1\]
    /// * `engine_rpm`       – current engine speed (RPM)
    /// * `wheel_slip`       – maximum driven-wheel slip ratio
    /// * `vehicle_speed`    – vehicle longitudinal speed (m/s)
    /// * `dt`               – time step (s)
    ///
    /// Returns throttle fraction `[0, 1]` sent to the engine.
    pub fn compute_launch_control(
        &mut self,
        driver_throttle: f64,
        engine_rpm: f64,
        wheel_slip: f64,
        vehicle_speed: f64,
        dt: f64,
    ) -> f64 {
        let launch_speed_threshold = 2.0_f64;
        if vehicle_speed > launch_speed_threshold {
            self.throttle_output = (self.throttle_output + self.throttle_ramp_rate * dt)
                .min(driver_throttle)
                .min(1.0);
            return self.throttle_output;
        }
        let rpm_error = engine_rpm - self.target_rpm;
        let rpm_correction = (self.kp_rpm * rpm_error).clamp(-0.3, 0.3);
        let slip_penalty = if wheel_slip > self.max_slip {
            (wheel_slip - self.max_slip) * 2.0
        } else {
            0.0
        };
        let base_throttle = (driver_throttle - rpm_correction - slip_penalty).max(0.0);
        self.throttle_output = base_throttle.clamp(0.0, 1.0);
        self.throttle_output
    }
}
/// Legacy piecewise-linear engine curve using sorted `(rpm, torque)` tuples.
///
/// Prefer [`EngineCurve`] for new code.
#[derive(Debug, Clone)]
pub struct EngineCurveLegacy {
    /// Sorted list of (RPM, torque_Nm) sample points.
    pub(super) points: Vec<(Real, Real)>,
    /// Idle RPM (minimum operating RPM).
    pub idle_rpm: Real,
    /// Redline RPM (maximum operating RPM).
    pub redline_rpm: Real,
}
impl EngineCurveLegacy {
    /// Create a new engine curve from sample points.
    pub fn new(mut points: Vec<(Real, Real)>) -> Self {
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let idle_rpm = points.first().map(|p| p.0).unwrap_or(800.0);
        let redline_rpm = points.last().map(|p| p.0).unwrap_or(7000.0);
        Self {
            points,
            idle_rpm,
            redline_rpm,
        }
    }
    /// Create a simple flat-torque engine curve for testing.
    pub fn flat(torque: Real, idle_rpm: Real, redline_rpm: Real) -> Self {
        Self::new(vec![(idle_rpm, torque), (redline_rpm, torque)])
    }
    /// Look up torque at a given RPM using linear interpolation.
    pub fn torque_at_rpm(&self, rpm: Real) -> Real {
        if self.points.is_empty() {
            return 0.0;
        }
        if self.points.len() == 1 {
            return self.points[0].1;
        }
        let rpm_clamped = rpm.clamp(self.idle_rpm, self.redline_rpm);
        for i in 0..self.points.len() - 1 {
            let (r0, t0) = self.points[i];
            let (r1, t1) = self.points[i + 1];
            if rpm_clamped >= r0 && rpm_clamped <= r1 {
                if (r1 - r0).abs() < 1e-10 {
                    return t0;
                }
                let t = (rpm_clamped - r0) / (r1 - r0);
                return t0 + t * (t1 - t0);
            }
        }
        self.points
            .last()
            .expect("collection should not be empty")
            .1
    }
}
/// Legacy combined drivetrain.  Prefer [`Drivetrain`] for new code.
#[derive(Debug, Clone)]
pub struct DrivetrainLegacy {
    /// Engine torque curve.
    pub engine: EngineCurveLegacy,
    /// Gearbox with gear ratios.
    pub gearbox: GearboxLegacy,
    /// Front differential.
    pub front_diff: DifferentialMode,
    /// Rear differential.
    pub rear_diff: DifferentialMode,
    /// Drive layout (FWD, RWD, AWD).
    pub layout: DriveLayout,
    /// Front/rear torque split for AWD (0.0 = all rear, 1.0 = all front).
    pub awd_front_bias: Real,
    /// Maximum brake torque per wheel (N·m).
    pub max_brake_torque: Real,
    /// Current engine RPM.
    pub engine_rpm: Real,
}
impl DrivetrainLegacy {
    /// Create a new drivetrain with default components.
    pub fn new() -> Self {
        Self::default()
    }
    /// Compute individual wheel torques for a 4-wheel vehicle.
    pub fn compute_wheel_torques(
        &mut self,
        throttle: Real,
        brake: Real,
        wheel_speeds: &[Real; 4],
    ) -> [Real; 4] {
        let throttle = throttle.clamp(0.0, 1.0);
        let brake = brake.clamp(0.0, 1.0);
        let driven_speed = match self.layout {
            DriveLayout::FrontWheelDrive => (wheel_speeds[0].abs() + wheel_speeds[1].abs()) * 0.5,
            DriveLayout::RearWheelDrive => (wheel_speeds[2].abs() + wheel_speeds[3].abs()) * 0.5,
            DriveLayout::AllWheelDrive => {
                (wheel_speeds[0].abs()
                    + wheel_speeds[1].abs()
                    + wheel_speeds[2].abs()
                    + wheel_speeds[3].abs())
                    * 0.25
            }
        };
        let rpm = self.gearbox.wheel_speed_to_rpm(driven_speed);
        self.engine_rpm = rpm.clamp(self.engine.idle_rpm, self.engine.redline_rpm);
        let engine_torque = self.engine.torque_at_rpm(self.engine_rpm) * throttle;
        let ratio = self.gearbox.total_ratio();
        let driveshaft_torque = engine_torque * ratio * self.gearbox.efficiency;
        let (front_torque, rear_torque) = match self.layout {
            DriveLayout::FrontWheelDrive => (driveshaft_torque, 0.0),
            DriveLayout::RearWheelDrive => (0.0, driveshaft_torque),
            DriveLayout::AllWheelDrive => {
                let ft = driveshaft_torque * self.awd_front_bias;
                let rt = driveshaft_torque * (1.0 - self.awd_front_bias);
                (ft, rt)
            }
        };
        let (fl, fr) = self
            .front_diff
            .split_torque(front_torque, wheel_speeds[0], wheel_speeds[1]);
        let (rl, rr) = self
            .rear_diff
            .split_torque(rear_torque, wheel_speeds[2], wheel_speeds[3]);
        let brake_torque = brake * self.max_brake_torque;
        let apply_brake = |drive_torque: Real, wheel_speed: Real| -> Real {
            if brake > 0.0 {
                let braking = -brake_torque * wheel_speed.signum();
                if wheel_speed.abs() < 0.1 {
                    drive_torque - brake_torque * 0.5 * wheel_speed.signum()
                } else {
                    drive_torque + braking
                }
            } else {
                drive_torque
            }
        };
        [
            apply_brake(fl, wheel_speeds[0]),
            apply_brake(fr, wheel_speeds[1]),
            apply_brake(rl, wheel_speeds[2]),
            apply_brake(rr, wheel_speeds[3]),
        ]
    }
}
