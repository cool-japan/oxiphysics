//! Extended drivetrain types - Gearbox, Torsen, and LSD types.
//!
//! Split from types.rs to keep modules under 2000 lines.
use super::types::Gear;
use oxiphysics_core::math::Real;

/// Torsen (torque sensing) limited-slip differential.
///
/// Uses an internal worm-gear bias mechanism. The torque bias ratio (TBR)
/// determines how much more torque can go to the faster wheel before the
/// slower wheel starts to spin. Typical TBR: 2.5 -- 5.0.
#[derive(Debug, Clone)]
pub struct TorsenDifferential {
    /// Torque bias ratio: max torque to slow wheel / min torque to fast wheel.
    pub bias_ratio: f64,
}
impl TorsenDifferential {
    /// Create a Torsen differential with a given bias ratio.
    pub fn new(bias_ratio: f64) -> Self {
        Self {
            bias_ratio: bias_ratio.max(1.0),
        }
    }
    /// Typical Type-1 Torsen (bias ratio ~ 4.0).
    pub fn type1() -> Self {
        Self::new(4.0)
    }
    /// Compute torque split (left, right) from input torque and wheel speeds.
    ///
    /// When wheel speeds are equal, torque is split 50/50. As speed difference
    /// grows, more torque is directed to the slower wheel up to the bias ratio limit.
    pub fn torque_split(&self, input_torque: f64, left_rpm: f64, right_rpm: f64) -> (f64, f64) {
        let sum = left_rpm.abs() + right_rpm.abs();
        if sum < 1e-6 {
            return (input_torque * 0.5, input_torque * 0.5);
        }
        let speed_diff = (left_rpm - right_rpm).abs();
        let avg = sum * 0.5;
        let diff_ratio = (speed_diff / avg.max(1e-6)).min(1.0);
        let bias = 1.0 + diff_ratio * (self.bias_ratio - 1.0);
        let slow_frac = bias / (1.0 + bias);
        let fast_frac = 1.0 / (1.0 + bias);
        if left_rpm.abs() <= right_rpm.abs() {
            (input_torque * slow_frac, input_torque * fast_frac)
        } else {
            (input_torque * fast_frac, input_torque * slow_frac)
        }
    }
}
/// Multi-speed gearbox with a final drive and automatic-shift support.
#[derive(Debug, Clone)]
pub struct Gearbox {
    /// Forward gears (index 0 = 1st gear).
    pub gears: Vec<Gear>,
    /// Reverse gear.
    pub reverse: Gear,
    /// Currently selected gear: 0 = neutral, -1 = reverse, 1..N = forward.
    pub current_gear: i32,
    /// Final drive (axle) ratio.
    pub final_drive: f64,
    /// Time required to complete a gear change (seconds).
    pub shift_time: f64,
}
impl Gearbox {
    /// Create a new gearbox.
    pub fn new(gears: Vec<Gear>, reverse: Gear, final_drive: f64) -> Self {
        Self {
            gears,
            reverse,
            current_gear: 1,
            final_drive,
            shift_time: 0.3,
        }
    }
    /// Typical 6-speed manual gearbox with common sport-car ratios.
    pub fn typical_6speed() -> Self {
        let gears = vec![
            Gear::new(3.82, 0.97),
            Gear::new(2.20, 0.97),
            Gear::new(1.52, 0.97),
            Gear::new(1.15, 0.97),
            Gear::new(0.89, 0.97),
            Gear::new(0.73, 0.97),
        ];
        Self::new(gears, Gear::new(3.50, 0.95), 3.73)
    }
    /// Combined gear ratio for the currently selected gear x final drive.
    ///
    /// Returns 0.0 for neutral.
    pub fn total_ratio(&self) -> f64 {
        match self.current_gear {
            0 => 0.0,
            -1 => -self.reverse.ratio * self.final_drive,
            g if g >= 1 && (g as usize) <= self.gears.len() => {
                self.gears[(g - 1) as usize].ratio * self.final_drive
            }
            _ => 0.0,
        }
    }
    /// Gear efficiency for the current selection.
    fn current_efficiency(&self) -> f64 {
        match self.current_gear {
            0 => 1.0,
            -1 => self.reverse.efficiency,
            g if g >= 1 && (g as usize) <= self.gears.len() => {
                self.gears[(g - 1) as usize].efficiency
            }
            _ => 1.0,
        }
    }
    /// Output torque: T_out = T_in x ratio x efficiency.
    pub fn output_torque(&self, engine_torque: f64) -> f64 {
        engine_torque * self.total_ratio().abs() * self.current_efficiency()
    }
    /// Shift up one gear if not already at the top.
    pub fn shift_up(&mut self) {
        let max = self.gears.len() as i32;
        if self.current_gear < max {
            self.current_gear += 1;
        }
    }
    /// Shift down one gear.
    pub fn shift_down(&mut self) {
        if self.current_gear > -1 {
            self.current_gear -= 1;
        }
    }
    /// Simple automatic shift logic based on engine RPM and throttle.
    ///
    /// Upshifts at high RPM (>85% of redline) and downshifts at low RPM (<25%).
    pub fn auto_shift(&mut self, engine_rpm: f64, _throttle: f64) {
        let upshift_rpm = 6000.0;
        let downshift_rpm = 1500.0;
        if engine_rpm > upshift_rpm {
            self.shift_up();
        } else if engine_rpm < downshift_rpm && self.current_gear > 1 {
            self.shift_down();
        }
    }
}
/// Legacy multi-speed gearbox.  Prefer [`Gearbox`] for new code.
#[derive(Debug, Clone)]
pub struct GearboxLegacy {
    /// Forward gear ratios (index 0 = 1st gear).
    pub gear_ratios: Vec<Real>,
    /// Reverse gear ratio (positive value; direction handled externally).
    pub reverse_ratio: Real,
    /// Final drive (differential) ratio.
    pub final_drive_ratio: Real,
    /// Currently selected gear (0 = neutral, 1..N = forward, -1 = reverse).
    pub current_gear: i32,
    /// Transmission efficiency (0..1).
    pub efficiency: Real,
}
impl GearboxLegacy {
    /// Create a new gearbox.
    pub fn new(gear_ratios: Vec<Real>, reverse_ratio: Real, final_drive_ratio: Real) -> Self {
        Self {
            gear_ratios,
            reverse_ratio,
            final_drive_ratio,
            current_gear: 1,
            efficiency: 0.9,
        }
    }
    /// Number of forward gears.
    pub fn num_gears(&self) -> usize {
        self.gear_ratios.len()
    }
    /// Shift to the given gear.  Clamps to valid range.
    pub fn shift(&mut self, gear: i32) {
        let max_gear = self.gear_ratios.len() as i32;
        self.current_gear = gear.clamp(-1, max_gear);
    }
    /// Shift up one gear.
    pub fn shift_up(&mut self) {
        let max_gear = self.gear_ratios.len() as i32;
        if self.current_gear < max_gear {
            self.current_gear += 1;
        }
    }
    /// Shift down one gear.
    pub fn shift_down(&mut self) {
        if self.current_gear > -1 {
            self.current_gear -= 1;
        }
    }
    /// Get the total gear ratio for the current gear.  Returns 0.0 for neutral.
    pub fn total_ratio(&self) -> Real {
        match self.current_gear {
            0 => 0.0,
            -1 => -self.reverse_ratio * self.final_drive_ratio,
            g if g >= 1 && (g as usize) <= self.gear_ratios.len() => {
                self.gear_ratios[(g - 1) as usize] * self.final_drive_ratio
            }
            _ => 0.0,
        }
    }
    /// Compute engine RPM from wheel angular velocity and current gear.
    pub fn wheel_speed_to_rpm(&self, wheel_angular_velocity: Real) -> Real {
        let ratio = self.total_ratio();
        if ratio.abs() < 1e-10 {
            return 0.0;
        }
        (wheel_angular_velocity * ratio.abs() * 60.0) / (2.0 * std::f64::consts::PI)
    }
}
/// Clutch-pack limited-slip differential.
#[derive(Debug, Clone)]
pub struct LimitedSlipDifferential {
    /// Clutch preload (N-m).
    pub preload_torque: f64,
    /// Ramp angle of the clutch-pack actuating ramps (degrees).
    pub ramp_angle_deg: f64,
    /// Current angular velocity of the left output shaft (rad/s).
    pub left_wheel_speed: f64,
    /// Current angular velocity of the right output shaft (rad/s).
    pub right_wheel_speed: f64,
}
impl LimitedSlipDifferential {
    /// Create a new LSD with given preload and a default 45 deg ramp angle.
    pub fn new(preload: f64) -> Self {
        Self {
            preload_torque: preload,
            ramp_angle_deg: 45.0,
            left_wheel_speed: 0.0,
            right_wheel_speed: 0.0,
        }
    }
    /// Compute torque split with LSD behaviour.
    pub fn torque_split(&self, input_torque: f64) -> (f64, f64) {
        let ramp_rad = self.ramp_angle_deg.to_radians();
        let half_torque = input_torque.abs() * 0.5 + 1e-6;
        let bias = 1.0 + self.preload_torque / half_torque + ramp_rad.tan();
        let fast_frac = 1.0 / (1.0 + bias);
        let slow_frac = bias / (1.0 + bias);
        if self.left_wheel_speed.abs() <= self.right_wheel_speed.abs() {
            (input_torque * slow_frac, input_torque * fast_frac)
        } else {
            (input_torque * fast_frac, input_torque * slow_frac)
        }
    }
}
