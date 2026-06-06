// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::{Real, Vec3};

/// Tire acoustic resonance (cavity resonance) model.
///
/// The air cavity inside a pneumatic tire resonates at a frequency determined
/// by the cavity geometry.  This is a major source of road noise ("cavity boom").
///
/// f_cavity ≈ c_sound / (π * D_rim)
#[derive(Debug, Clone)]
pub struct TireAcousticResonance {
    /// Speed of sound in air (m/s). Default ≈ 343 m/s.
    pub speed_of_sound: f64,
    /// Rim diameter (m).
    pub rim_diameter: f64,
    /// Quality factor Q of the cavity resonance.
    pub quality_factor: f64,
    /// Excitation amplitude at current resonance.
    pub current_amplitude: f64,
}
impl TireAcousticResonance {
    /// Create a cavity resonance model for a given rim.
    pub fn new(rim_diameter: f64, quality_factor: f64) -> Self {
        Self {
            speed_of_sound: 343.0,
            rim_diameter,
            quality_factor,
            current_amplitude: 0.0,
        }
    }
    /// Cavity resonance frequency (Hz).
    ///
    /// Uses the half-wavelength formula: f = c / (π * D).
    pub fn cavity_frequency(&self) -> f64 {
        self.speed_of_sound / (std::f64::consts::PI * self.rim_diameter.max(1e-6))
    }
    /// Bandwidth of the resonance (Hz): BW = f_res / Q.
    pub fn bandwidth(&self) -> f64 {
        self.cavity_frequency() / self.quality_factor.max(1e-6)
    }
    /// Amplitude response at a given excitation frequency (dimensionless).
    ///
    /// Lorentzian resonance shape:
    /// A(f) = 1 / sqrt(1 + Q^2 * ((f/f0) - (f0/f))^2)
    pub fn amplitude_at(&self, freq: f64) -> f64 {
        let f0 = self.cavity_frequency();
        if f0 < 1e-6 || freq < 1e-6 {
            return 0.0;
        }
        let x = self.quality_factor * (freq / f0 - f0 / freq);
        1.0 / (1.0 + x * x).sqrt()
    }
    /// Excitation frequency from vehicle speed and tire circumference.
    pub fn excitation_frequency(&self, speed: f64, circumference: f64) -> f64 {
        if circumference < 1e-6 {
            return 0.0;
        }
        speed / circumference
    }
}
/// Simplified rim-bead interaction: checks if tire is seated correctly
/// and estimates bead unseating risk.
#[derive(Debug, Clone)]
pub struct RimBeadModel {
    /// Bead seating force required to keep tire on rim (N).
    pub bead_seating_force: f64,
    /// Rim diameter (m).
    pub rim_diameter: f64,
    /// Minimum inflation pressure to maintain bead seating (Pa).
    pub min_inflation_pressure: f64,
}
impl RimBeadModel {
    /// Create a new rim-bead model.
    pub fn new(bead_seating_force: f64, rim_diameter: f64, min_inflation_pressure: f64) -> Self {
        Self {
            bead_seating_force,
            rim_diameter,
            min_inflation_pressure,
        }
    }
    /// Check if the bead remains seated at the given inflation pressure.
    pub fn is_seated(&self, inflation_pressure: f64) -> bool {
        inflation_pressure >= self.min_inflation_pressure
    }
    /// Bead retention margin: ratio of inflation pressure to minimum.
    /// Values > 1.0 = safe; < 1.0 = risk of unseating.
    pub fn retention_margin(&self, inflation_pressure: f64) -> f64 {
        inflation_pressure / self.min_inflation_pressure.max(1.0)
    }
}
/// Simplified brake disc temperature model.
#[derive(Debug, Clone)]
pub struct BrakeDiscTemperature {
    /// Current disc temperature (C).
    pub temperature: f64,
    /// Disc mass (kg).
    pub mass: f64,
    /// Specific heat capacity (J/(kg*K)).
    pub specific_heat: f64,
    /// Cooling coefficient (W/K) — models convective cooling.
    pub cooling_coeff: f64,
    /// Ambient temperature (C).
    pub ambient: f64,
}
impl BrakeDiscTemperature {
    /// Create a new brake disc temperature model.
    pub fn new(mass: f64, specific_heat: f64, cooling_coeff: f64, ambient: f64) -> Self {
        Self {
            temperature: ambient,
            mass,
            specific_heat,
            cooling_coeff,
            ambient,
        }
    }
    /// Step the temperature model.
    ///
    /// * `brake_power` — power dissipated by braking (W).
    /// * `dt` — time step (s).
    pub fn step(&mut self, brake_power: f64, dt: f64) {
        let heat_in = brake_power * dt;
        let heat_out = self.cooling_coeff * (self.temperature - self.ambient) * dt;
        let net_heat = heat_in - heat_out;
        let delta_t = net_heat / (self.mass * self.specific_heat).max(1e-10);
        self.temperature += delta_t;
    }
    /// Brake fade factor: friction coefficient drops at high temperatures.
    ///
    /// Returns a multiplier in \[0.3, 1.0\].
    pub fn fade_factor(&self) -> f64 {
        let fade_start = 300.0;
        let fade_end = 700.0;
        if self.temperature <= fade_start {
            1.0
        } else if self.temperature >= fade_end {
            0.3
        } else {
            let t = (self.temperature - fade_start) / (fade_end - fade_start);
            1.0 - 0.7 * t
        }
    }
}
/// Degressive load sensitivity model for tire grip.
///
/// Tires show degressive behavior: doubling the load does not double the
/// maximum lateral force. The grip coefficient falls as load increases.
///
/// mu_peak(Fz) = mu0 - k_deg * (Fz - Fz_nom) / Fz_nom
#[derive(Debug, Clone)]
pub struct LoadSensitivity {
    /// Peak friction at nominal load.
    pub mu0: f64,
    /// Degressive slope coefficient (dimensionless).
    pub k_deg: f64,
    /// Nominal vertical load (N).
    pub fz_nominal: f64,
}
impl LoadSensitivity {
    /// Create a load sensitivity model.
    pub fn new(mu0: f64, k_deg: f64, fz_nominal: f64) -> Self {
        Self {
            mu0,
            k_deg,
            fz_nominal,
        }
    }
    /// Typical road tire: mu0=1.0, k_deg=0.15, Fz_nom=4000 N.
    pub fn road_tire() -> Self {
        Self {
            mu0: 1.0,
            k_deg: 0.15,
            fz_nominal: 4000.0,
        }
    }
    /// Racing slick: high mu0, lower degression.
    pub fn racing_slick() -> Self {
        Self {
            mu0: 1.6,
            k_deg: 0.10,
            fz_nominal: 5000.0,
        }
    }
    /// Peak friction coefficient at the given vertical load.
    pub fn mu_at(&self, fz: f64) -> f64 {
        let rel = (fz - self.fz_nominal) / self.fz_nominal.max(1.0);
        (self.mu0 - self.k_deg * rel).max(0.1)
    }
    /// Maximum lateral force (N) at the given vertical load.
    pub fn max_lateral_force(&self, fz: f64) -> f64 {
        self.mu_at(fz) * fz.max(0.0)
    }
}
/// Wheel alignment parameters (camber, caster, toe).
#[derive(Debug, Clone)]
pub struct WheelAlignment {
    /// Camber angle (radians). Negative = top of wheel tilted inward.
    pub camber: f64,
    /// Caster angle (radians). Positive = steering axis tilted rearward.
    pub caster: f64,
    /// Toe angle (radians). Positive = toe-in (front of wheel points inward).
    pub toe: f64,
}
impl WheelAlignment {
    /// Create a new alignment spec.
    pub fn new(camber: f64, caster: f64, toe: f64) -> Self {
        Self {
            camber,
            caster,
            toe,
        }
    }
    /// Zero alignment (all angles 0).
    pub fn zero() -> Self {
        Self {
            camber: 0.0,
            caster: 0.0,
            toe: 0.0,
        }
    }
    /// Typical sports car front alignment.
    pub fn sport_front() -> Self {
        Self {
            camber: (-1.5_f64).to_radians(),
            caster: 6.0_f64.to_radians(),
            toe: 0.1_f64.to_radians(),
        }
    }
    /// Effective camber change due to body roll.
    ///
    /// `roll_angle` is the chassis roll angle (radians).
    /// Returns the adjusted camber angle.
    pub fn effective_camber(&self, roll_angle: f64) -> f64 {
        self.camber - roll_angle
    }
    /// Mechanical trail from caster angle and kingpin offset.
    ///
    /// `kingpin_offset` is the lateral offset of the steering axis (m).
    pub fn mechanical_trail(&self, kingpin_offset: f64) -> f64 {
        kingpin_offset * self.caster.sin()
    }
    /// Total toe for an axle (sum of left and right toe).
    pub fn total_toe(&self, other: &WheelAlignment) -> f64 {
        self.toe + other.toe
    }
}
/// Records a rolling window of wheel vertical loads for analysis.
#[derive(Debug, Clone)]
pub struct WheelLoadHistory {
    pub(super) loads: Vec<f64>,
    pub(super) times: Vec<f64>,
    pub(super) max_samples: usize,
}
impl WheelLoadHistory {
    /// Create a new load history with the given window size.
    pub fn new(max_samples: usize) -> Self {
        Self {
            loads: Vec::new(),
            times: Vec::new(),
            max_samples: max_samples.max(2),
        }
    }
    /// Record a load sample at time `t` (s).
    pub fn record(&mut self, t: f64, load: f64) {
        if self.loads.len() >= self.max_samples {
            self.loads.remove(0);
            self.times.remove(0);
        }
        self.loads.push(load.max(0.0));
        self.times.push(t);
    }
    /// Mean vertical load (N) over the window.
    pub fn mean_load(&self) -> f64 {
        if self.loads.is_empty() {
            return 0.0;
        }
        self.loads.iter().sum::<f64>() / self.loads.len() as f64
    }
    /// Peak vertical load (N) in the window.
    pub fn peak_load(&self) -> f64 {
        self.loads.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Load variation coefficient: std_dev / mean.  Returns 0 if mean ≈ 0.
    pub fn variation_coefficient(&self) -> f64 {
        if self.loads.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_load();
        if mean < 1e-6 {
            return 0.0;
        }
        let variance =
            self.loads.iter().map(|&l| (l - mean).powi(2)).sum::<f64>() / self.loads.len() as f64;
        variance.sqrt() / mean
    }
    /// Rate of change of load at the most recent sample (N/s).
    pub fn load_rate(&self) -> f64 {
        let n = self.loads.len();
        if n < 2 {
            return 0.0;
        }
        let dt = self.times[n - 1] - self.times[n - 2];
        if dt.abs() < 1e-12 {
            return 0.0;
        }
        (self.loads[n - 1] - self.loads[n - 2]) / dt
    }
}
/// Simple raycast wheel that tests against an infinite flat plane (y = 0).
#[derive(Debug, Clone)]
pub struct RaycastWheel {
    /// World-space origin of the ray (suspension top).
    pub position: [f64; 3],
    /// Normalised ray direction (typically downward = \[0,-1,0\]).
    pub direction: [f64; 3],
    /// Maximum ray length (suspension rest + travel, metres).
    pub max_length: f64,
    /// Wheel radius (metres).
    pub radius: f64,
}
impl RaycastWheel {
    /// Create a new raycast wheel.
    pub fn new(position: [f64; 3], direction: [f64; 3], max_length: f64, radius: f64) -> Self {
        Self {
            position,
            direction,
            max_length,
            radius,
        }
    }
    /// Cast the ray against a flat plane at height `plane_y` (default 0.0).
    ///
    /// In a real engine `bodies` would be the physics world; here we mock with
    /// a simple infinite plane test.
    pub fn cast(&self, plane_y: f64) -> WheelContact {
        let target_y = plane_y + self.radius;
        let dir_y = self.direction[1];
        if dir_y.abs() < 1e-10 {
            return WheelContact::default();
        }
        let t = (target_y - self.position[1]) / dir_y;
        if t < 0.0 || t > self.max_length {
            return WheelContact::default();
        }
        let contact_pos = [
            self.position[0] + t * self.direction[0],
            self.position[1] + t * self.direction[1],
            self.position[2] + t * self.direction[2],
        ];
        let penetration = (self.max_length - t).max(0.0);
        WheelContact {
            hit: true,
            position: contact_pos,
            normal: [0.0, 1.0, 0.0],
            penetration,
            body_idx: None,
        }
    }
}
/// Pneumatic trail: the distance behind the contact centre where the resultant
/// lateral force acts, creating a self-aligning torque.
///
/// As slip angle increases, pneumatic trail decreases (tire saturates).
#[derive(Debug, Clone)]
pub struct PneumaticTrail {
    /// Maximum pneumatic trail at small slip angles (m).
    pub t_max: f64,
    /// Characteristic slip angle at which trail halves (rad).
    pub alpha_half: f64,
}
impl PneumaticTrail {
    /// Create a pneumatic trail model.
    pub fn new(t_max: f64, alpha_half: f64) -> Self {
        Self { t_max, alpha_half }
    }
    /// Trail length at the given slip angle magnitude (m).
    ///
    /// Uses a Lorentzian decay: t = t_max / (1 + (alpha/alpha_half)^2)
    pub fn trail_at(&self, alpha: f64) -> f64 {
        let r = alpha.abs() / self.alpha_half.max(1e-12);
        self.t_max / (1.0 + r * r)
    }
    /// Self-aligning torque (N·m) given lateral force (N) and slip angle.
    pub fn self_aligning_torque(&self, lateral_force: f64, alpha: f64) -> f64 {
        self.trail_at(alpha) * lateral_force
    }
}
/// Road roughness power spectral density (PSD) model per ISO 8608.
///
/// The road PSD is modeled as: G_d(n) = G_d(n_0) * (n / n_0)^(-w)
/// where n is spatial frequency (cycles/m) and w is the waviness exponent.
#[derive(Debug, Clone)]
pub struct RoadTexturePSD {
    /// Reference PSD at n_0 = 0.1 cycles/m (m³/cycle).
    pub g_d_ref: f64,
    /// Waviness exponent (typically 2..3).
    pub waviness: f64,
    /// Reference spatial frequency (cycles/m).
    pub n_ref: f64,
}
impl RoadTexturePSD {
    /// ISO 8608 Class A road (very good).
    pub fn class_a() -> Self {
        Self {
            g_d_ref: 1.0e-6,
            waviness: 2.0,
            n_ref: 0.1,
        }
    }
    /// ISO 8608 Class B road (good).
    pub fn class_b() -> Self {
        Self {
            g_d_ref: 4.0e-6,
            waviness: 2.0,
            n_ref: 0.1,
        }
    }
    /// ISO 8608 Class C road (average).
    pub fn class_c() -> Self {
        Self {
            g_d_ref: 16.0e-6,
            waviness: 2.0,
            n_ref: 0.1,
        }
    }
    /// ISO 8608 Class D road (poor).
    pub fn class_d() -> Self {
        Self {
            g_d_ref: 64.0e-6,
            waviness: 2.0,
            n_ref: 0.1,
        }
    }
    /// PSD at spatial frequency `n` (cycles/m).
    pub fn psd_at(&self, n: f64) -> f64 {
        if n < 1e-12 {
            return 0.0;
        }
        self.g_d_ref * (n / self.n_ref).powf(-self.waviness)
    }
    /// RMS roughness over a spatial frequency band \[n1, n2\].
    ///
    /// Uses trapezoidal integration with `steps` intervals.
    pub fn rms_roughness(&self, n1: f64, n2: f64, steps: usize) -> f64 {
        if n2 <= n1 || steps == 0 {
            return 0.0;
        }
        let dn = (n2 - n1) / steps as f64;
        let mut integral = 0.0;
        for k in 0..=steps {
            let n = n1 + k as f64 * dn;
            let w = if k == 0 || k == steps { 0.5 } else { 1.0 };
            integral += w * self.psd_at(n) * dn;
        }
        integral.sqrt()
    }
    /// Temporal PSD at frequency `f` (Hz) for vehicle speed `v` (m/s).
    ///
    /// Converts spatial PSD to temporal: G_d(f) = G_d(n) / v, where n = f / v.
    pub fn temporal_psd(&self, freq_hz: f64, speed: f64) -> f64 {
        if speed < 1e-6 {
            return 0.0;
        }
        let n = freq_hz / speed;
        self.psd_at(n) / speed
    }
}
/// Runtime state of a wheel during simulation.
#[derive(Debug, Clone)]
pub struct WheelState {
    /// Current suspension compression length in meters.
    pub suspension_length: Real,
    /// Suspension force magnitude (positive = pushing wheel down).
    pub suspension_force: Real,
    /// World-space contact point (valid only if `is_in_contact`).
    pub contact_point: Vec3,
    /// World-space contact surface normal (valid only if `is_in_contact`).
    pub contact_normal: Vec3,
    /// Whether the wheel is currently in contact with a surface.
    pub is_in_contact: bool,
    /// Current steering angle in radians.
    pub steering_angle: Real,
    /// Cumulative rotation angle of the wheel (visual only).
    pub rotation_angle: Real,
    /// Angular velocity of the wheel about its axle (rad/s).
    pub spin_velocity: Real,
    /// Lateral slip angle in radians.
    pub slip_angle: Real,
    /// Longitudinal slip ratio (dimensionless).
    pub slip_ratio: Real,
    /// Lateral tire force in Newtons.
    pub lateral_force: Real,
    /// Longitudinal tire force in Newtons.
    pub longitudinal_force: Real,
    /// World-space position of the wheel center (for rendering).
    pub world_position: Vec3,
    /// Forward direction of the wheel in world space.
    pub forward_direction: Vec3,
    /// Right (axle) direction of the wheel in world space.
    pub right_direction: Vec3,
    /// Current angular velocity (rad/s) used for spin dynamics.
    pub angular_velocity: Real,
    /// Current load (normal force) on the wheel (N).
    pub load: Real,
    /// Current tire temperature (°C) — affects grip.
    pub temperature: Real,
}
impl WheelState {
    /// Update the wheel's angular velocity and spin given applied torques and contact.
    ///
    /// * `torque`         — drive/brake torque applied to the wheel (N·m). Positive = accelerate.
    /// * `brake_torque`   — braking torque magnitude (always opposing spin, N·m).
    /// * `contact_force`  — normal contact force (N).
    /// * `dt`             — time step (s).
    pub fn update(&mut self, torque: Real, brake_torque: Real, contact_force: Real, dt: Real) {
        let radius = 0.3_f64;
        let mass = 15.0_f64;
        let inertia = 0.5 * mass * radius * radius;
        let brake_sign = if self.angular_velocity.abs() < 1e-6 {
            0.0
        } else {
            -self.angular_velocity.signum()
        };
        let effective_brake = brake_sign * brake_torque.abs();
        let rolling_resistance = -0.01 * contact_force * self.angular_velocity.signum();
        let net_torque = torque + effective_brake + rolling_resistance;
        let angular_accel = net_torque / inertia;
        self.angular_velocity += angular_accel * dt;
        self.spin_velocity = self.angular_velocity;
        self.rotation_angle += self.angular_velocity * dt;
        self.load = contact_force;
    }
}
/// Rolling resistance force model.
///
/// Models the resistive force at the contact patch due to tire deformation.
/// The effective coefficient is speed-dependent per ISO 8767.
#[derive(Debug, Clone)]
pub struct RollingResistance {
    /// Base rolling resistance coefficient (dimensionless).
    pub c_r0: f64,
    /// Speed-dependent coefficient (1/(m/s)²).
    pub c_r1: f64,
}
impl RollingResistance {
    /// Create a rolling resistance model.
    ///
    /// Typical values for a passenger car: c_r0 = 0.01, c_r1 = 6e-7.
    pub fn new(c_r0: f64, c_r1: f64) -> Self {
        Self { c_r0, c_r1 }
    }
    /// Standard road tire (ISO 8767 default).
    pub fn standard_road() -> Self {
        Self {
            c_r0: 0.013,
            c_r1: 6.5e-7,
        }
    }
    /// Low-rolling-resistance tire (EV-optimised).
    pub fn low_resistance() -> Self {
        Self {
            c_r0: 0.007,
            c_r1: 4.0e-7,
        }
    }
    /// Racing slick tire (low c_r due to compound hardness).
    pub fn racing_slick() -> Self {
        Self {
            c_r0: 0.018,
            c_r1: 8.0e-7,
        }
    }
    /// Effective rolling resistance coefficient at the given speed (m/s).
    ///
    /// c_r(v) = c_r0 + c_r1 * v²
    pub fn coefficient_at_speed(&self, speed: f64) -> f64 {
        self.c_r0 + self.c_r1 * speed * speed
    }
    /// Rolling resistance force magnitude (N) given normal load and speed.
    pub fn force(&self, normal_load: f64, speed: f64) -> f64 {
        self.coefficient_at_speed(speed) * normal_load.max(0.0)
    }
    /// Power dissipated by rolling resistance at given speed (W).
    pub fn power(&self, normal_load: f64, speed: f64) -> f64 {
        self.force(normal_load, speed) * speed.abs()
    }
}
/// Combined tire and wheel assembly.
#[derive(Debug, Clone)]
pub struct TireWheelAssembly {
    /// Wheel configuration.
    pub wheel: Wheel,
    /// Tire wear model.
    pub wear: TireWear,
    /// Wheel alignment.
    pub alignment: WheelAlignment,
    /// Brake disc temperature.
    pub brake_temp: BrakeDiscTemperature,
}
impl TireWheelAssembly {
    /// Create a new assembly with default components.
    pub fn new(wheel: Wheel, compound: TireCompound) -> Self {
        let wear_rate = match compound {
            TireCompound::Soft => 2e-5,
            TireCompound::Medium => 1e-5,
            TireCompound::Hard => 5e-6,
        };
        Self {
            wheel,
            wear: TireWear::new(compound, wear_rate),
            alignment: WheelAlignment::zero(),
            brake_temp: BrakeDiscTemperature::new(5.0, 500.0, 10.0, 25.0),
        }
    }
    /// Step the assembly: update wear and brake temperature.
    pub fn step(&mut self, slip: f64, load: f64, brake_power: f64, dt: f64) {
        self.wear.update(slip, load, dt);
        self.brake_temp.step(brake_power, dt);
    }
    /// Current grip factor accounting for wear and temperature.
    pub fn grip_factor(&self) -> f64 {
        self.wear.grip_factor() * self.brake_temp.fade_factor()
    }
}
/// Gyroscopic precession moment for a spinning wheel.
///
/// When the axle orientation changes (steering or chassis pitch/roll), the
/// gyroscopic effect produces a moment perpendicular to both the spin axis
/// and the precession axis.
///
/// M = I × ω × Ω_prec
#[derive(Debug, Clone)]
pub struct GyroscopicPrecession {
    /// Spin moment of inertia about the axle (kg·m²).
    pub spin_inertia: f64,
    /// Current spin rate (rad/s).
    pub spin_rate: f64,
}
impl GyroscopicPrecession {
    /// Create a new gyroscopic precession model.
    pub fn new(spin_inertia: f64) -> Self {
        Self {
            spin_inertia,
            spin_rate: 0.0,
        }
    }
    /// Set the spin rate from wheel angular velocity.
    pub fn set_spin_rate(&mut self, omega: f64) {
        self.spin_rate = omega;
    }
    /// Angular momentum magnitude (kg·m²/s).
    pub fn angular_momentum(&self) -> f64 {
        self.spin_inertia * self.spin_rate.abs()
    }
    /// Gyroscopic moment (N·m) given precession rate about a perpendicular axis.
    pub fn gyroscopic_moment(&self, precession_rate: f64) -> f64 {
        self.spin_inertia * self.spin_rate * precession_rate
    }
    /// Gyroscopic moment vector \[mx, my, mz\] given spin axis and precession axis.
    ///
    /// `spin_axis` — unit vector along wheel axle (e.g. \[1,0,0\]).
    /// `prec_axis` — precession axis (e.g. yaw = \[0,1,0\]).
    /// Returns the gyroscopic moment vector M = I * omega * (spin_axis × prec_axis).
    pub fn moment_vector(&self, spin_axis: [f64; 3], prec_axis: [f64; 3]) -> [f64; 3] {
        let cross = [
            spin_axis[1] * prec_axis[2] - spin_axis[2] * prec_axis[1],
            spin_axis[2] * prec_axis[0] - spin_axis[0] * prec_axis[2],
            spin_axis[0] * prec_axis[1] - spin_axis[1] * prec_axis[0],
        ];
        let scale = self.spin_inertia * self.spin_rate;
        [cross[0] * scale, cross[1] * scale, cross[2] * scale]
    }
}
/// Configuration parameters for a wheel.
#[derive(Debug, Clone)]
pub struct Wheel {
    /// Wheel radius in meters.
    pub radius: Real,
    /// Wheel width in meters.
    pub width: Real,
    /// Wheel mass in kilograms.
    pub mass: Real,
    /// Rest length of the suspension spring in meters.
    pub suspension_rest_length: Real,
    /// Suspension spring stiffness (N/m).
    pub suspension_stiffness: Real,
    /// Damping coefficient for compression (N*s/m).
    pub damping_compression: Real,
    /// Damping coefficient for relaxation (N*s/m).
    pub damping_relaxation: Real,
    /// Maximum suspension travel in meters.
    pub max_suspension_travel: Real,
    /// Friction slip coefficient (peak friction multiplier).
    pub friction_slip: Real,
    /// Local-space position of the suspension mount point relative to chassis.
    pub connection_point: Vec3,
    /// Suspension direction in local space (typically downward).
    pub suspension_direction: Vec3,
    /// Axle direction in local space (typically right for left wheel).
    pub axle_direction: Vec3,
}
/// Tire compound types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TireCompound {
    /// Soft compound — high grip, fast wear.
    Soft,
    /// Medium compound — balanced.
    Medium,
    /// Hard compound — low wear, slightly lower peak grip.
    Hard,
}
/// 4-wheel layout configuration.
///
/// Corner indices: 0 = front-left, 1 = front-right, 2 = rear-left, 3 = rear-right.
#[derive(Debug, Clone)]
pub struct WheelSetConfig {
    /// World-space position of each corner's wheel centre at rest.
    pub positions: [[f64; 3]; 4],
    /// Track width (left–right distance between wheel centres, m).
    pub track_width: f64,
    /// Wheelbase (front–rear axle distance, m).
    pub wheelbase: f64,
}
impl WheelSetConfig {
    /// Build a symmetric 4-wheel layout from wheelbase and track width.
    ///
    /// Wheels are placed at ground level (y = 0).
    pub fn from_wheelbase_track(wb: f64, tw: f64) -> Self {
        let hw = tw * 0.5;
        let hb = wb * 0.5;
        Self {
            positions: [
                [-hw, 0.0, hb],
                [hw, 0.0, hb],
                [-hw, 0.0, -hb],
                [hw, 0.0, -hb],
            ],
            track_width: tw,
            wheelbase: wb,
        }
    }
    /// Return the rest position of a corner wheel by index (0..3).
    ///
    /// Panics if `idx >= 4`.
    pub fn corner_position(&self, idx: usize) -> [f64; 3] {
        self.positions[idx]
    }
}
/// Tire contact patch geometry (elliptical approximation).
///
/// Based on Hertz contact theory for a toroidal tire on a flat surface.
#[derive(Debug, Clone)]
pub struct ContactPatch {
    /// Half-length of the contact patch in the rolling direction (m).
    pub half_length: f64,
    /// Half-width of the contact patch (m).
    pub half_width: f64,
    /// Peak contact pressure (Pa).
    pub peak_pressure: f64,
}
impl ContactPatch {
    /// Estimate the contact patch dimensions from tire and load parameters.
    ///
    /// Uses a simplified elliptical Hertz model.
    ///
    /// * `normal_load` — normal force (N)
    /// * `tire_radius` — undeformed outer radius (m)
    /// * `tire_width` — nominal tire width (m)
    /// * `rubber_stiffness` — effective radial stiffness (N/m)
    pub fn from_load(
        normal_load: f64,
        tire_radius: f64,
        tire_width: f64,
        rubber_stiffness: f64,
    ) -> Self {
        let f = normal_load.max(0.0);
        let deflection = f / rubber_stiffness.max(1.0);
        let half_length = (2.0 * tire_radius * deflection).sqrt().max(0.0);
        let half_width = (tire_width * 0.5 * (f / (rubber_stiffness * tire_width)).sqrt())
            .clamp(0.001, tire_width * 0.5);
        let area = std::f64::consts::PI * half_length * half_width;
        let peak_pressure = if area > 1e-12 { f / area } else { 0.0 };
        Self {
            half_length,
            half_width,
            peak_pressure,
        }
    }
    /// Contact area (m²) using elliptical approximation.
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.half_length * self.half_width
    }
    /// Check whether a point (x, y) relative to the contact centre falls inside the patch.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let nx = x / self.half_length.max(1e-12);
        let ny = y / self.half_width.max(1e-12);
        nx * nx + ny * ny <= 1.0
    }
}
/// Result of a wheel raycast contact query.
#[derive(Debug, Clone)]
pub struct WheelContact {
    /// Whether the ray hit a surface.
    pub hit: bool,
    /// World-space contact position.
    pub position: [f64; 3],
    /// Surface normal at contact point (points away from surface).
    pub normal: [f64; 3],
    /// Penetration depth (metres, positive = overlapping).
    pub penetration: f64,
    /// Index of the colliding body, if applicable.
    pub body_idx: Option<usize>,
}
/// Advanced wheel spin dynamics with inertia and friction torque.
#[derive(Debug, Clone)]
pub struct WheelSpinModel {
    /// Moment of inertia (kg*m^2).
    pub inertia: f64,
    /// Current angular velocity (rad/s).
    pub omega: f64,
    /// Wheel radius (m).
    pub radius: f64,
}
impl WheelSpinModel {
    /// Create a new wheel spin model.
    pub fn new(mass: f64, radius: f64) -> Self {
        Self {
            inertia: 0.5 * mass * radius * radius,
            omega: 0.0,
            radius,
        }
    }
    /// Step the spin dynamics.
    ///
    /// * `drive_torque` — torque from engine/motor (N*m).
    /// * `brake_torque` — braking torque magnitude (N*m).
    /// * `friction_torque` — road friction torque (N*m), opposes motion.
    /// * `dt` — time step (s).
    pub fn step(&mut self, drive_torque: f64, brake_torque: f64, friction_torque: f64, dt: f64) {
        let brake_sign = if self.omega.abs() < 1e-6 {
            0.0
        } else {
            -self.omega.signum()
        };
        let net =
            drive_torque + brake_sign * brake_torque.abs() - friction_torque * self.omega.signum();
        let alpha = net / self.inertia.max(1e-10);
        self.omega += alpha * dt;
    }
    /// Current linear speed at the contact patch (m/s).
    pub fn contact_speed(&self) -> f64 {
        self.omega * self.radius
    }
    /// Kinetic energy stored in wheel rotation (J).
    pub fn rotational_energy(&self) -> f64 {
        0.5 * self.inertia * self.omega * self.omega
    }
}
/// Resonant wobble modes of a wheel assembly.
///
/// A wheel-tire assembly has several resonant modes:
/// - First lateral bending (shimmy)
/// - First longitudinal (hop)
/// - Torsional (twist)
///
/// These are modeled as independent second-order systems.
#[derive(Debug, Clone)]
pub struct WheelWobbleModes {
    /// Lateral shimmy natural frequency (Hz).
    pub shimmy_frequency: f64,
    /// Lateral shimmy damping ratio.
    pub shimmy_damping: f64,
    /// Longitudinal hop natural frequency (Hz).
    pub hop_frequency: f64,
    /// Longitudinal hop damping ratio.
    pub hop_damping: f64,
    /// Torsional natural frequency (Hz).
    pub torsional_frequency: f64,
    /// Torsional damping ratio.
    pub torsional_damping: f64,
    /// Current shimmy amplitude (rad).
    pub shimmy_amplitude: f64,
    /// Current hop amplitude (m).
    pub hop_amplitude: f64,
    /// Current torsional amplitude (rad).
    pub torsional_amplitude: f64,
}
impl WheelWobbleModes {
    /// Typical passenger car wheel wobble modes.
    pub fn passenger_car() -> Self {
        Self {
            shimmy_frequency: 8.0,
            shimmy_damping: 0.05,
            hop_frequency: 14.0,
            hop_damping: 0.10,
            torsional_frequency: 40.0,
            torsional_damping: 0.08,
            shimmy_amplitude: 0.0,
            hop_amplitude: 0.0,
            torsional_amplitude: 0.0,
        }
    }
    /// Create with explicit parameters.
    pub fn new(
        shimmy_freq: f64,
        shimmy_damp: f64,
        hop_freq: f64,
        hop_damp: f64,
        torsional_freq: f64,
        torsional_damp: f64,
    ) -> Self {
        Self {
            shimmy_frequency: shimmy_freq,
            shimmy_damping: shimmy_damp,
            hop_frequency: hop_freq,
            hop_damping: hop_damp,
            torsional_frequency: torsional_freq,
            torsional_damping: torsional_damp,
            shimmy_amplitude: 0.0,
            hop_amplitude: 0.0,
            torsional_amplitude: 0.0,
        }
    }
    /// Damped natural frequency for shimmy (Hz).
    pub fn shimmy_damped_frequency(&self) -> f64 {
        let zeta = self.shimmy_damping.clamp(0.0, 1.0);
        self.shimmy_frequency * (1.0 - zeta * zeta).sqrt()
    }
    /// Damped natural frequency for hop (Hz).
    pub fn hop_damped_frequency(&self) -> f64 {
        let zeta = self.hop_damping.clamp(0.0, 1.0);
        self.hop_frequency * (1.0 - zeta * zeta).sqrt()
    }
    /// Logarithmic decrement for shimmy mode.
    pub fn shimmy_log_decrement(&self) -> f64 {
        let z = self.shimmy_damping;
        2.0 * std::f64::consts::PI * z / (1.0 - z * z).sqrt().max(1e-12)
    }
    /// Decay rate (1/s) for shimmy: α = ζ * ωn.
    pub fn shimmy_decay_rate(&self) -> f64 {
        self.shimmy_damping * 2.0 * std::f64::consts::PI * self.shimmy_frequency
    }
    /// Step the shimmy mode given an excitation force and dt.
    ///
    /// Uses Euler integration of the SDOF oscillator:
    /// ẍ + 2ζωn ẋ + ωn² x = F/m
    ///
    /// This simplified version accumulates a state variable for amplitude.
    pub fn step_shimmy(&mut self, excitation: f64, dt: f64) {
        let _wn = 2.0 * std::f64::consts::PI * self.shimmy_frequency;
        let decay = self.shimmy_decay_rate();
        self.shimmy_amplitude = self.shimmy_amplitude * (1.0 - decay * dt) + excitation * dt;
    }
    /// Step the hop mode given a road excitation and dt.
    pub fn step_hop(&mut self, excitation: f64, dt: f64) {
        let wn = 2.0 * std::f64::consts::PI * self.hop_frequency;
        let decay = self.hop_damping * wn;
        self.hop_amplitude = self.hop_amplitude * (1.0 - decay * dt) + excitation * dt;
    }
}
/// Builder for creating [`Wheel`] configurations with a fluent API.
#[derive(Debug, Clone)]
pub struct WheelConfig {
    pub(super) wheel: Wheel,
}
impl WheelConfig {
    /// Create a new wheel config builder with default values.
    pub fn new() -> Self {
        Self {
            wheel: Wheel::default(),
        }
    }
    /// Set the wheel radius.
    pub fn radius(mut self, radius: Real) -> Self {
        self.wheel.radius = radius;
        self
    }
    /// Set the wheel width.
    pub fn width(mut self, width: Real) -> Self {
        self.wheel.width = width;
        self
    }
    /// Set the wheel mass.
    pub fn mass(mut self, mass: Real) -> Self {
        self.wheel.mass = mass;
        self
    }
    /// Set the suspension rest length.
    pub fn suspension_rest_length(mut self, length: Real) -> Self {
        self.wheel.suspension_rest_length = length;
        self
    }
    /// Set the suspension stiffness.
    pub fn suspension_stiffness(mut self, stiffness: Real) -> Self {
        self.wheel.suspension_stiffness = stiffness;
        self
    }
    /// Set the damping coefficient for compression.
    pub fn damping_compression(mut self, damping: Real) -> Self {
        self.wheel.damping_compression = damping;
        self
    }
    /// Set the damping coefficient for relaxation.
    pub fn damping_relaxation(mut self, damping: Real) -> Self {
        self.wheel.damping_relaxation = damping;
        self
    }
    /// Set the maximum suspension travel.
    pub fn max_suspension_travel(mut self, travel: Real) -> Self {
        self.wheel.max_suspension_travel = travel;
        self
    }
    /// Set the friction slip coefficient.
    pub fn friction_slip(mut self, friction: Real) -> Self {
        self.wheel.friction_slip = friction;
        self
    }
    /// Set the local-space connection point.
    pub fn connection_point(mut self, point: Vec3) -> Self {
        self.wheel.connection_point = point;
        self
    }
    /// Set the suspension direction in local space.
    pub fn suspension_direction(mut self, dir: Vec3) -> Self {
        self.wheel.suspension_direction = dir;
        self
    }
    /// Set the axle direction in local space.
    pub fn axle_direction(mut self, dir: Vec3) -> Self {
        self.wheel.axle_direction = dir;
        self
    }
    /// Build the final [`Wheel`].
    pub fn build(self) -> Wheel {
        self.wheel
    }
}
/// Tire wear model: tracks cumulative wear and degrades grip.
#[derive(Debug, Clone)]
pub struct TireWear {
    /// Tire compound.
    pub compound: TireCompound,
    /// Base wear rate coefficient (wear units per (slip * load * dt)).
    pub wear_rate_coeff: f64,
    /// Current normalised wear (0.0 = new, 1.0 = fully worn).
    pub current_wear: f64,
}
impl TireWear {
    /// Create a new tire wear tracker.
    pub fn new(compound: TireCompound, wear_rate_coeff: f64) -> Self {
        Self {
            compound,
            wear_rate_coeff,
            current_wear: 0.0,
        }
    }
    /// Advance wear given slip ratio, load (N) and timestep (s).
    pub fn update(&mut self, slip: f64, load: f64, dt: f64) {
        let delta = self.wear_rate_coeff * slip.abs() * load.abs() * dt;
        self.current_wear = (self.current_wear + delta).clamp(0.0, 1.0);
    }
    /// Return a grip multiplier (1.0 = fresh tires, 0.0 = fully worn).
    ///
    /// A worn tire loses grip quadratically; compound affects the starting grip.
    pub fn grip_factor(&self) -> f64 {
        let base = match self.compound {
            TireCompound::Soft => 1.0,
            TireCompound::Medium => 0.95,
            TireCompound::Hard => 0.90,
        };
        base * (1.0 - self.current_wear * self.current_wear)
    }
}
/// Linear spring-damper for wheel suspension.
#[derive(Debug, Clone)]
pub struct WheelDamper {
    /// Spring stiffness (N/m).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Rest (unloaded) length (m).
    pub rest_length: f64,
}
impl WheelDamper {
    /// Create a new wheel damper.
    pub fn new(stiffness: f64, damping: f64, rest_length: f64) -> Self {
        Self {
            stiffness,
            damping,
            rest_length,
        }
    }
    /// Compute suspension force.
    ///
    /// * `compression` — how much the spring is compressed (positive = compressed, m).
    /// * `velocity`    — rate of change of compression (positive = compressing faster, m/s).
    ///
    /// Returns the upward force (N).  Positive = pushing chassis up.
    pub fn suspension_force(&self, compression: f64, velocity: f64) -> f64 {
        self.stiffness * compression + self.damping * velocity
    }
}
/// Tire drumming noise model.
///
/// Drumming noise is caused by tread block impact on the road surface.  The
/// dominant frequencies are related to the tire rotation speed and the number
/// of tread blocks (pitch sequence).
#[derive(Debug, Clone)]
pub struct DrummingNoiseModel {
    /// Number of tread pitches around the circumference.
    pub tread_pitches: usize,
    /// Tire outer radius (m).
    pub outer_radius: f64,
    /// Noise amplitude gain (dimensionless).
    pub gain: f64,
    /// Acoustic radiation efficiency.
    pub radiation_efficiency: f64,
}
impl DrummingNoiseModel {
    /// Create a new drumming noise model.
    pub fn new(tread_pitches: usize, outer_radius: f64, gain: f64) -> Self {
        Self {
            tread_pitches,
            outer_radius,
            gain,
            radiation_efficiency: 0.5,
        }
    }
    /// Fundamental drumming frequency (Hz) at vehicle speed `v` (m/s).
    ///
    /// f = v * N / (2π * r)
    pub fn drumming_frequency(&self, speed: f64) -> f64 {
        let circumference = 2.0 * std::f64::consts::PI * self.outer_radius;
        speed * self.tread_pitches as f64 / circumference.max(1e-6)
    }
    /// Harmonic series of drumming frequencies.
    ///
    /// Returns the first `n_harmonics` harmonic frequencies (Hz).
    pub fn harmonics(&self, speed: f64, n_harmonics: usize) -> Vec<f64> {
        let f1 = self.drumming_frequency(speed);
        (1..=n_harmonics).map(|k| f1 * k as f64).collect()
    }
    /// Estimated sound pressure level (dB) relative to `p_ref = 20 μPa`.
    ///
    /// SPL = 20 * log10(gain * radiation_efficiency * v / v_ref)
    pub fn sound_pressure_level(&self, speed: f64) -> f64 {
        if speed < 1e-6 {
            return 0.0;
        }
        let v_ref = 1.0;
        20.0 * (self.gain * self.radiation_efficiency * speed / v_ref)
            .log10()
            .max(0.0)
    }
    /// Drumming noise increases with the square of speed (empirical).
    pub fn noise_index(&self, speed: f64) -> f64 {
        self.gain * speed * speed
    }
}
/// Tire inflation pressure effects on handling and performance.
///
/// Inflation pressure affects rolling resistance, cornering stiffness,
/// contact patch size, and carcass deflection.
#[derive(Debug, Clone)]
pub struct InflationPressureEffects {
    /// Nominal inflation pressure (Pa). Typically 2.0–2.5 bar.
    pub nominal_pressure: f64,
    /// Rolling resistance sensitivity (change in Cr per Pa).
    pub rolling_resistance_sensitivity: f64,
    /// Cornering stiffness sensitivity (change in Cf per Pa).
    pub cornering_stiffness_sensitivity: f64,
    /// Contact patch area at nominal pressure (m²).
    pub nominal_contact_area: f64,
}
impl InflationPressureEffects {
    /// Create a model for a typical passenger car tire.
    pub fn passenger_car() -> Self {
        Self {
            nominal_pressure: 220_000.0,
            rolling_resistance_sensitivity: -2e-9,
            cornering_stiffness_sensitivity: 5.0,
            nominal_contact_area: 0.015,
        }
    }
    /// Rolling resistance coefficient at the given inflation pressure.
    pub fn rolling_resistance(&self, pressure: f64, base_cr: f64) -> f64 {
        let dp = pressure - self.nominal_pressure;
        (base_cr + self.rolling_resistance_sensitivity * dp).max(0.001)
    }
    /// Cornering stiffness at the given inflation pressure (N/rad).
    pub fn cornering_stiffness(&self, pressure: f64, base_cf: f64) -> f64 {
        let dp = pressure - self.nominal_pressure;
        (base_cf + self.cornering_stiffness_sensitivity * dp).max(0.0)
    }
    /// Contact patch area at the given inflation pressure (m²).
    ///
    /// Higher pressure → smaller contact patch (inverse relationship).
    pub fn contact_area(&self, pressure: f64, normal_load: f64) -> f64 {
        if pressure < 1e-6 {
            return self.nominal_contact_area;
        }
        (normal_load / pressure).max(0.001)
    }
    /// Tire load capacity ratio: load / (pressure * nominal_area).
    pub fn load_capacity_ratio(&self, normal_load: f64, pressure: f64) -> f64 {
        let max_load = pressure * self.nominal_contact_area;
        if max_load < 1e-6 {
            return f64::INFINITY;
        }
        normal_load / max_load
    }
    /// Under-inflation penalty on rolling resistance (multiplicative factor).
    pub fn under_inflation_penalty(&self, pressure: f64) -> f64 {
        if pressure >= self.nominal_pressure {
            return 1.0;
        }
        let ratio = pressure / self.nominal_pressure.max(1.0);
        1.0 + 0.5 * (1.0 - ratio)
    }
}
/// Radial and lateral carcass stiffness for a pneumatic tire.
///
/// The carcass behaves like a spring in both radial (vertical) and lateral
/// directions; these stiffnesses couple into the effective contact mechanics.
#[derive(Debug, Clone)]
pub struct CarcassStiffness {
    /// Radial stiffness (N/m) — resistance to vertical deflection.
    pub k_radial: f64,
    /// Lateral stiffness (N/m) — resistance to lateral deflection.
    pub k_lateral: f64,
    /// Torsional stiffness (N·m/rad) — resistance to twist about vertical axis.
    pub k_torsional: f64,
    /// Inflation pressure (Pa).
    pub inflation_pressure: f64,
}
impl CarcassStiffness {
    /// Typical stiffness for a 225/45R17 road tire at 2.2 bar.
    pub fn passenger_225_45_r17() -> Self {
        Self {
            k_radial: 200_000.0,
            k_lateral: 120_000.0,
            k_torsional: 800.0,
            inflation_pressure: 220_000.0,
        }
    }
    /// Formula racing tire (stiff, low-profile).
    pub fn formula_slick() -> Self {
        Self {
            k_radial: 280_000.0,
            k_lateral: 180_000.0,
            k_torsional: 1_200.0,
            inflation_pressure: 130_000.0,
        }
    }
    /// Effective radial stiffness including inflation pressure contribution.
    ///
    /// k_eff = k_radial + pressure * contact_area
    pub fn effective_radial_stiffness(&self, contact_area: f64) -> f64 {
        self.k_radial + self.inflation_pressure * contact_area
    }
    /// Radial deflection (m) for a given normal load (N).
    pub fn radial_deflection(&self, normal_load: f64, contact_area: f64) -> f64 {
        let k_eff = self.effective_radial_stiffness(contact_area);
        normal_load / k_eff.max(1.0)
    }
    /// Lateral deflection (m) for a given lateral force (N).
    pub fn lateral_deflection(&self, lateral_force: f64) -> f64 {
        lateral_force / self.k_lateral.max(1.0)
    }
}
