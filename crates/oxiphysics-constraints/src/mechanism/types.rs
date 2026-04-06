//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

use std::f64::consts::PI;

/// Force analysis result for a slider-crank mechanism.
#[allow(dead_code)]
pub struct SliderCrankForces {
    /// Gas/load force on piston (N).
    pub piston_force: f64,
    /// Connecting rod force (N).
    pub rod_force: f64,
    /// Crankpin tangential force (N) — the useful torque-producing component.
    pub crankpin_tangential: f64,
    /// Crankpin radial force (N).
    pub crankpin_radial: f64,
    /// Torque at crankshaft (N·m).
    pub crank_torque: f64,
}
impl SliderCrankForces {
    /// Compute forces in the slider-crank given piston load and crank geometry.
    ///
    /// `theta` is crank angle (rad), `phi` is connecting-rod angle (rad).
    #[allow(dead_code)]
    pub fn compute(crank_radius: f64, theta: f64, phi: f64, piston_force: f64) -> Self {
        let rod_force = piston_force / phi.cos().abs().max(1e-9);
        let crankpin_tangential = rod_force * (theta + phi).sin();
        let crankpin_radial = rod_force * (theta + phi).cos();
        let crank_torque = crankpin_tangential * crank_radius;
        Self {
            piston_force,
            rod_force,
            crankpin_tangential,
            crankpin_radial,
            crank_torque,
        }
    }
}
/// Extended toggle clamp with full kinematics.
#[allow(dead_code)]
pub struct ToggleClampExtended {
    /// Crank link length `a` (metres) — driven by input force.
    pub crank_length: f64,
    /// Toggle link length `b` (metres).
    pub toggle_length: f64,
    /// Pivot separation `c` (metres).
    pub pivot_sep: f64,
    /// Current input link angle from horizontal (radians).
    pub theta: f64,
    /// Input force (N).
    pub f_input: f64,
    /// Friction angle at crank pivot (radians).
    pub friction_angle: f64,
}
impl ToggleClampExtended {
    /// Create a toggle clamp with given geometry.
    pub fn new(crank_length: f64, toggle_length: f64, pivot_sep: f64) -> Self {
        Self {
            crank_length,
            toggle_length,
            pivot_sep,
            theta: 30.0_f64.to_radians(),
            f_input: 0.0,
            friction_angle: 2.0_f64.to_radians(),
        }
    }
    /// Coupler point position (output end of toggle link), x coordinate.
    pub fn output_x(&self) -> f64 {
        self.crank_length * self.theta.cos() + self.toggle_length * self.toggle_angle().cos()
    }
    /// Toggle link angle from horizontal (radians).
    pub fn toggle_angle(&self) -> f64 {
        let cx = self.crank_length * self.theta.cos();
        let cy = self.crank_length * self.theta.sin();
        let dx = self.pivot_sep - cx;
        let dy = -cy;
        dy.atan2(dx)
    }
    /// Mechanical advantage (ideal).
    pub fn ideal_mechanical_advantage(&self) -> f64 {
        let phi = self.toggle_angle();
        let sin_phi = phi.sin().abs().max(1e-9);
        self.toggle_length / (self.crank_length * self.theta.sin().abs().max(1e-9)) * sin_phi
    }
    /// Output clamping force (N).
    pub fn output_force(&self) -> f64 {
        self.f_input * self.ideal_mechanical_advantage()
    }
    /// Check if near dead-centre (toggle angle ≈ 0 or π).
    pub fn is_near_dead_centre(&self) -> bool {
        let phi = self.toggle_angle().abs();
        phi < 5.0_f64.to_radians() || (phi - PI).abs() < 5.0_f64.to_radians()
    }
}
/// A locking mechanism (ratchet or detent).
#[allow(dead_code)]
pub struct LockingMechanism {
    /// Kind of locking mechanism.
    pub kind: LockingKind,
    /// Current angular position (radians).
    pub angle: f64,
    /// Whether the mechanism is currently locked.
    pub locked: bool,
}
impl LockingMechanism {
    /// Create a ratchet locking mechanism.
    pub fn new_ratchet(teeth: u32, spring_force: f64) -> Self {
        Self {
            kind: LockingKind::Ratchet {
                teeth,
                pawl_spring_force: spring_force,
            },
            angle: 0.0,
            locked: false,
        }
    }
    /// Create a detent locking mechanism.
    pub fn new_detent(ball_radius: f64, spring_force: f64) -> Self {
        Self {
            kind: LockingKind::Detent {
                ball_radius,
                spring_force,
            },
            angle: 0.0,
            locked: false,
        }
    }
    /// Advance the mechanism by `delta_angle` radians.
    ///
    /// For a ratchet, negative rotation is blocked.
    pub fn try_advance(&mut self, delta_angle: f64) -> bool {
        match &self.kind {
            LockingKind::Ratchet { teeth, .. } => {
                if delta_angle < 0.0 {
                    self.locked = true;
                    false
                } else {
                    let tooth_angle = 2.0 * PI / *teeth as f64;
                    self.angle += delta_angle;
                    self.locked = false;
                    self.angle = (self.angle / tooth_angle).round() * tooth_angle;
                    true
                }
            }
            LockingKind::Detent { .. } => {
                self.angle += delta_angle;
                self.locked = false;
                true
            }
        }
    }
    /// Detent engagement force required to move past a detent.
    pub fn detent_force(&self) -> f64 {
        match &self.kind {
            LockingKind::Detent {
                ball_radius,
                spring_force,
            } => spring_force * ball_radius,
            _ => 0.0,
        }
    }
    /// Pawl engagement force for a ratchet.
    pub fn pawl_force(&self) -> f64 {
        match &self.kind {
            LockingKind::Ratchet {
                pawl_spring_force, ..
            } => *pawl_spring_force,
            _ => 0.0,
        }
    }
}
/// Kinematic pair (joint) classification for mobility calculation.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JointClass {
    /// Revolute joint — 1 DOF.
    Revolute,
    /// Prismatic joint — 1 DOF.
    Prismatic,
    /// Cylindrical joint — 2 DOF.
    Cylindrical,
    /// Universal (Hooke's) joint — 2 DOF.
    Universal,
    /// Spherical ball joint — 3 DOF.
    Spherical,
    /// Plane (sliding surface) — 3 DOF.
    Planar,
    /// Screw joint — 1 DOF.
    Screw,
    /// Rigid (welded) joint — 0 DOF (adds a constraint).
    Rigid,
}
impl JointClass {
    /// Degrees of freedom allowed by this joint type.
    pub fn dof(&self) -> u32 {
        match self {
            JointClass::Revolute | JointClass::Prismatic | JointClass::Screw => 1,
            JointClass::Cylindrical | JointClass::Universal => 2,
            JointClass::Spherical | JointClass::Planar => 3,
            JointClass::Rigid => 0,
        }
    }
}
/// A rotating mass for balancing analysis.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct RotatingMass {
    /// Mass (kg).
    pub mass: f64,
    /// Radius from rotation axis (metres).
    pub radius: f64,
    /// Angular position (radians) from reference.
    pub angle: f64,
    /// Axial position along rotation axis (metres).
    pub axial_pos: f64,
}
impl RotatingMass {
    /// Create a rotating mass.
    pub fn new(mass: f64, radius: f64, angle_deg: f64, axial_pos: f64) -> Self {
        Self {
            mass,
            radius,
            angle: angle_deg.to_radians(),
            axial_pos,
        }
    }
    /// Centrifugal force vector (m·r) in the x-y plane.
    pub fn mr_vector(&self) -> [f64; 2] {
        let mr = self.mass * self.radius;
        [mr * self.angle.cos(), mr * self.angle.sin()]
    }
    /// Moment vector (m·r·l) for dynamic balancing.
    pub fn mrl_vector(&self) -> [f64; 2] {
        let mr = self.mass * self.radius;
        let mrl = mr * self.axial_pos;
        [mrl * self.angle.cos(), mrl * self.angle.sin()]
    }
}
/// Spur gear pair.
#[allow(dead_code)]
pub struct GearPair {
    /// Number of teeth on drive gear.
    pub teeth_drive: u32,
    /// Number of teeth on driven gear.
    pub teeth_driven: u32,
    /// Module (tooth size parameter, metres).
    pub module: f64,
    /// Pressure angle (radians), typically 20°.
    pub pressure_angle: f64,
    /// Backlash (metres).
    pub backlash: f64,
    /// Angular velocity of drive gear (rad/s).
    pub omega_drive: f64,
}
impl GearPair {
    /// Create a gear pair.
    pub fn new(teeth_drive: u32, teeth_driven: u32, module: f64) -> Self {
        Self {
            teeth_drive,
            teeth_driven,
            module,
            pressure_angle: 20.0_f64.to_radians(),
            backlash: 0.0,
            omega_drive: 0.0,
        }
    }
    /// Gear ratio (driven / drive).
    pub fn gear_ratio(&self) -> f64 {
        self.teeth_driven as f64 / self.teeth_drive as f64
    }
    /// Angular velocity of driven gear.
    pub fn omega_driven(&self) -> f64 {
        -self.omega_drive / self.gear_ratio()
    }
    /// Pitch radius of drive gear (metres).
    pub fn pitch_radius_drive(&self) -> f64 {
        self.module * self.teeth_drive as f64 / 2.0
    }
    /// Pitch radius of driven gear (metres).
    pub fn pitch_radius_driven(&self) -> f64 {
        self.module * self.teeth_driven as f64 / 2.0
    }
    /// Centre distance between gear axes (metres).
    pub fn centre_distance(&self) -> f64 {
        self.pitch_radius_drive() + self.pitch_radius_driven()
    }
    /// Compute involute profile x-coordinate at angle `phi` for drive gear.
    pub fn involute_x(&self, phi: f64) -> f64 {
        let r = self.pitch_radius_drive();
        r * (phi.cos() + phi * phi.sin())
    }
    /// Compute involute profile y-coordinate at angle `phi` for drive gear.
    pub fn involute_y(&self, phi: f64) -> f64 {
        let r = self.pitch_radius_drive();
        r * (phi.sin() - phi * phi.cos())
    }
}
/// Disc cam follower types.
#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CamFollowerType {
    /// Flat-faced follower.
    FlatFaced,
    /// Roller follower.
    Roller,
    /// Knife-edge follower.
    KnifeEdge,
}
/// Hooke's joint (universal joint / Cardan joint) kinematics.
///
/// Relates the angular velocity ratio between input and output shafts.
#[allow(dead_code)]
pub struct HookesJoint {
    /// Angle between the two shaft axes (radians).
    pub shaft_angle: f64,
    /// Input shaft angular velocity (rad/s).
    pub omega_in: f64,
    /// Current input shaft angle (radians).
    pub theta: f64,
}
impl HookesJoint {
    /// Create a Hooke's joint model.
    pub fn new(shaft_angle_deg: f64, omega_in: f64) -> Self {
        Self {
            shaft_angle: shaft_angle_deg.to_radians(),
            omega_in,
            theta: 0.0,
        }
    }
    /// Velocity ratio ω₂/ω₁ at input angle `theta` (radians).
    ///
    /// `ω₂/ω₁ = cos(β) / (1 - sin²(β) * cos²(θ))`.
    pub fn velocity_ratio(&self, theta: f64) -> f64 {
        let beta = self.shaft_angle;
        let num = beta.cos();
        let den = 1.0 - beta.sin().powi(2) * theta.cos().powi(2);
        num / den.max(1e-15)
    }
    /// Output shaft angular velocity (rad/s) at current `theta`.
    pub fn omega_out(&self) -> f64 {
        self.omega_in * self.velocity_ratio(self.theta)
    }
    /// Angular acceleration of output shaft (rad/s²).
    ///
    /// `α₂ = ω₁² * cos(β) * sin²(β) * sin(2θ) / (1 - sin²(β) * cos²(θ))²`.
    pub fn alpha_out(&self) -> f64 {
        let beta = self.shaft_angle;
        let theta = self.theta;
        let den = 1.0 - beta.sin().powi(2) * theta.cos().powi(2);
        self.omega_in.powi(2) * beta.cos() * beta.sin().powi(2) * (2.0 * theta).sin()
            / den.powi(2).max(1e-15)
    }
    /// Maximum speed variation: ratio of max to min ω₂.
    pub fn speed_variation(&self) -> f64 {
        let beta = self.shaft_angle;
        let cos_b = beta.cos();
        if cos_b < 1e-15 {
            return f64::INFINITY;
        }
        1.0 / (cos_b * cos_b)
    }
    /// Check if a double Cardan joint arrangement is needed (large angle).
    pub fn needs_double_cardan(&self) -> bool {
        self.shaft_angle > 15.0_f64.to_radians()
    }
    /// Phase angle needed for equal shaft phases in a double Cardan shaft.
    pub fn double_cardan_phase(&self) -> f64 {
        self.shaft_angle
    }
    /// Advance input shaft by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.theta = (self.theta + self.omega_in * dt) % (2.0 * PI);
    }
}
/// Worm gear set.
#[allow(dead_code)]
pub struct WormGear {
    /// Number of worm starts (thread count).
    pub worm_starts: u32,
    /// Number of worm wheel teeth.
    pub wheel_teeth: u32,
    /// Lead angle (radians).
    pub lead_angle: f64,
    /// Coefficient of friction.
    pub friction_coeff: f64,
    /// Input angular velocity (rad/s).
    pub omega_in: f64,
}
impl WormGear {
    /// Create a worm gear set.
    pub fn new(worm_starts: u32, wheel_teeth: u32, lead_angle: f64, friction_coeff: f64) -> Self {
        Self {
            worm_starts,
            wheel_teeth,
            lead_angle,
            friction_coeff,
            omega_in: 0.0,
        }
    }
    /// Speed reduction ratio.
    pub fn reduction_ratio(&self) -> f64 {
        self.wheel_teeth as f64 / self.worm_starts as f64
    }
    /// Output angular velocity (rad/s).
    pub fn omega_out(&self) -> f64 {
        self.omega_in / self.reduction_ratio()
    }
    /// Check if self-locking (worm cannot be back-driven).
    pub fn is_self_locking(&self) -> bool {
        self.lead_angle.tan() < self.friction_coeff
    }
    /// Mechanical efficiency (forward direction).
    pub fn efficiency(&self) -> f64 {
        let lambda = self.lead_angle;
        let mu = self.friction_coeff;
        (lambda.tan()) / (lambda.tan() + mu / lambda.cos())
    }
}
/// A planar four-bar linkage (ground – crank – coupler – follower).
///
/// Link lengths: `l1` = ground, `l2` = crank, `l3` = coupler, `l4` = follower.
#[allow(dead_code)]
pub struct FourBarLinkage {
    /// Ground link length (metres).
    pub l1: f64,
    /// Crank link length (metres).
    pub l2: f64,
    /// Coupler link length (metres).
    pub l3: f64,
    /// Follower link length (metres).
    pub l4: f64,
    /// Current crank angle (radians).
    pub theta2: f64,
}
impl FourBarLinkage {
    /// Create a four-bar linkage.
    pub fn new(l1: f64, l2: f64, l3: f64, l4: f64) -> Self {
        Self {
            l1,
            l2,
            l3,
            l4,
            theta2: 0.0,
        }
    }
    /// Check Grashof condition.
    ///
    /// Returns `true` if at least one link can make a full rotation.
    pub fn grashof(&self) -> bool {
        let s = [self.l1, self.l2, self.l3, self.l4]
            .iter()
            .cloned()
            .fold(f64::MAX, f64::min);
        let l = [self.l1, self.l2, self.l3, self.l4]
            .iter()
            .cloned()
            .fold(f64::MIN, f64::max);
        let sum_rest: f64 = [self.l1, self.l2, self.l3, self.l4].iter().sum::<f64>() - s - l;
        s + l <= sum_rest
    }
    /// Returns `true` if the linkage is a crank-rocker (l2 = crank can rotate fully).
    pub fn is_crank_rocker(&self) -> bool {
        self.grashof() && self.l2 <= self.l1 && self.l2 <= self.l3 && self.l2 <= self.l4
    }
    /// Compute the follower angle `theta4` for the current `theta2`.
    ///
    /// Uses the Freudenstein equation. Returns `None` on assembly error.
    pub fn follower_angle(&self, theta2: f64) -> Option<f64> {
        let k1 = self.l1 / self.l4;
        let k2 = self.l1 / self.l2;
        let k3 =
            (sq(self.l2) - sq(self.l3) + sq(self.l4) + sq(self.l1)) / (2.0 * self.l2 * self.l4);
        let a = theta2.cos() - k1 - k2 * theta2.cos() + k3;
        let b = -2.0 * theta2.sin();
        let c = k1 - (k2 + 1.0) * theta2.cos() + k3;
        let disc = sq(b) - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let t = (-b - disc.sqrt()) / (2.0 * a);
        Some(2.0 * t.atan())
    }
    /// Compute the transmission angle (angle between coupler and follower, radians).
    pub fn transmission_angle(&self, theta2: f64) -> f64 {
        let cos_mu = (sq(self.l3) + sq(self.l4) - sq(self.l1) - sq(self.l2)
            + 2.0 * self.l1 * self.l2 * theta2.cos())
            / (2.0 * self.l3 * self.l4);
        cos_mu.clamp(-1.0, 1.0).acos()
    }
}
/// Planetary gear set (simple epicyclic).
#[allow(dead_code)]
pub struct PlanetaryGear {
    /// Number of teeth on sun gear.
    pub sun_teeth: u32,
    /// Number of teeth on ring gear.
    pub ring_teeth: u32,
    /// Number of planet gears.
    pub planet_count: u32,
    /// Angular velocity of sun gear (rad/s).
    pub omega_sun: f64,
    /// Angular velocity of ring gear (rad/s).
    pub omega_ring: f64,
    /// Angular velocity of carrier (rad/s).
    pub omega_carrier: f64,
}
impl PlanetaryGear {
    /// Create a planetary gear set.
    pub fn new(sun_teeth: u32, ring_teeth: u32, planet_count: u32) -> Self {
        Self {
            sun_teeth,
            ring_teeth,
            planet_count,
            omega_sun: 0.0,
            omega_ring: 0.0,
            omega_carrier: 0.0,
        }
    }
    /// Number of teeth on each planet gear.
    pub fn planet_teeth(&self) -> u32 {
        (self.ring_teeth - self.sun_teeth) / 2
    }
    /// Gear ratio: omega_sun / omega_ring when carrier is fixed.
    pub fn sun_ring_ratio(&self) -> f64 {
        -(self.ring_teeth as f64 / self.sun_teeth as f64)
    }
    /// Carrier velocity given sun and ring angular velocities (Willis equation).
    pub fn carrier_velocity(&self) -> f64 {
        let r = self.ring_teeth as f64;
        let s = self.sun_teeth as f64;
        (s * self.omega_sun + r * self.omega_ring) / (s + r)
    }
    /// Planet angular velocity relative to the carrier.
    pub fn planet_velocity_relative(&self) -> f64 {
        let n_s = self.sun_teeth as f64;
        let n_p = self.planet_teeth() as f64;
        (self.omega_sun - self.omega_carrier) * n_s / n_p
    }
}
/// Rack and pinion mechanism.
#[allow(dead_code)]
pub struct RackAndPinion {
    /// Pinion pitch radius (metres).
    pub pinion_radius: f64,
    /// Number of teeth on pinion.
    pub pinion_teeth: u32,
    /// Module.
    pub module: f64,
    /// Current pinion angle (radians).
    pub pinion_angle: f64,
}
impl RackAndPinion {
    /// Create a rack-and-pinion.
    pub fn new(pinion_teeth: u32, module: f64) -> Self {
        let pinion_radius = module * pinion_teeth as f64 / 2.0;
        Self {
            pinion_radius,
            pinion_teeth,
            module,
            pinion_angle: 0.0,
        }
    }
    /// Compute rack displacement (metres) for a given pinion rotation (radians).
    pub fn rack_displacement(&self, pinion_rotation: f64) -> f64 {
        self.pinion_radius * pinion_rotation
    }
    /// Compute pinion rotation (radians) for a given rack displacement (metres).
    pub fn pinion_rotation(&self, rack_displacement: f64) -> f64 {
        rack_displacement / self.pinion_radius
    }
    /// Velocity amplification: rack speed (m/s) from pinion angular speed (rad/s).
    pub fn rack_velocity(&self, omega: f64) -> f64 {
        self.pinion_radius * omega
    }
}
/// Slider-crank mechanism (piston and crank).
#[allow(dead_code)]
pub struct SliderCrankMechanism {
    /// Crank radius (metres).
    pub crank_radius: f64,
    /// Connecting rod length (metres).
    pub rod_length: f64,
    /// Crank angular velocity (rad/s).
    pub omega: f64,
    /// Current crank angle (radians).
    pub theta: f64,
}
impl SliderCrankMechanism {
    /// Create a slider-crank mechanism.
    pub fn new(crank_radius: f64, rod_length: f64, omega: f64) -> Self {
        Self {
            crank_radius,
            rod_length,
            omega,
            theta: 0.0,
        }
    }
    /// Compute piston position (displacement from crank centre) at crank angle `theta`.
    pub fn piston_position(&self, theta: f64) -> f64 {
        let r = self.crank_radius;
        let l = self.rod_length;
        r * theta.cos() + (sq(l) - sq(r * theta.sin())).sqrt()
    }
    /// Compute piston velocity (m/s) at crank angle `theta`.
    pub fn piston_velocity(&self, theta: f64) -> f64 {
        let r = self.crank_radius;
        let l = self.rod_length;
        let omega = self.omega;
        -r * omega
            * (theta.sin() + (r * theta.sin() * theta.cos()) / (sq(l) - sq(r * theta.sin())).sqrt())
    }
    /// Compute piston acceleration (m/s²) at crank angle `theta`.
    pub fn piston_acceleration(&self, theta: f64) -> f64 {
        let r = self.crank_radius;
        let l = self.rod_length;
        let omega = self.omega;
        let lambda = r / l;
        -r * sq(omega) * (theta.cos() + lambda * (2.0 * theta).cos())
    }
    /// Advance crank by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.theta = (self.theta + self.omega * dt) % (2.0 * PI);
    }
}
/// Over-centre toggle clamp.
#[allow(dead_code)]
pub struct ToggleClamp {
    /// Toggle link length (metres).
    pub link_length: f64,
    /// Pivot eccentricity from dead-centre (metres).
    pub eccentricity: f64,
    /// Input force (N).
    pub input_force: f64,
    /// Current toggle angle from dead-centre (radians).
    pub angle: f64,
}
impl ToggleClamp {
    /// Create a toggle clamp.
    pub fn new(link_length: f64, eccentricity: f64) -> Self {
        Self {
            link_length,
            eccentricity,
            input_force: 0.0,
            angle: PI / 4.0,
        }
    }
    /// Mechanical advantage at the current angle.
    pub fn mechanical_advantage(&self) -> f64 {
        if self.angle.abs() < 1e-9 {
            f64::MAX
        } else {
            1.0 / (2.0 * self.angle.tan().abs().max(1e-9))
        }
    }
    /// Output clamping force (N).
    pub fn clamp_force(&self) -> f64 {
        self.input_force * self.mechanical_advantage().min(1e6)
    }
    /// Returns `true` when the mechanism is past dead-centre (locked).
    pub fn is_locked(&self) -> bool {
        self.angle.abs() < self.eccentricity / self.link_length
    }
}
/// Diamond scissor lift linkage.
#[allow(dead_code)]
pub struct ScissorLift {
    /// Arm length (metres).
    pub arm_length: f64,
    /// Number of scissor stages.
    pub stages: u32,
    /// Current input stroke (horizontal displacement, metres).
    pub stroke: f64,
    /// Base width at zero stroke (metres).
    pub base_width: f64,
}
impl ScissorLift {
    /// Create a scissor lift.
    pub fn new(arm_length: f64, stages: u32) -> Self {
        let base_width = arm_length * 1.8;
        Self {
            arm_length,
            stages,
            stroke: 0.0,
            base_width,
        }
    }
    /// Platform height (metres) for a given horizontal stroke.
    pub fn platform_height(&self, stroke: f64) -> f64 {
        let half_base = (self.base_width - stroke).max(0.01) / 2.0;
        let arm = self.arm_length;
        let h_single = (sq(arm) - sq(half_base)).max(0.0).sqrt();
        h_single * self.stages as f64
    }
    /// Mechanical advantage (ratio of output force to input force) at current stroke.
    pub fn mechanical_advantage(&self, stroke: f64) -> f64 {
        let h = self.platform_height(stroke);
        if h < 1e-6 {
            return 0.0;
        }
        h / (self.base_width / 2.0 - stroke / 2.0).max(1e-6)
    }
}
/// Disc cam and follower model.
#[allow(dead_code)]
pub struct CamFollower {
    /// Base circle radius (metres).
    pub base_radius: f64,
    /// Total lift height (metres).
    pub lift: f64,
    /// Rise angle (radians).
    pub rise_angle: f64,
    /// Follower type.
    pub follower_type: CamFollowerType,
    /// Cam profile type.
    pub profile_type: CamProfileType,
    /// Current cam angle (radians).
    pub theta: f64,
    /// Angular velocity (rad/s).
    pub omega: f64,
}
impl CamFollower {
    /// Create a cam follower.
    pub fn new(base_radius: f64, lift: f64, rise_angle: f64, omega: f64) -> Self {
        Self {
            base_radius,
            lift,
            rise_angle,
            follower_type: CamFollowerType::Roller,
            profile_type: CamProfileType::Cycloidal,
            theta: 0.0,
            omega,
        }
    }
    /// Compute follower displacement for a normalised rise fraction `tau` ∈ \[0, 1\].
    pub fn displacement(&self, tau: f64) -> f64 {
        let tau = tau.clamp(0.0, 1.0);
        match self.profile_type {
            CamProfileType::Polynomial345 => {
                self.lift * (10.0 * tau.powi(3) - 15.0 * tau.powi(4) + 6.0 * tau.powi(5))
            }
            CamProfileType::Cycloidal => self.lift * (tau - (2.0 * PI * tau).sin() / (2.0 * PI)),
            CamProfileType::Harmonic => self.lift * 0.5 * (1.0 - (PI * tau).cos()),
        }
    }
    /// Compute follower position for the current cam angle.
    pub fn follower_position(&self) -> f64 {
        let tau = (self.theta % self.rise_angle) / self.rise_angle;
        self.base_radius + self.displacement(tau)
    }
    /// Compute pressure angle (radians) — angle between follower motion axis and cam normal.
    pub fn pressure_angle(&self, tau: f64) -> f64 {
        let ds = self.displacement_prime(tau) / self.rise_angle;
        let r = self.base_radius + self.displacement(tau);
        (ds / r).atan()
    }
    /// Numerical derivative of displacement w.r.t. tau.
    fn displacement_prime(&self, tau: f64) -> f64 {
        let h = 1e-5;
        (self.displacement(tau + h) - self.displacement(tau - h)) / (2.0 * h)
    }
    /// Advance cam by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.theta = (self.theta + self.omega * dt) % (2.0 * PI);
    }
}
/// Disc cam profile type.
#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CamProfileType {
    /// Polynomial (3-4-5) rise profile.
    Polynomial345,
    /// Cycloidal rise profile.
    Cycloidal,
    /// Harmonic rise profile.
    Harmonic,
}
/// A mechanism described by links and joints for mobility analysis.
#[allow(dead_code)]
pub struct MechanismGraph {
    /// Number of links (including ground link).
    pub n_links: u32,
    /// Joints in the mechanism.
    pub joints: Vec<JointClass>,
    /// Whether the analysis is planar (true) or spatial (false).
    pub planar: bool,
}
impl MechanismGraph {
    /// Create a new mechanism graph.
    pub fn new(n_links: u32, planar: bool) -> Self {
        Self {
            n_links,
            joints: Vec::new(),
            planar,
        }
    }
    /// Add a joint to the mechanism.
    pub fn add_joint(&mut self, joint: JointClass) {
        self.joints.push(joint);
    }
    /// Compute mobility (DOF) of the mechanism.
    pub fn mobility(&self) -> i32 {
        if self.planar {
            grubler_kutzbach_planar(self.n_links, &self.joints)
        } else {
            grubler_kutzbach_spatial(self.n_links, &self.joints)
        }
    }
    /// Number of joints (pairs).
    pub fn joint_count(&self) -> usize {
        self.joints.len()
    }
    /// Check if the mechanism is a structure (M ≤ 0).
    pub fn is_structure(&self) -> bool {
        self.mobility() <= 0
    }
    /// Check if the mechanism has a unique DOF input (M == 1).
    pub fn is_single_dof(&self) -> bool {
        self.mobility() == 1
    }
}
/// Grashof classification of a four-bar linkage.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GrashofClass {
    /// Grashof: at least one link can rotate fully — crank-rocker.
    CrankRocker,
    /// Grashof: rocker-crank (follower is the shortest, can fully rotate).
    RockerCrank,
    /// Grashof: double-crank (double-drag link).
    DoubleCrank,
    /// Grashof: double-rocker.
    DoubleRocker,
    /// Non-Grashof: no link can make a full rotation.
    NonGrashof,
    /// Change-point mechanism (special Grashof case with collinear links).
    ChangePoint,
}
/// Geneva drive (Maltese cross) mechanism.
///
/// Converts continuous rotation into intermittent indexing motion.
#[allow(dead_code)]
pub struct GenevaMechanism {
    /// Number of slots on the Geneva wheel.
    pub n_slots: u32,
    /// Crank radius (metres).
    pub crank_radius: f64,
    /// Angular velocity of crank (rad/s).
    pub omega_crank: f64,
    /// Current crank angle (radians).
    pub crank_angle: f64,
}
impl GenevaMechanism {
    /// Create a Geneva mechanism.
    pub fn new(n_slots: u32, crank_radius: f64, omega_crank: f64) -> Self {
        Self {
            n_slots,
            crank_radius,
            omega_crank,
            crank_angle: 0.0,
        }
    }
    /// Centre distance between crank and wheel pivots (metres).
    ///
    /// `a = r / sin(π/n)`.
    pub fn centre_distance(&self) -> f64 {
        self.crank_radius / (PI / self.n_slots as f64).sin()
    }
    /// Half-angle of wheel rotation per index (radians) = π/n.
    pub fn wheel_half_angle(&self) -> f64 {
        PI / self.n_slots as f64
    }
    /// Half-angle of crank during which the wheel is driven (radians).
    ///
    /// `α = π/2 - π/n`.
    pub fn crank_driving_half_angle(&self) -> f64 {
        PI / 2.0 - PI / self.n_slots as f64
    }
    /// Fraction of crank revolution during which the wheel is stationary.
    pub fn dwell_fraction(&self) -> f64 {
        let driving_angle = 2.0 * self.crank_driving_half_angle();
        1.0 - driving_angle / (2.0 * PI)
    }
    /// Geneva wheel angular velocity at crank angle `phi` (rad/s).
    ///
    /// Valid only during the driving phase; returns 0 during dwell.
    pub fn wheel_velocity(&self, phi: f64) -> f64 {
        let alpha = self.crank_driving_half_angle();
        if phi.abs() > alpha {
            return 0.0;
        }
        let a = self.centre_distance();
        let r = self.crank_radius;
        let num = r * (a * phi.cos() - r);
        let den = a * a + r * r - 2.0 * a * r * phi.cos();
        self.omega_crank * num / den.max(1e-15)
    }
    /// Geneva wheel angle at crank angle `phi` (radians from start of driving phase).
    pub fn wheel_angle(&self, phi: f64) -> f64 {
        let a = self.centre_distance();
        let r = self.crank_radius;
        let alpha = self.crank_driving_half_angle();
        if phi.abs() > alpha {
            return self.wheel_half_angle() * (phi / phi.abs().max(1e-15));
        }
        (r * phi.sin() / (a - r * phi.cos()).max(1e-15)).atan()
    }
    /// Advance crank by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.crank_angle = (self.crank_angle + self.omega_crank * dt) % (2.0 * PI);
    }
}
/// Involute gear geometry parameters.
#[allow(dead_code)]
pub struct InvoluteGear {
    /// Number of teeth.
    pub teeth: u32,
    /// Module m (metres).
    pub module: f64,
    /// Pressure angle φ (radians).
    pub pressure_angle: f64,
    /// Addendum coefficient (typically 1.0).
    pub addendum_coeff: f64,
    /// Dedendum coefficient (typically 1.25).
    pub dedendum_coeff: f64,
    /// Profile shift coefficient x (0 = standard).
    pub profile_shift: f64,
}
impl InvoluteGear {
    /// Create a standard involute gear.
    pub fn new(teeth: u32, module: f64) -> Self {
        Self {
            teeth,
            module,
            pressure_angle: 20.0_f64.to_radians(),
            addendum_coeff: 1.0,
            dedendum_coeff: 1.25,
            profile_shift: 0.0,
        }
    }
    /// Pitch circle radius (metres).
    pub fn pitch_radius(&self) -> f64 {
        self.module * self.teeth as f64 / 2.0
    }
    /// Base circle radius (metres): r_b = r * cos(φ).
    pub fn base_radius(&self) -> f64 {
        self.pitch_radius() * self.pressure_angle.cos()
    }
    /// Addendum circle radius (metres): r_a = r + m*(1 + x).
    pub fn addendum_radius(&self) -> f64 {
        self.pitch_radius() + self.module * (self.addendum_coeff + self.profile_shift)
    }
    /// Dedendum circle radius (metres): r_d = r - m*(1.25 - x).
    pub fn dedendum_radius(&self) -> f64 {
        self.pitch_radius() - self.module * (self.dedendum_coeff - self.profile_shift)
    }
    /// Involute function: inv(φ) = tan(φ) - φ.
    pub fn involute_function(&self) -> f64 {
        self.pressure_angle.tan() - self.pressure_angle
    }
    /// Circular tooth thickness at pitch circle (metres).
    pub fn tooth_thickness(&self) -> f64 {
        self.module * (PI / 2.0 + 2.0 * self.profile_shift * self.pressure_angle.tan())
    }
    /// Circular pitch (metres): p = π * m.
    pub fn circular_pitch(&self) -> f64 {
        PI * self.module
    }
    /// Base pitch (metres): p_b = p * cos(φ).
    pub fn base_pitch(&self) -> f64 {
        self.circular_pitch() * self.pressure_angle.cos()
    }
}
/// Ratchet or detent locking mechanism.
#[allow(dead_code)]
pub enum LockingKind {
    /// One-way ratchet with pawl.
    Ratchet {
        /// Number of teeth on ratchet wheel.
        teeth: u32,
        /// Pawl spring force (N).
        pawl_spring_force: f64,
    },
    /// Spring-loaded detent.
    Detent {
        /// Detent ball radius (metres).
        ball_radius: f64,
        /// Spring preload force (N).
        spring_force: f64,
    },
}
/// Extended analysis results for a four-bar linkage at a given crank angle.
#[allow(dead_code)]
pub struct FourBarAnalysis {
    /// Crank angle (radians).
    pub theta2: f64,
    /// Coupler angle (radians).
    pub theta3: f64,
    /// Follower angle (radians).
    pub theta4: f64,
    /// Transmission angle (radians).
    pub transmission_angle: f64,
    /// Coupler point x coordinate (metres).
    pub coupler_x: f64,
    /// Coupler point y coordinate (metres).
    pub coupler_y: f64,
    /// Angular velocity of coupler link (rad/s).
    pub omega3: f64,
    /// Angular velocity of follower link (rad/s).
    pub omega4: f64,
}
impl FourBarAnalysis {
    /// Analyse a four-bar linkage at crank angle `theta2` with crank speed `omega2`.
    ///
    /// Returns `None` if the linkage cannot assemble at this angle.
    #[allow(dead_code)]
    pub fn analyse(l1: f64, l2: f64, l3: f64, l4: f64, theta2: f64, omega2: f64) -> Option<Self> {
        let (theta3, theta4) = four_bar_angles(l1, l2, l3, l4, theta2)?;
        let omega3 =
            -omega2 * l2 * (theta2 - theta4).sin() / (l3 * (theta3 - theta4).sin()).max(1e-12);
        let omega4 =
            -omega2 * l2 * (theta2 - theta3).sin() / (l4 * (theta4 - theta3).sin()).max(1e-12);
        let a_x = l2 * theta2.cos();
        let a_y = l2 * theta2.sin();
        let coupler_x = a_x + 0.5 * l3 * theta3.cos();
        let coupler_y = a_y + 0.5 * l3 * theta3.sin();
        let transmission_angle = (theta4 - theta3).abs();
        Some(Self {
            theta2,
            theta3,
            theta4,
            transmission_angle,
            coupler_x,
            coupler_y,
            omega3,
            omega4,
        })
    }
}
/// Configuration of an epicyclic (planetary) gear train.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpicyclicConfig {
    /// Standard planetary: sun input, ring fixed, carrier output.
    PlanetarySunIn,
    /// Star drive: carrier fixed, sun input, ring output.
    StarDrive,
    /// Solar drive: ring input, sun fixed, carrier output.
    SolarDrive,
    /// Differential: all three members free.
    Differential,
}
/// Epicyclic gear train with sun, planet, ring gears.
#[allow(dead_code)]
pub struct EpicyclicGearTrain {
    /// Number of teeth on sun gear.
    pub sun_teeth: u32,
    /// Number of teeth on each planet gear.
    pub planet_teeth: u32,
    /// Number of teeth on ring gear.
    pub ring_teeth: u32,
    /// Number of planet gears.
    pub n_planets: u32,
    /// Module (metres).
    pub module: f64,
    /// Configuration.
    pub config: EpicyclicConfig,
}
impl EpicyclicGearTrain {
    /// Create a standard planetary gear train satisfying `r = s + 2*p`.
    pub fn new(sun_teeth: u32, planet_teeth: u32, n_planets: u32, module: f64) -> Self {
        let ring_teeth = sun_teeth + 2 * planet_teeth;
        Self {
            sun_teeth,
            planet_teeth,
            ring_teeth,
            n_planets,
            module,
            config: EpicyclicConfig::PlanetarySunIn,
        }
    }
    /// Overall gear ratio using the Willis equation.
    ///
    /// Returns `ω_output / ω_input` for the current `config`.
    pub fn gear_ratio(&self) -> f64 {
        let r = self.ring_teeth as f64;
        let s = self.sun_teeth as f64;
        match self.config {
            EpicyclicConfig::PlanetarySunIn => 1.0 + r / s,
            EpicyclicConfig::StarDrive => -r / s,
            EpicyclicConfig::SolarDrive => r / (r + s),
            EpicyclicConfig::Differential => 0.0,
        }
    }
    /// Output speed (rad/s) given input speed.
    pub fn output_speed(&self, input_speed: f64) -> f64 {
        input_speed / self.gear_ratio().max(1e-15)
    }
    /// Sun gear pitch radius (metres).
    pub fn sun_radius(&self) -> f64 {
        self.module * self.sun_teeth as f64 / 2.0
    }
    /// Ring gear pitch radius (metres).
    pub fn ring_radius(&self) -> f64 {
        self.module * self.ring_teeth as f64 / 2.0
    }
    /// Planet gear pitch radius (metres).
    pub fn planet_radius(&self) -> f64 {
        self.module * self.planet_teeth as f64 / 2.0
    }
    /// Carrier radius (centre distance from main axis to planet centre) (metres).
    pub fn carrier_radius(&self) -> f64 {
        self.sun_radius() + self.planet_radius()
    }
    /// Check the assembly (equal spacing) condition: `n_planets` divides `sun + ring`.
    pub fn assembly_condition(&self) -> bool {
        (self.sun_teeth + self.ring_teeth).is_multiple_of(self.n_planets)
    }
    /// Check the neighbour condition: planet circles don't overlap.
    pub fn neighbour_condition(&self) -> bool {
        let half_angle = PI / self.n_planets as f64;
        let carrier_r = self.carrier_radius();
        carrier_r * half_angle.sin() > self.planet_radius()
    }
    /// Torque ratio: output torque / input torque (ideal, no friction).
    pub fn torque_ratio(&self) -> f64 {
        self.gear_ratio()
    }
}
