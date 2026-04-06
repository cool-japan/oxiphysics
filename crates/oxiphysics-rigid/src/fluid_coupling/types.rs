//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{add, cross, length, normalize, scale, sub};

/// Hull resistance model using the Froude number.
///
/// Separates resistance into frictional, residual (wave-making), and
/// air-resistance components.
#[allow(dead_code)]
pub struct HullResistance {
    /// Wetted surface area (m²).
    pub s_wet: f64,
    /// Waterplane area (m²).
    pub a_wp: f64,
    /// Ship length at waterline (m).
    pub lwl: f64,
    /// Hull block coefficient C_b (0–1).
    pub c_block: f64,
    /// Form factor (1 + k) for frictional resistance correction.
    pub form_factor: f64,
}
impl HullResistance {
    /// Create a hull resistance model.
    pub fn new(s_wet: f64, a_wp: f64, lwl: f64, c_block: f64, form_factor: f64) -> Self {
        Self {
            s_wet,
            a_wp,
            lwl,
            c_block,
            form_factor: form_factor.max(1.0),
        }
    }
    /// ITTC-57 friction coefficient C_F for model-ship correlation.
    ///
    /// C_F = 0.075 / (log10(Re) - 2)²
    pub fn cf_ittc57(re: f64) -> f64 {
        if re <= 0.0 {
            return 0.0;
        }
        let log_re = re.log10();
        if (log_re - 2.0).abs() < 1e-15 {
            return f64::MAX;
        }
        0.075 / (log_re - 2.0).powi(2)
    }
    /// Frictional resistance (N).
    ///
    /// R_F = 0.5 * ρ * V² * S_wet * (1+k) * C_F
    pub fn frictional_resistance(&self, rho: f64, speed: f64, re: f64) -> f64 {
        let cf = Self::cf_ittc57(re);
        0.5 * rho * speed * speed * self.s_wet * self.form_factor * cf
    }
    /// Froude number Fn = V / sqrt(g * L).
    pub fn froude_number(speed: f64, lwl: f64, g: f64) -> f64 {
        if g <= 0.0 || lwl <= 0.0 {
            return 0.0;
        }
        speed / (g * lwl).sqrt()
    }
    /// Simplified wave-making resistance based on Froude number.
    ///
    /// Uses the Michell/Savitsky-type approximation:
    /// R_W = displacement_weight * C_w  where C_w = exp(-0.5 / Fn²) for Fn < 0.6.
    pub fn wave_resistance(&self, displacement_n: f64, speed: f64, g: f64) -> f64 {
        let fn_val = Self::froude_number(speed, self.lwl, g);
        if fn_val < 1e-15 {
            return 0.0;
        }
        let cw = (-0.5 / (fn_val * fn_val)).exp();
        displacement_n * cw * self.c_block
    }
    /// Total resistance = frictional + wave.
    pub fn total_resistance(
        &self,
        rho: f64,
        speed: f64,
        re: f64,
        displacement_n: f64,
        g: f64,
    ) -> f64 {
        self.frictional_resistance(rho, speed, re) + self.wave_resistance(displacement_n, speed, g)
    }
}
/// Buoyancy force calculation for a sphere.
pub struct SphereBuoyancy;
impl SphereBuoyancy {
    /// Buoyancy force on a sphere (positive = upward).
    ///
    /// F = rho_fluid * g * submerged_fraction * volume
    pub fn buoyancy_force(rho_fluid: f64, g: f64, submerged_fraction: f64, volume: f64) -> f64 {
        rho_fluid * g * submerged_fraction.clamp(0.0, 1.0) * volume
    }
}
/// Added (virtual) mass calculations.
pub struct AddedMass;
impl AddedMass {
    /// Added mass of a sphere: m_a = 0.5·ρ·(4/3·π·r³).
    pub fn sphere_added_mass(radius: f64, fluid_density: f64) -> f64 {
        let volume = (4.0 / 3.0) * PI * radius.powi(3);
        0.5 * fluid_density * volume
    }
    /// Added mass of a cylinder: m_a = ρ·π·r²·L.
    pub fn cylinder_added_mass(radius: f64, length_val: f64, fluid_density: f64) -> f64 {
        fluid_density * PI * radius * radius * length_val
    }
    /// Added mass of an ellipsoid along the a-axis.
    ///
    /// For a prolate spheroid (a > b = c):
    /// m_a = alpha_0 * rho * (4/3)*pi*a*b*c
    /// where alpha_0 depends on the eccentricity.
    pub fn ellipsoid_added_mass_a(a: f64, b: f64, fluid_density: f64) -> f64 {
        let ratio = (b / a).min(1.0);
        let e_sq = 1.0 - ratio * ratio;
        let volume = (4.0 / 3.0) * PI * a * b * b;
        if e_sq < 1e-6 {
            return 0.5 * fluid_density * volume;
        }
        let e = e_sq.sqrt();
        let ln_term = ((1.0 + e) / (1.0 - e)).ln();
        let alpha = (2.0 * (1.0 - e_sq)) / (2.0 * e_sq - (1.0 - e_sq) * ln_term);
        let cm = alpha / (2.0 - alpha);
        cm.abs() * fluid_density * volume
    }
    /// Effective (total) mass: body mass + added mass.
    pub fn effective_mass(body_mass: f64, added_mass: f64) -> f64 {
        body_mass + added_mass
    }
    /// Effective (total) inertia: body inertia + added inertia.
    pub fn effective_inertia(body_inertia: f64, added_inertia: f64) -> f64 {
        body_inertia + added_inertia
    }
    /// Added mass force: F_am = -m_a * a_body.
    ///
    /// The added mass resists the body's acceleration through the fluid.
    pub fn added_mass_force(added_mass: f64, body_acceleration: [f64; 3]) -> [f64; 3] {
        scale(body_acceleration, -added_mass)
    }
}
/// Buoyancy model that accounts for body orientation.
///
/// When a body tilts, the submerged volume changes and the buoyancy force
/// shifts laterally (metacentric height effect). This model uses a capsule/box
/// approximation to compute orientation-dependent buoyancy.
#[allow(dead_code)]
pub struct OrientedBuoyancy {
    /// Fluid density (kg/m³).
    pub rho: f64,
    /// Gravitational acceleration (m/s²).
    pub g: f64,
    /// Body half-extents \[lx, ly, lz\] in body frame.
    pub half_extents: [f64; 3],
}
impl OrientedBuoyancy {
    /// Create an oriented buoyancy model for a box-shaped body.
    pub fn new(rho: f64, g: f64, half_extents: [f64; 3]) -> Self {
        Self {
            rho,
            g,
            half_extents,
        }
    }
    /// Full volume of the box.
    pub fn full_volume(&self) -> f64 {
        8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
    }
    /// Compute the buoyancy force vector (world space, pointing up) given
    /// the depth of the bottom of the body below the free surface.
    ///
    /// `depth_bottom` = distance from bottom of box to free surface (positive = submerged).
    /// `up` = world-space up direction (normalised).
    pub fn buoyancy_force(&self, depth_bottom: f64, up: [f64; 3]) -> [f64; 3] {
        let height = 2.0 * self.half_extents[1];
        let submerged = depth_bottom.clamp(0.0, height);
        let base_area = 4.0 * self.half_extents[0] * self.half_extents[2];
        let vol = base_area * submerged;
        let f_mag = self.rho * self.g * vol;
        [up[0] * f_mag, up[1] * f_mag, up[2] * f_mag]
    }
    /// Metacentric height approximation for a box-shaped hull.
    ///
    /// GM = I_waterplane / V_displaced - BG
    ///
    /// `draught` = draught of the hull (m), `bg` = buoyancy centre to gravity distance (m).
    pub fn metacentric_height(&self, draught: f64, bg: f64) -> f64 {
        let v = 4.0 * self.half_extents[0] * self.half_extents[2] * draught;
        if v < 1e-30 {
            return 0.0;
        }
        let i_wp = (2.0 * self.half_extents[0]) * (2.0 * self.half_extents[2]).powi(3) / 12.0;
        i_wp / v - bg
    }
}
/// 3×3 diagonal added mass tensor for common body shapes.
///
/// The added mass tensor relates body acceleration to the reaction force from
/// the surrounding fluid:  F_am = -\[m_a\] * a_body.
#[allow(dead_code)]
pub struct AddedMassTensor {
    /// Diagonal entries \[m_a_x, m_a_y, m_a_z\] (kg).
    pub diagonal: [f64; 3],
}
impl AddedMassTensor {
    /// Added mass tensor for a sphere (isotropic: all diagonal entries equal).
    ///
    /// m_a = 0.5 * ρ * (4/3 π r³) for each axis.
    pub fn sphere(radius: f64, fluid_density: f64) -> Self {
        let m_a = 0.5 * fluid_density * (4.0 / 3.0) * PI * radius.powi(3);
        Self { diagonal: [m_a; 3] }
    }
    /// Added mass tensor for a finite cylinder (axis along z).
    ///
    /// Transverse axes (x, y): m_a = ρ π r² L  (Lamb's result)
    /// Axial direction (z):    m_a ≈ 0 for a long slender cylinder.
    pub fn cylinder(radius: f64, length: f64, fluid_density: f64) -> Self {
        let m_transverse = fluid_density * PI * radius * radius * length;
        Self {
            diagonal: [m_transverse, m_transverse, 0.0],
        }
    }
    /// Added mass tensor for a thin flat plate (normal along z, dimensions a × b).
    ///
    /// Normal direction: m_a = (π/4) ρ (min(a,b))² * max(a,b)   (approximate)
    /// Tangential directions: negligible (≈ 0).
    pub fn flat_plate(a: f64, b: f64, fluid_density: f64) -> Self {
        let minor = a.min(b);
        let major = a.max(b);
        let m_normal = (PI / 4.0) * fluid_density * minor * minor * major;
        Self {
            diagonal: [0.0, 0.0, m_normal],
        }
    }
    /// Compute the added mass force: F = -\[m_a\] * acceleration.
    pub fn force(&self, acceleration: [f64; 3]) -> [f64; 3] {
        [
            -self.diagonal[0] * acceleration[0],
            -self.diagonal[1] * acceleration[1],
            -self.diagonal[2] * acceleration[2],
        ]
    }
    /// Effective mass along axis `i` (0=x,1=y,2=z): m_eff = m_body + m_a_i.
    pub fn effective_mass(&self, body_mass: f64, axis: usize) -> f64 {
        body_mass + self.diagonal[axis.min(2)]
    }
}
/// Simple vortex-induced vibration parameters based on Strouhal number.
pub struct VivParams {
    /// Strouhal number (dimensionless).
    pub st: f64,
    /// Fluid density (kg/m³).
    pub rho: f64,
    /// Cylinder diameter (m).
    pub d: f64,
}
impl VivParams {
    /// Vortex shedding frequency: f_vs = St * U / D (Hz).
    pub fn vortex_shedding_frequency(&self, u: f64) -> f64 {
        if self.d <= 0.0 {
            return 0.0;
        }
        self.st * u / self.d
    }
    /// Check if vortex lock-in is occurring.
    ///
    /// Lock-in occurs when the shedding frequency is within ±20% of the
    /// structural natural frequency.
    pub fn lock_in_range(&self, fn_hz: f64, u: f64) -> bool {
        let f_vs = self.vortex_shedding_frequency(u);
        if fn_hz <= 0.0 {
            return false;
        }
        let ratio = f_vs / fn_hz;
        (0.8..=1.2).contains(&ratio)
    }
}
/// Froude-Krylov and diffraction wave force on a submerged body.
///
/// The total first-order wave force on a body is the sum of
/// the Froude-Krylov force (undisturbed pressure field) and
/// the diffraction force (scattered wave field).
///
/// For small bodies (ka ≪ 1, k = wave number, a = body radius) the
/// diffraction force is negligible and only Froude-Krylov matters.
#[allow(dead_code)]
pub struct FroudeKrylovForce {
    /// Body volume (m³).
    pub volume: f64,
    /// Added mass coefficient (dimensionless). Typically 1.0 for small bodies.
    pub cm: f64,
}
impl FroudeKrylovForce {
    /// Create a new Froude-Krylov force model.
    pub fn new(volume: f64, cm: f64) -> Self {
        Self { volume, cm }
    }
    /// Froude-Krylov force: F_FK = ρ_f * V * a_fluid.
    ///
    /// This is the force due to the undisturbed pressure gradient acting on
    /// the volume occupied by the body.
    pub fn froude_krylov(&self, fluid_density: f64, fluid_acceleration: [f64; 3]) -> [f64; 3] {
        scale(fluid_acceleration, fluid_density * self.volume)
    }
    /// Total wave excitation force including diffraction:
    /// F = ρ_f * V * (1 + C_m) * a_fluid  (for small ka).
    ///
    /// C_m = 1 gives the classic Morison inertia coefficient of 2.
    pub fn total_excitation(&self, fluid_density: f64, fluid_acceleration: [f64; 3]) -> [f64; 3] {
        scale(
            fluid_acceleration,
            fluid_density * self.volume * (1.0 + self.cm),
        )
    }
    /// Diffraction force (scattered wave contribution):
    /// F_diff = ρ_f * V * C_m * a_fluid.
    pub fn diffraction(&self, fluid_density: f64, fluid_acceleration: [f64; 3]) -> [f64; 3] {
        scale(fluid_acceleration, fluid_density * self.volume * self.cm)
    }
    /// Radiation force on a body oscillating in calm water.
    ///
    /// F_rad = -m_a * a_body - b_rad * v_body
    /// where `m_a` is the added mass and `b_rad` is the radiation damping coefficient.
    pub fn radiation_force(
        added_mass: f64,
        radiation_damping: f64,
        body_velocity: [f64; 3],
        body_acceleration: [f64; 3],
    ) -> [f64; 3] {
        let f_added_mass = scale(body_acceleration, -added_mass);
        let f_damping = scale(body_velocity, -radiation_damping);
        add(f_added_mass, f_damping)
    }
}
/// Combined fluid forces on a fully submerged rigid body.
pub struct SubmergedBodyCoupling {
    /// Density of the body (kg/m³).
    pub body_density: f64,
    /// Volume of the body (m³).
    pub body_volume: f64,
    /// Surrounding fluid properties.
    pub fluid: FluidProperties,
    /// Drag coefficient.
    pub cd: f64,
}
impl SubmergedBodyCoupling {
    /// Net force = gravity + buoyancy + drag.
    pub fn net_force(
        &self,
        _position: [f64; 3],
        velocity: [f64; 3],
        gravity: [f64; 3],
    ) -> [f64; 3] {
        let body_mass = self.body_density * self.body_volume;
        let weight = scale(gravity, body_mass);
        let buoyancy = BuoyancyForce::archimedes(self.fluid.density, self.body_volume, gravity);
        let radius = (self.body_volume * 3.0 / (4.0 * PI)).cbrt();
        let cross_section = PI * radius * radius;
        let drag =
            HydrodynamicDrag::drag_force(velocity, cross_section, self.cd, self.fluid.density);
        add(add(weight, buoyancy), drag)
    }
    /// Net force including added mass effect.
    pub fn net_force_with_added_mass(
        &self,
        _position: [f64; 3],
        velocity: [f64; 3],
        acceleration: [f64; 3],
        gravity: [f64; 3],
    ) -> [f64; 3] {
        let base = self.net_force(_position, velocity, gravity);
        let radius = (self.body_volume * 3.0 / (4.0 * PI)).cbrt();
        let m_a = AddedMass::sphere_added_mass(radius, self.fluid.density);
        let f_am = AddedMass::added_mass_force(m_a, acceleration);
        add(base, f_am)
    }
    /// Terminal velocity magnitude where weight = buoyancy + drag.
    pub fn terminal_velocity(&self, gravity: [f64; 3], cross_section: f64) -> f64 {
        let g = length(gravity);
        let body_mass = self.body_density * self.body_volume;
        let buoyancy_mag = self.fluid.density * self.body_volume * g;
        let net_weight = body_mass * g - buoyancy_mag;
        if net_weight <= 0.0 {
            return 0.0;
        }
        (2.0 * net_weight / (self.fluid.density * self.cd * cross_section)).sqrt()
    }
    /// Whether the body will float (body density < fluid density).
    pub fn will_float(&self) -> bool {
        self.body_density < self.fluid.density
    }
    /// Equilibrium submerged fraction for a floating body.
    pub fn equilibrium_submerged_fraction(&self) -> f64 {
        if self.fluid.density <= 0.0 {
            return 1.0;
        }
        (self.body_density / self.fluid.density).min(1.0)
    }
}
/// Buoyancy force calculations.
pub struct BuoyancyForce;
impl BuoyancyForce {
    /// Archimedes buoyancy force: F = -ρ_f · V · g.
    pub fn archimedes(fluid_density: f64, displaced_volume: f64, gravity: [f64; 3]) -> [f64; 3] {
        scale(gravity, -fluid_density * displaced_volume)
    }
    /// Partial submersion: scales displaced volume by submerged fraction.
    pub fn partial_submersion(
        fluid_density: f64,
        submerged_volume: f64,
        total_volume: f64,
        gravity: [f64; 3],
    ) -> [f64; 3] {
        let fraction = if total_volume > 0.0 {
            submerged_volume / total_volume
        } else {
            0.0
        };
        Self::archimedes(fluid_density, submerged_volume * fraction, gravity)
    }
    /// Buoyancy torque: τ = (cob − com) × F.
    pub fn buoyancy_torque(
        center_of_buoyancy: [f64; 3],
        center_of_mass: [f64; 3],
        buoyancy_force: [f64; 3],
    ) -> [f64; 3] {
        cross(sub(center_of_buoyancy, center_of_mass), buoyancy_force)
    }
    /// Buoyancy force on a sphere partially submerged to depth `h` below the
    /// water surface (sphere of radius `r`).
    ///
    /// Displaced volume of a spherical cap: V = π*h²*(3r - h)/3.
    pub fn sphere_partial_buoyancy(
        radius: f64,
        submerged_depth: f64,
        fluid_density: f64,
        gravity: [f64; 3],
    ) -> [f64; 3] {
        let h = submerged_depth.clamp(0.0, 2.0 * radius);
        let vol = PI * h * h * (3.0 * radius - h) / 3.0;
        Self::archimedes(fluid_density, vol, gravity)
    }
    /// Buoyancy force on a horizontal cylinder partially submerged.
    ///
    /// The submerged cross-section area is computed from the chord geometry.
    pub fn cylinder_partial_buoyancy(
        radius: f64,
        submerged_depth: f64,
        cyl_length: f64,
        fluid_density: f64,
        gravity: [f64; 3],
    ) -> [f64; 3] {
        let h = submerged_depth.clamp(0.0, 2.0 * radius);
        let r = radius;
        let area = r * r * ((r - h + r).max(0.0) / r).min(1.0).acos()
            - (r - h) * (2.0 * r * h - h * h).max(0.0).sqrt();
        let vol = area.abs() * cyl_length;
        Self::archimedes(fluid_density, vol, gravity)
    }
}
/// A Morison equation element (cylinder segment) for wave force computation.
pub struct MorisonElement {
    /// Diameter of the element (m).
    pub diameter: f64,
    /// Length of the element (m).
    pub length: f64,
    /// Inertia coefficient (dimensionless, typically 2.0).
    pub cm: f64,
    /// Drag coefficient (dimensionless, typically 1.0).
    pub cd: f64,
}
impl MorisonElement {
    /// Morison force on the element.
    ///
    /// F = ρ·V·Cm·a_fluid + 0.5·ρ·Cd·A·|v_rel|·v_rel
    ///
    /// where v_rel = fluid_vel - body_vel.
    pub fn morison_force(
        &self,
        fluid_accel: [f64; 3],
        fluid_vel: [f64; 3],
        body_vel: [f64; 3],
        rho: f64,
    ) -> [f64; 3] {
        let cross_section = 0.25 * PI * self.diameter * self.diameter;
        let volume = cross_section * self.length;
        let v_rel = sub(fluid_vel, body_vel);
        let speed = length(v_rel);
        let drag = if speed < 1e-300 {
            [0.0; 3]
        } else {
            scale(v_rel, 0.5 * rho * self.cd * cross_section * speed)
        };
        let inertia = scale(fluid_accel, rho * volume * self.cm);
        add(drag, inertia)
    }
}
/// Pendulum analogy model for liquid sloshing in a rectangular tank.
///
/// The sloshing is modelled as a simple pendulum whose effective mass and
/// length are derived from the fluid fill level.  See Abramson (1966).
#[allow(dead_code)]
pub struct SloshingModel {
    /// Tank half-length in the sloshing direction (m).
    pub half_length: f64,
    /// Fluid fill height (m).
    pub fill_height: f64,
    /// Total fluid mass (kg).
    pub fluid_mass: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: f64,
    /// Current pendulum angle (rad).
    pub angle: f64,
    /// Current pendulum angular velocity (rad/s).
    pub angular_vel: f64,
}
impl SloshingModel {
    /// Create a sloshing model for a rectangular tank.
    pub fn new(half_length: f64, fill_height: f64, fluid_mass: f64, gravity: f64) -> Self {
        Self {
            half_length,
            fill_height,
            fluid_mass,
            gravity,
            angle: 0.0,
            angular_vel: 0.0,
        }
    }
    /// Effective pendulum length for the first sloshing mode.
    ///
    /// L_eff = h / tanh(π h / (2 a))  where h = fill height, a = half-length.
    ///
    /// For deep fill (h ≫ a): L_eff → h.
    pub fn effective_pendulum_length(&self) -> f64 {
        if self.half_length <= 0.0 {
            return self.fill_height;
        }
        let arg = PI * self.fill_height / (2.0 * self.half_length);
        if arg < 1e-12 {
            return self.fill_height;
        }
        self.fill_height / arg.tanh()
    }
    /// Natural sloshing frequency (rad/s) for the first mode.
    ///
    /// ω_n = sqrt(g / L_eff)
    pub fn natural_frequency_rad(&self) -> f64 {
        let l_eff = self.effective_pendulum_length();
        if l_eff <= 0.0 {
            return 0.0;
        }
        (self.gravity / l_eff).sqrt()
    }
    /// Natural sloshing period (s).
    pub fn natural_period(&self) -> f64 {
        let omega = self.natural_frequency_rad();
        if omega <= 0.0 {
            return f64::INFINITY;
        }
        2.0 * PI / omega
    }
    /// Effective sloshing mass (fraction of total fluid that participates).
    ///
    /// m_eff = m_fluid * tanh(π h / (2 a)) * 8 / π²  (first mode, approx.)
    pub fn effective_mass(&self) -> f64 {
        if self.half_length <= 0.0 {
            return self.fluid_mass;
        }
        let arg = PI * self.fill_height / (2.0 * self.half_length);
        let frac = 8.0 / (PI * PI) * arg.tanh().min(1.0);
        self.fluid_mass * frac
    }
    /// Step the pendulum analogy using a semi-implicit Euler scheme.
    ///
    /// `tank_accel_x` is the horizontal tank acceleration (m/s²).
    pub fn step(&mut self, tank_accel_x: f64, dt: f64) {
        let omega_n = self.natural_frequency_rad();
        let l_eff = self.effective_pendulum_length().max(1e-10);
        let alpha = -omega_n * omega_n * self.angle - tank_accel_x / l_eff;
        self.angular_vel += alpha * dt;
        self.angle += self.angular_vel * dt;
    }
    /// Sloshing force exerted on the tank walls (reaction force from the fluid).
    ///
    /// F_slosh = m_eff * g * θ  (linearised pendulum).
    pub fn sloshing_force(&self) -> f64 {
        self.effective_mass() * self.gravity * self.angle
    }
}
/// Basic vortex-induced vibration model.
///
/// Predicts whether VIV lock-in occurs and estimates the oscillation amplitude.
pub struct VortexInducedVibration {
    /// Diameter of the cylinder (m).
    pub diameter: f64,
    /// Natural frequency of the structure (Hz).
    pub fn_structure: f64,
    /// Structural damping ratio (dimensionless).
    pub damping_ratio: f64,
    /// Mass per unit length of the cylinder (kg/m).
    pub mass_per_length: f64,
}
impl VortexInducedVibration {
    /// Create a new VIV model.
    pub fn new(diameter: f64, fn_structure: f64, damping_ratio: f64, mass_per_length: f64) -> Self {
        Self {
            diameter,
            fn_structure,
            damping_ratio,
            mass_per_length,
        }
    }
    /// Strouhal number for a circular cylinder (typically ~0.2).
    pub fn strouhal_number() -> f64 {
        0.20
    }
    /// Vortex shedding frequency at the given flow speed (Hz).
    pub fn shedding_frequency(&self, flow_speed: f64) -> f64 {
        if self.diameter <= 0.0 {
            return 0.0;
        }
        Self::strouhal_number() * flow_speed / self.diameter
    }
    /// Reduced velocity: V_r = U / (f_n * D).
    pub fn reduced_velocity(&self, flow_speed: f64) -> f64 {
        if self.fn_structure <= 0.0 || self.diameter <= 0.0 {
            return 0.0;
        }
        flow_speed / (self.fn_structure * self.diameter)
    }
    /// Whether lock-in is expected (shedding frequency near the natural frequency).
    ///
    /// Lock-in typically occurs when 0.8 < f_s/f_n < 1.2 (i.e. V_r ≈ 4..7).
    pub fn is_lock_in(&self, flow_speed: f64) -> bool {
        let vr = self.reduced_velocity(flow_speed);
        (4.0..=7.0).contains(&vr)
    }
    /// Estimated peak amplitude ratio A/D during lock-in.
    ///
    /// Uses the Skop-Griffin correlation:
    /// A/D ≈ 1.29 / (1 + 0.43 * (2π * S_G)^0.35)^(0.85)
    /// where S_G = 2π * zeta * m_ratio is the Skop-Griffin parameter,
    /// and m_ratio = m / (rho * D^2).
    pub fn peak_amplitude_ratio(&self, fluid_density: f64) -> f64 {
        let m_ratio = if fluid_density > 0.0 && self.diameter > 0.0 {
            self.mass_per_length / (fluid_density * self.diameter * self.diameter)
        } else {
            1.0
        };
        let sg = 2.0 * PI * self.damping_ratio * m_ratio;
        let term = 2.0 * PI * sg;
        1.29 / (1.0 + 0.43 * term.powf(0.35)).powf(0.85)
    }
    /// Lift force per unit length during lock-in.
    ///
    /// F_L = 0.5 * rho * U^2 * D * C_L * sin(2π * f_s * t)
    pub fn lift_force_per_length(
        &self,
        flow_speed: f64,
        fluid_density: f64,
        cl: f64,
        time: f64,
    ) -> f64 {
        let f_s = self.shedding_frequency(flow_speed);
        let dynamic_pressure = 0.5 * fluid_density * flow_speed * flow_speed;
        dynamic_pressure * self.diameter * cl * (2.0 * PI * f_s * time).sin()
    }
}
/// A simplified floating body for stability analysis.
pub struct FloatingBody {
    /// Body mass (kg).
    pub mass: f64,
    /// Total displaced volume (m³).
    pub volume: f64,
    /// Center of mass in world frame.
    pub center_of_mass: [f64; 3],
    /// Geometric center of displaced volume in world frame.
    pub center_of_buoyancy: [f64; 3],
    /// Metacentric height GM (m). Positive means stable.
    pub metacentric_height: f64,
    /// Position of body origin.
    pub position: [f64; 3],
    /// Roll angle (simplified 2-D, radians).
    pub orientation_angle: f64,
    /// Linear velocity (m/s).
    pub linear_vel: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_vel: f64,
}
impl FloatingBody {
    /// Stability check: stable when GM > 0.
    pub fn is_stable(&self) -> bool {
        self.metacentric_height > 0.0
    }
    /// Compute GM = I_water/V − (z_com − z_cob).
    pub fn metacenter_height(&self, second_moment_waterplane: f64) -> f64 {
        let kb = self.center_of_buoyancy[2];
        let kg = self.center_of_mass[2];
        second_moment_waterplane / self.volume - (kg - kb)
    }
    /// Righting moment: τ = mass·g·GM·sin(θ).
    pub fn righting_moment(&self) -> f64 {
        const G: f64 = 9.81;
        self.mass * G * self.metacentric_height * self.orientation_angle.sin()
    }
    /// Natural roll period (s) for small oscillations.
    ///
    /// T = 2π * sqrt(I_xx / (mass * g * GM))
    ///
    /// `i_xx` is the roll moment of inertia (kg·m²).
    pub fn roll_period(&self, i_xx: f64) -> f64 {
        const G: f64 = 9.81;
        let denom = self.mass * G * self.metacentric_height;
        if denom <= 0.0 {
            return f64::INFINITY;
        }
        2.0 * PI * (i_xx / denom).sqrt()
    }
    /// Step the simplified 2-D roll dynamics (Euler).
    ///
    /// Uses the single-DOF equation: I_xx * theta_ddot = -m*g*GM*sin(theta)
    /// - c * omega, where c is a linear damping coefficient.
    pub fn step_roll(&mut self, i_xx: f64, damping: f64, dt: f64) {
        const G: f64 = 9.81;
        if i_xx <= 0.0 {
            return;
        }
        let restoring = -self.mass * G * self.metacentric_height * self.orientation_angle.sin();
        let damping_torque = -damping * self.angular_vel;
        let alpha = (restoring + damping_torque) / i_xx;
        self.angular_vel += alpha * dt;
        self.orientation_angle += self.angular_vel * dt;
    }
}
/// Propeller thrust model based on actuator-disk (momentum) theory.
///
/// The thrust of a propeller in open water is:
///
/// ```text
/// T = K_T(J) · ρ · n² · D⁴
/// J = V_a / (n · D)   (advance ratio)
/// ```
///
/// For the simplified model we use a linear K_T(J) fit:
/// K_T(J) = k_t0 - k_t1 * J
#[allow(dead_code)]
pub struct PropellerThrust {
    /// Propeller diameter (m).
    pub diameter: f64,
    /// Zero-advance thrust coefficient K_T0 (dimensionless).
    pub kt0: f64,
    /// Slope of K_T vs. J (dimensionless).
    pub kt1: f64,
    /// Torque coefficient K_Q0 (dimensionless).
    pub kq0: f64,
    /// Slope of K_Q vs. J.
    pub kq1: f64,
}
impl PropellerThrust {
    /// Create a propeller with given geometric and hydrodynamic coefficients.
    pub fn new(diameter: f64, kt0: f64, kt1: f64, kq0: f64, kq1: f64) -> Self {
        Self {
            diameter,
            kt0,
            kt1,
            kq0,
            kq1,
        }
    }
    /// Advance ratio J = V_a / (n * D).
    ///
    /// `va` = advance velocity (m/s), `n` = rotational speed (rev/s).
    pub fn advance_ratio(&self, va: f64, n: f64) -> f64 {
        if n.abs() < 1e-15 || self.diameter < 1e-15 {
            return 0.0;
        }
        va / (n * self.diameter)
    }
    /// Thrust coefficient K_T(J).
    pub fn kt(&self, j: f64) -> f64 {
        (self.kt0 - self.kt1 * j).max(0.0)
    }
    /// Torque coefficient K_Q(J).
    pub fn kq(&self, j: f64) -> f64 {
        (self.kq0 - self.kq1 * j).max(0.0)
    }
    /// Thrust force T = K_T * ρ * n² * D⁴ (N).
    ///
    /// `rho` = fluid density (kg/m³), `n` = rev/s.
    pub fn thrust(&self, rho: f64, n: f64, va: f64) -> f64 {
        let j = self.advance_ratio(va, n);
        let kt = self.kt(j);
        kt * rho * n * n * self.diameter.powi(4)
    }
    /// Propeller torque Q = K_Q * ρ * n² * D⁵ (N·m).
    pub fn torque(&self, rho: f64, n: f64, va: f64) -> f64 {
        let j = self.advance_ratio(va, n);
        let kq = self.kq(j);
        kq * rho * n * n * self.diameter.powi(5)
    }
    /// Open-water efficiency η = J / (2π) * K_T / K_Q.
    pub fn open_water_efficiency(&self, va: f64, n: f64) -> f64 {
        let j = self.advance_ratio(va, n);
        let kq = self.kq(j);
        if kq < 1e-30 {
            return 0.0;
        }
        let kt = self.kt(j);
        j * kt / (2.0 * PI * kq)
    }
}
/// Unified interface for computing common hydrodynamic forces on rigid bodies.
///
/// This struct consolidates Morison equation forces, added-mass acceleration
/// corrections, and vortex-induced vibration (VIV) lock-in forces into a
/// single API.
pub struct FluidCoupling {
    /// Fluid density (kg/m³).
    pub fluid_density: f64,
    /// Gravitational acceleration magnitude (m/s²).
    pub gravity: f64,
}
impl FluidCoupling {
    /// Create a new `FluidCoupling` with the specified fluid and gravity.
    pub fn new(fluid_density: f64, gravity: f64) -> Self {
        Self {
            fluid_density,
            gravity,
        }
    }
    /// Compute the Morison equation force on a circular cylinder element.
    ///
    /// The Morison force combines an inertia (added-mass) term and a drag
    /// term:
    ///
    /// ```text
    /// F = ρ · V · Cm · a_fluid  +  0.5 · ρ · Cd · A · |v_rel| · v_rel
    /// ```
    ///
    /// where `v_rel = fluid_vel − body_vel`.
    ///
    /// - `diameter`    : cylinder diameter (m)
    /// - `length`      : cylinder length (m)
    /// - `cm`          : inertia coefficient (dimensionless, typically 2.0)
    /// - `cd`          : drag coefficient (dimensionless, typically 1.0)
    /// - `fluid_accel` : fluid acceleration vector (m/s²)
    /// - `fluid_vel`   : fluid velocity vector (m/s)
    /// - `body_vel`    : body velocity vector (m/s)
    #[allow(clippy::too_many_arguments)]
    pub fn compute_morrison_force(
        &self,
        diameter: f64,
        elem_length: f64,
        cm: f64,
        cd: f64,
        fluid_accel: [f64; 3],
        fluid_vel: [f64; 3],
        body_vel: [f64; 3],
    ) -> [f64; 3] {
        let cross_section = 0.25 * PI * diameter * diameter;
        let volume = cross_section * elem_length;
        let rho = self.fluid_density;
        let v_rel = sub(fluid_vel, body_vel);
        let speed = length(v_rel);
        let drag = if speed < 1e-300 {
            [0.0; 3]
        } else {
            scale(v_rel, 0.5 * rho * cd * cross_section * speed)
        };
        let inertia = scale(fluid_accel, rho * volume * cm);
        add(drag, inertia)
    }
    /// Compute the added-mass force component on a sphere undergoing
    /// acceleration in a fluid.
    ///
    /// Added-mass force: F_a = −m_a · (a_body − a_fluid)
    ///
    /// where m_a = C_a · ρ_f · V is the added mass, and the net force
    /// opposes the body's acceleration relative to the fluid.
    ///
    /// - `radius`     : sphere radius (m)
    /// - `ca`         : added-mass coefficient (dimensionless, 0.5 for a sphere)
    /// - `body_accel` : body acceleration vector (m/s²)
    /// - `fluid_accel`: fluid acceleration vector (m/s²)
    pub fn compute_added_mass_force(
        &self,
        radius: f64,
        ca: f64,
        body_accel: [f64; 3],
        fluid_accel: [f64; 3],
    ) -> [f64; 3] {
        let volume = 4.0 / 3.0 * PI * radius * radius * radius;
        let added_mass = ca * self.fluid_density * volume;
        let rel_accel = sub(body_accel, fluid_accel);
        scale(rel_accel, -added_mass)
    }
    /// Compute the vortex-induced vibration (VIV) lock-in force per unit
    /// length on a cylinder.
    ///
    /// When the vortex shedding frequency is close to the structural natural
    /// frequency (lock-in condition), the oscillatory lift force is:
    ///
    /// ```text
    /// F_VIV = 0.5 · ρ · U² · D · C_L · sin(2π · f_s · t)
    /// ```
    ///
    /// Returns the force vector perpendicular to the flow in the given
    /// `lift_dir` direction; returns `[0; 3]` when not in lock-in.
    ///
    /// - `diameter`      : cylinder diameter (m)
    /// - `fn_structure`  : structural natural frequency (Hz)
    /// - `flow_speed`    : flow velocity magnitude (m/s)
    /// - `cl`            : lift coefficient (dimensionless)
    /// - `time`          : current simulation time (s)
    /// - `lift_dir`      : unit vector perpendicular to flow (lift direction)
    pub fn compute_vortex_induced_vibration(
        &self,
        diameter: f64,
        fn_structure: f64,
        flow_speed: f64,
        cl: f64,
        time: f64,
        lift_dir: [f64; 3],
    ) -> [f64; 3] {
        if diameter <= 0.0 || fn_structure <= 0.0 || flow_speed < 1e-12 {
            return [0.0; 3];
        }
        let st = 0.20_f64;
        let f_shedding = st * flow_speed / diameter;
        let ratio = f_shedding / fn_structure;
        if !(0.8..=1.2).contains(&ratio) {
            return [0.0; 3];
        }
        let dynamic_pressure = 0.5 * self.fluid_density * flow_speed * flow_speed;
        let force_per_length =
            dynamic_pressure * diameter * cl * (2.0 * PI * f_shedding * time).sin();
        scale(lift_dir, force_per_length)
    }
}
/// Flow regime classification based on the Froude number.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// Sub-critical flow: Fr < 1.
    SubCritical,
    /// Critical flow: Fr ≈ 1.
    Critical,
    /// Super-critical flow: Fr > 1.
    SuperCritical,
}
/// Reynolds number and related non-dimensional parameters.
#[allow(dead_code)]
pub struct ReynoldsUtils;
impl ReynoldsUtils {
    /// Reynolds number: Re = ρ * U * L / μ.
    pub fn reynolds(rho: f64, u: f64, l: f64, mu: f64) -> f64 {
        if mu <= 0.0 {
            return f64::INFINITY;
        }
        rho * u * l / mu
    }
    /// Weber number: We = ρ * U² * L / σ  (surface tension ratio).
    pub fn weber(rho: f64, u: f64, l: f64, surface_tension: f64) -> f64 {
        if surface_tension <= 0.0 {
            return f64::INFINITY;
        }
        rho * u * u * l / surface_tension
    }
    /// Euler number: Eu = ΔP / (0.5 ρ U²).
    pub fn euler(delta_p: f64, rho: f64, u: f64) -> f64 {
        let dyn_p = 0.5 * rho * u * u;
        if dyn_p <= 0.0 {
            return f64::INFINITY;
        }
        delta_p / dyn_p
    }
    /// Strouhal number: St = f * L / U.
    pub fn strouhal(freq: f64, l: f64, u: f64) -> f64 {
        if u <= 0.0 {
            return f64::INFINITY;
        }
        freq * l / u
    }
    /// Keulegan-Carpenter number: KC = U_m * T / D.
    ///
    /// Important for wave-structure interaction (oscillatory flow).
    pub fn keulegan_carpenter(u_max: f64, period: f64, diameter: f64) -> f64 {
        if diameter <= 0.0 {
            return 0.0;
        }
        u_max * period / diameter
    }
    /// Classify Reynolds number into laminar / transitional / turbulent regime.
    ///
    /// For pipe flow: Re < 2300 → laminar, 2300..4000 → transitional, > 4000 → turbulent.
    pub fn pipe_regime(re: f64) -> &'static str {
        if re < 2300.0 {
            "laminar"
        } else if re < 4000.0 {
            "transitional"
        } else {
            "turbulent"
        }
    }
}
/// Physical properties of a fluid.
pub struct FluidProperties {
    /// Fluid density in kg/m³.
    pub density: f64,
    /// Dynamic viscosity in Pa·s.
    pub dynamic_viscosity: f64,
}
impl FluidProperties {
    /// Kinematic viscosity ν = μ/ρ (m²/s).
    pub fn kinematic_viscosity(&self) -> f64 {
        self.dynamic_viscosity / self.density
    }
    /// Create water at 20 °C.
    pub fn water() -> Self {
        Self {
            density: 998.2,
            dynamic_viscosity: 1.002e-3,
        }
    }
    /// Create air at sea level, 20 °C.
    pub fn air() -> Self {
        Self {
            density: 1.225,
            dynamic_viscosity: 1.81e-5,
        }
    }
}
/// Hydrodynamic drag calculations.
pub struct HydrodynamicDrag;
impl HydrodynamicDrag {
    /// Quadratic drag force: F = -0.5·ρ·|v|²·Cd·A·v̂.
    pub fn drag_force(
        relative_vel: [f64; 3],
        cross_section: f64,
        cd: f64,
        fluid_density: f64,
    ) -> [f64; 3] {
        let speed = length(relative_vel);
        if speed < 1e-300 {
            return [0.0; 3];
        }
        let v_hat = normalize(relative_vel);
        let magnitude = 0.5 * fluid_density * speed * speed * cd * cross_section;
        scale(v_hat, -magnitude)
    }
    /// Stokes drag (low Reynolds number): F = -6·π·μ·r·v.
    pub fn stokes_drag(radius: f64, relative_vel: [f64; 3], fluid: &FluidProperties) -> [f64; 3] {
        scale(relative_vel, -6.0 * PI * fluid.dynamic_viscosity * radius)
    }
    /// Oseen drag correction for moderate Re: F = -6πμrv(1 + 3Re/16).
    pub fn oseen_drag(radius: f64, relative_vel: [f64; 3], fluid: &FluidProperties) -> [f64; 3] {
        let speed = length(relative_vel);
        let re = Self::reynolds_number(speed, 2.0 * radius, fluid);
        let correction = 1.0 + 3.0 * re / 16.0;
        scale(
            relative_vel,
            -6.0 * PI * fluid.dynamic_viscosity * radius * correction,
        )
    }
    /// Reynolds number: Re = v·L/ν.
    pub fn reynolds_number(velocity: f64, length_scale: f64, fluid: &FluidProperties) -> f64 {
        velocity * length_scale / fluid.kinematic_viscosity()
    }
    /// Drag coefficient for a sphere (Schlichting approximation).
    pub fn drag_coefficient_sphere(re: f64) -> f64 {
        if re < 1e-10 {
            return f64::INFINITY;
        }
        24.0 / re + 6.0 / (1.0 + re.sqrt()) + 0.4
    }
    /// Drag coefficient for a long cylinder (cross-flow).
    ///
    /// Empirical fit: Cd ≈ 1.0 + 10/Re^(2/3) for Re > 1.
    pub fn drag_coefficient_cylinder(re: f64) -> f64 {
        if re < 1e-10 {
            return f64::INFINITY;
        }
        if re < 1.0 {
            8.0 * PI / (re * (2.002 - (re / 4.0).max(0.001).ln()))
        } else {
            1.0 + 10.0 / re.powf(2.0 / 3.0)
        }
    }
    /// Drag coefficient for a flat plate normal to flow.
    pub fn drag_coefficient_flat_plate() -> f64 {
        1.98
    }
    /// Drag force on a sphere using the Schlichting Cd model.
    pub fn sphere_drag(radius: f64, relative_vel: [f64; 3], fluid: &FluidProperties) -> [f64; 3] {
        let speed = length(relative_vel);
        let re = Self::reynolds_number(speed, 2.0 * radius, fluid);
        let cd = Self::drag_coefficient_sphere(re);
        let area = PI * radius * radius;
        Self::drag_force(relative_vel, area, cd, fluid.density)
    }
    /// Drag force on a cylinder in cross-flow.
    pub fn cylinder_drag(
        radius: f64,
        cyl_length: f64,
        relative_vel: [f64; 3],
        fluid: &FluidProperties,
    ) -> [f64; 3] {
        let speed = length(relative_vel);
        let re = Self::reynolds_number(speed, 2.0 * radius, fluid);
        let cd = Self::drag_coefficient_cylinder(re);
        let area = 2.0 * radius * cyl_length;
        Self::drag_force(relative_vel, area, cd, fluid.density)
    }
    /// Drag force with an arbitrary Cd.
    pub fn arbitrary_drag(
        relative_vel: [f64; 3],
        reference_area: f64,
        cd: f64,
        fluid: &FluidProperties,
    ) -> [f64; 3] {
        Self::drag_force(relative_vel, reference_area, cd, fluid.density)
    }
}
/// Wave loading via the Morison equation.
pub struct WaveForce;
impl WaveForce {
    /// Morison force: F = 0.5·ρ·Cd·A·(u_rel|u_rel) + ρ·V·Cm·a.
    #[allow(clippy::too_many_arguments)]
    pub fn morison_force(
        diameter: f64,
        length_val: f64,
        fluid: &FluidProperties,
        wave_vel: [f64; 3],
        wave_acc: [f64; 3],
        body_vel: [f64; 3],
        cd: f64,
        cm: f64,
    ) -> [f64; 3] {
        let cross_section = 0.25 * PI * diameter * diameter;
        let volume = cross_section * length_val;
        let u_rel = sub(wave_vel, body_vel);
        let speed = length(u_rel);
        let drag = if speed < 1e-300 {
            [0.0; 3]
        } else {
            scale(u_rel, 0.5 * fluid.density * cd * cross_section * speed)
        };
        let inertia = scale(wave_acc, fluid.density * volume * cm);
        add(drag, inertia)
    }
    /// Linear wave theory particle velocity (surface approximation).
    pub fn wave_particle_velocity(
        amplitude: f64,
        omega: f64,
        k: f64,
        x: f64,
        _z: f64,
        t: f64,
    ) -> [f64; 3] {
        let phase = k * x - omega * t;
        [
            amplitude * omega * phase.cos(),
            0.0,
            -amplitude * omega * phase.sin(),
        ]
    }
    /// Linear wave theory particle acceleration (surface approximation).
    ///
    /// a_x = A*ω²*sin(kx - ωt), a_z = -A*ω²*cos(kx - ωt)
    pub fn wave_particle_acceleration(
        amplitude: f64,
        omega: f64,
        k: f64,
        x: f64,
        _z: f64,
        t: f64,
    ) -> [f64; 3] {
        let phase = k * x - omega * t;
        [
            amplitude * omega * omega * phase.sin(),
            0.0,
            -amplitude * omega * omega * phase.cos(),
        ]
    }
    /// Wave particle velocity with depth decay (Airy wave theory).
    ///
    /// u = A*ω * cosh(k*(z+d)) / sinh(k*d) * cos(kx - ωt)
    /// w = A*ω * sinh(k*(z+d)) / sinh(k*d) * sin(kx - ωt)  (note: sign from convention)
    #[allow(clippy::too_many_arguments)]
    pub fn wave_velocity_with_depth(
        amplitude: f64,
        omega: f64,
        k: f64,
        x: f64,
        z: f64,
        t: f64,
        water_depth: f64,
    ) -> [f64; 3] {
        let phase = k * x - omega * t;
        let kd = k * water_depth;
        let sinh_kd = kd.sinh();
        if sinh_kd.abs() < 1e-15 {
            return [0.0; 3];
        }
        let kzd = k * (z + water_depth);
        let u = amplitude * omega * kzd.cosh() / sinh_kd * phase.cos();
        let w = amplitude * omega * kzd.sinh() / sinh_kd * (-phase.sin());
        [u, 0.0, w]
    }
    /// Deep-water dispersion relation: ω² = g * k.
    pub fn deep_water_omega(k: f64, gravity: f64) -> f64 {
        (gravity * k).sqrt()
    }
    /// Wave length from period: L = g*T²/(2π) (deep water).
    pub fn deep_water_wavelength(period: f64, gravity: f64) -> f64 {
        gravity * period * period / (2.0 * PI)
    }
    /// Compute the maximum wave force on a vertical cylinder using
    /// Morison over one wave cycle.
    ///
    /// Returns (max_drag_force, max_inertia_force).
    pub fn morison_max_components(
        diameter: f64,
        length_val: f64,
        fluid_density: f64,
        amplitude: f64,
        omega: f64,
        cd: f64,
        cm: f64,
    ) -> (f64, f64) {
        let area = 0.25 * PI * diameter * diameter;
        let volume = area * length_val;
        let u_max = amplitude * omega;
        let a_max = amplitude * omega * omega;
        let f_drag_max = 0.5 * fluid_density * cd * area * u_max * u_max;
        let f_inertia_max = fluid_density * volume * cm * a_max;
        (f_drag_max, f_inertia_max)
    }
}
/// Froude number utilities and flow regime classification.
#[allow(dead_code)]
pub struct FroudeNumber;
impl FroudeNumber {
    /// Froude number: Fr = U / sqrt(g * L).
    ///
    /// * `u` – characteristic velocity (m/s)
    /// * `l` – characteristic length (m)
    /// * `g` – gravitational acceleration magnitude (m/s²)
    pub fn froude(u: f64, l: f64, g: f64) -> f64 {
        if l <= 0.0 || g <= 0.0 {
            return 0.0;
        }
        u / (g * l).sqrt()
    }
    /// Classify the flow regime from the Froude number.
    ///
    /// * Fr < 0.95 → SubCritical
    /// * 0.95 ≤ Fr ≤ 1.05 → Critical
    /// * Fr > 1.05 → SuperCritical
    pub fn classify(fr: f64) -> FlowRegime {
        if fr > 1.05 {
            FlowRegime::SuperCritical
        } else if fr >= 0.95 {
            FlowRegime::Critical
        } else {
            FlowRegime::SubCritical
        }
    }
    /// Ship Froude number (uses waterline length as the length scale).
    pub fn ship_froude(speed_knots: f64, waterline_length_m: f64) -> f64 {
        let u = speed_knots * 0.5144;
        Self::froude(u, waterline_length_m, 9.81)
    }
    /// Wave celerity (phase speed) at the Froude critical point: c = sqrt(g * h).
    pub fn critical_celerity(depth: f64, gravity: f64) -> f64 {
        (gravity * depth).sqrt()
    }
    /// Hydraulic jump indicator: if sub-critical → super-critical transition exists.
    pub fn has_hydraulic_jump(fr_upstream: f64, fr_downstream: f64) -> bool {
        fr_upstream < 1.0 && fr_downstream > 1.0 || fr_upstream > 1.0 && fr_downstream < 1.0
    }
    /// Conjugate depth ratio for a hydraulic jump (Bélanger equation).
    ///
    /// y2/y1 = 0.5 * (sqrt(1 + 8*Fr1²) - 1)
    pub fn conjugate_depth_ratio(fr1: f64) -> f64 {
        0.5 * ((1.0 + 8.0 * fr1 * fr1).sqrt() - 1.0)
    }
}
/// Tools to compute buoyancy forces from submerged volume fractions.
///
/// Supports sphere, cylinder, and rectangular box geometries.
#[allow(dead_code)]
pub struct SubmergedVolumeFraction;
impl SubmergedVolumeFraction {
    /// Submerged volume fraction for a sphere whose centre is at height `z_c`
    /// above the free surface (negative = below surface).
    ///
    /// Returns 0.0 (fully above) .. 1.0 (fully submerged).
    pub fn sphere(radius: f64, z_centre: f64) -> f64 {
        let h = (radius - z_centre).clamp(0.0, 2.0 * radius);
        let total_vol = (4.0 / 3.0) * PI * radius.powi(3);
        if total_vol <= 0.0 {
            return 0.0;
        }
        let cap_vol = PI * h * h * (3.0 * radius - h) / 3.0;
        (cap_vol / total_vol).clamp(0.0, 1.0)
    }
    /// Submerged volume fraction for an upright rectangular box
    /// (width × depth × height = `w × d × h`), whose bottom face is at `z_bottom`.
    ///
    /// Free surface at z = 0.
    pub fn box_fraction(box_height: f64, z_bottom: f64) -> f64 {
        let z_top = z_bottom + box_height;
        if z_bottom >= 0.0 {
            return 0.0;
        }
        if z_top <= 0.0 {
            return 1.0;
        }
        let h_sub = (-z_bottom).clamp(0.0, box_height);
        h_sub / box_height
    }
    /// Buoyancy force from a submerged volume fraction.
    ///
    /// F_b = ρ_f * g * fraction * total_volume  (upward).
    pub fn buoyancy_force(
        fluid_density: f64,
        gravity_mag: f64,
        fraction: f64,
        total_volume: f64,
    ) -> f64 {
        fluid_density * gravity_mag * fraction.clamp(0.0, 1.0) * total_volume
    }
}
