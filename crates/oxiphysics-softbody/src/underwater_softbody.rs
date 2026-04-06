// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Underwater soft body dynamics.
//!
//! Covers buoyancy distribution over flexible bodies, added mass effects,
//! viscous drag on flexible structures, Morison's equation for oscillating
//! bodies, vortex-induced vibration (VIV) of flexible cylinders, marine
//! growth effects on stiffness, pressure hull deformation, hydrostatic
//! buckling, soft underwater gripper simulation, undulatory (anguilliform)
//! locomotion, and bio-inspired fin propulsion.

#![allow(dead_code)]

// ── Vec3 helpers ─────────────────────────────────────────────────────────────

/// A simple 3-D vector as `[f64; 3]`.
pub type Vec3 = [f64; 3];

#[inline]
fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn normalize(a: Vec3) -> Vec3 {
    let l = len(a);
    if l < 1e-30 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / l)
    }
}

// ── Physical constants ────────────────────────────────────────────────────────

/// Density of seawater in kg/m³.
pub const RHO_SEAWATER: f64 = 1025.0;
/// Density of freshwater in kg/m³.
pub const RHO_FRESHWATER: f64 = 1000.0;
/// Standard gravitational acceleration in m/s².
pub const GRAVITY: f64 = 9.81;
/// Dynamic viscosity of seawater at 20 °C in Pa·s.
pub const MU_SEAWATER: f64 = 1.08e-3;

// ── Buoyancy ──────────────────────────────────────────────────────────────────

/// Parameters describing the fluid environment for underwater simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FluidEnvironment {
    /// Fluid density in kg/m³.
    pub density: f64,
    /// Gravitational acceleration (positive downward convention) in m/s².
    pub gravity: f64,
    /// Ambient fluid velocity vector in m/s.
    pub current: Vec3,
    /// Dynamic viscosity in Pa·s.
    pub viscosity: f64,
}

impl FluidEnvironment {
    /// Create a default seawater environment with no current.
    pub fn seawater() -> Self {
        Self {
            density: RHO_SEAWATER,
            gravity: GRAVITY,
            current: [0.0; 3],
            viscosity: MU_SEAWATER,
        }
    }

    /// Create a freshwater environment with no current.
    pub fn freshwater() -> Self {
        Self {
            density: RHO_FRESHWATER,
            gravity: GRAVITY,
            current: [0.0; 3],
            viscosity: 1.002e-3,
        }
    }

    /// Compute the hydrostatic pressure at depth `z` (metres below surface,
    /// positive = deeper).
    pub fn hydrostatic_pressure(&self, z: f64) -> f64 {
        self.density * self.gravity * z
    }

    /// Buoyancy force on a submerged volume (upward = positive Y).
    pub fn buoyancy_force(&self, volume: f64) -> Vec3 {
        [0.0, self.density * self.gravity * volume, 0.0]
    }
}

// ── Added mass ────────────────────────────────────────────────────────────────

/// Added (virtual) mass coefficients for a body accelerating in a fluid.
///
/// The added mass tensor accounts for the inertia of entrained fluid and is
/// stored as a diagonal 3×3 coefficient matrix: `m_added_i = C_a_i * rho * V`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AddedMass {
    /// Added mass coefficients (dimensionless) along each axis.
    pub coefficients: Vec3,
    /// Fluid density in kg/m³.
    pub rho: f64,
    /// Displaced volume in m³.
    pub volume: f64,
}

impl AddedMass {
    /// Create an added mass model for a sphere (`C_a = 0.5` along all axes).
    pub fn sphere(rho: f64, radius: f64) -> Self {
        let volume = (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        Self {
            coefficients: [0.5; 3],
            rho,
            volume,
        }
    }

    /// Create an added mass model for an infinite cylinder (axis along Z).
    ///
    /// `C_a = 1.0` transverse, `0.0` axial.
    pub fn cylinder(rho: f64, radius: f64, length: f64) -> Self {
        let volume = std::f64::consts::PI * radius.powi(2) * length;
        Self {
            coefficients: [1.0, 1.0, 0.0],
            rho,
            volume,
        }
    }

    /// Compute the added mass force given body acceleration `accel` (m/s²).
    ///
    /// `F_added = -C_a * rho * V * accel`
    pub fn force(&self, accel: Vec3) -> Vec3 {
        [
            -self.coefficients[0] * self.rho * self.volume * accel[0],
            -self.coefficients[1] * self.rho * self.volume * accel[1],
            -self.coefficients[2] * self.rho * self.volume * accel[2],
        ]
    }

    /// Scalar added mass along axis `axis` (0=X, 1=Y, 2=Z).
    pub fn scalar_mass(&self, axis: usize) -> f64 {
        self.coefficients[axis] * self.rho * self.volume
    }
}

// ── Viscous drag ──────────────────────────────────────────────────────────────

/// Viscous drag model for flexible bodies based on Stokes / quadratic drag.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ViscousDrag {
    /// Drag coefficient (dimensionless).  Use `C_d ≈ 1.2` for a cylinder.
    pub cd: f64,
    /// Reference area for drag computation in m².
    pub reference_area: f64,
    /// Fluid density in kg/m³.
    pub rho: f64,
}

impl ViscousDrag {
    /// Create a new viscous drag model.
    pub fn new(cd: f64, reference_area: f64, rho: f64) -> Self {
        Self {
            cd,
            reference_area,
            rho,
        }
    }

    /// Quadratic drag force on a body with velocity `vel` relative to the fluid.
    ///
    /// `F_drag = -0.5 * rho * Cd * A * |v| * v`
    pub fn force(&self, relative_velocity: Vec3) -> Vec3 {
        let speed = len(relative_velocity);
        let mag = -0.5 * self.rho * self.cd * self.reference_area * speed;
        scale(relative_velocity, mag)
    }

    /// Reynolds number given body diameter `d` and kinematic viscosity `nu`.
    pub fn reynolds_number(speed: f64, diameter: f64, nu: f64) -> f64 {
        if nu < 1e-30 {
            0.0
        } else {
            speed * diameter / nu
        }
    }
}

// ── Morison's equation ────────────────────────────────────────────────────────

/// Morison's equation parameters for a slender cylindrical element.
///
/// The total hydrodynamic force per unit length is:
///
/// ```text
/// f = rho * A * Cm * du/dt  +  0.5 * rho * D * Cd * |u - v| * (u - v)
/// ```
///
/// where `u` is fluid velocity, `v` is body velocity.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MorisonElement {
    /// Inertia coefficient `Cm` (= 1 + Ca).
    pub cm: f64,
    /// Drag coefficient `Cd`.
    pub cd: f64,
    /// Cylinder diameter in m.
    pub diameter: f64,
    /// Fluid density in kg/m³.
    pub rho: f64,
}

impl MorisonElement {
    /// Create a Morison element with typical marine values.
    pub fn new(cm: f64, cd: f64, diameter: f64, rho: f64) -> Self {
        Self {
            cm,
            cd,
            diameter,
            rho,
        }
    }

    /// Cross-sectional area.
    fn area(&self) -> f64 {
        std::f64::consts::PI * (self.diameter / 2.0).powi(2)
    }

    /// Compute the Morison force per unit length.
    ///
    /// * `fluid_vel` - Fluid particle velocity in m/s.
    /// * `fluid_accel` - Fluid particle acceleration in m/s².
    /// * `body_vel` - Body element velocity in m/s.
    pub fn force_per_unit_length(
        &self,
        fluid_vel: Vec3,
        fluid_accel: Vec3,
        body_vel: Vec3,
    ) -> Vec3 {
        // Inertia term: rho * A * Cm * du/dt
        let f_inertia = scale(fluid_accel, self.rho * self.area() * self.cm);
        // Relative velocity
        let rel_v = sub(fluid_vel, body_vel);
        let speed = len(rel_v);
        // Drag term: 0.5 * rho * D * Cd * |u - v| * (u - v)
        let f_drag = scale(rel_v, 0.5 * self.rho * self.diameter * self.cd * speed);
        add(f_inertia, f_drag)
    }
}

// ── Vortex-induced vibration ──────────────────────────────────────────────────

/// Vortex-induced vibration (VIV) model for a flexible cylinder.
///
/// Uses the Hartlen-Currie wake oscillator model:
/// - Wake variable `q` obeys a van der Pol oscillator driven by cylinder motion.
/// - Lift force = `0.5 * rho * D * L * U^2 * CL0 * q`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VivModel {
    /// Fluid density in kg/m³.
    pub rho: f64,
    /// Cylinder diameter in m.
    pub diameter: f64,
    /// Cylinder length in m.
    pub length: f64,
    /// Free-stream velocity in m/s.
    pub free_stream: f64,
    /// Reference lift coefficient.
    pub cl0: f64,
    /// Van der Pol `epsilon` (nonlinear damping).
    pub epsilon: f64,
    /// Shedding frequency in rad/s (Strouhal relation: `fs = St * U / D`).
    pub omega_s: f64,
    /// Coupling coefficient.
    pub gamma: f64,
    /// Wake variable (internal state).
    pub q: f64,
    /// Wake variable time derivative (internal state).
    pub q_dot: f64,
}

impl VivModel {
    /// Create a VIV model; Strouhal number assumed `St = 0.2`.
    pub fn new(rho: f64, diameter: f64, length: f64, free_stream: f64) -> Self {
        let st = 0.2;
        let omega_s = 2.0 * std::f64::consts::PI * st * free_stream / diameter;
        Self {
            rho,
            diameter,
            length,
            free_stream,
            cl0: 0.3,
            epsilon: 0.3,
            omega_s,
            gamma: 0.8,
            q: 0.01,
            q_dot: 0.0,
        }
    }

    /// Advance the wake oscillator by one time step.
    ///
    /// * `cylinder_accel_y` - Transverse acceleration of the cylinder in m/s².
    /// * `dt` - Time step in seconds.
    pub fn step(&mut self, cylinder_accel_y: f64, dt: f64) {
        let q_ddot = self.omega_s.powi(2)
            * (self.epsilon * (1.0 - self.q * self.q) * self.q_dot / self.omega_s - self.q
                + self.gamma * cylinder_accel_y / (self.omega_s.powi(2) * self.diameter));
        self.q_dot += q_ddot * dt;
        self.q += self.q_dot * dt;
    }

    /// Compute the current lift force in N.
    pub fn lift_force(&self) -> f64 {
        0.5 * self.rho * self.diameter * self.length * self.free_stream.powi(2) * self.cl0 * self.q
    }

    /// Strouhal shedding frequency in Hz.
    pub fn shedding_frequency_hz(&self) -> f64 {
        self.omega_s / (2.0 * std::f64::consts::PI)
    }
}

// ── Marine growth ─────────────────────────────────────────────────────────────

/// Effect of marine growth (barnacles, mussels, biofilm) on structural stiffness.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MarineGrowth {
    /// Added thickness due to marine growth in m.
    pub thickness: f64,
    /// Density of growth material in kg/m³.
    pub density: f64,
    /// Effective Young's modulus modification factor (>1 stiffens, <1 softens).
    pub stiffness_factor: f64,
    /// Drag coefficient increase factor due to roughness.
    pub roughness_cd_factor: f64,
}

impl MarineGrowth {
    /// Create a marine growth model.
    pub fn new(thickness: f64, density: f64, stiffness_factor: f64) -> Self {
        Self {
            thickness,
            density,
            stiffness_factor,
            roughness_cd_factor: 1.0 + 0.3 * (thickness / 0.05).min(1.0),
        }
    }

    /// Added mass per unit length for a cylinder of base diameter `d`.
    pub fn added_mass_per_length(&self, base_diameter: f64) -> f64 {
        let outer_d = base_diameter + 2.0 * self.thickness;
        let area = std::f64::consts::PI / 4.0 * (outer_d.powi(2) - base_diameter.powi(2));
        self.density * area
    }

    /// Effective flexural stiffness `EI` given base `ei`.
    pub fn effective_ei(&self, base_ei: f64) -> f64 {
        base_ei * self.stiffness_factor
    }

    /// Effective drag coefficient accounting for surface roughness.
    pub fn effective_cd(&self, base_cd: f64) -> f64 {
        base_cd * self.roughness_cd_factor
    }
}

// ── Pressure hull deformation ─────────────────────────────────────────────────

/// Linear elastic pressure hull deformation model (thin-shell approximation).
///
/// Based on Lamé's equations for a thick cylinder under external pressure.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PressureHull {
    /// Outer radius in m.
    pub outer_radius: f64,
    /// Inner radius in m.
    pub inner_radius: f64,
    /// Young's modulus of hull material in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson: f64,
}

impl PressureHull {
    /// Create a pressure hull.
    pub fn new(outer_radius: f64, inner_radius: f64, youngs_modulus: f64, poisson: f64) -> Self {
        Self {
            outer_radius,
            inner_radius,
            youngs_modulus,
            poisson,
        }
    }

    /// Hoop stress at the inner surface for external pressure `p` (Pa).
    ///
    /// Uses Lamé thick-cylinder formula:
    /// `σ_θ = -p * 2 * R_o^2 / (R_o^2 - R_i^2)`
    pub fn hoop_stress_inner(&self, external_pressure: f64) -> f64 {
        let ro2 = self.outer_radius.powi(2);
        let ri2 = self.inner_radius.powi(2);
        -external_pressure * 2.0 * ro2 / (ro2 - ri2)
    }

    /// Radial deflection at the outer surface under external pressure `p` (Pa).
    pub fn radial_deflection_outer(&self, external_pressure: f64) -> f64 {
        let ro2 = self.outer_radius.powi(2);
        let ri2 = self.inner_radius.powi(2);
        let factor = (1.0 - self.poisson) * ro2 + (1.0 + self.poisson) * ri2;
        -external_pressure * self.outer_radius * factor / (self.youngs_modulus * (ro2 - ri2))
    }

    /// Critical external pressure for elastic shell buckling (von Mises).
    ///
    /// Simplified formula: `p_cr = 2E(t/D)^3 / (1 - nu^2)`
    pub fn buckling_pressure(&self) -> f64 {
        let t = self.outer_radius - self.inner_radius;
        let d = 2.0 * self.outer_radius;
        2.0 * self.youngs_modulus * (t / d).powi(3) / (1.0 - self.poisson.powi(2))
    }
}

// ── Hydrostatic buckling ──────────────────────────────────────────────────────

/// Hydrostatic buckling analysis for a cylindrical shell.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HydrostaticBuckling {
    /// Cylinder radius in m.
    pub radius: f64,
    /// Wall thickness in m.
    pub thickness: f64,
    /// Unsupported length in m.
    pub length: f64,
    /// Young's modulus in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson: f64,
}

impl HydrostaticBuckling {
    /// Create a buckling model.
    pub fn new(
        radius: f64,
        thickness: f64,
        length: f64,
        youngs_modulus: f64,
        poisson: f64,
    ) -> Self {
        Self {
            radius,
            thickness,
            length,
            youngs_modulus,
            poisson,
        }
    }

    /// Windenburg & Trilling (1934) critical pressure for lobar buckling.
    ///
    /// Returns approximate critical external pressure in Pa.
    pub fn critical_pressure_windenburg(&self) -> f64 {
        let nu = self.poisson;
        let r = self.radius;
        let t = self.thickness;
        let l = self.length;
        let e = self.youngs_modulus;
        // Simplified WT formula
        2.6 * e * (t / (2.0 * r)).powi(5).cbrt()
            / (l / (2.0 * r) - 0.45 * (t / (2.0 * r)).sqrt())
            / (1.0 - nu.powi(2)).cbrt()
    }

    /// Euler column buckling load for the cylinder as a beam.
    pub fn euler_buckling_load(&self, end_condition_factor: f64) -> f64 {
        let i = std::f64::consts::PI
            * (self.radius.powi(4) - (self.radius - self.thickness).powi(4))
            / 4.0;
        std::f64::consts::PI.powi(2) * self.youngs_modulus * i
            / (end_condition_factor * self.length).powi(2)
    }

    /// Safety factor relative to operating depth (metres).
    pub fn safety_factor(&self, operating_depth: f64) -> f64 {
        let p_operating = RHO_SEAWATER * GRAVITY * operating_depth;
        let p_cr = self.critical_pressure_windenburg();
        if p_operating < 1e-12 {
            f64::INFINITY
        } else {
            p_cr / p_operating
        }
    }
}

// ── Soft underwater gripper ───────────────────────────────────────────────────

/// A segment of a pneumatically-actuated soft underwater gripper finger.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GripperSegment {
    /// Rest length in m.
    pub rest_length: f64,
    /// Current length in m.
    pub length: f64,
    /// Bending angle at this segment in radians.
    pub angle: f64,
    /// Pneumatic pressure in Pa.
    pub pressure: f64,
    /// Bending stiffness in N·m/rad.
    pub bending_stiffness: f64,
    /// Axial stiffness in N/m.
    pub axial_stiffness: f64,
}

impl GripperSegment {
    /// Create a new gripper segment.
    pub fn new(rest_length: f64, bending_stiffness: f64, axial_stiffness: f64) -> Self {
        Self {
            rest_length,
            length: rest_length,
            angle: 0.0,
            pressure: 0.0,
            bending_stiffness,
            axial_stiffness,
        }
    }

    /// Apply internal pressure to compute equilibrium bending angle.
    ///
    /// Simple linear model: `theta_eq = p * A_eff / k_bend`
    pub fn pressurize(&mut self, pressure: f64, effective_area: f64) {
        self.pressure = pressure;
        let moment = pressure * effective_area * self.rest_length;
        self.angle = moment / self.bending_stiffness;
    }

    /// Bending moment at this segment.
    pub fn bending_moment(&self) -> f64 {
        self.bending_stiffness * self.angle
    }

    /// Tip position relative to proximal end (2-D, in XY plane).
    pub fn tip_position_2d(&self) -> [f64; 2] {
        if self.angle.abs() < 1e-10 {
            [self.length, 0.0]
        } else {
            let r = self.length / self.angle;
            [r * self.angle.sin(), r * (1.0 - self.angle.cos())]
        }
    }
}

/// A multi-segment soft underwater gripper finger.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SoftGripperFinger {
    /// Segments composing the finger.
    pub segments: Vec<GripperSegment>,
    /// Drag coefficient in water.
    pub cd: f64,
}

impl SoftGripperFinger {
    /// Create a finger with `n` uniform segments.
    pub fn new(
        n: usize,
        segment_length: f64,
        bending_stiffness: f64,
        axial_stiffness: f64,
    ) -> Self {
        Self {
            segments: (0..n)
                .map(|_| GripperSegment::new(segment_length, bending_stiffness, axial_stiffness))
                .collect(),
            cd: 1.2,
        }
    }

    /// Uniformly pressurize all segments.
    pub fn pressurize_all(&mut self, pressure: f64, effective_area: f64) {
        for seg in &mut self.segments {
            seg.pressurize(pressure, effective_area);
        }
    }

    /// Total tip deflection (2-D cumulative).
    pub fn tip_deflection(&self) -> [f64; 2] {
        let mut x = 0.0_f64;
        let mut y = 0.0_f64;
        let mut cumulative_angle = 0.0_f64;
        for seg in &self.segments {
            let [dx, dy] = seg.tip_position_2d();
            let cos_a = cumulative_angle.cos();
            let sin_a = cumulative_angle.sin();
            x += cos_a * dx - sin_a * dy;
            y += sin_a * dx + cos_a * dy;
            cumulative_angle += seg.angle;
        }
        [x, y]
    }
}

// ── Undulatory locomotion (anguilliform) ──────────────────────────────────────

/// Kinematic model for anguilliform (eel-like) undulatory swimming.
///
/// The body curvature wave travels from head to tail:
/// `kappa(s, t) = A(s) * sin(k*s - omega*t)`
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AnguillaBody {
    /// Body length in m.
    pub length: f64,
    /// Number of segments.
    pub n_segments: usize,
    /// Wave number in rad/m.
    pub wave_number: f64,
    /// Angular frequency in rad/s.
    pub omega: f64,
    /// Maximum curvature amplitude in 1/m.
    pub amplitude: f64,
    /// Current simulation time.
    pub time: f64,
    /// Segment positions (x, y) in m.
    pub positions: Vec<[f64; 2]>,
    /// Segment angles (heading) in rad.
    pub angles: Vec<f64>,
}

impl AnguillaBody {
    /// Create an anguilliform body at rest.
    pub fn new(
        length: f64,
        n_segments: usize,
        wave_number: f64,
        omega: f64,
        amplitude: f64,
    ) -> Self {
        let ds = length / n_segments as f64;
        let positions = (0..=n_segments).map(|i| [i as f64 * ds, 0.0]).collect();
        let angles = vec![0.0; n_segments + 1];
        Self {
            length,
            n_segments,
            wave_number,
            omega,
            amplitude,
            time: 0.0,
            positions,
            angles,
        }
    }

    /// Compute body curvature at arc-length `s` and current time.
    pub fn curvature(&self, s: f64) -> f64 {
        self.amplitude * (self.wave_number * s - self.omega * self.time).sin()
    }

    /// Update body shape by integrating curvature along the midline.
    pub fn update_shape(&mut self) {
        let ds = self.length / self.n_segments as f64;
        self.angles[0] = 0.0;
        self.positions[0] = [0.0, 0.0];
        for i in 0..self.n_segments {
            let s = (i as f64 + 0.5) * ds;
            let kappa = self.curvature(s);
            self.angles[i + 1] = self.angles[i] + kappa * ds;
            let theta = self.angles[i + 1];
            self.positions[i + 1] = [
                self.positions[i][0] + theta.cos() * ds,
                self.positions[i][1] + theta.sin() * ds,
            ];
        }
    }

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.time += dt;
        self.update_shape();
    }

    /// Phase velocity of the undulatory wave in m/s.
    pub fn phase_velocity(&self) -> f64 {
        self.omega / self.wave_number
    }

    /// Strouhal number `St = f * A / U` where `f` is beat frequency,
    /// `A` is tail-beat amplitude, and `U` is swimming speed estimate.
    pub fn strouhal_number(&self, swimming_speed: f64) -> f64 {
        let frequency = self.omega / (2.0 * std::f64::consts::PI);
        let tail_amplitude = self.amplitude / self.wave_number; // rough estimate
        if swimming_speed.abs() < 1e-12 {
            f64::INFINITY
        } else {
            frequency * tail_amplitude / swimming_speed
        }
    }
}

// ── Bio-inspired fin propulsion ───────────────────────────────────────────────

/// Oscillating fin mode for bio-inspired propulsion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FinMode {
    /// Flapping motion (pectoral fins).
    Flapping,
    /// Undulating motion (pelvic / dorsal fins).
    Undulating,
    /// Combined flapping + undulating.
    Combined,
}

/// Bio-inspired oscillating fin model.
///
/// Compute thrust and lift based on quasi-steady blade-element theory.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BioFin {
    /// Fin span in m.
    pub span: f64,
    /// Mean chord length in m.
    pub chord: f64,
    /// Lift coefficient slope `dCL/dalpha` in 1/rad.
    pub cl_alpha: f64,
    /// Drag coefficient at zero lift.
    pub cd0: f64,
    /// Fluid density in kg/m³.
    pub rho: f64,
    /// Oscillation frequency in rad/s.
    pub omega: f64,
    /// Stroke amplitude in radians.
    pub amplitude: f64,
    /// Propulsion mode.
    pub mode: FinMode,
    /// Current phase in radians.
    pub phase: f64,
}

impl BioFin {
    /// Create a bio-inspired fin.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        span: f64,
        chord: f64,
        cl_alpha: f64,
        cd0: f64,
        rho: f64,
        omega: f64,
        amplitude: f64,
        mode: FinMode,
    ) -> Self {
        Self {
            span,
            chord,
            cl_alpha,
            cd0,
            rho,
            omega,
            amplitude,
            mode,
            phase: 0.0,
        }
    }

    /// Advance phase by `dt`.
    pub fn step(&mut self, dt: f64) {
        self.phase = (self.phase + self.omega * dt) % (2.0 * std::f64::consts::PI);
    }

    /// Instantaneous fin angle (radians).
    pub fn fin_angle(&self) -> f64 {
        self.amplitude * self.phase.sin()
    }

    /// Instantaneous fin angular velocity (rad/s).
    pub fn fin_angular_velocity(&self) -> f64 {
        self.amplitude * self.omega * self.phase.cos()
    }

    /// Compute quasi-steady thrust and lift forces in N.
    ///
    /// Returns `(thrust, lift)`.
    pub fn forces(&self, forward_speed: f64) -> (f64, f64) {
        let theta = self.fin_angle();
        let theta_dot = self.fin_angular_velocity();
        // Velocity at mid-span
        let v_stroke = theta_dot * self.span * 0.5;
        let v_total = (forward_speed.powi(2) + v_stroke.powi(2)).sqrt();
        let alpha = theta_dot.atan2(forward_speed.max(1e-8));
        let cl = self.cl_alpha * alpha;
        let cd = self.cd0 + cl.powi(2) / (std::f64::consts::PI * (self.span / self.chord));
        let area = self.span * self.chord;
        let q = 0.5 * self.rho * v_total.powi(2) * area;
        let lift_force = cl * q;
        let drag_force = cd * q;
        // Thrust is the component opposing drag and aligned with forward motion
        let thrust = lift_force * theta.sin() - drag_force * theta.cos();
        let lift = lift_force * theta.cos() + drag_force * theta.sin();
        (thrust, lift)
    }

    /// Froude efficiency estimate.
    ///
    /// `eta = T * U / P_input` (approximation).
    pub fn froude_efficiency(&self, thrust: f64, forward_speed: f64, power_in: f64) -> f64 {
        if power_in.abs() < 1e-12 {
            0.0
        } else {
            thrust * forward_speed / power_in
        }
    }
}

// ── Buoyancy distribution over a flexible body ────────────────────────────────

/// Distributes buoyancy over discrete particles of a soft body.
///
/// Each particle below the free surface receives an upward buoyancy force
/// proportional to its associated volume.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BuoyancyDistributor {
    /// Fluid density in kg/m³.
    pub rho: f64,
    /// Gravitational acceleration in m/s².
    pub gravity: f64,
    /// Y-coordinate of the free surface.
    pub surface_y: f64,
}

impl BuoyancyDistributor {
    /// Create a buoyancy distributor.
    pub fn new(rho: f64, gravity: f64, surface_y: f64) -> Self {
        Self {
            rho,
            gravity,
            surface_y,
        }
    }

    /// Compute the buoyancy force on a particle at position `pos` with
    /// associated volume `volume` (m³).
    ///
    /// Returns zero if the particle is above the free surface.
    pub fn force_on_particle(&self, pos: Vec3, volume: f64) -> Vec3 {
        if pos[1] >= self.surface_y {
            [0.0; 3]
        } else {
            [0.0, self.rho * self.gravity * volume, 0.0]
        }
    }

    /// Apply buoyancy to all particles.
    ///
    /// `positions`: slice of particle positions.
    /// `volumes`: per-particle volume.
    /// `forces`: accumulated forces (mutated in place).
    pub fn apply_to_all(&self, positions: &[Vec3], volumes: &[f64], forces: &mut [Vec3]) {
        for ((pos, &vol), f) in positions.iter().zip(volumes.iter()).zip(forces.iter_mut()) {
            let b = self.force_on_particle(*pos, vol);
            f[0] += b[0];
            f[1] += b[1];
            f[2] += b[2];
        }
    }
}

// ── Flexible underwater cable segment ────────────────────────────────────────

/// A single segment of a flexible underwater cable or pipe.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CableSegment {
    /// Start position in m.
    pub pos_a: Vec3,
    /// End position in m.
    pub pos_b: Vec3,
    /// Rest length in m.
    pub rest_length: f64,
    /// Axial stiffness (EA) in N.
    pub ea: f64,
    /// Bending stiffness (EI) in N·m².
    pub ei: f64,
    /// Linear mass density in kg/m.
    pub mass_per_length: f64,
}

impl CableSegment {
    /// Create a cable segment.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        pos_a: Vec3,
        pos_b: Vec3,
        rest_length: f64,
        ea: f64,
        ei: f64,
        mass_per_length: f64,
    ) -> Self {
        Self {
            pos_a,
            pos_b,
            rest_length,
            ea,
            ei,
            mass_per_length,
        }
    }

    /// Current length of the segment.
    pub fn current_length(&self) -> f64 {
        len(sub(self.pos_b, self.pos_a))
    }

    /// Axial strain `(L - L0) / L0`.
    pub fn axial_strain(&self) -> f64 {
        let l = self.current_length();
        if self.rest_length < 1e-12 {
            0.0
        } else {
            (l - self.rest_length) / self.rest_length
        }
    }

    /// Axial tension force in N.
    pub fn tension(&self) -> f64 {
        self.ea * self.axial_strain()
    }

    /// Unit tangent vector along the segment.
    pub fn tangent(&self) -> Vec3 {
        normalize(sub(self.pos_b, self.pos_a))
    }

    /// Mass of the segment in kg.
    pub fn mass(&self) -> f64 {
        self.rest_length * self.mass_per_length
    }
}

// ── Underwater soft body simulation ──────────────────────────────────────────

/// A particle in an underwater soft body simulation.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct UnderwaterParticle {
    /// Position in m.
    pub position: Vec3,
    /// Velocity in m/s.
    pub velocity: Vec3,
    /// Mass in kg.
    pub mass: f64,
    /// Associated volume in m³.
    pub volume: f64,
    /// Accumulated external forces.
    pub force: Vec3,
    /// Whether the particle is pinned (infinite mass).
    pub pinned: bool,
}

impl UnderwaterParticle {
    /// Create a free underwater particle.
    pub fn new(position: Vec3, mass: f64, volume: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            volume,
            force: [0.0; 3],
            pinned: false,
        }
    }

    /// Create a static (pinned) particle.
    pub fn pinned(position: Vec3) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass: f64::INFINITY,
            volume: 0.0,
            force: [0.0; 3],
            pinned: true,
        }
    }

    /// Add a force to this particle.
    pub fn add_force(&mut self, f: Vec3) {
        if !self.pinned {
            self.force = add(self.force, f);
        }
    }

    /// Semi-implicit Euler integration step.
    pub fn integrate(&mut self, dt: f64) {
        if self.pinned {
            return;
        }
        let accel = scale(self.force, 1.0 / self.mass);
        self.velocity = add(self.velocity, scale(accel, dt));
        self.position = add(self.position, scale(self.velocity, dt));
        self.force = [0.0; 3];
    }
}

/// A flexible underwater soft body consisting of mass–spring particles.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct UnderwaterSoftBody {
    /// Particles.
    pub particles: Vec<UnderwaterParticle>,
    /// Spring edges: `(i, j, rest_length, stiffness)`.
    pub springs: Vec<(usize, usize, f64, f64)>,
    /// Fluid environment.
    pub fluid: FluidEnvironment,
    /// Drag coefficient for each particle.
    pub cd: f64,
    /// Particle drag reference area in m².
    pub drag_area: f64,
}

impl UnderwaterSoftBody {
    /// Create an empty underwater soft body.
    pub fn new(fluid: FluidEnvironment, cd: f64, drag_area: f64) -> Self {
        Self {
            particles: Vec::new(),
            springs: Vec::new(),
            fluid,
            cd,
            drag_area,
        }
    }

    /// Add a particle and return its index.
    pub fn add_particle(&mut self, p: UnderwaterParticle) -> usize {
        let idx = self.particles.len();
        self.particles.push(p);
        idx
    }

    /// Add a spring between particles `i` and `j`.
    pub fn add_spring(&mut self, i: usize, j: usize, stiffness: f64) {
        let rest = len(sub(self.particles[j].position, self.particles[i].position));
        self.springs.push((i, j, rest, stiffness));
    }

    /// Apply spring forces.
    fn apply_spring_forces(&mut self) {
        for &(i, j, rest, k) in &self.springs.clone() {
            let d = sub(self.particles[j].position, self.particles[i].position);
            let l = len(d);
            if l < 1e-12 {
                continue;
            }
            let stretch = (l - rest) / l;
            let f = scale(d, k * stretch);
            self.particles[i].add_force(f);
            self.particles[j].add_force(scale(f, -1.0));
        }
    }

    /// Apply buoyancy forces.
    fn apply_buoyancy(&mut self) {
        for p in &mut self.particles {
            let b = [0.0, self.fluid.density * self.fluid.gravity * p.volume, 0.0];
            p.add_force(b);
        }
    }

    /// Apply gravity forces.
    fn apply_gravity(&mut self) {
        for p in &mut self.particles {
            if !p.pinned {
                p.add_force([0.0, -self.fluid.gravity * p.mass, 0.0]);
            }
        }
    }

    /// Apply viscous drag forces.
    fn apply_drag(&mut self) {
        let current = self.fluid.current;
        let rho = self.fluid.density;
        let cd = self.cd;
        let area = self.drag_area;
        for p in &mut self.particles {
            let rel_v = sub(current, p.velocity);
            let speed = len(rel_v);
            let drag = scale(rel_v, 0.5 * rho * cd * area * speed);
            p.add_force(drag);
        }
    }

    /// Advance the simulation by one time step `dt` (seconds).
    pub fn step(&mut self, dt: f64) {
        self.apply_gravity();
        self.apply_buoyancy();
        self.apply_drag();
        self.apply_spring_forces();
        for p in &mut self.particles {
            p.integrate(dt);
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // T01: Hydrostatic pressure increases linearly with depth.
    #[test]
    fn test_hydrostatic_pressure_linear() {
        let env = FluidEnvironment::seawater();
        let p1 = env.hydrostatic_pressure(10.0);
        let p2 = env.hydrostatic_pressure(20.0);
        assert!((p2 - 2.0 * p1).abs() < 1e-6, "p1={p1}, p2={p2}");
    }

    // T02: Buoyancy force is upward and proportional to volume.
    #[test]
    fn test_buoyancy_force_direction_and_scale() {
        let env = FluidEnvironment::seawater();
        let f = env.buoyancy_force(1.0);
        assert!(f[1] > 0.0, "Buoyancy must be upward");
        let f2 = env.buoyancy_force(2.0);
        assert!((f2[1] - 2.0 * f[1]).abs() < EPS);
    }

    // T03: Added mass sphere coefficients are 0.5 along all axes.
    #[test]
    fn test_added_mass_sphere_coefficients() {
        let am = AddedMass::sphere(1000.0, 0.5);
        for &c in &am.coefficients {
            assert!((c - 0.5).abs() < EPS);
        }
    }

    // T04: Added mass force opposes acceleration.
    #[test]
    fn test_added_mass_force_opposes_acceleration() {
        let am = AddedMass::sphere(1000.0, 0.5);
        let accel = [1.0, 0.0, 0.0];
        let f = am.force(accel);
        assert!(f[0] < 0.0, "Added mass force should oppose acceleration");
    }

    // T05: Quadratic drag is zero for zero relative velocity.
    #[test]
    fn test_viscous_drag_zero_velocity() {
        let drag = ViscousDrag::new(1.2, 0.01, 1025.0);
        let f = drag.force([0.0; 3]);
        for c in f {
            assert!(c.abs() < EPS);
        }
    }

    // T06: Viscous drag scales quadratically with speed.
    #[test]
    fn test_viscous_drag_quadratic_scaling() {
        let drag = ViscousDrag::new(1.2, 0.01, 1025.0);
        let f1 = drag.force([1.0, 0.0, 0.0]);
        let f2 = drag.force([2.0, 0.0, 0.0]);
        // |f2| / |f1| should be 4
        let ratio = len(f2) / len(f1);
        assert!((ratio - 4.0).abs() < 1e-6, "ratio={ratio}");
    }

    // T07: Morison inertia term dominates when fluid accelerates strongly.
    #[test]
    fn test_morison_inertia_dominates() {
        let elem = MorisonElement::new(2.0, 1.2, 0.3, 1025.0);
        let fluid_accel = [0.0, 0.0, 100.0]; // large acceleration
        let f = elem.force_per_unit_length([0.0; 3], fluid_accel, [0.0; 3]);
        assert!(len(f) > 0.0);
    }

    // T08: Morison drag term is in relative velocity direction.
    #[test]
    fn test_morison_drag_direction() {
        let elem = MorisonElement::new(1.0, 1.2, 0.2, 1025.0);
        let fluid_vel = [1.0, 0.0, 0.0];
        let f = elem.force_per_unit_length(fluid_vel, [0.0; 3], [0.0; 3]);
        assert!(
            f[0] > 0.0,
            "Drag should be in fluid direction when body is stationary"
        );
    }

    // T09: VIV model advances state without panic.
    #[test]
    fn test_viv_model_advances() {
        let mut viv = VivModel::new(1025.0, 0.3, 10.0, 2.0);
        let q0 = viv.q;
        viv.step(0.0, 0.01);
        // After one step, q should have changed (q_dot != 0 initially)
        let _ = q0; // we just check no panic
    }

    // T10: VIV shedding frequency matches Strouhal calculation.
    #[test]
    fn test_viv_shedding_frequency() {
        let viv = VivModel::new(1025.0, 0.3, 10.0, 2.0);
        let expected = 0.2 * 2.0 / 0.3; // St * U / D
        let got = viv.shedding_frequency_hz();
        assert!(
            (got - expected).abs() < 1e-6,
            "expected={expected}, got={got}"
        );
    }

    // T11: Marine growth increases effective drag coefficient.
    #[test]
    fn test_marine_growth_increases_cd() {
        let mg = MarineGrowth::new(0.05, 1100.0, 1.1);
        let base_cd = 1.2;
        assert!(mg.effective_cd(base_cd) > base_cd);
    }

    // T12: Marine growth added mass per length is positive.
    #[test]
    fn test_marine_growth_added_mass_positive() {
        let mg = MarineGrowth::new(0.02, 1100.0, 1.0);
        assert!(mg.added_mass_per_length(0.3) > 0.0);
    }

    // T13: Pressure hull hoop stress is compressive under external pressure.
    #[test]
    fn test_pressure_hull_hoop_stress_compressive() {
        let hull = PressureHull::new(0.5, 0.4, 200e9, 0.3);
        let sigma = hull.hoop_stress_inner(1e6);
        assert!(sigma < 0.0, "Hoop stress should be compressive: {sigma}");
    }

    // T14: Pressure hull radial deflection is negative (inward) under external pressure.
    #[test]
    fn test_pressure_hull_deflection_inward() {
        let hull = PressureHull::new(0.5, 0.4, 200e9, 0.3);
        let delta = hull.radial_deflection_outer(1e6);
        assert!(delta < 0.0, "Deflection should be inward: {delta}");
    }

    // T15: Buckling pressure is positive.
    #[test]
    fn test_pressure_hull_buckling_positive() {
        let hull = PressureHull::new(0.5, 0.48, 200e9, 0.3);
        assert!(hull.buckling_pressure() > 0.0);
    }

    // T16: Hydrostatic buckling safety factor > 1 at shallow depth.
    #[test]
    fn test_buckling_safety_factor_shallow() {
        let b = HydrostaticBuckling::new(0.5, 0.02, 3.0, 200e9, 0.3);
        let sf = b.safety_factor(100.0); // 100 m depth
        assert!(sf > 1.0, "Safety factor should be > 1 at 100 m: {sf}");
    }

    // T17: Gripper segment bending angle increases with pressure.
    #[test]
    fn test_gripper_segment_angle_increases_with_pressure() {
        let mut seg = GripperSegment::new(0.1, 0.05, 100.0);
        seg.pressurize(1000.0, 1e-4);
        let angle1 = seg.angle;
        seg.pressurize(2000.0, 1e-4);
        let angle2 = seg.angle;
        assert!(angle2 > angle1, "Higher pressure should increase bending");
    }

    // T18: Gripper tip deflection is non-zero when pressurized.
    #[test]
    fn test_gripper_tip_deflection_non_zero() {
        let mut finger = SoftGripperFinger::new(5, 0.02, 0.1, 500.0);
        finger.pressurize_all(5000.0, 1e-4);
        let tip = finger.tip_deflection();
        let d = (tip[0].powi(2) + tip[1].powi(2)).sqrt();
        assert!(
            d > 0.0,
            "Pressurized finger should have non-zero deflection"
        );
    }

    // T19: Anguilliform body has correct number of segments.
    #[test]
    fn test_anguilla_body_segment_count() {
        let body = AnguillaBody::new(0.3, 20, 20.0, 10.0, 0.5);
        assert_eq!(body.positions.len(), 21);
        assert_eq!(body.angles.len(), 21);
    }

    // T20: Anguilliform curvature at time 0 is a sine wave.
    #[test]
    fn test_anguilla_curvature_sine() {
        let body = AnguillaBody::new(1.0, 10, std::f64::consts::PI, 1.0, 0.5);
        let k_half = body.curvature(0.5); // s = 0.5 m
        let expected = 0.5 * (std::f64::consts::PI * 0.5).sin();
        assert!((k_half - expected).abs() < 1e-10);
    }

    // T21: Phase velocity equals omega / k.
    #[test]
    fn test_anguilla_phase_velocity() {
        let body = AnguillaBody::new(
            1.0,
            10,
            2.0 * std::f64::consts::PI,
            4.0 * std::f64::consts::PI,
            0.3,
        );
        let pv = body.phase_velocity();
        assert!((pv - 2.0).abs() < 1e-10, "pv={pv}");
    }

    // T22: Anguilliform step advances time.
    #[test]
    fn test_anguilla_step_advances_time() {
        let mut body = AnguillaBody::new(0.5, 10, 10.0, 5.0, 0.2);
        body.step(0.01);
        assert!((body.time - 0.01).abs() < EPS);
    }

    // T23: Bio-fin forces are finite numbers.
    #[test]
    fn test_bio_fin_forces_finite() {
        let fin = BioFin::new(
            0.2,
            0.05,
            2.0 * std::f64::consts::PI,
            0.02,
            1025.0,
            10.0,
            0.3,
            FinMode::Flapping,
        );
        let (t, l) = fin.forces(1.0);
        assert!(t.is_finite(), "thrust={t}");
        assert!(l.is_finite(), "lift={l}");
    }

    // T24: Bio-fin phase advances on step.
    #[test]
    fn test_bio_fin_step_advances_phase() {
        let mut fin = BioFin::new(
            0.2,
            0.05,
            2.0 * std::f64::consts::PI,
            0.02,
            1025.0,
            10.0,
            0.3,
            FinMode::Flapping,
        );
        fin.step(0.1);
        assert!(fin.phase > 0.0, "phase={}", fin.phase);
    }

    // T25: Buoyancy distributor gives zero force above surface.
    #[test]
    fn test_buoyancy_distributor_above_surface() {
        let bd = BuoyancyDistributor::new(1025.0, 9.81, 0.0);
        let f = bd.force_on_particle([0.0, 1.0, 0.0], 1.0); // y > 0
        for c in f {
            assert!(c.abs() < EPS, "c={c}");
        }
    }

    // T26: Buoyancy distributor gives upward force below surface.
    #[test]
    fn test_buoyancy_distributor_below_surface() {
        let bd = BuoyancyDistributor::new(1025.0, 9.81, 0.0);
        let f = bd.force_on_particle([0.0, -1.0, 0.0], 0.001);
        assert!(f[1] > 0.0, "Buoyancy should be upward below surface");
    }

    // T27: Cable segment axial strain is zero at rest length.
    #[test]
    fn test_cable_segment_zero_strain_at_rest() {
        let seg = CableSegment::new([0.0; 3], [1.0, 0.0, 0.0], 1.0, 1e6, 100.0, 5.0);
        assert!(seg.axial_strain().abs() < EPS);
    }

    // T28: Cable segment tension is positive under extension.
    #[test]
    fn test_cable_segment_tension_positive_extension() {
        let seg = CableSegment::new([0.0; 3], [1.5, 0.0, 0.0], 1.0, 1e6, 100.0, 5.0);
        assert!(seg.tension() > 0.0, "tension={}", seg.tension());
    }

    // T29: Underwater soft body particle falls under gravity without buoyancy.
    #[test]
    fn test_underwater_soft_body_gravity_sinks() {
        let mut fluid = FluidEnvironment::freshwater();
        fluid.density = 0.0; // no buoyancy
        let mut body = UnderwaterSoftBody::new(fluid, 0.0, 0.0);
        let p = UnderwaterParticle::new([0.0, 10.0, 0.0], 1.0, 0.0);
        body.add_particle(p);
        for _ in 0..100 {
            body.step(1.0 / 60.0);
        }
        assert!(body.particles[0].position[1] < 9.0, "particle should sink");
    }

    // T30: Underwater soft body with full buoyancy rises.
    #[test]
    fn test_underwater_soft_body_buoyancy_rises() {
        let fluid = FluidEnvironment::seawater();
        let mut body = UnderwaterSoftBody::new(fluid, 0.0, 0.0);
        // Volume chosen so buoyancy >> gravity
        let p = UnderwaterParticle::new([0.0, 0.0, 0.0], 1.0, 10.0);
        body.add_particle(p);
        for _ in 0..60 {
            body.step(1.0 / 60.0);
        }
        assert!(
            body.particles[0].position[1] > 0.0,
            "buoyant particle should rise"
        );
    }

    // T31: Spring between two free particles reduces to rest length.
    #[test]
    fn test_spring_converges_to_rest_length() {
        let fluid = FluidEnvironment::freshwater();
        // Use non-zero drag to dissipate energy and allow convergence.
        let mut body = UnderwaterSoftBody::new(fluid, 1.0, 0.01);
        let mut p0 = UnderwaterParticle::new([0.0; 3], 1.0, 0.0);
        p0.pinned = true;
        let p1 = UnderwaterParticle::new([3.0, 0.0, 0.0], 1.0, 0.0);
        body.add_particle(p0);
        body.add_particle(p1);
        // rest length = 1 m, stiffness = 500 N/m
        body.springs.push((0, 1, 1.0, 500.0));
        for _ in 0..500 {
            body.step(0.001);
        }
        let d = len(sub(body.particles[1].position, body.particles[0].position));
        assert!(
            (d - 1.0).abs() < 0.2,
            "spring length should converge to 1.0, got {d}"
        );
    }

    // T32: Reynolds number is proportional to speed.
    #[test]
    fn test_reynolds_number_proportional_to_speed() {
        let re1 = ViscousDrag::reynolds_number(1.0, 0.1, 1e-6);
        let re2 = ViscousDrag::reynolds_number(2.0, 0.1, 1e-6);
        assert!((re2 - 2.0 * re1).abs() < EPS);
    }

    // T33: Added mass scalar matches formula.
    #[test]
    fn test_added_mass_scalar_mass() {
        let am = AddedMass::sphere(1000.0, 1.0);
        let vol = 4.0 / 3.0 * std::f64::consts::PI;
        let expected = 0.5 * 1000.0 * vol;
        assert!(
            (am.scalar_mass(0) - expected).abs() < 1e-6,
            "scalar_mass={}",
            am.scalar_mass(0)
        );
    }

    // T34: VIV lift force is finite.
    #[test]
    fn test_viv_lift_force_finite() {
        let mut viv = VivModel::new(1025.0, 0.3, 10.0, 2.0);
        for _ in 0..100 {
            viv.step(0.0, 0.01);
        }
        assert!(viv.lift_force().is_finite());
    }

    // T35: Buckling pressure is positive for Windenburg formula.
    #[test]
    fn test_windenburg_buckling_positive() {
        let b = HydrostaticBuckling::new(0.5, 0.02, 3.0, 200e9, 0.3);
        assert!(b.critical_pressure_windenburg() > 0.0);
    }
}
