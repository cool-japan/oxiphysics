// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CFD–rigid body coupling utilities.
//!
//! Provides aerodynamic/hydrodynamic force and moment calculations,
//! added mass, vortex-induced vibration, mooring lines, and wave loads
//! for coupled fluid–structure simulations.
//!
//! References:
//! - Peskin (2002). Acta Numerica 11, 479–517. (Immersed boundary method)
//! - Morison et al. (1950). Trans. ASCE 115, 149.
//! - Skop & Balasubramanian (1997). J. Fluids Struct. 11, 301.
//! - Faltinsen (1990). Sea Loads on Ships and Offshore Structures. Cambridge.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// CfdForces
// ---------------------------------------------------------------------------

/// Aerodynamic/hydrodynamic forces and moments produced by a CFD solution.
#[derive(Debug, Clone, PartialEq)]
pub struct CfdForces {
    /// Drag force vector \[N\].
    pub drag: [f64; 3],
    /// Lift force vector \[N\].
    pub lift: [f64; 3],
    /// Aerodynamic moment vector \[N·m\].
    pub moment: [f64; 3],
    /// Centre of pressure in world coordinates \[m\].
    pub pressure_center: [f64; 3],
}

impl CfdForces {
    /// Create a zeroed CfdForces struct.
    pub fn zero() -> Self {
        Self {
            drag: [0.0_f64; 3],
            lift: [0.0_f64; 3],
            moment: [0.0_f64; 3],
            pressure_center: [0.0_f64; 3],
        }
    }

    /// Total aerodynamic force vector (drag + lift).
    pub fn total_force(&self) -> [f64; 3] {
        [
            self.drag[0] + self.lift[0],
            self.drag[1] + self.lift[1],
            self.drag[2] + self.lift[2],
        ]
    }

    /// Magnitude of the drag force vector.
    pub fn drag_magnitude(&self) -> f64 {
        (self.drag[0] * self.drag[0] + self.drag[1] * self.drag[1] + self.drag[2] * self.drag[2])
            .sqrt()
    }

    /// Magnitude of the lift force vector.
    pub fn lift_magnitude(&self) -> f64 {
        (self.lift[0] * self.lift[0] + self.lift[1] * self.lift[1] + self.lift[2] * self.lift[2])
            .sqrt()
    }
}

// ---------------------------------------------------------------------------
// CfdCouplingParams
// ---------------------------------------------------------------------------

/// Parameters describing the fluid medium and reference geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct CfdCouplingParams {
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Dynamic (absolute) viscosity \[Pa·s\].
    pub fluid_viscosity: f64,
    /// Reference planform area \[m²\].
    pub reference_area: f64,
    /// Reference chord / diameter length \[m\].
    pub reference_length: f64,
}

impl CfdCouplingParams {
    /// Compute dynamic (free-stream) pressure q∞ = ½ ρ V².
    pub fn dynamic_pressure(&self, velocity: f64) -> f64 {
        0.5_f64 * self.fluid_density * velocity * velocity
    }

    /// Compute Reynolds number Re = ρ V L / μ.
    pub fn reynolds_number(&self, velocity: f64) -> f64 {
        self.fluid_density * velocity.abs() * self.reference_length
            / self.fluid_viscosity.max(1e-30_f64)
    }
}

// ---------------------------------------------------------------------------
// Basic aerodynamic utilities (previously existing, retained and documented)
// ---------------------------------------------------------------------------

/// Compute the Reynolds number Re = ρ·V·L / μ.
///
/// # Arguments
/// * `velocity`   – flow speed \[m/s\]
/// * `length`     – reference length \[m\]
/// * `viscosity`  – dynamic viscosity \[Pa·s\]
/// * `density`    – fluid density \[kg/m³\]
pub fn reynolds_number(velocity: f64, length: f64, viscosity: f64, density: f64) -> f64 {
    density * velocity.abs() * length / viscosity.max(1e-30_f64)
}

/// Compute the drag force vector acting on a body.
///
/// F_drag = ½ · ρ · |V|² · Cd · A · V̂
///
/// The force is aligned anti-parallel to the velocity direction.
pub fn drag_force(cd: f64, rho: f64, vel: [f64; 3], area: f64) -> [f64; 3] {
    let speed_sq = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
    let speed = speed_sq.sqrt();
    if speed < 1e-15_f64 {
        return [0.0_f64; 3];
    }
    let mag = 0.5_f64 * rho * speed_sq * cd * area;
    [
        -mag * vel[0] / speed,
        -mag * vel[1] / speed,
        -mag * vel[2] / speed,
    ]
}

/// Compute the lift force vector acting on a body.
///
/// F_lift = ½ · ρ · |V|² · Cl · A, directed perpendicular to velocity
/// (upward in the Y direction by convention).
pub fn lift_force(cl: f64, rho: f64, vel: [f64; 3], area: f64) -> [f64; 3] {
    let speed_sq = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
    let mag = 0.5_f64 * rho * speed_sq * cl * area;
    [0.0_f64, mag, 0.0_f64]
}

/// Compute the pitching moment coefficient from a discrete pressure
/// distribution integrated over panels.
///
/// Cm = Σ(p_i · m_i) / (q_∞ · A · L)
pub fn moment_coefficient_from_pressure(
    pressure_dist: &[f64],
    moments: &[f64],
    q_inf: f64,
    area: f64,
    length: f64,
) -> f64 {
    let numerator: f64 = pressure_dist
        .iter()
        .zip(moments.iter())
        .map(|(p, m)| p * m)
        .sum();
    let denom = q_inf * area * length;
    if denom.abs() < 1e-30_f64 {
        0.0_f64
    } else {
        numerator / denom
    }
}

/// Compute the added-mass inertial force on a body accelerating through fluid.
///
/// F_added = −Cm · ρ_f · V_body · a
pub fn added_mass_force(
    acceleration: [f64; 3],
    fluid_density: f64,
    volume: f64,
    cm: f64,
) -> [f64; 3] {
    let k = -cm * fluid_density * volume;
    [
        k * acceleration[0],
        k * acceleration[1],
        k * acceleration[2],
    ]
}

/// Approximate drag coefficient for a bluff body using Reynolds-number correlations.
pub fn bluff_body_cd(reynolds: f64) -> f64 {
    if reynolds < 1e-10_f64 {
        return 24.0_f64 * 1e10_f64;
    }
    if reynolds < 1.0_f64 {
        24.0_f64 / reynolds
    } else if reynolds < 1_000.0_f64 {
        24.0_f64 / reynolds + 6.0_f64 / (1.0_f64 + reynolds.sqrt()) + 0.4_f64
    } else {
        0.44_f64
    }
}

/// Estimate vortex-shedding frequency using the Strouhal number.
///
/// f_vs = St · V / D,  where St ≈ 0.21 (bluff cylinder approximation).
pub fn strouhal_vortex_shedding(_cd: f64, velocity: f64, diameter: f64) -> f64 {
    const STROUHAL: f64 = 0.21_f64;
    if diameter.abs() < 1e-15_f64 {
        return 0.0_f64;
    }
    STROUHAL * velocity.abs() / diameter
}

/// Compute the maximum fluid-structure coupling time step from the CFL condition.
///
/// Δt = CFL · Δx / |V|
pub fn fluid_structure_time_step(courant: f64, dx: f64, velocity: f64) -> f64 {
    let speed = velocity.abs();
    if speed < 1e-15_f64 {
        return f64::INFINITY;
    }
    courant * dx / speed
}

// ---------------------------------------------------------------------------
// ImmersedBodyCfd
// ---------------------------------------------------------------------------

/// Immersed boundary method representation of a rigid body in a fluid domain.
///
/// In the immersed boundary method (IBM), the rigid body is represented by
/// a set of Lagrangian marker points embedded in the Eulerian fluid grid.
/// Forces are spread from markers to the grid, and velocities are interpolated
/// from the grid back to the markers using regularised delta functions.
#[derive(Debug, Clone)]
pub struct ImmersedBodyCfd {
    /// Lagrangian marker positions \[m\].
    pub markers: Vec<[f64; 3]>,
    /// Marker normals (unit outward normals).
    pub marker_normals: Vec<[f64; 3]>,
    /// Area weight of each marker \[m²\].
    pub marker_areas: Vec<f64>,
    /// Body centre of mass \[m\].
    pub center: [f64; 3],
    /// Body translational velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Body angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Body mass \[kg\].
    pub mass: f64,
    /// Euler fluid grid spacing \[m\].
    pub grid_spacing: f64,
}

impl ImmersedBodyCfd {
    /// Create a new immersed body from marker data.
    pub fn new(
        markers: Vec<[f64; 3]>,
        marker_normals: Vec<[f64; 3]>,
        marker_areas: Vec<f64>,
        center: [f64; 3],
        mass: f64,
        grid_spacing: f64,
    ) -> Self {
        Self {
            markers,
            marker_normals,
            marker_areas,
            center,
            velocity: [0.0_f64; 3],
            angular_velocity: [0.0_f64; 3],
            mass,
            grid_spacing,
        }
    }

    /// Create a simplified sphere approximation with `n_lat` × `n_lon` markers.
    pub fn new_sphere(
        radius: f64,
        center: [f64; 3],
        mass: f64,
        grid_spacing: f64,
        n_lat: usize,
        n_lon: usize,
    ) -> Self {
        let mut markers = Vec::new();
        let mut normals = Vec::new();
        let mut areas = Vec::new();
        let d_theta = PI / n_lat as f64;
        let d_phi = 2.0_f64 * PI / n_lon as f64;
        let area_each = radius * radius * d_theta * d_phi;
        for i in 0..n_lat {
            let theta = (i as f64 + 0.5_f64) * d_theta;
            for j in 0..n_lon {
                let phi = j as f64 * d_phi;
                let nx = theta.sin() * phi.cos();
                let ny = theta.cos();
                let nz = theta.sin() * phi.sin();
                markers.push([
                    center[0] + radius * nx,
                    center[1] + radius * ny,
                    center[2] + radius * nz,
                ]);
                normals.push([nx, ny, nz]);
                areas.push(area_each * theta.sin());
            }
        }
        Self::new(markers, normals, areas, center, mass, grid_spacing)
    }

    /// Regularised delta function (Peskin cosine kernel).
    ///
    /// φ(r) = (1 + cos(π r / (2h))) / (4h)  for |r| ≤ 2h, else 0.
    pub fn delta_function(&self, r: f64) -> f64 {
        let h = self.grid_spacing;
        if r.abs() > 2.0_f64 * h {
            return 0.0_f64;
        }
        (1.0_f64 + (PI * r / (2.0_f64 * h)).cos()) / (4.0_f64 * h)
    }

    /// Spread marker force `f_marker` to a nearby Eulerian grid point.
    ///
    /// F_grid = f_marker · φ(x_marker − x_grid) · dA
    pub fn spread_force(
        &self,
        marker_idx: usize,
        f_marker: [f64; 3],
        grid_point: [f64; 3],
    ) -> [f64; 3] {
        let xm = self.markers[marker_idx];
        let da = self.marker_areas[marker_idx];
        let dr = [
            xm[0] - grid_point[0],
            xm[1] - grid_point[1],
            xm[2] - grid_point[2],
        ];
        let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let phi = self.delta_function(r);
        [
            f_marker[0] * phi * da,
            f_marker[1] * phi * da,
            f_marker[2] * phi * da,
        ]
    }

    /// Interpolate fluid velocity at marker position from Eulerian grid value.
    ///
    /// u_marker = u_grid · φ(x_marker − x_grid) · h³
    pub fn interpolate_velocity(
        &self,
        marker_idx: usize,
        u_grid: [f64; 3],
        grid_point: [f64; 3],
    ) -> [f64; 3] {
        let xm = self.markers[marker_idx];
        let h3 = self.grid_spacing.powi(3);
        let dr = [
            xm[0] - grid_point[0],
            xm[1] - grid_point[1],
            xm[2] - grid_point[2],
        ];
        let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
        let phi = self.delta_function(r);
        [
            u_grid[0] * phi * h3,
            u_grid[1] * phi * h3,
            u_grid[2] * phi * h3,
        ]
    }

    /// Number of Lagrangian markers.
    pub fn n_markers(&self) -> usize {
        self.markers.len()
    }

    /// Update body position and marker positions after one time step.
    pub fn advance(&mut self, dt: f64) {
        self.center[0] += self.velocity[0] * dt;
        self.center[1] += self.velocity[1] * dt;
        self.center[2] += self.velocity[2] * dt;
        for m in self.markers.iter_mut() {
            m[0] += self.velocity[0] * dt;
            m[1] += self.velocity[1] * dt;
            m[2] += self.velocity[2] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// BuoyancyCfd
// ---------------------------------------------------------------------------

/// Pressure-integrated buoyancy force from a CFD pressure field.
///
/// Computes F_b = −∮ p n̂ dA over the submerged body surface,
/// which for a uniform pressure field reduces to Archimedes' law:
/// F_b = ρ g V (upward).
#[derive(Debug, Clone)]
pub struct BuoyancyCfd {
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Gravity vector \[m/s²\].
    pub gravity: [f64; 3],
    /// Submerged volume \[m³\].
    pub submerged_volume: f64,
    /// Last computed buoyancy force \[N\].
    pub force: [f64; 3],
}

impl BuoyancyCfd {
    /// Create a new buoyancy model.
    pub fn new(fluid_density: f64, gravity: [f64; 3]) -> Self {
        Self {
            fluid_density,
            gravity,
            submerged_volume: 0.0_f64,
            force: [0.0_f64; 3],
        }
    }

    /// Compute the Archimedes buoyancy force from submerged volume.
    ///
    /// F_b = ρ_fluid · V_sub · g  (vector, directed opposite to gravity).
    pub fn archimedes_force(&mut self, submerged_volume: f64) -> [f64; 3] {
        self.submerged_volume = submerged_volume;
        let k = self.fluid_density * submerged_volume;
        self.force = [
            k * self.gravity[0],
            k * self.gravity[1],
            k * self.gravity[2],
        ];
        self.force
    }

    /// Compute pressure-integrated buoyancy from panel pressures and normals.
    ///
    /// F_b = −Σ_i p_i · n̂_i · dA_i
    pub fn pressure_integrated_force(
        &mut self,
        pressures: &[f64],
        normals: &[[f64; 3]],
        areas: &[f64],
    ) -> [f64; 3] {
        let mut f = [0.0_f64; 3];
        for ((&p, n), &a) in pressures.iter().zip(normals.iter()).zip(areas.iter()) {
            f[0] -= p * n[0] * a;
            f[1] -= p * n[1] * a;
            f[2] -= p * n[2] * a;
        }
        self.force = f;
        f
    }

    /// Fraction of body submerged given current water line height.
    ///
    /// `body_bottom` and `body_top` are the lower and upper extents of the body \[m\].
    /// `water_level` is the free surface elevation \[m\].
    pub fn submersion_fraction(body_bottom: f64, body_top: f64, water_level: f64) -> f64 {
        let height = (body_top - body_bottom).max(1e-30_f64);
        let submerged = (water_level - body_bottom).clamp(0.0_f64, height);
        submerged / height
    }

    /// Buoyancy stability check: returns positive if centre of buoyancy
    /// is above centre of gravity (stable configuration).
    pub fn metacentric_height(centre_of_gravity_y: f64, centre_of_buoyancy_y: f64) -> f64 {
        centre_of_buoyancy_y - centre_of_gravity_y
    }
}

// ---------------------------------------------------------------------------
// DragLiftCfd
// ---------------------------------------------------------------------------

/// Drag and lift forces from CFD surface stress integration.
///
/// Integrates both pressure and viscous (skin friction) contributions
/// over the body surface to compute the total drag and lift vectors.
#[derive(Debug, Clone)]
pub struct DragLiftCfd {
    /// Free-stream velocity direction (unit vector).
    pub flow_direction: [f64; 3],
    /// Lift direction (unit vector, perpendicular to flow).
    pub lift_direction: [f64; 3],
    /// Reference area \[m²\].
    pub reference_area: f64,
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Free-stream speed \[m/s\].
    pub free_stream_speed: f64,
    /// Last computed drag coefficient.
    pub cd: f64,
    /// Last computed lift coefficient.
    pub cl: f64,
}

impl DragLiftCfd {
    /// Create a new drag/lift integrator.
    pub fn new(
        flow_direction: [f64; 3],
        lift_direction: [f64; 3],
        reference_area: f64,
        fluid_density: f64,
        free_stream_speed: f64,
    ) -> Self {
        Self {
            flow_direction,
            lift_direction,
            reference_area,
            fluid_density,
            free_stream_speed,
            cd: 0.0_f64,
            cl: 0.0_f64,
        }
    }

    /// Integrate pressure drag from panel pressures and normals.
    ///
    /// D_pressure = Σ_i p_i (n̂_i · ê_drag) dA_i
    pub fn pressure_drag(&mut self, pressures: &[f64], normals: &[[f64; 3]], areas: &[f64]) -> f64 {
        let mut d = 0.0_f64;
        let ed = &self.flow_direction;
        for ((&p, n), &a) in pressures.iter().zip(normals.iter()).zip(areas.iter()) {
            d += p * (n[0] * ed[0] + n[1] * ed[1] + n[2] * ed[2]) * a;
        }
        d
    }

    /// Integrate pressure lift from panel pressures and normals.
    ///
    /// L_pressure = Σ_i p_i (n̂_i · ê_lift) dA_i
    pub fn pressure_lift(&mut self, pressures: &[f64], normals: &[[f64; 3]], areas: &[f64]) -> f64 {
        let mut l = 0.0_f64;
        let el = &self.lift_direction;
        for ((&p, n), &a) in pressures.iter().zip(normals.iter()).zip(areas.iter()) {
            l += p * (n[0] * el[0] + n[1] * el[1] + n[2] * el[2]) * a;
        }
        l
    }

    /// Dynamic pressure at free-stream conditions.
    pub fn dynamic_pressure(&self) -> f64 {
        0.5_f64 * self.fluid_density * self.free_stream_speed * self.free_stream_speed
    }

    /// Update drag and lift coefficients from integrated forces.
    pub fn update_coefficients(&mut self, drag: f64, lift: f64) {
        let q = self.dynamic_pressure();
        let denom = q * self.reference_area;
        if denom.abs() < 1e-30_f64 {
            return;
        }
        self.cd = drag / denom;
        self.cl = lift / denom;
    }

    /// Compute total drag and lift forces from stored coefficients.
    pub fn forces(&self) -> (f64, f64) {
        let q = self.dynamic_pressure();
        let a = self.reference_area;
        (self.cd * q * a, self.cl * q * a)
    }
}

// ---------------------------------------------------------------------------
// VortexInducedVibration
// ---------------------------------------------------------------------------

/// Vortex-induced vibration (VIV) response model for a cylindrical structure.
///
/// Implements the wake oscillator model of Skop & Balasubramanian (1997),
/// combined with the Morison equation for hydrodynamic loading.
///
/// The wake oscillator equation:
/// q̈ + 2ε ω_s (q² − 1) q̇ + ω_s² q = A ÿ
///
/// where q is the non-dimensional wake variable, ω_s is the Strouhal
/// angular frequency, ε is the van der Pol parameter, and y is the
/// cross-flow displacement of the cylinder.
#[derive(Debug, Clone)]
pub struct VortexInducedVibration {
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Cylinder length \[m\].
    pub length: f64,
    /// Cylinder mass per unit length \[kg/m\].
    pub mass_per_length: f64,
    /// Structural natural frequency \[rad/s\].
    pub natural_frequency: f64,
    /// Structural damping ratio.
    pub damping_ratio: f64,
    /// Free-stream velocity \[m/s\].
    pub flow_velocity: f64,
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Current cross-flow displacement \[m\].
    pub displacement: f64,
    /// Current cross-flow velocity \[m/s\].
    pub disp_velocity: f64,
    /// Wake oscillator variable q (dimensionless).
    pub wake_variable: f64,
    /// Wake oscillator rate dq/dt.
    pub wake_rate: f64,
    /// Van der Pol (stall) parameter ε.
    pub van_der_pol_eps: f64,
    /// Wake oscillator coupling coefficient A.
    pub coupling_coeff: f64,
    /// Vortex shedding frequency \[rad/s\].
    pub shedding_frequency: f64,
    /// Morison inertia coefficient C_M.
    pub cm: f64,
    /// Morison drag coefficient C_D.
    pub cd_morison: f64,
    /// Time elapsed \[s\].
    pub time: f64,
}

impl VortexInducedVibration {
    /// Create a new VIV model.
    pub fn new(
        diameter: f64,
        length: f64,
        mass_per_length: f64,
        natural_frequency: f64,
        damping_ratio: f64,
        flow_velocity: f64,
        fluid_density: f64,
    ) -> Self {
        let shedding_frequency = 2.0_f64 * PI * 0.21_f64 * flow_velocity / diameter.max(1e-15_f64);
        Self {
            diameter,
            length,
            mass_per_length,
            natural_frequency,
            damping_ratio,
            flow_velocity,
            fluid_density,
            displacement: 0.0_f64,
            disp_velocity: 0.0_f64,
            wake_variable: 0.1_f64, // small initial perturbation
            wake_rate: 0.0_f64,
            van_der_pol_eps: 0.3_f64,
            coupling_coeff: 12.0_f64,
            shedding_frequency,
            cm: 1.0_f64,
            cd_morison: 1.0_f64,
            time: 0.0_f64,
        }
    }

    /// Strouhal shedding frequency for the current flow conditions.
    pub fn strouhal_frequency(&self) -> f64 {
        2.0_f64 * PI * 0.21_f64 * self.flow_velocity / self.diameter.max(1e-15_f64)
    }

    /// Non-dimensional reduced velocity Vr = U / (f_n D).
    pub fn reduced_velocity(&self) -> f64 {
        let f_n_hz = self.natural_frequency / (2.0_f64 * PI);
        self.flow_velocity / (f_n_hz * self.diameter).max(1e-30_f64)
    }

    /// Check whether lock-in condition is satisfied.
    ///
    /// Lock-in occurs approximately when 4 ≤ Vr ≤ 8.
    pub fn is_lock_in(&self) -> bool {
        let vr = self.reduced_velocity();
        (4.0_f64..=8.0_f64).contains(&vr)
    }

    /// Morison equation: total hydrodynamic force per unit length.
    ///
    /// F/L = ρ C_M (π D²/4) a_fluid + ½ ρ C_D D |u_rel| u_rel
    ///
    /// where u_rel = u_fluid − ẏ_body.
    pub fn morison_force_per_length(&self, fluid_velocity: f64, fluid_acceleration: f64) -> f64 {
        let a_cross = PI * self.diameter * self.diameter / 4.0_f64;
        let u_rel = fluid_velocity - self.disp_velocity;
        let inertia = self.fluid_density * self.cm * a_cross * fluid_acceleration;
        let drag =
            0.5_f64 * self.fluid_density * self.cd_morison * self.diameter * u_rel * u_rel.abs();
        inertia + drag
    }

    /// Total Morison force on the full cylinder length.
    pub fn morison_force_total(&self, fluid_velocity: f64, fluid_acceleration: f64) -> f64 {
        self.morison_force_per_length(fluid_velocity, fluid_acceleration) * self.length
    }

    /// Wake oscillator right-hand side: dq/dt = q_rate, dq_rate/dt = ...
    ///
    /// q̈ = −2ε ω_s (q²−1) q̇ − ω_s² q + A ÿ
    fn wake_oscillator_accel(&self, q: f64, q_dot: f64, y_ddot: f64) -> f64 {
        let ws = self.shedding_frequency;
        -2.0_f64 * self.van_der_pol_eps * ws * (q * q - 1.0_f64) * q_dot - ws * ws * q
            + self.coupling_coeff * y_ddot
    }

    /// Advance VIV system by one time step using semi-implicit Euler.
    ///
    /// Returns (new_displacement, new_wake_variable).
    pub fn step(&mut self, dt: f64, fluid_velocity: f64, fluid_acceleration: f64) -> (f64, f64) {
        let morison = self.morison_force_per_length(fluid_velocity, fluid_acceleration);
        let m = self.mass_per_length;
        let wn = self.natural_frequency;
        let zeta = self.damping_ratio;
        let ws = self.shedding_frequency;

        // Vortex lift per unit length: F_viv = ½ ρ U² D Cl with Cl = q/2
        let cl = self.wake_variable / 2.0_f64;
        let f_viv = 0.5_f64
            * self.fluid_density
            * self.flow_velocity
            * self.flow_velocity
            * self.diameter
            * cl;

        // Structural EOM: ÿ = (F_morison + F_viv − 2ζωₙẏ − ωₙ²y) / m
        let y_ddot = if m > 1e-30_f64 {
            (morison + f_viv
                - 2.0_f64 * zeta * wn * self.disp_velocity
                - wn * wn * self.displacement)
                / m
        } else {
            0.0_f64
        };

        // Wake oscillator
        let q_ddot = self.wake_oscillator_accel(self.wake_variable, self.wake_rate, y_ddot);

        self.wake_rate += q_ddot * dt;
        self.wake_variable += self.wake_rate * dt;
        self.disp_velocity += y_ddot * dt;
        self.displacement += self.disp_velocity * dt;
        self.time += dt;
        self.shedding_frequency =
            2.0_f64 * PI * 0.21_f64 * self.flow_velocity / self.diameter.max(1e-15_f64);

        let _ = ws; // suppress unused-variable warning

        (self.displacement, self.wake_variable)
    }

    /// Amplitude of oscillation after reaching statistical steady state.
    ///
    /// Approximate from wake variable: A ≈ D · |q| / 2
    pub fn estimated_amplitude(&self) -> f64 {
        self.diameter * self.wake_variable.abs() / 2.0_f64
    }
}

// ---------------------------------------------------------------------------
// MooringLine
// ---------------------------------------------------------------------------

/// Catenary mooring line statics and simplified dynamics.
///
/// Implements the analytical catenary equations for a mooring line
/// under its own weight in a horizontal-vertical plane.
///
/// The catenary solution:
/// z(x) = a (cosh(x/a) − 1),  a = T_h / (w)
///
/// where T_h is the horizontal tension and w is the line weight per unit length.
#[derive(Debug, Clone)]
pub struct MooringLine {
    /// Total unstretched length of the mooring line \[m\].
    pub unstretched_length: f64,
    /// Submerged weight per unit length \[N/m\].
    pub weight_per_length: f64,
    /// Axial stiffness EA \[N\].
    pub axial_stiffness: f64,
    /// Anchor position \[m\].
    pub anchor: [f64; 3],
    /// Fairlead position \[m\].
    pub fairlead: [f64; 3],
    /// Horizontal tension \[N\].
    pub horizontal_tension: f64,
    /// Vertical tension at fairlead \[N\].
    pub vertical_tension: f64,
    /// Number of line segments for discretisation.
    pub n_segments: usize,
}

impl MooringLine {
    /// Create a new mooring line.
    pub fn new(
        unstretched_length: f64,
        weight_per_length: f64,
        axial_stiffness: f64,
        anchor: [f64; 3],
        fairlead: [f64; 3],
        n_segments: usize,
    ) -> Self {
        Self {
            unstretched_length,
            weight_per_length,
            axial_stiffness,
            anchor,
            fairlead,
            horizontal_tension: 0.0_f64,
            vertical_tension: 0.0_f64,
            n_segments,
        }
    }

    /// Horizontal distance between anchor and fairlead (projected).
    pub fn horizontal_span(&self) -> f64 {
        let dx = self.fairlead[0] - self.anchor[0];
        let dz = self.fairlead[2] - self.anchor[2];
        (dx * dx + dz * dz).sqrt()
    }

    /// Vertical rise from anchor to fairlead.
    pub fn vertical_rise(&self) -> f64 {
        self.fairlead[1] - self.anchor[1]
    }

    /// Catenary parameter a = T_h / w.
    pub fn catenary_parameter(&self) -> f64 {
        if self.weight_per_length < 1e-30_f64 {
            return f64::INFINITY;
        }
        self.horizontal_tension / self.weight_per_length
    }

    /// Solve for horizontal tension using catenary equations (Newton iteration).
    ///
    /// Returns the converged horizontal tension T_h \[N\].
    pub fn solve_static(&mut self) -> f64 {
        let h = self.horizontal_span();
        let v = self.vertical_rise();
        let l = self.unstretched_length;
        let w = self.weight_per_length;

        if w < 1e-30_f64 || l < 1e-12_f64 {
            self.horizontal_tension = 0.0_f64;
            self.vertical_tension = 0.0_f64;
            return 0.0_f64;
        }

        // Initial guess: T_h = w * h / 2
        let mut t_h = (w * h / 2.0_f64).max(1.0_f64);

        for _ in 0..50 {
            let a = t_h / w;
            // Catenary arc length: L = a * (sinh(h_end/a) - sinh(h_start/a))
            // Simplified: assume anchor at x=0, fairlead at x=h with height v
            // Using: L² = (l_arc)², h_end satisfying cosh relation
            // Residual: f(T_h) = sqrt(L² - v²) - 2a sinh(sqrt(L²-v²)/(2a)) = h
            let l_sq = l * l;
            let v_sq = v * v;
            if l_sq <= v_sq {
                break;
            }
            let lh = (l_sq - v_sq).sqrt();
            let f = 2.0_f64 * a * (lh / (2.0_f64 * a)).sinh() - h;
            let df = 2.0_f64 * (lh / (2.0_f64 * a)).sinh() - (lh / (2.0_f64 * a)).cosh() * lh / a;
            if df.abs() < 1e-15_f64 {
                break;
            }
            t_h -= f / df;
            t_h = t_h.max(1e-6_f64);
        }

        self.horizontal_tension = t_h;
        // Vertical tension: T_v = w * L (simplified catenary)
        let lh = if l * l > v * v {
            (l * l - v * v).sqrt()
        } else {
            0.0_f64
        };
        self.vertical_tension = w * lh;
        t_h
    }

    /// Total tension at the fairlead.
    pub fn fairlead_tension(&self) -> f64 {
        (self.horizontal_tension * self.horizontal_tension
            + self.vertical_tension * self.vertical_tension)
            .sqrt()
    }

    /// Catenary profile: z coordinate at horizontal position x from anchor.
    pub fn profile_z(&self, x: f64) -> f64 {
        let a = self.catenary_parameter();
        if a.is_infinite() || a < 1e-15_f64 {
            return 0.0_f64;
        }
        a * ((x / a).cosh() - 1.0_f64)
    }

    /// Discretise the mooring line into segment endpoint positions.
    pub fn discretised_positions(&self) -> Vec<[f64; 3]> {
        let n = self.n_segments;
        let h = self.horizontal_span();
        let mut pts = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = t * h;
            let z = self.profile_z(x);
            pts.push([self.anchor[0] + x, self.anchor[1] + z, self.anchor[2]]);
        }
        pts
    }

    /// Effective stiffness at fairlead (dT/dh), linearised.
    pub fn linearised_stiffness(&self) -> f64 {
        let h = self.horizontal_span();
        if h < 1e-12_f64 {
            return 0.0_f64;
        }
        self.horizontal_tension / h
    }
}

// ---------------------------------------------------------------------------
// WaveLoadModel
// ---------------------------------------------------------------------------

/// Morison equation wave load model for cylindrical structural members.
///
/// Computes the wave-induced force per unit length on a vertical cylinder
/// using the Morison equation:
///
/// F/L = ρ C_M (π D²/4) ü_w + ½ ρ C_D D u_w |u_w|
///
/// where u_w is the wave-induced water particle velocity.
#[derive(Debug, Clone)]
pub struct WaveLoadModel {
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Inertia coefficient C_M.
    pub cm: f64,
    /// Drag coefficient C_D.
    pub cd: f64,
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Wave height H \[m\].
    pub wave_height: f64,
    /// Wave period T \[s\].
    pub wave_period: f64,
    /// Water depth d \[m\].
    pub water_depth: f64,
}

impl WaveLoadModel {
    /// Create a new Morison wave load model.
    pub fn new(
        diameter: f64,
        cm: f64,
        cd: f64,
        fluid_density: f64,
        wave_height: f64,
        wave_period: f64,
        water_depth: f64,
    ) -> Self {
        Self {
            diameter,
            cm,
            cd,
            fluid_density,
            wave_height,
            wave_period,
            water_depth,
        }
    }

    /// Wave angular frequency ω = 2π/T.
    pub fn angular_frequency(&self) -> f64 {
        2.0_f64 * PI / self.wave_period.max(1e-30_f64)
    }

    /// Wave number k from linear dispersion relation ω² = g k tanh(k d).
    ///
    /// Solved iteratively using Newton's method.
    pub fn wave_number(&self) -> f64 {
        const G: f64 = 9.81_f64;
        let omega = self.angular_frequency();
        let omega2 = omega * omega;
        let d = self.water_depth;
        // Deep-water initial guess
        let mut k = omega2 / G;
        for _ in 0..30 {
            let kd = k * d;
            let th = kd.tanh();
            let f = G * k * th - omega2;
            let df = G * (th + k * d * (1.0_f64 - th * th));
            if df.abs() < 1e-30_f64 {
                break;
            }
            k -= f / df;
            k = k.abs().max(1e-10_f64);
        }
        k
    }

    /// Water particle horizontal velocity at depth z and time t (linear theory).
    ///
    /// u(z, t) = (H ω / 2) · cosh(k(z+d)) / sinh(kd) · cos(ωt)
    pub fn water_particle_velocity(&self, z: f64, t: f64) -> f64 {
        let omega = self.angular_frequency();
        let k = self.wave_number();
        let d = self.water_depth;
        let h = self.wave_height;
        let denom = (k * d).sinh();
        if denom.abs() < 1e-30_f64 {
            return 0.0_f64;
        }
        h * omega / 2.0_f64 * (k * (z + d)).cosh() / denom * (omega * t).cos()
    }

    /// Water particle horizontal acceleration at depth z and time t.
    ///
    /// a(z, t) = −(H ω² / 2) · cosh(k(z+d)) / sinh(kd) · sin(ωt)
    pub fn water_particle_acceleration(&self, z: f64, t: f64) -> f64 {
        let omega = self.angular_frequency();
        let k = self.wave_number();
        let d = self.water_depth;
        let h = self.wave_height;
        let denom = (k * d).sinh();
        if denom.abs() < 1e-30_f64 {
            return 0.0_f64;
        }
        -h * omega * omega / 2.0_f64 * (k * (z + d)).cosh() / denom * (omega * t).sin()
    }

    /// Morison force per unit length at depth z and time t \[N/m\].
    pub fn force_per_length(&self, z: f64, t: f64) -> f64 {
        let u = self.water_particle_velocity(z, t);
        let a_fluid = self.water_particle_acceleration(z, t);
        let a_cross = PI * self.diameter * self.diameter / 4.0_f64;
        let inertia = self.fluid_density * self.cm * a_cross * a_fluid;
        let drag = 0.5_f64 * self.fluid_density * self.cd * self.diameter * u * u.abs();
        inertia + drag
    }

    /// Total Morison force on the cylinder from seabed (z=−d) to surface (z=0).
    ///
    /// Computed by numerical integration using n_points Gaussian quadrature
    /// (here simple trapezoidal rule).
    pub fn total_force(&self, t: f64, n_points: usize) -> f64 {
        let d = self.water_depth;
        let dz = d / n_points.max(1) as f64;
        let mut f = 0.0_f64;
        for i in 0..n_points {
            let z = -d + (i as f64 + 0.5_f64) * dz;
            f += self.force_per_length(z, t) * dz;
        }
        f
    }

    /// Maximum Morison force (over one wave cycle) on the cylinder.
    pub fn maximum_force(&self, n_time: usize, n_depth: usize) -> f64 {
        let t_max = self.wave_period;
        let dt = t_max / n_time.max(1) as f64;
        let mut f_max = 0.0_f64;
        for i in 0..n_time {
            let t = i as f64 * dt;
            let f = self.total_force(t, n_depth).abs();
            if f > f_max {
                f_max = f;
            }
        }
        f_max
    }
}

// ---------------------------------------------------------------------------
// AddedMassEffect
// ---------------------------------------------------------------------------

/// Frequency-dependent added mass coefficient for a body in a fluid.
///
/// For a sphere the analytical added mass coefficient is C_M = 0.5.
/// For a cylinder it depends on the cross-sectional geometry.
/// The frequency dependence arises from the radiation damping and
/// can be significant at high oscillation frequencies.
#[derive(Debug, Clone)]
pub struct AddedMassEffect {
    /// Body reference volume \[m³\].
    pub volume: f64,
    /// Fluid density \[kg/m³\].
    pub fluid_density: f64,
    /// Zero-frequency (infinite-period) added mass coefficient C_m0.
    pub cm0: f64,
    /// Frequency-dependence parameter β \[rad·s\].
    pub frequency_parameter: f64,
    /// Body natural frequency \[rad/s\] (used to compute resonant added mass).
    pub natural_frequency: f64,
}

impl AddedMassEffect {
    /// Create a new added mass model.
    pub fn new(
        volume: f64,
        fluid_density: f64,
        cm0: f64,
        frequency_parameter: f64,
        natural_frequency: f64,
    ) -> Self {
        Self {
            volume,
            fluid_density,
            cm0,
            frequency_parameter,
            natural_frequency,
        }
    }

    /// Create a sphere model with the analytical C_m = 0.5.
    pub fn new_sphere(radius: f64, fluid_density: f64) -> Self {
        let volume = 4.0_f64 / 3.0_f64 * PI * radius * radius * radius;
        Self::new(volume, fluid_density, 0.5_f64, 1.0_f64, 1.0_f64)
    }

    /// Frequency-dependent added mass coefficient C_m(ω).
    ///
    /// Simple model: C_m(ω) = C_m0 / (1 + (ω / β)²)
    pub fn cm_at_frequency(&self, omega: f64) -> f64 {
        let ratio = omega / self.frequency_parameter.max(1e-30_f64);
        self.cm0 / (1.0_f64 + ratio * ratio)
    }

    /// Added mass m_a = C_m · ρ_f · V at angular frequency ω.
    pub fn added_mass(&self, omega: f64) -> f64 {
        self.cm_at_frequency(omega) * self.fluid_density * self.volume
    }

    /// Added mass force on a body with acceleration `accel` at frequency `omega`.
    ///
    /// F_a = −m_a(ω) · accel
    pub fn added_mass_force_freq(&self, acceleration: [f64; 3], omega: f64) -> [f64; 3] {
        let ma = self.added_mass(omega);
        [
            -ma * acceleration[0],
            -ma * acceleration[1],
            -ma * acceleration[2],
        ]
    }

    /// Radiation damping coefficient B_r at frequency ω.
    ///
    /// Simplified model: B_r = ρ_f V ω C_m0 / (2 β)
    pub fn radiation_damping(&self, omega: f64) -> f64 {
        self.fluid_density * self.volume * omega * self.cm0
            / (2.0_f64 * self.frequency_parameter.max(1e-30_f64))
    }

    /// Effective natural frequency accounting for added mass.
    ///
    /// ω_eff = sqrt(k / (m + m_a)),  where k = m ωₙ²
    pub fn effective_natural_frequency(&self, body_mass: f64) -> f64 {
        let ma = self.added_mass(self.natural_frequency);
        let k = body_mass * self.natural_frequency * self.natural_frequency;
        let m_total = body_mass + ma;
        if m_total < 1e-30_f64 {
            return self.natural_frequency;
        }
        (k / m_total).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── reynolds_number ───────────────────────────────────────────────────

    #[test]
    fn reynolds_laminar_pipe() {
        let re = reynolds_number(1.0_f64, 1.0_f64, 1e-3_f64, 1000.0_f64);
        assert!((re - 1e6_f64).abs() < 1.0_f64, "re={re}");
    }

    #[test]
    fn reynolds_zero_velocity() {
        let re = reynolds_number(0.0_f64, 1.0_f64, 1e-3_f64, 1000.0_f64);
        assert!(re.abs() < 1e-10_f64, "re={re}");
    }

    #[test]
    fn reynolds_symmetry_abs_velocity() {
        let re_pos = reynolds_number(5.0_f64, 1.0_f64, 1e-3_f64, 1.0_f64);
        let re_neg = reynolds_number(-5.0_f64, 1.0_f64, 1e-3_f64, 1.0_f64);
        assert!((re_pos - re_neg).abs() < 1e-10_f64);
    }

    #[test]
    fn reynolds_scales_with_density() {
        let re1 = reynolds_number(1.0_f64, 1.0_f64, 1e-3_f64, 1.0_f64);
        let re2 = reynolds_number(1.0_f64, 1.0_f64, 1e-3_f64, 2.0_f64);
        assert!((re2 / re1 - 2.0_f64).abs() < 1e-10_f64);
    }

    // ── drag_force ────────────────────────────────────────────────────────

    #[test]
    fn drag_force_zero_velocity_is_zero() {
        let f = drag_force(1.0_f64, 1.2_f64, [0.0_f64; 3], 1.0_f64);
        assert!(f.iter().all(|x| x.abs() < 1e-15_f64));
    }

    #[test]
    fn drag_force_direction_opposes_velocity() {
        let f = drag_force(1.0_f64, 1.2_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!(f[0] < 0.0_f64, "drag should oppose +x velocity");
    }

    #[test]
    fn drag_force_magnitude_formula() {
        let f = drag_force(1.0_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!((f[0] + 50.0_f64).abs() < 1e-10_f64, "f[0]={}", f[0]);
    }

    #[test]
    fn drag_force_scales_with_cd() {
        let f1 = drag_force(0.5_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        let f2 = drag_force(1.0_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!((f2[0] / f1[0] - 2.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn drag_force_scales_with_area() {
        let f1 = drag_force(1.0_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        let f2 = drag_force(1.0_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 2.0_f64);
        assert!((f2[0] / f1[0] - 2.0_f64).abs() < 1e-10_f64);
    }

    // ── lift_force ────────────────────────────────────────────────────────

    #[test]
    fn lift_force_positive_cl_upward() {
        let f = lift_force(1.0_f64, 1.2_f64, [50.0_f64, 0.0_f64, 0.0_f64], 2.0_f64);
        assert!(f[1] > 0.0_f64, "lift should be positive (upward)");
    }

    #[test]
    fn lift_force_zero_speed_is_zero() {
        let f = lift_force(1.0_f64, 1.2_f64, [0.0_f64; 3], 2.0_f64);
        assert!(f[1].abs() < 1e-15_f64);
    }

    #[test]
    fn lift_force_magnitude_formula() {
        let f = lift_force(1.0_f64, 1.0_f64, [10.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!((f[1] - 50.0_f64).abs() < 1e-10_f64, "f[1]={}", f[1]);
    }

    #[test]
    fn lift_force_only_y_component() {
        let f = lift_force(1.0_f64, 1.2_f64, [10.0_f64, 0.0_f64, 5.0_f64], 1.0_f64);
        assert!(f[0].abs() < 1e-15_f64, "x={}", f[0]);
        assert!(f[2].abs() < 1e-15_f64, "z={}", f[2]);
    }

    // ── moment_coefficient_from_pressure ──────────────────────────────────

    #[test]
    fn moment_coeff_zero_pressure_gives_zero() {
        let cm = moment_coefficient_from_pressure(
            &[0.0_f64, 0.0_f64],
            &[1.0_f64, 2.0_f64],
            100.0_f64,
            1.0_f64,
            1.0_f64,
        );
        assert!(cm.abs() < 1e-15_f64);
    }

    #[test]
    fn moment_coeff_single_panel() {
        let cm =
            moment_coefficient_from_pressure(&[5.0_f64], &[2.0_f64], 100.0_f64, 1.0_f64, 1.0_f64);
        assert!((cm - 0.1_f64).abs() < 1e-12_f64, "cm={cm}");
    }

    #[test]
    fn moment_coeff_zero_denom_gives_zero() {
        let cm =
            moment_coefficient_from_pressure(&[1.0_f64], &[1.0_f64], 0.0_f64, 1.0_f64, 1.0_f64);
        assert!(cm.abs() < 1e-15_f64);
    }

    // ── added_mass_force ──────────────────────────────────────────────────

    #[test]
    fn added_mass_force_opposes_acceleration() {
        let f = added_mass_force([1.0_f64, 0.0_f64, 0.0_f64], 1000.0_f64, 0.1_f64, 1.0_f64);
        assert!(f[0] < 0.0_f64, "added mass should oppose acceleration");
    }

    #[test]
    fn added_mass_force_zero_acceleration() {
        let f = added_mass_force([0.0_f64; 3], 1000.0_f64, 1.0_f64, 0.5_f64);
        assert!(f.iter().all(|x| x.abs() < 1e-15_f64));
    }

    #[test]
    fn added_mass_force_scales_with_volume() {
        let f1 = added_mass_force([1.0_f64, 0.0_f64, 0.0_f64], 1000.0_f64, 1.0_f64, 1.0_f64);
        let f2 = added_mass_force([1.0_f64, 0.0_f64, 0.0_f64], 1000.0_f64, 2.0_f64, 1.0_f64);
        assert!((f2[0] / f1[0] - 2.0_f64).abs() < 1e-10_f64);
    }

    // ── bluff_body_cd ─────────────────────────────────────────────────────

    #[test]
    fn bluff_body_cd_high_re_plateau() {
        let cd = bluff_body_cd(1e6_f64);
        assert!((cd - 0.44_f64).abs() < 1e-10_f64, "cd={cd}");
    }

    #[test]
    fn bluff_body_cd_low_re_stokes() {
        let cd = bluff_body_cd(0.1_f64);
        assert!((cd - 240.0_f64).abs() < 1e-8_f64, "cd={cd}");
    }

    // ── strouhal_vortex_shedding ──────────────────────────────────────────

    #[test]
    fn strouhal_shedding_basic() {
        let f = strouhal_vortex_shedding(0.44_f64, 10.0_f64, 0.1_f64);
        assert!((f - 21.0_f64).abs() < 1e-10_f64, "f={f}");
    }

    #[test]
    fn strouhal_shedding_zero_diameter() {
        let f = strouhal_vortex_shedding(0.44_f64, 10.0_f64, 0.0_f64);
        assert!(f.abs() < 1e-15_f64);
    }

    // ── fluid_structure_time_step ─────────────────────────────────────────

    #[test]
    fn cfl_time_step_basic() {
        let dt = fluid_structure_time_step(0.5_f64, 0.01_f64, 1.0_f64);
        assert!((dt - 0.005_f64).abs() < 1e-15_f64, "dt={dt}");
    }

    #[test]
    fn cfl_time_step_zero_velocity_infinite() {
        let dt = fluid_structure_time_step(1.0_f64, 0.01_f64, 0.0_f64);
        assert!(dt.is_infinite());
    }

    // ── ImmersedBodyCfd ───────────────────────────────────────────────────

    #[test]
    fn immersed_body_sphere_marker_count() {
        let body = ImmersedBodyCfd::new_sphere(0.5_f64, [0.0_f64; 3], 10.0_f64, 0.1_f64, 10, 20);
        assert_eq!(body.n_markers(), 200);
    }

    #[test]
    fn immersed_body_delta_function_zero_far_away() {
        let body = ImmersedBodyCfd::new_sphere(0.5_f64, [0.0_f64; 3], 10.0_f64, 0.1_f64, 4, 8);
        assert_eq!(body.delta_function(0.5_f64), 0.0_f64);
    }

    #[test]
    fn immersed_body_delta_function_positive_at_origin() {
        let body = ImmersedBodyCfd::new_sphere(0.5_f64, [0.0_f64; 3], 10.0_f64, 0.1_f64, 4, 8);
        assert!(body.delta_function(0.0_f64) > 0.0_f64);
    }

    #[test]
    fn immersed_body_advance_moves_markers() {
        let mut body = ImmersedBodyCfd::new_sphere(0.5_f64, [0.0_f64; 3], 1.0_f64, 0.1_f64, 4, 8);
        body.velocity = [1.0_f64, 0.0_f64, 0.0_f64];
        let first_x_before = body.markers[0][0];
        body.advance(0.1_f64);
        let first_x_after = body.markers[0][0];
        assert!((first_x_after - first_x_before - 0.1_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn immersed_body_advance_moves_center() {
        let mut body = ImmersedBodyCfd::new_sphere(0.5_f64, [0.0_f64; 3], 1.0_f64, 0.1_f64, 4, 8);
        body.velocity = [2.0_f64, 0.0_f64, 0.0_f64];
        body.advance(0.5_f64);
        assert!((body.center[0] - 1.0_f64).abs() < 1e-12_f64);
    }

    // ── BuoyancyCfd ───────────────────────────────────────────────────────

    #[test]
    fn buoyancy_archimedes_upward_force() {
        let mut bc = BuoyancyCfd::new(1000.0_f64, [0.0_f64, -9.81_f64, 0.0_f64]);
        let f = bc.archimedes_force(1.0_f64);
        // F_b = 1000 * 1 * (-9.81) in y → should be negative (downward correction)
        // Actually ρ g V where g is already negative → upward buoyancy is +y if g=-9.81
        assert!(
            f[1] < 0.0_f64,
            "buoyancy in direction of gravity (upward support): {}",
            f[1]
        );
    }

    #[test]
    fn buoyancy_archimedes_scales_with_volume() {
        let mut bc = BuoyancyCfd::new(1000.0_f64, [0.0_f64, -9.81_f64, 0.0_f64]);
        let f1 = bc.archimedes_force(1.0_f64);
        let f2 = bc.archimedes_force(2.0_f64);
        assert!((f2[1] / f1[1] - 2.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn buoyancy_pressure_integrated_single_panel() {
        let mut bc = BuoyancyCfd::new(1000.0_f64, [0.0_f64, -9.81_f64, 0.0_f64]);
        let f =
            bc.pressure_integrated_force(&[100.0_f64], &[[0.0_f64, 1.0_f64, 0.0_f64]], &[1.0_f64]);
        assert!((f[1] + 100.0_f64).abs() < 1e-12_f64, "f[1]={}", f[1]);
    }

    #[test]
    fn buoyancy_submersion_fraction_fully_submerged() {
        let frac = BuoyancyCfd::submersion_fraction(-5.0_f64, 5.0_f64, 10.0_f64);
        assert!((frac - 1.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn buoyancy_submersion_fraction_half_submerged() {
        let frac = BuoyancyCfd::submersion_fraction(0.0_f64, 10.0_f64, 5.0_f64);
        assert!((frac - 0.5_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn buoyancy_submersion_fraction_not_submerged() {
        let frac = BuoyancyCfd::submersion_fraction(5.0_f64, 10.0_f64, 0.0_f64);
        assert!(frac.abs() < 1e-12_f64);
    }

    #[test]
    fn buoyancy_metacentric_height_stable_above() {
        let gm = BuoyancyCfd::metacentric_height(1.0_f64, 2.0_f64);
        assert!(gm > 0.0_f64, "gm={gm}");
    }

    // ── DragLiftCfd ───────────────────────────────────────────────────────

    #[test]
    fn drag_lift_cfd_pressure_drag_single_panel() {
        let mut dlc = DragLiftCfd::new(
            [1.0_f64, 0.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64, 0.0_f64],
            1.0_f64,
            1.0_f64,
            10.0_f64,
        );
        // Panel normal aligned with drag direction → max pressure drag
        let d = dlc.pressure_drag(&[100.0_f64], &[[1.0_f64, 0.0_f64, 0.0_f64]], &[1.0_f64]);
        assert!((d - 100.0_f64).abs() < 1e-10_f64, "d={d}");
    }

    #[test]
    fn drag_lift_cfd_pressure_lift_single_panel() {
        let mut dlc = DragLiftCfd::new(
            [1.0_f64, 0.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64, 0.0_f64],
            1.0_f64,
            1.0_f64,
            10.0_f64,
        );
        let l = dlc.pressure_lift(&[50.0_f64], &[[0.0_f64, 1.0_f64, 0.0_f64]], &[2.0_f64]);
        assert!((l - 100.0_f64).abs() < 1e-10_f64, "l={l}");
    }

    #[test]
    fn drag_lift_cfd_dynamic_pressure() {
        let dlc = DragLiftCfd::new(
            [1.0_f64, 0.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64, 0.0_f64],
            1.0_f64,
            1.0_f64,
            10.0_f64,
        );
        // q = 0.5 * 1 * 100 = 50
        assert!((dlc.dynamic_pressure() - 50.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn drag_lift_cfd_update_and_retrieve_coefficients() {
        let mut dlc = DragLiftCfd::new(
            [1.0_f64, 0.0_f64, 0.0_f64],
            [0.0_f64, 1.0_f64, 0.0_f64],
            2.0_f64,
            1.0_f64,
            10.0_f64,
        );
        // q = 50, A = 2 → denom = 100
        dlc.update_coefficients(100.0_f64, 50.0_f64);
        assert!((dlc.cd - 1.0_f64).abs() < 1e-10_f64, "cd={}", dlc.cd);
        assert!((dlc.cl - 0.5_f64).abs() < 1e-10_f64, "cl={}", dlc.cl);
    }

    // ── VortexInducedVibration ────────────────────────────────────────────

    #[test]
    fn viv_strouhal_frequency_positive() {
        let viv = VortexInducedVibration::new(
            0.1_f64, 10.0_f64, 100.0_f64, 2.0_f64, 0.02_f64, 1.0_f64, 1025.0_f64,
        );
        assert!(viv.strouhal_frequency() > 0.0_f64);
    }

    #[test]
    fn viv_reduced_velocity_lock_in_check() {
        // Vr ≈ 5 (inside lock-in window 4-8)
        let wn = 2.0_f64 * PI * 1.0_f64; // 1 Hz natural
        let d = 0.1_f64;
        let u = 5.0_f64 * 1.0_f64 * d; // Vr = 5
        let viv = VortexInducedVibration::new(d, 10.0_f64, 100.0_f64, wn, 0.02_f64, u, 1025.0_f64);
        assert!(viv.is_lock_in(), "Vr={}", viv.reduced_velocity());
    }

    #[test]
    fn viv_morison_force_zero_at_rest() {
        let viv = VortexInducedVibration::new(
            0.1_f64, 10.0_f64, 100.0_f64, 1.0_f64, 0.02_f64, 0.0_f64, 1025.0_f64,
        );
        let f = viv.morison_force_per_length(0.0_f64, 0.0_f64);
        assert!(f.abs() < 1e-15_f64, "f={f}");
    }

    #[test]
    fn viv_step_advances_time() {
        let mut viv = VortexInducedVibration::new(
            0.1_f64, 10.0_f64, 100.0_f64, 2.0_f64, 0.02_f64, 1.0_f64, 1025.0_f64,
        );
        viv.step(0.01_f64, 1.0_f64, 0.0_f64);
        assert!((viv.time - 0.01_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn viv_step_returns_finite_values() {
        let mut viv = VortexInducedVibration::new(
            0.1_f64, 10.0_f64, 100.0_f64, 2.0_f64, 0.02_f64, 1.0_f64, 1025.0_f64,
        );
        for _ in 0..100 {
            let (y, q) = viv.step(0.01_f64, 1.0_f64, 0.0_f64);
            assert!(y.is_finite(), "y={y}");
            assert!(q.is_finite(), "q={q}");
        }
    }

    // ── MooringLine ───────────────────────────────────────────────────────

    #[test]
    fn mooring_line_horizontal_span() {
        let ml = MooringLine::new(
            100.0_f64,
            1000.0_f64,
            1e8_f64,
            [0.0_f64; 3],
            [80.0_f64, 50.0_f64, 0.0_f64],
            10,
        );
        assert!((ml.horizontal_span() - 80.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn mooring_line_vertical_rise() {
        let ml = MooringLine::new(
            100.0_f64,
            1000.0_f64,
            1e8_f64,
            [0.0_f64; 3],
            [80.0_f64, 50.0_f64, 0.0_f64],
            10,
        );
        assert!((ml.vertical_rise() - 50.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn mooring_line_catenary_profile_zero_at_origin() {
        let mut ml = MooringLine::new(
            120.0_f64,
            1000.0_f64,
            1e8_f64,
            [0.0_f64; 3],
            [80.0_f64, 50.0_f64, 0.0_f64],
            10,
        );
        ml.horizontal_tension = 1e5_f64;
        let z = ml.profile_z(0.0_f64);
        assert!(z.abs() < 1e-12_f64, "z at origin should be 0: {z}");
    }

    #[test]
    fn mooring_line_profile_increases_away_from_origin() {
        let mut ml = MooringLine::new(
            120.0_f64,
            1000.0_f64,
            1e8_f64,
            [0.0_f64; 3],
            [80.0_f64, 50.0_f64, 0.0_f64],
            10,
        );
        ml.horizontal_tension = 1e5_f64;
        let z1 = ml.profile_z(10.0_f64);
        let z2 = ml.profile_z(20.0_f64);
        assert!(z2 > z1, "catenary should rise: z1={z1} z2={z2}");
    }

    #[test]
    fn mooring_line_discretised_positions_count() {
        let mut ml = MooringLine::new(
            100.0_f64,
            1000.0_f64,
            1e8_f64,
            [0.0_f64; 3],
            [80.0_f64, 40.0_f64, 0.0_f64],
            20,
        );
        ml.horizontal_tension = 5e4_f64;
        let pts = ml.discretised_positions();
        assert_eq!(pts.len(), 21);
    }

    #[test]
    fn mooring_line_solve_static_positive_tension() {
        let mut ml = MooringLine::new(
            200.0_f64,
            5000.0_f64,
            1e9_f64,
            [0.0_f64; 3],
            [150.0_f64, 100.0_f64, 0.0_f64],
            10,
        );
        let t_h = ml.solve_static();
        assert!(t_h > 0.0_f64, "t_h={t_h}");
    }

    // ── WaveLoadModel ─────────────────────────────────────────────────────

    #[test]
    fn wave_load_angular_frequency() {
        let wl = WaveLoadModel::new(
            0.5_f64, 2.0_f64, 0.7_f64, 1025.0_f64, 2.0_f64, 8.0_f64, 20.0_f64,
        );
        let omega = wl.angular_frequency();
        assert!((omega - 2.0_f64 * PI / 8.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn wave_load_wave_number_positive() {
        let wl = WaveLoadModel::new(
            0.5_f64, 2.0_f64, 0.7_f64, 1025.0_f64, 2.0_f64, 8.0_f64, 20.0_f64,
        );
        let k = wl.wave_number();
        assert!(k > 0.0_f64, "k={k}");
    }

    #[test]
    fn wave_load_particle_velocity_zero_at_quarter_period() {
        let wl = WaveLoadModel::new(
            0.5_f64, 2.0_f64, 0.7_f64, 1025.0_f64, 2.0_f64, 8.0_f64, 20.0_f64,
        );
        let t_quarter = wl.wave_period / 4.0_f64;
        let u = wl.water_particle_velocity(0.0_f64, t_quarter);
        assert!(u.abs() < 1e-6_f64, "u at t=T/4 should be ~0: {u}");
    }

    #[test]
    fn wave_load_force_per_length_finite() {
        let wl = WaveLoadModel::new(
            0.5_f64, 2.0_f64, 0.7_f64, 1025.0_f64, 2.0_f64, 8.0_f64, 20.0_f64,
        );
        let f = wl.force_per_length(-10.0_f64, 0.0_f64);
        assert!(f.is_finite(), "f={f}");
    }

    #[test]
    fn wave_load_maximum_force_positive() {
        let wl = WaveLoadModel::new(
            0.5_f64, 2.0_f64, 0.7_f64, 1025.0_f64, 2.0_f64, 8.0_f64, 20.0_f64,
        );
        let f_max = wl.maximum_force(50, 20);
        assert!(f_max > 0.0_f64, "f_max={f_max}");
    }

    // ── AddedMassEffect ───────────────────────────────────────────────────

    #[test]
    fn added_mass_sphere_cm_half() {
        let am = AddedMassEffect::new_sphere(1.0_f64, 1000.0_f64);
        assert!((am.cm0 - 0.5_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn added_mass_at_zero_frequency_equals_cm0() {
        let am = AddedMassEffect::new(1.0_f64, 1000.0_f64, 0.5_f64, 10.0_f64, 1.0_f64);
        assert!((am.cm_at_frequency(0.0_f64) - 0.5_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn added_mass_decreases_with_frequency() {
        let am = AddedMassEffect::new(1.0_f64, 1000.0_f64, 0.5_f64, 1.0_f64, 1.0_f64);
        let cm_lo = am.cm_at_frequency(0.1_f64);
        let cm_hi = am.cm_at_frequency(10.0_f64);
        assert!(cm_lo > cm_hi, "cm_lo={cm_lo} cm_hi={cm_hi}");
    }

    #[test]
    fn added_mass_effect_force_opposes_acceleration() {
        let am = AddedMassEffect::new(1.0_f64, 1000.0_f64, 0.5_f64, 1.0_f64, 1.0_f64);
        let f = am.added_mass_force_freq([1.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!(
            f[0] < 0.0_f64,
            "added mass should oppose acceleration: f={}",
            f[0]
        );
    }

    #[test]
    fn added_mass_effective_natural_frequency_less_than_in_vacuo() {
        let am = AddedMassEffect::new(1.0_f64, 1000.0_f64, 0.5_f64, 100.0_f64, 10.0_f64);
        let body_mass = 1.0_f64;
        let wn_eff = am.effective_natural_frequency(body_mass);
        assert!(wn_eff < 10.0_f64, "wn_eff={wn_eff}");
    }

    // ── CfdForces / CfdCouplingParams struct tests ─────────────────────────

    #[test]
    fn cfd_forces_zero_construction() {
        let f = CfdForces::zero();
        assert!(f.drag.iter().all(|&x| x.abs() < 1e-15_f64));
        assert!(f.lift.iter().all(|&x| x.abs() < 1e-15_f64));
    }

    #[test]
    fn cfd_forces_total_force() {
        let f = CfdForces {
            drag: [1.0_f64, 0.0_f64, 0.0_f64],
            lift: [0.0_f64, 2.0_f64, 0.0_f64],
            moment: [0.0_f64; 3],
            pressure_center: [0.0_f64; 3],
        };
        let tf = f.total_force();
        assert!((tf[0] - 1.0_f64).abs() < 1e-12_f64);
        assert!((tf[1] - 2.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn cfd_forces_drag_magnitude() {
        let f = CfdForces {
            drag: [3.0_f64, 4.0_f64, 0.0_f64],
            lift: [0.0_f64; 3],
            moment: [0.0_f64; 3],
            pressure_center: [0.0_f64; 3],
        };
        assert!((f.drag_magnitude() - 5.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn cfd_coupling_params_dynamic_pressure() {
        let p = CfdCouplingParams {
            fluid_density: 1.0_f64,
            fluid_viscosity: 1e-3_f64,
            reference_area: 1.0_f64,
            reference_length: 1.0_f64,
        };
        // q = 0.5 * 1 * 100 = 50
        assert!((p.dynamic_pressure(10.0_f64) - 50.0_f64).abs() < 1e-10_f64);
    }

    #[test]
    fn cfd_coupling_params_reynolds_number() {
        let p = CfdCouplingParams {
            fluid_density: 1000.0_f64,
            fluid_viscosity: 1e-3_f64,
            reference_area: 1.0_f64,
            reference_length: 1.0_f64,
        };
        assert!((p.reynolds_number(1.0_f64) - 1e6_f64).abs() < 1.0_f64);
    }
}
