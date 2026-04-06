// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Railway and tracked-vehicle dynamics.
//!
//! Provides wheel-rail contact mechanics (Hertz/Kalker), track geometry,
//! bogie hunting stability, and train resistance (Davis equation).
//!
//! # Overview
//!
//! - [`WheelRailContact`] — Hertz contact patch and Kalker creep forces
//! - [`TrackAlignment`] — tabulated track geometry with interpolation
//! - [`BogieDynamics`] — bogie hunting, yaw, and lateral equilibrium
//! - [`TrainResistance`] — Davis equation rolling resistance
//! - [`rail_bending_stress`] — simple beam bending stress in rail
//! - [`hertz_rail_contact_width`] — Hertz contact half-width
//! - [`grade_resistance`] — gravitational grade resistance

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Standard gravitational acceleration in m/s².
const G: f64 = 9.80665;

// ---------------------------------------------------------------------------
// WheelRailContact
// ---------------------------------------------------------------------------

/// Wheel-rail contact model based on Hertz contact theory and Kalker's linear
/// creep theory.
///
/// The contact patch is an ellipse with semi-axes a (longitudinal) and b
/// (lateral), computed from the equivalent radius and reduced elastic modulus.
pub struct WheelRailContact {
    /// Wheel rolling radius \[m\].
    pub wheel_radius: f64,
    /// Rail crown (transverse) radius \[m\].
    pub rail_crown_radius: f64,
    /// Rail inclination / contact angle \[rad\].
    pub contact_angle: f64,
    /// Kalker linear creep coefficient C₁₁ \[-\].
    pub creep_coeff: f64,
}

impl WheelRailContact {
    /// Create a standard wheel-rail contact for 1435 mm gauge railways.
    ///
    /// Typical values: wheel radius 0.46 m, rail crown 0.3 m, 1:20 rail tilt,
    /// Kalker C₁₁ = 200 000 N.
    pub fn standard() -> Self {
        Self {
            wheel_radius: 0.46,
            rail_crown_radius: 0.30,
            contact_angle: (1.0_f64 / 20.0_f64).atan(),
            creep_coeff: 200_000.0,
        }
    }

    /// Hertz contact patch semi-axes (a, b) for a given normal force.
    ///
    /// Uses the simplified Hertz formula for an elliptical contact between
    /// two cylinders.  The equivalent radius and reduced modulus E* are:
    ///
    /// ```text
    /// 1/R_eq = 1/R_w + 1/R_r
    /// (a, b) ≈ (3 F R_eq / (4 E*))^(1/3) · (scale factors)
    /// ```
    ///
    /// Returns `(a, b)` where a is the longitudinal and b the lateral
    /// semi-axis \[m\].
    ///
    /// # Parameters
    /// - `normal_force` – contact normal force N \[N\]
    /// - `e_star`       – reduced elastic modulus E* = E / (2(1−ν²)) \[Pa\]
    pub fn contact_patch_size(&self, normal_force: f64, e_star: f64) -> (f64, f64) {
        if normal_force <= 0.0 || e_star <= 0.0 {
            return (0.0, 0.0);
        }
        // Equivalent radius (series of two radii)
        let r_eq = 1.0 / (1.0 / self.wheel_radius + 1.0 / self.rail_crown_radius);
        // Hertz half-width for a cylinder-on-cylinder contact
        // a = sqrt(4 F R_eq / (π E* L)), but for wheel-rail we use sphere–flat
        // simplified: contact radius a = (3 F R_eq / (4 E*))^(1/3)
        let a = (3.0 * normal_force * r_eq / (4.0 * e_star)).cbrt();
        // Aspect ratio from rail crown vs wheel crown
        let aspect = (self.rail_crown_radius / self.wheel_radius)
            .sqrt()
            .clamp(0.5, 2.0);
        let b = a * aspect;
        (a, b)
    }

    /// Longitudinal (tractive) creep force using Kalker's linear theory.
    ///
    /// ```text
    /// F_x = −C₁₁ · G · a · b · ε_x
    /// ```
    ///
    /// where ε_x is the longitudinal creepage (slip / rolling velocity ratio).
    ///
    /// # Parameters
    /// - `creepage`     – longitudinal creepage ε_x \[-\]
    /// - `normal_force` – contact normal force \[N\]
    pub fn longitudinal_creep_force(&self, creepage: f64, normal_force: f64) -> f64 {
        // Kalker linear theory: F = C₁₁ * a*b * creepage (in units where G = shear mod)
        // Simplified normalisation: F_max ≈ μ * N, apply saturation via tanh
        let f_linear = self.creep_coeff * creepage;
        let f_max = 0.3 * normal_force; // nominal friction limit
        f_max * (f_linear / f_max.max(1e-6)).tanh()
    }

    /// Lateral creep force using Kalker's linear theory.
    ///
    /// ```text
    /// F_y = −C₂₂ · G · a · b · ε_y
    /// ```
    ///
    /// C₂₂ is typically slightly larger than C₁₁; we use C₂₂ = 1.1 · C₁₁.
    ///
    /// # Parameters
    /// - `creepage`     – lateral creepage ε_y \[-\]
    /// - `normal_force` – contact normal force \[N\]
    pub fn lateral_creep_force(&self, creepage: f64, normal_force: f64) -> f64 {
        let c22 = self.creep_coeff * 1.1;
        let f_linear = c22 * creepage;
        let f_max = 0.3 * normal_force;
        f_max * (f_linear / f_max.max(1e-6)).tanh()
    }

    /// Spin creep moment using Kalker's linear theory.
    ///
    /// ```text
    /// M_z = −C₂₃ · G · a · b² · φ
    /// ```
    ///
    /// where φ is the spin creepage \[1/m\].  C₂₃ ≈ 0.1 · C₁₁ (approximate).
    ///
    /// # Parameters
    /// - `spin`         – spin creepage φ \[1/m\]
    /// - `normal_force` – contact normal force \[N\]
    pub fn spin_creep_moment(&self, spin: f64, normal_force: f64) -> f64 {
        let c23 = self.creep_coeff * 0.1;
        let m_linear = c23 * spin;
        let m_max = 0.3 * normal_force * self.wheel_radius;
        m_max * (m_linear / m_max.max(1e-6)).tanh()
    }
}

// ---------------------------------------------------------------------------
// TrackAlignment
// ---------------------------------------------------------------------------

/// Tabulated track geometry along the route.
///
/// All arrays are indexed by chainage (running distance along the track) and
/// must all have the same length.
pub struct TrackAlignment {
    /// Chainage (arc length along track) stations \[m\].
    pub x: Vec<f64>,
    /// Track gauge at each station \[m\].
    pub gauge: Vec<f64>,
    /// Plan curvature 1/R at each station \[1/m\] (positive = left turn).
    pub curvature: Vec<f64>,
    /// Superelevation (cant) at each station \[m\].  Positive = left rail high.
    pub superelevation: Vec<f64>,
}

impl TrackAlignment {
    /// Construct a straight, level tangent track.
    ///
    /// # Parameters
    /// - `length`   – total track length \[m\]
    /// - `gauge`    – nominal track gauge \[m\] (1.435 m for standard gauge)
    /// - `n_points` – number of discretisation points
    pub fn from_tangent_track(length: f64, gauge: f64, n_points: usize) -> Self {
        let n = n_points.max(2);
        let x: Vec<f64> = (0..n).map(|i| i as f64 * length / (n - 1) as f64).collect();
        Self {
            gauge: vec![gauge; n],
            curvature: vec![0.0; n],
            superelevation: vec![0.0; n],
            x,
        }
    }

    /// Interpolate the track gauge at chainage `s` \[m\].
    pub fn gauge_at(&self, s: f64) -> f64 {
        self.interpolate(&self.gauge, s)
    }

    /// Interpolate curvature at chainage `s` \[m\].
    pub fn curvature_at(&self, s: f64) -> f64 {
        self.interpolate(&self.curvature, s)
    }

    /// Interpolate superelevation at chainage `s` \[m\].
    pub fn superelevation_at(&self, s: f64) -> f64 {
        self.interpolate(&self.superelevation, s)
    }

    /// Cant deficiency at chainage `s` for a vehicle travelling at `speed`
    /// with wheelbase `wheelbase`.
    ///
    /// Cant deficiency D = h_equil − h_actual, where the equilibrium cant for
    /// curve radius R at speed v is:
    ///
    /// ```text
    /// h_equil = v² · gauge / (g · R)
    /// ```
    ///
    /// Positive deficiency means the vehicle leans outward (more cant needed).
    ///
    /// # Parameters
    /// - `s`         – chainage \[m\]
    /// - `speed`     – vehicle speed \[m/s\]
    /// - `wheelbase` – rigid wheelbase \[m\] (used for cant gradient check — not
    ///   used in the scalar formula but kept for API completeness)
    pub fn cant_deficiency(&self, s: f64, speed: f64, _wheelbase: f64) -> f64 {
        let kappa = self.curvature_at(s);
        if kappa.abs() < 1e-12 {
            return 0.0; // straight track
        }
        let r = 1.0 / kappa.abs();
        let gauge = self.gauge_at(s);
        let h_equil = speed * speed * gauge / (G * r);
        let h_actual = self.superelevation_at(s);
        h_equil - h_actual
    }

    /// Linear interpolation helper for a data array indexed by `self.x`.
    fn interpolate(&self, data: &[f64], s: f64) -> f64 {
        let n = self.x.len();
        if n == 0 {
            return 0.0;
        }
        if s <= self.x[0] {
            return data[0];
        }
        if s >= self.x[n - 1] {
            return data[n - 1];
        }
        // Binary search for the bracket
        let idx = self
            .x
            .partition_point(|&xi| xi <= s)
            .saturating_sub(1)
            .min(n - 2);
        let x0 = self.x[idx];
        let x1 = self.x[idx + 1];
        let t = if (x1 - x0).abs() < 1e-14 {
            0.0
        } else {
            (s - x0) / (x1 - x0)
        };
        data[idx] * (1.0 - t) + data[idx + 1] * t
    }
}

// ---------------------------------------------------------------------------
// BogieDynamics
// ---------------------------------------------------------------------------

/// Two-axle bogie (truck) dynamics model.
///
/// Models lateral stability (hunting), steady-state yaw rate in curves,
/// and lateral force balance with superelevation.
pub struct BogieDynamics {
    /// Bogie mass \[kg\].
    pub mass: f64,
    /// Yaw moment of inertia about the vertical axis \[kg·m²\].
    pub inertia: f64,
    /// Axle spacing (rigid wheelbase) \[m\].
    pub wheelbase: f64,
    /// Lateral spring stiffness of primary suspension \[N/m\].
    pub spring_stiffness: f64,
    /// Lateral damping of primary suspension \[N·s/m\].
    pub damping: f64,
}

impl BogieDynamics {
    /// Critical (hunting instability) speed using the simplified Klingel formula.
    ///
    /// For a rigid two-axle bogie the critical speed is:
    ///
    /// ```text
    /// v_c = √(k_lat · a² / (m · λ))
    /// ```
    ///
    /// where λ is the conicity (contact angle) and a = wheelbase/2.
    ///
    /// This implementation uses the provided `lateral_stiffness` to refine
    /// the estimate via the effective suspension stiffness.
    ///
    /// Returns the critical speed in m/s (and f64::INFINITY for infinite
    /// stiffness).
    ///
    /// # Parameters
    /// - `lateral_stiffness` – effective lateral creep stiffness \[N/m\]
    pub fn critical_speed_hunting(&self, lateral_stiffness: f64) -> f64 {
        // Klingel formula: v_c = (L/2) * sqrt(k / (m * λ))
        // Using conicity λ = contact_angle ≈ 1/40 as a nominal value
        let conicity = 1.0_f64 / 40.0;
        let half_wheelbase = self.wheelbase * 0.5;
        let k_eff = self.spring_stiffness + lateral_stiffness;
        if k_eff <= 0.0 || self.mass <= 0.0 {
            return 0.0;
        }
        half_wheelbase * (k_eff / (self.mass * conicity)).sqrt()
    }

    /// Steady-state yaw rate of the bogie in a curve.
    ///
    /// For a vehicle traversing a curve of curvature κ at speed v:
    ///
    /// ```text
    /// ψ̇ = v · κ
    /// ```
    ///
    /// # Parameters
    /// - `speed`     – forward speed \[m/s\]
    /// - `curvature` – track curvature 1/R \[1/m\]
    pub fn yaw_rate(&self, speed: f64, curvature: f64) -> f64 {
        speed * curvature
    }

    /// Lateral force required to balance the bogie in a curve.
    ///
    /// Quasi-static lateral force from centripetal acceleration minus the
    /// compensating effect of superelevation:
    ///
    /// ```text
    /// F_lat = m · (v² · κ − g · sin(α))
    /// ```
    ///
    /// where α is the superelevation angle ≈ h_cant / gauge.
    ///
    /// # Parameters
    /// - `speed`           – forward speed \[m/s\]
    /// - `curvature`       – track curvature 1/R \[1/m\]
    /// - `superelevation`  – cant h \[m\] (assumes 1.435 m gauge)
    pub fn lateral_force_balance(&self, speed: f64, curvature: f64, superelevation: f64) -> f64 {
        let gauge = 1.435;
        let alpha = (superelevation / gauge).asin().min(PI / 6.0);
        let centripetal = speed * speed * curvature;
        self.mass * (centripetal - G * alpha.sin())
    }
}

// ---------------------------------------------------------------------------
// TrainResistance (Davis equation)
// ---------------------------------------------------------------------------

/// Train rolling resistance model based on the Davis equation.
///
/// The total running resistance is:
///
/// ```text
/// F = A + B·v + C·v²
/// ```
///
/// where A \[N\] is a constant term (journal bearing + wheel-rail deformation),
/// B \[N·s/m\] is a speed-proportional term, and C \[N·s²/m²\] is an
/// aerodynamic drag term.
pub struct TrainResistance {
    /// Train mass \[kg\].
    pub mass: f64,
    /// Constant resistance coefficient A \[N\].
    pub a: f64,
    /// Velocity-proportional resistance coefficient B \[N·s/m\].
    pub b: f64,
    /// Aerodynamic drag coefficient C \[N·s²/m²\].
    pub c: f64,
}

impl TrainResistance {
    /// Compute the total rolling resistance at speed `v` \[m/s\].
    ///
    /// ```text
    /// F = A + B·v + C·v²
    /// ```
    pub fn total_resistance(&self, speed: f64) -> f64 {
        self.a + self.b * speed + self.c * speed * speed
    }

    /// Resistance due to track gradient.
    ///
    /// ```text
    /// F_grade = m · g · (grade_per_mille / 1000)
    /// ```
    ///
    /// # Parameters
    /// - `grade_per_mille` – track gradient in ‰ (e.g. 10 means 1% or 1:100)
    pub fn gradient_resistance(&self, grade_per_mille: f64) -> f64 {
        self.mass * G * grade_per_mille / 1000.0
    }

    /// Construct a Davis model from specific resistance parameters.
    ///
    /// # Parameters
    /// - `mass` – train mass \[kg\]
    /// - `a`    – constant term \[N\]
    /// - `b`    – velocity term \[N·s/m\]
    /// - `c`    – quadratic term \[N·s²/m²\]
    pub fn new(mass: f64, a: f64, b: f64, c: f64) -> Self {
        Self { mass, a, b, c }
    }
}

// ---------------------------------------------------------------------------
// Stand-alone utility functions
// ---------------------------------------------------------------------------

/// Rail bending stress at the neutral axis under a point load.
///
/// For a simply supported rail of length L loaded at mid-span with force F:
///
/// ```text
/// σ = F · L / (4 · Z)
/// ```
///
/// where Z is the section modulus \[m³\].
///
/// # Parameters
/// - `force`            – wheel load (point force) \[N\]
/// - `section_modulus`  – rail section modulus Z \[m³\]
pub fn rail_bending_stress(force: f64, section_modulus: f64) -> f64 {
    // Simplified: σ = M / Z, with M ≈ F/2 (one-half of total, for single load)
    force / (2.0 * section_modulus)
}

/// Hertz contact half-width between a cylindrical wheel and a cylindrical rail.
///
/// For two parallel cylinders in contact:
///
/// ```text
/// b = √(4 · F · R_eq / (π · E* · L))
/// ```
///
/// Simplified to the sphere-on-flat 2-D case:
///
/// ```text
/// b = √(8 · F · R_eq / (π · E*))
/// ```
///
/// where R_eq = R_w · R_r / (R_w + R_r) and E* is the reduced modulus.
///
/// # Parameters
/// - `force`   – normal contact force \[N\]
/// - `r_wheel` – wheel rolling radius \[m\]
/// - `r_rail`  – rail crown radius \[m\]
/// - `e_star`  – reduced elastic modulus \[Pa\]
pub fn hertz_rail_contact_width(force: f64, r_wheel: f64, r_rail: f64, e_star: f64) -> f64 {
    if force <= 0.0 || e_star <= 0.0 {
        return 0.0;
    }
    let r_eq = r_wheel * r_rail / (r_wheel + r_rail);
    (8.0 * force * r_eq / (PI * e_star)).sqrt()
}

/// Grade resistance force on a vehicle climbing a slope.
///
/// ```text
/// F = m · g · sin(α)
/// ```
///
/// where α is the angle of the slope in radians.
///
/// For small grades, sin(α) ≈ grade (rise/run).
///
/// # Parameters
/// - `mass`  – vehicle mass \[kg\]
/// - `grade` – slope angle \[rad\]
pub fn grade_resistance(mass: f64, grade: f64) -> f64 {
    mass * G * grade.sin()
}

// ---------------------------------------------------------------------------
// WheelProfile
// ---------------------------------------------------------------------------

/// Wheel transverse profile parameters.
///
/// Defines the key geometric parameters of a railway wheel profile that
/// influence wheel-rail contact behaviour.
pub struct WheelProfile {
    /// Nominal rolling radius at the contact point \[m\].
    pub rolling_radius: f64,
    /// Flange angle measured from the tread tangent \[rad\].
    pub flange_angle: f64,
    /// Tread conicity λ (effective conicity, dimensionless).
    pub tread_conicity: f64,
    /// Flange height \[m\].
    pub flange_height: f64,
    /// Back-to-back distance between wheel flanges \[m\] (≈ gauge − 2×flange_thickness).
    pub back_to_back: f64,
}

impl WheelProfile {
    /// Create a standard UIC 60 wheel profile (approximate parameters).
    pub fn standard_uic60() -> Self {
        Self {
            rolling_radius: 0.46,
            flange_angle: (70.0_f64).to_radians(),
            tread_conicity: 0.025, // 1:40 conicity
            flange_height: 0.028,
            back_to_back: 1.360,
        }
    }

    /// Effective rolling radius at a lateral displacement `y` from nominal.
    ///
    /// Approximated by a linear taper: r(y) = r₀ − λ·y (right wheel).
    ///
    /// # Parameters
    /// - `y` – lateral displacement from nominal contact point \[m\]
    pub fn rolling_radius_at(&self, y: f64) -> f64 {
        (self.rolling_radius - self.tread_conicity * y).max(0.01)
    }

    /// Rolling radius difference between the two wheels of a wheelset at
    /// lateral displacement `y`.
    ///
    /// ΔR = 2λy (positive for rightward shift).
    pub fn rolling_radius_difference(&self, y: f64) -> f64 {
        2.0 * self.tread_conicity * y
    }
}

// ---------------------------------------------------------------------------
// RailProfile
// ---------------------------------------------------------------------------

/// Rail head profile parameters.
pub struct RailProfile {
    /// Rail head transverse crown radius (lateral radius of curvature) \[m\].
    pub head_radius: f64,
    /// Nominal track gauge (rail head to rail head inner face) \[m\].
    pub gauge: f64,
    /// Rail cant (inclination of the rail from vertical) \[rad\].
    pub cant: f64,
}

impl RailProfile {
    /// Standard UIC 60 rail on ballast (1:20 cant, 60E1 profile).
    pub fn standard_uic60() -> Self {
        Self {
            head_radius: 0.30,
            gauge: 1.435,
            cant: (1.0_f64 / 20.0_f64).atan(),
        }
    }

    /// Effective contact angle at a given lateral position on the rail head.
    ///
    /// For a cylindrical head, the contact angle varies as:
    ///
    /// ```text
    /// δ = cant + y / R_head
    /// ```
    ///
    /// # Parameters
    /// - `y` – lateral offset from rail head centre \[m\]
    pub fn contact_angle_at(&self, y: f64) -> f64 {
        self.cant + y / self.head_radius.max(1e-6)
    }
}

// ---------------------------------------------------------------------------
// HertzContact
// ---------------------------------------------------------------------------

/// Hertzian contact patch geometry for a wheel-rail contact pair.
///
/// Stores the semi-axes and maximum pressure of the Hertz elliptical contact.
pub struct HertzContact {
    /// Semi-axis in the rolling direction \[m\].
    pub a: f64,
    /// Semi-axis in the lateral direction \[m\].
    pub b: f64,
    /// Maximum contact pressure \[Pa\].
    pub p_max: f64,
}

impl HertzContact {
    /// Compute Hertz contact geometry for wheel-rail contact.
    ///
    /// Uses the simplified sphere-on-sphere model:
    ///
    /// ```text
    /// a = (3 F R_eq / (4 E*))^(1/3) × aspect
    /// b = a / aspect
    /// p_max = 3F / (2π a b)
    /// ```
    ///
    /// # Parameters
    /// - `normal_force`   – wheel load N \[N\]
    /// - `r_wheel`        – wheel rolling radius \[m\]
    /// - `r_rail_crown`   – rail crown radius \[m\]
    /// - `e_star`         – reduced elastic modulus E* \[Pa\]
    pub fn compute(normal_force: f64, r_wheel: f64, r_rail_crown: f64, e_star: f64) -> Self {
        if normal_force <= 0.0 || e_star <= 0.0 {
            return Self {
                a: 0.0,
                b: 0.0,
                p_max: 0.0,
            };
        }
        let r_eq = r_wheel * r_rail_crown / (r_wheel + r_rail_crown);
        let a_sphere = (3.0 * normal_force * r_eq / (4.0 * e_star)).cbrt();
        // Aspect ratio: wheel is longer in rolling direction
        let aspect = (r_rail_crown / r_wheel).sqrt().clamp(0.5, 2.0);
        let a = a_sphere * aspect;
        let b = a_sphere / aspect;
        let contact_area = PI * a * b;
        let p_max = if contact_area > 1e-20 {
            1.5 * normal_force / contact_area
        } else {
            0.0
        };
        Self { a, b, p_max }
    }

    /// Contact ellipse area \[m²\].
    pub fn area(&self) -> f64 {
        PI * self.a * self.b
    }
}

// ---------------------------------------------------------------------------
// KalkerContact
// ---------------------------------------------------------------------------

/// Kalker linear theory creep force model.
///
/// Computes tangential contact forces from creepages using Kalker's
/// linear creep coefficients C₁₁, C₂₂, C₂₃.
pub struct KalkerContact {
    /// Kalker coefficient C₁₁ for longitudinal creep \[-\].
    pub c11: f64,
    /// Kalker coefficient C₂₂ for lateral creep \[-\].
    pub c22: f64,
    /// Kalker coefficient C₂₃ for spin-lateral coupling \[-\].
    pub c23: f64,
    /// Shear modulus of contacting bodies G \[Pa\].
    pub shear_modulus: f64,
}

impl KalkerContact {
    /// Create a Kalker model with typical steel wheel-rail coefficients.
    ///
    /// For a contact ellipse with a/b ≈ 1.5, typical Kalker coefficients are
    /// C₁₁ ≈ 3.4, C₂₂ ≈ 5.0, C₂₃ ≈ 1.5 (non-dimensional).
    /// Shear modulus of steel G ≈ 80 GPa.
    pub fn standard_steel() -> Self {
        Self {
            c11: 3.4,
            c22: 5.0,
            c23: 1.5,
            shear_modulus: 80.0e9,
        }
    }

    /// Longitudinal creep force using Kalker's linear theory.
    ///
    /// ```text
    /// F_x = -G · a · b · C₁₁ · ε_x
    /// ```
    ///
    /// Saturated by tanh to respect the friction limit μ·N.
    ///
    /// # Parameters
    /// - `contact`    – Hertz contact geometry
    /// - `creepage_x` – longitudinal creepage ε_x \[-\]
    /// - `normal`     – contact normal force \[N\]
    /// - `mu`         – friction coefficient \[-\]
    pub fn longitudinal_force(
        &self,
        contact: &HertzContact,
        creepage_x: f64,
        normal: f64,
        mu: f64,
    ) -> f64 {
        let f_linear = self.shear_modulus * contact.a * contact.b * self.c11 * creepage_x;
        let f_max = mu * normal;
        if f_max < 1e-12 {
            return 0.0;
        }
        f_max * (f_linear / f_max).tanh()
    }

    /// Lateral creep force using Kalker's linear theory.
    ///
    /// ```text
    /// F_y = -G · a · b · (C₂₂ · ε_y + C₂₃ · √(a·b) · φ)
    /// ```
    ///
    /// # Parameters
    /// - `contact`    – Hertz contact geometry
    /// - `creepage_y` – lateral creepage ε_y \[-\]
    /// - `spin`       – spin creepage φ \[1/m\]
    /// - `normal`     – contact normal force \[N\]
    /// - `mu`         – friction coefficient \[-\]
    pub fn lateral_force(
        &self,
        contact: &HertzContact,
        creepage_y: f64,
        spin: f64,
        normal: f64,
        mu: f64,
    ) -> f64 {
        let ab = contact.a * contact.b;
        let f_linear =
            self.shear_modulus * ab * (self.c22 * creepage_y + self.c23 * ab.sqrt() * spin);
        let f_max = mu * normal;
        if f_max < 1e-12 {
            return 0.0;
        }
        f_max * (f_linear / f_max).tanh()
    }
}

// ---------------------------------------------------------------------------
// BogieFrame
// ---------------------------------------------------------------------------

/// Two-axle bogie (truck) frame parameters.
///
/// Defines the geometry and suspension properties of a railway bogie.
pub struct BogieFrame {
    /// Wheelset spacing (rigid wheelbase) \[m\].
    pub wheelset_spacing: f64,
    /// Primary suspension lateral stiffness \[N/m\].
    pub primary_stiffness: f64,
    /// Primary suspension lateral damping \[N·s/m\].
    pub primary_damping: f64,
    /// Secondary suspension lateral stiffness \[N/m\].
    pub secondary_stiffness: f64,
    /// Secondary suspension lateral damping \[N·s/m\].
    pub secondary_damping: f64,
    /// Bogie mass (frame only, not including wheelsets) \[kg\].
    pub frame_mass: f64,
    /// Bogie yaw moment of inertia \[kg·m²\].
    pub yaw_inertia: f64,
}

impl BogieFrame {
    /// Create a typical ICE-style bogie.
    pub fn standard_ice() -> Self {
        Self {
            wheelset_spacing: 2.5,
            primary_stiffness: 1_200_000.0,
            primary_damping: 30_000.0,
            secondary_stiffness: 250_000.0,
            secondary_damping: 80_000.0,
            frame_mass: 3_500.0,
            yaw_inertia: 2_500.0,
        }
    }

    /// Effective lateral stiffness seen by the car body (primary + secondary in series).
    ///
    /// ```text
    /// 1/k_eff = 1/k_primary + 1/k_secondary
    /// ```
    pub fn effective_lateral_stiffness(&self) -> f64 {
        let k1 = self.primary_stiffness;
        let k2 = self.secondary_stiffness;
        1.0 / (1.0 / k1 + 1.0 / k2)
    }
}

// ---------------------------------------------------------------------------
// TrackGeometry
// ---------------------------------------------------------------------------

/// Parameterised track geometry as a function of chainage.
///
/// All quantities vary linearly between the stored sample points.
pub struct TrackGeometry {
    /// Chainage stations \[m\].
    pub x: Vec<f64>,
    /// Curvature κ = 1/R at each station \[1/m\].
    pub curvature: Vec<f64>,
    /// Superelevation (cant) at each station \[m\].
    pub superelevation: Vec<f64>,
    /// Gauge widening at each station \[m\] (positive = wider).
    pub gauge_widening: Vec<f64>,
}

impl TrackGeometry {
    /// Create a circular curve with constant curvature and cant.
    ///
    /// # Parameters
    /// - `length`          – total arc length \[m\]
    /// - `radius`          – curve radius \[m\] (f64::INFINITY for straight)
    /// - `cant`            – superelevation \[m\]
    /// - `gauge_widening`  – extra gauge width \[m\]
    /// - `n_points`        – discretisation points
    pub fn circular_curve(
        length: f64,
        radius: f64,
        cant: f64,
        gauge_widening: f64,
        n_points: usize,
    ) -> Self {
        let n = n_points.max(2);
        let x: Vec<f64> = (0..n).map(|i| i as f64 * length / (n - 1) as f64).collect();
        let kappa = if radius.abs() > 1e-6 {
            1.0 / radius
        } else {
            0.0
        };
        Self {
            curvature: vec![kappa; n],
            superelevation: vec![cant; n],
            gauge_widening: vec![gauge_widening; n],
            x,
        }
    }

    /// Interpolate curvature at chainage `s`.
    pub fn curvature_at(&self, s: f64) -> f64 {
        Self::lerp(&self.x, &self.curvature, s)
    }

    /// Interpolate superelevation at chainage `s`.
    pub fn superelevation_at(&self, s: f64) -> f64 {
        Self::lerp(&self.x, &self.superelevation, s)
    }

    /// Interpolate gauge widening at chainage `s`.
    pub fn gauge_widening_at(&self, s: f64) -> f64 {
        Self::lerp(&self.x, &self.gauge_widening, s)
    }

    /// Linear interpolation helper.
    fn lerp(x: &[f64], y: &[f64], s: f64) -> f64 {
        let n = x.len();
        if n == 0 {
            return 0.0;
        }
        if s <= x[0] {
            return y[0];
        }
        if s >= x[n - 1] {
            return y[n - 1];
        }
        let idx = x
            .partition_point(|&xi| xi <= s)
            .saturating_sub(1)
            .min(n - 2);
        let t = (s - x[idx]) / (x[idx + 1] - x[idx]).max(1e-14);
        y[idx] * (1.0 - t) + y[idx + 1] * t
    }
}

// ---------------------------------------------------------------------------
// DerailmentCriteria
// ---------------------------------------------------------------------------

/// Wheel derailment criteria.
///
/// Implements the Nadal criterion for flange climb and provides detection of
/// wheel lift-off.
pub struct DerailmentCriteria {
    /// Flange angle \[rad\].
    pub flange_angle: f64,
    /// Wheel-rail friction coefficient.
    pub mu: f64,
}

impl DerailmentCriteria {
    /// Create a derailment criteria checker with default parameters.
    ///
    /// # Parameters
    /// - `flange_angle` – flange angle in radians
    /// - `mu`           – wheel-rail friction coefficient
    pub fn new(flange_angle: f64, mu: f64) -> Self {
        Self { flange_angle, mu }
    }

    /// Nadal limit value (L/V ratio at which flange climb begins).
    ///
    /// ```text
    /// (L/V)_Nadal = (tan δ − μ) / (1 + μ tan δ)
    /// ```
    ///
    /// where δ is the flange angle and μ is the friction coefficient.
    pub fn nadal_limit(&self) -> f64 {
        let tan_d = self.flange_angle.tan();
        (tan_d - self.mu) / (1.0 + self.mu * tan_d).max(1e-12)
    }

    /// Check if the wheel is at risk of flange climb given the lateral (L)
    /// and vertical (V) wheel-rail forces.
    ///
    /// Returns `true` if the L/V ratio exceeds the Nadal limit.
    ///
    /// # Parameters
    /// - `lateral_force`   – lateral guiding force L \[N\]
    /// - `vertical_force`  – vertical wheel load V \[N\]
    pub fn is_flange_climb(&self, lateral_force: f64, vertical_force: f64) -> bool {
        if vertical_force < 1.0 {
            return true; // wheel lift-off
        }
        let lv = lateral_force.abs() / vertical_force;
        lv >= self.nadal_limit()
    }

    /// Check for wheel lift-off (zero or negative vertical force).
    pub fn is_lift_off(&self, vertical_force: f64) -> bool {
        vertical_force < 1.0
    }
}

// ---------------------------------------------------------------------------
// CreepageCalc
// ---------------------------------------------------------------------------

/// Compute wheel-rail creepages from relative velocities.
///
/// Creepages are the normalised relative velocities at the contact patch.
pub struct CreepageCalc;

impl CreepageCalc {
    /// Longitudinal creepage ε_x (slip ratio).
    ///
    /// ```text
    /// ε_x = (v_wheel_rolling − V) / V
    /// ```
    ///
    /// where v_wheel_rolling = ω × r is the peripheral velocity of the wheel.
    ///
    /// # Parameters
    /// - `omega`        – wheel angular velocity \[rad/s\]
    /// - `wheel_radius` – rolling radius \[m\]
    /// - `vehicle_speed`– forward speed of the wheelset \[m/s\]
    pub fn longitudinal(omega: f64, wheel_radius: f64, vehicle_speed: f64) -> f64 {
        let v_rolling = omega * wheel_radius;
        let v_ref = vehicle_speed.abs().max(0.01);
        (v_rolling - vehicle_speed) / v_ref
    }

    /// Lateral creepage ε_y.
    ///
    /// ```text
    /// ε_y = V_lateral / V
    /// ```
    ///
    /// # Parameters
    /// - `lateral_velocity` – lateral velocity of the wheelset \[m/s\]
    /// - `vehicle_speed`    – forward speed \[m/s\]
    pub fn lateral(lateral_velocity: f64, vehicle_speed: f64) -> f64 {
        let v_ref = vehicle_speed.abs().max(0.01);
        lateral_velocity / v_ref
    }

    /// Spin creepage φ \[1/m\].
    ///
    /// ```text
    /// φ = (ψ̇ sin δ + Ω cos δ) / V
    /// ```
    ///
    /// Simplified to the yaw-dominated term:
    ///
    /// ```text
    /// φ ≈ ψ̇ / V
    /// ```
    ///
    /// # Parameters
    /// - `yaw_rate`      – wheelset yaw rate \[rad/s\]
    /// - `vehicle_speed` – forward speed \[m/s\]
    pub fn spin(yaw_rate: f64, vehicle_speed: f64) -> f64 {
        let v_ref = vehicle_speed.abs().max(0.01);
        yaw_rate / v_ref
    }
}

// ---------------------------------------------------------------------------
// TrackIrregularity
// ---------------------------------------------------------------------------

/// Random track irregularity generator based on a simplified PSD model.
///
/// Track irregularities (vertical, lateral, alignment, cross-level) are
/// modelled using power spectral density functions typical of Class 4/5
/// European railway track.
pub struct TrackIrregularity {
    /// Track quality class (1 = worst, 6 = best).
    pub quality_class: usize,
    /// Cutoff spatial frequency Ω_c \[rad/m\].
    pub omega_c: f64,
    /// Roughness level parameter A \[m²·rad/m\].
    pub roughness_a: f64,
}

impl TrackIrregularity {
    /// Create a new track irregularity model.
    ///
    /// # Parameters
    /// - `quality_class` – track quality class (1–6)
    pub fn new(quality_class: usize) -> Self {
        let qc = quality_class.clamp(1, 6);
        // Roughness parameter A scales with quality
        let roughness_a = match qc {
            1 => 4.032e-7,
            2 => 2.119e-7,
            3 => 1.080e-7,
            4 => 5.376e-8,
            5 => 2.688e-8,
            _ => 1.344e-8, // class 6
        };
        Self {
            quality_class: qc,
            omega_c: 0.8246,
            roughness_a,
        }
    }

    /// PSD of vertical track irregularity at spatial frequency Ω \[rad/m\].
    ///
    /// ```text
    /// S(Ω) = A · Ω_c² / (Ω² + Ω_c²)²
    /// ```
    pub fn psd_vertical(&self, omega: f64) -> f64 {
        let wc2 = self.omega_c * self.omega_c;
        let w2 = omega * omega;
        self.roughness_a * wc2 / (w2 + wc2).powi(2)
    }

    /// PSD of lateral track irregularity at spatial frequency Ω \[rad/m\].
    ///
    /// Lateral PSD is typically ≈ 0.8× vertical PSD.
    pub fn psd_lateral(&self, omega: f64) -> f64 {
        0.8 * self.psd_vertical(omega)
    }

    /// RMS amplitude of vertical irregularity over spatial frequency range.
    ///
    /// Integrates PSD from `omega_lo` to `omega_hi` with `n_pts` points.
    pub fn rms_vertical(&self, omega_lo: f64, omega_hi: f64, n_pts: usize) -> f64 {
        let n = n_pts.max(2);
        let dw = (omega_hi - omega_lo) / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let w = omega_lo + (i as f64 + 0.5) * dw;
            sum += self.psd_vertical(w) * dw;
        }
        sum.sqrt()
    }

    /// Standard deviation of cross-level irregularity \[rad\].
    ///
    /// Cross-level uses the same PSD but scaled by the reciprocal of gauge.
    ///
    /// # Parameters
    /// - `gauge` – track gauge \[m\]
    pub fn rms_cross_level(&self, gauge: f64, omega_lo: f64, omega_hi: f64, n_pts: usize) -> f64 {
        let rms_v = self.rms_vertical(omega_lo, omega_hi, n_pts);
        rms_v / gauge.max(0.1)
    }
}

// ---------------------------------------------------------------------------
// BogieDynamics::step  (extension of existing struct)
// ---------------------------------------------------------------------------

impl BogieDynamics {
    /// Integrate the bogie equations of motion by one time step `dt`.
    ///
    /// State vector: `[y, ẏ, ψ, ψ̇]` where y is lateral displacement \[m\]
    /// and ψ is yaw angle \[rad\].
    ///
    /// Forces applied:
    /// - Wheel-rail lateral creep forces (linearised Kalker)
    /// - Primary suspension restoring force
    /// - Track curvature input as a lateral disturbance
    ///
    /// Returns the updated state `[y, dy_dt, psi, dpsi_dt]`.
    ///
    /// # Parameters
    /// - `state`      – `[y, ẏ, ψ, ψ̇]`
    /// - `dt`         – time step \[s\]
    /// - `speed`      – forward speed \[m/s\]
    /// - `curvature`  – track curvature κ \[1/m\]
    /// - `lateral_k`  – lateral creep stiffness \[N/m\] (Kalker contribution)
    pub fn step(
        &self,
        state: [f64; 4],
        dt: f64,
        speed: f64,
        curvature: f64,
        lateral_k: f64,
    ) -> [f64; 4] {
        let [y, dy, psi, dpsi] = state;

        // Lateral equation of motion: m·ÿ = F_creep + F_suspension + F_curve
        let f_creep_y = -lateral_k * y; // restoring creep force
        let f_suspension = -self.spring_stiffness * y - self.damping * dy;
        let f_curve = self.mass * speed * speed * curvature; // centripetal input
        let ay = (f_creep_y + f_suspension + f_curve) / self.mass.max(1e-6);

        // Yaw equation: I·ψ̈ = M_yaw
        // Yaw moment from wheelset offset: M = -k_lateral * psi * (wheelbase/2)^2
        let half_wb = self.wheelbase * 0.5;
        let m_yaw = -lateral_k * psi * half_wb * half_wb
            - self.spring_stiffness * psi * half_wb
            - self.damping * dpsi * half_wb;
        let alpha_psi = m_yaw / self.inertia.max(1e-6);

        // Semi-implicit Euler integration
        let dy_new = dy + ay * dt;
        let y_new = y + dy_new * dt;
        let dpsi_new = dpsi + alpha_psi * dt;
        let psi_new = psi + dpsi_new * dt;

        [y_new, dy_new, psi_new, dpsi_new]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── WheelRailContact ─────────────────────────────────────────────────────

    #[test]
    fn wheel_rail_standard_positive_radii() {
        let c = WheelRailContact::standard();
        assert!(c.wheel_radius > 0.0);
        assert!(c.rail_crown_radius > 0.0);
    }

    #[test]
    fn contact_patch_size_positive_under_load() {
        let c = WheelRailContact::standard();
        let e_star = 1.0e11; // steel reduced modulus
        let (a, b) = c.contact_patch_size(100_000.0, e_star);
        assert!(a > 0.0, "a={a}");
        assert!(b > 0.0, "b={b}");
    }

    #[test]
    fn contact_patch_size_zero_for_zero_force() {
        let c = WheelRailContact::standard();
        let (a, b) = c.contact_patch_size(0.0, 1e11);
        assert_eq!(a, 0.0);
        assert_eq!(b, 0.0);
    }

    #[test]
    fn contact_patch_increases_with_force() {
        let c = WheelRailContact::standard();
        let e_star = 1e11;
        let (a1, _) = c.contact_patch_size(50_000.0, e_star);
        let (a2, _) = c.contact_patch_size(200_000.0, e_star);
        assert!(a2 > a1, "patch must grow with force");
    }

    #[test]
    fn longitudinal_creep_force_zero_at_zero_creepage() {
        let c = WheelRailContact::standard();
        let f = c.longitudinal_creep_force(0.0, 100_000.0);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    #[test]
    fn longitudinal_creep_force_sign_matches_creepage() {
        let c = WheelRailContact::standard();
        let f_pos = c.longitudinal_creep_force(0.01, 100_000.0);
        let f_neg = c.longitudinal_creep_force(-0.01, 100_000.0);
        assert!(f_pos > 0.0, "f_pos={f_pos}");
        assert!(f_neg < 0.0, "f_neg={f_neg}");
    }

    #[test]
    fn longitudinal_creep_force_saturates_at_friction_limit() {
        let c = WheelRailContact::standard();
        let n = 100_000.0;
        let f_small = c.longitudinal_creep_force(0.001, n);
        let f_large = c.longitudinal_creep_force(10.0, n);
        // Large creepage should saturate near μ·N = 30 000 N
        assert!(f_large < 0.35 * n, "should saturate: f_large={f_large}");
        assert!(
            f_large > f_small,
            "saturated force must exceed small creep force"
        );
    }

    #[test]
    fn lateral_creep_force_zero_at_zero_creepage() {
        let c = WheelRailContact::standard();
        let f = c.lateral_creep_force(0.0, 80_000.0);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    #[test]
    fn spin_creep_moment_zero_at_zero_spin() {
        let c = WheelRailContact::standard();
        let m = c.spin_creep_moment(0.0, 80_000.0);
        assert!(m.abs() < 1e-10, "m={m}");
    }

    #[test]
    fn spin_creep_moment_sign() {
        let c = WheelRailContact::standard();
        let m = c.spin_creep_moment(1.0, 80_000.0);
        assert!(m > 0.0, "m={m}");
    }

    // ── TrackAlignment ───────────────────────────────────────────────────────

    #[test]
    fn tangent_track_has_zero_curvature() {
        let t = TrackAlignment::from_tangent_track(1000.0, 1.435, 100);
        for k in &t.curvature {
            assert!(k.abs() < 1e-12, "curvature={k} on tangent track");
        }
    }

    #[test]
    fn tangent_track_gauge_is_nominal() {
        let t = TrackAlignment::from_tangent_track(500.0, 1.435, 50);
        assert!((t.gauge_at(250.0) - 1.435).abs() < 1e-10);
    }

    #[test]
    fn track_gauge_at_interpolates_correctly() {
        let mut t = TrackAlignment::from_tangent_track(100.0, 1.435, 11);
        // Widen the gauge linearly from 1.435 to 1.535
        for (i, g) in t.gauge.iter_mut().enumerate() {
            *g = 1.435 + i as f64 * 0.01;
        }
        let g_mid = t.gauge_at(50.0);
        assert!((g_mid - 1.485).abs() < 0.005, "g_mid={g_mid}");
    }

    #[test]
    fn cant_deficiency_zero_on_tangent() {
        let t = TrackAlignment::from_tangent_track(1000.0, 1.435, 100);
        let d = t.cant_deficiency(500.0, 30.0, 2.5);
        assert!(
            d.abs() < 1e-10,
            "cant deficiency on tangent must be zero: d={d}"
        );
    }

    #[test]
    fn cant_deficiency_positive_for_high_speed_curve() {
        let mut t = TrackAlignment::from_tangent_track(1000.0, 1.435, 100);
        // Set curvature to 1/500 (500 m radius)
        for k in t.curvature.iter_mut() {
            *k = 1.0 / 500.0;
        }
        // No superelevation
        let d = t.cant_deficiency(500.0, 40.0, 2.5); // 40 m/s ≈ 144 km/h
        assert!(d > 0.0, "cant deficiency should be positive: d={d}");
    }

    #[test]
    fn track_interpolation_clamps_at_boundaries() {
        let t = TrackAlignment::from_tangent_track(100.0, 1.435, 10);
        assert!((t.gauge_at(-10.0) - 1.435).abs() < 1e-10);
        assert!((t.gauge_at(200.0) - 1.435).abs() < 1e-10);
    }

    #[test]
    fn track_alignment_n_points_minimum_two() {
        let t = TrackAlignment::from_tangent_track(100.0, 1.435, 1);
        // n is clamped to max(n, 2)
        assert!(t.x.len() >= 2);
    }

    // ── BogieDynamics ────────────────────────────────────────────────────────

    #[test]
    fn bogie_critical_speed_positive() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1_000_000.0,
            damping: 20_000.0,
        };
        let v_c = bogie.critical_speed_hunting(500_000.0);
        assert!(v_c > 0.0, "v_c={v_c}");
    }

    #[test]
    fn bogie_critical_speed_increases_with_stiffness() {
        let make_bogie = |k: f64| BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: k,
            damping: 20_000.0,
        };
        let v_lo = make_bogie(100_000.0).critical_speed_hunting(0.0);
        let v_hi = make_bogie(4_000_000.0).critical_speed_hunting(0.0);
        assert!(v_hi > v_lo, "v_hi={v_hi} v_lo={v_lo}");
    }

    #[test]
    fn bogie_yaw_rate_proportional_to_speed() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        let kappa = 1.0 / 500.0;
        let psi1 = bogie.yaw_rate(20.0, kappa);
        let psi2 = bogie.yaw_rate(40.0, kappa);
        assert!(
            (psi2 - 2.0 * psi1).abs() < 1e-10,
            "yaw rate not proportional"
        );
    }

    #[test]
    fn bogie_yaw_rate_zero_on_tangent() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        assert!(bogie.yaw_rate(30.0, 0.0).abs() < 1e-12);
    }

    #[test]
    fn lateral_force_balance_zero_at_equilibrium_cant() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        // Equilibrium: h = v² * gauge / (g * R)
        let speed = 30.0;
        let r = 500.0;
        let kappa = 1.0 / r;
        let gauge = 1.435;
        let h_equil = speed * speed * gauge / (G * r);
        let f = bogie.lateral_force_balance(speed, kappa, h_equil);
        assert!(
            f.abs() < 10.0,
            "lateral force at equilibrium cant should be near zero: f={f}"
        );
    }

    // ── TrainResistance ──────────────────────────────────────────────────────

    #[test]
    fn train_resistance_at_zero_speed_equals_a() {
        let r = TrainResistance::new(100_000.0, 5000.0, 10.0, 0.5);
        assert!((r.total_resistance(0.0) - 5000.0).abs() < 1e-8);
    }

    #[test]
    fn train_resistance_increases_with_speed() {
        let r = TrainResistance::new(100_000.0, 5000.0, 10.0, 0.5);
        let f1 = r.total_resistance(10.0);
        let f2 = r.total_resistance(30.0);
        assert!(
            f2 > f1,
            "resistance should increase with speed: f1={f1} f2={f2}"
        );
    }

    #[test]
    fn train_resistance_quadratic_term_dominates_at_high_speed() {
        // a=0, b=0, c=5 → pure quadratic
        let r = TrainResistance::new(100_000.0, 0.0, 0.0, 5.0);
        let f_slow = r.total_resistance(1.0); // 5 N
        let f_fast = r.total_resistance(100.0); // 50 000 N
        // f_fast should be 10 000 × f_slow
        assert!(
            f_fast > 1000.0 * f_slow,
            "aerodynamic drag should dominate at high speed: slow={f_slow} fast={f_fast}"
        );
    }

    #[test]
    fn gradient_resistance_positive_on_uphill() {
        let r = TrainResistance::new(100_000.0, 0.0, 0.0, 0.0);
        let f = r.gradient_resistance(10.0); // 10‰ = 1%
        assert!(f > 0.0, "uphill grade resistance must be positive: f={f}");
    }

    #[test]
    fn gradient_resistance_proportional_to_mass() {
        let r1 = TrainResistance::new(50_000.0, 0.0, 0.0, 0.0);
        let r2 = TrainResistance::new(100_000.0, 0.0, 0.0, 0.0);
        let f1 = r1.gradient_resistance(10.0);
        let f2 = r2.gradient_resistance(10.0);
        assert!((f2 - 2.0 * f1).abs() < 1e-6, "f2={f2} 2*f1={}", 2.0 * f1);
    }

    #[test]
    fn gradient_resistance_zero_on_level() {
        let r = TrainResistance::new(100_000.0, 0.0, 0.0, 0.0);
        assert!(r.gradient_resistance(0.0).abs() < 1e-12);
    }

    // ── Utility functions ────────────────────────────────────────────────────

    #[test]
    fn rail_bending_stress_positive_under_load() {
        let sigma = rail_bending_stress(100_000.0, 5e-4);
        assert!(sigma > 0.0, "sigma={sigma}");
    }

    #[test]
    fn rail_bending_stress_proportional_to_force() {
        let s1 = rail_bending_stress(100_000.0, 5e-4);
        let s2 = rail_bending_stress(200_000.0, 5e-4);
        assert!((s2 - 2.0 * s1).abs() < 1e-6, "s2={s2} 2*s1={}", 2.0 * s1);
    }

    #[test]
    fn hertz_contact_width_positive_under_load() {
        let b = hertz_rail_contact_width(100_000.0, 0.46, 0.3, 1e11);
        assert!(b > 0.0, "b={b}");
    }

    #[test]
    fn hertz_contact_width_zero_for_zero_force() {
        let b = hertz_rail_contact_width(0.0, 0.46, 0.3, 1e11);
        assert_eq!(b, 0.0);
    }

    #[test]
    fn hertz_contact_width_increases_with_force() {
        let b1 = hertz_rail_contact_width(50_000.0, 0.46, 0.3, 1e11);
        let b2 = hertz_rail_contact_width(200_000.0, 0.46, 0.3, 1e11);
        assert!(b2 > b1, "contact width must grow with force");
    }

    #[test]
    fn grade_resistance_zero_on_level() {
        assert!(grade_resistance(100_000.0, 0.0).abs() < 1e-12);
    }

    #[test]
    fn grade_resistance_positive_uphill() {
        let f = grade_resistance(100_000.0, 0.01); // ~0.57° slope
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn grade_resistance_proportional_to_mass() {
        let f1 = grade_resistance(50_000.0, 0.05);
        let f2 = grade_resistance(100_000.0, 0.05);
        assert!((f2 - 2.0 * f1).abs() < 1e-6);
    }

    #[test]
    fn grade_resistance_known_value() {
        // m=1000 kg, grade=PI/2 (vertical) → F = 1000 * 9.80665 * 1.0
        let f = grade_resistance(1000.0, PI / 2.0);
        assert!((f - 1000.0 * G).abs() < 1e-6, "f={f}");
    }

    #[test]
    fn bogie_critical_speed_zero_mass_returns_zero() {
        let bogie = BogieDynamics {
            mass: 0.0,
            inertia: 1.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        let v = bogie.critical_speed_hunting(0.0);
        assert!(!v.is_nan(), "should not be NaN");
    }

    #[test]
    fn contact_patch_zero_for_zero_modulus() {
        let c = WheelRailContact::standard();
        let (a, b) = c.contact_patch_size(100_000.0, 0.0);
        assert_eq!(a, 0.0);
        assert_eq!(b, 0.0);
    }

    #[test]
    fn train_resistance_new_stores_params() {
        let r = TrainResistance::new(500_000.0, 3000.0, 20.0, 2.0);
        assert_eq!(r.mass, 500_000.0);
        assert_eq!(r.a, 3000.0);
        assert_eq!(r.b, 20.0);
        assert_eq!(r.c, 2.0);
    }

    #[test]
    fn lateral_creep_greater_than_longitudinal_for_same_creepage() {
        // C22 = 1.1 * C11, so lateral force should be slightly larger
        let c = WheelRailContact::standard();
        let f_lon = c.longitudinal_creep_force(0.001, 80_000.0);
        let f_lat = c.lateral_creep_force(0.001, 80_000.0);
        assert!(
            f_lat >= f_lon,
            "lateral force should be >= longitudinal: lat={f_lat} lon={f_lon}"
        );
    }

    #[test]
    fn hertz_contact_width_formula_check() {
        // F=1 N, R_w=R_r=1 m → R_eq=0.5, b = sqrt(8*1*0.5/(pi*1e11)) ≈ very small
        let b = hertz_rail_contact_width(1.0, 1.0, 1.0, 1e11);
        let expected = (8.0 * 1.0 * 0.5 / (PI * 1e11_f64)).sqrt();
        assert!((b - expected).abs() < 1e-15, "b={b} expected={expected}");
    }

    // ── WheelProfile ─────────────────────────────────────────────────────────

    #[test]
    fn wheel_profile_standard_positive_radius() {
        let w = WheelProfile::standard_uic60();
        assert!(w.rolling_radius > 0.0);
        assert!(w.tread_conicity > 0.0);
    }

    #[test]
    fn wheel_profile_rolling_radius_at_nominal_zero_offset() {
        let w = WheelProfile::standard_uic60();
        assert!((w.rolling_radius_at(0.0) - w.rolling_radius).abs() < 1e-10);
    }

    #[test]
    fn wheel_profile_rolling_radius_decreases_rightward() {
        let w = WheelProfile::standard_uic60();
        let r_left = w.rolling_radius_at(-0.01);
        let r_right = w.rolling_radius_at(0.01);
        assert!(
            r_right < r_left,
            "rightward displacement reduces right wheel radius"
        );
    }

    #[test]
    fn wheel_profile_radius_difference_proportional_to_offset() {
        let w = WheelProfile::standard_uic60();
        let dr1 = w.rolling_radius_difference(0.005);
        let dr2 = w.rolling_radius_difference(0.010);
        assert!(
            (dr2 - 2.0 * dr1).abs() < 1e-12,
            "dr2={dr2} 2*dr1={}",
            2.0 * dr1
        );
    }

    // ── RailProfile ───────────────────────────────────────────────────────────

    #[test]
    fn rail_profile_standard_gauge() {
        let r = RailProfile::standard_uic60();
        assert!((r.gauge - 1.435).abs() < 1e-10);
    }

    #[test]
    fn rail_profile_contact_angle_at_centre_equals_cant() {
        let r = RailProfile::standard_uic60();
        let angle = r.contact_angle_at(0.0);
        assert!((angle - r.cant).abs() < 1e-10, "angle at centre={angle}");
    }

    #[test]
    fn rail_profile_contact_angle_increases_with_offset() {
        let r = RailProfile::standard_uic60();
        let a1 = r.contact_angle_at(0.0);
        let a2 = r.contact_angle_at(0.01);
        assert!(a2 > a1, "angle should increase with lateral offset");
    }

    // ── HertzContact ──────────────────────────────────────────────────────────

    #[test]
    fn hertz_contact_positive_under_load() {
        let hc = HertzContact::compute(100_000.0, 0.46, 0.3, 1e11);
        assert!(hc.a > 0.0, "a={}", hc.a);
        assert!(hc.b > 0.0, "b={}", hc.b);
        assert!(hc.p_max > 0.0, "p_max={}", hc.p_max);
    }

    #[test]
    fn hertz_contact_zero_for_zero_force() {
        let hc = HertzContact::compute(0.0, 0.46, 0.3, 1e11);
        assert_eq!(hc.a, 0.0);
        assert_eq!(hc.b, 0.0);
    }

    #[test]
    fn hertz_contact_area_positive() {
        let hc = HertzContact::compute(100_000.0, 0.46, 0.3, 1e11);
        assert!(hc.area() > 0.0);
    }

    #[test]
    fn hertz_contact_a_greater_than_b_for_wheel_on_flat_rail() {
        // When r_rail_crown > r_wheel, aspect = sqrt(r_rail/r_wheel) > 1, so a > b
        // Use r_rail_crown (0.6) > r_wheel (0.46)
        let hc = HertzContact::compute(100_000.0, 0.46, 0.6, 1e11);
        assert!(hc.a > hc.b, "a={} b={}", hc.a, hc.b);
    }

    #[test]
    fn hertz_contact_grows_with_force() {
        let hc1 = HertzContact::compute(50_000.0, 0.46, 0.3, 1e11);
        let hc2 = HertzContact::compute(200_000.0, 0.46, 0.3, 1e11);
        assert!(hc2.a > hc1.a, "contact size must grow with force");
    }

    // ── KalkerContact ─────────────────────────────────────────────────────────

    #[test]
    fn kalker_longitudinal_force_zero_at_zero_creepage() {
        let k = KalkerContact::standard_steel();
        let hc = HertzContact::compute(100_000.0, 0.46, 0.3, 1e11);
        let f = k.longitudinal_force(&hc, 0.0, 100_000.0, 0.3);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    #[test]
    fn kalker_longitudinal_force_sign_follows_creepage() {
        let k = KalkerContact::standard_steel();
        let hc = HertzContact::compute(100_000.0, 0.46, 0.3, 1e11);
        let f_pos = k.longitudinal_force(&hc, 0.01, 100_000.0, 0.3);
        let f_neg = k.longitudinal_force(&hc, -0.01, 100_000.0, 0.3);
        assert!(f_pos > 0.0, "f_pos={f_pos}");
        assert!(f_neg < 0.0, "f_neg={f_neg}");
    }

    #[test]
    fn kalker_lateral_force_zero_at_zero_inputs() {
        let k = KalkerContact::standard_steel();
        let hc = HertzContact::compute(100_000.0, 0.46, 0.3, 1e11);
        let f = k.lateral_force(&hc, 0.0, 0.0, 100_000.0, 0.3);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    // ── BogieFrame ────────────────────────────────────────────────────────────

    #[test]
    fn bogie_frame_standard_positive_properties() {
        let bf = BogieFrame::standard_ice();
        assert!(bf.frame_mass > 0.0);
        assert!(bf.wheelset_spacing > 0.0);
        assert!(bf.yaw_inertia > 0.0);
    }

    #[test]
    fn bogie_frame_effective_stiffness_less_than_minimum() {
        let bf = BogieFrame::standard_ice();
        let k_eff = bf.effective_lateral_stiffness();
        let k_min = bf.primary_stiffness.min(bf.secondary_stiffness);
        assert!(
            k_eff < k_min,
            "series stiffness must be less than minimum: k_eff={k_eff} k_min={k_min}"
        );
    }

    // ── TrackGeometry ─────────────────────────────────────────────────────────

    #[test]
    fn track_geometry_circular_curve_constant_curvature() {
        let tg = TrackGeometry::circular_curve(1000.0, 500.0, 0.1, 0.0, 50);
        let kappa = tg.curvature_at(500.0);
        assert!((kappa - 1.0 / 500.0).abs() < 1e-12, "kappa={kappa}");
    }

    #[test]
    fn track_geometry_straight_track_zero_curvature() {
        let tg = TrackGeometry::circular_curve(1000.0, f64::INFINITY, 0.0, 0.0, 10);
        assert!(tg.curvature_at(500.0).abs() < 1e-12);
    }

    #[test]
    fn track_geometry_superelevation_constant() {
        let tg = TrackGeometry::circular_curve(500.0, 1000.0, 0.15, 0.003, 20);
        let cant = tg.superelevation_at(250.0);
        assert!((cant - 0.15).abs() < 1e-10, "cant={cant}");
    }

    #[test]
    fn track_geometry_gauge_widening_correct() {
        let tg = TrackGeometry::circular_curve(500.0, 300.0, 0.1, 0.005, 20);
        let gw = tg.gauge_widening_at(250.0);
        assert!((gw - 0.005).abs() < 1e-10, "gw={gw}");
    }

    // ── DerailmentCriteria ────────────────────────────────────────────────────

    #[test]
    fn derailment_nadal_limit_positive() {
        let d = DerailmentCriteria::new((70.0_f64).to_radians(), 0.3);
        assert!(d.nadal_limit() > 0.0);
    }

    #[test]
    fn derailment_nadal_limit_increases_with_flange_angle() {
        let d1 = DerailmentCriteria::new((60.0_f64).to_radians(), 0.3);
        let d2 = DerailmentCriteria::new((75.0_f64).to_radians(), 0.3);
        assert!(
            d2.nadal_limit() > d1.nadal_limit(),
            "steeper flange → higher Nadal limit"
        );
    }

    #[test]
    fn derailment_not_flange_climb_under_normal_forces() {
        let d = DerailmentCriteria::new((70.0_f64).to_radians(), 0.3);
        // L/V = 0.2 should be safe (Nadal limit ≈ 0.9 for 70° flange)
        assert!(
            !d.is_flange_climb(20_000.0, 100_000.0),
            "safe L/V should not trigger flange climb"
        );
    }

    #[test]
    fn derailment_flags_lift_off_at_zero_load() {
        let d = DerailmentCriteria::new((70.0_f64).to_radians(), 0.3);
        assert!(d.is_lift_off(0.0));
        assert!(!d.is_lift_off(50_000.0));
    }

    // ── CreepageCalc ──────────────────────────────────────────────────────────

    #[test]
    fn creepage_longitudinal_zero_free_rolling() {
        // Free rolling: omega * r = v
        let v = 30.0;
        let r = 0.46;
        let omega = v / r;
        let eps = CreepageCalc::longitudinal(omega, r, v);
        assert!(eps.abs() < 1e-10, "free rolling creepage={eps}");
    }

    #[test]
    fn creepage_longitudinal_positive_for_driven_wheel() {
        // Driven: omega * r > v
        let eps = CreepageCalc::longitudinal(70.0, 0.46, 30.0);
        assert!(eps > 0.0, "driven wheel creepage={eps}");
    }

    #[test]
    fn creepage_lateral_zero_when_no_lateral_velocity() {
        let eps = CreepageCalc::lateral(0.0, 30.0);
        assert!(eps.abs() < 1e-10, "zero lateral creepage={eps}");
    }

    #[test]
    fn creepage_spin_zero_for_straight_track() {
        let phi = CreepageCalc::spin(0.0, 30.0);
        assert!(phi.abs() < 1e-10, "spin={phi}");
    }

    #[test]
    fn creepage_spin_positive_for_positive_yaw_rate() {
        let phi = CreepageCalc::spin(0.5, 30.0);
        assert!(phi > 0.0, "phi={phi}");
    }

    // ── TrackIrregularity ─────────────────────────────────────────────────────

    #[test]
    fn track_irregularity_psd_positive() {
        let ti = TrackIrregularity::new(4);
        let psd = ti.psd_vertical(1.0);
        assert!(psd > 0.0, "PSD must be positive: psd={psd}");
    }

    #[test]
    fn track_irregularity_psd_decreases_at_high_freq() {
        let ti = TrackIrregularity::new(4);
        let psd_lo = ti.psd_vertical(0.1);
        let psd_hi = ti.psd_vertical(10.0);
        assert!(psd_hi < psd_lo, "PSD must decrease at high frequency");
    }

    #[test]
    fn track_irregularity_worse_class_has_higher_psd() {
        let ti_bad = TrackIrregularity::new(1);
        let ti_good = TrackIrregularity::new(5);
        assert!(
            ti_bad.psd_vertical(1.0) > ti_good.psd_vertical(1.0),
            "worse class → higher PSD"
        );
    }

    #[test]
    fn track_irregularity_rms_positive() {
        let ti = TrackIrregularity::new(3);
        let rms = ti.rms_vertical(0.1, 10.0, 100);
        assert!(rms > 0.0, "rms={rms}");
    }

    #[test]
    fn track_irregularity_lateral_less_than_vertical() {
        let ti = TrackIrregularity::new(4);
        let psd_v = ti.psd_vertical(1.0);
        let psd_l = ti.psd_lateral(1.0);
        assert!(psd_l < psd_v, "lateral PSD should be less than vertical");
    }

    // ── BogieDynamics::step ───────────────────────────────────────────────────

    #[test]
    fn bogie_step_zero_dt_returns_same_state() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        let state = [0.01, 0.0, 0.001, 0.0];
        let new_state = bogie.step(state, 0.0, 30.0, 0.002, 5e5);
        for (a, b) in state.iter().zip(new_state.iter()) {
            assert!((a - b).abs() < 1e-12, "state should not change with dt=0");
        }
    }

    #[test]
    fn bogie_step_state_finite() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 1e6,
            damping: 2e4,
        };
        let state = [0.01, 0.0, 0.001, 0.0];
        let new_state = bogie.step(state, 0.001, 30.0, 0.002, 5e5);
        for v in new_state {
            assert!(v.is_finite(), "state component not finite: {v}");
        }
    }

    #[test]
    fn bogie_step_lateral_velocity_changes_with_restoring_force() {
        let bogie = BogieDynamics {
            mass: 5000.0,
            inertia: 3000.0,
            wheelbase: 2.5,
            spring_stiffness: 2e6,
            damping: 5e4,
        };
        // Start with positive lateral displacement; restoring force should reduce it
        let state = [0.01, 0.0, 0.0, 0.0];
        let new_state = bogie.step(state, 0.01, 30.0, 0.0, 1e6);
        // Lateral velocity should now be negative (restoring)
        assert!(
            new_state[1] < 0.0,
            "velocity should be negative (restoring): {}",
            new_state[1]
        );
    }
}
