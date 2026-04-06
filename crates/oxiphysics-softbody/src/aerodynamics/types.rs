//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{biot_savart_filament, biot_savart_velocity};
use oxiphysics_core::math::{Real, Vec3};

use super::helpers::*;

/// Added mass tensor for a rigid body moving through fluid.
///
/// When an object accelerates through a fluid, it must accelerate some of
/// the surrounding fluid with it — effectively increasing its inertia.
///
/// For simple shapes the added mass matrix is diagonal and can be computed
/// analytically.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AddedMassTensor {
    /// Diagonal added mass coefficients `[m11, m22, m33]` in kg.
    pub diagonal: [f64; 3],
}
#[allow(dead_code)]
impl AddedMassTensor {
    /// Added mass for a sphere of radius `r` in fluid of density `rho`.
    ///
    /// `m_a = 0.5 * rho * (4/3) * pi * r^3`
    pub fn sphere(radius: f64, density: f64) -> Self {
        let m = 0.5 * density * (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        Self {
            diagonal: [m, m, m],
        }
    }
    /// Added mass for a prolate spheroid with semi-axes `a ≥ b = c`.
    ///
    /// Axial added mass ≈ `C_a * rho * V` where `V` is volume and `C_a` ≈
    /// depends on aspect ratio.  Transverse ≈ `0.5 * rho * V`.
    pub fn spheroid(a: f64, b: f64, density: f64) -> Self {
        let volume = (4.0 / 3.0) * std::f64::consts::PI * a * b * b;
        let ecc = (1.0 - (b / a).powi(2)).sqrt().max(0.0);
        let alpha0 = if ecc > 1e-6 {
            2.0 * (1.0 - ecc * ecc) / (ecc * ecc * ecc)
                * (0.5 * ((1.0 + ecc) / (1.0 - ecc)).ln() - ecc)
        } else {
            2.0 / 3.0
        };
        let c_ax = alpha0 / (2.0 - alpha0);
        let m_ax = c_ax * density * volume;
        let m_tr = 0.5 * density * volume;
        Self {
            diagonal: [m_ax, m_tr, m_tr],
        }
    }
    /// Force due to added mass: `F = -m_a * a`.
    ///
    /// `acceleration` is the body acceleration in its local frame.
    pub fn force(&self, acceleration: [f64; 3]) -> [f64; 3] {
        [
            -self.diagonal[0] * acceleration[0],
            -self.diagonal[1] * acceleration[1],
            -self.diagonal[2] * acceleration[2],
        ]
    }
    /// Total added mass (sum of diagonal elements).
    pub fn total_mass(&self) -> f64 {
        self.diagonal.iter().sum()
    }
}
/// A single aerodynamic panel (e.g. one triangle of a cloth mesh).
///
/// Panels are used by [`AeroSurface`] to compute distributed aerodynamic loads.
#[derive(Debug, Clone)]
pub struct AeroPanel {
    /// Outward unit normal of the panel.
    pub normal: [f64; 3],
    /// Panel area (m^2).
    pub area: f64,
    /// Centre point of the panel.
    pub center: [f64; 3],
    /// Lift coefficient (dimensionless).
    pub cl: f64,
    /// Drag coefficient (dimensionless).
    pub cd: f64,
    /// Panel velocity (m/s) -- set each frame from the underlying particles.
    pub velocity: [f64; 3],
}
impl AeroPanel {
    /// Create a new aerodynamic panel.
    #[allow(dead_code)]
    pub fn new(normal: [f64; 3], area: f64, center: [f64; 3], cl: f64, cd: f64) -> Self {
        Self {
            normal,
            area,
            center,
            cl,
            cd,
            velocity: [0.0; 3],
        }
    }
    /// Build a panel from three vertex positions (computes normal and area).
    #[allow(dead_code)]
    pub fn from_triangle(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], cl: f64, cd: f64) -> Self {
        let e1 = v3_sub(v1, v0);
        let e2 = v3_sub(v2, v0);
        let cross = v3_cross(e1, e2);
        let area_x2 = v3_norm(cross);
        let area = area_x2 * 0.5;
        let normal = v3_normalize(cross).unwrap_or([0.0, 1.0, 0.0]);
        let center = [
            (v0[0] + v1[0] + v2[0]) / 3.0,
            (v0[1] + v1[1] + v2[1]) / 3.0,
            (v0[2] + v1[2] + v2[2]) / 3.0,
        ];
        Self {
            normal,
            area,
            center,
            cl,
            cd,
            velocity: [0.0; 3],
        }
    }
    /// Compute lift and drag forces on this panel.
    ///
    /// Returns `(lift, drag)` as `[f64; 3]` vectors.
    #[allow(dead_code)]
    pub fn lift_drag(&self, wind_velocity: [f64; 3], density: f64) -> ([f64; 3], [f64; 3]) {
        let v_rel = v3_sub(wind_velocity, self.velocity);
        let half_rho_a = 0.5 * density * self.area;
        let vn = v3_dot(v_rel, self.normal);
        let drag = v3_scale(self.normal, half_rho_a * self.cd * vn);
        let cross1 = v3_cross(v_rel, self.normal);
        let lift = v3_scale(v3_cross(cross1, self.normal), half_rho_a * self.cl);
        (lift, drag)
    }
    /// Compute pressure coefficient Cp at this panel using Bernoulli equation.
    ///
    /// `Cp = 1 - (v_local / v_inf)^2`
    ///
    /// where `v_local` is the tangential velocity component at the panel surface.
    #[allow(dead_code)]
    pub fn pressure_coefficient(&self, wind_velocity: [f64; 3]) -> f64 {
        let v_inf_sq = v3_dot(wind_velocity, wind_velocity);
        if v_inf_sq < 1e-30 {
            return 0.0;
        }
        let v_rel = v3_sub(wind_velocity, self.velocity);
        let vn = v3_dot(v_rel, self.normal);
        let v_tang = v3_sub(v_rel, v3_scale(self.normal, vn));
        let v_tang_sq = v3_dot(v_tang, v_tang);
        1.0 - v_tang_sq / v_inf_sq
    }
    /// Compute the angle of attack (radians) between the wind and the panel plane.
    ///
    /// Returns a value in `[0, pi/2]`.
    #[allow(dead_code)]
    pub fn angle_of_attack(&self, wind_velocity: [f64; 3]) -> f64 {
        let v_rel = v3_sub(wind_velocity, self.velocity);
        let speed = v3_norm(v_rel);
        if speed < 1e-12 {
            return 0.0;
        }
        let cos_alpha = (v3_dot(v_rel, self.normal) / speed).abs();
        cos_alpha.clamp(0.0, 1.0).asin()
    }
}
/// Model for a pressurised membrane (e.g. an inflatable wing or parachute
/// canopy under internal pressure).
///
/// The membrane is assumed circular (simplified) and its deformation is
/// governed by the Young–Laplace equation:
///
/// `p = 2 * T / R`
///
/// where `T` is the tension, `R` is the radius of curvature, and `p` is the
/// pressure differential.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MembraneInflation {
    /// Unstressed membrane radius (m).
    pub unstressed_radius: f64,
    /// Membrane thickness (m).
    pub thickness: f64,
    /// Elastic modulus of the membrane material (Pa).
    pub elastic_modulus: f64,
    /// Poisson's ratio of the membrane material.
    pub poisson_ratio: f64,
}
#[allow(dead_code)]
impl MembraneInflation {
    /// Create a model for a typical textile membrane (low-porosity nylon).
    pub fn nylon_canopy(radius: f64) -> Self {
        Self {
            unstressed_radius: radius,
            thickness: 0.0005,
            elastic_modulus: 3e9,
            poisson_ratio: 0.35,
        }
    }
    /// Membrane stiffness per unit area (N/m²).
    ///
    /// `k = E * t / R²` (membrane stiffness in the spherical cap approx.)
    pub fn membrane_stiffness(&self) -> f64 {
        self.elastic_modulus * self.thickness / (self.unstressed_radius * self.unstressed_radius)
    }
    /// Equilibrium radius under pressure differential `dp` (Pa).
    ///
    /// From Young–Laplace: `R = 2T / p`, `T = E * t * ε` (hoop stress).
    /// Solved iteratively; here a linearised approximation is used.
    pub fn inflated_radius(&self, dp: f64) -> f64 {
        if dp <= 0.0 {
            return self.unstressed_radius;
        }
        let strain = dp * self.unstressed_radius / (2.0 * self.elastic_modulus * self.thickness);
        self.unstressed_radius * (1.0 + strain)
    }
    /// Inflation force (N) normal to the membrane for given pressure and area.
    pub fn inflation_force(&self, dp: f64, area: f64) -> f64 {
        dp * area
    }
    /// Hoop stress (Pa) in the membrane under pressure `dp`.
    ///
    /// `σ_hoop = dp * R / (2 * t)`
    pub fn hoop_stress(&self, dp: f64) -> f64 {
        dp * self.inflated_radius(dp) / (2.0 * self.thickness)
    }
    /// Safety factor: ratio of elastic limit to current hoop stress.
    ///
    /// Elastic limit ≈ `E / 1000` (0.1% strain limit).
    pub fn safety_factor(&self, dp: f64) -> f64 {
        let limit = self.elastic_modulus / 1000.0;
        let sigma = self.hoop_stress(dp);
        if sigma < 1e-10 {
            f64::INFINITY
        } else {
            limit / sigma
        }
    }
}
/// Prandtl lifting-line theory for a finite-span wing.
///
/// Computes the spanwise-averaged lift and induced drag using an elliptic
/// load distribution approximation.
///
/// - Lift:  L = 0.5 * rho * V^2 * S * CL
/// - CL    = CL_2D / (1 + 2/AR)  (finite-wing correction)
/// - CL_2D = 2π * alpha (thin-airfoil)
/// - Induced drag: CDi = CL^2 / (π * AR * e)  (e = Oswald efficiency ≈ 1)
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct LiftingLineTheory {
    /// Wing span (m).
    pub span: f64,
    /// Mean aerodynamic chord (m).
    pub chord: f64,
    /// Angle of attack (rad).
    pub angle_of_attack: f64,
    /// Aspect ratio AR = span^2 / S.
    pub aspect_ratio: f64,
    /// Oswald span efficiency factor (0..1].
    pub oswald_e: f64,
}
#[allow(dead_code)]
#[allow(non_snake_case)]
impl LiftingLineTheory {
    /// Create a new lifting-line wing.
    ///
    /// `aspect_ratio` can be computed as `span / chord` for rectangular wings
    /// or `span^2 / area` in general.
    pub fn new(span: f64, chord: f64, angle_of_attack: f64, aspect_ratio: f64) -> Self {
        Self {
            span,
            chord,
            angle_of_attack,
            aspect_ratio,
            oswald_e: 1.0,
        }
    }
    /// Create with explicit Oswald efficiency.
    pub fn with_oswald(mut self, e: f64) -> Self {
        self.oswald_e = e.max(1e-6);
        self
    }
    /// Reference wing area S = span * chord.
    pub fn ref_area(&self) -> f64 {
        self.span * self.chord
    }
    /// Finite-wing lift coefficient CL.
    ///
    /// Corrects the 2D thin-airfoil slope (2π/rad) for finite span.
    pub fn lift_coefficient(&self) -> f64 {
        let cl_2d = 2.0 * std::f64::consts::PI * self.angle_of_attack;
        cl_2d / (1.0 + 2.0 / self.aspect_ratio.max(1e-6))
    }
    /// Induced drag coefficient CDi = CL^2 / (π AR e).
    pub fn induced_drag_coefficient(&self) -> f64 {
        let cl = self.lift_coefficient();
        cl * cl / (std::f64::consts::PI * self.aspect_ratio.max(1e-6) * self.oswald_e)
    }
    /// Compute (lift, induced_drag) forces given free-stream speed and density.
    ///
    /// Returns forces in Newtons.
    pub fn compute_lift_drag(&self, velocity: f64, density: f64) -> (f64, f64) {
        let q = 0.5 * density * velocity * velocity;
        let s = self.ref_area();
        let lift = q * s * self.lift_coefficient();
        let drag = q * s * self.induced_drag_coefficient();
        (lift, drag)
    }
}
/// A simple vortex lattice method (VLM) panel solver for rectangular wings.
///
/// Discretises the wing into `n_span × n_chord` panels, places a horseshoe
/// vortex on each, and solves for the bound circulation strengths via a
/// collocation-point boundary condition (no-penetration).
///
/// The implementation uses an influence-coefficient matrix (AIC) approach:
/// each panel's collocation point sees induced velocities from all panels'
/// horseshoe vortices.  The system Γ = AIC⁻¹ * RHS is solved by Gaussian
/// elimination.
#[derive(Debug, Clone)]
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct VortexLatticeMethod {
    /// Number of spanwise panels.
    pub n_span: usize,
    /// Number of chordwise panels.
    pub n_chord: usize,
    /// Total wing span (m).
    pub span: f64,
    /// Mean chord length (m).
    pub chord: f64,
    /// Angle of attack (rad).
    pub angle_of_attack: f64,
    /// Panel centroid positions (n_span × n_chord), world space.
    pub panel_centers: Vec<[f64; 3]>,
    /// Bound vortex quarter-chord points.
    pub bound_vortex_pts: Vec<([f64; 3], [f64; 3])>,
}
#[allow(dead_code)]
#[allow(non_snake_case)]
impl VortexLatticeMethod {
    /// Create a flat, rectangular wing centred at the origin.
    ///
    /// Span is along X, chord along Y (freestream in +Y direction after AoA rotation).
    pub fn new(n_span: usize, n_chord: usize, span: f64, chord: f64, angle_of_attack: f64) -> Self {
        let n_span = n_span.max(1);
        let n_chord = n_chord.max(1);
        let dx = span / n_span as f64;
        let dy = chord / n_chord as f64;
        let mut panel_centers = Vec::with_capacity(n_span * n_chord);
        let mut bound_vortex_pts = Vec::with_capacity(n_span * n_chord);
        for i in 0..n_span {
            for j in 0..n_chord {
                let x_c = -span * 0.5 + (i as f64 + 0.5) * dx;
                let y_c = (j as f64 + 0.25) * dy;
                panel_centers.push([x_c, y_c, 0.0]);
                let x_l = -span * 0.5 + i as f64 * dx;
                let x_r = x_l + dx;
                let y_bv = (j as f64 + 0.25) * dy;
                bound_vortex_pts.push(([x_l, y_bv, 0.0], [x_r, y_bv, 0.0]));
            }
        }
        Self {
            n_span,
            n_chord,
            span,
            chord,
            angle_of_attack,
            panel_centers,
            bound_vortex_pts,
        }
    }
    /// Total number of panels.
    pub fn n_panels(&self) -> usize {
        self.n_span * self.n_chord
    }
    /// Solve for bound circulation strengths Γ using the AIC matrix.
    ///
    /// Returns a vector of length `n_panels()` with each panel's Γ (m²/s).
    pub fn solve_gammas(&self, velocity: f64, density: f64) -> Vec<f64> {
        let np = self.n_panels();
        if np == 0 {
            return vec![];
        }
        let _density = density;
        let aoa = self.angle_of_attack;
        let wing_normal = [0.0_f64, -aoa.sin(), aoa.cos()];
        let v_inf = [0.0_f64, velocity * aoa.cos(), -velocity * aoa.sin()];
        let mut aic = vec![0.0_f64; np * np];
        let trail_len = 100.0 * self.span;
        for (row, &cp) in self.panel_centers.iter().enumerate() {
            for (col, &(bv_l, bv_r)) in self.bound_vortex_pts.iter().enumerate() {
                let v_bound = biot_savart_velocity(bv_l, bv_r, 1.0, cp);
                let trail_l_end = [bv_l[0], bv_l[1] + trail_len, bv_l[2]];
                let v_trail_l = biot_savart_velocity(bv_r, trail_l_end, 1.0, cp);
                let trail_r_start = [bv_r[0], bv_r[1] + trail_len, bv_r[2]];
                let v_trail_r = biot_savart_velocity(trail_r_start, bv_l, 1.0, cp);
                let v_tot = [
                    v_bound[0] + v_trail_l[0] + v_trail_r[0],
                    v_bound[1] + v_trail_l[1] + v_trail_r[1],
                    v_bound[2] + v_trail_l[2] + v_trail_r[2],
                ];
                let w_n = v_tot[0] * wing_normal[0]
                    + v_tot[1] * wing_normal[1]
                    + v_tot[2] * wing_normal[2];
                aic[row * np + col] = w_n;
            }
        }
        let rhs: Vec<f64> = (0..np)
            .map(|_| {
                -(v_inf[0] * wing_normal[0] + v_inf[1] * wing_normal[1] + v_inf[2] * wing_normal[2])
            })
            .collect();
        gauss_solve_opt(np, &mut { aic }, &rhs).unwrap_or_else(|| vec![0.0; np])
    }
    /// Compute total aerodynamic lift using Kutta-Joukowski (Γ-based).
    ///
    /// `L = rho * V * sum(Gamma_i * span_i)`
    pub fn total_lift(&self, velocity: f64, density: f64) -> f64 {
        let gammas = self.solve_gammas(velocity, density);
        let dx = self.span / self.n_span as f64;
        let lift: f64 = gammas.iter().take(self.n_span).sum::<f64>() * dx;
        density * velocity * lift.abs()
    }
}
/// Standard drag polar: `CD = CD0 + k·CL²`
///
/// Models the total drag of a wing or aircraft as the sum of:
/// - `CD0`: parasitic (zero-lift) drag coefficient
/// - `k·CL²`: induced drag, where `k = 1 / (π·AR·e)`
///
/// Supports the computation of:
/// - Minimum-drag CL (best L/D)
/// - Maximum lift-to-drag ratio
/// - Stall-bounded envelope
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct DragPolar {
    /// Parasitic drag coefficient (zero-lift drag).
    pub cd0: f64,
    /// Induced drag factor k = 1/(π·AR·e).
    pub k: f64,
    /// Maximum lift coefficient (stall limit).
    pub cl_max: f64,
}
#[allow(dead_code)]
impl DragPolar {
    /// Create a drag polar for a wing with given `AR` and Oswald efficiency `e`.
    ///
    /// `cd0` is the zero-lift drag coefficient.
    pub fn new(cd0: f64, aspect_ratio: f64, oswald_e: f64) -> Self {
        let k = 1.0 / (std::f64::consts::PI * aspect_ratio.max(1e-6) * oswald_e.max(1e-6));
        Self {
            cd0,
            k,
            cl_max: 1.6,
        }
    }
    /// Set a custom maximum lift coefficient.
    pub fn with_cl_max(mut self, cl_max: f64) -> Self {
        self.cl_max = cl_max;
        self
    }
    /// Total drag coefficient for a given lift coefficient.
    ///
    /// `CD = CD0 + k·CL²`
    pub fn drag_coefficient(&self, cl: f64) -> f64 {
        self.cd0 + self.k * cl * cl
    }
    /// Lift-to-drag ratio `L/D = CL / CD`.
    pub fn lift_to_drag_ratio(&self, cl: f64) -> f64 {
        let cd = self.drag_coefficient(cl);
        if cd < 1e-30 { 0.0 } else { cl / cd }
    }
    /// CL at maximum L/D (best glide ratio):
    ///
    /// `CL_best = sqrt(CD0 / k)`
    pub fn cl_best_ld(&self) -> f64 {
        (self.cd0 / self.k).sqrt()
    }
    /// Maximum lift-to-drag ratio:
    ///
    /// `(L/D)_max = 1 / (2 * sqrt(CD0 * k))`
    pub fn max_lift_to_drag(&self) -> f64 {
        let denom = 2.0 * (self.cd0 * self.k).sqrt();
        if denom < 1e-30 { 0.0 } else { 1.0 / denom }
    }
    /// Drag coefficient at minimum total drag (CD minimum = CD0 + k·CL_best²).
    pub fn cd_min_drag(&self) -> f64 {
        self.drag_coefficient(self.cl_best_ld())
    }
    /// Build the drag polar as a series of (CL, CD) points.
    ///
    /// Samples `n` points uniformly from CL = 0 to `cl_max`.
    pub fn polar_curve(&self, n: usize) -> Vec<(f64, f64)> {
        if n == 0 {
            return Vec::new();
        }
        (0..n)
            .map(|i| {
                let cl = self.cl_max * i as f64 / (n - 1).max(1) as f64;
                (cl, self.drag_coefficient(cl))
            })
            .collect()
    }
    /// Aerodynamic forces for the given `CL`, flight conditions.
    ///
    /// Returns `(lift_N, drag_N)`.
    pub fn forces(&self, cl: f64, velocity: f64, density: f64, ref_area: f64) -> (f64, f64) {
        let q = 0.5 * density * velocity * velocity;
        let lift = cl * q * ref_area;
        let drag = self.drag_coefficient(cl) * q * ref_area;
        (lift, drag)
    }
    /// Range parameter `E = sqrt(CL^3 / CD^2)` (proportional to Breguet range).
    ///
    /// Maximum occurs at CL = sqrt(3 * CD0 / k).
    pub fn range_parameter(&self, cl: f64) -> f64 {
        let cd = self.drag_coefficient(cl);
        if cd < 1e-30 {
            0.0
        } else {
            (cl.powi(3) / (cd * cd)).sqrt()
        }
    }
    /// CL for maximum range (endurance):
    ///
    /// `CL_range = sqrt(3 * CD0 / k)`
    pub fn cl_max_range(&self) -> f64 {
        (3.0 * self.cd0 / self.k).sqrt()
    }
}
/// Aerodynamic force and moment state for a wing or body.
#[derive(Debug, Clone, Copy, Default)]
#[allow(dead_code)]
pub struct AeroForces {
    /// Total lift force (N).
    pub lift: f64,
    /// Total drag force (N).
    pub drag: f64,
    /// Pitching moment about the aerodynamic centre (N·m).
    /// Positive nose-up.
    pub pitching_moment: f64,
}
#[allow(dead_code)]
impl AeroForces {
    /// Compute from non-dimensional coefficients and flight conditions.
    ///
    /// `cm` is the pitching-moment coefficient about the quarter-chord.
    pub fn from_coefficients(
        cl: f64,
        cd: f64,
        cm: f64,
        dynamic_pressure: f64,
        ref_area: f64,
        ref_chord: f64,
    ) -> Self {
        Self {
            lift: cl * dynamic_pressure * ref_area,
            drag: cd * dynamic_pressure * ref_area,
            pitching_moment: cm * dynamic_pressure * ref_area * ref_chord,
        }
    }
    /// Lift-to-drag ratio.
    pub fn lift_to_drag(&self) -> f64 {
        if self.drag.abs() < 1e-30 {
            0.0
        } else {
            self.lift / self.drag
        }
    }
}
/// A horseshoe vortex element used in the vortex lattice method.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct HorseshoeVortex {
    /// Left bound vortex endpoint.
    pub p1: [f64; 3],
    /// Right bound vortex endpoint.
    pub p2: [f64; 3],
    /// Control point where the boundary condition is enforced.
    pub control_point: [f64; 3],
    /// Normal at the control point.
    pub normal: [f64; 3],
    /// Circulation strength (Gamma).
    pub gamma: f64,
}
/// Unsteady aerodynamic model for a thin oscillating airfoil (quasi-steady
/// plus added-mass correction).
///
/// For a thin airfoil pitching at rate `αdot` and plunging at rate `hdot`, the
/// unsteady lift correction adds the "added mass" (non-circulatory) term:
///
/// `ΔCl = π · (c/2) · (αdot / V + hdot / V²)` (per unit chord)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct UnsteadyAeroPanel {
    /// Chord length (m).
    pub chord: f64,
    /// Quasi-steady lift coefficient slope (1/rad).
    pub cl_alpha: f64,
    /// Air density (kg/m³).
    pub density: f64,
}
#[allow(dead_code)]
impl UnsteadyAeroPanel {
    /// Create a thin-airfoil panel with standard lift slope (2π).
    pub fn thin_airfoil(chord: f64, density: f64) -> Self {
        Self {
            chord,
            cl_alpha: 2.0 * std::f64::consts::PI,
            density,
        }
    }
    /// Quasi-steady lift force per unit span (N/m).
    ///
    /// `L_qs = 0.5 * rho * V^2 * c * cl_alpha * alpha`
    pub fn quasi_steady_lift(&self, speed: f64, alpha_rad: f64) -> f64 {
        0.5 * self.density * speed * speed * self.chord * self.cl_alpha * alpha_rad
    }
    /// Added-mass (non-circulatory) lift force per unit span.
    ///
    /// `L_am = rho * pi * (c/2)^2 * (V * alpha_dot + h_ddot)`
    ///
    /// where `alpha_dot` is pitch rate (rad/s) and `h_ddot` is heave
    /// acceleration (m/s²).
    pub fn added_mass_lift(&self, speed: f64, alpha_dot: f64, h_ddot: f64) -> f64 {
        let r = self.chord * 0.5;
        self.density * std::f64::consts::PI * r * r * (speed * alpha_dot + h_ddot)
    }
    /// Total unsteady lift per unit span.
    pub fn total_lift(&self, speed: f64, alpha_rad: f64, alpha_dot: f64, h_ddot: f64) -> f64 {
        self.quasi_steady_lift(speed, alpha_rad) + self.added_mass_lift(speed, alpha_dot, h_ddot)
    }
    /// Pitching moment about the mid-chord per unit span.
    ///
    /// For the added-mass term only (quasi-steady moment is zero for thin airfoil
    /// about the quarter-chord):
    ///
    /// `M_am = -rho * pi * (c/2)^2 * V * c/2 * alpha_dot / 8`
    pub fn pitching_moment(&self, speed: f64, alpha_dot: f64) -> f64 {
        let r = self.chord * 0.5;
        -self.density * std::f64::consts::PI * r * r * speed * self.chord * alpha_dot / 8.0
    }
}
/// Aerodynamic stall model based on a critical angle of attack.
///
/// Below `alpha_stall` the lift coefficient follows a linear slope.
/// Above it, lift drops off and drag rises sharply — modelling flow
/// separation.  The transition is smoothed with a logistic-style blending.
///
/// The post-stall lift uses a Kirchhoff-style recovery:
///
/// `CL_post = CL_max * sin(α) * cos(α)  (≈ separated-flow limit)`
///
/// The drag polar above stall adds a separation penalty.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct StallModel {
    /// Critical (stall) angle of attack (rad).  Typically 15–20°.
    pub alpha_stall: f64,
    /// Pre-stall lift slope dCL/dα (1/rad). Thin airfoil = 2π.
    pub cl_alpha: f64,
    /// Maximum lift coefficient at stall.
    pub cl_max: f64,
    /// Zero-lift angle of attack (rad). Typically 0 for symmetric profiles.
    pub alpha_zero_lift: f64,
    /// Profile drag coefficient (parasitic drag, constant term).
    pub cd0: f64,
    /// Post-stall drag increment at 90° (fully separated).
    pub cd_max_stall: f64,
}
#[allow(dead_code)]
impl StallModel {
    /// Create a NACA-4 symmetric profile model (e.g. NACA 0012).
    pub fn naca_symmetric() -> Self {
        Self {
            alpha_stall: 16.0_f64.to_radians(),
            cl_alpha: 2.0 * std::f64::consts::PI,
            cl_max: 1.4,
            alpha_zero_lift: 0.0,
            cd0: 0.006,
            cd_max_stall: 1.8,
        }
    }
    /// Create a cambered airfoil model (e.g. NACA 2412).
    pub fn naca_cambered() -> Self {
        Self {
            alpha_stall: 14.0_f64.to_radians(),
            cl_alpha: 2.0 * std::f64::consts::PI,
            cl_max: 1.6,
            alpha_zero_lift: -2.0_f64.to_radians(),
            cd0: 0.007,
            cd_max_stall: 1.9,
        }
    }
    /// Lift coefficient at angle of attack `alpha` (rad).
    ///
    /// Blends linearly between the linear (attached-flow) and separated-flow
    /// models using a smooth sigmoid transition at the stall angle.
    pub fn lift_coefficient(&self, alpha: f64) -> f64 {
        let eff = alpha - self.alpha_zero_lift;
        let cl_linear = self.cl_alpha * eff;
        let cl_post = self.cl_max * eff.sin() * eff.cos();
        let sigma = stall_blend(eff.abs(), self.alpha_stall, 2.0_f64.to_radians());
        (1.0 - sigma) * cl_linear + sigma * cl_post
    }
    /// Drag coefficient at angle of attack `alpha` (rad).
    ///
    /// Pre-stall: induced drag via polar (k·CL²).
    /// Post-stall: additional separation drag up to `cd_max_stall` at 90°.
    pub fn drag_coefficient(&self, alpha: f64) -> f64 {
        let cl = self.lift_coefficient(alpha);
        let k = 1.0 / (std::f64::consts::PI * 8.0);
        let cd_induced = k * cl * cl;
        let eff = (alpha - self.alpha_zero_lift).abs();
        let sigma = stall_blend(eff, self.alpha_stall, 2.0_f64.to_radians());
        let cd_sep = sigma * self.cd_max_stall * eff.sin().powi(2);
        self.cd0 + cd_induced + cd_sep
    }
    /// Pitching moment coefficient about the quarter-chord.
    ///
    /// Pre-stall: Cm₀ (constant, set to 0 for symmetric).
    /// Post-stall: Kirchhoff model gives Cm → −CL/4 approximately.
    pub fn moment_coefficient(&self, alpha: f64) -> f64 {
        let eff = alpha - self.alpha_zero_lift;
        let sigma = stall_blend(eff.abs(), self.alpha_stall, 2.0_f64.to_radians());
        let cl_post = self.cl_max * eff.sin() * eff.cos();
        -sigma * cl_post * 0.25
    }
    /// Full aerodynamic forces from flight conditions.
    ///
    /// Returns `(lift_N, drag_N)` for the given reference area, velocity, and density.
    #[allow(clippy::too_many_arguments)]
    pub fn forces(&self, alpha: f64, velocity: f64, density: f64, ref_area: f64) -> (f64, f64) {
        let q = 0.5 * density * velocity * velocity;
        let lift = self.lift_coefficient(alpha) * q * ref_area;
        let drag = self.drag_coefficient(alpha) * q * ref_area;
        (lift, drag)
    }
    /// Returns `true` if the given angle of attack exceeds the stall angle.
    pub fn is_stalled(&self, alpha: f64) -> bool {
        (alpha - self.alpha_zero_lift).abs() > self.alpha_stall
    }
    /// Stall onset fraction: 0 = attached flow, 1 = fully stalled.
    pub fn stall_fraction(&self, alpha: f64) -> f64 {
        let eff = (alpha - self.alpha_zero_lift).abs();
        stall_blend(eff, self.alpha_stall, 2.0_f64.to_radians())
    }
}
/// Simple look-up table aerodynamic model mapping angle of attack (degrees)
/// to CL and CD coefficients, with linear interpolation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AeroCoeffTable {
    /// Alpha values in degrees (must be sorted ascending).
    pub alpha_deg: Vec<f64>,
    /// Lift coefficient at each alpha.
    pub cl: Vec<f64>,
    /// Drag coefficient at each alpha.
    pub cd: Vec<f64>,
}
#[allow(dead_code)]
impl AeroCoeffTable {
    /// Construct from paired data.
    pub fn new(alpha_deg: Vec<f64>, cl: Vec<f64>, cd: Vec<f64>) -> Self {
        Self { alpha_deg, cl, cd }
    }
    /// Interpolate CL at a given angle of attack (degrees).
    pub fn cl_at(&self, alpha: f64) -> f64 {
        interp_linear(&self.alpha_deg, &self.cl, alpha)
    }
    /// Interpolate CD at a given angle of attack (degrees).
    pub fn cd_at(&self, alpha: f64) -> f64 {
        interp_linear(&self.alpha_deg, &self.cd, alpha)
    }
    /// Compute lift and drag forces (N).
    ///
    /// `alpha_deg`: angle of attack in degrees.
    /// `q`: dynamic pressure (Pa).
    /// `s_ref`: reference area (m²).
    pub fn forces(&self, alpha_deg: f64, q: f64, s_ref: f64) -> (f64, f64) {
        let cl = self.cl_at(alpha_deg);
        let cd = self.cd_at(alpha_deg);
        (cl * q * s_ref, cd * q * s_ref)
    }
}
/// A collection of [`AeroPanel`]s representing a complete aerodynamic surface.
#[derive(Debug, Clone, Default)]
pub struct AeroSurface {
    /// The panels composing this surface.
    pub panels: Vec<AeroPanel>,
}
impl AeroSurface {
    /// Create an empty surface.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a panel to the surface.
    #[allow(dead_code)]
    pub fn add_panel(&mut self, panel: AeroPanel) {
        self.panels.push(panel);
    }
    /// Compute the total aerodynamic force on the surface given wind and air
    /// density.
    ///
    /// Sums `lift + drag` contributions from all panels.
    #[allow(dead_code)]
    pub fn total_force(&self, wind_velocity: [f64; 3], density: f64) -> [f64; 3] {
        let mut total = [0.0_f64; 3];
        for panel in &self.panels {
            let (lift, drag) = panel.lift_drag(wind_velocity, density);
            total = v3_add(total, v3_add(lift, drag));
        }
        total
    }
    /// Compute lift and drag coefficients for the whole surface.
    ///
    /// Returns `(CL, CD)` based on total lift/drag forces and reference area.
    #[allow(dead_code)]
    pub fn lift_drag_coefficients(
        &self,
        wind_velocity: [f64; 3],
        density: f64,
        ref_area: f64,
    ) -> (f64, f64) {
        if ref_area < 1e-30 {
            return (0.0, 0.0);
        }
        let v_sq = v3_dot(wind_velocity, wind_velocity);
        if v_sq < 1e-30 {
            return (0.0, 0.0);
        }
        let q = 0.5 * density * v_sq;
        let wind_dir = v3_normalize(wind_velocity).unwrap_or([1.0, 0.0, 0.0]);
        let mut total_lift = [0.0_f64; 3];
        let mut total_drag = [0.0_f64; 3];
        for panel in &self.panels {
            let (lift, drag) = panel.lift_drag(wind_velocity, density);
            total_lift = v3_add(total_lift, lift);
            total_drag = v3_add(total_drag, drag);
        }
        let drag_mag = v3_dot(v3_add(total_lift, total_drag), wind_dir);
        let total_force = v3_add(total_lift, total_drag);
        let lift_vec = v3_sub(total_force, v3_scale(wind_dir, drag_mag));
        let lift_mag = v3_norm(lift_vec);
        let cd_val = drag_mag / (q * ref_area);
        let cl_val = lift_mag / (q * ref_area);
        (cl_val, cd_val)
    }
    /// Pressure distribution across all panels as Cp values.
    #[allow(dead_code)]
    pub fn pressure_distribution(&self, wind_velocity: [f64; 3]) -> Vec<f64> {
        self.panels
            .iter()
            .map(|p| p.pressure_coefficient(wind_velocity))
            .collect()
    }
    /// Total wetted area of the surface.
    #[allow(dead_code)]
    pub fn total_area(&self) -> f64 {
        self.panels.iter().map(|p| p.area).sum()
    }
    /// Compute the aerodynamic center (area-weighted centroid).
    #[allow(dead_code)]
    pub fn aero_center(&self) -> [f64; 3] {
        let total_a = self.total_area();
        if total_a < 1e-30 {
            return [0.0; 3];
        }
        let mut center = [0.0_f64; 3];
        for p in &self.panels {
            center = v3_add(center, v3_scale(p.center, p.area));
        }
        v3_scale(center, 1.0 / total_a)
    }
}
/// Physical aerodynamics parameters.
#[derive(Debug, Clone)]
pub struct AerodynamicsModel {
    /// Air density (kg/m^3), typical value: 1.225.
    pub rho_air: Real,
    /// Drag coefficient (dimensionless), typical range: 1.0--2.0.
    pub cd: Real,
    /// Lift coefficient (dimensionless), typical range: 0.5--1.0.
    pub cl: Real,
}
impl AerodynamicsModel {
    /// Create a new aerodynamics model.
    pub fn new(rho_air: Real, cd: Real, cl: Real) -> Self {
        Self { rho_air, cd, cl }
    }
    /// Compute aerodynamic forces on a triangle given the ambient wind velocity.
    ///
    /// Returns the force contribution `(f0, f1, f2)` for each vertex (equal
    /// thirds of the total triangle force).
    #[allow(clippy::too_many_arguments)]
    pub fn triangle_force(
        &self,
        v0: Vec3,
        v1: Vec3,
        v2: Vec3,
        vel0: Vec3,
        vel1: Vec3,
        vel2: Vec3,
        wind_velocity: Vec3,
    ) -> (Vec3, Vec3, Vec3) {
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let normal = e1.cross(&e2);
        let area_x2 = normal.norm();
        if area_x2 < 1e-12 {
            let zero = Vec3::zeros();
            return (zero, zero, zero);
        }
        let n_hat = normal / area_x2;
        let area = area_x2 * 0.5;
        let v_panel = (vel0 + vel1 + vel2) / 3.0;
        let v_rel = wind_velocity - v_panel;
        let half_rho_a = 0.5 * self.rho_air * area;
        let v_rel_dot_n = v_rel.dot(&n_hat);
        let f_drag = n_hat * (half_rho_a * self.cd * v_rel_dot_n);
        let cross1 = v_rel.cross(&n_hat);
        let f_lift = cross1.cross(&n_hat) * (half_rho_a * self.cl);
        let f_total = f_drag + f_lift;
        let f_each = f_total / 3.0;
        (f_each, f_each, f_each)
    }
}
/// Simple one-way fluid–structure interaction (FSI) coupling.
///
/// The fluid exerts aero forces on the structure; the structural deformation
/// modifies the panel geometry.  This struct tracks the panel pressures and
/// deformation history for one update cycle.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FsiCoupling {
    /// Aero pressure at each panel (Pa).
    pub pressures: Vec<f64>,
    /// Structural deformation magnitudes per panel (m, normal direction).
    pub deformations: Vec<f64>,
    /// Structural stiffness per panel (N/m).
    pub stiffness: Vec<f64>,
    /// Damping coefficient per panel (N·s/m).
    pub damping: Vec<f64>,
    /// Velocity of each panel in the deformation direction (m/s).
    pub velocities: Vec<f64>,
}
#[allow(dead_code)]
impl FsiCoupling {
    /// Create an FSI coupling with `n` panels.
    pub fn new(n: usize, stiffness: f64, damping: f64) -> Self {
        Self {
            pressures: vec![0.0; n],
            deformations: vec![0.0; n],
            stiffness: vec![stiffness; n],
            damping: vec![damping; n],
            velocities: vec![0.0; n],
        }
    }
    /// Update panel deformations given aero pressures and panel areas.
    ///
    /// Uses a simple spring-damper: `m * a = p * A - k * x - c * ẋ`.
    /// With panel mass `m` per panel.
    pub fn step(&mut self, areas: &[f64], panel_mass: f64, dt: f64) {
        let n = self
            .deformations
            .len()
            .min(areas.len())
            .min(self.pressures.len());
        for i in 0..n {
            let force = self.pressures[i] * areas[i]
                - self.stiffness[i] * self.deformations[i]
                - self.damping[i] * self.velocities[i];
            let accel = if panel_mass > 1e-20 {
                force / panel_mass
            } else {
                0.0
            };
            self.velocities[i] += accel * dt;
            self.deformations[i] += self.velocities[i] * dt;
        }
    }
    /// Set aero pressures from a Cp distribution.
    pub fn set_pressures_from_cp(&mut self, cp: &[f64], dynamic_pressure: f64) {
        let n = self.pressures.len().min(cp.len());
        for i in 0..n {
            self.pressures[i] = cp[i] * dynamic_pressure;
        }
    }
    /// Maximum deformation in the structure.
    pub fn max_deformation(&self) -> f64 {
        self.deformations.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Strain energy stored in the structural springs.
    ///
    /// `U = 0.5 * sum_i( k_i * x_i^2 )`
    pub fn strain_energy(&self) -> f64 {
        self.deformations
            .iter()
            .zip(self.stiffness.iter())
            .map(|(&x, &k)| 0.5 * k * x * x)
            .sum()
    }
}
/// Ground effect: increase in lift and decrease in induced drag when a wing
/// flies close to the ground.
///
/// Uses the classic out-of-ground-effect (OGE) / in-ground-effect (IGE) ratio
/// method by Eq. (7.16) from Raymer (2018) / NASA-TM-1998-208400:
///
/// `CLi_IGE / CLi_OGE ≈ 1 − σ_ground`
///
/// where σ_ground is derived from the wing height above the ground.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct GroundEffect {
    /// Wing semi-span (m).
    pub semi_span: f64,
    /// Wing aspect ratio.
    pub aspect_ratio: f64,
}
#[allow(dead_code)]
impl GroundEffect {
    /// Create a ground-effect model for a given wing geometry.
    pub fn new(semi_span: f64, aspect_ratio: f64) -> Self {
        Self {
            semi_span,
            aspect_ratio,
        }
    }
    /// Ground-effect factor φ (ratio of induced drag coefficient IGE/OGE).
    ///
    /// `φ = (h / (16 * b))^2 / (1 + (h / (16 * b))^2)`
    ///
    /// where `h` is the height above ground (m) and `b = 2 * semi_span`.
    ///
    /// Returns a value in \[0, 1\].  At h → 0, φ → 0 (no induced drag).
    /// At h → ∞, φ → 1 (OGE, no effect).
    pub fn induced_drag_factor(&self, height_above_ground: f64) -> f64 {
        let b = 2.0 * self.semi_span;
        let ratio = height_above_ground / (16.0 * b);
        let ratio2 = ratio * ratio;
        ratio2 / (1.0 + ratio2)
    }
    /// Ground-effect lift bonus: ratio of CL_IGE / CL_OGE.
    ///
    /// Lift increases near the ground because the downwash angle is reduced.
    /// A simplified expression from Torenbeek:
    ///
    /// `CL_IGE / CL_OGE ≈ 1 + δ_ge`  where  `δ_ge ≈ (b / (2h))^2 / AR`
    pub fn lift_ratio(&self, height_above_ground: f64) -> f64 {
        let h = height_above_ground.max(1e-6);
        let b = 2.0 * self.semi_span;
        let delta_ge = (b / (2.0 * h)).powi(2) / self.aspect_ratio.max(1.0);
        (1.0 + delta_ge).min(2.0)
    }
    /// Compute effective lift and induced-drag coefficients in ground effect.
    ///
    /// `cl_oge` and `cdi_oge` are the coefficients out-of-ground-effect.
    pub fn apply(&self, cl_oge: f64, cdi_oge: f64, height: f64) -> (f64, f64) {
        let phi = self.induced_drag_factor(height);
        let lr = self.lift_ratio(height);
        let cl_ige = cl_oge * lr;
        let cdi_ige = cdi_oge * phi;
        (cl_ige, cdi_ige)
    }
    /// Height at which ground effect becomes negligible (< 1 % correction).
    ///
    /// Returns the height where `1 − φ < 0.01` (i.e. φ > 0.99).
    pub fn negligible_height(&self) -> f64 {
        let b = 2.0 * self.semi_span;
        9.95 * 16.0 * b
    }
}
/// A rectangular vortex-ring panel for the doublet-lattice / full VLM wake
/// representation.
///
/// The panel has four corner nodes (counterclockwise viewed from above) and a
/// bound-vortex segment running span-wise along the quarter-chord line.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct VortexRingPanel {
    /// Four corner nodes of the panel \[front-left, front-right, rear-right, rear-left\].
    pub corners: [[f64; 3]; 4],
    /// Circulation strength Γ.
    pub gamma: f64,
}
#[allow(dead_code)]
impl VortexRingPanel {
    /// Construct a planar rectangular panel from leading-edge midpoint, span,
    /// chord and sweep angle (radians).
    pub fn new_rect(le_mid: [f64; 3], half_span: f64, chord: f64, sweep_rad: f64) -> Self {
        let dx_sweep = chord * sweep_rad.sin();
        let dz_chord = chord * sweep_rad.cos();
        let fl = [le_mid[0] - half_span, le_mid[1], le_mid[2]];
        let fr = [le_mid[0] + half_span, le_mid[1], le_mid[2]];
        let rr = [fr[0] + dx_sweep, fr[1] + dz_chord, fr[2]];
        let rl = [fl[0] + dx_sweep, fl[1] + dz_chord, fl[2]];
        Self {
            corners: [fl, fr, rr, rl],
            gamma: 0.0,
        }
    }
    /// Compute the collocation (control) point at the three-quarter chord.
    pub fn collocation_point(&self) -> [f64; 3] {
        let [c0, c1, c2, c3] = self.corners;
        let mid_le = v3_scale(v3_add(c0, c1), 0.5);
        let mid_te = v3_scale(v3_add(c2, c3), 0.5);
        v3_add(v3_scale(mid_le, 0.25), v3_scale(mid_te, 0.75))
    }
    /// Unit normal to the panel (using cross product of diagonals).
    pub fn normal(&self) -> [f64; 3] {
        let [c0, c1, c2, c3] = self.corners;
        let d1 = v3_sub(c2, c0);
        let d2 = v3_sub(c3, c1);
        let n = v3_cross(d1, d2);
        v3_normalize(n).unwrap_or([0.0, 1.0, 0.0])
    }
    /// Panel area (split into two triangles).
    pub fn area(&self) -> f64 {
        let [c0, c1, c2, c3] = self.corners;
        let t1_area = {
            let a = v3_sub(c1, c0);
            let b = v3_sub(c2, c0);
            v3_norm(v3_cross(a, b)) * 0.5
        };
        let t2_area = {
            let a = v3_sub(c2, c0);
            let b = v3_sub(c3, c0);
            v3_norm(v3_cross(a, b)) * 0.5
        };
        t1_area + t2_area
    }
    /// Induced velocity at point `p` due to a vortex ring of circulation `gamma`.
    ///
    /// Uses the Biot-Savart law applied to all four edge segments.
    pub fn induced_velocity_at(&self, p: [f64; 3], gamma: f64) -> [f64; 3] {
        let [c0, c1, c2, c3] = self.corners;
        let edges = [(c0, c1), (c1, c2), (c2, c3), (c3, c0)];
        let mut vel = [0.0_f64; 3];
        for (a, b) in &edges {
            let dv = biot_savart_filament(*a, *b, p, gamma);
            for d in 0..3 {
                vel[d] += dv[d];
            }
        }
        vel
    }
}
/// Simplified flexible wing model that couples aeroelastic bending with
/// aero-load feedback.
///
/// The wing is discretised into `n_sections` along the span.  Each section
/// has an independent twist angle that changes due to aerodynamic moments.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FlexibleWing {
    /// Number of spanwise sections.
    pub n_sections: usize,
    /// Semi-span length (m).
    pub semi_span: f64,
    /// Chord length (m, assumed constant).
    pub chord: f64,
    /// Torsional stiffness per section (N·m/rad).
    pub torsional_stiffness: f64,
    /// Zero-lift angle of attack per section (rad).
    pub alpha_zero: Vec<f64>,
    /// Current elastic twist angles per section (rad).
    pub twist: Vec<f64>,
}
#[allow(dead_code)]
impl FlexibleWing {
    /// Create a flexible wing with `n_sections` spanwise sections.
    pub fn new(n_sections: usize, semi_span: f64, chord: f64, torsional_stiffness: f64) -> Self {
        Self {
            n_sections,
            semi_span,
            chord,
            torsional_stiffness,
            alpha_zero: vec![0.0; n_sections],
            twist: vec![0.0; n_sections],
        }
    }
    /// Section area (m²) per section.
    pub fn section_area(&self) -> f64 {
        self.chord * self.semi_span / self.n_sections as f64
    }
    /// Compute aerodynamic moment at each section.
    ///
    /// Uses thin-airfoil theory: `M_aero = 0.5 * rho * V^2 * c^2 * Cm_alpha * alpha_eff`
    /// where `Cm_alpha ≈ -π/4` for a thin airfoil about the quarter-chord.
    pub fn aero_moments(&self, speed: f64, density: f64, base_alpha: f64) -> Vec<f64> {
        let cm_alpha = -std::f64::consts::PI / 4.0;
        let q = 0.5 * density * speed * speed;
        self.twist
            .iter()
            .zip(self.alpha_zero.iter())
            .map(|(&twist, &az)| {
                let alpha_eff = base_alpha + twist - az;
                q * self.chord * self.chord * cm_alpha * alpha_eff
            })
            .collect()
    }
    /// Update the elastic twist angles by one time step.
    ///
    /// Section inertia `J` (kg·m²) is assumed equal for all sections.
    pub fn step(&mut self, speed: f64, density: f64, base_alpha: f64, j_section: f64, dt: f64) {
        let moments = self.aero_moments(speed, density, base_alpha);
        for (i, &m_aero) in moments.iter().enumerate() {
            let restoring = -self.torsional_stiffness * self.twist[i];
            let alpha_ddot = if j_section > 1e-20 {
                (m_aero + restoring) / j_section
            } else {
                0.0
            };
            self.twist[i] += alpha_ddot * dt * dt;
        }
    }
    /// Effective lift coefficient for the wing integrated over all sections.
    pub fn total_lift_coefficient(&self, base_alpha: f64) -> f64 {
        let cl_alpha = 2.0 * std::f64::consts::PI;
        let cl_sum: f64 = self
            .twist
            .iter()
            .zip(self.alpha_zero.iter())
            .map(|(&twist, &az)| cl_alpha * (base_alpha + twist - az))
            .sum();
        cl_sum / self.n_sections as f64
    }
}
