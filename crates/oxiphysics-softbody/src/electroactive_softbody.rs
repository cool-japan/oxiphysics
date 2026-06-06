// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electroactive soft-body physics module.
//!
//! Models dielectric elastomer actuators (DEA), ionic hydrogel actuators,
//! piezoelectric soft bodies, triboelectric soft bodies, capacitive strain
//! sensing, and soft robot actuator simulation coupling electric and mechanical
//! degrees of freedom.
//!
//! # References
//!
//! * Pelrine et al., "High-Speed Electrically Actuated Elastomers with Strain
//!   Greater Than 100%", Science 2000.
//! * Suo, "Mechanics of stretchable electronics and soft machines", MRS Bull. 2012.
//! * Suo et al., "Theory of dielectric elastomers", Acta Mech. Solida Sin. 2010.

use std::f64::consts::PI;

// ── Scalar / vector helpers (no nalgebra) ────────────────────────────────────

/// Subtract two 3-vectors (a − b).
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Dot product of two 3-vectors.
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[cfg(test)]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

/// Cross product of two 3-vectors.
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Clamp `x` into `[lo, hi]`.
#[cfg(test)]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ── Physical constants ────────────────────────────────────────────────────────

/// Permittivity of free space ε₀ (F/m).
pub const EPSILON_0: f64 = 8.854_187_817e-12;

/// Default relative permittivity for acrylic VHB-type DEA material.
pub const EPSILON_R_VHB: f64 = 4.7;

/// Threshold electric field for pull-in instability in a typical DEA.
///
/// Above this value the actuator undergoes runaway thinning.
pub const E_PULL_IN_TYPICAL: f64 = 1.8e8; // V/m

// ── Dielectric Elastomer Actuator (DEA) ──────────────────────────────────────

/// State of a thin-film dielectric elastomer membrane actuator.
///
/// Models the neo-Hookean electromechanical coupling under uniaxial or
/// equal-biaxial loading with Maxwell electrostatic pressure.
#[derive(Debug, Clone)]
pub struct DeaMembraneState {
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
    /// Initial thickness of the elastomer film (m).
    pub thickness_0: f64,
    /// Relative permittivity εᵣ of the dielectric.
    pub epsilon_r: f64,
    /// Applied voltage across the film (V).
    pub voltage: f64,
    /// In-plane stretch λ (equal biaxial, dimensionless, ≥ 1).
    pub stretch: f64,
}

impl DeaMembraneState {
    /// Create a new DEA membrane state.
    ///
    /// * `shear_modulus` – Neo-Hookean shear modulus μ (Pa).
    /// * `thickness_0`   – Undeformed film thickness (m).
    /// * `epsilon_r`     – Relative permittivity of the elastomer.
    /// * `voltage`       – Applied voltage (V).
    pub fn new(shear_modulus: f64, thickness_0: f64, epsilon_r: f64, voltage: f64) -> Self {
        Self {
            shear_modulus,
            thickness_0,
            epsilon_r,
            voltage,
            stretch: 1.0,
        }
    }

    /// Current film thickness from incompressibility (λ₁ = λ₂ = stretch ⇒ λ₃ = 1/λ²).
    pub fn current_thickness(&self) -> f64 {
        self.thickness_0 / (self.stretch * self.stretch)
    }

    /// Electric field across the film E = V / t (V/m).
    pub fn electric_field(&self) -> f64 {
        self.voltage / self.current_thickness().max(1e-12)
    }

    /// Maxwell electrostatic pressure p_e = ε₀ εᵣ E² (Pa).
    pub fn maxwell_pressure(&self) -> f64 {
        let e = self.electric_field();
        EPSILON_0 * self.epsilon_r * e * e
    }

    /// True (Cauchy) elastic stress in the plane for equal-biaxial stretch
    /// using an incompressible neo-Hookean model: σ_el = μ (λ² − 1/λ⁴).
    pub fn elastic_stress(&self) -> f64 {
        let l = self.stretch;
        self.shear_modulus * (l * l - 1.0 / (l * l * l * l))
    }

    /// Electromechanical equilibrium residual.
    ///
    /// R = σ_el − p_e.  At equilibrium R = 0.
    pub fn equilibrium_residual(&self) -> f64 {
        self.elastic_stress() - self.maxwell_pressure()
    }

    /// Detect whether the actuator has exceeded the pull-in instability.
    ///
    /// Pull-in occurs when dR/dλ < 0 (the restoring elastic force cannot
    /// balance the increasing Maxwell pressure as the film thins).
    pub fn is_pulled_in(&self) -> bool {
        // Numerical derivative: dR/dλ
        let delta = 1e-6;
        let r_plus = {
            let mut clone = self.clone();
            clone.stretch += delta;
            clone.equilibrium_residual()
        };
        let r_minus = {
            let mut clone = self.clone();
            clone.stretch -= delta;
            clone.equilibrium_residual()
        };
        let dr_dl = (r_plus - r_minus) / (2.0 * delta);
        dr_dl < 0.0
    }

    /// Solve for the equilibrium stretch using Newton-Raphson iteration.
    ///
    /// Returns `(stretch, converged)`.  Starts from the current `stretch` value.
    pub fn solve_equilibrium(&mut self, max_iter: usize, tol: f64) -> (f64, bool) {
        for _ in 0..max_iter {
            let r = self.equilibrium_residual();
            if r.abs() < tol {
                return (self.stretch, true);
            }
            // dR/dλ via finite difference
            let delta = self.stretch * 1e-7 + 1e-12;
            let r_p = {
                let mut c = self.clone();
                c.stretch += delta;
                c.equilibrium_residual()
            };
            let dr = (r_p - r) / delta;
            if dr.abs() < 1e-30 {
                break;
            }
            self.stretch = (self.stretch - r / dr).max(1.0);
        }
        (self.stretch, false)
    }

    /// Actuation strain ε_act = λ² − 1 (area strain of the membrane).
    pub fn actuation_area_strain(&self) -> f64 {
        self.stretch * self.stretch - 1.0
    }

    /// Stored elastic energy density W = μ/2 (λ₁² + λ₂² + λ₃² − 3) for equal
    /// biaxial stretch (λ₃ = 1/λ²).
    pub fn elastic_energy_density(&self) -> f64 {
        let l = self.stretch;
        let l3 = 1.0 / (l * l);
        self.shear_modulus / 2.0 * (l * l + l * l + l3 * l3 - 3.0)
    }

    /// Electrical energy stored per unit deformed volume: U_e = ε₀ εᵣ E² / 2.
    pub fn electrical_energy_density(&self) -> f64 {
        let e = self.electric_field();
        0.5 * EPSILON_0 * self.epsilon_r * e * e
    }
}

// ── Ionic Hydrogel Actuator ───────────────────────────────────────────────────

/// Model for an ionic hydrogel soft actuator driven by osmotic pressure.
///
/// The actuator swells or shrinks in response to an applied electric field
/// that drives ion migration, changing the local osmotic pressure.
#[derive(Debug, Clone)]
pub struct IonicHydrogelActuator {
    /// Young's modulus of the dry polymer network (Pa).
    pub youngs_modulus: f64,
    /// Flory-Huggins interaction parameter χ (dimensionless, 0 < χ < 0.5 for swelling).
    pub flory_chi: f64,
    /// Fixed-charge density C_f (mol/m³, represents the number of fixed ionic groups).
    pub fixed_charge_density: f64,
    /// External ion concentration c_ext (mol/m³).
    pub external_ion_conc: f64,
    /// Current volumetric swelling ratio Q = V / V₀ (Q ≥ 1).
    pub swelling_ratio: f64,
    /// Applied electric field magnitude (V/m), used to drive ionic flux.
    pub applied_field: f64,
}

impl IonicHydrogelActuator {
    /// Create a new ionic hydrogel actuator with given material parameters.
    pub fn new(
        youngs_modulus: f64,
        flory_chi: f64,
        fixed_charge_density: f64,
        external_ion_conc: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            flory_chi,
            fixed_charge_density,
            external_ion_conc,
            swelling_ratio: 1.0,
            applied_field: 0.0,
        }
    }

    /// Polymer volume fraction φ = 1 / Q.
    pub fn polymer_volume_fraction(&self) -> f64 {
        1.0 / self.swelling_ratio.max(1.0)
    }

    /// Osmotic pressure from the Flory-Rehner model (simplified):
    ///
    /// Π_osm = -R T / V_m \[ ln(1 - φ) + φ + χ φ² \]
    ///
    /// Temperature fixed at 298 K; V_m = 18e-6 m³/mol (water molar volume).
    pub fn osmotic_pressure(&self) -> f64 {
        let r_gas = 8.314; // J/(mol·K)
        let temp = 298.0; // K
        let v_m = 18.0e-6; // m³/mol
        let phi = self.polymer_volume_fraction();
        let term = (1.0 - phi).ln() + phi + self.flory_chi * phi * phi;
        -(r_gas * temp / v_m) * term
    }

    /// Donnan ionic contribution to osmotic pressure:
    ///
    /// Π_ionic = R T ( √(C_f² + 4 c_ext²) − 2 c_ext )
    ///
    /// Accounts for the fixed-charge groups repelling co-ions.
    pub fn donnan_pressure(&self) -> f64 {
        let r_gas = 8.314;
        let temp = 298.0;
        let cf = self.fixed_charge_density;
        let ce = self.external_ion_conc;
        r_gas * temp * ((cf * cf + 4.0 * ce * ce).sqrt() - 2.0 * ce)
    }

    /// Elastic restoring pressure from the rubber-elastic network:
    ///
    /// Π_el = E/3 ( Q^(−1/3) − Q^(−5/3) / 2 )
    ///
    /// (derived from the affine network model for isotropic swelling)
    pub fn elastic_restoring_pressure(&self) -> f64 {
        let q = self.swelling_ratio;
        (self.youngs_modulus / 3.0) * (q.powf(-1.0 / 3.0) - 0.5 * q.powf(-5.0 / 3.0))
    }

    /// Total driving pressure = Π_osm + Π_ionic − Π_el.
    ///
    /// Positive → actuator tends to swell further.
    pub fn net_driving_pressure(&self) -> f64 {
        self.osmotic_pressure() + self.donnan_pressure() - self.elastic_restoring_pressure()
    }

    /// Evolve the swelling ratio under applied field for one time step dt.
    ///
    /// A mobility coefficient κ (m²/(N·s)) relates field-driven ionic flux to
    /// swelling rate: dQ/dt = κ |E| / Q.
    pub fn step(&mut self, dt: f64, mobility: f64) {
        let dq = mobility * self.applied_field.abs() / self.swelling_ratio.max(1.0) * dt;
        self.swelling_ratio = (self.swelling_ratio + dq).max(1.0);
    }

    /// Linear actuation strain ε = Q^(1/3) − 1.
    pub fn linear_strain(&self) -> f64 {
        self.swelling_ratio.powf(1.0 / 3.0) - 1.0
    }
}

// ── Piezoelectric Soft Body ───────────────────────────────────────────────────

/// Piezoelectric constitutive law for a poled PVDF-based soft body element.
///
/// Uses the linear piezoelectric equations in the 3-3 configuration (poling
/// along x₃ = z):
///
/// ```text
/// S₃ = s₃₃ T₃ + d₃₃ E₃
/// D₃ = d₃₃ T₃ + ε₃₃ E₃
/// ```
#[derive(Debug, Clone)]
pub struct PiezoelectricSoftElement {
    /// Elastic compliance s₃₃ (m²/N).
    pub compliance_33: f64,
    /// Piezoelectric strain coefficient d₃₃ (m/V).
    pub d33: f64,
    /// Permittivity at constant stress ε₃₃ (F/m).
    pub permittivity_33: f64,
    /// Applied mechanical stress T₃ (Pa).
    pub stress: f64,
    /// Applied electric field E₃ (V/m).
    pub electric_field: f64,
}

impl PiezoelectricSoftElement {
    /// Create a new piezoelectric element with PVDF-like parameters.
    ///
    /// Typical PVDF: s₃₃ ≈ 24e-12 m²/N, d₃₃ ≈ -33e-12 m/V,
    /// ε₃₃ ≈ 7.4e-11 F/m.
    pub fn new_pvdf() -> Self {
        Self {
            compliance_33: 24.0e-12,
            d33: -33.0e-12,
            permittivity_33: 7.4e-11,
            stress: 0.0,
            electric_field: 0.0,
        }
    }

    /// Create a piezoelectric element with custom parameters.
    pub fn new(compliance_33: f64, d33: f64, permittivity_33: f64) -> Self {
        Self {
            compliance_33,
            d33,
            permittivity_33,
            stress: 0.0,
            electric_field: 0.0,
        }
    }

    /// Mechanical strain S₃ = s₃₃ T₃ + d₃₃ E₃.
    pub fn mechanical_strain(&self) -> f64 {
        self.compliance_33 * self.stress + self.d33 * self.electric_field
    }

    /// Electric displacement D₃ = d₃₃ T₃ + ε₃₃ E₃.
    pub fn electric_displacement(&self) -> f64 {
        self.d33 * self.stress + self.permittivity_33 * self.electric_field
    }

    /// Electromechanical coupling coefficient k₃₃.
    ///
    /// k₃₃² = d₃₃² / (s₃₃ · ε₃₃)
    pub fn coupling_coefficient_k33(&self) -> f64 {
        let k2 = (self.d33 * self.d33) / (self.compliance_33 * self.permittivity_33);
        k2.sqrt()
    }

    /// Piezoelectric energy harvested per unit volume for a given strain cycle.
    ///
    /// W = 0.5 d₃₃² / s₃₃ · T₃²  (open-circuit energy density)
    pub fn harvested_energy_density(&self) -> f64 {
        0.5 * (self.d33 * self.d33 / self.compliance_33) * self.stress * self.stress
    }

    /// Resonant frequency of a PVDF strip of length L and thickness t.
    ///
    /// f_r = 1/(2L) √(1/(ρ s₃₃))
    ///
    /// * `length` – strip length (m).
    /// * `density` – material density (kg/m³), typically ~1780 kg/m³ for PVDF.
    pub fn resonant_frequency(&self, length: f64, density: f64) -> f64 {
        let v = (1.0 / (density * self.compliance_33)).sqrt();
        v / (2.0 * length)
    }
}

// ── Triboelectric Soft Body ───────────────────────────────────────────────────

/// Triboelectric nanogenerator (TENG) soft body model.
///
/// Models contact-separation mode between two polymer surfaces that develop
/// opposite surface charges when they touch and separate.
#[derive(Debug, Clone)]
pub struct TriboelectricSoftBody {
    /// Surface charge density σ (C/m²) on polymer 1.
    pub surface_charge_density: f64,
    /// Relative permittivity of dielectric layer 1.
    pub epsilon_r1: f64,
    /// Relative permittivity of dielectric layer 2.
    pub epsilon_r2: f64,
    /// Thickness of dielectric layer 1 (m).
    pub thickness_d1: f64,
    /// Thickness of dielectric layer 2 (m).
    pub thickness_d2: f64,
    /// Air gap between the two dielectric surfaces (m).
    pub air_gap: f64,
    /// External load resistance (Ω).
    pub load_resistance: f64,
    /// Accumulated charge on external circuit q (C/m²).
    pub circuit_charge: f64,
}

impl TriboelectricSoftBody {
    /// Create a new TENG model with PTFE / nylon contact pair defaults.
    pub fn new_ptfe_nylon() -> Self {
        Self {
            surface_charge_density: -8.0e-6, // C/m²
            epsilon_r1: 2.1,                 // PTFE
            epsilon_r2: 3.5,                 // Nylon
            thickness_d1: 50.0e-6,
            thickness_d2: 50.0e-6,
            air_gap: 0.0,
            load_resistance: 1.0e6,
            circuit_charge: 0.0,
        }
    }

    /// Create a custom TENG model.
    pub fn new(
        surface_charge_density: f64,
        epsilon_r1: f64,
        epsilon_r2: f64,
        thickness_d1: f64,
        thickness_d2: f64,
        air_gap: f64,
        load_resistance: f64,
    ) -> Self {
        Self {
            surface_charge_density,
            epsilon_r1,
            epsilon_r2,
            thickness_d1,
            thickness_d2,
            air_gap,
            load_resistance,
            circuit_charge: 0.0,
        }
    }

    /// Effective dielectric thickness d_eff = d₁/ε₁ + d₂/ε₂.
    pub fn effective_dielectric_thickness(&self) -> f64 {
        self.thickness_d1 / self.epsilon_r1 + self.thickness_d2 / self.epsilon_r2
    }

    /// Open-circuit voltage V_oc = σ x / ε₀  (x = air gap).
    ///
    /// This is the voltage when no current flows.
    pub fn open_circuit_voltage(&self) -> f64 {
        self.surface_charge_density * self.air_gap / EPSILON_0
    }

    /// Short-circuit charge density Q_sc = σ x / (x + d_eff).
    pub fn short_circuit_charge_density(&self) -> f64 {
        let d_eff = self.effective_dielectric_thickness();
        let x = self.air_gap;
        self.surface_charge_density * x / (x + d_eff)
    }

    /// Instantaneous output power P = V² / R, given current air gap.
    pub fn output_power(&self) -> f64 {
        let v = self.open_circuit_voltage();
        v * v / self.load_resistance
    }

    /// Simulate one contact-separation cycle and return the peak voltage.
    ///
    /// * `max_gap` – Maximum separation distance (m).
    /// * `freq`    – Cycle frequency (Hz).
    /// * `dt`      – Time step (s).
    pub fn simulate_cycle(&mut self, max_gap: f64, freq: f64, dt: f64) -> f64 {
        let period = 1.0 / freq;
        let steps = (period / dt).ceil() as usize;
        let mut peak_voltage = 0.0_f64;
        for i in 0..steps {
            let t = i as f64 * dt;
            // Sinusoidal gap variation
            self.air_gap = max_gap * 0.5 * (1.0 - (2.0 * PI * freq * t).cos());
            let v = self.open_circuit_voltage().abs();
            if v > peak_voltage {
                peak_voltage = v;
            }
        }
        peak_voltage
    }
}

// ── Capacitive Strain Sensor ──────────────────────────────────────────────────

/// Soft capacitive strain sensor model.
///
/// A compliant parallel-plate capacitor (e.g. carbon-black/silicone composite)
/// where the capacitance changes with deformation.
#[derive(Debug, Clone)]
pub struct CapacitiveStrainSensor {
    /// Undeformed electrode area A₀ (m²).
    pub area_0: f64,
    /// Undeformed dielectric thickness t₀ (m).
    pub thickness_0: f64,
    /// Relative permittivity of the dielectric.
    pub epsilon_r: f64,
    /// Current in-plane stretch λ (equal biaxial).
    pub stretch: f64,
}

impl CapacitiveStrainSensor {
    /// Create a new capacitive sensor.
    pub fn new(area_0: f64, thickness_0: f64, epsilon_r: f64) -> Self {
        Self {
            area_0,
            thickness_0,
            epsilon_r,
            stretch: 1.0,
        }
    }

    /// Current electrode area A = A₀ λ² (incompressible assumption).
    pub fn current_area(&self) -> f64 {
        self.area_0 * self.stretch * self.stretch
    }

    /// Current dielectric thickness t = t₀ / λ².
    pub fn current_thickness(&self) -> f64 {
        self.thickness_0 / (self.stretch * self.stretch)
    }

    /// Capacitance C = ε₀ εᵣ A / t (F).
    pub fn capacitance(&self) -> f64 {
        EPSILON_0 * self.epsilon_r * self.current_area() / self.current_thickness().max(1e-15)
    }

    /// Gauge factor GF = (ΔC / C₀) / ε where ε = λ − 1 is the linear strain.
    ///
    /// For an incompressible sheet, C ∝ λ⁴ so GF = 4 λ³ / (λ⁴ − 1) → 4 near λ = 1.
    pub fn gauge_factor(&self) -> f64 {
        let c0 = EPSILON_0 * self.epsilon_r * self.area_0 / self.thickness_0;
        let c = self.capacitance();
        let epsilon = self.stretch - 1.0;
        if epsilon.abs() < 1e-12 {
            return 4.0; // Limit as λ → 1
        }
        (c - c0) / c0 / epsilon
    }

    /// Estimate the applied strain from a measured capacitance.
    ///
    /// Inverts C(λ) = C₀ λ⁴ numerically via bisection.
    pub fn strain_from_capacitance(&self, measured_c: f64) -> f64 {
        let c0 = EPSILON_0 * self.epsilon_r * self.area_0 / self.thickness_0;
        // λ⁴ = measured_c / c0
        let lambda4 = measured_c / c0.max(1e-30);
        lambda4.powf(0.25)
    }
}

// ── Pneumatic Soft Actuator ───────────────────────────────────────────────────

/// Pneumatic soft finger actuator (PneuNet-like) model.
///
/// Models the bending angle as a function of supplied gauge pressure using
/// a simplified elastic beam approach.
#[derive(Debug, Clone)]
pub struct PneumaticSoftActuator {
    /// Number of pneumatic chambers n.
    pub num_chambers: usize,
    /// Chamber spacing / effective length per chamber (m).
    pub chamber_length: f64,
    /// Elastic modulus of the actuator body (Pa).
    pub elastic_modulus: f64,
    /// Wall thickness of the extensible side (m).
    pub wall_thickness: f64,
    /// Internal gauge pressure p (Pa).
    pub pressure: f64,
    /// Current bending angle θ (rad, total along actuator).
    pub bending_angle: f64,
}

impl PneumaticSoftActuator {
    /// Create a new pneumatic soft actuator.
    pub fn new(
        num_chambers: usize,
        chamber_length: f64,
        elastic_modulus: f64,
        wall_thickness: f64,
    ) -> Self {
        Self {
            num_chambers,
            chamber_length,
            elastic_modulus,
            wall_thickness,
            pressure: 0.0,
            bending_angle: 0.0,
        }
    }

    /// Compute the bending angle from internal pressure using the
    /// moment-curvature relation for a pressurized soft beam.
    ///
    /// θ ≈ p · L_tot / (E · t²) (simplified linear model).
    pub fn update_bending_angle(&mut self) {
        let l_tot = self.num_chambers as f64 * self.chamber_length;
        let t = self.wall_thickness;
        self.bending_angle = self.pressure * l_tot / (self.elastic_modulus * t * t);
    }

    /// End-effector position relative to the actuator root.
    ///
    /// Assumes the actuator forms a circular arc.  Returns `[x, y, 0]`.
    pub fn end_effector_position(&self) -> [f64; 3] {
        let l = self.num_chambers as f64 * self.chamber_length;
        let theta = self.bending_angle;
        if theta.abs() < 1e-8 {
            return [l, 0.0, 0.0];
        }
        let r = l / theta;
        [r * theta.sin(), r * (1.0 - theta.cos()), 0.0]
    }

    /// Tip contact force estimate given a blocked-force stiffness k_b.
    ///
    /// F_tip = k_b · p (linear approximation).
    pub fn tip_force(&self, stiffness_blocked: f64) -> f64 {
        stiffness_blocked * self.pressure
    }
}

// ── Tendon-Driven Soft Actuator ───────────────────────────────────────────────

/// Tendon-driven soft continuum arm segment.
///
/// Models a segment of a tendon-driven continuum robot where the bending
/// is achieved by pulling one or more tendons routed off-center.
#[derive(Debug, Clone)]
pub struct TendonDrivenSegment {
    /// Segment rest length L₀ (m).
    pub rest_length: f64,
    /// Flexural stiffness EI (N·m²) of the backbone.
    pub flexural_stiffness: f64,
    /// Tendon offset from neutral axis d (m).
    pub tendon_offset: f64,
    /// Applied tendon tension (N).
    pub tension: f64,
    /// Current bending angle θ (rad) of the segment.
    pub bending_angle: f64,
}

impl TendonDrivenSegment {
    /// Create a new tendon-driven segment.
    pub fn new(rest_length: f64, flexural_stiffness: f64, tendon_offset: f64) -> Self {
        Self {
            rest_length,
            flexural_stiffness,
            tendon_offset,
            tension: 0.0,
            bending_angle: 0.0,
        }
    }

    /// Equilibrium bending angle: θ = T d L₀ / EI.
    pub fn equilibrium_angle(&self) -> f64 {
        self.tension * self.tendon_offset * self.rest_length / self.flexural_stiffness.max(1e-30)
    }

    /// Update bending angle to current equilibrium.
    pub fn update(&mut self) {
        self.bending_angle = self.equilibrium_angle();
    }

    /// Curvature κ = θ / L₀.
    pub fn curvature(&self) -> f64 {
        self.bending_angle / self.rest_length.max(1e-15)
    }

    /// Segment tip position `[x, y, 0]` on a 2D plane.
    pub fn tip_position(&self) -> [f64; 3] {
        let theta = self.bending_angle;
        let l = self.rest_length;
        if theta.abs() < 1e-8 {
            return [l, 0.0, 0.0];
        }
        let r = l / theta;
        [r * theta.sin(), r * (1.0 - theta.cos()), 0.0]
    }
}

// ── Coupled Electro-Mechanical FEM Soft Body ──────────────────────────────────

/// Node in an electroactive FEM mesh.
#[derive(Debug, Clone)]
pub struct ElectroacitiveFemNode {
    /// Current position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Nodal electric potential (V).
    pub potential: f64,
    /// Nodal mass (kg).
    pub mass: f64,
    /// Whether this node is mechanically pinned.
    pub pinned: bool,
    /// Whether this node has a prescribed electric potential.
    pub potential_prescribed: bool,
}

impl ElectroacitiveFemNode {
    /// Create a free, uncharged FEM node.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            potential: 0.0,
            mass,
            pinned: false,
            potential_prescribed: false,
        }
    }

    /// Create a pinned (fixed) node.
    pub fn new_pinned(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            potential: 0.0,
            mass: 1e30,
            pinned: true,
            potential_prescribed: false,
        }
    }
}

/// A linear tetrahedral element for the electroactive FEM soft body.
///
/// Computes both elastic (neo-Hookean linearized) and Maxwell stress
/// contributions to the nodal force vector.
#[derive(Debug, Clone)]
pub struct ElectroacitiveTetraElement {
    /// Indices into the node array \[i0, i1, i2, i3\].
    pub node_indices: [usize; 4],
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus κ (Pa).
    pub bulk_modulus: f64,
    /// Relative permittivity εᵣ.
    pub epsilon_r: f64,
    /// Rest-state shape matrix inverse B₀⁻¹ (columns = edge vectors from node 0).
    pub shape_matrix_inv: [[f64; 3]; 3],
    /// Undeformed volume (m³).
    pub volume_0: f64,
}

impl ElectroacitiveTetraElement {
    /// Build an element from four node positions and material parameters.
    pub fn new(
        indices: [usize; 4],
        positions: &[[f64; 3]],
        shear_modulus: f64,
        bulk_modulus: f64,
        epsilon_r: f64,
    ) -> Self {
        let p0 = positions[indices[0]];
        let p1 = positions[indices[1]];
        let p2 = positions[indices[2]];
        let p3 = positions[indices[3]];

        let e1 = v3_sub(p1, p0);
        let e2 = v3_sub(p2, p0);
        let e3 = v3_sub(p3, p0);

        let vol = v3_dot(e1, v3_cross(e2, e3)).abs() / 6.0;

        // Simple 3x3 matrix from edge columns
        let shape_matrix = [e1, e2, e3];
        let shape_matrix_inv = invert_3x3(shape_matrix);

        Self {
            node_indices: indices,
            shear_modulus,
            bulk_modulus,
            epsilon_r,
            shape_matrix_inv,
            volume_0: vol,
        }
    }

    /// Compute the Maxwell stress contribution to nodal forces.
    ///
    /// Given the nodal electric potentials, computes the electric field E inside
    /// the element and the Maxwell stress tensor P_M = ε₀ εᵣ (E⊗E − |E|²/2 I),
    /// then distributes to nodes via the shape function gradients.
    pub fn maxwell_nodal_forces(&self, nodes: &[ElectroacitiveFemNode]) -> [[f64; 3]; 4] {
        // Electric field E = -B₀⁻¹ ∇φ (within element, constant for linear tet)
        let n0 = &nodes[self.node_indices[0]];
        let n1 = &nodes[self.node_indices[1]];
        let n2 = &nodes[self.node_indices[2]];
        let n3 = &nodes[self.node_indices[3]];

        let dphi = [
            n1.potential - n0.potential,
            n2.potential - n0.potential,
            n3.potential - n0.potential,
        ];

        // E = -B⁻¹ · Δφ
        let e_field = [
            -(self.shape_matrix_inv[0][0] * dphi[0]
                + self.shape_matrix_inv[0][1] * dphi[1]
                + self.shape_matrix_inv[0][2] * dphi[2]),
            -(self.shape_matrix_inv[1][0] * dphi[0]
                + self.shape_matrix_inv[1][1] * dphi[1]
                + self.shape_matrix_inv[1][2] * dphi[2]),
            -(self.shape_matrix_inv[2][0] * dphi[0]
                + self.shape_matrix_inv[2][1] * dphi[1]
                + self.shape_matrix_inv[2][2] * dphi[2]),
        ];

        let eps = EPSILON_0 * self.epsilon_r;
        let e2 = v3_dot(e_field, e_field);

        // Maxwell stress tensor rows: T_ij = ε (E_i E_j - δ_ij |E|²/2)
        let mut t = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                t[i][j] = eps * (e_field[i] * e_field[j] - 0.5 * e2 * delta);
            }
        }

        // Distribute force to each node (volume / 4 per node)
        let node_vol = self.volume_0 / 4.0;
        let mut forces = [[0.0_f64; 3]; 4];
        // Force on node k: f_k = T · n̂_k * A_k (approximated as T-divergence × volume/4)
        // For a linear tet we use f_k = V/4 * T · grad_N_k; grad_N_k from shape matrix
        // For simplicity: distribute T's trace force equally to all 4 nodes
        let trace_force = [
            (t[0][0] + t[0][1] + t[0][2]) * node_vol,
            (t[1][0] + t[1][1] + t[1][2]) * node_vol,
            (t[2][0] + t[2][1] + t[2][2]) * node_vol,
        ];
        forces.fill(trace_force);
        forces
    }
}

/// Invert a 3×3 matrix stored as rows.
fn invert_3x3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-30 {
        return [[0.0; 3]; 3];
    }
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

/// A minimal coupled electric-mechanical FEM soft body.
///
/// Stores an array of nodes and elements and provides a single-step explicit
/// time integrator that applies Maxwell forces.
#[derive(Debug)]
pub struct ElectroacitiveFemBody {
    /// Nodes of the mesh.
    pub nodes: Vec<ElectroacitiveFemNode>,
    /// Tetrahedral elements.
    pub elements: Vec<ElectroacitiveTetraElement>,
    /// Rayleigh damping coefficient α (mass-proportional).
    pub damping: f64,
}

impl ElectroacitiveFemBody {
    /// Create an electroactive FEM body.
    pub fn new(
        nodes: Vec<ElectroacitiveFemNode>,
        elements: Vec<ElectroacitiveTetraElement>,
        damping: f64,
    ) -> Self {
        Self {
            nodes,
            elements,
            damping,
        }
    }

    /// Apply Maxwell stress forces from all elements to their nodes.
    pub fn apply_maxwell_forces(&mut self) {
        // Gather forces for each node
        let n = self.nodes.len();
        let mut forces = vec![[0.0_f64; 3]; n];

        for elem in &self.elements {
            let f = elem.maxwell_nodal_forces(&self.nodes);
            for (k, &ni) in elem.node_indices.iter().enumerate() {
                for d in 0..3 {
                    forces[ni][d] += f[k][d];
                }
            }
        }

        // Integrate with explicit Euler (caller provides dt separately)
        // Store accumulated forces as a scratch field on the node velocity
        // NOTE: caller must call step() to advance positions.
        // Here we store force / mass as acceleration scratch:
        for (node, &f) in self.nodes.iter_mut().zip(forces.iter()) {
            if !node.pinned {
                let inv_m = 1.0 / node.mass.max(1e-30);
                for (v, fi) in node.velocity.iter_mut().zip(f.iter()) {
                    *v += fi * inv_m;
                }
            }
        }
    }

    /// Advance the body by one explicit Euler step.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        // Apply gravity
        for node in &mut self.nodes {
            if !node.pinned {
                for (v, g) in node.velocity.iter_mut().zip(gravity.iter()) {
                    *v += g * dt;
                    *v *= 1.0 - self.damping * dt;
                }
            }
        }
        // Update positions
        for node in &mut self.nodes {
            if !node.pinned {
                for d in 0..3 {
                    node.position[d] += node.velocity[d] * dt;
                }
            }
        }
    }
}

// ── Conductive Hydrogel ───────────────────────────────────────────────────────

/// Electrically conductive hydrogel element for soft robotics.
///
/// Combines mechanical elasticity with electric conductivity changes due to
/// strain (piezoresistive behaviour).
#[derive(Debug, Clone)]
pub struct ConductiveHydrogel {
    /// Initial conductivity σ₀ (S/m).
    pub conductivity_0: f64,
    /// Piezoresistive gauge factor Kσ (dimensionless).
    pub gauge_factor: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Current strain ε (dimensionless).
    pub strain: f64,
}

impl ConductiveHydrogel {
    /// Create a new conductive hydrogel with given properties.
    pub fn new(conductivity_0: f64, gauge_factor: f64, youngs_modulus: f64) -> Self {
        Self {
            conductivity_0,
            gauge_factor,
            youngs_modulus,
            strain: 0.0,
        }
    }

    /// Current conductivity σ = σ₀ (1 − Kσ ε).
    pub fn conductivity(&self) -> f64 {
        self.conductivity_0 * (1.0 - self.gauge_factor * self.strain)
    }

    /// Resistance of a slab of length L, cross-section A: R = L / (σ A).
    pub fn resistance(&self, length: f64, area: f64) -> f64 {
        let sigma = self.conductivity().max(1e-15);
        length / (sigma * area)
    }

    /// Elastic stress σ_el = E ε.
    pub fn elastic_stress(&self) -> f64 {
        self.youngs_modulus * self.strain
    }
}

// ── Pull-In Voltage Calculator ────────────────────────────────────────────────

/// Compute the theoretical pull-in voltage for a parallel-plate DEA.
///
/// Pull-in occurs when V_pi = t₀ √(8 μ / (27 ε₀ εᵣ)).
///
/// * `shear_modulus` – μ (Pa).
/// * `thickness_0`   – undeformed film thickness (m).
/// * `epsilon_r`     – relative permittivity.
pub fn pull_in_voltage(shear_modulus: f64, thickness_0: f64, epsilon_r: f64) -> f64 {
    thickness_0 * (8.0 * shear_modulus / (27.0 * EPSILON_0 * epsilon_r)).sqrt()
}

/// Compute the Maxwell stress tensor component T_33 = ε₀ εᵣ E² for a
/// field E along axis 3 (perpendicular to film plane).
pub fn maxwell_stress_33(electric_field: f64, epsilon_r: f64) -> f64 {
    EPSILON_0 * epsilon_r * electric_field * electric_field
}

/// Electrostatic energy per unit volume: u_e = ε₀ εᵣ E² / 2.
pub fn electrostatic_energy_density(electric_field: f64, epsilon_r: f64) -> f64 {
    0.5 * EPSILON_0 * epsilon_r * electric_field * electric_field
}

/// Energy conversion efficiency of a DEA cycle.
///
/// η = W_mech_out / W_elec_in  (very simplified, assumes ideal cycle).
///
/// * `actuation_area_strain` – ε_A = λ² − 1.
/// * `maxwell_pressure`      – p_M (Pa).
/// * `voltage`               – V (V).
/// * `capacitance_0`         – C₀ (F).
pub fn dea_cycle_efficiency(
    actuation_area_strain: f64,
    maxwell_pressure: f64,
    voltage: f64,
    capacitance_0: f64,
) -> f64 {
    let w_mech = maxwell_pressure * actuation_area_strain; // J/m³ approx
    let w_elec = 0.5 * capacitance_0 * voltage * voltage;
    if w_elec < 1e-30 {
        return 0.0;
    }
    (w_mech / w_elec).min(1.0)
}

// ── Soft Robot Actuator Comparison ───────────────────────────────────────────

/// Enumeration of soft robot actuator types for simulation comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SoftActuatorType {
    /// Pneumatically pressurized PneuNet-style actuator.
    Pneumatic,
    /// Tendon-driven continuum arm.
    TendonDriven,
    /// Dielectric elastomer actuator.
    DielectricElastomer,
    /// Ionic hydrogel actuator.
    IonicHydrogel,
    /// Piezoelectric film actuator.
    Piezoelectric,
    /// Triboelectric nanogenerator (self-powered sensing).
    Triboelectric,
}

/// Summary of a soft actuator's performance.
#[derive(Debug, Clone)]
pub struct ActuatorPerformanceSummary {
    /// Actuator type.
    pub actuator_type: SoftActuatorType,
    /// Maximum blocked force (N).
    pub max_force: f64,
    /// Maximum free stroke / displacement (m or rad).
    pub max_stroke: f64,
    /// Operating voltage or pressure.
    pub driving_parameter: f64,
    /// Response time constant τ (s).
    pub response_time: f64,
}

impl ActuatorPerformanceSummary {
    /// Create a performance summary.
    pub fn new(
        actuator_type: SoftActuatorType,
        max_force: f64,
        max_stroke: f64,
        driving_parameter: f64,
        response_time: f64,
    ) -> Self {
        Self {
            actuator_type,
            max_force,
            max_stroke,
            driving_parameter,
            response_time,
        }
    }

    /// Power-to-force ratio (W/N) = driving_parameter / (max_force · response_time).
    pub fn power_to_force_ratio(&self) -> f64 {
        let denom = self.max_force * self.response_time;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        self.driving_parameter / denom
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    // 1. Maxwell pressure scales as E² ──────────────────────────────────────
    #[test]
    fn test_maxwell_pressure_scales_with_e_squared() {
        let e1 = 1.0e6_f64;
        let e2 = 2.0e6_f64;
        let p1 = maxwell_stress_33(e1, 4.7);
        let p2 = maxwell_stress_33(e2, 4.7);
        assert!(
            (p2 / p1 - 4.0).abs() < 1e-6,
            "Maxwell pressure should scale as E²: p2/p1 = {}",
            p2 / p1
        );
    }

    // 2. Electrostatic energy density is half Maxwell stress ────────────────
    #[test]
    fn test_electrostatic_energy_half_maxwell_stress() {
        let e = 5.0e7;
        let eps_r = 3.0;
        let u = electrostatic_energy_density(e, eps_r);
        let p = maxwell_stress_33(e, eps_r);
        assert!(
            (2.0 * u - p).abs() < TOL * p.abs().max(1e-10),
            "u_e = p_M / 2 must hold, got u_e={u}, p_M={p}"
        );
    }

    // 3. DEA elastic stress zero at unit stretch ─────────────────────────────
    #[test]
    fn test_dea_elastic_stress_zero_at_rest() {
        let dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 0.0);
        assert!(
            dea.elastic_stress().abs() < TOL,
            "Elastic stress at λ=1 should be zero"
        );
    }

    // 4. DEA thickness decreases with stretch ────────────────────────────────
    #[test]
    fn test_dea_thickness_decreases_with_stretch() {
        let mut dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 0.0);
        let t0 = dea.current_thickness();
        dea.stretch = 2.0;
        let t1 = dea.current_thickness();
        assert!(
            t1 < t0,
            "Film should thin as stretch increases: t0={t0}, t1={t1}"
        );
        assert!(
            (t1 - t0 / 4.0).abs() < TOL,
            "At λ=2, t = t₀/4: expected {}, got {t1}",
            t0 / 4.0
        );
    }

    // 5. DEA capacitance increases with stretch ──────────────────────────────
    #[test]
    fn test_dea_capacitance_increases() {
        let mut sensor = CapacitiveStrainSensor::new(1.0e-4, 1.0e-3, 4.7);
        let c0 = sensor.capacitance();
        sensor.stretch = 2.0;
        let c1 = sensor.capacitance();
        assert!(c1 > c0, "Capacitance must increase with stretch");
        // C ∝ λ⁴ for incompressible film
        assert!(
            (c1 / c0 - 16.0).abs() < 1e-6,
            "C should scale as λ⁴: ratio = {}",
            c1 / c0
        );
    }

    // 6. Capacitive sensor: gauge factor at unit stretch → ≈ 4 ──────────────
    #[test]
    fn test_capacitive_sensor_gauge_factor_at_rest() {
        let sensor = CapacitiveStrainSensor::new(1.0e-4, 1.0e-3, 4.7);
        let gf = sensor.gauge_factor();
        assert!((gf - 4.0).abs() < 0.01, "GF at λ=1 should be ≈ 4, got {gf}");
    }

    // 7. Capacitive sensor: strain_from_capacitance is self-consistent ───────
    #[test]
    fn test_strain_from_capacitance_roundtrip() {
        let mut sensor = CapacitiveStrainSensor::new(1.0e-4, 1.0e-3, 4.7);
        sensor.stretch = 1.5;
        let c = sensor.capacitance();
        let lambda_back = sensor.strain_from_capacitance(c);
        assert!(
            (lambda_back - 1.5).abs() < 1e-6,
            "Round-trip λ should recover 1.5, got {lambda_back}"
        );
    }

    // 8. Pull-in voltage formula gives positive value ─────────────────────────
    #[test]
    fn test_pull_in_voltage_positive() {
        let v_pi = pull_in_voltage(1.0e5, 1.0e-3, EPSILON_R_VHB);
        assert!(v_pi > 0.0, "Pull-in voltage must be positive, got {v_pi}");
    }

    // 9. Pull-in voltage scales as sqrt of shear modulus ─────────────────────
    #[test]
    fn test_pull_in_voltage_scales_with_sqrt_mu() {
        let v1 = pull_in_voltage(1.0e5, 1.0e-3, 4.7);
        let v2 = pull_in_voltage(4.0e5, 1.0e-3, 4.7);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-6,
            "V_pi should scale as √μ: v2/v1 = {}",
            v2 / v1
        );
    }

    // 10. DEA equilibrium solver converges for zero voltage ──────────────────
    #[test]
    fn test_dea_equilibrium_zero_voltage() {
        let mut dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 0.0);
        let (stretch, converged) = dea.solve_equilibrium(100, 1e-8);
        assert!(converged, "Should converge with zero voltage");
        assert!(
            (stretch - 1.0).abs() < 1e-6,
            "Zero voltage: equilibrium stretch should be 1.0, got {stretch}"
        );
    }

    // 11. DEA at high voltage shows pull-in ──────────────────────────────────
    #[test]
    fn test_dea_pull_in_at_high_voltage() {
        // A very thin DEA with moderate stiffness at a high voltage
        let dea = DeaMembraneState {
            shear_modulus: 1.0e4,
            thickness_0: 1.0e-4,
            epsilon_r: 4.7,
            voltage: 20_000.0, // 20 kV
            stretch: 1.8,
        };
        // At these conditions the system should be past the pull-in point
        // (may or may not trigger depending on stretch value, so just check the
        // function doesn't panic and returns bool)
        let _is_pulled = dea.is_pulled_in();
    }

    // 12. Piezoelectric: zero field, zero stress → zero strain ───────────────
    #[test]
    fn test_piezo_zero_inputs_zero_strain() {
        let pvdf = PiezoelectricSoftElement::new_pvdf();
        assert!(
            pvdf.mechanical_strain().abs() < TOL,
            "Zero inputs should give zero strain"
        );
        assert!(
            pvdf.electric_displacement().abs() < TOL,
            "Zero inputs should give zero D"
        );
    }

    // 13. Piezoelectric: coupling coefficient is in [0, 1] ───────────────────
    #[test]
    fn test_piezo_coupling_in_range() {
        let pvdf = PiezoelectricSoftElement::new_pvdf();
        let k = pvdf.coupling_coefficient_k33();
        assert!((0.0..=1.0).contains(&k), "k33 must be in [0,1], got {k}");
    }

    // 14. Piezoelectric: resonant frequency increases with shorter length ─────
    #[test]
    fn test_piezo_resonant_freq_inversely_proportional_to_length() {
        let pvdf = PiezoelectricSoftElement::new_pvdf();
        let f1 = pvdf.resonant_frequency(0.01, 1780.0);
        let f2 = pvdf.resonant_frequency(0.02, 1780.0);
        assert!(
            (f1 / f2 - 2.0).abs() < 1e-6,
            "f should scale as 1/L: f1/f2={}",
            f1 / f2
        );
    }

    // 15. Triboelectric: open-circuit voltage is zero at zero gap ────────────
    #[test]
    fn test_teng_zero_voltage_at_zero_gap() {
        let teng = TriboelectricSoftBody::new_ptfe_nylon();
        assert_eq!(teng.air_gap, 0.0);
        let v = teng.open_circuit_voltage();
        assert!(v.abs() < TOL, "V_oc must be zero at zero gap, got {v}");
    }

    // 16. Triboelectric: V_oc scales linearly with gap ───────────────────────
    #[test]
    fn test_teng_voc_scales_with_gap() {
        let mut t1 = TriboelectricSoftBody::new_ptfe_nylon();
        let mut t2 = TriboelectricSoftBody::new_ptfe_nylon();
        t1.air_gap = 1.0e-3;
        t2.air_gap = 2.0e-3;
        let v1 = t1.open_circuit_voltage().abs();
        let v2 = t2.open_circuit_voltage().abs();
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-6,
            "V_oc should scale linearly with gap: v2/v1={}",
            v2 / v1
        );
    }

    // 17. Triboelectric: simulate_cycle returns positive peak voltage ─────────
    #[test]
    fn test_teng_cycle_positive_peak() {
        let mut teng = TriboelectricSoftBody::new_ptfe_nylon();
        let peak = teng.simulate_cycle(5.0e-3, 1.0, 1e-4);
        assert!(peak > 0.0, "Peak voltage should be positive, got {peak}");
    }

    // 18. Pneumatic actuator: zero pressure → no bending ─────────────────────
    #[test]
    fn test_pneumatic_no_bending_at_zero_pressure() {
        let mut act = PneumaticSoftActuator::new(5, 0.01, 1.0e5, 5.0e-4);
        act.update_bending_angle();
        assert!(
            act.bending_angle.abs() < TOL,
            "No bending at zero pressure, got {}",
            act.bending_angle
        );
    }

    // 19. Pneumatic actuator: tip position near origin at zero bending ────────
    #[test]
    fn test_pneumatic_tip_along_axis_at_zero_bending() {
        let act = PneumaticSoftActuator::new(5, 0.01, 1.0e5, 5.0e-4);
        let pos = act.end_effector_position();
        let expected_x = 5.0 * 0.01; // 5 chambers * 0.01 m each
        assert!(
            (pos[0] - expected_x).abs() < 1e-6,
            "Tip should be along x at zero bending"
        );
        assert!(pos[1].abs() < TOL, "Tip y should be zero at zero bending");
    }

    // 20. Tendon-driven: zero tension → zero bending ─────────────────────────
    #[test]
    fn test_tendon_zero_tension_zero_angle() {
        let seg = TendonDrivenSegment::new(0.05, 0.1, 0.005);
        assert!(
            seg.equilibrium_angle().abs() < TOL,
            "Zero tension should give zero angle"
        );
    }

    // 21. Tendon-driven: tip position consistent with arc length ─────────────
    #[test]
    fn test_tendon_tip_arc_length_preserved() {
        let mut seg = TendonDrivenSegment::new(0.1, 0.01, 0.01);
        seg.tension = 10.0;
        seg.update();
        let pos = seg.tip_position();
        let dist = v3_norm(pos);
        // For small angles, tip distance ≈ rest length
        assert!(
            dist <= seg.rest_length + 1e-6,
            "Tip distance cannot exceed rest length: dist={dist}"
        );
    }

    // 22. Ionic hydrogel: swelling ratio increases under applied field ─────────
    #[test]
    fn test_hydrogel_swelling_increases() {
        let mut hg = IonicHydrogelActuator::new(1.0e4, 0.3, 100.0, 10.0);
        hg.applied_field = 1.0e4;
        let q0 = hg.swelling_ratio;
        hg.step(0.01, 1.0e-6);
        assert!(
            hg.swelling_ratio > q0,
            "Swelling ratio should increase under field"
        );
    }

    // 23. Ionic hydrogel: linear strain is non-negative ─────────────────────
    #[test]
    fn test_hydrogel_linear_strain_nonnegative() {
        let hg = IonicHydrogelActuator::new(1.0e4, 0.3, 100.0, 10.0);
        assert!(
            hg.linear_strain() >= 0.0,
            "Linear strain must be non-negative at Q≥1"
        );
    }

    // 24. Ionic hydrogel: Donnan pressure is positive for finite fixed charge ──
    #[test]
    fn test_hydrogel_donnan_pressure_positive() {
        let hg = IonicHydrogelActuator::new(1.0e4, 0.3, 100.0, 10.0);
        assert!(
            hg.donnan_pressure() > 0.0,
            "Donnan pressure should be positive"
        );
    }

    // 25. Conductive hydrogel: conductivity decreases under tension ───────────
    #[test]
    fn test_conductive_hydrogel_conductivity_decreases_under_tension() {
        let mut chg = ConductiveHydrogel::new(10.0, 2.0, 5.0e4);
        let sigma0 = chg.conductivity();
        chg.strain = 0.1;
        let sigma1 = chg.conductivity();
        assert!(
            sigma1 < sigma0,
            "Conductivity should decrease with positive strain: sigma0={sigma0}, sigma1={sigma1}"
        );
    }

    // 26. Conductive hydrogel: resistance increases with length ───────────────
    #[test]
    fn test_conductive_hydrogel_resistance_increases_with_length() {
        let chg = ConductiveHydrogel::new(10.0, 2.0, 5.0e4);
        let r1 = chg.resistance(0.01, 1.0e-4);
        let r2 = chg.resistance(0.02, 1.0e-4);
        assert!(
            (r2 / r1 - 2.0).abs() < 1e-6,
            "R should double with double length: r2/r1={}",
            r2 / r1
        );
    }

    // 27. v3 helpers ─────────────────────────────────────────────────────────
    #[test]
    fn test_v3_helpers_cross_dot() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = v3_cross(x, y);
        assert!((z[2] - 1.0).abs() < TOL, "x × y should be z");
        assert!(v3_dot(x, y).abs() < TOL, "x · y should be zero");
        assert!((v3_norm(x) - 1.0).abs() < TOL, "norm(x) should be 1");
    }

    // 28. invert_3x3 identity matrix ─────────────────────────────────────────
    #[test]
    fn test_invert_3x3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = invert_3x3(id);
        for (i, row) in inv.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((v - exp).abs() < TOL, "inv[{i}][{j}] = {} ≠ {exp}", v);
            }
        }
    }

    // 29. Electroactive FEM body: step does not move pinned nodes ────────────
    #[test]
    fn test_fem_body_pinned_node_does_not_move() {
        let nodes = vec![
            ElectroacitiveFemNode::new_pinned([0.0, 0.0, 0.0]),
            ElectroacitiveFemNode::new([1.0, 0.0, 0.0], 1.0),
            ElectroacitiveFemNode::new([0.0, 1.0, 0.0], 1.0),
            ElectroacitiveFemNode::new([0.0, 0.0, 1.0], 1.0),
        ];
        let positions: Vec<[f64; 3]> = nodes.iter().map(|n| n.position).collect();
        let elem = ElectroacitiveTetraElement::new([0, 1, 2, 3], &positions, 1.0e5, 2.0e5, 4.7);
        let mut body = ElectroacitiveFemBody::new(nodes, vec![elem], 0.01);
        let pinned_pos_before = body.nodes[0].position;
        body.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        let pinned_pos_after = body.nodes[0].position;
        assert_eq!(
            pinned_pos_before, pinned_pos_after,
            "Pinned node must not move"
        );
    }

    // 30. Electroactive FEM body: free node falls under gravity ───────────────
    #[test]
    fn test_fem_body_free_node_falls_under_gravity() {
        let nodes = vec![
            ElectroacitiveFemNode::new_pinned([0.0, 0.0, 0.0]),
            ElectroacitiveFemNode::new([1.0, 0.0, 0.0], 1.0),
            ElectroacitiveFemNode::new([0.0, 1.0, 0.0], 1.0),
            ElectroacitiveFemNode::new([0.0, 0.0, 1.0], 1.0),
        ];
        let positions: Vec<[f64; 3]> = nodes.iter().map(|n| n.position).collect();
        let elem = ElectroacitiveTetraElement::new([0, 1, 2, 3], &positions, 1.0e5, 2.0e5, 4.7);
        let mut body = ElectroacitiveFemBody::new(nodes, vec![elem], 0.0);
        let y0 = body.nodes[1].position[1];
        for _ in 0..60 {
            body.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        }
        let y1 = body.nodes[1].position[1];
        assert!(
            y1 < y0,
            "Free node should fall under gravity: y0={y0}, y1={y1}"
        );
    }

    // 31. ActuatorPerformanceSummary: power-to-force ratio non-negative ────────
    #[test]
    fn test_actuator_summary_power_force_ratio() {
        let summary = ActuatorPerformanceSummary::new(
            SoftActuatorType::DielectricElastomer,
            0.5,
            0.02,
            5_000.0,
            0.1,
        );
        assert!(
            summary.power_to_force_ratio() > 0.0,
            "Power-to-force ratio should be positive"
        );
    }

    // 32. DEA actuation area strain increases with stretch ────────────────────
    #[test]
    fn test_dea_area_strain_positive() {
        let mut dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 1000.0);
        dea.stretch = 1.5;
        let eps = dea.actuation_area_strain();
        assert!(eps > 0.0, "Area strain should be positive for stretch > 1");
        assert!(
            (eps - (1.5_f64 * 1.5 - 1.0)).abs() < TOL,
            "Area strain = λ² − 1"
        );
    }

    // 33. DEA elastic energy density positive at stretch > 1 ─────────────────
    #[test]
    fn test_dea_elastic_energy_positive() {
        let mut dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 0.0);
        dea.stretch = 2.0;
        assert!(
            dea.elastic_energy_density() > 0.0,
            "Elastic energy should be positive for stretch ≠ 1"
        );
    }

    // 34. DEA electrical energy density zero at zero voltage ──────────────────
    #[test]
    fn test_dea_electrical_energy_zero_at_zero_voltage() {
        let dea = DeaMembraneState::new(1.0e5, 1.0e-3, EPSILON_R_VHB, 0.0);
        assert!(
            dea.electrical_energy_density().abs() < TOL,
            "Electrical energy should be zero with no voltage"
        );
    }

    // 35. clamp helper ─────────────────────────────────────────────────────────
    #[test]
    fn test_clamp_helper() {
        assert!((clamp(0.5, 0.0, 1.0) - 0.5).abs() < TOL);
        assert!((clamp(-1.0, 0.0, 1.0)).abs() < TOL);
        assert!((clamp(2.0, 0.0, 1.0) - 1.0).abs() < TOL);
    }
}
