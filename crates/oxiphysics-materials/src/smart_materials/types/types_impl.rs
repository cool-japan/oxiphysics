//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types_3::{EapType, ShapeMemoryAlloy};

/// Piezoelectric stack actuator model.
#[derive(Debug, Clone, Copy)]
pub struct PiezoActuator {
    /// Free stroke per unit voltage \[m/V\].
    pub stroke_per_volt: f64,
    /// Blocking force per unit voltage \[N/V\].
    pub force_per_volt: f64,
    /// Stiffness \[N/m\].
    pub stiffness: f64,
    /// Resonance frequency \[Hz\].
    pub resonance_freq: f64,
    /// Capacitance \[F\].
    pub capacitance: f64,
    /// Hysteresis coefficient (0 = linear, >0 = hysteresis).
    pub hysteresis: f64,
}
impl PiezoActuator {
    /// Create a piezoelectric stack actuator.
    pub fn new(
        stroke_per_volt: f64,
        force_per_volt: f64,
        stiffness: f64,
        resonance_freq: f64,
        capacitance: f64,
    ) -> Self {
        Self {
            stroke_per_volt,
            force_per_volt,
            stiffness,
            resonance_freq,
            capacitance,
            hysteresis: 0.1,
        }
    }
    /// Free displacement at voltage V \[m\].
    pub fn free_stroke(&self, voltage: f64) -> f64 {
        self.stroke_per_volt * voltage
    }
    /// Blocking force at voltage V \[N\].
    pub fn blocking_force(&self, voltage: f64) -> f64 {
        self.force_per_volt * voltage
    }
    /// Output force at voltage V with external load displacement δ \[N\].
    pub fn output_force(&self, voltage: f64, displacement: f64) -> f64 {
        let free = self.free_stroke(voltage);
        self.stiffness * (free - displacement)
    }
    /// Electrical energy input per cycle \[J\].
    pub fn energy_per_cycle(&self, voltage: f64, freq: f64) -> f64 {
        let _ = freq;
        0.5 * self.capacitance * voltage * voltage
    }
    /// Mechanical work output \[J\].
    pub fn mechanical_work(&self, voltage: f64) -> f64 {
        let f = self.blocking_force(voltage);
        let d = self.free_stroke(voltage);
        0.5 * f * d
    }
}
/// Smart piezoelectric actuator/sensor model.
///
/// Implements the full piezoelectric constitutive equations:
/// - Strain from electric field (converse effect): S = d · E
/// - Charge from stress (direct effect): D = d · T
/// - Resonance frequency from elastic compliance and geometry
/// - Actuation stroke and coupling coefficient k
#[derive(Debug, Clone, Copy)]
pub struct PiezoelectricSmart {
    /// d_33 piezoelectric charge coefficient \[C/N = m/V\].
    pub d33: f64,
    /// d_31 piezoelectric charge coefficient \[C/N = m/V\].
    pub d31: f64,
    /// d_15 (shear) piezoelectric charge coefficient \[C/N = m/V\].
    pub d15: f64,
    /// Elastic compliance s_33 at constant E field \[m²/N\].
    pub s33_e: f64,
    /// Elastic compliance s_11 at constant E field \[m²/N\].
    pub s11_e: f64,
    /// Permittivity ε_33 at constant stress \[F/m\].
    pub eps33_t: f64,
    /// Material density \[kg/m³\].
    pub density: f64,
    /// Geometry: length along poling axis \[m\].
    pub length: f64,
    /// Cross-sectional area \[m²\].
    pub area: f64,
}
impl PiezoelectricSmart {
    /// Create a piezoelectric smart material model.
    pub fn new(
        d33: f64,
        d31: f64,
        d15: f64,
        s33_e: f64,
        s11_e: f64,
        eps33_t: f64,
        density: f64,
        length: f64,
        area: f64,
    ) -> Self {
        Self {
            d33,
            d31,
            d15,
            s33_e,
            s11_e,
            eps33_t,
            density,
            length,
            area,
        }
    }
    /// Standard PZT-5H model.
    pub fn pzt5h() -> Self {
        const EPS0: f64 = 8.854e-12;
        Self::new(
            593e-12,
            -274e-12,
            741e-12,
            20.7e-12,
            16.5e-12,
            3400.0 * EPS0,
            7500.0,
            10e-3,
            1e-4,
        )
    }
    /// Standard PZT-4 model.
    pub fn pzt4() -> Self {
        const EPS0: f64 = 8.854e-12;
        Self::new(
            289e-12,
            -123e-12,
            496e-12,
            15.5e-12,
            12.3e-12,
            1300.0 * EPS0,
            7600.0,
            10e-3,
            1e-4,
        )
    }
    /// Electromechanical coupling coefficient k_33.
    ///
    /// k²_33 = d²_33 / (s³³^E · ε³³^T)
    pub fn k33(&self) -> f64 {
        let k2 = self.d33 * self.d33 / (self.s33_e * self.eps33_t);
        k2.sqrt().clamp(0.0, 1.0)
    }
    /// Planar coupling coefficient k_p (thin disc resonator).
    pub fn kp(&self) -> f64 {
        let k31_2 = self.d31 * self.d31 / (self.s11_e * self.eps33_t);
        let kp2 = 2.0 * k31_2 / 0.7_f64;
        kp2.sqrt().clamp(0.0, 1.0)
    }
    /// Resonance frequency of longitudinal mode \[Hz\].
    ///
    /// f_r = 1 / (2·L) · √(1 / (ρ · s³³^E))
    pub fn resonance_frequency(&self) -> f64 {
        let v_sound = (1.0 / (self.density * self.s33_e)).sqrt();
        v_sound / (2.0 * self.length)
    }
    /// Anti-resonance frequency f_a from resonance via k33.
    pub fn antiresonance_frequency(&self) -> f64 {
        let fr = self.resonance_frequency();
        let k = self.k33();
        fr * (1.0 / (1.0 - k * k)).sqrt()
    }
    /// Free stroke (actuation displacement) Δl \[m\] for voltage V.
    ///
    /// Δl = d_33 · V
    pub fn actuation_stroke(&self, voltage: f64) -> f64 {
        self.d33 * voltage
    }
    /// Generalized actuation stroke in 31 mode (transverse) \[m\].
    pub fn actuation_stroke_31(&self, voltage: f64) -> f64 {
        let thickness = self.length;
        let width = self.area.sqrt();
        self.d31 * voltage / thickness * width
    }
    /// Blocking force \[N\] for 33-mode actuator.
    ///
    /// F_block = d_33 · E · A / s³³^E = d_33 · V · A / (s³³^E · L)
    pub fn blocking_force(&self, voltage: f64) -> f64 {
        self.d33 * voltage * self.area / (self.s33_e * self.length)
    }
    /// Charge generated by applied stress (direct piezoelectric effect) \[C\].
    ///
    /// Q = d_33 · F (force along poling axis)
    pub fn charge_from_force(&self, force: f64) -> f64 {
        self.d33 * force
    }
    /// Open-circuit voltage from applied force \[V\].
    pub fn voltage_from_force(&self, force: f64) -> f64 {
        let stress = force / self.area;
        self.d33 * stress * self.length / self.eps33_t
    }
    /// Mechanical quality factor Q_m estimate from bandwidth.
    ///
    /// Q_m = f_r / (f_a - f_r) (approximate).
    pub fn mechanical_q(&self) -> f64 {
        let fr = self.resonance_frequency();
        let fa = self.antiresonance_frequency();
        if (fa - fr).abs() < 1e-3 {
            1000.0
        } else {
            fr / (fa - fr)
        }
    }
}
/// Phase of a shape memory alloy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SmaPhase {
    /// Fully martensitic (low temperature, high martensite fraction).
    Martensite,
    /// Fully austenitic (high temperature, low martensite fraction).
    Austenite,
    /// Mixed phase transformation region.
    Mixed,
}
/// Hydrogel volume change from solvent absorption using Flory-Rehner theory.
///
/// Models equilibrium swelling: chemical potential of mixing + elastic penalty = 0.
#[derive(Debug, Clone, Copy)]
pub struct HydrogelSwelling {
    /// Flory-Huggins interaction parameter χ.
    pub chi: f64,
    /// Polymer volume fraction at synthesis (reference state) φ_0.
    pub phi0: f64,
    /// Number of monomers per network strand N_c (crosslink-related).
    pub n_c: f64,
    /// Molar volume of solvent \[m³/mol\].
    pub v_s: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}
impl HydrogelSwelling {
    /// Create a hydrogel swelling model.
    pub fn new(chi: f64, phi0: f64, n_c: f64, v_s: f64, temperature: f64) -> Self {
        Self {
            chi,
            phi0,
            n_c,
            v_s,
            temperature,
        }
    }
    /// Standard poly(N-isopropylacrylamide) (PNIPAM) hydrogel model.
    pub fn pnipam() -> Self {
        Self::new(0.45, 0.05, 100.0, 18e-6, 298.0)
    }
    /// Flory-Huggins mixing free energy density \[J/m³\].
    ///
    /// f_mix = (RT/V_s) * \[φ ln φ + (1-φ) ln(1-φ) + χ φ (1-φ)\]
    pub fn mixing_free_energy(&self, phi: f64) -> f64 {
        const R: f64 = 8.314;
        let phi = phi.clamp(1e-6, 1.0 - 1e-6);
        let psi = 1.0 - phi;
        (R * self.temperature / self.v_s) * (phi * phi.ln() + psi * psi.ln() + self.chi * phi * psi)
    }
    /// Equilibrium polymer volume fraction φ_eq.
    ///
    /// Solved iteratively: chemical potential = 0.
    pub fn equilibrium_phi(&self) -> f64 {
        let phi_guess = self.phi0;
        let mut phi = phi_guess.clamp(0.01, 0.99);
        for _ in 0..50 {
            let psi = 1.0 - phi;
            let mu_mix = phi.ln() + psi + self.chi * psi * psi;
            let mu_el = (phi / (2.0 * self.n_c)) * (phi / self.phi0).powf(1.0 / 3.0);
            let mu = mu_mix + mu_el;
            let d_mu_mix = 1.0 / phi - 1.0 - 2.0 * self.chi * psi;
            let d_mu_el = (1.0 / (2.0 * self.n_c))
                * (phi / self.phi0).powf(1.0 / 3.0)
                * (1.0 + phi / (3.0 * self.phi0));
            let d_mu = d_mu_mix + d_mu_el;
            if d_mu.abs() < 1e-15 {
                break;
            }
            phi -= mu / d_mu;
            phi = phi.clamp(0.001, 0.999);
        }
        phi
    }
    /// Equilibrium swelling ratio Q = V_swollen / V_dry = 1/φ_eq.
    pub fn swelling_ratio(&self) -> f64 {
        1.0 / self.equilibrium_phi().max(1e-6)
    }
    /// Linear swelling strain ε = Q^(1/3) - 1.
    pub fn linear_strain(&self) -> f64 {
        self.swelling_ratio().powf(1.0 / 3.0) - 1.0
    }
    /// Osmotic pressure \[Pa\] at volume fraction φ.
    pub fn osmotic_pressure(&self, phi: f64) -> f64 {
        const R: f64 = 8.314;
        let phi = phi.clamp(1e-6, 1.0 - 1e-6);
        let psi = 1.0 - phi;
        -(R * self.temperature / self.v_s) * (psi + phi.ln() + self.chi * phi * phi)
    }
}
/// Electroactive polymer (EAP) model.
///
/// Covers both ionic (low voltage, large bending) and electronic (high voltage,
/// large area strain) EAP types.
#[derive(Debug, Clone)]
pub struct ElectroactivePoly {
    /// Type of EAP.
    pub eap_type: EapType,
    /// Maximum actuation strain (dimensionless).
    pub max_strain: f64,
    /// Half-wave voltage V_50 at which strain = 50 % max \[V\].
    pub v_half: f64,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
    /// Bandwidth (time constant) \[s\].
    pub time_constant: f64,
    /// Current actuation state (0..1).
    pub state: f64,
}
impl ElectroactivePoly {
    /// Create an EAP model.
    pub fn new(
        eap_type: EapType,
        max_strain: f64,
        v_half: f64,
        elastic_modulus: f64,
        time_constant: f64,
    ) -> Self {
        Self {
            eap_type,
            max_strain,
            v_half,
            elastic_modulus,
            time_constant,
            state: 0.0,
        }
    }
    /// Standard IPMC model (ionic).
    pub fn ipmc() -> Self {
        Self::new(EapType::Ionic, 0.03, 1.5, 1e6, 0.1)
    }
    /// Standard dielectric elastomer model (electronic).
    pub fn dielectric_elastomer() -> Self {
        Self::new(EapType::Electronic, 0.45, 1000.0, 1e5, 0.001)
    }
    /// Steady-state strain from applied voltage (sigmoid model).
    pub fn steady_state_strain(&self, voltage: f64) -> f64 {
        let k = 4.0 / self.v_half;
        let s = 1.0 / (1.0 + (-k * (voltage - self.v_half)).exp());
        self.max_strain * s
    }
    /// Transient strain at time t after step voltage.
    pub fn transient_strain(&self, voltage: f64, t: f64) -> f64 {
        let ss = self.steady_state_strain(voltage);
        ss * (1.0 - (-t / self.time_constant).exp())
    }
    /// Blocking stress (force per area) \[Pa\].
    pub fn blocking_stress(&self, voltage: f64) -> f64 {
        self.elastic_modulus * self.steady_state_strain(voltage)
    }
    /// Update state (first-order lag) by dt.
    pub fn update(&mut self, voltage: f64, dt: f64) {
        let target = self.steady_state_strain(voltage) / self.max_strain.max(1e-12);
        let tau = self.time_constant;
        self.state += (target - self.state) * (1.0 - (-dt / tau).exp());
        self.state = self.state.clamp(0.0, 1.0);
    }
    /// Current strain from state.
    pub fn current_strain(&self) -> f64 {
        self.state * self.max_strain
    }
}
/// Piezoelectric material coupling: voltage → strain (converse effect)
/// and stress → charge (direct effect).
///
/// Uses the d33 (longitudinal) and d31 (transverse) piezo coefficients.
#[derive(Debug, Clone, Copy)]
pub struct PiezoelectricCoupling {
    /// d33 piezoelectric coefficient \[C/N = m/V\].
    pub d33: f64,
    /// d31 piezoelectric coefficient \[C/N = m/V\] (typically negative).
    pub d31: f64,
    /// Elastic modulus along poling axis \[Pa\].
    pub elastic_modulus_33: f64,
    /// Elastic modulus transverse \[Pa\].
    pub elastic_modulus_11: f64,
    /// Relative permittivity ε_r along poling axis (free).
    pub epsilon_r33: f64,
    /// Coupling factor k33 (electromechanical).
    pub k33: f64,
}
impl PiezoelectricCoupling {
    /// Create a piezoelectric coupling model.
    pub fn new(
        d33: f64,
        d31: f64,
        elastic_modulus_33: f64,
        elastic_modulus_11: f64,
        epsilon_r33: f64,
    ) -> Self {
        const EPS0: f64 = 8.854e-12;
        let k33 = d33 * (elastic_modulus_33 / (EPS0 * epsilon_r33)).sqrt();
        Self {
            d33,
            d31,
            elastic_modulus_33,
            elastic_modulus_11,
            epsilon_r33,
            k33,
        }
    }
    /// Standard PZT-5A model.
    pub fn pzt5a() -> Self {
        Self::new(374e-12, -171e-12, 61e9, 61e9, 1700.0)
    }
    /// Strain along poling axis (33 direction) from applied voltage V over
    /// thickness t \[m\].  ε33 = d33 * E_field = d33 * V / t.
    pub fn strain_from_voltage_33(&self, voltage: f64, thickness: f64) -> f64 {
        self.d33 * voltage / thickness
    }
    /// Transverse strain (31 direction) from applied voltage.
    pub fn strain_from_voltage_31(&self, voltage: f64, thickness: f64) -> f64 {
        self.d31 * voltage / thickness
    }
    /// Charge density (polarization) from applied stress σ33 \[C/m²\].
    pub fn charge_from_stress_33(&self, stress: f64) -> f64 {
        self.d33 * stress
    }
    /// Open-circuit voltage developed from applied stress σ33 \[V\].
    pub fn voltage_from_stress(&self, stress: f64, thickness: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        let p = self.charge_from_stress_33(stress);
        p * thickness / (EPS0 * self.epsilon_r33)
    }
    /// Mechanical energy harvested from stress cycle \[J/m³\].
    pub fn harvested_energy_density(&self, stress_amplitude: f64) -> f64 {
        0.5 * self.d33 * self.d33 * stress_amplitude * stress_amplitude
            / (self.epsilon_r33 * 8.854e-12)
    }
}
/// Dielectric elastomer actuator using Maxwell stress.
#[derive(Debug, Clone, Copy)]
pub struct DielectricElastomer {
    /// Relative permittivity of the elastomer.
    pub epsilon_r: f64,
    /// Undeformed thickness \[m\].
    pub thickness: f64,
    /// Initial area \[m²\].
    pub initial_area: f64,
    /// Elastic modulus of elastomer \[Pa\].
    pub elastic_modulus: f64,
    /// Maximum electric field (breakdown) \[V/m\].
    pub breakdown_field: f64,
}
impl DielectricElastomer {
    /// Create a dielectric elastomer actuator.
    pub fn new(
        epsilon_r: f64,
        thickness: f64,
        initial_area: f64,
        elastic_modulus: f64,
        breakdown_field: f64,
    ) -> Self {
        Self {
            epsilon_r,
            thickness,
            initial_area,
            elastic_modulus,
            breakdown_field,
        }
    }
    /// Maxwell stress (electrostatic pressure) at field E \[Pa\].
    pub fn maxwell_stress(&self, e_field: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        self.epsilon_r * EPS0 * e_field * e_field
    }
    /// Actuation strain (thickness reduction fraction) at voltage V \[V\].
    pub fn actuation_strain(&self, voltage: f64) -> f64 {
        let e_field = voltage / self.thickness;
        let e_eff = e_field.min(self.breakdown_field);
        let p = self.maxwell_stress(e_eff);
        (p / self.elastic_modulus).min(0.9)
    }
    /// Area strain at voltage V (incompressible: λ_z²*λ_A = 1).
    pub fn area_strain(&self, voltage: f64) -> f64 {
        2.0 * self.actuation_strain(voltage)
    }
    /// Actuation pressure \[Pa\] at voltage V.
    pub fn actuation_pressure(&self, voltage: f64) -> f64 {
        let e_field = (voltage / self.thickness).min(self.breakdown_field);
        self.maxwell_stress(e_field)
    }
    /// Electric field for given voltage \[V/m\].
    pub fn electric_field(&self, voltage: f64) -> f64 {
        voltage / self.thickness
    }
}
/// Layered smart composite: SMA layer embedded in elastic matrix.
///
/// Computes effective thermo-mechanical properties using rule-of-mixtures
/// (Voigt/Reuss bounds) for a unidirectional SMA fiber composite.
#[derive(Debug, Clone)]
pub struct SmartComposite {
    /// SMA fiber volume fraction (0..1).
    pub fiber_volume_fraction: f64,
    /// SMA constitutive model.
    pub sma: ShapeMemoryAlloy,
    /// Matrix elastic modulus \[Pa\].
    pub matrix_modulus: f64,
    /// Matrix density \[kg/m³\].
    pub matrix_density: f64,
    /// SMA density \[kg/m³\].
    pub sma_density: f64,
    /// Total thickness \[m\].
    pub thickness: f64,
}
impl SmartComposite {
    /// Create a smart composite.
    pub fn new(
        fiber_volume_fraction: f64,
        sma: ShapeMemoryAlloy,
        matrix_modulus: f64,
        matrix_density: f64,
        sma_density: f64,
        thickness: f64,
    ) -> Self {
        Self {
            fiber_volume_fraction,
            sma,
            matrix_modulus,
            matrix_density,
            sma_density,
            thickness,
        }
    }
    /// Effective longitudinal modulus (Voigt rule of mixtures) \[Pa\].
    pub fn longitudinal_modulus(&self) -> f64 {
        let v_f = self.fiber_volume_fraction;
        let v_m = 1.0 - v_f;
        v_f * self.sma.elastic_modulus() + v_m * self.matrix_modulus
    }
    /// Effective transverse modulus (Reuss rule of mixtures) \[Pa\].
    pub fn transverse_modulus(&self) -> f64 {
        let v_f = self.fiber_volume_fraction;
        let v_m = 1.0 - v_f;
        1.0 / (v_f / self.sma.elastic_modulus() + v_m / self.matrix_modulus)
    }
    /// Effective density \[kg/m³\].
    pub fn density(&self) -> f64 {
        let v_f = self.fiber_volume_fraction;
        let v_m = 1.0 - v_f;
        v_f * self.sma_density + v_m * self.matrix_density
    }
    /// Recovery stress of the composite from SMA activation \[Pa\].
    pub fn composite_recovery_stress(&self, strain: f64) -> f64 {
        self.fiber_volume_fraction * self.sma.recovery_stress(strain)
    }
    /// Actuation curvature of an asymmetric laminate from thermal loading \[1/m\].
    ///
    /// Simplified: assumes SMA layer on one side and matrix on the other.
    pub fn actuation_curvature(&self, temperature: f64) -> f64 {
        let e_sma = self.sma.elastic_modulus();
        let e_mat = self.matrix_modulus;
        let t_sma = self.thickness * self.fiber_volume_fraction;
        let t_mat = self.thickness * (1.0 - self.fiber_volume_fraction);
        let eps_sma = self.sma.max_strain * (1.0 - self.sma.xi);
        let _ = temperature;
        6.0 * e_sma * e_mat * t_sma * t_mat * (t_sma + t_mat) * eps_sma
            / (e_sma * t_sma.powi(2) * (4.0 * t_sma + 3.0 * t_mat)
                + e_mat * t_mat.powi(2) * (4.0 * t_mat + 3.0 * t_sma)
                + 6.0 * e_sma * e_mat * t_sma * t_mat * (t_sma + t_mat))
                .max(1e-30)
    }
}
/// Electrostrictive material model.
///
/// Electrostriction is a quadratic electromechanical coupling:
/// ε_ij = Q_ijkl P_k P_l
///
/// where ε is strain, Q is the electrostriction coefficient tensor,
/// and P is electric polarization.
///
/// Implements:
/// - Maxwell stress tensor
/// - P-E relationship (nonlinear dielectric)
/// - Strain from polarization via Q coefficient
#[derive(Debug, Clone, Copy)]
pub struct Electrostrictive {
    /// Relative permittivity ε_r (linear regime).
    pub epsilon_r: f64,
    /// Electrostriction coefficient Q_33 \[m⁴/C²\].
    pub q33: f64,
    /// Electrostriction coefficient Q_31 \[m⁴/C²\].
    pub q31: f64,
    /// Nonlinear permittivity saturation coefficient \[V²/m²\].
    pub sat_field: f64,
    /// Young's modulus \[Pa\].
    pub elastic_modulus: f64,
    /// Material thickness \[m\].
    pub thickness: f64,
}
impl Electrostrictive {
    /// Create an electrostrictive material model.
    pub fn new(
        epsilon_r: f64,
        q33: f64,
        q31: f64,
        sat_field: f64,
        elastic_modulus: f64,
        thickness: f64,
    ) -> Self {
        Self {
            epsilon_r,
            q33,
            q31,
            sat_field,
            elastic_modulus,
            thickness,
        }
    }
    /// Standard PMN-PT (lead magnesium niobate - lead titanate) model.
    pub fn pmn_pt() -> Self {
        Self::new(5000.0, 0.026, -0.010, 2e7, 6e9, 1e-3)
    }
    /// Electric polarization P at applied field E \[C/m²\].
    ///
    /// Uses nonlinear dielectric: P = ε₀ (ε_r − 1) E / (1 + |E|/E_sat)
    pub fn polarization(&self, e_field: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        let chi = self.epsilon_r - 1.0;
        let denom = 1.0 + e_field.abs() / self.sat_field;
        EPS0 * chi * e_field / denom
    }
    /// Maxwell stress tensor component T_33 \[Pa\].
    ///
    /// T_33 = ε₀ ε_r E² / 2
    pub fn maxwell_stress(&self, e_field: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        0.5 * EPS0 * self.epsilon_r * e_field * e_field
    }
    /// Electrostrictive strain ε_33 from polarization P.
    ///
    /// ε_33 = Q_33 · P²
    pub fn strain_33(&self, e_field: f64) -> f64 {
        let p = self.polarization(e_field);
        self.q33 * p * p
    }
    /// Transverse electrostrictive strain ε_31.
    pub fn strain_31(&self, e_field: f64) -> f64 {
        let p = self.polarization(e_field);
        self.q31 * p * p
    }
    /// Actuation displacement Δt \[m\] (thickness change).
    pub fn thickness_change(&self, voltage: f64) -> f64 {
        let e_field = voltage / self.thickness;
        self.strain_33(e_field) * self.thickness
    }
    /// Blocking stress \[Pa\] for zero-displacement condition.
    pub fn blocking_stress(&self, voltage: f64) -> f64 {
        let e_field = voltage / self.thickness;
        self.elastic_modulus * self.strain_33(e_field)
    }
    /// Electromechanical coupling coefficient k_33.
    ///
    /// k²_33 = d²_33 / (s_33 ε_33) where d_33 = 2 Q_33 P ε_0 ε_r
    pub fn coupling_coefficient(&self, e_field: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        let p = self.polarization(e_field);
        let d33 = 2.0 * self.q33 * p * EPS0 * self.epsilon_r;
        let s33 = 1.0 / self.elastic_modulus;
        let eps33 = EPS0 * self.epsilon_r;
        let k2 = (d33 * d33) / (s33 * eps33);
        k2.sqrt().clamp(0.0, 1.0)
    }
}
/// Hydrogel actuator using Flory-Rehner swelling theory.
#[derive(Debug, Clone, Copy)]
pub struct HydrogelActuator {
    /// Flory-Huggins interaction parameter χ (temperature/pH dependent).
    pub chi: f64,
    /// Cross-link density \[mol/m³\].
    pub crosslink_density: f64,
    /// Reference swelling ratio Q_ref.
    pub q_ref: f64,
    /// Temperature coefficient dχ/dT \[1/K\].
    pub chi_temp_coeff: f64,
    /// pH at which χ has reference value.
    pub ph_ref: f64,
    /// pH sensitivity of χ.
    pub chi_ph_coeff: f64,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
}
impl HydrogelActuator {
    /// Create a hydrogel actuator.
    pub fn new(
        chi: f64,
        crosslink_density: f64,
        q_ref: f64,
        chi_temp_coeff: f64,
        ph_ref: f64,
        chi_ph_coeff: f64,
        elastic_modulus: f64,
    ) -> Self {
        Self {
            chi,
            crosslink_density,
            q_ref,
            chi_temp_coeff,
            ph_ref,
            chi_ph_coeff,
            elastic_modulus,
        }
    }
    /// Effective Flory parameter at temperature T \[K\] and pH.
    pub fn effective_chi(&self, temperature: f64, ph: f64) -> f64 {
        self.chi
            + self.chi_temp_coeff * (temperature - 298.0)
            + self.chi_ph_coeff * (ph - self.ph_ref)
    }
    /// Swelling ratio Q at given conditions (Flory-Rehner equilibrium, simplified).
    /// Q = V_swollen / V_dry.
    pub fn swelling_ratio(&self, temperature: f64, ph: f64) -> f64 {
        let chi_eff = self.effective_chi(temperature, ph);
        let q = self.q_ref * (-chi_eff * 0.5).exp();
        q.max(1.0)
    }
    /// Linear strain (relative to reference) from swelling ratio.
    pub fn linear_strain(&self, temperature: f64, ph: f64) -> f64 {
        let q = self.swelling_ratio(temperature, ph);
        q.powf(1.0 / 3.0) - 1.0
    }
    /// Swelling pressure (osmotic minus elastic) \[Pa\].
    pub fn swelling_pressure(&self, temperature: f64, ph: f64) -> f64 {
        let strain = self.linear_strain(temperature, ph);
        self.elastic_modulus * strain
    }
    /// Force generated by hydrogel strip of cross-section area A \[m²\].
    pub fn force(&self, temperature: f64, ph: f64, area: f64) -> f64 {
        self.swelling_pressure(temperature, ph) * area
    }
}
/// Electrorheological (ER) fluid model based on dielectric polarization.
#[derive(Debug, Clone, Copy)]
pub struct ElectrorheologicalFluid {
    /// Base (off-state) viscosity \[Pa·s\].
    pub base_viscosity: f64,
    /// Electric field coefficient for yield stress \[Pa·(m/V)²\].
    pub field_coeff: f64,
    /// Saturation field strength \[V/m\].
    pub saturation_field: f64,
    /// Dielectric constant of particles.
    pub epsilon_p: f64,
    /// Dielectric constant of carrier fluid.
    pub epsilon_f: f64,
}
impl ElectrorheologicalFluid {
    /// Create an ER fluid model.
    pub fn new(
        base_viscosity: f64,
        field_coeff: f64,
        saturation_field: f64,
        epsilon_p: f64,
        epsilon_f: f64,
    ) -> Self {
        Self {
            base_viscosity,
            field_coeff,
            saturation_field,
            epsilon_p,
            epsilon_f,
        }
    }
    /// Yield stress at electric field E \[V/m\].
    pub fn yield_stress(&self, e_field: f64) -> f64 {
        let e_eff = e_field.min(self.saturation_field);
        self.field_coeff * e_eff * e_eff
    }
    /// Shear stress using Bingham model with electric field.
    pub fn shear_stress(&self, shear_rate: f64, e_field: f64) -> f64 {
        let tau_y = self.yield_stress(e_field);
        if shear_rate.abs() < 1e-12 {
            return tau_y;
        }
        tau_y + self.base_viscosity * shear_rate
    }
    /// Dielectric loss factor (simplified).
    pub fn beta_parameter(&self) -> f64 {
        (self.epsilon_p - self.epsilon_f) / (self.epsilon_p + 2.0 * self.epsilon_f)
    }
}
/// Bimaterial beam (bimorph) actuated by differential thermal expansion.
#[derive(Debug, Clone, Copy)]
pub struct BimorphActuator {
    /// Length of beam \[m\].
    pub length: f64,
    /// Thickness of layer 1 \[m\].
    pub t1: f64,
    /// Thickness of layer 2 \[m\].
    pub t2: f64,
    /// Elastic modulus of layer 1 \[Pa\].
    pub e1: f64,
    /// Elastic modulus of layer 2 \[Pa\].
    pub e2: f64,
    /// Thermal expansion of layer 1 \[1/K\].
    pub alpha1: f64,
    /// Thermal expansion of layer 2 \[1/K\].
    pub alpha2: f64,
}
impl BimorphActuator {
    /// Create a bimorph actuator.
    pub fn new(length: f64, t1: f64, t2: f64, e1: f64, e2: f64, alpha1: f64, alpha2: f64) -> Self {
        Self {
            length,
            t1,
            t2,
            e1,
            e2,
            alpha1,
            alpha2,
        }
    }
    /// Tip deflection for temperature change ΔT \[m\].
    /// Uses Timoshenko bimetal formula.
    pub fn tip_deflection(&self, delta_t: f64) -> f64 {
        let m = self.e1 * self.t1 / (self.e2 * self.t2);
        let n = self.t1 / self.t2;
        let t_total = self.t1 + self.t2;
        let denom = 3.0 * (1.0 + m) * (1.0 + m) + (1.0 + m * n) * (m * m + 1.0 / (m * n));
        if denom.abs() < 1e-20 {
            return 0.0;
        }
        let curvature = 6.0 * (self.alpha1 - self.alpha2) * delta_t * (1.0 + m) / (t_total * denom);
        curvature * self.length * self.length / 2.0
    }
    /// Neutral axis position from layer 1 bottom \[m\].
    pub fn neutral_axis(&self) -> f64 {
        let n1 = self.e1 * self.t1;
        let n2 = self.e2 * self.t2;
        (n1 * self.t1 / 2.0 + n2 * (self.t1 + self.t2 / 2.0)) / (n1 + n2)
    }
    /// Curvature radius at ΔT \[m\].
    pub fn curvature_radius(&self, delta_t: f64) -> f64 {
        let deflection = self.tip_deflection(delta_t);
        if deflection.abs() < 1e-20 {
            return f64::INFINITY;
        }
        self.length * self.length / (2.0 * deflection)
    }
}
/// Ferroelectric P-E hysteresis loop using a simplified Preisach model.
///
/// The Preisach model represents the polarization as a superposition of
/// rectangular hysterons. Here we use a discretized 1D grid of hysterons.
#[derive(Debug, Clone)]
pub struct FerroelectricHysteresis {
    /// Saturation polarization P_sat \[C/m²\].
    pub p_sat: f64,
    /// Coercive field E_c \[V/m\].
    pub e_coercive: f64,
    /// Remanent polarization P_r \[C/m²\].
    pub p_remanent: f64,
    /// Number of Preisach hysteron bins.
    pub n_bins: usize,
    /// Hysteron states: +1 or -1.
    pub(super) hysteron_states: Vec<i8>,
    /// Hysteron switching fields (up/down) pairs.
    pub(super) hysteron_fields: Vec<[f64; 2]>,
    /// Current polarization \[C/m²\].
    pub polarization: f64,
}
impl FerroelectricHysteresis {
    /// Create a ferroelectric hysteresis model.
    pub fn new(p_sat: f64, e_coercive: f64, p_remanent: f64, n_bins: usize) -> Self {
        let mut hysteron_fields = Vec::with_capacity(n_bins);
        let mut hysteron_states = Vec::with_capacity(n_bins);
        for i in 0..n_bins {
            let frac = (i as f64 + 0.5) / n_bins as f64;
            let e_up = e_coercive * (0.5 + 1.5 * frac);
            let e_down = -e_up * (p_remanent / p_sat).max(0.1);
            hysteron_fields.push([e_down, e_up]);
            hysteron_states.push(-1i8);
        }
        Self {
            p_sat,
            e_coercive,
            p_remanent,
            n_bins,
            hysteron_states,
            hysteron_fields,
            polarization: -p_remanent,
        }
    }
    /// Standard BaTiO3 model.
    pub fn barium_titanate() -> Self {
        Self::new(0.26, 0.4e6, 0.18, 32)
    }
    /// Update polarization for applied field E \[V/m\].
    pub fn update(&mut self, e_field: f64) {
        for (i, &[e_down, e_up]) in self.hysteron_fields.iter().enumerate() {
            if e_field > e_up {
                self.hysteron_states[i] = 1;
            } else if e_field < e_down {
                self.hysteron_states[i] = -1;
            }
        }
        let sum: i32 = self.hysteron_states.iter().map(|&s| s as i32).sum();
        self.polarization = self.p_sat * sum as f64 / self.n_bins as f64;
    }
    /// Susceptibility at current state (finite difference).
    pub fn susceptibility(&self, e_field: f64, delta_e: f64) -> f64 {
        let mut twin = self.clone();
        twin.update(e_field + delta_e);
        let p1 = twin.polarization;
        let mut twin2 = self.clone();
        twin2.update(e_field - delta_e);
        let p2 = twin2.polarization;
        (p1 - p2) / (2.0 * delta_e)
    }
    /// Dielectric energy stored \[J/m³\].
    pub fn stored_energy(&self, e_field: f64) -> f64 {
        const EPS0: f64 = 8.854e-12;
        EPS0 * e_field * e_field / 2.0 + self.polarization * e_field
    }
}
