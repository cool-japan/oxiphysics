//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::f64::consts::PI;

/// Magnetorheological (MR) fluid modeled as a Bingham plastic.
#[derive(Debug, Clone, Copy)]
pub struct MagnetorheologicalFluid {
    /// Off-state viscosity \[Pa·s\].
    pub base_viscosity: f64,
    /// Field-independent yield stress offset \[Pa\].
    pub tau0_base: f64,
    /// Field-dependent yield stress coefficient \[Pa·(m/A)^alpha\].
    pub tau_h_coeff: f64,
    /// Field exponent (typically ~1-2).
    pub field_exponent: f64,
    /// Particle volume fraction.
    pub volume_fraction: f64,
}
impl MagnetorheologicalFluid {
    /// Create an MR fluid model.
    pub fn new(
        base_viscosity: f64,
        tau0_base: f64,
        tau_h_coeff: f64,
        field_exponent: f64,
        volume_fraction: f64,
    ) -> Self {
        Self {
            base_viscosity,
            tau0_base,
            tau_h_coeff,
            field_exponent,
            volume_fraction,
        }
    }
    /// Yield stress at magnetic field H \[A/m\].
    pub fn yield_stress(&self, h_field: f64) -> f64 {
        self.tau0_base + self.tau_h_coeff * h_field.powf(self.field_exponent)
    }
    /// Shear stress for given shear rate γ̇ \[1/s\] and field H \[A/m\] (Bingham model).
    pub fn shear_stress(&self, shear_rate: f64, h_field: f64) -> f64 {
        let tau_y = self.yield_stress(h_field);
        if shear_rate.abs() < 1e-12 {
            return tau_y;
        }
        tau_y + self.base_viscosity * shear_rate
    }
    /// Mason number: ratio of viscous to magnetic forces.
    /// Mn = η * γ̇ / (μ₀ * H²).
    pub fn mason_number(&self, shear_rate: f64, h_field: f64) -> f64 {
        const MU0: f64 = 4.0 * PI * 1e-7;
        if h_field < 1e-10 {
            return f64::INFINITY;
        }
        self.base_viscosity * shear_rate / (MU0 * h_field * h_field)
    }
    /// Apparent viscosity at given field and shear rate.
    pub fn apparent_viscosity(&self, shear_rate: f64, h_field: f64) -> f64 {
        if shear_rate.abs() < 1e-12 {
            return f64::INFINITY;
        }
        self.shear_stress(shear_rate, h_field) / shear_rate
    }
}
/// Shape memory alloy model using the Tanaka-Liang-Rogers cosine model.
/// Tracks martensite fraction ξ ∈ \[0, 1\], where 1 = full martensite.
#[derive(Debug, Clone)]
pub struct ShapeMemoryAlloy {
    /// Martensite start temperature \[K\].
    pub ms: f64,
    /// Martensite finish temperature \[K\].
    pub mf: f64,
    /// Austenite start temperature \[K\].
    pub a_s: f64,
    /// Austenite finish temperature \[K\].
    pub af: f64,
    /// Current martensite fraction ξ ∈ \[0, 1\].
    pub xi: f64,
    /// Austenite elastic modulus \[Pa\].
    pub e_a: f64,
    /// Martensite elastic modulus \[Pa\].
    pub e_m: f64,
    /// Maximum recoverable strain ε_L.
    pub max_strain: f64,
    /// Stress influence coefficient for martensite transformation \[Pa/K\].
    pub cm: f64,
    /// Stress influence coefficient for austenite transformation \[Pa/K\].
    pub ca: f64,
}
impl ShapeMemoryAlloy {
    /// Create an SMA model with transformation temperatures and elastic moduli.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ms: f64,
        mf: f64,
        a_s: f64,
        af: f64,
        e_a: f64,
        e_m: f64,
        max_strain: f64,
        cm: f64,
        ca: f64,
    ) -> Self {
        Self {
            ms,
            mf,
            a_s,
            af,
            xi: 1.0,
            e_a,
            e_m,
            max_strain,
            cm,
            ca,
        }
    }
    /// Create a standard NiTi Nitinol model with typical parameters.
    pub fn nitinol() -> Self {
        Self::new(291.0, 273.0, 307.0, 325.0, 75e9, 28e9, 0.08, 8e6, 13e6)
    }
    /// Elastic modulus at current martensite fraction.
    pub fn elastic_modulus(&self) -> f64 {
        self.e_a + self.xi * (self.e_m - self.e_a)
    }
    /// Update martensite fraction for given temperature T \[K\] and stress σ \[Pa\].
    /// Returns updated ξ.
    pub fn update_phase(&mut self, temperature: f64, stress: f64) -> f64 {
        let ms_stress = self.ms + stress / self.cm;
        let mf_stress = self.mf + stress / self.cm;
        if temperature <= ms_stress && temperature >= mf_stress {
            let xi_m =
                (1.0 + (PI * (temperature - ms_stress) / (mf_stress - ms_stress)).cos()) / 2.0;
            self.xi = xi_m.clamp(self.xi, 1.0);
        }
        let as_stress = self.a_s + stress / self.ca;
        let af_stress = self.af + stress / self.ca;
        if temperature >= as_stress && temperature <= af_stress {
            let xi_a = self.xi
                * (1.0 + (PI * (temperature - as_stress) / (af_stress - as_stress)).cos())
                / 2.0;
            self.xi = xi_a.clamp(0.0, self.xi);
        }
        if temperature < mf_stress {
            self.xi = 1.0;
        }
        if temperature > af_stress {
            self.xi = 0.0;
        }
        self.xi
    }
    /// Current phase classification.
    pub fn phase(&self) -> SmaPhase {
        if self.xi > 0.99 {
            SmaPhase::Martensite
        } else if self.xi < 0.01 {
            SmaPhase::Austenite
        } else {
            SmaPhase::Mixed
        }
    }
    /// Recovery stress for given strain \[Pa\].
    pub fn recovery_stress(&self, strain: f64) -> f64 {
        self.elastic_modulus() * (strain - self.xi * self.max_strain)
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
/// Ionic polymer-metal composite (IPMC) actuator model.
#[derive(Debug, Clone, Copy)]
pub struct PolymerActuator {
    /// Bending curvature constant \[1/(V·m)\].
    pub curvature_gain: f64,
    /// Time constant for bending response \[s\].
    pub time_constant: f64,
    /// Length of IPMC strip \[m\].
    pub length: f64,
    /// Thickness \[m\].
    pub thickness: f64,
    /// Blocking force per volt \[N/V\].
    pub blocking_force_gain: f64,
    /// Current bending state (time integral).
    pub bending_state: f64,
}
impl PolymerActuator {
    /// Create an IPMC actuator.
    pub fn new(
        curvature_gain: f64,
        time_constant: f64,
        length: f64,
        thickness: f64,
        blocking_force_gain: f64,
    ) -> Self {
        Self {
            curvature_gain,
            time_constant,
            length,
            thickness,
            blocking_force_gain,
            bending_state: 0.0,
        }
    }
    /// Steady-state tip deflection for step voltage V \[m\].
    pub fn steady_state_deflection(&self, voltage: f64) -> f64 {
        let kappa = self.curvature_gain * voltage;
        kappa * self.length * self.length / 2.0
    }
    /// Transient tip deflection at time t \[s\] after step voltage V.
    pub fn transient_deflection(&self, voltage: f64, t: f64) -> f64 {
        let delta_ss = self.steady_state_deflection(voltage);
        delta_ss * (1.0 - (-t / self.time_constant).exp())
    }
    /// Blocking force at given voltage \[N\].
    pub fn blocking_force(&self, voltage: f64) -> f64 {
        self.blocking_force_gain * voltage
    }
    /// Curvature for given voltage \[1/m\].
    pub fn curvature(&self, voltage: f64) -> f64 {
        self.curvature_gain * voltage
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
/// Bang-bang + proportional controller for SMA wire actuator temperature.
#[derive(Debug, Clone)]
pub struct SmaController {
    /// Reference SMA model.
    pub sma: ShapeMemoryAlloy,
    /// Target martensite fraction.
    pub target_xi: f64,
    /// Temperature setpoint for full austenite \[K\].
    pub t_hot: f64,
    /// Temperature setpoint for full martensite \[K\].
    pub t_cold: f64,
    /// Proportional gain for temperature error.
    pub kp: f64,
    /// Current temperature \[K\].
    pub temperature: f64,
    /// Maximum heating power \[W\].
    pub max_power: f64,
}
impl SmaController {
    /// Create an SMA controller.
    pub fn new(sma: ShapeMemoryAlloy, t_hot: f64, t_cold: f64, kp: f64, max_power: f64) -> Self {
        let temperature = t_cold;
        Self {
            sma,
            target_xi: 0.0,
            t_hot,
            t_cold,
            kp,
            temperature,
            max_power,
        }
    }
    /// Set target martensite fraction.
    pub fn set_target(&mut self, xi_target: f64) {
        self.target_xi = xi_target.clamp(0.0, 1.0);
    }
    /// Compute heating power to achieve target (bang-bang + proportional).
    pub fn compute_power(&self) -> f64 {
        let xi_err = self.target_xi - self.sma.xi;
        if xi_err > 0.1 {
            0.0
        } else if xi_err < -0.1 {
            self.max_power
        } else {
            (-self.kp * xi_err * self.max_power).clamp(0.0, self.max_power)
        }
    }
    /// Step the controller by dt \[s\] with thermal resistance R_th \[K/W\].
    /// `t_ambient` is environment temperature \[K\].
    pub fn step(&mut self, dt: f64, r_thermal: f64, t_ambient: f64) {
        let power = self.compute_power();
        let thermal_mass = 1e-4;
        let q_loss = (self.temperature - t_ambient) / r_thermal;
        let d_temp = (power - q_loss) / thermal_mass * dt;
        self.temperature = (self.temperature + d_temp).clamp(self.t_cold, self.t_hot);
        self.sma.update_phase(self.temperature, 0.0);
    }
}
/// Distributed actuator and sensor layout for a smart structure.
#[derive(Debug, Clone)]
pub struct SmartStructure {
    /// Number of structural modes.
    pub n_modes: usize,
    /// Modal frequencies \[Hz\].
    pub modal_frequencies: Vec<f64>,
    /// Modal damping ratios.
    pub modal_damping: Vec<f64>,
    /// Actuator positions (normalized 0..1 along structure).
    pub actuator_positions: Vec<f64>,
    /// Sensor positions.
    pub sensor_positions: Vec<f64>,
    /// Modal state: (displacement, velocity) per mode.
    pub modal_state: Vec<[f64; 2]>,
    /// Control gains per mode.
    pub control_gains: Vec<f64>,
}
impl SmartStructure {
    /// Create a smart structure model.
    pub fn new(
        n_modes: usize,
        modal_frequencies: Vec<f64>,
        modal_damping: Vec<f64>,
        actuator_positions: Vec<f64>,
        sensor_positions: Vec<f64>,
    ) -> Self {
        let control_gains = vec![1.0; n_modes];
        let modal_state = vec![[0.0; 2]; n_modes];
        Self {
            n_modes,
            modal_frequencies,
            modal_damping,
            actuator_positions,
            sensor_positions,
            modal_state,
            control_gains,
        }
    }
    /// Mode shape (assumed sinusoidal basis): φ_n(x) = sin(nπx).
    fn mode_shape(&self, mode: usize, position: f64) -> f64 {
        ((mode + 1) as f64 * PI * position).sin()
    }
    /// Compute modal control forces from sensor readings.
    pub fn modal_control(&self, sensor_readings: &[f64]) -> Vec<f64> {
        let mut u = vec![0.0; self.actuator_positions.len()];
        for (k, &x_act) in self.actuator_positions.iter().enumerate() {
            for m in 0..self.n_modes {
                let obs: f64 = sensor_readings
                    .iter()
                    .zip(self.sensor_positions.iter())
                    .map(|(&y, &x_s)| y * self.mode_shape(m, x_s))
                    .sum();
                let vel = self.modal_state[m][1];
                let phi_act = self.mode_shape(m, x_act);
                u[k] -= self.control_gains[m] * (obs + vel) * phi_act;
            }
        }
        u
    }
    /// Integrate modal equations of motion under distributed forcing.
    pub fn step(&mut self, forcing: &[f64], dt: f64) {
        for m in 0..self.n_modes {
            let omega = 2.0 * PI * self.modal_frequencies[m];
            let zeta = self.modal_damping[m];
            let q = self.modal_state[m][0];
            let qdot = self.modal_state[m][1];
            let f_gen: f64 = forcing
                .iter()
                .zip(self.actuator_positions.iter())
                .map(|(&f, &x)| f * self.mode_shape(m, x))
                .sum();
            let qddot = f_gen - 2.0 * zeta * omega * qdot - omega * omega * q;
            self.modal_state[m][0] = q + qdot * dt;
            self.modal_state[m][1] = qdot + qddot * dt;
        }
    }
    /// Physical displacement at position x from modal superposition.
    pub fn displacement_at(&self, x: f64) -> f64 {
        (0..self.n_modes)
            .map(|m| self.modal_state[m][0] * self.mode_shape(m, x))
            .sum()
    }
    /// Set control gains to avoid spillover (zero out high modes).
    pub fn set_control_gains(&mut self, gains: Vec<f64>) {
        let n = gains.len().min(self.n_modes);
        self.control_gains[..n].copy_from_slice(&gains[..n]);
    }
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
/// Type of electroactive polymer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EapType {
    /// Ionic EAP: actuation by ion migration (IPMC, conductive polymer).
    Ionic,
    /// Electronic EAP: actuation by electrostatic forces (dielectric elastomer, PVDF).
    Electronic,
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
/// Brinson constitutive model for shape memory alloys.
///
/// Explicitly tracks stress-induced martensite (ξ_s) and
/// temperature-induced martensite (ξ_T) fractions, enabling full
/// simulation of:
/// - Shape Memory Effect (SME): thermally driven full recovery
/// - Superelasticity: stress-driven reversible strain
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BrinsonModel {
    /// Martensite start temperature \[K\].
    pub ms: f64,
    /// Martensite finish temperature \[K\].
    pub mf: f64,
    /// Austenite start temperature \[K\].
    pub a_s: f64,
    /// Austenite finish temperature \[K\].
    pub af: f64,
    /// Stress-induced martensite fraction ξ_s ∈ \[0, 1\].
    pub xi_s: f64,
    /// Temperature-induced martensite fraction ξ_T ∈ \[0, 1\].
    pub xi_t: f64,
    /// Austenite elastic modulus \[Pa\].
    pub e_a: f64,
    /// Martensite elastic modulus \[Pa\].
    pub e_m: f64,
    /// Maximum recoverable transformation strain ε_L.
    pub eps_l: f64,
    /// Stress influence coefficient for M transformation \[Pa/K\].
    pub cm: f64,
    /// Stress influence coefficient for A transformation \[Pa/K\].
    pub ca: f64,
    /// Thermoelastic coefficient Θ \[Pa/K\].
    pub theta: f64,
    /// Current stress \[Pa\].
    pub stress: f64,
    /// Current strain.
    pub strain: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}
impl BrinsonModel {
    /// Create a Brinson model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ms: f64,
        mf: f64,
        a_s: f64,
        af: f64,
        e_a: f64,
        e_m: f64,
        eps_l: f64,
        cm: f64,
        ca: f64,
        theta: f64,
    ) -> Self {
        Self {
            ms,
            mf,
            a_s,
            af,
            xi_s: 0.0,
            xi_t: 1.0,
            e_a,
            e_m,
            eps_l,
            cm,
            ca,
            theta,
            stress: 0.0,
            strain: 0.0,
            temperature: mf,
        }
    }
    /// Create a standard NiTi model (Brinson parameters).
    pub fn nitinol() -> Self {
        Self::new(
            291.0, 273.0, 307.0, 325.0, 75e9, 28e9, 0.08, 8e6, 13e6, 0.55e6,
        )
    }
    /// Total martensite fraction ξ = ξ_s + ξ_T.
    pub fn total_xi(&self) -> f64 {
        (self.xi_s + self.xi_t).clamp(0.0, 1.0)
    }
    /// Effective elastic modulus E(ξ) = E_A + ξ(E_M − E_A).
    pub fn elastic_modulus(&self) -> f64 {
        self.e_a + self.total_xi() * (self.e_m - self.e_a)
    }
    /// Update phase fractions from temperature and stress using Brinson kinetics.
    ///
    /// Returns `(ξ_s, ξ_T, ξ_total)` after update.
    pub fn update(&mut self, temperature: f64, stress: f64) -> (f64, f64, f64) {
        self.temperature = temperature;
        self.stress = stress;
        let ms_s = self.ms + stress / self.cm;
        let mf_s = self.mf + stress / self.cm;
        if temperature <= ms_s && temperature >= mf_s {
            let cos_arg = PI * (temperature - ms_s) / (mf_s - ms_s);
            let xi_new = (1.0 + cos_arg.cos()) / 2.0;
            if xi_new > self.xi_t {
                self.xi_t = xi_new.min(1.0 - self.xi_s);
            }
        } else if temperature < mf_s {
            self.xi_t = 1.0 - self.xi_s;
        }
        let sigma_ms = self.cm * (temperature - self.ms).max(0.0);
        let sigma_mf = self.cm * (temperature - self.mf).max(0.0);
        if stress >= sigma_ms && stress <= sigma_mf + 1e-3 && temperature > self.mf {
            let cos_arg = PI * (stress - sigma_ms) / (sigma_mf - sigma_ms + 1e-6);
            let xi_s_new = (1.0 - self.xi_s) / 2.0 * (1.0 - cos_arg.cos());
            self.xi_s = (self.xi_s + xi_s_new).clamp(0.0, 1.0 - self.xi_t);
        }
        let as_s = self.a_s + stress / self.ca;
        let af_s = self.af + stress / self.ca;
        if temperature >= as_s && temperature <= af_s {
            let cos_arg = PI * (temperature - as_s) / (af_s - as_s);
            let xi_new = self.total_xi() / 2.0 * (1.0 + cos_arg.cos());
            let xi_total = xi_new.clamp(0.0, self.total_xi());
            let frac = if self.total_xi() > 1e-10 {
                xi_total / self.total_xi()
            } else {
                0.0
            };
            self.xi_s = (self.xi_s * frac).clamp(0.0, 1.0);
            self.xi_t = (self.xi_t * frac).clamp(0.0, 1.0);
        } else if temperature > af_s {
            self.xi_s = 0.0;
            self.xi_t = 0.0;
        }
        (self.xi_s, self.xi_t, self.total_xi())
    }
    /// Compute stress from strain using Brinson constitutive law.
    ///
    /// σ = E(ξ)·ε − E(ξ)·ε_L·ξ_s + Θ·(T − T_ref)
    pub fn stress_from_strain(&self, strain: f64, t_ref: f64) -> f64 {
        let e = self.elastic_modulus();
        e * (strain - self.eps_l * self.xi_s) + self.theta * (self.temperature - t_ref)
    }
    /// Recovery stress when cooled from austenite to martensite at fixed strain \[Pa\].
    pub fn recovery_stress(&self, strain: f64) -> f64 {
        self.elastic_modulus() * (strain - self.eps_l * self.xi_s)
    }
    /// SME stroke: change in length from xi_initial to xi_final \[m\] for wire of length L.
    pub fn sme_stroke(&self, xi_initial: f64, xi_final: f64, wire_length: f64) -> f64 {
        let delta_xi_s = xi_initial - xi_final;
        delta_xi_s * self.eps_l * wire_length
    }
}
/// Extended magnetorheological fluid model including particle chain formation.
///
/// Uses Bingham plastic with field-dependent yield stress derived from:
/// 1. Mason number analysis
/// 2. Chain formation model (dipole-dipole interactions)
/// 3. Off-state viscosity correction with volume fraction (Krieger-Dougherty)
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct Magnetorheological {
    /// Off-state dynamic viscosity η_0 \[Pa·s\].
    pub eta0: f64,
    /// Magnetic permeability of carrier fluid \[H/m\].
    pub mu_c: f64,
    /// Saturation magnetization of particles M_s \[A/m\].
    pub m_sat: f64,
    /// Particle volume fraction φ.
    pub phi: f64,
    /// Field-exponent for yield stress correlation.
    pub n_exp: f64,
    /// Empirical yield stress coefficient c_y \[Pa·(m/A)^n\].
    pub c_yield: f64,
    /// Maximum packing fraction φ_m.
    pub phi_max: f64,
    /// Intrinsic viscosity \[η\] for K-D equation.
    pub intrinsic_viscosity: f64,
}
impl Magnetorheological {
    /// Create an MR fluid model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        eta0: f64,
        mu_c: f64,
        m_sat: f64,
        phi: f64,
        n_exp: f64,
        c_yield: f64,
        phi_max: f64,
        intrinsic_viscosity: f64,
    ) -> Self {
        Self {
            eta0,
            mu_c,
            m_sat,
            phi,
            n_exp,
            c_yield,
            phi_max,
            intrinsic_viscosity,
        }
    }
    /// Standard carbonyl iron / silicone oil MR fluid.
    pub fn carbonyl_iron() -> Self {
        Self::new(0.1, 4.0 * PI * 1e-7, 1.36e6, 0.30, 1.5, 0.3e3, 0.74, 2.5)
    }
    /// Effective off-state viscosity using Krieger-Dougherty model.
    pub fn effective_viscosity_off(&self) -> f64 {
        let ratio = self.phi / self.phi_max;
        let kd = (1.0 - ratio).powf(-self.intrinsic_viscosity * self.phi_max);
        self.eta0 * kd.max(1.0)
    }
    /// Field-dependent yield stress τ_y(H) \[Pa\].
    pub fn yield_stress(&self, h_field: f64) -> f64 {
        self.c_yield * h_field.powf(self.n_exp)
    }
    /// Mason number Mn = η·γ̇ / (μ_0 H²).
    pub fn mason_number(&self, shear_rate: f64, h_field: f64) -> f64 {
        const MU0: f64 = 1.2566370614e-6;
        if h_field < 1e-10 {
            return f64::INFINITY;
        }
        self.eta0 * shear_rate / (MU0 * h_field * h_field)
    }
    /// Shear stress using Bingham model: τ = τ_y(H) + η_eff · γ̇.
    pub fn shear_stress(&self, shear_rate: f64, h_field: f64) -> f64 {
        let tau_y = self.yield_stress(h_field);
        let eta_eff = self.effective_viscosity_off();
        if shear_rate.abs() < 1e-12 {
            return tau_y;
        }
        tau_y + eta_eff * shear_rate
    }
    /// Apparent viscosity η_app = τ / γ̇ \[Pa·s\].
    pub fn apparent_viscosity(&self, shear_rate: f64, h_field: f64) -> f64 {
        if shear_rate.abs() < 1e-12 {
            return f64::INFINITY;
        }
        self.shear_stress(shear_rate, h_field) / shear_rate
    }
    /// Number of chains formed per unit volume (simplified model).
    ///
    /// N_chains ∝ φ / d³ where d is particle diameter.
    /// Returns a dimensionless chain fraction estimate.
    pub fn chain_fraction(&self, h_field: f64) -> f64 {
        const MU0: f64 = 1.2566370614e-6;
        let m = MU0 * self.m_sat * h_field;
        let lambda = m / (self.eta0.max(1e-10) + 1.0);
        (lambda * self.phi).tanh().clamp(0.0, 1.0)
    }
    /// Relative viscosity enhancement due to chain formation.
    ///
    /// Returns η_eff / η_0.
    pub fn viscosity_ratio(&self, h_field: f64) -> f64 {
        let cf = self.chain_fraction(h_field);
        1.0 + 100.0 * cf * self.phi
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
/// Extended hydrogel swelling model using full Flory-Rehner theory.
///
/// Models the equilibrium between:
/// - Mixing free energy (Flory-Huggins): ΔF_mix = kT\[φ ln φ + (1-φ) ln(1-φ) + χφ(1-φ)\]
/// - Elastic free energy (affine network): ΔF_el = (3νkT/2)(λ² − 1 − 2 ln λ)
///
/// where φ is polymer volume fraction and λ = Q^(1/3) is the stretch ratio.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct HydrogelFloryRehner {
    /// Flory-Huggins interaction parameter χ (dimensionless).
    pub chi: f64,
    /// Effective cross-link density ν \[mol/m³\].
    pub crosslink_density: f64,
    /// Polymer volume fraction at synthesis φ₀.
    pub phi0: f64,
    /// Molar volume of solvent \[m³/mol\].
    pub v_s: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Elastic modulus of the dry network \[Pa\].
    pub dry_modulus: f64,
}
impl HydrogelFloryRehner {
    /// Create a Flory-Rehner hydrogel model.
    pub fn new(
        chi: f64,
        crosslink_density: f64,
        phi0: f64,
        v_s: f64,
        temperature: f64,
        dry_modulus: f64,
    ) -> Self {
        Self {
            chi,
            crosslink_density,
            phi0,
            v_s,
            temperature,
            dry_modulus,
        }
    }
    /// Standard polyacrylamide hydrogel.
    pub fn polyacrylamide() -> Self {
        Self::new(0.45, 10.0, 0.05, 18e-6, 298.0, 1e4)
    }
    /// Mixing free energy density \[J/m³\] at volume fraction φ.
    pub fn mixing_free_energy(&self, phi: f64) -> f64 {
        const R: f64 = 8.314;
        let phi = phi.clamp(1e-6, 1.0 - 1e-6);
        let psi = 1.0 - phi;
        R * self.temperature / self.v_s * (phi * phi.ln() + psi * psi.ln() + self.chi * phi * psi)
    }
    /// Elastic free energy density \[J/m³\] at volume fraction φ.
    ///
    /// Uses affine network model: f_el = ν kT / 2 · (3λ² − 3 − 2 ln λ³)
    pub fn elastic_free_energy(&self, phi: f64) -> f64 {
        const R: f64 = 8.314;
        let phi = phi.clamp(1e-6, 1.0);
        let lambda = (self.phi0 / phi).powf(1.0 / 3.0);
        let lambda2 = lambda * lambda;
        let ln_lambda3 = 3.0 * lambda.ln();
        self.crosslink_density * R * self.temperature / 2.0
            * (3.0 * lambda2 - 3.0 - 2.0 * ln_lambda3)
    }
    /// Total free energy density \[J/m³\].
    pub fn total_free_energy(&self, phi: f64) -> f64 {
        self.mixing_free_energy(phi) + self.elastic_free_energy(phi)
    }
    /// Chemical potential of solvent μ (dimensionless, in kT units).
    ///
    /// μ = ∂(V_s · f_mix) / ∂(1-φ)
    pub fn chemical_potential(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-6, 1.0 - 1e-6);
        let psi = 1.0 - phi;
        psi.ln() + phi + self.chi * phi * phi
    }
    /// Elastic chemical potential contribution from network elasticity.
    pub fn elastic_chemical_potential(&self, phi: f64) -> f64 {
        let phi = phi.clamp(1e-6, 1.0);
        let n = (self.crosslink_density * self.v_s).max(1e-10);
        (phi / (2.0 * n)) * (phi / self.phi0).powf(1.0 / 3.0)
    }
    /// Find equilibrium swelling ratio Q = V_wet / V_dry by Newton iteration.
    pub fn equilibrium_swelling(&self) -> f64 {
        let mut phi = self.phi0.clamp(0.01, 0.99);
        for _ in 0..100 {
            let mu_mix = self.chemical_potential(phi);
            let mu_el = self.elastic_chemical_potential(phi);
            let mu = mu_mix + mu_el;
            let dphi = 1e-6;
            let d_mu = (self.chemical_potential(phi + dphi)
                + self.elastic_chemical_potential(phi + dphi)
                - self.chemical_potential(phi - dphi)
                - self.elastic_chemical_potential(phi - dphi))
                / (2.0 * dphi);
            if d_mu.abs() < 1e-15 {
                break;
            }
            phi -= mu / d_mu;
            phi = phi.clamp(0.001, 0.999);
        }
        1.0 / phi.max(1e-6)
    }
    /// Linear swelling strain ε = Q^(1/3) − 1.
    pub fn linear_strain(&self) -> f64 {
        self.equilibrium_swelling().powf(1.0 / 3.0) - 1.0
    }
    /// Osmotic pressure Π = −∂f_total/∂V at given φ \[Pa\].
    pub fn osmotic_pressure(&self, phi: f64) -> f64 {
        const R: f64 = 8.314;
        let phi = phi.clamp(1e-6, 1.0 - 1e-6);
        let psi = 1.0 - phi;
        let mu1 = psi.ln() + phi + self.chi * phi * phi;
        -(R * self.temperature / self.v_s) * mu1
    }
}
/// Self-healing material model.
///
/// Models:
/// - Damage-triggered healing kinetics (first-order reaction)
/// - Healing efficiency η_h as a function of time and temperature
/// - Strength recovery: σ_recovery(t) = σ₀ · η_h(t)
/// - Activation energy dependence (Arrhenius kinetics)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SelfHealingMaterial {
    /// Undamaged tensile strength σ₀ \[Pa\].
    pub virgin_strength: f64,
    /// Current damage fraction D ∈ \[0, 1\].
    pub damage: f64,
    /// Cumulative healing fraction η_h ∈ \[0, 1\].
    pub healing_efficiency: f64,
    /// Healing kinetics rate constant k_h \[1/s\] at reference temperature.
    pub k_healing: f64,
    /// Activation energy for healing reaction \[J/mol\].
    pub activation_energy: f64,
    /// Reference temperature for k_h \[K\].
    pub t_ref: f64,
    /// Healing agent volume fraction (0..1).
    pub agent_fraction: f64,
    /// Critical damage for healing activation (threshold D > D_c triggers healing).
    pub damage_threshold: f64,
    /// Time elapsed since damage \[s\].
    pub time_since_damage: f64,
    /// Maximum achievable healing efficiency (intrinsic limit).
    pub max_efficiency: f64,
}
impl SelfHealingMaterial {
    /// Create a self-healing material model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        virgin_strength: f64,
        k_healing: f64,
        activation_energy: f64,
        t_ref: f64,
        agent_fraction: f64,
        damage_threshold: f64,
        max_efficiency: f64,
    ) -> Self {
        Self {
            virgin_strength,
            damage: 0.0,
            healing_efficiency: 0.0,
            k_healing,
            activation_energy,
            t_ref,
            agent_fraction,
            damage_threshold,
            time_since_damage: 0.0,
            max_efficiency,
        }
    }
    /// Standard microencapsulated epoxy self-healing composite.
    pub fn microencapsulated_epoxy() -> Self {
        Self::new(60e6, 1e-4, 50e3, 298.0, 0.10, 0.05, 0.90)
    }
    /// Arrhenius rate constant at temperature T \[K\].
    pub fn healing_rate(&self, temperature: f64) -> f64 {
        const R: f64 = 8.314;
        self.k_healing
            * (-self.activation_energy / (R * temperature)).exp()
            * (-self.activation_energy / (R * self.t_ref)).exp().recip()
    }
    /// Apply damage D_new (D ∈ \[0, 1\]) to the material.
    pub fn apply_damage(&mut self, damage: f64) {
        self.damage = (self.damage + damage).clamp(0.0, 1.0);
        if self.damage > self.damage_threshold {
            self.time_since_damage = 0.0;
            self.healing_efficiency = 0.0;
        }
    }
    /// Advance healing by time `dt` at temperature `temperature`.
    ///
    /// Updates healing efficiency using first-order kinetics:
    /// dη/dt = k_h(T) · (η_max − η) · agent_fraction
    pub fn heal(&mut self, dt: f64, temperature: f64) {
        if self.damage < self.damage_threshold {
            return;
        }
        let k = self.healing_rate(temperature);
        let rate = k * (self.max_efficiency - self.healing_efficiency) * self.agent_fraction;
        self.healing_efficiency =
            (self.healing_efficiency + rate * dt).clamp(0.0, self.max_efficiency);
        self.time_since_damage += dt;
    }
    /// Current tensile strength \[Pa\] accounting for damage and healing.
    ///
    /// σ = σ₀ · (1 − D) · (1 + η_h · D)
    pub fn current_strength(&self) -> f64 {
        let intact_fraction = 1.0 - self.damage;
        let healed_contribution = self.healing_efficiency * self.damage;
        self.virgin_strength * (intact_fraction + healed_contribution).clamp(0.0, 1.0)
    }
    /// Effective stiffness reduction factor (CDM approach).
    ///
    /// E_eff / E₀ = 1 − D · (1 − η_h)
    pub fn stiffness_reduction(&self) -> f64 {
        (1.0 - self.damage * (1.0 - self.healing_efficiency)).clamp(0.0, 1.0)
    }
    /// Fracture toughness recovery ratio K_IC / K_IC0.
    ///
    /// Approximated as sqrt of stiffness reduction (Irwin relation).
    pub fn toughness_recovery(&self) -> f64 {
        self.stiffness_reduction().sqrt()
    }
    /// Healing completion fraction (η_h / η_max).
    pub fn healing_progress(&self) -> f64 {
        if self.max_efficiency < 1e-10 {
            0.0
        } else {
            self.healing_efficiency / self.max_efficiency
        }
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
/// Magnetostrictive material: strain vs applied magnetic field H.
///
/// Uses a simplified Langevin-based saturation model for Joule magnetostriction.
/// λ(H) = λ_sat * \[coth(H/H_sat) - H_sat/H\]³  (modified Langevin)
#[derive(Debug, Clone, Copy)]
pub struct MagnetostrictiveMaterial {
    /// Saturation magnetostriction λ_sat (dimensionless, e.g. 1000e-6 for Terfenol-D).
    pub lambda_sat: f64,
    /// Characteristic (saturation-scale) field H_sat \[A/m\].
    pub h_sat: f64,
    /// Young's modulus \[Pa\].
    pub elastic_modulus: f64,
    /// Magnetomechanical coupling factor d_33 \[m/A\].
    pub d33: f64,
    /// Relative permeability.
    pub mu_r: f64,
}
impl MagnetostrictiveMaterial {
    /// Create a magnetostrictive material model.
    pub fn new(lambda_sat: f64, h_sat: f64, elastic_modulus: f64, d33: f64, mu_r: f64) -> Self {
        Self {
            lambda_sat,
            h_sat,
            elastic_modulus,
            d33,
            mu_r,
        }
    }
    /// Standard Terfenol-D model.
    pub fn terfenol_d() -> Self {
        Self::new(1600e-6, 60e3, 50e9, 20e-9, 3.0)
    }
    /// Magnetostrictive strain λ at field H \[A/m\] using Langevin function.
    pub fn strain(&self, h_field: f64) -> f64 {
        if h_field.abs() < 1e-10 {
            return 0.0;
        }
        let x = h_field / self.h_sat;
        let l = if x.abs() > 30.0 {
            x.signum() * (1.0 - 1.0 / x.abs())
        } else {
            1.0 / x.tanh() - 1.0 / x
        };
        self.lambda_sat * l * l * h_field.signum().max(0.0)
            + self.lambda_sat * l * l * (-h_field.signum()).max(0.0)
    }
    /// Magnetostrictive strain with correct even-symmetry (λ ≥ 0 always).
    pub fn strain_magnitude(&self, h_field: f64) -> f64 {
        if h_field.abs() < 1e-10 {
            return 0.0;
        }
        let x = h_field.abs() / self.h_sat;
        let l = if x > 30.0 {
            1.0 - 1.0 / x
        } else {
            1.0 / x.tanh() - 1.0 / x
        };
        self.lambda_sat * l * l
    }
    /// Blocked stress (zero strain) at field H \[Pa\].
    pub fn blocked_stress(&self, h_field: f64) -> f64 {
        self.elastic_modulus * self.strain_magnitude(h_field)
    }
    /// Piezomagnetic coefficient dλ/dH \[1/(A/m)\] at field H.
    pub fn piezomagnetic_coeff(&self, h_field: f64) -> f64 {
        let dh = h_field * 1e-4 + 1.0;
        (self.strain_magnitude(h_field + dh) - self.strain_magnitude(h_field)) / dh
    }
}
/// Magnetic shape memory effect in Heusler alloys (Ni-Mn-Ga).
///
/// Models variant reorientation driven by magnetic field via energy minimization.
/// Strain output up to ~12 % for Ni2MnGa.
#[derive(Debug, Clone, Copy)]
pub struct MagneticShape {
    /// Maximum transformation strain ε_max.
    pub max_strain: f64,
    /// Magnetic anisotropy energy density K_u \[J/m³\].
    pub k_u: f64,
    /// Saturation magnetization M_s \[A/m\].
    pub m_sat: f64,
    /// Critical field for reorientation H_cr \[A/m\].
    pub h_critical: f64,
    /// Current variant fraction η ∈ \[0, 1\].
    pub eta: f64,
    /// Elastic modulus \[Pa\].
    pub elastic_modulus: f64,
}
impl MagneticShape {
    /// Create a magnetic shape memory material.
    pub fn new(
        max_strain: f64,
        k_u: f64,
        m_sat: f64,
        h_critical: f64,
        elastic_modulus: f64,
    ) -> Self {
        Self {
            max_strain,
            k_u,
            m_sat,
            h_critical,
            eta: 0.0,
            elastic_modulus,
        }
    }
    /// Standard Ni2MnGa model.
    pub fn ni2mnga() -> Self {
        Self::new(0.06, 1.65e5, 600e3, 300e3, 2.0e9)
    }
    /// Zeeman energy density at field H (field along easy axis) \[J/m³\].
    pub fn zeeman_energy(&self, h_field: f64, eta: f64) -> f64 {
        const MU0: f64 = 4.0 * PI * 1e-7;
        -MU0 * self.m_sat * h_field * eta
    }
    /// Anisotropy energy density \[J/m³\].
    pub fn anisotropy_energy(&self, eta: f64) -> f64 {
        self.k_u * eta * (1.0 - eta)
    }
    /// Update variant fraction η from applied field H \[A/m\].
    pub fn update(&mut self, h_field: f64) {
        if h_field.abs() > self.h_critical {
            self.eta = if h_field > 0.0 { 1.0 } else { 0.0 };
        } else {
            let fraction = (h_field / self.h_critical).clamp(-1.0, 1.0);
            self.eta = (0.5 + 0.5 * fraction).clamp(0.0, 1.0);
        }
    }
    /// Macroscopic strain at current η.
    pub fn strain(&self) -> f64 {
        self.eta * self.max_strain
    }
    /// Blocking stress (stress required to prevent strain) \[Pa\].
    pub fn blocking_stress(&self) -> f64 {
        self.elastic_modulus * self.strain()
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
#[allow(dead_code)]
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
/// Thermochromic material: color/reflectivity changes with temperature.
///
/// Models an RGB reflectance spectrum shift through a transition temperature
/// using a smooth sigmoid function.
#[derive(Debug, Clone)]
pub struct ThermochromicResponse {
    /// Transition temperature \[K\].
    pub t_transition: f64,
    /// Transition width (steepness) \[K\].
    pub delta_t: f64,
    /// Low-temperature RGB reflectance (0..1 each channel).
    pub color_low: [f64; 3],
    /// High-temperature RGB reflectance (0..1 each channel).
    pub color_high: [f64; 3],
    /// Latent heat of transition \[J/kg\].
    pub latent_heat: f64,
}
impl ThermochromicResponse {
    /// Create a thermochromic material.
    pub fn new(
        t_transition: f64,
        delta_t: f64,
        color_low: [f64; 3],
        color_high: [f64; 3],
        latent_heat: f64,
    ) -> Self {
        Self {
            t_transition,
            delta_t,
            color_low,
            color_high,
            latent_heat,
        }
    }
    /// Transition fraction s ∈ \[0, 1\] at temperature T.
    ///
    /// s = 0 → low-T color, s = 1 → high-T color.
    pub fn transition_fraction(&self, temperature: f64) -> f64 {
        let x = (temperature - self.t_transition) / self.delta_t;
        1.0 / (1.0 + (-x).exp())
    }
    /// RGB reflectance \[r, g, b\] at temperature T.
    pub fn reflectance(&self, temperature: f64) -> [f64; 3] {
        let s = self.transition_fraction(temperature);
        [
            self.color_low[0] + s * (self.color_high[0] - self.color_low[0]),
            self.color_low[1] + s * (self.color_high[1] - self.color_low[1]),
            self.color_low[2] + s * (self.color_high[2] - self.color_low[2]),
        ]
    }
    /// Grayscale luminance at temperature T.
    pub fn luminance(&self, temperature: f64) -> f64 {
        let [r, g, b] = self.reflectance(temperature);
        0.2126 * r + 0.7152 * g + 0.0722 * b
    }
    /// Specific heat capacity including latent heat peak \[J/(kg·K)\].
    ///
    /// c_eff(T) = c_base + L * ds/dT
    pub fn effective_heat_capacity(&self, temperature: f64, c_base: f64) -> f64 {
        let s = self.transition_fraction(temperature);
        let ds_dt = s * (1.0 - s) / self.delta_t;
        c_base + self.latent_heat * ds_dt
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
#[allow(dead_code)]
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
    #[allow(clippy::too_many_arguments)]
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
/// SMA wire actuator using Brinson constitutive model.
///
/// Models recovery stress and stroke output during a temperature cycle.
/// The Brinson model separates martensite into stress-induced (ξs)
/// and temperature-induced (ξT) fractions.
#[derive(Debug, Clone)]
pub struct SmaActuator {
    /// Wire cross-sectional area \[m²\].
    pub area: f64,
    /// Wire rest length \[m\].
    pub length: f64,
    /// Stress-induced martensite fraction ξ_s ∈ \[0, 1\].
    pub xi_s: f64,
    /// Temperature-induced martensite fraction ξ_T ∈ \[0, 1\].
    pub xi_t: f64,
    /// Underlying SMA constitutive model.
    pub sma: ShapeMemoryAlloy,
    /// Current applied strain.
    pub strain: f64,
    /// Current stress \[Pa\].
    pub stress: f64,
}
impl SmaActuator {
    /// Create an SMA wire actuator from SMA properties and geometry.
    pub fn new(sma: ShapeMemoryAlloy, area: f64, length: f64) -> Self {
        Self {
            area,
            length,
            xi_s: 0.0,
            xi_t: sma.xi,
            sma,
            strain: 0.0,
            stress: 0.0,
        }
    }
    /// Total martensite fraction ξ = ξ_s + ξ_T.
    pub fn total_xi(&self) -> f64 {
        (self.xi_s + self.xi_t).clamp(0.0, 1.0)
    }
    /// Recovery stress for current strain and temperature \[Pa\].
    ///
    /// Uses Brinson: σ = E(ξ) * (ε - ε_L * ξ_s) where ε_L = max_strain.
    pub fn recovery_stress(&self) -> f64 {
        let e = self.sma.e_a + self.total_xi() * (self.sma.e_m - self.sma.e_a);
        e * (self.strain - self.sma.max_strain * self.xi_s)
    }
    /// Stroke (shortening) of the actuator wire \[m\].
    ///
    /// Full stroke = max_strain * length when going from full martensite to austenite.
    pub fn stroke(&self, xi_final: f64) -> f64 {
        let delta_xi = self.total_xi() - xi_final.clamp(0.0, 1.0);
        delta_xi * self.sma.max_strain * self.length
    }
    /// Force output of the actuator \[N\].
    pub fn force(&self) -> f64 {
        self.recovery_stress() * self.area
    }
    /// Update actuator for given temperature \[K\] and external strain.
    pub fn update(&mut self, temperature: f64, applied_strain: f64) {
        self.strain = applied_strain;
        let xi_new = self.sma.update_phase(temperature, self.stress);
        let xi_s_new = (applied_strain / self.sma.max_strain.max(1e-12)).clamp(0.0, xi_new);
        self.xi_s = xi_s_new;
        self.xi_t = (xi_new - xi_s_new).max(0.0);
        self.stress = self.recovery_stress();
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
