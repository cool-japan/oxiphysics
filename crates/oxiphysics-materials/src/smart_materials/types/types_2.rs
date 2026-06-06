//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types_3::ShapeMemoryAlloy;

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
/// Self-healing material model.
///
/// Models:
/// - Damage-triggered healing kinetics (first-order reaction)
/// - Healing efficiency η_h as a function of time and temperature
/// - Strength recovery: σ_recovery(t) = σ₀ · η_h(t)
/// - Activation energy dependence (Arrhenius kinetics)
#[derive(Debug, Clone)]
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
/// Brinson constitutive model for shape memory alloys.
///
/// Explicitly tracks stress-induced martensite (ξ_s) and
/// temperature-induced martensite (ξ_T) fractions, enabling full
/// simulation of:
/// - Shape Memory Effect (SME): thermally driven full recovery
/// - Superelasticity: stress-driven reversible strain
#[derive(Debug, Clone)]
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
/// Extended hydrogel swelling model using full Flory-Rehner theory.
///
/// Models the equilibrium between:
/// - Mixing free energy (Flory-Huggins): ΔF_mix = kT\[φ ln φ + (1-φ) ln(1-φ) + χφ(1-φ)\]
/// - Elastic free energy (affine network): ΔF_el = (3νkT/2)(λ² − 1 − 2 ln λ)
///
/// where φ is polymer volume fraction and λ = Q^(1/3) is the stretch ratio.
#[derive(Debug, Clone, Copy)]
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
