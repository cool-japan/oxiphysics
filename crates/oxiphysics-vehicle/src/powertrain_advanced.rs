// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced powertrain simulation for OxiPhysics vehicle crate.
//!
//! Provides detailed combustion engine, fuel injection, exhaust/turbo, electric
//! motor, battery, hybrid drivetrain, CVT, DCT, and throttle-by-wire models.

use std::f64::consts::PI;

// ─── helpers ────────────────────────────────────────────────────────────────

/// Clamp value to \[lo, hi\].
#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// Linear interpolation of tabulated data (x-sorted).
fn interp1d(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    assert!(xs.len() == ys.len() && !xs.is_empty());
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[xs.len() - 1] {
        return ys[ys.len() - 1];
    }
    for i in 0..xs.len() - 1 {
        if x >= xs[i] && x <= xs[i + 1] {
            let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
            return ys[i] + t * (ys[i + 1] - ys[i]);
        }
    }
    ys[ys.len() - 1]
}

/// Lower heating value of gasoline \[J/kg\].
const GASOLINE_LHV_J_PER_KG: f64 = 43.0e6;

// ─── CombustionEngine ───────────────────────────────────────────────────────

/// Wiebe combustion heat release parameters.
#[derive(Debug, Clone)]
pub struct WiebeParams {
    /// Combustion efficiency factor (typically 0.90–0.99).
    pub efficiency: f64,
    /// Shape parameter m (typically 2–3).
    pub shape: f64,
    /// Crank angle at start of combustion \[deg\].
    pub theta_start: f64,
    /// Combustion duration \[deg\].
    pub delta_theta: f64,
}

impl Default for WiebeParams {
    fn default() -> Self {
        Self {
            efficiency: 0.95,
            shape: 2.5,
            theta_start: -10.0,
            delta_theta: 60.0,
        }
    }
}

/// Advanced combustion engine model with crank-angle-resolved simulation.
#[derive(Debug, Clone)]
pub struct CombustionEngine {
    /// Engine displacement \[m³\].
    pub displacement: f64,
    /// Number of cylinders.
    pub cylinders: u32,
    /// Compression ratio.
    pub compression_ratio: f64,
    /// Maximum RPM.
    pub max_rpm: f64,
    /// Idle RPM.
    pub idle_rpm: f64,
    /// Peak torque \[N·m\].
    pub peak_torque: f64,
    /// RPM at peak torque.
    pub peak_torque_rpm: f64,
    /// Peak power \[W\].
    pub peak_power: f64,
    /// RPM at peak power.
    pub peak_power_rpm: f64,
    /// Wiebe combustion parameters.
    pub wiebe: WiebeParams,
    /// Current RPM (runtime state).
    pub rpm: f64,
    /// Current throttle position \[0..1\] (runtime state).
    pub throttle: f64,
    /// Accumulated crank angle \[deg\] (runtime state).
    pub crank_angle: f64,
    /// Coolant temperature \[°C\] (runtime state).
    pub coolant_temp: f64,
    /// Oil temperature \[°C\] (runtime state).
    pub oil_temp: f64,
    /// Knock intensity \[0..1\] (runtime state).
    pub knock: f64,
}

impl CombustionEngine {
    /// Create a new combustion engine with typical parameters.
    pub fn new(displacement: f64, cylinders: u32, compression_ratio: f64) -> Self {
        Self {
            displacement,
            cylinders,
            compression_ratio,
            max_rpm: 7000.0,
            idle_rpm: 700.0,
            peak_torque: 350.0,
            peak_torque_rpm: 3500.0,
            peak_power: 147_000.0, // ~200 hp
            peak_power_rpm: 6000.0,
            wiebe: WiebeParams::default(),
            rpm: 700.0,
            throttle: 0.0,
            crank_angle: 0.0,
            coolant_temp: 90.0,
            oil_temp: 85.0,
            knock: 0.0,
        }
    }

    /// Volumetric efficiency as a function of RPM and throttle position.
    ///
    /// Returns a value in \[0..1\]. Peak VE ~0.85 at mid-range RPM.
    pub fn volumetric_efficiency(&self, rpm: f64, throttle: f64) -> f64 {
        // Simple bell-curve model peaking at 3500 RPM
        let rpm_factor = {
            let x = (rpm - 3500.0) / 2000.0;
            (-x * x).exp()
        };
        let base = 0.65 + 0.25 * rpm_factor;
        clamp(base * throttle.sqrt().max(0.1), 0.05, 1.0)
    }

    /// Brake mean effective pressure \[Pa\] at given RPM and load.
    ///
    /// BMEP = 2π * torque / displacement_per_cylinder
    pub fn bmep(&self, torque_nm: f64) -> f64 {
        let displacement_per_cyl = self.displacement / self.cylinders as f64;
        2.0 * PI * torque_nm / displacement_per_cyl
    }

    /// Wiebe function: fraction of fuel burned at crank angle theta \[deg\].
    ///
    /// Returns x_b in \[0, 1\].
    pub fn wiebe_burned_fraction(&self, theta: f64) -> f64 {
        let w = &self.wiebe;
        if theta < w.theta_start {
            return 0.0;
        }
        let x = (theta - w.theta_start) / w.delta_theta;
        if x > 1.0 {
            return w.efficiency;
        }
        w.efficiency * (1.0 - (-6.908 * x.powf(w.shape + 1.0)).exp())
    }

    /// Instantaneous heat release rate dQb/dtheta at crank angle \[J/deg\].
    ///
    /// `fuel_energy_j` is the total fuel energy released per cycle.
    pub fn heat_release_rate(&self, theta: f64, fuel_energy_j: f64) -> f64 {
        let w = &self.wiebe;
        if theta < w.theta_start || theta > w.theta_start + w.delta_theta * 1.2 {
            return 0.0;
        }
        let x = ((theta - w.theta_start) / w.delta_theta).clamp(0.0, 1.0);
        let m1 = w.shape + 1.0;
        let val = 6.908 * w.efficiency * m1 / w.delta_theta
            * x.powf(w.shape)
            * (-6.908 * x.powf(m1)).exp();
        fuel_energy_j * val
    }

    /// Engine torque \[N·m\] at given RPM and throttle using empirical curve.
    pub fn torque(&self, rpm: f64, throttle: f64) -> f64 {
        let rpm_norm = rpm / self.peak_torque_rpm;
        // Torque curve: quadratic bell around peak
        let shape = (-((rpm_norm - 1.0) * 2.0).powi(2)).exp();
        let base = self.peak_torque * shape;
        clamp(base * throttle, 0.0, self.peak_torque)
    }

    /// Engine power \[W\] at given RPM and throttle.
    pub fn power(&self, rpm: f64, throttle: f64) -> f64 {
        self.torque(rpm, throttle) * rpm * PI / 30.0
    }

    /// Indicated fuel consumption rate \[kg/s\].
    ///
    /// Uses simplified BSFC map: ~250 g/kWh at best efficiency point.
    pub fn fuel_flow_rate(&self, rpm: f64, throttle: f64) -> f64 {
        let power_kw = self.power(rpm, throttle) / 1000.0;
        let bsfc = self.bsfc(rpm, throttle); // g/kWh
        power_kw * bsfc / (3_600_000.0) // kg/s
    }

    /// Fuel mass flow \[kg/s\] needed to deliver `power_w` of brake power at `rpm`.
    ///
    /// [`CombustionEngine::fuel_flow_rate`] evaluates the engine map at a given
    /// (speed, throttle) point and therefore burns fuel for whatever power the torque
    /// curve yields there.  A powertrain controller works the other way round: it
    /// asks for a power and has to know what that costs.  This returns
    /// `P · BSFC / 3.6·10⁹` with the BSFC read from the map at the relative load
    /// `power_w / peak_power`.  Zero, negative and NaN power cost no fuel.
    pub fn fuel_flow_for_power(&self, rpm: f64, power_w: f64) -> f64 {
        let power_w = power_w.max(0.0);
        let relative_load = clamp(power_w / self.peak_power, 0.0, 1.0);
        power_w / 1000.0 * self.bsfc(rpm, relative_load) / 3_600_000.0
    }

    /// Brake-specific fuel consumption \[g/kWh\].
    pub fn bsfc(&self, rpm: f64, throttle: f64) -> f64 {
        // Island map: best ~240 g/kWh at ~2000 RPM, 80% load
        let rpm_norm = (rpm / 3000.0 - 1.0).powi(2);
        let load_norm = (throttle - 0.8).powi(2) * 4.0;
        240.0 + 80.0 * rpm_norm + 60.0 * load_norm
    }

    /// Knock probability based on compression ratio and RPM.
    pub fn knock_probability(&self, rpm: f64, throttle: f64) -> f64 {
        let base = ((self.compression_ratio - 10.0) / 4.0).max(0.0);
        let rpm_factor = (rpm / self.max_rpm).powi(2);
        clamp(base * rpm_factor * throttle, 0.0, 1.0)
    }

    /// Step the engine simulation by dt \[s\].
    pub fn step(&mut self, dt: f64, load_torque: f64) {
        // Simple angular momentum model: I * d(omega)/dt = T_engine - T_load
        let inertia = 0.3; // kg*m²
        let torque = self.torque(self.rpm, self.throttle);
        let net_torque = torque - load_torque;
        let omega = self.rpm * PI / 30.0;
        let d_omega = net_torque / inertia * dt;
        let new_omega = (omega + d_omega).max(self.idle_rpm * PI / 30.0);
        self.rpm = new_omega * 30.0 / PI;
        self.crank_angle = (self.crank_angle + self.rpm / 60.0 * 360.0 * dt) % 720.0;

        // Thermal model: warm-up
        let target_coolant = 90.0;
        self.coolant_temp += (target_coolant - self.coolant_temp) * 0.001 * dt;
        self.oil_temp += (85.0 - self.oil_temp) * 0.002 * dt;

        // Knock
        self.knock = self.knock_probability(self.rpm, self.throttle);
    }

    /// Set throttle position \[0..1\].
    pub fn set_throttle(&mut self, throttle: f64) {
        self.throttle = clamp(throttle, 0.0, 1.0);
    }
}

// ─── FuelInjection ──────────────────────────────────────────────────────────

/// Fuel injector specification.
#[derive(Debug, Clone)]
pub struct InjectorSpec {
    /// Maximum fuel flow rate \[kg/s\] per injector.
    pub max_flow: f64,
    /// Injection pressure \[Pa\].
    pub rail_pressure: f64,
    /// Number of injectors.
    pub count: u32,
    /// Injector spray angle \[deg\].
    pub spray_angle: f64,
}

/// Fuel injection system model.
#[derive(Debug, Clone)]
pub struct FuelInjection {
    /// Injector specification.
    pub spec: InjectorSpec,
    /// Rail pressure \[Pa\] (can vary with fuel pump speed).
    pub rail_pressure: f64,
    /// Start of injection \[deg BTDC\].
    pub soi: f64,
    /// Current injection duration \[μs\] (runtime).
    pub duration_us: f64,
    /// Lambda (air/fuel ratio normalized to stoichiometric).
    pub lambda: f64,
    /// Total fuel injected this cycle \[kg\] (runtime).
    pub fuel_mass_cycle: f64,
}

impl FuelInjection {
    /// Create a typical port-injection system.
    pub fn new_port_injection(injectors: u32) -> Self {
        Self {
            spec: InjectorSpec {
                max_flow: 3e-4,
                rail_pressure: 400_000.0, // 4 bar
                count: injectors,
                spray_angle: 15.0,
            },
            rail_pressure: 400_000.0,
            soi: 300.0, // 300° BTDC (intake)
            duration_us: 2000.0,
            lambda: 1.0,
            fuel_mass_cycle: 0.0,
        }
    }

    /// Create a GDI (direct injection) system.
    pub fn new_gdi(injectors: u32) -> Self {
        Self {
            spec: InjectorSpec {
                max_flow: 2.5e-4,
                rail_pressure: 20_000_000.0, // 200 bar
                count: injectors,
                spray_angle: 70.0,
            },
            rail_pressure: 20_000_000.0,
            soi: 60.0, // 60° BTDC (compression)
            duration_us: 1500.0,
            lambda: 1.0,
            fuel_mass_cycle: 0.0,
        }
    }

    /// Spray tip penetration \[m\] at time t \[s\] after injection start.
    ///
    /// Naber-Siebers correlation (simplified).
    pub fn spray_penetration(&self, t: f64, ambient_density: f64) -> f64 {
        let delta_p = self.rail_pressure - 101_325.0;
        if delta_p <= 0.0 || t <= 0.0 {
            return 0.0;
        }
        // Breakup time approximation
        let rho_fuel = 750.0;
        let v_inject = (2.0 * delta_p / rho_fuel).sqrt();
        let d_nozzle = 200e-6; // 200 μm nozzle
        let t_break = 40.0 * d_nozzle / v_inject;
        if t < t_break {
            v_inject * t
        } else {
            let coeff = (2.0 * delta_p / ambient_density).powf(0.25);
            coeff * (d_nozzle * t).sqrt()
        }
    }

    /// Compute injection duration \[μs\] for target fuel mass \[kg\].
    pub fn compute_duration(&self, fuel_mass_kg: f64) -> f64 {
        let flow_per_inj =
            self.spec.max_flow * (self.rail_pressure / self.spec.rail_pressure).sqrt();
        let total_flow = flow_per_inj * self.spec.count as f64;
        if total_flow < 1e-12 {
            return 0.0;
        }
        fuel_mass_kg / total_flow * 1e6 // convert to μs
    }

    /// Set target lambda and compute injection duration.
    pub fn set_lambda_and_compute(&mut self, air_mass_kg: f64, target_lambda: f64) {
        let stoich_afr = 14.7;
        let fuel_mass = air_mass_kg / (stoich_afr * target_lambda);
        self.fuel_mass_cycle = fuel_mass;
        self.duration_us = self.compute_duration(fuel_mass);
        self.lambda = target_lambda;
    }

    /// Actual injected fuel mass \[kg\] per cycle (clipped by duration).
    pub fn injected_mass(&self) -> f64 {
        let max_duration_us =
            self.compute_duration(self.spec.max_flow * self.spec.count as f64 * 0.02); // 20 ms max
        let actual_dur = self.duration_us.min(max_duration_us);
        let flow = self.spec.max_flow
            * self.spec.count as f64
            * (self.rail_pressure / self.spec.rail_pressure).sqrt();
        flow * actual_dur * 1e-6
    }
}

// ─── ExhaustSystem ──────────────────────────────────────────────────────────

/// Exhaust system back-pressure model.
#[derive(Debug, Clone)]
pub struct ExhaustSystem {
    /// Pipe diameter \[m\].
    pub pipe_diameter: f64,
    /// Catalytic converter pressure drop coefficient.
    pub cat_dp_coeff: f64,
    /// Muffler pressure drop coefficient.
    pub muffler_dp_coeff: f64,
    /// Back-pressure RPM breakpoints.
    bp_rpms: Vec<f64>,
    /// Back-pressure values \[Pa\] at breakpoints.
    bp_values: Vec<f64>,
}

impl ExhaustSystem {
    /// Create a typical exhaust system.
    pub fn new() -> Self {
        let bp_rpms = vec![0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0];
        let bp_values = vec![
            101325.0, 101600.0, 102000.0, 103000.0, 105000.0, 108000.0, 112000.0, 118000.0,
        ];
        Self {
            pipe_diameter: 0.063,
            cat_dp_coeff: 0.002,
            muffler_dp_coeff: 0.001,
            bp_rpms,
            bp_values,
        }
    }

    /// Back pressure \[Pa\] as a function of engine RPM.
    pub fn back_pressure(&self, rpm: f64) -> f64 {
        interp1d(&self.bp_rpms, &self.bp_values, rpm)
    }

    /// Exhaust mass flow rate \[kg/s\] at given RPM and throttle.
    pub fn mass_flow(&self, rpm: f64, throttle: f64, air_density: f64) -> f64 {
        // Approximate: proportional to RPM and throttle
        let displacement = 2.0e-3; // 2L engine
        let ve = 0.85;
        displacement * rpm / 120.0 * air_density * ve * throttle
    }

    /// Exhaust gas temperature \[°C\] at given load.
    pub fn exhaust_temp(&self, rpm: f64, throttle: f64) -> f64 {
        // Typical EGT: 400-900°C depending on load
        400.0 + 500.0 * throttle * (rpm / 6000.0).sqrt()
    }
}

impl Default for ExhaustSystem {
    fn default() -> Self {
        Self::new()
    }
}

// ─── TurboCharger ───────────────────────────────────────────────────────────

/// Turbocharger compressor map point.
#[derive(Debug, Clone)]
pub struct CompressorMapPoint {
    /// Corrected mass flow \[kg/s\].
    pub mass_flow: f64,
    /// Pressure ratio.
    pub pressure_ratio: f64,
    /// Efficiency.
    pub efficiency: f64,
    /// Corrected speed \[rpm/√K\].
    pub corrected_speed: f64,
}

/// Turbocharger model with compressor/turbine maps and boost control.
#[derive(Debug, Clone)]
pub struct TurboCharger {
    /// Turbine wheel inertia \[kg·m²\].
    pub inertia: f64,
    /// Compressor wheel diameter \[m\].
    pub compressor_diameter: f64,
    /// Turbine wheel diameter \[m\].
    pub turbine_diameter: f64,
    /// Maximum boost pressure \[Pa gauge\].
    pub max_boost: f64,
    /// Wastegate actuator pressure setpoint \[Pa gauge\].
    pub wastegate_setpoint: f64,
    /// Current turbo speed \[rpm\] (runtime).
    pub turbo_rpm: f64,
    /// Current boost pressure \[Pa gauge\] (runtime).
    pub boost_pressure: f64,
    /// Wastegate duty cycle \[0..1\] (runtime).
    pub wastegate_dc: f64,
    /// Turbo lag time constant \[s\].
    pub lag_tau: f64,
}

impl TurboCharger {
    /// Create a typical turbocharger.
    pub fn new() -> Self {
        Self {
            inertia: 1.5e-5,
            compressor_diameter: 0.052,
            turbine_diameter: 0.048,
            max_boost: 120_000.0,         // 1.2 bar gauge
            wastegate_setpoint: 80_000.0, // 0.8 bar
            turbo_rpm: 0.0,
            boost_pressure: 0.0,
            lag_tau: 0.3,
            wastegate_dc: 0.0,
        }
    }

    /// Compressor pressure ratio at given mass flow and speed.
    ///
    /// Simplified map: PR peaks at design point.
    pub fn compressor_pr(&self, mass_flow: f64, speed_rpm: f64) -> f64 {
        let design_flow = 0.15;
        let design_speed = 100_000.0;
        let design_pr = 2.5;
        let flow_factor = 1.0 - ((mass_flow - design_flow) / design_flow).powi(2);
        let speed_factor = (speed_rpm / design_speed).sqrt().min(1.5);
        (1.0 + (design_pr - 1.0) * flow_factor * speed_factor).max(1.0)
    }

    /// Compressor efficiency at given operating point.
    pub fn compressor_efficiency(&self, mass_flow: f64, speed_rpm: f64) -> f64 {
        let design_flow = 0.15;
        let flow_dev = ((mass_flow - design_flow) / design_flow).abs();
        let speed_factor = (speed_rpm / 100_000.0).min(1.0);
        clamp(0.78 * (1.0 - flow_dev) * speed_factor.sqrt(), 0.3, 0.82)
    }

    /// Turbine expansion ratio at given exhaust mass flow.
    pub fn turbine_er(&self, exhaust_flow: f64, speed_rpm: f64) -> f64 {
        let design_flow = 0.18;
        let flow_ratio = exhaust_flow / design_flow;
        let speed_factor = (speed_rpm / 80_000.0).sqrt().min(1.2);
        (1.0 + 1.8 * flow_ratio * speed_factor).min(4.0)
    }

    /// Step turbocharger dynamics by dt \[s\].
    ///
    /// `exhaust_power_w` is turbine input power from exhaust.
    pub fn step(&mut self, dt: f64, exhaust_power_w: f64) {
        // Simple angular momentum on turbo shaft
        let turbo_omega = self.turbo_rpm * PI / 30.0;
        let turbine_torque = if turbo_omega > 1.0 {
            exhaust_power_w / turbo_omega
        } else {
            exhaust_power_w
        };
        let compressor_power = self.boost_pressure * 0.15; // approximate
        let compressor_torque = if turbo_omega > 1.0 {
            compressor_power / turbo_omega
        } else {
            0.0
        };
        let friction_torque = 0.001 * turbo_omega;
        let net_torque = turbine_torque - compressor_torque - friction_torque;
        let d_omega = net_torque / self.inertia * dt;
        let new_omega = (turbo_omega + d_omega).max(0.0);
        self.turbo_rpm = new_omega * 30.0 / PI;

        // Boost pressure with lag
        let target_boost = self.compressor_pr(0.15, self.turbo_rpm) * 101_325.0 - 101_325.0;
        let target_boost = target_boost.min(self.max_boost);
        self.boost_pressure += (target_boost - self.boost_pressure) * dt / self.lag_tau;
        self.boost_pressure = self.boost_pressure.max(0.0);

        // Wastegate control
        self.wastegate_dc = if self.boost_pressure > self.wastegate_setpoint {
            clamp(
                (self.boost_pressure - self.wastegate_setpoint) / 20_000.0,
                0.0,
                1.0,
            )
        } else {
            0.0
        };
    }
}

impl Default for TurboCharger {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Intercooler ────────────────────────────────────────────────────────────

/// Intercooler / charge-air-cooler model.
#[derive(Debug, Clone)]
pub struct Intercooler {
    /// Core volume \[m³\].
    pub volume: f64,
    /// Pressure drop coefficient.
    pub dp_coeff: f64,
    /// Cooling effectiveness \[0..1\].
    pub effectiveness: f64,
    /// Ambient air temperature \[°C\].
    pub ambient_temp: f64,
    /// Current charge air temperature \[°C\] (runtime).
    pub charge_temp: f64,
}

impl Intercooler {
    /// Create a typical intercooler.
    pub fn new() -> Self {
        Self {
            volume: 0.004,
            dp_coeff: 0.03,
            effectiveness: 0.85,
            ambient_temp: 25.0,
            charge_temp: 25.0,
        }
    }

    /// Cooled charge air temperature \[°C\] after intercooler.
    pub fn outlet_temp(&self, inlet_temp: f64) -> f64 {
        inlet_temp - self.effectiveness * (inlet_temp - self.ambient_temp)
    }

    /// Pressure drop \[Pa\] across intercooler.
    pub fn pressure_drop(&self, mass_flow: f64, density: f64) -> f64 {
        if density < 1e-6 {
            return 0.0;
        }
        let velocity = mass_flow / (density * self.volume.cbrt().powi(2));
        0.5 * density * velocity * velocity * self.dp_coeff
    }
}

impl Default for Intercooler {
    fn default() -> Self {
        Self::new()
    }
}

// ─── ElectricMotor ──────────────────────────────────────────────────────────

/// Electric motor type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MotorType {
    /// Permanent magnet synchronous motor.
    Pmsm,
    /// Brushless DC motor.
    Bldc,
    /// Induction motor.
    InductionMotor,
}

/// Electric motor model (PMSM/BLDC).
#[derive(Debug, Clone)]
pub struct ElectricMotor {
    /// Motor type.
    pub motor_type: MotorType,
    /// Peak torque \[N·m\].
    pub peak_torque: f64,
    /// Peak power \[W\].
    pub peak_power: f64,
    /// Maximum motor speed \[RPM\].
    pub max_rpm: f64,
    /// Base speed (end of constant torque region) \[RPM\].
    pub base_rpm: f64,
    /// Motor resistance \[Ω\] (per-phase).
    pub resistance: f64,
    /// Back-EMF constant \[V/(rad/s)\].
    pub ke: f64,
    /// Thermal resistance \[K/W\].
    pub thermal_resistance: f64,
    /// Thermal mass \[J/K\].
    pub thermal_mass: f64,
    /// Motor winding temperature \[°C\] (runtime).
    pub winding_temp: f64,
    /// Current motor speed \[RPM\] (runtime).
    pub rpm: f64,
    /// Current torque command \[N·m\] (runtime).
    pub torque_cmd: f64,
}

impl ElectricMotor {
    /// Create a typical PMSM motor.
    pub fn new_pmsm(peak_torque: f64, peak_power: f64, max_rpm: f64) -> Self {
        let ke = peak_power / (max_rpm * PI / 30.0 * peak_torque).sqrt().max(1.0);
        Self {
            motor_type: MotorType::Pmsm,
            peak_torque,
            peak_power,
            max_rpm,
            base_rpm: max_rpm * 0.3,
            resistance: 0.015,
            ke: ke.max(0.01),
            thermal_resistance: 0.1,
            thermal_mass: 5000.0,
            winding_temp: 25.0,
            rpm: 0.0,
            torque_cmd: 0.0,
        }
    }

    /// Available torque at given speed \[N·m\].
    pub fn available_torque(&self, rpm: f64) -> f64 {
        if rpm <= self.base_rpm {
            self.peak_torque
        } else {
            (self.peak_power / (rpm * PI / 30.0)).min(self.peak_torque)
        }
    }

    /// Motor efficiency at given torque and speed.
    pub fn efficiency(&self, torque: f64, rpm: f64) -> f64 {
        if rpm < 1.0 || torque.abs() < 0.1 {
            return 0.0;
        }
        let omega = rpm * PI / 30.0;
        let mech_power = torque * omega;
        // Copper losses: I²R, where T ≈ ke * I
        let current = torque / self.ke.max(0.001);
        let copper_loss = current * current * self.resistance * 3.0;
        // Iron losses (proportional to speed²)
        let iron_loss = 0.001 * omega * omega;
        let total_input = mech_power + copper_loss + iron_loss;
        if total_input < 1e-6 {
            return 0.0;
        }
        clamp(mech_power / total_input, 0.0, 0.99)
    }

    /// Thermal derating factor \[0..1\] based on winding temperature.
    pub fn derating_factor(&self) -> f64 {
        let max_temp = 150.0;
        let derate_start = 120.0;
        if self.winding_temp < derate_start {
            1.0
        } else {
            clamp(
                1.0 - (self.winding_temp - derate_start) / (max_temp - derate_start),
                0.0,
                1.0,
            )
        }
    }

    /// Actual torque limited by thermal derating.
    pub fn actual_torque(&self, torque_cmd: f64) -> f64 {
        let available = self.available_torque(self.rpm);
        let derated = available * self.derating_factor();
        clamp(torque_cmd, -derated, derated)
    }

    /// Step motor thermal model by dt \[s\].
    pub fn step_thermal(&mut self, dt: f64) {
        let actual_t = self.actual_torque(self.torque_cmd);
        let omega = self.rpm * PI / 30.0;
        let current = actual_t.abs() / self.ke.max(0.001);
        let heat = current * current * self.resistance * 3.0 + 0.001 * omega * omega;
        let cooling = (self.winding_temp - 25.0) / self.thermal_resistance;
        let d_temp = (heat - cooling) / self.thermal_mass * dt;
        self.winding_temp += d_temp;
    }
}

// ─── Battery ────────────────────────────────────────────────────────────────

/// Battery cell chemistry type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellChemistry {
    /// NMC (Nickel Manganese Cobalt).
    Nmc,
    /// LFP (Lithium Iron Phosphate).
    Lfp,
    /// NCA (Nickel Cobalt Aluminum).
    Nca,
}

/// Battery equivalent circuit model (Thevenin 1st order).
#[derive(Debug, Clone)]
pub struct Battery {
    /// Cell chemistry.
    pub chemistry: CellChemistry,
    /// Total capacity \[Ah\].
    pub capacity_ah: f64,
    /// Nominal voltage \[V\].
    pub nominal_voltage: f64,
    /// Internal resistance \[Ω\].
    pub r0: f64,
    /// RC network resistance \[Ω\].
    pub r1: f64,
    /// RC network capacitance \[F\].
    pub c1: f64,
    /// Current SOC \[0..1\] (runtime).
    pub soc: f64,
    /// RC voltage (overpotential) \[V\] (runtime).
    pub v_rc: f64,
    /// Temperature \[°C\] (runtime).
    pub temperature: f64,
    /// Accumulated charge \[Ah\] (runtime).
    pub charge_ah: f64,
}

impl Battery {
    /// Create an NMC battery pack.
    pub fn new_nmc(capacity_ah: f64, nominal_voltage: f64) -> Self {
        Self {
            chemistry: CellChemistry::Nmc,
            capacity_ah,
            nominal_voltage,
            r0: 0.002 * nominal_voltage / (capacity_ah * 3.6),
            r1: 0.001 * nominal_voltage / (capacity_ah * 3.6),
            c1: 2000.0,
            soc: 0.8,
            v_rc: 0.0,
            temperature: 25.0,
            charge_ah: capacity_ah * 0.8,
        }
    }

    /// Open circuit voltage \[V\] as a function of SOC.
    pub fn ocv(&self, soc: f64) -> f64 {
        let soc = clamp(soc, 0.0, 1.0);
        match self.chemistry {
            CellChemistry::Nmc => {
                // NMC OCV curve approximation
                self.nominal_voltage * (0.85 + 0.25 * soc - 0.1 * soc * soc)
            }
            CellChemistry::Lfp => {
                // LFP has a flatter curve
                self.nominal_voltage * (0.9 + 0.1 * soc)
            }
            CellChemistry::Nca => self.nominal_voltage * (0.82 + 0.28 * soc - 0.12 * soc * soc),
        }
    }

    /// Terminal voltage \[V\] under current load I \[A\] (positive = discharge).
    pub fn terminal_voltage(&self, current_a: f64) -> f64 {
        let ocv = self.ocv(self.soc);
        let v_r0 = current_a * self.r0;
        ocv - v_r0 - self.v_rc
    }

    /// Instantaneous power capability \[W\] (positive = discharge).
    pub fn max_discharge_power(&self) -> f64 {
        let ocv = self.ocv(self.soc);
        // Max power at matched load: P_max = V_oc² / (4*R0)
        ocv * ocv / (4.0 * self.r0.max(1e-6))
    }

    /// SOC estimation using Coulomb counting.
    pub fn step(&mut self, current_a: f64, dt: f64) {
        // SOC change
        let d_soc = -current_a * dt / 3600.0 / self.capacity_ah;
        self.soc = clamp(self.soc + d_soc, 0.0, 1.0);
        self.charge_ah = self.soc * self.capacity_ah;

        // RC dynamics
        let tau = self.r1 * self.c1;
        self.v_rc += dt / tau * (current_a * self.r1 - self.v_rc);

        // Simple thermal model
        let heat = current_a * current_a * (self.r0 + self.r1) * dt;
        let mass = self.capacity_ah * 0.5; // rough kg
        let cp = 900.0;
        let d_temp = heat / (mass * cp);
        let cooling = (self.temperature - 25.0) * 0.01 * dt;
        self.temperature += d_temp - cooling;
    }

    /// Estimated remaining range factor \[0..1\] based on SOC.
    pub fn range_factor(&self) -> f64 {
        self.soc
    }
}

// ─── HybridPowertrain ───────────────────────────────────────────────────────

/// Hybrid powertrain mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HybridMode {
    /// Pure electric drive.
    ElectricOnly,
    /// Pure combustion engine.
    EngineOnly,
    /// Parallel hybrid: both engine and motor contribute.
    Parallel,
    /// Series hybrid: engine charges battery, motor drives wheels.
    Series,
    /// Regenerative braking.
    RegenBraking,
    /// Charging stationary.
    Charging,
}

/// Energy management strategy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EmsStrategy {
    /// Rule-based (simple threshold).
    RuleBased,
    /// Charge-depleting/charge-sustaining.
    Cdcs,
    /// ECMS (Equivalent Consumption Minimization Strategy): the engine/motor split
    /// that minimises fuel flow plus the fuel-equivalent of the battery power at
    /// every instant, see [`HybridPowertrain::ecms_split`].
    Ecms,
}

/// Parallel/series hybrid powertrain.
#[derive(Debug, Clone)]
pub struct HybridPowertrain {
    /// Combustion engine.
    pub engine: CombustionEngine,
    /// Electric motor.
    pub motor: ElectricMotor,
    /// Battery.
    pub battery: Battery,
    /// Current hybrid mode (runtime).
    pub mode: HybridMode,
    /// Energy management strategy.
    pub ems: EmsStrategy,
    /// SOC target for charge-sustaining mode; also the SOC at which the ECMS
    /// equivalence factor equals `lambda_ecms`.
    pub soc_target: f64,
    /// ECMS equivalence factor s₀ (dimensionless): joules of fuel energy that one
    /// joule of battery energy is worth when the SOC is at `soc_target`.
    pub lambda_ecms: f64,
    /// Motor speed per unit engine speed (the motor is geared to the crankshaft by
    /// this ratio; 1.0 puts it directly on the crankshaft).
    pub motor_ratio: f64,
    /// Total regenerated energy \[J\] (runtime).
    pub regen_energy_j: f64,
    /// Total fuel consumed \[kg\] (runtime).
    pub total_fuel_kg: f64,
}

impl HybridPowertrain {
    /// Create a new hybrid powertrain.
    pub fn new(engine: CombustionEngine, motor: ElectricMotor, battery: Battery) -> Self {
        Self {
            engine,
            motor,
            battery,
            mode: HybridMode::Parallel,
            ems: EmsStrategy::RuleBased,
            soc_target: 0.6,
            lambda_ecms: 2.5,
            motor_ratio: 1.0,
            regen_energy_j: 0.0,
            total_fuel_kg: 0.0,
        }
    }

    /// Lowest battery SOC at which ECMS still lets the motor discharge the battery.
    const ECMS_SOC_MIN: f64 = 0.2;
    /// Highest battery SOC at which ECMS still lets the motor charge the battery.
    const ECMS_SOC_MAX: f64 = 0.9;
    /// Number of motor-power candidates evaluated by the ECMS sweep.
    const ECMS_SWEEP_POINTS: usize = 41;
    /// Slope of the SOC correction of the ECMS equivalence factor.
    const ECMS_SOC_GAIN: f64 = 1.0;

    /// Determine hybrid mode from driver demand and battery SOC.
    pub fn determine_mode(&self, demanded_power_w: f64, braking: bool) -> HybridMode {
        if braking {
            return HybridMode::RegenBraking;
        }
        match self.ems {
            EmsStrategy::RuleBased => {
                let soc = self.battery.soc;
                if demanded_power_w < 10_000.0 && soc > 0.4 {
                    HybridMode::ElectricOnly
                } else if soc < 0.2 {
                    HybridMode::EngineOnly
                } else {
                    HybridMode::Parallel
                }
            }
            EmsStrategy::Cdcs => {
                if self.battery.soc > self.soc_target {
                    HybridMode::ElectricOnly
                } else {
                    HybridMode::EngineOnly
                }
            }
            EmsStrategy::Ecms => {
                if demanded_power_w < 0.0 {
                    return HybridMode::RegenBraking;
                }
                let (engine_power, motor_power) = self.ecms_split(demanded_power_w);
                if engine_power <= 0.0 {
                    HybridMode::ElectricOnly
                } else if motor_power == 0.0 {
                    HybridMode::EngineOnly
                } else {
                    HybridMode::Parallel
                }
            }
        }
    }

    /// ECMS equivalence factor `s(SOC) = λ · (1 + (SOC_target − SOC))`
    /// (dimensionless, never negative).
    ///
    /// Battery energy is dearer than `lambda_ecms` when the pack is below its target
    /// and cheaper when it is above, which steers the SOC towards `soc_target`.
    fn ecms_equivalence_factor(&self) -> f64 {
        let correction = 1.0 + Self::ECMS_SOC_GAIN * (self.soc_target - self.battery.soc);
        self.lambda_ecms * correction.max(0.0)
    }

    /// Engine/motor power split chosen by the Equivalent Consumption Minimization
    /// Strategy for a wheel power demand `p_dem` \[W\].
    ///
    /// Returns `(engine_power, motor_power)` like [`HybridPowertrain::power_split`];
    /// a negative motor power means that the motor generates (the engine then runs
    /// above the demand and charges the battery, or the wheels recuperate when
    /// `p_dem` is negative).
    ///
    /// The motor power `P_m` is swept over 41 equally spaced points of its feasible
    /// interval and the point with the least instantaneous cost
    ///
    /// `J(P_m) = ṁ_f(P_e) + s(SOC) · P_m / H_u`  with  `P_e = p_dem − P_m`
    ///
    /// is returned.  `ṁ_f` is the fuel flow that delivers the engine power at the
    /// current engine speed ([`CombustionEngine::fuel_flow_for_power`]), `H_u` the lower
    /// heating value of gasoline (43 MJ/kg) and `s` the equivalence factor
    /// `lambda_ecms · (1 + (soc_target − soc))`.  The battery power is the mechanical
    /// motor power, the same lossless convention as [`HybridPowertrain::step`].
    ///
    /// The feasible interval of `P_m` keeps the engine power in
    /// `[0, engine.peak_power]`, the motor power within
    /// `±motor.peak_power · derating_factor()`, forbids discharging at or below an
    /// SOC of 0.2 and charging at or above 0.9.  When the demand exceeds what engine
    /// and motor can deliver together, both run flat out.
    pub fn ecms_split(&self, p_dem: f64) -> (f64, f64) {
        if p_dem.is_nan() {
            return (0.0, 0.0);
        }
        let engine_max = self.engine.peak_power.max(0.0);
        let motor_max = (self.motor.peak_power * self.motor.derating_factor()).max(0.0);
        let soc = self.battery.soc;
        let discharge_max = if soc > Self::ECMS_SOC_MIN {
            motor_max
        } else {
            0.0
        };
        let charge_max = if soc < Self::ECMS_SOC_MAX {
            motor_max
        } else {
            0.0
        };

        if p_dem <= 0.0 {
            // Nothing to propel: the wheels return energy to the motor, as much as
            // the motor and the battery accept.
            return (0.0, p_dem.max(-charge_max));
        }

        let motor_hi = discharge_max.min(p_dem);
        let motor_lo = (-charge_max).max(p_dem - engine_max);
        if motor_lo > motor_hi {
            return (clamp(p_dem - motor_hi, 0.0, engine_max), motor_hi);
        }

        let s = self.ecms_equivalence_factor();
        let cost = |p_m: f64| -> f64 {
            let p_e = clamp(p_dem - p_m, 0.0, engine_max);
            self.engine.fuel_flow_for_power(self.engine.rpm, p_e) + s * p_m / GASOLINE_LHV_J_PER_KG
        };
        let last = (Self::ECMS_SWEEP_POINTS - 1) as f64;
        let mut best_motor = motor_lo;
        let mut best_cost = f64::INFINITY;
        for i in 0..Self::ECMS_SWEEP_POINTS {
            let t = i as f64 / last;
            let p_m = motor_lo * (1.0 - t) + motor_hi * t;
            let c = cost(p_m);
            if c < best_cost {
                best_cost = c;
                best_motor = p_m;
            }
        }
        (clamp(p_dem - best_motor, 0.0, engine_max), best_motor)
    }

    /// Power split between engine and motor \[W each\].
    ///
    /// Returns (engine_power, motor_power).  With the [`EmsStrategy::Ecms`]
    /// strategy the propulsion modes ([`HybridMode::ElectricOnly`],
    /// [`HybridMode::EngineOnly`], [`HybridMode::Parallel`]) and regenerative
    /// braking take the split from [`HybridPowertrain::ecms_split`].
    pub fn power_split(&self, demanded_power_w: f64) -> (f64, f64) {
        if self.ems == EmsStrategy::Ecms {
            match self.mode {
                HybridMode::ElectricOnly | HybridMode::EngineOnly | HybridMode::Parallel => {
                    return self.ecms_split(demanded_power_w);
                }
                HybridMode::RegenBraking => {
                    return self.ecms_split(-demanded_power_w.abs());
                }
                HybridMode::Series | HybridMode::Charging => {}
            }
        }
        match self.mode {
            HybridMode::ElectricOnly => (0.0, demanded_power_w),
            HybridMode::EngineOnly => (demanded_power_w, 0.0),
            HybridMode::Parallel => {
                let motor_max = self.motor.peak_power * self.motor.derating_factor();
                let motor_power = demanded_power_w.min(motor_max) * 0.4;
                let engine_power = demanded_power_w - motor_power;
                (engine_power, motor_power)
            }
            HybridMode::Series => {
                let gen_power = demanded_power_w.min(self.engine.peak_power);
                (gen_power, demanded_power_w)
            }
            HybridMode::RegenBraking => (0.0, -demanded_power_w.abs()),
            HybridMode::Charging => (demanded_power_w.min(self.engine.peak_power), 0.0),
        }
    }

    /// Electrical power drawn from the battery \[W\] (positive = discharge) when the
    /// engine delivers `engine_power_w` and the motor `motor_power_w` in the current
    /// mode.
    ///
    /// The motor always exchanges its power with the battery.  In
    /// [`HybridMode::Series`] the engine drives a generator that feeds the motor
    /// directly, so only the difference reaches the battery; in
    /// [`HybridMode::Charging`] the whole engine power charges it.
    pub fn battery_power(&self, engine_power_w: f64, motor_power_w: f64) -> f64 {
        match self.mode {
            HybridMode::Series => motor_power_w - engine_power_w,
            HybridMode::Charging => -engine_power_w,
            HybridMode::ElectricOnly
            | HybridMode::EngineOnly
            | HybridMode::Parallel
            | HybridMode::RegenBraking => motor_power_w,
        }
    }

    /// Step the hybrid system by dt \[s\].
    pub fn step(&mut self, dt: f64, demanded_power_w: f64, braking: bool) {
        self.mode = self.determine_mode(demanded_power_w, braking);
        self.advance(dt, demanded_power_w);
    }

    /// Advance engine, battery and motor by `dt` seconds in the current mode.
    fn advance(&mut self, dt: f64, demanded_power_w: f64) {
        let (eng_pow, mot_pow) = self.power_split(demanded_power_w);

        // Engine
        let rpm = self.engine.rpm;
        let throttle = if self.engine.peak_power > 0.0 {
            clamp(eng_pow / self.engine.peak_power, 0.0, 1.0)
        } else {
            0.0
        };
        self.engine.set_throttle(throttle);
        self.total_fuel_kg += self.engine.fuel_flow_for_power(rpm, eng_pow) * dt;
        self.engine.step(dt, eng_pow / (rpm * PI / 30.0 + 1.0));

        // Battery: `Battery::step` takes a positive current for discharge.
        let current = if self.battery.nominal_voltage > 0.0 {
            self.battery_power(eng_pow, mot_pow) / self.battery.nominal_voltage
        } else {
            0.0
        };
        self.battery.step(current, dt);
        if mot_pow < 0.0 && (self.mode == HybridMode::RegenBraking || demanded_power_w < 0.0) {
            self.regen_energy_j += mot_pow.abs() * dt;
        }

        // Motor: it turns with the engine, so it has a speed to develop torque at
        // and to dissipate heat at.
        self.motor.rpm = rpm * self.motor_ratio;
        let motor_omega = self.motor.rpm * PI / 30.0;
        self.motor.torque_cmd = if motor_omega.abs() > 1.0 {
            mot_pow / motor_omega
        } else {
            0.0
        };
        self.motor.step_thermal(dt);
    }

    /// Fuel economy \[L/100km\] given average power and speed.
    pub fn fuel_economy(&self, avg_power_w: f64, avg_speed_ms: f64) -> f64 {
        if avg_speed_ms < 0.01 {
            return 0.0;
        }
        let fuel_density = 0.75; // kg/L for gasoline
        let bsfc = self.engine.bsfc(
            self.engine.peak_torque_rpm * 0.6,
            avg_power_w / self.engine.peak_power,
        );
        let fuel_per_kwh = bsfc / 1000.0 / fuel_density; // L/kWh
        let power_kw = avg_power_w / 1000.0;
        let fuel_per_s = power_kw * fuel_per_kwh / 3600.0; // L/s
        fuel_per_s / avg_speed_ms * 100_000.0 // L/100km
    }
}

// ─── CvtTransmission ────────────────────────────────────────────────────────

/// Continuously Variable Transmission (CVT) model.
#[derive(Debug, Clone)]
pub struct CvtTransmission {
    /// Minimum transmission ratio (overdrive).
    pub ratio_min: f64,
    /// Maximum transmission ratio (underdrive).
    pub ratio_max: f64,
    /// Mechanical efficiency at peak \[0..1\].
    pub peak_efficiency: f64,
    /// Belt/chain slip coefficient.
    pub slip_coeff: f64,
    /// Current ratio (runtime).
    pub ratio: f64,
    /// Ratio change rate limit \[1/s\].
    pub ratio_rate_limit: f64,
    /// Primary pulley clamping force \[N\] (runtime).
    pub clamp_force: f64,
}

impl CvtTransmission {
    /// Create a typical CVT.
    pub fn new() -> Self {
        Self {
            ratio_min: 0.4,
            ratio_max: 2.5,
            peak_efficiency: 0.92,
            slip_coeff: 0.002,
            ratio: 2.0,
            ratio_rate_limit: 0.5,
            clamp_force: 10_000.0,
        }
    }

    /// Target ratio for given engine RPM and vehicle speed.
    pub fn target_ratio(&self, engine_rpm: f64, vehicle_speed_ms: f64) -> f64 {
        if vehicle_speed_ms < 0.1 {
            return self.ratio_max;
        }
        let wheel_rpm = vehicle_speed_ms / (2.0 * PI * 0.3) * 60.0;
        let final_drive = 3.5;
        let target = engine_rpm / (wheel_rpm * final_drive);
        clamp(target, self.ratio_min, self.ratio_max)
    }

    /// CVT efficiency as a function of torque and ratio.
    pub fn efficiency(&self, torque_nm: f64, ratio: f64) -> f64 {
        let ratio_norm = (ratio - self.ratio_min) / (self.ratio_max - self.ratio_min);
        // Efficiency drops at extremes and at high torque (slip)
        let ratio_factor = 1.0 - 0.05 * (2.0 * ratio_norm - 1.0).powi(2);
        let torque_factor = 1.0 - self.slip_coeff * (torque_nm / 200.0).abs();
        clamp(
            self.peak_efficiency * ratio_factor * torque_factor,
            0.7,
            0.95,
        )
    }

    /// Output torque considering efficiency and ratio.
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        let eff = self.efficiency(input_torque, self.ratio);
        input_torque * self.ratio * eff
    }

    /// Step CVT ratio toward target by dt \[s\].
    pub fn step(&mut self, dt: f64, target: f64) {
        let target = clamp(target, self.ratio_min, self.ratio_max);
        let delta = target - self.ratio;
        let max_change = self.ratio_rate_limit * dt;
        if delta.abs() > max_change {
            self.ratio += delta.signum() * max_change;
        } else {
            self.ratio = target;
        }
    }

    /// Belt slip velocity \[m/s\].
    pub fn belt_slip(&self, input_torque: f64) -> f64 {
        let clamping_friction = self.clamp_force * 0.1;
        let excess_torque = (input_torque.abs() - clamping_friction * 0.05).max(0.0);
        excess_torque * self.slip_coeff
    }
}

impl Default for CvtTransmission {
    fn default() -> Self {
        Self::new()
    }
}

// ─── DualClutchTransmission ─────────────────────────────────────────────────

/// Clutch engagement state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClutchState {
    /// Fully open (no torque transfer).
    Open,
    /// Slipping (partial engagement).
    Slipping,
    /// Fully locked.
    Locked,
}

/// DCT transmission model with odd/even gear shafts.
#[derive(Debug, Clone)]
pub struct DualClutchTransmission {
    /// Gear ratios for all gears (index 0 = 1st gear).
    pub ratios: Vec<f64>,
    /// Final drive ratio.
    pub final_drive: f64,
    /// Current gear (1-indexed, 0 = neutral).
    pub current_gear: usize,
    /// Pre-selected gear on other shaft.
    pub preselected_gear: usize,
    /// Odd shaft clutch state.
    pub clutch_odd: ClutchState,
    /// Even shaft clutch state.
    pub clutch_even: ClutchState,
    /// Clutch engagement fraction for odd shaft \[0..1\].
    pub engagement_odd: f64,
    /// Clutch engagement fraction for even shaft \[0..1\].
    pub engagement_even: f64,
    /// Shift in progress flag.
    pub shifting: bool,
    /// Shift time elapsed \[s\].
    pub shift_time: f64,
    /// Total shift duration \[s\].
    pub shift_duration: f64,
    /// Mechanical efficiency.
    pub efficiency: f64,
}

impl DualClutchTransmission {
    /// Create a 7-speed DCT.
    pub fn new_7speed() -> Self {
        Self {
            ratios: vec![3.82, 2.36, 1.69, 1.31, 1.00, 0.78, 0.64],
            final_drive: 3.77,
            current_gear: 1,
            preselected_gear: 2,
            clutch_odd: ClutchState::Locked,
            clutch_even: ClutchState::Open,
            engagement_odd: 1.0,
            engagement_even: 0.0,
            shifting: false,
            shift_time: 0.0,
            shift_duration: 0.15,
            efficiency: 0.97,
        }
    }

    /// Current total ratio (including final drive).
    pub fn total_ratio(&self) -> f64 {
        if self.current_gear == 0 || self.current_gear > self.ratios.len() {
            return 1.0;
        }
        self.ratios[self.current_gear - 1] * self.final_drive
    }

    /// Output torque at wheels \[N·m\].
    pub fn output_torque(&self, engine_torque: f64) -> f64 {
        let eng = if self.shifting {
            let shift_engagement = if self.current_gear % 2 == 1 {
                self.engagement_odd
            } else {
                self.engagement_even
            };
            engine_torque * shift_engagement
        } else {
            engine_torque
        };
        eng * self.total_ratio() * self.efficiency
    }

    /// Initiate an upshift.
    pub fn upshift(&mut self) {
        if !self.shifting && self.current_gear < self.ratios.len() {
            self.shifting = true;
            self.shift_time = 0.0;
            self.preselected_gear = self.current_gear + 1;
        }
    }

    /// Initiate a downshift.
    pub fn downshift(&mut self) {
        if !self.shifting && self.current_gear > 1 {
            self.shifting = true;
            self.shift_time = 0.0;
            self.preselected_gear = self.current_gear - 1;
        }
    }

    /// Step DCT state machine by dt \[s\].
    pub fn step(&mut self, dt: f64) {
        if !self.shifting {
            return;
        }
        self.shift_time += dt;
        let progress = (self.shift_time / self.shift_duration).min(1.0);
        // Overlap clutch engagement
        let outgoing_odd = self.current_gear % 2 == 1;
        if outgoing_odd {
            self.engagement_odd = 1.0 - progress;
            self.engagement_even = progress;
        } else {
            self.engagement_even = 1.0 - progress;
            self.engagement_odd = progress;
        }
        if progress >= 1.0 {
            self.current_gear = self.preselected_gear;
            self.shifting = false;
            if self.current_gear % 2 == 1 {
                self.clutch_odd = ClutchState::Locked;
                self.clutch_even = ClutchState::Open;
                self.engagement_odd = 1.0;
                self.engagement_even = 0.0;
            } else {
                self.clutch_even = ClutchState::Locked;
                self.clutch_odd = ClutchState::Open;
                self.engagement_even = 1.0;
                self.engagement_odd = 0.0;
            }
            // Preselect next gear
            if self.current_gear < self.ratios.len() {
                self.preselected_gear = self.current_gear + 1;
            }
        }
    }

    /// Recommended gear for given speed and throttle (simple rule).
    pub fn recommended_gear(&self, vehicle_speed_ms: f64, throttle: f64) -> usize {
        if vehicle_speed_ms < 2.0 {
            return 1;
        }
        let wheel_rpm = vehicle_speed_ms / (2.0 * PI * 0.3) * 60.0;
        // Target engine RPM based on throttle
        let target_rpm = 1500.0 + throttle * 4000.0;
        let mut best_gear = 1usize;
        let mut best_err = f64::MAX;
        for (i, &r) in self.ratios.iter().enumerate() {
            let eng_rpm = wheel_rpm * r * self.final_drive;
            let err = (eng_rpm - target_rpm).abs();
            if err < best_err {
                best_err = err;
                best_gear = i + 1;
            }
        }
        best_gear
    }
}

// ─── ThrottleByWire ─────────────────────────────────────────────────────────

/// Throttle pedal characteristic map.
#[derive(Debug, Clone)]
pub struct PedalMap {
    /// Pedal positions \[0..1\].
    pub pedal_pos: Vec<f64>,
    /// Corresponding throttle openings \[0..1\].
    pub throttle: Vec<f64>,
}

impl PedalMap {
    /// Linear pedal map.
    pub fn linear() -> Self {
        Self {
            pedal_pos: vec![0.0, 0.1, 0.5, 1.0],
            throttle: vec![0.0, 0.05, 0.5, 1.0],
        }
    }

    /// Eco mode pedal map (soft response).
    pub fn eco() -> Self {
        Self {
            pedal_pos: vec![0.0, 0.2, 0.5, 0.8, 1.0],
            throttle: vec![0.0, 0.03, 0.2, 0.6, 1.0],
        }
    }

    /// Sport mode pedal map (sharp initial response).
    pub fn sport() -> Self {
        Self {
            pedal_pos: vec![0.0, 0.1, 0.3, 0.6, 1.0],
            throttle: vec![0.0, 0.2, 0.6, 0.9, 1.0],
        }
    }

    /// Evaluate throttle for given pedal position.
    pub fn evaluate(&self, pedal: f64) -> f64 {
        interp1d(&self.pedal_pos, &self.throttle, pedal)
    }
}

/// Electronic throttle body (fly-by-wire) model.
#[derive(Debug, Clone)]
pub struct ThrottleByWire {
    /// Pedal map.
    pub pedal_map: PedalMap,
    /// Throttle actuator bandwidth \[Hz\].
    pub bandwidth: f64,
    /// Safety limp-home throttle limit \[0..1\].
    pub limp_home_limit: f64,
    /// Fault active (e.g., sensor failure).
    pub fault: bool,
    /// Current throttle position \[0..1\] (runtime).
    pub throttle_pos: f64,
    /// Target throttle from pedal \[0..1\] (runtime).
    pub target_throttle: f64,
    /// TCS override throttle limit \[0..1\] (runtime).
    pub tcs_limit: f64,
}

impl ThrottleByWire {
    /// Create a new throttle-by-wire system.
    pub fn new() -> Self {
        Self {
            pedal_map: PedalMap::linear(),
            bandwidth: 10.0,
            limp_home_limit: 0.15,
            fault: false,
            throttle_pos: 0.0,
            target_throttle: 0.0,
            tcs_limit: 1.0,
        }
    }

    /// Process driver pedal input.
    pub fn process_pedal(&mut self, pedal_pos: f64) {
        self.target_throttle = self.pedal_map.evaluate(clamp(pedal_pos, 0.0, 1.0));
        if self.fault {
            self.target_throttle = self.target_throttle.min(self.limp_home_limit);
        }
        self.target_throttle = self.target_throttle.min(self.tcs_limit);
    }

    /// Step actuator dynamics by dt \[s\].
    pub fn step(&mut self, dt: f64) {
        let tau = 1.0 / (2.0 * PI * self.bandwidth);
        let error = self.target_throttle - self.throttle_pos;
        self.throttle_pos += error * dt / (tau + dt);
        self.throttle_pos = clamp(self.throttle_pos, 0.0, 1.0);
    }

    /// Get actual throttle for engine control.
    pub fn actual_throttle(&self) -> f64 {
        self.throttle_pos
    }

    /// Set TCS override limit.
    pub fn set_tcs_limit(&mut self, limit: f64) {
        self.tcs_limit = clamp(limit, 0.0, 1.0);
    }

    /// Trigger fault condition (e.g., pedal sensor failure).
    pub fn trigger_fault(&mut self) {
        self.fault = true;
        self.target_throttle = self.target_throttle.min(self.limp_home_limit);
    }

    /// Clear fault condition.
    pub fn clear_fault(&mut self) {
        self.fault = false;
    }
}

impl Default for ThrottleByWire {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
