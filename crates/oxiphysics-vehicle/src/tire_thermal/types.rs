//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
/// Three-layer tire thermal model: tread, carcass, sidewall.
///
/// Heat is generated at the tread from friction and conducted inward through
/// the carcass to the sidewall. Convective cooling acts on the tread (road-side
/// air flow) and the sidewall (inner fender air).
pub struct TireLayerModel {
    /// Model parameters.
    pub params: TireLayerParams,
    /// Current thermal state.
    pub state: TireLayerState,
}
impl TireLayerModel {
    /// Create a new model initialised at ambient temperature.
    pub fn new(params: TireLayerParams) -> Self {
        let t = params.ambient_temperature;
        Self {
            state: TireLayerState {
                temps: [t, t, t],
                heat_generation_rate: 0.0,
            },
            params,
        }
    }
    /// Heat generation from slip: Q = mu * Fn * |v_slip|.
    pub fn compute_heat_generation(
        &mut self,
        normal_force: f64,
        slip_velocity: f64,
        friction_coef: f64,
    ) {
        self.state.heat_generation_rate = friction_coef * normal_force * slip_velocity.abs();
    }
    /// Heat generation from combined slip components.
    ///
    /// Uses the total slip speed: v_slip = sqrt(v_long² + v_lat²).
    pub fn compute_heat_generation_combined(
        &mut self,
        normal_force: f64,
        long_slip_vel: f64,
        lat_slip_vel: f64,
        friction_coef: f64,
    ) {
        let v_slip = (long_slip_vel * long_slip_vel + lat_slip_vel * lat_slip_vel).sqrt();
        self.state.heat_generation_rate = friction_coef * normal_force * v_slip;
    }
    /// Advance the thermal state by `dt` seconds (forward Euler).
    ///
    /// Energy balance per layer:
    /// - Tread: receives friction heat, loses heat to air convection and to carcass.
    /// - Carcass: gains from tread, loses to sidewall.
    /// - Sidewall: gains from carcass, loses to ambient air.
    pub fn step(&mut self, dt: f64) {
        let p = &self.params;
        let t = self.state.temps;
        let q = self.state.heat_generation_rate;
        let t_amb = p.ambient_temperature;
        let k01 = p.inter_layer_conductivity[0] * p.contact_area;
        let k12 = p.inter_layer_conductivity[1] * p.contact_area;
        let cap0 = p.layer_mass[0] * p.layer_cp[0];
        let dt0 =
            (q - p.tread_convection * p.contact_area * (t[0] - t_amb) - k01 * (t[0] - t[1])) / cap0;
        let cap1 = p.layer_mass[1] * p.layer_cp[1];
        let dt1 = (k01 * (t[0] - t[1]) - k12 * (t[1] - t[2])) / cap1;
        let cap2 = p.layer_mass[2] * p.layer_cp[2];
        let dt2 =
            (k12 * (t[1] - t[2]) - p.sidewall_convection * p.sidewall_area * (t[2] - t_amb)) / cap2;
        self.state.temps[0] += dt0 * dt;
        self.state.temps[1] += dt1 * dt;
        self.state.temps[2] += dt2 * dt;
    }
    /// Mean temperature across all three layers (K).
    pub fn mean_temperature(&self) -> f64 {
        (self.state.temps[0] + self.state.temps[1] + self.state.temps[2]) / 3.0
    }
    /// Temperature gradient between tread and sidewall (K).
    pub fn tread_sidewall_gradient(&self) -> f64 {
        self.state.temps[0] - self.state.temps[2]
    }
}
impl TireLayerModel {
    /// Compute heat generation from lateral (cornering) force at the tread.
    ///
    /// `Q_corn = mu_lat * |F_lat| * |v_lat|`
    ///
    /// Stores result in `state.heat_generation_rate`.
    pub fn compute_heat_generation_cornering(
        &mut self,
        lateral_force: f64,
        lat_slip_vel: f64,
        mu_lat: f64,
    ) {
        let q = mu_lat * lateral_force.abs() * lat_slip_vel.abs();
        self.state.heat_generation_rate = q.max(0.0);
    }
    /// Compute the peak-grip temperature range.
    ///
    /// Uses the layer model ambient temperature as the cold baseline.
    /// Peak grip is assumed at `T_amb + 80 K` (typical racing tire).
    /// Returns `(T_optimal_K, T_window_K)`.
    pub fn compute_peak_grip_temperature(&self) -> (f64, f64) {
        let t_opt = self.params.ambient_temperature + 80.0;
        let window = 30.0_f64;
        (t_opt, window)
    }
    /// Compute blister risk for the three-layer model.
    ///
    /// Uses the tread (layer 0) temperature.
    ///
    /// Returns a risk index in `[0, 1]`.
    pub fn compute_blister_risk(&self, t_blister_limit: f64, q_blister_ref: f64) -> f64 {
        let t_tread = self.state.temps[0];
        let q = self.state.heat_generation_rate;
        let t_range = 60.0_f64;
        let excess = (t_tread - t_blister_limit).max(0.0);
        let thermal_risk = (excess / t_range).min(1.0);
        let flux_risk = if q_blister_ref > 1e-9 {
            (q / q_blister_ref).min(1.0)
        } else {
            0.0
        };
        ((thermal_risk + flux_risk) * 0.5).clamp(0.0, 1.0)
    }
}
/// Temperature-dependent grip model using a Gaussian curve.
///
/// mu(T) = mu_cold + (mu_peak - mu_cold) * exp(-0.5 * ((T - T_opt) / sigma)^2)
pub struct TemperatureGrip {
    /// Peak friction coefficient at optimal temperature.
    pub mu_peak: f64,
    /// Friction coefficient at ambient / cold temperature.
    pub mu_cold: f64,
    /// Optimal operating temperature (K).
    pub t_optimal: f64,
    /// Width of the Gaussian window (K).
    pub sigma: f64,
}
impl TemperatureGrip {
    /// Compute grip coefficient at the given temperature.
    pub fn grip_at(&self, temperature: f64) -> f64 {
        let delta = (temperature - self.t_optimal) / self.sigma;
        self.mu_cold + (self.mu_peak - self.mu_cold) * (-0.5 * delta * delta).exp()
    }
    /// Compute the derivative d(mu)/dT at the given temperature.
    pub fn grip_derivative(&self, temperature: f64) -> f64 {
        let delta = (temperature - self.t_optimal) / self.sigma;
        let gauss = (-0.5 * delta * delta).exp();
        -(self.mu_peak - self.mu_cold) * gauss * delta / self.sigma
    }
    /// Whether the temperature is within one sigma of optimal.
    pub fn in_window(&self, temperature: f64) -> bool {
        (temperature - self.t_optimal).abs() <= self.sigma
    }
}
/// Two-node tire thermal model combining grip and temperature dynamics.
pub struct TireThermalModel {
    /// Model parameters.
    pub params: TireThermalParams,
    /// Current thermal state.
    pub state: TireThermalState,
}
impl TireThermalModel {
    /// Create a new model initialised at ambient temperature.
    pub fn new(params: TireThermalParams) -> Self {
        let t_amb = params.ambient_temperature;
        Self {
            state: TireThermalState {
                surface_temp: t_amb,
                core_temp: t_amb,
                road_temp: t_amb,
                heat_generation_rate: 0.0,
            },
            params,
        }
    }
    /// Current grip coefficient based on surface temperature (Gaussian curve).
    pub fn grip_coefficient(&self) -> f64 {
        let p = &self.params;
        let delta = (self.state.surface_temp - p.optimal_temperature) / p.temperature_width;
        p.cold_grip + (p.max_grip - p.cold_grip) * (-0.5 * delta * delta).exp()
    }
    /// Compute frictional heat generation and store in state.
    ///
    /// Q = friction_coef * F_n * |v_slip|
    pub fn compute_heat_generation(
        &mut self,
        normal_force: f64,
        slip_velocity: f64,
        friction_coef: f64,
    ) {
        self.state.heat_generation_rate = friction_coef * normal_force * slip_velocity.abs();
    }
    /// Advance thermal state by `dt` seconds using forward Euler.
    ///
    /// Two-node model:
    /// - dT_surface/dt = (Q/2 - h_conv·A·(T_s - T_amb) - k·(T_s - T_c)) / (m/2·c_p)
    /// - dT_core/dt   = (Q/2 + k·(T_s - T_c)         - k·(T_c - T_amb)) / (m/2·c_p)
    pub fn step(&mut self, dt: f64) {
        let p = &self.params;
        let t_s = self.state.surface_temp;
        let t_c = self.state.core_temp;
        let t_amb = p.ambient_temperature;
        let q = self.state.heat_generation_rate;
        let half_mass = p.mass * 0.5;
        let thermal_cap = half_mass * p.specific_heat;
        let k_tread = p.thermal_conductivity * p.contact_area;
        let k_rim = p.thermal_conductivity * p.contact_area;
        let dt_s = (q * 0.5
            - p.convection_coefficient * p.contact_area * (t_s - t_amb)
            - k_tread * (t_s - t_c))
            / thermal_cap;
        let dt_c = (q * 0.5 + k_tread * (t_s - t_c) - k_rim * (t_c - t_amb)) / thermal_cap;
        self.state.surface_temp += dt_s * dt;
        self.state.core_temp += dt_c * dt;
        self.state.road_temp = 0.5 * (self.state.surface_temp + t_amb);
    }
}
impl TireThermalModel {
    /// Compute heat generation from cornering lateral force.
    ///
    /// During cornering the tire contact patch scrubs laterally, generating
    /// frictional heat proportional to the lateral force magnitude and the
    /// lateral slip velocity.
    ///
    /// `Q_corn = mu_lat * F_lat * |v_lat_slip|`
    ///
    /// where `mu_lat` is the lateral friction coefficient (typically ~0.9–1.0
    /// of the longitudinal value), `F_lat` is the lateral tire force (N), and
    /// `v_lat_slip` is the lateral slip velocity at the contact patch (m/s).
    ///
    /// The computed heat is stored in `state.heat_generation_rate` (W).
    ///
    /// # Arguments
    /// * `lateral_force`  – lateral tire force magnitude (N)
    /// * `lat_slip_vel`   – lateral slip velocity (m/s)
    /// * `mu_lat`         – lateral friction coefficient (dimensionless)
    pub fn compute_heat_generation_cornering(
        &mut self,
        lateral_force: f64,
        lat_slip_vel: f64,
        mu_lat: f64,
    ) {
        let q = mu_lat * lateral_force.abs() * lat_slip_vel.abs();
        self.state.heat_generation_rate = q.max(0.0);
    }
    /// Compute the peak-grip temperature for this tire compound.
    ///
    /// The optimal temperature where the friction coefficient peaks is derived
    /// from the Gaussian grip model parameters.  This method returns the
    /// `optimal_temperature` field directly, but also computes the temperature
    /// bandwidth at 95 % of peak grip so the caller can assess the operating
    /// window.
    ///
    /// Returns `(T_optimal_K, T_window_95pct_K)`.
    ///
    /// The 95 % bandwidth is `2 * sigma * sqrt(2 * ln(20))`, derived from the
    /// Gaussian: `exp(-0.5 * (ΔT/sigma)^2) = 0.95` → `ΔT = sigma * sqrt(2 * 0.051)`.
    pub fn compute_peak_grip_temperature(&self) -> (f64, f64) {
        let p = &self.params;
        let delta = p.temperature_width * (-2.0_f64 * 0.95_f64.ln()).sqrt();
        (p.optimal_temperature, 2.0 * delta)
    }
    /// Compute the blister risk index (dimensionless, 0–1).
    ///
    /// Blistering occurs when the tread surface temperature remains far above
    /// the optimal temperature for an extended period, causing surface
    /// delamination.  This index combines:
    ///
    /// 1. Thermal excess: `excess = max(0, T_surface - T_blister_limit)`
    /// 2. Sustained heat flux: `flux_norm = heat_rate / Q_blister_ref`
    ///
    /// `risk = clamp((excess / T_range + flux_norm) / 2, 0, 1)`
    ///
    /// # Arguments
    /// * `t_blister_limit` – temperature above which blistering may start (K);
    ///   typically `T_optimal + 2 * sigma` (~40–60 K above opt.)
    /// * `q_blister_ref`   – reference heat flux above which blistering is
    ///   likely (W); typical ~8000 W for racing tires.
    ///
    /// Returns a risk index in `[0, 1]` where 1 = imminent blistering.
    pub fn compute_blister_risk(&self, t_blister_limit: f64, q_blister_ref: f64) -> f64 {
        let t_surface = self.state.surface_temp;
        let q = self.state.heat_generation_rate;
        let t_range = 60.0_f64;
        let excess = (t_surface - t_blister_limit).max(0.0);
        let thermal_risk = (excess / t_range).min(1.0);
        let flux_risk = if q_blister_ref > 1e-9 {
            (q / q_blister_ref).min(1.0)
        } else {
            0.0
        };
        ((thermal_risk + flux_risk) * 0.5).clamp(0.0, 1.0)
    }
}
/// Convective cooling model for a surface exposed to airflow.
///
/// Q_conv = h * A * (T_surface - T_air)
pub struct ConvectiveCooling {
    /// Base heat transfer coefficient (W/m²/K).
    pub h_base: f64,
    /// Exposed surface area (m²).
    pub area: f64,
    /// Speed-dependent coefficient: h_eff = h_base + h_speed * sqrt(v_air).
    pub h_speed: f64,
}
impl ConvectiveCooling {
    /// Effective heat transfer coefficient at the given air speed (m/s).
    pub fn effective_h(&self, air_speed: f64) -> f64 {
        self.h_base + self.h_speed * air_speed.abs().sqrt()
    }
    /// Heat flux removed from the surface (W). Positive = heat leaving surface.
    pub fn heat_flux(&self, surface_temp: f64, air_temp: f64, air_speed: f64) -> f64 {
        let h = self.effective_h(air_speed);
        h * self.area * (surface_temp - air_temp)
    }
    /// Apply cooling to a surface for `dt` seconds. Returns new temperature.
    pub fn cool(
        &self,
        surface_temp: f64,
        air_temp: f64,
        air_speed: f64,
        dt: f64,
        thermal_mass: f64,
    ) -> f64 {
        if thermal_mass <= 0.0 {
            return surface_temp;
        }
        let q = self.heat_flux(surface_temp, air_temp, air_speed);
        surface_temp - q * dt / thermal_mass
    }
}
/// Tracks a tire warming up to its optimal operating temperature.
pub struct TireWarmup {
    /// The underlying thermal model.
    pub model: TireThermalModel,
    /// Simulated time elapsed (s).
    pub time_elapsed: f64,
}
impl TireWarmup {
    /// Create a new warmup tracker.
    pub fn new(model: TireThermalModel) -> Self {
        Self {
            model,
            time_elapsed: 0.0,
        }
    }
    /// Returns `true` if the surface temperature is within `temperature_width` of optimal.
    pub fn is_in_window(&self) -> bool {
        let delta = (self.model.state.surface_temp - self.model.params.optimal_temperature).abs();
        delta <= self.model.params.temperature_width
    }
    /// Simulate warmup with constant inputs.
    ///
    /// Returns the time (s) at which the tire reaches optimal temperature,
    /// or `f64::MAX` if it is not reached within `n_steps`.
    pub fn simulate_warmup(&mut self, n_force: f64, v_slip: f64, dt: f64, n_steps: usize) -> f64 {
        for _ in 0..n_steps {
            let mu = self.model.grip_coefficient();
            self.model.compute_heat_generation(n_force, v_slip, mu);
            self.model.step(dt);
            self.time_elapsed += dt;
            let delta =
                (self.model.state.surface_temp - self.model.params.optimal_temperature).abs();
            if delta <= self.model.params.temperature_width {
                return self.time_elapsed;
            }
        }
        f64::MAX
    }
}
/// Cumulative rubber degradation model.
///
/// Tracks thermal stress cycles and computes residual strength fraction.
/// Based on an Arrhenius-style degradation: high temperatures accelerate aging.
#[derive(Debug, Clone)]
pub struct RubberDegradation {
    /// Reference temperature for degradation onset (K).
    pub t_ref: f64,
    /// Arrhenius activation parameter (K).
    pub activation_k: f64,
    /// Cumulative degradation index \[0, 1\]. 0 = new, 1 = fully degraded.
    pub degradation: f64,
}
impl RubberDegradation {
    /// Create a new degradation tracker.
    pub fn new(t_ref: f64, activation_k: f64) -> Self {
        Self {
            t_ref,
            activation_k,
            degradation: 0.0,
        }
    }
    /// Default parameters for a racing compound rubber.
    pub fn racing_compound() -> Self {
        Self {
            t_ref: 383.15,
            activation_k: 6000.0,
            degradation: 0.0,
        }
    }
    /// Instantaneous degradation rate at temperature `t` (K).
    ///
    /// rate = exp((T - T_ref) / activation_k)
    pub fn degradation_rate(&self, t: f64) -> f64 {
        ((t - self.t_ref) / self.activation_k).exp().max(0.0)
    }
    /// Step the model by `dt` at temperature `t`.
    pub fn step(&mut self, t: f64, dt: f64) {
        let rate = self.degradation_rate(t);
        self.degradation = (self.degradation + rate * dt).clamp(0.0, 1.0);
    }
    /// Residual strength fraction (1.0 = new, 0.0 = fully degraded).
    pub fn residual_strength(&self) -> f64 {
        1.0 - self.degradation
    }
    /// Whether the tire has exceeded the thermal degradation threshold.
    pub fn is_thermally_degraded(&self, threshold: f64) -> bool {
        self.degradation >= threshold
    }
}
/// One-dimensional thermal diffusion across N equally spaced nodes.
///
/// Uses forward Euler with constant diffusivity alpha = k / (rho * cp).
pub struct ThermalDiffusion1D {
    /// Thermal diffusivity (m²/s).
    pub alpha: f64,
    /// Node spacing (m).
    pub dx: f64,
    /// Temperatures at each node.
    pub temps: Vec<f64>,
}
impl ThermalDiffusion1D {
    /// Create a new 1-D diffusion model initialised to uniform temperature.
    pub fn new(n_nodes: usize, alpha: f64, dx: f64, initial_temp: f64) -> Self {
        Self {
            alpha,
            dx,
            temps: vec![initial_temp; n_nodes],
        }
    }
    /// Advance by `dt` using forward Euler with Neumann (insulated) boundaries.
    pub fn step(&mut self, dt: f64) {
        let n = self.temps.len();
        if n < 3 {
            return;
        }
        let coeff = self.alpha * dt / (self.dx * self.dx);
        let old = self.temps.clone();
        for i in 1..n - 1 {
            self.temps[i] = old[i] + coeff * (old[i + 1] - 2.0 * old[i] + old[i - 1]);
        }
        self.temps[0] = old[0] + coeff * (old[1] - old[0]);
        self.temps[n - 1] = old[n - 1] + coeff * (old[n - 2] - old[n - 1]);
    }
    /// Mean temperature.
    pub fn mean_temperature(&self) -> f64 {
        if self.temps.is_empty() {
            return 0.0;
        }
        self.temps.iter().sum::<f64>() / self.temps.len() as f64
    }
    /// Maximum temperature.
    pub fn max_temperature(&self) -> f64 {
        self.temps.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
    /// Set heat source at a specific node (adds energy directly).
    pub fn add_heat_at_node(&mut self, node: usize, heat: f64, thermal_mass_per_node: f64) {
        if node < self.temps.len() && thermal_mass_per_node > 0.0 {
            self.temps[node] += heat / thermal_mass_per_node;
        }
    }
}
/// Parameters for the tire thermal model.
pub struct TireThermalParams {
    /// Tire mass (kg).
    pub mass: f64,
    /// Specific heat capacity (J/kg/K), typical ~1000 for rubber.
    pub specific_heat: f64,
    /// Thermal conductivity between surface and core nodes (W/m/K).
    pub thermal_conductivity: f64,
    /// Tire-road contact area (m^2).
    pub contact_area: f64,
    /// Ambient (air) temperature (K).
    pub ambient_temperature: f64,
    /// Convective heat transfer coefficient (W/m^2/K).
    pub convection_coefficient: f64,
    /// Optimal operating temperature where grip peaks (K).
    pub optimal_temperature: f64,
    /// Width of the Gaussian grip curve (K).
    pub temperature_width: f64,
    /// Peak friction coefficient at optimal temperature.
    pub max_grip: f64,
    /// Friction coefficient at ambient temperature (cold tire).
    pub cold_grip: f64,
}
/// Simple one-step-ahead tire temperature predictor.
///
/// Uses the current state and instantaneous heat rate to predict
/// the surface temperature after `dt` seconds.
pub struct TireTemperaturePredictor;
impl TireTemperaturePredictor {
    /// Predict next tread surface temperature (K).
    ///
    /// * `current_temp` — current surface temperature (K).
    /// * `q_in` — heat generation rate (W).
    /// * `q_out` — cooling rate (W).
    /// * `thermal_mass` — effective thermal mass (J/K).
    /// * `dt` — prediction horizon (s).
    pub fn predict_surface_temp(
        current_temp: f64,
        q_in: f64,
        q_out: f64,
        thermal_mass: f64,
        dt: f64,
    ) -> f64 {
        if thermal_mass < 1e-12 {
            return current_temp;
        }
        let d_t = (q_in - q_out) / thermal_mass;
        current_temp + d_t * dt
    }
    /// Predict time to reach target temperature (s) at constant heat rate.
    ///
    /// Returns `f64::INFINITY` if target cannot be reached (net heat rate ≤ 0 when heating).
    pub fn time_to_target(
        current_temp: f64,
        target_temp: f64,
        q_in: f64,
        q_out: f64,
        thermal_mass: f64,
    ) -> f64 {
        let net = q_in - q_out;
        let delta_t = target_temp - current_temp;
        if delta_t.abs() < 1e-6 {
            return 0.0;
        }
        if (delta_t > 0.0 && net <= 0.0) || (delta_t < 0.0 && net >= 0.0) {
            return f64::INFINITY;
        }
        if thermal_mass < 1e-12 {
            return f64::INFINITY;
        }
        (delta_t * thermal_mass) / net
    }
}
/// Three-layer thermal model with layers explicitly named:
/// tread (outer), belt (structural), carcass (inner).
///
/// This is an alias-style wrapper over [`TireLayerModel`] that enforces
/// the physical naming convention described in the tire thermal literature.
#[derive(Debug, Clone)]
pub struct TireBeltModel {
    /// Tread temperature (K) — outermost rubber layer.
    pub tread_temp: f64,
    /// Belt temperature (K) — steel/Kevlar belt pack.
    pub belt_temp: f64,
    /// Carcass temperature (K) — inner liner and structural cords.
    pub carcass_temp: f64,
    /// Thermal mass of tread layer (J/K).
    pub tread_mass_cp: f64,
    /// Thermal mass of belt layer (J/K).
    pub belt_mass_cp: f64,
    /// Thermal mass of carcass layer (J/K).
    pub carcass_mass_cp: f64,
    /// Conductance tread→belt (W/K).
    pub k_tread_belt: f64,
    /// Conductance belt→carcass (W/K).
    pub k_belt_carcass: f64,
    /// Convection coefficient of tread to road/air (W/K).
    pub h_tread: f64,
    /// Convection coefficient of carcass to inner air (W/K).
    pub h_carcass: f64,
    /// Ambient temperature (K).
    pub t_ambient: f64,
}
impl TireBeltModel {
    /// Create a tread/belt/carcass model at ambient temperature.
    pub fn new(
        tread_mass_cp: f64,
        belt_mass_cp: f64,
        carcass_mass_cp: f64,
        k_tread_belt: f64,
        k_belt_carcass: f64,
        h_tread: f64,
        h_carcass: f64,
        t_ambient: f64,
    ) -> Self {
        Self {
            tread_temp: t_ambient,
            belt_temp: t_ambient,
            carcass_temp: t_ambient,
            tread_mass_cp,
            belt_mass_cp,
            carcass_mass_cp,
            k_tread_belt,
            k_belt_carcass,
            h_tread,
            h_carcass,
            t_ambient,
        }
    }
    /// Medium-compound racing tire defaults.
    pub fn medium_racing() -> Self {
        Self::new(
            3.0 * 1100.0,
            4.0 * 900.0,
            2.5 * 950.0,
            0.25 * 0.02,
            0.18 * 0.02,
            30.0 * 0.02,
            8.0 * 0.04,
            293.15,
        )
    }
    /// Advance thermal state by `dt` seconds.
    ///
    /// `q_slip` (W) — frictional heat input at tread surface.
    pub fn step(&mut self, q_slip: f64, dt: f64) {
        let t0 = self.tread_temp;
        let t1 = self.belt_temp;
        let t2 = self.carcass_temp;
        let t_amb = self.t_ambient;
        let dt0 = (q_slip - self.h_tread * (t0 - t_amb) - self.k_tread_belt * (t0 - t1))
            / self.tread_mass_cp;
        let dt1 =
            (self.k_tread_belt * (t0 - t1) - self.k_belt_carcass * (t1 - t2)) / self.belt_mass_cp;
        let dt2 = (self.k_belt_carcass * (t1 - t2) - self.h_carcass * (t2 - t_amb))
            / self.carcass_mass_cp;
        self.tread_temp += dt0 * dt;
        self.belt_temp += dt1 * dt;
        self.carcass_temp += dt2 * dt;
    }
    /// Mean temperature of the three layers (K).
    pub fn mean_temp(&self) -> f64 {
        (self.tread_temp + self.belt_temp + self.carcass_temp) / 3.0
    }
    /// Gradient from tread to carcass (K).
    pub fn gradient(&self) -> f64 {
        self.tread_temp - self.carcass_temp
    }
}
/// Tracks the cumulative thermal energy balance of a tire:
/// energy in from friction, energy out through cooling.
#[derive(Debug, Clone)]
pub struct TireHeatBalance {
    /// Total frictional heat energy generated (J).
    pub total_heat_in: f64,
    /// Total heat energy removed by convection (J).
    pub total_heat_out: f64,
    /// Peak instantaneous heat rate seen (W).
    pub peak_heat_rate: f64,
}
impl TireHeatBalance {
    /// Create a new zeroed heat balance.
    pub fn new() -> Self {
        Self {
            total_heat_in: 0.0,
            total_heat_out: 0.0,
            peak_heat_rate: 0.0,
        }
    }
    /// Record one step: `q_in` (W) heat generated, `q_out` (W) heat removed, `dt` (s).
    pub fn record_step(&mut self, q_in: f64, q_out: f64, dt: f64) {
        let q_in = q_in.max(0.0);
        let q_out = q_out.max(0.0);
        self.total_heat_in += q_in * dt;
        self.total_heat_out += q_out * dt;
        if q_in > self.peak_heat_rate {
            self.peak_heat_rate = q_in;
        }
    }
    /// Net heat stored in the tire (J). Positive = tire is heating up overall.
    pub fn net_heat_stored(&self) -> f64 {
        self.total_heat_in - self.total_heat_out
    }
    /// Cooling efficiency: fraction of generated heat removed.
    pub fn cooling_efficiency(&self) -> f64 {
        if self.total_heat_in < 1e-12 {
            return 0.0;
        }
        (self.total_heat_out / self.total_heat_in).clamp(0.0, 1.0)
    }
}
/// Inner liner (bead area) temperature model for TPMS simulation.
///
/// Estimates the temperature of the inner rubber liner — what a TPMS sensor
/// mounted on the wheel would read. The liner is a few K cooler than the carcass.
#[derive(Debug, Clone)]
pub struct InnerLinerTemperature {
    /// Thermal lag coefficient (1/s): how quickly liner tracks carcass temperature.
    pub lag_coeff: f64,
    /// Current inner liner temperature (K).
    pub temp: f64,
    /// Ambient temperature (K).
    pub t_ambient: f64,
}
impl InnerLinerTemperature {
    /// Create an inner liner model.
    pub fn new(lag_coeff: f64, t_ambient: f64) -> Self {
        Self {
            lag_coeff: lag_coeff.max(0.01),
            temp: t_ambient,
            t_ambient,
        }
    }
    /// Default TPMS sensor model (moderate thermal lag).
    pub fn default_tpms() -> Self {
        Self::new(0.05, 293.15)
    }
    /// Update the inner liner temperature toward the carcass temperature.
    ///
    /// First-order lag: `dT_liner/dt = lag * (T_carcass - T_liner)`
    pub fn update(&mut self, carcass_temp: f64, dt: f64) {
        self.temp += self.lag_coeff * (carcass_temp - self.temp) * dt;
    }
    /// TPMS-reported pressure using Gay-Lussac's law relative to cold state.
    ///
    /// `P_hot = P_cold * T_liner / T_cold`
    pub fn tpms_pressure(&self, p_cold: f64, t_cold: f64) -> f64 {
        if t_cold < 1.0 {
            return p_cold;
        }
        p_cold * self.temp / t_cold
    }
    /// Whether the tire TPMS would flag a pressure warning.
    ///
    /// Warning if pressure exceeds `max_pressure` or drops below `min_pressure`.
    pub fn pressure_warning(&self, p_cold: f64, t_cold: f64, min_p: f64, max_p: f64) -> bool {
        let p = self.tpms_pressure(p_cold, t_cold);
        p < min_p || p > max_p
    }
}
/// Heat exchange model at the tire–road contact patch.
///
/// Separately tracks:
/// 1. Heat generated by friction at the contact patch.
/// 2. Road surface temperature influence on tire tread.
#[derive(Debug, Clone)]
pub struct TireContactHeat {
    /// Fraction of slip heat absorbed by the tire (remainder goes to road).
    pub tire_fraction: f64,
    /// Road surface temperature (K).
    pub road_temp: f64,
    /// Contact conductance (W/K): heat flow between tread and road.
    pub contact_conductance: f64,
}
impl TireContactHeat {
    /// Create a contact heat model.
    pub fn new(tire_fraction: f64, road_temp: f64, contact_conductance: f64) -> Self {
        Self {
            tire_fraction: tire_fraction.clamp(0.0, 1.0),
            road_temp,
            contact_conductance,
        }
    }
    /// Default values for a dry asphalt circuit (road at 40 °C).
    pub fn dry_asphalt() -> Self {
        Self::new(0.6, 313.15, 40.0)
    }
    /// Heat power flowing into the tire tread from slip (W).
    pub fn tire_slip_heat(&self, total_slip_power: f64) -> f64 {
        self.tire_fraction * total_slip_power.max(0.0)
    }
    /// Heat power flowing to the road from slip (W).
    pub fn road_slip_heat(&self, total_slip_power: f64) -> f64 {
        (1.0 - self.tire_fraction) * total_slip_power.max(0.0)
    }
    /// Net tread heat from friction and conduction (W).
    ///
    /// Positive = tread heats up, negative = tread cools (e.g. cool road).
    pub fn net_tread_heat(&self, tread_temp: f64, total_slip_power: f64) -> f64 {
        let friction_in = self.tire_slip_heat(total_slip_power);
        let conduction = self.contact_conductance * (tread_temp - self.road_temp);
        friction_in - conduction
    }
}
/// Models the change in tire inflation pressure with temperature.
///
/// Based on the ideal gas law: P * V = n * R * T.
/// At constant volume: P2 / P1 = T2 / T1 (temperatures in K).
#[derive(Debug, Clone)]
pub struct TireInflationPressure {
    /// Reference (cold) inflation pressure (Pa).
    pub p_cold: f64,
    /// Cold (reference) temperature (K).
    pub t_cold: f64,
    /// Tire inner air volume (m³).
    pub volume: f64,
}
impl TireInflationPressure {
    /// Create a new pressure model.
    pub fn new(p_cold: f64, t_cold: f64, volume: f64) -> Self {
        Self {
            p_cold,
            t_cold,
            volume,
        }
    }
    /// Standard road tire at 2.2 bar cold (293 K, 8 L).
    pub fn road_tire_cold() -> Self {
        Self {
            p_cold: 220_000.0,
            t_cold: 293.15,
            volume: 0.008,
        }
    }
    /// Pressure at the current temperature (Pa) using Gay-Lussac's law.
    pub fn pressure_at(&self, temperature: f64) -> f64 {
        self.p_cold * temperature / self.t_cold.max(1.0)
    }
    /// Pressure rise from cold (Pa) at given temperature.
    pub fn pressure_rise(&self, temperature: f64) -> f64 {
        self.pressure_at(temperature) - self.p_cold
    }
    /// Hot pressure for a typical race operation temperature (373 K ≈ 100 °C).
    pub fn hot_pressure(&self) -> f64 {
        self.pressure_at(373.15)
    }
    /// Whether the tire is over-inflated at the given temperature.
    pub fn is_over_inflated(&self, temperature: f64, max_pressure: f64) -> bool {
        self.pressure_at(temperature) > max_pressure
    }
}
/// Tire compound presets.
pub enum TireCompound {
    /// Soft compound: low optimal temperature, narrow grip window, highest peak grip.
    Soft,
    /// Medium compound: balanced temperature window and grip.
    Medium,
    /// Hard compound: high optimal temperature, wide window, lower peak grip.
    Hard,
}
impl TireCompound {
    /// Return representative `TireThermalParams` for this compound.
    pub fn default_params(&self) -> TireThermalParams {
        let ambient = 293.15;
        match self {
            TireCompound::Soft => TireThermalParams {
                mass: 9.0,
                specific_heat: 1000.0,
                thermal_conductivity: 0.25,
                contact_area: 0.02,
                ambient_temperature: ambient,
                convection_coefficient: 30.0,
                optimal_temperature: 353.15,
                temperature_width: 15.0,
                max_grip: 1.6,
                cold_grip: 0.9,
            },
            TireCompound::Medium => TireThermalParams {
                mass: 9.5,
                specific_heat: 1000.0,
                thermal_conductivity: 0.22,
                contact_area: 0.02,
                ambient_temperature: ambient,
                convection_coefficient: 28.0,
                optimal_temperature: 373.15,
                temperature_width: 25.0,
                max_grip: 1.45,
                cold_grip: 0.85,
            },
            TireCompound::Hard => TireThermalParams {
                mass: 10.0,
                specific_heat: 1000.0,
                thermal_conductivity: 0.20,
                contact_area: 0.02,
                ambient_temperature: ambient,
                convection_coefficient: 25.0,
                optimal_temperature: 393.15,
                temperature_width: 40.0,
                max_grip: 1.3,
                cold_grip: 0.80,
            },
        }
    }
}
/// Heat exchange between the tire tread and the road surface.
///
/// Models two mechanisms:
/// 1. Conduction through the contact patch
/// 2. Frictional heat split between tire and road
#[derive(Debug, Clone)]
pub struct TireGroundHeatExchange {
    /// Thermal conductance of the contact patch (W/K).
    pub contact_conductance: f64,
    /// Fraction of frictional heat absorbed by the tire (0–1).
    /// The remaining fraction goes into the road surface.
    pub tire_heat_fraction: f64,
    /// Road surface temperature (K).
    pub road_temperature: f64,
}
impl TireGroundHeatExchange {
    /// Create a new ground heat exchange model.
    pub fn new(contact_conductance: f64, tire_heat_fraction: f64, road_temperature: f64) -> Self {
        Self {
            contact_conductance,
            tire_heat_fraction,
            road_temperature,
        }
    }
    /// Conductive heat flux from tire tread to road (W). Positive = tire loses heat.
    pub fn conductive_flux(&self, tread_temperature: f64) -> f64 {
        self.contact_conductance * (tread_temperature - self.road_temperature)
    }
    /// Frictional heat absorbed by the tire tread (W) from total friction power.
    pub fn tire_friction_heat(&self, total_friction_power: f64) -> f64 {
        self.tire_heat_fraction * total_friction_power.max(0.0)
    }
    /// Frictional heat deposited into road surface (W).
    pub fn road_friction_heat(&self, total_friction_power: f64) -> f64 {
        (1.0 - self.tire_heat_fraction) * total_friction_power.max(0.0)
    }
    /// Net heat added to tire tread in this step (W): friction heat minus conductive loss.
    pub fn net_tread_heat(&self, tread_temperature: f64, total_friction_power: f64) -> f64 {
        self.tire_friction_heat(total_friction_power) - self.conductive_flux(tread_temperature)
    }
}
/// Circumferential temperature distribution across tire sectors.
pub struct TireTemperatureMap {
    /// Number of sectors around the circumference.
    pub n_sectors: usize,
    /// Temperature of each sector (K).
    pub temps: Vec<f64>,
    /// Index of the sector currently in contact with the road.
    pub contact_sector: usize,
}
impl TireTemperatureMap {
    /// Create a new map initialised to uniform ambient temperature.
    pub fn new(n_sectors: usize, ambient_temperature: f64) -> Self {
        Self {
            n_sectors,
            temps: vec![ambient_temperature; n_sectors],
            contact_sector: 0,
        }
    }
    /// Rotate the contact sector based on angular velocity (rad/s) and timestep.
    pub fn rotate(&mut self, angular_velocity: f64, dt: f64) {
        if self.n_sectors == 0 {
            return;
        }
        let sector_angle = std::f64::consts::TAU / self.n_sectors as f64;
        let angle_turned = angular_velocity.abs() * dt;
        let sectors_moved = (angle_turned / sector_angle).floor() as usize;
        if sectors_moved > 0 {
            self.contact_sector = (self.contact_sector + sectors_moved) % self.n_sectors;
        }
    }
    /// Add frictional heat to the contact sector.
    pub fn heat_contact_sector(&mut self, q: f64, dt: f64, c_p: f64, mass_per_sector: f64) {
        if self.n_sectors == 0 || mass_per_sector <= 0.0 {
            return;
        }
        let delta_t = q * dt / (mass_per_sector * c_p);
        self.temps[self.contact_sector] += delta_t;
    }
    /// Apply convective cooling to all sectors.
    pub fn cool_all(
        &mut self,
        h: f64,
        area: f64,
        t_amb: f64,
        dt: f64,
        c_p: f64,
        mass_per_sector: f64,
    ) {
        if mass_per_sector <= 0.0 {
            return;
        }
        let thermal_cap = mass_per_sector * c_p;
        for t in self.temps.iter_mut() {
            let q_conv = h * area * (*t - t_amb);
            *t -= q_conv * dt / thermal_cap;
        }
    }
    /// Diffuse heat between adjacent sectors (circumferential conduction).
    ///
    /// Uses forward Euler with conductance `k_circ` (W/K) between neighbours.
    pub fn diffuse_circumferential(&mut self, k_circ: f64, dt: f64, thermal_mass_per_sector: f64) {
        let n = self.n_sectors;
        if n < 3 || thermal_mass_per_sector <= 0.0 {
            return;
        }
        let old = self.temps.clone();
        for i in 0..n {
            let prev = if i == 0 { n - 1 } else { i - 1 };
            let next = (i + 1) % n;
            let flux = k_circ * (old[prev] - 2.0 * old[i] + old[next]);
            self.temps[i] += flux * dt / thermal_mass_per_sector;
        }
    }
    /// Mean temperature across all sectors (K).
    pub fn mean_temperature(&self) -> f64 {
        if self.n_sectors == 0 {
            return 0.0;
        }
        self.temps.iter().sum::<f64>() / self.n_sectors as f64
    }
    /// Maximum temperature across all sectors (K).
    pub fn max_temperature(&self) -> f64 {
        self.temps.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }
    /// Temperature variance across sectors.
    pub fn temperature_variance(&self) -> f64 {
        if self.n_sectors == 0 {
            return 0.0;
        }
        let mean = self.mean_temperature();
        self.temps
            .iter()
            .map(|t| (t - mean) * (t - mean))
            .sum::<f64>()
            / self.n_sectors as f64
    }
}
/// Runtime thermal state of a tire (two-node: surface + core).
pub struct TireThermalState {
    /// Tread surface temperature (K).
    pub surface_temp: f64,
    /// Inner core temperature (K).
    pub core_temp: f64,
    /// Road interface temperature (K).
    pub road_temp: f64,
    /// Current heat generation rate (W).
    pub heat_generation_rate: f64,
}
/// Determines whether a tire is operating within its thermal performance window.
#[derive(Debug, Clone)]
pub struct TireOperatingWindow {
    /// Lower temperature bound of optimal window (K).
    pub t_min: f64,
    /// Upper temperature bound of optimal window (K).
    pub t_max: f64,
    /// Critical over-temperature threshold (K). Above this, degradation accelerates.
    pub t_critical: f64,
}
impl TireOperatingWindow {
    /// Create a new operating window.
    pub fn new(t_min: f64, t_max: f64, t_critical: f64) -> Self {
        Self {
            t_min,
            t_max,
            t_critical,
        }
    }
    /// Soft compound window: 70–100 °C optimal, 130 °C critical.
    pub fn soft_compound() -> Self {
        Self {
            t_min: 343.15,
            t_max: 373.15,
            t_critical: 403.15,
        }
    }
    /// Hard compound window: 90–130 °C optimal, 160 °C critical.
    pub fn hard_compound() -> Self {
        Self {
            t_min: 363.15,
            t_max: 403.15,
            t_critical: 433.15,
        }
    }
    /// Returns `true` if temperature is in the optimal window.
    pub fn in_optimal_window(&self, temperature: f64) -> bool {
        temperature >= self.t_min && temperature <= self.t_max
    }
    /// Returns `true` if temperature exceeds the critical threshold.
    pub fn is_critical(&self, temperature: f64) -> bool {
        temperature >= self.t_critical
    }
    /// Temperature margin to optimal window (K).
    /// Positive = needs more heat, negative = overheating (above t_max).
    pub fn warm_up_margin(&self, temperature: f64) -> f64 {
        if temperature < self.t_min {
            temperature - self.t_min
        } else if temperature > self.t_max {
            temperature - self.t_max
        } else {
            0.0
        }
    }
    /// Normalised temperature position within window \[0, 1\].
    /// 0.0 = at t_min, 1.0 = at t_max, clamped outside window.
    pub fn normalised_position(&self, temperature: f64) -> f64 {
        let range = self.t_max - self.t_min;
        if range < 1e-6 {
            return 0.5;
        }
        ((temperature - self.t_min) / range).clamp(0.0, 1.0)
    }
}
/// Runtime thermal state of the three layers.
pub struct TireLayerState {
    /// Temperature of each layer (K): \[tread, carcass, sidewall\].
    pub temps: [f64; 3],
    /// Current heat generation rate applied to the tread (W).
    pub heat_generation_rate: f64,
}
/// Vulcanization-based thermal degradation model.
///
/// Tracks the state of cure `S` of the rubber compound.
/// Under-cure (S < 1): tire is too soft / not fully vulcanized.
/// Over-cure (S >> 1): thermal degradation (reversion).
///
/// Rate follows an Arrhenius model: dS/dt = k₀ · exp(-Ea/(R·T)).
#[derive(Debug, Clone)]
pub struct VulcanizationDegradation {
    /// Pre-exponential frequency factor (1/s).
    pub k0: f64,
    /// Activation energy divided by gas constant (K): Ea/R.
    pub ea_over_r: f64,
    /// Optimal cure state (S = 1.0 → fully cured, no degradation yet).
    pub s_optimal: f64,
    /// Current cure state (dimensionless, starts at 0 for green rubber).
    pub cure_state: f64,
}
impl VulcanizationDegradation {
    /// Create a vulcanization model.
    pub fn new(k0: f64, ea_over_r: f64) -> Self {
        Self {
            k0,
            ea_over_r,
            s_optimal: 1.0,
            cure_state: 0.0,
        }
    }
    /// Typical racing slick compound parameters.
    pub fn racing_slick() -> Self {
        Self::new(1e8, 10000.0)
    }
    /// Instantaneous cure rate at temperature `t` (K).
    pub fn cure_rate(&self, t: f64) -> f64 {
        self.k0 * (-self.ea_over_r / t.max(1.0)).exp()
    }
    /// Advance cure state by `dt` seconds at temperature `t` (K).
    pub fn step(&mut self, t: f64, dt: f64) {
        self.cure_state += self.cure_rate(t) * dt;
    }
    /// Degradation factor \[0, 1\].
    ///
    /// Returns 0 when under-cured, rises above cure_optimal.
    /// Value of 1 indicates severe over-cure (full thermal degradation).
    pub fn degradation_factor(&self) -> f64 {
        if self.cure_state <= self.s_optimal {
            0.0
        } else {
            (1.0 - (-(self.cure_state - self.s_optimal)).exp()).clamp(0.0, 1.0)
        }
    }
    /// Effective grip multiplier accounting for cure state.
    ///
    /// Under-cured: grip is reduced.
    /// Optimal: full grip.
    /// Over-cured: grip decreases due to degradation.
    pub fn grip_multiplier(&self) -> f64 {
        if self.cure_state < self.s_optimal {
            self.cure_state / self.s_optimal
        } else {
            1.0 - self.degradation_factor()
        }
    }
}
/// Parameters for the three-layer tire thermal model (tread, carcass, sidewall).
pub struct TireLayerParams {
    /// Mass of each layer (kg): \[tread, carcass, sidewall\].
    pub layer_mass: [f64; 3],
    /// Specific heat capacity of each layer (J/kg/K).
    pub layer_cp: [f64; 3],
    /// Thermal conductivity between adjacent layers (W/m/K):
    /// \[tread↔carcass, carcass↔sidewall\].
    pub inter_layer_conductivity: [f64; 2],
    /// Contact area used for coupling terms (m²).
    pub contact_area: f64,
    /// Ambient temperature (K).
    pub ambient_temperature: f64,
    /// Convective heat transfer coefficient for tread surface (W/m²/K).
    pub tread_convection: f64,
    /// Convective heat transfer coefficient for sidewall surface (W/m²/K).
    pub sidewall_convection: f64,
    /// Sidewall exposed area (m²).
    pub sidewall_area: f64,
}
impl TireLayerParams {
    /// Default parameters representative of a medium-compound racing tire.
    pub fn default_medium() -> Self {
        Self {
            layer_mass: [3.0, 4.0, 2.5],
            layer_cp: [1100.0, 1000.0, 950.0],
            inter_layer_conductivity: [0.25, 0.18],
            contact_area: 0.02,
            ambient_temperature: 293.15,
            tread_convection: 30.0,
            sidewall_convection: 15.0,
            sidewall_area: 0.06,
        }
    }
}
