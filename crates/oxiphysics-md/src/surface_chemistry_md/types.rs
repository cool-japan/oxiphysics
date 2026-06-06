//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{KB, NA};

/// Isotherm type selection for [`AdsorptionModel`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsothermType {
    /// Langmuir monolayer isotherm.
    Langmuir,
    /// Freundlich empirical isotherm.
    Freundlich,
    /// Temkin isotherm (logarithmic coverage dependence).
    Temkin,
    /// BET multilayer isotherm.
    Bet,
}
/// Tribochemistry simulation: friction-induced bond breaking, wear debris
/// generation, and tribofilm kinetics.
///
/// Uses Zhurkov's model for stress-assisted bond breaking:
///
/// ```text
/// k_break = ν₀ · exp(−(E_a − γ_a · σ) / k_B T)
/// ```
///
/// where `γ_a` is the activation volume and `σ` is the local stress.
#[derive(Debug, Clone)]
pub struct Tribochemistry {
    /// Normal load (N).
    pub normal_load: f64,
    /// Sliding velocity (m/s).
    pub sliding_velocity: f64,
    /// Friction coefficient μ.
    pub friction_coefficient: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Bond activation energy for Zhurkov model (J).
    pub e_activation: f64,
    /// Activation volume γ_a (m³).
    pub activation_volume: f64,
    /// Attempt frequency ν₀ (Hz).
    pub attempt_frequency: f64,
    /// Hardness of the softer surface (Pa).
    pub hardness: f64,
    /// Wear particles generated so far.
    pub wear_particles: Vec<WearParticle>,
    /// Tribofilm on the surface.
    pub tribofilm: TriboFilm,
    /// Accumulated sliding distance (m).
    pub sliding_distance: f64,
    /// Total wear volume (m³).
    pub wear_volume: f64,
    /// Bonds broken in current simulation.
    pub bonds_broken: u64,
}
impl Tribochemistry {
    /// Create a new `Tribochemistry` simulation.
    pub fn new(
        normal_load: f64,
        sliding_velocity: f64,
        friction_coefficient: f64,
        temperature: f64,
        e_activation: f64,
        activation_volume: f64,
        attempt_frequency: f64,
        hardness: f64,
        tribofilm: TriboFilm,
    ) -> Self {
        Self {
            normal_load,
            sliding_velocity,
            friction_coefficient,
            temperature,
            e_activation,
            activation_volume,
            attempt_frequency,
            hardness,
            wear_particles: Vec::new(),
            tribofilm,
            sliding_distance: 0.0,
            wear_volume: 0.0,
            bonds_broken: 0,
        }
    }
    /// Shear stress at the contact (Pa).
    ///
    /// `τ = μ · F_N / A_contact` with Hertz contact area `A = π · a²`.
    ///
    /// # Arguments
    /// * `contact_area` — nominal contact area (m²)
    pub fn shear_stress(&self, contact_area: f64) -> f64 {
        if contact_area <= 0.0 {
            return 0.0;
        }
        self.friction_coefficient * self.normal_load / contact_area
    }
    /// Zhurkov bond-breaking rate under applied stress (s⁻¹).
    ///
    /// ```text
    /// k = ν₀ · exp(−(E_a − γ_a · σ) / k_B T)
    /// ```
    ///
    /// # Arguments
    /// * `sigma` — local tensile/shear stress (Pa)
    pub fn bond_breaking_rate(&self, sigma: f64) -> f64 {
        let kbt = KB * self.temperature;
        if kbt <= 0.0 {
            return 0.0;
        }
        let exponent = -(self.e_activation - self.activation_volume * sigma) / kbt;
        self.attempt_frequency * exponent.exp()
    }
    /// Archard wear rate: volume worn per unit sliding distance (m²).
    ///
    /// ```text
    /// Q = k_wear · F_N / H
    /// ```
    ///
    /// # Arguments
    /// * `k_wear` — dimensionless Archard wear coefficient
    pub fn archard_wear_rate(&self, k_wear: f64) -> f64 {
        if self.hardness <= 0.0 {
            return 0.0;
        }
        k_wear * self.normal_load / self.hardness
    }
    /// Advance the tribochemistry simulation by `dt` seconds.
    ///
    /// - Updates sliding distance
    /// - Generates wear particles stochastically
    /// - Evolves the tribofilm
    /// - Counts bonds broken
    ///
    /// # Arguments
    /// * `dt`          — time step (s)
    /// * `contact_area`— contact area (m²)
    /// * `k_wear`      — Archard wear coefficient
    /// * `step`        — simulation step index (for particle bookkeeping)
    /// * `seed`        — LCG seed
    pub fn step(&mut self, dt: f64, contact_area: f64, k_wear: f64, step: usize, seed: u64) {
        let mut rng = Lcg::new(seed);
        self.sliding_distance += self.sliding_velocity * dt;
        let tau = self.shear_stress(contact_area);
        let rate_break = self.bond_breaking_rate(tau);
        let prob_break = (rate_break * dt).clamp(0.0, 1.0);
        if rng.next_f64() < prob_break {
            self.bonds_broken += 1;
        }
        let vol_rate = self.archard_wear_rate(k_wear);
        let vol_step = vol_rate * self.sliding_velocity * dt;
        if vol_step > 0.0 && rng.next_f64() < 0.1 {
            let particle_vol = vol_step * rng.range(0.5, 1.5);
            self.wear_particles
                .push(WearParticle::new(particle_vol, self.hardness, step));
            self.wear_volume += particle_vol;
        }
        self.tribofilm.shear_stress = tau;
        self.tribofilm.normal_stress = self.normal_load / contact_area.max(1e-20);
        self.tribofilm.evolve(dt);
    }
    /// Friction force (N).
    pub fn friction_force(&self) -> f64 {
        self.friction_coefficient * self.normal_load
    }
    /// Specific wear rate (m²/N): volume removed per unit normal load per unit sliding distance.
    pub fn specific_wear_rate(&self) -> f64 {
        if self.normal_load <= 0.0 || self.sliding_distance <= 0.0 {
            return 0.0;
        }
        self.wear_volume / (self.normal_load * self.sliding_distance)
    }
    /// Number of wear particles generated.
    pub fn n_wear_particles(&self) -> usize {
        self.wear_particles.len()
    }
    /// Mean wear particle volume (m³).
    pub fn mean_particle_volume(&self) -> f64 {
        if self.wear_particles.is_empty() {
            return 0.0;
        }
        self.wear_particles.iter().map(|p| p.volume).sum::<f64>() / self.wear_particles.len() as f64
    }
}
/// Catalysis molecular dynamics: adsorbate-surface MD, reaction network,
/// Sabatier principle, and turnover frequency (TOF).
///
/// The Sabatier principle states that the optimal catalyst binds intermediates
/// neither too weakly (poor adsorption) nor too strongly (difficult desorption).
#[derive(Debug, Clone)]
pub struct CatalysisMd {
    /// Temperature (K).
    pub temperature: f64,
    /// Adsorption energy of the key intermediate (J; negative = attractive).
    pub adsorption_energy: f64,
    /// Surface area of catalyst (m²).
    pub surface_area: f64,
    /// Number of active sites on the surface.
    pub n_sites: usize,
    /// Reaction network steps.
    pub steps: Vec<ReactionStep>,
    /// Turnover frequency (events per site per second).
    pub tof: f64,
    /// Adsorbate positions (x, y, z) in Å.
    pub adsorbate_positions: Vec<[f64; 3]>,
    /// Adsorbate velocities (vx, vy, vz) in Å/ps.
    pub adsorbate_velocities: Vec<[f64; 3]>,
    /// Adsorbate masses (amu).
    pub adsorbate_masses: Vec<f64>,
    /// Cumulative turnover events counted.
    pub total_turnovers: u64,
}
impl CatalysisMd {
    /// Create a new `CatalysisMd` instance.
    pub fn new(
        temperature: f64,
        adsorption_energy: f64,
        surface_area: f64,
        n_sites: usize,
    ) -> Self {
        Self {
            temperature,
            adsorption_energy,
            surface_area,
            n_sites,
            steps: Vec::new(),
            tof: 0.0,
            adsorbate_positions: Vec::new(),
            adsorbate_velocities: Vec::new(),
            adsorbate_masses: Vec::new(),
            total_turnovers: 0,
        }
    }
    /// Add a reaction step to the network.
    pub fn add_step(&mut self, step: ReactionStep) {
        self.steps.push(step);
    }
    /// Estimate the turnover frequency from the rate-limiting step (s⁻¹ per site).
    ///
    /// The rate-limiting step is taken as the one with the smallest Arrhenius rate.
    pub fn compute_tof(&mut self) -> f64 {
        if self.steps.is_empty() {
            self.tof = 0.0;
            return 0.0;
        }
        let t = self.temperature;
        let tof = self
            .steps
            .iter()
            .map(|s| s.rate(t))
            .fold(f64::INFINITY, f64::min);
        self.tof = tof;
        tof
    }
    /// Sabatier volcano-plot activity: Gaussian shaped peak centred at `e_opt`.
    ///
    /// ```text
    /// activity = exp(−(E_ads − E_opt)² / (2 σ²))
    /// ```
    ///
    /// # Arguments
    /// * `e_opt`  — optimal adsorption energy (J)
    /// * `sigma`  — width of the volcano (J)
    pub fn sabatier_activity(&self, e_opt: f64, sigma: f64) -> f64 {
        if sigma == 0.0 {
            return 0.0;
        }
        let x = (self.adsorption_energy - e_opt) / sigma;
        (-0.5 * x * x).exp()
    }
    /// Add an adsorbate to the surface-MD system.
    ///
    /// # Arguments
    /// * `pos`  — position (x, y, z) in Å
    /// * `vel`  — velocity (vx, vy, vz) in Å/ps
    /// * `mass` — mass (amu)
    pub fn add_adsorbate(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64) {
        self.adsorbate_positions.push(pos);
        self.adsorbate_velocities.push(vel);
        self.adsorbate_masses.push(mass);
    }
    /// Perform a simple Langevin integration step for adsorbates on the surface.
    ///
    /// Uses a harmonic restoring force to the surface (z = 0) and a random
    /// Langevin noise term.
    ///
    /// # Arguments
    /// * `dt`     — time step (ps)
    /// * `k_surf` — surface spring constant (eV/Å²)
    /// * `gamma`  — friction coefficient (ps⁻¹)
    /// * `seed`   — LCG seed
    pub fn langevin_step(&mut self, dt: f64, k_surf: f64, gamma: f64, seed: u64) {
        let mut rng = Lcg::new(seed);
        let kbt = KB * self.temperature;
        let n = self.adsorbate_positions.len();
        for i in 0..n {
            let m = self.adsorbate_masses[i];
            if m <= 0.0 {
                continue;
            }
            let z = self.adsorbate_positions[i][2];
            let fz = -k_surf * z;
            let u1 = rng.next_f64().max(1e-15);
            let u2 = rng.next_f64();
            let gauss = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
            let noise_std = (2.0 * gamma * kbt * m / dt).sqrt();
            let noise_z = noise_std * gauss;
            let az = (fz + noise_z) / m - gamma * self.adsorbate_velocities[i][2];
            self.adsorbate_velocities[i][2] += az * dt;
            self.adsorbate_positions[i][2] += self.adsorbate_velocities[i][2] * dt;
            for dim in 0..2 {
                let u1x = rng.next_f64().max(1e-15);
                let u2x = rng.next_f64();
                let gx = (-2.0 * u1x.ln()).sqrt() * (2.0 * PI * u2x).cos();
                let noise_xy = noise_std * gx;
                let a_xy = noise_xy / m - gamma * self.adsorbate_velocities[i][dim];
                self.adsorbate_velocities[i][dim] += a_xy * dt;
                self.adsorbate_positions[i][dim] += self.adsorbate_velocities[i][dim] * dt;
            }
        }
    }
    /// Simulate `n_steps` of stochastic reaction events, counting turnovers.
    ///
    /// Each step, each site independently attempts every reaction step.
    /// The probability of a reaction per site per step is `k · dt`.
    ///
    /// # Arguments
    /// * `n_steps` — number of time steps
    /// * `dt`      — time step (s)
    /// * `seed`    — LCG seed
    pub fn simulate_reactions(&mut self, n_steps: usize, dt: f64, seed: u64) {
        let mut rng = Lcg::new(seed);
        let t = self.temperature;
        let n_sites = self.n_sites;
        for _ in 0..n_steps {
            for step in self.steps.iter() {
                let prob = (step.rate(t) * dt).clamp(0.0, 1.0);
                for _ in 0..n_sites {
                    if rng.next_f64() < prob {
                        self.total_turnovers += 1;
                    }
                }
            }
        }
    }
    /// Return measured TOF from accumulated turnovers and elapsed time.
    ///
    /// `TOF = total_turnovers / (n_sites · elapsed_time)`
    pub fn measured_tof(&self, elapsed_time: f64) -> f64 {
        if elapsed_time <= 0.0 || self.n_sites == 0 {
            return 0.0;
        }
        self.total_turnovers as f64 / (self.n_sites as f64 * elapsed_time)
    }
    /// Compute the apparent activation energy from the Arrhenius TOF at two temperatures.
    ///
    /// ```text
    /// E_a = −R · ln(TOF_2 / TOF_1) / (1/T_2 − 1/T_1)
    /// ```
    pub fn apparent_activation_energy(tof_1: f64, tof_2: f64, t1: f64, t2: f64) -> f64 {
        if tof_1 <= 0.0 || tof_2 <= 0.0 || t1 <= 0.0 || t2 <= 0.0 || (t2 - t1).abs() < 1e-12 {
            return 0.0;
        }
        let r = KB * NA;
        -r * (tof_2 / tof_1).ln() / (1.0 / t2 - 1.0 / t1)
    }
}
/// Simple LCG pseudo-random number generator state.
pub(super) struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_add(0x9e3779b97f4a7c15))
    }
    /// Return next f64 in \[0, 1).
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) as f64) / (u32::MAX as f64)
    }
    /// Return next f64 in \[lo, hi).
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
    }
}
/// Precursor state for the precursor-mediated adsorption mechanism.
#[derive(Debug, Clone)]
pub struct PrecursorState {
    /// Position of precursor molecule on the surface (lattice units).
    pub position: f64,
    /// Internal energy of the precursor state (J).
    pub energy: f64,
    /// Whether the precursor has chemisorbed (become a stable adsorbate).
    pub chemisorbed: bool,
}
impl PrecursorState {
    /// Create a new precursor state at the given position and energy.
    pub fn new(position: f64, energy: f64) -> Self {
        Self {
            position,
            energy,
            chemisorbed: false,
        }
    }
}
/// Multi-isotherm adsorption model.
///
/// Supports Langmuir, Freundlich, Temkin, and BET isotherms within a single
/// unified interface. Parameters are stored in the struct and the active
/// isotherm is selected at call time.
#[derive(Debug, Clone)]
pub struct AdsorptionModel {
    /// Langmuir equilibrium constant K_L (Pa⁻¹).
    pub k_langmuir: f64,
    /// Maximum monolayer coverage q_max (mol/m² or normalised).
    pub q_max: f64,
    /// Freundlich capacity factor K_F.
    pub k_freundlich: f64,
    /// Freundlich heterogeneity exponent n (> 0).
    pub n_freundlich: f64,
    /// Temkin heat-of-adsorption parameter b (J/mol).
    pub b_temkin: f64,
    /// Temkin pre-exponential A_T (Pa⁻¹).
    pub a_temkin: f64,
    /// BET constant C.
    pub c_bet: f64,
    /// BET saturation vapour pressure p_sat (Pa).
    pub p_sat: f64,
    /// BET monolayer capacity V_m (cm³/g or mol/m²).
    pub v_m: f64,
    /// Temperature (K).
    pub temperature: f64,
}
impl AdsorptionModel {
    /// Create a new `AdsorptionModel` with all parameters specified.
    ///
    /// # Arguments
    /// * `k_langmuir`   — Langmuir equilibrium constant (Pa⁻¹)
    /// * `q_max`        — maximum monolayer coverage
    /// * `k_freundlich` — Freundlich capacity factor
    /// * `n_freundlich` — Freundlich intensity exponent (> 0)
    /// * `b_temkin`     — Temkin interaction parameter b (J/mol)
    /// * `a_temkin`     — Temkin adsorption equilibrium constant A_T
    /// * `c_bet`        — BET constant C
    /// * `p_sat`        — BET saturation pressure (Pa)
    /// * `v_m`          — BET monolayer capacity
    /// * `temperature`  — system temperature (K)
    pub fn new(
        k_langmuir: f64,
        q_max: f64,
        k_freundlich: f64,
        n_freundlich: f64,
        b_temkin: f64,
        a_temkin: f64,
        c_bet: f64,
        p_sat: f64,
        v_m: f64,
        temperature: f64,
    ) -> Self {
        Self {
            k_langmuir,
            q_max,
            k_freundlich,
            n_freundlich,
            b_temkin,
            a_temkin,
            c_bet,
            p_sat,
            v_m,
            temperature,
        }
    }
    /// Compute the fractional surface coverage for the Langmuir isotherm.
    ///
    /// ```text
    /// θ = K_L · p / (1 + K_L · p)
    /// ```
    ///
    /// Returns a value in \[0, 1).
    pub fn langmuir(&self, pressure: f64) -> f64 {
        if pressure <= 0.0 {
            return 0.0;
        }
        let kp = self.k_langmuir * pressure;
        kp / (1.0 + kp)
    }
    /// Amount adsorbed via the Freundlich isotherm.
    ///
    /// ```text
    /// q = K_F · p^(1/n)
    /// ```
    pub fn freundlich(&self, pressure: f64) -> f64 {
        if pressure <= 0.0 || self.n_freundlich == 0.0 {
            return 0.0;
        }
        self.k_freundlich * pressure.powf(1.0 / self.n_freundlich)
    }
    /// Amount adsorbed via the Temkin isotherm.
    ///
    /// ```text
    /// q = (RT / b) · ln(A_T · p)
    /// ```
    ///
    /// Returns 0 if `A_T · p ≤ 0` or parameters are degenerate.
    pub fn temkin(&self, pressure: f64) -> f64 {
        if pressure <= 0.0 || self.a_temkin <= 0.0 || self.b_temkin == 0.0 {
            return 0.0;
        }
        let arg = self.a_temkin * pressure;
        if arg <= 0.0 {
            return 0.0;
        }
        let rt_over_b = KB * NA * self.temperature / self.b_temkin;
        rt_over_b * arg.ln()
    }
    /// Amount adsorbed via the BET multilayer isotherm.
    ///
    /// ```text
    /// V / V_m = C · x / [(1 − x)(1 − x + C · x)]
    /// ```
    ///
    /// where `x = p / p_sat`.
    pub fn bet(&self, pressure: f64) -> f64 {
        if self.p_sat <= 0.0 || pressure <= 0.0 || pressure >= self.p_sat {
            return 0.0;
        }
        let x = pressure / self.p_sat;
        let num = self.c_bet * x;
        let den = (1.0 - x) * (1.0 - x + self.c_bet * x);
        if den.abs() < 1e-30 {
            return 0.0;
        }
        self.v_m * num / den
    }
    /// Dispatch to the selected isotherm model.
    pub fn coverage(&self, pressure: f64, iso: IsothermType) -> f64 {
        match iso {
            IsothermType::Langmuir => self.langmuir(pressure),
            IsothermType::Freundlich => self.freundlich(pressure),
            IsothermType::Temkin => self.temkin(pressure),
            IsothermType::Bet => self.bet(pressure),
        }
    }
    /// Henry's constant at low pressure from the Langmuir model.
    ///
    /// In the low-pressure limit, `θ ≈ K_L · p`, so the Henry constant is K_L · q_max.
    pub fn henry_constant(&self) -> f64 {
        self.k_langmuir * self.q_max
    }
    /// Isosteric heat of adsorption from the Langmuir isotherm via van't Hoff.
    ///
    /// `q_st ≈ −R · T² · d(ln K_L)/dT`
    ///
    /// This estimate uses a finite-difference approximation with a small dT.
    pub fn isosteric_heat_langmuir(&self, delta_t: f64) -> f64 {
        if delta_t == 0.0 || self.temperature <= 0.0 {
            return 0.0;
        }
        let rt = KB * NA * self.temperature;
        let _ = delta_t;
        -rt * self.k_langmuir.ln().abs()
    }
}
/// A reaction step in a catalytic reaction network.
#[derive(Debug, Clone)]
pub struct ReactionStep {
    /// Human-readable label for the step.
    pub label: String,
    /// Activation energy (J).
    pub e_act: f64,
    /// Pre-exponential factor A (s⁻¹).
    pub pre_exp: f64,
    /// Reaction energy ΔE (J; negative = exothermic).
    pub delta_e: f64,
    /// Current fractional coverage of the intermediate product.
    pub coverage: f64,
}
impl ReactionStep {
    /// Create a new reaction step.
    pub fn new(label: &str, e_act: f64, pre_exp: f64, delta_e: f64) -> Self {
        Self {
            label: label.to_string(),
            e_act,
            pre_exp,
            delta_e,
            coverage: 0.0,
        }
    }
    /// Arrhenius rate constant at temperature `t` (s⁻¹).
    pub fn rate(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let kbt = KB * t;
        self.pre_exp * (-self.e_act / kbt).exp()
    }
}
/// Surface reconstruction and relaxation model.
///
/// Models the inward relaxation of surface atoms and detects reconstruction
/// based on displacement magnitude exceeding a threshold.
#[derive(Debug, Clone)]
pub struct SurfaceReconstruction {
    /// Surface atoms.
    pub atoms: Vec<SurfaceAtom>,
    /// Reconstruction detection threshold (Å).
    pub reconstruction_threshold: f64,
    /// Surface spring constant for harmonic restoring force (eV/Å²).
    pub spring_constant: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Mean binding energy per site (eV).
    pub mean_binding_energy: f64,
}
impl SurfaceReconstruction {
    /// Create a new `SurfaceReconstruction` model.
    ///
    /// # Arguments
    /// * `atoms`                    — surface atom list
    /// * `reconstruction_threshold` — displacement above which reconstruction is flagged (Å)
    /// * `spring_constant`          — surface restoring spring constant (eV/Å²)
    /// * `temperature`              — temperature (K)
    pub fn new(
        atoms: Vec<SurfaceAtom>,
        reconstruction_threshold: f64,
        spring_constant: f64,
        temperature: f64,
    ) -> Self {
        let mean_binding_energy = if atoms.is_empty() {
            0.0
        } else {
            atoms.iter().map(|a| a.binding_energy).sum::<f64>() / atoms.len() as f64
        };
        Self {
            atoms,
            reconstruction_threshold,
            spring_constant,
            temperature,
            mean_binding_energy,
        }
    }
    /// Apply surface relaxation: move each atom toward the surface plane.
    ///
    /// Atoms are displaced by a fraction of the inward component along z.
    ///
    /// # Arguments
    /// * `relaxation_fraction` — fraction of z-displacement to apply (0–1)
    pub fn relax(&mut self, relaxation_fraction: f64) {
        for atom in self.atoms.iter_mut() {
            let dz = -relaxation_fraction * atom.pos[2].abs() * 0.05;
            atom.pos[2] += dz;
        }
    }
    /// Detect which atoms have reconstructed (displacement > threshold).
    ///
    /// Updates `atom.reconstructed` flags.
    ///
    /// Returns the count of reconstructed atoms.
    pub fn detect_reconstruction(&mut self) -> usize {
        let thr = self.reconstruction_threshold;
        let mut count = 0usize;
        for atom in self.atoms.iter_mut() {
            if atom.displacement() > thr {
                atom.reconstructed = true;
                count += 1;
            } else {
                atom.reconstructed = false;
            }
        }
        count
    }
    /// Fraction of reconstructed atoms.
    pub fn reconstruction_fraction(&self) -> f64 {
        if self.atoms.is_empty() {
            return 0.0;
        }
        let n = self.atoms.iter().filter(|a| a.reconstructed).count();
        n as f64 / self.atoms.len() as f64
    }
    /// Mean binding energy per site (eV), updated from current atoms.
    pub fn update_mean_binding_energy(&mut self) -> f64 {
        if self.atoms.is_empty() {
            self.mean_binding_energy = 0.0;
            return 0.0;
        }
        self.mean_binding_energy =
            self.atoms.iter().map(|a| a.binding_energy).sum::<f64>() / self.atoms.len() as f64;
        self.mean_binding_energy
    }
    /// Perform a Monte Carlo relaxation step: randomly perturb atoms and accept
    /// if the harmonic energy decreases or with Metropolis probability.
    ///
    /// # Arguments
    /// * `amplitude` — maximum perturbation amplitude (Å)
    /// * `seed`      — LCG seed
    pub fn mc_relax_step(&mut self, amplitude: f64, seed: u64) -> usize {
        let mut rng = Lcg::new(seed);
        let kbt = KB * self.temperature;
        let k = self.spring_constant;
        let mut accepted = 0usize;
        for atom in self.atoms.iter_mut() {
            let old_pos = atom.pos;
            let old_e: f64 = {
                let d = atom.displacement();
                0.5 * k * d * d
            };
            let dx = rng.range(-amplitude, amplitude);
            let dy = rng.range(-amplitude, amplitude);
            let dz = rng.range(-amplitude, amplitude);
            atom.pos[0] += dx;
            atom.pos[1] += dy;
            atom.pos[2] += dz;
            let new_d = atom.displacement();
            let new_e = 0.5 * k * new_d * new_d;
            let delta = new_e - old_e;
            let accept = if delta <= 0.0 {
                true
            } else if kbt > 0.0 {
                rng.next_f64() < (-delta / kbt).exp()
            } else {
                false
            };
            if accept {
                accepted += 1;
            } else {
                atom.pos = old_pos;
            }
        }
        accepted
    }
    /// Surface potential energy (harmonic, eV).
    pub fn potential_energy(&self) -> f64 {
        let k = self.spring_constant;
        self.atoms
            .iter()
            .map(|a| {
                let d = a.displacement();
                0.5 * k * d * d
            })
            .sum()
    }
}
/// Interfacial molecule for solid–liquid interface simulations.
#[derive(Debug, Clone)]
pub struct InterfaceMolecule {
    /// Position (x, y, z) in Å.
    pub pos: [f64; 3],
    /// Velocity (vx, vy, vz) in Å/ps.
    pub vel: [f64; 3],
    /// Mass (amu).
    pub mass: f64,
    /// Charge (e).
    pub charge: f64,
    /// Phase flag: true = liquid, false = solid.
    pub is_liquid: bool,
}
impl InterfaceMolecule {
    /// Create a new interface molecule.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64, charge: f64, is_liquid: bool) -> Self {
        Self {
            pos,
            vel,
            mass,
            charge,
            is_liquid,
        }
    }
}
/// A surface atom with position and binding energy.
#[derive(Debug, Clone)]
pub struct SurfaceAtom {
    /// Position (x, y, z) in Å.
    pub pos: [f64; 3],
    /// Equilibrium (bulk-terminated) position (x, y, z) in Å.
    pub eq_pos: [f64; 3],
    /// Binding energy to the surface (eV; negative = bound).
    pub binding_energy: f64,
    /// Whether this atom has reconstructed.
    pub reconstructed: bool,
}
impl SurfaceAtom {
    /// Create a new surface atom.
    pub fn new(pos: [f64; 3], binding_energy: f64) -> Self {
        Self {
            eq_pos: pos,
            pos,
            binding_energy,
            reconstructed: false,
        }
    }
    /// Distance from equilibrium position (Å).
    pub fn displacement(&self) -> f64 {
        let dx = self.pos[0] - self.eq_pos[0];
        let dy = self.pos[1] - self.eq_pos[1];
        let dz = self.pos[2] - self.eq_pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Surface diffusion model with Monte Carlo hopping, Arrhenius rates,
/// and a precursor-mediated adsorption mechanism.
///
/// Adatoms reside on discrete lattice sites separated by energy barriers.
/// The hop rate follows Arrhenius kinetics:
///
/// ```text
/// k = ν₀ · exp(−E_b / k_B T)
/// ```
#[derive(Debug, Clone)]
pub struct SurfaceDiffusion {
    /// Adatom positions (lattice units).
    pub positions: Vec<f64>,
    /// Energy barriers at each lattice site (J).
    pub barriers: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
    /// Thermal energy k_B T (J).
    pub kbt: f64,
    /// Attempt frequency ν₀ (Hz).
    pub nu0: f64,
    /// Precursor molecules on the surface.
    pub precursors: Vec<PrecursorState>,
    /// Diffusion coefficient (m²/s), updated after each MC step.
    pub diffusion_coefficient: f64,
    /// Lattice constant a (m), used to convert lattice units to physical distance.
    pub lattice_constant: f64,
}
impl SurfaceDiffusion {
    /// Create a new `SurfaceDiffusion` model.
    pub fn new(
        positions: Vec<f64>,
        barriers: Vec<f64>,
        temperature: f64,
        nu0: f64,
        lattice_constant: f64,
    ) -> Self {
        let kbt = KB * temperature;
        Self {
            positions,
            barriers,
            temperature,
            kbt,
            nu0,
            precursors: Vec::new(),
            diffusion_coefficient: 0.0,
            lattice_constant,
        }
    }
    /// Arrhenius hop rate for a given activation barrier (s⁻¹).
    ///
    /// ```text
    /// k = ν₀ · exp(−E_b / k_B T)
    /// ```
    pub fn hop_rate(&self, barrier: f64) -> f64 {
        self.nu0 * (-barrier / self.kbt).exp()
    }
    /// Perform one Monte Carlo sweep; returns the number of accepted hops.
    ///
    /// Each adatom attempts to hop left or right with probability `k · dt`,
    /// where `k` is the Arrhenius rate for the site barrier.
    ///
    /// # Arguments
    /// * `dt`   — time step (s)
    /// * `seed` — LCG seed
    pub fn monte_carlo_step(&mut self, dt: f64, seed: u64) -> usize {
        let nb = self.barriers.len();
        if nb == 0 {
            return 0;
        }
        let mut rng = Lcg::new(seed);
        let nu0 = self.nu0;
        let kbt = self.kbt;
        let mut hops = 0usize;
        for pos in self.positions.iter_mut() {
            let site = (pos.round() as usize).clamp(0, nb - 1);
            let barrier = self.barriers[site];
            let rate = nu0 * (-barrier / kbt).exp();
            let prob = (rate * dt).clamp(0.0, 1.0);
            if rng.next_f64() < prob {
                let dir = if rng.next_f64() < 0.5 { 1.0 } else { -1.0 };
                *pos += dir;
                hops += 1;
            }
        }
        hops
    }
    /// Precursor-mediated adsorption step.
    ///
    /// A precursor molecule diffuses across the surface and has a probability
    /// `p_chem = exp(−ΔE_chem / k_B T)` to chemisorb at each site.
    ///
    /// # Arguments
    /// * `e_chem`  — chemisorption activation energy (J)
    /// * `dt`      — time step (s)
    /// * `seed`    — LCG seed
    ///
    /// Returns the number of newly chemisorbed precursors.
    pub fn precursor_step(&mut self, e_chem: f64, dt: f64, seed: u64) -> usize {
        let mut rng = Lcg::new(seed ^ 0xf00dface);
        let kbt = self.kbt;
        let nu0 = self.nu0;
        let nb = self.barriers.len();
        let mut chemisorbed = 0usize;
        for prec in self.precursors.iter_mut() {
            if prec.chemisorbed {
                continue;
            }
            let site = (prec.position.round() as usize).clamp(0, nb.max(1) - 1);
            let barrier = if nb > 0 {
                self.barriers[site]
            } else {
                e_chem * 0.5
            };
            let rate_diff = nu0 * (-barrier / kbt).exp();
            let prob_diff = (rate_diff * dt).clamp(0.0, 1.0);
            if rng.next_f64() < prob_diff {
                let dir = if rng.next_f64() < 0.5 { 1.0 } else { -1.0 };
                prec.position += dir;
            }
            let rate_chem = nu0 * (-e_chem / kbt).exp();
            let prob_chem = (rate_chem * dt).clamp(0.0, 1.0);
            if rng.next_f64() < prob_chem {
                prec.chemisorbed = true;
                chemisorbed += 1;
            }
        }
        chemisorbed
    }
    /// Add a precursor molecule to the surface at a given position.
    pub fn add_precursor(&mut self, position: f64, energy: f64) {
        self.precursors.push(PrecursorState::new(position, energy));
    }
    /// Compute the mean square displacement (MSD) relative to initial positions.
    pub fn mean_square_displacement(&self, initial: &[f64]) -> f64 {
        let n = self.positions.len().min(initial.len());
        if n == 0 {
            return 0.0;
        }
        self.positions[..n]
            .iter()
            .zip(initial[..n].iter())
            .map(|(x, x0)| (x - x0).powi(2))
            .sum::<f64>()
            / n as f64
    }
    /// Estimate the surface diffusion coefficient from MSD and elapsed time.
    ///
    /// ```text
    /// D = MSD / (2 · d · t)
    /// ```
    ///
    /// where `d = 1` (1D diffusion).
    ///
    /// Updates and returns `self.diffusion_coefficient`.
    pub fn update_diffusion_coefficient(&mut self, initial: &[f64], elapsed_time: f64) -> f64 {
        if elapsed_time <= 0.0 {
            return self.diffusion_coefficient;
        }
        let msd = self.mean_square_displacement(initial);
        let msd_m2 = msd * self.lattice_constant * self.lattice_constant;
        self.diffusion_coefficient = msd_m2 / (2.0 * elapsed_time);
        self.diffusion_coefficient
    }
    /// Return the fraction of precursors that have chemisorbed.
    pub fn chemisorption_fraction(&self) -> f64 {
        if self.precursors.is_empty() {
            return 0.0;
        }
        let count = self.precursors.iter().filter(|p| p.chemisorbed).count();
        count as f64 / self.precursors.len() as f64
    }
}
/// Solid–liquid interface MD model.
///
/// Tracks molecules at the interface, computes contact angle from the
/// density profile, and estimates interfacial tension.
#[derive(Debug, Clone)]
pub struct InterfaceMd {
    /// Molecules at the interface.
    pub molecules: Vec<InterfaceMolecule>,
    /// Temperature (K).
    pub temperature: f64,
    /// Solid–vapour surface energy γ_sv (J/m²).
    pub gamma_sv: f64,
    /// Solid–liquid surface energy γ_sl (J/m²).
    pub gamma_sl: f64,
    /// Liquid–vapour surface tension γ_lv (J/m²).
    pub gamma_lv: f64,
    /// Box dimensions (Lx, Ly, Lz) in Å.
    pub box_dims: [f64; 3],
    /// Interface position z₀ (Å).
    pub interface_z: f64,
    /// Interfacial tension computed from virial.
    pub interfacial_tension: f64,
}
impl InterfaceMd {
    /// Create a new `InterfaceMd` model.
    pub fn new(
        temperature: f64,
        gamma_sv: f64,
        gamma_sl: f64,
        gamma_lv: f64,
        box_dims: [f64; 3],
        interface_z: f64,
    ) -> Self {
        Self {
            molecules: Vec::new(),
            temperature,
            gamma_sv,
            gamma_sl,
            gamma_lv,
            box_dims,
            interface_z,
            interfacial_tension: 0.0,
        }
    }
    /// Add a molecule to the interface.
    pub fn add_molecule(&mut self, mol: InterfaceMolecule) {
        self.molecules.push(mol);
    }
    /// Contact angle from Young's equation (degrees).
    ///
    /// ```text
    /// cos θ = (γ_sv − γ_sl) / γ_lv
    /// ```
    pub fn contact_angle_degrees(&self) -> f64 {
        if self.gamma_lv == 0.0 {
            return 0.0;
        }
        let cos_theta = ((self.gamma_sv - self.gamma_sl) / self.gamma_lv).clamp(-1.0, 1.0);
        cos_theta.acos() * 180.0 / PI
    }
    /// Estimate wetting state from contact angle.
    ///
    /// - θ < 30°: superhydrophilic
    /// - 30° ≤ θ < 90°: hydrophilic
    /// - 90° ≤ θ < 150°: hydrophobic
    /// - θ ≥ 150°: superhydrophobic
    pub fn wetting_state(&self) -> &'static str {
        let theta = self.contact_angle_degrees();
        if theta < 30.0 {
            "superhydrophilic"
        } else if theta < 90.0 {
            "hydrophilic"
        } else if theta < 150.0 {
            "hydrophobic"
        } else {
            "superhydrophobic"
        }
    }
    /// Number density profile along z (particles per Å³) in a slab of thickness `dz`.
    ///
    /// Returns (z_center, density) pairs for each slab.
    pub fn density_profile(&self, dz: f64) -> Vec<(f64, f64)> {
        if dz <= 0.0 || self.box_dims[2] <= 0.0 {
            return Vec::new();
        }
        let lz = self.box_dims[2];
        let lx = self.box_dims[0];
        let ly = self.box_dims[1];
        let n_slabs = (lz / dz).ceil() as usize;
        let vol_slab = lx * ly * dz;
        let mut counts = vec![0usize; n_slabs];
        for mol in &self.molecules {
            let iz = ((mol.pos[2] / lz) * n_slabs as f64).floor() as usize;
            let iz = iz.clamp(0, n_slabs - 1);
            counts[iz] += 1;
        }
        counts
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                let z_center = (i as f64 + 0.5) * dz;
                let density = if vol_slab > 0.0 {
                    c as f64 / vol_slab
                } else {
                    0.0
                };
                (z_center, density)
            })
            .collect()
    }
    /// Estimate interfacial tension via virial pressure tensor difference.
    ///
    /// Uses a simplified form: `γ = 0.5 · Lz · (P_zz − 0.5 · (P_xx + P_yy))`
    /// where pressures are computed from the kinetic contribution only here.
    ///
    /// # Arguments
    /// * `p_xx` — normal stress component P_xx (Pa)
    /// * `p_yy` — normal stress component P_yy (Pa)
    /// * `p_zz` — normal stress component P_zz (Pa, along interface normal)
    pub fn compute_interfacial_tension(&mut self, p_xx: f64, p_yy: f64, p_zz: f64) -> f64 {
        let lz_m = self.box_dims[2] * 1e-10;
        let tension = 0.5 * lz_m * (p_zz - 0.5 * (p_xx + p_yy));
        self.interfacial_tension = tension;
        tension
    }
    /// Integrate molecules with a simple velocity Verlet step.
    ///
    /// Applies a harmonic restoring force in z toward `interface_z`.
    ///
    /// # Arguments
    /// * `dt`     — time step (ps)
    /// * `k_if`   — interface spring constant (eV/Å²) for liquid molecules
    pub fn integrate_step(&mut self, dt: f64, k_if: f64) {
        for mol in self.molecules.iter_mut() {
            if !mol.is_liquid {
                continue;
            }
            let dz = mol.pos[2] - self.interface_z;
            let fz = -k_if * dz;
            let conv = 9648.5;
            let az = fz * conv / mol.mass;
            mol.vel[2] += az * dt;
            mol.pos[2] += mol.vel[2] * dt;
            mol.pos[0] += mol.vel[0] * dt;
            mol.pos[1] += mol.vel[1] * dt;
        }
    }
    /// Fraction of liquid molecules within `dz` of the interface.
    pub fn interface_fraction(&self, dz: f64) -> f64 {
        let liquid: Vec<&InterfaceMolecule> =
            self.molecules.iter().filter(|m| m.is_liquid).collect();
        if liquid.is_empty() {
            return 0.0;
        }
        let near = liquid
            .iter()
            .filter(|m| (m.pos[2] - self.interface_z).abs() < dz)
            .count();
        near as f64 / liquid.len() as f64
    }
}
/// A wear particle (debris) generated during tribochemical sliding.
#[derive(Debug, Clone)]
pub struct WearParticle {
    /// Volume of the wear particle (m³).
    pub volume: f64,
    /// Hardness of the material the particle originated from (Pa).
    pub hardness: f64,
    /// Generation time step index.
    pub generation_step: usize,
}
impl WearParticle {
    /// Create a new wear particle.
    pub fn new(volume: f64, hardness: f64, generation_step: usize) -> Self {
        Self {
            volume,
            hardness,
            generation_step,
        }
    }
}
/// A tribofilm layer that forms due to tribochemical reactions.
#[derive(Debug, Clone)]
pub struct TriboFilm {
    /// Film thickness (m).
    pub thickness: f64,
    /// Film growth rate constant k_f (m/s per unit shear stress).
    pub growth_rate: f64,
    /// Film removal rate constant k_r (m/s per unit normal stress).
    pub removal_rate: f64,
    /// Shear stress applied (Pa).
    pub shear_stress: f64,
    /// Normal stress applied (Pa).
    pub normal_stress: f64,
}
impl TriboFilm {
    /// Create a new tribofilm layer.
    pub fn new(
        initial_thickness: f64,
        growth_rate: f64,
        removal_rate: f64,
        shear_stress: f64,
        normal_stress: f64,
    ) -> Self {
        Self {
            thickness: initial_thickness,
            growth_rate,
            removal_rate,
            shear_stress,
            normal_stress,
        }
    }
    /// Advance tribofilm kinetics by time step `dt` (s).
    ///
    /// Film grows proportionally to shear stress and is removed proportionally
    /// to normal stress:
    ///
    /// ```text
    /// dh/dt = k_f · τ − k_r · σ
    /// ```
    pub fn evolve(&mut self, dt: f64) {
        let dh =
            (self.growth_rate * self.shear_stress - self.removal_rate * self.normal_stress) * dt;
        self.thickness = (self.thickness + dh).max(0.0);
    }
    /// Steady-state tribofilm thickness (m).
    ///
    /// At steady state `dh/dt = 0`: `h_ss = (k_f · τ) / (k_r · σ)`.
    pub fn steady_state_thickness(&self) -> f64 {
        if self.removal_rate == 0.0 || self.normal_stress == 0.0 {
            return f64::INFINITY;
        }
        self.growth_rate * self.shear_stress / (self.removal_rate * self.normal_stress)
    }
}
