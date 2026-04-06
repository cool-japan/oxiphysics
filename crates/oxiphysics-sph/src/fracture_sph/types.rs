//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use std::f64::consts::PI;

use super::functions::Mat3;
#[allow(unused_imports)]
use super::functions::*;

/// Continuum damage model for brittle solids.
///
/// Implements both isotropic damage (scalar D) and anisotropic damage
/// (second-order damage tensor D_ij) with gradient regularization.
#[derive(Debug, Clone)]
pub struct ContinuumDamageModel {
    /// Critical strain ε_c above which damage begins.
    pub critical_strain: f64,
    /// Fracture energy per unit area G_f (J/m²).
    pub fracture_energy: f64,
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Characteristic length l_c for regularization (m).
    pub char_length: f64,
    /// Damage rate exponent m in Grady–Kipp model.
    pub damage_exponent: f64,
    /// Anisotropic damage tensor D_ij (currently stored per-model for reference).
    pub aniso_damage: Mat3,
    /// Model type: true = isotropic, false = anisotropic.
    pub isotropic: bool,
}
impl ContinuumDamageModel {
    /// Create a new isotropic damage model.
    ///
    /// # Arguments
    /// * `critical_strain` — strain threshold for damage onset
    /// * `fracture_energy` — energy release rate G_f (J/m²)
    /// * `youngs_modulus`  — Young's modulus E (Pa)
    /// * `poisson_ratio`   — Poisson's ratio
    /// * `char_length`     — regularization length (m)
    pub fn new_isotropic(
        critical_strain: f64,
        fracture_energy: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
        char_length: f64,
    ) -> Self {
        Self {
            critical_strain,
            fracture_energy,
            youngs_modulus,
            poisson_ratio,
            char_length,
            damage_exponent: 3.0,
            aniso_damage: [[0.0; 3]; 3],
            isotropic: true,
        }
    }
    /// Create a new anisotropic damage model.
    ///
    /// Uses a second-order damage tensor to capture directional degradation.
    pub fn new_anisotropic(
        critical_strain: f64,
        fracture_energy: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
        char_length: f64,
        damage_exponent: f64,
    ) -> Self {
        Self {
            critical_strain,
            fracture_energy,
            youngs_modulus,
            poisson_ratio,
            char_length,
            damage_exponent,
            aniso_damage: [[0.0; 3]; 3],
            isotropic: false,
        }
    }
    /// Compute damage increment dD/dt for isotropic model.
    ///
    /// Uses the Benz–Asphaug (1995) activation function:
    /// `dD/dt = m·D^{(m-1)/m} · ε̇ / ε_g` where ε_g is the activation strain.
    ///
    /// # Arguments
    /// * `current_damage` — current scalar damage D ∈ \[0,1)
    /// * `strain_rate`    — current strain rate ε̇ (1/s)
    /// * `activation_strain` — activation strain for flaw
    pub fn damage_rate_isotropic(
        &self,
        current_damage: f64,
        strain_rate: f64,
        activation_strain: f64,
    ) -> f64 {
        if current_damage >= 1.0 || strain_rate <= 0.0 || activation_strain < 1e-30 {
            return 0.0;
        }
        let m = self.damage_exponent;
        let d_clamped = current_damage.max(1e-10);
        m * d_clamped.powf((m - 1.0) / m) * strain_rate / activation_strain
    }
    /// Update damage for a single particle.
    ///
    /// Returns new damage value clipped to \[0, 1\].
    ///
    /// # Arguments
    /// * `particle`  — current particle state
    /// * `dt`        — timestep (s)
    pub fn update_damage(&self, particle: &SphFractureParticle, dt: f64) -> f64 {
        if particle.damage >= 1.0 {
            return 1.0;
        }
        let sigma_max = particle.max_principal_stress().max(0.0);
        let eff_strain = sigma_max / self.youngs_modulus.max(1e-30);
        if eff_strain < self.critical_strain {
            return particle.damage;
        }
        let strain_rate = (eff_strain - self.critical_strain)
            / (particle.activation_strain.max(1e-30) * dt.max(1e-30));
        let dd =
            self.damage_rate_isotropic(particle.damage, strain_rate, particle.activation_strain)
                * dt;
        (particle.damage + dd).clamp(0.0, 1.0)
    }
    /// Effective secant modulus E_eff = (1 - D) * E.
    pub fn effective_modulus(&self, damage: f64) -> f64 {
        (1.0 - damage.clamp(0.0, 1.0)) * self.youngs_modulus
    }
    /// Regularized damage: smooth using the Gaussian kernel over distance.
    ///
    /// Returns regularized damage at a point given neighbor damages.
    /// The regularization prevents mesh-size dependence.
    ///
    /// # Arguments
    /// * `damages`   — damage values at neighbor positions
    /// * `distances` — distances from the point to neighbors
    pub fn regularized_damage(&self, damages: &[f64], distances: &[f64]) -> f64 {
        let lc = self.char_length.max(1e-30);
        let mut sum_w = 0.0_f64;
        let mut sum_wd = 0.0_f64;
        for (&d, &r) in damages.iter().zip(distances.iter()) {
            let w = (-r * r / (2.0 * lc * lc)).exp();
            sum_w += w;
            sum_wd += w * d;
        }
        if sum_w < 1e-30 { 0.0 } else { sum_wd / sum_w }
    }
    /// Anisotropic damage tensor update.
    ///
    /// Adds damage increment in the crack normal direction n:
    /// `ΔD_ij = (D_new - D_old) * n_i n_j`.
    ///
    /// Returns updated anisotropic damage tensor.
    pub fn update_aniso_damage(
        &self,
        current: &Mat3,
        delta_d: f64,
        crack_normal: [f64; 3],
    ) -> Mat3 {
        let nn = outer3(crack_normal, crack_normal);
        let increment = mat3_scale(&nn, delta_d);
        mat3_add(current, &increment)
    }
    /// Bulk modulus K from E and ν.
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio).max(1e-10))
    }
    /// Shear modulus G from E and ν.
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio).max(1e-10))
    }
    /// Critical strain energy density for fracture initiation.
    ///
    /// `W_c = G_f / l_c`
    pub fn critical_energy_density(&self) -> f64 {
        self.fracture_energy / self.char_length.max(1e-30)
    }
}
/// Grady–Kipp model for strain-rate–dependent impact fragmentation.
///
/// Predicts fragment size distribution from a uniaxial strain pulse.
#[derive(Debug, Clone)]
pub struct ImpactFragmentation {
    /// Grady–Kipp coefficient K_gk (m³·s⁻ⁿ).
    pub k_gk: f64,
    /// Strain-rate exponent n (typically 3 for Grady 1982).
    pub n_exp: f64,
    /// Material density ρ (kg/m³).
    pub density: f64,
    /// Longitudinal wave speed c_l (m/s).
    pub wave_speed: f64,
    /// Weibull modulus m_w for flaw distribution.
    pub weibull_m: f64,
    /// Weibull scale parameter k_w (m⁻³).
    pub weibull_k: f64,
}
impl ImpactFragmentation {
    /// Create a new Grady–Kipp fragmentation model.
    ///
    /// # Arguments
    /// * `k_gk`       — Grady–Kipp coefficient
    /// * `n_exp`      — strain-rate exponent
    /// * `density`    — material density (kg/m³)
    /// * `wave_speed` — longitudinal wave speed (m/s)
    /// * `weibull_m`  — Weibull shape parameter (modulus)
    /// * `weibull_k`  — Weibull scale parameter
    pub fn new(
        k_gk: f64,
        n_exp: f64,
        density: f64,
        wave_speed: f64,
        weibull_m: f64,
        weibull_k: f64,
    ) -> Self {
        Self {
            k_gk,
            n_exp,
            density,
            wave_speed,
            weibull_m,
            weibull_k,
        }
    }
    /// Mean fragment size from Grady (1982) model.
    ///
    /// `s = (K_gk / ε̇ⁿ)^{1/3}` in metres.
    pub fn mean_fragment_size(&self, strain_rate: f64) -> f64 {
        if strain_rate < 1e-30 {
            return f64::INFINITY;
        }
        (self.k_gk / strain_rate.powf(self.n_exp)).cbrt()
    }
    /// Grady energy-based fragment size.
    ///
    /// `s_E = (24 G_f / (ρ ε̇²))^{1/3}` where G_f is fracture energy.
    pub fn grady_energy_fragment_size(&self, strain_rate: f64, fracture_energy: f64) -> f64 {
        if strain_rate < 1e-30 || self.density < 1e-30 {
            return f64::INFINITY;
        }
        let num = 24.0 * fracture_energy;
        let den = self.density * strain_rate * strain_rate;
        (num / den).cbrt()
    }
    /// Number of fragments per unit volume at strain rate ε̇.
    pub fn fragment_number_density(&self, strain_rate: f64) -> f64 {
        let s = self.mean_fragment_size(strain_rate);
        if s <= 0.0 || !s.is_finite() {
            return 0.0;
        }
        1.0 / (s * s * s)
    }
    /// Weibull cumulative flaw density up to strain ε.
    ///
    /// `n(ε) = k_w * ε^m_w` (flaws per unit volume).
    pub fn weibull_flaw_density(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            return 0.0;
        }
        self.weibull_k * strain.powf(self.weibull_m)
    }
    /// Activation strain for a given flaw (inverse Weibull).
    ///
    /// `ε_flaw(n) = (n / k_w)^{1/m_w}`.
    pub fn flaw_activation_strain(&self, flaw_density: f64) -> f64 {
        if flaw_density <= 0.0 || self.weibull_k <= 0.0 {
            return f64::INFINITY;
        }
        (flaw_density / self.weibull_k).powf(1.0 / self.weibull_m.max(1e-10))
    }
    /// Cumulative fragment size distribution (Mott distribution CDF).
    ///
    /// `P(s ≤ x) = 1 - exp(-x / s_mean)` (exponential distribution model).
    pub fn fragment_cdf(&self, x: f64, strain_rate: f64) -> f64 {
        let s_mean = self.mean_fragment_size(strain_rate);
        if !s_mean.is_finite() || s_mean < 1e-30 {
            return 0.0;
        }
        1.0 - (-x / s_mean).exp()
    }
    /// Kipp–Grady mean number of fragments from a sphere of radius R.
    ///
    /// `N = (4π/3 R³) * n_v(ε̇)`.
    pub fn fragment_count_sphere(&self, radius: f64, strain_rate: f64) -> f64 {
        let volume = 4.0 * PI / 3.0 * radius * radius * radius;
        volume * self.fragment_number_density(strain_rate)
    }
    /// Characteristic fragmentation time t_f = s_mean / c_l.
    pub fn fragmentation_time(&self, strain_rate: f64) -> f64 {
        let s = self.mean_fragment_size(strain_rate);
        s / self.wave_speed.max(1e-30)
    }
}
/// A single SPH particle carrying fracture-mechanics state variables.
///
/// Each particle tracks its continuum stress tensor, scalar damage variable,
/// and crack normal direction for anisotropic fracture.
#[derive(Debug, Clone)]
pub struct SphFractureParticle {
    /// Position \[x, y, z\] in metres.
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] in m/s.
    pub velocity: [f64; 3],
    /// Total Cauchy stress tensor σ_ij (3×3, symmetric).
    pub stress: Mat3,
    /// Deviatoric stress tensor S_ij = σ_ij - (1/3)tr(σ) δ_ij.
    pub deviatoric_stress: Mat3,
    /// Scalar damage variable D ∈ \[0, 1\] (0 = intact, 1 = fully damaged).
    pub damage: f64,
    /// Crack normal direction (unit vector in crack plane normal).
    pub crack_normal: [f64; 3],
    /// Mass of the particle (kg).
    pub mass: f64,
    /// Smoothing length h (m).
    pub smooth_h: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Pressure p (Pa).
    pub pressure: f64,
    /// Internal energy e (J/kg).
    pub energy: f64,
    /// Fragment label (−1 = unassigned).
    pub fragment_id: i32,
    /// Activation strain (for flaw-based damage models).
    pub activation_strain: f64,
    /// Whether this particle is on a crack surface.
    pub on_crack_surface: bool,
}
impl SphFractureParticle {
    /// Create a new intact particle at the given position.
    ///
    /// # Arguments
    /// * `position`    — initial position \[x, y, z\]
    /// * `velocity`    — initial velocity \[vx, vy, vz\]
    /// * `mass`        — particle mass (kg)
    /// * `smooth_h`    — smoothing length (m)
    /// * `density`     — initial density (kg/m³)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        smooth_h: f64,
        density: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            stress: [[0.0; 3]; 3],
            deviatoric_stress: [[0.0; 3]; 3],
            damage: 0.0,
            crack_normal: [0.0, 0.0, 1.0],
            mass,
            smooth_h,
            density,
            pressure: 0.0,
            energy: 0.0,
            fragment_id: -1,
            activation_strain: 1e-4,
            on_crack_surface: false,
        }
    }
    /// Return effective stress = (1 - D) * stress (damage-reduced).
    pub fn effective_stress(&self) -> Mat3 {
        mat3_scale(&self.stress, 1.0 - self.damage.clamp(0.0, 1.0))
    }
    /// Maximum principal stress (largest eigenvalue of stress tensor).
    pub fn max_principal_stress(&self) -> f64 {
        let eigs = symmetric_eigenvalues_3x3(&self.stress);
        eigs[2]
    }
    /// Minimum principal stress (smallest eigenvalue).
    pub fn min_principal_stress(&self) -> f64 {
        let eigs = symmetric_eigenvalues_3x3(&self.stress);
        eigs[0]
    }
    /// Von Mises equivalent stress sqrt(3/2 * S:S).
    pub fn von_mises_stress(&self) -> f64 {
        let s = &self.deviatoric_stress;
        let j2: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| s[i][j] * s[i][j]))
            .sum::<f64>()
            * 0.5;
        (3.0 * j2).sqrt()
    }
    /// Pressure from trace: p = -tr(σ)/3.
    pub fn hydrostatic_pressure(&self) -> f64 {
        -mat3_trace(&self.stress) / 3.0
    }
    /// Update deviatoric stress from total stress and pressure.
    pub fn update_deviatoric(&mut self) {
        let p = self.hydrostatic_pressure();
        for i in 0..3 {
            for j in 0..3 {
                self.deviatoric_stress[i][j] = self.stress[i][j] + if i == j { p } else { 0.0 };
            }
        }
    }
    /// Check if particle is intact (D < threshold).
    pub fn is_intact(&self, threshold: f64) -> bool {
        self.damage < threshold
    }
    /// Distance to another particle.
    pub fn distance_to(&self, other: &SphFractureParticle) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Tracks fragments formed by fracture using connected-component analysis.
///
/// Two particles are "connected" if they are within interaction range and
/// both have damage below a threshold.
#[derive(Debug, Clone)]
pub struct FragmentTracking {
    /// Maximum inter-particle distance for connectivity (m).
    pub interaction_radius: f64,
    /// Damage threshold above which particles are considered separated.
    pub damage_threshold: f64,
    /// Fragment labels for each particle (index → fragment ID).
    pub labels: Vec<i32>,
    /// Number of distinct fragments found.
    pub num_fragments: usize,
}
impl FragmentTracking {
    /// Create a new fragment tracker.
    ///
    /// # Arguments
    /// * `interaction_radius` — smoothing length × factor for connectivity
    /// * `damage_threshold`   — D above this value means the particle is "broken off"
    pub fn new(interaction_radius: f64, damage_threshold: f64) -> Self {
        Self {
            interaction_radius,
            damage_threshold,
            labels: Vec::new(),
            num_fragments: 0,
        }
    }
    /// Run connected-component labeling (union-find) on the particle array.
    ///
    /// Particles with `damage ≥ damage_threshold` are isolated nodes.
    /// Returns the number of fragments found.
    pub fn label_fragments(&mut self, particles: &[SphFractureParticle]) -> usize {
        let n = particles.len();
        self.labels = vec![-1; n];
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
        let r2 = self.interaction_radius * self.interaction_radius;
        for i in 0..n {
            if particles[i].damage >= self.damage_threshold {
                continue;
            }
            for j in (i + 1)..n {
                if particles[j].damage >= self.damage_threshold {
                    continue;
                }
                let dx = particles[i].position[0] - particles[j].position[0];
                let dy = particles[i].position[1] - particles[j].position[1];
                let dz = particles[i].position[2] - particles[j].position[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    union(&mut parent, i, j);
                }
            }
        }
        let mut label_map = std::collections::HashMap::new();
        let mut next_label = 0i32;
        for i in 0..n {
            if particles[i].damage >= self.damage_threshold {
                self.labels[i] = -1;
            } else {
                let root = find(&mut parent, i);
                let label = *label_map.entry(root).or_insert_with(|| {
                    let l = next_label;
                    next_label += 1;
                    l
                });
                self.labels[i] = label;
            }
        }
        self.num_fragments = next_label as usize;
        self.num_fragments
    }
    /// Return the number of particles in each fragment.
    pub fn fragment_sizes(&self) -> Vec<usize> {
        if self.num_fragments == 0 {
            return vec![];
        }
        let mut sizes = vec![0usize; self.num_fragments];
        for &lbl in &self.labels {
            if lbl >= 0 && (lbl as usize) < self.num_fragments {
                sizes[lbl as usize] += 1;
            }
        }
        sizes
    }
    /// Return the mass of each fragment given particle masses.
    pub fn fragment_masses(&self, particles: &[SphFractureParticle]) -> Vec<f64> {
        if self.num_fragments == 0 {
            return vec![];
        }
        let mut masses = vec![0.0f64; self.num_fragments];
        for (i, &lbl) in self.labels.iter().enumerate() {
            if lbl >= 0 && (lbl as usize) < self.num_fragments && i < particles.len() {
                masses[lbl as usize] += particles[i].mass;
            }
        }
        masses
    }
    /// Return the centroid of each fragment.
    pub fn fragment_centroids(&self, particles: &[SphFractureParticle]) -> Vec<[f64; 3]> {
        if self.num_fragments == 0 {
            return vec![];
        }
        let mut sums = vec![[0.0f64; 3]; self.num_fragments];
        let mut counts = vec![0usize; self.num_fragments];
        for (i, &lbl) in self.labels.iter().enumerate() {
            if lbl >= 0 && (lbl as usize) < self.num_fragments && i < particles.len() {
                let fid = lbl as usize;
                for k in 0..3 {
                    sums[fid][k] += particles[i].position[k];
                }
                counts[fid] += 1;
            }
        }
        sums.iter()
            .zip(counts.iter())
            .map(|(s, &c)| {
                if c == 0 {
                    [0.0; 3]
                } else {
                    [s[0] / c as f64, s[1] / c as f64, s[2] / c as f64]
                }
            })
            .collect()
    }
    /// Largest fragment index (by particle count).
    pub fn largest_fragment(&self) -> Option<usize> {
        let sizes = self.fragment_sizes();
        sizes
            .iter()
            .enumerate()
            .max_by_key(|&(_, s)| s)
            .map(|(i, _)| i)
    }
}
/// Full SPH fracture simulation loop.
///
/// Orchestrates: SPH stress update → continuum damage update →
/// crack propagation → fragment tracking.
#[derive(Debug, Clone)]
pub struct SphFractureSimulation {
    /// Particle array.
    pub particles: Vec<SphFractureParticle>,
    /// Damage model.
    pub damage_model: ContinuumDamageModel,
    /// Crack surfaces being tracked.
    pub cracks: Vec<CrackSurface>,
    /// Fracture energy model.
    pub fracture_energy: FractureEnergy,
    /// Fragment tracker.
    pub fragment_tracker: FragmentTracking,
    /// Current simulation time (s).
    pub time: f64,
    /// Timestep Δt (s).
    pub dt: f64,
    /// Total number of steps taken.
    pub step_count: u64,
    /// Smoothing length used by SPH kernels (m).
    pub smooth_h: f64,
}
impl SphFractureSimulation {
    /// Create a new SPH fracture simulation.
    ///
    /// # Arguments
    /// * `particles`       — initial particle configuration
    /// * `damage_model`    — continuum damage model
    /// * `fracture_energy` — fracture energy parameters
    /// * `dt`              — timestep (s)
    /// * `smooth_h`        — SPH smoothing length (m)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        particles: Vec<SphFractureParticle>,
        damage_model: ContinuumDamageModel,
        fracture_energy: FractureEnergy,
        dt: f64,
        smooth_h: f64,
    ) -> Self {
        let interaction_radius = 2.0 * smooth_h;
        let fragment_tracker = FragmentTracking::new(interaction_radius, 0.99);
        Self {
            particles,
            damage_model,
            cracks: Vec::new(),
            fracture_energy,
            fragment_tracker,
            time: 0.0,
            dt,
            step_count: 0,
            smooth_h,
        }
    }
    /// Step 1: update stress tensors using Hooke's law approximation.
    ///
    /// This simplified elastic constitutive update uses the current velocity
    /// gradient to compute the strain rate and update the Cauchy stress.
    pub fn update_stress(&mut self) {
        let n = self.particles.len();
        let e = self.damage_model.youngs_modulus;
        let nu = self.damage_model.poisson_ratio;
        let lame_lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)).max(1e-10);
        let lame_mu = e / (2.0 * (1.0 + nu).max(1e-10));
        let dt = self.dt;
        let h = self.smooth_h;
        for i in 0..n {
            let mut strain_rate = [[0.0f64; 3]; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r_vec = [
                    self.particles[j].position[0] - self.particles[i].position[0],
                    self.particles[j].position[1] - self.particles[i].position[1],
                    self.particles[j].position[2] - self.particles[i].position[2],
                ];
                let r = norm3(r_vec);
                if r < 1e-30 || r > 2.0 * h {
                    continue;
                }
                let dw_dr = wendland_c2_3d_grad(r, h);
                let vol_j = self.particles[j].mass / self.particles[j].density.max(1e-30);
                for alpha in 0..3 {
                    let dv_alpha =
                        self.particles[j].velocity[alpha] - self.particles[i].velocity[alpha];
                    for beta in 0..3 {
                        let grad_w = -dw_dr * r_vec[beta] / r;
                        strain_rate[alpha][beta] += vol_j * dv_alpha * grad_w;
                    }
                }
            }
            let mut eps_dot = [[0.0f64; 3]; 3];
            for alpha in 0..3 {
                for beta in 0..3 {
                    eps_dot[alpha][beta] =
                        0.5 * (strain_rate[alpha][beta] + strain_rate[beta][alpha]);
                }
            }
            let tr_eps = eps_dot[0][0] + eps_dot[1][1] + eps_dot[2][2];
            let d = self.particles[i].damage;
            let eff = 1.0 - d.clamp(0.0, 1.0);
            for alpha in 0..3 {
                for beta in 0..3 {
                    let delta = if alpha == beta { 1.0 } else { 0.0 };
                    self.particles[i].stress[alpha][beta] += eff
                        * (lame_lambda * tr_eps * delta + 2.0 * lame_mu * eps_dot[alpha][beta])
                        * dt;
                }
            }
            self.particles[i].update_deviatoric();
        }
    }
    /// Step 2: update damage for all particles.
    pub fn update_damage(&mut self) {
        let dt = self.dt;
        let n = self.particles.len();
        let mut new_damages = vec![0.0f64; n];
        for i in 0..n {
            new_damages[i] = self.damage_model.update_damage(&self.particles[i], dt);
        }
        for i in 0..n {
            self.particles[i].damage = new_damages[i];
            if self.particles[i].damage > 0.01 {
                let n_dir = max_eigenvector_3x3(&self.particles[i].stress);
                self.particles[i].crack_normal = n_dir;
                self.particles[i].on_crack_surface = self.particles[i].damage > 0.5;
            }
        }
    }
    /// Step 3: propagate all tracked crack surfaces.
    pub fn propagate_cracks(&mut self) {
        let e = self.damage_model.youngs_modulus;
        let nu_val = self.damage_model.poisson_ratio;
        let k_ic = self.fracture_energy.k_ic(nu_val);
        for crack in self.cracks.iter_mut() {
            let tip_idx = crack.nearest_particle_index(&self.particles);
            if tip_idx < self.particles.len() {
                let p = &self.particles[tip_idx];
                let sigma_yy = p.stress[1][1].max(0.0);
                let sigma_xy = p.stress[0][1];
                let r_tip = crack.smooth_h_or_default(self.smooth_h);
                crack.compute_k_factors(sigma_yy, sigma_xy, r_tip);
                crack.update_cod(e, r_tip);
                if crack.should_propagate() {
                    crack.update_propagation_dir();
                    crack.advance(self.smooth_h * 0.5);
                } else if (crack.k_i * crack.k_i + crack.k_ii * crack.k_ii).sqrt() < k_ic * 0.1 {
                    crack.arrest();
                }
            }
        }
    }
    /// Full timestep: stress → damage → cracks → time advance.
    pub fn step(&mut self) {
        self.update_stress();
        self.update_damage();
        self.propagate_cracks();
        let dt = self.dt;
        for p in self.particles.iter_mut() {
            for k in 0..3 {
                p.position[k] += p.velocity[k] * dt;
            }
        }
        self.time += dt;
        self.step_count += 1;
    }
    /// Run `n_steps` timesteps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }
    /// Identify fragments after simulation.
    pub fn identify_fragments(&mut self) -> usize {
        self.fragment_tracker.label_fragments(&self.particles)
    }
    /// Add a crack surface to track.
    pub fn add_crack(&mut self, crack: CrackSurface) {
        self.cracks.push(crack);
    }
    /// Maximum damage across all particles.
    pub fn max_damage(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.damage)
            .fold(0.0_f64, f64::max)
    }
    /// Mean damage across all particles.
    pub fn mean_damage(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let s: f64 = self.particles.iter().map(|p| p.damage).sum();
        s / self.particles.len() as f64
    }
    /// Number of fully damaged particles (D ≥ 0.99).
    pub fn fully_damaged_count(&self) -> usize {
        self.particles.iter().filter(|p| p.damage >= 0.99).count()
    }
    /// Total kinetic energy of all particles.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| 0.5 * p.mass * dot3(p.velocity, p.velocity))
            .sum()
    }
}
/// Fracture energy computations: Griffith criterion and cohesive zone.
///
/// Implements the energy balance for crack propagation and the
/// cohesive traction–separation law for the process zone.
#[derive(Debug, Clone)]
pub struct FractureEnergy {
    /// Critical energy release rate G_c (J/m²).
    pub g_c: f64,
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Cohesive zone peak traction t_0 (Pa).
    pub cohesive_traction: f64,
    /// Cohesive zone critical separation δ_c (m).
    pub critical_separation: f64,
    /// Plane stress (true) or plane strain (false).
    pub plane_stress: bool,
}
impl FractureEnergy {
    /// Create a new FractureEnergy from Griffith parameters.
    ///
    /// # Arguments
    /// * `g_c`                — critical energy release rate (J/m²)
    /// * `youngs_modulus`     — E (Pa)
    /// * `cohesive_traction`  — peak cohesive traction t_0 (Pa)
    /// * `critical_separation`— δ_c at which traction vanishes (m)
    /// * `plane_stress`       — true for plane stress, false for plane strain
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        g_c: f64,
        youngs_modulus: f64,
        cohesive_traction: f64,
        critical_separation: f64,
        plane_stress: bool,
    ) -> Self {
        Self {
            g_c,
            youngs_modulus,
            cohesive_traction,
            critical_separation,
            plane_stress,
        }
    }
    /// Critical stress intensity factor K_Ic from G_c and E.
    ///
    /// Plane stress: `K_Ic = sqrt(G_c * E)`
    /// Plane strain: `K_Ic = sqrt(G_c * E / (1 - ν²))`
    pub fn k_ic(&self, poisson_ratio: f64) -> f64 {
        if self.plane_stress {
            (self.g_c * self.youngs_modulus).sqrt()
        } else {
            let factor = 1.0 - poisson_ratio * poisson_ratio;
            (self.g_c * self.youngs_modulus / factor.max(1e-10)).sqrt()
        }
    }
    /// Griffith crack length: half-length `a` for a given applied stress σ.
    ///
    /// `a = G_c E / (π σ²)` (plane stress approximation).
    pub fn griffith_crack_length(&self, applied_stress: f64) -> f64 {
        if applied_stress.abs() < 1e-20 {
            return f64::INFINITY;
        }
        self.g_c * self.youngs_modulus / (PI * applied_stress * applied_stress)
    }
    /// Linear cohesive traction–separation law.
    ///
    /// `t(δ) = t_0 * (1 - δ/δ_c)` for `δ ≤ δ_c`, else 0.
    pub fn cohesive_traction_linear(&self, separation: f64) -> f64 {
        if separation >= self.critical_separation || separation < 0.0 {
            return 0.0;
        }
        self.cohesive_traction * (1.0 - separation / self.critical_separation)
    }
    /// Exponential cohesive traction–separation law.
    ///
    /// `t(δ) = t_0 * (δ/δ_c) * exp(1 - δ/δ_c)`.
    pub fn cohesive_traction_exponential(&self, separation: f64) -> f64 {
        if separation < 0.0 {
            return 0.0;
        }
        let xi = separation / self.critical_separation.max(1e-30);
        self.cohesive_traction * xi * (1.0 - xi).exp()
    }
    /// Energy dissipated in cohesive zone (area under T–δ curve, linear law).
    ///
    /// `W = 0.5 * t_0 * δ_c`
    pub fn cohesive_energy(&self) -> f64 {
        0.5 * self.cohesive_traction * self.critical_separation
    }
    /// Check Griffith criterion: G ≥ G_c.
    pub fn griffith_criterion_met(&self, g: f64) -> bool {
        g >= self.g_c
    }
    /// Cohesive zone length (Dugdale model).
    ///
    /// `r_p = (π/8) * (K_Ic / t_0)²`
    pub fn dugdale_process_zone(&self, poisson_ratio: f64) -> f64 {
        let kic = self.k_ic(poisson_ratio);
        PI / 8.0 * (kic / self.cohesive_traction.max(1e-30)).powi(2)
    }
    /// SPH cohesive force between two particles bridging a crack.
    ///
    /// Returns cohesive force magnitude based on opening displacement.
    pub fn sph_cohesive_force(&self, separation: f64, area: f64, use_exponential: bool) -> f64 {
        let traction = if use_exponential {
            self.cohesive_traction_exponential(separation)
        } else {
            self.cohesive_traction_linear(separation)
        };
        traction * area
    }
}
/// Tracks the crack tip position and propagation direction.
///
/// Uses the maximum-tensile-stress (MTS) criterion to determine the
/// crack propagation angle.
#[derive(Debug, Clone)]
pub struct CrackSurface {
    /// Crack tip position.
    pub tip_position: [f64; 3],
    /// Crack propagation direction (unit vector).
    pub propagation_dir: [f64; 3],
    /// Crack opening displacement (mode I) in metres.
    pub cod: f64,
    /// Crack length (accumulated arc length in metres).
    pub crack_length: f64,
    /// Stress intensity factor K_I (Pa·√m).
    pub k_i: f64,
    /// Stress intensity factor K_II (Pa·√m).
    pub k_ii: f64,
    /// Critical stress intensity factor K_Ic (Pa·√m).
    pub k_ic: f64,
    /// List of crack path points.
    pub crack_path: Vec<[f64; 3]>,
    /// Whether the crack has arrested.
    pub arrested: bool,
}
impl CrackSurface {
    /// Create a new crack surface at a given tip position.
    ///
    /// # Arguments
    /// * `tip_position`  — initial crack tip coordinates
    /// * `k_ic`          — fracture toughness K_Ic (Pa·√m)
    pub fn new(tip_position: [f64; 3], k_ic: f64) -> Self {
        Self {
            tip_position,
            propagation_dir: [1.0, 0.0, 0.0],
            cod: 0.0,
            crack_length: 0.0,
            k_i: 0.0,
            k_ii: 0.0,
            k_ic,
            crack_path: vec![tip_position],
            arrested: false,
        }
    }
    /// Compute stress intensity factors from the near-tip stress field.
    ///
    /// Uses the Irwin K-factor formulas for a mode-I / mode-II crack.
    ///
    /// # Arguments
    /// * `sigma_yy` — tensile stress perpendicular to crack (Pa)
    /// * `sigma_xy` — shear stress along crack plane (Pa)
    /// * `r_tip`    — distance from tip to evaluation point (m)
    pub fn compute_k_factors(&mut self, sigma_yy: f64, sigma_xy: f64, r_tip: f64) {
        if r_tip < 1e-15 {
            return;
        }
        let factor = (2.0 * PI * r_tip).sqrt();
        self.k_i = sigma_yy * factor;
        self.k_ii = sigma_xy * factor;
    }
    /// Maximum-tensile-stress criterion: compute crack propagation angle.
    ///
    /// Returns angle θ in radians relative to the current crack plane.
    /// Formula: θ = 2 * atan((K_I - sqrt(K_I² + 8 K_II²)) / (4 K_II)) for K_II ≠ 0.
    pub fn mts_propagation_angle(&self) -> f64 {
        let ki = self.k_i;
        let kii = self.k_ii;
        if kii.abs() < 1e-15 {
            return 0.0;
        }
        let num = ki - (ki * ki + 8.0 * kii * kii).sqrt();
        let den = 4.0 * kii;
        2.0 * (num / den).atan()
    }
    /// Check if crack should propagate (K_eff ≥ K_Ic).
    ///
    /// K_eff = sqrt(K_I² + K_II²) (simplified mixed-mode criterion).
    pub fn should_propagate(&self) -> bool {
        if self.arrested {
            return false;
        }
        let k_eff = (self.k_i * self.k_i + self.k_ii * self.k_ii).sqrt();
        k_eff >= self.k_ic
    }
    /// Advance crack tip by `delta_a` along the propagation direction.
    ///
    /// Updates tip position and crack path.
    pub fn advance(&mut self, delta_a: f64) {
        if self.arrested || delta_a <= 0.0 {
            return;
        }
        let dir = normalize3(self.propagation_dir);
        self.tip_position[0] += dir[0] * delta_a;
        self.tip_position[1] += dir[1] * delta_a;
        self.tip_position[2] += dir[2] * delta_a;
        self.crack_length += delta_a;
        self.crack_path.push(self.tip_position);
    }
    /// Update propagation direction from the stress field at the tip.
    ///
    /// Rotates the current direction by the MTS angle.
    pub fn update_propagation_dir(&mut self) {
        let theta = self.mts_propagation_angle();
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        let dx = self.propagation_dir[0];
        let dy = self.propagation_dir[1];
        self.propagation_dir[0] = cos_t * dx - sin_t * dy;
        self.propagation_dir[1] = sin_t * dx + cos_t * dy;
        self.propagation_dir = normalize3(self.propagation_dir);
    }
    /// Arrest the crack (sets `arrested = true`).
    pub fn arrest(&mut self) {
        self.arrested = true;
    }
    /// Crack opening displacement from K_I and material constants.
    ///
    /// COD = 4 K_I / (E * sqrt(2π/r)).
    pub fn update_cod(&mut self, youngs_modulus: f64, r_tip: f64) {
        if r_tip < 1e-15 || youngs_modulus < 1e-10 {
            return;
        }
        self.cod = 4.0 * self.k_i * (r_tip / (2.0 * PI)).sqrt() / youngs_modulus;
    }
    /// Nearest particle index to the crack tip.
    pub fn nearest_particle_index(&self, particles: &[SphFractureParticle]) -> usize {
        let mut best = 0usize;
        let mut best_dist = f64::INFINITY;
        for (k, p) in particles.iter().enumerate() {
            let dx = p.position[0] - self.tip_position[0];
            let dy = p.position[1] - self.tip_position[1];
            let dz = p.position[2] - self.tip_position[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_dist {
                best_dist = d2;
                best = k;
            }
        }
        best
    }
}
