#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::*;

/// Polytropic equation of state: P = K * ρ^Γ.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PolytropicEos {
    /// Polytropic constant K.
    pub k: f64,
    /// Polytropic exponent Γ = 1 + 1/n where n is the polytropic index.
    pub gamma_poly: f64,
}
#[allow(dead_code)]
impl PolytropicEos {
    /// Create a new polytropic EOS.
    pub fn new(k: f64, gamma_poly: f64) -> Self {
        Self { k, gamma_poly }
    }
    /// Create from polytropic index n: Γ = 1 + 1/n.
    pub fn from_index(k: f64, n: f64) -> Self {
        Self {
            k,
            gamma_poly: 1.0 + 1.0 / n,
        }
    }
    /// Isothermal case (Γ = 1, P = K * ρ).
    pub fn isothermal(k: f64) -> Self {
        Self { k, gamma_poly: 1.0 }
    }
    /// Pressure: P = K * ρ^Γ.
    pub fn pressure(&self, density: f64) -> f64 {
        self.k * density.powf(self.gamma_poly)
    }
    /// Sound speed: c_s = sqrt(Γ * P / ρ) = sqrt(Γ * K * ρ^{Γ-1}).
    pub fn sound_speed(&self, density: f64) -> f64 {
        (self.gamma_poly * self.k * density.powf(self.gamma_poly - 1.0)).sqrt()
    }
    /// Internal energy per unit mass: u = K * ρ^{Γ-1} / (Γ - 1).
    /// Only valid for Γ ≠ 1 (not isothermal).
    pub fn internal_energy(&self, density: f64) -> f64 {
        if (self.gamma_poly - 1.0).abs() < 1e-15 {
            self.k * density.ln().max(0.0)
        } else {
            self.k * density.powf(self.gamma_poly - 1.0) / (self.gamma_poly - 1.0)
        }
    }
    /// Density from pressure: ρ = (P / K)^{1/Γ}.
    pub fn density_from_pressure(&self, pressure: f64) -> f64 {
        (pressure / self.k).powf(1.0 / self.gamma_poly)
    }
    /// Lane-Emden characteristic radius scale: r_n = sqrt((n+1) K / (4π G)) * ρ_c^{(1-n)/(2n)}.
    pub fn lane_emden_scale(&self, central_density: f64, g: f64) -> f64 {
        let n = 1.0 / (self.gamma_poly - 1.0);
        let factor = (n + 1.0) * self.k / (4.0 * PI * g);
        factor.sqrt() * central_density.powf((1.0 - n) / (2.0 * n))
    }
}
/// Self-gravitating SPH with Plummer softening.
///
/// Computes gravitational acceleration between all particle pairs
/// using direct O(N^2) summation with a softened potential.
#[derive(Debug, Clone)]
pub struct GravitationalSph {
    /// Gravitational constant (in code units).
    pub g_const: f64,
    /// Softening length to prevent singularities.
    pub softening: f64,
    /// Opening angle for optional tree code (not used in direct sum).
    pub _theta: f64,
}
impl GravitationalSph {
    /// Create a new gravitational SPH solver.
    pub fn new(g_const: f64, softening: f64) -> Self {
        Self {
            g_const,
            softening,
            _theta: 0.5,
        }
    }
    /// Softened gravitational potential between two particles.
    ///
    /// phi = -G * m / sqrt(r^2 + eps^2)
    pub fn potential_pair(&self, mass: f64, r: f64) -> f64 {
        let eps2 = self.softening * self.softening;
        -self.g_const * mass / (r * r + eps2).sqrt()
    }
    /// Softened gravitational acceleration on particle i due to particle j.
    ///
    /// a_i = -G * m_j * (r_i - r_j) / (|r_ij|^2 + eps^2)^{3/2}
    pub fn acceleration_pair(&self, ri: [f64; 3], rj: [f64; 3], mj: f64) -> [f64; 3] {
        let dx = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
        let r2 = dot3(dx, dx);
        let eps2 = self.softening * self.softening;
        let denom = (r2 + eps2).powf(1.5);
        if denom < 1e-30 {
            return [0.0; 3];
        }
        let factor = -self.g_const * mj / denom;
        [factor * dx[0], factor * dx[1], factor * dx[2]]
    }
    /// Compute gravitational acceleration for all particles (direct O(N^2)).
    pub fn compute_gravity(&self, particles: &mut [AstroParticle]) {
        let n = particles.len();
        let positions: Vec<[f64; 3]> = particles.iter().map(|p| p.pos).collect();
        let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
        for i in 0..n {
            let mut acc = [0.0; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let a = self.acceleration_pair(positions[i], positions[j], masses[j]);
                acc[0] += a[0];
                acc[1] += a[1];
                acc[2] += a[2];
            }
            particles[i].acc[0] += acc[0];
            particles[i].acc[1] += acc[1];
            particles[i].acc[2] += acc[2];
        }
    }
    /// Total gravitational potential energy of the system.
    pub fn total_potential_energy(&self, particles: &[AstroParticle]) -> f64 {
        let n = particles.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(particles[i].pos, particles[j].pos);
                energy += self
                    .potential_pair(particles[i].mass * particles[j].mass / particles[j].mass, r)
                    * particles[j].mass;
            }
        }
        energy
    }
    /// Compute gravitational potential energy properly.
    pub fn gravitational_potential_energy(&self, particles: &[AstroParticle]) -> f64 {
        let n = particles.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(particles[i].pos, particles[j].pos);
                let eps2 = self.softening * self.softening;
                energy -=
                    self.g_const * particles[i].mass * particles[j].mass / (r * r + eps2).sqrt();
            }
        }
        energy
    }
    /// Tidal tensor (second derivative of potential) at a point.
    ///
    /// Returns the 3x3 tidal tensor T_ij = d^2 phi / dx_i dx_j.
    pub fn tidal_tensor(&self, pos: [f64; 3], particles: &[AstroParticle]) -> [[f64; 3]; 3] {
        let mut t = [[0.0; 3]; 3];
        let eps2 = self.softening * self.softening;
        for p in particles {
            let dx = sub3(pos, p.pos);
            let r2 = dot3(dx, dx) + eps2;
            let r5 = r2.powf(2.5);
            let r3 = r2.powf(1.5);
            if r3 < 1e-30 {
                continue;
            }
            for a in 0..3 {
                for b in 0..3 {
                    let delta_ab = if a == b { 1.0 } else { 0.0 };
                    t[a][b] += self.g_const * p.mass * (3.0 * dx[a] * dx[b] / r5 - delta_ab / r3);
                }
            }
        }
        t
    }
}
/// Jeans instability analysis for gravitational collapse.
///
/// Provides calculations for Jeans mass, Jeans length, and free-fall time.
#[derive(Debug, Clone, Copy)]
pub struct JeansInstability {
    /// Sound speed c_s.
    pub sound_speed: f64,
    /// Background density rho_0.
    pub density: f64,
    /// Gravitational constant.
    pub g_const: f64,
}
impl JeansInstability {
    /// Create a new Jeans instability analyzer.
    pub fn new(sound_speed: f64, density: f64, g_const: f64) -> Self {
        Self {
            sound_speed,
            density,
            g_const,
        }
    }
    /// Jeans length: lambda_J = c_s * sqrt(pi / (G * rho)).
    pub fn jeans_length(&self) -> f64 {
        self.sound_speed * (PI / (self.g_const * self.density)).sqrt()
    }
    /// Jeans mass: M_J = (pi/6) * rho * lambda_J^3.
    pub fn jeans_mass(&self) -> f64 {
        let lj = self.jeans_length();
        (PI / 6.0) * self.density * lj * lj * lj
    }
    /// Jeans wavenumber: k_J = sqrt(4 * pi * G * rho) / c_s.
    pub fn jeans_wavenumber(&self) -> f64 {
        (4.0 * PI * self.g_const * self.density).sqrt() / self.sound_speed
    }
    /// Free-fall time: t_ff = sqrt(3 * pi / (32 * G * rho)).
    pub fn free_fall_time(&self) -> f64 {
        (3.0 * PI / (32.0 * self.g_const * self.density)).sqrt()
    }
    /// Growth rate of the Jeans instability for a given wavenumber k.
    ///
    /// omega^2 = c_s^2 * k^2 - 4*pi*G*rho
    /// Returns real part of omega (imaginary if unstable).
    pub fn growth_rate(&self, k: f64) -> f64 {
        let omega2 =
            self.sound_speed * self.sound_speed * k * k - 4.0 * PI * self.g_const * self.density;
        if omega2 >= 0.0 {
            omega2.sqrt()
        } else {
            -(-omega2).sqrt()
        }
    }
    /// Check if a perturbation with wavenumber k is Jeans-unstable.
    pub fn is_unstable(&self, k: f64) -> bool {
        k < self.jeans_wavenumber()
    }
    /// Bonnor-Ebert mass (critical mass for isothermal sphere collapse).
    ///
    /// M_BE ~ 1.18 * c_s^4 / (G^{3/2} * rho^{1/2} * P_ext^{1/2})
    /// Simplified version using external pressure P_ext.
    pub fn bonnor_ebert_mass(&self, p_ext: f64) -> f64 {
        if p_ext <= 0.0 || self.g_const <= 0.0 {
            return 0.0;
        }
        1.18 * self.sound_speed.powi(4) / (self.g_const.powf(1.5) * p_ext.sqrt())
    }
}
/// Star formation model using sink particles.
///
/// Creates sink particles when gas exceeds density threshold
/// and satisfies additional collapse criteria.
#[derive(Debug, Clone)]
pub struct StarFormation {
    /// Density threshold for sink creation.
    pub density_threshold: f64,
    /// Accretion radius for sink particles.
    pub accretion_radius: f64,
    /// Minimum Jeans number (resolution criterion).
    pub min_jeans_number: f64,
    /// Star formation efficiency (0 to 1).
    pub efficiency: f64,
}
impl StarFormation {
    /// Create a new star formation model.
    pub fn new(
        density_threshold: f64,
        accretion_radius: f64,
        min_jeans_number: f64,
        efficiency: f64,
    ) -> Self {
        Self {
            density_threshold,
            accretion_radius,
            min_jeans_number,
            efficiency,
        }
    }
    /// Check if a particle should form a sink.
    ///
    /// Criteria: density > threshold, converging flow, gravitationally bound.
    pub fn should_form_sink(
        &self,
        particle: &AstroParticle,
        div_v: f64,
        neighbors: &[AstroParticle],
        g_const: f64,
    ) -> bool {
        if particle.density < self.density_threshold {
            return false;
        }
        if div_v >= 0.0 {
            return false;
        }
        let e_kin = particle.kinetic_energy();
        let e_therm = particle.thermal_energy();
        let mut e_grav = 0.0;
        for nb in neighbors {
            let r = dist3(particle.pos, nb.pos);
            if r > 1e-30 {
                e_grav -= g_const * particle.mass * nb.mass / r;
            }
        }
        e_grav.abs() > e_kin + e_therm
    }
    /// Create a sink particle from a gas particle.
    pub fn create_sink(&self, particle: &AstroParticle) -> AstroParticle {
        AstroParticle::new_sink(
            particle.mass * self.efficiency,
            particle.pos,
            particle.vel,
            self.accretion_radius,
        )
    }
    /// Accrete gas particles onto a sink particle.
    ///
    /// Returns the updated sink mass and momentum, and indices of accreted particles.
    pub fn accrete(&self, sink: &mut AstroParticle, particles: &[AstroParticle]) -> Vec<usize> {
        let mut accreted = Vec::new();
        for (i, p) in particles.iter().enumerate() {
            if p.is_sink {
                continue;
            }
            let r = dist3(sink.pos, p.pos);
            if r < sink.accretion_radius {
                let v_rel = sub3(p.vel, sink.vel);
                let e_kin = 0.5 * p.mass * dot3(v_rel, v_rel);
                let e_grav = if r > 1e-30 {
                    sink.mass * p.mass / r
                } else {
                    f64::MAX
                };
                if e_grav > e_kin {
                    let total_mass = sink.mass + p.mass;
                    for d in 0..3 {
                        sink.vel[d] = (sink.mass * sink.vel[d] + p.mass * p.vel[d]) / total_mass;
                    }
                    sink.mass = total_mass;
                    accreted.push(i);
                }
            }
        }
        accreted
    }
    /// Star formation rate (Schmidt law): dM*/dt = efficiency * rho / t_ff.
    pub fn star_formation_rate(&self, density: f64, g_const: f64) -> f64 {
        let t_ff = (3.0 * PI / (32.0 * g_const * density)).sqrt();
        if t_ff > 1e-30 {
            self.efficiency * density / t_ff
        } else {
            0.0
        }
    }
}
/// Adiabatic (ideal gas) equation of state: P = (gamma - 1) * rho * u.
///
/// For monatomic gas, gamma = 5/3; for diatomic, gamma = 7/5.
#[derive(Debug, Clone, Copy)]
pub struct AdiabaticEos {
    /// Adiabatic index (ratio of specific heats).
    pub gamma: f64,
}
impl AdiabaticEos {
    /// Create a new adiabatic EOS with given gamma.
    pub fn new(gamma: f64) -> Self {
        assert!(gamma > 1.0, "gamma must be > 1");
        Self { gamma }
    }
    /// Monatomic ideal gas (gamma = 5/3).
    pub fn monatomic() -> Self {
        Self::new(5.0 / 3.0)
    }
    /// Diatomic ideal gas (gamma = 7/5).
    pub fn diatomic() -> Self {
        Self::new(7.0 / 5.0)
    }
    /// Compute pressure from density and specific internal energy.
    pub fn pressure(&self, density: f64, internal_energy: f64) -> f64 {
        (self.gamma - 1.0) * density * internal_energy
    }
    /// Compute sound speed: c_s = sqrt(gamma * (gamma - 1) * u).
    pub fn sound_speed(&self, internal_energy: f64) -> f64 {
        (self.gamma * (self.gamma - 1.0) * internal_energy).sqrt()
    }
    /// Compute sound speed from density and pressure.
    pub fn sound_speed_from_pressure(&self, density: f64, pressure: f64) -> f64 {
        if density > 1e-30 {
            (self.gamma * pressure / density).sqrt()
        } else {
            0.0
        }
    }
    /// Compute temperature from internal energy (ideal gas: u = k_B T / ((gamma-1) mu m_p)).
    pub fn temperature(&self, internal_energy: f64, mean_molecular_weight: f64) -> f64 {
        (self.gamma - 1.0) * mean_molecular_weight * M_PROTON_CGS * internal_energy / K_BOLTZ_CGS
    }
    /// Compute internal energy from temperature.
    pub fn internal_energy_from_temperature(
        &self,
        temperature: f64,
        mean_molecular_weight: f64,
    ) -> f64 {
        K_BOLTZ_CGS * temperature / ((self.gamma - 1.0) * mean_molecular_weight * M_PROTON_CGS)
    }
    /// Adiabatic density-pressure relation: P = K * rho^gamma.
    pub fn polytropic_pressure(&self, density: f64, entropy_constant: f64) -> f64 {
        entropy_constant * density.powf(self.gamma)
    }
    /// Apply EOS to all particles, setting pressure.
    pub fn apply(&self, particles: &mut [AstroParticle]) {
        for p in particles.iter_mut() {
            p.pressure = self.pressure(p.density, p.internal_energy);
        }
    }
}
/// Supernova feedback model parameters.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SupernovaFeedback {
    /// Energy per supernova event \[erg\] (canonical: 1e51).
    pub energy_per_sn: f64,
    /// Mass fraction returned to the ISM per SN.
    pub mass_return_fraction: f64,
    /// Supernova rate per solar mass of stars formed per year.
    pub sn_rate_per_msun: f64,
    /// Coupling efficiency (fraction of energy actually deposited).
    pub coupling_efficiency: f64,
    /// Minimum number of neighbours to distribute energy.
    pub min_neighbours: usize,
}
#[allow(dead_code)]
impl SupernovaFeedback {
    /// Create a standard supernova feedback model.
    pub fn standard() -> Self {
        Self {
            energy_per_sn: 1e51,
            mass_return_fraction: 0.1,
            sn_rate_per_msun: 0.01,
            coupling_efficiency: 0.1,
            min_neighbours: 32,
        }
    }
    /// Compute the total energy injected for a given stellar mass formed.
    pub fn energy_injected(&self, stellar_mass_formed: f64) -> f64 {
        let n_sn = self.sn_rate_per_msun * stellar_mass_formed / M_SUN_CGS;
        n_sn * self.energy_per_sn * self.coupling_efficiency
    }
    /// Compute the mass returned to the ISM.
    pub fn mass_returned(&self, stellar_mass_formed: f64) -> f64 {
        self.mass_return_fraction * stellar_mass_formed
    }
    /// Distribute supernova energy thermally among neighbour particles.
    ///
    /// Returns the energy increment per neighbour particle.
    pub fn thermal_energy_per_neighbour(
        &self,
        stellar_mass_formed: f64,
        n_neighbours: usize,
    ) -> f64 {
        if n_neighbours == 0 {
            return 0.0;
        }
        self.energy_injected(stellar_mass_formed) / n_neighbours as f64
    }
    /// Compute the Sedov-Taylor blast wave radius for the given
    /// energy and ambient density at time t.
    pub fn sedov_taylor_radius(&self, energy: f64, ambient_density: f64, time: f64) -> f64 {
        let xi = 1.15;
        xi * (energy * time * time / ambient_density).powf(0.2)
    }
    /// Compute the momentum injection from a supernova in the
    /// snowplough (momentum-conserving) phase.
    pub fn momentum_injection(&self, energy: f64, ambient_density: f64) -> f64 {
        let e_ratio = energy / 1e51;
        let n_h = ambient_density / M_PROTON_CGS;
        3e5 * M_SUN_CGS * 1e5 * e_ratio.powf(0.8) * n_h.powf(-0.2)
    }
}
/// Cosmological parameters for comoving SPH.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Cosmology {
    /// Current scale factor a(t).
    pub scale_factor: f64,
    /// Hubble constant H_0 \[km/s/Mpc\] (in code units: 1/time).
    pub h0: f64,
    /// Matter density parameter Ω_m.
    pub omega_m: f64,
    /// Dark energy density parameter Ω_Λ.
    pub omega_lambda: f64,
    /// Radiation density parameter Ω_r (usually small).
    pub omega_r: f64,
}
#[allow(dead_code)]
impl Cosmology {
    /// Create a standard ΛCDM cosmology (Planck 2018 values).
    pub fn planck2018() -> Self {
        Self {
            scale_factor: 1.0,
            h0: 67.66,
            omega_m: 0.3111,
            omega_lambda: 0.6889,
            omega_r: 9.1e-5,
        }
    }
    /// Create a cosmology with given parameters.
    pub fn new(h0: f64, omega_m: f64, omega_lambda: f64) -> Self {
        Self {
            scale_factor: 1.0,
            h0,
            omega_m,
            omega_lambda,
            omega_r: 0.0,
        }
    }
    /// Hubble parameter H(a) = H_0 * E(a).
    pub fn hubble(&self, a: f64) -> f64 {
        self.h0 * self.expansion_function(a)
    }
    /// Dimensionless expansion function E(a) = sqrt(Ω_r/a^4 + Ω_m/a^3 + Ω_Λ).
    pub fn expansion_function(&self, a: f64) -> f64 {
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a2 * a2;
        (self.omega_r / a4 + self.omega_m / a3 + self.omega_lambda).sqrt()
    }
    /// Redshift from scale factor: z = 1/a - 1.
    pub fn redshift(&self, a: f64) -> f64 {
        1.0 / a - 1.0
    }
    /// Scale factor from redshift: a = 1/(1+z).
    pub fn scale_factor_from_z(&self, z: f64) -> f64 {
        1.0 / (1.0 + z)
    }
    /// Convert physical coordinates to comoving: r_comov = r_phys / a.
    pub fn to_comoving(&self, physical: [f64; 3], a: f64) -> [f64; 3] {
        [physical[0] / a, physical[1] / a, physical[2] / a]
    }
    /// Convert comoving coordinates to physical: r_phys = a * r_comov.
    pub fn to_physical(&self, comoving: [f64; 3], a: f64) -> [f64; 3] {
        [comoving[0] * a, comoving[1] * a, comoving[2] * a]
    }
    /// Comoving velocity: v_pec = v_phys - H*a * r_comov (peculiar velocity).
    pub fn peculiar_velocity(
        &self,
        physical_velocity: [f64; 3],
        comoving_position: [f64; 3],
        a: f64,
    ) -> [f64; 3] {
        let ha = self.hubble(a) * a;
        [
            physical_velocity[0] - ha * comoving_position[0],
            physical_velocity[1] - ha * comoving_position[1],
            physical_velocity[2] - ha * comoving_position[2],
        ]
    }
    /// Critical density ρ_crit = 3 H^2 / (8π G).
    pub fn critical_density(&self, a: f64) -> f64 {
        let h = self.hubble(a);
        let h_cgs = h * 1e5 / (3.086e24);
        3.0 * h_cgs * h_cgs / (8.0 * PI * G_CGS)
    }
    /// Comoving distance for a flat universe (simplified trapezoidal integration).
    pub fn comoving_distance(&self, z: f64, n_steps: usize) -> f64 {
        let dz = z / n_steps as f64;
        let c_over_h0 = 2.998e5 / self.h0;
        let mut integral = 0.0;
        for i in 0..n_steps {
            let z1 = i as f64 * dz;
            let z2 = (i + 1) as f64 * dz;
            let a1 = 1.0 / (1.0 + z1);
            let a2 = 1.0 / (1.0 + z2);
            integral +=
                0.5 * dz * (1.0 / self.expansion_function(a1) + 1.0 / self.expansion_function(a2));
        }
        c_over_h0 * integral
    }
    /// Advance the scale factor by dt using da/dt = a * H(a).
    pub fn advance_scale_factor(&mut self, dt: f64) {
        let h = self.hubble(self.scale_factor);
        let h_per_s = h * 1e5 / 3.086e24;
        self.scale_factor += self.scale_factor * h_per_s * dt;
    }
}
/// Radiative cooling function type.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub enum CoolingFunction {
    /// No cooling.
    None,
    /// Optically thin cooling with a piecewise power-law Λ(T).
    OpticallyThin,
    /// Bremsstrahlung (free-free) cooling: Λ ∝ T^{1/2}.
    Bremsstrahlung,
    /// Atomic line cooling with a peak around 1e5 K.
    AtomicLine,
    /// Tabulated cooling (simplified: Sutherland & Dopita style).
    Tabulated,
}
/// A node in the Barnes-Hut octree.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OctreeNode {
    /// Bounding box of this node.
    pub bbox: BBox3,
    /// Total mass in this node.
    pub total_mass: f64,
    /// Centre of mass.
    pub com: [f64; 3],
    /// Children (None = leaf or empty).
    pub children: [Option<Box<OctreeNode>>; 8],
    /// Particle index (only set for leaf nodes with exactly one particle).
    pub particle_idx: Option<usize>,
}
#[allow(dead_code)]
impl OctreeNode {
    /// Create an empty octree node.
    pub fn new(bbox: BBox3) -> Self {
        Self {
            bbox,
            total_mass: 0.0,
            com: [0.0; 3],
            children: Default::default(),
            particle_idx: None,
        }
    }
    /// Insert a particle into the tree.
    pub fn insert(&mut self, idx: usize, pos: [f64; 3], mass: f64, max_depth: usize) {
        if max_depth == 0 {
            let tm = self.total_mass + mass;
            if tm > 0.0 {
                for d in 0..3 {
                    self.com[d] = (self.com[d] * self.total_mass + pos[d] * mass) / tm;
                }
            }
            self.total_mass = tm;
            return;
        }
        if self.total_mass == 0.0 && self.particle_idx.is_none() {
            self.total_mass = mass;
            self.com = pos;
            self.particle_idx = Some(idx);
            return;
        }
        if let Some(old_idx) = self.particle_idx.take() {
            let old_pos = self.com;
            let old_mass = self.total_mass;
            self.total_mass = 0.0;
            self.com = [0.0; 3];
            let oct = self.bbox.octant(old_pos);
            let sub = self.bbox.sub_box(oct);
            let child = self.children[oct].get_or_insert_with(|| Box::new(OctreeNode::new(sub)));
            child.insert(old_idx, old_pos, old_mass, max_depth - 1);
            self.total_mass = old_mass;
            self.com = old_pos;
        }
        let oct = self.bbox.octant(pos);
        let sub = self.bbox.sub_box(oct);
        let child = self.children[oct].get_or_insert_with(|| Box::new(OctreeNode::new(sub)));
        child.insert(idx, pos, mass, max_depth - 1);
        let tm = self.total_mass + mass;
        if tm > 0.0 {
            for d in 0..3 {
                self.com[d] = (self.com[d] * self.total_mass + pos[d] * mass) / tm;
            }
        }
        self.total_mass = tm;
    }
    /// Compute the gravitational acceleration at position `pos` using the
    /// Barnes-Hut approximation with opening angle `theta`.
    pub fn gravity_at(&self, pos: [f64; 3], theta: f64, softening: f64, g: f64) -> [f64; 3] {
        if self.total_mass == 0.0 {
            return [0.0; 3];
        }
        let dx = [
            self.com[0] - pos[0],
            self.com[1] - pos[1],
            self.com[2] - pos[2],
        ];
        let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + softening * softening;
        let r = r2.sqrt();
        let width = 2.0 * self.bbox.half_width();
        let is_leaf = self.children.iter().all(|c| c.is_none());
        if is_leaf || (width / r < theta) {
            let inv_r3 = g * self.total_mass / (r2 * r);
            return [dx[0] * inv_r3, dx[1] * inv_r3, dx[2] * inv_r3];
        }
        let mut acc = [0.0; 3];
        for c in self.children.iter().flatten() {
            let a = c.gravity_at(pos, theta, softening, g);
            acc[0] += a[0];
            acc[1] += a[1];
            acc[2] += a[2];
        }
        acc
    }
}
/// Artificial viscosity for SPH shock capturing.
///
/// Implements the standard Monaghan viscosity with Balsara limiter.
#[derive(Debug, Clone, Copy)]
pub struct ArtificialViscosity {
    /// Linear coefficient alpha (typically ~1.0).
    pub alpha: f64,
    /// Quadratic coefficient beta (typically ~2.0).
    pub beta: f64,
    /// Signal velocity coefficient for time stepping.
    pub _eta: f64,
}
impl ArtificialViscosity {
    /// Create standard artificial viscosity (alpha=1, beta=2).
    pub fn standard() -> Self {
        Self {
            alpha: 1.0,
            beta: 2.0,
            _eta: 0.01,
        }
    }
    /// Create with custom coefficients.
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self {
            alpha,
            beta,
            _eta: 0.01,
        }
    }
    /// Compute the viscosity term Pi_ij between particles i and j.
    ///
    /// Uses the Monaghan (1992) formulation:
    /// Pi_ij = (-alpha * c_bar * mu_ij + beta * mu_ij^2) / rho_bar
    /// where mu_ij = h_bar * v_ij . r_ij / (|r_ij|^2 + eta^2)
    pub fn compute_pi(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        vi: [f64; 3],
        vj: [f64; 3],
        rho_i: f64,
        rho_j: f64,
        cs_i: f64,
        cs_j: f64,
        hi: f64,
        hj: f64,
    ) -> f64 {
        let rij = sub3(ri, rj);
        let vij = sub3(vi, vj);
        let v_dot_r = dot3(vij, rij);
        if v_dot_r >= 0.0 {
            return 0.0;
        }
        let r2 = dot3(rij, rij);
        let h_bar = 0.5 * (hi + hj);
        let eta2 = 0.01 * h_bar * h_bar;
        let mu = h_bar * v_dot_r / (r2 + eta2);
        let c_bar = 0.5 * (cs_i + cs_j);
        let rho_bar = 0.5 * (rho_i + rho_j);
        if rho_bar < 1e-30 {
            return 0.0;
        }
        (-self.alpha * c_bar * mu + self.beta * mu * mu) / rho_bar
    }
    /// Von Neumann-Richtmyer viscosity (1D): q = C^2 * rho * (du/dx)^2 * dx^2.
    pub fn von_neumann_richtmyer(&self, density: f64, velocity_divergence: f64, dx: f64) -> f64 {
        if velocity_divergence >= 0.0 {
            return 0.0;
        }
        self.beta * density * velocity_divergence * velocity_divergence * dx * dx
    }
    /// Compute the Balsara switch for a particle.
    ///
    /// f_i = |div v|_i / (|div v|_i + |curl v|_i + eta * c_s / h)
    pub fn balsara_switch(&self, div_v: f64, curl_v_mag: f64, sound_speed: f64, h: f64) -> f64 {
        let eta_cs_h = 1e-4 * sound_speed / h;
        let denom = div_v.abs() + curl_v_mag + eta_cs_h;
        if denom < 1e-30 {
            return 0.0;
        }
        div_v.abs() / denom
    }
    /// Apply Balsara-limited viscosity.
    pub fn compute_pi_balsara(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        vi: [f64; 3],
        vj: [f64; 3],
        rho_i: f64,
        rho_j: f64,
        cs_i: f64,
        cs_j: f64,
        hi: f64,
        hj: f64,
        fi: f64,
        fj: f64,
    ) -> f64 {
        let pi_ij = self.compute_pi(ri, rj, vi, vj, rho_i, rho_j, cs_i, cs_j, hi, hj);
        let f_bar = 0.5 * (fi + fj);
        pi_ij * f_bar
    }
    /// Signal velocity for time step estimation.
    pub fn signal_velocity(&self, cs_i: f64, cs_j: f64, v_approach: f64) -> f64 {
        0.5 * (cs_i + cs_j) + self.beta * v_approach.abs()
    }
}
/// Astrophysical analysis tools: virial theorem, energy conservation.
#[derive(Debug, Clone)]
pub struct AstroAnalysis {
    /// Gravitational constant.
    pub g_const: f64,
    /// Softening length.
    pub softening: f64,
}
impl AstroAnalysis {
    /// Create a new analysis tool.
    pub fn new(g_const: f64, softening: f64) -> Self {
        Self { g_const, softening }
    }
    /// Total kinetic energy of the system.
    pub fn total_kinetic_energy(&self, particles: &[AstroParticle]) -> f64 {
        particles.iter().map(|p| p.kinetic_energy()).sum()
    }
    /// Total thermal energy.
    pub fn total_thermal_energy(&self, particles: &[AstroParticle]) -> f64 {
        particles.iter().map(|p| p.thermal_energy()).sum()
    }
    /// Total gravitational potential energy.
    pub fn total_gravitational_energy(&self, particles: &[AstroParticle]) -> f64 {
        let grav = GravitationalSph::new(self.g_const, self.softening);
        grav.gravitational_potential_energy(particles)
    }
    /// Total energy: E = E_kin + E_therm + E_grav.
    pub fn total_energy(&self, particles: &[AstroParticle]) -> f64 {
        self.total_kinetic_energy(particles)
            + self.total_thermal_energy(particles)
            + self.total_gravitational_energy(particles)
    }
    /// Virial ratio: 2*E_kin / |E_grav|.
    ///
    /// For a virialized system, this should be approximately 1.
    pub fn virial_ratio(&self, particles: &[AstroParticle]) -> f64 {
        let e_kin = self.total_kinetic_energy(particles);
        let e_grav = self.total_gravitational_energy(particles);
        if e_grav.abs() < 1e-30 {
            return 0.0;
        }
        2.0 * e_kin / e_grav.abs()
    }
    /// Virial parameter alpha_vir = 2*E_kin / |E_grav| (same as virial_ratio).
    pub fn virial_parameter(&self, particles: &[AstroParticle]) -> f64 {
        self.virial_ratio(particles)
    }
    /// Center of mass position.
    pub fn center_of_mass(&self, particles: &[AstroParticle]) -> [f64; 3] {
        let total_mass: f64 = particles.iter().map(|p| p.mass).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com = [0.0; 3];
        for p in particles {
            for d in 0..3 {
                com[d] += p.mass * p.pos[d];
            }
        }
        for d in 0..3 {
            com[d] /= total_mass;
        }
        com
    }
    /// Center of mass velocity.
    pub fn center_of_mass_velocity(&self, particles: &[AstroParticle]) -> [f64; 3] {
        let total_mass: f64 = particles.iter().map(|p| p.mass).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut com_vel = [0.0; 3];
        for p in particles {
            for d in 0..3 {
                com_vel[d] += p.mass * p.vel[d];
            }
        }
        for d in 0..3 {
            com_vel[d] /= total_mass;
        }
        com_vel
    }
    /// Total angular momentum.
    pub fn total_angular_momentum(&self, particles: &[AstroParticle]) -> [f64; 3] {
        let mut l = [0.0; 3];
        for p in particles {
            let li = cross3(p.pos, p.vel);
            l[0] += p.mass * li[0];
            l[1] += p.mass * li[1];
            l[2] += p.mass * li[2];
        }
        l
    }
    /// Total mass.
    pub fn total_mass(&self, particles: &[AstroParticle]) -> f64 {
        particles.iter().map(|p| p.mass).sum()
    }
    /// Half-mass radius: radius containing half the total mass.
    pub fn half_mass_radius(&self, particles: &[AstroParticle]) -> f64 {
        let com = self.center_of_mass(particles);
        let total_mass = self.total_mass(particles);
        let half_mass = total_mass / 2.0;
        let mut radii: Vec<(f64, f64)> = particles
            .iter()
            .map(|p| (dist3(p.pos, com), p.mass))
            .collect();
        radii.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut cumulative = 0.0;
        for (r, m) in &radii {
            cumulative += m;
            if cumulative >= half_mass {
                return *r;
            }
        }
        radii.last().map_or(0.0, |(r, _)| *r)
    }
    /// Density profile: compute average density in radial bins.
    pub fn density_profile(
        &self,
        particles: &[AstroParticle],
        n_bins: usize,
        r_max: f64,
    ) -> Vec<(f64, f64)> {
        let com = self.center_of_mass(particles);
        let dr = r_max / n_bins as f64;
        let mut bins = vec![0.0; n_bins];
        for p in particles {
            let r = dist3(p.pos, com);
            let bin = (r / dr) as usize;
            if bin < n_bins {
                bins[bin] += p.mass;
            }
        }
        (0..n_bins)
            .map(|i| {
                let r_inner = i as f64 * dr;
                let r_outer = (i + 1) as f64 * dr;
                let r_mid = 0.5 * (r_inner + r_outer);
                let volume = 4.0 / 3.0 * PI * (r_outer.powi(3) - r_inner.powi(3));
                let density = if volume > 1e-30 {
                    bins[i] / volume
                } else {
                    0.0
                };
                (r_mid, density)
            })
            .collect()
    }
    /// Energy conservation check: returns relative energy change.
    pub fn energy_conservation(&self, e_initial: f64, e_current: f64) -> f64 {
        if e_initial.abs() < 1e-30 {
            return 0.0;
        }
        (e_current - e_initial).abs() / e_initial.abs()
    }
}
/// Disk dynamics for Keplerian rotation and stability analysis.
///
/// Provides tools for analyzing accretion disks and galactic disks.
#[derive(Debug, Clone)]
pub struct DiskDynamics {
    /// Central mass (for Keplerian rotation).
    pub central_mass: f64,
    /// Gravitational constant.
    pub g_const: f64,
}
impl DiskDynamics {
    /// Create a new disk dynamics analyzer.
    pub fn new(central_mass: f64, g_const: f64) -> Self {
        Self {
            central_mass,
            g_const,
        }
    }
    /// Keplerian orbital velocity at radius R: v_K = sqrt(G*M/R).
    pub fn keplerian_velocity(&self, radius: f64) -> f64 {
        if radius <= 0.0 {
            return 0.0;
        }
        (self.g_const * self.central_mass / radius).sqrt()
    }
    /// Keplerian angular velocity: Omega_K = sqrt(G*M/R^3).
    pub fn keplerian_omega(&self, radius: f64) -> f64 {
        if radius <= 0.0 {
            return 0.0;
        }
        (self.g_const * self.central_mass / (radius * radius * radius)).sqrt()
    }
    /// Keplerian orbital period: T = 2*pi / Omega_K.
    pub fn orbital_period(&self, radius: f64) -> f64 {
        let omega = self.keplerian_omega(radius);
        if omega > 1e-30 {
            2.0 * PI / omega
        } else {
            f64::MAX
        }
    }
    /// Epicyclic frequency kappa for a Keplerian disk: kappa = Omega_K.
    pub fn epicyclic_frequency(&self, radius: f64) -> f64 {
        self.keplerian_omega(radius)
    }
    /// Toomre Q parameter: Q = c_s * kappa / (pi * G * Sigma).
    ///
    /// Disk is stable if Q > 1, unstable if Q < 1.
    pub fn toomre_q(&self, sound_speed: f64, surface_density: f64, radius: f64) -> f64 {
        let kappa = self.epicyclic_frequency(radius);
        if surface_density.abs() < 1e-30 {
            return f64::MAX;
        }
        sound_speed * kappa / (PI * self.g_const * surface_density)
    }
    /// Check if the disk is Toomre-stable at a given radius.
    pub fn is_toomre_stable(&self, sound_speed: f64, surface_density: f64, radius: f64) -> bool {
        self.toomre_q(sound_speed, surface_density, radius) > 1.0
    }
    /// Disk scale height: H = c_s / Omega_K.
    pub fn scale_height(&self, sound_speed: f64, radius: f64) -> f64 {
        let omega = self.keplerian_omega(radius);
        if omega > 1e-30 {
            sound_speed / omega
        } else {
            0.0
        }
    }
    /// Viscous accretion time: t_visc = R^2 / nu.
    pub fn viscous_time(&self, radius: f64, viscosity: f64) -> f64 {
        if viscosity > 1e-30 {
            radius * radius / viscosity
        } else {
            f64::MAX
        }
    }
    /// Shakura-Sunyaev alpha viscosity: nu = alpha * c_s * H.
    pub fn alpha_viscosity(&self, alpha_ss: f64, sound_speed: f64, radius: f64) -> f64 {
        let h = self.scale_height(sound_speed, radius);
        alpha_ss * sound_speed * h
    }
    /// Set up a Keplerian disk of particles.
    ///
    /// Creates particles distributed in rings with Keplerian velocities.
    pub fn create_keplerian_disk(
        &self,
        r_inner: f64,
        r_outer: f64,
        n_particles: usize,
        particle_mass: f64,
        internal_energy: f64,
        h: f64,
    ) -> Vec<AstroParticle> {
        let mut particles = Vec::with_capacity(n_particles);
        let mut rng = rand::rng();
        for _ in 0..n_particles {
            let u: f64 = rand::RngExt::random_range(&mut rng, 0.0..1.0);
            let r = (r_inner * r_inner + u * (r_outer * r_outer - r_inner * r_inner)).sqrt();
            let theta: f64 = rand::RngExt::random_range(&mut rng, 0.0..(2.0 * PI));
            let pos = [r * theta.cos(), r * theta.sin(), 0.0];
            let v_k = self.keplerian_velocity(r);
            let vel = [-v_k * theta.sin(), v_k * theta.cos(), 0.0];
            particles.push(AstroParticle::new(
                particle_mass,
                pos,
                vel,
                internal_energy,
                h,
            ));
        }
        particles
    }
    /// Specific angular momentum at radius R: j = R * v_K.
    pub fn specific_angular_momentum(&self, radius: f64) -> f64 {
        radius * self.keplerian_velocity(radius)
    }
    /// ISCO (innermost stable circular orbit) for Schwarzschild: R_ISCO = 6 * G * M / c^2.
    ///
    /// Uses a simplified Newtonian analog; `speed_of_light` must be in code units.
    pub fn isco_radius(&self, speed_of_light: f64) -> f64 {
        6.0 * self.g_const * self.central_mass / (speed_of_light * speed_of_light)
    }
}
/// A single astrophysical SPH particle.
///
/// Stores mass, position, velocity, thermodynamic state, and smoothing length.
#[derive(Debug, Clone)]
pub struct AstroParticle {
    /// Particle mass \[code units\].
    pub mass: f64,
    /// Position vector \[x, y, z\].
    pub pos: [f64; 3],
    /// Velocity vector \[vx, vy, vz\].
    pub vel: [f64; 3],
    /// SPH density.
    pub density: f64,
    /// Specific internal energy u (energy per unit mass).
    pub internal_energy: f64,
    /// Smoothing length h.
    pub h: f64,
    /// Pressure (computed from EOS).
    pub pressure: f64,
    /// Acceleration vector.
    pub acc: [f64; 3],
    /// Rate of change of internal energy du/dt.
    pub du_dt: f64,
    /// Balsara switch value (0 to 1).
    pub balsara_switch: f64,
    /// Whether this particle is a sink particle.
    pub is_sink: bool,
    /// Accretion radius (for sink particles).
    pub accretion_radius: f64,
}
impl AstroParticle {
    /// Create a new astrophysical particle.
    pub fn new(mass: f64, pos: [f64; 3], vel: [f64; 3], internal_energy: f64, h: f64) -> Self {
        Self {
            mass,
            pos,
            vel,
            density: 0.0,
            internal_energy,
            h,
            pressure: 0.0,
            acc: [0.0; 3],
            du_dt: 0.0,
            balsara_switch: 1.0,
            is_sink: false,
            accretion_radius: 0.0,
        }
    }
    /// Create a sink particle at the given position.
    pub fn new_sink(mass: f64, pos: [f64; 3], vel: [f64; 3], accretion_radius: f64) -> Self {
        Self {
            mass,
            pos,
            vel,
            density: 0.0,
            internal_energy: 0.0,
            h: accretion_radius,
            pressure: 0.0,
            acc: [0.0; 3],
            du_dt: 0.0,
            balsara_switch: 0.0,
            is_sink: true,
            accretion_radius,
        }
    }
    /// Kinetic energy: 0.5 * m * v^2.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
    /// Thermal energy: m * u.
    pub fn thermal_energy(&self) -> f64 {
        self.mass * self.internal_energy
    }
    /// Specific angular momentum L = r x v.
    pub fn specific_angular_momentum(&self) -> [f64; 3] {
        cross3(self.pos, self.vel)
    }
    /// Cylindrical radius (distance from z-axis).
    pub fn cylindrical_radius(&self) -> f64 {
        (self.pos[0] * self.pos[0] + self.pos[1] * self.pos[1]).sqrt()
    }
    /// Speed (magnitude of velocity).
    pub fn speed(&self) -> f64 {
        norm3(self.vel)
    }
    /// Sound speed for ideal gas: c_s = sqrt(gamma * P / rho).
    pub fn sound_speed(&self, gamma: f64) -> f64 {
        if self.density > 1e-30 {
            (gamma * self.pressure / self.density).sqrt()
        } else {
            0.0
        }
    }
}
/// An axis-aligned bounding box in 3-D.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct BBox3 {
    /// Minimum corner.
    pub lo: [f64; 3],
    /// Maximum corner.
    pub hi: [f64; 3],
}
#[allow(dead_code)]
impl BBox3 {
    /// Create a new bounding box.
    pub fn new(lo: [f64; 3], hi: [f64; 3]) -> Self {
        Self { lo, hi }
    }
    /// Centre of the box.
    pub fn centre(&self) -> [f64; 3] {
        [
            0.5 * (self.lo[0] + self.hi[0]),
            0.5 * (self.lo[1] + self.hi[1]),
            0.5 * (self.lo[2] + self.hi[2]),
        ]
    }
    /// Half-width of the box (max of the three half-extents).
    pub fn half_width(&self) -> f64 {
        let dx = 0.5 * (self.hi[0] - self.lo[0]);
        let dy = 0.5 * (self.hi[1] - self.lo[1]);
        let dz = 0.5 * (self.hi[2] - self.lo[2]);
        dx.max(dy).max(dz)
    }
    /// Check whether a point lies inside.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        p[0] >= self.lo[0]
            && p[0] <= self.hi[0]
            && p[1] >= self.lo[1]
            && p[1] <= self.hi[1]
            && p[2] >= self.lo[2]
            && p[2] <= self.hi[2]
    }
    /// Return the octant index (0..7) for a point relative to the centre.
    pub fn octant(&self, p: [f64; 3]) -> usize {
        let c = self.centre();
        let mut idx = 0;
        if p[0] >= c[0] {
            idx |= 1;
        }
        if p[1] >= c[1] {
            idx |= 2;
        }
        if p[2] >= c[2] {
            idx |= 4;
        }
        idx
    }
    /// Return the sub-box for a given octant index.
    pub fn sub_box(&self, octant: usize) -> Self {
        let c = self.centre();
        let mut lo = self.lo;
        let mut hi = self.hi;
        if octant & 1 != 0 {
            lo[0] = c[0];
        } else {
            hi[0] = c[0];
        }
        if octant & 2 != 0 {
            lo[1] = c[1];
        } else {
            hi[1] = c[1];
        }
        if octant & 4 != 0 {
            lo[2] = c[2];
        } else {
            hi[2] = c[2];
        }
        Self { lo, hi }
    }
}
