//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    ReactionFn, cubic_kernel, cubic_kernel_grad, dot3, heat_exchange_rate, len3, scale3, sub3,
};

/// Advection-dominated scalar transport with stabilisation.
///
/// Adds an upwind or flux-limited term to suppress oscillations.
#[derive(Debug, Clone)]
pub struct AdvectionScheme {
    /// The advection scheme type.
    pub kind: AdvectionSchemeKind,
    /// Numerical diffusion coefficient for upwind scheme.
    pub epsilon_upwind: f64,
}
impl AdvectionScheme {
    /// Create a new advection scheme.
    pub fn new(kind: AdvectionSchemeKind, epsilon_upwind: f64) -> Self {
        Self {
            kind,
            epsilon_upwind,
        }
    }
    /// Compute the stabilised advection rate dφ/dt for particle `i`.
    ///
    /// For upwind: adds ε·∑_j (m_j/ρ_j)·|v·∇W|·(φ_j − φ_i) if v·r_ij < 0.
    pub fn compute_rate_i(&self, particles: &[TransportParticle], i: usize) -> f64 {
        let pi = &particles[i];
        let mut rate = 0.0f64;
        for (j, pj) in particles.iter().enumerate() {
            if i == j {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r = len3(r_ij);
            if r < 1e-12 {
                continue;
            }
            let h_ij = 0.5 * (pi.h + pj.h);
            let grad_w = cubic_kernel_grad(r_ij, h_ij);
            let vj = pj.volume();
            let v_dot_r = dot3(pi.velocity, r_ij);
            match self.kind {
                AdvectionSchemeKind::Standard => {
                    let v_dot_gradw = dot3(pi.velocity, grad_w);
                    rate -= vj * (pj.conc() - pi.conc()) * v_dot_gradw;
                }
                AdvectionSchemeKind::Upwind => {
                    let v_dot_gradw = dot3(pi.velocity, grad_w);
                    rate -= vj * (pj.conc() - pi.conc()) * v_dot_gradw;
                    if v_dot_r < 0.0 {
                        let e_ij = scale3(r_ij, 1.0 / r);
                        let term = 2.0 * (pj.conc() - pi.conc()) / r * dot3(e_ij, grad_w);
                        rate += self.epsilon_upwind * vj * term;
                    }
                }
                AdvectionSchemeKind::WenoSph => {
                    let v_dot_gradw = dot3(pi.velocity, grad_w);
                    let dphi = pj.conc() - pi.conc();
                    let r_ratio = if dphi.abs() > 1e-15 {
                        dphi / (dphi + 1e-6 * dphi.signum())
                    } else {
                        1.0
                    };
                    let phi_lim = 2.0 * r_ratio / (1.0 + r_ratio);
                    rate -= vj * dphi * phi_lim * v_dot_gradw;
                }
            }
        }
        rate
    }
    /// Compute rates for all particles.
    pub fn compute_all_rates(&self, particles: &[TransportParticle]) -> Vec<f64> {
        (0..particles.len())
            .map(|i| self.compute_rate_i(particles, i))
            .collect()
    }
    /// Apply one explicit Euler step.
    pub fn step(&self, particles: &mut [TransportParticle], dt: f64) {
        let rates = self.compute_all_rates(particles);
        for (p, &rate) in particles.iter_mut().zip(rates.iter()) {
            if let Some(c) = p.concentration.first_mut() {
                *c += rate * dt;
            }
        }
    }
}
/// Multi-species solute transport with cross-diffusion.
///
/// For species s: ∂c_s/∂t = ∑_{s'} D_{ss'} ∇²c_{s'} + advection
#[derive(Debug, Clone)]
pub struct SoluteTransport {
    /// Number of species.
    pub n_species: usize,
    /// Diffusion matrix D_{ss'} (flat row-major, n_species × n_species).
    pub diffusion_matrix: Vec<f64>,
}
impl SoluteTransport {
    /// Create a new multi-species solute transport.
    ///
    /// `diffusion_matrix` is row-major n_species × n_species.
    pub fn new(n_species: usize, diffusion_matrix: Vec<f64>) -> Self {
        assert_eq!(diffusion_matrix.len(), n_species * n_species);
        Self {
            n_species,
            diffusion_matrix,
        }
    }
    /// Create with diagonal diffusion (no cross-diffusion).
    pub fn diagonal(diffusivities: &[f64]) -> Self {
        let n = diffusivities.len();
        let mut mat = vec![0.0f64; n * n];
        for (i, &d) in diffusivities.iter().enumerate() {
            mat[i * n + i] = d;
        }
        Self {
            n_species: n,
            diffusion_matrix: mat,
        }
    }
    /// Get diffusion coefficient D_{ss'}.
    pub fn d(&self, s: usize, sp: usize) -> f64 {
        if s < self.n_species && sp < self.n_species {
            self.diffusion_matrix[s * self.n_species + sp]
        } else {
            0.0
        }
    }
    /// Compute dc_s/dt for species `s` for all particles using SPH Laplacian.
    pub fn compute_species_rate(&self, particles: &[TransportParticle], s: usize) -> Vec<f64> {
        let n = particles.len();
        let mut rate = vec![0.0f64; n];
        for sp in 0..self.n_species {
            let d_ssp = self.d(s, sp);
            if d_ssp.abs() < 1e-300 {
                continue;
            }
            let diff_op = DiffusionSph::new(d_ssp);
            let mut tmp: Vec<TransportParticle> = particles.to_vec();
            for (i, p) in tmp.iter_mut().enumerate() {
                let c_sp = if sp < particles[i].concentration.len() {
                    particles[i].concentration[sp]
                } else {
                    0.0
                };
                if p.concentration.is_empty() {
                    p.concentration.push(c_sp);
                } else {
                    p.concentration[0] = c_sp;
                }
            }
            let lap = diff_op.compute_dphidt(&tmp);
            for i in 0..n {
                rate[i] += lap[i];
            }
        }
        rate
    }
    /// Compute dc/dt for all species and all particles.
    ///
    /// Returns a `Vec<Vec`f64`>` of shape \[n_species\]\[n_particles\].
    pub fn compute_all_rates(&self, particles: &[TransportParticle]) -> Vec<Vec<f64>> {
        (0..self.n_species)
            .map(|s| self.compute_species_rate(particles, s))
            .collect()
    }
    /// Apply one explicit Euler step for all species.
    pub fn step(&self, particles: &mut [TransportParticle], dt: f64) {
        let rates = self.compute_all_rates(particles);
        for p in particles.iter_mut() {
            while p.concentration.len() < self.n_species {
                p.concentration.push(0.0);
            }
        }
        for (s, rate_s) in rates.iter().enumerate().take(self.n_species) {
            for (i, p) in particles.iter_mut().enumerate() {
                if i < rate_s.len() {
                    p.concentration[s] += rate_s[i] * dt;
                }
            }
        }
    }
}
/// Chemical reaction model for source/sink terms in transport equations.
///
/// Implements Arrhenius kinetics and first-/second-order reaction terms.
#[derive(Debug, Clone)]
pub struct ChemicalReaction {
    /// Pre-exponential factor A \[1/s or m³/(mol·s) for second-order\].
    pub pre_exponential: f64,
    /// Activation energy E_a \[J/mol\].
    pub activation_energy: f64,
    /// Universal gas constant R \[J/(mol·K)\].
    pub gas_constant: f64,
    /// Reaction order (1 = first-order, 2 = second-order bimolecular).
    pub order: u32,
    /// Stoichiometric coefficient (positive = production, negative = consumption).
    pub stoich: f64,
}
impl ChemicalReaction {
    /// Create a first-order Arrhenius reaction.
    pub fn first_order(pre_exponential: f64, activation_energy: f64) -> Self {
        Self {
            pre_exponential,
            activation_energy,
            gas_constant: 8.314_462_618,
            order: 1,
            stoich: -1.0,
        }
    }
    /// Create a second-order bimolecular reaction.
    pub fn second_order(pre_exponential: f64, activation_energy: f64) -> Self {
        Self {
            pre_exponential,
            activation_energy,
            gas_constant: 8.314_462_618,
            order: 2,
            stoich: -1.0,
        }
    }
    /// Arrhenius rate constant k(T) = A exp(−Ea / (R T)).
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        if temperature < 1e-300 {
            return 0.0;
        }
        self.pre_exponential * (-self.activation_energy / (self.gas_constant * temperature)).exp()
    }
    /// Source/sink term R(c, T) = stoich · k(T) · c^order \[mol/(m³·s)\].
    pub fn source_term(&self, concentration: f64, temperature: f64) -> f64 {
        let k = self.rate_constant(temperature);
        let c = concentration.max(0.0);
        self.stoich * k * c.powi(self.order as i32)
    }
    /// Integrate concentration forward by Euler: c ← c + Δt · R(c, T).
    pub fn euler_step(&self, concentration: f64, temperature: f64, dt: f64) -> f64 {
        (concentration + dt * self.source_term(concentration, temperature)).max(0.0)
    }
    /// Half-life t_{1/2} = ln(2) / k(T) (first-order only).
    pub fn half_life(&self, temperature: f64) -> f64 {
        let k = self.rate_constant(temperature);
        if k < 1e-300 {
            return f64::INFINITY;
        }
        2.0_f64.ln() / k
    }
}
/// SPH electroosmotic flow model based on the Helmholtz–Smoluchowski equation.
///
/// Models the electroosmotic velocity `v_eo = -ε ζ E / η` and Debye screening.
#[derive(Debug, Clone)]
pub struct SPHElectroosmosis {
    /// Zeta potential ζ \[V\] (negative for most surfaces in water).
    pub zeta_potential: f64,
    /// Dynamic viscosity η \[Pa·s\].
    pub viscosity: f64,
    /// Permittivity ε \[F/m\].
    pub permittivity: f64,
    /// Debye screening length λ_D \[m\].
    pub debye_length: f64,
}
impl SPHElectroosmosis {
    /// Create a new `SPHElectroosmosis` model.
    pub fn new(zeta_potential: f64, viscosity: f64, permittivity: f64, debye_length: f64) -> Self {
        Self {
            zeta_potential,
            viscosity,
            permittivity,
            debye_length,
        }
    }
    /// Helmholtz–Smoluchowski electroosmotic velocity.
    ///
    /// `v_eo = -ε ζ E / η`
    ///
    /// # Arguments
    /// * `e_field` – Applied electric field magnitude \[V/m\].
    pub fn smoluchowski_velocity(&self, e_field: f64) -> f64 {
        if self.viscosity <= 0.0_f64 {
            return 0.0_f64;
        }
        -self.permittivity * self.zeta_potential * e_field / self.viscosity
    }
    /// Electroosmotic mobility µ_eo = -ε ζ / η \[m²/(V·s)\].
    pub fn electroosmotic_mobility(&self) -> f64 {
        if self.viscosity <= 0.0_f64 {
            return 0.0_f64;
        }
        -self.permittivity * self.zeta_potential / self.viscosity
    }
    /// Debye–Hückel exponential profile of electric potential at distance y.
    ///
    /// `φ(y) = ζ exp(-y / λ_D)`
    pub fn potential_profile(&self, y: f64) -> f64 {
        if self.debye_length <= 0.0_f64 {
            return 0.0_f64;
        }
        self.zeta_potential * (-y / self.debye_length).exp()
    }
    /// Body force on a particle with charge density `rho_e` in electric field `E`.
    ///
    /// `f = rho_e E`
    pub fn body_force(&self, rho_e: f64, e_field: f64) -> f64 {
        rho_e * e_field
    }
}
/// Verifies global conservation of mass and energy in a transport step.
#[derive(Debug, Clone)]
pub struct ConservationCheck {
    /// Reference total scalar quantity (recorded at start of step).
    pub ref_scalar: f64,
    /// Reference total thermal energy (recorded at start of step).
    pub ref_thermal: f64,
}
impl ConservationCheck {
    /// Create a new conservation checker, recording initial state.
    pub fn new(particles: &[TransportParticle], rho_cv: f64) -> Self {
        let ref_scalar: f64 = particles.iter().map(|p| p.volume() * p.conc()).sum();
        let ref_thermal: f64 = particles
            .iter()
            .map(|p| rho_cv * p.volume() * p.temperature)
            .sum();
        Self {
            ref_scalar,
            ref_thermal,
        }
    }
    /// Compute current total scalar.
    pub fn current_scalar(particles: &[TransportParticle]) -> f64 {
        particles.iter().map(|p| p.volume() * p.conc()).sum()
    }
    /// Compute current total thermal energy.
    pub fn current_thermal(particles: &[TransportParticle], rho_cv: f64) -> f64 {
        particles
            .iter()
            .map(|p| rho_cv * p.volume() * p.temperature)
            .sum()
    }
    /// Relative scalar conservation error: |Q - Q0| / |Q0|.
    pub fn scalar_error(&self, particles: &[TransportParticle]) -> f64 {
        let q = Self::current_scalar(particles);
        if self.ref_scalar.abs() < 1e-300 {
            q.abs()
        } else {
            (q - self.ref_scalar).abs() / self.ref_scalar.abs()
        }
    }
    /// Relative thermal energy conservation error.
    pub fn thermal_error(&self, particles: &[TransportParticle], rho_cv: f64) -> f64 {
        let e = Self::current_thermal(particles, rho_cv);
        if self.ref_thermal.abs() < 1e-300 {
            e.abs()
        } else {
            (e - self.ref_thermal).abs() / self.ref_thermal.abs()
        }
    }
    /// Check that scalar is conserved to within `tol` (relative).
    pub fn is_scalar_conserved(&self, particles: &[TransportParticle], tol: f64) -> bool {
        self.scalar_error(particles) < tol
    }
}
/// SPH heat transfer solver: ρ cₚ DT/Dt = ∇·(k ∇T) + Q.
#[derive(Debug, Clone)]
pub struct HeatTransferSPH {
    /// Thermal conductivity k \[W/(m·K)\].
    pub thermal_conductivity: f64,
    /// Specific heat capacity cₚ \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Volumetric heat source Q \[W/m³\].
    pub heat_source: f64,
    /// Reference temperature \[K\].
    pub t_ref: f64,
}
impl HeatTransferSPH {
    /// Create a new heat transfer SPH model.
    pub fn new(thermal_conductivity: f64, specific_heat: f64) -> Self {
        Self {
            thermal_conductivity,
            specific_heat,
            heat_source: 0.0,
            t_ref: 293.15,
        }
    }
    /// Thermal diffusivity α = k / (ρ cₚ).
    pub fn alpha(&self, density: f64) -> f64 {
        if density < 1e-300 || self.specific_heat < 1e-300 {
            return 0.0;
        }
        self.thermal_conductivity / (density * self.specific_heat)
    }
    /// SPH heat conduction rate DT/Dt at particle `i`.
    ///
    /// Uses the SPH approximation:
    ///   ρ_i cₚ DT_i/Dt = Σ_j m_j 4 k_i k_j/(k_i+k_j) (T_j−T_i)/(r²_{ij}) (r·∇W)
    pub fn heat_conduction_rate(
        &self,
        particles: &[TransportParticle],
        i: usize,
        kernel_grad_fn: impl Fn([f64; 3], f64) -> [f64; 3],
    ) -> f64 {
        let pi = &particles[i];
        let ki = self.thermal_conductivity;
        let rhocpi = pi.density.max(1e-300) * self.specific_heat.max(1e-300);
        let mut heat_rate = 0.0_f64;
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r = len3(r_ij);
            if r < 1e-12 {
                continue;
            }
            let grad_w = kernel_grad_fn(r_ij, pi.h);
            let rdotgw = dot3(r_ij, grad_w);
            let kj = ki;
            let k_ij = 4.0 * ki * kj / (ki + kj + 1e-300);
            heat_rate += pj.mass / pj.density.max(1e-300)
                * k_ij
                * (pj.temperature - pi.temperature)
                * 2.0
                * rdotgw
                / (r * r + 1e-20);
        }
        heat_rate / rhocpi + self.heat_source / rhocpi
    }
    /// Advance temperature of all particles by Euler integration.
    pub fn euler_step(&self, particles: &mut [TransportParticle], rates: &[f64], dt: f64) {
        for (p, &rate) in particles.iter_mut().zip(rates.iter()) {
            p.temperature += rate * dt;
        }
    }
    /// Nusselt number for forced convection over a cylinder (Churchill–Bernstein).
    ///
    /// Nu = 0.3 + 0.62 Re^{1/2} Pr^{1/3} / (1 + (0.4/Pr)^{2/3})^{1/4}
    ///          · (1 + (Re/282000)^{5/8})^{4/5}
    pub fn nusselt_cylinder(reynolds: f64, prandtl: f64) -> f64 {
        if prandtl < 1e-300 || reynolds < 0.0 {
            return 0.3;
        }
        let t1 = 0.62 * reynolds.sqrt() * prandtl.powf(1.0 / 3.0);
        let denom = (1.0 + (0.4 / prandtl).powf(2.0 / 3.0)).powf(0.25);
        let t2 = 1.0 + (reynolds / 282000.0).powf(5.0 / 8.0);
        0.3 + (t1 / denom) * t2.powf(4.0 / 5.0)
    }
    /// Prandtl number Pr = cₚ μ / k.
    pub fn prandtl_number(&self, dynamic_viscosity: f64) -> f64 {
        if self.thermal_conductivity < 1e-300 {
            return 0.0;
        }
        self.specific_heat * dynamic_viscosity / self.thermal_conductivity
    }
}
/// SPH heat transfer using Fourier's law of conduction.
///
/// Computes inter-particle heat exchange: `dQ = 2 k (T_j - T_i) m_j / ρ_j · ∇W`.
#[derive(Debug, Clone)]
pub struct SPHHeatTransfer {
    /// Thermal conductivity k \[W/(m·K)\].
    pub thermal_conductivity: f64,
    /// Specific heat capacity cₚ \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Particle density ρ \[kg/m³\].
    pub density: f64,
}
impl SPHHeatTransfer {
    /// Create a new `SPHHeatTransfer` model.
    pub fn new(thermal_conductivity: f64, specific_heat: f64, density: f64) -> Self {
        Self {
            thermal_conductivity,
            specific_heat,
            density,
        }
    }
    /// Thermal diffusivity α = k / (ρ cₚ).
    pub fn thermal_diffusivity(&self) -> f64 {
        let rho_cp = self.density * self.specific_heat;
        if rho_cp <= 0.0_f64 {
            return 0.0_f64;
        }
        self.thermal_conductivity / rho_cp
    }
    /// Heat exchange from particle i to particle j.
    ///
    /// Simplified SPH conduction: `dQ_i = -2k (T_i - T_j) m_j W(r, h) / ρ_j`.
    ///
    /// # Arguments
    /// * `t_i`, `t_j` – Temperatures of particles i and j \[K\].
    /// * `mass_j`     – Mass of particle j \[kg\].
    /// * `rho_j`      – Density of particle j \[kg/m³\].
    /// * `w`          – Kernel weight W(r, h).
    pub fn heat_exchange(&self, t_i: f64, t_j: f64, mass_j: f64, rho_j: f64, w: f64) -> f64 {
        if rho_j <= 0.0_f64 {
            return 0.0_f64;
        }
        -2.0_f64 * self.thermal_conductivity * (t_i - t_j) * mass_j / rho_j * w
    }
    /// Rate of temperature change dT/dt for particle i from all neighbours.
    ///
    /// # Arguments
    /// * `t_i`         – Temperature of particle i \[K\].
    /// * `neighbours`  – Slice of `(T_j, r_ij, mass_j, rho_j)` tuples.
    /// * `h`           – Smoothing length.
    pub fn dtemp_dt(&self, t_i: f64, neighbours: &[(f64, [f64; 3], f64, f64)], h: f64) -> f64 {
        let alpha = self.thermal_diffusivity();
        let mut sum = 0.0_f64;
        for &(t_j, r_ij, mass_j, rho_j) in neighbours {
            let r = len3(r_ij);
            let w = cubic_kernel(r, h);
            sum += heat_exchange_rate(alpha, t_i, t_j, mass_j, rho_j, w);
        }
        sum
    }
}
/// SPH concentration diffusion via Fick's second law.
///
/// Implements `∂c/∂t = D ∇²c` using an SPH Laplacian approximation.
#[derive(Debug, Clone)]
pub struct SPHDiffusion {
    /// Diffusion coefficient D \[m²/s\].
    pub diffusivity: f64,
    /// Smoothing length h \[m\].
    pub smoothing_length: f64,
}
impl SPHDiffusion {
    /// Create a new `SPHDiffusion`.
    pub fn new(diffusivity: f64, smoothing_length: f64) -> Self {
        Self {
            diffusivity,
            smoothing_length,
        }
    }
    /// SPH Laplacian approximation: `∇²c_i ≈ Σ_j 2 (c_j - c_i) m_j W / ρ_j`.
    ///
    /// # Arguments
    /// * `c_i`        – Concentration at particle i.
    /// * `neighbours` – Slice of `(c_j, r_ij, mass_j, rho_j)`.
    pub fn concentration_laplacian(
        &self,
        c_i: f64,
        neighbours: &[(f64, [f64; 3], f64, f64)],
    ) -> f64 {
        let h = self.smoothing_length;
        let mut lap = 0.0_f64;
        for &(c_j, r_ij, mass_j, rho_j) in neighbours {
            let r = len3(r_ij);
            let w = cubic_kernel(r, h);
            if rho_j > 0.0_f64 {
                lap += 2.0_f64 * (c_j - c_i) * mass_j / rho_j * w;
            }
        }
        lap
    }
    /// Rate of change `dc/dt = D ∇²c`.
    pub fn dcdt(&self, c_i: f64, neighbours: &[(f64, [f64; 3], f64, f64)]) -> f64 {
        self.diffusivity * self.concentration_laplacian(c_i, neighbours)
    }
}
/// A Neumann flux boundary condition for transport equations.
///
/// Specifies dφ/dn = g at a boundary surface, implemented via
/// boundary ghost particles or direct flux injection.
#[derive(Debug, Clone)]
pub struct BoundaryFlux {
    /// Prescribed normal flux g = dφ/dn \[1/m × species units\].
    pub flux_value: f64,
    /// Outward normal direction of the boundary.
    pub normal: [f64; 3],
    /// List of particle indices subject to this boundary condition.
    pub particle_ids: Vec<usize>,
    /// Area associated with each boundary particle \[m²\].
    pub area_weights: Vec<f64>,
}
impl BoundaryFlux {
    /// Create a new Neumann flux boundary condition.
    pub fn new(
        flux_value: f64,
        normal: [f64; 3],
        particle_ids: Vec<usize>,
        area_weights: Vec<f64>,
    ) -> Self {
        Self {
            flux_value,
            normal,
            particle_ids,
            area_weights,
        }
    }
    /// Apply the boundary flux to particle concentrations.
    ///
    /// dφ_i/dt += g * A_i / V_i  (for each boundary particle i)
    pub fn apply(&self, particles: &mut [TransportParticle], dt: f64, diffusivity: f64) {
        for (idx, &pid) in self.particle_ids.iter().enumerate() {
            if pid >= particles.len() {
                continue;
            }
            let area = if idx < self.area_weights.len() {
                self.area_weights[idx]
            } else {
                1.0
            };
            let vol = particles[pid].volume();
            if vol.abs() < 1e-300 {
                continue;
            }
            let delta_phi = diffusivity * self.flux_value * area / vol * dt;
            if let Some(c) = particles[pid].concentration.first_mut() {
                *c += delta_phi;
            }
        }
    }
    /// Total flux injected: ∑ g * A_i.
    pub fn total_flux(&self) -> f64 {
        self.flux_value * self.area_weights.iter().sum::<f64>()
    }
}
/// A scalar concentration field stored per-particle.
///
/// Tracks molar concentration \[mol/m³\], mole fraction, and local fluid density.
#[derive(Debug, Clone)]
pub struct ConcentrationField {
    /// Molar concentration of each particle \[mol/m³\].
    pub molar_concentration: Vec<f64>,
    /// Mole fraction of each particle (0–1).
    pub mole_fraction: Vec<f64>,
    /// Local fluid density at each particle \[kg/m³\].
    pub fluid_density: Vec<f64>,
    /// Molar mass of the species \[kg/mol\].
    pub molar_mass: f64,
    /// Number of particles.
    pub n_particles: usize,
}
impl ConcentrationField {
    /// Create a new concentration field with `n` particles.
    ///
    /// Initial concentrations and fractions are set to zero.
    pub fn new(n_particles: usize, molar_mass: f64) -> Self {
        Self {
            molar_concentration: vec![0.0; n_particles],
            mole_fraction: vec![0.0; n_particles],
            fluid_density: vec![1000.0; n_particles],
            molar_mass,
            n_particles,
        }
    }
    /// Set the molar concentration of particle `i` and update mole fraction.
    ///
    /// `total_molar_conc` is the sum of all species concentrations at that particle.
    pub fn set_concentration(&mut self, i: usize, c: f64, total_molar_conc: f64) {
        self.molar_concentration[i] = c.max(0.0);
        if total_molar_conc > 1e-300 {
            self.mole_fraction[i] = c / total_molar_conc;
        } else {
            self.mole_fraction[i] = 0.0;
        }
    }
    /// Mass concentration \[kg/m³\] of particle `i`.
    pub fn mass_concentration(&self, i: usize) -> f64 {
        self.molar_concentration[i] * self.molar_mass
    }
    /// Mass fraction at particle `i`.
    pub fn mass_fraction(&self, i: usize) -> f64 {
        if self.fluid_density[i] > 1e-300 {
            self.mass_concentration(i) / self.fluid_density[i]
        } else {
            0.0
        }
    }
    /// Total amount of substance \[mol\] (trapezoidal integration with uniform spacing `dx`).
    pub fn total_moles(&self, dx: f64) -> f64 {
        self.molar_concentration.iter().sum::<f64>() * dx
    }
    /// Apply a uniform decay c_i ← c_i * exp(−λ Δt).
    pub fn apply_decay(&mut self, lambda: f64, dt: f64) {
        let factor = (-lambda * dt).exp();
        for c in &mut self.molar_concentration {
            *c *= factor;
        }
    }
}
/// Multi-species SPH transport with stoichiometric reactions.
///
/// Tracks `nspecies` species. Each species has its own diffusivity.
/// A single reversible reaction `Σ ν_α A_α → 0` proceeds with rate k_r.
#[derive(Debug, Clone)]
pub struct SPHMultispeciesTransport {
    /// Number of species.
    pub nspecies: usize,
    /// Diffusivities for each species \[m²/s\].
    pub diffusivities: Vec<f64>,
    /// Stoichiometric coefficients ν_α (negative = reactant, positive = product).
    pub stoichiometry: Vec<f64>,
    /// Molecular weights \[g/mol\].
    pub molecular_weights: Vec<f64>,
}
impl SPHMultispeciesTransport {
    /// Create a new `SPHMultispeciesTransport`.
    ///
    /// # Arguments
    /// * `diffusivities`  – Diffusivity for each species.
    /// * `stoichiometry`  – Stoichiometric coefficients (negative = reactant).
    /// * `molecular_weights` – Molecular weights \[g/mol\].
    pub fn new(
        diffusivities: Vec<f64>,
        stoichiometry: Vec<f64>,
        molecular_weights: Vec<f64>,
    ) -> Self {
        let nspecies = diffusivities.len();
        assert_eq!(stoichiometry.len(), nspecies);
        assert_eq!(molecular_weights.len(), nspecies);
        Self {
            nspecies,
            diffusivities,
            stoichiometry,
            molecular_weights,
        }
    }
    /// Reaction source term for species α at reaction rate k_r \[mol/(m³·s)\].
    ///
    /// `S_α = ν_α · k_r`
    pub fn reaction_source(&self, alpha: usize, k_r: f64) -> f64 {
        self.stoichiometry[alpha] * k_r
    }
    /// Mole fraction of species α given a concentration vector.
    ///
    /// Returns 0 when total molar concentration is zero.
    pub fn mole_fraction(&self, alpha: usize, concentrations: &[f64]) -> f64 {
        let total: f64 = concentrations.iter().sum();
        if total <= 0.0_f64 {
            return 0.0_f64;
        }
        concentrations[alpha] / total
    }
    /// Mass fraction of species α given a concentration vector.
    ///
    /// Uses molecular weights to convert mole fractions to mass fractions.
    pub fn mass_fraction(&self, alpha: usize, concentrations: &[f64]) -> f64 {
        let n = self.nspecies.min(concentrations.len());
        let mass_total: f64 = (0..n)
            .map(|a| concentrations[a] * self.molecular_weights[a])
            .sum();
        if mass_total <= 0.0_f64 {
            return 0.0_f64;
        }
        concentrations[alpha] * self.molecular_weights[alpha] / mass_total
    }
    /// SPH diffusion flux for species α at particle i.
    ///
    /// `dc_α/dt = D_α Σ_j 2(c_αj - c_αi) m_j W / ρ_j`
    pub fn dcdt(
        &self,
        alpha: usize,
        c_i: f64,
        neighbours: &[(f64, [f64; 3], f64, f64)],
        h: f64,
        source: f64,
    ) -> f64 {
        let mut lap = 0.0_f64;
        for &(c_j, r_ij, mass_j, rho_j) in neighbours {
            let r = len3(r_ij);
            let w = cubic_kernel(r, h);
            if rho_j > 0.0_f64 {
                lap += 2.0_f64 * (c_j - c_i) * mass_j / rho_j * w;
            }
        }
        self.diffusivities[alpha] * lap + source
    }
}
/// Advects a scalar field carried by SPH particles.
///
/// Uses the standard SPH advection: Dφ/Dt = 0 (material derivative).
/// In practice, updates particle positions and preserves φ values.
#[derive(Debug, Clone)]
pub struct ScalarTransport {
    /// Time step \[s\].
    pub dt: f64,
}
impl ScalarTransport {
    /// Create a new scalar transport operator.
    pub fn new(dt: f64) -> Self {
        Self { dt }
    }
    /// Advect particles: update positions by x += v * dt.
    ///
    /// Concentrations are carried passively (Lagrangian frame).
    pub fn advect(&self, particles: &mut [TransportParticle]) {
        for p in particles.iter_mut() {
            p.position[0] += p.velocity[0] * self.dt;
            p.position[1] += p.velocity[1] * self.dt;
            p.position[2] += p.velocity[2] * self.dt;
        }
    }
    /// Compute the SPH interpolated concentration at a probe point `x`.
    pub fn interpolate_at(&self, particles: &[TransportParticle], x: [f64; 3]) -> f64 {
        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for p in particles {
            let r_ij = sub3(x, p.position);
            let r = len3(r_ij);
            let w = cubic_kernel(r, p.h);
            let vj = p.volume();
            num += vj * p.conc() * w;
            denom += vj * w;
        }
        if denom.abs() > 1e-300 {
            num / denom
        } else {
            0.0
        }
    }
    /// Total scalar quantity: ∑ m_i * φ_i / ρ_i.
    pub fn total_scalar(&self, particles: &[TransportParticle]) -> f64 {
        particles.iter().map(|p| p.volume() * p.conc()).sum()
    }
}
/// SPH Laplacian diffusion operator.
///
/// dφ_i/dt = D · ∑_j (m_j/ρ_j) · (φ_j − φ_i) / r_ij² · (r_ij · ∇W_ij)
/// (Brookshaw / Morris form for the second derivative).
#[derive(Debug, Clone)]
pub struct DiffusionSph {
    /// Diffusion coefficient D \[m²/s\].
    pub diffusivity: f64,
}
impl DiffusionSph {
    /// Create a new diffusion operator.
    pub fn new(diffusivity: f64) -> Self {
        Self { diffusivity }
    }
    /// Compute concentration time derivatives for all particles.
    ///
    /// Returns a vector of dφ/dt values (one per particle).
    pub fn compute_dphidt(&self, particles: &[TransportParticle]) -> Vec<f64> {
        let n = particles.len();
        let mut dphidt = vec![0.0f64; n];
        for i in 0..n {
            let pi = &particles[i];
            for (j, pj) in particles.iter().enumerate().take(n) {
                if i == j {
                    continue;
                }
                let r_ij = sub3(pi.position, pj.position);
                let r = len3(r_ij);
                if r < 1e-12 {
                    continue;
                }
                let h_ij = 0.5 * (pi.h + pj.h);
                let grad_w = cubic_kernel_grad(r_ij, h_ij);
                let e_ij = scale3(r_ij, 1.0 / r);
                let dphi = pj.conc() - pi.conc();
                let vj = pj.volume();
                let term = 2.0 * dphi / r * dot3(e_ij, grad_w);
                dphidt[i] += self.diffusivity * vj * term;
            }
        }
        dphidt
    }
    /// Apply one explicit Euler diffusion step in-place.
    pub fn step(&self, particles: &mut [TransportParticle], dt: f64) {
        let dphidt = self.compute_dphidt(particles);
        for (p, &rate) in particles.iter_mut().zip(dphidt.iter()) {
            if let Some(c) = p.concentration.first_mut() {
                *c += rate * dt;
            }
        }
    }
}
/// Advection scheme for SPH scalar transport.
///
/// Supports pure Lagrangian advection (particles move with fluid) and a
/// semi-Lagrangian correction for upwind stabilisation.
#[derive(Debug, Clone)]
pub struct AdvectionSPH {
    /// CFL-limited time step \[s\].
    pub dt: f64,
    /// Artificial diffusion coefficient for upwind correction.
    pub artificial_diffusion: f64,
    /// Use semi-Lagrangian back-tracing instead of pure Lagrangian.
    pub semi_lagrangian: bool,
}
impl AdvectionSPH {
    /// Create a new [`AdvectionSPH`] scheme.
    pub fn new(dt: f64, artificial_diffusion: f64, semi_lagrangian: bool) -> Self {
        Self {
            dt,
            artificial_diffusion,
            semi_lagrangian,
        }
    }
    /// Advect a single scalar `c` at particle `i` using the material derivative.
    ///
    /// Dc/Dt = ∂c/∂t + u·∇c = 0 (in Lagrangian frame = 0).
    /// In practice, adds an upwind diffusion term if requested.
    pub fn advect_scalar(
        &self,
        particles: &[TransportParticle],
        i: usize,
        kernel_fn: impl Fn(f64, f64) -> f64,
        kernel_grad_fn: impl Fn([f64; 3], f64) -> [f64; 3],
    ) -> f64 {
        let pi = &particles[i];
        let ci = pi.concentration.first().copied().unwrap_or(0.0);
        let mut grad_c = [0.0_f64; 3];
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r = len3(r_ij);
            let grad_w = kernel_grad_fn(r_ij, pi.h);
            let cj = pj.concentration.first().copied().unwrap_or(0.0);
            let vol_j = pj.mass / pj.density.max(1e-300);
            grad_c[0] += vol_j * (cj - ci) * grad_w[0];
            grad_c[1] += vol_j * (cj - ci) * grad_w[1];
            grad_c[2] += vol_j * (cj - ci) * grad_w[2];
            let _ = r;
        }
        let speed = len3(pi.velocity);
        let upwind_diff = self.artificial_diffusion * speed * pi.h;
        let mut laplacian_c = 0.0_f64;
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r = len3(r_ij);
            if r < 1e-12 {
                continue;
            }
            let grad_w = kernel_grad_fn(r_ij, pi.h);
            let rdotgw = dot3(r_ij, grad_w);
            let cj = pj.concentration.first().copied().unwrap_or(0.0);
            laplacian_c +=
                pj.mass / pj.density.max(1e-300) * (cj - ci) * 2.0 * rdotgw / (r * r + 1e-20);
        }
        let conv = if self.semi_lagrangian {
            -(pi.velocity[0] * grad_c[0] + pi.velocity[1] * grad_c[1] + pi.velocity[2] * grad_c[2])
        } else {
            0.0
        };
        let _ = kernel_fn;
        conv + upwind_diff * laplacian_c
    }
    /// Update particle positions by Lagrangian advection: x ← x + v Δt.
    pub fn update_positions(&self, particles: &mut [TransportParticle]) {
        for p in particles.iter_mut() {
            p.position[0] += p.velocity[0] * self.dt;
            p.position[1] += p.velocity[1] * self.dt;
            p.position[2] += p.velocity[2] * self.dt;
        }
    }
    /// CFL number for a maximum velocity `v_max` and smoothing length `h`.
    pub fn cfl_number(v_max: f64, h: f64) -> f64 {
        if h < 1e-300 {
            return 0.0;
        }
        v_max * 1.0 / h
    }
}
/// A particle carrying scalar concentration, temperature, and velocity.
#[derive(Debug, Clone)]
pub struct TransportParticle {
    /// Position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Density \[kg/m³\].
    pub density: f64,
    /// Smoothing length \[m\].
    pub h: f64,
    /// Scalar concentration(s) — supports multi-species (index 0 = primary).
    pub concentration: Vec<f64>,
    /// Temperature \[K\].
    pub temperature: f64,
}
impl TransportParticle {
    /// Create a new transport particle with a single species.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        density: f64,
        h: f64,
        concentration: f64,
        temperature: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            mass,
            density,
            h,
            concentration: vec![concentration],
            temperature,
        }
    }
    /// Volume of this particle: V = m / ρ.
    pub fn volume(&self) -> f64 {
        if self.density.abs() < 1e-300 {
            0.0
        } else {
            self.mass / self.density
        }
    }
    /// Kinetic energy: ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }
    /// Primary species concentration (index 0).
    pub fn conc(&self) -> f64 {
        self.concentration.first().copied().unwrap_or(0.0)
    }
}
/// Advection scheme selector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AdvectionSchemeKind {
    /// Standard SPH (central, no upwinding).
    Standard,
    /// Upwind-biased: add numerical diffusion in the upwind direction.
    Upwind,
    /// WENO-inspired: use a flux limiter based on concentration gradient.
    WenoSph,
}
/// Multi-species transport particle for reactive-flow SPH.
///
/// Carries species concentrations and temperature alongside SPH kinematics.
#[derive(Clone, Debug)]
pub struct MultiTransportParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Density (kg/m³).
    pub rho: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Species volume/mass fractions.
    pub concentration: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
}
impl MultiTransportParticle {
    /// Create a new multi-species transport particle with `n_species` species.
    pub fn new(pos: [f64; 3], vel: [f64; 3], n_species: usize) -> Self {
        Self {
            pos,
            vel,
            rho: 1000.0,
            mass: 1.0,
            h: 0.1,
            concentration: vec![0.0; n_species],
            temperature: 293.15,
        }
    }
    /// Number of tracked species.
    pub fn species_count(&self) -> usize {
        self.concentration.len()
    }
    /// Sum of all species concentrations (should be ≤ 1 for volume fractions).
    pub fn total_mass_fraction(&self) -> f64 {
        self.concentration.iter().sum()
    }
}
/// Reactive SPH simulation driver.
pub struct ReactiveFlowSPH {
    /// Particle list.
    pub particles: Vec<MultiTransportParticle>,
    /// Number of tracked species.
    pub n_species: usize,
    /// Global reaction rate constant.
    pub reaction_rate: f64,
}
impl ReactiveFlowSPH {
    /// Create a new reactive flow SPH with `n_species` species.
    pub fn new(n_species: usize, reaction_rate: f64) -> Self {
        Self {
            particles: Vec::new(),
            n_species,
            reaction_rate,
        }
    }
    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: MultiTransportParticle) {
        self.particles.push(p);
    }
    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }
    /// Total mass of species `species` across all particles.
    pub fn total_species_mass(&self, species: usize) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let c = if species < p.concentration.len() {
                    p.concentration[species]
                } else {
                    0.0
                };
                c * p.mass
            })
            .sum()
    }
}
/// Phase change state for a particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhaseState {
    /// Fully solid.
    Solid,
    /// In the mushy (partially solid) zone.
    Mushy,
    /// Fully liquid.
    Liquid,
}
/// Reaction-diffusion system: ∂c/∂t = D∇²c + R(c).
///
/// Combines SPH diffusion with a pointwise reaction term.
#[derive(Debug, Clone)]
pub struct ReactionDiffusion {
    /// SPH diffusion operator.
    pub diffusion: DiffusionSph,
    /// Reaction rate coefficient k \[1/s\].
    pub reaction_k: f64,
}
impl ReactionDiffusion {
    /// Create a new reaction-diffusion system.
    pub fn new(diffusivity: f64, reaction_k: f64) -> Self {
        Self {
            diffusion: DiffusionSph::new(diffusivity),
            reaction_k,
        }
    }
    /// Compute total dφ/dt = D∇²φ + k·φ·(1 − φ) (Fisher-KPP default reaction).
    pub fn compute_dphidt_fisher(&self, particles: &[TransportParticle]) -> Vec<f64> {
        let mut dphidt = self.diffusion.compute_dphidt(particles);
        for (i, p) in particles.iter().enumerate() {
            let c = p.conc();
            dphidt[i] += self.reaction_k * c * (1.0 - c);
        }
        dphidt
    }
    /// One explicit Euler step with Fisher-KPP reaction.
    pub fn step_fisher(&self, particles: &mut [TransportParticle], dt: f64) {
        let dphidt = self.compute_dphidt_fisher(particles);
        for (p, &rate) in particles.iter_mut().zip(dphidt.iter()) {
            if let Some(c) = p.concentration.first_mut() {
                *c = (*c + rate * dt).clamp(0.0, 1.0);
            }
        }
    }
    /// One explicit Euler step with a custom reaction function.
    pub fn step_custom(&self, particles: &mut [TransportParticle], dt: f64, reaction: ReactionFn) {
        let mut dphidt = self.diffusion.compute_dphidt(particles);
        for (i, p) in particles.iter().enumerate() {
            dphidt[i] += reaction(p.conc());
        }
        for (p, &rate) in particles.iter_mut().zip(dphidt.iter()) {
            if let Some(c) = p.concentration.first_mut() {
                *c += rate * dt;
            }
        }
    }
}
/// SPH heat conduction: ρ cv ∂T/∂t = ∇·(k∇T).
///
/// Uses the Brookshaw/Morris SPH form for the Laplacian of temperature.
#[derive(Debug, Clone)]
pub struct ThermalConduction {
    /// Thermal conductivity k \[W/(m·K)\].
    pub conductivity: f64,
    /// Volumetric heat capacity ρ·cv \[J/(m³·K)\].
    pub rho_cv: f64,
}
impl ThermalConduction {
    /// Create a new thermal conduction operator.
    pub fn new(conductivity: f64, rho_cv: f64) -> Self {
        Self {
            conductivity,
            rho_cv,
        }
    }
    /// Thermal diffusivity α = k / (ρ·cv) \[m²/s\].
    pub fn diffusivity(&self) -> f64 {
        if self.rho_cv.abs() < 1e-300 {
            0.0
        } else {
            self.conductivity / self.rho_cv
        }
    }
    /// Compute temperature time derivatives dT_i/dt for all particles.
    pub fn compute_dtdt(&self, particles: &[TransportParticle]) -> Vec<f64> {
        let n = particles.len();
        let alpha = self.diffusivity();
        let mut dtdt = vec![0.0f64; n];
        for i in 0..n {
            let pi = &particles[i];
            for (j, pj) in particles.iter().enumerate().take(n) {
                if i == j {
                    continue;
                }
                let r_ij = sub3(pi.position, pj.position);
                let r = len3(r_ij);
                if r < 1e-12 {
                    continue;
                }
                let h_ij = 0.5 * (pi.h + pj.h);
                let grad_w = cubic_kernel_grad(r_ij, h_ij);
                let e_ij = scale3(r_ij, 1.0 / r);
                let dt_ij = pj.temperature - pi.temperature;
                let vj = pj.volume();
                let term = 2.0 * dt_ij / r * dot3(e_ij, grad_w);
                dtdt[i] += alpha * vj * term;
            }
        }
        dtdt
    }
    /// Apply one explicit Euler thermal conduction step in-place.
    pub fn step(&self, particles: &mut [TransportParticle], dt: f64) {
        let dtdt = self.compute_dtdt(particles);
        for (p, &rate) in particles.iter_mut().zip(dtdt.iter()) {
            p.temperature += rate * dt;
        }
    }
    /// Total thermal energy: ∑ (ρ cv) V_i T_i.
    pub fn total_thermal_energy(&self, particles: &[TransportParticle]) -> f64 {
        particles
            .iter()
            .map(|p| self.rho_cv * p.volume() * p.temperature)
            .sum()
    }
}
/// Fickian diffusion model J = −D ∇c.
///
/// Supports a scalar diffusivity or a full 3×3 anisotropic diffusivity tensor.
#[derive(Debug, Clone)]
pub struct FickianDiffusion {
    /// Scalar diffusivity D \[m²/s\] (isotropic case).
    pub diffusivity: f64,
    /// Optional 3×3 diffusivity tensor stored row-major \[D_xx, D_xy, D_xz, …\].
    pub tensor: Option<[f64; 9]>,
    /// Tortuosity factor τ (effective D_eff = D / τ).
    pub tortuosity: f64,
}
impl FickianDiffusion {
    /// Create an isotropic Fickian diffusion model.
    pub fn isotropic(diffusivity: f64) -> Self {
        Self {
            diffusivity,
            tensor: None,
            tortuosity: 1.0,
        }
    }
    /// Create an anisotropic model with a full diffusivity tensor.
    pub fn anisotropic(tensor: [f64; 9]) -> Self {
        let d = (tensor[0] + tensor[4] + tensor[8]) / 3.0;
        Self {
            diffusivity: d,
            tensor: Some(tensor),
            tortuosity: 1.0,
        }
    }
    /// Effective scalar diffusivity accounting for tortuosity.
    pub fn effective_diffusivity(&self) -> f64 {
        self.diffusivity / self.tortuosity.max(1e-300)
    }
    /// Compute the diffusive flux vector J = −D_eff ∇c at a point.
    ///
    /// `grad_c` is the concentration gradient \[mol/m⁴\].
    pub fn flux(&self, grad_c: [f64; 3]) -> [f64; 3] {
        let d = self.effective_diffusivity();
        match self.tensor {
            None => [-d * grad_c[0], -d * grad_c[1], -d * grad_c[2]],
            Some(t) => {
                let s = 1.0 / self.tortuosity.max(1e-300);
                [
                    -(t[0] * grad_c[0] + t[1] * grad_c[1] + t[2] * grad_c[2]) * s,
                    -(t[3] * grad_c[0] + t[4] * grad_c[1] + t[5] * grad_c[2]) * s,
                    -(t[6] * grad_c[0] + t[7] * grad_c[1] + t[8] * grad_c[2]) * s,
                ]
            }
        }
    }
    /// SPH Laplacian diffusion operator: Σ_j m_j/ρ_j (c_j − c_i) ∇²W_{ij}.
    ///
    /// Uses the second-order SPH Laplacian  2D Σ_j m_j(c_j−c_i)/(ρ_j r_{ij}²) (r·∇W).
    pub fn sph_laplacian(
        &self,
        particles: &[TransportParticle],
        i: usize,
        kernel_grad_fn: impl Fn([f64; 3], f64) -> [f64; 3],
    ) -> f64 {
        let pi = &particles[i];
        let d = self.effective_diffusivity();
        let mut laplacian = 0.0_f64;
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let r_ij = sub3(pi.position, pj.position);
            let r = len3(r_ij);
            if r < 1e-12 {
                continue;
            }
            let grad_w = kernel_grad_fn(r_ij, pi.h);
            let rdotgw = dot3(r_ij, grad_w);
            let ci = pi.concentration.first().copied().unwrap_or(0.0);
            let cj = pj.concentration.first().copied().unwrap_or(0.0);
            laplacian +=
                pj.mass / pj.density.max(1e-300) * (cj - ci) * 2.0 * rdotgw / (r * r + 1e-20);
        }
        d * laplacian
    }
}
/// SPH model for solidification and melting with a mushy zone.
///
/// Uses the enthalpy–porosity method to track the liquid fraction.
#[derive(Debug, Clone)]
pub struct PhaseChangeSPH {
    /// Solidus temperature T_s \[K\].
    pub t_solidus: f64,
    /// Liquidus temperature T_l \[K\].
    pub t_liquidus: f64,
    /// Latent heat of fusion L \[J/kg\].
    pub latent_heat: f64,
    /// Specific heat capacity cₚ \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Kozeny–Carman constant for mushy-zone resistance.
    pub kozeny_constant: f64,
}
impl PhaseChangeSPH {
    /// Create a new phase-change model.
    pub fn new(t_solidus: f64, t_liquidus: f64, latent_heat: f64, specific_heat: f64) -> Self {
        Self {
            t_solidus,
            t_liquidus,
            latent_heat,
            specific_heat,
            kozeny_constant: 1.6e5,
        }
    }
    /// Liquid fraction f_L(T) (0 = solid, 1 = liquid).
    pub fn liquid_fraction(&self, temperature: f64) -> f64 {
        if temperature <= self.t_solidus {
            return 0.0;
        }
        if temperature >= self.t_liquidus {
            return 1.0;
        }
        let range = self.t_liquidus - self.t_solidus;
        if range < 1e-300 {
            return 0.5;
        }
        (temperature - self.t_solidus) / range
    }
    /// Phase state of a particle at temperature T.
    pub fn phase_state(&self, temperature: f64) -> PhaseState {
        if temperature <= self.t_solidus {
            PhaseState::Solid
        } else if temperature >= self.t_liquidus {
            PhaseState::Liquid
        } else {
            PhaseState::Mushy
        }
    }
    /// Enthalpy h(T) = cₚ T + f_L · L \[J/kg\].
    pub fn enthalpy(&self, temperature: f64) -> f64 {
        self.specific_heat * temperature + self.liquid_fraction(temperature) * self.latent_heat
    }
    /// Temperature from enthalpy (inverse of enthalpy).
    ///
    /// Uses a bisection solve for the mushy zone.
    pub fn temperature_from_enthalpy(&self, h: f64) -> f64 {
        let h_s = self.specific_heat * self.t_solidus;
        let h_l = self.specific_heat * self.t_liquidus + self.latent_heat;
        if h <= h_s {
            return h / self.specific_heat.max(1e-300);
        }
        if h >= h_l {
            return self.t_liquidus + (h - h_l) / self.specific_heat.max(1e-300);
        }
        let mut lo = self.t_solidus;
        let mut hi = self.t_liquidus;
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.enthalpy(mid) < h {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
    /// Darcy mushy-zone resistance: A(f_L) = −C(1−f_L)² / (f_L³ + ε).
    ///
    /// Used in the momentum equation as a source term.
    pub fn mushy_zone_resistance(&self, liquid_fraction: f64) -> f64 {
        let fl = liquid_fraction.clamp(0.0, 1.0);
        let eps = 1e-6;
        -self.kozeny_constant * (1.0 - fl).powi(2) / (fl.powi(3) + eps)
    }
    /// Effective thermal capacity including latent heat:
    ///   c_eff = cₚ + L · df_L/dT.
    pub fn effective_heat_capacity(&self, temperature: f64) -> f64 {
        let range = self.t_liquidus - self.t_solidus;
        let dfl_dt =
            if temperature > self.t_solidus && temperature < self.t_liquidus && range > 1e-300 {
                1.0 / range
            } else {
                0.0
            };
        self.specific_heat + self.latent_heat * dfl_dt
    }
    /// Update enthalpy by Euler step: h ← h + Δt · DT/Dt · c_eff.
    pub fn enthalpy_euler_step(&self, temperature: f64, dt_dt: f64, dt: f64) -> f64 {
        let c_eff = self.effective_heat_capacity(temperature);
        let h0 = self.enthalpy(temperature);
        h0 + c_eff * dt_dt * dt
    }
}
/// SPH advection-diffusion-reaction transport.
///
/// Models `∂c/∂t = D ∇²c - k c + S` where k is a first-order decay rate
/// and S is an external source term.
#[derive(Debug, Clone)]
pub struct SPHReactiveTransport {
    /// Diffusion coefficient D \[m²/s\].
    pub diffusivity: f64,
    /// First-order reaction rate k \[1/s\] (positive = decay).
    pub reaction_rate: f64,
    /// Smoothing length h \[m\].
    pub smoothing_length: f64,
}
impl SPHReactiveTransport {
    /// Create a new `SPHReactiveTransport`.
    pub fn new(diffusivity: f64, reaction_rate: f64, smoothing_length: f64) -> Self {
        Self {
            diffusivity,
            reaction_rate,
            smoothing_length,
        }
    }
    /// Rate of change `dc/dt = D ∇²c - k c + source`.
    ///
    /// # Arguments
    /// * `c_i`        – Concentration at particle i.
    /// * `neighbours` – Slice of `(c_j, r_ij, mass_j, rho_j)`.
    /// * `source`     – External source S \[mol/(m³·s)\].
    pub fn dcdt(&self, c_i: f64, neighbours: &[(f64, [f64; 3], f64, f64)], source: f64) -> f64 {
        let h = self.smoothing_length;
        let mut lap = 0.0_f64;
        for &(c_j, r_ij, mass_j, rho_j) in neighbours {
            let r = len3(r_ij);
            let w = cubic_kernel(r, h);
            if rho_j > 0.0_f64 {
                lap += 2.0_f64 * (c_j - c_i) * mass_j / rho_j * w;
            }
        }
        self.diffusivity * lap - self.reaction_rate * c_i + source
    }
    /// Damköhler number Da = k L² / D (ratio of reaction to diffusion).
    ///
    /// Returns infinity when D = 0.
    pub fn damkohler(&self, length_scale: f64) -> f64 {
        if self.diffusivity <= 0.0_f64 {
            return f64::INFINITY;
        }
        self.reaction_rate * length_scale * length_scale / self.diffusivity
    }
}
/// Gray-Scott reaction-diffusion system via SPH.
///
/// ∂u/∂t = Du ∇²u − u·v² + f·(1 − u)
/// ∂v/∂t = Dv ∇²v + u·v² − (f + k)·v
///
/// Particles carry two species: `concentration[0]` = u, `concentration[1]` = v.
#[derive(Debug, Clone)]
pub struct GrayScottSph {
    /// Diffusion coefficient for species u.
    pub du: f64,
    /// Diffusion coefficient for species v.
    pub dv: f64,
    /// Feed rate f.
    pub feed: f64,
    /// Kill rate k.
    pub kill: f64,
}
impl GrayScottSph {
    /// Create a new Gray-Scott model.
    pub fn new(du: f64, dv: f64, feed: f64, kill: f64) -> Self {
        Self { du, dv, feed, kill }
    }
    /// Ensure all particles have two concentration slots.
    pub fn init_particles(particles: &mut [TransportParticle]) {
        for p in particles.iter_mut() {
            while p.concentration.len() < 2 {
                p.concentration.push(0.0);
            }
        }
    }
    /// Compute (du/dt, dv/dt) for all particles.
    pub fn compute_rates(&self, particles: &[TransportParticle]) -> (Vec<f64>, Vec<f64>) {
        let n = particles.len();
        let diff_u = DiffusionSph::new(self.du);
        let diff_v = DiffusionSph::new(self.dv);
        let mut particles_u: Vec<TransportParticle> = particles.to_vec();
        let mut particles_v: Vec<TransportParticle> = particles.to_vec();
        for i in 0..n {
            if particles_u[i].concentration.len() >= 2 {
                let v_val = particles_u[i].concentration[1];
                particles_u[i].concentration[0] = particles[i].concentration[0];
                particles_v[i].concentration[0] = v_val;
            }
        }
        let lap_u = diff_u.compute_dphidt(&particles_u);
        let lap_v = diff_v.compute_dphidt(&particles_v);
        let mut du_dt = vec![0.0f64; n];
        let mut dv_dt = vec![0.0f64; n];
        for i in 0..n {
            let u = if !particles[i].concentration.is_empty() {
                particles[i].concentration[0]
            } else {
                0.0
            };
            let v = if particles[i].concentration.len() >= 2 {
                particles[i].concentration[1]
            } else {
                0.0
            };
            let uvv = u * v * v;
            du_dt[i] = lap_u[i] - uvv + self.feed * (1.0 - u);
            dv_dt[i] = lap_v[i] + uvv - (self.feed + self.kill) * v;
        }
        (du_dt, dv_dt)
    }
    /// One explicit Euler step.
    pub fn step(&self, particles: &mut [TransportParticle], dt: f64) {
        let (du_dt, dv_dt) = self.compute_rates(particles);
        for (i, p) in particles.iter_mut().enumerate() {
            while p.concentration.len() < 2 {
                p.concentration.push(0.0);
            }
            p.concentration[0] = (p.concentration[0] + du_dt[i] * dt).clamp(0.0, 1.0);
            p.concentration[1] = (p.concentration[1] + dv_dt[i] * dt).clamp(0.0, 1.0);
        }
    }
}
