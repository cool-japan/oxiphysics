//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

/// Elastic SPH solver managing a collection of [`ElasticParticle`]s.
pub struct ElasticSolver {
    /// All elastic particles managed by this solver.
    pub particles: Vec<ElasticParticle>,
    /// Young's modulus (Pa) applied to newly added particles.
    pub e_mod: f64,
    /// Poisson's ratio applied to newly added particles.
    pub nu: f64,
}
impl ElasticSolver {
    /// Create a new empty [`ElasticSolver`].
    ///
    /// # Arguments
    /// * `e_mod` – Young's modulus (Pa)
    /// * `nu`    – Poisson's ratio
    pub fn new(e_mod: f64, nu: f64) -> Self {
        Self {
            particles: Vec::new(),
            e_mod,
            nu,
        }
    }
    /// Add a particle to the solver.
    pub fn add_particle(&mut self, p: ElasticParticle) {
        self.particles.push(p);
    }
    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }
    /// Total strain energy (J) = Σᵢ (mᵢ / ρᵢ) · (½ σᵢ : εᵢ).
    ///
    /// Uses the infinitesimal strain assumption.
    pub fn strain_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let vol = p.mass / p.density.max(1e-30);
                let strain = green_lagrange_strain(&p.deformation_gradient);
                let sigma = linear_elastic_stress(&strain, p.youngs_modulus, p.poissons_ratio);
                let double_contract: f64 = (0..3)
                    .flat_map(|i| (0..3).map(move |j| sigma[i][j] * strain[i][j]))
                    .sum();
                vol * 0.5 * double_contract
            })
            .sum()
    }
}
/// Extended SPH particle for solid/elastic mechanics with stress tensor,
/// strain rate, and plasticity state.
#[allow(dead_code)]
pub struct ElasticSphParticle {
    /// Current position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Acceleration (m/s²).
    pub acceleration: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Cauchy stress tensor σ (3×3, Pa), row-major.
    pub stress: [[f64; 3]; 3],
    /// Deviatoric stress deviator s (3×3, Pa).
    pub stress_dev: [[f64; 3]; 3],
    /// Strain rate tensor ε̇ (3×3, 1/s), symmetric.
    pub strain_rate: [[f64; 3]; 3],
    /// Equivalent plastic strain (dimensionless, accumulated).
    pub plastic_strain: f64,
    /// Plastic strain rate (1/s).
    pub plastic_strain_rate: f64,
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν.
    pub poissons_ratio: f64,
    /// Fracture flag.
    pub fractured: bool,
}
impl ElasticSphParticle {
    /// Create a new [`ElasticSphParticle`] with zero stress and strain rate.
    ///
    /// # Arguments
    /// * `position`       – initial position (m)
    /// * `mass`           – particle mass (kg)
    /// * `density`        – initial density (kg/m³)
    /// * `h`              – smoothing length (m)
    /// * `youngs_modulus` – E (Pa)
    /// * `poissons_ratio` – ν (dimensionless)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: [f64; 3],
        mass: f64,
        density: f64,
        h: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
    ) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            mass,
            density,
            h,
            stress: [[0.0; 3]; 3],
            stress_dev: [[0.0; 3]; 3],
            strain_rate: [[0.0; 3]; 3],
            plastic_strain: 0.0,
            plastic_strain_rate: 0.0,
            youngs_modulus,
            poissons_ratio,
            fractured: false,
        }
    }
    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poissons_ratio))
    }
    /// Bulk modulus K = E / (3(1−2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poissons_ratio))
    }
    /// Lamé first parameter λ.
    pub fn lame_lambda(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poissons_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }
    /// Hydrostatic pressure from trace of stress: p = −tr(σ)/3.
    pub fn pressure(&self) -> f64 {
        -(self.stress[0][0] + self.stress[1][1] + self.stress[2][2]) / 3.0
    }
    /// Equivalent (von Mises) stress: σ_eq = sqrt(3/2 s:s).
    pub fn von_mises_eq(&self) -> f64 {
        let s = &self.stress_dev;
        let s_sq: f64 = s.iter().flatten().map(|v| v * v).sum();
        (1.5 * s_sq).sqrt()
    }
}
/// Johnson-Cook constitutive model for high-strain-rate plasticity.
///
/// Yield stress: σ_y = (A + B ε_p^n)(1 + C ln(ε̇*/ε̇₀))(1 − T*^m)
///
/// where:
/// - A = initial yield stress (Pa)
/// - B = strain hardening coefficient (Pa)
/// - n = strain hardening exponent
/// - C = strain rate sensitivity coefficient
/// - m = thermal softening exponent
/// - ε_p = equivalent plastic strain
/// - ε̇* = ε̇/ε̇₀ = normalised strain rate
/// - T* = (T − T_ref)/(T_melt − T_ref) = homologous temperature
#[allow(dead_code)]
#[derive(Clone)]
pub struct JohnsonCookModel {
    /// Initial yield stress A (Pa).
    pub a: f64,
    /// Strain hardening coefficient B (Pa).
    pub b: f64,
    /// Strain hardening exponent n.
    pub n: f64,
    /// Strain rate sensitivity C.
    pub c: f64,
    /// Thermal softening exponent m.
    pub m: f64,
    /// Reference strain rate ε̇₀ (1/s).
    pub eps_dot_0: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Melt temperature T_melt (K).
    pub t_melt: f64,
}
impl JohnsonCookModel {
    /// Create a new Johnson-Cook model for steel (default parameters).
    pub fn new_steel() -> Self {
        Self {
            a: 792e6,
            b: 510e6,
            n: 0.26,
            c: 0.014,
            m: 1.03,
            eps_dot_0: 1.0,
            t_ref: 293.0,
            t_melt: 1793.0,
        }
    }
    /// Create a Johnson-Cook model with explicit parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a: f64,
        b: f64,
        n: f64,
        c: f64,
        m: f64,
        eps_dot_0: f64,
        t_ref: f64,
        t_melt: f64,
    ) -> Self {
        Self {
            a,
            b,
            n,
            c,
            m,
            eps_dot_0,
            t_ref,
            t_melt,
        }
    }
    /// Evaluate Johnson-Cook yield stress.
    ///
    /// # Arguments
    /// * `eps_p`     – equivalent plastic strain
    /// * `eps_dot`   – equivalent plastic strain rate (1/s)
    /// * `temp`      – current temperature (K)
    pub fn yield_stress(&self, eps_p: f64, eps_dot: f64, temp: f64) -> f64 {
        let strain_hardening = self.a + self.b * eps_p.powf(self.n);
        let rate_factor = if self.c.abs() < 1e-30 {
            1.0
        } else {
            let eps_star = (eps_dot / self.eps_dot_0).max(1.0);
            1.0 + self.c * eps_star.ln()
        };
        let t_star = if (self.t_melt - self.t_ref).abs() < 1.0 {
            0.0
        } else {
            ((temp - self.t_ref) / (self.t_melt - self.t_ref)).clamp(0.0, 1.0)
        };
        let thermal_factor = 1.0 - t_star.powf(self.m);
        strain_hardening * rate_factor * thermal_factor
    }
    /// Radial return mapping: project stress deviator to yield surface.
    ///
    /// Returns updated deviatoric stress and incremental plastic strain.
    /// Uses the associative flow rule with von Mises yield criterion.
    pub fn radial_return(
        &self,
        s_trial: [[f64; 3]; 3],
        eps_p: f64,
        eps_dot: f64,
        temp: f64,
        shear_mod: f64,
    ) -> ([[f64; 3]; 3], f64) {
        let s_norm_sq: f64 = s_trial.iter().flatten().map(|v| v * v).sum();
        let s_norm = s_norm_sq.sqrt();
        let yield_stress = self.yield_stress(eps_p, eps_dot, temp);
        let yield_limit = (2.0 / 3.0_f64).sqrt() * yield_stress;
        if s_norm <= yield_limit + 1e-30 {
            return (s_trial, 0.0);
        }
        let scale = yield_limit / s_norm;
        let mut s_ret = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                s_ret[i][j] = scale * s_trial[i][j];
            }
        }
        let d_eps_p = (s_norm - yield_limit) / (2.0 * shear_mod * (2.0 / 3.0_f64).sqrt());
        (s_ret, d_eps_p.max(0.0))
    }
}
/// A single elastic SPH particle.
pub struct ElasticParticle {
    /// Current position (m).
    pub position: [f64; 3],
    /// Reference (initial) position (m).
    pub ref_position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Acceleration (m/s²) – updated each step.
    pub acceleration: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Deformation gradient F (3×3), row-major.
    pub deformation_gradient: [[f64; 3]; 3],
    /// Cauchy stress tensor σ (3×3), row-major (Pa).
    pub stress: [[f64; 3]; 3],
    /// First Piola-Kirchhoff stress P (3×3), row-major (Pa).
    pub pk1_stress: [[f64; 3]; 3],
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν.
    pub poissons_ratio: f64,
    /// Fracture flag: true if the particle has failed.
    pub fractured: bool,
}
impl ElasticParticle {
    /// Create a new elastic particle.  `deformation_gradient` is initialised to
    /// the identity and `stress` to zero.
    ///
    /// # Arguments
    /// * `position`       – initial (reference) position (m)
    /// * `mass`           – particle mass (kg)
    /// * `density`        – initial density (kg/m³)
    /// * `h`              – smoothing length (m)
    /// * `youngs_modulus` – E (Pa)
    /// * `poissons_ratio` – ν (dimensionless)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        position: [f64; 3],
        mass: f64,
        density: f64,
        h: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
    ) -> Self {
        Self {
            position,
            ref_position: position,
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            mass,
            density,
            h,
            deformation_gradient: IDENTITY3,
            stress: [[0.0; 3]; 3],
            pk1_stress: [[0.0; 3]; 3],
            youngs_modulus,
            poissons_ratio,
            fractured: false,
        }
    }
    /// Lamé parameter μ (shear modulus).
    pub fn shear_modulus(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poissons_ratio;
        e / (2.0 * (1.0 + nu))
    }
    /// Lamé parameter λ.
    pub fn lame_lambda(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poissons_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }
}
/// Symmetric velocity gradient (strain rate) tensor ε̇ = ½(L + Lᵀ)
/// where L = ∂vᵢ/∂xⱼ is the velocity gradient.
///
/// Computed via SPH summation over neighbours.
#[allow(dead_code)]
pub struct StrainRateTensor {
    /// The 3×3 symmetric strain-rate tensor (1/s), row-major.
    pub tensor: [[f64; 3]; 3],
    /// Volumetric strain rate: tr(ε̇) = ε̇₁₁ + ε̇₂₂ + ε̇₃₃.
    pub volumetric: f64,
    /// Second invariant J₂ = ½ dev(ε̇):dev(ε̇).
    pub j2: f64,
}
impl StrainRateTensor {
    /// Construct from a raw 3×3 tensor, computing derived quantities.
    pub fn from_tensor(raw: [[f64; 3]; 3]) -> Self {
        let mut sym = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sym[i][j] = 0.5 * (raw[i][j] + raw[j][i]);
            }
        }
        let vol = sym[0][0] + sym[1][1] + sym[2][2];
        let vol3 = vol / 3.0;
        let mut dev = sym;
        for i in 0..3 {
            dev[i][i] -= vol3;
        }
        let j2 = 0.5 * dev.iter().flatten().map(|v| v * v).sum::<f64>();
        Self {
            tensor: sym,
            volumetric: vol,
            j2,
        }
    }
    /// Compute strain-rate tensor for particle `i` via SPH summation.
    ///
    /// ε̇ᵢⱼ = ½ Σₖ mₖ/ρₖ (vₖ − vᵢ)ᵢ ∂Wᵢₖ/∂xⱼ + (symmetric)
    pub fn compute_sph(particles: &[ElasticSphParticle], idx: usize) -> Self {
        if idx >= particles.len() {
            return Self::from_tensor([[0.0; 3]; 3]);
        }
        let pi = &particles[idx];
        let mut l = [[0.0f64; 3]; 3];
        for (k, pk) in particles.iter().enumerate() {
            if k == idx {
                continue;
            }
            let rij: [f64; 3] = [
                pi.position[0] - pk.position[0],
                pi.position[1] - pk.position[1],
                pi.position[2] - pk.position[2],
            ];
            let r = (rij[0].powi(2) + rij[1].powi(2) + rij[2].powi(2)).sqrt();
            if r < 1e-300 {
                continue;
            }
            let h_avg = 0.5 * (pi.h + pk.h);
            let dw_dr = {
                let q = r / h_avg;
                let alpha = 1.0 / (std::f64::consts::PI * h_avg.powi(3));
                if q < 1.0 {
                    alpha / h_avg * (-3.0 * q + 2.25 * q * q)
                } else if q < 2.0 {
                    -alpha / h_avg * 0.75 * (2.0 - q).powi(2)
                } else {
                    0.0
                }
            };
            let grad_w: [f64; 3] = [dw_dr * rij[0] / r, dw_dr * rij[1] / r, dw_dr * rij[2] / r];
            let vol_k = pk.mass / pk.density.max(1e-30);
            let dv: [f64; 3] = [
                pk.velocity[0] - pi.velocity[0],
                pk.velocity[1] - pi.velocity[1],
                pk.velocity[2] - pi.velocity[2],
            ];
            for a in 0..3 {
                for b in 0..3 {
                    l[a][b] += vol_k * dv[a] * grad_w[b];
                }
            }
        }
        let mut sym = [[0.0f64; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                sym[a][b] = 0.5 * (l[a][b] + l[b][a]);
            }
        }
        let vol = sym[0][0] + sym[1][1] + sym[2][2];
        let vol3 = vol / 3.0;
        let mut dev = sym;
        for a in 0..3 {
            dev[a][a] -= vol3;
        }
        let j2 = 0.5 * dev.iter().flatten().map(|v| v * v).sum::<f64>();
        Self {
            tensor: sym,
            volumetric: vol,
            j2,
        }
    }
    /// Equivalent strain rate: ε̇_eq = sqrt(2/3 · ε̇:ε̇).
    pub fn equivalent(&self) -> f64 {
        let sq: f64 = self.tensor.iter().flatten().map(|v| v * v).sum();
        (2.0 / 3.0 * sq).sqrt()
    }
}
/// Wave speed utilities for elastic media.
#[allow(dead_code)]
pub struct ElasticWaveSpeed;
impl ElasticWaveSpeed {
    /// Longitudinal (P-wave) speed: c_P = sqrt((λ + 2G)/ρ).
    pub fn p_wave(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let nu = poissons_ratio;
        let e = youngs_modulus;
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let g = e / (2.0 * (1.0 + nu));
        ((lambda + 2.0 * g) / density).sqrt()
    }
    /// Shear (S-wave) speed: c_S = sqrt(G/ρ).
    pub fn s_wave(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let g = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
        (g / density).sqrt()
    }
    /// Both wave speeds: returns (c_P, c_S).
    pub fn both(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> (f64, f64) {
        (
            Self::p_wave(youngs_modulus, poissons_ratio, density),
            Self::s_wave(youngs_modulus, poissons_ratio, density),
        )
    }
    /// SPH Courant time step: Δt_CFL = CFL * h / c_P.
    pub fn courant_dt(
        h: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
        density: f64,
        cfl: f64,
    ) -> f64 {
        let cp = Self::p_wave(youngs_modulus, poissons_ratio, density);
        if cp < 1e-30 {
            return f64::INFINITY;
        }
        cfl * h / cp
    }
    /// Rayleigh wave speed approximation: c_R ≈ c_S * (0.862 + 1.14ν)/(1 + ν).
    pub fn rayleigh_wave(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
        let cs = Self::s_wave(youngs_modulus, poissons_ratio, density);
        let nu = poissons_ratio;
        cs * (0.862 + 1.14 * nu) / (1.0 + nu)
    }
    /// Bulk wave speed: c_B = sqrt(K/ρ) = sqrt(λ + 2G/3)/ρ.
    pub fn bulk_wave(youngs_modulus: f64, poissons_ratio: f64, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let k = youngs_modulus / (3.0 * (1.0 - 2.0 * poissons_ratio));
        (k / density).sqrt()
    }
}
/// Grady-Kipp fragmentation model for SPH particles.
///
/// Implements the scalar damage variable D ∈ \[0, 1\] where D = 0 is intact
/// and D = 1 is fully fractured.  Based on Grady-Kipp (1980).
#[allow(dead_code)]
#[derive(Clone)]
pub struct SphFracture {
    /// Weibull modulus m (typically 5–10 for rocks, higher for metals).
    pub weibull_m: f64,
    /// Weibull scale coefficient k (material-dependent, Pa^-m · s^-1).
    pub weibull_k: f64,
    /// Critical stress threshold below which no damage accumulates (Pa).
    pub sigma_threshold: f64,
    /// Fracture energy G_f (J/m²).
    pub fracture_energy: f64,
}
impl SphFracture {
    /// Create a Grady-Kipp model for granite (representative parameters).
    pub fn new_granite() -> Self {
        Self {
            weibull_m: 6.0,
            weibull_k: 1e24,
            sigma_threshold: 50e6,
            fracture_energy: 100.0,
        }
    }
    /// Create with explicit parameters.
    pub fn new(weibull_m: f64, weibull_k: f64, sigma_threshold: f64, fracture_energy: f64) -> Self {
        Self {
            weibull_m,
            weibull_k,
            sigma_threshold,
            fracture_energy,
        }
    }
    /// Compute damage rate Ḋ from current tensile pressure `p_tens` (Pa).
    ///
    /// Ḋ = k * p_tens^m  when p_tens > σ_threshold,
    /// Ḋ = 0             otherwise.
    pub fn damage_rate(&self, p_tens: f64) -> f64 {
        if p_tens <= self.sigma_threshold {
            return 0.0;
        }
        self.weibull_k * p_tens.powf(self.weibull_m)
    }
    /// Advance damage variable by one time step.
    ///
    /// D(t+dt) = min(D + Ḋ * dt, 1).
    pub fn advance_damage(&self, damage: f64, p_tens: f64, dt: f64) -> f64 {
        (damage + self.damage_rate(p_tens) * dt).min(1.0)
    }
    /// Apply damage to stress tensor: σ_eff = (1 − D) σ.
    ///
    /// Fully damaged (D=1) material cannot carry tensile stress.
    pub fn effective_stress(&self, stress: [[f64; 3]; 3], damage: f64) -> [[f64; 3]; 3] {
        let scale = (1.0 - damage).max(0.0);
        let mut eff = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                eff[i][j] = scale * stress[i][j];
            }
        }
        eff
    }
    /// Fragment size estimate (Grady): s ≈ sqrt(20 K_Ic² / (ρ ε̇²))
    ///
    /// # Arguments
    /// * `k_ic`       – fracture toughness (Pa·m^0.5)
    /// * `density`    – material density (kg/m³)
    /// * `strain_rate`– strain rate at fracture (1/s)
    pub fn fragment_size(&self, k_ic: f64, density: f64, strain_rate: f64) -> f64 {
        if density <= 0.0 || strain_rate <= 0.0 {
            return f64::INFINITY;
        }
        (20.0 * k_ic * k_ic / (density * strain_rate * strain_rate)).sqrt()
    }
    /// Grady fragment velocity: v_f = (4/3 * G_f * ε̇ / ρ)^{1/3}.
    pub fn fragment_velocity(&self, density: f64, strain_rate: f64) -> f64 {
        if density <= 0.0 || strain_rate <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0 * self.fracture_energy * strain_rate / density).powf(1.0 / 3.0)
    }
}
/// Elastic stress update for SPH particles using Hooke's law in a
/// co-rotational (Jaumann) formulation to handle finite rotations.
#[allow(dead_code)]
pub struct SphElasticStress {
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
}
impl SphElasticStress {
    /// Create a new [`SphElasticStress`] from Young's modulus and Poisson's ratio.
    pub fn from_elastic_params(youngs_modulus: f64, poissons_ratio: f64) -> Self {
        let g = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
        let k = youngs_modulus / (3.0 * (1.0 - 2.0 * poissons_ratio));
        Self {
            shear_modulus: g,
            bulk_modulus: k,
        }
    }
    /// Jaumann stress-rate: Ṡᵢⱼ = 2G ε̇ᵢⱼᵈᵉᵛ + K ε̇ₖₖ δᵢⱼ + σᵢₖ Ωₖⱼ − Ωᵢₖ σₖⱼ
    ///
    /// where Ω = ½(L − Lᵀ) is the spin tensor.
    /// Returns the incremental stress dσ = Ṡ * dt.
    pub fn jaumann_stress_increment(
        &self,
        stress: [[f64; 3]; 3],
        strain_rate: [[f64; 3]; 3],
        spin: [[f64; 3]; 3],
        dt: f64,
    ) -> [[f64; 3]; 3] {
        let vol_rate = strain_rate[0][0] + strain_rate[1][1] + strain_rate[2][2];
        let mut ds = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let eps_dev = strain_rate[i][j] - if i == j { vol_rate / 3.0 } else { 0.0 };
                ds[i][j] += 2.0 * self.shear_modulus * eps_dev;
                if i == j {
                    ds[i][j] += self.bulk_modulus * vol_rate;
                }
                for k in 0..3 {
                    ds[i][j] += stress[i][k] * spin[k][j] - spin[i][k] * stress[k][j];
                }
                ds[i][j] *= dt;
            }
        }
        ds
    }
    /// Update particle stress using Jaumann formulation for one time step.
    ///
    /// Modifies `particle.stress` and recomputes `particle.stress_dev`.
    pub fn update_particle_stress(
        &self,
        particle: &mut ElasticSphParticle,
        velocity_gradient: [[f64; 3]; 3],
        dt: f64,
    ) {
        let mut eps_dot = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                eps_dot[i][j] = 0.5 * (velocity_gradient[i][j] + velocity_gradient[j][i]);
            }
        }
        let mut omega = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                omega[i][j] = 0.5 * (velocity_gradient[i][j] - velocity_gradient[j][i]);
            }
        }
        let ds = self.jaumann_stress_increment(particle.stress, eps_dot, omega, dt);
        for i in 0..3 {
            for j in 0..3 {
                particle.stress[i][j] += ds[i][j];
            }
        }
        let tr = particle.stress[0][0] + particle.stress[1][1] + particle.stress[2][2];
        let p = tr / 3.0;
        for i in 0..3 {
            for j in 0..3 {
                particle.stress_dev[i][j] = particle.stress[i][j];
            }
            particle.stress_dev[i][i] -= p;
        }
    }
}
/// Choice of elastic constitutive model.
pub enum ConstitutiveModel {
    /// Linear St. Venant-Kirchhoff model (suitable for small strains).
    StVenantKirchhoff,
    /// Compressible Neo-Hookean model (valid for large deformations).
    NeoHookean,
}
