//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    EX9, EY9, W9, drag_coeff_particle, feq_d2q9, mag3, mixture_viscosity, richardson_zaki_exponent,
    stokes_settling,
};

/// Bed formation tracker.
///
/// Tracks deposited particles, bed height evolution, and bed shear stress.
#[derive(Debug, Clone)]
pub struct BedFormation {
    /// Grid width.
    pub nx: usize,
    /// Bed height profile h\[x\] in meters.
    pub bed_height: Vec<f64>,
    /// Deposited mass per cell in kg.
    pub deposited_mass: Vec<f64>,
    /// Particle diameter in meters.
    pub particle_diameter: f64,
    /// Fluid density kg/m³.
    pub rho_fluid: f64,
    /// Particle density kg/m³.
    pub rho_solid: f64,
}
impl BedFormation {
    /// Create a new bed formation tracker.
    pub fn new(nx: usize, particle_diameter: f64, rho_fluid: f64, rho_solid: f64) -> Self {
        Self {
            nx,
            bed_height: vec![0.0; nx],
            deposited_mass: vec![0.0; nx],
            particle_diameter,
            rho_fluid,
            rho_solid,
        }
    }
    /// Deposit a particle at grid column x.
    pub fn deposit(&mut self, x: usize, mass: f64) {
        if x < self.nx {
            let vol = mass / self.rho_solid;
            let packing = 0.64;
            let dh = vol / (self.particle_diameter * packing);
            self.bed_height[x] += dh;
            self.deposited_mass[x] += mass;
        }
    }
    /// Compute bed shear stress at column x given fluid velocity u.
    pub fn shear_stress(&self, x: usize, u_fluid: f64, mu: f64) -> f64 {
        if x >= self.nx {
            return 0.0;
        }
        let h = self.bed_height[x];
        if h < 1e-10 {
            return 0.0;
        }
        mu * u_fluid / (h / 2.0 + 1e-10)
    }
    /// Shields parameter (non-dimensional shear stress).
    pub fn shields_parameter(&self, x: usize, u_fluid: f64, mu: f64) -> f64 {
        let tau = self.shear_stress(x, u_fluid, mu);
        let d = self.particle_diameter;
        tau / ((self.rho_solid - self.rho_fluid) * 9.81 * d)
    }
    /// Average bed height.
    pub fn avg_height(&self) -> f64 {
        self.bed_height.iter().sum::<f64>() / self.nx as f64
    }
    /// Maximum bed height.
    pub fn max_height(&self) -> f64 {
        self.bed_height
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}
/// Darcy-Brinkman-Forchheimer porous flow in LBM.
///
/// Implements permeability, porosity, and drag terms.
#[derive(Debug, Clone)]
pub struct DarcyFlowLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Porosity field epsilon\[y\]\[x\].
    pub epsilon: Vec<Vec<f64>>,
    /// Permeability field K\[y\]\[x\] in m².
    pub permeability: Vec<Vec<f64>>,
    /// Distribution functions.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Relaxation parameter.
    pub omega: f64,
    /// Fluid density.
    pub rho: Vec<Vec<f64>>,
    /// Velocity field.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Forchheimer coefficient.
    pub beta: f64,
}
impl DarcyFlowLbm {
    /// Create a new Darcy flow LBM.
    pub fn new(nx: usize, ny: usize, omega: f64, epsilon_val: f64, k_perm: f64, beta: f64) -> Self {
        let epsilon = vec![vec![epsilon_val; nx]; ny];
        let permeability = vec![vec![k_perm; nx]; ny];
        let f = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let rho = vec![vec![1.0f64; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let mut sim = Self {
            nx,
            ny,
            epsilon,
            permeability,
            f,
            omega,
            rho,
            u,
            beta,
        };
        sim.initialize();
        sim
    }
    /// Initialize equilibrium distributions.
    pub fn initialize(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let feq = feq_d2q9(1.0, 0.0, 0.0);
                self.f[y][x] = feq.to_vec();
            }
        }
    }
    /// Compute Darcy drag force term.
    pub fn darcy_force(&self, x: usize, y: usize) -> [f64; 2] {
        let k_perm = self.permeability[y][x];
        let eps = self.epsilon[y][x];
        let u = self.u[y][x];
        let umag = (u[0] * u[0] + u[1] * u[1]).sqrt();
        let nu = (1.0 / self.omega - 0.5) / 3.0;
        let coeff = nu / k_perm + self.beta / eps * umag;
        [-coeff * u[0], -coeff * u[1]]
    }
    /// One BGK collision step with Darcy body force.
    pub fn collide(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.rho[y][x];
                let [ux, uy] = self.u[y][x];
                let feq = feq_d2q9(rho, ux, uy);
                let fd = self.darcy_force(x, y);
                for q in 0..9 {
                    let fq = W9[q] * 3.0 * (EX9[q] as f64 * fd[0] + EY9[q] as f64 * fd[1]);
                    self.f[y][x][q] =
                        self.f[y][x][q] * (1.0 - self.omega) + self.omega * feq[q] + fq;
                }
            }
        }
    }
    /// One streaming step.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let xn = (x as i32 + EX9[q]).rem_euclid(nx as i32) as usize;
                    let yn = (y as i32 + EY9[q]).rem_euclid(ny as i32) as usize;
                    f_new[yn][xn][q] = self.f[y][x][q];
                }
            }
        }
        self.f = f_new;
    }
    /// Update macroscopic variables.
    pub fn update_macros(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for q in 0..9 {
                    rho += self.f[y][x][q];
                    ux += EX9[q] as f64 * self.f[y][x][q];
                    uy += EY9[q] as f64 * self.f[y][x][q];
                }
                self.rho[y][x] = rho;
                if rho > 1e-10 {
                    self.u[y][x] = [ux / rho, uy / rho];
                }
            }
        }
    }
    /// Full time step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }
    /// Average Darcy velocity magnitude.
    pub fn avg_darcy_velocity(&self) -> f64 {
        let mut sum = 0.0;
        let n = self.nx * self.ny;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let u = self.u[y][x];
                sum += (u[0] * u[0] + u[1] * u[1]).sqrt() * self.epsilon[y][x];
            }
        }
        sum / n as f64
    }
}
/// Floc aggregation–breakup dynamics model.
///
/// Uses a simplified population balance approach with aggregation kernel
/// and turbulent breakup.
#[derive(Debug, Clone)]
pub struct FlocDynamics {
    /// Fluid dynamic viscosity (Pa·s).
    pub mu: f64,
    /// Aggregation efficiency α ∈ \[0, 1\].
    pub aggregation_efficiency: f64,
    /// Breakup rate constant B_0 (1/s).
    pub breakup_rate: f64,
    /// Turbulent energy dissipation rate ε (m²/s³).
    pub epsilon_turb: f64,
    /// Turbulent breakup exponent p.
    pub breakup_exponent: f64,
}
impl FlocDynamics {
    /// Create a new floc dynamics model.
    pub fn new(mu: f64, aggregation_efficiency: f64, breakup_rate: f64, epsilon_turb: f64) -> Self {
        Self {
            mu,
            aggregation_efficiency,
            breakup_rate,
            epsilon_turb,
            breakup_exponent: 0.5,
        }
    }
    /// Kolmogorov microscale η = (ν³/ε)^(1/4) (m).
    pub fn kolmogorov_scale(&self, rho: f64) -> f64 {
        let nu = self.mu / rho;
        if self.epsilon_turb < 1e-30 {
            return f64::INFINITY;
        }
        (nu.powi(3) / self.epsilon_turb).powf(0.25)
    }
    /// Orthokinetic (turbulent shear) aggregation kernel β (m³/s).
    ///
    /// β = (4/3) * G * (r_i + r_j)³  where G = sqrt(ε/ν).
    pub fn aggregation_kernel(&self, r_i: f64, r_j: f64, rho: f64) -> f64 {
        let nu = self.mu / rho.max(1e-30);
        let g_shear = (self.epsilon_turb / nu.max(1e-30)).sqrt();
        4.0 / 3.0 * g_shear * (r_i + r_j).powi(3)
    }
    /// Aggregation rate (flocs/s) for a population with number density n.
    pub fn aggregation_rate(&self, r_i: f64, r_j: f64, n_i: f64, n_j: f64, rho: f64) -> f64 {
        let beta = self.aggregation_kernel(r_i, r_j, rho);
        self.aggregation_efficiency * beta * n_i * n_j
    }
    /// Turbulent breakup rate (1/s) for a floc of diameter d_f.
    ///
    /// B = B_0 * (ε * d_f² / ν)^p
    pub fn turbulent_breakup_rate(&self, d_floc: f64, rho: f64) -> f64 {
        let nu = self.mu / rho.max(1e-30);
        let arg = self.epsilon_turb * d_floc * d_floc / nu.max(1e-30);
        self.breakup_rate * arg.powf(self.breakup_exponent)
    }
    /// Net growth rate of mean floc size (simplified).
    pub fn net_growth_rate(&self, floc: &Floc, n_floc: f64, rho_fluid: f64) -> f64 {
        let r_f = floc.diameter() / 2.0;
        let agg = self.aggregation_rate(r_f, r_f, n_floc, n_floc, rho_fluid);
        let brk = self.turbulent_breakup_rate(floc.diameter(), rho_fluid) * n_floc;
        agg - brk
    }
}
/// Turbulent dispersion of particles using the eddy lifetime model.
///
/// Accounts for preferential concentration (St > 1) and turbophoresis.
#[derive(Debug, Clone)]
pub struct TurbulentDispersion {
    /// Turbulent kinetic energy k (m²/s²).
    pub k_turb: f64,
    /// Energy dissipation rate ε (m²/s³).
    pub epsilon_turb: f64,
    /// Fluid kinematic viscosity ν (m²/s).
    pub nu: f64,
    /// Particle response time τ_p (s).
    pub tau_p: f64,
}
impl TurbulentDispersion {
    /// Create a new turbulent dispersion model.
    pub fn new(k_turb: f64, epsilon_turb: f64, nu: f64, tau_p: f64) -> Self {
        Self {
            k_turb,
            epsilon_turb,
            nu,
            tau_p,
        }
    }
    /// Integral time scale T_L = 0.6 k / ε (s).
    pub fn integral_time_scale(&self) -> f64 {
        if self.epsilon_turb < 1e-30 {
            return f64::INFINITY;
        }
        0.6 * self.k_turb / self.epsilon_turb
    }
    /// Kolmogorov time scale τ_η = (ν/ε)^(1/2) (s).
    pub fn kolmogorov_time(&self) -> f64 {
        if self.epsilon_turb < 1e-30 {
            return f64::INFINITY;
        }
        (self.nu / self.epsilon_turb).sqrt()
    }
    /// Stokes number based on Kolmogorov time scale.
    pub fn stokes_number_kol(&self) -> f64 {
        let tau_eta = self.kolmogorov_time();
        if tau_eta < 1e-30 {
            return f64::INFINITY;
        }
        self.tau_p / tau_eta
    }
    /// Stokes number based on integral time scale.
    pub fn stokes_number_int(&self) -> f64 {
        let t_l = self.integral_time_scale();
        if t_l < 1e-30 {
            return f64::INFINITY;
        }
        self.tau_p / t_l
    }
    /// Turbulent diffusivity D_t = k * τ_L (simplified mixing length).
    pub fn turbulent_diffusivity(&self) -> f64 {
        let st_int = self.stokes_number_int();
        let t_l = self.integral_time_scale();
        let d_t_fluid = self.k_turb * t_l / 3.0;
        d_t_fluid / (1.0 + st_int)
    }
    /// Turbophoretic velocity toward low-TKE region (m/s).
    ///
    /// v_tp = -τ_p * dK/dy (1D approximation with linear TKE gradient).
    pub fn turbophoretic_velocity(&self, dk_dy: f64) -> f64 {
        -self.tau_p * dk_dy
    }
    /// Preferential concentration factor C_phi (Maxey theory).
    ///
    /// Particles with St ~ 1 concentrate in high-strain regions.
    pub fn preferential_concentration(&self) -> f64 {
        let st = self.stokes_number_kol();
        (-0.5 * (st.ln()).powi(2)).exp()
    }
}
/// Sedimentation front tracker using the Kynch theory.
///
/// Tracks the upward-propagating clear-liquid/suspension interface
/// and the downward-propagating sediment/suspension interface.
#[derive(Debug, Clone)]
pub struct SedimentationFront {
    /// Container height H (m).
    pub height: f64,
    /// Initial uniform concentration φ_0.
    pub phi_0: f64,
    /// Richardson-Zaki exponent n.
    pub rz_exp: f64,
    /// Terminal velocity v_t (m/s, magnitude, downward positive).
    pub v_t: f64,
    /// Current time t (s).
    pub time: f64,
    /// Height of upper clear-liquid interface (measured from bottom).
    pub upper_front: f64,
    /// Height of lower sediment bed top.
    pub lower_front: f64,
}
impl SedimentationFront {
    /// Create a new sedimentation front tracker.
    pub fn new(height: f64, phi_0: f64, rz_exp: f64, v_t: f64) -> Self {
        Self {
            height,
            phi_0,
            rz_exp,
            v_t,
            time: 0.0,
            upper_front: height,
            lower_front: 0.0,
        }
    }
    /// Batch flux function q(φ) = v_t φ (1-φ)^n.
    pub fn batch_flux(&self, phi: f64) -> f64 {
        self.v_t * phi * (1.0 - phi).powf(self.rz_exp)
    }
    /// Propagation speed of a concentration wave at φ.
    ///
    /// dq/dφ = v_t \[(1-φ)^n - n φ (1-φ)^(n-1)\]
    pub fn wave_speed(&self, phi: f64) -> f64 {
        let n = self.rz_exp;
        self.v_t * ((1.0 - phi).powf(n) - n * phi * (1.0 - phi).powf(n - 1.0))
    }
    /// Advance the front tracker by one time step.
    pub fn step(&mut self, dt: f64) {
        let v_upper = self.wave_speed(self.phi_0);
        self.upper_front -= v_upper * dt;
        self.upper_front = self.upper_front.max(self.lower_front);
        let phi_bed = 0.64;
        let flux_bed = self.batch_flux(self.phi_0);
        let dh = flux_bed * dt / (phi_bed - self.phi_0).max(1e-10);
        self.lower_front += dh;
        self.lower_front = self.lower_front.min(self.upper_front);
        self.time += dt;
    }
    /// Has sedimentation completed (fronts have met)?
    pub fn is_complete(&self) -> bool {
        self.upper_front <= self.lower_front + 1e-10
    }
    /// Fraction of domain that is clear liquid.
    pub fn clear_fraction(&self) -> f64 {
        (self.height - self.upper_front) / self.height
    }
    /// Fraction of domain that is settled bed.
    pub fn bed_fraction(&self) -> f64 {
        self.lower_front / self.height
    }
}
/// Turbulent suspension model.
///
/// Models turbophoresis: particle migration in turbulent boundary layer.
#[derive(Debug, Clone)]
pub struct TurbulentSuspension {
    /// Number of grid points in wall-normal direction.
    pub ny: usize,
    /// Wall-normal positions y+.
    pub y_plus: Vec<f64>,
    /// Particle concentration profile.
    pub concentration: Vec<f64>,
    /// Turbophoresis velocity profile.
    pub v_turbo: Vec<f64>,
    /// Stokes number of particles.
    pub st: f64,
    /// Friction velocity u_tau.
    pub u_tau: f64,
    /// Kinematic viscosity.
    pub nu: f64,
    /// Time step (wall units).
    pub dt_plus: f64,
}
impl TurbulentSuspension {
    /// Create a new turbulent suspension model.
    pub fn new(ny: usize, st: f64, u_tau: f64, nu: f64, dt_plus: f64) -> Self {
        let y_plus: Vec<f64> = (0..ny)
            .map(|i| 0.1 + i as f64 * 200.0 / ny as f64)
            .collect();
        let concentration = vec![1.0f64; ny];
        let v_turbo = vec![0.0f64; ny];
        Self {
            ny,
            y_plus,
            concentration,
            v_turbo,
            st,
            u_tau,
            nu,
            dt_plus,
        }
    }
    /// Compute turbophoresis velocity from TKE gradient.
    pub fn compute_turbophoresis(&mut self) {
        for i in 0..self.ny {
            let yp = self.y_plus[i];
            let k = 4.5 * (-0.0018 * yp * yp).exp() * yp * yp + 0.2;
            let dkdy = if i + 1 < self.ny {
                let yp1 = self.y_plus[i + 1];
                let k1 = 4.5 * (-0.0018 * yp1 * yp1).exp() * yp1 * yp1 + 0.2;
                (k1 - k) / (yp1 - yp)
            } else {
                0.0
            };
            self.v_turbo[i] = -self.st * dkdy;
        }
    }
    /// Advance concentration by one step.
    pub fn step(&mut self) {
        self.compute_turbophoresis();
        let c_old = self.concentration.clone();
        for i in 1..self.ny - 1 {
            let dy = self.y_plus[i + 1] - self.y_plus[i - 1];
            let flux_p = self.v_turbo[i] * c_old[i];
            let flux_m = self.v_turbo[i - 1] * c_old[i - 1];
            self.concentration[i] = (c_old[i] - self.dt_plus * (flux_p - flux_m) / dy).max(0.0);
        }
        self.concentration[0] = self.concentration[1];
        self.concentration[self.ny - 1] = self.concentration[self.ny - 2];
    }
    /// Average concentration (should be conserved).
    pub fn avg_concentration(&self) -> f64 {
        self.concentration.iter().sum::<f64>() / self.ny as f64
    }
    /// Near-wall accumulation factor (c_wall / c_avg).
    pub fn near_wall_accumulation(&self) -> f64 {
        let avg = self.avg_concentration();
        if avg < 1e-15 {
            return 1.0;
        }
        self.concentration[0] / avg
    }
}
/// Extended Krieger-Dougherty viscosity model for concentrated suspensions.
///
/// Computes relative viscosity, yield stress, shear thinning and thickening
/// for dense particle suspensions.
#[derive(Debug, Clone)]
pub struct KriegerDougherty {
    /// Maximum packing fraction φ_m.
    pub phi_max: f64,
    /// Intrinsic viscosity \[η\] (typically 2.5 for spheres).
    pub intrinsic_viscosity: f64,
    /// Reference viscosity of continuous phase μ_0 (Pa·s).
    pub mu_0: f64,
}
impl KriegerDougherty {
    /// Create a Krieger-Dougherty model.
    pub fn new(phi_max: f64, intrinsic_viscosity: f64, mu_0: f64) -> Self {
        Self {
            phi_max,
            intrinsic_viscosity,
            mu_0,
        }
    }
    /// Standard spheres in random close packing.
    pub fn spheres(mu_0: f64) -> Self {
        Self::new(0.64, 2.5, mu_0)
    }
    /// Relative viscosity μ_r = (1 - φ/φ_m)^(-\[η\] φ_m).
    pub fn relative_viscosity(&self, phi: f64) -> f64 {
        let phi = phi.clamp(0.0, self.phi_max * 0.999);
        (1.0 - phi / self.phi_max).powf(-self.intrinsic_viscosity * self.phi_max)
    }
    /// Absolute mixture viscosity (Pa·s).
    pub fn mixture_viscosity(&self, phi: f64) -> f64 {
        self.mu_0 * self.relative_viscosity(phi)
    }
    /// Differential viscosity ∂μ_mix/∂φ.
    pub fn viscosity_derivative(&self, phi: f64) -> f64 {
        let phi = phi.clamp(0.0, self.phi_max * 0.999);
        let exp = -self.intrinsic_viscosity * self.phi_max;
        let base = 1.0 - phi / self.phi_max;
        self.mu_0 * exp * base.powf(exp - 1.0) * (-1.0 / self.phi_max)
    }
    /// Onset of shear thickening (critical φ where mu_r doubles).
    ///
    /// Solves (1 - φ/φ_m)^(-\[η\]φ_m) = 2 → φ = φ_m (1 - 2^(-1/(\[η\]φ_m))).
    pub fn shear_thickening_onset(&self) -> f64 {
        let exp = -1.0 / (self.intrinsic_viscosity * self.phi_max);
        self.phi_max * (1.0 - 2.0_f64.powf(exp))
    }
    /// Einstein approximation (dilute limit): μ_r ≈ 1 + 2.5 φ.
    pub fn einstein_viscosity(&self, phi: f64) -> f64 {
        self.mu_0 * (1.0 + self.intrinsic_viscosity * phi)
    }
}
/// Coupled settling and convection simulation.
///
/// Tracks particles settling through a fluid while accounting for
/// drag, buoyancy, and fluid feedback forces.
#[derive(Debug, Clone)]
pub struct ForcedConvectionSedi {
    /// Collection of sediment particles.
    pub particles: Vec<SedimentParticle>,
    /// Fluid density in kg/m³.
    pub rho_fluid: f64,
    /// Dynamic viscosity in Pa·s.
    pub mu_fluid: f64,
    /// Gravity vector \[gx, gy, gz\] in m/s².
    pub gravity: [f64; 3],
    /// Domain size \[Lx, Ly, Lz\] in meters.
    pub domain: [f64; 3],
    /// Time step in seconds.
    pub dt: f64,
}
impl ForcedConvectionSedi {
    /// Create a new forced convection sedimentation setup.
    pub fn new(
        rho_fluid: f64,
        mu_fluid: f64,
        gravity: [f64; 3],
        domain: [f64; 3],
        dt: f64,
    ) -> Self {
        Self {
            particles: Vec::new(),
            rho_fluid,
            mu_fluid,
            gravity,
            domain,
            dt,
        }
    }
    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, particle: SedimentParticle) {
        self.particles.push(particle);
    }
    /// Compute drag force on particle.
    pub fn drag_force(&self, p: &SedimentParticle) -> [f64; 3] {
        let re = 2.0 * p.radius * self.rho_fluid * mag3(p.velocity) / self.mu_fluid;
        let cd = drag_coeff_particle(re);
        let area = std::f64::consts::PI * p.radius * p.radius;
        let fmag = 0.5 * self.rho_fluid * cd * area * mag3(p.velocity);
        let v = p.velocity;
        let vmag = mag3(v).max(1e-30);
        [
            -fmag * v[0] / vmag,
            -fmag * v[1] / vmag,
            -fmag * v[2] / vmag,
        ]
    }
    /// Compute buoyancy force on particle.
    pub fn buoyancy_force(&self, p: &SedimentParticle) -> [f64; 3] {
        let vol = p.volume();
        [
            -self.rho_fluid * self.gravity[0] * vol,
            -self.rho_fluid * self.gravity[1] * vol,
            -self.rho_fluid * self.gravity[2] * vol,
        ]
    }
    /// Compute gravity force on particle.
    pub fn gravity_force(&self, p: &SedimentParticle) -> [f64; 3] {
        let mass = p.mass();
        [
            mass * self.gravity[0],
            mass * self.gravity[1],
            mass * self.gravity[2],
        ]
    }
    /// Compute net force on particle (gravity + buoyancy + drag).
    pub fn net_force(&self, p: &SedimentParticle) -> [f64; 3] {
        let fg = self.gravity_force(p);
        let fb = self.buoyancy_force(p);
        let fd = self.drag_force(p);
        [
            fg[0] + fb[0] + fd[0],
            fg[1] + fb[1] + fd[1],
            fg[2] + fb[2] + fd[2],
        ]
    }
    /// Advance all particles by one time step.
    pub fn step(&mut self) {
        let dt = self.dt;
        let forces: Vec<[f64; 3]> = self.particles.iter().map(|p| self.net_force(p)).collect();
        for (i, p) in self.particles.iter_mut().enumerate() {
            let mass = p.mass();
            let f = forces[i];
            p.velocity[0] += f[0] / mass * dt;
            p.velocity[1] += f[1] / mass * dt;
            p.velocity[2] += f[2] / mass * dt;
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
        }
    }
    /// Get terminal velocity for a particle (iterative).
    pub fn terminal_velocity(&self, p: &SedimentParticle) -> f64 {
        stokes_settling(p.radius, p.density, self.rho_fluid, self.mu_fluid, 9.81)
    }
}
/// Fluidization regime enumeration.
#[derive(Debug, Clone, PartialEq)]
pub enum FluidizationRegime {
    /// Packed bed (u < umf).
    PackedBed,
    /// Bubbling bed.
    BubblingBed,
    /// Slug flow.
    SlugFlow,
    /// Fast fluidization.
    FastFluidization,
}
/// Volume fraction field with concentration gradient and diffusion.
///
/// Tracks the evolution of φ(x,t) under settling, diffusion and
/// hindered settling corrections.
#[derive(Debug, Clone)]
pub struct VolumeFractionField {
    /// 1D volume fraction array (vertical column).
    pub phi: Vec<f64>,
    /// Number of grid cells.
    pub n: usize,
    /// Cell height Δy (m).
    pub dy: f64,
    /// Richardson-Zaki exponent.
    pub rz_exp: f64,
    /// Terminal settling velocity (m/s).
    pub v_t: f64,
    /// Particle diffusivity D_p (m²/s).
    pub diffusivity: f64,
}
impl VolumeFractionField {
    /// Create a new volume fraction field.
    pub fn new(n: usize, dy: f64, v_t: f64, rz_exp: f64, diffusivity: f64) -> Self {
        Self {
            phi: vec![0.0; n],
            n,
            dy,
            rz_exp,
            v_t,
            diffusivity,
        }
    }
    /// Set uniform initial concentration.
    pub fn set_uniform(&mut self, phi_0: f64) {
        for x in self.phi.iter_mut() {
            *x = phi_0;
        }
    }
    /// Hindered settling flux at cell i.
    pub fn flux(&self, i: usize) -> f64 {
        let phi_i = self.phi[i];
        phi_i * self.v_t * (1.0 - phi_i).powf(self.rz_exp)
    }
    /// Advance by one time step using upwind finite-volume scheme.
    pub fn step(&mut self, dt: f64) {
        let phi_old = self.phi.clone();
        for i in 1..self.n - 1 {
            let f_i = phi_old[i] * self.v_t * (1.0 - phi_old[i]).powf(self.rz_exp);
            let f_im1 = phi_old[i - 1] * self.v_t * (1.0 - phi_old[i - 1]).powf(self.rz_exp);
            let advection = (f_i - f_im1) / self.dy;
            let diffusion = self.diffusivity * (phi_old[i + 1] - 2.0 * phi_old[i] + phi_old[i - 1])
                / (self.dy * self.dy);
            self.phi[i] = (phi_old[i] - dt * advection + dt * diffusion).clamp(0.0, 1.0);
        }
        self.phi[0] = self.phi[1];
        self.phi[self.n - 1] = self.phi[self.n - 2];
    }
    /// Total solid volume (sum φ_i Δy).
    pub fn total_solid_volume(&self) -> f64 {
        self.phi.iter().sum::<f64>() * self.dy
    }
    /// Sedimentation front position: highest cell with φ > threshold.
    pub fn front_position(&self, threshold: f64) -> f64 {
        for i in (0..self.n).rev() {
            if self.phi[i] > threshold {
                return i as f64 * self.dy;
            }
        }
        0.0
    }
    /// Bed height: lowest column where φ > packing threshold.
    pub fn bed_height(&self, phi_pack: f64) -> f64 {
        let mut h = 0.0;
        for i in 0..self.n {
            if self.phi[i] > phi_pack {
                h = (i + 1) as f64 * self.dy;
            }
        }
        h
    }
}
/// Statistics for sedimentation simulations.
///
/// Computes settling velocity distributions, polydisperse size effects, and Stokes numbers.
#[derive(Debug, Clone)]
pub struct SediStatistics {
    /// Particle radii array.
    pub radii: Vec<f64>,
    /// Particle densities array.
    pub densities: Vec<f64>,
    /// Fluid density.
    pub rho_fluid: f64,
    /// Fluid viscosity.
    pub mu: f64,
    /// Gravity.
    pub gravity: f64,
    /// Characteristic length scale for Stokes number.
    pub l_char: f64,
    /// Characteristic velocity for Stokes number.
    pub u_char: f64,
}
impl SediStatistics {
    /// Create new sedimentation statistics.
    pub fn new(rho_fluid: f64, mu: f64, gravity: f64, l_char: f64, u_char: f64) -> Self {
        Self {
            radii: Vec::new(),
            densities: Vec::new(),
            rho_fluid,
            mu,
            gravity,
            l_char,
            u_char,
        }
    }
    /// Add a particle size/density.
    pub fn add_particle(&mut self, radius: f64, density: f64) {
        self.radii.push(radius);
        self.densities.push(density);
    }
    /// Compute settling velocities for all particles.
    pub fn settling_velocities(&self) -> Vec<f64> {
        self.radii
            .iter()
            .zip(self.densities.iter())
            .map(|(&r, &rho)| stokes_settling(r, rho, self.rho_fluid, self.mu, self.gravity))
            .collect()
    }
    /// Compute Stokes numbers for all particles.
    pub fn stokes_numbers(&self) -> Vec<f64> {
        let vs = self.settling_velocities();
        vs.iter()
            .zip(self.radii.iter())
            .map(|(&_v, &r)| {
                let tau_p = 2.0 * self.densities[0] * r * r / (9.0 * self.mu);
                tau_p * self.u_char / self.l_char
            })
            .collect()
    }
    /// Mean settling velocity.
    pub fn mean_settling_velocity(&self) -> f64 {
        let vs = self.settling_velocities();
        vs.iter().sum::<f64>() / vs.len().max(1) as f64
    }
    /// Standard deviation of settling velocities.
    pub fn std_settling_velocity(&self) -> f64 {
        let vs = self.settling_velocities();
        if vs.is_empty() {
            return 0.0;
        }
        let mean = vs.iter().sum::<f64>() / vs.len() as f64;
        let var = vs.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vs.len() as f64;
        var.sqrt()
    }
    /// Segregation index: ratio of std to mean.
    pub fn segregation_index(&self) -> f64 {
        let mean = self.mean_settling_velocity();
        if mean.abs() < 1e-15 {
            return 0.0;
        }
        self.std_settling_velocity() / mean.abs()
    }
}
/// Coupled LBM-DEM sedimentation solver (2D slice representation).
///
/// Couples a D2Q9 LBM fluid field with an ensemble of DEM particles
/// using an immersed boundary–type momentum exchange.
#[derive(Debug, Clone)]
pub struct LbmDemSolver {
    /// D2Q9 LBM fluid solver.
    pub lbm: SediLbm,
    /// DEM particles.
    pub particles: Vec<DemParticle>,
    /// Lattice spacing Δx (m).
    pub dx: f64,
    /// Time step Δt (s).
    pub dt: f64,
    /// Gravity vector \[gx, gy\] in lattice units.
    pub gravity_latt: [f64; 2],
    /// Fluid density (physical).
    pub rho_fluid: f64,
}
impl LbmDemSolver {
    /// Create a new LBM-DEM solver.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        rho_fluid: f64,
        rho_solid: f64,
        dx: f64,
        dt: f64,
    ) -> Self {
        let lbm = SediLbm::new(nx, ny, omega, rho_fluid, rho_solid, 9.81);
        Self {
            lbm,
            particles: Vec::new(),
            dx,
            dt,
            gravity_latt: [0.0, -9.81 * dt * dt / dx],
            rho_fluid,
        }
    }
    /// Add a DEM particle.
    pub fn add_particle(&mut self, p: DemParticle) {
        self.particles.push(p);
    }
    /// Map particle position to LBM grid cell.
    pub fn grid_cell(&self, pos: [f64; 3]) -> (usize, usize) {
        let ix = (pos[0] / self.dx).floor() as i64;
        let iy = (pos[1] / self.dx).floor() as i64;
        let ix = ix.clamp(0, self.lbm.nx as i64 - 1) as usize;
        let iy = iy.clamp(0, self.lbm.ny as i64 - 1) as usize;
        (ix, iy)
    }
    /// Compute fluid-particle interaction force (simplified direct forcing).
    ///
    /// Returns hydrodynamic force \[fx, fy\] on particle in physical units.
    pub fn hydro_force(&self, p: &DemParticle) -> [f64; 2] {
        let (ix, iy) = self.grid_cell(p.position);
        let u_f = self.lbm.u[iy][ix];
        let rho_f = self.lbm.rho[iy][ix];
        let nu = (1.0 / self.lbm.omega - 0.5) / 3.0;
        let mu = rho_f * nu * self.dx * self.dx / self.dt;
        let factor = 6.0 * std::f64::consts::PI * mu * p.radius;
        [
            factor * (u_f[0] * self.dx / self.dt - p.velocity[0]),
            factor * (u_f[1] * self.dx / self.dt - p.velocity[1]),
        ]
    }
    /// Update volume fraction field from particle positions.
    pub fn update_phi(&mut self) {
        for y in 0..self.lbm.ny {
            for x in 0..self.lbm.nx {
                self.lbm.phi[y][x] = 0.0;
            }
        }
        let dx = self.dx;
        let cell_vol = dx * dx;
        for p in &self.particles {
            let (ix, iy) = self.grid_cell(p.position);
            let pvol = std::f64::consts::PI * p.radius * p.radius;
            self.lbm.phi[iy][ix] = (self.lbm.phi[iy][ix] + pvol / cell_vol).min(0.99);
        }
    }
    /// Advance one time step.
    pub fn step(&mut self) {
        self.update_phi();
        self.lbm.step();
        let forces: Vec<[f64; 2]> = self.particles.iter().map(|p| self.hydro_force(p)).collect();
        for (i, p) in self.particles.iter_mut().enumerate() {
            let fh = forces[i];
            p.force[0] += fh[0];
            p.force[1] += fh[1];
            p.integrate([0.0, -9.81, 0.0], self.rho_fluid, self.dt);
        }
    }
}
/// Gravity current (density current) driven by suspended sediment.
///
/// Models the front speed, head height, and mixing of a turbidity current.
#[derive(Debug, Clone)]
pub struct DensityCurrent {
    /// Ambient fluid density ρ_a (kg/m³).
    pub rho_ambient: f64,
    /// Current density ρ_c (kg/m³).
    pub rho_current: f64,
    /// Current head height h (m).
    pub head_height: f64,
    /// Current length L (m).
    pub length: f64,
    /// Gravitational acceleration g (m/s²).
    pub g: f64,
    /// Friction factor f_r.
    pub friction: f64,
    /// Current time t (s).
    pub time: f64,
}
impl DensityCurrent {
    /// Create a new density current.
    pub fn new(rho_ambient: f64, rho_current: f64, head_height: f64, g: f64) -> Self {
        Self {
            rho_ambient,
            rho_current,
            head_height,
            length: 0.0,
            g,
            friction: 0.01,
            time: 0.0,
        }
    }
    /// Reduced gravity g' = g (ρ_c - ρ_a) / ρ_a (m/s²).
    pub fn reduced_gravity(&self) -> f64 {
        let dr = self.rho_current - self.rho_ambient;
        self.g * dr / self.rho_ambient.max(1e-30)
    }
    /// Froude number Fr = u_f / sqrt(g' h).
    pub fn front_froude_number(&self) -> f64 {
        0.5_f64.sqrt()
    }
    /// Front speed u_f = Fr * sqrt(g' h) (m/s).
    pub fn front_speed(&self) -> f64 {
        let gp = self.reduced_gravity();
        if gp < 0.0 {
            return 0.0;
        }
        self.front_froude_number() * (gp * self.head_height).sqrt()
    }
    /// Advance the current head by one time step.
    pub fn step(&mut self, dt: f64) {
        let u_f = self.front_speed();
        self.length += u_f * dt;
        self.time += dt;
        let dilution_rate = self.friction * u_f / self.head_height.max(1e-30);
        self.rho_current -= (self.rho_current - self.rho_ambient) * dilution_rate * dt;
    }
    /// Richardson number Ri = g' h / u_f² (bulk).
    pub fn richardson_number(&self) -> f64 {
        let u = self.front_speed();
        if u < 1e-30 {
            return f64::INFINITY;
        }
        self.reduced_gravity() * self.head_height / (u * u)
    }
    /// Deposit mass flux from turbidity current (simplified).
    ///
    /// Ws = settling velocity, C = volume fraction.
    pub fn deposit_flux(&self, settling_vel: f64, concentration: f64) -> f64 {
        settling_vel * concentration
    }
}
/// DEM (Discrete Element Method) particle for LBM-DEM coupling.
///
/// Stores full 6-DOF dynamics state for a spherical particle.
#[derive(Debug, Clone)]
pub struct DemParticle {
    /// Position \[x, y, z\] m.
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] m/s.
    pub velocity: [f64; 3],
    /// Angular velocity \[ωx, ωy, ωz\] rad/s.
    pub omega: [f64; 3],
    /// Radius m.
    pub radius: f64,
    /// Density kg/m³.
    pub density: f64,
    /// Accumulated force \[fx, fy, fz\] N.
    pub force: [f64; 3],
    /// Accumulated torque \[tx, ty, tz\] N·m.
    pub torque: [f64; 3],
    /// Particle identifier.
    pub id: usize,
}
impl DemParticle {
    /// Create a new DEM particle.
    pub fn new(position: [f64; 3], radius: f64, density: f64, id: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            omega: [0.0; 3],
            radius,
            density,
            force: [0.0; 3],
            torque: [0.0; 3],
            id,
        }
    }
    /// Particle mass kg.
    pub fn mass(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3) * self.density
    }
    /// Moment of inertia (solid sphere) kg·m².
    pub fn inertia(&self) -> f64 {
        2.0 / 5.0 * self.mass() * self.radius * self.radius
    }
    /// Add hydrodynamic force and torque from LBM.
    pub fn add_hydro_force(&mut self, force: [f64; 3], torque: [f64; 3]) {
        for k in 0..3 {
            self.force[k] += force[k];
            self.torque[k] += torque[k];
        }
    }
    /// Integrate equations of motion using Verlet (one step).
    pub fn integrate(&mut self, gravity: [f64; 3], rho_fluid: f64, dt: f64) {
        let mass = self.mass();
        let vol = 4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3);
        let buoy_frac = rho_fluid / self.density;
        let _ = vol;
        for (k, (vel_k, pos_k)) in self
            .velocity
            .iter_mut()
            .zip(self.position.iter_mut())
            .enumerate()
        {
            let net_f = self.force[k] + mass * gravity[k] * (1.0 - buoy_frac);
            let acc = net_f / mass;
            *vel_k += acc * dt;
            *pos_k += *vel_k * dt;
        }
        let inertia = self.inertia();
        for (omega_k, &torque_k) in self.omega.iter_mut().zip(self.torque.iter()) {
            let alpha = torque_k / inertia;
            *omega_k += alpha * dt;
        }
        self.force = [0.0; 3];
        self.torque = [0.0; 3];
    }
    /// Detect collision with another particle; returns overlap distance.
    pub fn overlap(&self, other: &DemParticle) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let min_dist = self.radius + other.radius;
        (min_dist - dist).max(0.0)
    }
}
/// Fluidization LBM simulation.
///
/// Models minimum fluidization velocity, bubbling bed, and slug flow.
#[derive(Debug, Clone)]
pub struct FluidizationLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Solid volume fraction field.
    pub phi: Vec<Vec<f64>>,
    /// Fluid velocity field u\[y\]\[x\]\[2\].
    pub u: Vec<Vec<[f64; 2]>>,
    /// Particle diameter in meters.
    pub dp: f64,
    /// Particle density kg/m³.
    pub rho_solid: f64,
    /// Fluid density kg/m³.
    pub rho_fluid: f64,
    /// Fluid viscosity Pa·s.
    pub mu: f64,
    /// Inlet superficial velocity m/s.
    pub u_inlet: f64,
    /// Minimum fluidization velocity.
    pub u_mf: f64,
}
impl FluidizationLbm {
    /// Create a new fluidization LBM simulation.
    pub fn new(
        nx: usize,
        ny: usize,
        dp: f64,
        rho_solid: f64,
        rho_fluid: f64,
        mu: f64,
        u_inlet: f64,
    ) -> Self {
        let phi_init = 0.6;
        let phi = vec![vec![phi_init; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let u_mf = Self::calc_umf(dp, rho_solid, rho_fluid, mu, phi_init);
        Self {
            nx,
            ny,
            phi,
            u,
            dp,
            rho_solid,
            rho_fluid,
            mu,
            u_inlet,
            u_mf,
        }
    }
    /// Compute minimum fluidization velocity using Ergun equation.
    pub fn calc_umf(dp: f64, rho_solid: f64, rho_fluid: f64, mu: f64, phi_mf: f64) -> f64 {
        let eps = 1.0 - phi_mf;
        let drho = rho_solid - rho_fluid;
        let g = 9.81;
        let ar = rho_fluid * drho * g * dp.powi(3) / mu.powi(2);
        let re_mf = (33.7_f64.powi(2) + 0.0408 * ar).sqrt() - 33.7;
        let _ = eps;
        re_mf * mu / (rho_fluid * dp)
    }
    /// Check if bed is fluidized.
    pub fn is_fluidized(&self) -> bool {
        self.u_inlet > self.u_mf
    }
    /// Fluidization regime classification.
    pub fn regime(&self) -> FluidizationRegime {
        let ratio = self.u_inlet / self.u_mf;
        if ratio < 1.0 {
            FluidizationRegime::PackedBed
        } else if ratio < 5.0 {
            FluidizationRegime::BubblingBed
        } else if ratio < 20.0 {
            FluidizationRegime::SlugFlow
        } else {
            FluidizationRegime::FastFluidization
        }
    }
    /// Void fraction at fluidization.
    pub fn epsilon_mf(&self) -> f64 {
        let re = self.rho_fluid * self.u_inlet * self.dp / self.mu;
        let ar = self.rho_fluid * (self.rho_solid - self.rho_fluid) * 9.81 * self.dp.powi(3)
            / self.mu.powi(2);
        let _ = re;
        (0.0408 * ar).powf(1.0 / 3.0) / (self.dp * 100.0 + 1.0).max(0.4)
    }
    /// Update phi field (simple settling model).
    pub fn update_phi(&mut self) {
        let eps_mf = 1.0 - 0.4;
        for y in 0..self.ny {
            for x in 0..self.nx {
                if self.is_fluidized() {
                    self.phi[y][x] = (self.phi[y][x] * 0.99).max(0.1);
                } else {
                    self.phi[y][x] = (self.phi[y][x] + 0.001).min(1.0 - eps_mf);
                }
            }
        }
    }
}
/// Richardson-Zaki hindered settling factor.
///
/// Computes the effective settling velocity using the Richardson-Zaki correlation:
/// Ut_eff = Ut * (1 - phi)^n
#[derive(Debug, Clone)]
pub struct HindranceFactor {
    /// Terminal settling velocity of a single particle in m/s.
    pub terminal_velocity: f64,
    /// Richardson-Zaki exponent (typically 4.65 for Re < 0.2).
    pub exponent: f64,
}
impl HindranceFactor {
    /// Create a new hindrance factor calculator.
    pub fn new(terminal_velocity: f64, exponent: f64) -> Self {
        Self {
            terminal_velocity,
            exponent,
        }
    }
    /// Compute hindered settling velocity at volume fraction phi.
    pub fn hindered_velocity(&self, phi: f64) -> f64 {
        let phi = phi.clamp(0.0, 1.0);
        self.terminal_velocity * (1.0 - phi).powf(self.exponent)
    }
    /// Compute flux density at volume fraction phi.
    pub fn batch_flux(&self, phi: f64) -> f64 {
        phi * self.hindered_velocity(phi)
    }
    /// Compute maximum flux (Kynch flux).
    pub fn max_flux_phi(&self) -> f64 {
        1.0 / (1.0 + self.exponent)
    }
}
/// Free-settling velocity correction model.
///
/// Extends Stokes settling to account for particle shape (Corey shape factor),
/// wall effects, and non-sphericity.
#[derive(Debug, Clone)]
pub struct FreeSettlingCorrection {
    /// Corey shape factor ψ (= 1 for sphere, < 1 for irregular).
    pub corey_shape: f64,
    /// Particle sphericity Φ (= 1 for sphere).
    pub sphericity: f64,
    /// Wall correction factor Fw (1 for no wall effects).
    pub wall_correction: f64,
}
impl FreeSettlingCorrection {
    /// Create a correction model (sphere defaults).
    pub fn new(corey_shape: f64, sphericity: f64, wall_correction: f64) -> Self {
        Self {
            corey_shape,
            sphericity,
            wall_correction,
        }
    }
    /// Spherical particle (no corrections).
    pub fn sphere() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }
    /// Drag coefficient for non-spherical particle (Haider-Levenspiel).
    ///
    /// Uses the modified correlation that accounts for sphericity φ.
    pub fn drag_coeff_nonspherical(&self, re: f64) -> f64 {
        if re < 1e-10 {
            return 1e6;
        }
        let phi = self.sphericity;
        let a = (-2.3288 + 6.4581 * phi - 2.4486 * phi * phi).exp();
        let b = 0.0964 + 0.5565 * phi;
        let c = (4.905 - 13.8944 * phi + 18.4222 * phi * phi - 10.2599 * phi.powi(3)).exp();
        let d = (1.4681 + 12.2584 * phi - 20.7322 * phi * phi + 15.8855 * phi.powi(3)).exp();
        24.0 / re * (1.0 + a * re.powf(b)) + c / (1.0 + d / re)
    }
    /// Corrected terminal velocity ratio v_t/v_t_sphere.
    pub fn velocity_correction_factor(&self, re: f64) -> f64 {
        let cd_sphere = drag_coeff_particle(re);
        let cd_ns = self.drag_coeff_nonspherical(re);
        if cd_ns < 1e-30 {
            return 1.0;
        }
        (cd_sphere / cd_ns).sqrt() * self.wall_correction * self.corey_shape
    }
    /// Corrected settling velocity.
    pub fn corrected_velocity(&self, stokes_vel: f64, re: f64) -> f64 {
        stokes_vel * self.velocity_correction_factor(re)
    }
}
/// Full 3D settling tank simulation.
///
/// Handles creaming, sedimentation, and cake compression.
#[derive(Debug, Clone)]
pub struct Sedimentation3D {
    /// Grid dimensions \[nx, ny, nz\].
    pub dims: [usize; 3],
    /// Volume fraction field.
    pub phi: Vec<f64>,
    /// Settling velocity field.
    pub vs: Vec<f64>,
    /// Particle diameter.
    pub dp: f64,
    /// Particle density.
    pub rho_p: f64,
    /// Fluid density.
    pub rho_f: f64,
    /// Fluid viscosity.
    pub mu: f64,
    /// Hindered settling exponent.
    pub rz_exp: f64,
    /// Time step.
    pub dt: f64,
}
impl Sedimentation3D {
    /// Create a new 3D sedimentation simulation.
    pub fn new(dims: [usize; 3], dp: f64, rho_p: f64, rho_f: f64, mu: f64, dt: f64) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        let phi = vec![0.1f64; n];
        let vs = vec![0.0f64; n];
        let rz_exp = richardson_zaki_exponent(
            2.0 * rho_f * stokes_settling(dp / 2.0, rho_p, rho_f, mu, 9.81) * (dp / 2.0) / mu,
        );
        Self {
            dims,
            phi,
            vs,
            dp,
            rho_p,
            rho_f,
            mu,
            rz_exp,
            dt,
        }
    }
    /// Index into flat array.
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.dims[1] * self.dims[0] + y * self.dims[0] + x
    }
    /// Compute terminal settling velocity (single particle).
    pub fn vt_single(&self) -> f64 {
        stokes_settling(self.dp / 2.0, self.rho_p, self.rho_f, self.mu, 9.81)
    }
    /// Update settling velocities using Richardson-Zaki.
    pub fn update_settling(&mut self) {
        let vt = self.vt_single();
        for i in 0..self.phi.len() {
            let phi = self.phi[i];
            self.vs[i] = vt * (1.0 - phi).powf(self.rz_exp);
        }
    }
    /// Advance sedimentation by one time step (1D vertical settling in y).
    pub fn step(&mut self) {
        self.update_settling();
        let [nx, ny, _nz] = self.dims;
        let dt = self.dt;
        let _nx = nx;
        let phi_old = self.phi.clone();
        for z in 0..self.dims[2] {
            for x in 0..self.dims[0] {
                for y in 1..ny {
                    let i = self.idx(x, y, z);
                    let i_up = self.idx(x, y - 1, z);
                    let flux = phi_old[i] * self.vs[i];
                    self.phi[i] = (phi_old[i] + dt * flux / 0.01).clamp(0.0, 0.99);
                    let _ = i_up;
                }
            }
        }
    }
    /// Total volume of settled particles.
    pub fn settled_volume(&self, domain_cell_vol: f64) -> f64 {
        self.phi.iter().sum::<f64>() * domain_cell_vol
    }
    /// Check if creaming (upward floating) is occurring.
    pub fn is_creaming(&self) -> bool {
        self.rho_p < self.rho_f
    }
}
/// Floc (aggregate) particle for cohesive sediment simulation.
///
/// A floc is a porous aggregate of primary particles held together
/// by inter-particle forces. Its effective density and drag differ
/// from solid spheres.
#[derive(Debug, Clone)]
pub struct Floc {
    /// Number of primary particles in floc.
    pub n_primary: f64,
    /// Primary particle diameter in meters.
    pub d_primary: f64,
    /// Fractal dimension D_f (typically 1.7–2.5).
    pub fractal_dim: f64,
    /// Floc mass in kg.
    pub mass: f64,
    /// Floc velocity \[vx, vy, vz\] m/s.
    pub velocity: [f64; 3],
    /// Floc position \[x, y, z\] m.
    pub position: [f64; 3],
}
impl Floc {
    /// Create a new floc from primary particle properties.
    pub fn new(n_primary: f64, d_primary: f64, fractal_dim: f64, rho_primary: f64) -> Self {
        let r_p = d_primary / 2.0;
        let vol_primary = 4.0 / 3.0 * std::f64::consts::PI * r_p.powi(3);
        let mass = n_primary * rho_primary * vol_primary;
        Self {
            n_primary,
            d_primary,
            fractal_dim,
            mass,
            velocity: [0.0; 3],
            position: [0.0; 3],
        }
    }
    /// Floc diameter via fractal scaling: d_f = d_p * N^(1/D_f).
    pub fn diameter(&self) -> f64 {
        self.d_primary * self.n_primary.powf(1.0 / self.fractal_dim)
    }
    /// Effective floc density (kg/m³).
    ///
    /// ρ_eff = ρ_fluid + (ρ_p - ρ_fluid) * (d_f/d_p)^(D_f - 3)
    pub fn effective_density(&self, rho_fluid: f64, rho_primary: f64) -> f64 {
        let ratio = self.diameter() / self.d_primary;
        let exponent = self.fractal_dim - 3.0;
        rho_fluid + (rho_primary - rho_fluid) * ratio.powf(exponent)
    }
    /// Floc volume (m³) from diameter.
    pub fn volume(&self) -> f64 {
        let r = self.diameter() / 2.0;
        4.0 / 3.0 * std::f64::consts::PI * r.powi(3)
    }
    /// Settling velocity of floc using effective density and Stokes.
    pub fn settling_velocity(&self, rho_fluid: f64, rho_primary: f64, mu: f64, g: f64) -> f64 {
        let rho_eff = self.effective_density(rho_fluid, rho_primary);
        stokes_settling(self.diameter() / 2.0, rho_eff, rho_fluid, mu, g)
    }
}
/// A particle in a sedimentation simulation.
///
/// Stores position, velocity, radius, density and local volume fraction.
#[derive(Debug, Clone)]
pub struct SedimentParticle {
    /// Position \[x, y, z\] in meters.
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] in m/s.
    pub velocity: [f64; 3],
    /// Particle radius in meters.
    pub radius: f64,
    /// Particle density in kg/m³.
    pub density: f64,
    /// Local solid volume fraction (0..1).
    pub volume_fraction: f64,
}
impl SedimentParticle {
    /// Create a new sediment particle.
    pub fn new(position: [f64; 3], radius: f64, density: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            radius,
            density,
            volume_fraction: 0.0,
        }
    }
    /// Particle mass in kg.
    pub fn mass(&self) -> f64 {
        let vol = 4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3);
        self.density * vol
    }
    /// Particle volume in m³.
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3)
    }
}
/// Sediment bed compaction model (self-weight consolidation).
///
/// Models the time-dependent compaction of a freshly deposited sediment
/// bed under its own weight, with permeability that depends on void ratio.
#[derive(Debug, Clone)]
pub struct BedCompaction {
    /// Initial void ratio e_0.
    pub void_ratio: f64,
    /// Compression index C_c (from consolidation test).
    pub compression_index: f64,
    /// Initial permeability K_0 (m²) at void ratio e_0.
    pub permeability_0: f64,
    /// Kozeny-Carman exponent for permeability-void ratio relation.
    pub kozeny_exp: f64,
    /// Effective stress σ' (Pa).
    pub effective_stress: f64,
    /// Particle density ρ_p (kg/m³).
    pub rho_p: f64,
    /// Fluid density ρ_f (kg/m³).
    pub rho_f: f64,
    /// Bed height H (m).
    pub bed_height: f64,
    /// Gravity g (m/s²).
    pub g: f64,
}
impl BedCompaction {
    /// Create a new bed compaction model.
    pub fn new(
        void_ratio: f64,
        compression_index: f64,
        permeability_0: f64,
        rho_p: f64,
        rho_f: f64,
        bed_height: f64,
        g: f64,
    ) -> Self {
        let buoyant_weight = (rho_p - rho_f) / (1.0 + void_ratio) * g * bed_height;
        Self {
            void_ratio,
            compression_index,
            permeability_0,
            kozeny_exp: 3.0,
            effective_stress: buoyant_weight,
            rho_p,
            rho_f,
            bed_height,
            g,
        }
    }
    /// Current porosity from void ratio: n = e / (1 + e).
    pub fn porosity(&self) -> f64 {
        self.void_ratio / (1.0 + self.void_ratio)
    }
    /// Permeability via Kozeny-Carman (void ratio form).
    ///
    /// K = K_0 * (e/e_0)^kozeny_exp
    pub fn permeability(&self, e_0: f64) -> f64 {
        if e_0 < 1e-30 {
            return 0.0;
        }
        self.permeability_0 * (self.void_ratio / e_0).powf(self.kozeny_exp)
    }
    /// Compress bed by increasing effective stress and reducing void ratio.
    ///
    /// Uses Terzaghi 1D consolidation: Δe = -C_c log10(σ_new/σ_old).
    pub fn consolidate(&mut self, delta_stress: f64) {
        if delta_stress <= 0.0 || self.effective_stress < 1e-10 {
            return;
        }
        let sigma_new = self.effective_stress + delta_stress;
        let de = -self.compression_index * (sigma_new / self.effective_stress).log10();
        self.void_ratio = (self.void_ratio + de).max(0.01);
        self.effective_stress = sigma_new;
        self.bed_height *= (1.0 + self.void_ratio) / (1.0 + self.void_ratio - de);
    }
    /// Time to consolidate 90% (Terzaghi): T_v90 = 0.848 H²/c_v.
    ///
    /// `c_v` coefficient of consolidation (m²/s).
    pub fn time_to_90_consolidation(&self, c_v: f64) -> f64 {
        if c_v < 1e-30 {
            return f64::INFINITY;
        }
        0.848 * self.bed_height.powi(2) / c_v
    }
    /// Solid volume fraction in the bed.
    pub fn solid_fraction(&self) -> f64 {
        1.0 / (1.0 + self.void_ratio)
    }
}
/// Two-phase LBM with mixture model for solid-liquid flow.
///
/// Uses a single-fluid approach where the mixture density and viscosity
/// are computed from the local solid volume fraction.
#[derive(Debug, Clone)]
pub struct SediLbm {
    /// Grid width (number of cells in x).
    pub nx: usize,
    /// Grid height (number of cells in y).
    pub ny: usize,
    /// Relaxation parameter omega for fluid phase.
    pub omega: f64,
    /// Distribution functions f\[y\]\[x\]\[q\] for D2Q9.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Solid volume fraction field phi\[y\]\[x\].
    pub phi: Vec<Vec<f64>>,
    /// Density field rho\[y\]\[x\].
    pub rho: Vec<Vec<f64>>,
    /// Velocity field u\[y\]\[x\]\[2\].
    pub u: Vec<Vec<[f64; 2]>>,
    /// Particle density (solid phase density).
    pub rho_solid: f64,
    /// Fluid phase density.
    pub rho_fluid: f64,
    /// Gravity magnitude (m/s²).
    pub gravity: f64,
}
impl SediLbm {
    /// Create a new SediLbm simulation.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        rho_fluid: f64,
        rho_solid: f64,
        gravity: f64,
    ) -> Self {
        let f = vec![vec![vec![0.0f64; 9]; nx]; ny];
        let phi = vec![vec![0.0f64; nx]; ny];
        let rho = vec![vec![rho_fluid; nx]; ny];
        let u = vec![vec![[0.0f64; 2]; nx]; ny];
        let mut sim = Self {
            nx,
            ny,
            omega,
            f,
            phi,
            rho,
            u,
            rho_solid,
            rho_fluid,
            gravity,
        };
        sim.initialize_equilibrium();
        sim
    }
    /// Initialize all distributions to equilibrium with zero velocity.
    pub fn initialize_equilibrium(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let rho = self.rho[y][x];
                let feq = feq_d2q9(rho, 0.0, 0.0);
                self.f[y][x] = feq.to_vec();
            }
        }
    }
    /// Compute mixture density and viscosity from volume fraction.
    pub fn mixture_density(&self, phi: f64) -> f64 {
        self.rho_fluid * (1.0 - phi) + self.rho_solid * phi
    }
    /// Compute mixture relaxation from volume fraction using Krieger-Dougherty.
    pub fn mixture_omega(&self, phi: f64) -> f64 {
        let mu_rel = mixture_viscosity(phi, 0.64);
        let nu_mix = (1.0 / self.omega - 0.5) / mu_rel;
        1.0 / (nu_mix + 0.5)
    }
    /// Perform one BGK collision step with sedimentation body force.
    pub fn collide(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let phi = self.phi[y][x];
                let rho = self.rho[y][x];
                let [ux, uy] = self.u[y][x];
                let omega = self.mixture_omega(phi);
                let feq = feq_d2q9(rho, ux, uy);
                let buoy = (self.rho_solid - self.rho_fluid) * phi * self.gravity;
                for q in 0..9 {
                    let force_q = W9[q] * 3.0 * (-buoy) * EY9[q] as f64;
                    self.f[y][x][q] = self.f[y][x][q] * (1.0 - omega) + omega * feq[q] + force_q;
                }
            }
        }
    }
    /// Perform one streaming step with periodic boundaries.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![vec![0.0f64; 9]; nx]; ny];
        for y in 0..ny {
            for x in 0..nx {
                for q in 0..9 {
                    let xn = (x as i32 + EX9[q]).rem_euclid(nx as i32) as usize;
                    let yn = (y as i32 + EY9[q]).rem_euclid(ny as i32) as usize;
                    f_new[yn][xn][q] = self.f[y][x][q];
                }
            }
        }
        self.f = f_new;
    }
    /// Update macroscopic variables from distribution functions.
    pub fn update_macros(&mut self) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for q in 0..9 {
                    rho += self.f[y][x][q];
                    ux += EX9[q] as f64 * self.f[y][x][q];
                    uy += EY9[q] as f64 * self.f[y][x][q];
                }
                self.rho[y][x] = rho;
                if rho > 1e-10 {
                    self.u[y][x] = [ux / rho, uy / rho];
                }
            }
        }
    }
    /// Perform one full LBM time step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }
    /// Set uniform solid volume fraction.
    pub fn set_uniform_phi(&mut self, phi: f64) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                self.phi[y][x] = phi;
            }
        }
    }
    /// Get average vertical velocity.
    pub fn avg_vy(&self) -> f64 {
        let mut sum = 0.0;
        let n = self.nx * self.ny;
        for y in 0..self.ny {
            for x in 0..self.nx {
                sum += self.u[y][x][1];
            }
        }
        sum / n as f64
    }
}
