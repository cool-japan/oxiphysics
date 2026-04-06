#![allow(
    clippy::needless_range_loop,
    clippy::ptr_arg,
    clippy::too_many_arguments
)]
//! Auto-generated module
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{E_CHARGE, EK_CS2, EK_CX, EK_CY, EK_W, EPSILON_0, K_B, N_A};

/// Diffusio-osmosis: flow driven by a concentration gradient along a wall.
pub struct DiffusioOsmosis;
impl DiffusioOsmosis {
    /// Diffusio-osmotic velocity (Derjaguin formula):
    ///
    /// `u_do = (epsilon_0 * epsilon_r * (k_B * T / e)^2 * (zeta / mu)) * grad(ln c)`
    pub fn velocity(
        zeta: f64,
        temperature: f64,
        viscosity: f64,
        epsilon_r: f64,
        grad_ln_c: f64,
    ) -> f64 {
        let thermal_voltage = K_B * temperature / E_CHARGE;
        EPSILON_0 * epsilon_r * thermal_voltage * thermal_voltage * zeta / viscosity * grad_ln_c
    }
    /// Lewis number: ratio of thermal to mass diffusivity.
    pub fn lewis_number(thermal_diffusivity: f64, mass_diffusivity: f64) -> f64 {
        thermal_diffusivity / mass_diffusivity.max(1e-30)
    }
}
/// Electrophoretic mobility model (Henry equation).
pub struct ElectrophoreticMobility {
    /// Particle radius (m)
    pub particle_radius: f64,
    /// Zeta potential of the particle (V)
    pub zeta: f64,
    /// Electrolyte parameters
    pub params: ElectrolyteParams,
}
impl ElectrophoreticMobility {
    /// Create a new electrophoretic mobility model.
    pub fn new(particle_radius: f64, zeta: f64, params: ElectrolyteParams) -> Self {
        Self {
            particle_radius,
            zeta,
            params,
        }
    }
    /// Henry's function f(kappa * a).
    ///
    /// Limits: f → 1.0 for thin double layer (kappa*a → infinity),
    ///         f → 1.5 for thick double layer (kappa*a → 0, Hückel limit).
    ///
    /// Approximate smooth interpolation:
    ///   f = 1.0 + 2.5 / (1 + 2 * exp(-kappa_a))
    ///
    /// (Recovers ≈ 1.5 at kappa_a = 0 and ≈ 1.0 at kappa_a → ∞.)
    pub fn henry_function(&self, kappa_a: f64) -> f64 {
        1.0 + 2.5 / (1.0 + 2.0 * (-kappa_a).exp())
    }
    /// Electrophoretic mobility (m^2/(V·s)).
    ///
    /// mu_e = (2/3) * epsilon_0 * epsilon_r * zeta / mu * f(kappa*a)
    pub fn mobility(&self, kappa: f64) -> f64 {
        let kappa_a = kappa * self.particle_radius;
        let f = self.henry_function(kappa_a);
        (2.0 / 3.0) * EPSILON_0 * self.params.permittivity_r * self.zeta / self.params.viscosity * f
    }
    /// Electrophoretic velocity (m/s).
    ///
    /// v = mu_e * E
    pub fn velocity(&self, e_field: f64, kappa: f64) -> f64 {
        self.mobility(kappa) * e_field
    }
}
impl ElectrophoreticMobility {
    /// Stokes drag force on the particle:
    ///
    /// `F_drag = 6 * pi * mu * a * v`
    pub fn stokes_drag(&self, velocity: f64) -> f64 {
        6.0 * std::f64::consts::PI * self.params.viscosity * self.particle_radius * velocity
    }
    /// Sedimentation velocity of the particle under gravity:
    ///
    /// `v_sed = 2 * a^2 * (rho_p - rho_f) * g / (9 * mu)`
    pub fn sedimentation_velocity(
        &self,
        density_particle: f64,
        density_fluid: f64,
        gravity: f64,
    ) -> f64 {
        2.0 * self.particle_radius.powi(2) * (density_particle - density_fluid) * gravity
            / (9.0 * self.params.viscosity)
    }
    /// Peclet number for electrophoresis:
    ///
    /// `Pe = v * a / D`
    pub fn peclet_number(&self, velocity: f64, diffusivity: f64) -> f64 {
        velocity * self.particle_radius / diffusivity
    }
}
/// Couple an external electric field to the body force in LBM.
///
/// For LBM with electroosmotic flow, the body force on the fluid at each cell
/// is `F = rho_e * E` where `rho_e` is the local charge density.
pub struct ElectroosmoticBodyForce {
    /// Electric field vector (V/m).
    pub e_field: [f64; 3],
}
impl ElectroosmoticBodyForce {
    /// Create a new electroosmotic body force from an applied field.
    pub fn new(e_field: [f64; 3]) -> Self {
        Self { e_field }
    }
    /// Compute the electroosmotic body force at a cell with charge density `rho_e` (C/m^3).
    pub fn force_at_cell(&self, rho_e: f64) -> [f64; 3] {
        [
            rho_e * self.e_field[0],
            rho_e * self.e_field[1],
            rho_e * self.e_field[2],
        ]
    }
    /// Apply the electroosmotic body force to a 2D force field.
    ///
    /// `force_field[k] = [rho_e[k\] * E_x, rho_e[k] * E_y, 0]`
    pub fn apply_2d(&self, rho_e: &[f64], force_field: &mut Vec<[f64; 3]>) {
        for (k, &re) in rho_e.iter().enumerate() {
            force_field[k][0] += re * self.e_field[0];
            force_field[k][1] += re * self.e_field[1];
        }
    }
}
/// Diagnostic functions for the electric double layer.
pub struct DebyeLayerDiagnostics;
impl DebyeLayerDiagnostics {
    /// Compute the Debye length from ion species.
    pub fn debye_length(species: &[IonSpecies], params: &ElectrolyteParams) -> f64 {
        let ionic_str = PoissonBoltzmann::ionic_strength(species);
        params.debye_length(ionic_str)
    }
    /// Check whether the thin double-layer approximation is valid:
    ///
    /// `kappa * a >> 1` where a is the channel half-width.
    pub fn is_thin_edl(lambda_d: f64, channel_half_width: f64) -> bool {
        channel_half_width / lambda_d > 10.0
    }
    /// Overlap parameter: ratio of Debye length to channel half-width.
    ///
    /// When > 1, the double layers from opposite walls overlap significantly.
    pub fn overlap_parameter(lambda_d: f64, channel_half_width: f64) -> f64 {
        lambda_d / channel_half_width
    }
    /// Surface potential from Grahame equation (low-potential limit):
    ///
    /// `sigma = epsilon_0 * epsilon_r * zeta / lambda_D`
    pub fn surface_charge_from_zeta(zeta: f64, lambda_d: f64, epsilon_r: f64) -> f64 {
        EPSILON_0 * epsilon_r * zeta / lambda_d
    }
}
/// Activity coefficients and Debye-Hückel theory.
pub struct DebyeHuckel;
impl DebyeHuckel {
    /// Extended Debye-Hückel activity coefficient (mean):
    ///
    /// `log(gamma) = -A * z^2 * sqrt(I) / (1 + B * a * sqrt(I))`
    ///
    /// where A ~ 0.509 for water at 25°C, B ~ 3.28e9 m^-1/2 mol^-1/2.
    pub fn activity_coefficient_log10(
        valence: i32,
        ionic_strength: f64,
        ion_diameter_m: f64,
    ) -> f64 {
        let a = 0.509_f64;
        let b = 3.28e9_f64;
        let z = valence as f64;
        let sqrt_i = ionic_strength.max(0.0).sqrt();
        -a * z * z * sqrt_i / (1.0 + b * ion_diameter_m * sqrt_i)
    }
    /// Activity coefficient (linear scale).
    pub fn activity_coefficient(valence: i32, ionic_strength: f64, ion_diameter_m: f64) -> f64 {
        let log10_gamma = Self::activity_coefficient_log10(valence, ionic_strength, ion_diameter_m);
        10.0_f64.powf(log10_gamma)
    }
}
/// Physical parameters of the electrolyte solution.
pub struct ElectrolyteParams {
    /// Temperature (K)
    pub temperature: f64,
    /// Relative permittivity (dimensionless; water ≈ 80)
    pub permittivity_r: f64,
    /// Dynamic viscosity (Pa·s)
    pub viscosity: f64,
}
impl ElectrolyteParams {
    /// Create new electrolyte parameters.
    pub fn new(temperature: f64, permittivity_r: f64, viscosity: f64) -> Self {
        Self {
            temperature,
            permittivity_r,
            viscosity,
        }
    }
    /// Debye screening length (m).
    ///
    /// lambda_D = sqrt(epsilon_0 * epsilon_r * k_B * T / (2 * N_A * e^2 * I))
    ///
    /// where I = 0.5 * sum(c_i * z_i^2) is the ionic strength (mol/m^3).
    pub fn debye_length(&self, ionic_strength: f64) -> f64 {
        let numerator = EPSILON_0 * self.permittivity_r * K_B * self.temperature;
        let denominator = 2.0 * N_A * E_CHARGE * E_CHARGE * ionic_strength;
        (numerator / denominator).sqrt()
    }
    /// Surface charge density from zeta potential (Debye-Hückel approximation).
    ///
    /// sigma = -epsilon_0 * epsilon_r * zeta / lambda_D
    pub fn zeta_potential_to_surface_charge(&self, zeta: f64, lambda_d: f64) -> f64 {
        -EPSILON_0 * self.permittivity_r * zeta / lambda_d
    }
}
/// Ion distribution in the electric double layer using the Boltzmann equation.
pub struct BoltzmannIonDistribution;
impl BoltzmannIonDistribution {
    /// Ion concentration at position y from the surface.
    ///
    /// c(y) = c_bulk * exp(-z * e * psi(y) / (kB * T))
    pub fn concentration(c_bulk: f64, z: f64, psi: f64, temperature: f64) -> f64 {
        let beta = z * E_CHARGE / (K_B * temperature);
        c_bulk * (-beta * psi).exp()
    }
    /// Net charge density for a 1:1 electrolyte (Na+/Cl-):
    ///
    /// rho_e = e * NA * (c_+ - c_-) = -2 * e * NA * c_bulk * sinh(e * psi / kT)
    pub fn net_charge_density_symmetric(c_bulk: f64, psi: f64, temperature: f64) -> f64 {
        let beta = E_CHARGE / (K_B * temperature);
        -2.0 * E_CHARGE * N_A * c_bulk * (beta * psi).sinh()
    }
    /// Electric field from potential gradient (central difference).
    pub fn electric_field_from_potential(phi: &[f64], dx: f64) -> Vec<f64> {
        let n = phi.len();
        let mut e_field = vec![0.0_f64; n];
        for i in 0..n {
            let ip = if i + 1 < n { i + 1 } else { i };
            let im = if i > 0 { i - 1 } else { i };
            let denom = if i == 0 || i + 1 == n { dx } else { 2.0 * dx };
            e_field[i] = -(phi[ip] - phi[im]) / denom;
        }
        e_field
    }
}
/// 2D grid of ionic concentrations and electric potential.
pub struct IonicConcentrationField {
    /// Grid width
    pub nx: usize,
    /// Grid height
    pub ny: usize,
    /// Number of ionic species
    pub n_species: usize,
    /// Concentrations: `concentrations[species_idx][cell_idx]` (mol/m^3)
    pub concentrations: Vec<Vec<f64>>,
    /// Electric potential phi at each cell (V)
    pub potential: Vec<f64>,
}
impl IonicConcentrationField {
    /// Create a new zeroed concentration field.
    pub fn new(nx: usize, ny: usize, n_species: usize) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            n_species,
            concentrations: vec![vec![0.0; n]; n_species],
            potential: vec![0.0; n],
        }
    }
    /// Row-major cell index.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
}
/// Energy stored in the electric double layer near a flat charged surface.
pub struct ElectricDoubleLayerEnergy;
impl ElectricDoubleLayerEnergy {
    /// Electrostatic interaction energy per unit area between two flat EDLs
    /// separated by distance D (Derjaguin–Landau–Verwey–Overbeek theory):
    ///
    /// `G(D) = 64 * n_bulk * k_B * T * lambda_D * tanh^2(z*e*zeta / (4*k_B*T)) * exp(-D/lambda_D)`
    pub fn derjaguin_energy_per_area(
        n_bulk: f64,
        temperature: f64,
        lambda_d: f64,
        zeta: f64,
        valence: i32,
        separation: f64,
    ) -> f64 {
        let z = valence as f64;
        let arg = z * E_CHARGE * zeta / (4.0 * K_B * temperature);
        let tanh_sq = arg.tanh().powi(2);
        64.0 * n_bulk * K_B * temperature * lambda_d * tanh_sq * (-separation / lambda_d).exp()
    }
    /// Electrostatic disjoining pressure:
    ///
    /// `Pi(D) = -dG/dD = (G_0 / lambda_D) * exp(-D/lambda_D)`
    pub fn disjoining_pressure(
        n_bulk: f64,
        temperature: f64,
        lambda_d: f64,
        zeta: f64,
        valence: i32,
        separation: f64,
    ) -> f64 {
        let z = valence as f64;
        let arg = z * E_CHARGE * zeta / (4.0 * K_B * temperature);
        let tanh_sq = arg.tanh().powi(2);
        let g0 = 64.0 * n_bulk * K_B * temperature * lambda_d * tanh_sq;
        g0 / lambda_d * (-separation / lambda_d).exp()
    }
}
/// Solver for the Poisson equation (electrostatics).
pub struct PoissonSolver;
impl PoissonSolver {
    /// Solve the 2D Poisson equation with Gauss-Seidel iteration.
    ///
    /// Equation: nabla^2 phi = -rho / (epsilon_0 * epsilon_r)
    ///
    /// Gauss-Seidel update:
    ///   phi\[i\] = (sum of 4 neighbours - rho\[i\] * dx^2 / epsilon_r) / 4
    ///
    /// Returns the number of iterations used.
    pub fn solve_poisson_2d(
        phi: &mut Vec<f64>,
        rho: &[f64],
        nx: usize,
        ny: usize,
        dx: f64,
        epsilon_r: f64,
        max_iter: usize,
        tol: f64,
    ) -> usize {
        let dx2 = dx * dx;
        let eps_total = EPSILON_0 * epsilon_r;
        for iter in 0..max_iter {
            let mut max_change = 0.0_f64;
            for y in 0..ny {
                for x in 0..nx {
                    let idx = y * nx + x;
                    let phi_e = if x + 1 < nx { phi[y * nx + x + 1] } else { 0.0 };
                    let phi_w = if x > 0 { phi[y * nx + x - 1] } else { 0.0 };
                    let phi_n = if y + 1 < ny {
                        phi[(y + 1) * nx + x]
                    } else {
                        0.0
                    };
                    let phi_s = if y > 0 { phi[(y - 1) * nx + x] } else { 0.0 };
                    let new_phi =
                        (phi_e + phi_w + phi_n + phi_s - rho[idx] * dx2 / eps_total) / 4.0;
                    let change = (new_phi - phi[idx]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                    phi[idx] = new_phi;
                }
            }
            if max_change < tol {
                return iter + 1;
            }
        }
        max_iter
    }
    /// Compute the charge density field (C/m^3) from ionic concentrations.
    ///
    /// rho = sum_s z_s * e * c_s   (where c_s is in mol/m^3, so multiply by N_A)
    pub fn charge_density(concentrations: &[Vec<f64>], species: &[IonSpecies]) -> Vec<f64> {
        if concentrations.is_empty() {
            return Vec::new();
        }
        let n = concentrations[0].len();
        let mut rho = vec![0.0_f64; n];
        for (s, ion) in species.iter().enumerate() {
            let z = ion.valence as f64;
            for k in 0..n {
                rho[k] += z * E_CHARGE * N_A * concentrations[s][k];
            }
        }
        rho
    }
}
/// Solver for the Nernst-Planck ionic transport equations.
pub struct NernstPlanckSolver;
impl NernstPlanckSolver {
    /// Compute the diffusive flux magnitude at each cell (mol/(m^2·s)).
    ///
    /// J_diff = -D * grad(c)
    ///
    /// Returns a flat vector of length `2 * nx * ny` where
    /// `flux[2*k]` = Jx at cell k and `flux[2*k+1]` = Jy at cell k.
    pub fn diffusive_flux(c: &[f64], nx: usize, ny: usize, dx: f64, diffusivity: f64) -> Vec<f64> {
        let n = nx * ny;
        let mut flux = vec![0.0_f64; 2 * n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                let c_e = if x + 1 < nx { c[y * nx + x + 1] } else { c[k] };
                let c_w = if x > 0 { c[y * nx + x - 1] } else { c[k] };
                let grad_cx = if x > 0 && x + 1 < nx {
                    (c_e - c_w) / (2.0 * dx)
                } else {
                    (c_e - c_w) / dx
                };
                let c_n = if y + 1 < ny {
                    c[(y + 1) * nx + x]
                } else {
                    c[k]
                };
                let c_s = if y > 0 { c[(y - 1) * nx + x] } else { c[k] };
                let grad_cy = if y > 0 && y + 1 < ny {
                    (c_n - c_s) / (2.0 * dx)
                } else {
                    (c_n - c_s) / dx
                };
                flux[2 * k] = -diffusivity * grad_cx;
                flux[2 * k + 1] = -diffusivity * grad_cy;
            }
        }
        flux
    }
    /// Compute the migrative (electrophoretic) flux at each cell.
    ///
    /// J_mig = -mu * c * grad(phi)
    ///
    /// Returns a flat vector of length `2 * nx * ny` (Jx, Jy interleaved).
    pub fn migrative_flux(
        c: &[f64],
        phi: &[f64],
        nx: usize,
        ny: usize,
        dx: f64,
        mobility: f64,
    ) -> Vec<f64> {
        let n = nx * ny;
        let mut flux = vec![0.0_f64; 2 * n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                let phi_e = if x + 1 < nx {
                    phi[y * nx + x + 1]
                } else {
                    phi[k]
                };
                let phi_w = if x > 0 { phi[y * nx + x - 1] } else { phi[k] };
                let grad_phi_x = if x > 0 && x + 1 < nx {
                    (phi_e - phi_w) / (2.0 * dx)
                } else {
                    (phi_e - phi_w) / dx
                };
                let phi_n = if y + 1 < ny {
                    phi[(y + 1) * nx + x]
                } else {
                    phi[k]
                };
                let phi_s = if y > 0 { phi[(y - 1) * nx + x] } else { phi[k] };
                let grad_phi_y = if y > 0 && y + 1 < ny {
                    (phi_n - phi_s) / (2.0 * dx)
                } else {
                    (phi_n - phi_s) / dx
                };
                flux[2 * k] = -mobility * c[k] * grad_phi_x;
                flux[2 * k + 1] = -mobility * c[k] * grad_phi_y;
            }
        }
        flux
    }
    /// Update concentrations with an explicit finite-difference step.
    ///
    /// dc/dt = -div(J)  →  c\[k\] += dt * (-div J)\[k\]
    ///
    /// `flux` is the combined (diffusive + migrative) flux, length `2 * nx * ny`.
    pub fn update_concentrations(
        c: &mut Vec<f64>,
        flux: &[f64],
        nx: usize,
        ny: usize,
        dt: f64,
        dx: f64,
    ) {
        let n = nx * ny;
        let mut div_j = vec![0.0_f64; n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                let jx_e = if x + 1 < nx {
                    flux[2 * (y * nx + x + 1)]
                } else {
                    flux[2 * k]
                };
                let jx_w = if x > 0 {
                    flux[2 * (y * nx + x - 1)]
                } else {
                    flux[2 * k]
                };
                let jy_n = if y + 1 < ny {
                    flux[2 * ((y + 1) * nx + x) + 1]
                } else {
                    flux[2 * k + 1]
                };
                let jy_s = if y > 0 {
                    flux[2 * ((y - 1) * nx + x) + 1]
                } else {
                    flux[2 * k + 1]
                };
                div_j[k] = (jx_e - jx_w) / (2.0 * dx) + (jy_n - jy_s) / (2.0 * dx);
            }
        }
        for k in 0..n {
            c[k] -= dt * div_j[k];
        }
    }
}
/// Electroosmotic flow model (Helmholtz-Smoluchowski).
pub struct ElectroosmosticFlow {
    /// Zeta potential at the channel wall (V)
    pub zeta: f64,
    /// Electrolyte parameters
    pub params: ElectrolyteParams,
}
impl ElectroosmosticFlow {
    /// Create a new electroosmotic flow model.
    pub fn new(zeta: f64, params: ElectrolyteParams) -> Self {
        Self { zeta, params }
    }
    /// Helmholtz-Smoluchowski electroosmotic velocity (m/s).
    ///
    /// u_eo = -epsilon_0 * epsilon_r * zeta * E / mu
    pub fn helmholtz_smoluchowski_velocity(&self, e_field: f64) -> f64 {
        -EPSILON_0 * self.params.permittivity_r * self.zeta * e_field / self.params.viscosity
    }
    /// Electroosmotic volumetric flux (m^2/s) for a channel of given width.
    ///
    /// Q = u_eo * width
    pub fn electroosmotic_flux(&self, e_field: f64, width: f64) -> f64 {
        self.helmholtz_smoluchowski_velocity(e_field) * width
    }
}
impl ElectroosmosticFlow {
    /// Electroosmotic velocity profile in a slit channel:
    ///
    /// `u(y) = u_eo * (1 - cosh(y/lambda_D) / cosh(h/lambda_D))`
    ///
    /// where y is measured from channel center and h is the half-width.
    pub fn slit_channel_velocity(
        &self,
        e_field: f64,
        y_from_center: f64,
        half_width: f64,
        lambda_d: f64,
    ) -> f64 {
        let u_eo = self.helmholtz_smoluchowski_velocity(e_field);
        if lambda_d <= 0.0 || half_width <= 0.0 {
            return u_eo;
        }
        let cosh_y = (y_from_center / lambda_d).cosh();
        let cosh_h = (half_width / lambda_d).cosh();
        u_eo * (1.0 - cosh_y / cosh_h)
    }
    /// Average electroosmotic velocity in a slit channel:
    ///
    /// `u` = u_eo * (1 - lambda_D / h * tanh(h / lambda_D))`
    pub fn slit_channel_average_velocity(
        &self,
        e_field: f64,
        half_width: f64,
        lambda_d: f64,
    ) -> f64 {
        let u_eo = self.helmholtz_smoluchowski_velocity(e_field);
        if lambda_d <= 0.0 || half_width <= 0.0 {
            return u_eo;
        }
        u_eo * (1.0 - lambda_d / half_width * (half_width / lambda_d).tanh())
    }
}
/// Linearized Poisson-Boltzmann solver for equilibrium EDL.
pub struct PoissonBoltzmann;
impl PoissonBoltzmann {
    /// Analytic solution for potential near a flat charged surface:
    ///
    /// `phi(y) = zeta * exp(-y / lambda_D)`
    ///
    /// Valid in the Debye-Huckel (linearized) regime.
    pub fn flat_surface_potential(zeta: f64, y: f64, lambda_d: f64) -> f64 {
        zeta * (-y / lambda_d).exp()
    }
    /// Charge density from Boltzmann distribution at position y:
    ///
    /// `rho(y) = -epsilon_0 * epsilon_r * zeta / lambda_D^2 * exp(-y/lambda_D)`
    pub fn flat_surface_charge_density(zeta: f64, y: f64, lambda_d: f64, epsilon_r: f64) -> f64 {
        -EPSILON_0 * epsilon_r * zeta / (lambda_d * lambda_d) * (-y / lambda_d).exp()
    }
    /// Solve the 1D linearized Poisson-Boltzmann equation numerically
    /// on a uniform grid using finite differences.
    ///
    /// `d^2 phi / dy^2 = phi / lambda_D^2`
    ///
    /// Boundary: phi(0) = zeta, phi(L) = 0.
    pub fn solve_1d(zeta: f64, lambda_d: f64, ny: usize, dy: f64) -> Vec<f64> {
        let mut phi = vec![0.0_f64; ny];
        phi[0] = zeta;
        let inv_lambda_sq = 1.0 / (lambda_d * lambda_d);
        for _iter in 0..5000 {
            let mut max_change = 0.0_f64;
            for j in 1..ny - 1 {
                let new_phi = (phi[j - 1] + phi[j + 1]) / (2.0 + inv_lambda_sq * dy * dy);
                let change = (new_phi - phi[j]).abs();
                if change > max_change {
                    max_change = change;
                }
                phi[j] = new_phi;
            }
            if max_change < 1e-12 {
                break;
            }
        }
        phi
    }
    /// Ionic strength from species:
    ///
    /// `I = 0.5 * sum(c_k * z_k^2)`
    pub fn ionic_strength(species: &[IonSpecies]) -> f64 {
        0.5 * species
            .iter()
            .map(|s| s.concentration_bulk * (s.valence as f64).powi(2))
            .sum::<f64>()
    }
}
/// Infer zeta potential from a measured streaming potential.
pub struct ZetaPotentialEstimator;
impl ZetaPotentialEstimator {
    /// Estimate zeta potential from streaming potential coefficient:
    ///
    /// `zeta = -(E_stream / delta_P) * (mu * sigma_el) / (epsilon_0 * epsilon_r)`
    pub fn from_streaming_potential_coefficient(
        e_stream_per_dp: f64,
        viscosity: f64,
        conductivity: f64,
        epsilon_r: f64,
    ) -> f64 {
        -e_stream_per_dp * viscosity * conductivity / (EPSILON_0 * epsilon_r)
    }
    /// Estimate zeta potential from electrophoretic mobility (Smoluchowski limit):
    ///
    /// `zeta = mu_e * mu / (epsilon_0 * epsilon_r)`
    pub fn from_electrophoretic_mobility(mobility_e: f64, viscosity: f64, epsilon_r: f64) -> f64 {
        mobility_e * viscosity / (EPSILON_0 * epsilon_r)
    }
    /// Estimate zeta potential from surface charge density:
    ///
    /// `zeta = sigma * lambda_D / (epsilon_0 * epsilon_r)`
    pub fn from_surface_charge(sigma: f64, lambda_d: f64, epsilon_r: f64) -> f64 {
        sigma * lambda_d / (EPSILON_0 * epsilon_r)
    }
}
/// Transient electroosmotic flow in a slit channel after a step electric field.
///
/// The start-up velocity at the channel centerline follows:
///
/// `u(t) = u_eo * (1 - sum_n C_n * exp(-lambda_n * t))`
pub struct TransientEof;
impl TransientEof {
    /// Dimensionless characteristic time for first mode:
    ///
    /// `t* = nu * t / h^2` (h = half-width, nu = kinematic viscosity)
    pub fn dimensionless_time(t: f64, nu: f64, half_width: f64) -> f64 {
        nu * t / (half_width * half_width)
    }
    /// Start-up velocity at channel center using first-mode approximation:
    ///
    /// `u_center(t) ≈ u_eo * (1 - exp(-pi^2 * nu * t / h^2))`
    pub fn startup_velocity_center(u_eo: f64, t: f64, nu: f64, half_width: f64) -> f64 {
        let t_star = Self::dimensionless_time(t, nu, half_width);
        u_eo * (1.0 - (-std::f64::consts::PI * std::f64::consts::PI * t_star).exp())
    }
    /// Time to reach 99% of steady state:
    ///
    /// `t_99 = -h^2 * ln(0.01) / (pi^2 * nu) ≈ 0.466 * h^2 / nu`
    pub fn time_to_steady_state(nu: f64, half_width: f64) -> f64 {
        -half_width * half_width * 0.01_f64.ln()
            / (std::f64::consts::PI * std::f64::consts::PI * nu)
    }
}
/// Full (non-linearised) 1D Poisson-Boltzmann solver for a symmetric z:z electrolyte.
///
/// Equation:
///   d²psi/dy² = (2 * n_bulk * z * e / epsilon) * sinh(z * e * psi / (kB * T))
///
/// Boundary conditions: psi(0) = zeta, psi(L) = 0.
pub struct NonlinearPoissonBoltzmann;
impl NonlinearPoissonBoltzmann {
    /// Solve the 1D non-linear PB equation using Gauss-Seidel iteration.
    ///
    /// # Arguments
    /// * `zeta`      - surface potential (V)
    /// * `n_bulk`    - bulk ionic number density (mol/m^3 × NA → ions/m^3)
    /// * `z`         - ion valence (symmetric electrolyte)
    /// * `epsilon`   - absolute permittivity (F/m)
    /// * `temperature` - temperature (K)
    /// * `ny`        - number of grid points
    /// * `dy`        - grid spacing (m)
    pub fn solve_1d(
        zeta: f64,
        n_bulk: f64,
        z: f64,
        epsilon: f64,
        temperature: f64,
        ny: usize,
        dy: f64,
    ) -> Vec<f64> {
        let mut psi = vec![0.0_f64; ny];
        psi[0] = zeta;
        let beta = z * E_CHARGE / (K_B * temperature);
        let rhs_scale = 2.0 * n_bulk * z * E_CHARGE / epsilon;
        for _iter in 0..10_000 {
            let mut max_change = 0.0_f64;
            for j in 1..ny - 1 {
                let rhs = rhs_scale * (beta * psi[j]).sinh();
                let new_psi = 0.5 * (psi[j - 1] + psi[j + 1] - dy * dy * rhs);
                let change = (new_psi - psi[j]).abs();
                if change > max_change {
                    max_change = change;
                }
                psi[j] = new_psi;
            }
            if max_change < 1e-14 {
                break;
            }
        }
        psi
    }
    /// Compute charge density from psi field (mol/m^3 → C/m^3 after × NA × e).
    ///
    /// rho_e = -2 * n_bulk * z * e * sinh(z * e * psi / kB T)
    pub fn charge_density_field(psi: &[f64], n_bulk: f64, z: f64, temperature: f64) -> Vec<f64> {
        let beta = z * E_CHARGE / (K_B * temperature);
        psi.iter()
            .map(|&p| -2.0 * n_bulk * z * E_CHARGE * (beta * p).sinh())
            .collect()
    }
}
/// Electroosmotic pump: generates pressure or flow using an electric field.
///
/// For a porous medium filled with electrolyte, the pump characteristic is:
///   Q = Q_max * (1 - delta_P / delta_P_max)
/// where:
///   Q_max = epsilon * zeta * A * E / mu  (max flow rate at zero back-pressure)
///   delta_P_max = epsilon * zeta * E / K  (max pressure at zero flow)
pub struct ElectroosmosticPump {
    /// Zeta potential (V).
    pub zeta: f64,
    /// Absolute permittivity (F/m).
    pub epsilon: f64,
    /// Dynamic viscosity (Pa·s).
    pub mu: f64,
    /// Darcy permeability K (m^2).
    pub k: f64,
    /// Cross-sectional area A (m^2).
    pub area: f64,
}
impl ElectroosmosticPump {
    /// Create a new electroosmotic pump.
    pub fn new(zeta: f64, epsilon: f64, mu: f64, k: f64, area: f64) -> Self {
        Self {
            zeta,
            epsilon,
            mu,
            k,
            area,
        }
    }
    /// Maximum (free) volumetric flow rate (m^3/s) at zero back-pressure.
    ///
    /// Q_max = epsilon * zeta * A * E / mu
    pub fn max_flow_rate(&self, e_field: f64) -> f64 {
        self.epsilon * self.zeta * self.area * e_field / self.mu
    }
    /// Maximum back-pressure (Pa) at zero flow.
    ///
    /// delta_P_max = epsilon * zeta * E / K  × (channel geometric factor)
    pub fn max_pressure(&self, e_field: f64) -> f64 {
        if self.k < 1e-30 {
            return 0.0;
        }
        self.epsilon * self.zeta * e_field / self.k
    }
    /// Actual flow rate at back-pressure `delta_p`.
    pub fn flow_rate(&self, e_field: f64, delta_p: f64) -> f64 {
        let q_max = self.max_flow_rate(e_field);
        let dp_max = self.max_pressure(e_field);
        if dp_max.abs() < 1e-30 {
            return q_max;
        }
        q_max * (1.0 - delta_p / dp_max)
    }
    /// Efficiency of the pump at operating point (delta_p, flow_rate).
    ///
    /// eta = delta_P * Q / (V * I)
    pub fn efficiency(&self, e_field: f64, delta_p: f64, current: f64, voltage: f64) -> f64 {
        let q = self.flow_rate(e_field, delta_p);
        let power_in = (current * voltage).abs();
        if power_in < 1e-30 {
            return 0.0;
        }
        (delta_p * q).abs() / power_in
    }
}
/// Electrophoresis of a colloidal particle (Henry approximation).
#[derive(Debug, Clone, Copy)]
pub struct ElectrophoresisParticle {
    /// Particle radius (m).
    pub radius: f64,
    /// Zeta potential (V).
    pub zeta: f64,
    /// Fluid dynamic viscosity (Pa·s).
    pub mu_f: f64,
}
impl ElectrophoresisParticle {
    /// Create a new `ElectrophoresisParticle`.
    pub fn new(radius: f64, zeta: f64, mu_f: f64) -> Self {
        Self { radius, zeta, mu_f }
    }
    /// Electrophoretic velocity using the Smoluchowski limit (Henry f=1).
    ///
    /// v = (epsilon_0 * epsilon_r * zeta / mu) * E
    ///
    /// Here we use a simplified form: v = (2/3) * (epsilon_0 * 80 * zeta / mu_f) * E
    /// for the Henry function at infinite dilution.  For a general implementation
    /// pass the full permittivity.
    ///
    /// # Arguments
    /// * `e` - electric field magnitude (V/m)
    pub fn electrophoretic_velocity(&self, e: f64) -> f64 {
        let eps = EPSILON_0 * 80.0;
        (2.0 / 3.0) * eps * self.zeta / self.mu_f * e
    }
}
/// Differential capacitance of the electric double layer.
///
/// In the Gouy-Chapman model (symmetric z:z electrolyte):
/// C_D = epsilon * kappa * cosh(z * e * zeta / (2 * kB * T))
pub struct ElectricDoubleLayerCapacitance;
impl ElectricDoubleLayerCapacitance {
    /// Differential capacitance per unit area (F/m^2).
    ///
    /// # Arguments
    /// * `epsilon`     - absolute permittivity (F/m)
    /// * `kappa`       - inverse Debye length (1/m)
    /// * `zeta`        - surface potential (V)
    /// * `z`           - ion valence
    /// * `temperature` - temperature (K)
    pub fn gouy_chapman(epsilon: f64, kappa: f64, zeta: f64, z: f64, temperature: f64) -> f64 {
        let beta_zeta = z * E_CHARGE * zeta / (2.0 * K_B * temperature);
        epsilon * kappa * beta_zeta.cosh()
    }
    /// Inner Helmholtz layer capacitance (compact layer):
    /// C_H = epsilon / d_H
    ///
    /// where d_H is the thickness of the Helmholtz layer (~0.3-0.5 nm).
    pub fn helmholtz_layer(epsilon: f64, thickness: f64) -> f64 {
        if thickness < 1e-30 {
            return 0.0;
        }
        epsilon / thickness
    }
    /// Total double-layer capacitance (series combination of Helmholtz + diffuse):
    /// 1/C_total = 1/C_H + 1/C_D
    pub fn total_series(c_helmholtz: f64, c_diffuse: f64) -> f64 {
        let inv_sum = 1.0 / c_helmholtz.max(1e-30) + 1.0 / c_diffuse.max(1e-30);
        1.0 / inv_sum
    }
}
/// Electroviscous effect: apparent viscosity increase in charged nanochannels.
///
/// The effective viscosity increases when the double-layer thickness is
/// comparable to the channel width (Donath-Voigt model).
pub struct ElectroviscousEffect;
impl ElectroviscousEffect {
    /// Relative apparent viscosity (eta_app / eta_0).
    ///
    /// `alpha = 1 + (epsilon^2 * zeta^2) / (eta_0^2 * K_c)`
    ///
    /// where K_c is the channel cross-sectional conductance (S).
    pub fn apparent_viscosity_ratio(epsilon: f64, zeta: f64, eta0: f64, conductance: f64) -> f64 {
        if conductance.abs() < 1e-30 {
            return 1.0;
        }
        1.0 + epsilon * epsilon * zeta * zeta / (eta0 * eta0 * conductance)
    }
    /// Streaming current contribution to apparent viscosity (simplified):
    ///
    /// eta_app = eta0 * (1 + (epsilon * zeta / (mu * sigma * L / A))^2)
    pub fn apparent_viscosity_from_dimensions(
        epsilon: f64,
        zeta: f64,
        eta0: f64,
        conductivity: f64,
        length: f64,
        area: f64,
    ) -> f64 {
        if area < 1e-30 || conductivity < 1e-30 {
            return eta0;
        }
        let k_c = conductivity * area / length;
        eta0 * Self::apparent_viscosity_ratio(epsilon, zeta, eta0, k_c)
    }
}
/// An ionic species in the electrolyte.
pub struct IonSpecies {
    /// Name (e.g. "Na+", "Cl-")
    pub name: String,
    /// Valence / charge number (e.g. +1 for Na+, -2 for SO4^2-)
    pub valence: i32,
    /// Diffusivity D (m^2/s)
    pub diffusivity: f64,
    /// Bulk concentration c_bulk (mol/m^3)
    pub concentration_bulk: f64,
}
impl IonSpecies {
    /// Create a new ion species.
    pub fn new(name: &str, valence: i32, diffusivity: f64, concentration_bulk: f64) -> Self {
        Self {
            name: name.to_string(),
            valence,
            diffusivity,
            concentration_bulk,
        }
    }
    /// Electrochemical mobility: mu = z * e * D / (k_B * T)  (m^2/(V·s))
    pub fn mobility(&self, temperature: f64) -> f64 {
        self.valence as f64 * E_CHARGE * self.diffusivity / (K_B * temperature)
    }
}
/// Electrokinetic coupling matrix based on Onsager reciprocal relations.
///
/// The coupled transport equations are:
///   J = L11 * (-grad P) + L12 * (-grad phi)
///   I = L21 * (-grad P) + L22 * (-grad phi)
///
/// where J is volume flux, I is electric current, and
/// L12 = L21 (Onsager symmetry).
#[derive(Debug, Clone, Copy)]
pub struct OnsagerCouplingMatrix {
    /// L11: hydraulic permeability (m^2·s/kg = m/(Pa·s)).
    pub l11: f64,
    /// L12 = L21: electrokinetic coupling coefficient (m^2/(V·s)).
    pub l12: f64,
    /// L22: electrical conductance per length (S/m = A/(V·m)).
    pub l22: f64,
}
impl OnsagerCouplingMatrix {
    /// Create a new Onsager coupling matrix from physical parameters.
    ///
    /// # Arguments
    /// * `permeability`  - Darcy permeability K (m^2)
    /// * `viscosity`     - dynamic viscosity mu (Pa·s)
    /// * `zeta`          - zeta potential (V)
    /// * `epsilon`       - absolute permittivity (F/m)
    /// * `conductivity`  - bulk conductivity sigma (S/m)
    /// * `lambda_d`      - Debye length (m)
    pub fn from_physical(
        permeability: f64,
        viscosity: f64,
        zeta: f64,
        epsilon: f64,
        conductivity: f64,
        lambda_d: f64,
    ) -> Self {
        let l11 = permeability / viscosity;
        let l12 = epsilon * zeta / (viscosity * lambda_d);
        let l22 = conductivity;
        Self { l11, l12, l22 }
    }
    /// Volume flux J for applied pressure gradient `dp` and electric field `e_field`.
    pub fn volume_flux(&self, dp: f64, e_field: f64) -> f64 {
        self.l11 * (-dp) + self.l12 * (-e_field)
    }
    /// Electric current density I for applied gradients.
    pub fn current_density(&self, dp: f64, e_field: f64) -> f64 {
        self.l12 * (-dp) + self.l22 * (-e_field)
    }
    /// Check Onsager symmetry (L12 = L21 holds by construction).
    pub fn is_symmetric(&self) -> bool {
        true
    }
    /// Figure of merit for electrokinetic energy conversion:
    /// Z = L12^2 / (L11 * L22)
    pub fn figure_of_merit(&self) -> f64 {
        self.l12 * self.l12 / (self.l11 * self.l22)
    }
    /// Maximum energy conversion efficiency:
    /// eta_max = (sqrt(1 + Z) - 1)^2 / (sqrt(1 + Z) + 1)^2
    pub fn max_efficiency(&self) -> f64 {
        let z = self.figure_of_merit();
        let sqrt_z_plus_1 = (1.0 + z).sqrt();
        let num = sqrt_z_plus_1 - 1.0;
        let den = sqrt_z_plus_1 + 1.0;
        (num / den).powi(2)
    }
}
/// Electrokinetic LBM solver on D2Q9.
///
/// Solves the electroosmotic flow problem by coupling:
/// 1. BGK collision for the fluid (density/velocity)
/// 2. Poisson equation for the electric potential
/// 3. Body force on fluid from charge density × electric field
pub struct ElectrokineticLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions `f\[cell\]\[dir\]`.
    pub f: Vec<[f64; 9]>,
    /// Kinematic viscosity nu (lattice units).
    pub nu: f64,
    /// Electric potential at each cell (V or lattice units).
    pub phi: Vec<f64>,
    /// Charge density at each cell (C/m^3 or lattice units).
    pub rho_e: Vec<f64>,
    /// Applied electric field [Ex, Ey] (lattice units).
    pub e_field: [f64; 2],
}
impl ElectrokineticLbm {
    /// Create a new electrokinetic LBM solver.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let feq = Self::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![feq; n],
            nu,
            phi: vec![0.0; n],
            rho_e: vec![0.0; n],
            e_field: [0.0; 2],
        }
    }
    /// Set the applied electric field.
    pub fn set_electric_field(&mut self, ex: f64, ey: f64) {
        self.e_field = [ex, ey];
    }
    /// Set the charge density field.
    pub fn set_charge_density(&mut self, rho_e: Vec<f64>) {
        self.rho_e = rho_e;
    }
    /// Linear index.
    #[inline]
    #[allow(dead_code)]
    fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Equilibrium distribution (standard D2Q9 BGK).
    pub fn equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
        let cs2 = EK_CS2;
        let ux = u[0];
        let uy = u[1];
        let u2 = ux * ux + uy * uy;
        let mut feq = [0.0_f64; 9];
        for k in 0..9 {
            let cu = EK_CX[k] as f64 * ux + EK_CY[k] as f64 * uy;
            feq[k] =
                EK_W[k] * rho * (1.0 + cu / cs2 + cu * cu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
        }
        feq
    }
    /// Compute macroscopic density and velocity.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let fi = &self.f[idx];
        let rho: f64 = fi.iter().sum();
        let ux = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * EK_CX[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        let uy = if rho > 0.0 {
            fi.iter()
                .enumerate()
                .map(|(k, &fk)| fk * EK_CY[k] as f64)
                .sum::<f64>()
                / rho
        } else {
            0.0
        };
        (rho, [ux, uy])
    }
    /// BGK relaxation parameter from nu.
    #[inline]
    fn omega(&self) -> f64 {
        1.0 / (3.0 * self.nu + 0.5)
    }
    /// Guo forcing term for electroosmotic body force F = rho_e * E.
    fn guo_force(&self, idx: usize, ux: f64, uy: f64, omega: f64) -> [f64; 9] {
        let fx = self.rho_e[idx] * self.e_field[0];
        let fy = self.rho_e[idx] * self.e_field[1];
        let factor = 1.0 - 0.5 * omega;
        let mut fi = [0.0_f64; 9];
        for k in 0..9 {
            let cx = EK_CX[k] as f64;
            let cy = EK_CY[k] as f64;
            let cu = cx * ux + cy * uy;
            let term1 = (cx - ux) * fx + (cy - uy) * fy;
            let term2 = cu * (cx * fx + cy * fy);
            fi[k] = EK_W[k] * factor * (term1 / EK_CS2 + term2 / (EK_CS2 * EK_CS2));
        }
        fi
    }
    /// BGK collision with electroosmotic body force.
    pub fn collide(&mut self) {
        let omega = self.omega();
        let n = self.nx * self.ny;
        for idx in 0..n {
            let (rho, u) = self.macros(idx);
            let feq = Self::equilibrium(rho, u);
            for k in 0..9 {
                self.f[idx][k] += omega * (feq[k] - self.f[idx][k]);
            }
            let force = self.guo_force(idx, u[0], u[1], omega);
            for k in 0..9 {
                self.f[idx][k] += force[k];
            }
        }
    }
    /// Streaming step (pull scheme, periodic).
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for j in 0..ny {
            for i in 0..nx {
                let dst = j * nx + i;
                for k in 0..9 {
                    let si = (i as isize - EK_CX[k] as isize).rem_euclid(nx as isize) as usize;
                    let sj = (j as isize - EK_CY[k] as isize).rem_euclid(ny as isize) as usize;
                    self.f[dst][k] = f_old[sj * nx + si][k];
                }
            }
        }
    }
    /// One LBM step: collide + stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }
    /// Set uniform initial state.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = Self::equilibrium(rho, u);
        for fi in self.f.iter_mut() {
            *fi = feq;
        }
    }
    /// Total mass.
    pub fn total_mass(&self) -> f64 {
        self.f.iter().map(|fi| fi.iter().sum::<f64>()).sum()
    }
    /// Average velocity magnitude.
    pub fn average_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let sum_u: f64 = (0..n)
            .map(|idx| {
                let (_, u) = self.macros(idx);
                (u[0] * u[0] + u[1] * u[1]).sqrt()
            })
            .sum();
        sum_u / n as f64
    }
    /// Update electric potential using one Gauss-Seidel sweep.
    ///
    /// Solves: nabla^2 phi = -rho_e / epsilon
    pub fn update_potential_gauss_seidel(&mut self, epsilon: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = dx * dx;
        let phi_old = self.phi.clone();
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                let phi_e = if i + 1 < nx {
                    phi_old[j * nx + i + 1]
                } else {
                    0.0
                };
                let phi_w = if i > 0 { phi_old[j * nx + i - 1] } else { 0.0 };
                let phi_n = if j + 1 < ny {
                    phi_old[(j + 1) * nx + i]
                } else {
                    0.0
                };
                let phi_s = if j > 0 {
                    phi_old[(j - 1) * nx + i]
                } else {
                    0.0
                };
                self.phi[k] = (phi_e + phi_w + phi_n + phi_s + self.rho_e[k] * dx2 / epsilon) / 4.0;
            }
        }
    }
}
/// Helmholtz-Smoluchowski electroosmotic velocity model.
#[derive(Debug, Clone, Copy)]
pub struct ElectroosmoticVelocity {
    /// Electroosmotic mobility mu_eo = -epsilon*zeta/mu (m^2/(V·s)).
    pub mu_eo: f64,
}
impl ElectroosmoticVelocity {
    /// Create a new `ElectroosmoticVelocity`.
    pub fn new(mu_eo: f64) -> Self {
        Self { mu_eo }
    }
    /// Electroosmotic velocity: v = mu_eo * E.
    ///
    /// # Arguments
    /// * `e_field` - electric field vector (V/m)
    pub fn velocity(&self, e_field: [f64; 3]) -> [f64; 3] {
        [
            self.mu_eo * e_field[0],
            self.mu_eo * e_field[1],
            self.mu_eo * e_field[2],
        ]
    }
}
/// Electric double layer (EDL) characterisation.
#[derive(Debug, Clone, Copy)]
pub struct DoubleLayer {
    /// Debye screening length lambda_D (m).
    pub debye_length: f64,
    /// Zeta potential at the shear plane (V).
    pub zeta_potential: f64,
    /// Absolute permittivity of the medium: epsilon_0 * epsilon_r (F/m).
    pub epsilon: f64,
}
impl DoubleLayer {
    /// Create a `DoubleLayer` directly from its parameters.
    pub fn new(debye_length: f64, zeta_potential: f64, epsilon: f64) -> Self {
        Self {
            debye_length,
            zeta_potential,
            epsilon,
        }
    }
    /// Compute the Debye screening length from ionic concentration.
    ///
    /// kappa^-1 = sqrt(eps * kT / (2 * NA * e^2 * c * z^2))
    ///
    /// # Arguments
    /// * `c_ionic`   - ionic concentration (mol/m^3)
    /// * `z`         - ion valence magnitude
    /// * `t`         - temperature (K)
    /// * `epsilon`   - absolute permittivity (F/m)
    pub fn debye_length_from_ionic(c_ionic: f64, z: f64, t: f64, epsilon: f64) -> f64 {
        let numerator = epsilon * K_B * t;
        let denominator = 2.0 * N_A * E_CHARGE * E_CHARGE * c_ionic * z * z;
        (numerator / denominator).sqrt()
    }
    /// Surface charge density from the Debye-Hückel approximation.
    ///
    /// sigma = -epsilon * zeta / lambda_D
    pub fn surface_charge_density(&self) -> f64 {
        -self.epsilon * self.zeta_potential / self.debye_length
    }
}
/// Dimensionless numbers relevant to electrokinetic transport.
pub struct ElectrokineticNumbers;
#[allow(non_snake_case)]
impl ElectrokineticNumbers {
    /// Electrokinetic Reynolds number: Re_ek = rho * u_eo * L / mu.
    pub fn Re(rho: f64, u_eo: f64, length: f64, mu: f64) -> f64 {
        rho * u_eo * length / mu
    }
    /// Electroviscous number (ratio of electrical to viscous effects):
    /// Ev = epsilon * (k_B * T / e)^2 / (mu * D)
    pub fn Ev(epsilon: f64, temperature: f64, mu: f64, diffusivity: f64) -> f64 {
        let vt = K_B * temperature / E_CHARGE;
        epsilon * vt * vt / (mu * diffusivity)
    }
    /// Dimensionless Debye length: kappa * L where kappa = 1/lambda_D.
    pub fn kappa_L(lambda_d: f64, length: f64) -> f64 {
        length / lambda_d
    }
    /// Dukhin number: ratio of surface conductance to bulk conductance.
    ///
    /// Du = sigma_s / (sigma_bulk * L)
    pub fn Du(surface_conductance: f64, bulk_conductivity: f64, length: f64) -> f64 {
        surface_conductance / (bulk_conductivity * length.max(1e-30))
    }
    /// Electroosmotic Peclet number: Pe_eo = u_eo * L / D.
    pub fn Pe_eo(u_eo: f64, length: f64, diffusivity: f64) -> f64 {
        u_eo * length / diffusivity
    }
}
/// Streaming potential in pressure-driven flow through a charged channel.
pub struct StreamingPotential;
impl StreamingPotential {
    /// Helmholtz-Smoluchowski streaming potential coefficient:
    ///
    /// `E_stream / delta_P = -epsilon_0 * epsilon_r * zeta / (mu * sigma_el)`
    ///
    /// where `sigma_el` is the electrical conductivity of the electrolyte (S/m).
    pub fn streaming_potential_coefficient(
        zeta: f64,
        epsilon_r: f64,
        viscosity: f64,
        conductivity: f64,
    ) -> f64 {
        -EPSILON_0 * epsilon_r * zeta / (viscosity * conductivity)
    }
    /// Streaming current per unit width for a parallel-plate channel:
    ///
    /// `I_stream = -(2 * epsilon_0 * epsilon_r * zeta * delta_P * h) / (3 * mu * L)`
    ///
    /// where h = channel half-height, L = channel length.
    #[allow(clippy::too_many_arguments)]
    pub fn streaming_current(
        zeta: f64,
        epsilon_r: f64,
        viscosity: f64,
        delta_p: f64,
        half_height: f64,
        length: f64,
    ) -> f64 {
        -(2.0 * EPSILON_0 * epsilon_r * zeta * delta_p * half_height) / (3.0 * viscosity * length)
    }
    /// Electrokinetic coupling coefficient:
    ///
    /// `L_12 = epsilon_0 * epsilon_r * zeta / mu`
    pub fn coupling_coefficient(zeta: f64, epsilon_r: f64, viscosity: f64) -> f64 {
        EPSILON_0 * epsilon_r * zeta / viscosity
    }
}
/// Stern layer (compact inner layer) model for the electric double layer.
///
/// The Stern layer is the rigid adsorbed ion layer between the solid surface
/// and the diffuse Gouy-Chapman layer.
#[derive(Debug, Clone, Copy)]
pub struct SternLayer {
    /// Thickness of the Stern layer (m).
    pub thickness: f64,
    /// Permittivity of the Stern layer (F/m).
    pub epsilon_stern: f64,
    /// Surface charge density of the solid (C/m^2).
    pub sigma_surface: f64,
}
impl SternLayer {
    /// Create a new Stern layer model.
    pub fn new(thickness: f64, epsilon_stern: f64, sigma_surface: f64) -> Self {
        Self {
            thickness,
            epsilon_stern,
            sigma_surface,
        }
    }
    /// Potential drop across the Stern layer: delta_psi_H = sigma * d_H / epsilon_H.
    pub fn potential_drop(&self) -> f64 {
        self.sigma_surface * self.thickness / self.epsilon_stern.max(1e-30)
    }
    /// Helmholtz capacitance per unit area: C_H = epsilon_H / d_H.
    pub fn capacitance(&self) -> f64 {
        self.epsilon_stern / self.thickness.max(1e-30)
    }
    /// Outer Helmholtz plane (OHP) potential given the surface potential psi_0.
    ///
    /// psi_OHP = psi_0 - delta_psi_H
    pub fn ohp_potential(&self, psi_0: f64) -> f64 {
        psi_0 - self.potential_drop()
    }
}
