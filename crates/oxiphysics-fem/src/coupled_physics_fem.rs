#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled multi-physics FEM: thermo-mechanical, thermo-electric,
//! electro-magnetic, chemo-mechanical, hydro-mechanical, and staggered /
//! monolithic coupled solvers with Aitken relaxation.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Dot product of two slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Euclidean norm.
fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// Element-wise vector addition.
fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

/// Scale a vector.
fn vec_scale(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|x| x * s).collect()
}

/// Dense matrix–vector product (row-major).
fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

// ---------------------------------------------------------------------------
// § 1  ThermoMechanical — thermal expansion / Duhamel-Neumann
// ---------------------------------------------------------------------------

/// Isotropic linear thermo-mechanical material.
///
/// The Duhamel-Neumann constitutive law gives:
/// σ = C : ε − β (T − T_ref) I
/// where β = α (3λ + 2μ) = α · E / (1 − 2ν).
#[derive(Debug, Clone)]
pub struct ThermoMechanical {
    /// Young's modulus \[Pa\].
    pub e_mod: f64,
    /// Poisson ratio.
    pub nu: f64,
    /// Linear coefficient of thermal expansion \[1/K\].
    pub alpha: f64,
    /// Thermal conductivity \[W/(m·K)\].
    pub kappa: f64,
    /// Density \[kg/m³\].
    pub rho: f64,
    /// Specific heat capacity \[J/(kg·K)\].
    pub cp: f64,
    /// Reference temperature \[K\].
    pub t_ref: f64,
}

impl ThermoMechanical {
    /// Create a new `ThermoMechanical` material.
    pub fn new(e_mod: f64, nu: f64, alpha: f64, kappa: f64, rho: f64, cp: f64, t_ref: f64) -> Self {
        Self {
            e_mod,
            nu,
            alpha,
            kappa,
            rho,
            cp,
            t_ref,
        }
    }

    /// Lamé constant λ.
    pub fn lambda(&self) -> f64 {
        self.e_mod * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
    }

    /// Shear modulus μ = G.
    pub fn mu(&self) -> f64 {
        self.e_mod / (2.0 * (1.0 + self.nu))
    }

    /// Thermal stress modulus β = α · E / (1 − 2ν).
    pub fn beta(&self) -> f64 {
        self.alpha * self.e_mod / (1.0 - 2.0 * self.nu)
    }

    /// Build 3-D isotropic 6×6 elastic stiffness tensor (Voigt).
    pub fn stiffness_tensor(&self) -> [[f64; 6]; 6] {
        let lam = self.lambda();
        let mu = self.mu();
        let mut c = [[0.0f64; 6]; 6];
        c[0][0] = lam + 2.0 * mu;
        c[1][1] = lam + 2.0 * mu;
        c[2][2] = lam + 2.0 * mu;
        c[0][1] = lam;
        c[1][0] = lam;
        c[0][2] = lam;
        c[2][0] = lam;
        c[1][2] = lam;
        c[2][1] = lam;
        c[3][3] = mu;
        c[4][4] = mu;
        c[5][5] = mu;
        c
    }

    /// Thermal eigenstrain vector (Voigt): ε_th = α (T − T_ref) \[1,1,1,0,0,0\].
    pub fn thermal_eigenstrain(&self, temp: f64) -> [f64; 6] {
        let delta_t = temp - self.t_ref;
        [
            self.alpha * delta_t,
            self.alpha * delta_t,
            self.alpha * delta_t,
            0.0,
            0.0,
            0.0,
        ]
    }

    /// Duhamel-Neumann stress: σ = C:(ε − ε_th).
    pub fn duhamel_neumann_stress(&self, strain: &[f64; 6], temp: f64) -> [f64; 6] {
        let eps_th = self.thermal_eigenstrain(temp);
        let c = self.stiffness_tensor();
        let mut sigma = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                sigma[i] += c[i][j] * (strain[j] - eps_th[j]);
            }
        }
        sigma
    }

    /// Heat equation residual: ρ c_p ∂T/∂t − ∇·(κ ∇T) − Q_mech.
    ///
    /// Returns the thermal coupling source Q_mech = β T ∇·u̇
    /// (volumetric heating from mechanical dissipation, linearized).
    pub fn mechanical_heat_source(&self, temp: f64, div_u_dot: f64) -> f64 {
        self.beta() * temp * div_u_dot
    }
}

// ---------------------------------------------------------------------------
// § 2  ThermoElectric — Seebeck / Peltier / Thomson / Joule
// ---------------------------------------------------------------------------

/// Thermo-electric material encapsulating Seebeck, Peltier and Thomson effects.
///
/// The coupled equations are:
///   J = σ_e (−∇φ − S ∇T)       (Seebeck, current density)
///   J_q = Π J − κ ∇T            (Peltier + Fourier heat flux)
///   Q_joule = J · E = J²/σ_e    (Joule heating)
///
/// where Π = S T is the Peltier coefficient and S is Seebeck coefficient.
#[derive(Debug, Clone)]
pub struct ThermoElectric {
    /// Electrical conductivity σ_e \[S/m\].
    pub sigma_e: f64,
    /// Seebeck coefficient S \[V/K\].
    pub seebeck: f64,
    /// Thermal conductivity κ \[W/(m·K)\].
    pub kappa: f64,
    /// Reference temperature \[K\].
    pub t_ref: f64,
    /// Thomson coefficient τ \[V/K\] (= T · dS/dT; zero for constant S).
    pub thomson: f64,
}

impl ThermoElectric {
    /// Create a new `ThermoElectric` material.
    pub fn new(sigma_e: f64, seebeck: f64, kappa: f64, t_ref: f64) -> Self {
        Self {
            sigma_e,
            seebeck,
            kappa,
            t_ref,
            thomson: 0.0,
        }
    }

    /// Peltier coefficient Π = S · T \[V\].
    pub fn peltier(&self, temp: f64) -> f64 {
        self.seebeck * temp
    }

    /// Figure of merit ZT = S² σ_e T / κ (dimensionless).
    pub fn zt(&self, temp: f64) -> f64 {
        self.seebeck * self.seebeck * self.sigma_e * temp / self.kappa
    }

    /// Seebeck voltage ΔV = S · ΔT.
    pub fn seebeck_voltage(&self, delta_t: f64) -> f64 {
        self.seebeck * delta_t
    }

    /// Current density J from electric field E and temperature gradient ∇T.
    ///
    /// J = σ_e (E − S ∇T)
    pub fn current_density(&self, e_field: f64, grad_t: f64) -> f64 {
        self.sigma_e * (e_field - self.seebeck * grad_t)
    }

    /// Joule heating per unit volume Q = J² / σ_e.
    pub fn joule_heating(&self, current: f64) -> f64 {
        current * current / self.sigma_e
    }

    /// Thomson heat source per unit volume Q_t = τ · J · ∇T.
    pub fn thomson_heat(&self, current: f64, grad_t: f64) -> f64 {
        self.thomson * current * grad_t
    }

    /// Total heat flux q = Π J − κ ∇T.
    pub fn heat_flux(&self, temp: f64, current: f64, grad_t: f64) -> f64 {
        self.peltier(temp) * current - self.kappa * grad_t
    }
}

// ---------------------------------------------------------------------------
// § 3  ElectroMagnetic — eddy currents, Lorentz body force
// ---------------------------------------------------------------------------

/// Electromagnetic material / field data.
///
/// Implements a simplified low-frequency (eddy-current) model and
/// the Lorentz body force f = J × B.
#[derive(Debug, Clone)]
pub struct ElectroMagnetic {
    /// Electrical conductivity \[S/m\].
    pub sigma_e: f64,
    /// Magnetic permeability μ = μ_r μ_0 \[H/m\].
    pub mu_r: f64,
    /// Vacuum permeability \[H/m\].
    pub mu_0: f64,
    /// Characteristic frequency \[Hz\].
    pub freq: f64,
}

impl ElectroMagnetic {
    /// Create a new `ElectroMagnetic` material.
    pub fn new(sigma_e: f64, mu_r: f64, freq: f64) -> Self {
        Self {
            sigma_e,
            mu_r,
            mu_0: 4.0 * std::f64::consts::PI * 1e-7,
            freq,
        }
    }

    /// Magnetic permeability μ = μ_r μ_0 \[H/m\].
    pub fn mu(&self) -> f64 {
        self.mu_r * self.mu_0
    }

    /// Skin depth δ = sqrt(2 / (ω μ σ)) \[m\].
    pub fn skin_depth(&self) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * self.freq;
        (2.0 / (omega * self.mu() * self.sigma_e)).sqrt()
    }

    /// Eddy current density magnitude from induced EMF.
    ///
    /// J_eddy ≈ σ_e · E_induced (simplified).
    pub fn eddy_current(&self, e_induced: f64) -> f64 {
        self.sigma_e * e_induced
    }

    /// Lorentz body force f = J × B (scalar for aligned J, B).
    pub fn lorentz_force(&self, j: f64, b: f64) -> f64 {
        j * b
    }

    /// Magnetic Reynolds number Re_m = μ σ_e v L.
    pub fn magnetic_reynolds(&self, velocity: f64, length: f64) -> f64 {
        self.mu() * self.sigma_e * velocity * length
    }

    /// Joule dissipation density P = J²/σ_e \[W/m³\].
    pub fn joule_dissipation(&self, j: f64) -> f64 {
        j * j / self.sigma_e
    }

    /// Maxwell stress tensor component T_ij = (1/μ)(B_i B_j − δ_ij B²/2).
    ///
    /// For 1-D: T_11 = B²/(2μ).
    pub fn maxwell_stress_1d(&self, b: f64) -> f64 {
        b * b / (2.0 * self.mu())
    }
}

// ---------------------------------------------------------------------------
// § 4  ChemoMechanical — concentration-induced stress / diffusion coupling
// ---------------------------------------------------------------------------

/// Chemo-mechanical coupling material.
///
/// Models concentration-induced eigenstrain: ε_c = β_c (c − c_ref) I
/// and deformation-enhanced diffusion via mechanical driving force.
#[derive(Debug, Clone)]
pub struct ChemoMechanical {
    /// Young's modulus \[Pa\].
    pub e_mod: f64,
    /// Poisson ratio.
    pub nu: f64,
    /// Chemical expansion coefficient β_c \[m³/mol\] (or dimensionless).
    pub beta_c: f64,
    /// Reference concentration c_ref \[mol/m³\].
    pub c_ref: f64,
    /// Diffusivity D \[m²/s\].
    pub diffusivity: f64,
    /// Partial molar volume Ω \[m³/mol\].
    pub partial_molar_vol: f64,
    /// Temperature \[K\] (for chemical potential).
    pub temp: f64,
}

impl ChemoMechanical {
    /// Create a new `ChemoMechanical` material.
    pub fn new(
        e_mod: f64,
        nu: f64,
        beta_c: f64,
        c_ref: f64,
        diffusivity: f64,
        partial_molar_vol: f64,
        temp: f64,
    ) -> Self {
        Self {
            e_mod,
            nu,
            beta_c,
            c_ref,
            diffusivity,
            partial_molar_vol,
            temp,
        }
    }

    /// Concentration eigenstrain ε_c = β_c (c − c_ref) \[1,1,1,0,0,0\].
    pub fn concentration_eigenstrain(&self, conc: f64) -> [f64; 6] {
        let ec = self.beta_c * (conc - self.c_ref);
        [ec, ec, ec, 0.0, 0.0, 0.0]
    }

    /// Concentration-induced stress: σ = C:(ε − ε_c).
    pub fn chemo_stress(&self, strain: &[f64; 6], conc: f64) -> [f64; 6] {
        let eps_c = self.concentration_eigenstrain(conc);
        let lam = self.e_mod * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu));
        let mu = self.e_mod / (2.0 * (1.0 + self.nu));
        let mut c_mat = [[0.0f64; 6]; 6];
        c_mat[0][0] = lam + 2.0 * mu;
        c_mat[1][1] = lam + 2.0 * mu;
        c_mat[2][2] = lam + 2.0 * mu;
        c_mat[0][1] = lam;
        c_mat[1][0] = lam;
        c_mat[0][2] = lam;
        c_mat[2][0] = lam;
        c_mat[1][2] = lam;
        c_mat[2][1] = lam;
        c_mat[3][3] = mu;
        c_mat[4][4] = mu;
        c_mat[5][5] = mu;
        let mut sigma = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                sigma[i] += c_mat[i][j] * (strain[j] - eps_c[j]);
            }
        }
        sigma
    }

    /// Chemical potential μ_c = μ_0 + Ω σ_h (mechanically modified).
    ///
    /// σ_h = (σ_11 + σ_22 + σ_33) / 3 is the hydrostatic stress.
    pub fn chemical_potential(&self, stress: &[f64; 6]) -> f64 {
        let sigma_h = (stress[0] + stress[1] + stress[2]) / 3.0;
        // μ_c = R T ln(c/c_ref) + Ω σ_h  (linearized around c_ref)
        self.partial_molar_vol * sigma_h
    }

    /// Effective diffusivity with mechanical coupling D_eff = D exp(-Ω σ_h / RT).
    pub fn effective_diffusivity(&self, stress: &[f64; 6]) -> f64 {
        let sigma_h = (stress[0] + stress[1] + stress[2]) / 3.0;
        let r_gas = 8.314_462_618;
        let exponent = -self.partial_molar_vol * sigma_h / (r_gas * self.temp);
        self.diffusivity * exponent.exp().clamp(0.01, 100.0)
    }

    /// Flux J = −D_eff ∇c.
    pub fn diffusion_flux(&self, grad_c: f64, stress: &[f64; 6]) -> f64 {
        -self.effective_diffusivity(stress) * grad_c
    }
}

// ---------------------------------------------------------------------------
// § 5  HydroMechanical — Biot consolidation / poromechanics
// ---------------------------------------------------------------------------

/// Biot poroelastic material for coupled hydro-mechanical analysis.
///
/// The Biot consolidation equations are:
///   ∇·σ' + ∇·(α p I) + f = 0        (mechanical equilibrium, effective stress)
///   (α ∇·u̇ + p̈/M) + ∇·q = 0         (mass conservation / continuity)
///   q = −(k/μ_f) (∇p − ρ_f g)        (Darcy flow)
#[derive(Debug, Clone)]
pub struct HydroMechanical {
    /// Drained Young's modulus \[Pa\].
    pub e_drained: f64,
    /// Drained Poisson ratio.
    pub nu_drained: f64,
    /// Biot-Willis coefficient α (0 < α ≤ 1).
    pub biot_alpha: f64,
    /// Biot modulus M \[Pa\].
    pub biot_m: f64,
    /// Intrinsic permeability k \[m²\].
    pub permeability: f64,
    /// Fluid dynamic viscosity μ_f \[Pa·s\].
    pub fluid_viscosity: f64,
    /// Fluid density ρ_f \[kg/m³\].
    pub fluid_density: f64,
}

impl HydroMechanical {
    /// Create a new `HydroMechanical` material.
    pub fn new(
        e_drained: f64,
        nu_drained: f64,
        biot_alpha: f64,
        biot_m: f64,
        permeability: f64,
        fluid_viscosity: f64,
        fluid_density: f64,
    ) -> Self {
        Self {
            e_drained,
            nu_drained,
            biot_alpha,
            biot_m,
            permeability,
            fluid_viscosity,
            fluid_density,
        }
    }

    /// Drained bulk modulus K_d = E / (3(1-2ν)).
    pub fn k_drained(&self) -> f64 {
        self.e_drained / (3.0 * (1.0 - 2.0 * self.nu_drained))
    }

    /// Drained shear modulus G = E / (2(1+ν)).
    pub fn g_drained(&self) -> f64 {
        self.e_drained / (2.0 * (1.0 + self.nu_drained))
    }

    /// Undrained bulk modulus K_u = K_d + α² M.
    pub fn k_undrained(&self) -> f64 {
        self.k_drained() + self.biot_alpha * self.biot_alpha * self.biot_m
    }

    /// Skempton coefficient B = α M / K_u.
    pub fn skempton_b(&self) -> f64 {
        self.biot_alpha * self.biot_m / self.k_undrained()
    }

    /// Effective stress σ' = σ − α p I.
    ///
    /// Returns the effective stress vector (Voigt) from total stress and pore pressure.
    pub fn effective_stress(&self, total_stress: &[f64; 6], pore_pressure: f64) -> [f64; 6] {
        let mut sigma_eff = *total_stress;
        let ap = self.biot_alpha * pore_pressure;
        sigma_eff[0] -= ap;
        sigma_eff[1] -= ap;
        sigma_eff[2] -= ap;
        sigma_eff
    }

    /// Total stress from effective stress and pore pressure.
    pub fn total_stress(&self, sigma_eff: &[f64; 6], pore_pressure: f64) -> [f64; 6] {
        let mut sigma = *sigma_eff;
        let ap = self.biot_alpha * pore_pressure;
        sigma[0] += ap;
        sigma[1] += ap;
        sigma[2] += ap;
        sigma
    }

    /// Darcy flux q = −(k/μ) ∇p \[m/s\].
    pub fn darcy_flux(&self, grad_p: f64) -> f64 {
        -(self.permeability / self.fluid_viscosity) * grad_p
    }

    /// Hydraulic conductivity K_h = k ρ_f g / μ_f \[m/s\].
    pub fn hydraulic_conductivity(&self) -> f64 {
        const G: f64 = 9.81;
        self.permeability * self.fluid_density * G / self.fluid_viscosity
    }

    /// Consolidation coefficient c_v = k M_v / μ_f where M_v is the oedometric modulus.
    pub fn consolidation_coeff(&self) -> f64 {
        let m_v = self.k_drained() + 4.0 * self.g_drained() / 3.0;
        self.permeability * m_v / self.fluid_viscosity
    }
}

// ---------------------------------------------------------------------------
// § 6  CoupledSolver — staggered and monolithic schemes
// ---------------------------------------------------------------------------

/// Convergence criteria for coupled solvers.
#[derive(Debug, Clone)]
pub struct ConvergenceCriteria {
    /// Maximum allowed relative displacement residual.
    pub tol_disp: f64,
    /// Maximum allowed relative pressure / temperature residual.
    pub tol_scalar: f64,
    /// Maximum number of staggered iterations.
    pub max_iter: usize,
}

impl ConvergenceCriteria {
    /// Create convergence criteria with standard tolerances.
    pub fn new(tol_disp: f64, tol_scalar: f64, max_iter: usize) -> Self {
        Self {
            tol_disp,
            tol_scalar,
            max_iter,
        }
    }
}

impl Default for ConvergenceCriteria {
    fn default() -> Self {
        Self {
            tol_disp: 1e-6,
            tol_scalar: 1e-6,
            max_iter: 100,
        }
    }
}

/// Staggered scheme type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverScheme {
    /// Staggered (operator-splitting) approach.
    Staggered,
    /// Monolithic (fully-coupled) approach.
    Monolithic,
    /// Iteratively-staggered with convergence check.
    IterativeStaggered,
}

/// Convergence result.
#[derive(Debug, Clone)]
pub struct CoupledResult {
    /// Displacement field.
    pub displacement: Vec<f64>,
    /// Scalar field (temperature, pressure, or concentration).
    pub scalar_field: Vec<f64>,
    /// Number of iterations taken.
    pub n_iter: usize,
    /// Final displacement residual.
    pub res_disp: f64,
    /// Final scalar residual.
    pub res_scalar: f64,
    /// Whether the solve converged.
    pub converged: bool,
}

/// Aitken dynamic relaxation accelerator.
///
/// Computes the relaxation parameter ω_{n+1} via Aitken's Δ² method.
#[derive(Debug, Clone)]
pub struct AitkenRelaxation {
    /// Current relaxation factor ω.
    pub omega: f64,
    /// Previous residual vector.
    pub res_prev: Vec<f64>,
    /// Minimum allowed ω.
    pub omega_min: f64,
    /// Maximum allowed ω.
    pub omega_max: f64,
}

impl AitkenRelaxation {
    /// Create a new Aitken relaxation with initial ω.
    pub fn new(omega_init: f64) -> Self {
        Self {
            omega: omega_init,
            res_prev: Vec::new(),
            omega_min: 1e-4,
            omega_max: 1.0,
        }
    }

    /// Update ω using Aitken's formula and return the new relaxation factor.
    ///
    /// ω_{n+1} = −ω_n · (r_n^T (r_{n+1} − r_n)) / |r_{n+1} − r_n|²
    pub fn update(&mut self, res_new: &[f64]) -> f64 {
        if self.res_prev.is_empty() {
            self.res_prev = res_new.to_vec();
            return self.omega;
        }
        let dr: Vec<f64> = res_new
            .iter()
            .zip(self.res_prev.iter())
            .map(|(a, b)| a - b)
            .collect();
        let dr2 = dot(&dr, &dr);
        if dr2 < 1.0e-60 {
            self.res_prev = res_new.to_vec();
            return self.omega;
        }
        let num = dot(&self.res_prev, &dr);
        self.omega = (-self.omega * num / dr2)
            .abs()
            .clamp(self.omega_min, self.omega_max);
        self.res_prev = res_new.to_vec();
        self.omega
    }

    /// Reset to initial omega.
    pub fn reset(&mut self, omega_init: f64) {
        self.omega = omega_init;
        self.res_prev.clear();
    }
}

/// Coupled physics solver supporting staggered and monolithic schemes.
///
/// Manages the coupling between a displacement field and a scalar field
/// (temperature, pressure, or concentration), applying the chosen
/// solution scheme with optional Aitken relaxation.
#[derive(Debug, Clone)]
pub struct CoupledSolver {
    /// Number of displacement degrees of freedom.
    pub n_disp_dof: usize,
    /// Number of scalar degrees of freedom.
    pub n_scalar_dof: usize,
    /// Solution scheme.
    pub scheme: SolverScheme,
    /// Convergence criteria.
    pub criteria: ConvergenceCriteria,
    /// Aitken relaxation.
    pub relaxation: AitkenRelaxation,
}

impl CoupledSolver {
    /// Create a new `CoupledSolver`.
    pub fn new(
        n_disp_dof: usize,
        n_scalar_dof: usize,
        scheme: SolverScheme,
        criteria: ConvergenceCriteria,
    ) -> Self {
        Self {
            n_disp_dof,
            n_scalar_dof,
            scheme,
            criteria,
            relaxation: AitkenRelaxation::new(1.0),
        }
    }

    /// Solve the coupled system.
    ///
    /// Accepts a right-hand side for displacements `f_u` and a right-hand
    /// side for the scalar field `f_s`.  The coupling matrices C_us and
    /// C_su are provided as dense row-major matrices.
    ///
    /// Returns a `CoupledResult` with the converged solution.
    pub fn solve(
        &mut self,
        k_u: &[Vec<f64>],
        k_s: &[Vec<f64>],
        c_us: &[Vec<f64>],
        c_su: &[Vec<f64>],
        f_u: &[f64],
        f_s: &[f64],
    ) -> CoupledResult {
        let mut u = vec![0.0f64; self.n_disp_dof];
        let mut s = vec![0.0f64; self.n_scalar_dof];
        let mut n_iter = 0;
        let mut res_u = f64::MAX;
        let mut res_s = f64::MAX;
        let mut converged = false;

        match self.scheme {
            SolverScheme::Monolithic => {
                // Build the coupled block system and solve once
                let (u_sol, s_sol) = self.monolithic_solve(k_u, k_s, c_us, c_su, f_u, f_s);
                u = u_sol;
                s = s_sol;
                res_u = 0.0;
                res_s = 0.0;
                n_iter = 1;
                converged = true;
            }
            SolverScheme::Staggered => {
                // One-pass staggered: first displacement, then scalar
                s = self.solve_scalar_field(k_s, f_s);
                let f_u_corrected = self.coupling_rhs(c_us, &s, f_u);
                u = self.solve_disp_field(k_u, &f_u_corrected);
                n_iter = 1;
                converged = true;
                res_u = 0.0;
                res_s = 0.0;
            }
            SolverScheme::IterativeStaggered => {
                self.relaxation.reset(1.0);
                for iter in 0..self.criteria.max_iter {
                    let u_prev = u.clone();
                    let s_prev = s.clone();

                    // Sub-step 1: solve scalar field with fixed u
                    let f_s_c = self.coupling_rhs(c_su, &u, f_s);
                    let s_new = self.solve_scalar_field(k_s, &f_s_c);

                    // Sub-step 2: solve displacement with updated s
                    let f_u_c = self.coupling_rhs(c_us, &s_new, f_u);
                    let u_new = self.solve_disp_field(k_u, &f_u_c);

                    // Compute residuals
                    let du: Vec<f64> = u_new
                        .iter()
                        .zip(u_prev.iter())
                        .map(|(a, b)| a - b)
                        .collect();
                    let ds: Vec<f64> = s_new
                        .iter()
                        .zip(s_prev.iter())
                        .map(|(a, b)| a - b)
                        .collect();
                    res_u = norm(&du) / (norm(&u_new).max(1.0e-30));
                    res_s = norm(&ds) / (norm(&s_new).max(1.0e-30));

                    // Aitken relaxation on scalar residual
                    let omega = self.relaxation.update(&ds);
                    s = s_new
                        .iter()
                        .zip(s_prev.iter())
                        .map(|(sn, sp)| sp + omega * (sn - sp))
                        .collect();
                    u = u_new;

                    n_iter = iter + 1;
                    if res_u < self.criteria.tol_disp && res_s < self.criteria.tol_scalar {
                        converged = true;
                        break;
                    }
                }
            }
        }

        CoupledResult {
            displacement: u,
            scalar_field: s,
            n_iter,
            res_disp: res_u,
            res_scalar: res_s,
            converged,
        }
    }

    /// Solve a linear system K x = f using conjugate gradient (no preconditioner).
    fn cg_solve(k: &[Vec<f64>], f: &[f64]) -> Vec<f64> {
        let n = f.len();
        if n == 0 {
            return vec![];
        }
        let mut x = vec![0.0; n];
        let mut r = f.to_vec();
        let mut p = r.clone();
        let mut rs = dot(&r, &r);
        for _ in 0..n * 10 {
            let kp = mat_vec(k, &p);
            let alpha = rs / dot(&p, &kp).max(1e-300);
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * kp[i];
            }
            let rs_new = dot(&r, &r);
            if rs_new.sqrt() < 1e-12 * norm(f).max(1e-30) {
                break;
            }
            let beta = rs_new / rs.max(1e-300);
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }
            rs = rs_new;
        }
        x
    }

    /// Solve displacement sub-system K_u u = f_u.
    fn solve_disp_field(&self, k_u: &[Vec<f64>], f_u: &[f64]) -> Vec<f64> {
        Self::cg_solve(k_u, f_u)
    }

    /// Solve scalar sub-system K_s s = f_s.
    fn solve_scalar_field(&self, k_s: &[Vec<f64>], f_s: &[f64]) -> Vec<f64> {
        Self::cg_solve(k_s, f_s)
    }

    /// Compute coupling RHS: f_coupled = f + C x_other.
    fn coupling_rhs(&self, c: &[Vec<f64>], x_other: &[f64], f: &[f64]) -> Vec<f64> {
        if c.is_empty() || c[0].is_empty() {
            return f.to_vec();
        }
        let cx = mat_vec(c, x_other);
        let nf = f.len();
        (0..nf)
            .map(|i| f[i] + if i < cx.len() { cx[i] } else { 0.0 })
            .collect()
    }

    /// Monolithic block solve.
    ///
    /// Builds the 2×2 block system \[K_u, C_us; C_su, K_s\] \[u; s\] = \[f_u; f_s\]
    /// and solves via block-Gauss elimination.
    fn monolithic_solve(
        &self,
        k_u: &[Vec<f64>],
        k_s: &[Vec<f64>],
        c_us: &[Vec<f64>],
        c_su: &[Vec<f64>],
        f_u: &[f64],
        f_s: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        let nu = self.n_disp_dof;
        let ns = self.n_scalar_dof;
        let ntot = nu + ns;
        // Build full block matrix
        let mut a = vec![vec![0.0; ntot]; ntot];
        for i in 0..nu {
            for j in 0..nu {
                if i < k_u.len() && j < k_u[i].len() {
                    a[i][j] = k_u[i][j];
                }
            }
        }
        for i in 0..nu {
            for j in 0..ns {
                if i < c_us.len() && j < c_us[i].len() {
                    a[i][nu + j] = c_us[i][j];
                }
            }
        }
        for i in 0..ns {
            for j in 0..nu {
                if i < c_su.len() && j < c_su[i].len() {
                    a[nu + i][j] = c_su[i][j];
                }
            }
        }
        for i in 0..ns {
            for j in 0..ns {
                if i < k_s.len() && j < k_s[i].len() {
                    a[nu + i][nu + j] = k_s[i][j];
                }
            }
        }
        let mut rhs = vec![0.0; ntot];
        for i in 0..nu {
            rhs[i] = if i < f_u.len() { f_u[i] } else { 0.0 };
        }
        for i in 0..ns {
            rhs[nu + i] = if i < f_s.len() { f_s[i] } else { 0.0 };
        }

        let x = Self::cg_solve(&a, &rhs);
        let u: Vec<f64> = x[..nu].to_vec();
        let s: Vec<f64> = x[nu..].to_vec();
        (u, s)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ──────────────────────────────────────────────────────────────

    fn steel_tm() -> ThermoMechanical {
        ThermoMechanical::new(210e9, 0.3, 12e-6, 50.0, 7850.0, 500.0, 293.15)
    }

    fn bi2te3() -> ThermoElectric {
        ThermoElectric::new(1e5, 200e-6, 1.5, 300.0)
    }

    fn copper_em() -> ElectroMagnetic {
        ElectroMagnetic::new(5.96e7, 1.0, 50.0)
    }

    fn li_ion_cm() -> ChemoMechanical {
        ChemoMechanical::new(100e9, 0.3, 3.64e-6, 10e3, 1e-14, 8.9e-6, 300.0)
    }

    fn clay_hm() -> HydroMechanical {
        // Physically consistent soft clay:
        // K_d = 10e6/(3*(1-0.6)) = 8.33e6 Pa
        // For Skempton B ≤ 1: M ≤ K_d / (α(1−α)) = 8.33e6/(0.8*0.2) = 52e6 Pa
        HydroMechanical::new(10e6, 0.3, 0.8, 40e6, 1e-17, 1e-3, 1000.0)
    }

    fn diag_k(n: usize, val: f64) -> Vec<Vec<f64>> {
        let mut k = vec![vec![0.0; n]; n];
        for i in 0..n {
            k[i][i] = val;
        }
        k
    }

    // ── § 1  ThermoMechanical ─────────────────────────────────────────────────

    #[test]
    fn test_tm_lambda() {
        let m = steel_tm();
        let expected = 210e9 * 0.3 / ((1.0 + 0.3) * (1.0 - 0.6));
        assert!((m.lambda() - expected).abs() < 1e3);
    }

    #[test]
    fn test_tm_mu() {
        let m = steel_tm();
        let expected = 210e9 / (2.0 * 1.3);
        assert!((m.mu() - expected).abs() < 1e3);
    }

    #[test]
    fn test_tm_beta_positive() {
        assert!(steel_tm().beta() > 0.0);
    }

    #[test]
    fn test_tm_stiffness_symmetry() {
        let c = steel_tm().stiffness_tensor();
        for i in 0..6 {
            for j in 0..6 {
                assert!((c[i][j] - c[j][i]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_tm_stiffness_positive_definite_diagonal() {
        let c = steel_tm().stiffness_tensor();
        for i in 0..6 {
            assert!(c[i][i] > 0.0, "C[{i}][{i}] should be positive");
        }
    }

    #[test]
    fn test_tm_thermal_eigenstrain_zero_at_ref() {
        let m = steel_tm();
        let eps_th = m.thermal_eigenstrain(m.t_ref);
        for v in &eps_th {
            assert!(v.abs() < 1e-20);
        }
    }

    #[test]
    fn test_tm_thermal_eigenstrain_positive_above_ref() {
        let m = steel_tm();
        let eps_th = m.thermal_eigenstrain(m.t_ref + 100.0);
        assert!(eps_th[0] > 0.0);
        assert_eq!(eps_th[3], 0.0); // no shear
    }

    #[test]
    fn test_tm_duhamel_neumann_no_thermal_strain() {
        let m = steel_tm();
        let strain = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma = m.duhamel_neumann_stress(&strain, m.t_ref);
        // At reference temp, no thermal effect
        let c = m.stiffness_tensor();
        let sigma_check: Vec<f64> = (0..6)
            .map(|i| (0..6).map(|j| c[i][j] * strain[j]).sum())
            .collect();
        for k in 0..6 {
            assert!((sigma[k] - sigma_check[k]).abs() < 1e3);
        }
    }

    #[test]
    fn test_tm_mechanical_heat_source_zero_divu() {
        let m = steel_tm();
        assert_eq!(m.mechanical_heat_source(300.0, 0.0), 0.0);
    }

    // ── § 2  ThermoElectric ───────────────────────────────────────────────────

    #[test]
    fn test_te_peltier_coefficient() {
        let te = bi2te3();
        assert!((te.peltier(300.0) - 200e-6 * 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_te_zt_positive() {
        assert!(bi2te3().zt(300.0) > 0.0);
    }

    #[test]
    fn test_te_seebeck_voltage_linear() {
        let te = bi2te3();
        let v1 = te.seebeck_voltage(10.0);
        let v2 = te.seebeck_voltage(20.0);
        assert!((v2 - 2.0 * v1).abs() < 1e-15);
    }

    #[test]
    fn test_te_current_density_zero_grad() {
        let te = bi2te3();
        assert!((te.current_density(1.0, 0.0) - te.sigma_e).abs() < 1e3);
    }

    #[test]
    fn test_te_joule_heating_positive() {
        let te = bi2te3();
        assert!(te.joule_heating(1e3) > 0.0);
    }

    #[test]
    fn test_te_heat_flux_sign() {
        let te = bi2te3();
        // Heat flows from hot to cold: negative grad_t, positive flux contribution
        let q = te.heat_flux(300.0, 0.0, -10.0);
        assert!(q > 0.0); // -κ * (-10) > 0
    }

    #[test]
    fn test_te_thomson_zero_when_coefficient_zero() {
        let te = bi2te3();
        assert_eq!(te.thomson_heat(1e3, 5.0), 0.0);
    }

    // ── § 3  ElectroMagnetic ──────────────────────────────────────────────────

    #[test]
    fn test_em_skin_depth_positive() {
        assert!(copper_em().skin_depth() > 0.0);
    }

    #[test]
    fn test_em_skin_depth_decreases_with_freq() {
        let em1 = ElectroMagnetic::new(5.96e7, 1.0, 50.0);
        let em2 = ElectroMagnetic::new(5.96e7, 1.0, 50000.0);
        assert!(em1.skin_depth() > em2.skin_depth());
    }

    #[test]
    fn test_em_lorentz_force_proportional() {
        let em = copper_em();
        let f1 = em.lorentz_force(1e3, 0.5);
        let f2 = em.lorentz_force(2e3, 0.5);
        assert!((f2 - 2.0 * f1).abs() < 1e-10);
    }

    #[test]
    fn test_em_magnetic_reynolds_positive() {
        assert!(copper_em().magnetic_reynolds(1.0, 0.1) > 0.0);
    }

    #[test]
    fn test_em_joule_dissipation_positive() {
        assert!(copper_em().joule_dissipation(1e3) > 0.0);
    }

    #[test]
    fn test_em_maxwell_stress_positive() {
        assert!(copper_em().maxwell_stress_1d(1.0) > 0.0);
    }

    #[test]
    fn test_em_mu_gt_vacuum() {
        let em = copper_em();
        assert!(em.mu() >= em.mu_0);
    }

    // ── § 4  ChemoMechanical ──────────────────────────────────────────────────

    #[test]
    fn test_cm_eigenstrain_zero_at_ref() {
        let m = li_ion_cm();
        let eps = m.concentration_eigenstrain(m.c_ref);
        for v in &eps {
            assert!(v.abs() < 1e-20);
        }
    }

    #[test]
    fn test_cm_eigenstrain_positive_above_ref() {
        let m = li_ion_cm();
        let eps = m.concentration_eigenstrain(m.c_ref + 1e3);
        assert!(eps[0] > 0.0);
        assert_eq!(eps[3], 0.0);
    }

    #[test]
    fn test_cm_chemo_stress_no_concentration_change() {
        let m = li_ion_cm();
        let strain = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma = m.chemo_stress(&strain, m.c_ref);
        // No concentration change: stress from strain alone
        let lam = m.e_mod * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
        let mu = m.e_mod / (2.0 * (1.0 + m.nu));
        let sigma11_expected = (lam + 2.0 * mu) * strain[0];
        assert!((sigma[0] - sigma11_expected).abs() < 1e3);
    }

    #[test]
    fn test_cm_chemical_potential_zero_hydrostatic() {
        let m = li_ion_cm();
        let sigma = [0.0; 6];
        assert_eq!(m.chemical_potential(&sigma), 0.0);
    }

    #[test]
    fn test_cm_effective_diffusivity_near_one_without_stress() {
        let m = li_ion_cm();
        let sigma = [0.0; 6];
        let d_eff = m.effective_diffusivity(&sigma);
        assert!((d_eff - m.diffusivity).abs() < 1e-25);
    }

    #[test]
    fn test_cm_diffusion_flux_sign() {
        let m = li_ion_cm();
        // Flux should be negative for positive grad_c
        let flux = m.diffusion_flux(100.0, &[0.0; 6]);
        assert!(flux < 0.0);
    }

    // ── § 5  HydroMechanical ─────────────────────────────────────────────────

    #[test]
    fn test_hm_k_drained_positive() {
        assert!(clay_hm().k_drained() > 0.0);
    }

    #[test]
    fn test_hm_k_undrained_gt_drained() {
        let m = clay_hm();
        assert!(m.k_undrained() > m.k_drained());
    }

    #[test]
    fn test_hm_skempton_b_range() {
        let b = clay_hm().skempton_b();
        assert!(b > 0.0 && b <= 1.0);
    }

    #[test]
    fn test_hm_effective_stress() {
        let m = clay_hm();
        let sigma_total = [-1e6f64, -1e6, -1e6, 0.0, 0.0, 0.0];
        let p = 5e5;
        let sigma_eff = m.effective_stress(&sigma_total, p);
        let expected = -1e6 - m.biot_alpha * p;
        assert!((sigma_eff[0] - expected).abs() < 1e3);
    }

    #[test]
    fn test_hm_total_from_effective_roundtrip() {
        let m = clay_hm();
        let sigma_eff = [-2e6f64, -2e6, -2e6, 0.0, 0.0, 0.0];
        let p = 3e5;
        let sigma_tot = m.total_stress(&sigma_eff, p);
        let sigma_eff_back = m.effective_stress(&sigma_tot, p);
        for k in 0..6 {
            assert!((sigma_eff[k] - sigma_eff_back[k]).abs() < 1e3);
        }
    }

    #[test]
    fn test_hm_darcy_flux_sign() {
        let m = clay_hm();
        // Positive pressure gradient → negative flux (flow from high to low)
        let q = m.darcy_flux(1e6);
        assert!(q < 0.0);
    }

    #[test]
    fn test_hm_hydraulic_conductivity_positive() {
        assert!(clay_hm().hydraulic_conductivity() > 0.0);
    }

    #[test]
    fn test_hm_consolidation_coeff_positive() {
        assert!(clay_hm().consolidation_coeff() > 0.0);
    }

    // ── § 6  CoupledSolver ────────────────────────────────────────────────────

    #[test]
    fn test_aitken_first_step_returns_init_omega() {
        let mut ar = AitkenRelaxation::new(1.0);
        let res = vec![0.1, 0.2, 0.3];
        let omega = ar.update(&res);
        assert!((omega - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_aitken_update_stays_in_bounds() {
        let mut ar = AitkenRelaxation::new(0.5);
        ar.update(&[1.0, 0.0]);
        let omega = ar.update(&[0.0, 1.0]);
        assert!(omega >= ar.omega_min && omega <= ar.omega_max);
    }

    #[test]
    fn test_aitken_reset() {
        let mut ar = AitkenRelaxation::new(0.5);
        ar.update(&[1.0, 2.0]);
        ar.reset(0.8);
        assert!((ar.omega - 0.8).abs() < 1e-10);
        assert!(ar.res_prev.is_empty());
    }

    #[test]
    fn test_coupled_solver_staggered_solves() {
        let n = 3;
        let k_u = diag_k(n, 1e9);
        let k_s = diag_k(n, 1.0);
        let c_us = vec![vec![0.0; n]; n];
        let c_su = vec![vec![0.0; n]; n];
        let f_u = vec![1e6; n];
        let f_s = vec![300.0; n];
        let mut solver = CoupledSolver::new(
            n,
            n,
            SolverScheme::Staggered,
            ConvergenceCriteria::default(),
        );
        let result = solver.solve(&k_u, &k_s, &c_us, &c_su, &f_u, &f_s);
        assert!(result.converged);
        for &u in &result.displacement {
            assert!(u.abs() < 1e-2);
        }
    }

    #[test]
    fn test_coupled_solver_monolithic_solves() {
        let n = 2;
        let k_u = diag_k(n, 2e9);
        let k_s = diag_k(n, 1.0);
        let c_us = vec![vec![0.0; n]; n];
        let c_su = vec![vec![0.0; n]; n];
        let f_u = vec![1e6; n];
        let f_s = vec![350.0; n];
        let mut solver = CoupledSolver::new(
            n,
            n,
            SolverScheme::Monolithic,
            ConvergenceCriteria::default(),
        );
        let result = solver.solve(&k_u, &k_s, &c_us, &c_su, &f_u, &f_s);
        assert!(result.converged);
        assert_eq!(result.n_iter, 1);
    }

    #[test]
    fn test_coupled_solver_iterative_staggered() {
        let n = 2;
        let k_u = diag_k(n, 1e9);
        let k_s = diag_k(n, 10.0);
        let c_us = vec![vec![0.0; n]; n];
        let c_su = vec![vec![0.0; n]; n];
        let f_u = vec![1e6; n];
        let f_s = vec![100.0; n];
        let mut solver = CoupledSolver::new(
            n,
            n,
            SolverScheme::IterativeStaggered,
            ConvergenceCriteria::new(1e-8, 1e-8, 50),
        );
        let result = solver.solve(&k_u, &k_s, &c_us, &c_su, &f_u, &f_s);
        // Should converge since uncoupled
        assert!(result.converged);
    }

    #[test]
    fn test_coupled_result_displacement_length() {
        let n = 4;
        let k_u = diag_k(n, 1e9);
        let k_s = diag_k(n, 1.0);
        let c = vec![vec![0.0; n]; n];
        let f_u = vec![1e6; n];
        let f_s = vec![300.0; n];
        let mut solver = CoupledSolver::new(
            n,
            n,
            SolverScheme::Staggered,
            ConvergenceCriteria::default(),
        );
        let result = solver.solve(&k_u, &k_s, &c, &c, &f_u, &f_s);
        assert_eq!(result.displacement.len(), n);
        assert_eq!(result.scalar_field.len(), n);
    }

    #[test]
    fn test_convergence_criteria_default() {
        let cc = ConvergenceCriteria::default();
        assert_eq!(cc.max_iter, 100);
        assert!(cc.tol_disp < 1e-4);
    }

    #[test]
    fn test_solver_scheme_variants_distinct() {
        assert_ne!(SolverScheme::Staggered, SolverScheme::Monolithic);
        assert_ne!(SolverScheme::Monolithic, SolverScheme::IterativeStaggered);
    }

    #[test]
    fn test_hm_biot_alpha_range() {
        let m = clay_hm();
        assert!(m.biot_alpha > 0.0 && m.biot_alpha <= 1.0);
    }

    #[test]
    fn test_em_eddy_current() {
        let em = copper_em();
        let j = em.eddy_current(1.0);
        assert!((j - em.sigma_e).abs() < 1.0);
    }

    #[test]
    fn test_tm_duhamel_compression_from_heating() {
        // Free thermal expansion → zero stress (but constrained → compressive stress)
        let m = steel_tm();
        // Simulate constrained heating: strain = 0, T > T_ref
        let strain = [0.0; 6];
        let sigma = m.duhamel_neumann_stress(&strain, m.t_ref + 50.0);
        // Should be compressive (negative normal stress)
        assert!(sigma[0] < 0.0);
    }
}
