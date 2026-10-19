// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled multi-physics FEM: thermo-mechanical, thermo-electric,
//! electro-magnetic, chemo-mechanical, hydro-mechanical, and staggered /
//! monolithic coupled solvers with Aitken relaxation.

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

/// Dense matrix–vector product y = A x for a flat row-major `n_rows × n_cols` matrix.
fn mat_vec(a: &[f64], n_rows: usize, n_cols: usize, x: &[f64]) -> Vec<f64> {
    (0..n_rows)
        .map(|i| dot(&a[i * n_cols..(i + 1) * n_cols], x))
        .collect()
}

/// Dense row-major `n_rows × n_cols` copy of a (possibly ragged or short) row list.
///
/// Missing rows and entries are zero; entries outside the requested shape are ignored.
fn dense_block(rows: &[Vec<f64>], n_rows: usize, n_cols: usize) -> Vec<f64> {
    let mut out = vec![0.0; n_rows * n_cols];
    for (i, row) in rows.iter().take(n_rows).enumerate() {
        for (j, &v) in row.iter().take(n_cols).enumerate() {
            out[i * n_cols + j] = v;
        }
    }
    out
}

/// Copy of `v` zero-padded or truncated to length `n`.
fn padded(v: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    let m = v.len().min(n);
    out[..m].copy_from_slice(&v[..m]);
    out
}

/// LU factorisation with scaled partial pivoting of a dense row-major `n × n` matrix.
struct DenseLu {
    n: usize,
    /// Unit lower factor below the diagonal, upper factor on and above it.
    lu: Vec<f64>,
    /// Row `i` of the factorised matrix is row `perm[i]` of the original one.
    perm: Vec<usize>,
}

impl DenseLu {
    /// Factorise `a`.  Returns `None` if `a` has non-finite entries, a zero row or
    /// is (numerically) singular, i.e. a row-scaled pivot drops below `n·ε`.
    fn factor(mut lu: Vec<f64>, n: usize) -> Option<Self> {
        if lu.len() != n * n || lu.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let mut row_scale = vec![0.0; n];
        for (i, scale) in row_scale.iter_mut().enumerate() {
            let big = lu[i * n..(i + 1) * n]
                .iter()
                .fold(0.0_f64, |m, v| m.max(v.abs()));
            if big == 0.0 {
                return None;
            }
            *scale = 1.0 / big;
        }
        let mut perm: Vec<usize> = (0..n).collect();
        let singular_threshold = n as f64 * f64::EPSILON;
        for k in 0..n {
            let mut pivot_row = k;
            let mut pivot_size = 0.0;
            for (i, &scale) in row_scale.iter().enumerate().skip(k) {
                let size = scale * lu[i * n + k].abs();
                if size > pivot_size {
                    pivot_size = size;
                    pivot_row = i;
                }
            }
            if pivot_size <= singular_threshold {
                return None;
            }
            if pivot_row != k {
                for j in 0..n {
                    lu.swap(k * n + j, pivot_row * n + j);
                }
                row_scale.swap(k, pivot_row);
                perm.swap(k, pivot_row);
            }
            let pivot = lu[k * n + k];
            for i in (k + 1)..n {
                let factor = lu[i * n + k] / pivot;
                lu[i * n + k] = factor;
                if factor != 0.0 {
                    for j in (k + 1)..n {
                        let upper = lu[k * n + j];
                        lu[i * n + j] -= factor * upper;
                    }
                }
            }
        }
        Some(Self { n, lu, perm })
    }

    /// Solve `A x = b` with the stored factors (`b` is zero-padded / truncated to `n`).
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let b = padded(b, n);
        let mut x: Vec<f64> = self.perm.iter().map(|&p| b[p]).collect();
        for i in 0..n {
            let row = &self.lu[i * n..i * n + i];
            let sum: f64 = row.iter().zip(x.iter()).map(|(l, xj)| l * xj).sum();
            x[i] -= sum;
        }
        for i in (0..n).rev() {
            let row = &self.lu[i * n + i + 1..(i + 1) * n];
            let sum: f64 = row
                .iter()
                .zip(x[i + 1..].iter())
                .map(|(u, xj)| u * xj)
                .sum();
            x[i] = (x[i] - sum) / self.lu[i * n + i];
        }
        x
    }
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
/// Manages the coupling between a displacement field `u` and a scalar field `s`
/// (temperature, pressure, or concentration) governed by the block system
///
/// ```text
/// K_u u + C_us s = f_u
/// C_su u + K_s s = f_s
/// ```
///
/// and applies the chosen solution scheme (with Aitken relaxation of the scalar
/// field for the iterative staggered scheme).  All three schemes solve this same
/// system; the reported residuals are the true block residuals of the returned
/// fields, so `converged` reflects the actual accuracy of the solution.
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

/// The assembled dense block system `[K_u C_us; C_su K_s] [u; s] = [f_u; f_s]`.
struct BlockSystem {
    n_u: usize,
    n_s: usize,
    k_u: Vec<f64>,
    k_s: Vec<f64>,
    c_us: Vec<f64>,
    c_su: Vec<f64>,
    f_u: Vec<f64>,
    f_s: Vec<f64>,
}

/// `‖r‖` relative to `scale`, with a non-finite result mapped to infinity.
fn relative_residual(r: &[f64], scale: f64) -> f64 {
    let r_norm = norm(r);
    let value = if scale > 0.0 { r_norm / scale } else { r_norm };
    if value.is_finite() {
        value
    } else {
        f64::INFINITY
    }
}

impl BlockSystem {
    /// Relative residuals `(res_u, res_s)` of the fields `(u, s)`.
    ///
    /// `r_u = f_u − K_u u − C_us s` and `r_s = f_s − C_su u − K_s s` are measured
    /// relative to the magnitude of the terms of their own equation,
    /// `‖f‖ + ‖K x‖ + ‖C x_other‖`, so that a block whose load vanishes (for
    /// example a purely thermally driven displacement field) is still judged
    /// meaningfully.
    fn residuals(&self, u: &[f64], s: &[f64]) -> (f64, f64) {
        let ku_u = mat_vec(&self.k_u, self.n_u, self.n_u, u);
        let cus_s = mat_vec(&self.c_us, self.n_u, self.n_s, s);
        let r_u: Vec<f64> = (0..self.n_u)
            .map(|i| self.f_u[i] - ku_u[i] - cus_s[i])
            .collect();
        let scale_u = norm(&self.f_u) + norm(&ku_u) + norm(&cus_s);

        let csu_u = mat_vec(&self.c_su, self.n_s, self.n_u, u);
        let ks_s = mat_vec(&self.k_s, self.n_s, self.n_s, s);
        let r_s: Vec<f64> = (0..self.n_s)
            .map(|i| self.f_s[i] - csu_u[i] - ks_s[i])
            .collect();
        let scale_s = norm(&self.f_s) + norm(&csu_u) + norm(&ks_s);

        (
            relative_residual(&r_u, scale_u),
            relative_residual(&r_s, scale_s),
        )
    }

    /// `f_u − C_us s`.
    fn displacement_rhs(&self, s: &[f64]) -> Vec<f64> {
        let cus_s = mat_vec(&self.c_us, self.n_u, self.n_s, s);
        self.f_u
            .iter()
            .zip(cus_s.iter())
            .map(|(f, c)| f - c)
            .collect()
    }

    /// `f_s − C_su u`.
    fn scalar_rhs(&self, u: &[f64]) -> Vec<f64> {
        let csu_u = mat_vec(&self.c_su, self.n_s, self.n_u, u);
        self.f_s
            .iter()
            .zip(csu_u.iter())
            .map(|(f, c)| f - c)
            .collect()
    }
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
    /// side for the scalar field `f_s`.  The stiffness-like blocks `k_u`
    /// (`n_disp_dof × n_disp_dof`), `k_s` (`n_scalar_dof × n_scalar_dof`) and the
    /// coupling blocks `c_us` (`n_disp_dof × n_scalar_dof`), `c_su`
    /// (`n_scalar_dof × n_disp_dof`) are dense row-major matrices; missing rows,
    /// entries or load components are taken as zero.  The system solved is
    /// `K_u u + C_us s = f_u`, `C_su u + K_s s = f_s`.
    ///
    /// Returns a [`CoupledResult`] holding the fields, the true relative block
    /// residuals `res_disp`, `res_scalar` and `converged`, which is set only if
    /// both residuals are below the tolerances of [`ConvergenceCriteria`].  A
    /// singular system yields `converged == false` with infinite residuals
    /// (see [`CoupledSolver::try_solve`] for a typed error instead).
    pub fn solve(
        &mut self,
        k_u: &[Vec<f64>],
        k_s: &[Vec<f64>],
        c_us: &[Vec<f64>],
        c_su: &[Vec<f64>],
        f_u: &[f64],
        f_s: &[f64],
    ) -> CoupledResult {
        match self.try_solve(k_u, k_s, c_us, c_su, f_u, f_s) {
            Ok(result) => result,
            Err(_) => CoupledResult {
                displacement: vec![0.0; self.n_disp_dof],
                scalar_field: vec![0.0; self.n_scalar_dof],
                n_iter: 0,
                res_disp: f64::INFINITY,
                res_scalar: f64::INFINITY,
                converged: false,
            },
        }
    }

    /// Fallible variant of [`CoupledSolver::solve`].
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::SingularSystem`] if the matrix that has to be
    /// inverted by the chosen scheme (the full block matrix for
    /// [`SolverScheme::Monolithic`], `K_u` and `K_s` for the staggered schemes)
    /// is singular or contains non-finite entries.
    pub fn try_solve(
        &mut self,
        k_u: &[Vec<f64>],
        k_s: &[Vec<f64>],
        c_us: &[Vec<f64>],
        c_su: &[Vec<f64>],
        f_u: &[f64],
        f_s: &[f64],
    ) -> crate::Result<CoupledResult> {
        let (n_u, n_s) = (self.n_disp_dof, self.n_scalar_dof);
        let system = BlockSystem {
            n_u,
            n_s,
            k_u: dense_block(k_u, n_u, n_u),
            k_s: dense_block(k_s, n_s, n_s),
            c_us: dense_block(c_us, n_u, n_s),
            c_su: dense_block(c_su, n_s, n_u),
            f_u: padded(f_u, n_u),
            f_s: padded(f_s, n_s),
        };
        match self.scheme {
            SolverScheme::Monolithic => self.solve_monolithic(&system),
            SolverScheme::Staggered => self.solve_staggered_one_pass(&system),
            SolverScheme::IterativeStaggered => self.solve_iterative_staggered(&system),
        }
    }

    /// Build the result of a scheme from its fields, measuring the true residuals.
    fn finish(
        &self,
        system: &BlockSystem,
        u: Vec<f64>,
        s: Vec<f64>,
        n_iter: usize,
    ) -> CoupledResult {
        let (res_disp, res_scalar) = system.residuals(&u, &s);
        let converged = res_disp < self.criteria.tol_disp && res_scalar < self.criteria.tol_scalar;
        CoupledResult {
            displacement: u,
            scalar_field: s,
            n_iter,
            res_disp,
            res_scalar,
            converged,
        }
    }

    /// Monolithic scheme: factorise the full block matrix once and solve directly.
    ///
    /// Gaussian elimination with scaled partial pivoting is used because the block
    /// matrix is generally nonsymmetric (`C_us ≠ C_suᵀ`).
    fn solve_monolithic(&self, system: &BlockSystem) -> crate::Result<CoupledResult> {
        let (n_u, n_s) = (system.n_u, system.n_s);
        let n_tot = n_u + n_s;
        let mut a = vec![0.0; n_tot * n_tot];
        for i in 0..n_u {
            for j in 0..n_u {
                a[i * n_tot + j] = system.k_u[i * n_u + j];
            }
            for j in 0..n_s {
                a[i * n_tot + n_u + j] = system.c_us[i * n_s + j];
            }
        }
        for i in 0..n_s {
            for j in 0..n_u {
                a[(n_u + i) * n_tot + j] = system.c_su[i * n_u + j];
            }
            for j in 0..n_s {
                a[(n_u + i) * n_tot + n_u + j] = system.k_s[i * n_s + j];
            }
        }
        let lu = DenseLu::factor(a, n_tot).ok_or(crate::Error::SingularSystem)?;
        let rhs: Vec<f64> = system
            .f_u
            .iter()
            .chain(system.f_s.iter())
            .copied()
            .collect();
        let x = lu.solve(&rhs);
        Ok(self.finish(system, x[..n_u].to_vec(), x[n_u..].to_vec(), 1))
    }

    /// One-pass staggered scheme: the scalar field is solved first with the
    /// displacement still at zero (the `C_su u` feedback is neglected), then the
    /// displacement is solved with the resulting scalar field.
    ///
    /// The pass is exact only for one-way coupling (`C_su = 0`); otherwise the
    /// reported residual of the scalar equation exposes the neglected feedback and
    /// `converged` is `false`.
    fn solve_staggered_one_pass(&self, system: &BlockSystem) -> crate::Result<CoupledResult> {
        let lu_u =
            DenseLu::factor(system.k_u.clone(), system.n_u).ok_or(crate::Error::SingularSystem)?;
        let lu_s =
            DenseLu::factor(system.k_s.clone(), system.n_s).ok_or(crate::Error::SingularSystem)?;
        let s = lu_s.solve(&system.f_s);
        let u = lu_u.solve(&system.displacement_rhs(&s));
        Ok(self.finish(system, u, s, 1))
    }

    /// Iteratively staggered (block Gauss–Seidel) scheme with Aitken relaxation
    /// of the scalar field.
    ///
    /// Each sweep solves the scalar equation with the displacement of the previous
    /// sweep, relaxes the scalar update with the Aitken factor ω, and then solves
    /// the displacement equation with the relaxed scalar field.  The iteration stops
    /// when the true block residuals fall below the tolerances or `max_iter` sweeps
    /// have been made.
    fn solve_iterative_staggered(&mut self, system: &BlockSystem) -> crate::Result<CoupledResult> {
        let lu_u =
            DenseLu::factor(system.k_u.clone(), system.n_u).ok_or(crate::Error::SingularSystem)?;
        let lu_s =
            DenseLu::factor(system.k_s.clone(), system.n_s).ok_or(crate::Error::SingularSystem)?;
        self.relaxation.reset(1.0);
        let mut u = vec![0.0; system.n_u];
        let mut s = vec![0.0; system.n_s];
        let mut n_iter = 0;
        for iter in 0..self.criteria.max_iter {
            let s_new = lu_s.solve(&system.scalar_rhs(&u));
            let ds: Vec<f64> = s_new.iter().zip(s.iter()).map(|(n, o)| n - o).collect();
            let omega = self.relaxation.update(&ds);
            for (s_i, d_i) in s.iter_mut().zip(ds.iter()) {
                *s_i += omega * d_i;
            }
            u = lu_u.solve(&system.displacement_rhs(&s));
            n_iter = iter + 1;
            let (res_disp, res_scalar) = system.residuals(&u, &s);
            if res_disp < self.criteria.tol_disp && res_scalar < self.criteria.tol_scalar {
                break;
            }
        }
        Ok(self.finish(system, u, s, n_iter))
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
        for (i, row) in k.iter_mut().enumerate() {
            row[i] = val;
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
        for (i, row) in c.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                assert!((v - c[j][i]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_tm_stiffness_positive_definite_diagonal() {
        let c = steel_tm().stiffness_tensor();
        for (i, row) in c.iter().enumerate() {
            assert!(row[i] > 0.0, "C[{i}][{i}] should be positive");
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
        for (&s, &sc) in sigma.iter().zip(sigma_check.iter()) {
            assert!((s - sc).abs() < 1e3);
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
        for (&eff, &eff_back) in sigma_eff.iter().zip(sigma_eff_back.iter()) {
            assert!((eff - eff_back).abs() < 1e3);
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
        assert!((ar.omega_min..=ar.omega_max).contains(&omega));
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

#[cfg(test)]
mod tests_coupled_solver {
    use super::*;

    /// Deterministic pseudo-random numbers in [0, 1).
    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn scaled_identity(n: usize, value: f64) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| (0..n).map(|j| if i == j { value } else { 0.0 }).collect())
            .collect()
    }

    fn constant_block(rows: usize, cols: usize, value: f64) -> Vec<Vec<f64>> {
        vec![vec![value; cols]; rows]
    }

    /// Diagonally dominant random matrix (nonsymmetric).
    fn random_dominant(n: usize, state: &mut u64) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            8.0 + 2.0 * lcg(state)
                        } else {
                            2.0 * lcg(state) - 1.0
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn random_block(rows: usize, cols: usize, amplitude: f64, state: &mut u64) -> Vec<Vec<f64>> {
        (0..rows)
            .map(|_| {
                (0..cols)
                    .map(|_| amplitude * (2.0 * lcg(state) - 1.0))
                    .collect()
            })
            .collect()
    }

    fn random_vector(n: usize, state: &mut u64) -> Vec<f64> {
        (0..n).map(|_| 2.0 * lcg(state) - 1.0).collect()
    }

    fn row_dot(row: &[f64], x: &[f64]) -> f64 {
        row.iter().zip(x).map(|(a, b)| a * b).sum()
    }

    /// A complete coupled problem K_u u + C_us s = f_u, C_su u + K_s s = f_s.
    struct Problem {
        k_u: Vec<Vec<f64>>,
        k_s: Vec<Vec<f64>>,
        c_us: Vec<Vec<f64>>,
        c_su: Vec<Vec<f64>>,
        f_u: Vec<f64>,
        f_s: Vec<f64>,
    }

    impl Problem {
        /// Random diagonally dominant problem with `n_u` + `n_s` unknowns.
        fn random(n_u: usize, n_s: usize, coupling: f64, one_way: bool, seed: u64) -> Self {
            let mut state = seed;
            let k_u = random_dominant(n_u, &mut state);
            let k_s = random_dominant(n_s, &mut state);
            let c_us = random_block(n_u, n_s, coupling, &mut state);
            let c_su = if one_way {
                constant_block(n_s, n_u, 0.0)
            } else {
                random_block(n_s, n_u, coupling, &mut state)
            };
            let f_u = random_vector(n_u, &mut state);
            let f_s = random_vector(n_s, &mut state);
            Self {
                k_u,
                k_s,
                c_us,
                c_su,
                f_u,
                f_s,
            }
        }

        fn sizes(&self) -> (usize, usize) {
            (self.f_u.len(), self.f_s.len())
        }

        fn solver(&self, scheme: SolverScheme, criteria: ConvergenceCriteria) -> CoupledSolver {
            let (n_u, n_s) = self.sizes();
            CoupledSolver::new(n_u, n_s, scheme, criteria)
        }

        fn solve(&self, scheme: SolverScheme, criteria: ConvergenceCriteria) -> CoupledResult {
            self.solver(scheme, criteria).solve(
                &self.k_u, &self.k_s, &self.c_us, &self.c_su, &self.f_u, &self.f_s,
            )
        }

        fn try_solve(
            &self,
            scheme: SolverScheme,
            criteria: ConvergenceCriteria,
        ) -> crate::Result<CoupledResult> {
            self.solver(scheme, criteria).try_solve(
                &self.k_u, &self.k_s, &self.c_us, &self.c_su, &self.f_u, &self.f_s,
            )
        }

        /// Full block matrix [[K_u, C_us], [C_su, K_s]].
        fn full_matrix(&self) -> Vec<Vec<f64>> {
            let mut rows: Vec<Vec<f64>> = Vec::new();
            for (ku_row, cus_row) in self.k_u.iter().zip(self.c_us.iter()) {
                rows.push(ku_row.iter().chain(cus_row.iter()).copied().collect());
            }
            for (csu_row, ks_row) in self.c_su.iter().zip(self.k_s.iter()) {
                rows.push(csu_row.iter().chain(ks_row.iter()).copied().collect());
            }
            rows
        }

        /// Solution of the full system by plain Gaussian elimination (independent of `DenseLu`).
        fn direct_solution(&self) -> Vec<f64> {
            let mut a = self.full_matrix();
            let mut b: Vec<f64> = self.f_u.iter().chain(self.f_s.iter()).copied().collect();
            let n = b.len();
            for k in 0..n {
                let mut p = k;
                for (i, row) in a.iter().enumerate().skip(k + 1) {
                    if row[k].abs() > a[p][k].abs() {
                        p = i;
                    }
                }
                a.swap(k, p);
                b.swap(k, p);
                let pivot_row = a[k].clone();
                let pivot_rhs = b[k];
                for (row, rhs) in a.iter_mut().zip(b.iter_mut()).skip(k + 1) {
                    let f = row[k] / pivot_row[k];
                    for (entry, pivot_entry) in row.iter_mut().zip(pivot_row.iter()).skip(k) {
                        *entry -= f * pivot_entry;
                    }
                    *rhs -= f * pivot_rhs;
                }
            }
            let mut x = vec![0.0; n];
            for i in (0..n).rev() {
                let tail: f64 = a[i][i + 1..]
                    .iter()
                    .zip(x[i + 1..].iter())
                    .map(|(aij, xj)| aij * xj)
                    .sum();
                x[i] = (b[i] - tail) / a[i][i];
            }
            x
        }

        /// Independent absolute block residual norms for the fields `(u, s)`.
        fn residual_norms(&self, u: &[f64], s: &[f64]) -> (f64, f64) {
            let ru: f64 = (0..self.f_u.len())
                .map(|i| self.f_u[i] - row_dot(&self.k_u[i], u) - row_dot(&self.c_us[i], s))
                .map(|v| v * v)
                .sum();
            let rs: f64 = (0..self.f_s.len())
                .map(|i| self.f_s[i] - row_dot(&self.c_su[i], u) - row_dot(&self.k_s[i], s))
                .map(|v| v * v)
                .sum();
            (ru.sqrt(), rs.sqrt())
        }
    }

    fn all_schemes() -> [SolverScheme; 3] {
        [
            SolverScheme::Monolithic,
            SolverScheme::Staggered,
            SolverScheme::IterativeStaggered,
        ]
    }

    fn uniform_problem(n: usize, k: f64, c_us: f64, c_su: f64, f_u: f64, f_s: f64) -> Problem {
        Problem {
            k_u: scaled_identity(n, k),
            k_s: scaled_identity(n, k),
            c_us: scaled_identity(n, c_us),
            c_su: scaled_identity(n, c_su),
            f_u: vec![f_u; n],
            f_s: vec![f_s; n],
        }
    }

    #[test]
    fn test_coupled_sign_convention_is_k_u_plus_c_s_equals_f() {
        // K_u = K_s = 2, C_us = 0.5, one-way coupling: s = 1/2, u = (1 - 0.5 * s) / 2 = 0.375.
        let problem = uniform_problem(1, 2.0, 0.5, 0.0, 1.0, 1.0);
        let r = problem.solve(SolverScheme::Monolithic, ConvergenceCriteria::default());
        assert!(r.converged);
        assert!((r.scalar_field[0] - 0.5).abs() < 1e-14);
        assert!((r.displacement[0] - 0.375).abs() < 1e-14);
    }

    #[test]
    fn test_coupled_schemes_agree_for_symmetric_coupling() {
        // Regression: K = 2 I, C = 0.5 I, f = 1 gave 0.4 (monolithic) vs 0.667 (iterative).
        let n = 3;
        let problem = uniform_problem(n, 2.0, 0.5, 0.5, 1.0, 1.0);
        let criteria = ConvergenceCriteria::new(1e-12, 1e-12, 200);
        let r_mono = problem.solve(SolverScheme::Monolithic, criteria.clone());
        let r_iter = problem.solve(SolverScheme::IterativeStaggered, criteria);
        assert!(r_mono.converged && r_iter.converged);
        for &v in r_mono.displacement.iter().chain(r_mono.scalar_field.iter()) {
            assert!((v - 0.4).abs() < 1e-12, "monolithic value {v}");
        }
        for &v in r_iter.displacement.iter().chain(r_iter.scalar_field.iter()) {
            assert!((v - 0.4).abs() < 1e-9, "iterative value {v}");
        }
    }

    #[test]
    fn test_one_pass_staggered_with_feedback_is_not_reported_converged() {
        let problem = uniform_problem(2, 2.0, 0.5, 0.5, 1.0, 1.0);
        let r = problem.solve(SolverScheme::Staggered, ConvergenceCriteria::default());
        // s = 0.5, u = 0.375: the scalar equation is violated by C_su u (the neglected feedback).
        assert!((r.scalar_field[0] - 0.5).abs() < 1e-14);
        assert!((r.displacement[0] - 0.375).abs() < 1e-14);
        assert!(
            !r.converged,
            "feedback was neglected: must not be converged"
        );
        assert_eq!(r.n_iter, 1);
        assert!(r.res_scalar > 0.05, "res_scalar = {}", r.res_scalar);
        assert!(r.res_disp < 1e-14, "displacement equation is satisfied");
    }

    #[test]
    fn test_one_pass_staggered_is_exact_for_one_way_coupling() {
        let problem = Problem::random(4, 3, 0.5, true, 7);
        let r = problem.solve(SolverScheme::Staggered, ConvergenceCriteria::default());
        assert!(r.converged, "residuals ({}, {})", r.res_disp, r.res_scalar);
        let reference = problem.direct_solution();
        for (i, &u) in r.displacement.iter().enumerate() {
            assert!((u - reference[i]).abs() < 1e-10);
        }
        for (i, &s) in r.scalar_field.iter().enumerate() {
            assert!((s - reference[4 + i]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_nonsymmetric_system_all_schemes_match_direct_solve() {
        let problem = Problem::random(4, 3, 0.6, false, 0xdead_beef);
        // The full block matrix really is nonsymmetric.
        let a = problem.full_matrix();
        let asym = (0..a.len())
            .flat_map(|i| (0..a.len()).map(move |j| (i, j)))
            .map(|(i, j)| (a[i][j] - a[j][i]).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            asym > 0.1,
            "test matrix must be nonsymmetric, asymmetry {asym}"
        );
        let reference = problem.direct_solution();
        let criteria = ConvergenceCriteria::new(1e-13, 1e-13, 200);
        for scheme in [SolverScheme::Monolithic, SolverScheme::IterativeStaggered] {
            let r = problem.solve(scheme, criteria.clone());
            assert!(
                r.converged,
                "{scheme:?}: residuals ({:e}, {:e}) after {} iterations",
                r.res_disp, r.res_scalar, r.n_iter
            );
            for (i, &u) in r.displacement.iter().enumerate() {
                assert!(
                    (u - reference[i]).abs() < 1e-8,
                    "{scheme:?}: u[{i}] = {u} vs direct {}",
                    reference[i]
                );
            }
            for (i, &s) in r.scalar_field.iter().enumerate() {
                assert!(
                    (s - reference[4 + i]).abs() < 1e-8,
                    "{scheme:?}: s[{i}] = {s} vs direct {}",
                    reference[4 + i]
                );
            }
            let (ru, rs) = problem.residual_norms(&r.displacement, &r.scalar_field);
            assert!(
                ru < 1e-10 && rs < 1e-10,
                "independent residuals ({ru:e}, {rs:e})"
            );
        }
        // The iterative scheme needs real sweeps (more than one) for genuinely coupled blocks.
        let r = problem.solve(SolverScheme::IterativeStaggered, criteria);
        assert!(r.n_iter > 1, "n_iter = {}", r.n_iter);
    }

    #[test]
    fn test_reported_residuals_are_the_true_block_residuals() {
        // A one-pass solve of a strongly coupled system leaves a real residual; res_* must reflect it.
        let problem = Problem::random(3, 2, 1.5, false, 99);
        let r = problem.solve(SolverScheme::Staggered, ConvergenceCriteria::default());
        let (ru, rs) = problem.residual_norms(&r.displacement, &r.scalar_field);
        assert!(rs > 1e-3, "independent scalar residual {rs}");
        assert!(ru < 1e-12, "displacement equation is solved exactly: {ru}");
        assert!(
            r.res_scalar > 1e-4,
            "reported res_scalar = {}",
            r.res_scalar
        );
        assert!(r.res_disp < 1e-12, "reported res_disp = {}", r.res_disp);
        assert!(!r.converged);
    }

    #[test]
    fn test_thermally_driven_displacement_with_zero_mechanical_load() {
        // f_u = 0: the displacement is driven only by the coupling term; residuals must stay finite
        // and all schemes must agree.
        let n = 2;
        let problem = Problem {
            k_u: scaled_identity(n, 4.0),
            k_s: scaled_identity(n, 2.0),
            c_us: constant_block(n, n, 0.5),
            c_su: constant_block(n, n, 0.0),
            f_u: vec![0.0; n],
            f_s: vec![1.0; n],
        };
        for scheme in all_schemes() {
            let r = problem.solve(scheme, ConvergenceCriteria::default());
            assert!(r.converged, "{scheme:?}");
            // s = 1/2 each; K_u u = -C_us s = -0.5 (0.5 + 0.5) = -0.5  ->  u = -0.125.
            for &u in &r.displacement {
                assert!((u + 0.125).abs() < 1e-10, "{scheme:?}: u = {u}");
            }
            assert!(r.res_disp.is_finite() && r.res_scalar.is_finite());
        }
    }

    #[test]
    fn test_singular_system_is_reported_not_solved() {
        let n = 2;
        let problem = Problem {
            k_u: constant_block(n, n, 0.0),
            k_s: scaled_identity(n, 1.0),
            c_us: constant_block(n, n, 0.0),
            c_su: constant_block(n, n, 0.0),
            f_u: vec![1.0; n],
            f_s: vec![1.0; n],
        };
        for scheme in all_schemes() {
            let r = problem.solve(scheme, ConvergenceCriteria::default());
            assert!(!r.converged, "{scheme:?}: singular K_u cannot converge");
            assert!(r.res_disp.is_infinite() && r.res_scalar.is_infinite());
            assert_eq!(r.displacement.len(), n);
            assert_eq!(r.scalar_field.len(), n);
            let err = problem
                .try_solve(scheme, ConvergenceCriteria::default())
                .unwrap_err();
            assert!(
                matches!(err, crate::Error::SingularSystem),
                "{scheme:?}: {err}"
            );
        }
        // The block matrix [[1, 1], [1, 1]] is singular although K_u and K_s are regular.
        let coupled = uniform_problem(1, 1.0, 1.0, 1.0, 1.0, 1.0);
        assert!(
            coupled
                .try_solve(SolverScheme::Monolithic, ConvergenceCriteria::default())
                .is_err()
        );
    }

    #[test]
    fn test_iterative_staggered_reports_exhausted_iterations() {
        // Strong feedback needs many sweeps; two are not enough for a 1e-14 tolerance.
        let problem = uniform_problem(2, 2.0, 1.5, 1.5, 1.0, 1.0);
        let r = problem.solve(
            SolverScheme::IterativeStaggered,
            ConvergenceCriteria::new(1e-14, 1e-14, 2),
        );
        assert_eq!(r.n_iter, 2);
        assert!(!r.converged);
        assert!(r.res_disp.is_finite() && r.res_scalar.is_finite());
    }

    #[test]
    fn test_dense_lu_pivoting_scaling_and_singularity() {
        // Needs a row exchange.
        let lu = DenseLu::factor(vec![0.0, 1.0, 1.0, 0.0], 2).unwrap();
        let x = lu.solve(&[3.0, 5.0]);
        assert!((x[0] - 5.0).abs() < 1e-14 && (x[1] - 3.0).abs() < 1e-14);
        // Badly scaled but regular.
        let lu = DenseLu::factor(vec![1.0e9, 0.0, 0.0, 1.0e-9], 2).unwrap();
        let x = lu.solve(&[2.0e9, 3.0e-9]);
        assert!((x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12);
        // Rank deficient / zero row / non-finite.
        assert!(DenseLu::factor(vec![1.0, 2.0, 2.0, 4.0], 2).is_none());
        assert!(DenseLu::factor(vec![1.0, 2.0, 0.0, 0.0], 2).is_none());
        assert!(DenseLu::factor(vec![f64::NAN, 0.0, 0.0, 1.0], 2).is_none());
        // Empty system is trivially solvable.
        assert!(
            DenseLu::factor(Vec::new(), 0)
                .unwrap()
                .solve(&[])
                .is_empty()
        );
    }

    mod agreement;
}
