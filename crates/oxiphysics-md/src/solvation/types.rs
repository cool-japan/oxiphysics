//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{ELEM_CHARGE, EPS0, J_TO_KJMOL};

/// Thermodynamic integration (TI) estimator for solvation free energy.
///
/// TI integrates ⟨∂U/∂λ⟩_λ over the coupling parameter λ ∈ \[0, 1\]:
///
/// ΔG = ∫₀¹ ⟨∂U/∂λ⟩_λ dλ
///
/// This struct stores (λ, ⟨∂U/∂λ⟩) pairs and integrates numerically
/// using the trapezoidal rule.
#[derive(Debug, Clone)]
pub struct ThermodynamicIntegrationSolvation {
    /// Coupling parameters λ_i ∈ \[0, 1\].
    pub lambda_values: Vec<f64>,
    /// Ensemble averages ⟨∂U/∂λ⟩ at each λ_i (kJ/mol).
    pub dhdl_averages: Vec<f64>,
    /// Statistical uncertainties σ(⟨∂U/∂λ⟩) at each λ_i (kJ/mol).
    pub dhdl_errors: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
}
impl ThermodynamicIntegrationSolvation {
    /// Construct a TI solvation free energy estimator.
    ///
    /// # Arguments
    /// * `temperature` – simulation temperature (K)
    pub fn new(temperature: f64) -> Self {
        Self {
            lambda_values: Vec::new(),
            dhdl_averages: Vec::new(),
            dhdl_errors: Vec::new(),
            temperature,
        }
    }
    /// Add a data point (λ, ⟨∂U/∂λ⟩, σ).
    ///
    /// # Arguments
    /// * `lambda`    – coupling parameter
    /// * `dhdl_avg`  – ensemble average of ∂U/∂λ (kJ/mol)
    /// * `dhdl_err`  – standard error of ∂U/∂λ (kJ/mol)
    pub fn add_point(&mut self, lambda: f64, dhdl_avg: f64, dhdl_err: f64) {
        self.lambda_values.push(lambda);
        self.dhdl_averages.push(dhdl_avg);
        self.dhdl_errors.push(dhdl_err);
    }
    /// Compute ΔG_solv via the trapezoidal rule (kJ/mol).
    ///
    /// Returns 0 if fewer than 2 data points are present.
    pub fn integrate(&self) -> f64 {
        let n = self.lambda_values.len();
        if n < 2 {
            return 0.0;
        }
        let mut sum = 0.0_f64;
        for i in 0..n - 1 {
            let dl = self.lambda_values[i + 1] - self.lambda_values[i];
            let avg = 0.5 * (self.dhdl_averages[i] + self.dhdl_averages[i + 1]);
            sum += avg * dl;
        }
        sum
    }
    /// Propagated uncertainty on ΔG via the trapezoidal rule (kJ/mol).
    ///
    /// σ²(ΔG) = Σ_i (δλ/2)² (σ_i² + σ_{i+1}²)
    pub fn uncertainty(&self) -> f64 {
        let n = self.lambda_values.len();
        if n < 2 {
            return f64::INFINITY;
        }
        let mut var = 0.0_f64;
        for i in 0..n - 1 {
            let dl = self.lambda_values[i + 1] - self.lambda_values[i];
            var += (0.5 * dl).powi(2)
                * (self.dhdl_errors[i].powi(2) + self.dhdl_errors[i + 1].powi(2));
        }
        var.sqrt()
    }
    /// Compute ΔG using Simpson's rule (requires odd number of points).
    ///
    /// Falls back to trapezoidal if point count is even.
    pub fn integrate_simpson(&self) -> f64 {
        let n = self.lambda_values.len();
        if n < 3 || n.is_multiple_of(2) {
            return self.integrate();
        }
        let mut sum = 0.0_f64;
        let i_max = n - 2;
        let mut i = 0;
        while i < i_max {
            let dl = self.lambda_values[i + 2] - self.lambda_values[i];
            sum += (dl / 6.0)
                * (self.dhdl_averages[i]
                    + 4.0 * self.dhdl_averages[i + 1]
                    + self.dhdl_averages[i + 2]);
            i += 2;
        }
        sum
    }
    /// Free energy difference between two λ windows (partial integration).
    ///
    /// # Arguments
    /// * `lam_a` – lower λ bound
    /// * `lam_b` – upper λ bound
    pub fn partial_free_energy(&self, lam_a: f64, lam_b: f64) -> f64 {
        let n = self.lambda_values.len();
        let mut sum = 0.0_f64;
        for i in 0..n.saturating_sub(1) {
            let la = self.lambda_values[i];
            let lb = self.lambda_values[i + 1];
            if lb <= lam_a || la >= lam_b {
                continue;
            }
            let dl = lb - la;
            let avg = 0.5 * (self.dhdl_averages[i] + self.dhdl_averages[i + 1]);
            sum += avg * dl;
        }
        sum
    }
}
/// COSMO-like continuum solvation model.
///
/// Implements a simplified COSMO (COnductor-like Screening MOdel) approach:
/// - Dielectric screening energy: ΔG_el = f(ε) · q² / (8π ε₀ r)
/// - f(ε) = (ε − 1) / (ε + 1/2) (COSMO factor)
/// - Cavity surface area with a dielectric penalty
#[derive(Debug, Clone)]
pub struct ContinuumSolvation {
    /// Solvent relative dielectric constant ε_r.
    pub dielectric: f64,
    /// Solute radius (Å).
    pub radius: f64,
    /// Surface tension coefficient γ (kJ/mol Å⁻²).
    pub gamma: f64,
    /// Ion charge (elementary charge units).
    pub charge: f64,
}
impl ContinuumSolvation {
    /// Create a continuum solvation model.
    ///
    /// # Arguments
    /// * `eps`    – solvent dielectric
    /// * `r`      – solute radius (Å)
    /// * `gamma`  – surface tension (kJ/mol Å⁻²)
    /// * `charge` – solute charge (e)
    pub fn new(eps: f64, r: f64, gamma: f64, charge: f64) -> Self {
        Self {
            dielectric: eps,
            radius: r,
            gamma,
            charge,
        }
    }
    /// COSMO dielectric scaling factor f(ε) = (ε−1)/(ε+1/2).
    pub fn cosmo_factor(&self) -> f64 {
        if self.dielectric <= 0.0 {
            return 0.0;
        }
        (self.dielectric - 1.0) / (self.dielectric + 0.5)
    }
    /// COSMO electrostatic solvation energy (kJ/mol).
    ///
    /// ΔG_el = −f(ε) · (q e)² N_A / (8π ε₀ r)
    pub fn electrostatic_energy(&self) -> f64 {
        let r_m = self.radius * 1e-10;
        if r_m <= 0.0 {
            return 0.0;
        }
        let q_c = self.charge * ELEM_CHARGE;
        let prefactor = q_c * q_c * J_TO_KJMOL / (8.0 * PI * EPS0 * r_m);
        -self.cosmo_factor() * prefactor
    }
    /// Cavity surface energy ΔG_cav = γ · 4π r² (kJ/mol).
    pub fn cavity_energy(&self) -> f64 {
        let sasa = 4.0 * PI * self.radius * self.radius;
        self.gamma * sasa
    }
    /// Total solvation energy ΔG = ΔG_el + ΔG_cav (kJ/mol).
    pub fn total_energy(&self) -> f64 {
        self.electrostatic_energy() + self.cavity_energy()
    }
    /// Effective dielectric at distance r from solute using sigmoidal switching.
    ///
    /// Switches from 1 (inside) to ε_r (bulk) over the surface.
    ///
    /// # Arguments
    /// * `r` – distance from solute centre (Å)
    pub fn effective_dielectric(&self, r: f64) -> f64 {
        let s = r - self.radius;
        1.0 + (self.dielectric - 1.0) / (1.0 + (-4.0 * s).exp())
    }
}
/// Solvation free energy from test particle insertion (Widom method).
///
/// The excess chemical potential is:
/// ```text
/// μ_ex = −k_B T ln⟨exp(−ΔU / k_B T)⟩
/// ```
/// where ΔU is the interaction energy of the inserted test particle.
#[derive(Debug, Clone)]
pub struct SolvationMD {
    /// Temperature (K).
    pub temperature: f64,
    /// Collected Boltzmann factors exp(−ΔU / k_B T) from trial insertions.
    pub boltzmann_factors: Vec<f64>,
    /// Lennard-Jones ε for the test particle (kJ/mol).
    pub eps_lj: f64,
    /// Lennard-Jones σ for the test particle (Å).
    pub sigma_lj: f64,
}
impl SolvationMD {
    /// Create a new Widom test particle insertion sampler.
    ///
    /// # Arguments
    /// * `t`        – temperature (K)
    /// * `eps_lj`   – LJ epsilon (kJ/mol)
    /// * `sigma_lj` – LJ sigma (Å)
    pub fn new(t: f64, eps_lj: f64, sigma_lj: f64) -> Self {
        Self {
            temperature: t,
            boltzmann_factors: Vec::new(),
            eps_lj,
            sigma_lj,
        }
    }
    /// Add a Boltzmann factor from one test particle insertion.
    ///
    /// # Arguments
    /// * `delta_u` – insertion energy ΔU (kJ/mol)
    pub fn add_insertion(&mut self, delta_u: f64) {
        let kbt_kjmol = 8.314_462_618e-3 * self.temperature;
        let bf = (-delta_u / kbt_kjmol).exp();
        self.boltzmann_factors.push(bf);
    }
    /// Excess chemical potential μ_ex (kJ/mol) from Widom estimator.
    pub fn excess_chemical_potential(&self) -> f64 {
        if self.boltzmann_factors.is_empty() {
            return 0.0;
        }
        let mean_bf: f64 =
            self.boltzmann_factors.iter().sum::<f64>() / self.boltzmann_factors.len() as f64;
        if mean_bf <= 0.0 {
            return f64::INFINITY;
        }
        let kbt_kjmol = 8.314_462_618e-3 * self.temperature;
        -kbt_kjmol * mean_bf.ln()
    }
    /// Hydration free energy estimate using a simple LJ insertion at a given distance.
    ///
    /// Inserts the test particle at distance r from a single solvent molecule and
    /// accumulates the Boltzmann factor.
    ///
    /// # Arguments
    /// * `r` – insertion distance (Å)
    pub fn lj_insertion(&mut self, r: f64) {
        if r <= 0.0 {
            self.add_insertion(1e6_f64);
            return;
        }
        let sr = self.sigma_lj / r;
        let lj = 4.0 * self.eps_lj * (sr.powi(12) - sr.powi(6));
        self.add_insertion(lj);
    }
    /// Number of trial insertions accumulated.
    pub fn n_insertions(&self) -> usize {
        self.boltzmann_factors.len()
    }
    /// Standard error of the mean Boltzmann factor.
    pub fn standard_error(&self) -> f64 {
        let n = self.boltzmann_factors.len();
        if n < 2 {
            return f64::INFINITY;
        }
        let mean: f64 = self.boltzmann_factors.iter().sum::<f64>() / n as f64;
        let var: f64 = self
            .boltzmann_factors
            .iter()
            .map(|&b| (b - mean).powi(2))
            .sum::<f64>()
            / (n - 1) as f64;
        (var / n as f64).sqrt()
    }
}
/// Combined solvation free energy (Born + cavity + non-polar).
pub struct SolvationFreeEnergy {
    /// Solute charge (elementary charge units).
    pub solute_charge: f64,
    /// Solvent relative dielectric constant ε_r.
    pub dielectric_epsilon: f64,
    /// Born radius (m).
    pub born_radius: f64,
}
impl SolvationFreeEnergy {
    /// Create a new solvation free energy calculator.
    pub fn new(solute_charge: f64, dielectric_epsilon: f64, born_radius: f64) -> Self {
        Self {
            solute_charge,
            dielectric_epsilon,
            born_radius,
        }
    }
    /// Compute the Born solvation free energy (kJ/mol).
    pub fn compute_born(&self, radius: f64, charge: f64, eps: f64) -> f64 {
        let r = if radius == 0.0 {
            self.born_radius
        } else {
            radius
        };
        let q = if charge == 0.0 {
            self.solute_charge
        } else {
            charge
        };
        let e = if eps == 0.0 {
            self.dielectric_epsilon
        } else {
            eps
        };
        born_energy(q, r, 1.0, e)
    }
    /// Total Born solvation free energy using stored parameters.
    pub fn total_born_energy(&self) -> f64 {
        born_energy(
            self.solute_charge,
            self.born_radius,
            1.0,
            self.dielectric_epsilon,
        )
    }
}
/// Geometric model of the solvation shell with RDF-derived properties.
///
/// Provides coordination numbers, shell radii from a radial distribution
/// function, and mean residence time estimates.
#[derive(Debug, Clone)]
pub struct SolvationShell {
    /// Number of solvent molecules in the first hydration shell.
    pub first_shell_count: usize,
    /// Number of solvent molecules in the second hydration shell.
    pub second_shell_count: usize,
    /// Radius of the first shell minimum in the RDF (Å).
    pub shell_radius: f64,
    /// Shell thickness Δr (Å).
    pub shell_thickness: f64,
    /// Radial distribution function g(r) values.
    pub rdf: Vec<f64>,
    /// Radial grid points r (Å).
    pub r_grid: Vec<f64>,
}
impl SolvationShell {
    /// Create a solvation shell with explicit shell parameters.
    ///
    /// # Arguments
    /// * `first`   – coordination number of first shell
    /// * `second`  – coordination number of second shell
    /// * `r_shell` – first-shell maximum radius (Å)
    /// * `dr`      – shell thickness (Å)
    pub fn new(first: usize, second: usize, r_shell: f64, dr: f64) -> Self {
        Self {
            first_shell_count: first,
            second_shell_count: second,
            shell_radius: r_shell,
            shell_thickness: dr,
            rdf: Vec::new(),
            r_grid: Vec::new(),
        }
    }
    /// Set the radial distribution function data.
    ///
    /// # Arguments
    /// * `r`    – radial distances (Å)
    /// * `gr`   – g(r) values
    pub fn set_rdf(&mut self, r: Vec<f64>, gr: Vec<f64>) {
        self.r_grid = r;
        self.rdf = gr;
    }
    /// Coordination number from RDF integration: n = 4π ρ ∫ g(r) r² dr.
    ///
    /// # Arguments
    /// * `rho`   – solvent number density (Å⁻³)
    /// * `r_min` – lower bound of integration (Å)
    /// * `r_max` – upper bound of integration (Å)
    pub fn coordination_number_from_rdf(&self, rho: f64, r_min: f64, r_max: f64) -> f64 {
        if self.r_grid.len() < 2 || self.rdf.len() != self.r_grid.len() {
            return 0.0;
        }
        let mut integral = 0.0_f64;
        for i in 0..self.r_grid.len() - 1 {
            let r = self.r_grid[i];
            let rp1 = self.r_grid[i + 1];
            if r >= r_max {
                break;
            }
            if rp1 < r_min {
                continue;
            }
            let dr = rp1 - r;
            let g_mid = 0.5 * (self.rdf[i] + self.rdf[i + 1]);
            let r_mid = 0.5 * (r + rp1);
            integral += g_mid * r_mid * r_mid * dr;
        }
        4.0 * PI * rho * integral
    }
    /// Total hydration number (first + second shell).
    pub fn total_hydration_number(&self) -> usize {
        self.first_shell_count + self.second_shell_count
    }
    /// Shell volume (Å³) for shell index 1 or 2.
    pub fn shell_volume(&self, shell_index: u8) -> f64 {
        let r_inner = match shell_index {
            1 => 2.0_f64,
            2 => 2.0 + self.shell_thickness,
            _ => return 0.0,
        };
        let r_outer = r_inner + self.shell_thickness;
        4.0 / 3.0 * PI * (r_outer.powi(3) - r_inner.powi(3))
    }
    /// Estimate mean residence time (ps) using a heuristic model.
    pub fn compute_residence_time(&self) -> f64 {
        let tau0 = 20.0_f64;
        let n1 = self.first_shell_count as f64;
        let n2 = self.second_shell_count as f64;
        let ratio = if n1 == 0.0 {
            1.0
        } else {
            (n2 / n1.max(1.0)).max(1.0)
        };
        tau0 * n1 / ratio
    }
    /// RDF peak value (maximum of g(r)).
    pub fn rdf_peak(&self) -> f64 {
        self.rdf.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Distance at which g(r) first exceeds a threshold.
    pub fn rdf_first_peak_position(&self, threshold: f64) -> f64 {
        for (i, &g) in self.rdf.iter().enumerate() {
            if g >= threshold {
                return self.r_grid.get(i).cloned().unwrap_or(0.0);
            }
        }
        0.0
    }
    /// First minimum of g(r) after the first peak (Å), using finite differences.
    ///
    /// The first minimum separates the first and second solvation shells.
    pub fn first_shell_boundary(&self) -> f64 {
        let n = self.rdf.len();
        if n < 3 {
            return self.shell_radius;
        }
        let mut found_peak = false;
        for i in 1..n - 1 {
            if !found_peak && self.rdf[i] > self.rdf[i - 1] && self.rdf[i] >= self.rdf[i + 1] {
                found_peak = true;
            }
            if found_peak && self.rdf[i] < self.rdf[i - 1] && self.rdf[i] <= self.rdf[i + 1] {
                return self.r_grid.get(i).cloned().unwrap_or(self.shell_radius);
            }
        }
        self.shell_radius
    }
    /// Running integral of g(r): cumulative coordination N(r) = 4π ρ ∫₀ʳ g(r') r'² dr'.
    ///
    /// # Arguments
    /// * `rho` – solvent number density (Å⁻³)
    ///
    /// Returns a vector of N values evaluated at each r_grid point.
    pub fn cumulative_coordination(&self, rho: f64) -> Vec<f64> {
        let n = self.r_grid.len();
        let mut result = vec![0.0_f64; n];
        let mut integral = 0.0_f64;
        for i in 0..n.saturating_sub(1) {
            let dr = self.r_grid[i + 1] - self.r_grid[i];
            let g_mid = 0.5 * (self.rdf[i] + self.rdf[i + 1]);
            let r_mid = 0.5 * (self.r_grid[i] + self.r_grid[i + 1]);
            integral += 4.0 * PI * rho * g_mid * r_mid * r_mid * dr;
            result[i + 1] = integral;
        }
        result
    }
}
/// Transfer free energy between two solvents with log P correlation.
pub struct TransferFreeEnergy {
    /// Source solvent name.
    pub from_solvent: String,
    /// Target solvent name.
    pub to_solvent: String,
    /// Per-fragment ΔΔG contributions.
    pub ddg: Vec<(String, f64)>,
}
impl TransferFreeEnergy {
    /// Create a transfer free energy container.
    pub fn new(from_solvent: impl Into<String>, to_solvent: impl Into<String>) -> Self {
        Self {
            from_solvent: from_solvent.into(),
            to_solvent: to_solvent.into(),
            ddg: Vec::new(),
        }
    }
    /// Add a fragment contribution.
    pub fn add_fragment(&mut self, name: impl Into<String>, ddg: f64) {
        self.ddg.push((name.into(), ddg));
    }
    /// Compute the total ΔΔG = Σ ddG_i (kJ/mol).
    pub fn total_ddg(&self) -> f64 {
        self.ddg.iter().map(|(_, v)| v).sum()
    }
    /// Estimate log P from the transfer free energy.
    pub fn estimate_log_p(&self) -> f64 {
        let rt = 2.479_f64;
        -self.total_ddg() / (rt * std::f64::consts::LN_2 * std::f64::consts::LOG2_E)
    }
    /// Pearson correlation of ddG values with supplied log P values.
    pub fn correlation_with_logp(&self, log_p_values: &[f64]) -> f64 {
        let n = self.ddg.len();
        if n < 2 || log_p_values.len() != n {
            return 0.0;
        }
        let x: Vec<f64> = self.ddg.iter().map(|(_, v)| *v).collect();
        let y = log_p_values;
        let mx = x.iter().sum::<f64>() / n as f64;
        let my = y.iter().sum::<f64>() / n as f64;
        let cov: f64 = x
            .iter()
            .zip(y.iter())
            .map(|(xi, yi)| (xi - mx) * (yi - my))
            .sum::<f64>();
        let sx: f64 = x.iter().map(|xi| (xi - mx).powi(2)).sum::<f64>().sqrt();
        let sy: f64 = y.iter().map(|yi| (yi - my).powi(2)).sum::<f64>().sqrt();
        if sx == 0.0 || sy == 0.0 {
            0.0
        } else {
            cov / (sx * sy)
        }
    }
}
/// Continuum solvent model (PB/SA).
pub struct ContinuumSolvent {
    /// Dielectric constant of the solute interior.
    pub dielectric_inner: f64,
    /// Dielectric constant of the solvent.
    pub dielectric_outer: f64,
    /// Ionic strength I (mol/m³).
    pub ionic_strength: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Surface tension coefficient for non-polar SASA term (kJ/mol/Å²).
    pub gamma_np: f64,
    /// Solvent-accessible surface area (Å²).
    pub sasa: f64,
}
impl ContinuumSolvent {
    /// Create a continuum solvent model.
    pub fn new(
        eps_in: f64,
        eps_out: f64,
        ionic_strength: f64,
        temperature: f64,
        gamma_np: f64,
        sasa: f64,
    ) -> Self {
        Self {
            dielectric_inner: eps_in,
            dielectric_outer: eps_out,
            ionic_strength,
            temperature,
            gamma_np,
            sasa,
        }
    }
    /// Compute the Debye screening length (m).
    pub fn debye_length(&self) -> f64 {
        debye_length(self.dielectric_outer, self.ionic_strength, self.temperature)
    }
    /// Compute the total PB/SA solvation energy (kJ/mol).
    pub fn compute_pbsa_energy(&self) -> f64 {
        self.gamma_np * self.sasa
    }
    /// Estimate the electrostatic screening factor relative to vacuum.
    pub fn dielectric_screening_factor(&self) -> f64 {
        if self.dielectric_outer == 0.0 {
            return 0.0;
        }
        self.dielectric_inner / self.dielectric_outer
    }
}
/// Born model electrostatic solvation energy for a single ion.
///
/// Models the ion as a charged sphere of radius R embedded in a continuum
/// solvent with relative permittivity ε_r.  The solvation energy is the work
/// of transferring the charged sphere from vacuum (ε = 1) into the solvent:
///
/// ΔG_Born = −(q e)² N_A (1 − 1/ε_r) / (8π ε₀ R)
///
/// This is equivalent to calling [`born_energy`] with `eps_in = 1`.
#[derive(Debug, Clone)]
pub struct BornSolvation {
    /// Ion charge (elementary charge units).
    pub charge: f64,
    /// Born radius (m).  Typically 1.2–3 Å for common ions.
    pub radius: f64,
    /// Solvent relative dielectric constant ε_r (water ≈ 78.5).
    pub dielectric: f64,
    /// Temperature (K).
    pub temperature: f64,
}
impl BornSolvation {
    /// Construct a new Born solvation model.
    ///
    /// # Arguments
    /// * `charge`      – ion charge in units of elementary charge
    /// * `radius`      – Born radius (m)
    /// * `dielectric`  – solvent dielectric constant ε_r
    /// * `temperature` – temperature (K)
    pub fn new(charge: f64, radius: f64, dielectric: f64, temperature: f64) -> Self {
        Self {
            charge,
            radius,
            dielectric,
            temperature,
        }
    }
    /// Electrostatic solvation free energy ΔG_Born (kJ/mol).
    ///
    /// Returns the Born solvation free energy for transfer from vacuum to
    /// solvent.  The result is negative for ions in high-ε solvents.
    pub fn solvation_energy(&self) -> f64 {
        born_energy(self.charge, self.radius, 1.0, self.dielectric)
    }
    /// Self-energy in vacuum: G_vac = q²N_A / (8π ε₀ R) (kJ/mol).
    pub fn self_energy_vacuum(&self) -> f64 {
        if self.radius <= 0.0 {
            return 0.0;
        }
        let q_c = self.charge * ELEM_CHARGE;
        q_c * q_c * J_TO_KJMOL / (8.0 * PI * EPS0 * self.radius)
    }
    /// Self-energy in solvent: G_sol = q²N_A / (8π ε₀ ε_r R) (kJ/mol).
    pub fn self_energy_solvent(&self) -> f64 {
        if self.radius <= 0.0 || self.dielectric <= 0.0 {
            return 0.0;
        }
        let q_c = self.charge * ELEM_CHARGE;
        q_c * q_c * J_TO_KJMOL / (8.0 * PI * EPS0 * self.dielectric * self.radius)
    }
    /// Entropy contribution −T ∂ΔG/∂T (kJ/mol), assuming ε_r(T) is constant.
    ///
    /// For a pure Born model with constant ε, the entropy term is zero.
    /// This method returns 0 as a placeholder.
    pub fn entropy_contribution(&self) -> f64 {
        0.0
    }
    /// Scale the solvation energy by a partial-charge Born radius product.
    ///
    /// ΔG_scaled = ΔG_Born · scale_factor
    ///
    /// Useful for parameterisation against experimental hydration data.
    pub fn scaled_energy(&self, scale_factor: f64) -> f64 {
        self.solvation_energy() * scale_factor
    }
}
/// Generalized Born (GB) model for multi-atom solutes.
///
/// Uses the Hawkins-Cramer-Truhlar (HCT) pairwise descreening approximation
/// to compute effective Born radii α_i for each atom, then evaluates the
/// GB electrostatic solvation energy:
///
/// ΔG_GB = −(1/2)(1 − 1/ε_r) Σ_{i,j} q_i q_j e² N_A / f_GB(r_ij, α_i, α_j)
///
/// where the GB function is:
///
/// f_GB = sqrt(r_ij² + α_i α_j exp(−r_ij²/(4 α_i α_j)))
///
/// Charges are in elementary charge units; positions in metres.
#[derive(Debug, Clone)]
pub struct GeneralizedBorn {
    /// Partial charges q_i (elementary charge units).
    pub charges: Vec<f64>,
    /// Atomic positions (m) stored flat as \[x0,y0,z0, x1,y1,z1, ...\].
    pub positions: Vec<f64>,
    /// Intrinsic (van der Waals) radii ρ_i (m).
    pub intrinsic_radii: Vec<f64>,
    /// Scaling factors S_i for HCT descreening (dimensionless).
    pub scaling_factors: Vec<f64>,
    /// Solvent dielectric constant ε_r.
    pub dielectric: f64,
}
impl GeneralizedBorn {
    /// Construct a new GB model with explicit data.
    ///
    /// All slices must have the same length `n_atoms`.
    /// `positions` has length `3 * n_atoms`.
    pub fn new(
        charges: Vec<f64>,
        positions: Vec<f64>,
        intrinsic_radii: Vec<f64>,
        scaling_factors: Vec<f64>,
        dielectric: f64,
    ) -> Self {
        Self {
            charges,
            positions,
            intrinsic_radii,
            scaling_factors,
            dielectric,
        }
    }
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.charges.len()
    }
    /// Compute the distance between atoms i and j (m).
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        let (ix, iy, iz) = (
            self.positions[3 * i],
            self.positions[3 * i + 1],
            self.positions[3 * i + 2],
        );
        let (jx, jy, jz) = (
            self.positions[3 * j],
            self.positions[3 * j + 1],
            self.positions[3 * j + 2],
        );
        let dx = ix - jx;
        let dy = iy - jy;
        let dz = iz - jz;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Compute the HCT effective Born radius α_i (m) for atom i.
    ///
    /// α_i = 1 / (1/ρ_i − Σ_{j≠i} S_j H(r_ij, ρ_j))
    ///
    /// where H is the HCT pairwise integral approximation.
    pub fn effective_born_radius(&self, i: usize) -> f64 {
        let rho_i = self.intrinsic_radii[i];
        if rho_i <= 0.0 {
            return 1e-10_f64;
        }
        let mut sum = 0.0_f64;
        for j in 0..self.n_atoms() {
            if j == i {
                continue;
            }
            let rij = self.distance(i, j);
            let rho_j = self.intrinsic_radii[j];
            let sj = self.scaling_factors[j];
            let h = if rij > 0.0 {
                (sj * rho_j / rij).tanh() / rij
            } else {
                0.0
            };
            sum += h;
        }
        let inv_alpha = 1.0 / rho_i - sum;
        if inv_alpha <= 0.0 {
            rho_i
        } else {
            1.0 / inv_alpha
        }
    }
    /// GB function f_GB(r, α_i, α_j).
    ///
    /// f = sqrt(r² + α_i α_j exp(−r²/(4 α_i α_j)))
    pub fn gb_function(r: f64, alpha_i: f64, alpha_j: f64) -> f64 {
        let ai_aj = alpha_i * alpha_j;
        if ai_aj <= 0.0 {
            return r;
        }
        let exp_term = (-r * r / (4.0 * ai_aj)).exp();
        (r * r + ai_aj * exp_term).sqrt()
    }
    /// Compute the GB electrostatic solvation energy (kJ/mol).
    ///
    /// ΔG_GB = −(1/2)(1 − 1/ε_r) N_A e² Σ_{i,j} q_i q_j / f_GB(r_ij, α_i, α_j)
    pub fn solvation_energy(&self) -> f64 {
        if self.dielectric <= 0.0 {
            return 0.0;
        }
        let factor = -(1.0 - 1.0 / self.dielectric) * 0.5 * J_TO_KJMOL * ELEM_CHARGE * ELEM_CHARGE
            / (4.0 * PI * EPS0);
        let n = self.n_atoms();
        let mut energy = 0.0_f64;
        let alphas: Vec<f64> = (0..n).map(|i| self.effective_born_radius(i)).collect();
        for i in 0..n {
            for j in 0..n {
                let rij = if i == j { 0.0 } else { self.distance(i, j) };
                let fgb = Self::gb_function(rij, alphas[i], alphas[j]);
                if fgb > 0.0 {
                    energy += self.charges[i] * self.charges[j] / fgb;
                }
            }
        }
        factor * energy
    }
    /// Salt-screening correction to GB energy via a Debye-Hückel term.
    ///
    /// Adds an ionic strength correction:
    /// ΔΔG_salt = −κ (1 − 1/ε_r) N_A e² Σ_{i,j} q_i q_j exp(−κ f_GB) / f_GB / 2
    ///
    /// # Arguments
    /// * `kappa` – inverse Debye length (m⁻¹)
    pub fn salt_correction(&self, kappa: f64) -> f64 {
        if self.dielectric <= 0.0 || kappa <= 0.0 {
            return 0.0;
        }
        let factor = -(1.0 - 1.0 / self.dielectric) * 0.5 * J_TO_KJMOL * ELEM_CHARGE * ELEM_CHARGE
            / (4.0 * PI * EPS0);
        let n = self.n_atoms();
        let mut corr = 0.0_f64;
        let alphas: Vec<f64> = (0..n).map(|i| self.effective_born_radius(i)).collect();
        for i in 0..n {
            for j in 0..n {
                let rij = if i == j { 0.0 } else { self.distance(i, j) };
                let fgb = Self::gb_function(rij, alphas[i], alphas[j]);
                if fgb > 0.0 {
                    let screen = (-kappa * fgb).exp();
                    corr += self.charges[i] * self.charges[j] * (screen - 1.0) / fgb;
                }
            }
        }
        factor * corr
    }
}
/// Solvation free energy decomposed into Born, cavity, dispersion, and
/// repulsion contributions.
///
/// Total: ΔG = ΔG_Born + ΔG_cav + ΔG_disp + ΔG_rep (kJ/mol)
#[derive(Debug, Clone)]
pub struct HydrationFreeEnergy {
    /// Solute charge (elementary charge units).
    pub charge: f64,
    /// Born radius (m).
    pub born_radius: f64,
    /// Solvent relative dielectric ε_r.
    pub dielectric: f64,
    /// Cavity surface tension γ_cav (kJ/mol per Å²).
    pub gamma_cavity: f64,
    /// Solute radius for cavity estimate (Å).
    pub solute_radius: f64,
    /// Dispersion scale ε_disp (kJ/mol).
    pub eps_disp: f64,
    /// Repulsion scale ε_rep (kJ/mol).
    pub eps_rep: f64,
    /// LJ sigma parameter (Å).
    pub sigma_lj: f64,
}
impl HydrationFreeEnergy {
    /// Create a hydration free energy calculator.
    ///
    /// # Arguments
    /// * `charge`       – solute charge (e)
    /// * `born_radius`  – Born radius (m)
    /// * `dielectric`   – solvent ε_r
    /// * `gamma_cavity` – cavity surface tension (kJ/mol Å⁻²)
    /// * `radius`       – solute radius (Å)
    /// * `eps_disp`     – dispersion energy scale (kJ/mol)
    /// * `eps_rep`      – repulsion energy scale (kJ/mol)
    /// * `sigma_lj`     – LJ σ (Å)
    pub fn new(
        charge: f64,
        born_radius: f64,
        dielectric: f64,
        gamma_cavity: f64,
        radius: f64,
        eps_disp: f64,
        eps_rep: f64,
        sigma_lj: f64,
    ) -> Self {
        Self {
            charge,
            born_radius,
            dielectric,
            gamma_cavity,
            solute_radius: radius,
            eps_disp,
            eps_rep,
            sigma_lj,
        }
    }
    /// Born electrostatic contribution ΔG_Born (kJ/mol).
    pub fn born_contribution(&self) -> f64 {
        born_energy(self.charge, self.born_radius, 1.0, self.dielectric)
    }
    /// Cavity formation energy ΔG_cav = γ · 4π r² (kJ/mol).
    pub fn cavity_contribution(&self) -> f64 {
        let sasa = 4.0 * PI * self.solute_radius * self.solute_radius;
        self.gamma_cavity * sasa
    }
    /// Dispersion contribution ΔG_disp (kJ/mol).
    pub fn dispersion_contribution(&self) -> f64 {
        dispersion_energy(self.eps_disp, self.sigma_lj, self.solute_radius)
    }
    /// Repulsion contribution ΔG_rep (kJ/mol).
    pub fn repulsion_contribution(&self) -> f64 {
        repulsion_energy(self.eps_rep, self.sigma_lj, self.solute_radius)
    }
    /// Total solvation free energy ΔG_total = ΔG_Born + ΔG_cav + ΔG_disp + ΔG_rep.
    pub fn total_solvation_energy(&self) -> f64 {
        self.born_contribution()
            + self.cavity_contribution()
            + self.dispersion_contribution()
            + self.repulsion_contribution()
    }
}
/// Hydrophobic solvation model based on solvent-accessible surface area.
///
/// ΔG_np = γ · SASA + ΔG_transfer
///
/// Also provides contact area (overlap) and transfer free energy between
/// two environments.
#[derive(Debug, Clone)]
pub struct HydrophobicEffect {
    /// Solvent-accessible surface area (Å²).
    pub solvent_accessible_surface_area: f64,
    /// Surface tension coefficient γ (kJ/mol Å⁻²).
    pub gamma: f64,
    /// Transfer free energy per unit area ΔG_tr (kJ/mol Å⁻²).
    pub transfer_coeff: f64,
    /// Contact area between two hydrophobic surfaces (Å²).
    pub contact_area: f64,
}
impl HydrophobicEffect {
    /// Create a hydrophobic effect model.
    ///
    /// # Arguments
    /// * `sasa`  – solvent-accessible surface area (Å²)
    /// * `gamma` – surface tension coefficient (kJ/mol Å⁻²)
    pub fn new(sasa: f64, gamma: f64) -> Self {
        Self {
            solvent_accessible_surface_area: sasa,
            gamma,
            transfer_coeff: 0.0,
            contact_area: 0.0,
        }
    }
    /// Set the transfer coefficient and contact area.
    pub fn with_transfer(mut self, tr_coeff: f64, contact_area: f64) -> Self {
        self.transfer_coeff = tr_coeff;
        self.contact_area = contact_area;
        self
    }
    /// Compute the hydrophobic free energy ΔG_np = γ · SASA (kJ/mol).
    pub fn compute_hydrophobic_energy(&self) -> f64 {
        self.gamma * self.solvent_accessible_surface_area
    }
    /// Transfer free energy ΔG_tr = transfer_coeff · contact_area (kJ/mol).
    pub fn transfer_free_energy(&self) -> f64 {
        self.transfer_coeff * self.contact_area
    }
    /// Total non-polar solvation energy including transfer (kJ/mol).
    pub fn total_nonpolar_energy(&self) -> f64 {
        self.compute_hydrophobic_energy() + self.transfer_free_energy()
    }
    /// Scale the SASA by burial fraction (0 = buried, 1 = exposed).
    pub fn buried_energy(&self, fraction_exposed: f64) -> f64 {
        self.gamma * self.solvent_accessible_surface_area * fraction_exposed
    }
    /// Estimate the hydrophobic collapse driving force (kJ/mol).
    ///
    /// Returns the energy reduction when contact area increases from 0 to
    /// `contact_area`.
    pub fn collapse_driving_force(&self) -> f64 {
        -self.transfer_free_energy().abs()
    }
    /// Compute SASA for a sphere given atomic radius and probe radius.
    ///
    /// # Arguments
    /// * `atom_radius`  – atom van der Waals radius (Å)
    /// * `probe_radius` – solvent probe radius (Å)
    pub fn sphere_sasa(atom_radius: f64, probe_radius: f64) -> f64 {
        solvent_accessible_area(atom_radius, probe_radius)
    }
    /// Lee-Richards SASA approximation for a set of spheres (no overlap correction).
    ///
    /// Returns the sum of individual SASAs (upper bound, ignores burial).
    ///
    /// # Arguments
    /// * `radii`        – per-atom radii (Å)
    /// * `probe_radius` – solvent probe radius (Å)
    pub fn sum_sasa(radii: &[f64], probe_radius: f64) -> f64 {
        radii
            .iter()
            .map(|&r| solvent_accessible_area(r, probe_radius))
            .sum()
    }
}
/// Poisson-Boltzmann / Debye-Hückel solvation model.
///
/// Implements:
/// - Ionic strength calculation
/// - Debye-Hückel activity coefficient
/// - Linearised PB Debye-Hückel electrostatic potential
/// - Mean-field electrostatic free energy
#[derive(Debug, Clone)]
pub struct PoissonBoltzmann {
    /// Solvent relative dielectric constant ε_r.
    pub dielectric: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Ion concentrations (mol/m³) for each species.
    pub ion_concentrations: Vec<f64>,
    /// Corresponding ion valences (signed integers as f64).
    pub ion_valences: Vec<f64>,
}
impl PoissonBoltzmann {
    /// Create a new Poisson-Boltzmann model.
    ///
    /// # Arguments
    /// * `eps`  – solvent dielectric
    /// * `t`    – temperature (K)
    pub fn new(eps: f64, t: f64) -> Self {
        Self {
            dielectric: eps,
            temperature: t,
            ion_concentrations: Vec::new(),
            ion_valences: Vec::new(),
        }
    }
    /// Add an ionic species.
    ///
    /// # Arguments
    /// * `concentration` – concentration (mol/m³)
    /// * `valence`       – signed valence (e.g. +1, -2)
    pub fn add_ion(&mut self, concentration: f64, valence: f64) {
        self.ion_concentrations.push(concentration);
        self.ion_valences.push(valence);
    }
    /// Compute the ionic strength I = 0.5 Σ c_i z_i² (mol/m³).
    pub fn ionic_strength(&self) -> f64 {
        0.5 * self
            .ion_concentrations
            .iter()
            .zip(self.ion_valences.iter())
            .map(|(c, z)| c * z * z)
            .sum::<f64>()
    }
    /// Debye screening length λ_D (m).
    pub fn debye_length(&self) -> f64 {
        debye_length(self.dielectric, self.ionic_strength(), self.temperature)
    }
    /// Debye-Hückel activity coefficient ln(γ±) = −A |z+z-| √I.
    ///
    /// Uses the extended Debye-Hückel formula at low ionic strength.
    /// A ≈ 1.172 (mol/m³)^(-1/2) for water at 298 K.
    ///
    /// # Arguments
    /// * `z_plus`  – cation valence
    /// * `z_minus` – anion valence (magnitude)
    pub fn debye_huckel_activity(&self, z_plus: f64, z_minus: f64) -> f64 {
        let i = self.ionic_strength();
        let a = 1.172_f64;
        -a * (z_plus * z_minus).abs() * i.sqrt()
    }
    /// Linearised Debye-Hückel electrostatic potential φ(r) (V).
    ///
    /// φ(r) = (q / 4π ε₀ ε_r r) · exp(−r/λ_D)
    ///
    /// # Arguments
    /// * `charge` – solute charge (C)
    /// * `r`      – distance from solute centre (m)
    pub fn dh_potential(&self, charge: f64, r: f64) -> f64 {
        if r <= 0.0 {
            return f64::INFINITY;
        }
        let ld = self.debye_length();
        let screened = if ld.is_infinite() {
            1.0
        } else {
            (-r / ld).exp()
        };
        charge / (4.0 * PI * EPS0 * self.dielectric * r) * screened
    }
    /// Electrostatic free energy of a point charge q in Debye-Hückel medium (kJ/mol).
    ///
    /// ΔG_el = −q² / (8π ε₀ ε_r a) · κ_D a / (1 + κ_D a) × N_A
    ///
    /// # Arguments
    /// * `q`    – charge (elementary charge units)
    /// * `a`    – ion radius (m)
    pub fn electrostatic_free_energy(&self, q: f64, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let ld = self.debye_length();
        let kappa = if ld.is_infinite() { 0.0 } else { 1.0 / ld };
        let q_c = q * ELEM_CHARGE;
        let prefactor = q_c * q_c / (8.0 * PI * EPS0 * self.dielectric * a);
        let screen_factor = kappa * a / (1.0 + kappa * a);
        -prefactor * screen_factor * J_TO_KJMOL
    }
}
/// Excess chemical potential estimator via Widom test particle insertion.
///
/// This struct wraps the general Widom method and provides analysis of the
/// convergence and uncertainty of the μ_ex estimate.  Each call to
/// [`Self::insert`] registers one random test insertion.
///
/// The estimator uses:
///   μ_ex = −k_B T ln⟨e^{−β ΔU}⟩
///
/// where the average is over all random insertion positions.
#[derive(Debug, Clone)]
pub struct ExcessChemicalPotential {
    /// Temperature (K).
    pub temperature: f64,
    /// Stored Boltzmann weights e^{−β ΔU}.
    pub(super) weights: Vec<f64>,
    /// Stored insertion energies ΔU (kJ/mol).
    pub(super) energies: Vec<f64>,
}
impl ExcessChemicalPotential {
    /// Create a new excess chemical potential sampler.
    ///
    /// # Arguments
    /// * `temperature` – temperature (K)
    pub fn new(temperature: f64) -> Self {
        Self {
            temperature,
            weights: Vec::new(),
            energies: Vec::new(),
        }
    }
    /// kBT in kJ/mol.
    fn kbt(&self) -> f64 {
        8.314_462_618e-3 * self.temperature
    }
    /// Register one test particle insertion with energy ΔU (kJ/mol).
    pub fn insert(&mut self, delta_u: f64) {
        let w = (-delta_u / self.kbt()).exp();
        self.weights.push(w);
        self.energies.push(delta_u);
    }
    /// Estimate μ_ex from all accumulated insertions (kJ/mol).
    ///
    /// Returns `f64::INFINITY` if no insertions are recorded.
    pub fn estimate(&self) -> f64 {
        let n = self.weights.len();
        if n == 0 {
            return f64::INFINITY;
        }
        let mean_w = self.weights.iter().sum::<f64>() / n as f64;
        if mean_w <= 0.0 {
            return f64::INFINITY;
        }
        -self.kbt() * mean_w.ln()
    }
    /// Number of insertions recorded.
    pub fn n_samples(&self) -> usize {
        self.weights.len()
    }
    /// Mean insertion energy ⟨ΔU⟩ (kJ/mol).
    pub fn mean_insertion_energy(&self) -> f64 {
        let n = self.energies.len();
        if n == 0 {
            return 0.0;
        }
        self.energies.iter().sum::<f64>() / n as f64
    }
    /// Statistical uncertainty (standard error) of the mean Boltzmann weight.
    pub fn standard_error(&self) -> f64 {
        let n = self.weights.len();
        if n < 2 {
            return f64::INFINITY;
        }
        let mean = self.weights.iter().sum::<f64>() / n as f64;
        let var = self
            .weights
            .iter()
            .map(|&w| (w - mean).powi(2))
            .sum::<f64>()
            / (n - 1) as f64;
        (var / n as f64).sqrt()
    }
    /// Block-average estimate of μ_ex for convergence analysis.
    ///
    /// Divides data into `n_blocks` blocks, estimates μ_ex from each,
    /// and returns the mean and standard deviation of block estimates.
    pub fn block_average(&self, n_blocks: usize) -> (f64, f64) {
        let n = self.weights.len();
        if n < n_blocks || n_blocks == 0 {
            return (self.estimate(), f64::INFINITY);
        }
        let block_size = n / n_blocks;
        let mus: Vec<f64> = (0..n_blocks)
            .map(|b| {
                let start = b * block_size;
                let end = start + block_size;
                let mean_w = self.weights[start..end].iter().sum::<f64>() / block_size as f64;
                if mean_w <= 0.0 {
                    f64::INFINITY
                } else {
                    -self.kbt() * mean_w.ln()
                }
            })
            .filter(|v| v.is_finite())
            .collect();
        let nb = mus.len() as f64;
        if nb == 0.0 {
            return (f64::INFINITY, f64::INFINITY);
        }
        let mean_mu = mus.iter().sum::<f64>() / nb;
        let std_mu = (mus.iter().map(|&m| (m - mean_mu).powi(2)).sum::<f64>() / nb).sqrt();
        (mean_mu, std_mu)
    }
}
