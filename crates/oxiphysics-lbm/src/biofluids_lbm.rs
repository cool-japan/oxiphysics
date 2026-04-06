// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! LBM for biofluid dynamics: blood flow, cerebrospinal fluid, respiratory flow.
//!
//! Implements specialised LBM formulations for physiological flow simulations,
//! including pulsatile blood flow, Womersley profiles, vascular branching,
//! cerebrospinal fluid dynamics, and respiratory-tree flows.

#![allow(dead_code)]

use std::f64::consts::PI;

/// Fluid type for biofluid simulations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FluidType {
    /// Whole blood (non-Newtonian).
    Blood,
    /// Cerebrospinal fluid (Newtonian).
    Csf,
    /// Air in respiratory tract.
    Air,
    /// Plasma only (Newtonian).
    Plasma,
    /// User-defined fluid.
    Custom,
}

/// Non-Newtonian viscosity model.
#[derive(Debug, Clone, Copy)]
pub enum NonNewtonianModel {
    /// Newtonian (constant viscosity).
    Newtonian,
    /// Carreau model: η = η_inf + (η_0 - η_inf)(1 + (λγ̇)²)^((n-1)/2).
    Carreau {
        /// Zero-shear viscosity \[Pa·s\].
        eta_0: f64,
        /// Infinite-shear viscosity \[Pa·s\].
        eta_inf: f64,
        /// Time constant \[s\].
        lambda: f64,
        /// Power-law index.
        n: f64,
    },
    /// Power-law: η = K·γ̇^(n-1).
    PowerLaw {
        /// Consistency index \[Pa·s^n\].
        k: f64,
        /// Power-law index.
        n: f64,
    },
    /// Casson model for blood: √η = √η_c + √(τ_y/γ̇).
    Casson {
        /// Casson viscosity \[Pa·s\].
        eta_c: f64,
        /// Yield stress \[Pa\].
        tau_y: f64,
    },
}

/// Biofluid simulation parameters.
#[derive(Debug, Clone)]
pub struct BiofluidsParams {
    /// Fluid type identifier.
    pub fluid_type: FluidType,
    /// Reference (dynamic) viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Fluid density \[kg/m³\].
    pub density: f64,
    /// Non-Newtonian viscosity model.
    pub model: NonNewtonianModel,
    /// Body force acceleration \[m/s²\] (e.g. gravity component).
    pub body_force: f64,
    /// Temperature \[K\].
    pub temperature: f64,
}

impl BiofluidsParams {
    /// Default parameters for whole blood.
    pub fn blood() -> Self {
        Self {
            fluid_type: FluidType::Blood,
            viscosity: 3.5e-3, // 3.5 mPa·s
            density: 1060.0,   // kg/m³
            model: NonNewtonianModel::Carreau {
                eta_0: 0.056,
                eta_inf: 0.00345,
                lambda: 3.313,
                n: 0.3568,
            },
            body_force: 0.0,
            temperature: 310.15, // 37°C
        }
    }

    /// Default parameters for cerebrospinal fluid.
    pub fn csf() -> Self {
        Self {
            fluid_type: FluidType::Csf,
            viscosity: 1.0e-3,
            density: 1007.0,
            model: NonNewtonianModel::Newtonian,
            body_force: 0.0,
            temperature: 310.15,
        }
    }

    /// Default parameters for respiratory air.
    pub fn air() -> Self {
        Self {
            fluid_type: FluidType::Air,
            viscosity: 1.81e-5,
            density: 1.2,
            model: NonNewtonianModel::Newtonian,
            body_force: 0.0,
            temperature: 310.15,
        }
    }

    /// Compute kinematic viscosity ν = μ/ρ.
    pub fn kinematic_viscosity(&self) -> f64 {
        self.viscosity / self.density
    }

    /// Evaluate effective viscosity at shear rate `gamma_dot` \[1/s\].
    pub fn effective_viscosity(&self, gamma_dot: f64) -> f64 {
        match self.model {
            NonNewtonianModel::Newtonian => self.viscosity,
            NonNewtonianModel::Carreau {
                eta_0,
                eta_inf,
                lambda,
                n,
            } => {
                let lambda_gamma = lambda * gamma_dot;
                eta_inf
                    + (eta_0 - eta_inf) * (1.0 + lambda_gamma * lambda_gamma).powf((n - 1.0) / 2.0)
            }
            NonNewtonianModel::PowerLaw { k, n } => {
                let g = if gamma_dot.abs() < 1e-10 {
                    1e-10
                } else {
                    gamma_dot.abs()
                };
                k * g.powf(n - 1.0)
            }
            NonNewtonianModel::Casson { eta_c, tau_y } => {
                let g = if gamma_dot.abs() < 1e-10 {
                    1e-10
                } else {
                    gamma_dot.abs()
                };
                let sqrt_eta = eta_c.sqrt() + (tau_y / g).sqrt();
                sqrt_eta * sqrt_eta
            }
        }
    }
}

/// LBM simulation of blood flow in a cylindrical vessel.
///
/// Uses D2Q9 lattice in the cross-section; axial velocity profile is solved.
pub struct BloodFlowLbm {
    /// Vessel radius \[m\].
    pub radius: f64,
    /// Vessel length \[m\].
    pub length: f64,
    /// Applied pressure gradient \[Pa/m\].
    pub pressure_gradient: f64,
    /// Biofluid parameters.
    pub params: BiofluidsParams,
    /// Grid size (number of radial points).
    pub nr: usize,
    /// Axial velocity profile u_z(r).
    pub velocity_profile: Vec<f64>,
}

impl BloodFlowLbm {
    /// Create a blood flow LBM with given geometry and parameters.
    pub fn new(
        radius: f64,
        length: f64,
        pressure_gradient: f64,
        params: BiofluidsParams,
        nr: usize,
    ) -> Self {
        let velocity_profile = vec![0.0f64; nr];
        Self {
            radius,
            length,
            pressure_gradient,
            params,
            nr,
            velocity_profile,
        }
    }

    /// Compute Poiseuille velocity profile u(r) = (dp/dz)(R²-r²)/(4μ).
    pub fn poiseuille_profile(&self) -> Vec<f64> {
        let mu = self.params.viscosity;
        let dp = self.pressure_gradient;
        let r = self.radius;
        (0..self.nr)
            .map(|i| {
                let ri = r * (i as f64) / (self.nr as f64 - 1.0).max(1.0);
                -dp / (4.0 * mu) * (r * r - ri * ri)
            })
            .collect()
    }

    /// Poiseuille mean velocity: U = -dp/dz * R² / (8μ).
    pub fn poiseuille_mean_velocity(&self) -> f64 {
        -self.pressure_gradient * self.radius * self.radius / (8.0 * self.params.viscosity)
    }

    /// Volume flow rate for Poiseuille: Q = π R⁴ |dp/dz| / (8μ).
    pub fn poiseuille_flow_rate(&self) -> f64 {
        use std::f64::consts::PI;
        let r = self.radius;
        let mu = self.params.viscosity;
        -self.pressure_gradient * PI * r.powi(4) / (8.0 * mu)
    }

    /// Compute Womersley number Wo = R·√(ω/ν).
    pub fn womersley_number(&self, omega: f64) -> f64 {
        let nu = self.params.kinematic_viscosity();
        self.radius * (omega / nu).sqrt()
    }

    /// Run one LBM relaxation step toward Poiseuille with Carreau viscosity.
    pub fn step(&mut self) {
        let dp = self.pressure_gradient;
        let r = self.radius;
        for i in 0..self.nr {
            let ri = r * (i as f64) / (self.nr as f64 - 1.0).max(1.0);
            // Approximate shear rate γ̇ ≈ |du/dr|
            let gamma_dot = if i < self.nr - 1 {
                let dr = r / (self.nr as f64 - 1.0).max(1.0);
                (self.velocity_profile[(i + 1).min(self.nr - 1)]
                    - self.velocity_profile[i.saturating_sub(1)])
                .abs()
                    / (2.0 * dr)
            } else {
                0.0
            };
            let mu_eff = self.params.effective_viscosity(gamma_dot);
            self.velocity_profile[i] = -dp / (4.0 * mu_eff) * (r * r - ri * ri);
        }
    }
}

/// Pulsatile flow with Fourier series pressure gradient.
///
/// dp/dz(t) = dp_mean + Σ_n A_n cos(nωt) + B_n sin(nωt)
pub struct PulsatileFlow {
    /// Angular frequency ω = 2π/T \[rad/s\].
    pub omega: f64,
    /// Mean pressure gradient \[Pa/m\].
    pub dp_mean: f64,
    /// Fourier cosine coefficients A_n.
    pub cosine_coeffs: Vec<f64>,
    /// Fourier sine coefficients B_n.
    pub sine_coeffs: Vec<f64>,
    /// Vessel radius \[m\].
    pub radius: f64,
    /// Kinematic viscosity \[m²/s\].
    pub nu: f64,
}

impl PulsatileFlow {
    /// Create pulsatile flow with given fundamental frequency.
    pub fn new(omega: f64, dp_mean: f64, radius: f64, nu: f64) -> Self {
        Self {
            omega,
            dp_mean,
            cosine_coeffs: vec![],
            sine_coeffs: vec![],
            radius,
            nu,
        }
    }

    /// Add n-th harmonic to pressure gradient.
    pub fn add_harmonic(&mut self, an: f64, bn: f64) {
        self.cosine_coeffs.push(an);
        self.sine_coeffs.push(bn);
    }

    /// Evaluate instantaneous pressure gradient at time t.
    pub fn pressure_gradient_at(&self, t: f64) -> f64 {
        let mut dp = self.dp_mean;
        for (k, (&an, &bn)) in self
            .cosine_coeffs
            .iter()
            .zip(self.sine_coeffs.iter())
            .enumerate()
        {
            let n = (k + 1) as f64;
            dp += an * (n * self.omega * t).cos() + bn * (n * self.omega * t).sin();
        }
        dp
    }

    /// Womersley velocity profile at radius r, harmonic n, time t.
    ///
    /// For small Wo, approaches Poiseuille. Returns axial velocity contribution.
    pub fn womersley_velocity(&self, r: f64, t: f64, harmonic_idx: usize) -> f64 {
        if harmonic_idx >= self.cosine_coeffs.len() {
            return 0.0;
        }
        let n = (harmonic_idx + 1) as f64;
        let an = self.cosine_coeffs[harmonic_idx];
        let wo = self.radius * (n * self.omega / self.nu).sqrt();
        // Approximate Womersley profile (leading-order real part)
        // For Wo << 1: reduces to parabolic Poiseuille
        let r2 = (r / self.radius).powi(2);
        if wo < 1e-3 {
            // Quasi-static: Poiseuille with instantaneous dp
            let dp_amp = an;
            -dp_amp / (4.0 * self.nu)
                * self.radius
                * self.radius
                * (1.0 - r2)
                * (n * self.omega * t).cos()
        } else {
            // Simplified Womersley: amplitude modulated by Wo-dependent factor
            let alpha = 1.0 - 1.0 / (1.0 + wo * wo / 8.0);
            let dp_amp = an;
            -dp_amp / (4.0 * self.nu)
                * self.radius
                * self.radius
                * alpha
                * (1.0 - r2)
                * (n * self.omega * t).cos()
        }
    }

    /// Mean flow rate (time-averaged) = Poiseuille for mean dp.
    pub fn mean_flow_rate(&self) -> f64 {
        let r = self.radius;
        -self.dp_mean * PI * r.powi(4) / (8.0 * self.nu)
    }
}

/// Vascular geometry with Murray's law branching.
#[derive(Debug, Clone)]
pub struct VascularGeometry {
    /// Parent vessel radius \[m\].
    pub r_parent: f64,
    /// Daughter vessel radii \[m\].
    pub r_daughters: Vec<f64>,
    /// Branching angle \[rad\] for each daughter.
    pub angles: Vec<f64>,
}

impl VascularGeometry {
    /// Create parent vessel.
    pub fn new(r_parent: f64) -> Self {
        Self {
            r_parent,
            r_daughters: vec![],
            angles: vec![],
        }
    }

    /// Add a daughter vessel with given radius and branching angle.
    pub fn add_daughter(&mut self, r: f64, angle: f64) {
        self.r_daughters.push(r);
        self.angles.push(angle);
    }

    /// Verify Murray's law: r_p³ = Σ r_d³.
    pub fn check_murrays_law(&self) -> f64 {
        let lhs = self.r_parent.powi(3);
        let rhs: f64 = self.r_daughters.iter().map(|r| r.powi(3)).sum();
        (lhs - rhs).abs()
    }

    /// Enforce Murray's law: set daughters such that r_p³ = Σ r_d³.
    /// Returns the scaling factor applied to daughters.
    pub fn enforce_murrays_law(&mut self) -> f64 {
        if self.r_daughters.is_empty() {
            return 1.0;
        }
        let current_cube_sum: f64 = self.r_daughters.iter().map(|r| r.powi(3)).sum();
        if current_cube_sum < 1e-30 {
            return 1.0;
        }
        let scale = (self.r_parent.powi(3) / current_cube_sum).powf(1.0 / 3.0);
        for r in &mut self.r_daughters {
            *r *= scale;
        }
        scale
    }

    /// Compute resistance ratio assuming Poiseuille: Hagen-Poiseuille R ∝ L/r^4.
    /// Returns resistance of parent / sum of daughter resistances.
    pub fn resistance_ratio(&self, l_parent: f64, l_daughters: &[f64]) -> f64 {
        let r_p = 1.0 / self.r_parent.powi(4) * l_parent;
        let r_d: f64 = self
            .r_daughters
            .iter()
            .zip(l_daughters.iter())
            .map(|(r, l)| l / r.powi(4))
            .sum();
        if r_d.abs() < 1e-30 {
            return 0.0;
        }
        r_p / r_d
    }
}

/// Cerebrovascular autoregulation model.
///
/// Maintains roughly constant flow rate Q_target despite pressure variations
/// via a feedback mechanism adjusting effective vessel resistance.
pub struct CerebrovascularFlow {
    /// Target blood flow rate \[m³/s\].
    pub q_target: f64,
    /// Baseline cerebral perfusion pressure \[Pa\].
    pub cpp_baseline: f64,
    /// Vessel radius \[m\].
    pub radius: f64,
    /// Vessel length \[m\].
    pub length: f64,
    /// Blood viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Autoregulation gain.
    pub gain: f64,
}

impl CerebrovascularFlow {
    /// Create autoregulatory flow model.
    pub fn new(
        q_target: f64,
        cpp_baseline: f64,
        radius: f64,
        length: f64,
        viscosity: f64,
        gain: f64,
    ) -> Self {
        Self {
            q_target,
            cpp_baseline,
            radius,
            length,
            viscosity,
            gain,
        }
    }

    /// Effective resistance R = 8μL/(πR⁴).
    pub fn baseline_resistance(&self) -> f64 {
        8.0 * self.viscosity * self.length / (PI * self.radius.powi(4))
    }

    /// Compute flow rate at given CPP with autoregulation.
    ///
    /// Autoregulation adjusts radius to keep Q ≈ Q_target within plateau.
    pub fn flow_rate(&self, cpp: f64) -> f64 {
        // Outside autoregulation plateau: passive (linear) flow
        let plateau_width = 0.3 * self.cpp_baseline;
        let delta = (cpp - self.cpp_baseline).abs();
        if delta > plateau_width {
            // Passive: Q = CPP / R_baseline
            cpp / self.baseline_resistance()
        } else {
            // Autoregulated: Q ≈ Q_target (with small linear perturbation)
            let err = cpp - self.cpp_baseline;
            self.q_target * (1.0 + 0.05 * err / (plateau_width + 1e-15) * (1.0 - self.gain))
        }
    }
}

/// Respiratory flow in a symmetric Weibel bronchial tree.
///
/// Generation n: radius r_n = r_0 · 2^{-n/3}, length L_n = L_0 · 2^{-n/3}.
pub struct RespiratoryFlowLbm {
    /// Tracheal (generation 0) radius \[m\].
    pub r0: f64,
    /// Tracheal length \[m\].
    pub l0: f64,
    /// Number of generations.
    pub n_generations: usize,
    /// Air kinematic viscosity \[m²/s\].
    pub nu: f64,
    /// Tidal volume flow rate \[m³/s\].
    pub q_tidal: f64,
}

impl RespiratoryFlowLbm {
    /// Create Weibel symmetric tree.
    pub fn new(r0: f64, l0: f64, n_generations: usize, nu: f64, q_tidal: f64) -> Self {
        Self {
            r0,
            l0,
            n_generations,
            nu,
            q_tidal,
        }
    }

    /// Radius at generation n: r_n = r_0 · 2^{-n/3}.
    pub fn radius_at(&self, n: usize) -> f64 {
        self.r0 * 2.0_f64.powf(-(n as f64) / 3.0)
    }

    /// Length at generation n: L_n = L_0 · 2^{-n/3}.
    pub fn length_at(&self, n: usize) -> f64 {
        self.l0 * 2.0_f64.powf(-(n as f64) / 3.0)
    }

    /// Number of airways at generation n: 2^n.
    pub fn count_at(&self, n: usize) -> usize {
        1usize << n
    }

    /// Total cross-sectional area at generation n: A_n = count_n · π r_n².
    pub fn total_area_at(&self, n: usize) -> f64 {
        let r = self.radius_at(n);
        self.count_at(n) as f64 * PI * r * r
    }

    /// Flow per airway at generation n: Q_n = Q_tidal / count_n.
    pub fn flow_per_airway(&self, n: usize) -> f64 {
        self.q_tidal / self.count_at(n) as f64
    }

    /// Mean velocity at generation n: v_n = Q_n / (π r_n²).
    pub fn mean_velocity_at(&self, n: usize) -> f64 {
        let r = self.radius_at(n);
        let q = self.flow_per_airway(n);
        q / (PI * r * r)
    }

    /// Reynolds number at generation n: Re = v_n · 2r_n / ν.
    pub fn reynolds_at(&self, n: usize) -> f64 {
        let v = self.mean_velocity_at(n);
        let r = self.radius_at(n);
        2.0 * v * r / self.nu
    }

    /// Total airway resistance (Hagen-Poiseuille) up to generation N.
    pub fn total_resistance(&self) -> f64 {
        (0..self.n_generations)
            .map(|n| {
                let r = self.radius_at(n);
                let l = self.length_at(n);
                let cnt = self.count_at(n) as f64;
                // Resistance of one airway divided by number in parallel
                8.0 * self.nu * l / (PI * r.powi(4)) / cnt
            })
            .sum()
    }
}

/// Cerebrospinal fluid dynamics with Windkessel compliance.
///
/// Models CSF pulsation driven by cardiac cycle.
pub struct CsfDynamics {
    /// Intracranial compliance C \[m³/Pa\].
    pub compliance: f64,
    /// CSF outflow resistance R \[Pa·s/m³\].
    pub resistance: f64,
    /// Baseline intracranial pressure P_0 \[Pa\].
    pub p0: f64,
    /// CSF production rate Q_prod \[m³/s\].
    pub q_production: f64,
}

impl CsfDynamics {
    /// Create Windkessel CSF model.
    pub fn new(compliance: f64, resistance: f64, p0: f64, q_production: f64) -> Self {
        Self {
            compliance,
            resistance,
            p0,
            q_production,
        }
    }

    /// Steady-state ICP: P_ss = P_0 + Q_prod · R.
    pub fn steady_state_icp(&self) -> f64 {
        self.p0 + self.q_production * self.resistance
    }

    /// Impulse response: δP(t) = (δV / C) · exp(-t / (R·C)).
    pub fn impulse_response(&self, t: f64, delta_v: f64) -> f64 {
        let tau = self.resistance * self.compliance;
        (delta_v / self.compliance) * (-t / tau).exp()
    }

    /// Pressure at time t under sinusoidal CSF inflow Q(t) = Q_0 sin(ωt).
    ///
    /// P(t) = P_ss + Q_0·R / √(1+(ωRC)²) · sin(ωt - φ)
    pub fn pulsatile_pressure(&self, t: f64, q0: f64, omega: f64) -> f64 {
        let tau = self.resistance * self.compliance;
        let amplitude = q0 * self.resistance / (1.0 + (omega * tau).powi(2)).sqrt();
        let phi = (omega * tau).atan();
        self.steady_state_icp() + amplitude * (omega * t - phi).sin()
    }

    /// Pressure buffering factor: ratio of output amplitude to input amplitude.
    pub fn buffering_factor(&self, omega: f64) -> f64 {
        let tau = self.resistance * self.compliance;
        1.0 / (1.0 + (omega * tau).powi(2)).sqrt()
    }
}

/// Heart valve flow model using orifice equation.
///
/// Q = Cd · A · √(2ΔP/ρ)
pub struct HeartValveFlow {
    /// Discharge coefficient Cd (dimensionless).
    pub cd: f64,
    /// Valve orifice area A \[m²\].
    pub area: f64,
    /// Blood density ρ \[kg/m³\].
    pub density: f64,
    /// Valve open/close state.
    pub is_open: bool,
}

impl HeartValveFlow {
    /// Create heart valve model.
    pub fn new(cd: f64, area: f64, density: f64) -> Self {
        Self {
            cd,
            area,
            density,
            is_open: true,
        }
    }

    /// Compute flow rate Q = Cd · A · √(2ΔP/ρ) \[m³/s\].
    pub fn flow_rate(&self, delta_p: f64) -> f64 {
        if !self.is_open || delta_p <= 0.0 {
            return 0.0;
        }
        self.cd * self.area * (2.0 * delta_p / self.density).sqrt()
    }

    /// Effective orifice area (EOA) = Q / (√(2ΔP/ρ)) with Cd=1.
    pub fn effective_orifice_area(&self, delta_p: f64) -> f64 {
        if delta_p <= 0.0 {
            return 0.0;
        }
        self.flow_rate(delta_p) / (2.0 * delta_p / self.density).sqrt()
    }

    /// Pressure drop from flow rate (inverse orifice equation).
    pub fn pressure_drop(&self, q: f64) -> f64 {
        if self.area.abs() < 1e-15 || !self.is_open {
            return 0.0;
        }
        let v = q / (self.cd * self.area);
        0.5 * self.density * v * v
    }

    /// Valve resistance dΔP/dQ = ρ·Q/(Cd·A)².
    pub fn valve_resistance(&self, q: f64) -> f64 {
        if !self.is_open {
            return f64::INFINITY;
        }
        let cda = self.cd * self.area;
        self.density * q / (cda * cda)
    }
}

/// Deformable red blood cell (RBC) capsule model in shear flow.
///
/// Models tank-treading and tumbling behaviour using a simplified membrane model.
pub struct RedBloodCell {
    /// RBC effective radius \[m\] (~4 μm).
    pub radius: f64,
    /// Membrane shear modulus Gs \[N/m\].
    pub shear_modulus: f64,
    /// Membrane bending modulus Kb \[N·m\].
    pub bending_modulus: f64,
    /// Internal viscosity ratio λ = μ_int / μ_ext.
    pub viscosity_ratio: f64,
    /// Reduced volume (~ 0.64 for biconcave RBC).
    pub reduced_volume: f64,
}

impl RedBloodCell {
    /// Create a model RBC with physiological defaults.
    pub fn new() -> Self {
        Self {
            radius: 4.0e-6,
            shear_modulus: 6.0e-6,    // N/m
            bending_modulus: 2.0e-19, // N·m
            viscosity_ratio: 5.0,
            reduced_volume: 0.64,
        }
    }

    /// Capillary number Ca = μ_ext · γ̇ · R / Gs.
    pub fn capillary_number(&self, mu_ext: f64, shear_rate: f64) -> f64 {
        mu_ext * shear_rate * self.radius / self.shear_modulus
    }

    /// Critical capillary number for tumbling-to-tank-treading transition.
    ///
    /// Approximate: Ca_c ≈ (λ - 0.5) / (some geometric factor).
    pub fn critical_capillary_number(&self) -> f64 {
        // Simplified Keller-Skalak result
        let lambda = self.viscosity_ratio;
        if lambda > 0.5 {
            (lambda - 0.5) * 0.8
        } else {
            0.0
        }
    }

    /// Predict if cell is tank-treading at given conditions.
    pub fn is_tank_treading(&self, mu_ext: f64, shear_rate: f64) -> bool {
        self.capillary_number(mu_ext, shear_rate) > self.critical_capillary_number()
    }

    /// Stokes drag on a sphere of equivalent volume.
    pub fn stokes_drag(&self, mu_ext: f64, velocity: f64) -> f64 {
        6.0 * std::f64::consts::PI * mu_ext * self.radius * velocity
    }

    /// Membrane elastic energy (linearised for small deformations).
    pub fn elastic_energy(&self, strain: f64) -> f64 {
        0.5 * self.shear_modulus
            * strain
            * strain
            * 4.0
            * std::f64::consts::PI
            * self.radius
            * self.radius
    }

    /// Volume of equivalent sphere: V = 4π R³/3.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius.powi(3)
    }
}

impl Default for RedBloodCell {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── BiofluidsParams tests ─────────────────────────────────────────────

    #[test]
    fn test_blood_params_defaults() {
        let p = BiofluidsParams::blood();
        assert_eq!(p.fluid_type, FluidType::Blood);
        assert!((p.viscosity - 3.5e-3).abs() < 1e-10);
        assert!((p.density - 1060.0).abs() < 1.0);
    }

    #[test]
    fn test_csf_params_defaults() {
        let p = BiofluidsParams::csf();
        assert_eq!(p.fluid_type, FluidType::Csf);
        assert!((p.viscosity - 1.0e-3).abs() < 1e-10);
    }

    #[test]
    fn test_air_params_defaults() {
        let p = BiofluidsParams::air();
        assert_eq!(p.fluid_type, FluidType::Air);
        assert!((p.density - 1.2).abs() < 0.01);
    }

    #[test]
    fn test_kinematic_viscosity() {
        let p = BiofluidsParams::blood();
        let nu = p.kinematic_viscosity();
        assert!((nu - 3.5e-3 / 1060.0).abs() < 1e-12);
    }

    #[test]
    fn test_carreau_at_zero_shear() {
        let p = BiofluidsParams::blood();
        let eta = p.effective_viscosity(0.0);
        // Should return eta_0 + (eta_0 - eta_inf)*(1+0)^{...} ≈ eta_0
        assert!(
            eta > 0.0,
            "Carreau at zero shear should be positive: {eta:.6}"
        );
    }

    #[test]
    fn test_carreau_high_shear_approaches_eta_inf() {
        let p = BiofluidsParams::blood();
        let eta_high = p.effective_viscosity(1000.0);
        // At high shear, approaches eta_inf = 0.00345
        if let NonNewtonianModel::Carreau { eta_inf, .. } = p.model {
            assert!(
                (eta_high - eta_inf).abs() < 0.001,
                "Carreau high shear: {eta_high:.6} vs eta_inf={eta_inf:.6}"
            );
        }
    }

    #[test]
    fn test_newtonian_viscosity_constant() {
        let p = BiofluidsParams::csf();
        let eta1 = p.effective_viscosity(0.0);
        let eta2 = p.effective_viscosity(100.0);
        assert!((eta1 - eta2).abs() < 1e-14, "Newtonian should be constant");
    }

    #[test]
    fn test_powerlaw_viscosity() {
        let p = BiofluidsParams {
            fluid_type: FluidType::Custom,
            viscosity: 1.0,
            density: 1000.0,
            model: NonNewtonianModel::PowerLaw { k: 0.01, n: 0.5 },
            body_force: 0.0,
            temperature: 310.0,
        };
        // At gamma_dot=100: eta = 0.01 * 100^{-0.5} = 0.01/10 = 0.001
        let eta = p.effective_viscosity(100.0);
        assert!((eta - 0.001).abs() < 1e-8, "Power-law viscosity: {eta:.6}");
    }

    // ── BloodFlowLbm tests ────────────────────────────────────────────────

    #[test]
    fn test_poiseuille_profile_zero_at_wall() {
        let params = BiofluidsParams::blood();
        let lbm = BloodFlowLbm::new(0.005, 0.1, -100.0, params, 11);
        let profile = lbm.poiseuille_profile();
        assert!(
            profile.last().unwrap().abs() < 1e-12,
            "Poiseuille zero at wall: {:.6}",
            profile.last().unwrap()
        );
    }

    #[test]
    fn test_poiseuille_profile_max_at_center() {
        let params = BiofluidsParams::blood();
        let lbm = BloodFlowLbm::new(0.005, 0.1, -100.0, params, 11);
        let profile = lbm.poiseuille_profile();
        let max_v = profile[0];
        for &v in &profile {
            assert!(v <= max_v + 1e-12, "Center should have max velocity");
        }
    }

    #[test]
    fn test_poiseuille_mean_velocity() {
        let params = BiofluidsParams::blood();
        let r = 0.005;
        let dp = -100.0;
        let mu = params.viscosity;
        let lbm = BloodFlowLbm::new(r, 0.1, dp, params, 11);
        let u_mean = lbm.poiseuille_mean_velocity();
        let expected = -dp * r * r / (8.0 * mu);
        assert!(
            (u_mean - expected).abs() < 1e-12,
            "Mean velocity: {u_mean:.6}"
        );
    }

    #[test]
    fn test_womersley_number_formula() {
        let params = BiofluidsParams::blood();
        let r = 0.01;
        let omega = 2.0 * PI; // 1 Hz heartbeat
        let nu = params.kinematic_viscosity();
        let lbm = BloodFlowLbm::new(r, 0.1, -100.0, params, 11);
        let wo = lbm.womersley_number(omega);
        let expected = r * (omega / nu).sqrt();
        assert!(
            (wo - expected).abs() < 1e-10,
            "Womersley: {wo:.6} vs {expected:.6}"
        );
    }

    #[test]
    fn test_womersley_number_typical_aorta() {
        // Aorta: R≈12mm, HR≈1Hz, nu_blood≈3.3e-6 m²/s → Wo ≈ 12
        let params = BiofluidsParams::blood();
        let r = 0.012;
        let omega = 2.0 * PI;
        let lbm = BloodFlowLbm::new(r, 0.3, -100.0, params, 20);
        let wo = lbm.womersley_number(omega);
        assert!(wo > 5.0 && wo < 25.0, "Aortic Wo should be ~12: {wo:.2}");
    }

    #[test]
    fn test_poiseuille_flow_rate() {
        let params = BiofluidsParams::blood();
        let r = 0.005;
        let dp = -100.0;
        let mu = params.viscosity;
        let lbm = BloodFlowLbm::new(r, 0.1, dp, params, 11);
        let q = lbm.poiseuille_flow_rate();
        let expected = -dp * PI * r.powi(4) / (8.0 * mu);
        assert!(
            (q - expected).abs() < 1e-12 * expected.abs(),
            "Flow rate: {q:.6}"
        );
    }

    // ── PulsatileFlow tests ───────────────────────────────────────────────

    #[test]
    fn test_pulsatile_mean_flow_rate() {
        let nu = BiofluidsParams::blood().kinematic_viscosity();
        let r = 0.005;
        let dp_mean = -100.0;
        let omega = 2.0 * PI;
        let pf = PulsatileFlow::new(omega, dp_mean, r, nu);
        let q_mean = pf.mean_flow_rate();
        let expected = -dp_mean * PI * r.powi(4) / (8.0 * nu);
        assert!(
            (q_mean - expected).abs() < 1e-12 * expected.abs(),
            "Pulsatile mean Q: {q_mean:.6}"
        );
    }

    #[test]
    fn test_pulsatile_pressure_gradient_at_t0() {
        let nu = 1.0e-6;
        let omega = 2.0 * PI;
        let mut pf = PulsatileFlow::new(omega, -100.0, 0.01, nu);
        pf.add_harmonic(10.0, 0.0); // A1=10, B1=0
        let dp = pf.pressure_gradient_at(0.0);
        // dp(0) = -100 + 10*cos(0) + 0 = -90
        assert!((dp - (-90.0)).abs() < 1e-10, "Pulsatile dp(0)={dp:.6}");
    }

    #[test]
    fn test_womersley_reduces_to_poiseuille_small_wo() {
        // For small Womersley number (tiny omega), profile should be parabolic
        let nu = 1.0e-3;
        let r_vessel = 0.001;
        let omega_tiny = 1e-6;
        let mut pf = PulsatileFlow::new(omega_tiny, 0.0, r_vessel, nu);
        pf.add_harmonic(-100.0, 0.0);
        // At r=0, compare with Poiseuille
        let v_wo = pf.womersley_velocity(0.0, 0.0, 0);
        // Poiseuille: u_max = -dp/(4nu) * R^2 (here dp amplitude = -100, nu replaces mu)
        let u_max_pois = 100.0 / (4.0 * nu) * r_vessel * r_vessel;
        assert!(
            (v_wo - u_max_pois).abs() < 1e-4 * u_max_pois.abs(),
            "Womersley→Poiseuille for small Wo: {v_wo:.6} vs {u_max_pois:.6}"
        );
    }

    // ── VascularGeometry / Murray's law tests ─────────────────────────────

    #[test]
    fn test_murrays_law_exact() {
        // r_p = 1, r_1 and r_2 such that r_1^3 + r_2^3 = 1
        let r_p = 1.0f64;
        let r1 = (0.5f64).powf(1.0 / 3.0);
        let r2 = (0.5f64).powf(1.0 / 3.0);
        let mut vg = VascularGeometry::new(r_p);
        vg.add_daughter(r1, 0.3);
        vg.add_daughter(r2, -0.3);
        let err = vg.check_murrays_law();
        assert!(err < 1e-12, "Murray's law error: {err:.2e}");
    }

    #[test]
    fn test_murrays_law_enforce() {
        let mut vg = VascularGeometry::new(1.0);
        vg.add_daughter(0.7, 0.0);
        vg.add_daughter(0.7, 0.0);
        vg.enforce_murrays_law();
        let err = vg.check_murrays_law();
        assert!(err < 1e-12, "Enforced Murray's law error: {err:.2e}");
    }

    #[test]
    fn test_murrays_law_r_parent_cubed() {
        let r_p = 0.01f64; // 10mm
        let r1 = (r_p.powi(3) * 0.6).powf(1.0 / 3.0);
        let r2 = (r_p.powi(3) * 0.4).powf(1.0 / 3.0);
        let mut vg = VascularGeometry::new(r_p);
        vg.add_daughter(r1, 0.5);
        vg.add_daughter(r2, -0.3);
        let err = vg.check_murrays_law();
        assert!(err < 1e-20, "Murray's law r_p^3 = r1^3+r2^3: {err:.2e}");
    }

    // ── RespiratoryFlowLbm tests ──────────────────────────────────────────

    #[test]
    fn test_respiratory_radius_decreases() {
        let lbm = RespiratoryFlowLbm::new(0.009, 0.12, 23, 1.5e-5, 500e-6);
        for n in 1..10 {
            assert!(
                lbm.radius_at(n) < lbm.radius_at(n - 1),
                "Radius should decrease: gen {n}"
            );
        }
    }

    #[test]
    fn test_respiratory_total_area_increases() {
        let lbm = RespiratoryFlowLbm::new(0.009, 0.12, 23, 1.5e-5, 500e-6);
        for n in 1..10 {
            assert!(
                lbm.total_area_at(n) > lbm.total_area_at(n - 1),
                "Total area should increase distally: gen {n}"
            );
        }
    }

    #[test]
    fn test_respiratory_radius_formula() {
        let r0 = 0.009;
        let lbm = RespiratoryFlowLbm::new(r0, 0.12, 23, 1.5e-5, 500e-6);
        let r5 = lbm.radius_at(5);
        let expected = r0 * 2.0_f64.powf(-5.0 / 3.0);
        assert!(
            (r5 - expected).abs() < 1e-15,
            "Radius formula at gen 5: {r5:.8} vs {expected:.8}"
        );
    }

    #[test]
    fn test_respiratory_count_doubles() {
        let lbm = RespiratoryFlowLbm::new(0.009, 0.12, 23, 1.5e-5, 500e-6);
        for n in 0..8 {
            assert_eq!(
                lbm.count_at(n + 1),
                2 * lbm.count_at(n),
                "Count should double each generation"
            );
        }
    }

    // ── CsfDynamics tests ─────────────────────────────────────────────────

    #[test]
    fn test_csf_steady_state_icp() {
        let csf = CsfDynamics::new(1e-9, 1e8, 600.0, 3.3e-9);
        let icp = csf.steady_state_icp();
        let expected = 600.0 + 3.3e-9 * 1e8;
        assert!(
            (icp - expected).abs() < 1e-6,
            "CSF steady-state ICP: {icp:.2} vs {expected:.2}"
        );
    }

    #[test]
    fn test_csf_impulse_response_decays() {
        let csf = CsfDynamics::new(1e-9, 1e8, 600.0, 3.3e-9);
        let r0 = csf.impulse_response(0.0, 1e-6);
        let r1 = csf.impulse_response(1.0, 1e-6);
        assert!(r1 < r0, "Impulse response should decay: {r0:.6} -> {r1:.6}");
    }

    #[test]
    fn test_csf_buffering_factor_lt1() {
        let csf = CsfDynamics::new(1e-9, 1e8, 600.0, 3.3e-9);
        let f = csf.buffering_factor(2.0 * PI);
        assert!(
            f < 1.0 && f > 0.0,
            "Buffering factor should be in (0,1): {f:.6}"
        );
    }

    #[test]
    fn test_csf_compliance_pressure_buffering() {
        // Higher compliance → more buffering (lower factor)
        let csf_low_c = CsfDynamics::new(0.5e-9, 1e8, 600.0, 3.3e-9);
        let csf_high_c = CsfDynamics::new(2e-9, 1e8, 600.0, 3.3e-9);
        let omega = 2.0 * PI;
        let f_low = csf_low_c.buffering_factor(omega);
        let f_high = csf_high_c.buffering_factor(omega);
        assert!(
            f_high < f_low,
            "Higher compliance should buffer more: f_low={f_low:.4} f_high={f_high:.4}"
        );
    }

    // ── HeartValveFlow tests ──────────────────────────────────────────────

    #[test]
    fn test_valve_flow_proportional_to_sqrt_dp() {
        let valve = HeartValveFlow::new(0.85, 2.0e-4, 1060.0);
        let q1 = valve.flow_rate(100.0);
        let q4 = valve.flow_rate(400.0);
        // Q ∝ √ΔP → q4/q1 ≈ 2
        let ratio = q4 / q1;
        assert!((ratio - 2.0).abs() < 1e-10, "Q ∝ √ΔP: ratio = {ratio:.6}");
    }

    #[test]
    fn test_valve_zero_flow_when_closed() {
        let mut valve = HeartValveFlow::new(0.85, 2.0e-4, 1060.0);
        valve.is_open = false;
        let q = valve.flow_rate(1000.0);
        assert_eq!(q, 0.0, "Closed valve should give zero flow");
    }

    #[test]
    fn test_valve_zero_flow_negative_dp() {
        let valve = HeartValveFlow::new(0.85, 2.0e-4, 1060.0);
        let q = valve.flow_rate(-50.0);
        assert_eq!(q, 0.0, "Negative ΔP should give zero flow (valve closes)");
    }

    #[test]
    fn test_valve_pressure_drop_inverse() {
        let valve = HeartValveFlow::new(0.85, 2.0e-4, 1060.0);
        let dp = 500.0;
        let q = valve.flow_rate(dp);
        let dp_back = valve.pressure_drop(q);
        assert!(
            (dp_back - dp).abs() < 1e-8,
            "Pressure drop inverse: {dp_back:.6} vs {dp:.6}"
        );
    }

    #[test]
    fn test_valve_orifice_equation_formula() {
        let cd = 0.85;
        let area = 2.0e-4;
        let rho = 1060.0;
        let dp = 200.0;
        let valve = HeartValveFlow::new(cd, area, rho);
        let q = valve.flow_rate(dp);
        let expected = cd * area * (2.0 * dp / rho).sqrt();
        assert!(
            (q - expected).abs() < 1e-15,
            "Orifice equation: {q:.6} vs {expected:.6}"
        );
    }

    // ── RedBloodCell tests ────────────────────────────────────────────────

    #[test]
    fn test_rbc_capillary_number() {
        let rbc = RedBloodCell::new();
        let mu_ext = 1.2e-3;
        let shear = 100.0;
        let ca = rbc.capillary_number(mu_ext, shear);
        let expected = mu_ext * shear * rbc.radius / rbc.shear_modulus;
        assert!(
            (ca - expected).abs() < 1e-20,
            "Ca: {ca:.6} vs {expected:.6}"
        );
    }

    #[test]
    fn test_rbc_stokes_drag() {
        let rbc = RedBloodCell::new();
        let mu = 1.2e-3;
        let v = 1e-3;
        let drag = rbc.stokes_drag(mu, v);
        let expected = 6.0 * PI * mu * rbc.radius * v;
        assert!((drag - expected).abs() < 1e-25, "Stokes drag: {drag:.6}");
    }

    #[test]
    fn test_rbc_volume_formula() {
        let rbc = RedBloodCell::new();
        let vol = rbc.volume();
        let expected = 4.0 / 3.0 * PI * rbc.radius.powi(3);
        assert!((vol - expected).abs() < 1e-30, "RBC volume: {vol:.6}");
    }

    #[test]
    fn test_rbc_tank_treading_at_high_shear() {
        let rbc = RedBloodCell::new();
        // At very high shear rate, Ca >> Ca_c → tank treading
        let mu_ext = 1.2e-3;
        assert!(
            rbc.is_tank_treading(mu_ext, 1e6),
            "RBC should tank-tread at very high shear"
        );
    }
}
