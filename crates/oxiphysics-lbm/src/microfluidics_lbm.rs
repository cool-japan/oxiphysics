#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Microfluidic Lattice Boltzmann Method module.
//!
//! Implements LBM-based simulation of microfluidic phenomena including:
//!
//! - Knudsen number effects and Maxwell slip boundary conditions
//! - Electroosmotic flow (Debye layer, zeta potential)
//! - Electrophoresis of charged particles
//! - Droplet microfluidics: T-junction and flow-focusing geometries
//! - Digital microfluidics via electrowetting-on-dielectric (EWOD)
//! - Capillary electrophoresis separation
//! - Diffusiophoresis driven by concentration gradients
//! - Inertial microfluidics and Dean flow in curved channels
//! - PDMS channel wall deformation coupling
//!
//! # References
//! - Karniadakis, G., Beskok, A., Aluru, N. (2005). *Microflows and Nanoflows*.
//! - Li, D. (2004). *Electrokinetics in Microfluidics*.
//! - Stone, H. A., Stroock, A. D., Ajdari, A. (2004). Annu. Rev. Fluid Mech.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Elementary charge *e* (C).
pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;

/// Vacuum permittivity ε₀ (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// Boltzmann constant k_B (J K⁻¹).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;

/// Avogadro constant N_A (mol⁻¹).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Universal gas constant R = k_B · N_A (J mol⁻¹ K⁻¹).
pub const GAS_CONSTANT: f64 = 8.314_462_618;

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 discrete velocity vectors \[cx, cy\].
const D2Q9_C: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];

/// D2Q9 equilibrium weights.
const D2Q9_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Speed of sound squared for D2Q9 (lattice units): cs² = 1/3.
pub const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// Free functions — dimensionless numbers and physical relationships
// ---------------------------------------------------------------------------

/// Knudsen number Kn = λ / L.
///
/// `lambda` — mean free path (m), `length` — characteristic length (m).
pub fn knudsen_number(lambda: f64, length: f64) -> f64 {
    lambda / length
}

/// Debye screening length λ_D (m).
///
/// λ_D = √(ε_r ε₀ k_B T / (2 n₀ z² e²))
///
/// `eps_r` — relative permittivity, `temp_k` — temperature (K),
/// `n0` — bulk ion number density (m⁻³), `z` — ion valence.
pub fn debye_length(eps_r: f64, temp_k: f64, n0: f64, z: f64) -> f64 {
    let num = eps_r * EPSILON_0 * K_BOLTZMANN * temp_k;
    let den = 2.0 * n0 * z * z * ELEM_CHARGE * ELEM_CHARGE;
    (num / den).sqrt()
}

/// Electro-osmotic velocity (Helmholtz–Smoluchowski formula).
///
/// u_eo = -ε_r ε₀ ζ E_x / μ
///
/// `zeta` — zeta potential (V), `e_field` — applied field (V m⁻¹), `mu` — dynamic viscosity (Pa·s).
pub fn electroosmotic_velocity(eps_r: f64, zeta: f64, e_field: f64, mu: f64) -> f64 {
    -eps_r * EPSILON_0 * zeta * e_field / mu
}

/// Electrophoretic mobility μ_ep = z e / (6π μ r_p).
///
/// `z` — particle charge number, `r_p` — particle radius (m), `mu` — dynamic viscosity.
pub fn electrophoretic_mobility(z: f64, r_p: f64, mu: f64) -> f64 {
    z * ELEM_CHARGE / (6.0 * PI * mu * r_p)
}

/// Maxwell slip velocity at a wall.
///
/// u_slip = (2 - σ_v) / σ_v · Kn · (∂u/∂n) · L
///
/// `sigma_v` — tangential momentum accommodation coefficient (0..1],
/// `kn` — Knudsen number, `du_dn` — velocity gradient at wall (s⁻¹), `length` — channel half-height (m).
pub fn maxwell_slip_velocity(sigma_v: f64, kn: f64, du_dn: f64, length: f64) -> f64 {
    (2.0 - sigma_v) / sigma_v * kn * du_dn * length
}

/// Poiseuille velocity profile with Maxwell slip.
///
/// For a channel of half-width `h`, pressure gradient `dp_dx` and viscosity `mu`.
/// Returns centreline velocity (m s⁻¹).
pub fn poiseuille_slip_centerline(dp_dx: f64, h: f64, mu: f64, kn: f64, sigma_v: f64) -> f64 {
    let slip_coeff = (2.0 - sigma_v) / sigma_v * kn;
    -dp_dx / (2.0 * mu) * h * h * (1.0 + 2.0 * slip_coeff)
}

/// Dean number De = Re √(D_h / (2 R_c)).
///
/// `re` — Reynolds number, `dh` — hydraulic diameter (m), `rc` — centreline radius of curvature (m).
pub fn dean_number(re: f64, dh: f64, rc: f64) -> f64 {
    re * (dh / (2.0 * rc)).sqrt()
}

/// Critical Dean number above which secondary Dean vortices develop.
/// Empirical: De_c ≈ 11.6 for circular cross-section.
pub const DEAN_CRITICAL: f64 = 11.6;

/// Capillary electrophoresis plate count (theoretical plates).
///
/// N = μ_ep V_app / (2 D_mol)
///
/// `mu_ep` — electrophoretic mobility (m² V⁻¹ s⁻¹), `v_app` — applied voltage (V),
/// `d_mol` — molecular diffusivity (m² s⁻¹).
pub fn capillary_electrophoresis_plates(mu_ep: f64, v_app: f64, d_mol: f64) -> f64 {
    mu_ep * v_app / (2.0 * d_mol)
}

/// Diffusiophoretic velocity driven by a concentration gradient.
///
/// u_dp = Γ · ∇ ln c
///
/// `gamma` — diffusiophoretic coefficient (m² s⁻¹), `dc_dx` — concentration gradient (mol m⁻⁴),
/// `c` — local concentration (mol m⁻³).
pub fn diffusiophoretic_velocity(gamma: f64, dc_dx: f64, c: f64) -> f64 {
    if c <= 0.0 {
        return 0.0;
    }
    gamma * dc_dx / c
}

/// Electrowetting contact angle via Lippmann-Young equation.
///
/// cos θ(V) = cos θ₀ + ε_r ε₀ V² / (2 γ d)
///
/// `theta0` — equilibrium contact angle (rad), `eps_r` — dielectric constant,
/// `voltage` — applied voltage (V), `gamma` — surface tension (N m⁻¹),
/// `thickness` — dielectric thickness (m).
pub fn ewod_contact_angle(
    theta0: f64,
    eps_r: f64,
    voltage: f64,
    gamma: f64,
    thickness: f64,
) -> f64 {
    let cos_theta =
        theta0.cos() + eps_r * EPSILON_0 * voltage * voltage / (2.0 * gamma * thickness);
    cos_theta.clamp(-1.0, 1.0).acos()
}

/// T-junction droplet formation: droplet length scaling.
///
/// L_drop / w_c ≈ 1 + q_d / q_c
///
/// `w_c` — continuous phase channel width (m), `q_d` — dispersed phase flow rate (m³ s⁻¹),
/// `q_c` — continuous phase flow rate (m³ s⁻¹).
pub fn t_junction_droplet_length(w_c: f64, q_d: f64, q_c: f64) -> f64 {
    w_c * (1.0 + q_d / q_c)
}

/// Flow focusing droplet diameter scaling.
///
/// d_drop ≈ w_c · (μ_d / μ_c)^(1/3) · (q_d / q_c)^(1/3)
///
/// `w_c` — orifice width (m), `mu_d`, `mu_c` — dispersed/continuous viscosities,
/// `q_d`, `q_c` — flow rates (m³ s⁻¹).
pub fn flow_focusing_droplet_diameter(w_c: f64, mu_d: f64, mu_c: f64, q_d: f64, q_c: f64) -> f64 {
    w_c * (mu_d / mu_c).cbrt() * (q_d / q_c).cbrt()
}

/// Poiseuille flow Darcy-Weisbach pressure drop over length L.
///
/// Δp = 128 μ L Q / (π D_h⁴)
///
/// `mu` — dynamic viscosity, `length` — channel length, `q` — volumetric flow rate (m³ s⁻¹),
/// `dh` — hydraulic diameter (m).
pub fn poiseuille_pressure_drop(mu: f64, length: f64, q: f64, dh: f64) -> f64 {
    128.0 * mu * length * q / (PI * dh.powi(4))
}

/// PDMS channel wall deflection (thin-plate model, rectangular channel).
///
/// δ_max ≈ α Δp w⁴ / (E h³)
///
/// `alpha` — geometry factor (~0.014 for a clamped rectangular plate),
/// `dp` — pressure difference (Pa), `w` — channel width (m),
/// `e_mod` — Young's modulus (Pa), `h_wall` — wall thickness (m).
pub fn pdms_wall_deflection(alpha: f64, dp: f64, w: f64, e_mod: f64, h_wall: f64) -> f64 {
    alpha * dp * w.powi(4) / (e_mod * h_wall.powi(3))
}

// ---------------------------------------------------------------------------
// KnudsenFlowRegime
// ---------------------------------------------------------------------------

/// Flow regime classification based on Knudsen number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KnudsenFlowRegime {
    /// Continuum (Kn < 0.001).
    Continuum,
    /// Slip flow (0.001 ≤ Kn < 0.1).
    SlipFlow,
    /// Transition regime (0.1 ≤ Kn < 10).
    Transition,
    /// Free molecular (Kn ≥ 10).
    FreeMolecular,
}

impl KnudsenFlowRegime {
    /// Classify the Knudsen number into a flow regime.
    pub fn classify(kn: f64) -> Self {
        if kn < 0.001 {
            Self::Continuum
        } else if kn < 0.1 {
            Self::SlipFlow
        } else if kn < 10.0 {
            Self::Transition
        } else {
            Self::FreeMolecular
        }
    }

    /// Return effective viscosity correction factor for BGK-LBM slip regime.
    ///
    /// Based on Beskok-Karniadakis unified model.
    pub fn effective_viscosity_factor(&self, kn: f64) -> f64 {
        let alpha = 1.0; // rarefaction coefficient (simplified)
        1.0 + alpha * kn
    }
}

// ---------------------------------------------------------------------------
// MaxwellSlipBC
// ---------------------------------------------------------------------------

/// Maxwell slip boundary condition parameters for rarefied gas LBM.
#[derive(Debug, Clone)]
pub struct MaxwellSlipBC {
    /// Tangential momentum accommodation coefficient σ_v ∈ (0, 1].
    pub sigma_v: f64,
    /// Thermal accommodation coefficient σ_T ∈ (0, 1].
    pub sigma_t: f64,
    /// Knudsen number at this boundary.
    pub knudsen: f64,
    /// Wall temperature (K).
    pub temp_wall: f64,
    /// Fluid temperature (K).
    pub temp_fluid: f64,
}

impl MaxwellSlipBC {
    /// Create a new Maxwell slip boundary condition.
    pub fn new(sigma_v: f64, sigma_t: f64, knudsen: f64, temp_wall: f64, temp_fluid: f64) -> Self {
        Self {
            sigma_v,
            sigma_t,
            knudsen,
            temp_wall,
            temp_fluid,
        }
    }

    /// First-order Maxwell slip velocity at the wall given bulk velocity gradient.
    ///
    /// u_slip = (2 - σ_v)/σ_v · λ · (∂u/∂y)|_wall
    pub fn slip_velocity(&self, du_dn: f64, mean_free_path: f64) -> f64 {
        (2.0 - self.sigma_v) / self.sigma_v * mean_free_path * du_dn
    }

    /// Second-order Maxwell slip (including thermal creep correction).
    pub fn second_order_slip(&self, du_dn: f64, d2u_dn2: f64, mean_free_path: f64) -> f64 {
        let first = (2.0 - self.sigma_v) / self.sigma_v * mean_free_path * du_dn;
        let second = -(2.0 - self.sigma_v) / self.sigma_v * mean_free_path.powi(2) * d2u_dn2 / 2.0;
        first + second
    }

    /// Temperature jump at the wall (Smoluchowski formula).
    pub fn temperature_jump(
        &self,
        dt_dn: f64,
        mean_free_path: f64,
        gamma: f64,
        prandtl: f64,
    ) -> f64 {
        (2.0 - self.sigma_t) / self.sigma_t * 2.0 * gamma / (gamma + 1.0) * mean_free_path / prandtl
            * dt_dn
    }

    /// LBM slip relaxation parameter ω_slip for effective-viscosity approach.
    ///
    /// τ_slip = τ_bulk + C_Kn · Kn where C_Kn is the Knudsen layer correction.
    pub fn slip_relaxation_omega(&self, omega_bulk: f64) -> f64 {
        let tau_bulk = 1.0 / omega_bulk;
        let c_kn = (2.0 - self.sigma_v) / self.sigma_v * self.knudsen;
        let tau_slip = tau_bulk + c_kn;
        1.0 / tau_slip
    }
}

// ---------------------------------------------------------------------------
// ElectroosmosisLBM
// ---------------------------------------------------------------------------

/// Electroosmotic flow simulation state for a 1D channel.
#[derive(Debug, Clone)]
pub struct ElectroosmosisLBM {
    /// Number of lattice nodes across the channel.
    pub ny: usize,
    /// Relative permittivity of the electrolyte.
    pub eps_r: f64,
    /// Zeta potential at the walls (V).
    pub zeta: f64,
    /// Applied external electric field E_x (V m⁻¹).
    pub e_field: f64,
    /// Dynamic viscosity (Pa·s).
    pub mu: f64,
    /// Debye length (lattice units, normalized by channel half-height).
    pub debye_lu: f64,
    /// Electric potential at each node (V), normalized.
    pub phi: Vec<f64>,
    /// Velocity at each node (lattice units).
    pub u: Vec<f64>,
    /// Distribution functions f\[node\]\[direction\].
    pub f: Vec<[f64; 9]>,
}

impl ElectroosmosisLBM {
    /// Create a new `ElectroosmosisLBM` simulation.
    ///
    /// `ny` — nodes across channel, `eps_r` — relative permittivity,
    /// `zeta` — wall zeta potential (V), `e_field` — applied field (V m⁻¹),
    /// `mu` — viscosity, `debye_lu` — Debye length in lattice units.
    pub fn new(ny: usize, eps_r: f64, zeta: f64, e_field: f64, mu: f64, debye_lu: f64) -> Self {
        let phi = Self::init_debye_potential(ny, zeta, debye_lu);
        let u = vec![0.0; ny];
        let f = vec![[0.0; 9]; ny];
        Self {
            ny,
            eps_r,
            zeta,
            e_field,
            mu,
            debye_lu,
            phi,
            u,
            f,
        }
    }

    /// Initialize Debye-Hückel potential profile across the channel.
    ///
    /// φ(y) = ζ · cosh((H/2 - y)/λ_D) / cosh(H/(2λ_D))
    fn init_debye_potential(ny: usize, zeta: f64, debye_lu: f64) -> Vec<f64> {
        let h = ny as f64;
        (0..ny)
            .map(|j| {
                let y = j as f64 + 0.5;
                let arg = (h / 2.0 - y) / debye_lu;
                let norm = (h / (2.0 * debye_lu)).cosh();
                zeta * arg.cosh() / norm
            })
            .collect()
    }

    /// Compute equilibrium distribution for D2Q9 (1D slice, ux only).
    fn feq_1d(&self, rho: f64, ux: f64, i: usize) -> f64 {
        let cx = D2Q9_C[i][0];
        let cu = cx * ux;
        D2Q9_W[i] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - ux * ux / (2.0 * CS2))
    }

    /// Compute electroosmotic body force at node j.
    ///
    /// F_eo = -ε_r ε₀ ∂²φ/∂y² · E_x (Poisson–Boltzmann source).
    /// Approximated as charge density ρ_e · E_x.
    pub fn eo_force(&self, j: usize) -> f64 {
        if j == 0 || j + 1 >= self.ny {
            return 0.0;
        }
        // ρ_e ≈ -ε_r ε₀ ∂²φ/∂y² (Poisson equation)
        let d2phi = self.phi[j + 1] - 2.0 * self.phi[j] + self.phi[j - 1];
        -self.eps_r * EPSILON_0 * d2phi * self.e_field
    }

    /// Run one BGK collision+streaming step with electroosmotic forcing.
    ///
    /// `omega` — relaxation frequency (1/τ).
    pub fn step(&mut self, omega: f64) {
        // Collision with forcing (Guo et al. scheme simplified to 1D)
        for j in 0..self.ny {
            let rho = self.f[j].iter().sum::<f64>().max(1e-10);
            let mut ux = 0.0;
            for i in 0..9 {
                ux += D2Q9_C[i][0] * self.f[j][i];
            }
            ux = (ux + 0.5 * self.eo_force(j)) / rho;
            self.u[j] = ux;
            for i in 0..9 {
                let feq = self.feq_1d(rho, ux, i);
                // Guo body force term
                let force_term =
                    (1.0 - omega * 0.5) * D2Q9_W[i] * D2Q9_C[i][0] * self.eo_force(j) / CS2;
                self.f[j][i] += -omega * (self.f[j][i] - feq) + force_term;
            }
        }
        // Periodic streaming in x (1D: we just update u directly)
        // For a full 2D LBM, streaming would propagate along cx, cy.
    }

    /// Analytic Helmholtz-Smoluchowski electro-osmotic velocity.
    pub fn hs_velocity(&self) -> f64 {
        -self.eps_r * EPSILON_0 * self.zeta * self.e_field / self.mu
    }

    /// Initialize with equilibrium distributions.
    pub fn init_equilibrium(&mut self) {
        for j in 0..self.ny {
            for i in 0..9 {
                self.f[j][i] = self.feq_1d(1.0, 0.0, i);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ElectrophoreticParticle
// ---------------------------------------------------------------------------

/// A charged particle undergoing electrophoresis in a microfluidic channel.
#[derive(Debug, Clone)]
pub struct ElectrophoreticParticle {
    /// Particle position \[x, y\] (m).
    pub pos: [f64; 2],
    /// Particle velocity \[vx, vy\] (m s⁻¹).
    pub vel: [f64; 2],
    /// Particle radius (m).
    pub radius: f64,
    /// Surface charge number z (signed).
    pub charge_z: f64,
    /// Zeta potential of the particle surface (V).
    pub zeta: f64,
    /// Dynamic viscosity of the carrier fluid.
    pub mu: f64,
}

impl ElectrophoreticParticle {
    /// Create a new `ElectrophoreticParticle`.
    pub fn new(pos: [f64; 2], radius: f64, charge_z: f64, zeta: f64, mu: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 2],
            radius,
            charge_z,
            zeta,
            mu,
        }
    }

    /// Electrophoretic mobility μ_ep = ε_r ε₀ ζ / μ (Smoluchowski limit, κ a → ∞).
    pub fn smoluchowski_mobility(&self, eps_r: f64) -> f64 {
        eps_r * EPSILON_0 * self.zeta / self.mu
    }

    /// Henry factor f(κ a) for intermediate κ a (Henry correction).
    ///
    /// Interpolation between Hückel (f=1) and Smoluchowski (f=1.5) limits.
    pub fn henry_factor(&self, kappa_a: f64) -> f64 {
        // Henry's function approximation: f(κa) ≈ 1 + (κa)² / (2(1 + (κa)))
        1.0 + kappa_a * kappa_a / (2.0 * (1.0 + kappa_a))
    }

    /// Electrophoretic velocity in field E_x (m s⁻¹).
    pub fn electrophoretic_velocity(&self, eps_r: f64, e_x: f64) -> f64 {
        self.smoluchowski_mobility(eps_r) * e_x
    }

    /// Stokes drag force on the particle given relative fluid velocity.
    pub fn stokes_drag(&self, fluid_vel: [f64; 2]) -> [f64; 2] {
        let factor = 6.0 * PI * self.mu * self.radius;
        [
            factor * (fluid_vel[0] - self.vel[0]),
            factor * (fluid_vel[1] - self.vel[1]),
        ]
    }

    /// Advance particle position by time step `dt` under applied field.
    pub fn advance(&mut self, eps_r: f64, e_x: f64, fluid_vel: [f64; 2], mass: f64, dt: f64) {
        let ep_vel = self.electrophoretic_velocity(eps_r, e_x);
        let drag = self.stokes_drag(fluid_vel);
        let ax = (drag[0]) / mass;
        let ay = (drag[1]) / mass;
        self.vel[0] += ax * dt;
        self.vel[1] += ay * dt;
        // Superimpose EP velocity directly (overdamped limit)
        self.pos[0] += (self.vel[0] + ep_vel) * dt;
        self.pos[1] += self.vel[1] * dt;
    }
}

// ---------------------------------------------------------------------------
// TJunctionDropletGenerator
// ---------------------------------------------------------------------------

/// T-junction droplet formation model.
#[derive(Debug, Clone)]
pub struct TJunctionDropletGenerator {
    /// Width of the continuous phase channel (m).
    pub width_c: f64,
    /// Width of the dispersed phase channel (m).
    pub width_d: f64,
    /// Dynamic viscosity of continuous phase (Pa·s).
    pub mu_c: f64,
    /// Dynamic viscosity of dispersed phase (Pa·s).
    pub mu_d: f64,
    /// Interfacial tension (N m⁻¹).
    pub sigma: f64,
    /// Continuous phase volumetric flow rate (m³ s⁻¹).
    pub q_c: f64,
    /// Dispersed phase volumetric flow rate (m³ s⁻¹).
    pub q_d: f64,
}

impl TJunctionDropletGenerator {
    /// Create a new `TJunctionDropletGenerator`.
    pub fn new(
        width_c: f64,
        width_d: f64,
        mu_c: f64,
        mu_d: f64,
        sigma: f64,
        q_c: f64,
        q_d: f64,
    ) -> Self {
        Self {
            width_c,
            width_d,
            mu_c,
            mu_d,
            sigma,
            q_c,
            q_d,
        }
    }

    /// Capillary number of the continuous phase Ca_c = μ_c U_c / σ.
    pub fn capillary_number(&self, height: f64) -> f64 {
        let u_c = self.q_c / (self.width_c * height);
        self.mu_c * u_c / self.sigma
    }

    /// Viscosity ratio λ = μ_d / μ_c.
    pub fn viscosity_ratio(&self) -> f64 {
        self.mu_d / self.mu_c
    }

    /// Flow rate ratio φ = q_d / q_c.
    pub fn flow_ratio(&self) -> f64 {
        self.q_d / self.q_c
    }

    /// Droplet length L_drop ≈ w_c (1 + q_d/q_c) — squeezing regime.
    pub fn droplet_length_squeezing(&self) -> f64 {
        self.width_c * (1.0 + self.flow_ratio())
    }

    /// Droplet length in dripping regime (Garstecki et al. scaling).
    ///
    /// L_drop / w_c = α + β · q_d / q_c, with empirical α≈1.0, β≈1.0.
    pub fn droplet_length_dripping(&self, alpha: f64, beta: f64) -> f64 {
        self.width_c * (alpha + beta * self.flow_ratio())
    }

    /// Droplet generation frequency f_drop = q_d / V_drop.
    pub fn generation_frequency(&self, droplet_volume: f64) -> f64 {
        self.q_d / droplet_volume
    }

    /// Droplet volume (sphere approximation from radius).
    pub fn droplet_volume_sphere(&self, radius: f64) -> f64 {
        4.0 / 3.0 * PI * radius.powi(3)
    }

    /// Weber number We = ρ U² L / σ for continuous phase.
    pub fn weber_number(&self, rho_c: f64, height: f64) -> f64 {
        let u_c = self.q_c / (self.width_c * height);
        rho_c * u_c * u_c * self.width_c / self.sigma
    }
}

// ---------------------------------------------------------------------------
// FlowFocusingDroplet
// ---------------------------------------------------------------------------

/// Flow-focusing droplet/bubble generator model.
#[derive(Debug, Clone)]
pub struct FlowFocusingDroplet {
    /// Orifice (nozzle) width w_o (m).
    pub width_orifice: f64,
    /// Continuous phase viscosity (Pa·s).
    pub mu_c: f64,
    /// Dispersed phase viscosity (Pa·s).
    pub mu_d: f64,
    /// Continuous phase flow rate (m³ s⁻¹).
    pub q_c: f64,
    /// Dispersed phase flow rate (m³ s⁻¹).
    pub q_d: f64,
    /// Interfacial tension (N m⁻¹).
    pub sigma: f64,
}

impl FlowFocusingDroplet {
    /// Create a new `FlowFocusingDroplet`.
    pub fn new(width_orifice: f64, mu_c: f64, mu_d: f64, q_c: f64, q_d: f64, sigma: f64) -> Self {
        Self {
            width_orifice,
            mu_c,
            mu_d,
            q_c,
            q_d,
            sigma,
        }
    }

    /// Droplet diameter scaling: d ≈ w_o (μ_d/μ_c)^(1/3) (q_d/q_c)^(1/3).
    pub fn droplet_diameter(&self) -> f64 {
        self.width_orifice * (self.mu_d / self.mu_c).cbrt() * (self.q_d / self.q_c).cbrt()
    }

    /// Capillary number Ca = μ_c U_c / σ at the orifice (h = w_o assumed).
    pub fn capillary_number(&self) -> f64 {
        let u_c = self.q_c / (self.width_orifice * self.width_orifice);
        self.mu_c * u_c / self.sigma
    }

    /// Flow ratio φ = q_d / q_c.
    pub fn flow_ratio(&self) -> f64 {
        self.q_d / self.q_c
    }

    /// Critical capillary number for transition dripping→jetting (~0.1).
    pub fn is_jetting(&self) -> bool {
        self.capillary_number() > 0.1
    }
}

// ---------------------------------------------------------------------------
// EwodDroplet (digital microfluidics)
// ---------------------------------------------------------------------------

/// Droplet in a digital microfluidics EWOD device.
#[derive(Debug, Clone)]
pub struct EwodDroplet {
    /// Droplet position \[x, y\] (m).
    pub pos: [f64; 2],
    /// Droplet radius (m).
    pub radius: f64,
    /// Contact angle (rad) — current value.
    pub contact_angle: f64,
    /// Equilibrium (zero-voltage) contact angle (rad).
    pub theta0: f64,
    /// Surface tension liquid-vapor (N m⁻¹).
    pub gamma_lv: f64,
}

impl EwodDroplet {
    /// Create a new `EwodDroplet`.
    pub fn new(pos: [f64; 2], radius: f64, theta0: f64, gamma_lv: f64) -> Self {
        Self {
            pos,
            radius,
            contact_angle: theta0,
            theta0,
            gamma_lv,
        }
    }

    /// Apply voltage via Lippmann-Young equation.
    ///
    /// cos θ(V) = cos θ₀ + ε_r ε₀ V² / (2 γ d)
    pub fn apply_voltage(&mut self, voltage: f64, eps_r: f64, dielectric_thickness: f64) {
        self.contact_angle = ewod_contact_angle(
            self.theta0,
            eps_r,
            voltage,
            self.gamma_lv,
            dielectric_thickness,
        );
    }

    /// EWOD driving pressure: Δp = 2 γ cos θ / R (Young-Laplace for thin gap).
    pub fn driving_pressure(&self, gap_height: f64) -> f64 {
        2.0 * self.gamma_lv * self.contact_angle.cos() / gap_height
    }

    /// Droplet volume (thin disk approximation, radius r, height h).
    pub fn volume_thin_disk(&self, gap_height: f64) -> f64 {
        PI * self.radius * self.radius * gap_height
    }

    /// Check if droplet will move toward a higher-voltage electrode.
    pub fn will_actuate(&self, theta_adjacent: f64) -> bool {
        theta_adjacent < self.contact_angle
    }
}

// ---------------------------------------------------------------------------
// CapillaryElectrophoresisSeparator
// ---------------------------------------------------------------------------

/// Capillary electrophoresis separator model.
#[derive(Debug, Clone)]
pub struct CapillaryElectrophoresisSeparator {
    /// Capillary length (m).
    pub length: f64,
    /// Applied voltage (V).
    pub voltage: f64,
    /// Electroosmotic mobility μ_eo (m² V⁻¹ s⁻¹).
    pub mu_eo: f64,
    /// Background electrolyte dynamic viscosity (Pa·s).
    pub mu_fluid: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Capillary inner diameter (m).
    pub diameter: f64,
}

impl CapillaryElectrophoresisSeparator {
    /// Create a new `CapillaryElectrophoresisSeparator`.
    pub fn new(
        length: f64,
        voltage: f64,
        mu_eo: f64,
        mu_fluid: f64,
        temperature: f64,
        diameter: f64,
    ) -> Self {
        Self {
            length,
            voltage,
            mu_eo,
            mu_fluid,
            temperature,
            diameter,
        }
    }

    /// Applied electric field E = V / L (V m⁻¹).
    pub fn electric_field(&self) -> f64 {
        self.voltage / self.length
    }

    /// Joule heating power per unit volume P = σ E² (W m⁻³).
    ///
    /// `conductivity` — electrolyte conductivity (S m⁻¹).
    pub fn joule_heating(&self, conductivity: f64) -> f64 {
        let e = self.electric_field();
        conductivity * e * e
    }

    /// Migration time for analyte with electrophoretic mobility μ_ep.
    ///
    /// t = L / ((μ_ep + μ_eo) · E)
    pub fn migration_time(&self, mu_ep: f64) -> f64 {
        let e = self.electric_field();
        self.length / ((mu_ep + self.mu_eo) * e)
    }

    /// Resolution between two analytes.
    ///
    /// R = (t2 - t1) / (0.5 (σ_t1 + σ_t2))
    ///
    /// `mu_ep1`, `mu_ep2` — electrophoretic mobilities, `d_mol` — diffusivity (m² s⁻¹).
    pub fn resolution(&self, mu_ep1: f64, mu_ep2: f64, d_mol: f64) -> f64 {
        let t1 = self.migration_time(mu_ep1);
        let t2 = self.migration_time(mu_ep2);
        let sigma_t = (2.0 * d_mol * t1).sqrt(); // width ≈ std dev of Gaussian zone
        (t2 - t1).abs() / (sigma_t)
    }

    /// Theoretical plate count N = μ_app V / (2 D_mol).
    ///
    /// `mu_app` — apparent mobility = μ_ep + μ_eo, `d_mol` — diffusivity.
    pub fn theoretical_plates(&self, mu_ep: f64, d_mol: f64) -> f64 {
        let mu_app = mu_ep + self.mu_eo;
        mu_app * self.voltage / (2.0 * d_mol)
    }

    /// Plate height H = L / N.
    pub fn plate_height(&self, mu_ep: f64, d_mol: f64) -> f64 {
        let n = self.theoretical_plates(mu_ep, d_mol);
        if n <= 0.0 {
            f64::INFINITY
        } else {
            self.length / n
        }
    }
}

// ---------------------------------------------------------------------------
// DiffusiophoresisSolver
// ---------------------------------------------------------------------------

/// Diffusiophoresis solver: particle motion driven by solute gradients.
#[derive(Debug, Clone)]
pub struct DiffusiophoresisSolver {
    /// Diffusiophoretic coefficient Γ (m² s⁻¹).
    pub gamma_dp: f64,
    /// Number of nodes in 1D domain.
    pub nx: usize,
    /// Grid spacing (m).
    pub dx: f64,
    /// Solute concentration at each node (mol m⁻³).
    pub concentration: Vec<f64>,
    /// Particle velocity at each node (m s⁻¹) from diffusiophoresis.
    pub dp_velocity: Vec<f64>,
}

impl DiffusiophoresisSolver {
    /// Create a new `DiffusiophoresisSolver` with a linear concentration profile.
    pub fn new(gamma_dp: f64, nx: usize, dx: f64, c_left: f64, c_right: f64) -> Self {
        let concentration = (0..nx)
            .map(|i| c_left + (c_right - c_left) * i as f64 / (nx as f64 - 1.0).max(1.0))
            .collect();
        let dp_velocity = vec![0.0; nx];
        Self {
            gamma_dp,
            nx,
            dx,
            concentration,
            dp_velocity,
        }
    }

    /// Update diffusiophoretic velocity field from concentration gradient.
    ///
    /// u_dp(x) = Γ · d ln c / dx ≈ Γ / c · dc/dx
    pub fn update_dp_velocity(&mut self) {
        for i in 0..self.nx {
            let c = self.concentration[i].max(1e-30);
            let dc_dx = if i == 0 {
                (self.concentration[1] - self.concentration[0]) / self.dx
            } else if i + 1 == self.nx {
                (self.concentration[self.nx - 1] - self.concentration[self.nx - 2]) / self.dx
            } else {
                (self.concentration[i + 1] - self.concentration[i - 1]) / (2.0 * self.dx)
            };
            self.dp_velocity[i] = self.gamma_dp * dc_dx / c;
        }
    }

    /// Diffuse concentration by one step (explicit finite difference).
    ///
    /// `d_solute` — solute diffusivity (m² s⁻¹), `dt` — time step (s).
    pub fn diffuse_concentration(&mut self, d_solute: f64, dt: f64) {
        let mut c_new = self.concentration.clone();
        let r = d_solute * dt / (self.dx * self.dx);
        for i in 1..self.nx - 1 {
            c_new[i] = self.concentration[i]
                + r * (self.concentration[i + 1] - 2.0 * self.concentration[i]
                    + self.concentration[i - 1]);
        }
        self.concentration = c_new;
    }

    /// Maximum diffusiophoretic speed in the domain.
    pub fn max_dp_speed(&self) -> f64 {
        self.dp_velocity
            .iter()
            .cloned()
            .fold(0.0_f64, |a, v| a.max(v.abs()))
    }
}

// ---------------------------------------------------------------------------
// DeanFlowChannel
// ---------------------------------------------------------------------------

/// Dean flow in a curved rectangular microchannel.
#[derive(Debug, Clone)]
pub struct DeanFlowChannel {
    /// Channel width (m).
    pub width: f64,
    /// Channel height (m).
    pub height: f64,
    /// Centreline radius of curvature (m).
    pub radius_curve: f64,
    /// Fluid dynamic viscosity (Pa·s).
    pub mu: f64,
    /// Fluid density (kg m⁻³).
    pub rho: f64,
    /// Mean axial velocity (m s⁻¹).
    pub u_mean: f64,
}

impl DeanFlowChannel {
    /// Create a new `DeanFlowChannel`.
    pub fn new(width: f64, height: f64, radius_curve: f64, mu: f64, rho: f64, u_mean: f64) -> Self {
        Self {
            width,
            height,
            radius_curve,
            mu,
            rho,
            u_mean,
        }
    }

    /// Hydraulic diameter D_h = 4A/P for rectangular cross-section.
    pub fn hydraulic_diameter(&self) -> f64 {
        4.0 * self.width * self.height / (2.0 * (self.width + self.height))
    }

    /// Reynolds number Re = ρ U D_h / μ.
    pub fn reynolds_number(&self) -> f64 {
        self.rho * self.u_mean * self.hydraulic_diameter() / self.mu
    }

    /// Dean number De = Re √(D_h / (2 R_c)).
    pub fn dean_number(&self) -> f64 {
        let re = self.reynolds_number();
        let dh = self.hydraulic_diameter();
        dean_number(re, dh, self.radius_curve)
    }

    /// Secondary Dean flow velocity (centrifugal estimate).
    ///
    /// U_Dean ≈ U_axial · (D_h / (2 R_c))^0.5
    pub fn dean_velocity(&self) -> f64 {
        self.u_mean * (self.hydraulic_diameter() / (2.0 * self.radius_curve)).sqrt()
    }

    /// Inertial focusing position: equilibrium distance from channel axis.
    ///
    /// x_eq / (D_h/2) ≈ 0.6 (empirical for rectangular channel).
    pub fn inertial_focusing_position(&self) -> f64 {
        0.6 * self.hydraulic_diameter() / 2.0
    }

    /// Dean flow vortex strength (dimensionless): De/De_c.
    pub fn vortex_strength(&self) -> f64 {
        self.dean_number() / DEAN_CRITICAL
    }

    /// Particle separation factor based on Dean number difference for two sizes.
    ///
    /// Particles separate when |De_1 - De_2| > threshold.
    pub fn separation_factor(&self, r_p1: f64, r_p2: f64) -> f64 {
        // Confinement ratio λ = r_p / D_h/2
        let dh = self.hydraulic_diameter();
        let lambda1 = r_p1 / (dh / 2.0);
        let lambda2 = r_p2 / (dh / 2.0);
        (lambda1 - lambda2).abs() * self.dean_number()
    }
}

// ---------------------------------------------------------------------------
// PdmsChannelDeformation
// ---------------------------------------------------------------------------

/// PDMS (polydimethylsiloxane) channel wall deformation coupling model.
#[derive(Debug, Clone)]
pub struct PdmsChannelDeformation {
    /// Channel width (m).
    pub width: f64,
    /// Channel height (undeformed, m).
    pub height: f64,
    /// Channel length (m).
    pub length: f64,
    /// PDMS Young's modulus (Pa, typically ~1–3 MPa).
    pub young_modulus: f64,
    /// Poisson's ratio of PDMS (~0.5).
    pub poisson: f64,
    /// PDMS wall thickness (m).
    pub wall_thickness: f64,
    /// Internal gauge pressure (Pa).
    pub pressure: f64,
}

impl PdmsChannelDeformation {
    /// Create a new `PdmsChannelDeformation`.
    pub fn new(
        width: f64,
        height: f64,
        length: f64,
        young_modulus: f64,
        poisson: f64,
        wall_thickness: f64,
        pressure: f64,
    ) -> Self {
        Self {
            width,
            height,
            length,
            young_modulus,
            poisson,
            wall_thickness,
            pressure,
        }
    }

    /// Maximum wall deflection (thin-plate clamped model).
    ///
    /// δ_max = 0.014 Δp w⁴ / (E h³) (rectangular, clamped edges).
    pub fn max_deflection(&self) -> f64 {
        pdms_wall_deflection(
            0.014,
            self.pressure,
            self.width,
            self.young_modulus,
            self.wall_thickness,
        )
    }

    /// Effective channel height after deformation (approx).
    pub fn effective_height(&self) -> f64 {
        self.height + self.max_deflection()
    }

    /// Flow resistance correction factor due to deformation.
    ///
    /// R_eff / R_0 ≈ (h_0 / h_eff)³  (Hagen-Poiseuille scaling).
    pub fn flow_resistance_factor(&self) -> f64 {
        let h_eff = self.effective_height();
        (self.height / h_eff).powi(3)
    }

    /// Axial strain in PDMS wall from internal pressure (hoop-like).
    ///
    /// ε_x ≈ (1 - ν²) Δp w / (2 E h_wall)
    pub fn axial_strain(&self) -> f64 {
        (1.0 - self.poisson * self.poisson) * self.pressure * self.width
            / (2.0 * self.young_modulus * self.wall_thickness)
    }

    /// Natural frequency of PDMS wall (clamped plate, fundamental mode).
    ///
    /// f_n = (π/2) √(E h²_wall / (12(1-ν²) ρ w⁴))
    ///
    /// `rho_pdms` — density of PDMS (kg m⁻³, ~970).
    pub fn natural_frequency(&self, rho_pdms: f64) -> f64 {
        let factor = self.young_modulus * self.wall_thickness.powi(2)
            / (12.0 * (1.0 - self.poisson.powi(2)) * rho_pdms * self.width.powi(4));
        (PI / 2.0) * factor.sqrt()
    }
}

// ---------------------------------------------------------------------------
// MicrofluidicLBMGrid
// ---------------------------------------------------------------------------

/// Simplified 2D LBM grid for microfluidic channel simulations.
///
/// Uses D2Q9 with BGK collision and basic electroosmotic forcing.
#[derive(Debug, Clone)]
pub struct MicrofluidicLBMGrid {
    /// Grid width (nodes).
    pub nx: usize,
    /// Grid height (nodes).
    pub ny: usize,
    /// BGK relaxation rate ω = 1/τ.
    pub omega: f64,
    /// External body force \[fx, fy\] (lattice units).
    pub body_force: [f64; 2],
    /// Distribution functions f\[x*ny + y\]\[dir\].
    pub f: Vec<[f64; 9]>,
    /// Density at each node.
    pub rho: Vec<f64>,
    /// Velocity at each node \[ux, uy\].
    pub vel: Vec<[f64; 2]>,
}

impl MicrofluidicLBMGrid {
    /// Create a new `MicrofluidicLBMGrid`.
    pub fn new(nx: usize, ny: usize, omega: f64, body_force: [f64; 2]) -> Self {
        let n = nx * ny;
        let mut f = vec![[0.0; 9]; n];
        let rho = vec![1.0; n];
        let vel = vec![[0.0; 2]; n];
        // Initialize to equilibrium at rest
        for cell in f.iter_mut() {
            for (i, fi) in cell.iter_mut().enumerate() {
                *fi = D2Q9_W[i];
            }
        }
        Self {
            nx,
            ny,
            omega,
            body_force,
            f,
            rho,
            vel,
        }
    }

    /// Index helper.
    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        x * self.ny + y
    }

    /// D2Q9 equilibrium distribution.
    fn feq(rho: f64, ux: f64, uy: f64, i: usize) -> f64 {
        let cx = D2Q9_C[i][0];
        let cy = D2Q9_C[i][1];
        let cu = cx * ux + cy * uy;
        D2Q9_W[i]
            * rho
            * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - (ux * ux + uy * uy) / (2.0 * CS2))
    }

    /// BGK collision step with Guo body-force scheme.
    pub fn collide(&mut self) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                let idx = self.idx(x, y);
                let rho = self.rho[idx];
                let ux = self.vel[idx][0];
                let uy = self.vel[idx][1];
                for i in 0..9 {
                    let feq = Self::feq(rho, ux, uy, i);
                    let cx = D2Q9_C[i][0];
                    let cy = D2Q9_C[i][1];
                    let fi_term = (1.0 - 0.5 * self.omega)
                        * D2Q9_W[i]
                        * ((cx - ux) * self.body_force[0] + (cy - uy) * self.body_force[1])
                        / CS2;
                    self.f[idx][i] += -self.omega * (self.f[idx][i] - feq) + fi_term;
                }
            }
        }
    }

    /// Streaming step with periodic boundary in x, no-slip (bounce-back) in y.
    pub fn stream(&mut self) {
        let mut f_new = vec![[0.0_f64; 9]; self.nx * self.ny];
        for x in 0..self.nx {
            for y in 0..self.ny {
                let idx_src = self.idx(x, y);
                for i in 0..9 {
                    let cx = D2Q9_C[i][0] as isize;
                    let cy = D2Q9_C[i][1] as isize;
                    let x_dst = ((x as isize + cx).rem_euclid(self.nx as isize)) as usize;
                    let y_dst_i = y as isize + cy;
                    if y_dst_i < 0 || y_dst_i >= self.ny as isize {
                        // Bounce-back: reverse direction, return to source cell
                        let opp = [0, 3, 4, 1, 2, 7, 8, 5, 6][i];
                        f_new[idx_src][opp] += self.f[idx_src][i];
                    } else {
                        let y_dst = y_dst_i as usize;
                        let idx_dst = self.idx(x_dst, y_dst);
                        f_new[idx_dst][i] = self.f[idx_src][i];
                    }
                }
            }
        }
        self.f = f_new;
    }

    /// Compute macroscopic density and velocity from distribution functions.
    pub fn update_macroscopic(&mut self) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                let idx = self.idx(x, y);
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for i in 0..9 {
                    rho += self.f[idx][i];
                    ux += D2Q9_C[i][0] * self.f[idx][i];
                    uy += D2Q9_C[i][1] * self.f[idx][i];
                }
                rho = rho.max(1e-10);
                self.rho[idx] = rho;
                // Add half body force (Guo correction)
                self.vel[idx] = [
                    (ux + 0.5 * self.body_force[0]) / rho,
                    (uy + 0.5 * self.body_force[1]) / rho,
                ];
            }
        }
    }

    /// Run `n_steps` simulation steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.collide();
            self.stream();
            self.update_macroscopic();
        }
    }

    /// Mean axial velocity (x-direction) across the channel.
    pub fn mean_axial_velocity(&self) -> f64 {
        let sum: f64 = self.vel.iter().map(|v| v[0]).sum();
        sum / (self.nx * self.ny) as f64
    }

    /// Maximum axial velocity in the channel.
    pub fn max_axial_velocity(&self) -> f64 {
        self.vel.iter().map(|v| v[0]).fold(0.0_f64, f64::max)
    }

    /// Kinematic viscosity from relaxation parameter: ν = cs² (τ - 0.5).
    pub fn kinematic_viscosity(&self) -> f64 {
        CS2 * (1.0 / self.omega - 0.5)
    }

    /// Viscous stress at a node (simplified, τ_xy component).
    pub fn stress_xy(&self, x: usize, y: usize) -> f64 {
        let idx = self.idx(x, y);
        let mut s_xy = 0.0;
        for i in 0..9 {
            let cx = D2Q9_C[i][0];
            let cy = D2Q9_C[i][1];
            let ux = self.vel[idx][0];
            let uy = self.vel[idx][1];
            let feq = Self::feq(self.rho[idx], ux, uy, i);
            s_xy += cx * cy * (self.f[idx][i] - feq);
        }
        -self.omega * s_xy
    }
}

// ---------------------------------------------------------------------------
// InertialFocusingChannel
// ---------------------------------------------------------------------------

/// Inertial microfluidics: particle focusing in straight rectangular channels.
#[derive(Debug, Clone)]
pub struct InertialFocusingChannel {
    /// Channel width (m).
    pub width: f64,
    /// Channel height (m).
    pub height: f64,
    /// Mean flow velocity (m s⁻¹).
    pub u_mean: f64,
    /// Fluid density (kg m⁻³).
    pub rho: f64,
    /// Dynamic viscosity (Pa·s).
    pub mu: f64,
}

impl InertialFocusingChannel {
    /// Create a new `InertialFocusingChannel`.
    pub fn new(width: f64, height: f64, u_mean: f64, rho: f64, mu: f64) -> Self {
        Self {
            width,
            height,
            u_mean,
            rho,
            mu,
        }
    }

    /// Channel Reynolds number.
    pub fn reynolds_number(&self) -> f64 {
        self.rho * self.u_mean * self.height / self.mu
    }

    /// Particle Reynolds number Re_p = Re (r_p / D_h)².
    pub fn particle_reynolds_number(&self, r_p: f64) -> f64 {
        let dh = 2.0 * self.width * self.height / (self.width + self.height);
        self.reynolds_number() * (r_p / dh).powi(2)
    }

    /// Inertial lift coefficient C_L (empirical, ~0.5 near channel centerline).
    pub fn lift_coefficient(&self, confinement_ratio: f64) -> f64 {
        // C_L ≈ 0.5 * λ² for confinement ratio λ = d_p / D_h
        0.5 * confinement_ratio * confinement_ratio
    }

    /// Net inertial lift force on a particle (N).
    ///
    /// F_L = ρ U² r_p⁴ / D_h² · C_L
    pub fn inertial_lift_force(&self, r_p: f64) -> f64 {
        let dh = 2.0 * self.width * self.height / (self.width + self.height);
        let cl = self.lift_coefficient(2.0 * r_p / dh);
        self.rho * self.u_mean * self.u_mean * (2.0 * r_p).powi(4) / (dh * dh) * cl
    }

    /// Equilibrium focusing position (fraction of channel height from center).
    ///
    /// Empirically, equilibrium at y_eq / (H/2) ≈ 0.6 for square channels.
    pub fn equilibrium_position_fraction(&self) -> f64 {
        0.6
    }

    /// Entry length for inertial focusing L_f = D_h³ / (Re * ν * r_p² * C_L).
    ///
    /// Equivalently: D_h² * ν / (U * r_p² * C_L).
    pub fn focusing_entry_length(&self, r_p: f64) -> f64 {
        let nu = self.mu / self.rho;
        let dh = 2.0 * self.width * self.height / (self.width + self.height);
        let cl = self.lift_coefficient(2.0 * r_p / dh).max(1e-10);
        dh * dh * nu / (self.u_mean * 4.0 * r_p * r_p * cl)
    }
}

// ---------------------------------------------------------------------------
// MicrofluidicsResult
// ---------------------------------------------------------------------------

/// Summary result from a microfluidic LBM simulation.
#[derive(Debug, Clone)]
pub struct MicrofluidicsResult {
    /// Mean velocity (lattice or SI units).
    pub mean_velocity: f64,
    /// Maximum velocity.
    pub max_velocity: f64,
    /// Flow regime classification.
    pub flow_regime: KnudsenFlowRegime,
    /// Knudsen number used.
    pub knudsen: f64,
    /// Electroosmotic contribution to mean velocity.
    pub eo_velocity: f64,
    /// Number of simulation steps performed.
    pub steps: usize,
    /// Residual (L2 norm change in velocity between last two steps).
    pub residual: f64,
}

impl MicrofluidicsResult {
    /// Create a zeroed `MicrofluidicsResult`.
    pub fn new() -> Self {
        Self {
            mean_velocity: 0.0,
            max_velocity: 0.0,
            flow_regime: KnudsenFlowRegime::Continuum,
            knudsen: 0.0,
            eo_velocity: 0.0,
            steps: 0,
            residual: 0.0,
        }
    }

    /// Is the flow converged (residual < tolerance)?
    pub fn is_converged(&self, tol: f64) -> bool {
        self.residual < tol
    }
}

impl Default for MicrofluidicsResult {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// DropletSizeDistribution
// ---------------------------------------------------------------------------

/// Statistical distribution of droplet sizes from a microfluidic generator.
#[derive(Debug, Clone)]
pub struct DropletSizeDistribution {
    /// Mean droplet diameter (m).
    pub mean_diameter: f64,
    /// Standard deviation of diameter (m).
    pub std_diameter: f64,
    /// Number of droplets sampled.
    pub count: usize,
    /// Coefficient of variation CV = σ / μ.
    pub cv: f64,
}

impl DropletSizeDistribution {
    /// Create a new `DropletSizeDistribution` from raw diameters.
    pub fn from_diameters(diameters: &[f64]) -> Self {
        let n = diameters.len();
        if n == 0 {
            return Self {
                mean_diameter: 0.0,
                std_diameter: 0.0,
                count: 0,
                cv: 0.0,
            };
        }
        let mean = diameters.iter().sum::<f64>() / n as f64;
        let var = diameters.iter().map(|&d| (d - mean).powi(2)).sum::<f64>() / n as f64;
        let std = var.sqrt();
        let cv = if mean > 0.0 { std / mean } else { 0.0 };
        Self {
            mean_diameter: mean,
            std_diameter: std,
            count: n,
            cv,
        }
    }

    /// Is the generation monodisperse? (CV < 5%).
    pub fn is_monodisperse(&self) -> bool {
        self.cv < 0.05
    }

    /// Sauter mean diameter D_32 = Σd³ / Σd².
    pub fn sauter_mean_diameter(diameters: &[f64]) -> f64 {
        let sum_d3: f64 = diameters.iter().map(|&d| d.powi(3)).sum();
        let sum_d2: f64 = diameters.iter().map(|&d| d.powi(2)).sum();
        if sum_d2 == 0.0 { 0.0 } else { sum_d3 / sum_d2 }
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- free functions ---

    #[test]
    fn test_knudsen_number() {
        let kn = knudsen_number(1e-7, 1e-5);
        assert!((kn - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_debye_length_water() {
        // 1 mM KCl in water at 298 K, ε_r ≈ 80
        let lambda_d = debye_length(80.0, 298.0, 6.022e23, 1.0);
        // Should be ~9.6 nm
        assert!(
            lambda_d > 5e-9 && lambda_d < 20e-9,
            "Debye length out of range: {}",
            lambda_d
        );
    }

    #[test]
    fn test_electroosmotic_velocity_sign() {
        // Negative zeta potential, positive field → positive EO velocity
        let u_eo = electroosmotic_velocity(80.0, -0.05, 1e4, 1e-3);
        assert!(
            u_eo > 0.0,
            "EO velocity should be positive for negative zeta + positive field"
        );
    }

    #[test]
    fn test_maxwell_slip_velocity() {
        let u_slip = maxwell_slip_velocity(1.0, 0.05, 1e4, 1e-6);
        assert!(u_slip > 0.0, "Slip velocity should be positive");
    }

    #[test]
    fn test_poiseuille_slip_centerline() {
        let u0 = poiseuille_slip_centerline(-100.0, 5e-5, 1e-3, 0.0, 1.0);
        let u_slip = poiseuille_slip_centerline(-100.0, 5e-5, 1e-3, 0.01, 1.0);
        assert!(u_slip > u0, "Slip flow should be faster than no-slip");
    }

    #[test]
    fn test_dean_number() {
        let de = dean_number(100.0, 1e-4, 1e-2);
        assert!(de > 0.0 && de < 1000.0, "Dean number out of range: {}", de);
    }

    #[test]
    fn test_ewod_contact_angle_decreases_with_voltage() {
        let theta0 = 110.0_f64.to_radians();
        let theta_v = ewod_contact_angle(theta0, 2.2, 100.0, 0.072, 1e-6);
        assert!(
            theta_v < theta0,
            "EWOD should reduce contact angle with applied voltage"
        );
    }

    #[test]
    fn test_t_junction_droplet_length() {
        let l = t_junction_droplet_length(100e-6, 0.5e-9, 2.0e-9);
        assert!(
            (l - 100e-6 * 1.25).abs() < 1e-15,
            "T-junction droplet length mismatch"
        );
    }

    #[test]
    fn test_flow_focusing_droplet_diameter() {
        let d = flow_focusing_droplet_diameter(50e-6, 1e-3, 1e-3, 1e-9, 1e-9);
        assert!(
            (d - 50e-6).abs() < 1e-15,
            "Equal viscosities and flows should give w_o"
        );
    }

    #[test]
    fn test_pdms_wall_deflection() {
        let delta = pdms_wall_deflection(0.014, 1000.0, 100e-6, 1e6, 50e-6);
        assert!(
            delta > 0.0 && delta < 1e-3,
            "Deflection out of physical range: {}",
            delta
        );
    }

    #[test]
    fn test_diffusiophoretic_velocity_zero_gradient() {
        let v = diffusiophoretic_velocity(1e-10, 0.0, 1.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_diffusiophoretic_velocity_positive() {
        let v = diffusiophoretic_velocity(1e-10, 100.0, 0.01);
        assert!(v > 0.0);
    }

    // --- KnudsenFlowRegime ---

    #[test]
    fn test_knudsen_regime_continuum() {
        assert_eq!(
            KnudsenFlowRegime::classify(0.0005),
            KnudsenFlowRegime::Continuum
        );
    }

    #[test]
    fn test_knudsen_regime_slip() {
        assert_eq!(
            KnudsenFlowRegime::classify(0.01),
            KnudsenFlowRegime::SlipFlow
        );
    }

    #[test]
    fn test_knudsen_regime_transition() {
        assert_eq!(
            KnudsenFlowRegime::classify(1.0),
            KnudsenFlowRegime::Transition
        );
    }

    #[test]
    fn test_knudsen_regime_free_molecular() {
        assert_eq!(
            KnudsenFlowRegime::classify(100.0),
            KnudsenFlowRegime::FreeMolecular
        );
    }

    #[test]
    fn test_knudsen_viscosity_factor_increases_with_kn() {
        let f1 = KnudsenFlowRegime::SlipFlow.effective_viscosity_factor(0.05);
        let f2 = KnudsenFlowRegime::Transition.effective_viscosity_factor(5.0);
        assert!(f2 > f1, "Higher Kn should give larger viscosity factor");
    }

    // --- MaxwellSlipBC ---

    #[test]
    fn test_maxwell_slip_bc_slip_velocity() {
        let bc = MaxwellSlipBC::new(1.0, 1.0, 0.01, 300.0, 300.0);
        let u_slip = bc.slip_velocity(1e4, 6.8e-8);
        assert!(u_slip > 0.0);
    }

    #[test]
    fn test_maxwell_slip_bc_temp_jump() {
        let bc = MaxwellSlipBC::new(1.0, 0.9, 0.01, 300.0, 300.0);
        let tj = bc.temperature_jump(1e5, 6.8e-8, 1.4, 0.71);
        assert!(tj > 0.0);
    }

    #[test]
    fn test_maxwell_slip_relaxation_omega() {
        let bc = MaxwellSlipBC::new(1.0, 1.0, 0.05, 300.0, 300.0);
        let omega_slip = bc.slip_relaxation_omega(1.5);
        assert!(omega_slip < 1.5, "Slip should reduce relaxation rate");
        assert!(omega_slip > 0.0);
    }

    // --- ElectroosmosisLBM ---

    #[test]
    fn test_electroosmosis_lbm_hs_velocity() {
        let eo = ElectroosmosisLBM::new(20, 80.0, -0.05, 1e4, 1e-3, 2.0);
        let u_hs = eo.hs_velocity();
        assert!(u_hs > 0.0, "HS velocity should be positive");
    }

    #[test]
    fn test_electroosmosis_lbm_init_equilibrium() {
        let mut eo = ElectroosmosisLBM::new(10, 80.0, -0.02, 5000.0, 1e-3, 1.5);
        eo.init_equilibrium();
        // Weights should sum to 1 at each node
        for j in 0..10 {
            let s: f64 = eo.f[j].iter().sum();
            assert!((s - 1.0).abs() < 1e-12, "f sum at node {j}: {s}");
        }
    }

    // --- ElectrophoreticParticle ---

    #[test]
    fn test_electrophoretic_particle_mobility() {
        let p = ElectrophoreticParticle::new([0.0, 0.0], 1e-7, 1.0, -0.05, 1e-3);
        let mob = p.smoluchowski_mobility(80.0);
        assert!(mob.abs() > 0.0);
    }

    #[test]
    fn test_electrophoretic_velocity_direction() {
        let p = ElectrophoreticParticle::new([0.0, 0.0], 1e-7, 1.0, -0.05, 1e-3);
        let vel = p.electrophoretic_velocity(80.0, 1e4);
        assert!(
            vel < 0.0,
            "Negative zeta → negative EP velocity for positive field"
        );
    }

    #[test]
    fn test_stokes_drag_opposing() {
        let p = ElectrophoreticParticle::new([0.0, 0.0], 1e-7, 1.0, -0.05, 1e-3);
        let drag = p.stokes_drag([1e-3, 0.0]); // fluid faster than particle
        assert!(
            drag[0] > 0.0,
            "Drag should push particle in fluid direction"
        );
    }

    // --- TJunctionDropletGenerator ---

    #[test]
    fn test_t_junction_flow_ratio() {
        let tjgen = TJunctionDropletGenerator::new(100e-6, 50e-6, 1e-3, 1e-2, 0.01, 2e-9, 5e-10);
        let phi = tjgen.flow_ratio();
        assert!((phi - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_t_junction_droplet_length_squeezing() {
        let tjgen = TJunctionDropletGenerator::new(100e-6, 50e-6, 1e-3, 1e-2, 0.01, 1e-9, 1e-9);
        let l = tjgen.droplet_length_squeezing();
        assert!((l - 200e-6).abs() < 1e-18, "Equal flows → L = 2 w_c: {}", l);
    }

    // --- DeanFlowChannel ---

    #[test]
    fn test_dean_flow_hydraulic_diameter() {
        let ch = DeanFlowChannel::new(100e-6, 100e-6, 5e-3, 1e-3, 1000.0, 0.01);
        let dh = ch.hydraulic_diameter();
        assert!(
            (dh - 100e-6).abs() < 1e-18,
            "Square channel D_h = side: {}",
            dh
        );
    }

    #[test]
    fn test_dean_number_positive() {
        let ch = DeanFlowChannel::new(100e-6, 100e-6, 5e-3, 1e-3, 1000.0, 0.01);
        let de = ch.dean_number();
        assert!(de > 0.0);
    }

    // --- PdmsChannelDeformation ---

    #[test]
    fn test_pdms_deflection_increases_with_pressure() {
        let ch1 = PdmsChannelDeformation::new(100e-6, 50e-6, 1e-2, 1e6, 0.5, 50e-6, 1000.0);
        let ch2 = PdmsChannelDeformation::new(100e-6, 50e-6, 1e-2, 1e6, 0.5, 50e-6, 5000.0);
        assert!(ch2.max_deflection() > ch1.max_deflection());
    }

    #[test]
    fn test_pdms_effective_height_larger() {
        let ch = PdmsChannelDeformation::new(100e-6, 50e-6, 1e-2, 1e6, 0.5, 50e-6, 1000.0);
        assert!(ch.effective_height() > ch.height);
    }

    // --- MicrofluidicLBMGrid ---

    #[test]
    fn test_lbm_grid_kinematic_viscosity() {
        let grid = MicrofluidicLBMGrid::new(10, 5, 1.0, [0.0, 0.0]);
        let nu = grid.kinematic_viscosity();
        // ω=1.0 → τ=1.0 → ν = 1/3 * 0.5 = 1/6
        assert!((nu - 1.0 / 6.0).abs() < 1e-12, "Viscosity mismatch: {}", nu);
    }

    #[test]
    fn test_lbm_grid_run_conserves_mass() {
        let mut grid = MicrofluidicLBMGrid::new(8, 8, 1.2, [1e-4, 0.0]);
        let mass_before: f64 = grid.rho.iter().sum();
        grid.run(5);
        let mass_after: f64 = grid.rho.iter().sum();
        assert!(
            (mass_after - mass_before).abs() / mass_before < 0.01,
            "Mass not conserved: before={mass_before}, after={mass_after}"
        );
    }

    #[test]
    fn test_lbm_grid_velocity_driven_by_force() {
        let mut grid = MicrofluidicLBMGrid::new(8, 8, 1.5, [1e-3, 0.0]);
        grid.run(100);
        let u_mean = grid.mean_axial_velocity();
        assert!(
            u_mean > 0.0,
            "Body force should drive positive mean velocity"
        );
    }

    // --- CapillaryElectrophoresisSeparator ---

    #[test]
    fn test_ce_electric_field() {
        let ce = CapillaryElectrophoresisSeparator::new(0.5, 25000.0, 5e-8, 1e-3, 298.0, 75e-6);
        let e = ce.electric_field();
        assert!((e - 50000.0).abs() < 0.1);
    }

    #[test]
    fn test_ce_migration_time_ordering() {
        let ce = CapillaryElectrophoresisSeparator::new(0.5, 25000.0, 5e-8, 1e-3, 298.0, 75e-6);
        let t1 = ce.migration_time(4e-8);
        let t2 = ce.migration_time(2e-8);
        assert!(t2 > t1, "Slower mobility → longer migration time");
    }

    #[test]
    fn test_ce_theoretical_plates_positive() {
        let ce = CapillaryElectrophoresisSeparator::new(0.5, 25000.0, 5e-8, 1e-3, 298.0, 75e-6);
        let n = ce.theoretical_plates(4e-8, 1e-10);
        assert!(n > 1000.0, "Should have >1000 theoretical plates: {n}");
    }

    // --- DiffusiophoresisSolver ---

    #[test]
    fn test_diffusiophoresis_gradient_direction() {
        let mut dp = DiffusiophoresisSolver::new(1e-10, 20, 1e-6, 1.0, 0.1);
        dp.update_dp_velocity();
        // Negative gradient (c_left > c_right) → negative DP velocity
        assert!(
            dp.dp_velocity[10] < 0.0,
            "DP velocity should be negative for decreasing c"
        );
    }

    #[test]
    fn test_diffusiophoresis_diffusion_step() {
        let mut dp = DiffusiophoresisSolver::new(1e-10, 10, 1e-6, 1.0, 0.0);
        let c_init = dp.concentration.clone();
        dp.diffuse_concentration(1e-12, 0.1);
        // Concentration should be smoothed
        let changed = dp
            .concentration
            .iter()
            .zip(c_init.iter())
            .any(|(a, b)| (a - b).abs() > 1e-20);
        assert!(changed, "Diffusion should change concentration profile");
    }

    // --- DropletSizeDistribution ---

    #[test]
    fn test_droplet_distribution_monodisperse() {
        let diameters = vec![100e-6; 20];
        let dist = DropletSizeDistribution::from_diameters(&diameters);
        assert!(
            dist.is_monodisperse(),
            "Identical diameters should be monodisperse"
        );
        assert!(
            dist.cv < 1e-10,
            "CV of identical diameters should be ~0, got {}",
            dist.cv
        );
    }

    #[test]
    fn test_droplet_distribution_polydisperse() {
        let diameters = vec![50e-6, 100e-6, 150e-6, 200e-6, 250e-6];
        let dist = DropletSizeDistribution::from_diameters(&diameters);
        assert!(
            !dist.is_monodisperse(),
            "Varying diameters should not be monodisperse"
        );
    }

    #[test]
    fn test_sauter_mean_diameter() {
        let d = vec![1.0, 2.0, 3.0];
        let d32 = DropletSizeDistribution::sauter_mean_diameter(&d);
        let expected = (1.0 + 8.0 + 27.0) / (1.0 + 4.0 + 9.0);
        assert!(
            (d32 - expected).abs() < 1e-12,
            "D32 mismatch: {d32} vs {expected}"
        );
    }

    // --- InertialFocusingChannel ---

    #[test]
    fn test_inertial_lift_force_positive() {
        let ch = InertialFocusingChannel::new(100e-6, 100e-6, 0.1, 1000.0, 1e-3);
        let fl = ch.inertial_lift_force(5e-6);
        assert!(fl > 0.0, "Lift force should be positive");
    }

    #[test]
    fn test_inertial_particle_re_positive() {
        let ch = InertialFocusingChannel::new(100e-6, 100e-6, 0.1, 1000.0, 1e-3);
        let re_p = ch.particle_reynolds_number(5e-6);
        assert!(
            re_p > 0.0 && re_p < 1.0,
            "Particle Re should be small: {re_p}"
        );
    }

    #[test]
    fn test_inertial_focusing_entry_length() {
        let ch = InertialFocusingChannel::new(100e-6, 100e-6, 0.1, 1000.0, 1e-3);
        let l = ch.focusing_entry_length(5e-6);
        assert!(
            l > 0.0 && l < 10.0,
            "Entry length should be in mm–cm range: {l}"
        );
    }

    // --- EwodDroplet ---

    #[test]
    fn test_ewod_droplet_actuate() {
        let mut drop = EwodDroplet::new([0.0, 0.0], 0.5e-3, 110.0_f64.to_radians(), 0.072);
        drop.apply_voltage(50.0, 2.2, 1e-6);
        let theta_adj = drop.contact_angle - 5.0_f64.to_radians();
        assert!(
            drop.will_actuate(theta_adj),
            "Should actuate toward lower angle"
        );
    }

    #[test]
    fn test_ewod_driving_pressure() {
        let drop = EwodDroplet::new([0.0, 0.0], 0.5e-3, 30.0_f64.to_radians(), 0.072);
        let p = drop.driving_pressure(100e-6);
        assert!(
            p > 0.0,
            "Driving pressure for acute angle should be positive"
        );
    }

    // --- MicrofluidicsResult ---

    #[test]
    fn test_microfluidics_result_convergence() {
        let mut res = MicrofluidicsResult::new();
        res.residual = 1e-8;
        assert!(res.is_converged(1e-6));
        assert!(!res.is_converged(1e-10));
    }

    // --- Capillary electrophoresis plates (free function) ---

    #[test]
    fn test_ce_plates_free_function() {
        let n = capillary_electrophoresis_plates(4e-8, 30000.0, 1e-10);
        assert!(n > 1000.0 && n < 1e8, "Plates out of range: {n}");
    }
}
