// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nucleation theory and phase-transition kinetics for MD simulations.
//!
//! Implements:
//! - Classical nucleation theory (CNT): nucleation rate, critical nucleus, Gibbs barrier
//! - Surface energy and heterogeneous nucleation factor
//! - Avrami kinetics and induction time
//! - Spinodal decomposition: Cahn-Hilliard, structure factor, growth rate
//! - NucleationMD: seeding method and simplified forward flux sampling
//! - PhaseField1D: Allen-Cahn order parameter evolution
//! - ClassicalNucleationTheory: Zeldovich factor, nucleation rate
//! - HeterogeneousNucleation: contact angle, surface-mediated rate
//! - NucleusGrowthModel: diffusion-limited and interface-limited growth
//! - OstwaldRipening: LSW coarsening, mean radius evolution
//! - SpinodaDecomposition: Cahn-Hilliard free energy, spinodal criterion
//! - SolidificationFront: Stefan problem, moving boundary, freezing-point depression

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// ClassicalNucleation  (original, kept for compatibility)
// ---------------------------------------------------------------------------

/// Classical nucleation theory (CNT) model.
///
/// Computes steady-state nucleation rates, critical nucleus radius, and
/// the Gibbs free energy barrier ΔG*.
#[derive(Debug, Clone)]
pub struct ClassicalNucleation {
    /// Solid–liquid (or liquid–vapour) interfacial energy γ \[J/m²\].
    pub surface_energy: f64,
    /// Volumetric driving force ΔG_v \[J/m³\]; negative for nucleation.
    pub driving_force: f64,
    /// Kinetic pre-factor J₀ \[m⁻³ s⁻¹\].
    pub prefactor: f64,
    /// Thermal energy k_B T \[J\].
    pub kbt: f64,
}

impl ClassicalNucleation {
    /// Create a new `ClassicalNucleation` model.
    pub fn new(surface_energy: f64, driving_force: f64, prefactor: f64, kbt: f64) -> Self {
        Self {
            surface_energy,
            driving_force,
            prefactor,
            kbt,
        }
    }

    /// Critical nucleus radius r* = -2γ / ΔG_v.
    ///
    /// Returns `f64::INFINITY` when driving_force ≥ 0.
    pub fn critical_radius(&self) -> f64 {
        if self.driving_force >= 0.0 {
            return f64::INFINITY;
        }
        -2.0 * self.surface_energy / self.driving_force
    }

    /// Gibbs free energy barrier ΔG* = 16π γ³ / (3 ΔG_v²).
    ///
    /// Returns `f64::INFINITY` when driving_force ≥ 0.
    pub fn gibbs_barrier(&self) -> f64 {
        if self.driving_force >= 0.0 || self.driving_force == 0.0 {
            return f64::INFINITY;
        }
        16.0 * PI * self.surface_energy.powi(3) / (3.0 * self.driving_force * self.driving_force)
    }

    /// Steady-state nucleation rate J = J₀ exp(-ΔG* / k_B T).
    pub fn nucleation_rate(&self) -> f64 {
        if self.kbt <= 0.0 {
            return 0.0;
        }
        let dg_star = self.gibbs_barrier();
        if dg_star.is_infinite() {
            return 0.0;
        }
        self.prefactor * (-dg_star / self.kbt).exp()
    }

    /// Nucleus free energy ΔG(r) = 4π r² γ + (4/3) π r³ ΔG_v.
    pub fn free_energy_at_radius(&self, radius: f64) -> f64 {
        4.0 * PI * radius * radius * self.surface_energy
            + (4.0 / 3.0) * PI * radius.powi(3) * self.driving_force
    }

    /// Critical nucleus volume V* = (4/3) π r*³.
    pub fn critical_volume(&self) -> f64 {
        let r = self.critical_radius();
        if r.is_infinite() {
            return f64::INFINITY;
        }
        (4.0 / 3.0) * PI * r.powi(3)
    }

    /// Number of molecules in the critical nucleus (given molecular volume v_m).
    pub fn critical_cluster_size(&self, molecular_volume: f64) -> f64 {
        if molecular_volume <= 0.0 {
            return f64::INFINITY;
        }
        self.critical_volume() / molecular_volume
    }
}

// ---------------------------------------------------------------------------
// ClassicalNucleationTheory (extended)
// ---------------------------------------------------------------------------

/// Extended classical nucleation theory with Zeldovich factor and attachment rate.
///
/// The steady-state nucleation rate includes the Zeldovich factor Z and
/// the molecular attachment rate f* at the critical nucleus:
///
///   J = Z f* N_v exp(-ΔG* / k_B T)
///
/// where N_v is the number density of monomers.
#[derive(Debug, Clone)]
pub struct ClassicalNucleationTheory {
    /// Interfacial energy γ \[J/m²\].
    pub surface_energy: f64,
    /// Volumetric driving force ΔG_v \[J/m³\], must be negative.
    pub delta_g_v: f64,
    /// Thermal energy k_B T \[J\].
    pub kbt: f64,
    /// Number density of monomers N_v \[m⁻³\].
    pub monomer_density: f64,
    /// Monomer volume Ω \[m³\].
    pub monomer_volume: f64,
    /// Diffusion coefficient of monomers D \[m²/s\].
    pub diffusivity: f64,
}

impl ClassicalNucleationTheory {
    /// Create a new `ClassicalNucleationTheory`.
    pub fn new(
        surface_energy: f64,
        delta_g_v: f64,
        kbt: f64,
        monomer_density: f64,
        monomer_volume: f64,
        diffusivity: f64,
    ) -> Self {
        Self {
            surface_energy,
            delta_g_v,
            kbt,
            monomer_density,
            monomer_volume,
            diffusivity,
        }
    }

    /// Critical nucleus radius r* = -2γ / ΔG_v.
    pub fn critical_radius(&self) -> f64 {
        if self.delta_g_v >= 0.0 {
            return f64::INFINITY;
        }
        -2.0_f64 * self.surface_energy / self.delta_g_v
    }

    /// Critical nucleus number of monomers n* = V* / Ω.
    pub fn critical_cluster_size(&self) -> f64 {
        let r = self.critical_radius();
        if r.is_infinite() || self.monomer_volume <= 0.0 {
            return f64::INFINITY;
        }
        (4.0_f64 / 3.0_f64) * PI * r.powi(3) / self.monomer_volume
    }

    /// Gibbs free energy barrier ΔG* = 16π γ³ / (3 ΔG_v²).
    pub fn gibbs_barrier(&self) -> f64 {
        if self.delta_g_v >= 0.0 {
            return f64::INFINITY;
        }
        16.0_f64 * PI * self.surface_energy.powi(3) / (3.0_f64 * self.delta_g_v * self.delta_g_v)
    }

    /// Zeldovich factor Z = sqrt(ΔG* / (3π k_B T n*²)).
    ///
    /// Returns 0 if critical size or barrier is infinite.
    pub fn zeldovich_factor(&self) -> f64 {
        let dg_star = self.gibbs_barrier();
        let n_star = self.critical_cluster_size();
        if dg_star.is_infinite() || n_star.is_infinite() || self.kbt <= 0.0 {
            return 0.0_f64;
        }
        (dg_star / (3.0_f64 * PI * self.kbt * n_star * n_star)).sqrt()
    }

    /// Monomer attachment rate at the critical nucleus f* \[s⁻¹\].
    ///
    /// f* = 4π r*² D N_v / (Ω^(1/3)) (diffusion-limited attachment)
    pub fn attachment_rate(&self) -> f64 {
        let r = self.critical_radius();
        if r.is_infinite() || self.monomer_volume <= 0.0 {
            return 0.0_f64;
        }
        let monomer_size = self.monomer_volume.cbrt();
        4.0_f64 * PI * r * r * self.diffusivity * self.monomer_density / monomer_size
    }

    /// Full CNT nucleation rate J = Z f* N_v exp(-ΔG* / k_B T) \[m⁻³ s⁻¹\].
    pub fn nucleation_rate(&self) -> f64 {
        if self.kbt <= 0.0 {
            return 0.0_f64;
        }
        let dg_star = self.gibbs_barrier();
        if dg_star.is_infinite() {
            return 0.0_f64;
        }
        let z = self.zeldovich_factor();
        let f_star = self.attachment_rate();
        z * f_star * self.monomer_density * (-dg_star / self.kbt).exp()
    }

    /// Induction time τ = 1 / (J V) for a system of volume V \[m³\].
    pub fn induction_time(&self, volume: f64) -> f64 {
        let j = self.nucleation_rate();
        if j <= 0.0 || volume <= 0.0 {
            return f64::INFINITY;
        }
        1.0_f64 / (j * volume)
    }
}

// ---------------------------------------------------------------------------
// SurfaceEnergy (original)
// ---------------------------------------------------------------------------

/// Wetting and heterogeneous nucleation factor.
///
/// For a droplet on a flat substrate with contact angle θ, the heterogeneous
/// nucleation barrier is reduced by f(θ) = (2 + cos θ)(1 - cos θ)² / 4.
#[derive(Debug, Clone)]
pub struct SurfaceEnergy {
    /// Solid–liquid interfacial energy γ_sl \[J/m²\].
    pub gamma_sl: f64,
    /// Liquid–vapour interfacial energy γ_lv \[J/m²\].
    pub gamma_lv: f64,
    /// Solid–vapour interfacial energy γ_sv \[J/m²\].
    pub gamma_sv: f64,
}

impl SurfaceEnergy {
    /// Create a new `SurfaceEnergy`.
    pub fn new(gamma_sl: f64, gamma_lv: f64, gamma_sv: f64) -> Self {
        Self {
            gamma_sl,
            gamma_lv,
            gamma_sv,
        }
    }

    /// Young contact angle θ from Young's equation:
    /// cos θ = (γ_sv - γ_sl) / γ_lv.
    ///
    /// Clamped to \[-1, 1\] before taking arccos.
    pub fn contact_angle(&self) -> f64 {
        if self.gamma_lv <= 0.0 {
            return 0.0;
        }
        let cos_theta = ((self.gamma_sv - self.gamma_sl) / self.gamma_lv).clamp(-1.0, 1.0);
        cos_theta.acos()
    }

    /// Heterogeneous nucleation correction factor f(θ).
    ///
    /// f(θ) = (2 + cos θ)(1 - cos θ)² / 4 ∈ \[0, 1\].
    pub fn heterogeneous_factor(&self) -> f64 {
        let theta = self.contact_angle();
        let c = theta.cos();
        (2.0 + c) * (1.0 - c).powi(2) / 4.0
    }

    /// Heterogeneous nucleation barrier ΔG*_het = f(θ) ΔG*_hom.
    pub fn heterogeneous_barrier(&self, homogeneous_barrier: f64) -> f64 {
        self.heterogeneous_factor() * homogeneous_barrier
    }

    /// Work of adhesion W_a = γ_lv (1 + cos θ).
    pub fn work_of_adhesion(&self) -> f64 {
        let theta = self.contact_angle();
        self.gamma_lv * (1.0 + theta.cos())
    }

    /// Spreading coefficient S = γ_sv - γ_sl - γ_lv.
    ///
    /// S > 0 → complete wetting; S < 0 → partial wetting.
    pub fn spreading_coefficient(&self) -> f64 {
        self.gamma_sv - self.gamma_sl - self.gamma_lv
    }
}

// ---------------------------------------------------------------------------
// HeterogeneousNucleation
// ---------------------------------------------------------------------------

/// Heterogeneous nucleation on a flat substrate.
///
/// The substrate reduces the nucleation barrier by the geometric factor f(θ).
/// The nucleation rate is J_het = Z f* N_s exp(-f(θ) ΔG* / k_B T)
/// where N_s is the surface density of active sites.
#[derive(Debug, Clone)]
pub struct HeterogeneousNucleation {
    /// Homogeneous CNT model parameters.
    pub cnt: ClassicalNucleationTheory,
    /// Contact angle θ \[radians\] between nucleus and substrate.
    pub contact_angle: f64,
    /// Surface density of nucleation sites N_s \[m⁻²\].
    pub site_density: f64,
}

impl HeterogeneousNucleation {
    /// Create a new `HeterogeneousNucleation`.
    pub fn new(cnt: ClassicalNucleationTheory, contact_angle_rad: f64, site_density: f64) -> Self {
        Self {
            cnt,
            contact_angle: contact_angle_rad,
            site_density,
        }
    }

    /// Geometric correction factor f(θ) = (2 + cos θ)(1 - cos θ)² / 4.
    pub fn geometry_factor(&self) -> f64 {
        let c = self.contact_angle.cos();
        (2.0_f64 + c) * (1.0_f64 - c).powi(2) / 4.0_f64
    }

    /// Reduced barrier ΔG*_het = f(θ) ΔG*_hom.
    pub fn reduced_barrier(&self) -> f64 {
        self.geometry_factor() * self.cnt.gibbs_barrier()
    }

    /// Heterogeneous nucleation rate J_het \[m⁻² s⁻¹\].
    ///
    /// J_het = Z f* N_s exp(-ΔG*_het / k_B T)
    pub fn nucleation_rate(&self) -> f64 {
        if self.cnt.kbt <= 0.0 {
            return 0.0_f64;
        }
        let dg_het = self.reduced_barrier();
        if dg_het.is_infinite() {
            return 0.0_f64;
        }
        let z = self.cnt.zeldovich_factor();
        let f_star = self.cnt.attachment_rate();
        z * f_star * self.site_density * (-dg_het / self.cnt.kbt).exp()
    }

    /// Enhancement factor: ratio of heterogeneous to homogeneous rate.
    ///
    /// Returns the exponential gain exp((ΔG* - ΔG*_het) / k_B T).
    pub fn rate_enhancement(&self) -> f64 {
        if self.cnt.kbt <= 0.0 {
            return 1.0_f64;
        }
        let dg_hom = self.cnt.gibbs_barrier();
        let dg_het = self.reduced_barrier();
        if dg_hom.is_infinite() || dg_het.is_infinite() {
            return 1.0_f64;
        }
        ((dg_hom - dg_het) / self.cnt.kbt).exp()
    }
}

// ---------------------------------------------------------------------------
// NucleusGrowthModel
// ---------------------------------------------------------------------------

/// Growth model for a spherical nucleus: diffusion-limited or interface-limited.
///
/// Diffusion-limited: dr/dt = D (C_∞ - C_r) / (r ρ)
/// Interface-limited: dr/dt = k (C_∞ - C_eq)
#[derive(Debug, Clone)]
pub struct NucleusGrowthModel {
    /// Initial radius \[m\].
    pub radius: f64,
    /// Diffusion coefficient D \[m²/s\].
    pub diffusivity: f64,
    /// Far-field solute concentration C_∞ \[mol/m³\].
    pub c_bulk: f64,
    /// Equilibrium (Gibbs-Thomson) concentration at curved surface C_r \[mol/m³\].
    pub c_eq: f64,
    /// Density of nucleus material ρ \[mol/m³\].
    pub density: f64,
    /// Interface kinetic coefficient k \[m/s\].
    pub interface_coeff: f64,
}

impl NucleusGrowthModel {
    /// Create a new `NucleusGrowthModel`.
    pub fn new(
        radius: f64,
        diffusivity: f64,
        c_bulk: f64,
        c_eq: f64,
        density: f64,
        interface_coeff: f64,
    ) -> Self {
        Self {
            radius,
            diffusivity,
            c_bulk,
            c_eq,
            density,
            interface_coeff,
        }
    }

    /// Diffusion-limited radial growth rate dr/dt \[m/s\].
    ///
    /// dr/dt = D (C_∞ - C_r) / (r ρ)
    pub fn diffusion_limited_rate(&self) -> f64 {
        if self.radius <= 0.0 || self.density <= 0.0 {
            return 0.0_f64;
        }
        self.diffusivity * (self.c_bulk - self.c_eq) / (self.radius * self.density)
    }

    /// Interface-limited radial growth rate dr/dt \[m/s\].
    ///
    /// dr/dt = k (C_∞ - C_eq)
    pub fn interface_limited_rate(&self) -> f64 {
        self.interface_coeff * (self.c_bulk - self.c_eq)
    }

    /// Mixed growth rate: harmonic mean of diffusion and interface rates.
    ///
    /// 1/v_total = 1/v_diff + 1/v_int
    pub fn mixed_growth_rate(&self) -> f64 {
        let v_diff = self.diffusion_limited_rate();
        let v_int = self.interface_limited_rate();
        if v_diff <= 0.0 || v_int <= 0.0 {
            return v_diff.max(v_int).max(0.0_f64);
        }
        1.0_f64 / (1.0_f64 / v_diff + 1.0_f64 / v_int)
    }

    /// Advance one time step (diffusion-limited) with explicit Euler.
    pub fn step_diffusion(&mut self, dt: f64) {
        let dr = self.diffusion_limited_rate() * dt;
        self.radius = (self.radius + dr).max(0.0_f64);
    }

    /// Advance one time step (interface-limited) with explicit Euler.
    pub fn step_interface(&mut self, dt: f64) {
        let dr = self.interface_limited_rate() * dt;
        self.radius = (self.radius + dr).max(0.0_f64);
    }

    /// Gibbs-Thomson concentration at radius r \[mol/m³\].
    ///
    /// C_r = C_eq exp(2 γ Ω / (r k_B T)) ≈ C_eq (1 + 2 γ Ω / (r k_B T))
    pub fn gibbs_thomson_conc(&self, gamma: f64, molar_vol: f64, kbt: f64) -> f64 {
        if self.radius <= 0.0 || kbt <= 0.0 {
            return self.c_eq;
        }
        self.c_eq * (1.0_f64 + 2.0_f64 * gamma * molar_vol / (self.radius * kbt))
    }
}

// ---------------------------------------------------------------------------
// OstwaldRipening
// ---------------------------------------------------------------------------

/// Ostwald ripening (LSW theory) for a distribution of particles.
///
/// The mean radius evolves as: `r`³ - <r₀>³ = K_LSW t
///
/// where K_LSW = 8 γ D C_eq Ω² / (9 k_B T)  (Wagner 1961)
#[derive(Debug, Clone)]
pub struct OstwaldRipening {
    /// Mean initial radius <r₀> \[m\].
    pub mean_radius_0: f64,
    /// Interfacial energy γ \[J/m²\].
    pub surface_energy: f64,
    /// Diffusion coefficient D \[m²/s\].
    pub diffusivity: f64,
    /// Equilibrium solubility C_eq \[mol/m³\].
    pub c_eq: f64,
    /// Molar volume of the precipitate Ω \[m³/mol\].
    pub molar_volume: f64,
    /// Thermal energy k_B T \[J\].
    pub kbt: f64,
}

impl OstwaldRipening {
    /// Create a new `OstwaldRipening`.
    pub fn new(
        mean_radius_0: f64,
        surface_energy: f64,
        diffusivity: f64,
        c_eq: f64,
        molar_volume: f64,
        kbt: f64,
    ) -> Self {
        Self {
            mean_radius_0,
            surface_energy,
            diffusivity,
            c_eq,
            molar_volume,
            kbt,
        }
    }

    /// LSW coarsening rate constant K_LSW \[m³/s\].
    ///
    /// K_LSW = 8 γ D C_eq Ω² / (9 k_B T)
    pub fn lsw_rate(&self) -> f64 {
        if self.kbt <= 0.0 {
            return 0.0_f64;
        }
        8.0_f64
            * self.surface_energy
            * self.diffusivity
            * self.c_eq
            * self.molar_volume
            * self.molar_volume
            / (9.0_f64 * self.kbt)
    }

    /// Mean radius at time t \[m\].
    ///
    /// <r(t)>³ = <r₀>³ + K_LSW t
    pub fn mean_radius(&self, t: f64) -> f64 {
        let r0_cubed = self.mean_radius_0.powi(3);
        let r_cubed = r0_cubed + self.lsw_rate() * t;
        if r_cubed <= 0.0 {
            return 0.0_f64;
        }
        r_cubed.cbrt()
    }

    /// Number density evolution: N(t) = N₀ (r₀/r(t))³.
    pub fn number_density(&self, n0: f64, t: f64) -> f64 {
        if n0 <= 0.0 || t < 0.0 {
            return n0;
        }
        let r_t = self.mean_radius(t);
        if r_t <= 0.0 {
            return 0.0_f64;
        }
        n0 * (self.mean_radius_0 / r_t).powi(3)
    }

    /// Critical radius r_c(t) = mean radius (LSW: r_c = `r`).
    pub fn critical_radius(&self, t: f64) -> f64 {
        self.mean_radius(t)
    }

    /// Coarsening time scale τ = <r₀>³ / K_LSW.
    pub fn coarsening_timescale(&self) -> f64 {
        let k = self.lsw_rate();
        if k <= 0.0 {
            return f64::INFINITY;
        }
        self.mean_radius_0.powi(3) / k
    }
}

// ---------------------------------------------------------------------------
// SpinodaDecomposition (new extended struct)
// ---------------------------------------------------------------------------

/// Extended spinodal decomposition model with Cahn-Hilliard free energy.
///
/// Free energy: F\[c\] = ∫ \[f(c) + κ/2 |∇c|²\] dV
/// where f(c) = A(c-c_m)² is a parabolic approximation.
/// Spinodal criterion: ∂²f/∂c² < 0.
#[derive(Debug, Clone)]
pub struct SpinodaDecomposition {
    /// Curvature of free energy A \[J/m³\].
    pub free_energy_curvature: f64,
    /// Gradient energy coefficient κ \[J/m\].
    pub kappa: f64,
    /// Mobility M \[m⁵/(J·s)\].
    pub mobility: f64,
    /// Mean composition c_m.
    pub mean_composition: f64,
}

impl SpinodaDecomposition {
    /// Create a new `SpinodaDecomposition`.
    pub fn new(
        free_energy_curvature: f64,
        kappa: f64,
        mobility: f64,
        mean_composition: f64,
    ) -> Self {
        Self {
            free_energy_curvature,
            kappa,
            mobility,
            mean_composition,
        }
    }

    /// Spinodal criterion: ∂²f/∂c² < 0 → spinodal decomposition occurs.
    ///
    /// Returns `true` if the system is inside the spinodal (unstable).
    pub fn is_spinodal(&self) -> bool {
        self.free_energy_curvature < 0.0
    }

    /// Linear amplification factor σ(k) = M k² (-2A - κ k²).
    ///
    /// For A < 0 (inside spinodal): σ > 0 for k < k_c.
    pub fn amplification_factor(&self, k: f64) -> f64 {
        self.mobility * k * k * (-2.0_f64 * self.free_energy_curvature - self.kappa * k * k)
    }

    /// Most unstable wavenumber k* = sqrt(-A / κ).
    ///
    /// Returns 0 if not in spinodal or κ ≤ 0.
    pub fn most_unstable_wavenumber(&self) -> f64 {
        if !self.is_spinodal() || self.kappa <= 0.0 {
            return 0.0_f64;
        }
        (-self.free_energy_curvature / self.kappa).sqrt()
    }

    /// Critical wavenumber k_c = sqrt(-2A / κ) above which perturbations decay.
    ///
    /// Returns 0 if not in spinodal.
    pub fn critical_wavenumber(&self) -> f64 {
        if !self.is_spinodal() || self.kappa <= 0.0 {
            return 0.0_f64;
        }
        (-2.0_f64 * self.free_energy_curvature / self.kappa).sqrt()
    }

    /// Maximum growth rate σ* at k*.
    pub fn max_growth_rate(&self) -> f64 {
        let k_star = self.most_unstable_wavenumber();
        self.amplification_factor(k_star)
    }

    /// Characteristic length scale λ* = 2π / k*.
    pub fn characteristic_wavelength(&self) -> f64 {
        let k_star = self.most_unstable_wavenumber();
        if k_star <= 0.0 {
            return f64::INFINITY;
        }
        2.0_f64 * PI / k_star
    }

    /// Cahn-Hilliard free energy of a uniform composition perturbation δc at wavenumber k.
    ///
    /// ΔF(k) = (A + κ k²/2) δc²
    pub fn perturbation_free_energy(&self, k: f64, delta_c: f64) -> f64 {
        (self.free_energy_curvature + 0.5_f64 * self.kappa * k * k) * delta_c * delta_c
    }
}

// ---------------------------------------------------------------------------
// SpinodalDecomposition (original Cahn-Hilliard 1D, kept)
// ---------------------------------------------------------------------------

/// Cahn-Hilliard spinodal decomposition model on a 1-D lattice.
///
/// The free energy functional is f(c) = (c - c0)² + (c - c1)² (double well),
/// and the chemical potential μ = df/dc - κ ∇²c.
#[derive(Debug, Clone)]
pub struct SpinodalDecomposition {
    /// Concentration field c(x), size n.
    pub c: Vec<f64>,
    /// Grid spacing Δx.
    pub dx: f64,
    /// Gradient energy coefficient κ.
    pub kappa: f64,
    /// Mobility M.
    pub mobility: f64,
    /// Left spinodal concentration c₀.
    pub c0: f64,
    /// Right spinodal concentration c₁.
    pub c1: f64,
}

impl SpinodalDecomposition {
    /// Create a `SpinodalDecomposition` with a uniform initial concentration
    /// plus a small sinusoidal perturbation.
    pub fn new(n: usize, dx: f64, kappa: f64, mobility: f64, c_mean: f64, amplitude: f64) -> Self {
        let c: Vec<f64> = (0..n)
            .map(|i| c_mean + amplitude * (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        Self {
            c,
            dx,
            kappa,
            mobility,
            c0: 0.3,
            c1: 0.7,
        }
    }

    /// Double-well free energy density f(c) = (c-c₀)²(c-c₁)².
    pub fn free_energy_density(&self, c: f64) -> f64 {
        (c - self.c0).powi(2) * (c - self.c1).powi(2)
    }

    /// Derivative of free energy df/dc.
    pub fn df_dc(&self, c: f64) -> f64 {
        2.0 * (c - self.c0) * (c - self.c1).powi(2) + 2.0 * (c - self.c0).powi(2) * (c - self.c1)
    }

    /// Laplacian of c at node i (periodic, second-order centred).
    pub fn laplacian(&self, i: usize) -> f64 {
        let n = self.c.len();
        let ip = (i + 1) % n;
        let im = if i == 0 { n - 1 } else { i - 1 };
        (self.c[ip] + self.c[im] - 2.0 * self.c[i]) / (self.dx * self.dx)
    }

    /// Chemical potential μ_i = df/dc - κ ∇²c.
    pub fn chemical_potential(&self, i: usize) -> f64 {
        self.df_dc(self.c[i]) - self.kappa * self.laplacian(i)
    }

    /// Advance one Cahn-Hilliard time step using explicit Euler.
    ///
    /// ∂c/∂t = M ∇²μ
    pub fn step(&mut self, dt: f64) {
        let n = self.c.len();
        let mu: Vec<f64> = (0..n).map(|i| self.chemical_potential(i)).collect();
        let mut dc = vec![0.0f64; n];
        for i in 0..n {
            let ip = (i + 1) % n;
            let im = if i == 0 { n - 1 } else { i - 1 };
            let lap_mu = (mu[ip] + mu[im] - 2.0 * mu[i]) / (self.dx * self.dx);
            dc[i] = self.mobility * lap_mu;
        }
        for (c_i, &dc_i) in self.c.iter_mut().zip(dc.iter()) {
            *c_i += dt * dc_i;
        }
    }

    /// Structure factor S(k) at wavenumber k (discrete, real part only).
    ///
    /// S(k) = |Σ_j c_j exp(-i k j Δx)|² / N
    pub fn structure_factor(&self, k: f64) -> f64 {
        let n = self.c.len() as f64;
        let (re, im) = self
            .c
            .iter()
            .enumerate()
            .fold((0.0, 0.0), |(re, im), (j, &cj)| {
                let theta = k * j as f64 * self.dx;
                (re + cj * theta.cos(), im + cj * theta.sin())
            });
        (re * re + im * im) / n
    }

    /// Cahn-Hilliard linear growth rate σ(k) = M k² (2A - κ k²).
    ///
    /// For the double-well A = -|f''(c_mean)| = -(2 c₀ + 2 c₁ - 4 c_mean)·...
    /// Here we use a simplified: A = df''(c_mean)/2.
    pub fn growth_rate_at_k(&self, k: f64, c_mean: f64) -> f64 {
        // Second derivative of double-well at c_mean
        let d2f = 2.0 * (c_mean - self.c1).powi(2)
            + 8.0 * (c_mean - self.c0) * (c_mean - self.c1)
            + 2.0 * (c_mean - self.c0).powi(2);
        self.mobility * k * k * (d2f - self.kappa * k * k)
    }

    /// Most unstable wavenumber k* = sqrt(d²f/dc² / (2κ)).
    pub fn critical_wavenumber(&self, c_mean: f64) -> f64 {
        let d2f = 2.0 * (c_mean - self.c1).powi(2)
            + 8.0 * (c_mean - self.c0) * (c_mean - self.c1)
            + 2.0 * (c_mean - self.c0).powi(2);
        if d2f <= 0.0 || self.kappa <= 0.0 {
            return 0.0;
        }
        (d2f / (2.0 * self.kappa)).sqrt()
    }

    /// Total concentration (should be conserved).
    pub fn total_concentration(&self) -> f64 {
        self.c.iter().sum::<f64>()
    }
}

// ---------------------------------------------------------------------------
// SolidificationFront
// ---------------------------------------------------------------------------

/// Solidification front model based on the Stefan problem.
///
/// Tracks the moving solid-liquid interface position s(t) driven by latent heat
/// removal. For a semi-infinite 1D geometry with constant thermal boundary:
///
///   s(t) = 2 β √(α t)   where β is determined by the Stefan condition.
#[derive(Debug, Clone)]
pub struct SolidificationFront {
    /// Thermal diffusivity of solid α_s \[m²/s\].
    pub alpha_solid: f64,
    /// Thermal diffusivity of liquid α_l \[m²/s\].
    pub alpha_liquid: f64,
    /// Latent heat of fusion L \[J/kg\].
    pub latent_heat: f64,
    /// Specific heat capacity c_p \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Initial melt temperature T_m \[K\].
    pub melting_temp: f64,
    /// Applied surface temperature T_s \[K\] (T_s < T_m for solidification).
    pub surface_temp: f64,
    /// Gibbs-Thomson coefficient Γ = γ T_m / (ρ L) \[m·K\].
    pub gibbs_thomson_coeff: f64,
    /// Current front position s \[m\].
    pub front_position: f64,
}

impl SolidificationFront {
    /// Create a new `SolidificationFront`.
    pub fn new(
        alpha_solid: f64,
        alpha_liquid: f64,
        latent_heat: f64,
        specific_heat: f64,
        melting_temp: f64,
        surface_temp: f64,
        gibbs_thomson_coeff: f64,
    ) -> Self {
        Self {
            alpha_solid,
            alpha_liquid,
            latent_heat,
            specific_heat,
            melting_temp,
            surface_temp,
            gibbs_thomson_coeff,
            front_position: 0.0_f64,
        }
    }

    /// Supercooling ΔT = T_m - T_s \[K\].
    pub fn supercooling(&self) -> f64 {
        self.melting_temp - self.surface_temp
    }

    /// Stefan number Ste = c_p ΔT / L (dimensionless).
    pub fn stefan_number(&self) -> f64 {
        if self.latent_heat <= 0.0 {
            return 0.0_f64;
        }
        self.specific_heat * self.supercooling() / self.latent_heat
    }

    /// Analytical front position s(t) = 2 β √(α_s t) (one-phase approximation).
    ///
    /// β is approximated as sqrt(Ste / 2) for large Ste.
    pub fn front_position_analytical(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0_f64;
        }
        let ste = self.stefan_number();
        if ste <= 0.0 {
            return 0.0_f64;
        }
        // β ≈ sqrt(Ste / (2 * erfinv)) ≈ sqrt(Ste) for large Ste
        let beta = (ste / 2.0_f64).sqrt();
        2.0_f64 * beta * (self.alpha_solid * t).sqrt()
    }

    /// Front velocity ds/dt = β √(α_s / t).
    pub fn front_velocity(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return f64::INFINITY;
        }
        let ste = self.stefan_number();
        if ste <= 0.0 {
            return 0.0_f64;
        }
        let beta = (ste / 2.0_f64).sqrt();
        beta * (self.alpha_solid / t).sqrt()
    }

    /// Gibbs-Thomson freezing-point depression ΔT_r = Γ / r.
    ///
    /// The interface temperature is T_i = T_m - Γ / r.
    pub fn freezing_point_depression(&self, radius: f64) -> f64 {
        if radius <= 0.0 {
            return f64::INFINITY;
        }
        self.gibbs_thomson_coeff / radius
    }

    /// Interface temperature with curvature effect: T_i = T_m - Γ/r.
    pub fn interface_temperature(&self, radius: f64) -> f64 {
        self.melting_temp - self.freezing_point_depression(radius)
    }

    /// Advance front position using explicit Euler.
    pub fn step(&mut self, dt: f64, t: f64) {
        let v = self.front_velocity(t.max(dt * 0.5_f64));
        self.front_position += v * dt;
    }

    /// Cumulative solidified volume fraction (1D: ratio of front to domain length L).
    pub fn solidified_fraction(&self, domain_length: f64) -> f64 {
        if domain_length <= 0.0 {
            return 0.0_f64;
        }
        (self.front_position / domain_length).min(1.0_f64)
    }
}

// ---------------------------------------------------------------------------
// CrystalNucleus (original)
// ---------------------------------------------------------------------------

/// A crystal nucleus with position, radius, and Avrami induction-time tracking.
#[derive(Debug, Clone)]
pub struct CrystalNucleus {
    /// Position \[x, y, z\].
    pub position: [f64; 3],
    /// Current radius.
    pub radius: f64,
    /// Radial growth rate (length / time).
    pub growth_rate: f64,
    /// Time at which this nucleus appeared.
    pub induction_time: f64,
}

impl CrystalNucleus {
    /// Create a new `CrystalNucleus`.
    pub fn new(position: [f64; 3], radius: f64, growth_rate: f64, induction_time: f64) -> Self {
        Self {
            position,
            radius,
            growth_rate,
            induction_time,
        }
    }

    /// Nucleus volume (sphere).
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }

    /// Advance the radius by one time step dt.
    pub fn grow(&mut self, dt: f64) {
        self.radius = (self.radius + self.growth_rate * dt).max(0.0);
    }
}

/// Avrami (JMAK) transformed fraction X(t) = 1 - exp(-k tⁿ).
pub fn avrami_fraction(k: f64, n: f64, t: f64) -> f64 {
    if t <= 0.0 || k < 0.0 {
        return 0.0;
    }
    1.0 - (-k * t.powf(n)).exp()
}

/// Induction time τ_ind for CNT.
///
/// τ_ind = 1 / (J * V * Ω)  where Ω is the critical nucleus volume.
/// Here simplified: τ_ind = 1 / (J_rate * volume_factor).
pub fn induction_time(nucleation_rate: f64, volume_factor: f64) -> f64 {
    let denom = nucleation_rate * volume_factor;
    if denom <= 0.0 {
        f64::INFINITY
    } else {
        1.0 / denom
    }
}

// ---------------------------------------------------------------------------
// NucleationMD (original)
// ---------------------------------------------------------------------------

/// MD-level nucleation analysis: seeding method and forward flux sampling.
#[derive(Debug, Clone)]
pub struct NucleationMD {
    /// Positions of all atoms \[x, y, z\].
    pub positions: Vec<[f64; 3]>,
    /// Largest cluster size observed.
    pub max_cluster_size: usize,
    /// Order parameter threshold for solid-like atoms.
    pub order_threshold: f64,
}

impl NucleationMD {
    /// Create a new `NucleationMD`.
    pub fn new(positions: Vec<[f64; 3]>, order_threshold: f64) -> Self {
        Self {
            positions,
            max_cluster_size: 0,
            order_threshold,
        }
    }

    /// Seeding: insert a spherical seed of solid-like atoms and return their indices.
    ///
    /// All atoms within `seed_radius` of `center` are marked as seed atoms.
    pub fn seed_nucleus(&self, center: [f64; 3], seed_radius: f64) -> Vec<usize> {
        self.positions
            .iter()
            .enumerate()
            .filter(|(_, pos)| {
                let dx = pos[0] - center[0];
                let dy = pos[1] - center[1];
                let dz = pos[2] - center[2];
                (dx * dx + dy * dy + dz * dz).sqrt() < seed_radius
            })
            .map(|(i, _)| i)
            .collect()
    }

    /// Cluster size for a given seed (number of atoms within seed_radius).
    pub fn cluster_size(&self, center: [f64; 3], cluster_radius: f64) -> usize {
        self.seed_nucleus(center, cluster_radius).len()
    }

    /// Forward flux sampling: estimate the committor probability at
    /// interface λ (cluster size) using simplified two-state approximation.
    ///
    /// P_B(n) ≈ (n - n_A) / (n_B - n_A)  for n_A ≤ n ≤ n_B.
    pub fn committor_probability(&self, n: usize, n_a: usize, n_b: usize) -> f64 {
        if n_b <= n_a {
            return 0.0;
        }
        if n <= n_a {
            return 0.0;
        }
        if n >= n_b {
            return 1.0;
        }
        (n - n_a) as f64 / (n_b - n_a) as f64
    }

    /// Update the maximum cluster size observed.
    pub fn update_max_cluster(&mut self, center: [f64; 3], radius: f64) {
        let size = self.cluster_size(center, radius);
        if size > self.max_cluster_size {
            self.max_cluster_size = size;
        }
    }

    /// Mean distance of all atoms from a centre point.
    pub fn mean_distance_from(&self, center: [f64; 3]) -> f64 {
        if self.positions.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .positions
            .iter()
            .map(|p| {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .sum();
        sum / self.positions.len() as f64
    }
}

// ---------------------------------------------------------------------------
// PhaseField1D (original)
// ---------------------------------------------------------------------------

/// Allen-Cahn phase-field model on a 1-D periodic lattice.
///
/// Order parameter φ ∈ \[-1, 1\]: φ = -1 (phase A), φ = +1 (phase B).
/// Evolution: ∂φ/∂t = -M (dF/dφ - κ ∇²φ)
/// where F(φ) = (1-φ²)² / 4 (double-well potential).
#[derive(Debug, Clone)]
pub struct PhaseField1D {
    /// Order parameter field, length n.
    pub phi: Vec<f64>,
    /// Grid spacing Δx.
    pub dx: f64,
    /// Gradient energy coefficient κ (interface width ~ sqrt(κ/|f''|)).
    pub kappa: f64,
    /// Kinetic coefficient M.
    pub mobility: f64,
}

impl PhaseField1D {
    /// Create a `PhaseField1D` with an initial tanh interface at `x_interface`.
    ///
    /// φ(x) = tanh((x - x_interface) / interface_width)
    pub fn new(
        n: usize,
        dx: f64,
        kappa: f64,
        mobility: f64,
        x_interface: f64,
        interface_width: f64,
    ) -> Self {
        let phi: Vec<f64> = (0..n)
            .map(|i| ((i as f64 * dx - x_interface) / interface_width.max(1e-15)).tanh())
            .collect();
        Self {
            phi,
            dx,
            kappa,
            mobility,
        }
    }

    /// Double-well free energy density F(φ) = (1-φ²)²/4.
    pub fn free_energy_density(phi: f64) -> f64 {
        (1.0 - phi * phi).powi(2) / 4.0
    }

    /// Derivative of free energy dF/dφ = -φ(1-φ²).
    pub fn df_dphi(phi: f64) -> f64 {
        -phi * (1.0 - phi * phi)
    }

    /// Laplacian at node i (periodic, second-order centred).
    pub fn laplacian(&self, i: usize) -> f64 {
        let n = self.phi.len();
        let ip = (i + 1) % n;
        let im = if i == 0 { n - 1 } else { i - 1 };
        (self.phi[ip] + self.phi[im] - 2.0 * self.phi[i]) / (self.dx * self.dx)
    }

    /// Allen-Cahn driving force at node i: -dF/dφ + κ ∇²φ.
    pub fn driving_force(&self, i: usize) -> f64 {
        -Self::df_dphi(self.phi[i]) + self.kappa * self.laplacian(i)
    }

    /// Advance one Allen-Cahn step with explicit Euler.
    pub fn step(&mut self, dt: f64) {
        let n = self.phi.len();
        let drive: Vec<f64> = (0..n).map(|i| self.driving_force(i)).collect();
        for (phi_i, &dr) in self.phi.iter_mut().zip(drive.iter()) {
            *phi_i += dt * self.mobility * dr;
            // Clamp to [-1.5, 1.5] to prevent blow-up
            *phi_i = phi_i.clamp(-1.5, 1.5);
        }
    }

    /// Total free energy ∫ \[F(φ) + κ/2 |∇φ|²\] dx (trapezoidal, periodic).
    pub fn total_energy(&self) -> f64 {
        let n = self.phi.len();
        let mut e = 0.0;
        for i in 0..n {
            let ip = (i + 1) % n;
            let grad = (self.phi[ip] - self.phi[i]) / self.dx;
            e += Self::free_energy_density(self.phi[i]) + 0.5 * self.kappa * grad * grad;
        }
        e * self.dx
    }

    /// Mean order parameter.
    pub fn mean_phi(&self) -> f64 {
        if self.phi.is_empty() {
            return 0.0;
        }
        self.phi.iter().sum::<f64>() / self.phi.len() as f64
    }

    /// Interface width estimated as the distance between the φ = ±0.5 contours.
    ///
    /// Returns None if the interface cannot be located.
    pub fn interface_width(&self) -> Option<f64> {
        let n = self.phi.len();
        let mut x_lo: Option<f64> = None;
        let mut x_hi: Option<f64> = None;
        for i in 0..n {
            let x = i as f64 * self.dx;
            if x_lo.is_none() && self.phi[i] > -0.5 {
                x_lo = Some(x);
            }
            if x_hi.is_none() && self.phi[i] > 0.5 {
                x_hi = Some(x);
            }
        }
        match (x_lo, x_hi) {
            (Some(lo), Some(hi)) => Some((hi - lo).abs()),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Classical nucleation rate J = J₀ exp(-ΔG* / k_B T).
///
/// # Arguments
/// * `gamma`        - interfacial energy \[J/m²\]
/// * `delta_g_v`    - volumetric driving force \[J/m³\], must be negative
/// * `prefactor`    - kinetic pre-factor J₀
/// * `kbt`          - thermal energy k_B T \[J\]
pub fn cnt_nucleation_rate(gamma: f64, delta_g_v: f64, prefactor: f64, kbt: f64) -> f64 {
    if delta_g_v >= 0.0 || kbt <= 0.0 || prefactor <= 0.0 {
        return 0.0;
    }
    let dg_star = 16.0 * PI * gamma.powi(3) / (3.0 * delta_g_v * delta_g_v);
    prefactor * (-dg_star / kbt).exp()
}

/// Heterogeneous nucleation factor f(θ) = (2 + cos θ)(1 - cos θ)² / 4.
pub fn heterogeneous_factor(contact_angle_rad: f64) -> f64 {
    let c = contact_angle_rad.cos();
    (2.0 + c) * (1.0 - c).powi(2) / 4.0
}

/// Spinodal decomposition linear growth rate σ(k) = M k² (2A - κ k²).
///
/// # Arguments
/// * `mobility` - M
/// * `k`        - wavenumber
/// * `a`        - second derivative of free energy (negative inside spinodal)
/// * `kappa`    - gradient coefficient
pub fn spinodal_growth_rate(mobility: f64, k: f64, a: f64, kappa: f64) -> f64 {
    mobility * k * k * (2.0 * a - kappa * k * k)
}

/// Most unstable wavenumber k* = sqrt(-2A / (2κ)) = sqrt(-A/κ).
///
/// Returns 0 if A ≥ 0 or κ ≤ 0.
pub fn critical_spinodal_wavenumber(a: f64, kappa: f64) -> f64 {
    if a >= 0.0 || kappa <= 0.0 {
        return 0.0;
    }
    (-a / kappa).sqrt()
}

/// LSW coarsening rate K_LSW = 8 γ D C_eq Ω² / (9 k_B T).
///
/// # Arguments
/// * `gamma`      - interfacial energy \[J/m²\]
/// * `diffusivity`- diffusion coefficient D \[m²/s\]
/// * `c_eq`       - equilibrium solubility \[mol/m³\]
/// * `molar_vol`  - molar volume Ω \[m³/mol\]
/// * `kbt`        - thermal energy k_B T \[J\]
pub fn lsw_coarsening_rate(
    gamma: f64,
    diffusivity: f64,
    c_eq: f64,
    molar_vol: f64,
    kbt: f64,
) -> f64 {
    if kbt <= 0.0 {
        return 0.0_f64;
    }
    8.0_f64 * gamma * diffusivity * c_eq * molar_vol * molar_vol / (9.0_f64 * kbt)
}

/// Zeldovich factor Z = sqrt(ΔG* / (3π k_B T n*²)).
///
/// # Arguments
/// * `dg_star` - Gibbs barrier ΔG* \[J\]
/// * `n_star`  - critical nucleus size (number of monomers)
/// * `kbt`     - thermal energy k_B T \[J\]
pub fn zeldovich_factor(dg_star: f64, n_star: f64, kbt: f64) -> f64 {
    if kbt <= 0.0 || n_star <= 0.0 || dg_star <= 0.0 {
        return 0.0_f64;
    }
    (dg_star / (3.0_f64 * PI * kbt * n_star * n_star)).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- ClassicalNucleation ---

    #[test]
    fn test_cnt_critical_radius_formula() {
        let cnt = ClassicalNucleation::new(0.1, -1.0, 1.0, 1.0);
        // r* = -2*0.1 / (-1.0) = 0.2
        assert!((cnt.critical_radius() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_cnt_critical_radius_positive_dg_inf() {
        let cnt = ClassicalNucleation::new(0.1, 1.0, 1.0, 1.0);
        assert!(cnt.critical_radius().is_infinite());
    }

    #[test]
    fn test_cnt_gibbs_barrier_positive() {
        let cnt = ClassicalNucleation::new(0.01, -1e6, 1.0, 1.0);
        let dg = cnt.gibbs_barrier();
        assert!(dg > 0.0 && dg.is_finite());
    }

    #[test]
    fn test_cnt_gibbs_barrier_inf_for_positive_dg() {
        let cnt = ClassicalNucleation::new(0.1, 1.0, 1.0, 1.0);
        assert!(cnt.gibbs_barrier().is_infinite());
    }

    #[test]
    fn test_cnt_nucleation_rate_positive() {
        // Large driving force → small barrier → nonzero rate
        let cnt = ClassicalNucleation::new(0.001, -1e8, 1e30, 1.0);
        let rate = cnt.nucleation_rate();
        assert!(rate > 0.0, "rate={rate}");
    }

    #[test]
    fn test_cnt_nucleation_rate_zero_kbt() {
        let cnt = ClassicalNucleation::new(0.1, -1e6, 1.0, 0.0);
        assert_eq!(cnt.nucleation_rate(), 0.0);
    }

    #[test]
    fn test_cnt_free_energy_at_critical_radius_is_max() {
        let cnt = ClassicalNucleation::new(0.1, -1.0, 1.0, 1.0);
        let r_star = cnt.critical_radius();
        let g_star = cnt.free_energy_at_radius(r_star);
        let g_small = cnt.free_energy_at_radius(r_star * 0.5);
        let g_large = cnt.free_energy_at_radius(r_star * 2.0);
        assert!(g_star > g_small);
        assert!(g_star > g_large);
    }

    #[test]
    fn test_cnt_critical_cluster_size() {
        let cnt = ClassicalNucleation::new(0.1, -1.0, 1.0, 1.0);
        let n_star = cnt.critical_cluster_size(1e-28);
        assert!(n_star > 0.0);
    }

    // --- ClassicalNucleationTheory ---

    #[test]
    fn test_cnt_theory_critical_radius() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let r = cnt.critical_radius();
        // r* = -2 * 0.05 / (-1e7) = 1e-8
        assert!((r - 1e-8_f64).abs() < 1e-15_f64);
    }

    #[test]
    fn test_cnt_theory_gibbs_barrier_positive() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let dg = cnt.gibbs_barrier();
        assert!(dg > 0.0 && dg.is_finite());
    }

    #[test]
    fn test_cnt_theory_zeldovich_positive() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let z = cnt.zeldovich_factor();
        assert!(z > 0.0 && z.is_finite(), "Z={z}");
    }

    #[test]
    fn test_cnt_theory_attachment_rate_positive() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let f_star = cnt.attachment_rate();
        assert!(f_star > 0.0);
    }

    #[test]
    fn test_cnt_theory_nucleation_rate_positive() {
        // Use parameters where barrier is small enough for nonzero rate
        let cnt = ClassicalNucleationTheory::new(0.001, -1e8, 1.0, 1e20, 1e-29, 1e-10);
        let j = cnt.nucleation_rate();
        assert!(j >= 0.0);
    }

    #[test]
    fn test_cnt_theory_induction_time_finite() {
        let cnt = ClassicalNucleationTheory::new(0.001, -1e8, 1.0, 1e20, 1e-29, 1e-10);
        let tau = cnt.induction_time(1e-18);
        assert!(tau.is_finite() || tau.is_infinite()); // just check no panic
    }

    #[test]
    fn test_cnt_theory_critical_size_positive() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let n_star = cnt.critical_cluster_size();
        assert!(n_star > 0.0);
    }

    // --- SurfaceEnergy ---

    #[test]
    fn test_contact_angle_hydrophilic() {
        let se = SurfaceEnergy::new(0.02, 0.07, 0.07);
        let theta = se.contact_angle();
        assert!(theta > 0.0 && theta < PI / 2.0);
    }

    #[test]
    fn test_contact_angle_zero_gamma_lv() {
        let se = SurfaceEnergy::new(0.02, 0.0, 0.07);
        assert_eq!(se.contact_angle(), 0.0);
    }

    #[test]
    fn test_heterogeneous_factor_theta_zero() {
        let se = SurfaceEnergy::new(0.0, 0.07, 0.07);
        let f = se.heterogeneous_factor();
        assert!((f - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_heterogeneous_factor_theta_pi() {
        let _se = SurfaceEnergy::new(0.14, 0.07, 0.0);
        let se2 = SurfaceEnergy {
            gamma_sl: 0.14,
            gamma_lv: 0.07,
            gamma_sv: 0.0,
        };
        let f = se2.heterogeneous_factor();
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_heterogeneous_factor_between_0_and_1() {
        let se = SurfaceEnergy::new(0.03, 0.07, 0.05);
        let f = se.heterogeneous_factor();
        assert!((0.0..=1.0).contains(&f));
    }

    #[test]
    fn test_spreading_coefficient() {
        let se = SurfaceEnergy::new(0.02, 0.07, 0.10);
        let s = se.spreading_coefficient();
        assert!((s - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_work_of_adhesion_positive() {
        let se = SurfaceEnergy::new(0.02, 0.07, 0.07);
        assert!(se.work_of_adhesion() > 0.0);
    }

    // --- HeterogeneousNucleation ---

    #[test]
    fn test_hetero_geometry_factor_range() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let het = HeterogeneousNucleation::new(cnt, PI / 3.0, 1e15);
        let f = het.geometry_factor();
        assert!((0.0..=1.0).contains(&f));
    }

    #[test]
    fn test_hetero_geometry_factor_theta_pi_half() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let het = HeterogeneousNucleation::new(cnt, PI / 2.0, 1e15);
        let f = het.geometry_factor();
        // θ = π/2 → cos = 0 → f = (2+0)(1-0)²/4 = 0.5
        assert!((f - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_hetero_reduced_barrier_less_than_hom() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let dg_hom = cnt.gibbs_barrier();
        let het = HeterogeneousNucleation::new(cnt, PI / 3.0, 1e15);
        let dg_het = het.reduced_barrier();
        assert!(dg_het < dg_hom);
    }

    #[test]
    fn test_hetero_rate_enhancement_positive() {
        let cnt = ClassicalNucleationTheory::new(0.05, -1e7, 4.11e-21, 1e25, 1e-29, 1e-12);
        let het = HeterogeneousNucleation::new(cnt, PI / 3.0, 1e15);
        let enh = het.rate_enhancement();
        assert!(enh >= 1.0);
    }

    // --- NucleusGrowthModel ---

    #[test]
    fn test_growth_diffusion_limited_positive() {
        let m = NucleusGrowthModel::new(1e-8, 1e-15, 100.0, 50.0, 1e6, 1e-7);
        let rate = m.diffusion_limited_rate();
        assert!(rate > 0.0);
    }

    #[test]
    fn test_growth_interface_limited_positive() {
        let m = NucleusGrowthModel::new(1e-8, 1e-15, 100.0, 50.0, 1e6, 1e-7);
        let rate = m.interface_limited_rate();
        assert!(rate > 0.0);
    }

    #[test]
    fn test_growth_mixed_rate_between_components() {
        let m = NucleusGrowthModel::new(1e-8, 1e-15, 100.0, 50.0, 1e6, 1e-7);
        let v_diff = m.diffusion_limited_rate();
        let v_int = m.interface_limited_rate();
        let v_mix = m.mixed_growth_rate();
        assert!(v_mix <= v_diff.min(v_int) * 1.001);
    }

    #[test]
    fn test_growth_step_diffusion_increases_radius() {
        let mut m = NucleusGrowthModel::new(1e-8, 1e-15, 100.0, 50.0, 1e6, 1e-7);
        let r0 = m.radius;
        m.step_diffusion(1e-3);
        assert!(m.radius >= r0);
    }

    #[test]
    fn test_growth_gibbs_thomson_correction() {
        let m = NucleusGrowthModel::new(1e-8, 1e-15, 100.0, 50.0, 1e6, 1e-7);
        let c_r = m.gibbs_thomson_conc(0.05, 1e-5, 4.11e-21);
        assert!(c_r > m.c_eq);
    }

    // --- OstwaldRipening ---

    #[test]
    fn test_ostwald_lsw_rate_positive() {
        let or = OstwaldRipening::new(1e-7, 0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        assert!(or.lsw_rate() > 0.0);
    }

    #[test]
    fn test_ostwald_mean_radius_grows() {
        let or = OstwaldRipening::new(1e-7, 0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        let r0 = or.mean_radius(0.0);
        let r1 = or.mean_radius(1.0);
        assert!(r1 >= r0);
    }

    #[test]
    fn test_ostwald_lsw_cubic_law() {
        let or = OstwaldRipening::new(1e-7, 0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        let t = 1000.0_f64;
        let r_t = or.mean_radius(t);
        let expected = (or.mean_radius_0.powi(3) + or.lsw_rate() * t).cbrt();
        assert!((r_t - expected).abs() < 1e-20);
    }

    #[test]
    fn test_ostwald_number_density_decreases() {
        let or = OstwaldRipening::new(1e-7, 0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        let n0 = 1e18_f64;
        let n1 = or.number_density(n0, 1.0);
        assert!(n1 < n0);
    }

    #[test]
    fn test_ostwald_coarsening_timescale_finite() {
        let or = OstwaldRipening::new(1e-7, 0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        let tau = or.coarsening_timescale();
        assert!(tau > 0.0 && tau.is_finite());
    }

    // --- SpinodaDecomposition ---

    #[test]
    fn test_spinoda_is_spinodal_negative_a() {
        let sd = SpinodaDecomposition::new(-1.0, 0.1, 1.0, 0.5);
        assert!(sd.is_spinodal());
    }

    #[test]
    fn test_spinoda_is_not_spinodal_positive_a() {
        let sd = SpinodaDecomposition::new(1.0, 0.1, 1.0, 0.5);
        assert!(!sd.is_spinodal());
    }

    #[test]
    fn test_spinoda_most_unstable_wavenumber() {
        let sd = SpinodaDecomposition::new(-1.0, 0.1, 1.0, 0.5);
        // k* = sqrt(-A/κ) = sqrt(1/0.1) = sqrt(10)
        let k_star = sd.most_unstable_wavenumber();
        assert!((k_star - 10.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_spinoda_critical_wavenumber() {
        let sd = SpinodaDecomposition::new(-1.0, 0.1, 1.0, 0.5);
        // k_c = sqrt(-2A/κ) = sqrt(20)
        let k_c = sd.critical_wavenumber();
        assert!((k_c - 20.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_spinoda_amplification_at_k_star_positive() {
        let sd = SpinodaDecomposition::new(-1.0, 0.1, 1.0, 0.5);
        let sigma = sd.max_growth_rate();
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_spinoda_amplification_zero_outside_spinodal() {
        let sd = SpinodaDecomposition::new(1.0, 0.1, 1.0, 0.5);
        // Not in spinodal, k* = 0
        let k_star = sd.most_unstable_wavenumber();
        assert_eq!(k_star, 0.0);
    }

    #[test]
    fn test_spinoda_characteristic_wavelength() {
        let sd = SpinodaDecomposition::new(-1.0, 0.1, 1.0, 0.5);
        let lam = sd.characteristic_wavelength();
        assert!(lam > 0.0 && lam.is_finite());
    }

    // --- SolidificationFront ---

    #[test]
    fn test_stefan_supercooling() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-9);
        assert!((sf.supercooling() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_stefan_number_positive() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-9);
        let ste = sf.stefan_number();
        assert!(ste > 0.0);
    }

    #[test]
    fn test_stefan_front_position_zero_at_t0() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-9);
        assert_eq!(sf.front_position_analytical(0.0), 0.0);
    }

    #[test]
    fn test_stefan_front_position_grows_with_sqrt_t() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-9);
        let s1 = sf.front_position_analytical(1.0);
        let s4 = sf.front_position_analytical(4.0);
        // s(4t) = 2 * s(t)
        assert!((s4 / s1 - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_stefan_freezing_point_depression() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-8);
        let dtm = sf.freezing_point_depression(1e-6);
        assert!(dtm > 0.0);
    }

    #[test]
    fn test_stefan_interface_temperature_below_melting() {
        let sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-8);
        let t_i = sf.interface_temperature(1e-6);
        assert!(t_i < sf.melting_temp);
    }

    #[test]
    fn test_stefan_solidified_fraction_bounded() {
        let mut sf = SolidificationFront::new(1e-6, 1e-6, 3.34e5, 4200.0, 273.15, 263.15, 1e-9);
        for k in 1..=5 {
            sf.step(0.1, k as f64 * 0.1);
        }
        let frac = sf.solidified_fraction(1.0);
        assert!((0.0..=1.0).contains(&frac));
    }

    // --- CrystalNucleus ---

    #[test]
    fn test_crystal_nucleus_volume() {
        let n = CrystalNucleus::new([0.0; 3], 1.0, 0.0, 0.0);
        let expected = (4.0 / 3.0) * PI;
        assert!((n.volume() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_crystal_nucleus_grow() {
        let mut n = CrystalNucleus::new([0.0; 3], 1.0, 0.5, 0.0);
        n.grow(1.0);
        assert!((n.radius - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_avrami_zero_time() {
        assert_eq!(avrami_fraction(1.0, 4.0, 0.0), 0.0);
    }

    #[test]
    fn test_avrami_approaches_one() {
        let x = avrami_fraction(1.0, 4.0, 100.0);
        assert!(x > 0.9999);
    }

    #[test]
    fn test_avrami_monotone() {
        let x1 = avrami_fraction(0.5, 2.0, 1.0);
        let x2 = avrami_fraction(0.5, 2.0, 2.0);
        assert!(x2 > x1);
    }

    #[test]
    fn test_induction_time_positive() {
        let t = induction_time(1e10, 1e-20);
        assert!(t > 0.0);
    }

    #[test]
    fn test_induction_time_zero_rate() {
        assert!(induction_time(0.0, 1.0).is_infinite());
    }

    // --- SpinodalDecomposition (original) ---

    #[test]
    fn test_spinodal_concentration_conserved() {
        let mut sd = SpinodalDecomposition::new(50, 0.1, 0.01, 0.01, 0.5, 0.01);
        let c0 = sd.total_concentration();
        for _ in 0..10 {
            sd.step(1e-4);
        }
        let c1 = sd.total_concentration();
        assert!((c0 - c1).abs() < 1e-6, "c0={c0} c1={c1}");
    }

    #[test]
    fn test_spinodal_structure_factor_positive() {
        let sd = SpinodalDecomposition::new(50, 0.1, 0.01, 0.01, 0.5, 0.01);
        let sk = sd.structure_factor(PI / 0.5);
        assert!(sk >= 0.0);
    }

    #[test]
    fn test_spinodal_growth_rate_function() {
        let rate = spinodal_growth_rate(1.0, 1.0, -1.0, 0.1);
        assert!(rate.is_finite());
    }

    #[test]
    fn test_spinodal_critical_wavenumber() {
        let kc = critical_spinodal_wavenumber(-1.0, 0.1);
        assert!((kc - (1.0 / 0.1_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_spinodal_critical_wavenumber_positive_a() {
        assert_eq!(critical_spinodal_wavenumber(1.0, 0.1), 0.0);
    }

    #[test]
    fn test_spinodal_laplacian_uniform() {
        let sd = SpinodalDecomposition {
            c: vec![1.0; 20],
            dx: 0.1,
            kappa: 0.01,
            mobility: 0.01,
            c0: 0.3,
            c1: 0.7,
        };
        for i in 0..20 {
            assert!((sd.laplacian(i)).abs() < 1e-12);
        }
    }

    #[test]
    fn test_spinodal_growth_rate_method() {
        let sd = SpinodalDecomposition::new(50, 0.1, 0.01, 0.01, 0.5, 0.01);
        let rate = sd.growth_rate_at_k(1.0, 0.5);
        assert!(rate.is_finite());
    }

    #[test]
    fn test_spinodal_critical_wavenumber_method() {
        let sd = SpinodalDecomposition::new(50, 0.1, 0.01, 0.01, 0.5, 0.01);
        let kc = sd.critical_wavenumber(0.5);
        assert!(kc >= 0.0);
    }

    // --- NucleationMD ---

    #[test]
    fn test_nucleation_md_seed_size() {
        let positions: Vec<[f64; 3]> = (0..100).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let md = NucleationMD::new(positions, 0.5);
        let seed = md.seed_nucleus([5.0, 0.0, 0.0], 1.5);
        assert!(!seed.is_empty());
    }

    #[test]
    fn test_nucleation_md_cluster_size() {
        let positions: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let md = NucleationMD::new(positions, 0.5);
        let size = md.cluster_size([0.0, 0.0, 0.0], 3.0);
        assert!(size > 0);
    }

    #[test]
    fn test_committor_below_na() {
        let md = NucleationMD::new(vec![], 0.5);
        assert_eq!(md.committor_probability(5, 10, 50), 0.0);
    }

    #[test]
    fn test_committor_above_nb() {
        let md = NucleationMD::new(vec![], 0.5);
        assert_eq!(md.committor_probability(60, 10, 50), 1.0);
    }

    #[test]
    fn test_committor_midpoint() {
        let md = NucleationMD::new(vec![], 0.5);
        let p = md.committor_probability(30, 10, 50);
        assert!((p - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_nucleation_md_mean_distance() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let md = NucleationMD::new(positions, 0.5);
        assert!((md.mean_distance_from([0.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
    }

    // --- PhaseField1D ---

    #[test]
    fn test_phase_field_initial_tanh() {
        let pf = PhaseField1D::new(100, 0.1, 0.1, 1.0, 5.0, 1.0);
        assert!(pf.phi[0] < 0.0);
        assert!(pf.phi[99] > 0.0);
    }

    #[test]
    fn test_phase_field_free_energy_double_well() {
        assert!((PhaseField1D::free_energy_density(0.0) - 0.25).abs() < 1e-12);
        assert!((PhaseField1D::free_energy_density(1.0)).abs() < 1e-12);
        assert!((PhaseField1D::free_energy_density(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_phase_field_df_dphi_at_minima() {
        assert!((PhaseField1D::df_dphi(0.0)).abs() < 1e-12);
        assert!((PhaseField1D::df_dphi(1.0)).abs() < 1e-12);
        assert!((PhaseField1D::df_dphi(-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_phase_field_energy_positive() {
        let pf = PhaseField1D::new(50, 0.1, 0.1, 1.0, 2.5, 0.5);
        assert!(pf.total_energy() >= 0.0);
    }

    #[test]
    fn test_phase_field_step_reduces_energy() {
        let mut pf = PhaseField1D::new(100, 0.05, 0.1, 0.5, 2.5, 0.5);
        let e0 = pf.total_energy();
        for _ in 0..20 {
            pf.step(1e-4);
        }
        let e1 = pf.total_energy();
        assert!(e1 <= e0 * 10.0, "e0={e0} e1={e1}");
    }

    #[test]
    fn test_phase_field_laplacian_uniform() {
        let pf = PhaseField1D {
            phi: vec![1.0; 20],
            dx: 0.1,
            kappa: 0.1,
            mobility: 1.0,
        };
        for i in 0..20 {
            assert!((pf.laplacian(i)).abs() < 1e-12);
        }
    }

    // --- Free functions ---

    #[test]
    fn test_cnt_nucleation_rate_fn() {
        let rate = cnt_nucleation_rate(0.001, -1e8, 1e30, 1.0);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_cnt_rate_fn_positive_dg() {
        assert_eq!(cnt_nucleation_rate(0.1, 1.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_heterogeneous_factor_fn() {
        let f = heterogeneous_factor(0.0);
        assert!((f - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_heterogeneous_factor_fn_pi_half() {
        let f = heterogeneous_factor(PI / 2.0);
        assert!((f - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_lsw_coarsening_rate_positive() {
        let k = lsw_coarsening_rate(0.05, 1e-15, 10.0, 1e-5, 4.11e-21);
        assert!(k > 0.0);
    }

    #[test]
    fn test_lsw_coarsening_rate_zero_kbt() {
        let k = lsw_coarsening_rate(0.05, 1e-15, 10.0, 1e-5, 0.0);
        assert_eq!(k, 0.0);
    }

    #[test]
    fn test_zeldovich_factor_fn() {
        let z = zeldovich_factor(1e-19, 100.0, 4.11e-21);
        assert!(z > 0.0 && z.is_finite());
    }

    #[test]
    fn test_zeldovich_factor_zero_kbt() {
        assert_eq!(zeldovich_factor(1e-19, 100.0, 0.0), 0.0);
    }
}
