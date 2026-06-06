//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{CS2, CX, CY, W};

/// A cell with anisotropic permeability and thermal state.
#[derive(Debug, Clone)]
pub struct AnisotropicPorousCell {
    /// Volume-fraction porosity.
    pub porosity: f64,
    /// Full permeability tensor.
    pub permeability: PermeabilityTensor,
    /// Fluid temperature (K).
    pub t_fluid: f64,
    /// Solid temperature (K).
    pub t_solid: f64,
    /// Darcy velocity.
    pub velocity: [f64; 2],
    /// Local pressure.
    pub pressure: f64,
}
impl AnisotropicPorousCell {
    /// Create a new cell with isotropic permeability.
    pub fn new(porosity: f64, perm: f64) -> Self {
        Self {
            porosity,
            permeability: PermeabilityTensor::isotropic(perm),
            t_fluid: 300.0,
            t_solid: 300.0,
            velocity: [0.0; 2],
            pressure: 0.0,
        }
    }
}
/// Forchheimer inertial correction to Darcy's law (3-D).
///
/// F = -(mu/K + Cf * rho * |u| / sqrt(K)) * u
#[derive(Debug, Clone, Copy)]
pub struct ForchhheimerTerm {
    /// Permeability K (m^2).
    pub permeability: f64,
    /// Forchheimer coefficient Cf (dimensionless).
    pub cf: f64,
    /// Fluid density rho (kg/m^3).
    pub rho: f64,
}
impl ForchhheimerTerm {
    /// Create a new `ForchhheimerTerm`.
    pub fn new(permeability: f64, cf: f64, rho: f64) -> Self {
        Self {
            permeability,
            cf,
            rho,
        }
    }
    /// Compute Forchheimer force vector.
    pub fn force(&self, u: [f64; 3], viscosity: f64) -> [f64; 3] {
        let speed = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let k_sqrt = self.permeability.max(1e-300).sqrt();
        let coeff = -(viscosity / self.permeability + self.cf * self.rho * speed / k_sqrt);
        [coeff * u[0], coeff * u[1], coeff * u[2]]
    }
}
/// Extended Kozeny-Carman model with explicit tortuosity.
///
/// `K = eps^3 / (k_kc * S_v^2 * (1 - eps)^2 * T^2)`
///
/// where `S_v = 6/d_p` for spheres and `T` is the tortuosity.
#[derive(Debug, Clone)]
pub struct KozenyCarmanExtended {
    /// Particle diameter (m).
    pub particle_diameter: f64,
    /// Kozeny constant (typically 5.0 for random packing).
    pub kozeny_constant: f64,
    /// Tortuosity (dimensionless, >= 1).
    pub tortuosity: f64,
}
impl KozenyCarmanExtended {
    /// Create a new extended Kozeny-Carman model.
    pub fn new(particle_diameter: f64, kozeny_constant: f64, tortuosity: f64) -> Self {
        Self {
            particle_diameter,
            kozeny_constant,
            tortuosity,
        }
    }
    /// Compute permeability using the extended Kozeny-Carman equation.
    pub fn permeability(&self, porosity: f64) -> f64 {
        let one_minus_eps = 1.0 - porosity;
        if one_minus_eps.abs() < 1e-12 {
            return f64::MAX;
        }
        let sv = 6.0 / self.particle_diameter;
        let t2 = self.tortuosity * self.tortuosity;
        porosity.powi(3) / (self.kozeny_constant * sv * sv * one_minus_eps * one_minus_eps * t2)
    }
    /// Hydraulic diameter: `d_h = 4 * eps / (S_v * (1 - eps))`
    pub fn hydraulic_diameter(&self, porosity: f64) -> f64 {
        let sv = 6.0 / self.particle_diameter;
        let one_minus_eps = (1.0 - porosity).max(1e-12);
        4.0 * porosity / (sv * one_minus_eps)
    }
}
/// Functions for computing effective (homogenised) transport properties
/// of a fluid-saturated porous medium.
pub struct EffectiveMediumProperties;
impl EffectiveMediumProperties {
    /// Effective thermal conductivity using the series (harmonic mean) model:
    ///
    /// `1/k_eff = epsilon/k_f + (1-epsilon)/k_s`
    ///
    /// # Arguments
    /// * `k_fluid`  - fluid thermal conductivity (W/m/K)
    /// * `k_solid`  - solid thermal conductivity (W/m/K)
    /// * `porosity` - volume fraction of fluid
    pub fn effective_thermal_conductivity(k_fluid: f64, k_solid: f64, porosity: f64) -> f64 {
        let inv = porosity / k_fluid + (1.0 - porosity) / k_solid;
        1.0 / inv
    }
    /// Effective diffusivity (mass diffusion) accounting for tortuosity:
    ///
    /// `D_eff = D_free * epsilon / tortuosity`
    ///
    /// # Arguments
    /// * `d_free`    - free-fluid diffusion coefficient (m^2/s)
    /// * `porosity`  - volume fraction of fluid
    /// * `tortuosity`- geometric tortuosity factor (>= 1)
    pub fn effective_diffusivity(d_free: f64, porosity: f64, tortuosity: f64) -> f64 {
        d_free * porosity / tortuosity
    }
    /// Estimate tortuosity using the Millington-Quirk approximation:
    ///
    /// `T = porosity^(-1/3)`
    ///
    /// # Arguments
    /// * `porosity` - volume fraction of fluid
    pub fn tortuosity_estimate(porosity: f64) -> f64 {
        porosity.powf(-1.0 / 3.0)
    }
}
/// Brooks-Corey capillary pressure model for two-phase flow in porous media.
///
/// `P_c = P_d * S_e^(-1/lambda)`
///
/// where `S_e = (S_w - S_wr) / (1 - S_wr)` is the effective water saturation.
#[derive(Debug, Clone)]
pub struct BrooksCoreyCapillary {
    /// Entry (displacement) pressure P_d (Pa).
    pub entry_pressure: f64,
    /// Pore-size distribution index lambda (0.2..5.0).
    pub lambda: f64,
    /// Residual water saturation S_wr (dimensionless, 0..1).
    pub residual_saturation: f64,
}
impl BrooksCoreyCapillary {
    /// Create a new Brooks-Corey model.
    pub fn new(entry_pressure: f64, lambda: f64, residual_saturation: f64) -> Self {
        Self {
            entry_pressure,
            lambda,
            residual_saturation,
        }
    }
    /// Capillary pressure as a function of water saturation S_w.
    ///
    /// Returns `f64::MAX` if S_w ≤ S_wr (irreducible saturation).
    pub fn capillary_pressure(&self, s_w: f64) -> f64 {
        let s_e = (s_w - self.residual_saturation) / (1.0 - self.residual_saturation);
        if s_e <= 0.0 {
            return f64::MAX;
        }
        self.entry_pressure * s_e.powf(-1.0 / self.lambda)
    }
    /// Water relative permeability (Mualem-Brooks-Corey):
    ///
    /// `k_rw = S_e^((2+3*lambda)/lambda)`
    pub fn water_relative_permeability(&self, s_w: f64) -> f64 {
        let s_e =
            ((s_w - self.residual_saturation) / (1.0 - self.residual_saturation)).clamp(0.0, 1.0);
        s_e.powf((2.0 + 3.0 * self.lambda) / self.lambda)
    }
    /// Non-wetting relative permeability:
    ///
    /// `k_rnw = (1 - S_e)^2 * (1 - S_e^((2+lambda)/lambda))`
    pub fn nonwetting_relative_permeability(&self, s_w: f64) -> f64 {
        let s_e =
            ((s_w - self.residual_saturation) / (1.0 - self.residual_saturation)).clamp(0.0, 1.0);
        let one_minus_se = 1.0 - s_e;
        one_minus_se * one_minus_se * (1.0 - s_e.powf((2.0 + self.lambda) / self.lambda))
    }
}
/// Unified Darcy-Brinkman-Forchheimer (DBF) body-force model for 2D LBM.
///
/// The total drag force per unit volume is:
///
/// `F = -eps * mu / K * u  +  eps * mu_eff * lap(u)  -  eps * F_c * rho / sqrt(K) * |u| * u`
///
/// In the LBM context the Laplacian term is handled implicitly via the
/// modified relaxation time; this struct provides the explicit Darcy and
/// Forchheimer terms only.
#[derive(Debug, Clone)]
pub struct DarcyBrinkmanForchheimer {
    /// Permeability K (m^2).
    pub permeability: f64,
    /// Porosity epsilon.
    pub porosity: f64,
    /// Dynamic viscosity mu (Pa·s).
    pub viscosity: f64,
    /// Forchheimer constant C_F (dimensionless, ~ 0.55 for sphere packs).
    pub forchheimer_cf: f64,
    /// Fluid density rho (kg/m^3).
    pub density: f64,
}
impl DarcyBrinkmanForchheimer {
    /// Create a new DBF model.
    pub fn new(
        permeability: f64,
        porosity: f64,
        viscosity: f64,
        forchheimer_cf: f64,
        density: f64,
    ) -> Self {
        Self {
            permeability,
            porosity,
            viscosity,
            forchheimer_cf,
            density,
        }
    }
    /// Compute the Darcy drag force per unit volume:
    ///
    /// `F_D = -eps * mu / K * u`
    pub fn darcy_force(&self, u: [f64; 2]) -> [f64; 2] {
        let coeff = -self.porosity * self.viscosity / self.permeability;
        [coeff * u[0], coeff * u[1]]
    }
    /// Compute the Forchheimer inertial force per unit volume:
    ///
    /// `F_F = -eps * C_F * rho / sqrt(K) * |u| * u`
    pub fn forchheimer_force(&self, u: [f64; 2]) -> [f64; 2] {
        let speed = (u[0] * u[0] + u[1] * u[1]).sqrt();
        let k_sqrt = self.permeability.max(1e-300).sqrt();
        let coeff = -self.porosity * self.forchheimer_cf * self.density * speed / k_sqrt;
        [coeff * u[0], coeff * u[1]]
    }
    /// Total DBF force (Darcy + Forchheimer) per unit volume.
    pub fn total_force(&self, u: [f64; 2]) -> [f64; 2] {
        let fd = self.darcy_force(u);
        let ff = self.forchheimer_force(u);
        [fd[0] + ff[0], fd[1] + ff[1]]
    }
    /// Modified relaxation frequency accounting for Brinkman viscosity.
    ///
    /// `tau_eff = nu_eff / cs^2 + 0.5`  where `nu_eff = mu / (rho * eps)`
    pub fn brinkman_omega(&self) -> f64 {
        let nu_eff = self.viscosity / (self.density * self.porosity.max(1e-12));
        let tau_eff = nu_eff / CS2 + 0.5;
        1.0 / tau_eff
    }
}
/// Full 2x2 symmetric permeability tensor for anisotropic porous media.
///
/// Stored as `[kxx, kxy, kyy]` (symmetric: kyx = kxy).
#[derive(Debug, Clone, Copy)]
pub struct PermeabilityTensor {
    /// K_xx component (m^2).
    pub kxx: f64,
    /// K_xy = K_yx component (m^2).
    pub kxy: f64,
    /// K_yy component (m^2).
    pub kyy: f64,
}
impl PermeabilityTensor {
    /// Create an isotropic permeability tensor with K = k * I.
    pub fn isotropic(k: f64) -> Self {
        Self {
            kxx: k,
            kxy: 0.0,
            kyy: k,
        }
    }
    /// Create an anisotropic permeability tensor.
    pub fn anisotropic(kxx: f64, kxy: f64, kyy: f64) -> Self {
        Self { kxx, kxy, kyy }
    }
    /// Determinant of the 2x2 tensor: det(K) = kxx*kyy - kxy^2.
    pub fn determinant(&self) -> f64 {
        self.kxx * self.kyy - self.kxy * self.kxy
    }
    /// Compute Darcy velocity using the full tensor:
    ///
    /// `u = -(K / mu) * grad(P)`
    ///
    /// Returns `[ux, uy]`.
    pub fn darcy_velocity(&self, viscosity: f64, pressure_gradient: [f64; 2]) -> [f64; 2] {
        let inv_mu = -1.0 / viscosity;
        [
            inv_mu * (self.kxx * pressure_gradient[0] + self.kxy * pressure_gradient[1]),
            inv_mu * (self.kxy * pressure_gradient[0] + self.kyy * pressure_gradient[1]),
        ]
    }
    /// Principal permeabilities: eigenvalues of the 2x2 symmetric tensor.
    ///
    /// Returns `(k_min, k_max)`.
    pub fn principal_permeabilities(&self) -> (f64, f64) {
        let trace = self.kxx + self.kyy;
        let det = self.determinant();
        let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
        let k1 = 0.5 * (trace - disc);
        let k2 = 0.5 * (trace + disc);
        (k1, k2)
    }
    /// Anisotropy ratio: k_max / k_min.
    pub fn anisotropy_ratio(&self) -> f64 {
        let (k_min, k_max) = self.principal_permeabilities();
        if k_min <= 0.0 {
            return f64::MAX;
        }
        k_max / k_min
    }
    /// Rotate the tensor by angle `theta` (radians):
    ///
    /// K' = R^T * K * R
    pub fn rotate(&self, theta: f64) -> Self {
        let c = theta.cos();
        let s = theta.sin();
        let kxx_new = c * c * self.kxx + 2.0 * c * s * self.kxy + s * s * self.kyy;
        let kyy_new = s * s * self.kxx - 2.0 * c * s * self.kxy + c * c * self.kyy;
        let kxy_new = c * s * (self.kyy - self.kxx) + (c * c - s * s) * self.kxy;
        Self {
            kxx: kxx_new,
            kxy: kxy_new,
            kyy: kyy_new,
        }
    }
}
/// A single cell in a porous medium, carrying porosity, permeability,
/// Forchheimer coefficient, and local macroscopic state.
#[derive(Debug, Clone)]
pub struct PorousCell {
    /// Volume fraction occupied by fluid (epsilon in \[0, 1\]).
    /// 0 = fully solid, 1 = fully fluid.
    pub porosity: f64,
    /// Darcy permeability K (m^2).
    pub permeability: f64,
    /// Forchheimer inertial resistance coefficient F (m^-1).
    pub forchheimer_coeff: f64,
    /// Solid fraction = 1 - porosity.
    pub solid_fraction: f64,
    /// 2D Darcy (superficial) velocity \[ux, uy\].
    pub fluid_velocity: [f64; 2],
    /// Local fluid pressure.
    pub pressure: f64,
}
impl PorousCell {
    /// Create a new `PorousCell` with given porosity and permeability.
    pub fn new(porosity: f64, permeability: f64) -> Self {
        Self {
            porosity,
            permeability,
            forchheimer_coeff: 0.0,
            solid_fraction: 1.0 - porosity,
            fluid_velocity: [0.0; 2],
            pressure: 0.0,
        }
    }
}
/// Brinkman extension to Darcy's law, which adds an effective viscous
/// (Laplacian) term to capture viscous boundary layers near solid surfaces.
pub struct BrinkmanExtension {
    /// Ratio of effective (Brinkman) viscosity to fluid viscosity:
    /// `mu_eff / mu`.
    pub effective_viscosity_ratio: f64,
}
impl BrinkmanExtension {
    /// Create a new `BrinkmanExtension`.
    pub fn new(effective_viscosity_ratio: f64) -> Self {
        Self {
            effective_viscosity_ratio,
        }
    }
    /// Compute the Brinkman body force per unit volume:
    ///
    /// `F = -mu/K * u  +  mu_eff * laplacian(u)`
    ///
    /// The Laplacian term requires spatial stencil information not available
    /// at the cell level, so this method returns only the Darcy drag term:
    ///
    /// `F ≈ -mu/K * u`
    ///
    /// The caller is responsible for adding the Laplacian contribution when
    /// spatial derivatives are available.
    ///
    /// # Arguments
    /// * `u`           - velocity \[ux, uy\]
    /// * `permeability`- K (m^2)
    /// * `viscosity`   - dynamic viscosity mu (Pa·s)
    pub fn brinkman_force(&self, u: [f64; 2], permeability: f64, viscosity: f64) -> [f64; 2] {
        let coeff = -viscosity / permeability;
        [coeff * u[0], coeff * u[1]]
    }
}
/// Functions implementing Darcy's law and the Forchheimer extension.
pub struct DarcyFlow;
impl DarcyFlow {
    /// Compute the Darcy velocity:
    ///
    /// `u = -(K / mu) * grad(P)`
    ///
    /// # Arguments
    /// * `permeability` - K (m^2)
    /// * `viscosity`    - dynamic viscosity mu (Pa·s)
    /// * `pressure_gradient` - grad(P) \[dP/dx, dP/dy\] (Pa/m)
    pub fn darcy_velocity(
        permeability: f64,
        viscosity: f64,
        pressure_gradient: [f64; 2],
    ) -> [f64; 2] {
        let coeff = -permeability / viscosity;
        [coeff * pressure_gradient[0], coeff * pressure_gradient[1]]
    }
    /// Forchheimer correction body force for high-Reynolds porous flow.
    ///
    /// Additional inertial drag (body force per unit volume):
    ///
    /// `F_F = -F * rho * |u| * u`
    ///
    /// # Arguments
    /// * `u`             - Darcy velocity \[ux, uy\]
    /// * `permeability`  - K (m^2) — unused in pure Forchheimer term but
    ///   retained for interface symmetry
    /// * `forchheimer`   - F coefficient (m^-1)
    /// * `density`       - fluid density rho (kg/m^3)
    pub fn forchheimer_correction(
        u: [f64; 2],
        _permeability: f64,
        forchheimer: f64,
        density: f64,
    ) -> [f64; 2] {
        let speed = (u[0] * u[0] + u[1] * u[1]).sqrt();
        let coeff = -forchheimer * density * speed;
        [coeff * u[0], coeff * u[1]]
    }
    /// Darcy resistance body force per unit mass:
    ///
    /// `F_resistance = -(nu / K) * u`
    ///
    /// This is the linearised drag term used as a source in the momentum
    /// equation (units: m/s^2).
    ///
    /// # Arguments
    /// * `porosity`     - epsilon (dimensionless, for interface consistency)
    /// * `permeability` - K (m^2)
    /// * `viscosity`    - dynamic viscosity mu (Pa·s)
    /// * `u`            - velocity \[ux, uy\]
    pub fn darcy_resistance_force(
        _porosity: f64,
        permeability: f64,
        viscosity: f64,
        u: [f64; 2],
    ) -> [f64; 2] {
        let coeff = -viscosity / permeability;
        [coeff * u[0], coeff * u[1]]
    }
}
impl DarcyFlow {
    /// Compute the Darcy velocity field from a pressure field on a 2D grid.
    ///
    /// Uses central differences for the pressure gradient.
    pub fn darcy_velocity_field(
        pressure: &[f64],
        permeability: &[f64],
        viscosity: f64,
        nx: usize,
        ny: usize,
        dx: f64,
    ) -> Vec<[f64; 2]> {
        let n = nx * ny;
        let mut vel = vec![[0.0_f64; 2]; n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                let perm = permeability[k];
                let p_e = if x + 1 < nx {
                    pressure[y * nx + x + 1]
                } else {
                    pressure[k]
                };
                let p_w = if x > 0 {
                    pressure[y * nx + x - 1]
                } else {
                    pressure[k]
                };
                let p_n = if y + 1 < ny {
                    pressure[(y + 1) * nx + x]
                } else {
                    pressure[k]
                };
                let p_s = if y > 0 {
                    pressure[(y - 1) * nx + x]
                } else {
                    pressure[k]
                };
                let div_x = if x > 0 && x + 1 < nx { 2.0 * dx } else { dx };
                let div_y = if y > 0 && y + 1 < ny { 2.0 * dx } else { dx };
                let grad_px = (p_e - p_w) / div_x;
                let grad_py = (p_n - p_s) / div_y;
                vel[k] = Self::darcy_velocity(perm, viscosity, [grad_px, grad_py]);
            }
        }
        vel
    }
    /// Compute the magnitude of a 2D velocity vector.
    pub fn velocity_magnitude(u: [f64; 2]) -> f64 {
        (u[0] * u[0] + u[1] * u[1]).sqrt()
    }
    /// Reynolds number for porous flow:
    ///
    /// `Re_p = rho * |u| * d_p / mu`
    pub fn porous_reynolds_number(
        velocity_mag: f64,
        particle_diameter: f64,
        density: f64,
        viscosity: f64,
    ) -> f64 {
        density * velocity_mag * particle_diameter / viscosity
    }
}
/// Kozeny-Carman model for computing permeability of a packed-bed or
/// granular porous medium from particle diameter and porosity.
#[derive(Debug, Clone)]
pub struct KozenyCarmanModel {
    /// Representative particle (grain) diameter d_p (m).
    pub particle_diameter: f64,
}
impl KozenyCarmanModel {
    /// Create a new `KozenyCarmanModel`.
    pub fn new(particle_diameter: f64) -> Self {
        Self { particle_diameter }
    }
    /// Kozeny-Carman permeability:
    ///
    /// `K = d_p^2 * epsilon^3 / (180 * (1 - epsilon)^2)`
    ///
    /// Returns very large values as porosity approaches 1.
    pub fn permeability(&self, porosity: f64) -> f64 {
        let dp = self.particle_diameter;
        let one_minus_eps = 1.0 - porosity;
        if one_minus_eps.abs() < 1e-12 {
            return f64::MAX;
        }
        dp * dp * porosity.powi(3) / (180.0 * one_minus_eps * one_minus_eps)
    }
    /// Return the Kozeny constant (approximately 5.0 for random packings).
    pub fn kozeny_constant(&self) -> f64 {
        5.0
    }
    /// Specific surface area per unit volume:
    ///
    /// `S = 6 * (1 - epsilon) / d_p`
    pub fn specific_surface_area(&self, porosity: f64) -> f64 {
        6.0 * (1.0 - porosity) / self.particle_diameter
    }
}
/// A single D2Q9 LBM cell in a porous medium.
#[derive(Debug, Clone)]
pub struct PorousLbmCell {
    /// Distribution functions f_i for D2Q9 (9 directions).
    pub f: [f64; 9],
    /// Local porosity epsilon (0 = solid, 1 = fully fluid).
    pub porosity: f64,
    /// Dimensionless resistance coefficient (sigma) for porous drag.
    pub resistance: f64,
    /// Macroscopic density.
    pub density: f64,
    /// Macroscopic velocity \[ux, uy\].
    pub velocity: [f64; 2],
}
impl PorousLbmCell {
    /// Create a new `PorousLbmCell` at rest with equilibrium distributions.
    pub fn new(porosity: f64) -> Self {
        let mut f = [0.0_f64; 9];
        f.copy_from_slice(&W);
        Self {
            f,
            porosity,
            resistance: 0.0,
            density: 1.0,
            velocity: [0.0; 2],
        }
    }
}
/// Classification of porous-medium flow model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PorousMediumType {
    /// Pure Darcy flow (linear drag only).
    Darcy,
    /// Brinkman extension (adds viscous Laplacian term).
    Brinkman,
    /// Forchheimer–Darcy (adds inertial correction).
    ForchhheimerDarcy,
}
/// Brinkman drag body force for 3-D porous flow.
///
/// F = -mu/K * u
#[derive(Debug, Clone, Copy)]
pub struct BrinkmanForce {
    /// Permeability K (m^2).
    pub permeability: f64,
    /// Dynamic viscosity mu (Pa·s).
    pub viscosity: f64,
}
impl BrinkmanForce {
    /// Create a new `BrinkmanForce`.
    pub fn new(permeability: f64, viscosity: f64) -> Self {
        Self {
            permeability,
            viscosity,
        }
    }
    /// Compute the Brinkman body force F = -mu/K * u (N/m^3).
    pub fn body_force(&self, u: [f64; 3]) -> [f64; 3] {
        let coeff = -self.viscosity / self.permeability;
        [coeff * u[0], coeff * u[1], coeff * u[2]]
    }
}
/// Darcy resistance for 1-D pressure drop calculation.
#[derive(Debug, Clone, Copy)]
pub struct DarcyResistance {
    /// Permeability k (m^2).
    pub k: f64,
    /// Dynamic viscosity mu (Pa·s).
    pub mu: f64,
}
impl DarcyResistance {
    /// Create a new `DarcyResistance`.
    pub fn new(k: f64, mu: f64) -> Self {
        Self { k, mu }
    }
    /// Darcy law pressure drop: ΔP = (mu / K) * u * L.
    pub fn pressure_drop(&self, l: f64, u: f64) -> f64 {
        (self.mu / self.k) * u * l
    }
}
/// Driver for a volume-averaged porous-media LBM simulation.
///
/// At each step the standard BGK collision is applied using a porous-corrected
/// relaxation time, then the Darcy-Brinkman-Forchheimer body force is injected
/// via the Guo forcing scheme.
pub struct PorousMediaDriver {
    /// The underlying D2Q9 LBM porous grid.
    pub grid: PorousLbmGrid,
    /// DBF model per cell (same length as grid.cells).
    pub dbf_models: Vec<DarcyBrinkmanForchheimer>,
    /// Kinematic viscosity (lattice units).
    pub nu: f64,
    /// External body force `[fx, fy]` (lattice units/step^2).
    pub body_force: [f64; 2],
    /// Current step counter.
    pub step_count: usize,
}
impl PorousMediaDriver {
    /// Create a new porous media driver.
    ///
    /// All cells are initialised with uniform porosity and permeability.
    pub fn new(nx: usize, ny: usize, nu: f64, porosity: f64, permeability: f64) -> Self {
        let n = nx * ny;
        let mut grid = PorousLbmGrid::new(nx, ny);
        grid.set_porosity_field(|_, _| porosity);
        let dbf_models =
            vec![DarcyBrinkmanForchheimer::new(permeability, porosity, nu, 0.55, 1.0); n];
        Self {
            grid,
            dbf_models,
            nu,
            body_force: [0.0; 2],
            step_count: 0,
        }
    }
    /// Set a uniform external body force applied to all fluid cells.
    pub fn set_body_force(&mut self, fx: f64, fy: f64) {
        self.body_force = [fx, fy];
    }
    /// BGK collision with Darcy+Forchheimer implicit source.
    fn collide(&mut self) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for y in 0..ny {
            for x in 0..nx {
                let k = self.grid.idx(x, y);
                let dbf = &self.dbf_models[k];
                let omega = dbf.brinkman_omega().min(1.99);
                let rho = self.grid.cells[k].density;
                let sigma = self.grid.cells[k].resistance;
                let tau = 1.0 / omega;
                let denom = 1.0 + sigma * tau;
                let u_eff = [
                    self.grid.cells[k].velocity[0] / denom,
                    self.grid.cells[k].velocity[1] / denom,
                ];
                for i in 0..9 {
                    let feq = PorousLbmGrid::equilibrium(rho, u_eff, i);
                    self.grid.cells[k].f[i] += omega * (feq - self.grid.cells[k].f[i]);
                }
                let total_f = dbf.darcy_force(u_eff);
                let fx_total = self.body_force[0] + total_f[0];
                let fy_total = self.body_force[1] + total_f[1];
                for i in 0..9 {
                    let cu = CX[i] * u_eff[0] + CY[i] * u_eff[1];
                    let fi_term = W[i]
                        * (1.0 - 0.5 * omega)
                        * ((CX[i] / CS2 + cu * CX[i] / (CS2 * CS2)) * fx_total
                            + (CY[i] / CS2 + cu * CY[i] / (CS2 * CS2)) * fy_total);
                    self.grid.cells[k].f[i] += fi_term;
                }
            }
        }
    }
    /// Perform one complete step: collide → stream → macroscopic.
    pub fn step(&mut self) {
        self.grid.compute_macroscopic();
        self.collide();
        self.grid.stream();
        self.grid.compute_macroscopic();
        self.step_count += 1;
    }
    /// Run for `n` steps.
    pub fn run(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }
    /// Compute volume-averaged superficial velocity.
    pub fn mean_superficial_velocity(&self) -> [f64; 2] {
        let n = self.grid.nx * self.grid.ny;
        let mut sum = [0.0_f64; 2];
        for k in 0..n {
            let eps = self.grid.cells[k].porosity;
            sum[0] += self.grid.cells[k].velocity[0] * eps;
            sum[1] += self.grid.cells[k].velocity[1] * eps;
        }
        [sum[0] / n as f64, sum[1] / n as f64]
    }
    /// Estimate the Darcy pressure gradient from the x-momentum balance.
    ///
    /// At steady state: `dP/dx = -mu/K * u_D + f_x`
    /// Returns the estimated `dP/dx` using volume averages.
    pub fn estimate_pressure_gradient_x(&self) -> f64 {
        let u_sup = self.mean_superficial_velocity();
        let n = self.grid.nx * self.grid.ny;
        let mut k_mean = 0.0_f64;
        let mut mu_mean = 0.0_f64;
        for k in 0..n {
            k_mean += self.dbf_models[k].permeability;
            mu_mean += self.dbf_models[k].viscosity;
        }
        k_mean /= n as f64;
        mu_mean /= n as f64;
        -(mu_mean / k_mean) * u_sup[0] + self.body_force[0]
    }
}
/// 2D D2Q9 LBM grid for porous media simulation.
pub struct PorousLbmGrid {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Flat array of porous LBM cells (row-major: idx = y * nx + x).
    pub cells: Vec<PorousLbmCell>,
}
impl PorousLbmGrid {
    /// Create a new `PorousLbmGrid` initialised with unit porosity (fully
    /// fluid) and equilibrium distributions at rest (rho=1, u=0).
    pub fn new(nx: usize, ny: usize) -> Self {
        let cells = vec![PorousLbmCell::new(1.0); nx * ny];
        Self { nx, ny, cells }
    }
    /// Flat index for cell at grid position (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Assign porosity values from a closure `porosity_fn(x, y) -> f64`.
    pub fn set_porosity_field(&mut self, porosity_fn: impl Fn(usize, usize) -> f64) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                let idx = self.idx(x, y);
                self.cells[idx].porosity = porosity_fn(x, y);
            }
        }
    }
    /// D2Q9 Maxwell-Boltzmann equilibrium distribution for direction `i`.
    ///
    /// `f_i^eq = w_i * rho * [1 + (c_i·u)/cs^2 + (c_i·u)^2/(2 cs^4) - u^2/(2 cs^2)]`
    pub fn equilibrium(rho: f64, u: [f64; 2], i: usize) -> f64 {
        let cu = CX[i] * u[0] + CY[i] * u[1];
        let u2 = u[0] * u[0] + u[1] * u[1];
        W[i] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }
    /// BGK collision with a simple porous resistance force.
    ///
    /// For each cell the distribution is relaxed toward the local equilibrium
    /// with relaxation time `tau = 1/omega`.  An implicit linear resistance
    /// force is applied via a velocity correction before equilibrium
    /// evaluation:
    ///
    /// `u_eff = u / (1 + resistance * tau)`
    ///
    /// where `resistance` is the dimensionless drag coefficient stored in
    /// each `PorousLbmCell`.
    pub fn collide(&mut self, tau: f64) {
        let omega = 1.0 / tau;
        let n = self.nx * self.ny;
        for k in 0..n {
            let rho = self.cells[k].density;
            let sigma = self.cells[k].resistance;
            let denom = 1.0 + sigma * tau;
            let u_eff = [
                self.cells[k].velocity[0] / denom,
                self.cells[k].velocity[1] / denom,
            ];
            for i in 0..9 {
                let feq = Self::equilibrium(rho, u_eff, i);
                self.cells[k].f[i] += omega * (feq - self.cells[k].f[i]);
            }
        }
    }
    /// Standard D2Q9 streaming with periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let old: Vec<[f64; 9]> = self.cells.iter().map(|c| c.f).collect();
        for y in 0..ny {
            for x in 0..nx {
                let dst = self.idx(x, y);
                for i in 0..9 {
                    let cx_i = CX[i] as isize;
                    let cy_i = CY[i] as isize;
                    let sx = ((x as isize - cx_i).rem_euclid(nx as isize)) as usize;
                    let sy = ((y as isize - cy_i).rem_euclid(ny as isize)) as usize;
                    let src = self.idx(sx, sy);
                    self.cells[dst].f[i] = old[src][i];
                }
            }
        }
    }
    /// Compute macroscopic density and velocity from distribution functions.
    ///
    /// For porous media the velocity is porosity-corrected:
    /// `u_macro = (1 / (rho * epsilon)) * sum_i f_i * c_i`
    pub fn compute_macroscopic(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let eps = self.cells[k].porosity.max(1e-12);
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            for i in 0..9 {
                let fi = self.cells[k].f[i];
                rho += fi;
                mx += fi * CX[i];
                my += fi * CY[i];
            }
            self.cells[k].density = rho;
            self.cells[k].velocity = [mx / (rho * eps), my / (rho * eps)];
        }
    }
}
/// Representative Elementary Volume (REV) averaging for porous media.
pub struct RevAveraging;
impl RevAveraging {
    /// Compute volume-averaged porosity over a rectangular sub-region.
    ///
    /// Sub-region: `[x0, x1) x [y0, y1)` in a grid of `nx x ny` cells.
    pub fn average_porosity(
        grid: &PorousLbmGrid,
        x0: usize,
        x1: usize,
        y0: usize,
        y1: usize,
    ) -> f64 {
        let mut sum = 0.0;
        let mut count = 0u64;
        for y in y0..x1.min(grid.ny) {
            let _ = y1;
            if y >= y1 {
                break;
            }
            for x in x0..x1.min(grid.nx) {
                let idx = grid.idx(x, y);
                sum += grid.cells[idx].porosity;
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
    /// Compute volume-averaged velocity over a rectangular sub-region.
    pub fn average_velocity(
        grid: &PorousLbmGrid,
        x0: usize,
        x1: usize,
        y0: usize,
        y1: usize,
    ) -> [f64; 2] {
        let mut sum = [0.0, 0.0];
        let mut count = 0u64;
        for y in y0..y1.min(grid.ny) {
            for x in x0..x1.min(grid.nx) {
                let idx = grid.idx(x, y);
                sum[0] += grid.cells[idx].velocity[0];
                sum[1] += grid.cells[idx].velocity[1];
                count += 1;
            }
        }
        if count == 0 {
            [0.0, 0.0]
        } else {
            [sum[0] / count as f64, sum[1] / count as f64]
        }
    }
    /// Compute volume-averaged density over a rectangular sub-region.
    pub fn average_density(
        grid: &PorousLbmGrid,
        x0: usize,
        x1: usize,
        y0: usize,
        y1: usize,
    ) -> f64 {
        let mut sum = 0.0;
        let mut count = 0u64;
        for y in y0..y1.min(grid.ny) {
            for x in x0..x1.min(grid.nx) {
                let idx = grid.idx(x, y);
                sum += grid.cells[idx].density;
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
}
/// Local thermal non-equilibrium (LTNE) model for porous media heat transfer.
///
/// Solves two-temperature model where fluid and solid temperatures differ.
#[derive(Debug, Clone)]
pub struct PorousHeatTransfer {
    /// Fluid thermal conductivity (W/m/K).
    pub k_fluid: f64,
    /// Solid thermal conductivity (W/m/K).
    pub k_solid: f64,
    /// Volumetric heat transfer coefficient between phases (W/m^3/K).
    pub h_sf: f64,
    /// Fluid specific heat * density (J/m^3/K).
    pub rho_cp_fluid: f64,
    /// Solid specific heat * density (J/m^3/K).
    pub rho_cp_solid: f64,
}
impl PorousHeatTransfer {
    /// Create a new heat transfer model.
    pub fn new(
        k_fluid: f64,
        k_solid: f64,
        h_sf: f64,
        rho_cp_fluid: f64,
        rho_cp_solid: f64,
    ) -> Self {
        Self {
            k_fluid,
            k_solid,
            h_sf,
            rho_cp_fluid,
            rho_cp_solid,
        }
    }
    /// Compute the inter-phase heat transfer rate (W/m^3):
    ///
    /// `Q_sf = h_sf * (T_solid - T_fluid)`
    pub fn interphase_heat_rate(&self, t_fluid: f64, t_solid: f64) -> f64 {
        self.h_sf * (t_solid - t_fluid)
    }
    /// Fluid temperature rate of change (K/s):
    ///
    /// `dT_f/dt = (k_f * eps * lap_T_f + h_sf * (T_s - T_f)) / (eps * rho_cp_f)`
    pub fn fluid_temp_rate(
        &self,
        porosity: f64,
        lap_t_fluid: f64,
        t_fluid: f64,
        t_solid: f64,
    ) -> f64 {
        let eps = porosity.max(1e-12);
        let cond = self.k_fluid * eps * lap_t_fluid;
        let exchange = self.h_sf * (t_solid - t_fluid);
        (cond + exchange) / (eps * self.rho_cp_fluid)
    }
    /// Solid temperature rate of change (K/s):
    ///
    /// `dT_s/dt = (k_s * (1-eps) * lap_T_s - h_sf * (T_s - T_f)) / ((1-eps) * rho_cp_s)`
    pub fn solid_temp_rate(
        &self,
        porosity: f64,
        lap_t_solid: f64,
        t_fluid: f64,
        t_solid: f64,
    ) -> f64 {
        let solid_frac = (1.0 - porosity).max(1e-12);
        let cond = self.k_solid * solid_frac * lap_t_solid;
        let exchange = self.h_sf * (t_solid - t_fluid);
        (cond - exchange) / (solid_frac * self.rho_cp_solid)
    }
    /// Nusselt number for a packed bed (Wakao-Kaguei correlation):
    ///
    /// `Nu = 2 + 1.1 * Re^0.6 * Pr^(1/3)`
    pub fn nusselt_wakao_kaguei(re: f64, pr: f64) -> f64 {
        2.0 + 1.1 * re.powf(0.6) * pr.powf(1.0 / 3.0)
    }
}
