//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

use crate::kernel::{CubicSplineKernel, SphKernel};
use crate::neighbor::SpatialHash;
use oxiphysics_core::math::Vec3;

/// Relative threshold on the colour-gradient magnitude `|∇c|` (in units of
/// `1/h`) below which a particle is regarded as lying away from every interface
/// and gets the zero unit normal.
const NORMAL_THRESHOLD: f64 = 0.01;

/// Reference density (kg/m³) assumed for particles of a phase that has not been registered.
const DEFAULT_DENSITY_REF: f64 = 1000.0;

/// Neighbour lists (self excluded) of all particles within `radius`.
fn particle_neighbours(particles: &[MultiphaseParticle], radius: f64) -> Vec<Vec<usize>> {
    let positions: Vec<Vec3> = particles
        .iter()
        .map(|p| Vec3::new(p.position[0], p.position[1], p.position[2]))
        .collect();
    SpatialHash::find_all_neighbors(&positions, radius)
}

/// Sum `Σ_{a,b} W(s √(a² + b²), h)` of the 3-D cubic-spline kernel over the
/// points of an infinite planar square lattice of spacing `spacing`
/// (including the origin).
///
/// A planar lattice of identical particles of mass `m` has the kernel-summed
/// density `m · lattice_kernel_sum`, so `m = ρ / lattice_kernel_sum`
/// reproduces the density `ρ` at interior sites.
pub(super) fn planar_lattice_kernel_sum(spacing: f64, h: f64) -> f64 {
    let reach = ((2.0 * h / spacing).ceil() as i64).clamp(1, 1024);
    let mut sum = 0.0_f64;
    for a in -reach..=reach {
        for b in -reach..=reach {
            let r = spacing * ((a * a + b * b) as f64).sqrt();
            sum += cubic_spline_w(r, h);
        }
    }
    sum
}

/// Density-ratio multiphase SPH solver (Hu & Adams 2006 variant).
///
/// Couples distinct phases (liquid, gas, …) via a modified symmetric pressure
/// force that accounts for the density discontinuity at the interface.  The
/// pressure term used is:
///
/// ```text
/// p_i / ρ_i²  +  p_j / ρ_j²
/// ```
///
/// which is the standard symmetric SPH pressure discretisation and remains
/// well-behaved across density jumps.
pub struct DensityRatioMultiphase {
    /// Per-particle phase data.
    pub phases: Vec<PhaseParticle>,
    /// Interface-smoothing parameter η (used to regularise normal vectors).
    pub interface_smoothing: f64,
}
impl DensityRatioMultiphase {
    /// Create a new empty multiphase system.
    pub fn new() -> Self {
        Self {
            phases: Vec::new(),
            interface_smoothing: 1e-4,
        }
    }
    /// Add a particle with the given phase properties.
    ///
    /// Returns the index of the newly added particle.
    pub fn add_particle(&mut self, phase: Phase, density_ref: f64, viscosity: f64) -> usize {
        let idx = self.phases.len();
        self.phases.push(PhaseParticle {
            phase,
            density_ref,
            viscosity,
            surface_tension: 0.0,
        });
        idx
    }
    /// Colour function for phase α at particle *i*.
    ///
    /// Computes  `C_α(i) = Σ_{j∈α} V_j · W(r_ij, h)`  where the sum runs
    /// over all particles that belong to the same phase as particle *i*.
    ///
    /// The result is a number in `[0, 1]` that measures the local volume
    /// fraction occupied by phase α.
    pub fn color_function(&self, i: usize, positions: &[[f64; 3]], volumes: &[f64], h: f64) -> f64 {
        let kernel = CubicSplineKernel;
        let phase_i = self.phases[i].phase;
        let pi = positions[i];
        let mut c = 0.0_f64;
        for j in 0..positions.len() {
            if self.phases[j].phase != phase_i {
                continue;
            }
            let dx = pi[0] - positions[j][0];
            let dy = pi[1] - positions[j][1];
            let dz = pi[2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            c += volumes[j] * kernel.w(r, h);
        }
        c
    }
    /// Interface normal at particle *i* (gradient of the colour function).
    ///
    /// Points from phase α into the complementary phase and has magnitude
    /// proportional to the local interface curvature.
    pub fn interface_normal(
        &self,
        i: usize,
        positions: &[[f64; 3]],
        volumes: &[f64],
        h: f64,
    ) -> [f64; 3] {
        let kernel = CubicSplineKernel;
        let phase_i = self.phases[i].phase;
        let pi = positions[i];
        let mut n = [0.0_f64; 3];
        for j in 0..positions.len() {
            if self.phases[j].phase != phase_i {
                continue;
            }
            let dx = pi[0] - positions[j][0];
            let dy = pi[1] - positions[j][1];
            let dz = pi[2] - positions[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-14 {
                continue;
            }
            let dw_dr = kernel.grad_w(r, h);
            let scale = volumes[j] * dw_dr / r;
            n[0] += scale * dx;
            n[1] += scale * dy;
            n[2] += scale * dz;
        }
        n
    }
    /// Symmetric SPH pressure correction factor across the interface.
    ///
    /// Returns the value `p_i / ρ_i² + p_j / ρ_j²`, which is the symmetric
    /// pressure term used in the standard SPH discretisation.  Multiplying
    /// this by `−m_j · ∇W` gives the pressure force contribution on particle
    /// *i* from neighbour *j*.
    pub fn pressure_correction(rho_i: f64, rho_j: f64, p_i: f64, p_j: f64) -> f64 {
        p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)
    }
    /// Viscous force on particle *i* due to neighbour *j* with phase-dependent
    /// viscosity.
    ///
    /// Uses the harmonic-mean viscosity `μ_ij = 2 μ_i μ_j / (μ_i + μ_j)` to
    /// give a smooth transition across the interface.  The force formula is:
    ///
    /// ```text
    /// F_visc = μ_ij · (m_i + m_j) · (v_ij · r_ij)
    ///          / (ρ_i · ρ_j · (r² + 0.01 h²)) · ∇W(r, h) · r̂
    /// ```
    pub fn viscous_force(
        &self,
        i: usize,
        j: usize,
        r_ij: [f64; 3],
        v_ij: [f64; 3],
        rho_i: f64,
        rho_j: f64,
        h: f64,
    ) -> [f64; 3] {
        let kernel = CubicSplineKernel;
        let mu_i = self.phases[i].viscosity;
        let mu_j = self.phases[j].viscosity;
        let mu_ij = if mu_i + mu_j > 1e-30 {
            2.0 * mu_i * mu_j / (mu_i + mu_j)
        } else {
            0.0
        };
        let r2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
        let r = r2.sqrt();
        if r < 1e-14 {
            return [0.0; 3];
        }
        let v_dot_r = v_ij[0] * r_ij[0] + v_ij[1] * r_ij[1] + v_ij[2] * r_ij[2];
        let m_sum = 1.0 + 1.0;
        let dw_dr = kernel.grad_w(r, h);
        let denom = rho_i * rho_j * (r2 + 0.01 * h * h);
        let scale = mu_ij * m_sum * v_dot_r / denom * dw_dr / r;
        [scale * r_ij[0], scale * r_ij[1], scale * r_ij[2]]
    }
}
/// Interface sharpening (re-initialization) for the SPH color function.
///
/// Applies a single step of the Allen–Cahn–type interface compression:
/// ```text
/// ∂c/∂τ = ε ∇·( c(1−c) n̂ )  −  ε ∇²c
/// ```
/// where `τ` is a pseudo-time, `ε` controls the interface width, and `n̂` is
/// the interface normal (unit gradient of the color field).
///
/// In practice we advance one pseudo-timestep `dtau` using a finite SPH
/// approximation, returning the new color values.
pub struct InterfaceSharpening {
    /// Interface width parameter ε (m).
    pub epsilon: f64,
    /// Smoothing length h (m).
    pub h: f64,
}
impl InterfaceSharpening {
    /// Create a new interface sharpening operator.
    pub fn new(epsilon: f64, h: f64) -> Self {
        Self { epsilon, h }
    }
    /// Perform one re-initialization pseudo-step.
    ///
    /// `color[i]` should be in `[0, 1]` before calling.  Returns updated color values.
    pub fn step(
        &self,
        positions: &[[f64; 3]],
        volumes: &[f64],
        color: &[f64],
        dtau: f64,
    ) -> Vec<f64> {
        let n = positions.len();
        let h = self.h;
        let mut grads = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for j in 0..n {
                let dx = [
                    positions[i][0] - positions[j][0],
                    positions[i][1] - positions[j][1],
                    positions[i][2] - positions[j][2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let r = r2.sqrt();
                if r < 1e-14 {
                    continue;
                }
                let dw = cubic_spline_dw(r, h);
                let scale = volumes[j] * (color[j] - color[i]) * dw / r;
                grads[i][0] += scale * dx[0];
                grads[i][1] += scale * dx[1];
                grads[i][2] += scale * dx[2];
            }
        }
        let n_hats: Vec<[f64; 3]> = grads
            .iter()
            .map(|g| {
                let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                if mag > 1e-14 {
                    [g[0] / mag, g[1] / mag, g[2] / mag]
                } else {
                    [0.0; 3]
                }
            })
            .collect();
        let flux: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let ci = color[i];
                let fac = ci * (1.0 - ci);
                [fac * n_hats[i][0], fac * n_hats[i][1], fac * n_hats[i][2]]
            })
            .collect();
        let mut div_flux = vec![0.0_f64; n];
        let mut lap_c = vec![0.0_f64; n];
        for i in 0..n {
            for j in 0..n {
                let dx = [
                    positions[i][0] - positions[j][0],
                    positions[i][1] - positions[j][1],
                    positions[i][2] - positions[j][2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let r = r2.sqrt();
                if r < 1e-14 {
                    continue;
                }
                let dw = cubic_spline_dw(r, h);
                let gw = [dw / r * dx[0], dw / r * dx[1], dw / r * dx[2]];
                let df = [
                    flux[j][0] - flux[i][0],
                    flux[j][1] - flux[i][1],
                    flux[j][2] - flux[i][2],
                ];
                div_flux[i] += volumes[j] * (df[0] * gw[0] + df[1] * gw[1] + df[2] * gw[2]);
                let lap_w = cubic_spline_laplacian(r, h);
                lap_c[i] += volumes[j] * (color[j] - color[i]) * lap_w;
            }
        }
        (0..n)
            .map(|i| {
                let dc = dtau * self.epsilon * (div_flux[i] - lap_c[i]);
                (color[i] + dc).clamp(0.0, 1.0)
            })
            .collect()
    }
}
/// Model for wetting and contact-angle behaviour at solid–liquid–gas triple lines.
///
/// Implements the Young–Dupré equation:
/// ```text
/// cos θ = (σ_SG − σ_SL) / σ_LG
/// ```
/// where θ is the equilibrium contact angle, σ_SG the solid–gas surface energy,
/// σ_SL the solid–liquid surface energy, and σ_LG the liquid–gas surface tension.
#[derive(Debug, Clone)]
pub struct ContactAngleModel {
    /// Solid–gas surface energy (J/m²).
    pub sigma_sg: f64,
    /// Solid–liquid surface energy (J/m²).
    pub sigma_sl: f64,
    /// Liquid–gas surface tension (J/m²).
    pub sigma_lg: f64,
}
impl ContactAngleModel {
    /// Create a new contact-angle model.
    ///
    /// Uses water on glass as default values: θ ≈ 20°.
    pub fn new(sigma_sg: f64, sigma_sl: f64, sigma_lg: f64) -> Self {
        Self {
            sigma_sg,
            sigma_sl,
            sigma_lg,
        }
    }
    /// Equilibrium contact angle θ (radians) from Young's equation.
    ///
    /// Returns `None` if `|cos θ| > 1` (non-physical wettability parameters)
    /// or if `sigma_lg ≈ 0`.
    pub fn contact_angle_rad(&self) -> Option<f64> {
        if self.sigma_lg.abs() < 1e-30 {
            return None;
        }
        let cos_theta = (self.sigma_sg - self.sigma_sl) / self.sigma_lg;
        if cos_theta.abs() > 1.0 {
            None
        } else {
            Some(cos_theta.acos())
        }
    }
    /// Equilibrium contact angle θ in degrees.
    pub fn contact_angle_deg(&self) -> Option<f64> {
        self.contact_angle_rad().map(|rad| rad.to_degrees())
    }
    /// Classify the wetting regime for the contact angle θ.
    ///
    /// - θ < 10°   → superhydrophilic
    /// - θ ≤ 90°   → hydrophilic (wetting)
    /// - θ ≤ 150°  → hydrophobic (non-wetting)
    /// - θ > 150°  → superhydrophobic
    pub fn wetting_regime(&self) -> WettingRegime {
        match self.contact_angle_deg() {
            None => WettingRegime::NonPhysical,
            Some(deg) if deg < 10.0 => WettingRegime::Superhydrophilic,
            Some(deg) if deg <= 90.0 => WettingRegime::Hydrophilic,
            Some(deg) if deg <= 150.0 => WettingRegime::Hydrophobic,
            Some(_) => WettingRegime::Superhydrophobic,
        }
    }
    /// Spreading coefficient S = σ_SG − σ_SL − σ_LG.
    ///
    /// S ≥ 0 implies spontaneous spreading (complete wetting).
    pub fn spreading_coefficient(&self) -> f64 {
        self.sigma_sg - self.sigma_sl - self.sigma_lg
    }
    /// Adhesion work W_A = σ_LG (1 + cos θ) per unit area (J/m²).
    ///
    /// Returns 0 if the contact angle is non-physical.
    pub fn adhesion_work(&self) -> f64 {
        match self.contact_angle_rad() {
            None => 0.0,
            Some(theta) => self.sigma_lg * (1.0 + theta.cos()),
        }
    }
    /// Capillary pressure across a spherical meniscus of radius `r` (Pa).
    ///
    /// Uses the Young–Laplace equation: ΔP = 2 σ_LG cos θ / r.
    /// Returns `None` for non-physical contact angles or zero radius.
    pub fn capillary_pressure(&self, r: f64) -> Option<f64> {
        if r < 1e-30 {
            return None;
        }
        let theta = self.contact_angle_rad()?;
        Some(2.0 * self.sigma_lg * theta.cos() / r)
    }
}
/// Thermodynamic phase change model (liquid ↔ gas via temperature).
///
/// Determines the equilibrium phase from temperature and computes the latent
/// heat required for transitions.
pub struct PhaseChange {
    /// Boiling temperature (K) — liquid/gas transition.
    pub boiling_temperature: f64,
    /// Melting temperature (K) — solid/liquid transition.
    pub melting_temperature: f64,
    /// Latent heat of vaporisation (J/kg).
    pub latent_heat_vaporization: f64,
    /// Latent heat of fusion (J/kg).
    pub latent_heat_fusion: f64,
}
impl PhaseChange {
    /// Create a new phase change model.
    pub fn new(t_boil: f64, t_melt: f64, lv: f64, lf: f64) -> Self {
        Self {
            boiling_temperature: t_boil,
            melting_temperature: t_melt,
            latent_heat_vaporization: lv,
            latent_heat_fusion: lf,
        }
    }
    /// Determine the equilibrium phase from temperature.
    ///
    /// - `T < T_melt`  → [`Phase::Solid`]
    /// - `T_melt ≤ T < T_boil` → [`Phase::Liquid`]
    /// - `T ≥ T_boil`  → [`Phase::Gas`]
    pub fn phase_from_temperature(&self, temperature: f64) -> Phase {
        if temperature < self.melting_temperature {
            Phase::Solid
        } else if temperature < self.boiling_temperature {
            Phase::Liquid
        } else {
            Phase::Gas
        }
    }
    /// Energy (J/kg) required to transition from `from` to `to`.
    ///
    /// Endothermic transitions return a positive value; exothermic transitions
    /// return a negative value.  Returns `0.0` for same-phase or multi-step
    /// transitions that are not directly modelled.
    pub fn latent_heat(&self, from: Phase, to: Phase) -> f64 {
        match (from, to) {
            (Phase::Liquid, Phase::Gas) => self.latent_heat_vaporization,
            (Phase::Gas, Phase::Liquid) => -self.latent_heat_vaporization,
            (Phase::Solid, Phase::Liquid) => self.latent_heat_fusion,
            (Phase::Liquid, Phase::Solid) => -self.latent_heat_fusion,
            _ => 0.0,
        }
    }
    /// Check whether a phase transition occurs given an energy change.
    ///
    /// Returns the new phase and the residual energy after the transition has
    /// consumed (or released) the latent heat.  If no transition is triggered
    /// the current phase is returned unchanged and the full `energy_delta` is
    /// returned as residual.
    pub fn check_transition(
        &self,
        current_phase: Phase,
        temperature: f64,
        energy_delta: f64,
    ) -> (Phase, f64) {
        let equilibrium = self.phase_from_temperature(temperature);
        if equilibrium == current_phase {
            return (current_phase, energy_delta);
        }
        let lh = self.latent_heat(current_phase, equilibrium);
        let residual = energy_delta - lh;
        (equilibrium, residual)
    }
}
/// Mixture SPH: each particle carries fractional phase composition.
///
/// The three fractions `[liquid, gas, solid]` always sum to 1.  Effective
/// material properties are computed as volume-fraction-weighted averages.
pub struct MixtureSph {
    /// Per-particle phase fractions `[f_liquid, f_gas, f_solid]` (sum = 1).
    pub phase_fractions: Vec<[f64; 3]>,
}
impl MixtureSph {
    /// Create a new mixture system with `n_particles` particles, all
    /// initialised as pure liquid (`[1, 0, 0]`).
    pub fn new(n_particles: usize) -> Self {
        Self {
            phase_fractions: vec![[1.0, 0.0, 0.0]; n_particles],
        }
    }
    /// Set the phase fractions for particle `i`.
    ///
    /// The supplied values are normalised so that they sum to 1.  If all three
    /// values are zero the particle remains as pure liquid.
    pub fn set_fraction(&mut self, i: usize, liquid: f64, gas: f64, solid: f64) {
        let sum = liquid + gas + solid;
        if sum < 1e-30 {
            self.phase_fractions[i] = [1.0, 0.0, 0.0];
        } else {
            self.phase_fractions[i] = [liquid / sum, gas / sum, solid / sum];
        }
    }
    /// Effective density of particle `i` using volume-fraction weighting:
    /// `ρ = f_l ρ_l + f_g ρ_g + f_s ρ_s`.
    pub fn effective_density(&self, i: usize, rho_l: f64, rho_g: f64, rho_s: f64) -> f64 {
        let [fl, fg, fs] = self.phase_fractions[i];
        fl * rho_l + fg * rho_g + fs * rho_s
    }
    /// Effective dynamic viscosity of particle `i` using linear mixing of
    /// liquid and gas fractions (solid is treated as rigid and excluded).
    ///
    /// `μ = (f_l · μ_l + f_g · μ_g) / (f_l + f_g)`
    ///
    /// Returns `μ_l` when the particle is pure liquid and `0` when it is pure
    /// solid.
    pub fn effective_viscosity(&self, i: usize, mu_l: f64, mu_g: f64) -> f64 {
        let [fl, fg, _fs] = self.phase_fractions[i];
        let total = fl + fg;
        if total < 1e-30 {
            0.0
        } else {
            (fl * mu_l + fg * mu_g) / total
        }
    }
}
/// Phase identifier for a particle or material region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Liquid phase.
    Liquid,
    /// Gas / vapour phase.
    Gas,
    /// Solid phase.
    Solid,
}
/// A single particle in the immiscible multiphase SPH system.
#[derive(Debug, Clone)]
pub struct MultiphaseParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Smoothed density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Index into `MultiphaseSystem::phases`.
    pub phase_id: usize,
    /// Particle mass (kg).
    pub mass: f64,
    /// SPH-smoothed colour function `c = Σ_k V_k W_ik [phase_k = phase_i]` of the
    /// particle's own phase (about 1 inside the phase, about 0.5 on an interface).
    pub color_field: f64,
    /// Colour gradient `∇c` (1/m) of the particle's own phase.
    ///
    /// Points into the particle's own phase; its magnitude is the inverse
    /// interface thickness, so `σ κ ∇c` is the CSF force density.
    pub normal: [f64; 3],
    /// Interface curvature `κ = −∇·n̂` (1/m), evaluated in the frame of the
    /// particle's own phase: `+2/R` on the drop side of a spherical interface
    /// of radius `R` and `−2/R` on the ambient side.
    pub curvature: f64,
}
/// A fluid phase (e.g. water, oil) with physical properties and a visualisation colour.
pub struct FluidPhase {
    /// Unique identifier for this phase.
    pub id: usize,
    /// Reference (rest) density (kg/m³).
    pub density_ref: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub dynamic_viscosity: f64,
    /// Surface-tension coefficient σ (N/m).
    pub surface_tension: f64,
    /// RGB colour for visualisation.
    pub color: [f32; 3],
}
/// Analyser for Rayleigh–Taylor instability diagnostics.
///
/// Computes mixing-layer width, Atwood number, and growth-rate estimates from
/// a set of multiphase particles.
///
/// The Atwood number is defined as:
/// ```text
/// A = (ρ_heavy − ρ_light) / (ρ_heavy + ρ_light)
/// ```
/// and appears in the early-time linear growth rate `σ = sqrt(A g k)`, where
/// `g` is gravitational acceleration and `k` is the perturbation wavenumber.
#[derive(Debug, Clone)]
pub struct RayleighTaylorAnalyzer {
    /// Density of the heavy fluid (kg/m³).
    pub rho_heavy: f64,
    /// Density of the light fluid (kg/m³).
    pub rho_light: f64,
    /// Gravitational acceleration magnitude (m/s²).
    pub g: f64,
}
impl RayleighTaylorAnalyzer {
    /// Create a new analyser.
    pub fn new(rho_heavy: f64, rho_light: f64, g: f64) -> Self {
        Self {
            rho_heavy,
            rho_light,
            g,
        }
    }
    /// Atwood number A = (ρ_h − ρ_l) / (ρ_h + ρ_l).
    pub fn atwood_number(&self) -> f64 {
        (self.rho_heavy - self.rho_light) / (self.rho_heavy + self.rho_light)
    }
    /// Early-time linear instability growth rate σ = √(A g k) (s⁻¹).
    ///
    /// `k` is the perturbation wavenumber (rad/m).
    pub fn linear_growth_rate(&self, k: f64) -> f64 {
        (self.atwood_number() * self.g * k).sqrt()
    }
    /// Estimate the mixing-layer half-width h(t) = α A g t².
    ///
    /// Uses the standard turbulent-bubble formula with the empirical
    /// Youngs constant α ≈ 0.06 for bubbles and 0.035 for spikes.
    pub fn mixing_layer_width(&self, t: f64, alpha: f64) -> f64 {
        alpha * self.atwood_number() * self.g * t * t
    }
    /// Compute the centre-of-mass height of each phase from particle positions.
    ///
    /// `phase_ids[i]` should be 0 for heavy, 1 for light.  Returns
    /// `(y_com_heavy, y_com_light)`.
    pub fn phase_centers_of_mass(
        positions: &[[f64; 3]],
        phase_ids: &[usize],
        masses: &[f64],
    ) -> (f64, f64) {
        let mut sum_y_heavy = 0.0_f64;
        let mut mass_heavy = 0.0_f64;
        let mut sum_y_light = 0.0_f64;
        let mut mass_light = 0.0_f64;
        for ((pos, &pid), &m) in positions.iter().zip(phase_ids).zip(masses) {
            if pid == 0 {
                sum_y_heavy += m * pos[1];
                mass_heavy += m;
            } else {
                sum_y_light += m * pos[1];
                mass_light += m;
            }
        }
        let y_heavy = if mass_heavy > 1e-30 {
            sum_y_heavy / mass_heavy
        } else {
            0.0
        };
        let y_light = if mass_light > 1e-30 {
            sum_y_light / mass_light
        } else {
            0.0
        };
        (y_heavy, y_light)
    }
    /// Measure the mixing-layer width as the y-extent of the interfacial region.
    ///
    /// The "interfacial region" is defined as the band where the local heavy-phase
    /// volume fraction φ satisfies `0.05 < φ < 0.95`.  This is a column-averaged
    /// metric using a 1D binning along the y-axis.
    pub fn measured_mixing_width(
        positions: &[[f64; 3]],
        phase_ids: &[usize],
        n_bins: usize,
        y_min: f64,
        y_max: f64,
    ) -> f64 {
        if n_bins == 0 || y_max <= y_min {
            return 0.0;
        }
        let dy = (y_max - y_min) / n_bins as f64;
        let mut heavy_count = vec![0usize; n_bins];
        let mut total_count = vec![0usize; n_bins];
        for (pos, &pid) in positions.iter().zip(phase_ids) {
            let bin = ((pos[1] - y_min) / dy).floor() as isize;
            if bin < 0 || bin >= n_bins as isize {
                continue;
            }
            let bin = bin as usize;
            total_count[bin] += 1;
            if pid == 0 {
                heavy_count[bin] += 1;
            }
        }
        let mut y_lo = y_max;
        let mut y_hi = y_min;
        for b in 0..n_bins {
            if total_count[b] == 0 {
                continue;
            }
            let phi = heavy_count[b] as f64 / total_count[b] as f64;
            if phi > 0.05 && phi < 0.95 {
                let y_mid = y_min + (b as f64 + 0.5) * dy;
                if y_mid < y_lo {
                    y_lo = y_mid;
                }
                if y_mid > y_hi {
                    y_hi = y_mid;
                }
            }
        }
        if y_hi >= y_lo { y_hi - y_lo } else { 0.0 }
    }
}
/// Physical properties associated with a single phase particle.
pub struct PhaseParticle {
    /// Which phase this particle belongs to.
    pub phase: Phase,
    /// Reference (rest) density for this phase (kg/m³).
    pub density_ref: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Surface-tension coefficient (N/m).
    pub surface_tension: f64,
}
/// High-level immiscible multiphase SPH simulation driver.
pub struct MultiphaseSystem {
    /// All particles in the simulation.
    pub particles: Vec<MultiphaseParticle>,
    /// Registered fluid phases.
    pub phases: Vec<FluidPhase>,
    /// Smoothing length h (m).
    pub h: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: [f64; 3],
    /// Inter-phase surface-tension matrix.
    pub tension: InterphaseTension,
}
impl MultiphaseSystem {
    /// Create an empty simulation with the given smoothing length and gravity.
    pub fn new(h: f64, gravity: [f64; 3]) -> Self {
        Self {
            particles: Vec::new(),
            phases: Vec::new(),
            h,
            gravity,
            tension: InterphaseTension::new(0),
        }
    }
    /// Register a fluid phase; returns its assigned index.
    pub fn add_phase(&mut self, phase: FluidPhase) -> usize {
        let idx = self.phases.len();
        self.phases.push(phase);
        let n = self.phases.len();
        self.tension = InterphaseTension::new(n);
        idx
    }
    /// Add a particle at `pos` with initial velocity `vel` and mass `mass`
    /// belonging to `phase_id`.
    ///
    /// The mass must be consistent with the kernel normalisation: a lattice of
    /// spacing `s` fills space with `mass ≈ density_ref · s³`.  The initial
    /// density is the phase's reference density.
    ///
    /// Returns the index of the new particle.
    pub fn add_particle(
        &mut self,
        pos: [f64; 3],
        vel: [f64; 3],
        phase_id: usize,
        mass: f64,
    ) -> usize {
        let idx = self.particles.len();
        let density_ref = if phase_id < self.phases.len() {
            self.phases[phase_id].density_ref
        } else {
            DEFAULT_DENSITY_REF
        };
        self.particles.push(MultiphaseParticle {
            position: pos,
            velocity: vel,
            density: density_ref,
            pressure: 0.0,
            phase_id,
            mass,
            color_field: 0.0,
            normal: [0.0; 3],
            curvature: 0.0,
        });
        idx
    }
    /// Compute smoothed density for every particle: ρ_i = Σ_j m_j W(|x_i − x_j|, h).
    pub fn compute_density(&mut self) {
        let n = self.particles.len();
        let h = self.h;
        let mut densities = vec![0.0_f64; n];
        for (i, d_i) in densities.iter_mut().enumerate() {
            for pj in self.particles.iter() {
                let dx = [
                    self.particles[i].position[0] - pj.position[0],
                    self.particles[i].position[1] - pj.position[1],
                    self.particles[i].position[2] - pj.position[2],
                ];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                *d_i += pj.mass * cubic_spline_w(r, h);
            }
        }
        for (p, d) in self.particles.iter_mut().zip(densities.iter()) {
            p.density = *d;
        }
    }
    /// Compute pressure via Tait EOS: p = (ρ₀ c₀² / 7) · ((ρ/ρ₀)⁷ − 1).
    ///
    /// Uses c₀ = 100 m/s (corresponds to ≈10 × max velocity of 10 m/s).
    pub fn compute_pressure(&mut self) {
        let c0 = 100.0_f64;
        for p in &mut self.particles {
            let rho0 = if p.phase_id < self.phases.len() {
                self.phases[p.phase_id].density_ref
            } else {
                1000.0
            };
            let ratio = p.density / rho0;
            p.pressure = rho0 * c0 * c0 / 7.0 * (ratio.powi(7) - 1.0);
        }
    }
    /// Compute interface normals and curvature with the colour-field (CSF) approach.
    ///
    /// For every phase `p` present in the system the SPH-smoothed colour field
    /// `c^p_j = Σ_k V_k W_jk [phase_k = p]` (`V_k = m_k / ρ_k`) is formed and
    /// differentiated with the SPH gradient
    /// `∇c^p_i = Σ_j V_j (c^p_j − c^p_i) ∇_i W_ij`.  A particle `i` of phase `p`
    /// then receives
    ///
    /// * `normal = ∇c^p_i`, pointing into its own phase;
    /// * the unit normal field `n̂^p = ∇c^p / |∇c^p|` (zero where
    ///   `|∇c^p| ≤ 0.01 / h`);
    /// * `curvature = −Σ_j V_j (n̂^p_j − n̂^p_i) · ∇_i W_ij`, i.e. `κ = −∇·n̂`
    ///   with the neighbour normals taken from the same phase field, so both
    ///   sides of an interface are evaluated in a consistent frame.
    ///
    /// The result depends only on positions, masses, densities and phases, not
    /// on the storage order of the particles.  Densities must be up to date
    /// (see [`compute_density`](Self::compute_density)).  The cost is
    /// `O(P · N · k)` for `P` phases, `N` particles and `k` neighbours.
    pub fn compute_interface_normals(&mut self) {
        let n = self.particles.len();
        let h = self.h;
        if n == 0 {
            return;
        }
        if !(h > 0.0 && h.is_finite()) {
            for p in self.particles.iter_mut() {
                p.color_field = 0.0;
                p.normal = [0.0; 3];
                p.curvature = 0.0;
            }
            return;
        }

        let volumes: Vec<f64> = self
            .particles
            .iter()
            .map(|p| p.mass / p.density.max(1e-10))
            .collect();
        let neighbours = particle_neighbours(&self.particles, 2.0 * h);

        // Kernel values and gradients of every neighbour pair, shared by all phases.
        let mut weights: Vec<Vec<f64>> = Vec::with_capacity(n);
        let mut gradients: Vec<Vec<[f64; 3]>> = Vec::with_capacity(n);
        for (pi, list) in self.particles.iter().zip(neighbours.iter()) {
            let mut w_row = Vec::with_capacity(list.len());
            let mut g_row = Vec::with_capacity(list.len());
            for &j in list {
                let pj = &self.particles[j];
                let dx = [
                    pi.position[0] - pj.position[0],
                    pi.position[1] - pj.position[1],
                    pi.position[2] - pj.position[2],
                ];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                w_row.push(cubic_spline_w(r, h));
                g_row.push(kernel_grad(dx, r, h));
            }
            weights.push(w_row);
            gradients.push(g_row);
        }
        let self_weight = cubic_spline_w(0.0, h);
        let threshold = NORMAL_THRESHOLD / h;

        let mut phase_ids: Vec<usize> = self.particles.iter().map(|p| p.phase_id).collect();
        phase_ids.sort_unstable();
        phase_ids.dedup();

        let mut colours = vec![0.0_f64; n];
        let mut normals = vec![[0.0_f64; 3]; n];
        let mut curvatures = vec![0.0_f64; n];
        for &phase in &phase_ids {
            let colour: Vec<f64> = (0..n)
                .map(|j| {
                    let mut c = if self.particles[j].phase_id == phase {
                        volumes[j] * self_weight
                    } else {
                        0.0
                    };
                    for (&k, &w) in neighbours[j].iter().zip(weights[j].iter()) {
                        if self.particles[k].phase_id == phase {
                            c += volumes[k] * w;
                        }
                    }
                    c
                })
                .collect();
            let colour_gradient: Vec<[f64; 3]> = (0..n)
                .map(|i| {
                    let mut g = [0.0_f64; 3];
                    for (&j, gw) in neighbours[i].iter().zip(gradients[i].iter()) {
                        let f = volumes[j] * (colour[j] - colour[i]);
                        g[0] += f * gw[0];
                        g[1] += f * gw[1];
                        g[2] += f * gw[2];
                    }
                    g
                })
                .collect();
            let unit_normal: Vec<[f64; 3]> = colour_gradient
                .iter()
                .map(|g| {
                    let magnitude = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    if magnitude > threshold {
                        [g[0] / magnitude, g[1] / magnitude, g[2] / magnitude]
                    } else {
                        [0.0; 3]
                    }
                })
                .collect();

            for i in (0..n).filter(|&i| self.particles[i].phase_id == phase) {
                let mut divergence = 0.0_f64;
                for (&j, gw) in neighbours[i].iter().zip(gradients[i].iter()) {
                    let dn = [
                        unit_normal[j][0] - unit_normal[i][0],
                        unit_normal[j][1] - unit_normal[i][1],
                        unit_normal[j][2] - unit_normal[i][2],
                    ];
                    divergence += volumes[j] * (dn[0] * gw[0] + dn[1] * gw[1] + dn[2] * gw[2]);
                }
                colours[i] = colour[i];
                normals[i] = colour_gradient[i];
                curvatures[i] = -divergence;
            }
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.color_field = colours[i];
            p.normal = normals[i];
            p.curvature = curvatures[i];
        }
    }
    /// Compute per-particle accelerations: pressure + viscosity + gravity + surface tension.
    ///
    /// The surface tension is the continuum surface force `σ κ ∇c / ρ` built from
    /// the normals and curvatures of [`compute_interface_normals`](Self::compute_interface_normals);
    /// `σ` is the interphase coefficient, averaged over the neighbours of other
    /// phases with the weights `V_j W_ij`.
    ///
    /// Returns a `Vec<[f64; 3]>` of accelerations, one per particle.
    pub fn compute_forces(&self) -> Vec<[f64; 3]> {
        let n = self.particles.len();
        let h = self.h;
        let mut acc = vec![[0.0_f64; 3]; n];
        for (i, acc_i) in acc.iter_mut().enumerate() {
            let pi = &self.particles[i];
            let rho_i = pi.density.max(1e-10);
            let p_i = pi.pressure;
            let mu_i = if pi.phase_id < self.phases.len() {
                self.phases[pi.phase_id].dynamic_viscosity
            } else {
                0.0
            };
            let mut sigma_weighted = 0.0_f64;
            let mut sigma_weight = 0.0_f64;
            for (j, pj) in self.particles.iter().enumerate() {
                if i == j {
                    continue;
                }
                let rho_j = pj.density.max(1e-10);
                let p_j = pj.pressure;
                let mu_j = if pj.phase_id < self.phases.len() {
                    self.phases[pj.phase_id].dynamic_viscosity
                } else {
                    0.0
                };
                let dx = [
                    pi.position[0] - pj.position[0],
                    pi.position[1] - pj.position[1],
                    pi.position[2] - pj.position[2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let r = r2.sqrt();
                if r >= 2.0 * h {
                    continue;
                }
                let grad = kernel_grad(dx, r, h);
                let pres_coeff = -pj.mass * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j));
                acc_i[0] += pres_coeff * grad[0];
                acc_i[1] += pres_coeff * grad[1];
                acc_i[2] += pres_coeff * grad[2];
                let mu_ij = if mu_i + mu_j > 1e-30 {
                    2.0 * mu_i * mu_j / (mu_i + mu_j)
                } else {
                    0.0
                };
                let dv = [
                    pi.velocity[0] - pj.velocity[0],
                    pi.velocity[1] - pj.velocity[1],
                    pi.velocity[2] - pj.velocity[2],
                ];
                let v_dot_r = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
                let visc_scale = mu_ij * pj.mass * v_dot_r / (rho_i * rho_j * (r2 + 0.01 * h * h));
                acc_i[0] += visc_scale * grad[0];
                acc_i[1] += visc_scale * grad[1];
                acc_i[2] += visc_scale * grad[2];
                if pi.phase_id != pj.phase_id {
                    let n_phases = self.phases.len();
                    if pi.phase_id < n_phases && pj.phase_id < n_phases {
                        let weight = pj.mass / rho_j * cubic_spline_w(r, h);
                        sigma_weighted += weight * self.tension.get(pi.phase_id, pj.phase_id);
                        sigma_weight += weight;
                    }
                }
            }
            // Continuum surface force, applied once per particle:
            // a_i = σ_i κ_i ∇c_i / ρ_i with σ_i the tension of the interface seen by i.
            if sigma_weight > 0.0 {
                let sigma_i = sigma_weighted / sigma_weight;
                let csf =
                    InterphaseTension::csf_acceleration(sigma_i, pi.curvature, pi.normal, rho_i);
                acc_i[0] += csf[0];
                acc_i[1] += csf[1];
                acc_i[2] += csf[2];
            }
            acc_i[0] += self.gravity[0];
            acc_i[1] += self.gravity[1];
            acc_i[2] += self.gravity[2];
        }
        acc
    }
    /// Advance the simulation by `dt` seconds (leapfrog integrator).
    ///
    /// Order: density → pressure → interface normals → forces → velocity → position.
    pub fn step(&mut self, dt: f64) {
        self.compute_density();
        self.compute_pressure();
        self.compute_interface_normals();
        let acc = self.compute_forces();
        for (p, a) in self.particles.iter_mut().zip(acc.iter()) {
            p.velocity[0] += a[0] * dt;
            p.velocity[1] += a[1] * dt;
            p.velocity[2] += a[2] * dt;
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
        }
    }
    /// Phase-separation index in \[0, 1\]: 0 = fully mixed, 1 = fully separated.
    ///
    /// Computed as 1 − (mean cross-phase neighbour fraction).
    pub fn phase_separation_index(&self) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        let h = self.h;
        let mut total_cross = 0.0_f64;
        let mut total_neighbours = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = [
                    self.particles[i].position[0] - self.particles[j].position[0],
                    self.particles[i].position[1] - self.particles[j].position[1],
                    self.particles[i].position[2] - self.particles[j].position[2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                if r2 < (2.0 * h) * (2.0 * h) {
                    total_neighbours += 1.0;
                    if self.particles[i].phase_id != self.particles[j].phase_id {
                        total_cross += 1.0;
                    }
                }
            }
        }
        if total_neighbours < 1.0 {
            return 1.0;
        }
        1.0 - total_cross / total_neighbours
    }
    /// Build a two-layer particle lattice: lower half = phase_a, upper half = phase_b.
    ///
    /// The `nx × ny` lattice lies in the plane `z = 0`.  Because the SPH kernel
    /// is three-dimensional, the particle masses are calibrated so that the
    /// kernel-summed density at an interior lattice site equals the reference
    /// density of its phase: `m = ρ_ref / Σ_{a,b} W(s √(a² + b²), h)` with the
    /// spacing `s` and this system's smoothing length `h`.  (A mass of
    /// `ρ_ref · s³` would only fill a three-dimensional lattice.)  Phases that
    /// have not been registered use a reference density of 1000 kg/m³.
    ///
    /// Returns a flat `Vec<MultiphaseParticle>` (does *not* add them to `self`),
    /// or an empty vector if `spacing` or the smoothing length is not a
    /// positive finite number.
    pub fn initialize_two_layer(
        &self,
        nx: usize,
        ny: usize,
        phase_a: usize,
        phase_b: usize,
        spacing: f64,
    ) -> Vec<MultiphaseParticle> {
        let valid = |x: f64| x > 0.0 && x.is_finite();
        if !valid(spacing) || !valid(self.h) {
            return Vec::new();
        }
        let kernel_sum = planar_lattice_kernel_sum(spacing, self.h);
        let density_ref = |phase_id: usize| {
            self.phases
                .get(phase_id)
                .map_or(DEFAULT_DENSITY_REF, |phase| phase.density_ref)
        };
        let mut particles = Vec::with_capacity(nx * ny);
        let half = ny / 2;
        for row in 0..ny {
            for col in 0..nx {
                let x = col as f64 * spacing;
                let y = row as f64 * spacing;
                let phase_id = if row < half { phase_a } else { phase_b };
                let rho = density_ref(phase_id);
                particles.push(MultiphaseParticle {
                    position: [x, y, 0.0],
                    velocity: [0.0; 3],
                    density: rho,
                    pressure: 0.0,
                    phase_id,
                    mass: rho / kernel_sum,
                    color_field: 0.0,
                    normal: [0.0; 3],
                    curvature: 0.0,
                });
            }
        }
        particles
    }
}
/// Symmetric matrix of inter-phase surface-tension coefficients σ_ij.
pub struct InterphaseTension {
    /// `sigma_ij[i][j]` holds σ between phase i and phase j.
    pub sigma_ij: Vec<Vec<f64>>,
}
impl InterphaseTension {
    /// Create a zero-initialised n×n tension matrix.
    pub fn new(n_phases: usize) -> Self {
        Self {
            sigma_ij: vec![vec![0.0; n_phases]; n_phases],
        }
    }
    /// Set the surface-tension coefficient between phases i and j (symmetric).
    pub fn set(&mut self, i: usize, j: usize, sigma: f64) {
        self.sigma_ij[i][j] = sigma;
        self.sigma_ij[j][i] = sigma;
    }
    /// Get the surface-tension coefficient between phases i and j.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.sigma_ij[i][j]
    }
    /// Continuum-surface-force (CSF) acceleration of a particle,
    /// `a = σ κ ∇c / ρ` (m/s²).
    ///
    /// # Arguments
    /// * `sigma`     – interface tension coefficient σ (N/m)
    /// * `curvature` – interface curvature κ (1/m) in the frame of the particle's own phase
    /// * `normal`    – colour gradient `∇c` (1/m) pointing into the particle's own phase
    /// * `density`   – particle density ρ (kg/m³)
    ///
    /// κ and `∇c` change sign together when the roles of the two phases are
    /// exchanged, so both sides of an interface are pulled towards the centre
    /// of curvature.
    pub fn csf_acceleration(
        sigma: f64,
        curvature: f64,
        normal: [f64; 3],
        density: f64,
    ) -> [f64; 3] {
        let scale = sigma * curvature / density.max(1e-10);
        [scale * normal[0], scale * normal[1], scale * normal[2]]
    }
}
/// Simple connected-component labelling for phase morphology analysis.
///
/// Groups particles of the same phase into spatially connected clusters using
/// a union-find data structure.  Two same-phase particles are considered
/// connected if their separation is ≤ `h`.
#[derive(Debug, Clone)]
pub struct PhaseMorphologyAnalyzer {
    /// Smoothing length used for connectivity (m).
    pub h: f64,
}
impl PhaseMorphologyAnalyzer {
    /// Create a new analyser with the given connectivity radius.
    pub fn new(h: f64) -> Self {
        Self { h }
    }
    /// Label connected components of same-phase particles.
    ///
    /// Returns a `labels` vector where `labels[i]` is the component ID of
    /// particle `i`.  Components are numbered 0, 1, 2, … in order of first
    /// appearance.
    pub fn label_components(&self, positions: &[[f64; 3]], phase_ids: &[usize]) -> Vec<usize> {
        let n = positions.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let h2 = self.h * self.h;
        for i in 0..n {
            for j in (i + 1)..n {
                if phase_ids[i] != phase_ids[j] {
                    continue;
                }
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                if dx * dx + dy * dy + dz * dz <= h2 {
                    let ri = find(&mut parent, i);
                    let rj = find(&mut parent, j);
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }
        let mut root_to_label: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        let mut next_label = 0usize;
        let mut labels = vec![0usize; n];
        for (i, label_i) in labels.iter_mut().enumerate() {
            let root = find(&mut parent, i);
            *label_i = *root_to_label.entry(root).or_insert_with(|| {
                let l = next_label;
                next_label += 1;
                l
            });
        }
        labels
    }
    /// Count the number of connected components for a given phase.
    pub fn component_count(
        &self,
        positions: &[[f64; 3]],
        phase_ids: &[usize],
        target_phase: usize,
    ) -> usize {
        let n = positions.len();
        let filtered_pos: Vec<[f64; 3]> = (0..n)
            .filter(|&i| phase_ids[i] == target_phase)
            .map(|i| positions[i])
            .collect();
        let filtered_phase: Vec<usize> = vec![target_phase; filtered_pos.len()];
        if filtered_pos.is_empty() {
            return 0;
        }
        let labels = self.label_components(&filtered_pos, &filtered_phase);
        let max_label = labels.iter().cloned().max().unwrap_or(0);
        max_label + 1
    }
    /// Compute the size (particle count) of each connected component.
    ///
    /// Returns a sorted `Vec`usize` of component sizes in descending order.
    pub fn component_sizes(&self, positions: &[[f64; 3]], phase_ids: &[usize]) -> Vec<usize> {
        let labels = self.label_components(positions, phase_ids);
        let max_label = labels.iter().cloned().max().unwrap_or(0);
        let mut sizes = vec![0usize; max_label + 1];
        for &l in &labels {
            sizes[l] += 1;
        }
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        sizes
    }
    /// Compute the radius of gyration of component `component_id`.
    ///
    /// `R_g = sqrt(1/N Σ |r_i − r_cm|²)` where `r_cm` is the centre of mass.
    /// Returns `None` if the component does not exist or has fewer than 2 particles.
    pub fn radius_of_gyration(
        &self,
        positions: &[[f64; 3]],
        phase_ids: &[usize],
        component_id: usize,
    ) -> Option<f64> {
        let labels = self.label_components(positions, phase_ids);
        let members: Vec<[f64; 3]> = positions
            .iter()
            .zip(labels.iter())
            .filter(|(_, l)| **l == component_id)
            .map(|(&p, _)| p)
            .collect();
        let n = members.len();
        if n < 2 {
            return None;
        }
        let cm = [
            members.iter().map(|p| p[0]).sum::<f64>() / n as f64,
            members.iter().map(|p| p[1]).sum::<f64>() / n as f64,
            members.iter().map(|p| p[2]).sum::<f64>() / n as f64,
        ];
        let rg2 = members
            .iter()
            .map(|p: &[f64; 3]| {
                let dx = p[0] - cm[0];
                let dy = p[1] - cm[1];
                let dz = p[2] - cm[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            / n as f64;
        Some(rg2.sqrt())
    }
}
/// Pairwise coalescence detector for a set of droplets.
///
/// Two droplets are deemed candidates for coalescence when the gap between
/// their surfaces is smaller than `threshold`.  The gap is
/// `|c_i − c_j| − r_i − r_j`.
pub struct CoalescenceDetector {
    /// Gap threshold below which coalescence is triggered (m).
    pub gap_threshold: f64,
}
impl CoalescenceDetector {
    /// Create a new detector with the given gap threshold.
    pub fn new(gap_threshold: f64) -> Self {
        Self { gap_threshold }
    }
    /// Find all coalescence pairs from a list of droplet descriptors.
    ///
    /// `centers\[i\]` is the 3D centre of droplet *i*; `radii\[i\]` its effective
    /// radius.  Returns a list of `(i, j)` pairs that should coalesce.
    pub fn find_pairs(&self, centers: &[[f64; 3]], radii: &[f64]) -> Vec<(usize, usize)> {
        let n = centers.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = centers[i][0] - centers[j][0];
                let dy = centers[i][1] - centers[j][1];
                let dz = centers[i][2] - centers[j][2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let gap = dist - radii[i] - radii[j];
                if gap < self.gap_threshold {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
    /// Merge two droplets conserving volume: R_merged = (R₁³ + R₂³)^(1/3).
    pub fn merged_radius(r1: f64, r2: f64) -> f64 {
        (r1.powi(3) + r2.powi(3)).cbrt()
    }
    /// Merged centre of mass from two droplets with masses m1, m2.
    pub fn merged_center(c1: [f64; 3], m1: f64, c2: [f64; 3], m2: f64) -> [f64; 3] {
        let total = m1 + m2;
        if total < 1e-30 {
            return c1;
        }
        [
            (m1 * c1[0] + m2 * c2[0]) / total,
            (m1 * c1[1] + m2 * c2[1]) / total,
            (m1 * c1[2] + m2 * c2[2]) / total,
        ]
    }
}
/// Tracks liquid droplets (connected components of phase 0) over time.
///
/// Each droplet is characterised by its centre of mass, radius of gyration,
/// and particle count.
#[derive(Debug, Clone)]
pub struct DropletInfo {
    /// Droplet label (component ID at the current time step).
    pub label: usize,
    /// Number of particles in this droplet.
    pub particle_count: usize,
    /// Centre of mass position.
    pub center: [f64; 3],
    /// Radius of gyration.
    pub radius_of_gyration: f64,
}
/// Simplified van der Waals equation of state for vapor–liquid equilibrium.
///
/// `(p + a/v²)(v − b) = R_specific T`
///
/// where `v = 1/ρ` is the specific volume.  This is solved for pressure given
/// density and temperature.
#[derive(Debug, Clone)]
pub struct VanDerWaalsEos {
    /// Attraction parameter a (Pa·m⁶/kg²).
    pub a: f64,
    /// Co-volume parameter b (m³/kg).
    pub b: f64,
    /// Specific gas constant R_s = R/M (J/(kg·K)).
    pub r_specific: f64,
}
impl VanDerWaalsEos {
    /// Construct a vdW EOS for a given substance.
    ///
    /// Critical properties: T_c, p_c.  The vdW parameters follow:
    /// `a = 27 R_s² T_c² / (64 p_c)`,  `b = R_s T_c / (8 p_c)`.
    pub fn from_critical(t_c: f64, p_c: f64, r_specific: f64) -> Self {
        let a = 27.0 * r_specific * r_specific * t_c * t_c / (64.0 * p_c);
        let b = r_specific * t_c / (8.0 * p_c);
        Self { a, b, r_specific }
    }
    /// Pressure from density and temperature.
    ///
    /// Returns `None` if `ρ b ≥ 1` (unphysical – compressed beyond close packing).
    pub fn pressure(&self, rho: f64, temperature: f64) -> Option<f64> {
        let v = 1.0 / rho;
        let denom = v - self.b;
        if denom <= 0.0 {
            return None;
        }
        Some(self.r_specific * temperature / denom - self.a * rho * rho)
    }
    /// Speed of sound c = √(dp/dρ) at constant temperature.
    pub fn speed_of_sound(&self, rho: f64, temperature: f64) -> Option<f64> {
        let v = 1.0 / rho;
        let denom = v - self.b;
        if denom <= 0.0 {
            return None;
        }
        let dp_dv =
            -self.r_specific * temperature / (denom * denom) + 2.0 * self.a * rho * rho * rho;
        let dp_drho = -dp_dv / (rho * rho);
        if dp_drho < 0.0 {
            None
        } else {
            Some(dp_drho.sqrt())
        }
    }
    /// Check if the state (ρ, T) is in the mechanically unstable spinodal region.
    ///
    /// `(∂p/∂ρ)_T < 0` implies spinodal decomposition.
    pub fn in_spinodal_region(&self, rho: f64, temperature: f64) -> bool {
        self.speed_of_sound(rho, temperature).is_none()
    }
    /// Reduced temperature T_r = T / T_c.
    pub fn reduced_temperature(&self, temperature: f64) -> f64 {
        let t_c = 8.0 * self.a / (27.0 * self.r_specific * self.b);
        if t_c < 1e-30 { 0.0 } else { temperature / t_c }
    }
}
/// Wetting regime classification based on contact angle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WettingRegime {
    /// θ < 10°: complete or near-complete wetting.
    Superhydrophilic,
    /// 10° ≤ θ ≤ 90°: partial wetting.
    Hydrophilic,
    /// 90° < θ ≤ 150°: non-wetting.
    Hydrophobic,
    /// θ > 150°: extreme non-wetting.
    Superhydrophobic,
    /// Contact angle not physically achievable with given parameters.
    NonPhysical,
}
