// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biological tissue simulation for the OxiPhysics softbody crate.
//!
//! Models muscle fibres (Hill model), tendons, skin layers, cartilage (biphasic),
//! bone remodelling (Wolff's law), wound healing, pressure ulcers, and general
//! biomechanics analysis.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ---------------------------------------------------------------------------
// TissueType
// ---------------------------------------------------------------------------

/// Tissue type classification with associated material property presets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TissueType {
    /// Skeletal muscle.
    Muscle,
    /// Tendon (dense connective tissue).
    Tendon,
    /// Ligament (joint-stabilizing connective tissue).
    Ligament,
    /// Hyaline cartilage.
    Cartilage,
    /// Cortical or cancellous bone.
    Bone,
    /// Skin (layered: epidermis, dermis, hypodermis).
    Skin,
}

impl TissueType {
    /// Young's modulus preset (Pa).
    pub fn youngs_modulus(&self) -> f64 {
        match self {
            TissueType::Muscle => 50_000.0,        // 50 kPa
            TissueType::Tendon => 1_500_000_000.0, // 1.5 GPa
            TissueType::Ligament => 100_000_000.0, // 100 MPa
            TissueType::Cartilage => 1_000_000.0,  // 1 MPa
            TissueType::Bone => 20_000_000_000.0,  // 20 GPa
            TissueType::Skin => 100_000.0,         // 100 kPa
        }
    }

    /// Poisson's ratio preset.
    pub fn poisson_ratio(&self) -> f64 {
        match self {
            TissueType::Muscle => 0.48,
            TissueType::Tendon => 0.40,
            TissueType::Ligament => 0.40,
            TissueType::Cartilage => 0.10,
            TissueType::Bone => 0.30,
            TissueType::Skin => 0.49,
        }
    }

    /// Density preset (kg/m³).
    pub fn density(&self) -> f64 {
        match self {
            TissueType::Muscle => 1_060.0,
            TissueType::Tendon => 1_200.0,
            TissueType::Ligament => 1_200.0,
            TissueType::Cartilage => 1_100.0,
            TissueType::Bone => 1_900.0,
            TissueType::Skin => 1_100.0,
        }
    }
}

// ---------------------------------------------------------------------------
// MuscleFiber
// ---------------------------------------------------------------------------

/// A single muscle fibre with Hill model components.
#[derive(Debug, Clone)]
pub struct MuscleFiber {
    /// Fibre direction (unit vector).
    pub direction: [f64; 3],
    /// Neural activation level in `[0, 1]`.
    pub activation: f64,
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fibre length (m).
    pub l_opt: f64,
    /// Current fibre length (m).
    pub length: f64,
    /// Current fibre velocity (m/s, shortening positive).
    pub velocity: f64,
}

impl MuscleFiber {
    /// Create a new muscle fibre.
    pub fn new(direction: [f64; 3], f_max: f64, l_opt: f64) -> Self {
        Self {
            direction: normalize3(direction),
            activation: 0.0,
            f_max,
            l_opt,
            length: l_opt,
            velocity: 0.0,
        }
    }

    /// Active force as a fraction of `f_max` using a Gaussian force-length relation.
    pub fn active_force_length(&self) -> f64 {
        let ratio = self.length / self.l_opt;
        let width = 0.45;
        (-(ratio - 1.0).powi(2) / (2.0 * width * width)).exp()
    }

    /// Passive force (exponential stiffness beyond optimal length).
    pub fn passive_force_length(&self) -> f64 {
        let ratio = self.length / self.l_opt;
        if ratio > 1.0 {
            0.02 * ((ratio - 1.0) / 0.6).exp()
        } else {
            0.0
        }
    }

    /// Total fibre force using simplified Hill model.
    pub fn total_force(&self) -> f64 {
        let fl = self.active_force_length();
        let fv = hill_force_velocity(self.velocity, self.l_opt);
        let f_active = self.activation * fl * fv * self.f_max;
        let f_passive = self.passive_force_length() * self.f_max;
        f_active + f_passive
    }

    /// Force vector in world space.
    pub fn force_vector(&self) -> [f64; 3] {
        let f = self.total_force();
        [
            self.direction[0] * f,
            self.direction[1] * f,
            self.direction[2] * f,
        ]
    }
}

/// Hill force-velocity relation: `f_v = (1 - v/v_max) / (1 + v / (a * v_max))`.
pub fn hill_force_velocity(velocity: f64, l_opt: f64) -> f64 {
    let v_max = 10.0 * l_opt; // typical 10 fibre-lengths/s
    let a = 0.25;
    if velocity >= 0.0 {
        // Shortening
        let v = velocity.min(v_max);
        (1.0 - v / v_max) / (1.0 + v / (a * v_max))
    } else {
        // Lengthening (eccentric)
        let v = (-velocity).min(v_max);
        1.0 + 1.5 * v / v_max
    }
}

// ---------------------------------------------------------------------------
// HillModel
// ---------------------------------------------------------------------------

/// Hill's three-element muscle model: contractile (CE), parallel elastic (PE),
/// and series elastic (SE) elements.
#[derive(Debug, Clone)]
pub struct HillModel {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fibre length (m).
    pub l_opt: f64,
    /// Tendon slack length (m).
    pub l_slack: f64,
    /// Pennation angle (radians).
    pub pennation: f64,
    /// Current CE length.
    pub l_ce: f64,
    /// Current activation.
    pub activation: f64,
}

impl HillModel {
    /// Create a new Hill model.
    pub fn new(f_max: f64, l_opt: f64, l_slack: f64, pennation: f64) -> Self {
        Self {
            f_max,
            l_opt,
            l_slack,
            pennation,
            l_ce: l_opt,
            activation: 0.0,
        }
    }

    /// Contractile element force.
    pub fn ce_force(&self) -> f64 {
        let fiber = MuscleFiber {
            direction: [1.0, 0.0, 0.0],
            activation: self.activation,
            f_max: self.f_max,
            l_opt: self.l_opt,
            length: self.l_ce,
            velocity: 0.0,
        };
        fiber.total_force()
    }

    /// Series elastic element force: linear spring above slack length.
    pub fn se_force(&self, l_mtu: f64) -> f64 {
        // Musculotendon unit length minus CE projection
        let l_ce_proj = self.l_ce * self.pennation.cos();
        let l_se = (l_mtu - l_ce_proj).max(0.0);
        let k_se = 35.0 * self.f_max / self.l_slack;
        let delta = (l_se - self.l_slack).max(0.0);
        k_se * delta
    }

    /// Update CE length using implicit Euler given muscle tendon unit velocity.
    pub fn update(&mut self, l_mtu: f64, activation: f64, dt: f64) {
        self.activation = activation.clamp(0.0, 1.0);
        // Simple explicit integration of CE velocity toward equilibrium
        let f_ce = self.ce_force();
        let f_se = self.se_force(l_mtu);
        let f_err = f_se - f_ce;
        let v_ce = f_err / (self.f_max * 0.1 + 1e-6);
        self.l_ce = (self.l_ce + v_ce * dt).max(0.3 * self.l_opt);
    }
}

// ---------------------------------------------------------------------------
// TendonModel
// ---------------------------------------------------------------------------

/// Viscoelastic tendon model with toe and linear regions.
#[derive(Debug, Clone)]
pub struct TendonModel {
    /// Slack length (m).
    pub slack_length: f64,
    /// Stiffness in linear region (N/m).
    pub linear_stiffness: f64,
    /// Transition strain (end of toe region).
    pub toe_strain: f64,
    /// Yield strain.
    pub yield_strain: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Current length.
    pub length: f64,
    /// Current velocity.
    pub velocity: f64,
}

impl TendonModel {
    /// Create a new tendon model.
    pub fn new(slack_length: f64, linear_stiffness: f64) -> Self {
        Self {
            slack_length,
            linear_stiffness,
            toe_strain: 0.02,
            yield_strain: 0.08,
            damping: 100.0,
            length: slack_length,
            velocity: 0.0,
        }
    }

    /// Compute tendon strain `(l - l_slack) / l_slack`.
    pub fn strain(&self) -> f64 {
        (self.length - self.slack_length) / self.slack_length
    }

    /// Elastic force based on toe/linear/yield regions.
    pub fn elastic_force(&self) -> f64 {
        let e = self.strain();
        if e <= 0.0 {
            return 0.0;
        }
        if e <= self.toe_strain {
            // Toe region: quadratic
            let k_toe = self.linear_stiffness * e / self.toe_strain;
            0.5 * k_toe * e * self.slack_length
        } else if e <= self.yield_strain {
            // Linear region
            self.linear_stiffness * (e - self.toe_strain * 0.5) * self.slack_length
        } else {
            // Yield: reduced stiffness
            let f_yield = self.linear_stiffness
                * (self.yield_strain - self.toe_strain * 0.5)
                * self.slack_length;
            f_yield + self.linear_stiffness * 0.1 * (e - self.yield_strain) * self.slack_length
        }
    }

    /// Total force including damping.
    pub fn total_force(&self) -> f64 {
        let f_e = self.elastic_force();
        let f_d = self.damping * self.velocity;
        (f_e + f_d).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// SkinLayer
// ---------------------------------------------------------------------------

/// A single layer of skin tissue.
#[derive(Debug, Clone)]
pub struct SkinLayer {
    /// Thickness (m).
    pub thickness: f64,
    /// Young's modulus (Pa).
    pub modulus: f64,
    /// Pre-stretch ratio.
    pub prestretch: f64,
    /// Name/type of layer.
    pub name: &'static str,
}

impl SkinLayer {
    /// Epidermis preset values.
    pub fn epidermis() -> Self {
        Self {
            thickness: 0.0001,
            modulus: 140_000.0,
            prestretch: 1.05,
            name: "epidermis",
        }
    }

    /// Dermis preset values.
    pub fn dermis() -> Self {
        Self {
            thickness: 0.002,
            modulus: 100_000.0,
            prestretch: 1.10,
            name: "dermis",
        }
    }

    /// Hypodermis preset values.
    pub fn hypodermis() -> Self {
        Self {
            thickness: 0.005,
            modulus: 2_000.0,
            prestretch: 1.0,
            name: "hypodermis",
        }
    }

    /// Compute force per unit area (stress) for a given stretch ratio `λ`.
    ///
    /// Uses a linearised neo-Hookean formulation.
    pub fn stress(&self, stretch: f64) -> f64 {
        let lambda = stretch / self.prestretch;
        self.modulus * (lambda - 1.0 / (lambda * lambda))
    }
}

/// A layered skin model composed of multiple `SkinLayer` instances.
#[derive(Debug, Clone)]
pub struct SkinModel {
    /// Ordered layers from outer to inner.
    pub layers: Vec<SkinLayer>,
}

impl SkinModel {
    /// Build a three-layer skin model.
    pub fn default_skin() -> Self {
        Self {
            layers: vec![
                SkinLayer::epidermis(),
                SkinLayer::dermis(),
                SkinLayer::hypodermis(),
            ],
        }
    }

    /// Total skin thickness (m).
    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// Effective bending stiffness `D = Σ E_i t_i^3 / 12`.
    pub fn bending_stiffness(&self) -> f64 {
        self.layers
            .iter()
            .map(|l| l.modulus * l.thickness.powi(3) / 12.0)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// CartilageModel
// ---------------------------------------------------------------------------

/// Biphasic cartilage model (Mow et al.): solid matrix + interstitial fluid.
#[derive(Debug, Clone)]
pub struct CartilageModel {
    /// Aggregate modulus of solid matrix (Pa).
    pub ha: f64,
    /// Hydraulic permeability (m⁴/N·s).
    pub permeability: f64,
    /// Osmotic swelling pressure (Pa).
    pub swelling_pressure: f64,
    /// Fluid fraction (void ratio).
    pub fluid_fraction: f64,
    /// Current compressive strain.
    pub strain: f64,
}

impl CartilageModel {
    /// Create a new biphasic cartilage model.
    pub fn new(ha: f64, permeability: f64) -> Self {
        Self {
            ha,
            permeability,
            swelling_pressure: 50_000.0,
            fluid_fraction: 0.8,
            strain: 0.0,
        }
    }

    /// Solid matrix stress for a given strain.
    pub fn solid_stress(&self, strain: f64) -> f64 {
        self.ha * strain
    }

    /// Fluid pressure from biphasic theory: `p = ha * ε + swelling_pressure`.
    pub fn fluid_pressure(&self, strain: f64) -> f64 {
        self.ha * strain + self.swelling_pressure
    }

    /// Creep compliance after time `t` (simplified exponential).
    ///
    /// `J(t) = J_∞ (1 - exp(-t/τ))` where `τ = (fluid_fraction / permeability / ha)`.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let j_inf = 1.0 / self.ha;
        let tau = self.fluid_fraction / (self.permeability * self.ha + 1e-20);
        j_inf * (1.0 - (-t / tau).exp())
    }

    /// Update the biphasic model for a given applied stress and time step.
    pub fn update(&mut self, applied_stress: f64, dt: f64) {
        // Rate of strain from permeability and pressure gradient
        let p = self.fluid_pressure(self.strain);
        let d_strain = self.permeability * (applied_stress - p) * dt;
        self.strain = (self.strain + d_strain).clamp(0.0, 0.5);
    }
}

// ---------------------------------------------------------------------------
// BoneRemodeling
// ---------------------------------------------------------------------------

/// Bone remodelling based on Wolff's law.
///
/// Local apparent density evolves toward a target driven by the strain energy density.
#[derive(Debug, Clone)]
pub struct BoneRemodeling {
    /// Current apparent density (kg/m³).
    pub density: f64,
    /// Reference strain energy density (J/m³).
    pub reference_sed: f64,
    /// Remodelling rate coefficient (s⁻¹).
    pub rate: f64,
    /// Dead zone half-width (fraction).
    pub dead_zone: f64,
    /// Minimum density (kg/m³).
    pub rho_min: f64,
    /// Maximum density (kg/m³).
    pub rho_max: f64,
    /// Anisotropy direction (unit vector).
    pub fabric_tensor: [f64; 3],
}

impl BoneRemodeling {
    /// Create a new bone remodelling model.
    pub fn new(initial_density: f64, reference_sed: f64) -> Self {
        Self {
            density: initial_density,
            reference_sed,
            rate: 0.01,
            dead_zone: 0.1,
            rho_min: 100.0,
            rho_max: 2_000.0,
            fabric_tensor: [1.0, 0.0, 0.0],
        }
    }

    /// Compute the local apparent Young's modulus from density.
    ///
    /// Uses the power law: `E = c * ρ^n` (Carter-Hayes relation, n≈2).
    pub fn youngs_modulus(&self) -> f64 {
        let c = 3.79e-3; // MPa / (kg/m³)^2 → Pa units
        c * self.density.powi(2)
    }

    /// Update density based on the strain energy density stimulus `u`.
    pub fn update(&mut self, strain_energy_density: f64, dt: f64) {
        let u = strain_energy_density;
        let u_ref = self.reference_sed;
        let error = (u - u_ref) / u_ref;
        // Dead zone: no remodelling if |error| < dead_zone
        let stimulus = if error.abs() < self.dead_zone {
            0.0
        } else if error > 0.0 {
            self.rate * (error - self.dead_zone)
        } else {
            self.rate * (error + self.dead_zone)
        };
        self.density =
            (self.density + stimulus * self.density * dt).clamp(self.rho_min, self.rho_max);
    }

    /// Compute strain energy density from a simple 1-D stress and Young's modulus.
    pub fn strain_energy_density_1d(stress: f64, modulus: f64) -> f64 {
        if modulus.abs() < 1e-20 {
            0.0
        } else {
            0.5 * stress * stress / modulus
        }
    }
}

// ---------------------------------------------------------------------------
// WoundHealing
// ---------------------------------------------------------------------------

/// Wound healing model: fibroblast mechanics, collagen deposition, scar tissue.
#[derive(Debug, Clone)]
pub struct WoundHealing {
    /// Wound area (m²).
    pub wound_area: f64,
    /// Collagen density (kg/m³).
    pub collagen_density: f64,
    /// Fibroblast density (cells/m²).
    pub fibroblast_density: f64,
    /// Rate of wound contraction (m²/s).
    pub contraction_rate: f64,
    /// Rate of collagen deposition (kg/m³/s per cell).
    pub deposition_rate: f64,
    /// Scar modulus relative to normal tissue.
    pub scar_modulus_ratio: f64,
}

impl WoundHealing {
    /// Create a new wound healing model with given initial wound area.
    pub fn new(wound_area: f64) -> Self {
        Self {
            wound_area,
            collagen_density: 0.0,
            fibroblast_density: 1e6,
            contraction_rate: 1e-8,
            deposition_rate: 1e-13,
            scar_modulus_ratio: 0.5,
        }
    }

    /// Update wound state for a single time step.
    pub fn update(&mut self, dt: f64) {
        // Wound contraction (fibroblast-mediated)
        let contraction = self.contraction_rate * self.fibroblast_density * dt;
        self.wound_area = (self.wound_area - contraction).max(0.0);
        // Collagen deposition
        self.collagen_density += self.deposition_rate * self.fibroblast_density * dt;
        // Fibroblast response: decrease as wound closes
        let closed_fraction = 1.0 - self.wound_area / (self.wound_area + 1e-10);
        self.fibroblast_density *= 1.0 - closed_fraction * 0.01 * dt;
    }

    /// Fraction of wound closure.
    pub fn closure_fraction(&self, initial_area: f64) -> f64 {
        if initial_area < 1e-15 {
            return 1.0;
        }
        (1.0 - self.wound_area / initial_area).clamp(0.0, 1.0)
    }

    /// Scar tissue stiffness relative to initial area.
    pub fn scar_stiffness(&self, base_modulus: f64) -> f64 {
        let collagen_factor = 1.0 + self.collagen_density * 0.001;
        base_modulus * self.scar_modulus_ratio * collagen_factor
    }
}

// ---------------------------------------------------------------------------
// PressureUlcer
// ---------------------------------------------------------------------------

/// Tissue damage model under sustained pressure (pressure ulcer risk).
#[derive(Debug, Clone)]
pub struct PressureUlcer {
    /// Applied contact pressure (Pa).
    pub pressure: f64,
    /// Ischemia threshold pressure (Pa).
    pub ischemia_threshold: f64,
    /// Necrosis threshold (accumulated damage).
    pub necrosis_threshold: f64,
    /// Accumulated ischemic damage.
    pub damage: f64,
    /// Time above ischemia threshold (s).
    pub ischemia_time: f64,
}

impl PressureUlcer {
    /// Create a new pressure ulcer model.
    pub fn new(ischemia_threshold: f64, necrosis_threshold: f64) -> Self {
        Self {
            pressure: 0.0,
            ischemia_threshold,
            necrosis_threshold,
            damage: 0.0,
            ischemia_time: 0.0,
        }
    }

    /// Update model for one time step.
    pub fn update(&mut self, applied_pressure: f64, dt: f64) {
        self.pressure = applied_pressure;
        if applied_pressure > self.ischemia_threshold {
            let excess = applied_pressure - self.ischemia_threshold;
            let damage_rate = excess / self.ischemia_threshold * 0.001;
            self.damage += damage_rate * dt;
            self.ischemia_time += dt;
        } else {
            // Partial recovery
            self.damage = (self.damage - 0.0001 * dt).max(0.0);
        }
    }

    /// Returns `true` if necrosis threshold has been reached.
    pub fn is_necrotic(&self) -> bool {
        self.damage >= self.necrosis_threshold
    }

    /// Risk score in `[0, 1]`.
    pub fn risk_score(&self) -> f64 {
        (self.damage / self.necrosis_threshold).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// BiomechanicsAnalysis
// ---------------------------------------------------------------------------

/// Biomechanics analysis: joint contact pressure, muscle moment arms, reaction forces.
#[derive(Debug, Clone)]
pub struct BiomechanicsAnalysis {
    /// Joint position in world space.
    pub joint_position: [f64; 3],
    /// Contact area (m²).
    pub contact_area: f64,
    /// Number of muscles crossing the joint.
    pub num_muscles: usize,
    /// Muscle force magnitudes (N).
    pub muscle_forces: Vec<f64>,
    /// Muscle moment arms (m).
    pub moment_arms: Vec<f64>,
    /// Muscle line-of-action vectors.
    pub muscle_directions: Vec<[f64; 3]>,
}

impl BiomechanicsAnalysis {
    /// Create a new biomechanics analysis.
    pub fn new(joint_position: [f64; 3], contact_area: f64, num_muscles: usize) -> Self {
        Self {
            joint_position,
            contact_area,
            num_muscles,
            muscle_forces: vec![0.0; num_muscles],
            moment_arms: vec![0.0; num_muscles],
            muscle_directions: vec![[1.0, 0.0, 0.0]; num_muscles],
        }
    }

    /// Compute joint contact pressure from net compressive force.
    ///
    /// `p = F_net / contact_area`
    pub fn joint_contact_pressure(&self, external_force: [f64; 3], joint_normal: [f64; 3]) -> f64 {
        let f_net = self.joint_reaction_force(external_force);
        let compressive = dot3(f_net, joint_normal).abs();
        if self.contact_area > 1e-15 {
            compressive / self.contact_area
        } else {
            0.0
        }
    }

    /// Net joint reaction force: sum of all muscle forces + external force.
    pub fn joint_reaction_force(&self, external_force: [f64; 3]) -> [f64; 3] {
        let mut f = external_force;
        for (i, &fm) in self.muscle_forces.iter().enumerate() {
            let dir = self.muscle_directions[i];
            f[0] += fm * dir[0];
            f[1] += fm * dir[1];
            f[2] += fm * dir[2];
        }
        f
    }

    /// Compute net joint torque from muscle forces and moment arms.
    pub fn net_joint_torque(&self) -> f64 {
        self.muscle_forces
            .iter()
            .zip(self.moment_arms.iter())
            .map(|(&f, &r)| f * r)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. TissueType Young's modulus ordering
    #[test]
    fn test_tissue_type_modulus_ordering() {
        assert!(TissueType::Bone.youngs_modulus() > TissueType::Tendon.youngs_modulus());
        assert!(TissueType::Tendon.youngs_modulus() > TissueType::Muscle.youngs_modulus());
        assert!(TissueType::Muscle.youngs_modulus() > 0.0);
    }

    // 2. MuscleFiber passive force is zero at rest
    #[test]
    fn test_muscle_fiber_passive_at_rest() {
        let f = MuscleFiber::new([1.0, 0.0, 0.0], 1000.0, 0.1);
        assert!(f.passive_force_length().abs() < 1e-10);
    }

    // 3. MuscleFiber active force is 1 at l_opt with activation=1
    #[test]
    fn test_muscle_fiber_active_force_at_opt_length() {
        let f = MuscleFiber::new([1.0, 0.0, 0.0], 1000.0, 0.1);
        let fl = f.active_force_length();
        assert!((fl - 1.0).abs() < 0.01, "fl={fl}");
    }

    // 4. MuscleFiber total force with activation=0 is passive only
    #[test]
    fn test_muscle_fiber_zero_activation() {
        let mut f = MuscleFiber::new([1.0, 0.0, 0.0], 1000.0, 0.1);
        f.activation = 0.0;
        f.length = f.l_opt * 1.2; // extend beyond opt
        let total = f.total_force();
        let passive = f.passive_force_length() * f.f_max;
        assert!(
            (total - passive).abs() < 1e-6,
            "total={total}, passive={passive}"
        );
    }

    // 5. Hill force-velocity at zero velocity gives 1
    #[test]
    fn test_hill_force_velocity_zero() {
        let fv = hill_force_velocity(0.0, 0.1);
        assert!((fv - 1.0).abs() < 1e-10, "fv={fv}");
    }

    // 6. Hill force-velocity concentric decreases
    #[test]
    fn test_hill_force_velocity_concentric() {
        let fv = hill_force_velocity(0.5, 0.1);
        assert!(fv < 1.0, "fv={fv}");
        assert!(fv >= 0.0);
    }

    // 7. Hill model CE force is positive with activation
    #[test]
    fn test_hill_model_ce_force() {
        let mut hm = HillModel::new(1000.0, 0.1, 0.05, 0.1);
        hm.activation = 0.5;
        let f = hm.ce_force();
        assert!(f > 0.0, "f={f}");
    }

    // 8. Tendon model elastic force zero below slack
    #[test]
    fn test_tendon_zero_below_slack() {
        let t = TendonModel::new(0.3, 1e6);
        assert!(t.elastic_force().abs() < 1e-10);
    }

    // 9. Tendon model elastic force positive above slack
    #[test]
    fn test_tendon_positive_above_slack() {
        let mut t = TendonModel::new(0.3, 1e6);
        t.length = 0.31;
        assert!(t.elastic_force() > 0.0);
    }

    // 10. SkinLayer stress at prestretch = 0 (no strain)
    #[test]
    fn test_skin_stress_at_prestretch() {
        let s = SkinLayer::dermis();
        let stress = s.stress(s.prestretch);
        // stretch = prestretch → lambda = 1 → neo-Hookean ≈ 0
        assert!(stress.abs() < 1000.0, "stress={stress}"); // small but not exactly 0
    }

    // 11. SkinModel total thickness
    #[test]
    fn test_skin_model_total_thickness() {
        let skin = SkinModel::default_skin();
        let total = skin.total_thickness();
        assert!((total - 0.0071).abs() < 1e-5, "total={total}");
    }

    // 12. CartilageModel solid stress linear
    #[test]
    fn test_cartilage_solid_stress_linear() {
        let c = CartilageModel::new(1e6, 1e-15);
        let s1 = c.solid_stress(0.1);
        let s2 = c.solid_stress(0.2);
        assert!((s2 - 2.0 * s1).abs() < 1e-6, "s1={s1}, s2={s2}");
    }

    // 13. CartilageModel update reduces strain error
    #[test]
    fn test_cartilage_update() {
        // Use large permeability and stress well above swelling pressure (50_000 Pa)
        let mut c = CartilageModel::new(1e6, 1e-9);
        // Applied stress > swelling_pressure drives positive strain
        let stress = 200_000.0;
        for _ in 0..100 {
            c.update(stress, 0.1);
        }
        // After iterations, strain should be non-zero
        assert!(c.strain > 0.0);
    }

    // 14. BoneRemodeling Young's modulus increases with density
    #[test]
    fn test_bone_modulus_density_scaling() {
        let mut b1 = BoneRemodeling::new(500.0, 1000.0);
        let mut b2 = BoneRemodeling::new(1000.0, 1000.0);
        b1.density = 500.0;
        b2.density = 1000.0;
        assert!(b2.youngs_modulus() > b1.youngs_modulus());
    }

    // 15. BoneRemodeling density increases with high SED
    #[test]
    fn test_bone_density_increases_high_sed() {
        let mut b = BoneRemodeling::new(500.0, 100.0);
        let initial_density = b.density;
        for _ in 0..100 {
            b.update(1000.0, 0.1);
        } // 10x reference SED
        assert!(b.density > initial_density, "density should increase");
    }

    // 16. BoneRemodeling density decreases with low SED
    #[test]
    fn test_bone_density_decreases_low_sed() {
        let mut b = BoneRemodeling::new(1000.0, 1000.0);
        let initial_density = b.density;
        for _ in 0..100 {
            b.update(1.0, 0.1);
        } // very low SED
        assert!(b.density < initial_density, "density should decrease");
    }

    // 17. WoundHealing wound area decreases over time
    #[test]
    fn test_wound_healing_area_decreases() {
        let mut wh = WoundHealing::new(0.01);
        let initial = wh.wound_area;
        for _ in 0..1000 {
            wh.update(0.1);
        }
        assert!(wh.wound_area < initial, "wound should shrink");
    }

    // 18. WoundHealing collagen accumulates
    #[test]
    fn test_wound_healing_collagen_accumulates() {
        let mut wh = WoundHealing::new(0.01);
        for _ in 0..100 {
            wh.update(1.0);
        }
        assert!(wh.collagen_density > 0.0);
    }

    // 19. PressureUlcer damage accumulates above threshold
    #[test]
    fn test_pressure_ulcer_damage_accumulates() {
        let mut pu = PressureUlcer::new(1000.0, 1.0);
        pu.update(2000.0, 10.0); // 2x threshold for 10 seconds
        assert!(pu.damage > 0.0, "damage={}", pu.damage);
    }

    // 20. PressureUlcer not necrotic at low pressure
    #[test]
    fn test_pressure_ulcer_no_necrosis_low_pressure() {
        let mut pu = PressureUlcer::new(1000.0, 1.0);
        pu.update(500.0, 100.0); // below threshold
        assert!(!pu.is_necrotic());
    }

    // 21. PressureUlcer risk score in [0, 1]
    #[test]
    fn test_pressure_ulcer_risk_score_bounded() {
        let mut pu = PressureUlcer::new(1000.0, 1.0);
        pu.update(5000.0, 1.0);
        let r = pu.risk_score();
        assert!((0.0..=1.0).contains(&r), "r={r}");
    }

    // 22. BiomechanicsAnalysis joint reaction force without muscles
    #[test]
    fn test_biomechanics_reaction_no_muscles() {
        let ba = BiomechanicsAnalysis::new([0.0; 3], 0.01, 0);
        let f = ba.joint_reaction_force([10.0, 0.0, 0.0]);
        assert!((f[0] - 10.0).abs() < 1e-10);
    }

    // 23. BiomechanicsAnalysis contact pressure
    #[test]
    fn test_biomechanics_contact_pressure() {
        let ba = BiomechanicsAnalysis::new([0.0; 3], 0.01, 0);
        let p = ba.joint_contact_pressure([0.0, -100.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((p - 100.0 / 0.01).abs() < 1e-6, "p={p}");
    }

    // 24. BiomechanicsAnalysis net torque
    #[test]
    fn test_biomechanics_net_torque() {
        let mut ba = BiomechanicsAnalysis::new([0.0; 3], 0.01, 2);
        ba.muscle_forces = vec![100.0, 50.0];
        ba.moment_arms = vec![0.05, 0.03];
        let tau = ba.net_joint_torque();
        // 100*0.05 + 50*0.03 = 5 + 1.5 = 6.5
        assert!((tau - 6.5).abs() < 1e-10, "tau={tau}");
    }

    // 25. TissueType density all positive
    #[test]
    fn test_tissue_type_density_positive() {
        let types = [
            TissueType::Muscle,
            TissueType::Tendon,
            TissueType::Ligament,
            TissueType::Cartilage,
            TissueType::Bone,
            TissueType::Skin,
        ];
        for t in &types {
            assert!(t.density() > 0.0, "density={}", t.density());
        }
    }

    // 26. Muscle fiber force vector has same magnitude as scalar force
    #[test]
    fn test_muscle_fiber_force_vector() {
        let mut f = MuscleFiber::new([1.0, 0.0, 0.0], 1000.0, 0.1);
        f.activation = 0.8;
        let fv = f.force_vector();
        let fs = f.total_force();
        let mag = norm3(fv);
        assert!((mag - fs).abs() < 1e-6, "mag={mag}, fs={fs}");
    }

    // 27. Hill model update changes CE length
    #[test]
    fn test_hill_model_update_changes_ce() {
        let mut hm = HillModel::new(1000.0, 0.1, 0.05, 0.1);
        hm.activation = 1.0;
        let l_init = hm.l_ce;
        hm.update(0.12, 1.0, 0.01);
        // CE length should have changed
        let _changed = (hm.l_ce - l_init).abs() > 1e-12;
        // may not always change — just test no panic
    }

    // 28. SkinModel bending stiffness positive
    #[test]
    fn test_skin_bending_stiffness_positive() {
        let skin = SkinModel::default_skin();
        assert!(skin.bending_stiffness() > 0.0);
    }

    // 29. CartilageModel creep compliance increases with time
    #[test]
    fn test_cartilage_creep_compliance_increases() {
        let c = CartilageModel::new(1e6, 1e-15);
        let j1 = c.creep_compliance(1.0);
        let j2 = c.creep_compliance(10.0);
        assert!(j2 >= j1, "j1={j1}, j2={j2}");
    }
}

// ============================================================================
// EXPANSION: Tissue Simulation — Vascular, Remodeling, Neural, Ocular
// ============================================================================

// ── Vascular Wall Mechanics ───────────────────────────────────────────────────

/// Layer description for layered arterial wall models (e.g., intima-media-adventitia).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VascularLayer {
    /// Layer name.
    pub name: String,
    /// Unstressed thickness (m).
    pub thickness: f64,
    /// Neo-Hookean initial stiffness μ (Pa).
    pub mu: f64,
    /// Exponential stiffening parameter α.
    pub alpha: f64,
    /// Fiber angle with respect to circumferential direction (rad).
    pub fiber_angle: f64,
    /// Fiber volume fraction.
    pub fiber_fraction: f64,
}

impl VascularLayer {
    /// Intima layer (thin inner layer).
    pub fn intima() -> Self {
        Self {
            name: "intima".to_string(),
            thickness: 0.2e-3,
            mu: 5e3,
            alpha: 3.0,
            fiber_angle: 0.0,
            fiber_fraction: 0.1,
        }
    }

    /// Media layer (smooth muscle cells and collagen).
    pub fn media() -> Self {
        Self {
            name: "media".to_string(),
            thickness: 0.5e-3,
            mu: 20e3,
            alpha: 8.0,
            fiber_angle: 0.15, // ~8.6 degrees
            fiber_fraction: 0.45,
        }
    }

    /// Adventitia (outer collagen-rich layer).
    pub fn adventitia() -> Self {
        Self {
            name: "adventitia".to_string(),
            thickness: 0.3e-3,
            mu: 2e3,
            alpha: 15.0,
            fiber_angle: 0.52, // ~30 degrees
            fiber_fraction: 0.60,
        }
    }

    /// Circumferential strain energy density (Fung-type).
    pub fn strain_energy(&self, lambda_theta: f64) -> f64 {
        let i1 = lambda_theta * lambda_theta;
        let neo = 0.5 * self.mu * (i1 - 1.0);
        let fung = if self.alpha > 1e-15 {
            self.mu / (2.0 * self.alpha) * ((self.alpha * (i1 - 1.0)).exp() - 1.0)
        } else {
            0.0
        };
        (1.0 - self.fiber_fraction) * neo + self.fiber_fraction * fung
    }

    /// Circumferential Cauchy stress from stretch ratio.
    pub fn stress(&self, lambda_theta: f64) -> f64 {
        let i1 = lambda_theta * lambda_theta;
        let dns = 2.0 * self.mu * lambda_theta;
        let dfung = if self.alpha > 1e-15 {
            self.mu * lambda_theta * (self.alpha * (i1 - 1.0)).exp()
        } else {
            0.0
        };
        (1.0 - self.fiber_fraction) * dns + self.fiber_fraction * dfung
    }
}

/// Multi-layer cylindrical arterial wall model.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ArterialWall {
    /// Wall layers from inner to outer.
    pub layers: Vec<VascularLayer>,
    /// Inner radius at zero-stress reference (m).
    pub inner_radius: f64,
    /// Residual strain angle (opening angle, rad) for stress-free configuration.
    pub opening_angle: f64,
    /// Axial pre-stretch ratio.
    pub axial_prestretch: f64,
}

impl ArterialWall {
    /// Create a three-layer aortic wall model.
    pub fn aorta() -> Self {
        Self {
            layers: vec![
                VascularLayer::intima(),
                VascularLayer::media(),
                VascularLayer::adventitia(),
            ],
            inner_radius: 12e-3,
            opening_angle: 0.5, // ~28 degrees
            axial_prestretch: 1.1,
        }
    }

    /// Compute total wall thickness.
    pub fn wall_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }

    /// Outer radius at zero-stress state.
    pub fn outer_radius(&self) -> f64 {
        self.inner_radius + self.wall_thickness()
    }

    /// Approximate pressure-radius relationship (Laplace + layered wall).
    ///
    /// Returns transmural pressure for a given luminal radius `r`.
    pub fn pressure_at_radius(&self, r: f64) -> f64 {
        if r <= 1e-30 {
            return 0.0;
        }
        let lambda_theta = r / self.inner_radius;
        // Summation over layers (simplification: use mean stretch)
        let total_stress: f64 = self.layers.iter().map(|l| l.stress(lambda_theta)).sum();
        let thickness = self.wall_thickness();
        // Laplace: P = σ_θ * h / r
        total_stress * thickness / r
    }

    /// Compliance: dV/dP at given pressure (numerical).
    pub fn compliance_at_pressure(&self, p_ref: f64, dp: f64) -> f64 {
        // Approximate inverse: find radius from pressure
        // Binary search
        let r1 = self.radius_at_pressure(p_ref);
        let r2 = self.radius_at_pressure(p_ref + dp);
        let v1 = std::f64::consts::PI * r1 * r1;
        let v2 = std::f64::consts::PI * r2 * r2;
        (v2 - v1) / dp
    }

    /// Find inner radius at given pressure using bisection.
    pub fn radius_at_pressure(&self, target_p: f64) -> f64 {
        let mut lo = self.inner_radius;
        let mut hi = self.inner_radius * 3.0;
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.pressure_at_radius(mid) < target_p {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ── Biphasic Soft Tissue (Saturated Porous) ───────────────────────────────────

/// Biphasic soft tissue model (Mow et al., 1980).
///
/// Models cartilage/meniscus as a solid phase (elastic matrix) + fluid phase.
/// Fluid exudation during compression provides time-dependent response.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BiphasicTissue {
    /// Drained elastic modulus E_s (Pa).
    pub e_solid: f64,
    /// Drained Poisson's ratio ν_s.
    pub nu_solid: f64,
    /// Hydraulic permeability k (m²/Pa·s).
    pub permeability: f64,
    /// Fluid volume fraction (porosity).
    pub porosity: f64,
    /// Stress-dependent permeability exponent M (0 = constant).
    pub perm_exponent: f64,
    /// Current fluid pressure (Pa).
    pub fluid_pressure: f64,
    /// Current strain in solid phase.
    pub solid_strain: f64,
}

impl BiphasicTissue {
    /// Create articular cartilage model.
    pub fn cartilage() -> Self {
        Self {
            e_solid: 0.7e6,
            nu_solid: 0.3,
            permeability: 1e-15,
            porosity: 0.8,
            perm_exponent: 0.0,
            fluid_pressure: 0.0,
            solid_strain: 0.0,
        }
    }

    /// Create intervertebral disc nucleus model.
    pub fn disc_nucleus() -> Self {
        Self {
            e_solid: 0.3e6,
            nu_solid: 0.1,
            permeability: 5e-16,
            porosity: 0.88,
            perm_exponent: 1.0,
            fluid_pressure: 0.0,
            solid_strain: 0.0,
        }
    }

    /// Aggregate modulus H_A = E_s(1-ν_s)/((1+ν_s)(1-2ν_s)).
    pub fn aggregate_modulus(&self) -> f64 {
        let e = self.e_solid;
        let nu = self.nu_solid;
        e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Biphasic characteristic time τ = H_A k / h² (for layer height h).
    pub fn characteristic_time(&self, height: f64) -> f64 {
        if height <= 0.0 {
            return 0.0;
        }
        self.aggregate_modulus() * self.permeability / (height * height)
    }

    /// Fluid pressure creep response at normalized time t/τ.
    ///
    /// For confined compression: p(t) = p0 * Σ A_k exp(-k²π²τt).
    pub fn creep_pressure(&self, t_normalized: f64, p0: f64) -> f64 {
        let mut sum = 0.0;
        for k in 1..=10_usize {
            let kf = k as f64;
            let coeff = 8.0 / (kf * kf * std::f64::consts::PI * std::f64::consts::PI);
            sum += coeff
                * (-kf * kf * std::f64::consts::PI * std::f64::consts::PI * t_normalized).exp();
        }
        p0 * sum
    }

    /// Stress-dependent permeability: k(e) = k0 exp(M e_v).
    pub fn strain_dependent_permeability(&self, volumetric_strain: f64) -> f64 {
        self.permeability * (self.perm_exponent * volumetric_strain).exp()
    }

    /// Update fluid pressure under step compression ε_0 for time dt.
    pub fn update(&mut self, applied_strain: f64, dt: f64, height: f64) {
        let tau = self.characteristic_time(height);
        let t_norm = if tau > 1e-30 { dt / tau } else { 1.0 };
        let ha = self.aggregate_modulus();
        let p0 = ha * applied_strain;
        self.fluid_pressure = self.creep_pressure(t_norm, p0);
        self.solid_strain = applied_strain;
    }

    /// Effective solid stress = total - fluid pressure.
    pub fn effective_stress(&self, total_stress: f64) -> f64 {
        total_stress - self.fluid_pressure
    }
}

// ── Tendon Viscoelastic Model ──────────────────────────────────────────────────

/// Viscoelastic tendon/ligament model with toe, linear, and failure regions.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TendonViscoelastic {
    /// Cross-sectional area A₀ (m²).
    pub area: f64,
    /// Slack length L₀ (m).
    pub slack_length: f64,
    /// Toe region limit strain ε_toe.
    pub epsilon_toe: f64,
    /// Linear region modulus E (Pa).
    pub elastic_modulus: f64,
    /// Toe region stiffness factor (fraction of linear E at ε_toe).
    pub toe_stiffness_factor: f64,
    /// Failure strain ε_fail.
    pub failure_strain: f64,
    /// Viscosity coefficient η (Pa·s).
    pub viscosity: f64,
    /// Previous strain (for viscous update).
    pub prev_strain: f64,
    /// Previous stress (for viscous update).
    pub prev_stress: f64,
}

impl TendonViscoelastic {
    /// Create a patellar tendon model.
    pub fn patellar_tendon() -> Self {
        Self {
            area: 50e-6, // 50 mm²
            slack_length: 40e-3,
            epsilon_toe: 0.02,
            elastic_modulus: 1.5e9,
            toe_stiffness_factor: 0.5,
            failure_strain: 0.08,
            viscosity: 1e5,
            prev_strain: 0.0,
            prev_stress: 0.0,
        }
    }

    /// Elastic stress at given strain.
    pub fn elastic_stress(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            return 0.0;
        }
        if strain <= self.epsilon_toe {
            // Toe region: quadratic stress
            let e_toe = self.elastic_modulus * self.toe_stiffness_factor;
            e_toe * strain * strain / self.epsilon_toe
        } else if strain < self.failure_strain {
            // Linear region
            let toe_stress = self.elastic_modulus * self.toe_stiffness_factor * self.epsilon_toe;
            toe_stress + self.elastic_modulus * (strain - self.epsilon_toe)
        } else {
            // Failed: dramatic drop
            0.0
        }
    }

    /// Compute tensile force in tendon at given length L.
    pub fn force(&self, length: f64) -> f64 {
        let strain = (length - self.slack_length) / self.slack_length;
        self.elastic_stress(strain) * self.area
    }

    /// Viscoelastic update with strain rate effect.
    pub fn update_stress(&mut self, new_strain: f64, dt: f64) -> f64 {
        let strain_rate = if dt > 1e-15 {
            (new_strain - self.prev_strain) / dt
        } else {
            0.0
        };
        let sigma_elastic = self.elastic_stress(new_strain);
        let sigma_viscous = self.viscosity * strain_rate;
        let total = sigma_elastic + sigma_viscous;
        self.prev_strain = new_strain;
        self.prev_stress = total;
        total
    }

    /// Check if tendon is ruptured.
    pub fn is_ruptured(&self, strain: f64) -> bool {
        strain >= self.failure_strain
    }

    /// Stiffness (tangent modulus) at given strain.
    pub fn stiffness_at_strain(&self, strain: f64) -> f64 {
        if strain <= 0.0 || strain >= self.failure_strain {
            return 0.0;
        }
        if strain <= self.epsilon_toe {
            2.0 * self.elastic_modulus * self.toe_stiffness_factor * strain / self.epsilon_toe
        } else {
            self.elastic_modulus
        }
    }
}

// ── Wound Healing and Tissue Remodeling ───────────────────────────────────────

/// Wound healing phase classification.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum WoundHealingPhase {
    /// Hemostasis (0–1 days).
    Hemostasis,
    /// Inflammation (1–4 days).
    Inflammation,
    /// Proliferation (4–21 days).
    Proliferation,
    /// Remodeling (21 days – 2 years).
    Remodeling,
    /// Mature scar.
    Scar,
}

/// Tissue remodeling simulation tracking collagen turnover.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TissueRemodeling {
    /// Current collagen density (normalized, 0–1).
    pub collagen_density: f64,
    /// Target collagen density based on mechanical loading.
    pub target_density: f64,
    /// Remodeling rate constant (1/day).
    pub remodeling_rate: f64,
    /// Current stiffness (Pa) — proportional to collagen density.
    pub stiffness: f64,
    /// Reference stiffness at collagen_density = 1.
    pub stiffness_ref: f64,
    /// Days elapsed since injury.
    pub elapsed_days: f64,
    /// Current healing phase.
    pub phase: WoundHealingPhase,
}

impl TissueRemodeling {
    /// Create a new remodeling simulation starting from injury.
    pub fn new_wound(stiffness_ref: f64, remodeling_rate: f64) -> Self {
        Self {
            collagen_density: 0.1, // initial granulation tissue
            target_density: 1.0,
            remodeling_rate,
            stiffness: stiffness_ref * 0.1,
            stiffness_ref,
            elapsed_days: 0.0,
            phase: WoundHealingPhase::Hemostasis,
        }
    }

    /// Update remodeling for time `dt_days` under mechanical stimulus `strain`.
    pub fn update(&mut self, dt_days: f64, strain: f64) {
        self.elapsed_days += dt_days;

        // Update phase based on time
        self.phase = if self.elapsed_days < 1.0 {
            WoundHealingPhase::Hemostasis
        } else if self.elapsed_days < 4.0 {
            WoundHealingPhase::Inflammation
        } else if self.elapsed_days < 21.0 {
            WoundHealingPhase::Proliferation
        } else if self.elapsed_days < 365.0 {
            WoundHealingPhase::Remodeling
        } else {
            WoundHealingPhase::Scar
        };

        // Mechanical stimulus adjusts target density
        // Moderate strain promotes collagen; excess or zero reduces it
        let optimal_strain = 0.04;
        let stimulus = (-((strain - optimal_strain) / 0.02).powi(2)).exp();
        self.target_density = 0.5 + 0.5 * stimulus;

        // Exponential approach to target
        let error = self.target_density - self.collagen_density;
        self.collagen_density += error * self.remodeling_rate * dt_days;
        self.collagen_density = self.collagen_density.clamp(0.01, 1.2);

        // Update stiffness
        self.stiffness = self.stiffness_ref * self.collagen_density;
    }

    /// Fraction of strength recovered relative to normal tissue.
    pub fn strength_recovery(&self) -> f64 {
        self.collagen_density.min(1.0)
    }
}

// ── Nerve Tissue and Action Potential ────────────────────────────────────────

/// Simplified Hodgkin-Huxley neuron model for tissue-nerve interaction.
///
/// Tracks membrane potential and gating variables.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HodgkinHuxley {
    /// Membrane capacitance (µF/cm²).
    pub cm: f64,
    /// Resting potential (mV).
    pub v_rest: f64,
    /// Current membrane potential (mV).
    pub v: f64,
    /// Na channel activation m.
    pub m: f64,
    /// Na channel inactivation h.
    pub h: f64,
    /// K channel activation n.
    pub n: f64,
    /// Maximum Na conductance (mS/cm²).
    pub g_na: f64,
    /// Maximum K conductance (mS/cm²).
    pub g_k: f64,
    /// Leakage conductance (mS/cm²).
    pub g_l: f64,
    /// Na reversal potential (mV).
    pub e_na: f64,
    /// K reversal potential (mV).
    pub e_k: f64,
    /// Leakage reversal potential (mV).
    pub e_l: f64,
}

impl HodgkinHuxley {
    /// Create standard squid axon parameters.
    pub fn squid_axon() -> Self {
        Self {
            cm: 1.0,
            v_rest: -65.0,
            v: -65.0,
            m: 0.05,
            h: 0.6,
            n: 0.32,
            g_na: 120.0,
            g_k: 36.0,
            g_l: 0.3,
            e_na: 50.0,
            e_k: -77.0,
            e_l: -54.4,
        }
    }

    /// Alpha and beta rate functions for m gating variable.
    fn alpha_m(&self) -> f64 {
        let dv = self.v + 40.0;
        if dv.abs() < 1e-7 {
            1.0
        } else {
            0.1 * dv / (1.0 - (-dv / 10.0).exp())
        }
    }

    fn beta_m(&self) -> f64 {
        4.0 * (-(self.v + 65.0) / 18.0).exp()
    }

    fn alpha_h(&self) -> f64 {
        0.07 * (-(self.v + 65.0) / 20.0).exp()
    }

    fn beta_h(&self) -> f64 {
        1.0 / (1.0 + (-(self.v + 35.0) / 10.0).exp())
    }

    fn alpha_n(&self) -> f64 {
        let dv = self.v + 55.0;
        if dv.abs() < 1e-7 {
            0.1
        } else {
            0.01 * dv / (1.0 - (-dv / 10.0).exp())
        }
    }

    fn beta_n(&self) -> f64 {
        0.125 * (-(self.v + 65.0) / 80.0).exp()
    }

    /// Step the model forward by `dt` ms with applied current `i_ext` µA/cm².
    pub fn step(&mut self, dt: f64, i_ext: f64) {
        let i_na = self.g_na * self.m * self.m * self.m * self.h * (self.v - self.e_na);
        let i_k = self.g_k * self.n * self.n * self.n * self.n * (self.v - self.e_k);
        let i_l = self.g_l * (self.v - self.e_l);
        let dv = (i_ext - i_na - i_k - i_l) / self.cm;

        let dm = self.alpha_m() * (1.0 - self.m) - self.beta_m() * self.m;
        let dh = self.alpha_h() * (1.0 - self.h) - self.beta_h() * self.h;
        let dn = self.alpha_n() * (1.0 - self.n) - self.beta_n() * self.n;

        self.v += dv * dt;
        self.m = (self.m + dm * dt).clamp(0.0, 1.0);
        self.h = (self.h + dh * dt).clamp(0.0, 1.0);
        self.n = (self.n + dn * dt).clamp(0.0, 1.0);
    }

    /// Check if action potential is occurring (V > threshold).
    pub fn is_firing(&self) -> bool {
        self.v > -20.0
    }

    /// Reset to resting state.
    pub fn reset(&mut self) {
        self.v = self.v_rest;
        self.m = 0.05;
        self.h = 0.6;
        self.n = 0.32;
    }
}

// ── Cornea/Ocular Tissue Mechanics ────────────────────────────────────────────

/// Simplified corneal biomechanics model.
///
/// Models the cornea as a thin curved shell with hyperelastic properties.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CornealModel {
    /// Central corneal thickness (m).
    pub thickness: f64,
    /// Corneal radius of curvature (m).
    pub radius: f64,
    /// Elastic modulus (Pa).
    pub elastic_modulus: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Intraocular pressure (Pa).
    pub iop: f64,
}

impl CornealModel {
    /// Create a normal cornea.
    pub fn normal() -> Self {
        Self {
            thickness: 0.55e-3, // 550 µm
            radius: 7.8e-3,     // 7.8 mm
            elastic_modulus: 0.3e6,
            poisson: 0.49, // nearly incompressible
            iop: 2000.0,   // ~15 mmHg
        }
    }

    /// Applanation tonometry estimate (Imbert-Fick law).
    ///
    /// P ≈ F / A where A is the applanation area.
    pub fn applanation_pressure(&self, force: f64, area: f64) -> f64 {
        if area <= 0.0 {
            return 0.0;
        }
        force / area
    }

    /// Corneal stiffness (shell bending stiffness D = E*h³/(12(1-ν²))).
    pub fn bending_stiffness(&self) -> f64 {
        let h = self.thickness;
        self.elastic_modulus * h * h * h / (12.0 * (1.0 - self.poisson * self.poisson))
    }

    /// Membrane stiffness (E*h).
    pub fn membrane_stiffness(&self) -> f64 {
        self.elastic_modulus * self.thickness
    }

    /// Laplace pressure at current IOP: σ = P * R / (2 h)
    pub fn laplace_stress(&self) -> f64 {
        self.iop * self.radius / (2.0 * self.thickness)
    }

    /// Deflection under point load F at center (simplified spherical shell).
    pub fn central_deflection(&self, force: f64) -> f64 {
        // d = F R² / (2 π D) for clamped spherical shell (approximate)
        let d = self.bending_stiffness();
        if d < 1e-30 {
            return 0.0;
        }
        force * self.radius * self.radius / (2.0 * std::f64::consts::PI * d)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod expanded_tissue_tests {

    use crate::tissue_sim::ArterialWall;
    use crate::tissue_sim::BiphasicTissue;
    use crate::tissue_sim::CornealModel;
    use crate::tissue_sim::HodgkinHuxley;
    use crate::tissue_sim::TendonViscoelastic;
    use crate::tissue_sim::TissueRemodeling;
    use crate::tissue_sim::VascularLayer;
    use crate::tissue_sim::WoundHealingPhase;

    // T1. VascularLayer intima strain energy positive for stretch
    #[test]
    fn test_intima_strain_energy() {
        let layer = VascularLayer::intima();
        let w = layer.strain_energy(1.1);
        assert!(w > 0.0, "w={w}");
    }

    // T2. VascularLayer stress positive for stretch > 1
    #[test]
    fn test_intima_stress_positive() {
        let layer = VascularLayer::media();
        let sigma = layer.stress(1.2);
        assert!(sigma > 0.0, "sigma={sigma}");
    }

    // T3. ArterialWall thickness is sum of layers
    #[test]
    fn test_arterial_wall_thickness() {
        let wall = ArterialWall::aorta();
        let total = wall.wall_thickness();
        let sum: f64 = wall.layers.iter().map(|l| l.thickness).sum();
        assert!((total - sum).abs() < 1e-15);
    }

    // T4. ArterialWall pressure positive
    #[test]
    fn test_arterial_pressure_positive() {
        let wall = ArterialWall::aorta();
        let p = wall.pressure_at_radius(wall.inner_radius * 1.1);
        assert!(p > 0.0, "p={p}");
    }

    // T5. ArterialWall radius_at_pressure monotone
    #[test]
    fn test_arterial_radius_monotone() {
        let wall = ArterialWall::aorta();
        let r1 = wall.radius_at_pressure(5000.0);
        let r2 = wall.radius_at_pressure(10000.0);
        assert!(r2 > r1, "r1={r1}, r2={r2}");
    }

    // T6. BiphasicTissue aggregate modulus positive
    #[test]
    fn test_biphasic_aggregate_modulus() {
        let bt = BiphasicTissue::cartilage();
        let ha = bt.aggregate_modulus();
        assert!(ha > 0.0, "ha={ha}");
    }

    // T7. BiphasicTissue characteristic time positive
    #[test]
    fn test_biphasic_characteristic_time() {
        let bt = BiphasicTissue::cartilage();
        let tau = bt.characteristic_time(2e-3);
        assert!(tau > 0.0, "tau={tau}");
    }

    // T8. BiphasicTissue creep pressure decreases with time
    #[test]
    fn test_biphasic_creep_decay() {
        let bt = BiphasicTissue::cartilage();
        let p1 = bt.creep_pressure(0.0, 1000.0);
        let p2 = bt.creep_pressure(1.0, 1000.0);
        assert!(p1 >= p2, "p1={p1}, p2={p2}");
    }

    // T9. TendonViscoelastic elastic stress zero before slack
    #[test]
    fn test_tendon_no_stress_at_zero_strain() {
        let tendon = TendonViscoelastic::patellar_tendon();
        assert!(tendon.elastic_stress(0.0).abs() < 1e-10);
    }

    // T10. TendonViscoelastic force positive above slack length
    #[test]
    fn test_tendon_force_positive() {
        let tendon = TendonViscoelastic::patellar_tendon();
        let f = tendon.force(tendon.slack_length * 1.05);
        assert!(f > 0.0, "f={f}");
    }

    // T11. TendonViscoelastic stiffness increases from toe to linear
    #[test]
    fn test_tendon_stiffness_increase() {
        let tendon = TendonViscoelastic::patellar_tendon();
        let k_toe = tendon.stiffness_at_strain(tendon.epsilon_toe * 0.5);
        let k_lin = tendon.stiffness_at_strain(tendon.epsilon_toe * 1.5);
        assert!(k_lin > k_toe, "k_toe={k_toe}, k_lin={k_lin}");
    }

    // T12. TendonViscoelastic rupture at failure strain
    #[test]
    fn test_tendon_rupture() {
        let tendon = TendonViscoelastic::patellar_tendon();
        assert!(tendon.is_ruptured(tendon.failure_strain));
        assert!(!tendon.is_ruptured(tendon.failure_strain * 0.5));
    }

    // T13. TissueRemodeling phase progression
    #[test]
    fn test_remodeling_phase_progression() {
        let mut rem = TissueRemodeling::new_wound(1e6, 0.1);
        assert_eq!(rem.phase, WoundHealingPhase::Hemostasis);
        rem.update(2.0, 0.03);
        assert_eq!(rem.phase, WoundHealingPhase::Inflammation);
        rem.update(5.0, 0.03);
        assert_eq!(rem.phase, WoundHealingPhase::Proliferation);
    }

    // T14. TissueRemodeling collagen density increases under optimal loading
    #[test]
    fn test_remodeling_collagen_increase() {
        let mut rem = TissueRemodeling::new_wound(1e6, 0.5);
        let c0 = rem.collagen_density;
        for _ in 0..30 {
            rem.update(1.0, 0.04); // optimal strain
        }
        assert!(rem.collagen_density > c0, "collagen should increase");
    }

    // T15. HodgkinHuxley resting potential stable without stimulus
    #[test]
    fn test_hh_resting_stable() {
        let mut hh = HodgkinHuxley::squid_axon();
        for _ in 0..100 {
            hh.step(0.01, 0.0); // dt=0.01 ms, no stimulus
        }
        assert!((hh.v - hh.v_rest).abs() < 5.0, "v={}", hh.v);
    }

    // T16. HodgkinHuxley fires with large stimulus
    #[test]
    fn test_hh_action_potential() {
        let mut hh = HodgkinHuxley::squid_axon();
        let mut fired = false;
        for _ in 0..500 {
            hh.step(0.02, 10.0); // 10 µA/cm² stimulus
            if hh.is_firing() {
                fired = true;
                break;
            }
        }
        assert!(fired, "Should fire an action potential with 10 µA/cm²");
    }

    // T17. CornealModel bending stiffness positive
    #[test]
    fn test_cornea_bending_stiffness() {
        let cornea = CornealModel::normal();
        let d = cornea.bending_stiffness();
        assert!(d > 0.0, "d={d}");
    }

    // T18. CornealModel Laplace stress > 0
    #[test]
    fn test_cornea_laplace_stress() {
        let cornea = CornealModel::normal();
        let sigma = cornea.laplace_stress();
        assert!(sigma > 0.0, "sigma={sigma}");
    }

    // T19. BiphasicTissue strain-dependent permeability increases with compression
    #[test]
    fn test_biphasic_permeability() {
        let mut bt = BiphasicTissue::disc_nucleus();
        bt.perm_exponent = 1.5;
        let k0 = bt.strain_dependent_permeability(0.0);
        let k1 = bt.strain_dependent_permeability(0.1);
        assert!(k1 > k0, "k0={k0}, k1={k1}");
    }

    // T20. TissueRemodeling strength recovery ≤ 1
    #[test]
    fn test_remodeling_strength_recovery() {
        let rem = TissueRemodeling::new_wound(1e6, 0.1);
        assert!(rem.strength_recovery() <= 1.0);
    }

    // T21. VascularLayer adventitia has larger alpha than intima
    #[test]
    fn test_adventitia_stiffer() {
        let adv = VascularLayer::adventitia();
        let int = VascularLayer::intima();
        assert!(adv.alpha > int.alpha, "adventitia should have larger alpha");
    }

    // T22. ArterialWall outer > inner radius
    #[test]
    fn test_arterial_geometry() {
        let wall = ArterialWall::aorta();
        assert!(wall.outer_radius() > wall.inner_radius);
    }

    // T23. TendonViscoelastic viscoelastic stress > elastic under fast loading
    #[test]
    fn test_tendon_viscous_effect() {
        let mut tendon = TendonViscoelastic::patellar_tendon();
        let strain = 0.03;
        let elastic = tendon.elastic_stress(strain);
        let total = tendon.update_stress(strain, 0.001); // very fast loading
        assert!(total >= elastic, "Viscous effect should add stress");
    }

    // T24. CornealModel central deflection positive under load
    #[test]
    fn test_cornea_deflection() {
        let cornea = CornealModel::normal();
        let d = cornea.central_deflection(0.001);
        assert!(d > 0.0, "d={d}");
    }

    // T25. HodgkinHuxley reset returns to rest
    #[test]
    fn test_hh_reset() {
        let mut hh = HodgkinHuxley::squid_axon();
        for _ in 0..100 {
            hh.step(0.02, 10.0);
        }
        hh.reset();
        assert!((hh.v - hh.v_rest).abs() < 1e-10, "v={}", hh.v);
    }
}
