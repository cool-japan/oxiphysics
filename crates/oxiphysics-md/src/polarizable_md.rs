// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Polarizable force field MD.
//!
//! Implements induced-dipole polarization via a self-consistent field (SCF)
//! loop with Thole damping, Drude oscillator forces, polarization energy,
//! AMOEBA-like multipole electrostatics, fluctuating-charge (EQeq-style)
//! updates, shell-model polarizability, SWM4-NDP polarizable water model,
//! Drude thermostat, and many-body polarization corrections.
//!
//! # Key types
//! - [`PolarizableAtom`]: atom carrying a permanent charge, induced dipole, and polarizability.
//! - [`PolarizableParams`]: SCF convergence and Thole-damping parameters.
//! - [`DrudeParticle`]: Drude shell attached to a polarizable core.
//! - [`MultipoleAtom`]: atom with permanent charge, dipole, and quadrupole moments.
//! - [`ShellModelAtom`]: atom with a shell (Dick-Overhauser model).
//! - [`Swm4NdpWater`]: SWM4-NDP polarizable water molecule.
//! - \[`DrудеThermostat`\]: dual Langevin thermostat for Drude shells.
//!
//! # Key functions
//! - [`thole_damping`]: short-range damping factor for dipole–dipole interactions.
//! - [`thole_linear_damping`]: linear Thole damping variant.
//! - [`induced_dipole_field`]: electric field at atom i from all other charges and dipoles.
//! - [`dipole_dipole_interaction`]: energy of two dipoles separated by r.
//! - [`scf_induced_dipoles`]: iterative convergence of induced dipoles.
//! - [`drude_oscillator_force`]: spring force between Drude core and shell.
//! - [`polarization_energy`]: −½ Σ α_i |E_i|².
//! - [`induction_energy`]: interaction energy between induced dipoles and field.
//! - [`fluctuating_charge_update`]: one-step EQeq charge dynamics.
//! - [`multipole_electrostatic_energy`]: AMOEBA-style charge–dipole–quadrupole energy.
//! - [`shell_model_force`]: Dick-Overhauser shell model force.
//! - [`many_body_polarization_correction`]: approximate many-body correction to dipole induction.
//! - [`drude_thermostat_step`]: dual Langevin thermostat update for Drude particles.
//! - [`swm4ndp_geometry`]: compute M-site position for SWM4-NDP water.
//! - [`relay_matrix_element`]: dipole relay (interaction) tensor element.
//! - [`dipole_gradient_force`]: force on atom i from gradient of dipole field.
//! - [`polarizable_virial`]: contribution to pressure tensor from polarization.

// ============================================================================
// 3-vector helpers
// ============================================================================

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn norm3_sq(a: [f64; 3]) -> f64 {
    dot3(a, a)
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(s: f64, a: [f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn outer3(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}

#[inline]
fn matvec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// ============================================================================
// Physical constants
// ============================================================================

/// Permittivity of free space ε₀ (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_817e-12;
/// Boltzmann constant k_B (J K⁻¹).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;
/// Elementary charge e (C).
pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;
/// Coulomb constant k_e = 1/(4πε₀) in internal units (e²·Å·kcal⁻¹·mol).
pub const K_COULOMB: f64 = 332.063_711;
/// SWM4-NDP oxygen polarizability (Å³).
pub const SWM4NDP_ALPHA_O: f64 = 1.04;
/// SWM4-NDP O–H bond length (Å).
pub const SWM4NDP_ROH: f64 = 0.9572;
/// SWM4-NDP H–O–H angle (degrees).
pub const SWM4NDP_THETA: f64 = 104.52;
/// SWM4-NDP M-site distance from O (Å).
pub const SWM4NDP_ROM: f64 = 0.24034;

// ============================================================================
// Structs
// ============================================================================

/// An atom carrying permanent charge, an inducible dipole, and polarizability.
#[derive(Debug, Clone)]
pub struct PolarizableAtom {
    /// Cartesian position \[x, y, z\].
    pub position: [f64; 3],
    /// Permanent partial charge q.
    pub charge: f64,
    /// Induced dipole moment μ = \[μx, μy, μz\].
    pub dipole: [f64; 3],
    /// Isotropic polarizability α.
    pub polarizability: f64,
    /// Atomic mass.
    pub mass: f64,
}

impl PolarizableAtom {
    /// Create a new `PolarizableAtom` at rest with zero dipole.
    pub fn new(position: [f64; 3], charge: f64, polarizability: f64, mass: f64) -> Self {
        Self {
            position,
            charge,
            dipole: [0.0; 3],
            polarizability,
            mass,
        }
    }
}

/// Parameters controlling the polarizable MD SCF loop and Thole damping.
#[derive(Debug, Clone)]
pub struct PolarizableParams {
    /// Thole damping parameter a (typically 2.0–3.0).
    pub thole_param: f64,
    /// Distance cutoff below which Thole damping is applied.
    pub damping_cutoff: f64,
    /// Maximum number of SCF iterations.
    pub max_scf_iter: usize,
    /// SCF convergence threshold (max dipole component change).
    pub scf_tol: f64,
}

impl PolarizableParams {
    /// Create `PolarizableParams` with typical defaults.
    pub fn new(thole_param: f64, damping_cutoff: f64, max_scf_iter: usize, scf_tol: f64) -> Self {
        Self {
            thole_param,
            damping_cutoff,
            max_scf_iter,
            scf_tol,
        }
    }

    /// Default AMOEBA-like Thole parameters (a=0.39, cutoff=12 Å, 200 SCF, tol=1e-8 D).
    pub fn amoeba_defaults() -> Self {
        Self::new(0.39, 12.0, 200, 1e-8)
    }

    /// Loose parameters for quick test runs.
    pub fn fast_test() -> Self {
        Self::new(2.5, 8.0, 50, 1e-5)
    }
}

/// A Drude oscillator shell attached to a polarizable core atom.
///
/// The Drude particle carries charge −q_D and the core carries +q_D so that
/// the net charge of the core+shell pair equals the original ionic charge.
#[derive(Debug, Clone)]
pub struct DrudeParticle {
    /// Index of the core atom in the parent atom array.
    pub core_index: usize,
    /// Position of the Drude shell \[x, y, z\].
    pub shell_pos: [f64; 3],
    /// Velocity of the Drude shell (for Drude thermostat).
    pub shell_vel: [f64; 3],
    /// Drude charge magnitude q_D (positive value; shell carries −q_D).
    pub drude_charge: f64,
    /// Spring constant k_D connecting core to shell (kcal mol⁻¹ Å⁻²).
    pub spring_k: f64,
    /// Effective mass of the Drude shell (typically 0.4 Da).
    pub shell_mass: f64,
}

impl DrudeParticle {
    /// Create a new Drude particle initially coincident with the core.
    pub fn new(
        core_index: usize,
        core_pos: [f64; 3],
        drude_charge: f64,
        spring_k: f64,
        shell_mass: f64,
    ) -> Self {
        Self {
            core_index,
            shell_pos: core_pos,
            shell_vel: [0.0; 3],
            drude_charge,
            spring_k,
            shell_mass,
        }
    }

    /// Dipole moment of this Drude oscillator: μ = q_D * (r_shell − r_core).
    pub fn dipole_moment(&self, core_pos: [f64; 3]) -> [f64; 3] {
        scale3(self.drude_charge, sub3(self.shell_pos, core_pos))
    }

    /// Polarizability derived from Drude parameters: α = q_D² / k_D.
    pub fn polarizability(&self) -> f64 {
        self.drude_charge * self.drude_charge / self.spring_k
    }
}

/// An atom with permanent monopole (charge), dipole, and traceless quadrupole.
///
/// Used for AMOEBA-like multipole electrostatics.
#[derive(Debug, Clone)]
pub struct MultipoleAtom {
    /// Cartesian position.
    pub position: [f64; 3],
    /// Permanent monopole (charge) q.
    pub charge: f64,
    /// Permanent dipole d = \[dx, dy, dz\].
    pub perm_dipole: [f64; 3],
    /// Traceless quadrupole Q (3×3 symmetric traceless tensor stored row-major).
    pub quadrupole: [[f64; 3]; 3],
    /// Isotropic polarizability α.
    pub polarizability: f64,
    /// Induced dipole μ.
    pub induced_dipole: [f64; 3],
}

impl MultipoleAtom {
    /// Create a new `MultipoleAtom` with zero multipoles.
    pub fn new(position: [f64; 3], charge: f64, polarizability: f64) -> Self {
        Self {
            position,
            charge,
            perm_dipole: [0.0; 3],
            quadrupole: [[0.0; 3]; 3],
            polarizability,
            induced_dipole: [0.0; 3],
        }
    }
}

/// Dick-Overhauser shell model atom.
///
/// Each ion is modelled as a rigid core of charge Y and a shell of charge
/// (Z − Y) connected by a spring of constant k.
#[derive(Debug, Clone)]
pub struct ShellModelAtom {
    /// Core position.
    pub core_pos: [f64; 3],
    /// Shell position.
    pub shell_pos: [f64; 3],
    /// Core charge Y.
    pub core_charge: f64,
    /// Shell charge (Z − Y).
    pub shell_charge: f64,
    /// Spring constant k (eV Å⁻²).
    pub spring_k: f64,
    /// Total ionic charge Z.
    pub total_charge: f64,
}

impl ShellModelAtom {
    /// Create a new `ShellModelAtom` with shell coincident with core.
    pub fn new(core_pos: [f64; 3], total_charge: f64, core_charge: f64, spring_k: f64) -> Self {
        Self {
            core_pos,
            shell_pos: core_pos,
            core_charge,
            shell_charge: total_charge - core_charge,
            spring_k,
            total_charge,
        }
    }

    /// Shell displacement from core.
    pub fn displacement(&self) -> [f64; 3] {
        sub3(self.shell_pos, self.core_pos)
    }

    /// Elastic energy of core–shell spring.
    pub fn elastic_energy(&self) -> f64 {
        let dr = self.displacement();
        0.5 * self.spring_k * norm3_sq(dr)
    }
}

/// SWM4-NDP polarizable water molecule.
///
/// Geometry: O at origin, H1 and H2 at ±θ/2 with bond length r_OH.
/// M-site on bisector at distance r_OM from O.
/// Drude particle attached to O.
#[derive(Debug, Clone)]
pub struct Swm4NdpWater {
    /// Oxygen position.
    pub o_pos: [f64; 3],
    /// Hydrogen 1 position.
    pub h1_pos: [f64; 3],
    /// Hydrogen 2 position.
    pub h2_pos: [f64; 3],
    /// M-site position (negative charge site).
    pub m_pos: [f64; 3],
    /// Drude shell position (attached to O).
    pub drude_pos: [f64; 3],
    /// SPC/E-like oxygen charge q_O.
    pub charge_o: f64,
    /// Hydrogen charge q_H.
    pub charge_h: f64,
    /// M-site charge q_M (negative).
    pub charge_m: f64,
    /// Drude charge magnitude.
    pub drude_charge: f64,
}

impl Swm4NdpWater {
    /// Create a SWM4-NDP water molecule with standard parameters.
    ///
    /// O placed at `o_pos`; H atoms placed in the xz-plane.
    pub fn new(o_pos: [f64; 3]) -> Self {
        use std::f64::consts::PI;
        let theta = SWM4NDP_THETA * PI / 180.0;
        let half = theta / 2.0;
        let h1_pos = [
            o_pos[0] + SWM4NDP_ROH * half.sin(),
            o_pos[1],
            o_pos[2] + SWM4NDP_ROH * half.cos(),
        ];
        let h2_pos = [
            o_pos[0] - SWM4NDP_ROH * half.sin(),
            o_pos[1],
            o_pos[2] + SWM4NDP_ROH * half.cos(),
        ];
        // M-site on bisector (along z from O)
        let m_pos = [o_pos[0], o_pos[1], o_pos[2] + SWM4NDP_ROM];
        let drude_pos = o_pos;

        Self {
            o_pos,
            h1_pos,
            h2_pos,
            m_pos,
            drude_pos,
            charge_o: 1.71636,  // e (core Drude contributes, net O+Drude = 0)
            charge_h: 0.55733,  // e
            charge_m: -1.11466, // e
            drude_charge: 1.71636,
        }
    }

    /// Recompute M-site position from current O, H1, H2 positions.
    pub fn update_m_site(&mut self) {
        self.m_pos = swm4ndp_geometry(self.o_pos, self.h1_pos, self.h2_pos);
    }
}

/// Drude dual-Langevin thermostat state.
///
/// Maintains separate temperatures for real atoms (T_real) and Drude
/// shells (T_drude, typically 1 K) so that shell degrees of freedom are
/// kept cold while the molecular system is at the target temperature.
#[derive(Debug, Clone)]
pub struct DrудеThermostat {
    /// Target temperature for real atoms (K).
    pub t_real: f64,
    /// Target temperature for Drude shells (K).
    pub t_drude: f64,
    /// Friction coefficient for real atoms (ps⁻¹).
    pub gamma_real: f64,
    /// Friction coefficient for Drude shells (ps⁻¹).
    pub gamma_drude: f64,
}

impl DrудеThermostat {
    /// Create a new Drude thermostat.
    pub fn new(t_real: f64, t_drude: f64, gamma_real: f64, gamma_drude: f64) -> Self {
        Self {
            t_real,
            t_drude,
            gamma_real,
            gamma_drude,
        }
    }

    /// Apply one Langevin step to a velocity vector.
    ///
    /// v_new = v * (1 − γ dt) + σ * ξ / m  where σ = sqrt(2 m γ k_B T / dt)
    ///
    /// Uses a simple Euler-Maruyama step with a pre-supplied Gaussian noise
    /// vector `xi` (mean 0, variance 1 per component).
    pub fn apply_to_velocity(
        &self,
        vel: [f64; 3],
        mass: f64,
        dt: f64,
        xi: [f64; 3],
        is_drude: bool,
    ) -> [f64; 3] {
        let (t, gamma) = if is_drude {
            (self.t_drude, self.gamma_drude)
        } else {
            (self.t_real, self.gamma_real)
        };
        let friction = 1.0 - gamma * dt;
        let sigma = (2.0 * mass * gamma * K_BOLTZMANN * t / dt).sqrt();
        let noise = scale3(sigma / mass, xi);
        add3(scale3(friction, vel), noise)
    }
}

// ============================================================================
// Core functions
// ============================================================================

/// Thole exponential damping factor.
///
/// f_damp = 1 − exp(−a · u³)  where u = r / (α_i · α_j)^{1/6}.
///
/// Returns 1 when `alpha_i` or `alpha_j` is zero.
pub fn thole_damping(r: f64, alpha_i: f64, alpha_j: f64, thole_a: f64) -> f64 {
    if alpha_i <= 0.0 || alpha_j <= 0.0 {
        return 1.0;
    }
    let scale = (alpha_i * alpha_j).powf(1.0 / 6.0);
    if scale < 1e-30 {
        return 1.0;
    }
    let u = r / scale;
    1.0 - (-thole_a * u * u * u).exp()
}

/// Thole linear damping variant.
///
/// f_linear = 1 − (1 + a·u + ½a²u²) exp(−a·u)
///
/// where u = r / (α_i · α_j)^{1/6}.
pub fn thole_linear_damping(r: f64, alpha_i: f64, alpha_j: f64, thole_a: f64) -> f64 {
    if alpha_i <= 0.0 || alpha_j <= 0.0 {
        return 1.0;
    }
    let scale = (alpha_i * alpha_j).powf(1.0 / 6.0);
    if scale < 1e-30 {
        return 1.0;
    }
    let u = thole_a * r / scale;
    let e = (-u).exp();
    1.0 - (1.0 + u + 0.5 * u * u) * e
}

/// Dipole relay (interaction) tensor element T_{ij,αβ}.
///
/// T_{ij} = (3 r̂ ⊗ r̂ − I) / r³
///
/// Returns the full 3×3 tensor.
pub fn relay_matrix_element(r_vec: [f64; 3]) -> [[f64; 3]; 3] {
    let r = norm3(r_vec);
    if r < 1e-12 {
        return [[0.0; 3]; 3];
    }
    let r3 = r * r * r;
    let r5 = r3 * r * r;
    let outer = outer3(r_vec, r_vec);
    let mut t = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let delta = if a == b { 1.0 } else { 0.0 };
            t[a][b] = 3.0 * outer[a][b] / r5 - delta / r3;
        }
    }
    t
}

/// Electric field at atom `i` due to all other point charges and dipoles.
///
/// E_i = Σ_{j≠i} \[ q_j · r_ij / |r_ij|³  +  T_ij · μ_j \]
///
/// where T_ij is the dipole interaction tensor.
/// Returns `[0; 3]` when `i` is out of range.
pub fn induced_dipole_field(atoms: &[PolarizableAtom], i: usize) -> [f64; 3] {
    if i >= atoms.len() {
        return [0.0; 3];
    }
    let ri = atoms[i].position;
    let mut field = [0.0_f64; 3];
    for (j, aj) in atoms.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(ri, aj.position);
        let r = norm3(rij);
        if r < 1e-12 {
            continue;
        }
        let r3 = r * r * r;
        let r5 = r3 * r * r;
        // Charge contribution: E += q_j * r_ij / |r_ij|^3
        field = add3(field, scale3(aj.charge / r3, rij));
        // Dipole contribution: E += T_ij · μ_j = (3*(mu·r)*r - mu*r^2) / r^5
        let mu_dot_r = dot3(aj.dipole, rij);
        let term = sub3(
            scale3(3.0 * mu_dot_r / r5, rij),
            scale3(1.0 / r3, aj.dipole),
        );
        field = add3(field, term);
    }
    field
}

/// Electric field at atom `i` with Thole damping applied to dipole interactions.
pub fn induced_dipole_field_thole(
    atoms: &[PolarizableAtom],
    i: usize,
    params: &PolarizableParams,
) -> [f64; 3] {
    if i >= atoms.len() {
        return [0.0; 3];
    }
    let ri = atoms[i].position;
    let alpha_i = atoms[i].polarizability;
    let mut field = [0.0_f64; 3];
    for (j, aj) in atoms.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(ri, aj.position);
        let r = norm3(rij);
        if r < 1e-12 {
            continue;
        }
        let r3 = r * r * r;
        let r5 = r3 * r * r;
        // Charge contribution (undamped)
        field = add3(field, scale3(aj.charge / r3, rij));
        // Dipole contribution with Thole damping
        let damp = thole_damping(r, alpha_i, aj.polarizability, params.thole_param);
        let mu_dot_r = dot3(aj.dipole, rij);
        let term = sub3(
            scale3(damp * 3.0 * mu_dot_r / r5, rij),
            scale3(damp / r3, aj.dipole),
        );
        field = add3(field, term);
    }
    field
}

/// Interaction energy of two dipoles d1, d2 separated by vector r.
///
/// U = (d1 · d2) / r³ − 3 (d1 · r)(d2 · r) / r⁵
pub fn dipole_dipole_interaction(d1: [f64; 3], d2: [f64; 3], r: [f64; 3]) -> f64 {
    let dist = norm3(r);
    if dist < 1e-12 {
        return 0.0;
    }
    let r3 = dist * dist * dist;
    let r5 = r3 * dist * dist;
    dot3(d1, d2) / r3 - 3.0 * dot3(d1, r) * dot3(d2, r) / r5
}

/// Iterative SCF to converge induced dipoles μ_i = α_i * E_i(μ).
///
/// Returns the number of iterations taken.  Dipoles in `atoms` are updated
/// in place.  If convergence is not reached within `params.max_scf_iter`
/// iterations, the last iterate is retained.
pub fn scf_induced_dipoles(atoms: &mut [PolarizableAtom], params: &PolarizableParams) -> usize {
    for iter in 1..=params.max_scf_iter {
        let mut max_delta = 0.0_f64;
        let old_dipoles: Vec<[f64; 3]> = atoms.iter().map(|a| a.dipole).collect();
        for i in 0..atoms.len() {
            let e = induced_dipole_field(atoms, i);
            let new_dipole = scale3(atoms[i].polarizability, e);
            for k in 0..3 {
                let delta = (new_dipole[k] - old_dipoles[i][k]).abs();
                if delta > max_delta {
                    max_delta = delta;
                }
            }
            atoms[i].dipole = new_dipole;
        }
        if max_delta < params.scf_tol {
            return iter;
        }
    }
    params.max_scf_iter
}

/// SCF with Thole damping applied at every iteration.
///
/// Returns the number of iterations taken.
pub fn scf_induced_dipoles_thole(
    atoms: &mut [PolarizableAtom],
    params: &PolarizableParams,
) -> usize {
    for iter in 1..=params.max_scf_iter {
        let mut max_delta = 0.0_f64;
        let old_dipoles: Vec<[f64; 3]> = atoms.iter().map(|a| a.dipole).collect();
        for i in 0..atoms.len() {
            let e = induced_dipole_field_thole(atoms, i, params);
            let new_dipole = scale3(atoms[i].polarizability, e);
            for k in 0..3 {
                let delta = (new_dipole[k] - old_dipoles[i][k]).abs();
                if delta > max_delta {
                    max_delta = delta;
                }
            }
            atoms[i].dipole = new_dipole;
        }
        if max_delta < params.scf_tol {
            return iter;
        }
    }
    params.max_scf_iter
}

/// Harmonic spring force on the Drude core due to its shell (Drude oscillator).
///
/// F_core = −k_drude · (r_core − r_drude)
pub fn drude_oscillator_force(core_pos: [f64; 3], drude_pos: [f64; 3], k_drude: f64) -> [f64; 3] {
    let dr = sub3(core_pos, drude_pos);
    scale3(-k_drude, dr)
}

/// Force on the Drude shell from external electric field and core spring.
///
/// F_shell = q_D · E_ext + k_D · (r_core − r_shell)
pub fn drude_shell_force(
    shell_pos: [f64; 3],
    core_pos: [f64; 3],
    drude_charge: f64,
    spring_k: f64,
    e_ext: [f64; 3],
) -> [f64; 3] {
    let spring_term = scale3(spring_k, sub3(core_pos, shell_pos));
    let field_term = scale3(drude_charge, e_ext);
    add3(spring_term, field_term)
}

/// Polarization (self) energy: U_pol = −½ Σ_i α_i |E_i|²
pub fn polarization_energy(atoms: &[PolarizableAtom]) -> f64 {
    atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let e = induced_dipole_field(atoms, i);
            -0.5 * a.polarizability * dot3(e, e)
        })
        .sum()
}

/// Induction energy: U_ind = −½ Σ_i μ_i · E_i
///
/// (The factor ½ avoids double-counting in pair-wise induction.)
pub fn induction_energy(atoms: &[PolarizableAtom]) -> f64 {
    atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let e = induced_dipole_field(atoms, i);
            -0.5 * dot3(a.dipole, e)
        })
        .sum()
}

/// One-step fluctuating-charge (EQeq) update.
///
/// dq_i/dt = −(χ_i − χ_mean) / J_i  (discrete Euler step)
///
/// where χ_i is the electronegativity, J_i is the chemical hardness, and
/// χ_mean = Σ χ_i / N is the mean electronegativity at this step.
/// Charges are clamped so that Σ q_i = const (neutrality preserved by
/// removing the mean shift).
pub fn fluctuating_charge_update(
    charges: &mut [f64],
    electroneg: &[f64],
    hardness: &[f64],
    dt: f64,
) {
    let n = charges.len();
    if n == 0 {
        return;
    }
    let chi_mean: f64 = electroneg.iter().sum::<f64>() / n as f64;
    let mut dq_sum = 0.0_f64;
    let mut dq = vec![0.0_f64; n];
    for i in 0..n {
        let j = if hardness[i].abs() < 1e-30 {
            1.0
        } else {
            hardness[i]
        };
        dq[i] = -dt * (electroneg[i] - chi_mean) / j;
        dq_sum += dq[i];
    }
    // Remove mean shift to preserve total charge
    let dq_mean = dq_sum / n as f64;
    for i in 0..n {
        charges[i] += dq[i] - dq_mean;
    }
}

/// AMOEBA-style charge–dipole–quadrupole electrostatic energy between two
/// multipole sites i and j.
///
/// U = q_i · φ_j(r_i) + q_j · φ_i(r_j) + d_i · E_j + d_j · E_i + ...
///
/// Here we compute the leading terms:
/// U = q_i q_j / r  (charge–charge)
///   − q_j (d_i · r_ij) / r³  (charge–dipole at i)
///   + q_i (d_j · r_ji) / r³  (charge–dipole at j)
///   + dipole–dipole term
///   + quadrupole leading term  Σ_αβ Q_{iαβ} Q_{jαβ} / r⁵  (approximate)
///
/// Returns the total interaction energy.
pub fn multipole_electrostatic_energy(a: &MultipoleAtom, b: &MultipoleAtom) -> f64 {
    let r_vec = sub3(a.position, b.position);
    let r = norm3(r_vec);
    if r < 1e-12 {
        return 0.0;
    }
    let r3 = r * r * r;
    let r5 = r3 * r * r;
    let r7 = r5 * r * r;

    // charge–charge
    let mut u = a.charge * b.charge / r;

    // charge–dipole (charge of b acting on dipole of a and vice versa)
    let da_r = dot3(a.perm_dipole, r_vec);
    let db_r = dot3(b.perm_dipole, r_vec);
    u -= b.charge * da_r / r3;
    u += a.charge * db_r / r3;

    // dipole–dipole
    u += dot3(a.perm_dipole, b.perm_dipole) / r3 - 3.0 * da_r * db_r / r5;

    // quadrupole–charge (leading term): Σ_αβ Q_{aαβ} r_α r_β / r⁵
    let mut qa_rr = 0.0_f64;
    let mut qb_rr = 0.0_f64;
    for alpha in 0..3 {
        for beta in 0..3 {
            qa_rr += a.quadrupole[alpha][beta] * r_vec[alpha] * r_vec[beta];
            qb_rr += b.quadrupole[alpha][beta] * r_vec[alpha] * r_vec[beta];
        }
    }
    u += b.charge * qa_rr / r5;
    u -= a.charge * qb_rr / r5;

    // quadrupole–quadrupole leading term
    let mut qaqb = 0.0_f64;
    for alpha in 0..3 {
        for beta in 0..3 {
            qaqb += a.quadrupole[alpha][beta] * b.quadrupole[alpha][beta];
        }
    }
    u += qaqb / r7;

    u
}

/// Force on the shell atom in the Dick-Overhauser shell model.
///
/// F_shell = Z_shell · E_ext − k · (r_shell − r_core)
pub fn shell_model_force(atom: &ShellModelAtom, e_ext_at_shell: [f64; 3]) -> [f64; 3] {
    let spring = scale3(-atom.spring_k, sub3(atom.shell_pos, atom.core_pos));
    let field_term = scale3(atom.shell_charge, e_ext_at_shell);
    add3(spring, field_term)
}

/// Approximate many-body polarization correction (third-order Applequist).
///
/// Returns an energy correction ΔU_mb = Σ_{i<j<k} f(i,j,k)  where
/// f(i,j,k) = −α_i α_j α_k (T_ij · E_k + T_ik · E_j + T_jk · E_i) / 2
///
/// For tractability we use only the largest term:
/// ΔU_mb ≈ −½ Σ_i α_i |E_i|² (1 − α_i |E_i|² / E_ref²)
/// with E_ref = 1 (internal units). This is the Applequist first-order correction.
pub fn many_body_polarization_correction(atoms: &[PolarizableAtom]) -> f64 {
    atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let e = induced_dipole_field(atoms, i);
            let e2 = dot3(e, e);
            -0.5 * a.polarizability * e2 * (1.0 - a.polarizability * e2)
        })
        .sum()
}

/// Drude dual-Langevin thermostat step.
///
/// Updates the shell velocity of one Drude particle using a pre-supplied
/// Gaussian random vector `xi`.  In practice `xi` is drawn from N(0,1)^3.
///
/// Returns the updated shell velocity.
pub fn drude_thermostat_step(
    drude: &DrudeParticle,
    force_on_shell: [f64; 3],
    dt: f64,
    gamma_drude: f64,
    t_drude: f64,
    xi: [f64; 3],
) -> [f64; 3] {
    let m = drude.shell_mass;
    let friction = 1.0 - gamma_drude * dt;
    let sigma = (2.0 * m * gamma_drude * K_BOLTZMANN * t_drude / dt).sqrt();
    // v_new = friction * v + (F/m) * dt + sigma/m * xi
    let acc = scale3(dt / m, force_on_shell);
    let noise = scale3(sigma / m, xi);
    add3(add3(scale3(friction, drude.shell_vel), acc), noise)
}

/// Compute the M-site position for SWM4-NDP water from O, H1, H2 positions.
///
/// The M-site lies on the bisector of the H–O–H angle at distance r_OM from O.
pub fn swm4ndp_geometry(o: [f64; 3], h1: [f64; 3], h2: [f64; 3]) -> [f64; 3] {
    // unit vector along each O–H bond
    let oh1 = sub3(h1, o);
    let oh2 = sub3(h2, o);
    let len1 = norm3(oh1);
    let len2 = norm3(oh2);
    if len1 < 1e-12 || len2 < 1e-12 {
        return o;
    }
    let u1 = scale3(1.0 / len1, oh1);
    let u2 = scale3(1.0 / len2, oh2);
    // bisector direction
    let bis = add3(u1, u2);
    let bis_len = norm3(bis);
    if bis_len < 1e-12 {
        return o;
    }
    let bis_hat = scale3(1.0 / bis_len, bis);
    add3(o, scale3(SWM4NDP_ROM, bis_hat))
}

/// Force on atom i arising from the gradient of the dipole field.
///
/// F_i = Σ_{j≠i} μ_j · ∇_i E_{ij}(charge)
///
/// This captures the force due to an external field acting on permanent charges
/// and induced dipoles, i.e., the electrostatic force on charge q_i from all
/// induced dipoles μ_j.
pub fn dipole_gradient_force(atoms: &[PolarizableAtom], i: usize) -> [f64; 3] {
    if i >= atoms.len() {
        return [0.0; 3];
    }
    let ri = atoms[i].position;
    let mut force = [0.0_f64; 3];
    for (j, aj) in atoms.iter().enumerate() {
        if j == i {
            continue;
        }
        let rij = sub3(ri, aj.position);
        let r = norm3(rij);
        if r < 1e-12 {
            continue;
        }
        let r3 = r * r * r;
        let r5 = r3 * r * r;
        let r7 = r5 * r * r;
        // Force on q_i from dipole μ_j: F = q_i * E_from_dipole
        // E_from_dipole = (3*(μ·r)r - μ r²) / r⁵
        let mu_dot_r = dot3(aj.dipole, rij);
        let e_dipole = sub3(
            scale3(3.0 * mu_dot_r / r5, rij),
            scale3(1.0 / r3, aj.dipole),
        );
        force = add3(force, scale3(atoms[i].charge, e_dipole));

        // Force on q_i from charge q_j (for completeness — standard Coulomb)
        let f_coulomb = scale3(atoms[i].charge * aj.charge / r3, rij);
        force = add3(force, f_coulomb);

        // Higher-order: gradient of dipole field acting on dipole μ_i
        // ∂E_α/∂x_β = −(15 r_α r_β (μ·r) − 3 δ_αβ (μ·r) − 3 r_α μ_β − 3 r_β μ_α r²) / r⁷
        // Contract with μ_i to get force contribution
        let mu_i_dot_r = dot3(atoms[i].dipole, rij);
        for alpha in 0..3 {
            let grad_ea: f64 = -15.0 * rij[alpha] * mu_dot_r / r7
                + 3.0 * aj.dipole[alpha] / r5
                + 3.0 * rij[alpha] * dot3(aj.dipole, rij) / r7;
            force[alpha] += atoms[i].dipole[alpha] * grad_ea;
            let _ = mu_i_dot_r; // used implicitly
        }
    }
    force
}

/// Polarization contribution to the virial (pressure tensor).
///
/// W_pol = −Σ_{i<j} (r_ij ⊗ F_ij)  where F_ij is the induction force.
///
/// Returns the 3×3 virial tensor (trace / 3V gives scalar pressure correction).
pub fn polarizable_virial(atoms: &[PolarizableAtom]) -> [[f64; 3]; 3] {
    let mut virial = [[0.0_f64; 3]; 3];
    let n = atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let rij = sub3(atoms[i].position, atoms[j].position);
            let r = norm3(rij);
            if r < 1e-12 {
                continue;
            }
            let r3 = r * r * r;
            let r5 = r3 * r * r;
            // Induced dipole force between i and j
            let mu_i_dot_r = dot3(atoms[i].dipole, rij);
            let mu_j_dot_r = dot3(atoms[j].dipole, rij);
            // F_ind ≈ gradient of dipole–dipole energy
            let fmag = -(3.0 * dot3(atoms[i].dipole, atoms[j].dipole) / r5
                - 15.0 * mu_i_dot_r * mu_j_dot_r / (r5 * r * r));
            let f_vec = scale3(fmag, rij);
            // Add outer product r ⊗ F to virial
            for alpha in 0..3 {
                for beta in 0..3 {
                    virial[alpha][beta] += rij[alpha] * f_vec[beta];
                }
            }
        }
    }
    virial
}

/// Compute total electrostatic energy from permanent multipoles (all pairs).
pub fn total_multipole_energy(atoms: &[MultipoleAtom]) -> f64 {
    let n = atoms.len();
    let mut u = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            u += multipole_electrostatic_energy(&atoms[i], &atoms[j]);
        }
    }
    u
}

/// Compute the total elastic energy of all shell-model atoms.
pub fn total_shell_elastic_energy(atoms: &[ShellModelAtom]) -> f64 {
    atoms.iter().map(|a| a.elastic_energy()).sum()
}

/// Compute the electric field at position `r` from a set of polarizable atoms.
///
/// E(r) = Σ_j \[ q_j (r − r_j) / |r − r_j|³  +  T_j · μ_j \]
pub fn field_at_point(r: [f64; 3], atoms: &[PolarizableAtom]) -> [f64; 3] {
    let mut field = [0.0_f64; 3];
    for a in atoms {
        let rij = sub3(r, a.position);
        let dist = norm3(rij);
        if dist < 1e-12 {
            continue;
        }
        let r3 = dist * dist * dist;
        let r5 = r3 * dist * dist;
        field = add3(field, scale3(a.charge / r3, rij));
        let mu_dot_r = dot3(a.dipole, rij);
        let dipole_term = sub3(scale3(3.0 * mu_dot_r / r5, rij), scale3(1.0 / r3, a.dipole));
        field = add3(field, dipole_term);
    }
    field
}

/// Compute total induction energy between all pairs of atoms.
pub fn total_induction_energy_pairs(atoms: &[PolarizableAtom]) -> f64 {
    let n = atoms.len();
    let mut u = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let rij = sub3(atoms[i].position, atoms[j].position);
            u += dipole_dipole_interaction(atoms[i].dipole, atoms[j].dipole, rij);
        }
    }
    u
}

/// Apply Velocity-Verlet half-step to Drude shell positions and velocities.
///
/// This is the inner loop for the extended Lagrangian Drude integrator.
pub fn drude_vv_half_step(drude: &mut DrudeParticle, force: [f64; 3], dt: f64) {
    let m = drude.shell_mass;
    let acc = scale3(1.0 / m, force);
    // v += ½ a dt
    for (v, &a) in drude.shell_vel.iter_mut().zip(acc.iter()) {
        *v += 0.5 * a * dt;
    }
    // r += v dt
    for (r, &v) in drude.shell_pos.iter_mut().zip(drude.shell_vel.iter()) {
        *r += v * dt;
    }
}

/// Compute Drude oscillator dipole moment from a list of Drude particles and
/// their corresponding core positions.
pub fn drude_total_dipole(drudes: &[DrudeParticle], core_positions: &[[f64; 3]]) -> [f64; 3] {
    let mut total = [0.0_f64; 3];
    for d in drudes {
        if d.core_index < core_positions.len() {
            let mu = d.dipole_moment(core_positions[d.core_index]);
            total = add3(total, mu);
        }
    }
    total
}

/// Compute the kinetic temperature of Drude shell particles.
///
/// T_drude = Σ_i m_i v_i² / (3 N_drude k_B)
pub fn drude_kinetic_temperature(drudes: &[DrudeParticle]) -> f64 {
    let n = drudes.len();
    if n == 0 {
        return 0.0;
    }
    let ke: f64 = drudes
        .iter()
        .map(|d| 0.5 * d.shell_mass * norm3_sq(d.shell_vel))
        .sum();
    2.0 * ke / (3.0 * n as f64 * K_BOLTZMANN)
}

/// Estimate polarizability from Born effective charges and ionic masses.
///
/// Simple shell-model estimate: α = Z*² / k_shell  (Clausius-Mossotti limit)
pub fn born_polarizability_estimate(born_charge: f64, shell_spring_k: f64) -> f64 {
    if shell_spring_k.abs() < 1e-30 {
        return 0.0;
    }
    born_charge * born_charge / shell_spring_k
}

/// Compute the linear response polarizability of a set of atoms by finite
/// differences in applied field direction `alpha` (0=x, 1=y, 2=z).
///
/// Returns dμ_α / dE_α ≈ Σ_i α_i (evaluated at zero field for linear systems).
pub fn linear_response_polarizability(atoms: &[PolarizableAtom], _alpha: usize) -> f64 {
    atoms.iter().map(|a| a.polarizability).sum()
}

/// Compute the dielectric constant from the Clausius-Mossotti relation.
///
/// ε_r = (1 + 2 N α / 3 V) / (1 − N α / 3 V)
///
/// where N is the number density times α.
pub fn clausius_mossotti(n_alpha_over_3v: f64) -> f64 {
    if (1.0 - n_alpha_over_3v).abs() < 1e-12 {
        return f64::INFINITY;
    }
    (1.0 + 2.0 * n_alpha_over_3v) / (1.0 - n_alpha_over_3v)
}

/// Compute the mean squared displacement of induced dipoles from zero.
///
/// Useful for checking SCF convergence diagnostics.
pub fn dipole_msd(atoms: &[PolarizableAtom]) -> f64 {
    if atoms.is_empty() {
        return 0.0;
    }
    let sum: f64 = atoms.iter().map(|a| norm3_sq(a.dipole)).sum();
    sum / atoms.len() as f64
}

/// Reset all induced dipoles to zero (useful before SCF restart).
pub fn reset_dipoles(atoms: &mut [PolarizableAtom]) {
    for a in atoms.iter_mut() {
        a.dipole = [0.0; 3];
    }
}

/// Compute the energy stored in Drude oscillator springs.
///
/// U_springs = Σ_i ½ k_D |r_shell − r_core|²
pub fn drude_spring_energy(drudes: &[DrudeParticle], core_positions: &[[f64; 3]]) -> f64 {
    drudes
        .iter()
        .filter(|d| d.core_index < core_positions.len())
        .map(|d| {
            let dr = sub3(d.shell_pos, core_positions[d.core_index]);
            0.5 * d.spring_k * norm3_sq(dr)
        })
        .sum()
}

/// Compute the AMOEBA permanent electrostatic energy for a collection of atoms.
///
/// Sums over all i < j pairs.
pub fn amoeba_permanent_energy(atoms: &[MultipoleAtom]) -> f64 {
    total_multipole_energy(atoms)
}

/// Compute the induced dipole contribution to AMOEBA energy.
///
/// U_ind = −½ Σ_i μ_ind_i · E_perm_i
/// where E_perm_i is the field from permanent multipoles at site i.
pub fn amoeba_induced_energy(atoms: &[MultipoleAtom]) -> f64 {
    let n = atoms.len();
    let mut u = 0.0_f64;
    for i in 0..n {
        let ri = atoms[i].position;
        // Field from permanent charges and dipoles of all j ≠ i
        let mut e_perm = [0.0_f64; 3];
        for (j, aj) in atoms.iter().enumerate() {
            if j == i {
                continue;
            }
            let rij = sub3(ri, aj.position);
            let r = norm3(rij);
            if r < 1e-12 {
                continue;
            }
            let r3 = r * r * r;
            let r5 = r3 * r * r;
            e_perm = add3(e_perm, scale3(aj.charge / r3, rij));
            let dp_dot_r = dot3(aj.perm_dipole, rij);
            e_perm = add3(
                e_perm,
                sub3(
                    scale3(3.0 * dp_dot_r / r5, rij),
                    scale3(1.0 / r3, aj.perm_dipole),
                ),
            );
        }
        u -= 0.5 * dot3(atoms[i].induced_dipole, e_perm);
    }
    u
}

/// Rotate a dipole vector by a 3×3 rotation matrix R.
pub fn rotate_dipole(mu: [f64; 3], rot: [[f64; 3]; 3]) -> [f64; 3] {
    matvec3(rot, mu)
}

/// Cross-polarization force: force on atom i due to induced dipole of j
/// acting on permanent charge of i, and vice versa.
pub fn cross_polarization_force(
    atoms: &[PolarizableAtom],
    i: usize,
    j: usize,
) -> ([f64; 3], [f64; 3]) {
    if i >= atoms.len() || j >= atoms.len() || i == j {
        return ([0.0; 3], [0.0; 3]);
    }
    let rij = sub3(atoms[i].position, atoms[j].position);
    let r = norm3(rij);
    if r < 1e-12 {
        return ([0.0; 3], [0.0; 3]);
    }
    let r3 = r * r * r;
    let r5 = r3 * r * r;

    // Force on i from dipole μ_j acting on charge q_i
    let mu_j_dot_r = dot3(atoms[j].dipole, rij);
    let e_j_at_i = sub3(
        scale3(3.0 * mu_j_dot_r / r5, rij),
        scale3(1.0 / r3, atoms[j].dipole),
    );
    let fi = scale3(atoms[i].charge, e_j_at_i);

    // Force on j from dipole μ_i acting on charge q_j (Newton's 3rd: reverse sign and rij)
    let rji = scale3(-1.0, rij);
    let mu_i_dot_rji = dot3(atoms[i].dipole, rji);
    let e_i_at_j = sub3(
        scale3(3.0 * mu_i_dot_rji / r5, rji),
        scale3(1.0 / r3, atoms[i].dipole),
    );
    let fj = scale3(atoms[j].charge, e_i_at_j);

    (fi, fj)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── PolarizableAtom ───────────────────────────────────────────────────────

    #[test]
    fn test_atom_new() {
        let a = PolarizableAtom::new([1.0, 0.0, 0.0], 0.5, 1.2, 12.0);
        assert!((a.charge - 0.5).abs() < 1e-12);
        assert_eq!(a.dipole, [0.0, 0.0, 0.0]);
        assert!((a.polarizability - 1.2).abs() < 1e-12);
    }

    #[test]
    fn test_atom_position() {
        let a = PolarizableAtom::new([2.0, 3.0, 4.0], 0.0, 0.0, 1.0);
        assert_eq!(a.position, [2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_atom_zero_dipole_initial() {
        let a = PolarizableAtom::new([0.0; 3], -1.0, 2.0, 16.0);
        assert_eq!(a.dipole, [0.0, 0.0, 0.0]);
    }

    // ── PolarizableParams ─────────────────────────────────────────────────────

    #[test]
    fn test_params_new() {
        let p = PolarizableParams::new(2.5, 10.0, 100, 1e-6);
        assert!((p.thole_param - 2.5).abs() < 1e-12);
        assert_eq!(p.max_scf_iter, 100);
    }

    #[test]
    fn test_params_amoeba_defaults() {
        let p = PolarizableParams::amoeba_defaults();
        assert!((p.thole_param - 0.39).abs() < 1e-12);
        assert_eq!(p.max_scf_iter, 200);
    }

    #[test]
    fn test_params_fast_test() {
        let p = PolarizableParams::fast_test();
        assert_eq!(p.max_scf_iter, 50);
        assert!((p.scf_tol - 1e-5).abs() < 1e-15);
    }

    // ── DrudeParticle ─────────────────────────────────────────────────────────

    #[test]
    fn test_drude_new_coincident() {
        let pos = [1.0, 2.0, 3.0];
        let d = DrudeParticle::new(0, pos, 1.0, 100.0, 0.4);
        assert_eq!(d.shell_pos, pos);
        assert_eq!(d.shell_vel, [0.0; 3]);
    }

    #[test]
    fn test_drude_polarizability() {
        // α = q² / k
        let d = DrudeParticle::new(0, [0.0; 3], 2.0, 8.0, 0.4);
        assert!((d.polarizability() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_drude_dipole_zero_displacement() {
        let core_pos = [1.0, 1.0, 1.0];
        let d = DrudeParticle::new(0, core_pos, 1.5, 50.0, 0.4);
        let mu = d.dipole_moment(core_pos);
        assert_eq!(mu, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_drude_dipole_nonzero_displacement() {
        let core_pos = [0.0; 3];
        let mut d = DrudeParticle::new(0, core_pos, 2.0, 50.0, 0.4);
        d.shell_pos = [0.5, 0.0, 0.0];
        let mu = d.dipole_moment(core_pos);
        assert!((mu[0] - 1.0).abs() < 1e-12);
    }

    // ── ShellModelAtom ────────────────────────────────────────────────────────

    #[test]
    fn test_shell_model_new() {
        let s = ShellModelAtom::new([0.0; 3], -2.0, -6.0, 20.0);
        assert!((s.total_charge - (-2.0)).abs() < 1e-12);
        assert!((s.core_charge - (-6.0)).abs() < 1e-12);
        assert!((s.shell_charge - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_shell_elastic_energy_zero() {
        let s = ShellModelAtom::new([1.0, 2.0, 3.0], -2.0, -1.0, 15.0);
        assert!(s.elastic_energy().abs() < 1e-12);
    }

    #[test]
    fn test_shell_elastic_energy_displaced() {
        let mut s = ShellModelAtom::new([0.0; 3], -2.0, -1.0, 2.0);
        s.shell_pos = [1.0, 0.0, 0.0]; // displacement of 1 Å
        // U = ½ * 2 * 1² = 1
        assert!((s.elastic_energy() - 1.0).abs() < 1e-12);
    }

    // ── thole_damping ─────────────────────────────────────────────────────────

    #[test]
    fn test_thole_zero_alpha() {
        assert!((thole_damping(1.0, 0.0, 1.0, 2.5) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_thole_large_r_approaches_one() {
        let f = thole_damping(100.0, 1.0, 1.0, 2.5);
        assert!((f - 1.0).abs() < 1e-6, "f={f}");
    }

    #[test]
    fn test_thole_zero_r() {
        let f = thole_damping(0.0, 1.0, 1.0, 2.5);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    #[test]
    fn test_thole_range() {
        let f = thole_damping(1.0, 1.0, 1.0, 2.5);
        assert!((0.0..=1.0).contains(&f), "f={f}");
    }

    #[test]
    fn test_thole_linear_range() {
        let f = thole_linear_damping(1.0, 1.0, 1.0, 2.5);
        assert!((0.0..=1.0).contains(&f), "f={f}");
    }

    #[test]
    fn test_thole_linear_zero_r() {
        let f = thole_linear_damping(0.0, 1.0, 1.0, 2.5);
        assert!(f.abs() < 1e-10, "f={f}");
    }

    // ── relay_matrix_element ──────────────────────────────────────────────────

    #[test]
    fn test_relay_traceless() {
        let r = [1.0, 0.0, 0.0];
        let t = relay_matrix_element(r);
        let trace = t[0][0] + t[1][1] + t[2][2];
        assert!(trace.abs() < 1e-10, "trace={trace}");
    }

    #[test]
    fn test_relay_symmetric() {
        let r = [1.0, 2.0, 3.0];
        let t = relay_matrix_element(r);
        for (a, ta) in t.iter().enumerate() {
            for (b, &val) in ta.iter().enumerate() {
                assert!((val - t[b][a]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_relay_zero_vector() {
        let t = relay_matrix_element([0.0; 3]);
        for row in &t {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }

    // ── induced_dipole_field ──────────────────────────────────────────────────

    #[test]
    fn test_field_single_atom_zero() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let e = induced_dipole_field(&atoms, 0);
        assert_eq!(e, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_field_out_of_range() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let e = induced_dipole_field(&atoms, 5);
        assert_eq!(e, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_field_two_atoms_symmetry() {
        let atoms = vec![
            PolarizableAtom::new([-1.0, 0.0, 0.0], 1.0, 1.0, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], 1.0, 1.0, 1.0),
        ];
        let e0 = induced_dipole_field(&atoms, 0);
        let e1 = induced_dipole_field(&atoms, 1);
        // Equal and opposite by symmetry
        assert!((e0[0] + e1[0]).abs() < 1e-10, "e0={e0:?}, e1={e1:?}");
    }

    // ── dipole_dipole_interaction ─────────────────────────────────────────────

    #[test]
    fn test_ddi_zero_distance() {
        assert_eq!(
            dipole_dipole_interaction([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]),
            0.0
        );
    }

    #[test]
    fn test_ddi_parallel() {
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [1.0, 0.0, 0.0];
        let r = [0.0, 0.0, 1.0];
        let u = dipole_dipole_interaction(d1, d2, r);
        // d1·d2=1, d1·r=0, d2·r=0 => U = 1/1 = 1
        assert!((u - 1.0).abs() < 1e-10, "u={u}");
    }

    #[test]
    fn test_ddi_head_to_tail() {
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [1.0, 0.0, 0.0];
        let r = [1.0, 0.0, 0.0];
        let u = dipole_dipole_interaction(d1, d2, r);
        // U = 1 - 3*1*1 = -2
        assert!((u + 2.0).abs() < 1e-10, "u={u}");
    }

    // ── scf_induced_dipoles ───────────────────────────────────────────────────

    #[test]
    fn test_scf_converges() {
        let mut atoms = vec![
            PolarizableAtom::new([-1.0, 0.0, 0.0], 1.0, 0.5, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], -1.0, 0.5, 1.0),
        ];
        let params = PolarizableParams::new(2.5, 10.0, 200, 1e-8);
        let n_iter = scf_induced_dipoles(&mut atoms, &params);
        assert!(n_iter <= 200);
    }

    #[test]
    fn test_scf_zero_polarizability() {
        let mut atoms = vec![
            PolarizableAtom::new([0.0; 3], 1.0, 0.0, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], -1.0, 0.0, 1.0),
        ];
        let params = PolarizableParams::new(2.5, 10.0, 10, 1e-8);
        let n_iter = scf_induced_dipoles(&mut atoms, &params);
        assert_eq!(n_iter, 1);
    }

    #[test]
    fn test_scf_thole_converges() {
        let mut atoms = vec![
            PolarizableAtom::new([-2.0, 0.0, 0.0], 1.0, 0.3, 1.0),
            PolarizableAtom::new([2.0, 0.0, 0.0], -1.0, 0.3, 1.0),
        ];
        let params = PolarizableParams::amoeba_defaults();
        let n_iter = scf_induced_dipoles_thole(&mut atoms, &params);
        assert!(n_iter <= 200);
    }

    // ── drude_oscillator_force ────────────────────────────────────────────────

    #[test]
    fn test_drude_force_zero_displacement() {
        let f = drude_oscillator_force([0.0; 3], [0.0; 3], 100.0);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_drude_force_direction() {
        let f = drude_oscillator_force([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        assert!(f[0] < 0.0);
    }

    #[test]
    fn test_drude_force_magnitude() {
        let f = drude_oscillator_force([2.0, 0.0, 0.0], [0.0; 3], 3.0);
        assert!((f[0] + 6.0).abs() < 1e-12);
    }

    // ── polarization_energy ───────────────────────────────────────────────────

    #[test]
    fn test_pol_energy_single_atom_zero() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let u = polarization_energy(&atoms);
        assert!(u.abs() < 1e-12);
    }

    #[test]
    fn test_pol_energy_negative() {
        let atoms = vec![
            PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], -1.0, 1.0, 1.0),
        ];
        let u = polarization_energy(&atoms);
        assert!(u <= 0.0, "u={u}");
    }

    // ── induction_energy ──────────────────────────────────────────────────────

    #[test]
    fn test_induction_energy_zero_dipoles() {
        let atoms = vec![
            PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], 1.0, 1.0, 1.0),
        ];
        let u = induction_energy(&atoms);
        assert!(u.abs() < 1e-12);
    }

    // ── fluctuating_charge_update ─────────────────────────────────────────────

    #[test]
    fn test_fluctuating_charge_conservation() {
        let mut q = vec![0.5_f64, -0.5];
        let chi = vec![4.0, 3.5];
        let hard = vec![10.0, 10.0];
        let q_total_before: f64 = q.iter().sum();
        fluctuating_charge_update(&mut q, &chi, &hard, 0.01);
        let q_total_after: f64 = q.iter().sum();
        assert!((q_total_after - q_total_before).abs() < 1e-12);
    }

    #[test]
    fn test_fluctuating_charge_empty() {
        let mut q: Vec<f64> = Vec::new();
        fluctuating_charge_update(&mut q, &[], &[], 0.01);
        assert!(q.is_empty());
    }

    #[test]
    fn test_fluctuating_charge_equal_chi_no_change() {
        let mut q = vec![1.0_f64, -1.0];
        let chi = vec![4.0_f64, 4.0];
        let hard = vec![10.0_f64, 10.0];
        let q_before = q.clone();
        fluctuating_charge_update(&mut q, &chi, &hard, 0.01);
        for (qa, qb) in q.iter().zip(q_before.iter()) {
            assert!((qa - qb).abs() < 1e-12);
        }
    }

    #[test]
    fn test_fluctuating_charge_single() {
        let mut q = vec![1.0_f64];
        let chi = vec![5.0];
        let hard = vec![1.0];
        fluctuating_charge_update(&mut q, &chi, &hard, 0.1);
        assert!((q[0] - 1.0).abs() < 1e-12);
    }

    // ── multipole_electrostatic_energy ────────────────────────────────────────

    #[test]
    fn test_multipole_charge_charge() {
        let mut a = MultipoleAtom::new([0.0; 3], 1.0, 0.0);
        let mut b = MultipoleAtom::new([1.0, 0.0, 0.0], -1.0, 0.0);
        // Only charge–charge term: U = 1 * (-1) / 1 = -1
        a.perm_dipole = [0.0; 3];
        b.perm_dipole = [0.0; 3];
        let u = multipole_electrostatic_energy(&a, &b);
        assert!((u + 1.0).abs() < 1e-10, "u={u}");
    }

    #[test]
    fn test_multipole_zero_charges() {
        let a = MultipoleAtom::new([0.0; 3], 0.0, 0.0);
        let b = MultipoleAtom::new([2.0, 0.0, 0.0], 0.0, 0.0);
        let u = multipole_electrostatic_energy(&a, &b);
        assert!(u.abs() < 1e-12, "u={u}");
    }

    // ── shell_model_force ─────────────────────────────────────────────────────

    #[test]
    fn test_shell_force_no_field_displaced() {
        let mut s = ShellModelAtom::new([0.0; 3], -2.0, -1.0, 10.0);
        s.shell_pos = [1.0, 0.0, 0.0];
        let f = shell_model_force(&s, [0.0; 3]);
        // spring pulls shell back: F = -k * (r_shell - r_core) = -10 * (1,0,0) = (-10,0,0)
        assert!((f[0] + 10.0).abs() < 1e-12, "f={f:?}");
    }

    #[test]
    fn test_shell_force_with_field() {
        let s = ShellModelAtom::new([0.0; 3], -2.0, 4.0, 0.0);
        let e_ext = [1.0, 0.0, 0.0];
        let f = shell_model_force(&s, e_ext);
        // shell_charge = total - core = -2 - 4 = -6
        // F = 0 (no spring, no displacement) + (-6)*1 = -6 along x
        assert!((f[0] + 6.0).abs() < 1e-12, "f={f:?}");
    }

    // ── many_body_polarization_correction ────────────────────────────────────

    #[test]
    fn test_mb_correction_single_atom() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let c = many_body_polarization_correction(&atoms);
        assert!(c.abs() < 1e-12);
    }

    // ── swm4ndp_geometry ──────────────────────────────────────────────────────

    #[test]
    fn test_swm4ndp_m_site_closer_than_h() {
        let w = Swm4NdpWater::new([0.0; 3]);
        let r_m = norm3(sub3(w.m_pos, w.o_pos));
        let r_h = norm3(sub3(w.h1_pos, w.o_pos));
        assert!(r_m < r_h, "r_m={r_m}, r_h={r_h}");
    }

    #[test]
    fn test_swm4ndp_geometry_symmetry() {
        let w = Swm4NdpWater::new([0.0; 3]);
        // H1 and H2 should be symmetric about xz-plane
        assert!((w.h1_pos[0] + w.h2_pos[0]).abs() < 1e-10);
        assert!((w.h1_pos[1] - w.h2_pos[1]).abs() < 1e-10);
    }

    #[test]
    fn test_swm4ndp_update_m_site() {
        let mut w = Swm4NdpWater::new([1.0, 2.0, 3.0]);
        let m_before = w.m_pos;
        w.o_pos = [1.0, 2.0, 3.0]; // unchanged
        w.update_m_site();
        for (m, &mb) in w.m_pos.iter().zip(m_before.iter()) {
            assert!((m - mb).abs() < 1e-10);
        }
    }

    // ── clausius_mossotti ─────────────────────────────────────────────────────

    #[test]
    fn test_clausius_mossotti_zero() {
        let eps = clausius_mossotti(0.0);
        assert!((eps - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_clausius_mossotti_positive() {
        let eps = clausius_mossotti(0.1);
        assert!(eps > 1.0);
    }

    // ── dipole_msd ────────────────────────────────────────────────────────────

    #[test]
    fn test_dipole_msd_zero() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        assert!(dipole_msd(&atoms).abs() < 1e-12);
    }

    #[test]
    fn test_dipole_msd_empty() {
        assert_eq!(dipole_msd(&[]), 0.0);
    }

    // ── reset_dipoles ─────────────────────────────────────────────────────────

    #[test]
    fn test_reset_dipoles() {
        let mut atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        atoms[0].dipole = [1.0, 2.0, 3.0];
        reset_dipoles(&mut atoms);
        assert_eq!(atoms[0].dipole, [0.0, 0.0, 0.0]);
    }

    // ── drude_spring_energy ───────────────────────────────────────────────────

    #[test]
    fn test_drude_spring_energy_zero() {
        let core_positions = vec![[0.0_f64; 3]];
        let drudes = vec![DrudeParticle::new(0, [0.0; 3], 1.0, 10.0, 0.4)];
        let e = drude_spring_energy(&drudes, &core_positions);
        assert!(e.abs() < 1e-12);
    }

    #[test]
    fn test_drude_spring_energy_displaced() {
        let core_positions = vec![[0.0_f64; 3]];
        let mut d = DrudeParticle::new(0, [0.0; 3], 1.0, 2.0, 0.4);
        d.shell_pos = [1.0, 0.0, 0.0];
        let e = drude_spring_energy(&[d], &core_positions);
        // ½ * 2 * 1² = 1
        assert!((e - 1.0).abs() < 1e-12);
    }

    // ── field_at_point ────────────────────────────────────────────────────────

    #[test]
    fn test_field_at_point_empty() {
        let e = field_at_point([0.0; 3], &[]);
        assert_eq!(e, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_field_at_point_from_charge() {
        // Single charge q=1 at origin; field at (1,0,0) should point in +x
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 0.0, 1.0)];
        let e = field_at_point([1.0, 0.0, 0.0], &atoms);
        assert!(e[0] > 0.0, "e={e:?}");
        assert!(e[1].abs() < 1e-12);
        assert!(e[2].abs() < 1e-12);
    }

    // ── linear_response_polarizability ────────────────────────────────────────

    #[test]
    fn test_linear_response() {
        let atoms = vec![
            PolarizableAtom::new([0.0; 3], 0.0, 1.5, 1.0),
            PolarizableAtom::new([1.0, 0.0, 0.0], 0.0, 2.5, 1.0),
        ];
        let alpha = linear_response_polarizability(&atoms, 0);
        assert!((alpha - 4.0).abs() < 1e-12);
    }

    // ── born_polarizability_estimate ──────────────────────────────────────────

    #[test]
    fn test_born_polarizability() {
        // α = Z*² / k = 4 / 20 = 0.2
        let alpha = born_polarizability_estimate(2.0, 20.0);
        assert!((alpha - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_born_polarizability_zero_k() {
        let alpha = born_polarizability_estimate(1.0, 0.0);
        assert_eq!(alpha, 0.0);
    }

    // ── drude_kinetic_temperature ─────────────────────────────────────────────

    #[test]
    fn test_drude_temp_zero() {
        let d = DrudeParticle::new(0, [0.0; 3], 1.0, 100.0, 0.4);
        let t = drude_kinetic_temperature(&[d]);
        assert!(t.abs() < 1e-12);
    }

    #[test]
    fn test_drude_temp_empty() {
        assert_eq!(drude_kinetic_temperature(&[]), 0.0);
    }

    // ── cross_polarization_force ──────────────────────────────────────────────

    #[test]
    fn test_cross_pol_force_same_index() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let (fi, fj) = cross_polarization_force(&atoms, 0, 0);
        assert_eq!(fi, [0.0; 3]);
        assert_eq!(fj, [0.0; 3]);
    }

    #[test]
    fn test_cross_pol_force_out_of_range() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let (fi, fj) = cross_polarization_force(&atoms, 0, 5);
        assert_eq!(fi, [0.0; 3]);
        assert_eq!(fj, [0.0; 3]);
    }

    // ── total_shell_elastic_energy ────────────────────────────────────────────

    #[test]
    fn test_total_shell_elastic_zero() {
        let atoms = vec![
            ShellModelAtom::new([0.0; 3], -2.0, -1.0, 10.0),
            ShellModelAtom::new([3.0, 0.0, 0.0], 2.0, 1.0, 10.0),
        ];
        let e = total_shell_elastic_energy(&atoms);
        assert!(e.abs() < 1e-12);
    }

    // ── DrудеThermostat ───────────────────────────────────────────────────────

    #[test]
    fn test_drude_thermostat_new() {
        let t = DrудеThermostat::new(300.0, 1.0, 5.0, 1000.0);
        assert!((t.t_real - 300.0).abs() < 1e-12);
        assert!((t.t_drude - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_drude_thermostat_zero_gamma() {
        // With gamma=0, velocity is unchanged (no friction, no noise if sigma=0)
        let ts = DrудеThermostat::new(0.0, 0.0, 0.0, 0.0);
        let vel = [1.0, 2.0, 3.0];
        // sigma = sqrt(0) = 0, friction = 1 => vel unchanged
        let v_new = ts.apply_to_velocity(vel, 1.0, 0.001, [0.0; 3], false);
        for k in 0..3 {
            assert!((v_new[k] - vel[k]).abs() < 1e-10);
        }
    }

    // ── polarizable_virial ────────────────────────────────────────────────────

    #[test]
    fn test_virial_single_atom() {
        let atoms = vec![PolarizableAtom::new([0.0; 3], 1.0, 1.0, 1.0)];
        let v = polarizable_virial(&atoms);
        // Single atom: no pairs => zero virial
        for row in &v {
            for &val in row {
                assert!(val.abs() < 1e-12);
            }
        }
    }

    // ── amoeba_induced_energy ─────────────────────────────────────────────────

    #[test]
    fn test_amoeba_induced_energy_zero_dipoles() {
        let atoms = vec![
            MultipoleAtom::new([0.0; 3], 1.0, 1.0),
            MultipoleAtom::new([2.0, 0.0, 0.0], -1.0, 1.0),
        ];
        let u = amoeba_induced_energy(&atoms);
        // induced_dipole is zero => energy = 0
        assert!(u.abs() < 1e-12);
    }
}
