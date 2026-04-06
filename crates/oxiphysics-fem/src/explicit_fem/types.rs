//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// Failure criterion for element erosion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FailureCriterion {
    /// Delete element when effective plastic strain exceeds threshold.
    PlasticStrain(f64),
    /// Delete element when maximum principal stress exceeds threshold.
    PrincipalStress(f64),
    /// Delete element when damage variable exceeds threshold.
    Damage(f64),
    /// Delete element when element volume drops below fraction of initial.
    VolumeFraction(f64),
}
/// Mass lumping strategy for the explicit solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MassLumping {
    /// Row-sum lumping: each diagonal entry is the sum of its row.
    RowSum,
    /// HRZ lumping (Hinton–Rock–Zienkiewicz): scale diagonal to preserve total mass.
    Hrz,
    /// Diagonal scaling based on nodal volume.
    NodalVolume,
}
/// Type of hourglass stabilization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HourglassType {
    /// Stiffness-based hourglass control (Flanagan–Belytschko).
    Stiffness,
    /// Viscosity-based hourglass control.
    Viscous,
    /// Combined stiffness + viscous.
    Combined,
}
/// Contact detection and force computation for FEM.
pub struct ContactDetectionFem {
    /// Enforcement method.
    pub enforcement: ContactEnforcement,
    /// Penalty stiffness (N/m³).
    pub penalty_stiffness: f64,
    /// Augmentation parameter for augmented Lagrangian.
    pub augmentation: f64,
    /// Coulomb friction coefficient.
    pub friction_coeff: f64,
    /// Active contact pairs.
    pub pairs: Vec<ContactPair>,
    /// Search radius for contact detection.
    pub search_radius: f64,
}
impl ContactDetectionFem {
    /// Create a new contact detection system.
    pub fn new(
        enforcement: ContactEnforcement,
        penalty_stiffness: f64,
        friction_coeff: f64,
    ) -> Self {
        Self {
            enforcement,
            penalty_stiffness,
            augmentation: penalty_stiffness,
            friction_coeff,
            pairs: Vec::new(),
            search_radius: 1e-3,
        }
    }
    /// Compute the gap between a slave node and a master segment in 2D.
    ///
    /// Returns (gap, parameter_t, normal) where `t` is the projection parameter.
    pub fn node_to_segment_gap_2d(
        slave: [f64; 2],
        master_a: [f64; 2],
        master_b: [f64; 2],
    ) -> (f64, f64, [f64; 2]) {
        let ab = [master_b[0] - master_a[0], master_b[1] - master_a[1]];
        let len2 = ab[0] * ab[0] + ab[1] * ab[1];
        let t = if len2 > 1e-30 {
            let as_ = [slave[0] - master_a[0], slave[1] - master_a[1]];
            (as_[0] * ab[0] + as_[1] * ab[1]) / len2
        } else {
            0.5
        };
        let t_clamp = t.clamp(0.0, 1.0);
        let proj = [master_a[0] + t_clamp * ab[0], master_a[1] + t_clamp * ab[1]];
        let diff = [slave[0] - proj[0], slave[1] - proj[1]];
        let _dist = (diff[0] * diff[0] + diff[1] * diff[1]).sqrt();
        let seg_len = len2.sqrt().max(1e-30);
        let normal = [-ab[1] / seg_len, ab[0] / seg_len];
        let gap = diff[0] * normal[0] + diff[1] * normal[1];
        (gap, t_clamp, normal)
    }
    /// Compute the gap between a slave node and a master segment in 3D.
    ///
    /// Returns `(gap, t, normal)` where:
    /// - `gap > 0` means separation (no contact)
    /// - `gap < 0` means penetration
    /// - `t` is the parameter of the closest point on the segment
    /// - `normal` points from closest point toward the slave node
    pub fn node_to_segment_gap_3d(
        slave: [f64; 3],
        master_a: [f64; 3],
        master_b: [f64; 3],
    ) -> (f64, f64, [f64; 3]) {
        let ab = [
            master_b[0] - master_a[0],
            master_b[1] - master_a[1],
            master_b[2] - master_a[2],
        ];
        let len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        let t = if len2 > 1e-30 {
            let as_ = [
                slave[0] - master_a[0],
                slave[1] - master_a[1],
                slave[2] - master_a[2],
            ];
            (as_[0] * ab[0] + as_[1] * ab[1] + as_[2] * ab[2]) / len2
        } else {
            0.5
        };
        let t_clamp = t.clamp(0.0, 1.0);
        let proj = [
            master_a[0] + t_clamp * ab[0],
            master_a[1] + t_clamp * ab[1],
            master_a[2] + t_clamp * ab[2],
        ];
        let diff = [slave[0] - proj[0], slave[1] - proj[1], slave[2] - proj[2]];
        let dist = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
        let normal = if dist > 1e-30 {
            [diff[0] / dist, diff[1] / dist, diff[2] / dist]
        } else {
            [0.0, 1.0, 0.0]
        };
        (dist, t_clamp, normal)
    }
    /// Compute the penalty contact force for a given gap.
    pub fn penalty_force(&self, gap: f64) -> f64 {
        if gap < 0.0 {
            -self.penalty_stiffness * gap
        } else {
            0.0
        }
    }
    /// Compute the augmented Lagrangian contact force.
    pub fn augmented_lagrangian_force(&self, gap: f64, lambda: f64) -> f64 {
        let trial = lambda + self.augmentation * gap;
        if trial < 0.0 { -trial } else { 0.0 }
    }
    /// Detect contacts using a node-to-segment search.
    ///
    /// `node_positions` — flat array of node positions (3 per node)
    /// `segments` — list of segment node pairs
    /// `slave_nodes` — indices of slave nodes to check
    pub fn detect_node_to_segment(
        &mut self,
        node_positions: &[f64],
        segments: &[[usize; 2]],
        slave_nodes: &[usize],
    ) {
        self.pairs.clear();
        for &sn in slave_nodes {
            let sp = [
                node_positions[sn * 3],
                node_positions[sn * 3 + 1],
                node_positions[sn * 3 + 2],
            ];
            for &seg in segments {
                let ma = [
                    node_positions[seg[0] * 3],
                    node_positions[seg[0] * 3 + 1],
                    node_positions[seg[0] * 3 + 2],
                ];
                let mb = [
                    node_positions[seg[1] * 3],
                    node_positions[seg[1] * 3 + 1],
                    node_positions[seg[1] * 3 + 2],
                ];
                let (gap, _t, normal) = Self::node_to_segment_gap_3d(sp, ma, mb);
                if gap < self.search_radius {
                    let mut pair = ContactPair::new(sn, seg);
                    pair.gap = gap;
                    pair.normal = normal;
                    self.pairs.push(pair);
                }
            }
        }
    }
    /// Assemble contact forces into the global force vector.
    pub fn assemble_contact_forces(&self, force: &mut [f64]) {
        for pair in &self.pairs {
            let f = match self.enforcement {
                ContactEnforcement::Penalty => self.penalty_force(pair.gap),
                ContactEnforcement::Lagrange => {
                    if pair.gap < 0.0 {
                        pair.lambda
                    } else {
                        0.0
                    }
                }
                ContactEnforcement::AugmentedLagrangian => {
                    self.augmented_lagrangian_force(pair.gap, pair.lambda)
                }
            };
            let sn = pair.slave_node;
            for dir in 0..3 {
                force[sn * 3 + dir] += f * pair.normal[dir];
            }
            for (k, mn) in pair.master_segment.iter().enumerate() {
                let _ = k;
                for dir in 0..3 {
                    force[mn * 3 + dir] -= 0.5 * f * pair.normal[dir];
                }
            }
        }
    }
    /// Update Lagrange multipliers for augmented Lagrangian iteration.
    pub fn update_lagrange_multipliers(&mut self) {
        for pair in &mut self.pairs {
            if pair.gap < 0.0 {
                pair.lambda += self.augmentation * pair.gap;
                pair.lambda = pair.lambda.max(0.0);
            }
        }
    }
}
/// Stress wave propagation in bars and plates.
///
/// Implements analytical and FEM-based wave propagation analysis
/// including dispersion curves, reflection/transmission at interfaces.
pub struct WavePropagationFem {
    /// Wave type being analyzed.
    pub wave_type: WaveType,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Shear modulus (Pa).
    pub shear_modulus: f64,
    /// Material density (kg/m³).
    pub density: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Bar/plate cross-sectional area or thickness (m or m²).
    pub characteristic_size: f64,
}
impl WavePropagationFem {
    /// Create a new wave propagation analyzer.
    pub fn new(
        wave_type: WaveType,
        youngs_modulus: f64,
        density: f64,
        poisson_ratio: f64,
        characteristic_size: f64,
    ) -> Self {
        let shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
        Self {
            wave_type,
            youngs_modulus,
            shear_modulus,
            density,
            poisson_ratio,
            characteristic_size,
        }
    }
    /// Dilatational (P-wave) speed in a 3D solid: `c_p = sqrt((K + 4G/3) / rho)`.
    pub fn p_wave_speed(&self) -> f64 {
        let k = self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio));
        let g = self.shear_modulus;
        ((k + 4.0 * g / 3.0) / self.density).sqrt()
    }
    /// Shear (S-wave) speed: `c_s = sqrt(G / rho)`.
    pub fn s_wave_speed(&self) -> f64 {
        (self.shear_modulus / self.density).sqrt()
    }
    /// Bar wave speed (1D): `c_bar = sqrt(E / rho)`.
    pub fn bar_wave_speed(&self) -> f64 {
        (self.youngs_modulus / self.density).sqrt()
    }
    /// Rayleigh wave speed (approximate): `c_R ≈ c_s * (0.862 + 1.14*nu) / (1 + nu)`.
    pub fn rayleigh_wave_speed(&self) -> f64 {
        let nu = self.poisson_ratio;
        let cs = self.s_wave_speed();
        cs * (0.862 + 1.14 * nu) / (1.0 + nu)
    }
    /// Compute dispersion curves for a range of frequencies.
    ///
    /// For longitudinal waves in a bar using Love (lateral inertia) correction:
    /// `c_ph = c_bar / sqrt(1 + nu^2 * k^2 * r^2)`
    ///
    /// where `r` is the cross-sectional radius of gyration.
    pub fn dispersion_curve(&self, frequencies: &[f64]) -> WaveResult {
        let c0 = match self.wave_type {
            WaveType::Longitudinal => self.bar_wave_speed(),
            WaveType::Transverse => self.s_wave_speed(),
            WaveType::Flexural => self.bar_wave_speed(),
            WaveType::Rayleigh => self.rayleigh_wave_speed(),
        };
        let nu = self.poisson_ratio;
        let r = self.characteristic_size;
        let n = frequencies.len();
        let mut wave_numbers = vec![0.0; n];
        let mut phase_velocities = vec![0.0; n];
        let mut group_velocities = vec![0.0; n];
        let mut dispersion_param = vec![0.0; n];
        for (i, &f) in frequencies.iter().enumerate() {
            let omega = 2.0 * std::f64::consts::PI * f;
            let k = match self.wave_type {
                WaveType::Longitudinal => {
                    let k0 = omega / c0;
                    let denom = (1.0 + nu * nu * k0 * k0 * r * r).sqrt();
                    omega / (c0 * denom)
                }
                WaveType::Transverse => omega / c0,
                WaveType::Flexural => {
                    let ei = self.youngs_modulus * r * r * r * r / 12.0;
                    (omega * omega * self.density / ei).sqrt().sqrt()
                }
                WaveType::Rayleigh => omega / c0,
            };
            let cp = if k > 1e-30 { omega / k } else { c0 };
            let cg = match self.wave_type {
                WaveType::Flexural => 2.0 * cp,
                _ => cp,
            };
            wave_numbers[i] = k;
            phase_velocities[i] = cp;
            group_velocities[i] = cg;
            dispersion_param[i] = (cp - c0).abs() / c0;
        }
        WaveResult {
            frequencies: frequencies.to_vec(),
            phase_velocities,
            group_velocities,
            wave_numbers,
            dispersion_parameter: dispersion_param,
        }
    }
    /// Compute reflection and transmission coefficients at an interface.
    ///
    /// For a normal-incidence wave at a boundary between material 1 (current)
    /// and material 2 (given by `z2`):
    ///
    /// `R = (Z2 - Z1) / (Z2 + Z1)`, `T = 2*Z2 / (Z2 + Z1)`
    ///
    /// where `Z = rho * c` is the acoustic impedance.
    pub fn reflection_transmission(&self, density2: f64, wave_speed2: f64) -> (f64, f64) {
        let z1 = self.density * self.bar_wave_speed();
        let z2 = density2 * wave_speed2;
        if (z1 + z2).abs() < 1e-30 {
            return (0.0, 1.0);
        }
        let r = (z2 - z1) / (z2 + z1);
        let t = 2.0 * z2 / (z2 + z1);
        (r, t)
    }
    /// Compute the stress wave arrival time for a given distance.
    pub fn arrival_time(&self, distance: f64) -> f64 {
        distance / self.bar_wave_speed()
    }
    /// Compute the Hopkinson bar stress pulse shape (rectangular pulse).
    ///
    /// Returns pressure at given time for a rectangular pulse of duration `duration`.
    pub fn hopkinson_pulse(&self, time: f64, amplitude: f64, duration: f64) -> f64 {
        if time >= 0.0 && time <= duration {
            amplitude
        } else {
            0.0
        }
    }
    /// Compute FEM numerical wave speed dispersion.
    ///
    /// For a 1D uniform mesh with element size `h`, the FEM disperses waves.
    /// `omega * h / c` is the non-dimensional frequency.
    pub fn fem_dispersion(element_size: f64, wave_speed: f64, frequency: f64) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * frequency;
        let xi = omega * element_size / wave_speed;
        if xi <= 2.0 {
            let k_fem = (2.0 / element_size) * (xi / 2.0).asin();
            let k_exact = omega / wave_speed;
            if k_exact > 1e-30 {
                k_fem / k_exact
            } else {
                1.0
            }
        } else {
            f64::NAN
        }
    }
}
/// Blast loading model type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BlastModel {
    /// User-defined pressure-time history (tabulated).
    UserDefined,
    /// Friedlander waveform: `p(t) = p0 * (1 - t/td) * exp(-b*t/td)`.
    Friedlander,
    /// ConWep empirical model (Baker–Westine–Dodge).
    ConWep,
    /// Simplified triangular pulse.
    Triangular,
}
/// Hourglass mode shapes for a single hexahedral element.
///
/// The standard 8-node hex has 12 hourglass modes (4 independent hourglass
/// patterns × 3 displacement components).
#[derive(Debug, Clone)]
pub struct HourglassControl {
    /// Type of stabilization.
    pub hg_type: HourglassType,
    /// Hourglass stiffness coefficient (fraction of element stiffness).
    pub stiffness_coeff: f64,
    /// Hourglass viscosity coefficient.
    pub viscosity_coeff: f64,
    /// Number of elements.
    pub n_elements: usize,
    /// Accumulated hourglass energy per element.
    pub hourglass_energy: Vec<f64>,
}
impl HourglassControl {
    /// Standard hourglass basis vectors for a unit hexahedron (Flanagan–Belytschko).
    /// Returns a 4×8 matrix of hourglass base vectors (rows = modes, cols = nodes).
    pub const HOURGLASS_GAMMA: [[f64; 8]; 4] = [
        [1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0],
        [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0],
    ];
    /// Create a new hourglass controller.
    pub fn new(
        hg_type: HourglassType,
        stiffness_coeff: f64,
        viscosity_coeff: f64,
        n_elements: usize,
    ) -> Self {
        Self {
            hg_type,
            stiffness_coeff,
            viscosity_coeff,
            n_elements,
            hourglass_energy: vec![0.0; n_elements],
        }
    }
    /// Compute the hourglass mode amplitudes for a single element.
    ///
    /// `node_displacements` is an 8×3 array of nodal displacements.
    /// Returns a 4×3 matrix of mode amplitudes.
    pub fn compute_mode_amplitudes(&self, node_displacements: &[[f64; 3]; 8]) -> [[f64; 3]; 4] {
        let mut amplitudes = [[0.0f64; 3]; 4];
        for (mode, gamma) in Self::HOURGLASS_GAMMA.iter().enumerate() {
            for dir in 0..3 {
                amplitudes[mode][dir] = gamma
                    .iter()
                    .zip(node_displacements.iter())
                    .map(|(&g, nd)| g * nd[dir])
                    .sum();
            }
        }
        amplitudes
    }
    /// Compute stiffness-based hourglass forces for a single element.
    ///
    /// `q_modes` — 4×3 mode amplitudes
    /// `element_modulus` — Young's modulus
    /// `element_volume` — element volume
    ///
    /// Returns nodal hourglass forces (8×3).
    pub fn stiffness_forces(
        &self,
        q_modes: &[[f64; 3]; 4],
        element_modulus: f64,
        element_volume: f64,
    ) -> [[f64; 3]; 8] {
        let k_hg = self.stiffness_coeff * element_modulus * element_volume.cbrt();
        let mut forces = [[0.0f64; 3]; 8];
        for (mode, gamma) in Self::HOURGLASS_GAMMA.iter().enumerate() {
            for dir in 0..3 {
                let f_mode = -k_hg * q_modes[mode][dir];
                for (node, &g) in gamma.iter().enumerate() {
                    forces[node][dir] += f_mode * g;
                }
            }
        }
        forces
    }
    /// Compute viscosity-based hourglass forces for a single element.
    ///
    /// `q_dot_modes` — 4×3 mode velocity amplitudes
    /// `density` — material density
    /// `wave_speed` — dilatational wave speed
    /// `element_volume` — element volume
    pub fn viscous_forces(
        &self,
        q_dot_modes: &[[f64; 3]; 4],
        density: f64,
        wave_speed: f64,
        element_volume: f64,
    ) -> [[f64; 3]; 8] {
        let c_hg = self.viscosity_coeff * density * wave_speed * element_volume.powf(1.0 / 3.0);
        let mut forces = [[0.0f64; 3]; 8];
        for (mode, gamma) in Self::HOURGLASS_GAMMA.iter().enumerate() {
            for dir in 0..3 {
                let f_mode = -c_hg * q_dot_modes[mode][dir];
                for (node, &g) in gamma.iter().enumerate() {
                    forces[node][dir] += f_mode * g;
                }
            }
        }
        forces
    }
    /// Compute combined hourglass forces.
    pub fn combined_forces(
        &self,
        q_modes: &[[f64; 3]; 4],
        q_dot_modes: &[[f64; 3]; 4],
        element_modulus: f64,
        density: f64,
        wave_speed: f64,
        element_volume: f64,
    ) -> [[f64; 3]; 8] {
        let f_stiff = self.stiffness_forces(q_modes, element_modulus, element_volume);
        let f_visc = self.viscous_forces(q_dot_modes, density, wave_speed, element_volume);
        let mut total = [[0.0f64; 3]; 8];
        for node in 0..8 {
            for dir in 0..3 {
                total[node][dir] = f_stiff[node][dir] + f_visc[node][dir];
            }
        }
        total
    }
    /// Compute and accumulate hourglass energy for an element.
    ///
    /// `element_idx` — index into `hourglass_energy`
    pub fn update_hourglass_energy(
        &mut self,
        element_idx: usize,
        q_modes: &[[f64; 3]; 4],
        k_hg: f64,
    ) {
        let energy: f64 = q_modes
            .iter()
            .flat_map(|q| q.iter())
            .map(|q| 0.5 * k_hg * q * q)
            .sum();
        self.hourglass_energy[element_idx] = energy;
    }
    /// Return total hourglass energy across all elements.
    pub fn total_hourglass_energy(&self) -> f64 {
        self.hourglass_energy.iter().sum()
    }
    /// Check whether hourglass energy is acceptable (< fraction of total strain energy).
    pub fn is_hourglass_acceptable(&self, total_strain_energy: f64, threshold: f64) -> bool {
        let hg_energy = self.total_hourglass_energy();
        if total_strain_energy > 1e-30 {
            hg_energy / total_strain_energy < threshold
        } else {
            true
        }
    }
}
/// Wave type for dispersion analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WaveType {
    /// Longitudinal (P-wave) in a bar.
    Longitudinal,
    /// Transverse (S-wave) in a bar.
    Transverse,
    /// Flexural wave in a beam.
    Flexural,
    /// Rayleigh surface wave.
    Rayleigh,
}
/// Contact enforcement method.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactEnforcement {
    /// Penalty method: `f = k_pen * gap` when gap < 0.
    Penalty,
    /// Lagrange multiplier method: exact enforcement.
    Lagrange,
    /// Augmented Lagrangian: combines accuracy and robustness.
    AugmentedLagrangian,
}
/// Explicit FEM solver using the central difference time integration scheme.
///
/// The central difference method is conditionally stable; the time step must
/// satisfy the CFL (Courant–Friedrichs–Lewy) condition:
///
/// `dt ≤ cfl * L_min / c`
///
/// where `L_min` is the smallest element length and `c` is the dilatational
/// wave speed `sqrt((K + 4G/3) / rho)`.
pub struct ExplicitFemSolver {
    /// Solver configuration.
    pub config: ExplicitFemConfig,
    /// Lumped mass vector (one entry per DOF).
    pub mass: Vec<f64>,
    /// External force vector.
    pub external_force: Vec<f64>,
    /// Internal force vector (assembled from element contributions).
    pub internal_force: Vec<f64>,
    /// Damping vector (for viscous damping `c * v`).
    pub damping: Vec<f64>,
    /// Current solver state.
    pub state: ExplicitState,
    /// Boundary condition mask: `true` means the DOF is fixed.
    pub fixed_dof: Vec<bool>,
    /// Energy history: (time, kinetic, internal, total).
    pub energy_history: Vec<(f64, f64, f64, f64)>,
}
impl ExplicitFemSolver {
    /// Create a new explicit FEM solver.
    ///
    /// # Arguments
    /// * `n_dof`  — total number of degrees of freedom
    /// * `config` — solver configuration
    pub fn new(n_dof: usize, config: ExplicitFemConfig) -> Self {
        Self {
            config,
            mass: vec![1.0; n_dof],
            external_force: vec![0.0; n_dof],
            internal_force: vec![0.0; n_dof],
            damping: vec![0.0; n_dof],
            state: ExplicitState::new(n_dof),
            fixed_dof: vec![false; n_dof],
            energy_history: Vec::new(),
        }
    }
    /// Compute the critical time step from the CFL condition.
    ///
    /// `dt_crit = cfl * L_min / c`
    pub fn critical_time_step(&self) -> f64 {
        self.config.cfl_factor * self.config.min_element_length / self.config.wave_speed
    }
    /// Set the lumped mass for all DOFs from element masses.
    ///
    /// `element_masses` is a list of `(dof_indices, element_mass)` pairs.
    pub fn assemble_lumped_mass(&mut self, element_masses: &[(Vec<usize>, f64)]) {
        for m in self.mass.iter_mut() {
            *m = 0.0;
        }
        for (dofs, em) in element_masses {
            let share = em / dofs.len() as f64;
            for &d in dofs {
                self.mass[d] += share;
            }
        }
    }
    /// Apply HRZ lumping correction to preserve total element mass.
    ///
    /// Scales the diagonal so that the total mass matches the consistent mass trace.
    pub fn apply_hrz_correction(&mut self, element_dofs: &[Vec<usize>], consistent_traces: &[f64]) {
        for (dofs, &trace) in element_dofs.iter().zip(consistent_traces.iter()) {
            let diag_sum: f64 = dofs.iter().map(|&d| self.mass[d]).sum();
            if diag_sum > 1e-30 {
                let scale = trace / diag_sum;
                for &d in dofs {
                    self.mass[d] *= scale;
                }
            }
        }
    }
    /// Perform one explicit central difference step.
    ///
    /// The algorithm:
    /// 1. Half-step velocity: `v_{n+1/2} = v_n + dt/2 * a_n`
    /// 2. Full displacement:   `u_{n+1} = u_n + dt * v_{n+1/2}`
    /// 3. Compute internal forces from `u_{n+1}`
    /// 4. New acceleration:    `a_{n+1} = M^{-1} (F_ext - F_int - C * v_{n+1/2})`
    /// 5. Full velocity:       `v_{n+1} = v_{n+1/2} + dt/2 * a_{n+1}`
    pub fn step(&mut self, dt: f64, internal_force_fn: impl Fn(&[f64]) -> Vec<f64>) {
        let n = self.state.n_dof();
        let mut vel_half = vec![0.0; n];
        for i in 0..n {
            if !self.fixed_dof[i] {
                vel_half[i] = self.state.velocity[i] + 0.5 * dt * self.state.acceleration[i];
            }
        }
        for i in 0..n {
            if !self.fixed_dof[i] {
                self.state.displacement[i] += dt * vel_half[i];
            }
        }
        self.internal_force = internal_force_fn(&self.state.displacement);
        for i in 0..n {
            if !self.fixed_dof[i] && self.mass[i] > 1e-30 {
                let f_net =
                    self.external_force[i] - self.internal_force[i] - self.damping[i] * vel_half[i];
                self.state.acceleration[i] = f_net / self.mass[i];
            } else {
                self.state.acceleration[i] = 0.0;
            }
        }
        for i in 0..n {
            if !self.fixed_dof[i] {
                self.state.velocity[i] = vel_half[i] + 0.5 * dt * self.state.acceleration[i];
            }
        }
        self.state.time += dt;
        self.state.step += 1;
        if self.config.check_energy {
            self.record_energy();
        }
    }
    /// Record energy at the current state.
    fn record_energy(&mut self) {
        let ke: f64 = self
            .state
            .velocity
            .iter()
            .zip(self.mass.iter())
            .map(|(v, m)| 0.5 * m * v * v)
            .sum();
        let ie: f64 = self
            .internal_force
            .iter()
            .zip(self.state.displacement.iter())
            .map(|(f, u)| f * u)
            .sum::<f64>()
            .abs()
            * 0.5;
        self.energy_history.push((self.state.time, ke, ie, ke + ie));
    }
    /// Compute kinetic energy at current state.
    pub fn kinetic_energy(&self) -> f64 {
        self.state
            .velocity
            .iter()
            .zip(self.mass.iter())
            .map(|(v, m)| 0.5 * m * v * v)
            .sum()
    }
    /// Fix a set of DOFs (Dirichlet boundary condition).
    pub fn fix_dofs(&mut self, dofs: &[usize]) {
        for &d in dofs {
            self.fixed_dof[d] = true;
            self.state.displacement[d] = 0.0;
            self.state.velocity[d] = 0.0;
            self.state.acceleration[d] = 0.0;
        }
    }
    /// Set initial velocity for specific DOFs.
    pub fn set_initial_velocity(&mut self, dof: usize, value: f64) {
        self.state.velocity[dof] = value;
    }
    /// Set external force on a DOF.
    pub fn set_external_force(&mut self, dof: usize, value: f64) {
        self.external_force[dof] = value;
    }
    /// Estimate the stable time step from element stiffness and mass.
    ///
    /// Uses the element eigenvalue estimate `omega_max^2 ≈ 2*k_max/m_min`.
    pub fn estimate_stable_dt_from_stiffness(
        &self,
        element_stiffness_max: f64,
        element_mass_min: f64,
    ) -> f64 {
        if element_stiffness_max > 0.0 && element_mass_min > 0.0 {
            self.config.cfl_factor * (2.0 * element_mass_min / element_stiffness_max).sqrt()
        } else {
            self.critical_time_step()
        }
    }
    /// Run multiple steps until `end_time`.
    pub fn run_until(&mut self, end_time: f64, internal_force_fn: impl Fn(&[f64]) -> Vec<f64>) {
        let dt = self.critical_time_step();
        while self.state.time < end_time {
            let dt_use = dt.min(end_time - self.state.time);
            self.step(dt_use, &internal_force_fn);
        }
    }
}
/// Blast loading FEM solver.
///
/// Applies time-varying pressure loads to surface nodes based on blast
/// wave models including Friedlander waveform and simplified ConWep.
pub struct BlastLoadingFem {
    /// Blast model type.
    pub model: BlastModel,
    /// Blast parameters.
    pub params: BlastParameters,
    /// Surface node indices.
    pub surface_nodes: Vec<usize>,
    /// Outward normals for each surface node (3 per node).
    pub node_normals: Vec<[f64; 3]>,
    /// Tributary areas for each surface node (m²).
    pub tributary_areas: Vec<f64>,
    /// Pressure history: (time, pressure) recorded at each step.
    pub pressure_history: Vec<(f64, f64)>,
    /// Structural response history: (time, max_displacement).
    pub response_history: Vec<(f64, f64)>,
}
impl BlastLoadingFem {
    /// Create a new blast loading system.
    pub fn new(model: BlastModel, params: BlastParameters) -> Self {
        Self {
            model,
            params,
            surface_nodes: Vec::new(),
            node_normals: Vec::new(),
            tributary_areas: Vec::new(),
            pressure_history: Vec::new(),
            response_history: Vec::new(),
        }
    }
    /// Evaluate the overpressure at time `t`.
    pub fn pressure_at(&self, t: f64) -> f64 {
        match self.model {
            BlastModel::Friedlander => {
                let p0 = self.params.peak_pressure;
                let td = self.params.duration;
                let b = self.params.decay_coeff;
                if t < 0.0 || t > td {
                    0.0
                } else {
                    p0 * (1.0 - t / td) * (-b * t / td).exp()
                }
            }
            BlastModel::Triangular => {
                let p0 = self.params.peak_pressure;
                let td = self.params.duration;
                if t < 0.0 || t > td {
                    0.0
                } else {
                    p0 * (1.0 - t / td)
                }
            }
            BlastModel::ConWep => {
                let r = self.params.standoff;
                let w = self.params.charge_mass;
                let z = r / w.cbrt();
                let p_mpa = if z < 0.5 {
                    67.36 / (z * z * z)
                } else {
                    0.84 + 3.35 / z + 0.66 / (z * z)
                };
                let p0 = p_mpa * 1e6;
                let td = self.params.duration;
                if t < 0.0 || t > td {
                    0.0
                } else {
                    p0 * (1.0 - t / td) * (-b_coeff(z) * t / td).exp()
                }
            }
            BlastModel::UserDefined => interpolate_table(&self.params.table, t),
        }
    }
    /// Assemble blast pressure forces into the global force vector.
    ///
    /// Each surface node receives `p * A * normal` force.
    pub fn assemble_blast_forces(&self, time: f64, force: &mut [f64]) {
        let p = self.pressure_at(time);
        for (idx, &node) in self.surface_nodes.iter().enumerate() {
            let area = if idx < self.tributary_areas.len() {
                self.tributary_areas[idx]
            } else {
                1.0
            };
            let normal = if idx < self.node_normals.len() {
                self.node_normals[idx]
            } else {
                [0.0, 1.0, 0.0]
            };
            for dir in 0..3 {
                force[node * 3 + dir] += p * area * normal[dir];
            }
        }
    }
    /// Record current pressure and structural response.
    pub fn record_state(&mut self, time: f64, max_displacement: f64) {
        let p = self.pressure_at(time);
        self.pressure_history.push((time, p));
        self.response_history.push((time, max_displacement));
    }
    /// Compute the impulse (integral of pressure over time) numerically.
    pub fn compute_impulse(&self) -> f64 {
        if self.pressure_history.len() < 2 {
            return 0.0;
        }
        self.pressure_history
            .windows(2)
            .map(|w| {
                let dt = w[1].0 - w[0].0;
                0.5 * (w[0].1 + w[1].1) * dt
            })
            .sum()
    }
    /// Compute the dynamic load factor (DLF) from the response history.
    ///
    /// DLF = max_dynamic_displacement / static_displacement
    pub fn dynamic_load_factor(&self, static_displacement: f64) -> f64 {
        let max_disp = self
            .response_history
            .iter()
            .map(|(_, d)| d.abs())
            .fold(0.0_f64, f64::max);
        if static_displacement.abs() > 1e-30 {
            max_disp / static_displacement.abs()
        } else {
            0.0
        }
    }
    /// Compute the peak reflected pressure using the reflection factor.
    ///
    /// `pr = cr * pi` where `cr ≈ 2 + (γ+1)/2 * (p/patm)` for strong shocks.
    pub fn reflected_pressure(&self, time: f64, p_atm: f64) -> f64 {
        let p_i = self.pressure_at(time);
        let gamma = 1.4;
        let cr = 2.0 + (gamma + 1.0) / 2.0 * (p_i / p_atm);
        cr * p_i
    }
}
/// State vector for explicit time integration.
#[derive(Debug, Clone)]
pub struct ExplicitState {
    /// Nodal displacements (3 DOF per node, flattened).
    pub displacement: Vec<f64>,
    /// Nodal velocities (3 DOF per node, flattened).
    pub velocity: Vec<f64>,
    /// Nodal accelerations (3 DOF per node, flattened).
    pub acceleration: Vec<f64>,
    /// Current simulation time.
    pub time: f64,
    /// Current time step index.
    pub step: u64,
}
impl ExplicitState {
    /// Create a zero-initialized state for `n_dof` degrees of freedom.
    pub fn new(n_dof: usize) -> Self {
        Self {
            displacement: vec![0.0; n_dof],
            velocity: vec![0.0; n_dof],
            acceleration: vec![0.0; n_dof],
            time: 0.0,
            step: 0,
        }
    }
    /// Returns the number of degrees of freedom.
    pub fn n_dof(&self) -> usize {
        self.displacement.len()
    }
}
/// Blast loading parameters.
#[derive(Debug, Clone)]
pub struct BlastParameters {
    /// Peak overpressure (Pa).
    pub peak_pressure: f64,
    /// Positive phase duration (s).
    pub duration: f64,
    /// Friedlander decay coefficient.
    pub decay_coeff: f64,
    /// Standoff distance (m).
    pub standoff: f64,
    /// Charge mass equivalent TNT (kg).
    pub charge_mass: f64,
    /// User-defined (time, pressure) table.
    pub table: Vec<(f64, f64)>,
}
/// State of a single element in erosion tracking.
#[derive(Debug, Clone)]
pub struct ElementErosionState {
    /// Whether the element has been eroded.
    pub eroded: bool,
    /// Current damage variable (0 = intact, 1 = fully damaged).
    pub damage: f64,
    /// Current effective plastic strain.
    pub plastic_strain: f64,
    /// Maximum principal stress.
    pub max_principal_stress: f64,
    /// Current volume relative to initial.
    pub volume_ratio: f64,
    /// Mass associated with this element.
    pub mass: f64,
    /// Erosion time step (when element was eroded, `u64::MAX` if not).
    pub eroded_at_step: u64,
}
impl ElementErosionState {
    /// Create a new element state.
    pub fn new(mass: f64) -> Self {
        Self {
            eroded: false,
            damage: 0.0,
            plastic_strain: 0.0,
            max_principal_stress: 0.0,
            volume_ratio: 1.0,
            mass,
            eroded_at_step: u64::MAX,
        }
    }
    /// Check if element should be eroded based on the given criterion.
    pub fn should_erode(&self, criterion: &FailureCriterion) -> bool {
        match criterion {
            FailureCriterion::PlasticStrain(eps_f) => self.plastic_strain >= *eps_f,
            FailureCriterion::PrincipalStress(sigma_f) => self.max_principal_stress >= *sigma_f,
            FailureCriterion::Damage(d_f) => self.damage >= *d_f,
            FailureCriterion::VolumeFraction(v_f) => self.volume_ratio <= *v_f,
        }
    }
}
/// A contact pair: slave node index and master segment (two node indices).
#[derive(Debug, Clone)]
pub struct ContactPair {
    /// Index of the slave node.
    pub slave_node: usize,
    /// Indices of the two master segment nodes.
    pub master_segment: [usize; 2],
    /// Current gap value (negative = penetration).
    pub gap: f64,
    /// Contact normal (pointing slave → master).
    pub normal: [f64; 3],
    /// Contact pressure.
    pub pressure: f64,
    /// Lagrange multiplier (for Lagrange enforcement).
    pub lambda: f64,
}
impl ContactPair {
    /// Create a new contact pair with zero gap.
    pub fn new(slave_node: usize, master_segment: [usize; 2]) -> Self {
        Self {
            slave_node,
            master_segment,
            gap: 0.0,
            normal: [0.0, 1.0, 0.0],
            pressure: 0.0,
            lambda: 0.0,
        }
    }
}
/// Configuration for the explicit central difference solver.
#[derive(Debug, Clone)]
pub struct ExplicitFemConfig {
    /// CFL safety factor (0 < cfl ≤ 1). Typical value: 0.9.
    pub cfl_factor: f64,
    /// Mass lumping strategy.
    pub lumping: MassLumping,
    /// Material wave speed (m/s), used for critical time step estimate.
    pub wave_speed: f64,
    /// Minimum characteristic element length (m).
    pub min_element_length: f64,
    /// Bulk viscosity coefficient (dimensionless), for shock smoothing.
    pub bulk_viscosity: f64,
    /// Whether to use energy balance checking.
    pub check_energy: bool,
}
/// Element erosion algorithm for high-deformation impact simulations.
///
/// When an element's damage/strain exceeds the failure criterion, it is
/// deleted from the mesh. Mass is redistributed to neighboring nodes to
/// conserve total mass.
pub struct ErosionAlgorithm {
    /// Failure criterion for element deletion.
    pub criterion: FailureCriterion,
    /// Per-element erosion states.
    pub states: Vec<ElementErosionState>,
    /// Total eroded mass.
    pub total_eroded_mass: f64,
    /// Number of elements eroded.
    pub n_eroded: usize,
    /// Whether to redistribute mass to neighboring nodes.
    pub conserve_mass: bool,
    /// Debris particle positions (for debris tracking).
    pub debris_positions: Vec<[f64; 3]>,
    /// Debris particle velocities.
    pub debris_velocities: Vec<[f64; 3]>,
    /// Debris particle masses.
    pub debris_masses: Vec<f64>,
}
impl ErosionAlgorithm {
    /// Create a new erosion algorithm.
    pub fn new(criterion: FailureCriterion, n_elements: usize, conserve_mass: bool) -> Self {
        Self {
            criterion,
            states: vec![ElementErosionState::new(0.0); n_elements],
            total_eroded_mass: 0.0,
            n_eroded: 0,
            conserve_mass,
            debris_positions: Vec::new(),
            debris_velocities: Vec::new(),
            debris_masses: Vec::new(),
        }
    }
    /// Initialize element masses.
    pub fn set_element_mass(&mut self, element_idx: usize, mass: f64) {
        self.states[element_idx].mass = mass;
    }
    /// Update damage state and check for erosion. Returns `true` if eroded.
    ///
    /// # Arguments
    /// * `element_idx` — element index
    /// * `damage` — new damage value
    /// * `plastic_strain` — new effective plastic strain
    /// * `max_principal_stress` — new maximum principal stress
    /// * `volume_ratio` — current volume / initial volume
    /// * `step` — current time step index
    pub fn update_element(
        &mut self,
        element_idx: usize,
        damage: f64,
        plastic_strain: f64,
        max_principal_stress: f64,
        volume_ratio: f64,
        step: u64,
    ) -> bool {
        let state = &mut self.states[element_idx];
        if state.eroded {
            return false;
        }
        state.damage = damage;
        state.plastic_strain = plastic_strain;
        state.max_principal_stress = max_principal_stress;
        state.volume_ratio = volume_ratio;
        if state.should_erode(&self.criterion) {
            state.eroded = true;
            state.eroded_at_step = step;
            self.total_eroded_mass += state.mass;
            self.n_eroded += 1;
            true
        } else {
            false
        }
    }
    /// Redistribute eroded element mass to neighboring nodes.
    ///
    /// `node_masses` — global lumped mass vector
    /// `element_nodes` — node indices for the eroded element (4 for tet)
    pub fn redistribute_mass_to_nodes(
        &self,
        element_idx: usize,
        node_masses: &mut [f64],
        element_nodes: &[usize],
    ) {
        if !self.conserve_mass {
            return;
        }
        let state = &self.states[element_idx];
        let share = state.mass / element_nodes.len() as f64;
        for &node in element_nodes {
            node_masses[node] += share;
        }
    }
    /// Create a debris particle at the centroid of the eroded element.
    ///
    /// `node_positions` — nodal position array (flat, 3 per node)
    /// `node_velocities` — nodal velocity array (flat, 3 per node)
    /// `element_nodes` — node indices for the element
    pub fn spawn_debris(
        &mut self,
        element_idx: usize,
        node_positions: &[f64],
        node_velocities: &[f64],
        element_nodes: &[usize],
    ) {
        if !self.states[element_idx].eroded {
            return;
        }
        let n = element_nodes.len() as f64;
        let pos = element_nodes.iter().fold([0.0f64; 3], |mut acc, &ni| {
            acc[0] += node_positions[ni * 3] / n;
            acc[1] += node_positions[ni * 3 + 1] / n;
            acc[2] += node_positions[ni * 3 + 2] / n;
            acc
        });
        let vel = element_nodes.iter().fold([0.0f64; 3], |mut acc, &ni| {
            acc[0] += node_velocities[ni * 3] / n;
            acc[1] += node_velocities[ni * 3 + 1] / n;
            acc[2] += node_velocities[ni * 3 + 2] / n;
            acc
        });
        self.debris_positions.push(pos);
        self.debris_velocities.push(vel);
        self.debris_masses.push(self.states[element_idx].mass);
    }
    /// Return fraction of elements that have been eroded.
    pub fn erosion_fraction(&self) -> f64 {
        if self.states.is_empty() {
            return 0.0;
        }
        self.n_eroded as f64 / self.states.len() as f64
    }
    /// Return indices of all eroded elements.
    pub fn eroded_elements(&self) -> Vec<usize> {
        self.states
            .iter()
            .enumerate()
            .filter(|(_, s)| s.eroded)
            .map(|(i, _)| i)
            .collect()
    }
    /// Return indices of all active (non-eroded) elements.
    pub fn active_elements(&self) -> Vec<usize> {
        self.states
            .iter()
            .enumerate()
            .filter(|(_, s)| !s.eroded)
            .map(|(i, _)| i)
            .collect()
    }
}
/// Result of a wave propagation analysis.
#[derive(Debug, Clone)]
pub struct WaveResult {
    /// Frequencies (Hz).
    pub frequencies: Vec<f64>,
    /// Phase velocities (m/s).
    pub phase_velocities: Vec<f64>,
    /// Group velocities (m/s).
    pub group_velocities: Vec<f64>,
    /// Wave numbers (rad/m).
    pub wave_numbers: Vec<f64>,
    /// Dispersion parameter (dimensionless).
    pub dispersion_parameter: Vec<f64>,
}
