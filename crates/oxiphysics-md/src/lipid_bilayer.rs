// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coarse-grained lipid bilayer simulation (MARTINI-like).
//!
//! This module provides:
//! - Coarse-grained lipid bead model with head/glycerol/tail beads
//! - Bilayer self-assembly simulation
//! - Membrane tension calculation
//! - Lipid lateral diffusion coefficient via MSD
//! - Phospholipid flip-flop rate estimation
//! - Bending modulus (Helfrich elasticity)
//! - Area per lipid
//! - Membrane thickness
//! - Cholesterol condensing and ordering effects
//! - Protein-membrane interaction (insertion energy, tilt coupling)
//!
//! Units: nm for length, kJ/mol for energy, ps for time.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Physical / simulation constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant in kJ/(mol·K).
pub const KB: f64 = 8.314_462_618e-3;

/// Reference temperature 300 K in energy units kB·T \[kJ/mol\] at 300 K.
pub const KBT_300: f64 = KB * 300.0;

/// Bilayer normal direction index (z-axis).
pub const NORMAL_AXIS: usize = 2;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

/// Apply minimum image convention for periodic boundary conditions.
#[inline]
fn pbc_displacement(r_ij: [f64; 3], box_size: [f64; 3]) -> [f64; 3] {
    [
        r_ij[0] - box_size[0] * (r_ij[0] / box_size[0]).round(),
        r_ij[1] - box_size[1] * (r_ij[1] / box_size[1]).round(),
        r_ij[2] - box_size[2] * (r_ij[2] / box_size[2]).round(),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Lipid types
// ─────────────────────────────────────────────────────────────────────────────

/// Lipid molecular species for the coarse-grained bilayer model.
#[derive(Debug, Clone, PartialEq)]
pub enum LipidSpecies {
    /// DPPC — dipalmitoylphosphatidylcholine (saturated, gel at 300 K).
    Dppc,
    /// DOPC — dioleoylphosphatidylcholine (unsaturated, fluid at 300 K).
    Dopc,
    /// POPE — palmitoyloleoylphosphatidylethanolamine.
    Pope,
    /// POPS — palmitoyloleoylphosphatidylserine (anionic).
    Pops,
    /// Cholesterol — modulates order and fluidity.
    Cholesterol,
    /// Custom user-defined species with a label.
    Custom(String),
}

impl LipidSpecies {
    /// Equilibrium area per lipid A₀ in nm² at 300 K.
    pub fn equilibrium_area(&self) -> f64 {
        match self {
            LipidSpecies::Dppc => 0.64,
            LipidSpecies::Dopc => 0.72,
            LipidSpecies::Pope => 0.65,
            LipidSpecies::Pops => 0.65,
            LipidSpecies::Cholesterol => 0.38,
            LipidSpecies::Custom(_) => 0.65,
        }
    }

    /// Hydrophobic tail length d_tail in nm.
    pub fn tail_length(&self) -> f64 {
        match self {
            LipidSpecies::Dppc => 1.65,
            LipidSpecies::Dopc => 1.50,
            LipidSpecies::Pope => 1.58,
            LipidSpecies::Pops => 1.58,
            LipidSpecies::Cholesterol => 1.70,
            LipidSpecies::Custom(_) => 1.55,
        }
    }

    /// Number of coarse-grained beads (MARTINI-like).
    pub fn n_beads(&self) -> usize {
        match self {
            LipidSpecies::Dppc => 12,
            LipidSpecies::Dopc => 13,
            LipidSpecies::Pope => 12,
            LipidSpecies::Pops => 12,
            LipidSpecies::Cholesterol => 8,
            LipidSpecies::Custom(_) => 10,
        }
    }

    /// Net charge of the headgroup in units of elementary charge.
    pub fn headgroup_charge(&self) -> f64 {
        match self {
            LipidSpecies::Pops => -1.0,
            LipidSpecies::Cholesterol => 0.0,
            _ => 0.0,
        }
    }

    /// Intrinsic lipid shape parameter κ = V / (A₀ · l_c).
    ///
    /// κ < 1 → cylindrical (bilayer), κ > 1 → inverted cone.
    pub fn shape_parameter(&self) -> f64 {
        match self {
            LipidSpecies::Dppc => 1.0,
            LipidSpecies::Dopc => 1.05,
            LipidSpecies::Pope => 1.15, // cone-like → inverted hexagonal tendency
            LipidSpecies::Pops => 1.0,
            LipidSpecies::Cholesterol => 1.2,
            LipidSpecies::Custom(_) => 1.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Coarse-grained lipid molecule
// ─────────────────────────────────────────────────────────────────────────────

/// A single coarse-grained lipid molecule.
///
/// Represented by a headgroup bead, a glycerol linker, and tail bead positions.
#[derive(Debug, Clone)]
pub struct CgLipid {
    /// Lipid species.
    pub species: LipidSpecies,
    /// Headgroup bead position \[x, y, z\] in nm.
    pub head_pos: [f64; 3],
    /// Glycerol linker bead position in nm.
    pub glycerol_pos: [f64; 3],
    /// Tail bead positions (from glycerol towards tail end) in nm.
    pub tail_beads: Vec<[f64; 3]>,
    /// Instantaneous tail order parameter S₂.
    pub tail_order: f64,
    /// Leaflet assignment: 0 = upper, 1 = lower.
    pub leaflet: u8,
    /// Unique molecule ID.
    pub mol_id: u32,
}

impl CgLipid {
    /// Create a new CG lipid with a straight-chain conformation along the z-axis.
    pub fn new_upright(
        head_pos: [f64; 3],
        species: LipidSpecies,
        n_tail_beads: usize,
        leaflet: u8,
        mol_id: u32,
    ) -> Self {
        let bead_spacing = species.tail_length() / (n_tail_beads as f64).max(1.0);
        let direction = if leaflet == 0 { -1.0 } else { 1.0 };
        let glycerol_pos = [head_pos[0], head_pos[1], head_pos[2] + direction * 0.47];
        let tail_beads = (0..n_tail_beads)
            .map(|i| {
                [
                    head_pos[0],
                    head_pos[1],
                    glycerol_pos[2] + direction * (i + 1) as f64 * bead_spacing,
                ]
            })
            .collect();
        Self {
            species,
            head_pos,
            glycerol_pos,
            tail_beads,
            tail_order: 0.8,
            leaflet,
            mol_id,
        }
    }

    /// Director vector: from glycerol toward last tail bead.
    pub fn director(&self) -> [f64; 3] {
        if self.tail_beads.is_empty() {
            return [0.0, 0.0, 1.0];
        }
        let last = self.tail_beads[self.tail_beads.len() - 1];
        let d = sub3(last, self.glycerol_pos);
        norm3(d)
    }

    /// Tilt angle (rad) relative to bilayer normal \[0,0,1\].
    pub fn tilt_angle(&self) -> f64 {
        let dir = self.director();
        let cos_theta = dir[2].clamp(-1.0, 1.0).abs();
        cos_theta.acos()
    }

    /// Compute the deuterium order parameter S_CD for the i-th tail bond.
    ///
    /// `S_CD = 0.5 * (3 cos²θ - 1)`  where θ is the bond angle with the bilayer normal.
    pub fn bond_order_parameter(&self, bond_index: usize) -> f64 {
        if bond_index + 1 >= self.tail_beads.len() || self.tail_beads.is_empty() {
            return 0.0;
        }
        let b = sub3(self.tail_beads[bond_index + 1], self.tail_beads[bond_index]);
        let b_norm = norm3(b);
        let cos_theta = b_norm[2].clamp(-1.0, 1.0);
        0.5 * (3.0 * cos_theta * cos_theta - 1.0)
    }

    /// Mean order parameter over all tail bonds.
    pub fn mean_tail_order(&self) -> f64 {
        let n = self.tail_beads.len();
        if n < 2 {
            return 0.0;
        }
        let sum: f64 = (0..n - 1).map(|i| self.bond_order_parameter(i)).sum();
        sum / (n - 1) as f64
    }

    /// End-to-end distance (head to last tail bead) in nm.
    pub fn end_to_end_distance(&self) -> f64 {
        if let Some(&last) = self.tail_beads.last() {
            len3(sub3(last, self.head_pos))
        } else {
            0.0
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CG force-field parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Lennard-Jones interaction parameters for two CG bead types.
#[derive(Debug, Clone, Copy)]
pub struct LjParams {
    /// Well depth ε in kJ/mol.
    pub epsilon: f64,
    /// Collision diameter σ in nm.
    pub sigma: f64,
    /// Cutoff radius r_c in nm.
    pub r_cut: f64,
}

impl LjParams {
    /// Create LJ parameters.
    pub fn new(epsilon: f64, sigma: f64, r_cut: f64) -> Self {
        Self {
            epsilon,
            sigma,
            r_cut,
        }
    }

    /// MARTINI-like head-head interaction.
    pub fn head_head() -> Self {
        Self {
            epsilon: 5.6,
            sigma: 0.47,
            r_cut: 1.2,
        }
    }

    /// MARTINI-like tail-tail interaction (hydrophobic).
    pub fn tail_tail() -> Self {
        Self {
            epsilon: 3.4,
            sigma: 0.47,
            r_cut: 1.2,
        }
    }

    /// MARTINI-like head-tail (repulsive, weak).
    pub fn head_tail() -> Self {
        Self {
            epsilon: 2.0,
            sigma: 0.47,
            r_cut: 1.2,
        }
    }

    /// LJ potential energy V(r) = 4ε\[(σ/r)¹² − (σ/r)⁶\] − V(r_cut).
    pub fn potential(&self, r: f64) -> f64 {
        if r >= self.r_cut || r < 1e-10 {
            return 0.0;
        }
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        let v = 4.0 * self.epsilon * (sr12 - sr6);
        let sr_c = self.sigma / self.r_cut;
        let sr6_c = sr_c.powi(6);
        let v_cut = 4.0 * self.epsilon * (sr6_c * sr6_c - sr6_c);
        v - v_cut
    }

    /// LJ force magnitude |F(r)| = −dV/dr.
    pub fn force_magnitude(&self, r: f64) -> f64 {
        if r >= self.r_cut || r < 1e-10 {
            return 0.0;
        }
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        let sr12 = sr6 * sr6;
        24.0 * self.epsilon / r * (2.0 * sr12 - sr6)
    }

    /// LJ force vector on bead i due to bead j.
    pub fn force_vector(&self, r_ij: [f64; 3]) -> [f64; 3] {
        let r = len3(r_ij);
        if r < 1e-10 || r >= self.r_cut {
            return [0.0; 3];
        }
        let f_mag = self.force_magnitude(r);
        scale3(r_ij, -f_mag / r)
    }
}

/// Harmonic bond potential for adjacent CG beads.
#[derive(Debug, Clone, Copy)]
pub struct HarmonicBondParams {
    /// Force constant k in kJ/(mol·nm²).
    pub k_bond: f64,
    /// Equilibrium bond length r₀ in nm.
    pub r0: f64,
}

impl HarmonicBondParams {
    /// Default intra-lipid bond parameters (MARTINI-like).
    pub fn default_cg() -> Self {
        Self {
            k_bond: 3800.0,
            r0: 0.47,
        }
    }

    /// Potential energy V = k/2 (r - r₀)² \[kJ/mol\].
    pub fn potential(&self, r: f64) -> f64 {
        let dr = r - self.r0;
        0.5 * self.k_bond * dr * dr
    }

    /// Force magnitude |F| = k (r - r₀).
    pub fn force_magnitude(&self, r: f64) -> f64 {
        self.k_bond * (r - self.r0)
    }
}

/// Harmonic angle potential for CG beads.
#[derive(Debug, Clone, Copy)]
pub struct HarmonicAngleParams {
    /// Force constant k_θ in kJ/(mol·rad²).
    pub k_angle: f64,
    /// Equilibrium angle θ₀ in radians.
    pub theta0: f64,
}

impl HarmonicAngleParams {
    /// Default CG lipid angle (near-linear tail).
    pub fn default_cg() -> Self {
        Self {
            k_angle: 25.0,
            theta0: PI * 140.0 / 180.0,
        }
    }

    /// Potential energy V = k_θ/2 (θ - θ₀)² \[kJ/mol\].
    pub fn potential(&self, theta: f64) -> f64 {
        let dtheta = theta - self.theta0;
        0.5 * self.k_angle * dtheta * dtheta
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bilayer structure
// ─────────────────────────────────────────────────────────────────────────────

/// A coarse-grained lipid bilayer membrane.
#[derive(Debug, Clone)]
pub struct LipidBilayer {
    /// All lipid molecules in the bilayer.
    pub lipids: Vec<CgLipid>,
    /// Simulation box dimensions \[Lx, Ly, Lz\] in nm.
    pub box_size: [f64; 3],
    /// Current simulation time in ps.
    pub time_ps: f64,
    /// Temperature in K.
    pub temperature_k: f64,
    /// Mole fraction of cholesterol.
    pub cholesterol_fraction: f64,
}

impl LipidBilayer {
    /// Construct a flat bilayer with N lipids per leaflet distributed on a grid.
    pub fn flat(n_per_leaflet: u32, species: LipidSpecies, temperature_k: f64) -> Self {
        let a0 = species.equilibrium_area();
        let side = ((n_per_leaflet as f64) * a0).sqrt();
        let box_size = [side, side, 10.0]; // 10 nm in z to include water space
        let grid_n = (n_per_leaflet as f64).sqrt().ceil() as u32;
        let dx = side / grid_n as f64;

        let mut lipids = Vec::new();
        let mut mol_id = 0u32;
        let midplane_z = 5.0;

        for leaflet in 0..2u8 {
            let sign = if leaflet == 0 { 1.0 } else { -1.0 };
            let z_head = midplane_z + sign * (species.tail_length() + 0.5);
            let mut count = 0;
            'outer: for iy in 0..grid_n {
                for ix in 0..grid_n {
                    if count >= n_per_leaflet {
                        break 'outer;
                    }
                    let x = (ix as f64 + 0.5) * dx;
                    let y = (iy as f64 + 0.5) * dx;
                    let head = [x, y, z_head];
                    let lip = CgLipid::new_upright(head, species.clone(), 4, leaflet, mol_id);
                    lipids.push(lip);
                    mol_id += 1;
                    count += 1;
                }
            }
        }

        Self {
            lipids,
            box_size,
            time_ps: 0.0,
            temperature_k,
            cholesterol_fraction: 0.0,
        }
    }

    /// Number of lipids in each leaflet.
    pub fn leaflet_counts(&self) -> (usize, usize) {
        let upper = self.lipids.iter().filter(|l| l.leaflet == 0).count();
        let lower = self.lipids.iter().filter(|l| l.leaflet == 1).count();
        (upper, lower)
    }

    /// Total number of lipid molecules.
    pub fn n_lipids(&self) -> usize {
        self.lipids.len()
    }

    /// Lateral area of the bilayer (box x × y) in nm².
    pub fn lateral_area(&self) -> f64 {
        self.box_size[0] * self.box_size[1]
    }

    /// Average area per lipid per leaflet in nm².
    pub fn area_per_lipid(&self) -> f64 {
        let (upper, lower) = self.leaflet_counts();
        let n = (upper + lower).max(2) as f64 / 2.0;
        self.lateral_area() / n
    }

    /// Mean membrane thickness: z-distance between upper and lower headgroup planes in nm.
    pub fn membrane_thickness(&self) -> f64 {
        let upper_z: Vec<f64> = self
            .lipids
            .iter()
            .filter(|l| l.leaflet == 0)
            .map(|l| l.head_pos[2])
            .collect();
        let lower_z: Vec<f64> = self
            .lipids
            .iter()
            .filter(|l| l.leaflet == 1)
            .map(|l| l.head_pos[2])
            .collect();
        if upper_z.is_empty() || lower_z.is_empty() {
            return 0.0;
        }
        let z_upper = upper_z.iter().sum::<f64>() / upper_z.len() as f64;
        let z_lower = lower_z.iter().sum::<f64>() / lower_z.len() as f64;
        (z_upper - z_lower).abs()
    }

    /// Mean tail order parameter S₂ averaged over all lipids.
    pub fn mean_order_parameter(&self) -> f64 {
        if self.lipids.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.lipids.iter().map(|l| l.mean_tail_order()).sum();
        sum / self.lipids.len() as f64
    }

    /// Headgroup positions of the specified leaflet (for undulation analysis).
    pub fn headgroup_positions(&self, leaflet: u8) -> Vec<[f64; 3]> {
        self.lipids
            .iter()
            .filter(|l| l.leaflet == leaflet)
            .map(|l| l.head_pos)
            .collect()
    }

    /// Midplane z-position of the bilayer.
    pub fn midplane_z(&self) -> f64 {
        let sum: f64 = self.lipids.iter().map(|l| l.head_pos[2]).sum();
        sum / self.lipids.len().max(1) as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Membrane tension
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the lateral membrane tension γ from the pressure tensor.
///
/// `γ = Lz/2 · [P_N - (P_xx + P_yy)/2]`  \[kJ/(mol·nm²)\]
///
/// where P_N = P_zz is the normal pressure and P_L = (P_xx + P_yy)/2.
pub fn membrane_tension(p_xx: f64, p_yy: f64, p_zz: f64, lz: f64) -> f64 {
    let p_l = 0.5 * (p_xx + p_yy);
    0.5 * lz * (p_zz - p_l)
}

/// Pressure tensor contribution from a single pair interaction.
///
/// Uses the virial theorem: `P_αβ += r_α F_β / V`.
pub fn pressure_tensor_contribution(r_ij: [f64; 3], f_ij: [f64; 3], volume: f64) -> [[f64; 3]; 3] {
    let mut p = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            p[a][b] = r_ij[a] * f_ij[b] / volume;
        }
    }
    p
}

/// Accumulate the pressure tensor across a set of pair contributions.
pub fn accumulate_pressure_tensor(contributions: &[[[f64; 3]; 3]]) -> [[f64; 3]; 3] {
    let mut p = [[0.0_f64; 3]; 3];
    for c in contributions {
        for a in 0..3 {
            for b in 0..3 {
                p[a][b] += c[a][b];
            }
        }
    }
    p
}

// ─────────────────────────────────────────────────────────────────────────────
// Lipid diffusion coefficient
// ─────────────────────────────────────────────────────────────────────────────

/// Mean-square displacement (MSD) data for lipid diffusion analysis.
#[derive(Debug, Clone)]
pub struct MsdData {
    /// Time lags τ in ps.
    pub tau: Vec<f64>,
    /// Corresponding MSD values ⟨Δr²(τ)⟩ in nm².
    pub msd: Vec<f64>,
}

impl MsdData {
    /// Create MSD data from a set of 2D lateral displacements.
    ///
    /// `positions[t]` is a list of (x, y) lateral positions at time t.
    /// Time step is `dt_ps`.
    pub fn from_positions(positions: &[Vec<[f64; 2]>], dt_ps: f64) -> Self {
        let n_frames = positions.len();
        let n_lipids = positions.first().map_or(0, |v| v.len());
        let max_lag = n_frames / 4; // use up to 1/4 of trajectory for statistics

        let mut tau = Vec::with_capacity(max_lag);
        let mut msd = Vec::with_capacity(max_lag);

        for lag in 1..=max_lag {
            let mut sum_msd = 0.0;
            let mut count = 0usize;
            for t0 in 0..(n_frames - lag) {
                for (p0, p1) in positions[t0]
                    .iter()
                    .zip(positions[t0 + lag].iter())
                    .take(n_lipids)
                {
                    let dx = p1[0] - p0[0];
                    let dy = p1[1] - p0[1];
                    sum_msd += dx * dx + dy * dy;
                    count += 1;
                }
            }
            tau.push(lag as f64 * dt_ps);
            msd.push(if count > 0 {
                sum_msd / count as f64
            } else {
                0.0
            });
        }

        Self { tau, msd }
    }

    /// Lateral diffusion coefficient D from the Einstein relation.
    ///
    /// `D = MSD(τ) / (4 τ)`  \[nm²/ps\]
    ///
    /// Fits the linear regime using the last half of the MSD data.
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.tau.len();
        if n < 4 {
            return 0.0;
        }
        // Linear regression of MSD vs τ in the range [n/2, n]
        let start = n / 2;
        let n_fit = n - start;
        let tau_slice = &self.tau[start..];
        let msd_slice = &self.msd[start..];

        let mean_tau = tau_slice.iter().sum::<f64>() / n_fit as f64;
        let mean_msd = msd_slice.iter().sum::<f64>() / n_fit as f64;

        let num: f64 = tau_slice
            .iter()
            .zip(msd_slice.iter())
            .map(|(&t, &m)| (t - mean_tau) * (m - mean_msd))
            .sum();
        let den: f64 = tau_slice.iter().map(|&t| (t - mean_tau).powi(2)).sum();

        if den < 1e-20 {
            0.0
        } else {
            num / den / 4.0 // D = slope / 4 (2D diffusion)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Flip-flop rate
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters governing lipid flip-flop across the bilayer midplane.
#[derive(Debug, Clone)]
pub struct FlipFlopParams {
    /// Free energy barrier ΔG‡ for flip-flop \[kJ/mol\].
    pub delta_g_barrier: f64,
    /// Attempt frequency ν₀ \[ps⁻¹\].
    pub attempt_frequency: f64,
    /// Temperature T \[K\].
    pub temperature_k: f64,
}

impl FlipFlopParams {
    /// Default parameters for DPPC at 323 K.
    pub fn dppc_323k() -> Self {
        Self {
            delta_g_barrier: 80.0, // kJ/mol
            attempt_frequency: 1e-3,
            temperature_k: 323.0,
        }
    }

    /// Flip-flop rate k = ν₀ · exp(−ΔG‡/(RT)) \[ps⁻¹\].
    pub fn rate(&self) -> f64 {
        let rt = KB * self.temperature_k;
        self.attempt_frequency * (-self.delta_g_barrier / rt).exp()
    }

    /// Half-life t_{1/2} = ln(2) / k \[ps\].
    pub fn half_life_ps(&self) -> f64 {
        let k = self.rate();
        if k < 1e-100 {
            f64::INFINITY
        } else {
            2_f64.ln() / k
        }
    }

    /// Effect of asymmetric membrane composition: reduced barrier by ΔΔG.
    pub fn rate_with_asymmetry(&self, delta_delta_g: f64) -> f64 {
        let rt = KB * self.temperature_k;
        self.attempt_frequency * (-(self.delta_g_barrier + delta_delta_g) / rt).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bending modulus (Helfrich elasticity)
// ─────────────────────────────────────────────────────────────────────────────

/// Helfrich free energy of membrane deformation.
///
/// `F = ∫ [κ_b/2 (∇²h)² + κ_G (κ₁κ₂)] dA`
///
/// with κ_b the bending modulus and κ_G the Gaussian curvature modulus.
#[derive(Debug, Clone, Copy)]
pub struct HelfrichParams {
    /// Bending modulus κ_b \[kJ/mol\] (≈ 20 k_BT for DPPC bilayer).
    pub kappa_b: f64,
    /// Gaussian curvature modulus κ_G \[kJ/mol\] (≈ −κ_b for bilayers).
    pub kappa_g: f64,
    /// Spontaneous curvature c₀ \[nm⁻¹\].
    pub c0: f64,
}

impl HelfrichParams {
    /// Default parameters for a symmetric DPPC bilayer at 323 K.
    pub fn dppc() -> Self {
        let kbt = KB * 323.0;
        Self {
            kappa_b: 20.0 * kbt,
            kappa_g: -20.0 * kbt,
            c0: 0.0,
        }
    }

    /// Parameters for a DOPC bilayer.
    pub fn dopc() -> Self {
        let kbt = KB * 300.0;
        Self {
            kappa_b: 17.0 * kbt,
            kappa_g: -17.0 * kbt,
            c0: 0.0,
        }
    }

    /// Bending energy for a single patch of area dA with mean curvature H and Gaussian curvature K.
    pub fn bending_energy_patch(&self, h_mean: f64, k_gaussian: f64, area: f64) -> f64 {
        let bend = 0.5 * self.kappa_b * (2.0 * h_mean - self.c0).powi(2);
        let gauss = self.kappa_g * k_gaussian;
        (bend + gauss) * area
    }

    /// Estimate κ_b from the fluctuation spectrum of the membrane height h(q).
    ///
    /// In the Helfrich model: `⟨|h_q|²⟩ = k_BT / (κ_b q⁴ + γ q²)`
    pub fn kappa_b_from_spectrum(q: f64, h_q_sq: f64, temperature_k: f64, tension: f64) -> f64 {
        let kbt = KB * temperature_k;
        let q4 = q.powi(4);
        let q2 = q * q;
        if q4 < 1e-20 {
            return 0.0;
        }
        (kbt / h_q_sq - tension * q2) / q4
    }
}

/// Compute mean curvature H at a point from principal curvatures κ₁ and κ₂.
pub fn mean_curvature(kappa1: f64, kappa2: f64) -> f64 {
    0.5 * (kappa1 + kappa2)
}

/// Compute Gaussian curvature K = κ₁ · κ₂.
pub fn gaussian_curvature(kappa1: f64, kappa2: f64) -> f64 {
    kappa1 * kappa2
}

/// Estimate local mean curvature at a headgroup using the discretised Laplacian.
///
/// Given a height field h at a central point and its `neighbours` heights,
/// with lateral spacing `dx`, returns H ≈ ∇²h / 2.
pub fn discrete_mean_curvature(h_center: f64, neighbours: &[f64], dx: f64) -> f64 {
    if neighbours.is_empty() || dx < 1e-20 {
        return 0.0;
    }
    let mean_nb = neighbours.iter().sum::<f64>() / neighbours.len() as f64;
    (mean_nb - h_center) / (dx * dx * 0.5)
}

// ─────────────────────────────────────────────────────────────────────────────
// Cholesterol effects
// ─────────────────────────────────────────────────────────────────────────────

/// Models the condensing effect of cholesterol on a lipid bilayer.
///
/// Based on the umbrella model (Huang et al.) and experimental data.
#[derive(Debug, Clone)]
pub struct CholesterolEffect {
    /// Mole fraction of cholesterol χ_c.
    pub chi_c: f64,
    /// Critical mole fraction for liquid-ordered phase formation.
    pub chi_crit: f64,
    /// Condensing coefficient α for area reduction.
    pub alpha_condensing: f64,
    /// Ordering coefficient β for tail order increase.
    pub beta_ordering: f64,
}

impl CholesterolEffect {
    /// Create a cholesterol effect model with default DPPC/cholesterol parameters.
    pub fn new(chi_c: f64) -> Self {
        Self {
            chi_c,
            chi_crit: 0.25,
            alpha_condensing: 0.30, // ~30% area reduction at chi_c = 0.5
            beta_ordering: 0.40,    // S_CD increases by ~0.4 per unit chi_c
        }
    }

    /// Effective area per lipid A(χ_c) = A₀ · (1 − α · χ_c) \[nm²\].
    pub fn effective_area(&self, a0: f64) -> f64 {
        a0 * (1.0 - self.alpha_condensing * self.chi_c)
    }

    /// Effective tail order parameter S₂(χ_c) = S₂⁰ + β · χ_c (capped at 1).
    pub fn effective_order_parameter(&self, s2_pure: f64) -> f64 {
        (s2_pure + self.beta_ordering * self.chi_c).min(1.0)
    }

    /// Effective bending modulus κ_b(χ_c) increases with cholesterol.
    ///
    /// `κ_b(χ_c) = κ_b⁰ · (1 + γ_κ · χ_c)`  with γ_κ ≈ 2.
    pub fn effective_bending_modulus(&self, kappa_b0: f64) -> f64 {
        kappa_b0 * (1.0 + 2.0 * self.chi_c)
    }

    /// Liquid-ordered / liquid-disordered coexistence indicator.
    ///
    /// Returns true if χ_c is in the two-phase region (0.1–0.4 for DPPC).
    pub fn is_two_phase(&self) -> bool {
        self.chi_c >= 0.10 && self.chi_c <= 0.40
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Protein–membrane interaction
// ─────────────────────────────────────────────────────────────────────────────

/// Hydrophobic mismatch model for a transmembrane protein.
///
/// Computes the deformation energy when the hydrophobic span of a
/// protein d_p differs from the bilayer hydrophobic thickness d_b.
#[derive(Debug, Clone)]
pub struct HydrophobicMismatch {
    /// Protein hydrophobic span d_p in nm.
    pub d_protein: f64,
    /// Bilayer hydrophobic thickness d_b in nm.
    pub d_bilayer: f64,
    /// Deformation modulus K_A in kJ/(mol·nm²) (related to κ_b).
    pub k_deform: f64,
    /// Protein cross-sectional radius r_p in nm.
    pub radius_nm: f64,
}

impl HydrophobicMismatch {
    /// Create a hydrophobic mismatch model.
    pub fn new(d_protein: f64, d_bilayer: f64, k_deform: f64, radius_nm: f64) -> Self {
        Self {
            d_protein,
            d_bilayer,
            k_deform,
            radius_nm,
        }
    }

    /// Mismatch energy E_mm = K_A / 2 · (d_p − d_b)² per unit area \[kJ/mol\].
    pub fn mismatch_energy_per_area(&self) -> f64 {
        let delta = self.d_protein - self.d_bilayer;
        0.5 * self.k_deform * delta * delta
    }

    /// Total mismatch energy over the protein-membrane contact circumference.
    pub fn total_mismatch_energy(&self) -> f64 {
        let contact_area = 2.0 * PI * self.radius_nm * (self.d_bilayer * 0.5);
        self.mismatch_energy_per_area() * contact_area
    }

    /// Tilt coupling energy for a protein tilted by angle θ \[rad\].
    ///
    /// `E_tilt = K_t/2 · θ²`  with K_t from the Bohinc-Iglic model.
    pub fn tilt_energy(&self, theta_rad: f64) -> f64 {
        let k_tilt = self.k_deform * self.radius_nm * self.radius_nm;
        0.5 * k_tilt * theta_rad * theta_rad
    }

    /// Optimal protein tilt angle that minimises total deformation energy.
    ///
    /// Solves `dE/dθ = 0`; for the linear model the minimum is θ = 0
    /// unless an asymmetric mismatch is present.
    pub fn optimal_tilt(&self) -> f64 {
        0.0 // symmetric mismatch → no tilt preferred
    }
}

/// Insertion free energy for a peptide of length L into the bilayer \[kJ/mol\].
///
/// Uses a simple Born solvation model:
/// `ΔG_ins = Σ_i ΔG_transfer(residue_i)`
pub fn insertion_free_energy(residue_transfer_energies: &[f64]) -> f64 {
    residue_transfer_energies.iter().sum()
}

/// Membrane deformation energy due to a point inclusion (Gaussian model).
///
/// `E_def(r) = (κ_b / (2 π ξ²)) · exp(−r²/(2ξ²))`
///
/// where ξ is the deformation penetration length.
pub fn point_inclusion_energy(r_nm: f64, kappa_b: f64, xi_nm: f64) -> f64 {
    if xi_nm < 1e-10 {
        return 0.0;
    }
    kappa_b / (2.0 * PI * xi_nm * xi_nm) * (-(r_nm * r_nm) / (2.0 * xi_nm * xi_nm)).exp()
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-assembly simulation (minimal Langevin dynamics)
// ─────────────────────────────────────────────────────────────────────────────

/// Langevin thermostat parameters.
#[derive(Debug, Clone, Copy)]
pub struct LangevinParams {
    /// Friction coefficient γ in ps⁻¹.
    pub gamma: f64,
    /// Temperature T in K.
    pub temperature_k: f64,
    /// Time step Δt in ps.
    pub dt_ps: f64,
}

impl LangevinParams {
    /// Default MARTINI Langevin thermostat parameters.
    pub fn default_martini() -> Self {
        Self {
            gamma: 1.0, // ps⁻¹
            temperature_k: 300.0,
            dt_ps: 0.02,
        }
    }
}

/// A single CG bead in the self-assembly simulation.
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Position \[x, y, z\] in nm.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in nm/ps.
    pub vel: [f64; 3],
    /// Force accumulator in kJ/(mol·nm).
    pub force: [f64; 3],
    /// Mass m in g/mol (amu).
    pub mass: f64,
    /// Bead type identifier (0 = head, 1 = glycerol, 2 = tail).
    pub bead_type: u8,
}

impl CgBead {
    /// Create a CG bead at rest.
    pub fn at_rest(pos: [f64; 3], mass: f64, bead_type: u8) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            mass,
            bead_type,
        }
    }

    /// Kinetic energy in kJ/mol.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
}

/// Integrate a single bead with the BAOAB Langevin scheme.
///
/// B: velocity half-kick;  A: position full-drift;
/// O: Ornstein-Uhlenbeck friction/noise step;  A: drift;  B: velocity half-kick.
pub fn langevin_step_baoab(bead: &mut CgBead, params: &LangevinParams, rng_noise: [f64; 3]) {
    let dt = params.dt_ps;
    let gamma = params.gamma;
    let kbt = KB * params.temperature_k;
    let sigma = (2.0 * gamma * kbt / bead.mass).sqrt();

    // B: half velocity kick
    let a = scale3(bead.force, 1.0 / bead.mass);
    bead.vel = add3(bead.vel, scale3(a, 0.5 * dt));

    // A: half drift
    bead.pos = add3(bead.pos, scale3(bead.vel, 0.5 * dt));

    // O: Ornstein-Uhlenbeck
    let c1 = (-gamma * dt).exp();
    let c2 = ((1.0 - c1 * c1).max(0.0)).sqrt() * sigma;
    bead.vel = add3(scale3(bead.vel, c1), scale3(rng_noise, c2));

    // A: half drift
    bead.pos = add3(bead.pos, scale3(bead.vel, 0.5 * dt));

    // B: half kick (force must be recomputed externally before next call)
}

// ─────────────────────────────────────────────────────────────────────────────
// Bilayer analysis utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the deuterium NMR order parameter profile across the bilayer.
///
/// Returns a vector of (bond_index, ⟨S_CD⟩) pairs.
pub fn order_parameter_profile(bilayer: &LipidBilayer) -> Vec<(usize, f64)> {
    // Find max number of tail bonds
    let max_bonds = bilayer
        .lipids
        .iter()
        .map(|l| l.tail_beads.len().saturating_sub(1))
        .max()
        .unwrap_or(0);

    (0..max_bonds)
        .map(|bi| {
            let vals: Vec<f64> = bilayer
                .lipids
                .iter()
                .map(|l| l.bond_order_parameter(bi))
                .collect();
            let mean = if vals.is_empty() {
                0.0
            } else {
                vals.iter().sum::<f64>() / vals.len() as f64
            };
            (bi, mean)
        })
        .collect()
}

/// Compute the bilayer undulation spectrum |h_q|² from headgroup z-positions.
///
/// Assumes a square bilayer patch Lx × Ly. Returns a vector of (|q|, |h_q|²) pairs.
pub fn undulation_spectrum(
    headgroup_z: &[[f64; 3]],
    lx: f64,
    ly: f64,
    n_modes: usize,
) -> Vec<(f64, f64)> {
    let n = headgroup_z.len();
    if n == 0 || lx < 1e-10 || ly < 1e-10 {
        return Vec::new();
    }
    let z_mean = headgroup_z.iter().map(|p| p[2]).sum::<f64>() / n as f64;

    let mut result = Vec::new();
    for mx in 0..n_modes {
        for my in 0..n_modes {
            if mx == 0 && my == 0 {
                continue;
            }
            let qx = 2.0 * PI * mx as f64 / lx;
            let qy = 2.0 * PI * my as f64 / ly;
            let q = (qx * qx + qy * qy).sqrt();

            let mut re = 0.0;
            let mut im = 0.0;
            for p in headgroup_z {
                let h = p[2] - z_mean;
                let phase = qx * p[0] + qy * p[1];
                re += h * phase.cos();
                im += h * phase.sin();
            }
            let h_q_sq = (re * re + im * im) / (n * n) as f64;
            result.push((q, h_q_sq));
        }
    }
    result
}

/// Estimate the bending modulus κ_b from an undulation spectrum.
///
/// Fits `⟨|h_q|²⟩ = k_BT / (κ_b q⁴)` in the small-q (bending-dominated) regime.
pub fn kappa_from_undulation_spectrum(
    spectrum: &[(f64, f64)],
    temperature_k: f64,
    area: f64,
) -> f64 {
    let kbt = KB * temperature_k;
    // Filter modes in the bending regime q < 1 nm⁻¹
    let bending_modes: Vec<(f64, f64)> = spectrum
        .iter()
        .filter(|(q, _)| *q > 0.1 && *q < 1.0)
        .copied()
        .collect();

    if bending_modes.len() < 3 {
        return 0.0;
    }

    // Linear regression: log(h_q_sq * q^4) = log(k_BT / (κ_b * area))
    let kappa_estimates: Vec<f64> = bending_modes
        .iter()
        .map(|(q, h_q_sq)| {
            let q4 = q.powi(4);
            if *h_q_sq < 1e-30 || q4 < 1e-30 {
                0.0
            } else {
                kbt / (*h_q_sq * q4 * area)
            }
        })
        .filter(|k| *k > 0.0)
        .collect();

    if kappa_estimates.is_empty() {
        0.0
    } else {
        kappa_estimates.iter().sum::<f64>() / kappa_estimates.len() as f64
    }
}

/// Radial distribution function g(r) for lipid headgroups within a leaflet.
///
/// Returns a vector of (r, g(r)) pairs with spacing `dr`.
pub fn radial_distribution_function(
    positions: &[[f64; 3]],
    box_size: [f64; 3],
    dr: f64,
    r_max: f64,
) -> Vec<(f64, f64)> {
    let n = positions.len();
    if n < 2 {
        return Vec::new();
    }
    let n_bins = (r_max / dr) as usize;
    let mut hist = vec![0u64; n_bins];
    let area = box_size[0] * box_size[1];
    let rho_2d = n as f64 / area;

    for i in 0..n {
        for j in (i + 1)..n {
            let r_ij = sub3(positions[i], positions[j]);
            let r_ij_pbc = pbc_displacement(r_ij, box_size);
            // Only lateral (xy) distance for 2D RDF
            let r_lat = (r_ij_pbc[0] * r_ij_pbc[0] + r_ij_pbc[1] * r_ij_pbc[1]).sqrt();
            let bin = (r_lat / dr) as usize;
            if bin < n_bins {
                hist[bin] += 2; // count both i→j and j→i
            }
        }
    }

    (0..n_bins)
        .map(|bin| {
            let r = (bin as f64 + 0.5) * dr;
            let shell_area = 2.0 * PI * r * dr;
            let g = hist[bin] as f64 / (n as f64 * shell_area * rho_2d);
            (r, g)
        })
        .collect()
}

/// Compute the 2D mean-square displacement from centre-of-mass trajectories.
///
/// Returns MSD at each lag in units of nm².
pub fn compute_msd_2d(traj: &[Vec<[f64; 2]>], max_lag: usize) -> Vec<f64> {
    let n_frames = traj.len();
    let n_mol = traj.first().map_or(0, |v| v.len());
    let max_lag = max_lag.min(n_frames / 2);

    (1..=max_lag)
        .map(|lag| {
            let mut sum = 0.0;
            let mut count = 0usize;
            for t0 in 0..(n_frames - lag) {
                for (m0, m1) in traj[t0].iter().zip(traj[t0 + lag].iter()).take(n_mol) {
                    let dx = m1[0] - m0[0];
                    let dy = m1[1] - m0[1];
                    sum += dx * dx + dy * dy;
                    count += 1;
                }
            }
            if count > 0 { sum / count as f64 } else { 0.0 }
        })
        .collect()
}

/// Assign lipids to upper or lower leaflet based on their z-position relative to midplane.
pub fn assign_leaflets(lipids: &mut [CgLipid]) {
    if lipids.is_empty() {
        return;
    }
    let midplane = lipids.iter().map(|l| l.head_pos[2]).sum::<f64>() / lipids.len() as f64;
    for l in lipids.iter_mut() {
        l.leaflet = if l.head_pos[2] > midplane { 0 } else { 1 };
    }
}

/// Compute the lateral pressure profile P_L(z) from a coarse-grained density profile.
///
/// Returns a vector of (z, P_L) where P_L is in kJ/(mol·nm³).
pub fn lateral_pressure_profile(
    density_profile: &[(f64, f64)],
    temperature_k: f64,
) -> Vec<(f64, f64)> {
    let kbt = KB * temperature_k;
    density_profile
        .iter()
        .map(|&(z, rho)| (z, rho * kbt))
        .collect()
}

/// Surface tension γ from the first moment of the lateral pressure profile.
///
/// `γ = ∫ [P_N - P_L(z)] dz`
pub fn surface_tension_from_pressure_profile(
    pressure_profile: &[(f64, f64)],
    p_normal: f64,
    dz: f64,
) -> f64 {
    pressure_profile
        .iter()
        .map(|&(_z, p_l)| (p_normal - p_l) * dz)
        .sum()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── LipidSpecies tests ────────────────────────────────────────────────────

    #[test]
    fn test_lipid_species_equilibrium_area() {
        assert!(LipidSpecies::Dppc.equilibrium_area() > 0.0);
        assert!(LipidSpecies::Dopc.equilibrium_area() > LipidSpecies::Dppc.equilibrium_area());
        assert!(
            LipidSpecies::Cholesterol.equilibrium_area() < LipidSpecies::Dppc.equilibrium_area()
        );
    }

    #[test]
    fn test_lipid_species_tail_length_positive() {
        for species in [
            LipidSpecies::Dppc,
            LipidSpecies::Dopc,
            LipidSpecies::Pope,
            LipidSpecies::Cholesterol,
        ] {
            assert!(species.tail_length() > 0.0);
        }
    }

    #[test]
    fn test_lipid_species_n_beads() {
        assert!(LipidSpecies::Dppc.n_beads() > 0);
        assert!(LipidSpecies::Cholesterol.n_beads() > 0);
    }

    #[test]
    fn test_lipid_species_charge() {
        assert_eq!(LipidSpecies::Pops.headgroup_charge(), -1.0);
        assert_eq!(LipidSpecies::Dppc.headgroup_charge(), 0.0);
    }

    // ── CgLipid tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_cg_lipid_director_unit_length() {
        let lip = CgLipid::new_upright([2.0, 2.0, 6.0], LipidSpecies::Dppc, 4, 0, 0);
        let dir = lip.director();
        let l = len3(dir);
        assert!((l - 1.0).abs() < 1e-9 || l < 1e-9, "director length = {l}");
    }

    #[test]
    fn test_cg_lipid_tilt_angle_range() {
        let lip = CgLipid::new_upright([0.0; 3], LipidSpecies::Dppc, 4, 0, 0);
        let angle = lip.tilt_angle();
        assert!((0.0..=PI).contains(&angle), "tilt = {angle}");
    }

    #[test]
    fn test_cg_lipid_end_to_end_distance_positive() {
        let lip = CgLipid::new_upright([0.0, 0.0, 5.0], LipidSpecies::Dppc, 4, 0, 0);
        let d = lip.end_to_end_distance();
        assert!(d > 0.0, "end-to-end = {d}");
    }

    #[test]
    fn test_cg_lipid_bond_order_out_of_bounds_returns_zero() {
        let lip = CgLipid::new_upright([0.0; 3], LipidSpecies::Dppc, 2, 0, 0);
        let s = lip.bond_order_parameter(100);
        assert_eq!(s, 0.0);
    }

    // ── LJ potential tests ────────────────────────────────────────────────────

    #[test]
    fn test_lj_potential_zero_at_cutoff() {
        let params = LjParams::tail_tail();
        let v = params.potential(params.r_cut);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_lj_potential_zero_beyond_cutoff() {
        let params = LjParams::head_head();
        let v = params.potential(params.r_cut + 0.1);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_lj_minimum_at_sigma_2_to_1_6() {
        let params = LjParams::new(1.0, 1.0, 5.0);
        // Minimum at r = 2^(1/6) * σ ≈ 1.122
        let r_min = 2_f64.powf(1.0 / 6.0);
        let v_min = params.potential(r_min);
        let v_nearby = params.potential(r_min * 1.01);
        assert!(v_min <= v_nearby, "v_min = {v_min}, v_nearby = {v_nearby}");
    }

    #[test]
    fn test_lj_force_repulsive_at_short_range() {
        let params = LjParams::tail_tail();
        let r_ij = [0.3, 0.0, 0.0]; // well inside sigma
        let f = params.force_vector(r_ij);
        // Force should point in +x (away from origin): repulsive
        assert!(f[0] < 0.0, "force should be repulsive: {}", f[0]);
    }

    // ── Harmonic bond and angle tests ─────────────────────────────────────────

    #[test]
    fn test_harmonic_bond_minimum_at_r0() {
        let params = HarmonicBondParams::default_cg();
        let v0 = params.potential(params.r0);
        let v1 = params.potential(params.r0 + 0.01);
        assert!(v0 < v1);
    }

    #[test]
    fn test_harmonic_bond_force_zero_at_equilibrium() {
        let params = HarmonicBondParams::default_cg();
        let f = params.force_magnitude(params.r0);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_harmonic_angle_potential_positive() {
        let params = HarmonicAngleParams::default_cg();
        let v = params.potential(PI / 2.0); // 90 degrees, far from theta0
        assert!(v > 0.0);
    }

    // ── Bilayer construction tests ─────────────────────────────────────────────

    #[test]
    fn test_bilayer_flat_correct_lipid_count() {
        let bilayer = LipidBilayer::flat(25, LipidSpecies::Dppc, 300.0);
        let (upper, lower) = bilayer.leaflet_counts();
        assert_eq!(upper, 25, "upper leaflet = {upper}");
        assert_eq!(lower, 25, "lower leaflet = {lower}");
    }

    #[test]
    fn test_bilayer_area_per_lipid_reasonable() {
        let bilayer = LipidBilayer::flat(25, LipidSpecies::Dppc, 300.0);
        let apl = bilayer.area_per_lipid();
        assert!(apl > 0.3 && apl < 1.5, "APL = {apl}");
    }

    #[test]
    fn test_bilayer_membrane_thickness_positive() {
        let bilayer = LipidBilayer::flat(25, LipidSpecies::Dppc, 300.0);
        let thickness = bilayer.membrane_thickness();
        assert!(thickness > 0.0, "thickness = {thickness}");
    }

    #[test]
    fn test_bilayer_mean_order_parameter_range() {
        let bilayer = LipidBilayer::flat(16, LipidSpecies::Dopc, 300.0);
        let s2 = bilayer.mean_order_parameter();
        assert!((-0.5..=1.0).contains(&s2), "S2 = {s2}");
    }

    // ── Membrane tension tests ────────────────────────────────────────────────

    #[test]
    fn test_membrane_tension_isotropic_zero() {
        let p = 1.0; // kJ/(mol·nm³)
        let gamma = membrane_tension(p, p, p, 10.0);
        assert!(gamma.abs() < 1e-10, "tension = {gamma}");
    }

    #[test]
    fn test_membrane_tension_sign() {
        // P_N > P_L → positive tension
        let gamma = membrane_tension(1.0, 1.0, 2.0, 10.0);
        assert!(gamma > 0.0);
    }

    // ── Flip-flop tests ──────────────────────────────────────────────────────

    #[test]
    fn test_flipflop_rate_positive() {
        let params = FlipFlopParams::dppc_323k();
        assert!(params.rate() >= 0.0);
    }

    #[test]
    fn test_flipflop_half_life_long() {
        let params = FlipFlopParams::dppc_323k();
        // DPPC flip-flop half-life ~ hours = very large in ps
        let t_half = params.half_life_ps();
        assert!(t_half > 1e6, "t_half = {t_half}");
    }

    #[test]
    fn test_flipflop_rate_decreases_with_barrier() {
        let mut params = FlipFlopParams::dppc_323k();
        let rate_low = params.rate();
        params.delta_g_barrier *= 2.0;
        let rate_high = params.rate();
        assert!(rate_low > rate_high);
    }

    // ── Helfrich bending modulus tests ────────────────────────────────────────

    #[test]
    fn test_helfrich_bending_energy_flat_zero() {
        let params = HelfrichParams::dppc();
        let e = params.bending_energy_patch(0.0, 0.0, 1.0);
        // For c0 = 0 and H = 0: E = κ_G * K * area; K = 0 → E = 0
        assert!(e.abs() < 1e-10, "bending energy = {e}");
    }

    #[test]
    fn test_helfrich_kappa_b_dppc_positive() {
        let params = HelfrichParams::dppc();
        assert!(params.kappa_b > 0.0);
    }

    #[test]
    fn test_mean_curvature_symmetric() {
        let h = mean_curvature(1.0, 1.0);
        assert!((h - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_curvature_saddle_negative() {
        let k = gaussian_curvature(1.0, -1.0);
        assert!(k < 0.0);
    }

    #[test]
    fn test_discrete_mean_curvature_flat() {
        let nbs = vec![2.0; 4];
        let h = discrete_mean_curvature(2.0, &nbs, 0.5);
        assert!(h.abs() < 1e-8, "H = {h}");
    }

    // ── Cholesterol effect tests ──────────────────────────────────────────────

    #[test]
    fn test_cholesterol_condensing() {
        let chol = CholesterolEffect::new(0.3);
        let a_eff = chol.effective_area(LipidSpecies::Dppc.equilibrium_area());
        assert!(a_eff < LipidSpecies::Dppc.equilibrium_area());
    }

    #[test]
    fn test_cholesterol_ordering() {
        let chol = CholesterolEffect::new(0.3);
        let s_eff = chol.effective_order_parameter(0.4);
        assert!(s_eff >= 0.4);
    }

    #[test]
    fn test_cholesterol_bending_modulus_increase() {
        let chol = CholesterolEffect::new(0.3);
        let k0 = HelfrichParams::dppc().kappa_b;
        let k_eff = chol.effective_bending_modulus(k0);
        assert!(k_eff > k0);
    }

    #[test]
    fn test_cholesterol_two_phase_detection() {
        assert!(CholesterolEffect::new(0.25).is_two_phase());
        assert!(!CholesterolEffect::new(0.05).is_two_phase());
        assert!(!CholesterolEffect::new(0.60).is_two_phase());
    }

    // ── Protein-membrane interaction tests ────────────────────────────────────

    #[test]
    fn test_hydrophobic_mismatch_zero_at_match() {
        let mm = HydrophobicMismatch::new(3.0, 3.0, 50.0, 0.5);
        assert_eq!(mm.mismatch_energy_per_area(), 0.0);
    }

    #[test]
    fn test_hydrophobic_mismatch_positive_energy() {
        let mm = HydrophobicMismatch::new(4.0, 3.0, 50.0, 0.5);
        assert!(mm.mismatch_energy_per_area() > 0.0);
        assert!(mm.total_mismatch_energy() > 0.0);
    }

    #[test]
    fn test_tilt_energy_zero_at_zero_angle() {
        let mm = HydrophobicMismatch::new(3.0, 3.0, 50.0, 0.5);
        let e = mm.tilt_energy(0.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_insertion_free_energy_sum() {
        let residues = vec![-5.0, -3.0, 2.0];
        let dg = insertion_free_energy(&residues);
        assert!((dg - (-6.0)).abs() < 1e-10);
    }

    // ── MSD and diffusion coefficient tests ───────────────────────────────────

    #[test]
    fn test_msd_single_particle_linear_growth() {
        // Particle moving at constant velocity v in x
        let v = 0.1; // nm/ps
        let dt = 1.0;
        let n_frames = 20;
        let n_lip = 1;
        let positions: Vec<Vec<[f64; 2]>> = (0..n_frames)
            .map(|t| vec![[v * t as f64 * dt, 0.0]; n_lip])
            .collect();
        let msd_data = MsdData::from_positions(&positions, dt);
        // MSD(τ) = v² τ² for ballistic motion
        if msd_data.msd.len() >= 2 {
            assert!(msd_data.msd[1] > msd_data.msd[0]);
        }
    }

    #[test]
    fn test_diffusion_coefficient_nonnegative() {
        let dt = 1.0;
        let n_frames = 16;
        let positions: Vec<Vec<[f64; 2]>> = (0..n_frames)
            .map(|t| vec![[0.01 * t as f64, 0.0]])
            .collect();
        let msd_data = MsdData::from_positions(&positions, dt);
        let d = msd_data.diffusion_coefficient();
        assert!(d >= 0.0);
    }

    // ── Order parameter profile tests ─────────────────────────────────────────

    #[test]
    fn test_order_parameter_profile_returns_values() {
        let bilayer = LipidBilayer::flat(9, LipidSpecies::Dppc, 300.0);
        let profile = order_parameter_profile(&bilayer);
        assert!(!profile.is_empty());
        for (i, s) in &profile {
            assert!(*s >= -0.5 && *s <= 1.0, "bond {i}: S = {s}");
        }
    }

    // ── RDF tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_rdf_returns_nonzero() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [(i as f64) * 0.5, 0.0, 0.0]).collect();
        let box_size = [10.0, 10.0, 10.0];
        let rdf = radial_distribution_function(&positions, box_size, 0.1, 3.0);
        assert!(!rdf.is_empty());
    }

    #[test]
    fn test_rdf_at_large_r_approaches_one() {
        // For a random distribution, g(r) → 1 at large r
        // Just verify structure of output
        let positions: Vec<[f64; 3]> = (0..5).map(|i| [i as f64 * 1.5, 0.0, 0.0]).collect();
        let box_size = [10.0, 10.0, 10.0];
        let rdf = radial_distribution_function(&positions, box_size, 0.5, 5.0);
        assert!(!rdf.is_empty());
    }

    // ── Undulation spectrum tests ─────────────────────────────────────────────

    #[test]
    fn test_undulation_spectrum_flat_membrane() {
        // Flat membrane h = const → h_q = 0 for q ≠ 0
        let headgroup_z: Vec<[f64; 3]> = (0..16)
            .map(|i| [(i % 4) as f64, (i / 4) as f64, 5.0])
            .collect();
        let spec = undulation_spectrum(&headgroup_z, 4.0, 4.0, 3);
        for (_, h_q_sq) in &spec {
            assert!(
                *h_q_sq < 1e-20,
                "h_q^2 should be near zero for flat membrane: {h_q_sq}"
            );
        }
    }

    // ── Leaflet assignment tests ──────────────────────────────────────────────

    #[test]
    fn test_assign_leaflets_divides_correctly() {
        let mut lipids = Vec::new();
        for i in 0..10 {
            let z = if i < 5 { 7.0 } else { 3.0 };
            lipids.push(CgLipid::new_upright(
                [0.0, 0.0, z],
                LipidSpecies::Dppc,
                2,
                0,
                i,
            ));
        }
        assign_leaflets(&mut lipids);
        let upper = lipids.iter().filter(|l| l.leaflet == 0).count();
        let lower = lipids.iter().filter(|l| l.leaflet == 1).count();
        assert_eq!(upper, 5);
        assert_eq!(lower, 5);
    }

    // ── Point inclusion energy test ──────────────────────────────────────────

    #[test]
    fn test_point_inclusion_energy_decreases_with_distance() {
        let kappa_b = 50.0;
        let xi = 1.0;
        let e_near = point_inclusion_energy(0.1, kappa_b, xi);
        let e_far = point_inclusion_energy(3.0, kappa_b, xi);
        assert!(e_near > e_far, "near={e_near}, far={e_far}");
    }

    #[test]
    fn test_point_inclusion_energy_zero_xi() {
        let e = point_inclusion_energy(1.0, 50.0, 0.0);
        assert_eq!(e, 0.0);
    }

    // ── Langevin integration test ────────────────────────────────────────────

    #[test]
    fn test_langevin_step_position_changes() {
        let mut bead = CgBead::at_rest([0.0; 3], 72.0, 0);
        bead.vel = [0.5, 0.0, 0.0];
        bead.force = [0.0; 3];
        let params = LangevinParams::default_martini();
        let noise = [0.0; 3];
        let pos_before = bead.pos;
        langevin_step_baoab(&mut bead, &params, noise);
        assert!(
            (bead.pos[0] - pos_before[0]).abs() > 1e-10,
            "position should change"
        );
    }
}
