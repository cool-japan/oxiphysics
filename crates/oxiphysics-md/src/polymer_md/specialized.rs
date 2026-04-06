// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Specialized polymer models: DendriticPolymer, PolymerAdsorption, GoModel,
//! ElasticNetworkModel, SelfAvoidingWalk, StarPolymer, CombPolymer.

use super::{dist3_vec, dot3_vec, sub3_vec};
use std::f64::consts::PI;

/// Dendritic polymer (dendrimer) generation model.
///
/// Computes branching, fractal dimension, and end-group statistics.
#[derive(Debug, Clone)]
pub struct DendriticPolymer {
    /// Number of generations g.
    pub generation: usize,
    /// Core branching multiplicity nc (e.g., 3).
    pub core_multiplicity: usize,
    /// Branch multiplicity nb (e.g., 2).
    pub branch_multiplicity: usize,
    /// Spacer length ns (monomers between branching points).
    pub spacer_length: usize,
    /// Segment length b.
    pub b: f64,
}

impl DendriticPolymer {
    /// Create a new dendrimer.
    pub fn new(
        generation: usize,
        core_multiplicity: usize,
        branch_multiplicity: usize,
        spacer_length: usize,
        b: f64,
    ) -> Self {
        Self {
            generation,
            core_multiplicity,
            branch_multiplicity,
            spacer_length,
            b,
        }
    }

    /// Number of end groups.
    pub fn n_end_groups(&self) -> usize {
        self.core_multiplicity * self.branch_multiplicity.pow(self.generation as u32)
    }

    /// Total number of branching points.
    pub fn n_branching_points(&self) -> usize {
        // Sum over generations
        (0..self.generation)
            .map(|g| self.core_multiplicity * self.branch_multiplicity.pow(g as u32))
            .sum::<usize>()
            + 1 // include core
    }

    /// Total number of monomers.
    pub fn n_total_monomers(&self) -> usize {
        let n_segments = self.n_end_groups() + self.n_branching_points() - 1;
        n_segments * self.spacer_length
    }

    /// Fractal dimension Df (Flory exponent based).
    pub fn fractal_dimension(&self) -> f64 {
        // For dendrimers: Df ~ 3 at high generation (dense packing)
        let g = self.generation as f64;
        let n = self.n_end_groups() as f64;
        if n < 1.0 {
            return 3.0;
        }
        3.0 * g.ln() / (n.ln() + 1.0)
    }

    /// Radius estimate: R ~ b * sqrt(g * ns).
    pub fn radius(&self) -> f64 {
        let g = self.generation as f64;
        let ns = self.spacer_length as f64;
        self.b * (g * ns).sqrt()
    }

    /// Density: rho ~ N / R^3.
    pub fn density(&self) -> f64 {
        let r = self.radius();
        if r < 1e-30 {
            return 0.0;
        }
        self.n_total_monomers() as f64 / (r * r * r)
    }
}

/// Polymer adsorption on a surface.
///
/// Tracks train/loop/tail fractions and adsorption energy.
#[derive(Debug, Clone)]
pub struct PolymerAdsorption {
    /// Number of beads.
    pub n_beads: usize,
    /// Adsorption energy per contact epsilon (in kT units).
    pub epsilon: f64,
    /// Fraction of beads in trains (adsorbed).
    pub f_train: f64,
    /// Fraction of beads in loops.
    pub f_loop: f64,
    /// Fraction of beads in tails.
    pub f_tail: f64,
    /// Surface coverage theta.
    pub theta: f64,
    /// Bulk concentration c.
    pub c_bulk: f64,
    /// Temperature kT.
    pub kt: f64,
}

impl PolymerAdsorption {
    /// Create a new polymer adsorption model.
    pub fn new(n_beads: usize, epsilon: f64, kt: f64, c_bulk: f64) -> Self {
        // Initial guess: de Gennes scaling
        let f_train = 0.1;
        let f_loop = 0.6;
        let f_tail = 0.3;
        let theta = c_bulk * n_beads as f64 * epsilon.exp();
        Self {
            n_beads,
            epsilon,
            f_train,
            f_loop,
            f_tail,
            theta,
            c_bulk,
            kt,
        }
    }

    /// Adsorption free energy per chain.
    pub fn adsorption_energy(&self) -> f64 {
        -self.epsilon * self.f_train * self.n_beads as f64 * self.kt
    }

    /// Number of contacts with surface.
    pub fn n_contacts(&self) -> f64 {
        self.f_train * self.n_beads as f64
    }

    /// Effective adsorption thickness delta.
    pub fn layer_thickness(&self, b: f64) -> f64 {
        // de Gennes: delta ~ b * (N_loop)^(3/5)
        let n_loop = self.f_loop * self.n_beads as f64;
        b * n_loop.powf(0.6)
    }

    /// Desorption probability at temperature T.
    pub fn desorption_probability(&self) -> f64 {
        (-self.adsorption_energy().abs() / self.kt).exp()
    }

    /// Update fractions based on temperature (simple mean-field).
    pub fn update_fractions(&mut self) {
        let exp_e = (-self.epsilon).exp();
        let total = 1.0 + exp_e;
        self.f_train = 1.0 / total;
        self.f_loop = exp_e / total * 0.7;
        self.f_tail = exp_e / total * 0.3;
    }

    /// Langmuir isotherm surface coverage.
    pub fn langmuir_coverage(&self) -> f64 {
        let k = self.epsilon.exp();
        k * self.c_bulk / (1.0 + k * self.c_bulk)
    }
}

/// Cα Go-model for protein folding simulations.
///
/// Each residue is represented as a single bead at the Cα position.
/// Native contacts are assigned based on crystal structure.
#[derive(Debug, Clone)]
pub struct GoModel {
    /// Number of residues.
    pub n_residues: usize,
    /// Cα bead positions in native state.
    pub native_positions: Vec<[f64; 3]>,
    /// Native contact list: (i, j, r_native).
    pub native_contacts: Vec<(usize, usize, f64)>,
    /// Go model energy parameter ε.
    pub epsilon: f64,
    /// Bond length.
    pub bond_length: f64,
    /// Current positions.
    pub positions: Vec<[f64; 3]>,
}

impl GoModel {
    /// Create a new Go model from native positions.
    pub fn new(native_positions: Vec<[f64; 3]>, epsilon: f64, cutoff: f64) -> Self {
        let n = native_positions.len();
        let bond_length = if n > 1 {
            dist3_vec(native_positions[0], native_positions[1])
        } else {
            3.8
        };
        // Find native contacts (non-bonded pairs within cutoff)
        let mut native_contacts = Vec::new();
        for i in 0..n {
            for j in (i + 2)..n {
                let r = dist3_vec(native_positions[i], native_positions[j]);
                if r < cutoff {
                    native_contacts.push((i, j, r));
                }
            }
        }
        let positions = native_positions.clone();
        GoModel {
            n_residues: n,
            native_positions,
            native_contacts,
            epsilon,
            bond_length,
            positions,
        }
    }

    /// Go model potential energy.
    ///
    /// Attractive 10-12 LJ for native contacts, repulsive for non-native.
    pub fn potential_energy(&self) -> f64 {
        let mut energy = 0.0f64;
        // Native contacts: Lennard-Jones 10-12
        for &(i, j, r0) in &self.native_contacts {
            let r = dist3_vec(self.positions[i], self.positions[j]).max(1e-10);
            let ratio = r0 / r;
            energy += self.epsilon * (5.0 * ratio.powi(12) - 6.0 * ratio.powi(10));
        }
        energy
    }

    /// Fraction of native contacts Q.
    pub fn fraction_native_contacts(&self) -> f64 {
        if self.native_contacts.is_empty() {
            return 1.0;
        }
        let threshold_factor = 1.2;
        let n_formed = self
            .native_contacts
            .iter()
            .filter(|&&(i, j, r0)| {
                dist3_vec(self.positions[i], self.positions[j]) < threshold_factor * r0
            })
            .count();
        n_formed as f64 / self.native_contacts.len() as f64
    }

    /// Root-mean-square deviation (RMSD) from native structure.
    pub fn rmsd_from_native(&self) -> f64 {
        if self.n_residues == 0 {
            return 0.0;
        }
        let msd = self
            .positions
            .iter()
            .zip(self.native_positions.iter())
            .map(|(&p, &n)| {
                let d = sub3_vec(p, n);
                dot3_vec(d, d)
            })
            .sum::<f64>()
            / self.n_residues as f64;
        msd.sqrt()
    }

    /// Is the protein folded (Q > 0.7)?
    pub fn is_folded(&self) -> bool {
        self.fraction_native_contacts() > 0.7
    }
}

/// Elastic Network Model (ENM) for normal mode analysis of proteins.
///
/// Each residue pair within a cutoff is connected by a spring of constant γ.
#[derive(Debug, Clone)]
pub struct ElasticNetworkModel {
    /// Number of nodes (residues).
    pub n: usize,
    /// Node positions.
    pub positions: Vec<[f64; 3]>,
    /// Spring constant γ.
    pub gamma: f64,
    /// Cutoff distance R_c.
    pub r_cutoff: f64,
    /// Contact pairs (i, j).
    pub contacts: Vec<(usize, usize)>,
}

impl ElasticNetworkModel {
    /// Create a new ENM from positions.
    pub fn new(positions: Vec<[f64; 3]>, gamma: f64, r_cutoff: f64) -> Self {
        let n = positions.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if dist3_vec(positions[i], positions[j]) < r_cutoff {
                    contacts.push((i, j));
                }
            }
        }
        ElasticNetworkModel {
            n,
            positions,
            gamma,
            r_cutoff,
            contacts,
        }
    }

    /// Number of contacts in the network.
    pub fn n_contacts(&self) -> usize {
        self.contacts.len()
    }

    /// Kirchhoff (Laplacian) matrix element (i, j) for 1D projection.
    ///
    /// K_ij = -γ if (i,j) in contacts, K_ii = sum of off-diagonal.
    pub fn kirchhoff_element(&self, i: usize, j: usize) -> f64 {
        if i == j {
            self.contacts
                .iter()
                .filter(|&&(a, b)| a == i || b == i)
                .map(|_| self.gamma)
                .sum()
        } else if self.contacts.contains(&(i.min(j), i.max(j))) {
            -self.gamma
        } else {
            0.0
        }
    }

    /// Mean-square fluctuation of residue i: <Δrᵢ²> ~ (K⁻¹)ᵢᵢ.
    ///
    /// Approximated using the connectivity (inverse of Kirchhoff diagonal).
    pub fn mean_sq_fluctuation(&self, i: usize) -> f64 {
        let k_ii = self.kirchhoff_element(i, i);
        if k_ii < 1e-30 {
            f64::INFINITY
        } else {
            1.0 / k_ii
        }
    }

    /// B-factor of residue i (proportional to mean-square fluctuation).
    pub fn b_factor(&self, i: usize) -> f64 {
        8.0 * PI * PI / 3.0 * self.mean_sq_fluctuation(i)
    }

    /// Average B-factor over all residues.
    pub fn average_b_factor(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        (0..self.n).map(|i| self.b_factor(i)).sum::<f64>() / self.n as f64
    }
}

/// Self-avoiding walk (SAW) model for excluded-volume polymers.
///
/// Uses Flory theory: R ~ N^ν b where ν = 3/5 (d=3).
#[derive(Debug, Clone)]
pub struct SelfAvoidingWalk {
    /// Number of steps N.
    pub n: usize,
    /// Step length b.
    pub b: f64,
    /// Flory exponent ν.
    pub nu: f64,
    /// Excluded volume parameter v.
    pub v: f64,
}

impl SelfAvoidingWalk {
    /// Create a new SAW model.
    pub fn new(n: usize, b: f64) -> Self {
        SelfAvoidingWalk {
            n,
            b,
            nu: 0.588,
            v: b * b * b,
        }
    }

    /// Flory radius: R_F = b * N^ν.
    pub fn flory_radius(&self) -> f64 {
        self.b * (self.n as f64).powf(self.nu)
    }

    /// End-to-end distance `R²` = b² N^(2ν).
    pub fn mean_r2(&self) -> f64 {
        let rf = self.flory_radius();
        rf * rf
    }

    /// Radius of gyration: Rg ~ Rf / sqrt(6) for SAW.
    pub fn radius_of_gyration(&self) -> f64 {
        self.flory_radius() / 6.0_f64.sqrt()
    }

    /// Second virial coefficient A₂ ~ N^(3-dν) in d=3.
    pub fn second_virial_coefficient(&self) -> f64 {
        let n = self.n as f64;
        self.v * n * n / self.flory_radius().powi(3)
    }

    /// Osmotic pressure from virial expansion π/ρ = 1 + A₂ρ + ...
    pub fn osmotic_pressure(&self, concentration: f64) -> f64 {
        let rho = concentration;
        rho * (1.0 + self.second_virial_coefficient() * rho)
    }
}

/// Star polymer with f arms.
///
/// Daoud-Cotton blob model for star polymers in good solvent.
#[derive(Debug, Clone)]
pub struct StarPolymer {
    /// Arm length N (beads per arm).
    pub n_arm: usize,
    /// Number of arms f.
    pub n_arms: usize,
    /// Segment length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl StarPolymer {
    /// Create a new star polymer.
    pub fn new(n_arm: usize, n_arms: usize, b: f64, kt: f64) -> Self {
        StarPolymer {
            n_arm,
            n_arms,
            b,
            kt,
        }
    }

    /// Total number of monomers.
    pub fn n_total(&self) -> usize {
        self.n_arm * self.n_arms
    }

    /// Star radius R_star ~ b * N^(3/5) * f^(1/5).
    pub fn star_radius(&self) -> f64 {
        let n = self.n_arm as f64;
        let f = self.n_arms as f64;
        self.b * n.powf(0.6) * f.powf(0.2)
    }

    /// Osmotic pressure of star polymer solution.
    ///
    /// π ~ f^(3/2) kT / R³.
    pub fn osmotic_pressure(&self) -> f64 {
        let f = self.n_arms as f64;
        let r = self.star_radius();
        f.powf(1.5) * self.kt / (r * r * r)
    }

    /// Radius of gyration for star (Zimm-Stockmayer).
    ///
    /// Rg² = N_arm b² / 3 (3f - 2) / f².
    pub fn radius_of_gyration_sq(&self) -> f64 {
        let n = self.n_arm as f64;
        let f = self.n_arms as f64;
        n * self.b * self.b / 3.0 * (3.0 * f - 2.0) / (f * f)
    }

    /// Second virial coefficient A₂ for star polymers.
    pub fn second_virial(&self) -> f64 {
        let r = self.star_radius();
        r * r * r * (self.n_arms as f64).powf(0.5)
    }
}

/// Comb polymer (backbone + side chains).
#[derive(Debug, Clone)]
pub struct CombPolymer {
    /// Number of backbone monomers.
    pub n_backbone: usize,
    /// Number of side chains.
    pub n_side_chains: usize,
    /// Length of each side chain.
    pub n_side: usize,
    /// Segment length b.
    pub b: f64,
    /// Spacing between grafting points s.
    pub spacing: usize,
}

impl CombPolymer {
    /// Create a new comb polymer.
    pub fn new(
        n_backbone: usize,
        n_side_chains: usize,
        n_side: usize,
        b: f64,
        spacing: usize,
    ) -> Self {
        CombPolymer {
            n_backbone,
            n_side_chains,
            n_side,
            b,
            spacing,
        }
    }

    /// Total number of monomers.
    pub fn n_total(&self) -> usize {
        self.n_backbone + self.n_side_chains * self.n_side
    }

    /// Radius of gyration estimate.
    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.n_total() as f64;
        self.b * n.sqrt()
    }

    /// Backbone stretching ratio α = L_actual / L_ideal.
    pub fn backbone_stretching(&self) -> f64 {
        let n_b = self.n_backbone as f64;
        let n_s = self.n_side as f64;
        let s = self.spacing as f64;
        // Stretch increases with grafting density σ_graft ~ 1/s
        1.0 + 0.5 * n_s / (n_b.sqrt() * s)
    }
}
