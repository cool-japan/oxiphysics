// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Species definitions, reaction kinetics, and multi-step reactions.

/// Universal gas constant R (J/(mol*K)).
pub(crate) const R_GAS: f64 = 8.314;

/// A chemical species with transport and thermodynamic properties.
#[derive(Debug, Clone)]
pub struct Species {
    /// Name of the species (e.g., "O2", "CO2").
    pub name: String,
    /// Molecular diffusivity D (m^2/s or lattice units).
    pub diffusivity: f64,
    /// Molar mass (g/mol).
    pub molar_mass: f64,
}

impl Species {
    /// Create a new species.
    pub fn new(name: String, diffusivity: f64, molar_mass: f64) -> Self {
        Self {
            name,
            diffusivity,
            molar_mass,
        }
    }
}

// ---------------------------------------------------------------------------
// ArrheniusRate
// ---------------------------------------------------------------------------

/// Arrhenius reaction-rate model: `k(T) = A * exp(-E_a / (R * T))`.
#[derive(Debug, Clone, Copy)]
pub struct ArrheniusRate {
    /// Pre-exponential factor A (1/s or appropriate units).
    pub pre_exponential: f64,
    /// Activation energy E_a (J/mol).
    pub activation_energy: f64,
}

impl ArrheniusRate {
    /// Create a new Arrhenius rate model.
    pub fn new(pre_exponential: f64, activation_energy: f64) -> Self {
        Self {
            pre_exponential,
            activation_energy,
        }
    }

    /// Compute the rate constant at temperature `t` (K).
    pub fn rate(&self, t: f64) -> f64 {
        self.pre_exponential * (-self.activation_energy / (R_GAS * t)).exp()
    }

    /// Compute the ratio of rates at two temperatures (Arrhenius ratio).
    pub fn rate_ratio(&self, t1: f64, t2: f64) -> f64 {
        ((self.activation_energy / R_GAS) * (1.0 / t1 - 1.0 / t2)).exp()
    }
}

// ---------------------------------------------------------------------------
// BimolecularRate
// ---------------------------------------------------------------------------

/// Simple bimolecular reaction rate: `rate = k * c_a * c_b`.
#[derive(Debug, Clone, Copy)]
pub struct BimolecularRate {
    /// Rate constant k.
    pub k: f64,
}

impl BimolecularRate {
    /// Create a new bimolecular rate model.
    pub fn new(k: f64) -> Self {
        Self { k }
    }

    /// Compute the reaction rate given concentrations of A and B.
    pub fn rate(&self, c_a: f64, c_b: f64) -> f64 {
        self.k * c_a * c_b
    }
}

// ---------------------------------------------------------------------------
// Reaction descriptor for ReactiveLattice
// ---------------------------------------------------------------------------

/// A reaction with stoichiometric coefficients and a rate constant.
#[derive(Debug, Clone)]
pub struct Reaction {
    /// Stoichiometric coefficients for each species (negative = consumed).
    pub stoichiometry: Vec<f64>,
    /// Rate constant for this reaction.
    pub rate: f64,
}

// ---------------------------------------------------------------------------
// MultiStepReaction
// ---------------------------------------------------------------------------

/// A multi-step reaction chain: A -> B -> C -> ...
///
/// Each step has its own rate constant. The reaction proceeds sequentially:
/// Step 0: species\[0\] -> species\[1\] with rate k\[0\]
/// Step 1: species\[1\] -> species\[2\] with rate k\[1\], etc.
#[derive(Debug, Clone)]
pub struct MultiStepReaction {
    /// Rate constants for each step.
    pub rate_constants: Vec<f64>,
    /// Species indices involved (length = rate_constants.len() + 1).
    pub species_indices: Vec<usize>,
}

impl MultiStepReaction {
    /// Create a new multi-step reaction chain.
    ///
    /// `species_indices` has length n+1, `rate_constants` has length n.
    pub fn new(species_indices: Vec<usize>, rate_constants: Vec<f64>) -> Self {
        assert_eq!(
            species_indices.len(),
            rate_constants.len() + 1,
            "species_indices must be one longer than rate_constants"
        );
        Self {
            rate_constants,
            species_indices,
        }
    }

    /// Apply all steps of this reaction chain to a reactive lattice.
    pub fn apply(&self, concentrations: &mut [Vec<f64>], dt: f64) {
        for step in 0..self.rate_constants.len() {
            let src = self.species_indices[step];
            let dst = self.species_indices[step + 1];
            let rate_k = self.rate_constants[step];
            // Collect deltas to avoid simultaneous mutable borrows of two sub-slices.
            let deltas: Vec<f64> = concentrations[src]
                .iter()
                .map(|&c| rate_k * c * dt)
                .collect();
            for (c_src, &delta) in concentrations[src].iter_mut().zip(deltas.iter()) {
                *c_src -= delta;
            }
            for (c_dst, &delta) in concentrations[dst].iter_mut().zip(deltas.iter()) {
                *c_dst += delta;
            }
        }
    }

    /// Number of reaction steps.
    pub fn num_steps(&self) -> usize {
        self.rate_constants.len()
    }
}
