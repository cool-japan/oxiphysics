// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Multi-species Nernst-Planck LBM for ion transport in an electrolyte.

use super::{
    STANDARD_MOVING_WEIGHT, constants, drift_equilibrium, neg_gradient, stream_periodic,
    weighted_lattice_diffusivity,
};
use crate::error::{Error, Result};

/// LBM solver for ion transport in an electrolyte (Nernst-Planck equations).
///
/// Each species s obeys ∂c/∂t + ∇·(c **v** − D ∇c) = 0 with the drift
/// velocity **v** = −z·μ·∇φ and the Einstein relation D = μ·V_T, where μ is
/// the mobility per elementary charge (m²/V·s) and z the valence.  The
/// potential φ is supplied by the caller in [`SpaceChargeLbm::phi`] (gradients
/// use central differences inside the grid and one-sided differences at its
/// edges).
///
/// All species share one physical clock and the user-supplied BGK rate ω:
/// the most mobile species uses the standard D2Q5 weights (moving weight 1/6),
/// which fixes `Δt = D_lu·Δx²/(μ_max·V_T)` with `D_lu = (1/ω − ½)/3`; every
/// other species uses the moving weight `w_s = (μ_s/μ_max)/6`, which
/// reproduces its own diffusivity for the same Δt (see
/// [`SpaceChargeLbm::species_lattice_weight`]).  Species with zero mobility
/// are immobile.  Drift enters through the exponentially fitted equilibrium
/// (see the module documentation), so Boltzmann distributions in strong
/// fields (e.g. electric double layers) are reproduced exactly.  Streaming is
/// periodic.
pub struct SpaceChargeLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Number of ionic species.
    pub n_species: usize,
    /// Distribution functions per species \[species\]\[cell * 5 + q\].
    pub f: Vec<Vec<f64>>,
    /// Ion concentration per species (mol/m³).
    pub conc: Vec<Vec<f64>>,
    /// Valence of each species (signed integer).
    pub valence: Vec<f64>,
    /// Mobility of each species per elementary charge (m²/V·s).
    pub mobility: Vec<f64>,
    /// Thermal voltage (V).
    pub vt: f64,
    /// Electric potential φ (V).
    pub phi: Vec<f64>,
    /// BGK relaxation rate shared by all species.
    pub omega: f64,
    /// Cell spacing (m).
    pub dx: f64,
}

impl SpaceChargeLbm {
    /// Create a new space-charge LBM solver with uniform unit concentrations
    /// (1 mol/m³ for every species) at equilibrium.
    ///
    /// # Arguments
    /// * `nx`, `ny` - grid dimensions
    /// * `valence` - valence of each species
    /// * `mobility` - mobility of each species (m²/V·s, ≥ 0)
    /// * `temperature` - temperature (K)
    /// * `omega` - BGK relaxation rate
    /// * `dx` - cell spacing (m)
    ///
    /// # Errors
    /// * [`Error::ComponentCountMismatch`] if `valence` and `mobility` differ
    ///   in length.
    /// * [`Error::InvalidSpeciesCount`] if no species is given.
    /// * [`Error::InvalidGridDimensions`] for an empty grid.
    /// * [`Error::ParameterOutOfBounds`] for a negative/non-finite mobility, a
    ///   non-positive temperature or spacing, or ω ∉ (0, 2).
    pub fn new(
        nx: usize,
        ny: usize,
        valence: Vec<f64>,
        mobility: Vec<f64>,
        temperature: f64,
        omega: f64,
        dx: f64,
    ) -> Result<Self> {
        if valence.len() != mobility.len() {
            return Err(Error::ComponentCountMismatch {
                expected: valence.len(),
                found: mobility.len(),
                field: "mobility".to_string(),
            });
        }
        if valence.is_empty() {
            return Err(Error::InvalidSpeciesCount {
                min_species: 1,
                actual: 0,
            });
        }
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidGridDimensions { nx, ny });
        }
        for &mu in &mobility {
            if !(mu >= 0.0 && mu.is_finite()) {
                return Err(Error::ParameterOutOfBounds {
                    name: "mobility".to_string(),
                    value: mu,
                    min: 0.0,
                    max: f64::MAX,
                });
            }
        }
        let scalar_checks = [
            ("temperature", temperature, 0.0, f64::MAX, temperature > 0.0),
            ("omega", omega, 0.0, 2.0, omega > 0.0 && omega < 2.0),
            ("dx", dx, 0.0, f64::MAX, dx > 0.0),
        ];
        for (name, value, min, max, ok) in scalar_checks {
            if !(ok && value.is_finite()) {
                return Err(Error::ParameterOutOfBounds {
                    name: name.to_string(),
                    value,
                    min,
                    max,
                });
            }
        }
        let ns = valence.len();
        let n = nx * ny;
        let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
        let mut solver = Self {
            nx,
            ny,
            n_species: ns,
            f: vec![vec![0.0; n * 5]; ns],
            conc: vec![vec![1.0; n]; ns],
            valence,
            mobility,
            vt,
            phi: vec![0.0; n],
            omega,
            dx,
        };
        solver.init();
        Ok(solver)
    }

    /// Re-initialise all distributions at equilibrium from
    /// [`SpaceChargeLbm::conc`]; call after editing the concentrations.
    pub fn init(&mut self) {
        for s in 0..self.n_species {
            let weight = self.species_lattice_weight(s);
            let omega = self.omega;
            for (i, &c) in self.conc[s].iter().enumerate() {
                let feq = drift_equilibrium(c, 0.0, 0.0, omega, weight);
                self.f[s][i * 5..i * 5 + 5].copy_from_slice(&feq);
            }
        }
    }

    /// Largest species mobility (m²/V·s).
    fn max_mobility(&self) -> f64 {
        self.mobility.iter().copied().fold(0.0_f64, f64::max)
    }

    /// Physical duration of one lattice step (s); infinite if every species
    /// is immobile.
    pub fn time_step(&self) -> f64 {
        let d_max = self.max_mobility() * self.vt;
        if d_max > 0.0 {
            weighted_lattice_diffusivity(self.omega, STANDARD_MOVING_WEIGHT) * self.dx * self.dx
                / d_max
        } else {
            f64::INFINITY
        }
    }

    /// Weight of each moving D2Q5 population of species `s` on the common
    /// clock, `(μ_s/μ_max)/6` (0 for an immobile species).
    pub fn species_lattice_weight(&self, s: usize) -> f64 {
        let mu_max = self.max_mobility();
        match self.mobility.get(s) {
            Some(&mu) if mu_max > 0.0 => STANDARD_MOVING_WEIGHT * mu / mu_max,
            _ => 0.0,
        }
    }

    /// Lattice diffusivity of species `s` on the common clock (0 for an
    /// immobile species).
    pub fn species_lattice_diffusivity(&self, s: usize) -> f64 {
        weighted_lattice_diffusivity(self.omega, self.species_lattice_weight(s))
    }

    /// Return the net charge density at cell index `i` (C/m³), ρ = F·Σ z·c.
    pub fn charge_density(&self, i: usize) -> f64 {
        let mut rho = 0.0;
        for s in 0..self.n_species {
            rho += self.valence[s] * self.conc[s][i] * constants::FARADAY;
        }
        rho
    }

    /// Check electroneutrality: average |ρ| over all cells.
    pub fn mean_charge_density(&self) -> f64 {
        let n = self.nx * self.ny;
        let sum: f64 = (0..n).map(|i| self.charge_density(i).abs()).sum();
        sum / n as f64
    }

    /// Perform one Nernst-Planck step for all species.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let cells = nx * ny;
        let mut field = vec![(0.0_f64, 0.0_f64); cells];
        for y in 0..ny {
            for x in 0..nx {
                field[y * nx + x] = neg_gradient(&self.phi, nx, ny, self.dx, x, y, false);
            }
        }
        let mut f_new = vec![0.0_f64; cells * 5];
        for s in 0..self.n_species {
            let weight = self.species_lattice_weight(s);
            if weight <= 0.0 {
                continue; // immobile species
            }
            let omega = self.omega;
            // drift strength per cell a = z·E·Δx/V_T; the weak-field lattice
            // velocity is u = z·μ·E·Δt/Δx = D_lu·a
            let k = self.valence[s] * self.dx / self.vt;
            let f_s = &mut self.f[s];
            for (i, &(ex, ey)) in field.iter().enumerate() {
                let base = i * 5;
                let c: f64 = f_s[base..base + 5].iter().sum();
                let feq = drift_equilibrium(c, k * ex, k * ey, omega, weight);
                for q in 0..5_usize {
                    f_s[base + q] += omega * (feq[q] - f_s[base + q]);
                }
            }
            stream_periodic(f_s, nx, ny, &mut f_new);
            std::mem::swap(f_s, &mut f_new);
            for (i, c) in self.conc[s].iter_mut().enumerate() {
                *c = self.f[s][i * 5..i * 5 + 5].iter().sum();
            }
        }
    }
}
