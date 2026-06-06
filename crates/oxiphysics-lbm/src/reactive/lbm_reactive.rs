// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM-based reactive transport: passive scalar solver, reactive lattice,
//! species transport, and reactive LBM coupling.

use super::species::{R_GAS, Reaction, Species};

// ---------------------------------------------------------------------------
// D2Q9 constants (same ordering as entropic.rs)
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights used for passive scalar equilibrium.
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 lattice velocities: (cx, cy) for each direction.
const C: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

/// Speed of sound squared for D2Q9: cs^2 = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// ReactiveLattice
// ---------------------------------------------------------------------------

/// A multi-species reactive lattice holding concentration fields and reactions.
#[derive(Debug, Clone)]
pub struct ReactiveLattice {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Species descriptors.
    pub species: Vec<Species>,
    /// Concentration fields: `concentrations[species_idx][cell_idx]`.
    pub concentrations: Vec<Vec<f64>>,
    /// List of reactions.
    pub reactions: Vec<Reaction>,
}

impl ReactiveLattice {
    /// Create an empty reactive lattice (no species yet).
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            species: Vec::new(),
            concentrations: Vec::new(),
            reactions: Vec::new(),
        }
    }

    /// Add a species with initially zero concentration everywhere.
    pub fn add_species(&mut self, species: Species) {
        let n = self.nx * self.ny;
        self.species.push(species);
        self.concentrations.push(vec![0.0; n]);
    }

    /// Advance diffusion for all species using explicit finite-difference
    /// with periodic boundary conditions.
    pub fn advance_diffusion(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for s in 0..self.species.len() {
            let d = self.species[s].diffusivity;
            let coeff = d * dt / (dx * dx);
            let old = self.concentrations[s].clone();
            for j in 0..ny {
                for i in 0..nx {
                    let ip = (i + 1) % nx;
                    let im = (i + nx - 1) % nx;
                    let jp = (j + 1) % ny;
                    let jm = (j + ny - 1) % ny;
                    let laplacian =
                        old[j * nx + ip] + old[j * nx + im] + old[jp * nx + i] + old[jm * nx + i]
                            - 4.0 * old[j * nx + i];
                    self.concentrations[s][j * nx + i] = old[j * nx + i] + coeff * laplacian;
                }
            }
        }
    }

    /// Apply all reactions using forward Euler source terms.
    ///
    /// For each reaction, the source term for species `s` at each cell is:
    ///
    /// ```text
    /// S_s = stoich_s * rate * prod(c_i)
    /// ```
    pub fn apply_reactions(&mut self, dt: f64) {
        let n = self.nx * self.ny;
        let n_species = self.species.len();
        for rxn in &self.reactions {
            assert_eq!(rxn.stoichiometry.len(), n_species);
            for k in 0..n {
                // Compute product of all concentrations.
                let mut prod = 1.0;
                for s in 0..n_species {
                    prod *= self.concentrations[s][k];
                }
                let rate_val = rxn.rate * prod;
                for s in 0..n_species {
                    self.concentrations[s][k] += rxn.stoichiometry[s] * rate_val * dt;
                }
            }
        }
    }

    /// Total concentration of a given species across all cells.
    pub fn total_concentration(&self, species_idx: usize) -> f64 {
        self.concentrations[species_idx].iter().sum()
    }

    /// Number of species in this lattice.
    pub fn num_species(&self) -> usize {
        self.species.len()
    }

    /// Set uniform concentration for a given species.
    pub fn set_uniform(&mut self, species_idx: usize, value: f64) {
        for v in self.concentrations[species_idx].iter_mut() {
            *v = value;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_reaction_source
// ---------------------------------------------------------------------------

/// Compute the reaction source term for each species given current
/// concentrations, stoichiometric coefficients, and a rate constant.
///
/// Returns a vector of source values (one per species).
///
/// ```text
/// source_i = stoich_i * rate * prod(concentrations)
/// ```
pub fn compute_reaction_source(
    concentrations: &[f64],
    stoichiometry: &[f64],
    rate: f64,
) -> Vec<f64> {
    let prod: f64 = concentrations.iter().product();
    stoichiometry.iter().map(|&s| s * rate * prod).collect()
}

// ---------------------------------------------------------------------------
// LbmPassiveScalar
// ---------------------------------------------------------------------------

/// LBM-based passive scalar (concentration) solver on D2Q9.
///
/// Advection-diffusion is modelled by a BGK collision on the scalar
/// distribution functions `g_i`, with equilibrium:
///
///   g_eq_i = w_i * c * (1 + (e_i * u) / cs^2)
///
/// The relaxation time for the scalar is tau_s = D / cs^2 + 0.5.
pub struct LbmPassiveScalar {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Scalar distribution functions stored as `g[cell_index][direction]`.
    pub g: Vec<[f64; 9]>,
    /// Scalar diffusivity D.
    pub diffusivity: f64,
}

impl LbmPassiveScalar {
    /// Create a new passive scalar solver initialised to zero concentration.
    pub fn new(nx: usize, ny: usize, diffusivity: f64) -> Self {
        Self {
            nx,
            ny,
            g: vec![[0.0; 9]; nx * ny],
            diffusivity,
        }
    }

    /// Compute the equilibrium scalar distribution for concentration `c`
    /// and fluid velocity `u`.
    ///
    /// g_eq_i = w_i * c * (1 + (e_i * u) / cs^2)
    pub fn equilibrium(c: f64, u: [f64; 2]) -> [f64; 9] {
        let mut geq = [0.0; 9];
        for i in 0..9 {
            let cx = C[i].0 as f64;
            let cy = C[i].1 as f64;
            let eu = cx * u[0] + cy * u[1];
            geq[i] = W[i] * c * (1.0 + eu / CS2);
        }
        geq
    }

    /// Perform BGK collision of the scalar distributions.
    ///
    /// The relaxation rate is omega_s = 1 / (D / cs^2 + 0.5).
    ///
    /// `velocities` must have length `nx * ny` and contain the fluid velocity
    /// `[ux, uy]` at each cell in row-major order.
    pub fn collide(&mut self, velocities: &[[f64; 2]]) {
        let omega_s = 1.0 / (self.diffusivity / CS2 + 0.5);
        for (k, g_k) in self.g.iter_mut().enumerate() {
            let c: f64 = g_k.iter().sum();
            let u = velocities[k];
            let geq = Self::equilibrium(c, u);
            for (i, g_ki) in g_k.iter_mut().enumerate() {
                *g_ki -= omega_s * (*g_ki - geq[i]);
            }
        }
    }

    /// Perform pull-scheme streaming with periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let g_old = self.g.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = y * nx + x;
                for i in 0..9 {
                    let cx = C[i].0;
                    let cy = C[i].1;
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let k_src = sy * nx + sx;
                    self.g[k_dst][i] = g_old[k_src][i];
                }
            }
        }
    }

    /// Return the concentration at cell index `idx` (sum of all g_i).
    pub fn concentration(&self, idx: usize) -> f64 {
        self.g[idx].iter().sum()
    }

    /// Initialise concentration to a given field with zero velocity.
    pub fn initialize_concentration(&mut self, concentrations: &[f64]) {
        let n = self.nx * self.ny;
        assert_eq!(concentrations.len(), n);
        for (g_k, &c_k) in self.g.iter_mut().zip(concentrations.iter()) {
            *g_k = Self::equilibrium(c_k, [0.0, 0.0]);
        }
    }

    /// Total concentration across the entire domain.
    pub fn total_concentration(&self) -> f64 {
        let n = self.nx * self.ny;
        (0..n).map(|k| self.concentration(k)).sum()
    }

    /// Apply a source term to the scalar field (adds delta_c at each cell).
    pub fn add_source(&mut self, sources: &[f64]) {
        // Distribute source equally across all directions (weighted)
        for (k, g_k) in self.g.iter_mut().enumerate() {
            for (i, g_ki) in g_k.iter_mut().enumerate() {
                *g_ki += W[i] * sources[k];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// DamkoehlerNumber helper
// ---------------------------------------------------------------------------

/// Compute the Damkoehler number: Da = k * L^2 / D.
///
/// Da >> 1: reaction-limited regime.
/// Da << 1: diffusion-limited regime.
pub fn damkoehler_number(rate_constant: f64, length_scale: f64, diffusivity: f64) -> f64 {
    rate_constant * length_scale * length_scale / diffusivity
}

/// Compute the Peclet number for species transport: Pe = u * L / D.
pub fn peclet_number(velocity: f64, length_scale: f64, diffusivity: f64) -> f64 {
    velocity * length_scale / diffusivity
}

// ---------------------------------------------------------------------------
// SpeciesTransport
// ---------------------------------------------------------------------------

/// Multi-species transport data: diffusivities, species count, and concentration fields.
///
/// Each species has its own molecular diffusivity.  Concentrations are stored
/// per-cell as `concentrations[species_index][cell_index]`.
pub struct SpeciesTransport {
    /// Molecular diffusivities for each species (lattice units).
    pub diffusivities: Vec<f64>,
    /// Number of species.
    pub n_species: usize,
    /// Concentration fields: `concentrations[s][cell]`.
    pub concentrations: Vec<Vec<f64>>,
}

impl SpeciesTransport {
    /// Create a new `SpeciesTransport` with `n_species` species on a grid of
    /// `n_cells` cells.  All concentrations are initialised to zero.
    pub fn new(n_species: usize, n_cells: usize, diffusivities: Vec<f64>) -> Self {
        assert_eq!(
            diffusivities.len(),
            n_species,
            "diffusivities length must equal n_species"
        );
        Self {
            diffusivities,
            n_species,
            concentrations: vec![vec![0.0; n_cells]; n_species],
        }
    }

    /// Set the concentration of species `s` at cell `k` to `v`.
    pub fn set(&mut self, s: usize, k: usize, v: f64) {
        self.concentrations[s][k] = v;
    }

    /// Get the concentration of species `s` at cell `k`.
    pub fn get(&self, s: usize, k: usize) -> f64 {
        self.concentrations[s][k]
    }

    /// Return the total mass (sum of all concentrations across cells) for species `s`.
    pub fn total_mass(&self, s: usize) -> f64 {
        self.concentrations[s].iter().sum()
    }
}

// ---------------------------------------------------------------------------
// ReactionRate enum
// ---------------------------------------------------------------------------

/// Reaction-rate model variants.
#[derive(Debug, Clone, Copy)]
pub enum ReactionRate {
    /// Arrhenius rate: `k = A * exp(-Ea / (R * T))`.
    Arrhenius {
        /// Pre-exponential factor A.
        a: f64,
        /// Activation energy Ea (J/mol).
        ea: f64,
        /// Temperature T (K).
        t: f64,
    },
    /// Power-law rate: `rate = k * gamma^n`.
    PowerLaw {
        /// Rate constant k.
        k: f64,
        /// Exponent n.
        n: f64,
    },
    /// Custom constant rate value.
    Custom(f64),
}

impl ReactionRate {
    /// Evaluate the reaction rate.
    ///
    /// For `Arrhenius`, uses `R = 8.314 J/(mol·K)`.
    /// For `PowerLaw(k, n)`, returns `k` (a shear-rate argument is not provided here;
    /// call `arrhenius_rate` directly for full control).
    /// For `Custom(v)`, returns `v`.
    pub fn evaluate(&self) -> f64 {
        match *self {
            ReactionRate::Arrhenius { a, ea, t } => arrhenius_rate(a, ea, t, R_GAS),
            ReactionRate::PowerLaw { k, n: _ } => k,
            ReactionRate::Custom(v) => v,
        }
    }
}

// ---------------------------------------------------------------------------
// arrhenius_rate — free function
// ---------------------------------------------------------------------------

/// Compute the Arrhenius reaction rate constant.
///
/// ```text
/// k(T) = A * exp(-Ea / (R * T))
/// ```
///
/// # Arguments
/// * `a`  — pre-exponential factor A
/// * `ea` — activation energy Ea (J/mol)
/// * `t`  — temperature T (K)
/// * `r`  — universal gas constant R (J/(mol·K)); pass `8.314` for SI.
pub fn arrhenius_rate(a: f64, ea: f64, t: f64, r: f64) -> f64 {
    a * (-ea / (r * t)).exp()
}

// ---------------------------------------------------------------------------
// ReactiveLbm
// ---------------------------------------------------------------------------

/// LBM reactive-flow solver: couples BGK hydrodynamics with multi-species
/// advection and reaction source terms.
///
/// Species concentrations are transported on the same D2Q9 grid as the fluid.
pub struct ReactiveLbm {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Base (background) density (lattice units).
    pub base_rho: f64,
    /// Species concentration fields: `species[s][cell]`.
    pub species: Vec<Vec<f64>>,
    /// Registered reactions: `(species_index, rate_model, stoichiometric_sign)`.
    ///
    /// The stoichiometric sign is `+1.0` for production and `-1.0` for consumption.
    pub reactions: Vec<(usize, ReactionRate, f64)>,
}

impl ReactiveLbm {
    /// Create a new `ReactiveLbm` with `n_species` species on an `nx × ny` grid.
    ///
    /// All concentrations are initialised to zero.
    pub fn new(nx: usize, ny: usize, n_species: usize, base_rho: f64) -> Self {
        Self {
            nx,
            ny,
            base_rho,
            species: vec![vec![0.0; nx * ny]; n_species],
            reactions: Vec::new(),
        }
    }

    /// Register a reaction: species index `s`, rate model `rate`, stoichiometric
    /// sign `sign` (`-1` consumed, `+1` produced).
    pub fn add_reaction(&mut self, s: usize, rate: ReactionRate, sign: f64) {
        self.reactions.push((s, rate, sign));
    }

    /// Number of species.
    pub fn n_species(&self) -> usize {
        self.species.len()
    }

    /// First-order upwind advection of all species concentrations.
    ///
    /// Uses a simple 1-D upwind scheme in x and y:
    ///
    /// ```text
    /// C_new(i,j) = C(i,j)
    ///   - dt * ux_+ * (C(i,j) - C(i-1,j)) / 1
    ///   - dt * ux_- * (C(i+1,j) - C(i,j)) / 1
    ///   - dt * uy_+ * (C(i,j) - C(i,j-1)) / 1
    ///   - dt * uy_- * (C(i,j+1) - C(i,j)) / 1
    /// ```
    ///
    /// where `ux_+ = max(ux, 0)`, `ux_- = min(ux, 0)`, etc.
    /// Periodic boundary conditions are applied.
    ///
    /// `u` must have length `nx * ny`.  Each element is `[ux, uy, _]`; the
    /// z-component is ignored.
    pub fn advect_species(&mut self, u: &[[f64; 3]], dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let _n = nx * ny;
        for s in 0..self.species.len() {
            let c_old = self.species[s].clone();
            for j in 0..ny {
                for i in 0..nx {
                    let k = j * nx + i;
                    let ux = u[k][0];
                    let uy = u[k][1];
                    // Upstream neighbours with periodic wrapping.
                    let im1 = if i == 0 { nx - 1 } else { i - 1 };
                    let ip1 = if i + 1 == nx { 0 } else { i + 1 };
                    let jm1 = if j == 0 { ny - 1 } else { j - 1 };
                    let jp1 = if j + 1 == ny { 0 } else { j + 1 };
                    let c_im1 = c_old[j * nx + im1];
                    let c_ip1 = c_old[j * nx + ip1];
                    let c_jm1 = c_old[jm1 * nx + i];
                    let c_jp1 = c_old[jp1 * nx + i];
                    let c_k = c_old[k];
                    // Upwind differences.
                    let flux_x = if ux >= 0.0 {
                        ux * (c_k - c_im1)
                    } else {
                        ux * (c_ip1 - c_k)
                    };
                    let flux_y = if uy >= 0.0 {
                        uy * (c_k - c_jm1)
                    } else {
                        uy * (c_jp1 - c_k)
                    };
                    self.species[s][k] = c_k - dt * (flux_x + flux_y);
                }
            }
            // Clamp negatives produced by numerical diffusion.
            for c in self.species[s].iter_mut() {
                if *c < 0.0 {
                    *c = 0.0;
                }
            }
        }
    }

    /// Apply reaction source terms to all species concentrations.
    ///
    /// For each registered reaction `(s, rate_model, sign)`, updates:
    ///
    /// ```text
    /// C_s += sign * rate.evaluate() * dt
    /// ```
    ///
    /// applied per cell, then clamped to non-negative.
    pub fn react_species(&mut self, dt: f64) {
        let _n = self.nx * self.ny;
        for &(s, rate, sign) in &self.reactions {
            let r = rate.evaluate();
            for c in self.species[s].iter_mut() {
                *c += sign * r * dt;
                if *c < 0.0 {
                    *c = 0.0;
                }
            }
        }
    }

    /// Compute the mixture density at each cell as the sum of partial densities.
    ///
    /// Partial density of species `s` at cell `k` is `base_rho * C_s[k]`.
    pub fn compute_mixture_density(&self) -> Vec<f64> {
        let n = self.nx * self.ny;
        let mut rho_mix = vec![0.0; n];
        for s_conc in &self.species {
            for (rho_k, &c_k) in rho_mix.iter_mut().zip(s_conc.iter()) {
                *rho_k += self.base_rho * c_k;
            }
        }
        rho_mix
    }
}
