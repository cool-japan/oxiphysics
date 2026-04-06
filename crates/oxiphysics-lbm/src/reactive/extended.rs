#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended reactive models: multi-component fields, Arrhenius kinetics,
//! species LBM solver, Gray-Scott pattern formation, and reactive flow solver.

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
// MultiComponentField — coupled species + temperature field
// ---------------------------------------------------------------------------

/// A grid that stores concentrations of multiple species plus a temperature field.
///
/// Supports explicit Fickian diffusion for each species and first-order reactions
/// with Arrhenius temperature dependence.
#[allow(dead_code)]
pub struct MultiComponentField {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Number of chemical species.
    pub n_species: usize,
    /// Concentrations: `conc[s][j * nx + i]`.
    pub conc: Vec<Vec<f64>>,
    /// Temperature field: `temp[j * nx + i]`.
    pub temp: Vec<f64>,
    /// Diffusivity for each species.
    pub diffusivity: Vec<f64>,
    /// Thermal diffusivity.
    pub alpha_thermal: f64,
}

impl MultiComponentField {
    /// Create a field with uniform initial concentrations and temperature.
    #[allow(dead_code)]
    pub fn new(
        nx: usize,
        ny: usize,
        n_species: usize,
        init_conc: &[f64],
        init_temp: f64,
        diffusivity: &[f64],
        alpha_thermal: f64,
    ) -> Self {
        let n = nx * ny;
        let conc = (0..n_species).map(|s| vec![init_conc[s]; n]).collect();
        Self {
            nx,
            ny,
            n_species,
            conc,
            temp: vec![init_temp; n],
            diffusivity: diffusivity.to_vec(),
            alpha_thermal,
        }
    }

    /// Linear index for cell `(i, j)`.
    #[inline]
    #[allow(dead_code)]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Diffuse all species for one step with periodic BC.
    #[allow(dead_code)]
    pub fn diffuse_species(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        for s in 0..self.n_species {
            let d = self.diffusivity[s];
            let src = self.conc[s].clone();
            for j in 0..ny {
                for i in 0..nx {
                    let ip = (i + 1) % nx;
                    let im = (i + nx - 1) % nx;
                    let jp = (j + 1) % ny;
                    let jm = (j + ny - 1) % ny;
                    let laplacian =
                        (src[j * nx + ip] + src[j * nx + im] + src[jp * nx + i] + src[jm * nx + i]
                            - 4.0 * src[j * nx + i])
                            / (dx * dx);
                    self.conc[s][j * nx + i] += d * dt * laplacian;
                }
            }
        }
    }

    /// Diffuse temperature for one step with periodic BC.
    #[allow(dead_code)]
    pub fn diffuse_temperature(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let src = self.temp.clone();
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                let lap =
                    (src[j * nx + ip] + src[j * nx + im] + src[jp * nx + i] + src[jm * nx + i]
                        - 4.0 * src[j * nx + i])
                        / (dx * dx);
                self.temp[j * nx + i] += self.alpha_thermal * dt * lap;
            }
        }
    }

    /// Apply Arrhenius first-order consumption of species `s_reactant` at each cell.
    ///
    /// Rate = A * exp(-Ea / (R * T)) * C_s
    #[allow(clippy::too_many_arguments)]
    #[allow(dead_code)]
    pub fn apply_arrhenius_reaction(
        &mut self,
        s_reactant: usize,
        pre_exp: f64,
        ea: f64,
        r_gas: f64,
        dt: f64,
    ) {
        for idx in 0..self.nx * self.ny {
            let t = self.temp[idx];
            if t > 0.0 {
                let k = pre_exp * (-(ea / (r_gas * t))).exp();
                let c = self.conc[s_reactant][idx];
                let dc = -k * c * dt;
                self.conc[s_reactant][idx] = (c + dc).max(0.0);
            }
        }
    }

    /// Apply bimolecular reaction A + B → products.
    ///
    /// dA/dt = -k * A * B,  dB/dt = -k * A * B
    #[allow(dead_code)]
    pub fn apply_bimolecular_reaction(&mut self, s_a: usize, s_b: usize, rate_k: f64, dt: f64) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            let a = self.conc[s_a][idx];
            let b = self.conc[s_b][idx];
            let delta = rate_k * a * b * dt;
            self.conc[s_a][idx] = (a - delta).max(0.0);
            self.conc[s_b][idx] = (b - delta).max(0.0);
        }
    }

    /// Apply heat release from reaction: T += Q_rxn * k * C_A * C_B * dt.
    #[allow(clippy::too_many_arguments)]
    #[allow(dead_code)]
    pub fn apply_heat_release_bimolecular(
        &mut self,
        s_a: usize,
        s_b: usize,
        rate_k: f64,
        heat_release: f64,
        dt: f64,
    ) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            let a = self.conc[s_a][idx];
            let b = self.conc[s_b][idx];
            self.temp[idx] += heat_release * rate_k * a * b * dt;
        }
    }

    /// Compute total amount of species `s`.
    #[allow(dead_code)]
    pub fn total_concentration(&self, s: usize) -> f64 {
        self.conc[s].iter().sum()
    }

    /// Compute mean temperature.
    #[allow(dead_code)]
    pub fn mean_temperature(&self) -> f64 {
        self.temp.iter().sum::<f64>() / self.temp.len() as f64
    }

    /// Find cell index with maximum temperature.
    #[allow(dead_code)]
    pub fn hotspot_index(&self) -> usize {
        self.temp
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// ArrheniusKinetics — temperature-dependent rate with Lindemann fall-off
// ---------------------------------------------------------------------------

/// Arrhenius kinetics with optional Lindemann pressure fall-off correction.
#[allow(dead_code)]
pub struct ArrheniusKinetics {
    /// Pre-exponential factor A \[1/s or m³/mol/s for bimolecular\].
    pub pre_exp: f64,
    /// Activation energy Ea \[J/mol\].
    pub ea: f64,
    /// Universal gas constant R \[J/mol/K\].
    pub r_gas: f64,
    /// Temperature exponent β (modified Arrhenius: A * T^β * exp(-Ea/RT)).
    pub beta: f64,
}

impl ArrheniusKinetics {
    /// Create a new Arrhenius kinetics descriptor.
    #[allow(dead_code)]
    pub fn new(pre_exp: f64, ea: f64, r_gas: f64, beta: f64) -> Self {
        Self {
            pre_exp,
            ea,
            r_gas,
            beta,
        }
    }

    /// Compute rate constant k(T) = A * T^β * exp(-Ea / (R T)).
    #[allow(dead_code)]
    pub fn rate(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        self.pre_exp * temperature.powf(self.beta) * (-(self.ea / (self.r_gas * temperature))).exp()
    }

    /// Sensitivity ∂k/∂T = k(T) * (β/T + Ea/(R*T²)).
    #[allow(dead_code)]
    pub fn d_rate_d_temp(&self, temperature: f64) -> f64 {
        if temperature <= 0.0 {
            return 0.0;
        }
        let k = self.rate(temperature);
        k * (self.beta / temperature + self.ea / (self.r_gas * temperature * temperature))
    }

    /// Activation temperature (Ea / R).
    #[allow(dead_code)]
    pub fn activation_temperature(&self) -> f64 {
        self.ea / self.r_gas
    }
}

// ---------------------------------------------------------------------------
// SpeciesLbm2D — multi-species advection-diffusion on D2Q9
// ---------------------------------------------------------------------------

/// LBM advection-diffusion solver for a single chemical species on D2Q9.
///
/// The species concentration is carried by a set of distribution functions `g`.
/// The equilibrium is `g_i^eq = w_i * C * (1 + e_i·u / cs²)`.
#[allow(dead_code)]
pub struct SpeciesLbm2D {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Distribution functions: `g[j * nx + i]` = `[g_0 .. g_8]`.
    pub g: Vec<[f64; 9]>,
    /// Relaxation frequency for species transport (ω_D = 1/τ_D).
    pub omega_d: f64,
}

impl SpeciesLbm2D {
    /// Create a new species LBM solver with zero initial concentration.
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, omega_d: f64) -> Self {
        Self {
            nx,
            ny,
            g: vec![[0.0; 9]; nx * ny],
            omega_d,
        }
    }

    /// Linear index for cell `(i, j)`.
    #[inline]
    #[allow(dead_code)]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Equilibrium distribution for species: g_i^eq = w_i * C * (1 + e·u/cs²).
    #[allow(dead_code)]
    pub fn equilibrium(c: f64, ux: f64, uy: f64) -> [f64; 9] {
        let mut geq = [0.0f64; 9];
        for q in 0..9 {
            let cx = C[q].0 as f64;
            let cy = C[q].1 as f64;
            let eu = cx * ux + cy * uy;
            geq[q] = W[q] * c * (1.0 + eu / CS2);
        }
        geq
    }

    /// Compute local concentration from distributions.
    #[allow(dead_code)]
    pub fn concentration_at(&self, i: usize, j: usize) -> f64 {
        self.g[self.idx(i, j)].iter().sum()
    }

    /// Total concentration across all cells.
    #[allow(dead_code)]
    pub fn total_concentration(&self) -> f64 {
        self.g.iter().flat_map(|n| n.iter()).sum()
    }

    /// Initialize uniform concentration.
    #[allow(dead_code)]
    pub fn initialize_uniform(&mut self, c: f64, ux: f64, uy: f64) {
        let geq = Self::equilibrium(c, ux, uy);
        for node in self.g.iter_mut() {
            *node = geq;
        }
    }

    /// BGK collision step.
    #[allow(dead_code)]
    pub fn collide(&mut self, ux: &[f64], uy: &[f64]) {
        let n = self.nx * self.ny;
        for idx in 0..n {
            let c: f64 = self.g[idx].iter().sum();
            let geq = Self::equilibrium(c, ux[idx], uy[idx]);
            for q in 0..9 {
                self.g[idx][q] += self.omega_d * (geq[q] - self.g[idx][q]);
            }
        }
    }

    /// Periodic streaming step.
    #[allow(dead_code)]
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let src = self.g.clone();
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let si = ((i as i32 - C[q].0).rem_euclid(nx as i32)) as usize;
                    let sj = ((j as i32 - C[q].1).rem_euclid(ny as i32)) as usize;
                    self.g[j * nx + i][q] = src[sj * nx + si][q];
                }
            }
        }
    }

    /// Add a source term (e.g., reaction source) S(x,y) to distributions.
    ///
    /// Distributes equally: g_i += w_i * S * dt.
    #[allow(dead_code)]
    pub fn add_source(&mut self, source: &[f64], dt: f64) {
        for (idx, &s) in source.iter().enumerate() {
            for q in 0..9 {
                self.g[idx][q] += W[q] * s * dt;
            }
        }
    }

    /// Execute one advection-diffusion step.
    #[allow(dead_code)]
    pub fn step(&mut self, ux: &[f64], uy: &[f64]) {
        self.collide(ux, uy);
        self.stream();
    }
}

// ---------------------------------------------------------------------------
// TwoSpeciesReaction — A + B → P with Arrhenius kinetics
// ---------------------------------------------------------------------------

/// Coupled solver for two reacting species A and B undergoing: A + B → P.
///
/// Uses separate `SpeciesLbm2D` instances for A and B, plus explicit
/// Arrhenius reaction sourcing.
#[allow(dead_code)]
pub struct TwoSpeciesReaction {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// LBM solver for species A.
    pub species_a: SpeciesLbm2D,
    /// LBM solver for species B.
    pub species_b: SpeciesLbm2D,
    /// Arrhenius kinetics descriptor.
    pub kinetics: ArrheniusKinetics,
    /// Temperature field (uniform for isothermal, spatially varying otherwise).
    pub temperature: Vec<f64>,
}

impl TwoSpeciesReaction {
    /// Create a two-species reactive solver.
    #[allow(dead_code)]
    pub fn new(
        nx: usize,
        ny: usize,
        omega_d: f64,
        kinetics: ArrheniusKinetics,
        init_temp: f64,
    ) -> Self {
        Self {
            nx,
            ny,
            species_a: SpeciesLbm2D::new(nx, ny, omega_d),
            species_b: SpeciesLbm2D::new(nx, ny, omega_d),
            kinetics,
            temperature: vec![init_temp; nx * ny],
        }
    }

    /// Apply reaction A + B → P for one time step.
    ///
    /// Returns the total reaction source (integrated rate).
    #[allow(dead_code)]
    pub fn react(&mut self, dt: f64) -> f64 {
        let mut total_rate = 0.0f64;
        let n = self.nx * self.ny;
        let mut sa = vec![0.0f64; n];
        let mut sb = vec![0.0f64; n];
        for idx in 0..n {
            let t = self.temperature[idx];
            let k = self.kinetics.rate(t);
            let ca: f64 = self.species_a.g[idx].iter().sum();
            let cb: f64 = self.species_b.g[idx].iter().sum();
            let rate = k * ca * cb;
            sa[idx] = -rate;
            sb[idx] = -rate;
            total_rate += rate;
        }
        self.species_a.add_source(&sa, dt);
        self.species_b.add_source(&sb, dt);
        total_rate * dt
    }

    /// Execute one full step: react then advect-diffuse.
    #[allow(dead_code)]
    pub fn step(&mut self, ux: &[f64], uy: &[f64], dt: f64) -> f64 {
        let consumed = self.react(dt);
        self.species_a.step(ux, uy);
        self.species_b.step(ux, uy);
        consumed
    }

    /// Total concentration of species A.
    #[allow(dead_code)]
    pub fn total_a(&self) -> f64 {
        self.species_a.total_concentration()
    }

    /// Total concentration of species B.
    #[allow(dead_code)]
    pub fn total_b(&self) -> f64 {
        self.species_b.total_concentration()
    }
}

// ---------------------------------------------------------------------------
// ReactionDiffusionSystem — Gray-Scott model
// ---------------------------------------------------------------------------

/// Gray-Scott reaction-diffusion model.
///
/// du/dt = Du * ∇²u − u*v² + f*(1−u)
/// dv/dt = Dv * ∇²v + u*v² − (f+k)*v
///
/// Generates Turing-pattern-like structures.
#[allow(dead_code)]
pub struct GrayScott {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Activator concentration field U.
    pub u: Vec<f64>,
    /// Inhibitor concentration field V.
    pub v: Vec<f64>,
    /// Diffusivity of U.
    pub du: f64,
    /// Diffusivity of V.
    pub dv: f64,
    /// Feed rate f.
    pub feed: f64,
    /// Kill rate k.
    pub kill: f64,
}

impl GrayScott {
    /// Create a new Gray-Scott system (U=1, V=0 everywhere).
    #[allow(dead_code)]
    pub fn new(nx: usize, ny: usize, du: f64, dv: f64, feed: f64, kill: f64) -> Self {
        Self {
            nx,
            ny,
            u: vec![1.0; nx * ny],
            v: vec![0.0; nx * ny],
            du,
            dv,
            feed,
            kill,
        }
    }

    /// Seed a square patch of V at the centre.
    #[allow(dead_code)]
    pub fn seed_centre(&mut self, half_size: usize, v_init: f64) {
        let cx = self.nx / 2;
        let cy = self.ny / 2;
        for j in cy.saturating_sub(half_size)..=(cy + half_size).min(self.ny - 1) {
            for i in cx.saturating_sub(half_size)..=(cx + half_size).min(self.nx - 1) {
                self.u[j * self.nx + i] = 0.5;
                self.v[j * self.nx + i] = v_init;
            }
        }
    }

    fn laplacian(field: &[f64], nx: usize, ny: usize) -> Vec<f64> {
        let mut lap = vec![0.0f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                lap[j * nx + i] = field[j * nx + ip]
                    + field[j * nx + im]
                    + field[jp * nx + i]
                    + field[jm * nx + i]
                    - 4.0 * field[j * nx + i];
            }
        }
        lap
    }

    /// Advance by one explicit Euler step with `dt` (assumes `dx = 1`).
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let lap_u = Self::laplacian(&self.u, nx, ny);
        let lap_v = Self::laplacian(&self.v, nx, ny);
        for idx in 0..nx * ny {
            let u = self.u[idx];
            let v = self.v[idx];
            let uvv = u * v * v;
            self.u[idx] += dt * (self.du * lap_u[idx] - uvv + self.feed * (1.0 - u));
            self.v[idx] += dt * (self.dv * lap_v[idx] + uvv - (self.feed + self.kill) * v);
            // clamp to [0,1]
            self.u[idx] = self.u[idx].clamp(0.0, 1.0);
            self.v[idx] = self.v[idx].clamp(0.0, 1.0);
        }
    }

    /// Total U concentration.
    #[allow(dead_code)]
    pub fn total_u(&self) -> f64 {
        self.u.iter().sum()
    }

    /// Total V concentration.
    #[allow(dead_code)]
    pub fn total_v(&self) -> f64 {
        self.v.iter().sum()
    }

    /// Mean U.
    #[allow(dead_code)]
    pub fn mean_u(&self) -> f64 {
        self.total_u() / (self.nx * self.ny) as f64
    }
}

// ---------------------------------------------------------------------------
// ReactiveFlowSolver — coupled LBM flow + reactive species transport
// ---------------------------------------------------------------------------

/// A simplified coupled solver: LBM flow (BGK D2Q9) + passive scalar species transport.
///
/// Species are advected by the flow and react according to a first-order rate.
#[allow(dead_code)]
pub struct ReactiveFlowSolver {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// LBM populations for flow: `f[j*nx+i]`.
    pub f: Vec<[f64; 9]>,
    /// Species LBM for each species.
    pub species: Vec<SpeciesLbm2D>,
    /// BGK flow relaxation frequency.
    pub omega_flow: f64,
    /// First-order rate constants for each species.
    pub rate_constants: Vec<f64>,
}

impl ReactiveFlowSolver {
    /// Create a new reactive flow solver.
    #[allow(dead_code)]
    pub fn new(
        nx: usize,
        ny: usize,
        omega_flow: f64,
        n_species: usize,
        omega_d: f64,
        rate_constants: Vec<f64>,
    ) -> Self {
        let feq = crate::lattice::equilibrium_d2q9(1.0, 0.0, 0.0);
        Self {
            nx,
            ny,
            f: vec![feq; nx * ny],
            species: (0..n_species)
                .map(|_| SpeciesLbm2D::new(nx, ny, omega_d))
                .collect(),
            omega_flow,
            rate_constants,
        }
    }

    /// Perform one coupled step: flow collide → stream → species step with reactions.
    #[allow(dead_code)]
    pub fn step(&mut self, dt: f64) {
        use crate::lattice::{bgk_d2q9, macros_from_d2q9, stream_d2q9_periodic};
        let n = self.nx * self.ny;

        // Flow collision
        for idx in 0..n {
            let (rho, ux, uy) = macros_from_d2q9(&self.f[idx]);
            bgk_d2q9(&mut self.f[idx], rho, ux, uy, self.omega_flow);
        }
        stream_d2q9_periodic(&mut self.f, self.nx, self.ny);

        // Extract velocity field
        let ux: Vec<f64> = (0..n).map(|i| macros_from_d2q9(&self.f[i]).1).collect();
        let uy: Vec<f64> = (0..n).map(|i| macros_from_d2q9(&self.f[i]).2).collect();

        // Species transport + reaction
        for (s, sp) in self.species.iter_mut().enumerate() {
            sp.step(&ux, &uy);
            // First-order reaction: dC/dt = -k * C
            let k = self.rate_constants[s];
            let src: Vec<f64> =
                sp.g.iter()
                    .map(|node| {
                        let c: f64 = node.iter().sum();
                        -k * c
                    })
                    .collect();
            sp.add_source(&src, dt);
        }
    }

    /// Total concentration of species `s`.
    #[allow(dead_code)]
    pub fn total_species(&self, s: usize) -> f64 {
        self.species[s].total_concentration()
    }
}

// ---------------------------------------------------------------------------
// Dimensionless numbers and auxiliary functions
// ---------------------------------------------------------------------------

/// Compute the Zeldovich number: Ze = Ea * (T_ad - T_0) / (R * T_ad²).
#[allow(dead_code)]
pub fn zeldovich_number(ea: f64, t_ad: f64, t_0: f64, r_gas: f64) -> f64 {
    ea * (t_ad - t_0) / (r_gas * t_ad * t_ad)
}

/// Compute the Lewis number: Le = α / D.
#[allow(dead_code)]
pub fn lewis_number(alpha_thermal: f64, diffusivity: f64) -> f64 {
    if diffusivity < 1e-30 {
        f64::INFINITY
    } else {
        alpha_thermal / diffusivity
    }
}

/// Compute the Schmidt number: Sc = ν / D.
#[allow(dead_code)]
pub fn schmidt_number(nu: f64, diffusivity: f64) -> f64 {
    if diffusivity < 1e-30 {
        f64::INFINITY
    } else {
        nu / diffusivity
    }
}

/// Compute the reaction progress variable χ = (Y - Y_0) / (Y_eq - Y_0).
#[allow(dead_code)]
pub fn progress_variable(y: f64, y_0: f64, y_eq: f64) -> f64 {
    let denom = y_eq - y_0;
    if denom.abs() < 1e-30 {
        0.0
    } else {
        (y - y_0) / denom
    }
}

/// Compute mixture fraction Z from fuel and oxidiser mass fractions.
#[allow(dead_code)]
pub fn mixture_fraction(y_f: f64, y_ox: f64, nu_ox: f64, y_f0: f64, y_ox0: f64) -> f64 {
    (y_f / y_f0 - y_ox / (nu_ox * y_ox0) + 1.0) / (1.0 + 1.0)
}

/// Equivalence ratio Φ = (F/O)_actual / (F/O)_stoichiometric.
#[allow(dead_code)]
pub fn equivalence_ratio(y_f: f64, y_ox: f64, stoich_ratio: f64) -> f64 {
    if y_ox < 1e-30 {
        return f64::INFINITY;
    }
    (y_f / y_ox) / stoich_ratio
}
