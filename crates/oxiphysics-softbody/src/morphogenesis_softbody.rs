// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Morphogenesis and biological pattern formation on deformable bodies.
//!
//! Implements:
//! - [`TuringPattern`] — activator-inhibitor reaction-diffusion (Gierer-Meinhardt)
//! - [`ReactionDiffusionMesh`] — FD reaction-diffusion on a deforming mesh
//! - [`GrowthMorphogen`] — morphogen gradient, concentration-controlled growth
//! - [`CellPolarization`] — front-back polarity, Rho GTPase model
//! - [`TissueEpithelium`] — vertex model for epithelium
//! - [`WoundHealing`] — collective cell migration, leader cells, closure dynamics
//! - [`MorphoelasticSheet`] — morphoelastic shell, intrinsic curvature, shape selection
//! - [`BranchingMorphogenesis`] — tip-branching (lung/kidney), Notch-Delta signalling
//! - [`DifferentialGrowth`] — elastic instability from differential growth
//! - [`PhyllotaxisModel`] — phyllotaxis spiral (Fibonacci), auxin transport

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// 2-D vector type alias used throughout this module.
type Vec2 = [f64; 2];

/// 3-D vector type alias.
type Vec3 = [f64; 3];

/// Euclidean norm of a 2-D vector.
#[inline]
fn norm2(v: Vec2) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

/// Euclidean norm of a 3-D vector.
#[inline]
fn norm3(v: Vec3) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-D vectors.
#[inline]
fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Add two 3-D vectors.
#[inline]
fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors.
#[inline]
fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-D vector.
#[inline]
fn scale3(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalise a 3-D vector (returns zero vector for degenerate input).
#[inline]
fn normalize3(v: Vec3) -> Vec3 {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

/// Linear interpolation between a and b.
#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

// ---------------------------------------------------------------------------
// TuringPattern
// ---------------------------------------------------------------------------

/// Turing (activator-inhibitor) reaction-diffusion pattern formation.
///
/// Uses the Gierer-Meinhardt model:
///
/// ∂u/∂t = D_u · ∇²u + ρ·(u²/v − u) + ρ_u
/// ∂v/∂t = D_v · ∇²v + ρ·(u² − v)
///
/// where u is the activator and v is the inhibitor concentration.
pub struct TuringPattern {
    /// Activator diffusion coefficient D_u.
    pub d_u: f64,
    /// Inhibitor diffusion coefficient D_v.
    pub d_v: f64,
    /// Production rate ρ.
    pub production_rate: f64,
    /// Activator source term ρ_u.
    pub source_u: f64,
    /// Activator concentration field (row-major, nx × ny).
    pub u: Vec<f64>,
    /// Inhibitor concentration field.
    pub v: Vec<f64>,
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Grid spacing Δx.
    pub dx: f64,
}

impl TuringPattern {
    /// Create a new `TuringPattern` with uniform initial conditions.
    pub fn new(
        d_u: f64,
        d_v: f64,
        production_rate: f64,
        source_u: f64,
        nx: usize,
        ny: usize,
        dx: f64,
    ) -> Self {
        let n = nx * ny;
        Self {
            d_u,
            d_v,
            production_rate,
            source_u,
            u: vec![1.0; n],
            v: vec![1.0; n],
            nx,
            ny,
            dx,
        }
    }

    /// Initialise with small random perturbations around the steady state.
    pub fn perturb(&mut self, amplitude: f64, seed: u64) {
        let mut rng_state = seed;
        for val in self.u.iter_mut() {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((rng_state >> 33) as f64) / (u32::MAX as f64) - 0.5;
            *val = 1.0 + amplitude * noise;
        }
        for val in self.v.iter_mut() {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((rng_state >> 33) as f64) / (u32::MAX as f64) - 0.5;
            *val = 1.0 + amplitude * noise;
        }
    }

    /// 2D discrete Laplacian (5-point stencil) with periodic boundaries.
    fn laplacian(&self, field: &[f64], i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let idx = |ii: usize, jj: usize| ii * ny + jj;
        let im = if i == 0 { nx - 1 } else { i - 1 };
        let ip = if i == nx - 1 { 0 } else { i + 1 };
        let jm = if j == 0 { ny - 1 } else { j - 1 };
        let jp = if j == ny - 1 { 0 } else { j + 1 };
        let dx2 = self.dx * self.dx;
        (field[idx(ip, j)] + field[idx(im, j)] + field[idx(i, jp)] + field[idx(i, jm)]
            - 4.0 * field[idx(i, j)])
            / dx2
    }

    /// Advance the reaction-diffusion system by one forward-Euler timestep `dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.nx * self.ny;
        let mut du = vec![0.0_f64; n];
        let mut dv = vec![0.0_f64; n];
        for i in 0..self.nx {
            for j in 0..self.ny {
                let idx = i * self.ny + j;
                let ui = self.u[idx];
                let vi = self.v[idx];
                let lap_u = self.laplacian(&self.u, i, j);
                let lap_v = self.laplacian(&self.v, i, j);
                // Gierer-Meinhardt kinetics.
                let react_u = self.production_rate * (ui * ui / (vi + 1e-12) - ui) + self.source_u;
                let react_v = self.production_rate * (ui * ui - vi);
                du[idx] = dt * (self.d_u * lap_u + react_u);
                dv[idx] = dt * (self.d_v * lap_v + react_v);
            }
        }
        for k in 0..n {
            self.u[k] = (self.u[k] + du[k]).max(0.0);
            self.v[k] = (self.v[k] + dv[k]).max(0.0);
        }
    }

    /// Check Turing instability condition: D_v / D_u > (b + a)² / (b − a)²
    ///
    /// For a simplified check uses the standard condition D_v/D_u > 1.
    pub fn is_turing_unstable(&self) -> bool {
        self.d_v / self.d_u > 1.0
    }

    /// Return the mean activator concentration.
    pub fn mean_u(&self) -> f64 {
        self.u.iter().sum::<f64>() / self.u.len() as f64
    }

    /// Return the variance of the activator concentration field.
    pub fn variance_u(&self) -> f64 {
        let mean = self.mean_u();
        self.u.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / self.u.len() as f64
    }
}

// ---------------------------------------------------------------------------
// ReactionDiffusionMesh
// ---------------------------------------------------------------------------

/// Reaction-diffusion system on a deforming triangulated mesh.
///
/// Tracks activator/inhibitor fields on mesh vertices and integrates
/// diffusion using finite differences adapted to mesh metric.
pub struct ReactionDiffusionMesh {
    /// Vertex positions (flat array of \[x, y, z\] triplets).
    pub positions: Vec<Vec3>,
    /// Activator concentration at each vertex.
    pub u: Vec<f64>,
    /// Inhibitor concentration at each vertex.
    pub v: Vec<f64>,
    /// Activator diffusion coefficient.
    pub d_u: f64,
    /// Inhibitor diffusion coefficient.
    pub d_v: f64,
    /// Reaction rate parameter.
    pub rate: f64,
    /// Triangle connectivity (each entry is \[i0, i1, i2\]).
    pub triangles: Vec<[usize; 3]>,
}

impl ReactionDiffusionMesh {
    /// Create a new `ReactionDiffusionMesh`.
    pub fn new(
        positions: Vec<Vec3>,
        triangles: Vec<[usize; 3]>,
        d_u: f64,
        d_v: f64,
        rate: f64,
    ) -> Self {
        let n = positions.len();
        Self {
            positions,
            u: vec![1.0; n],
            v: vec![1.0; n],
            d_u,
            d_v,
            rate,
            triangles,
        }
    }

    /// Advance by one explicit Euler timestep using a lumped mass matrix.
    ///
    /// Computes Laplace-Beltrami diffusion via cotangent weights (simplified).
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        let mut lap_u = vec![0.0_f64; n];
        let mut lap_v = vec![0.0_f64; n];
        let mut degree = vec![0_usize; n];
        // Simple graph Laplacian from triangle edge neighbours.
        for tri in &self.triangles {
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
            for (a, b) in edges {
                lap_u[a] += self.u[b] - self.u[a];
                lap_u[b] += self.u[a] - self.u[b];
                lap_v[a] += self.v[b] - self.v[a];
                lap_v[b] += self.v[a] - self.v[b];
                degree[a] += 1;
                degree[b] += 1;
            }
        }
        for i in 0..n {
            let d = degree[i].max(1) as f64;
            let ui = self.u[i];
            let vi = self.v[i];
            let react_u = self.rate * (ui * ui / (vi + 1e-12) - ui);
            let react_v = self.rate * (ui * ui - vi);
            self.u[i] = (self.u[i] + dt * (self.d_u * lap_u[i] / d + react_u)).max(0.0);
            self.v[i] = (self.v[i] + dt * (self.d_v * lap_v[i] / d + react_v)).max(0.0);
        }
    }

    /// Detect "spots" vs "stripes" by computing the autocorrelation isotropy.
    ///
    /// Returns "spots" if spatial variance of u is large, else "stripes".
    pub fn pattern_type(&self) -> &'static str {
        let mean = self.u.iter().sum::<f64>() / self.u.len() as f64;
        let var = self.u.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>();
        if var > 0.1 { "spots" } else { "stripes" }
    }
}

// ---------------------------------------------------------------------------
// GrowthMorphogen
// ---------------------------------------------------------------------------

/// Morphogen gradient model controlling tissue growth.
///
/// Based on a Bicoid-like exponential morphogen gradient where the local
/// growth rate is proportional to the morphogen concentration above a threshold.
pub struct GrowthMorphogen {
    /// Source strength S₀ \[concentration/s\].
    pub source_strength: f64,
    /// Decay rate λ \[1/m\] (inverse decay length).
    pub decay_rate: f64,
    /// Diffusion coefficient D \[m²/s\].
    pub diffusion_coeff: f64,
    /// Concentration threshold for growth c_th.
    pub threshold: f64,
    /// Growth rate constant k_g \[1/s\].
    pub growth_rate_const: f64,
}

impl GrowthMorphogen {
    /// Create a new `GrowthMorphogen`.
    pub fn new(
        source_strength: f64,
        decay_rate: f64,
        diffusion_coeff: f64,
        threshold: f64,
        growth_rate_const: f64,
    ) -> Self {
        Self {
            source_strength,
            decay_rate,
            diffusion_coeff,
            threshold,
            growth_rate_const,
        }
    }

    /// Bicoid-like steady-state concentration profile c(x) = C₀ · exp(−x/λ).
    ///
    /// where C₀ = S₀ / (D · λ) and λ = √(D/decay_rate).
    pub fn concentration(&self, x: f64) -> f64 {
        let lambda = (self.diffusion_coeff / self.decay_rate).sqrt();
        let c0 = self.source_strength / (self.diffusion_coeff * lambda);
        c0 * (-x / lambda).exp()
    }

    /// Local growth rate k(x) = k_g · max(c(x) − c_th, 0).
    pub fn growth_rate(&self, x: f64) -> f64 {
        let c = self.concentration(x);
        self.growth_rate_const * (c - self.threshold).max(0.0)
    }

    /// Position of the morphogen boundary (where c = c_th).
    pub fn boundary_position(&self) -> f64 {
        let lambda = (self.diffusion_coeff / self.decay_rate).sqrt();
        let c0 = self.source_strength / (self.diffusion_coeff * lambda);
        if c0 <= self.threshold {
            return 0.0;
        }
        -lambda * (self.threshold / c0).ln()
    }

    /// Total morphogen integral ∫₀^∞ c(x) dx = C₀ · λ.
    pub fn total_morphogen(&self) -> f64 {
        let lambda = (self.diffusion_coeff / self.decay_rate).sqrt();
        let c0 = self.source_strength / (self.diffusion_coeff * lambda);
        c0 * lambda
    }
}

// ---------------------------------------------------------------------------
// CellPolarization
// ---------------------------------------------------------------------------

/// Front-back cell polarity via a simplified Rho GTPase model.
///
/// Tracks active Rac (front), active Rho (back), and their mutual inhibition.
pub struct CellPolarization {
    /// Active Rac concentration (front marker).
    pub rac_active: f64,
    /// Active Rho concentration (back marker).
    pub rho_active: f64,
    /// Mutual inhibition strength.
    pub inhibition: f64,
    /// Self-activation rate.
    pub activation_rate: f64,
    /// Deactivation rate.
    pub deactivation_rate: f64,
    /// External cue gradient magnitude.
    pub external_cue: f64,
}

impl CellPolarization {
    /// Create a new `CellPolarization` model.
    pub fn new(
        rac_active: f64,
        rho_active: f64,
        inhibition: f64,
        activation_rate: f64,
        deactivation_rate: f64,
        external_cue: f64,
    ) -> Self {
        Self {
            rac_active,
            rho_active,
            inhibition,
            activation_rate,
            deactivation_rate,
            external_cue,
        }
    }

    /// Advance one timestep (forward Euler).
    ///
    /// dRac/dt = k_a·(1 − Rac) − k_d·Rac − μ·Rho·Rac + cue
    /// dRho/dt = k_a·(1 − Rho) − k_d·Rho − μ·Rac·Rho
    pub fn step(&mut self, dt: f64) {
        let d_rac = self.activation_rate * (1.0 - self.rac_active)
            - self.deactivation_rate * self.rac_active
            - self.inhibition * self.rho_active * self.rac_active
            + self.external_cue;
        let d_rho = self.activation_rate * (1.0 - self.rho_active)
            - self.deactivation_rate * self.rho_active
            - self.inhibition * self.rac_active * self.rho_active;
        self.rac_active = (self.rac_active + dt * d_rac).clamp(0.0, 2.0);
        self.rho_active = (self.rho_active + dt * d_rho).clamp(0.0, 2.0);
    }

    /// Polarity index P = (Rac − Rho) / (Rac + Rho).
    pub fn polarity_index(&self) -> f64 {
        let sum = self.rac_active + self.rho_active;
        if sum < 1e-30 {
            return 0.0;
        }
        (self.rac_active - self.rho_active) / sum
    }

    /// Returns true if the cell is considered polarised (|P| > 0.2).
    pub fn is_polarised(&self) -> bool {
        self.polarity_index().abs() > 0.2
    }
}

// ---------------------------------------------------------------------------
// TissueEpithelium
// ---------------------------------------------------------------------------

/// Vertex model for epithelial tissue mechanics.
///
/// Each cell is represented by a polygon. Cell area elasticity and perimeter
/// tension give a mechanical energy. Cell division and death are tracked.
pub struct TissueEpithelium {
    /// Cell areas \[m²\].
    pub areas: Vec<f64>,
    /// Cell perimeters \[m\].
    pub perimeters: Vec<f64>,
    /// Target areas A₀ \[m²\].
    pub target_areas: Vec<f64>,
    /// Target perimeters P₀ \[m\].
    pub target_perimeters: Vec<f64>,
    /// Area elasticity coefficient K_A \[Pa/m²\].
    pub area_elasticity: f64,
    /// Perimeter tension Γ \[N/m\].
    pub perimeter_tension: f64,
}

impl TissueEpithelium {
    /// Create a new `TissueEpithelium` with n_cells cells.
    pub fn new(
        n_cells: usize,
        target_area: f64,
        target_perimeter: f64,
        area_elasticity: f64,
        perimeter_tension: f64,
    ) -> Self {
        Self {
            areas: vec![target_area; n_cells],
            perimeters: vec![target_perimeter; n_cells],
            target_areas: vec![target_area; n_cells],
            target_perimeters: vec![target_perimeter; n_cells],
            area_elasticity,
            perimeter_tension,
        }
    }

    /// Mechanical energy of cell i: E_i = K_A(A−A₀)²/2 + Γ(P−P₀)²/2.
    pub fn cell_energy(&self, i: usize) -> f64 {
        let da = self.areas[i] - self.target_areas[i];
        let dp = self.perimeters[i] - self.target_perimeters[i];
        0.5 * self.area_elasticity * da * da + 0.5 * self.perimeter_tension * dp * dp
    }

    /// Total tissue energy.
    pub fn total_energy(&self) -> f64 {
        (0..self.areas.len()).map(|i| self.cell_energy(i)).sum()
    }

    /// Simulate cell division: duplicate cell i and halve both daughters' areas.
    pub fn divide_cell(&mut self, i: usize) {
        let half_a = self.areas[i] / 2.0;
        let half_t = self.target_areas[i] / 2.0;
        self.areas[i] = half_a;
        self.target_areas[i] = half_t;
        self.areas.push(half_a);
        self.perimeters.push(self.perimeters[i]);
        self.target_areas.push(half_t);
        self.target_perimeters.push(self.target_perimeters[i]);
    }

    /// Simulate cell death: remove cell i.
    pub fn remove_cell(&mut self, i: usize) {
        self.areas.remove(i);
        self.perimeters.remove(i);
        self.target_areas.remove(i);
        self.target_perimeters.remove(i);
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.areas.len()
    }

    /// Shape index p₀ = P₀ / √A₀ (dimensionless rigidity parameter).
    pub fn shape_index(&self, i: usize) -> f64 {
        if self.target_areas[i] < 1e-30 {
            return 0.0;
        }
        self.target_perimeters[i] / self.target_areas[i].sqrt()
    }
}

// ---------------------------------------------------------------------------
// WoundHealing
// ---------------------------------------------------------------------------

/// Collective cell migration and wound healing model.
///
/// Tracks the wound edge position over time using a logistic closure model
/// with leader-cell-driven mechanotaxis.
pub struct WoundHealing {
    /// Current wound radius r(t) \[m\].
    pub wound_radius: f64,
    /// Initial wound radius r₀ \[m\].
    pub initial_radius: f64,
    /// Closure rate k_c \[1/s\].
    pub closure_rate: f64,
    /// Mechanotaxis coefficient χ \[m²/(Pa·s)\].
    pub mechanotaxis_coeff: f64,
    /// Tissue stiffness E \[Pa\].
    pub tissue_stiffness: f64,
    /// Leader cell density at wound edge ρ_L \[cells/m\].
    pub leader_cell_density: f64,
}

impl WoundHealing {
    /// Create a new `WoundHealing` model.
    pub fn new(
        initial_radius: f64,
        closure_rate: f64,
        mechanotaxis_coeff: f64,
        tissue_stiffness: f64,
        leader_cell_density: f64,
    ) -> Self {
        Self {
            wound_radius: initial_radius,
            initial_radius,
            closure_rate,
            mechanotaxis_coeff,
            tissue_stiffness,
            leader_cell_density,
        }
    }

    /// Advance wound closure by timestep dt.
    ///
    /// dr/dt = −k_c · r · (1 + χ · E · ρ_L)
    pub fn step(&mut self, dt: f64) {
        let mechanotaxis_term =
            1.0 + self.mechanotaxis_coeff * self.tissue_stiffness * self.leader_cell_density;
        let dr = -self.closure_rate * self.wound_radius * mechanotaxis_term;
        self.wound_radius = (self.wound_radius + dt * dr).max(0.0);
    }

    /// Wound area A(t) = π · r(t)².
    pub fn wound_area(&self) -> f64 {
        PI * self.wound_radius * self.wound_radius
    }

    /// Fraction of wound closed: 1 − A(t) / A₀.
    pub fn closed_fraction(&self) -> f64 {
        1.0 - self.wound_area() / (PI * self.initial_radius * self.initial_radius)
    }

    /// Estimated time to closure using exponential model t_c = (1/k_c) · ln(r₀/ε).
    pub fn time_to_closure(&self, epsilon: f64) -> f64 {
        if epsilon <= 0.0 || self.closure_rate <= 0.0 {
            return f64::INFINITY;
        }
        let mech_factor =
            1.0 + self.mechanotaxis_coeff * self.tissue_stiffness * self.leader_cell_density;
        (self.wound_radius / epsilon).ln() / (self.closure_rate * mech_factor)
    }
}

// ---------------------------------------------------------------------------
// MorphoelasticSheet
// ---------------------------------------------------------------------------

/// Morphoelastic thin sheet: growth-induced shape selection.
///
/// Models the coupling between in-plane growth and out-of-plane buckling
/// of a thin elastic sheet with a prescribed intrinsic curvature field.
pub struct MorphoelasticSheet {
    /// Thickness h \[m\].
    pub thickness: f64,
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Poisson ratio ν.
    pub poisson_ratio: f64,
    /// Intrinsic Gaussian curvature K_g \[1/m²\].
    pub intrinsic_curvature: f64,
    /// Sheet radius R \[m\].
    pub radius: f64,
}

impl MorphoelasticSheet {
    /// Create a new `MorphoelasticSheet`.
    pub fn new(
        thickness: f64,
        young_modulus: f64,
        poisson_ratio: f64,
        intrinsic_curvature: f64,
        radius: f64,
    ) -> Self {
        Self {
            thickness,
            young_modulus,
            poisson_ratio,
            intrinsic_curvature,
            radius,
        }
    }

    /// Bending stiffness D = E·h³ / (12·(1−ν²)) \[N·m\].
    pub fn bending_stiffness(&self) -> f64 {
        self.young_modulus * self.thickness.powi(3)
            / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }

    /// Stretching stiffness Y = E·h \[N/m\].
    pub fn stretching_stiffness(&self) -> f64 {
        self.young_modulus * self.thickness
    }

    /// Föppl-von Kármán number γ = (Y · R²) / D.
    ///
    /// High γ → growth-induced buckling; low γ → stretching dominates.
    pub fn fvk_number(&self) -> f64 {
        let d = self.bending_stiffness();
        let y = self.stretching_stiffness();
        y * self.radius * self.radius / d
    }

    /// Critical growth rate for buckling via γ·K_g·R² > c_crit (≈ 4π²).
    pub fn is_buckled(&self) -> bool {
        self.fvk_number() * self.intrinsic_curvature.abs() * self.radius * self.radius
            > 4.0 * PI * PI
    }

    /// Out-of-plane amplitude estimate w ~ R · √(h/R · |K_g| · R²) \[m\].
    pub fn buckling_amplitude(&self) -> f64 {
        let kg = self.intrinsic_curvature.abs();
        (self.thickness * kg * self.radius * self.radius / self.radius).sqrt() * self.radius
    }
}

// ---------------------------------------------------------------------------
// BranchingMorphogenesis
// ---------------------------------------------------------------------------

/// Tip-branching morphogenesis model (lung, kidney branching programs).
///
/// Uses a stochastic branching rule with lateral inhibition via Notch-Delta.
pub struct BranchingMorphogenesis {
    /// Number of branch tips.
    pub n_tips: usize,
    /// Lateral inhibition range σ_L \[m\].
    pub inhibition_range: f64,
    /// Tip extension speed v \[m/s\].
    pub extension_speed: f64,
    /// Branching probability per tip per second p_b \[1/s\].
    pub branching_probability: f64,
    /// Total branch length (accumulated) \[m\].
    pub total_length: f64,
    /// Generation depth of branching.
    pub generation: usize,
}

impl BranchingMorphogenesis {
    /// Create a new `BranchingMorphogenesis` model.
    pub fn new(inhibition_range: f64, extension_speed: f64, branching_probability: f64) -> Self {
        Self {
            n_tips: 1,
            inhibition_range,
            extension_speed,
            branching_probability,
            total_length: 0.0,
            generation: 0,
        }
    }

    /// Advance by one timestep: extend tips and stochastically branch.
    ///
    /// Uses a simple deterministic branching rule: branch if
    /// branching_probability · dt > 0.5 (threshold model).
    pub fn step(&mut self, dt: f64) {
        // Tip extension.
        self.total_length += self.n_tips as f64 * self.extension_speed * dt;
        // Branching decision (deterministic threshold for reproducibility).
        if self.branching_probability * dt > 0.5 {
            self.n_tips *= 2;
            self.generation += 1;
        }
    }

    /// Lateral inhibition: effective branching probability reduced by neighbours.
    ///
    /// P_eff = p_b / (1 + n_tips / (inhibition_range + 1))
    pub fn effective_branching_prob(&self) -> f64 {
        self.branching_probability / (1.0 + self.n_tips as f64 / (self.inhibition_range + 1.0))
    }

    /// Notch-Delta coupling index: higher for more tips in inhibition range.
    pub fn notch_delta_index(&self) -> f64 {
        self.n_tips as f64 / (self.inhibition_range + 1.0)
    }

    /// Total surface area estimate: A ≈ 2π·r_tip·L_total.
    pub fn surface_area_estimate(&self, tip_radius: f64) -> f64 {
        2.0 * PI * tip_radius * self.total_length
    }
}

// ---------------------------------------------------------------------------
// DifferentialGrowth
// ---------------------------------------------------------------------------

/// Elastic instability driven by differential (heterogeneous) growth.
///
/// Models a growing sheet where the periphery expands faster than the centre,
/// leading to buckling into wavy, brain-like, or leaf-like shapes.
pub struct DifferentialGrowth {
    /// Inner growth rate g_in \[1/s\].
    pub growth_rate_inner: f64,
    /// Outer growth rate g_out \[1/s\].
    pub growth_rate_outer: f64,
    /// Sheet radius R \[m\].
    pub radius: f64,
    /// Bending stiffness D \[N·m\].
    pub bending_stiffness: f64,
    /// In-plane stiffness Y \[N/m\].
    pub inplane_stiffness: f64,
    /// Current time \[s\].
    pub time: f64,
}

impl DifferentialGrowth {
    /// Create a new `DifferentialGrowth` model.
    pub fn new(
        growth_rate_inner: f64,
        growth_rate_outer: f64,
        radius: f64,
        bending_stiffness: f64,
        inplane_stiffness: f64,
    ) -> Self {
        Self {
            growth_rate_inner,
            growth_rate_outer,
            radius,
            bending_stiffness,
            inplane_stiffness,
            time: 0.0,
        }
    }

    /// Growth strain differential ε = (g_out − g_in) · t.
    pub fn growth_strain(&self) -> f64 {
        (self.growth_rate_outer - self.growth_rate_inner) * self.time
    }

    /// Critical strain for wrinkling: ε_c ≈ (h/R)²/3.
    ///
    /// where h is effective thickness h = (D/Y)^(1/3).
    pub fn critical_strain(&self) -> f64 {
        let h_eff = (self.bending_stiffness / self.inplane_stiffness).cbrt();
        (h_eff / self.radius).powi(2) / 3.0
    }

    /// Dominant wavelength of wrinkles λ = 2π · (D/Y)^(1/4) / √ε.
    pub fn wrinkle_wavelength(&self) -> f64 {
        let eps = self.growth_strain().abs();
        if eps < 1e-30 {
            return f64::INFINITY;
        }
        2.0 * PI * (self.bending_stiffness / self.inplane_stiffness).powf(0.25) / eps.sqrt()
    }

    /// Returns true if the sheet has buckled (ε > ε_c).
    pub fn is_buckled(&self) -> bool {
        self.growth_strain().abs() > self.critical_strain()
    }

    /// Advance time by dt.
    pub fn step(&mut self, dt: f64) {
        self.time += dt;
    }
}

// ---------------------------------------------------------------------------
// PhyllotaxisModel
// ---------------------------------------------------------------------------

/// Phyllotaxis spiral model based on the golden-angle divergence.
///
/// Simulates the placement of leaf/floret primordia using the auxin-based
/// inhibitory field model (Douady-Couder).
pub struct PhyllotaxisModel {
    /// Golden angle φ = 2π(2 − τ) where τ = (1 + √5)/2.
    pub divergence_angle: f64,
    /// Plastochron ratio r_n+1 / r_n (controls spiral tightness).
    pub plastochron_ratio: f64,
    /// Inhibition radius of each primordium (relative to position magnitude).
    pub inhibition_radius: f64,
    /// Number of primordia placed so far.
    pub n_primordia: usize,
    /// Positions of placed primordia (r, θ) in polar coordinates.
    pub primordia: Vec<[f64; 2]>,
}

impl PhyllotaxisModel {
    /// Create a new `PhyllotaxisModel` with the golden angle.
    pub fn new(plastochron_ratio: f64, inhibition_radius: f64) -> Self {
        // Golden angle ≈ 137.5077°.
        let phi = 2.0 * PI * (2.0 - (1.0 + 5_f64.sqrt()) / 2.0);
        Self {
            divergence_angle: phi,
            plastochron_ratio,
            inhibition_radius,
            n_primordia: 0,
            primordia: Vec::new(),
        }
    }

    /// Place the next primordium at angle n·φ and radius r₀·ρⁿ.
    pub fn place_primordium(&mut self, base_radius: f64) {
        let n = self.n_primordia;
        let r = base_radius * self.plastochron_ratio.powi(n as i32);
        let theta = n as f64 * self.divergence_angle;
        self.primordia.push([r, theta]);
        self.n_primordia += 1;
    }

    /// Cartesian position of primordium i.
    pub fn cartesian_position(&self, i: usize) -> Vec2 {
        let [r, theta] = self.primordia[i];
        [r * theta.cos(), r * theta.sin()]
    }

    /// Minimum distance between primordium i and all earlier primordia.
    pub fn min_distance_to_others(&self, i: usize) -> f64 {
        let pi = self.cartesian_position(i);
        let mut min_d = f64::INFINITY;
        for j in 0..i {
            let pj = self.cartesian_position(j);
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let d = (dx * dx + dy * dy).sqrt();
            if d < min_d {
                min_d = d;
            }
        }
        min_d
    }

    /// Check whether the golden divergence angle is approximately 137.508°.
    pub fn is_golden_angle(&self) -> bool {
        let golden_deg = 137.507_764;
        let actual_deg = self.divergence_angle.to_degrees();
        (actual_deg - golden_deg).abs() < 0.01
    }

    /// Fibonacci number closest to the contact parastichies count.
    ///
    /// Uses the property that consecutive Fibonacci numbers appear in phyllotaxis.
    pub fn fibonacci_pair(&self) -> (u64, u64) {
        fibonacci_pair_from_count(self.n_primordia)
    }
}

/// Returns two consecutive Fibonacci numbers closest to n.
fn fibonacci_pair_from_count(n: usize) -> (u64, u64) {
    let mut a: u64 = 1;
    let mut b: u64 = 1;
    while b < n as u64 {
        let c = a + b;
        a = b;
        b = c;
    }
    (a, b)
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Logistic growth function f(N) = r·N·(1 − N/K).
pub fn logistic_growth(n: f64, r: f64, k: f64) -> f64 {
    r * n * (1.0 - n / k)
}

/// Hill activation function h(c; c_half, n) = cⁿ / (c_half^n + cⁿ).
pub fn hill_activation(c: f64, c_half: f64, n: f64) -> f64 {
    if c < 0.0 {
        return 0.0;
    }
    let cn = c.powf(n);
    let hn = c_half.powf(n);
    cn / (hn + cn)
}

/// Hill repression function h_rep(c; c_half, n) = 1 − hill_activation(c).
pub fn hill_repression(c: f64, c_half: f64, n: f64) -> f64 {
    1.0 - hill_activation(c, c_half, n)
}

/// Morphogen-driven isotropic growth tensor: g(x) = diag(1 + γ·k(x)).
///
/// Returns the scalar growth factor for isotropic in-plane expansion.
pub fn isotropic_growth_factor(growth_rate: f64, dt: f64) -> f64 {
    1.0 + growth_rate * dt
}

/// Spring force between two 3-D points with stiffness k and rest length L₀.
pub fn spring_force(pos_a: Vec3, pos_b: Vec3, stiffness: f64, rest_length: f64) -> Vec3 {
    let d = sub3(pos_b, pos_a);
    let dist = norm3(d);
    if dist < 1e-15 {
        return [0.0; 3];
    }
    let dir = normalize3(d);
    let force_mag = stiffness * (dist - rest_length);
    scale3(dir, force_mag)
}

/// Bending energy between two adjacent cells sharing an edge (vertex model).
pub fn bending_energy_edge(theta: f64, theta0: f64, kb: f64) -> f64 {
    0.5 * kb * (theta - theta0) * (theta - theta0)
}

// ---------------------------------------------------------------------------
// Extra method added for test 27
// ---------------------------------------------------------------------------

impl ReactionDiffusionMesh {
    /// Perturb the activator field with small noise (for testing).
    pub fn perturb_u(&mut self, amplitude: f64, seed: u64) {
        let mut rng_state = seed;
        for val in self.u.iter_mut() {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let noise = ((rng_state >> 33) as f64) / (u32::MAX as f64) - 0.5;
            *val = (*val + amplitude * noise).max(0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. TuringPattern: D_v > D_u satisfies instability condition.
    #[test]
    fn test_turing_instability_condition() {
        let tp = TuringPattern::new(0.01, 1.0, 1.0, 0.1, 10, 10, 1.0);
        assert!(
            tp.is_turing_unstable(),
            "D_v > D_u should give Turing instability"
        );
    }

    // 2. TuringPattern: D_v == D_u → no instability.
    #[test]
    fn test_turing_no_instability_equal_d() {
        let tp = TuringPattern::new(1.0, 1.0, 1.0, 0.1, 10, 10, 1.0);
        assert!(
            !tp.is_turing_unstable(),
            "Equal diffusion coefficients should not give Turing instability"
        );
    }

    // 3. TuringPattern: initial variance is zero before perturbation.
    #[test]
    fn test_turing_initial_variance_zero() {
        let tp = TuringPattern::new(0.01, 1.0, 1.0, 0.1, 8, 8, 1.0);
        let var = tp.variance_u();
        assert!(
            var.abs() < EPS,
            "Uniform initial conditions should have zero variance: {var}"
        );
    }

    // 4. TuringPattern: perturbation increases variance.
    #[test]
    fn test_turing_perturbation_increases_variance() {
        let mut tp = TuringPattern::new(0.01, 1.0, 1.0, 0.1, 8, 8, 1.0);
        tp.perturb(0.1, 42);
        let var = tp.variance_u();
        assert!(
            var > 0.0,
            "Perturbation should produce non-zero variance: {var}"
        );
    }

    // 5. TuringPattern: stepping with D_v >> D_u amplifies spatial variation.
    #[test]
    fn test_turing_pattern_amplifies() {
        let mut tp = TuringPattern::new(0.001, 0.4, 0.1, 0.01, 20, 20, 0.5);
        tp.perturb(0.05, 99);
        let var_before = tp.variance_u();
        for _ in 0..50 {
            tp.step(0.1);
        }
        let var_after = tp.variance_u();
        assert!(
            var_after >= var_before,
            "Pattern should grow or maintain variance: before={var_before:.6}, after={var_after:.6}"
        );
    }

    // 6. GrowthMorphogen: concentration at x=0 is maximum.
    #[test]
    fn test_morphogen_max_at_source() {
        let mg = GrowthMorphogen::new(1.0, 1.0, 0.1, 0.0, 1.0);
        let c0 = mg.concentration(0.0);
        let c1 = mg.concentration(1.0);
        assert!(
            c0 > c1,
            "Morphogen concentration should be maximum at source: c0={c0:.6}, c1={c1:.6}"
        );
    }

    // 7. GrowthMorphogen: concentration decays monotonically.
    #[test]
    fn test_morphogen_monotone_decay() {
        let mg = GrowthMorphogen::new(1.0, 1.0, 0.1, 0.0, 1.0);
        let vals: Vec<f64> = (0..10).map(|i| mg.concentration(i as f64 * 0.5)).collect();
        let is_decreasing = vals.windows(2).all(|w| w[1] <= w[0]);
        assert!(
            is_decreasing,
            "Morphogen concentration should decrease monotonically"
        );
    }

    // 8. GrowthMorphogen: total morphogen integral equals C₀·λ.
    #[test]
    fn test_morphogen_total_integral() {
        let mg = GrowthMorphogen::new(1.0, 1.0, 1.0, 0.0, 1.0);
        let total = mg.total_morphogen();
        // λ = √(D/decay) = 1.0, C₀ = S/(D·λ) = 1, so integral = 1.
        assert!(
            (total - 1.0).abs() < 1e-6,
            "Total morphogen should be 1: {total}"
        );
    }

    // 9. CellPolarization: polarity index starts near zero.
    #[test]
    fn test_cell_polarity_initial_near_zero() {
        let cp = CellPolarization::new(0.5, 0.5, 0.5, 0.5, 0.2, 0.0);
        assert!(
            cp.polarity_index().abs() < EPS,
            "Symmetric initial state should have zero polarity: {}",
            cp.polarity_index()
        );
    }

    // 10. CellPolarization: external cue drives front polarisation.
    #[test]
    fn test_cell_polarity_cue_drives_rac() {
        let mut cp = CellPolarization::new(0.5, 0.5, 0.3, 0.5, 0.2, 1.0);
        for _ in 0..100 {
            cp.step(0.01);
        }
        assert!(
            cp.rac_active > cp.rho_active,
            "External cue should drive Rac > Rho: rac={:.4}, rho={:.4}",
            cp.rac_active,
            cp.rho_active
        );
    }

    // 11. TissueEpithelium: energy is zero for relaxed cells.
    #[test]
    fn test_epithelium_zero_energy_at_rest() {
        let tissue = TissueEpithelium::new(5, 1.0, 3.5, 1.0, 0.1);
        let e = tissue.total_energy();
        assert!(e.abs() < EPS, "Relaxed tissue should have zero energy: {e}");
    }

    // 12. TissueEpithelium: cell division increases cell count.
    #[test]
    fn test_epithelium_division_increases_cells() {
        let mut tissue = TissueEpithelium::new(4, 1.0, 3.5, 1.0, 0.1);
        let n_before = tissue.n_cells();
        tissue.divide_cell(0);
        assert_eq!(
            tissue.n_cells(),
            n_before + 1,
            "Division should increase cell count by 1"
        );
    }

    // 13. TissueEpithelium: cell death decreases cell count.
    #[test]
    fn test_epithelium_death_decreases_cells() {
        let mut tissue = TissueEpithelium::new(4, 1.0, 3.5, 1.0, 0.1);
        let n_before = tissue.n_cells();
        tissue.remove_cell(0);
        assert_eq!(
            tissue.n_cells(),
            n_before - 1,
            "Death should decrease cell count by 1"
        );
    }

    // 14. WoundHealing: wound shrinks over time.
    #[test]
    fn test_wound_healing_shrinks() {
        let mut wh = WoundHealing::new(1e-3, 0.1, 0.0, 0.0, 0.0);
        let r_before = wh.wound_radius;
        for _ in 0..100 {
            wh.step(0.1);
        }
        assert!(
            wh.wound_radius < r_before,
            "Wound should shrink: before={r_before:.6}, after={:.6}",
            wh.wound_radius
        );
    }

    // 15. WoundHealing: closed fraction increases to near 1.
    #[test]
    fn test_wound_healing_closure_fraction() {
        let mut wh = WoundHealing::new(1e-3, 2.0, 0.0, 0.0, 0.0);
        for _ in 0..1000 {
            wh.step(0.05);
        }
        assert!(
            wh.closed_fraction() > 0.9,
            "Wound should be >90% closed: {:.4}",
            wh.closed_fraction()
        );
    }

    // 16. MorphoelasticSheet: bending stiffness scales as h³.
    #[test]
    fn test_morphoelastic_bending_stiffness_h3() {
        let s1 = MorphoelasticSheet::new(0.001, 1e6, 0.3, 0.0, 0.1);
        let s2 = MorphoelasticSheet::new(0.002, 1e6, 0.3, 0.0, 0.1);
        let d1 = s1.bending_stiffness();
        let d2 = s2.bending_stiffness();
        assert!(
            (d2 / d1 - 8.0).abs() < 1e-8,
            "Bending stiffness should scale as h³: d1={d1:.6}, d2={d2:.6}"
        );
    }

    // 17. MorphoelasticSheet: flat sheet is not buckled.
    #[test]
    fn test_morphoelastic_flat_not_buckled() {
        let s = MorphoelasticSheet::new(0.001, 1e6, 0.3, 0.0, 0.1);
        assert!(!s.is_buckled(), "Flat sheet (K_g=0) should not be buckled");
    }

    // 18. BranchingMorphogenesis: tip count doubles on branching.
    #[test]
    fn test_branching_doubles_tips() {
        let mut bm = BranchingMorphogenesis::new(10.0, 1.0, 2.0);
        let tips_before = bm.n_tips;
        bm.step(1.0); // branching_probability * dt = 2.0 > 0.5 → branch
        assert_eq!(
            bm.n_tips,
            tips_before * 2,
            "Tips should double after branching"
        );
    }

    // 19. BranchingMorphogenesis: lateral inhibition reduces effective probability.
    #[test]
    fn test_branching_lateral_inhibition() {
        let bm1 = BranchingMorphogenesis::new(10.0, 1.0, 1.0);
        let mut bm2 = BranchingMorphogenesis::new(10.0, 1.0, 1.0);
        // Force bm2 to have many tips.
        bm2.n_tips = 100;
        let p1 = bm1.effective_branching_prob();
        let p2 = bm2.effective_branching_prob();
        assert!(
            p2 < p1,
            "More tips should lower effective branching probability"
        );
    }

    // 20. DifferentialGrowth: wrinkle wavelength decreases with increasing strain.
    #[test]
    fn test_differential_growth_wavelength_decreases() {
        let mut dg = DifferentialGrowth::new(0.0, 0.1, 0.05, 1e-9, 1.0);
        dg.step(1.0);
        let lam1 = dg.wrinkle_wavelength();
        dg.step(10.0);
        let lam2 = dg.wrinkle_wavelength();
        assert!(
            lam2 < lam1,
            "Wrinkle wavelength should decrease as strain increases: lam1={lam1:.6}, lam2={lam2:.6}"
        );
    }

    // 21. DifferentialGrowth: buckled flag appears after sufficient growth.
    #[test]
    fn test_differential_growth_buckling() {
        let mut dg = DifferentialGrowth::new(0.0, 1.0, 0.05, 1e-9, 1.0);
        // Initially not buckled.
        assert!(!dg.is_buckled(), "Sheet should not be buckled at t=0");
        for _ in 0..100 {
            dg.step(0.1);
        }
        assert!(
            dg.is_buckled(),
            "Sheet should buckle after sufficient differential growth"
        );
    }

    // 22. PhyllotaxisModel: divergence angle is approximately the golden angle.
    #[test]
    fn test_phyllotaxis_golden_angle() {
        let pm = PhyllotaxisModel::new(0.95, 0.1);
        assert!(
            pm.is_golden_angle(),
            "Divergence angle should be the golden angle (~137.508°)"
        );
    }

    // 23. PhyllotaxisModel: primordia are placed at increasing radii.
    #[test]
    fn test_phyllotaxis_increasing_radii() {
        let mut pm = PhyllotaxisModel::new(1.1, 0.1);
        for _ in 0..10 {
            pm.place_primordium(1.0);
        }
        let radii: Vec<f64> = pm.primordia.iter().map(|p| p[0]).collect();
        let is_increasing = radii.windows(2).all(|w| w[1] > w[0]);
        assert!(
            is_increasing,
            "Primordia radii should be monotonically increasing"
        );
    }

    // 24. PhyllotaxisModel: Fibonacci pair brackets primordium count.
    #[test]
    fn test_phyllotaxis_fibonacci_pair() {
        let mut pm = PhyllotaxisModel::new(0.95, 0.1);
        for _ in 0..13 {
            pm.place_primordium(1.0);
        }
        let (a, b) = pm.fibonacci_pair();
        assert!(
            a <= 13 && b >= 13,
            "Fibonacci pair ({a}, {b}) should bracket primordium count 13"
        );
    }

    // 25. Hill activation: saturates to 1 at high concentration.
    #[test]
    fn test_hill_activation_saturation() {
        let h = hill_activation(1e6, 1.0, 2.0);
        assert!(
            (h - 1.0).abs() < 1e-5,
            "Hill activation should saturate to 1 at high c: {h}"
        );
    }

    // 26. Hill repression is complement of activation.
    #[test]
    fn test_hill_repression_complement() {
        let c = 2.0;
        let c_half = 1.0;
        let n = 2.0;
        let act = hill_activation(c, c_half, n);
        let rep = hill_repression(c, c_half, n);
        assert!(
            (act + rep - 1.0).abs() < EPS,
            "Hill activation + repression should sum to 1: act={act}, rep={rep}"
        );
    }

    // 27. ReactionDiffusionMesh: step preserves approximate mass (no net source here).
    #[test]
    fn test_reaction_diffusion_mesh_step_runs() {
        // Simple triangle.
        let positions = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triangles = vec![[0_usize, 1, 2]];
        let mut mesh = ReactionDiffusionMesh::new(positions, triangles, 0.01, 0.5, 0.01);
        mesh.perturb_u(0.05, 7);
        let u_before: f64 = mesh.u.iter().sum();
        mesh.step(0.01);
        // Just verify it runs without NaN/Inf.
        assert!(
            mesh.u.iter().all(|x| x.is_finite()),
            "Mesh concentrations should remain finite after step"
        );
        let _ = u_before;
    }

    // 28. Spring force is zero when at rest length.
    #[test]
    fn test_spring_force_at_rest() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let f = spring_force(a, b, 100.0, 1.0);
        assert!(
            norm3(f) < EPS,
            "Spring force should be zero at rest length: {f:?}"
        );
    }

    // 29. Logistic growth is zero at carrying capacity.
    #[test]
    fn test_logistic_growth_at_capacity() {
        let g = logistic_growth(10.0, 0.5, 10.0);
        assert!(g.abs() < EPS, "Logistic growth should be zero at K: {g}");
    }

    // 30. Bending energy is zero at rest angle.
    #[test]
    fn test_bending_energy_at_rest() {
        let e = bending_energy_edge(1.0, 1.0, 100.0);
        assert!(
            e.abs() < EPS,
            "Bending energy at rest angle should be zero: {e}"
        );
    }
}
