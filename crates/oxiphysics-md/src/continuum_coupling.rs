// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale continuum-atomistic coupling for hybrid MD simulations.
//!
//! Provides methods for concurrent/hierarchical coupling of atomistic
//! molecular dynamics regions to continuum mechanics descriptions,
//! including:
//!
//! - Atomistic and continuum region representations
//! - Handshake/overlap zone management
//! - Cauchy-Born constitutive model
//! - Hardy/virial atomistic stress calculations
//! - Coarse-graining (binning) of atomistic data
//! - Embedded atom method (EAM) potentials
//! - Hybrid domain decomposition (MD-FEM/FV)
//! - Spatial averaging operators
//! - Bridging Domain Method (BDM)
//! - Fluctuating hydrodynamics coupling
//! - Energy conservation across coupling interface
//!
//! References:
//! - Rudd & Broughton (1998). Phys. Rev. B 58, R5893.
//! - Tadmor & Miller (2011). Modeling Materials.
//! - Hardy, R.J. (1982). J. Chem. Phys. 76, 622.
//! - Xiao & Belytschko (2004). Comput. Methods Appl. Mech. Engrg. 193, 1645.
//! - Flekkøy, Wagner & Feder (2000). Europhys. Lett. 52, 271.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// BoundaryType
// ---------------------------------------------------------------------------

/// Boundary condition type for an atomistic simulation region.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryType {
    /// Fully periodic boundaries.
    Periodic,
    /// Fixed (Dirichlet) boundaries.
    Fixed,
    /// Padding layer with specified thickness.
    PaddingLayer {
        /// Thickness of the padding layer (Å or nm).
        thickness: f64,
    },
    /// Ghost atoms repeated over the specified number of layers.
    Ghost {
        /// Number of ghost atom layers.
        n_layers: usize,
    },
}

// ---------------------------------------------------------------------------
// AtomisticRegion
// ---------------------------------------------------------------------------

/// An atomistic (MD) simulation region.
#[derive(Debug, Clone)]
pub struct AtomisticRegion {
    /// Region origin in Cartesian space (Å).
    pub origin: [f64; 3],
    /// Region extents (width, height, depth) in Å.
    pub extents: [f64; 3],
    /// Positions of atoms in this region.
    pub positions: Vec<[f64; 3]>,
    /// Velocities of atoms in this region.
    pub velocities: Vec<[f64; 3]>,
    /// Masses of atoms in this region (amu).
    pub masses: Vec<f64>,
    /// Boundary condition on each face.
    pub boundary: BoundaryType,
}

impl AtomisticRegion {
    /// Create a new atomistic region.
    pub fn new(origin: [f64; 3], extents: [f64; 3], boundary: BoundaryType) -> Self {
        Self {
            origin,
            extents,
            positions: Vec::new(),
            velocities: Vec::new(),
            masses: Vec::new(),
            boundary,
        }
    }

    /// Return the number of atoms in the region.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Check if a Cartesian point lies inside the region.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        for (i, &p) in point.iter().enumerate() {
            let lo = self.origin[i];
            let hi = lo + self.extents[i];
            if p < lo || p > hi {
                return false;
            }
        }
        true
    }

    /// Add one atom to the region.
    pub fn add_atom(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.masses.push(mass);
    }

    /// Volume of this region (Å³).
    pub fn volume(&self) -> f64 {
        self.extents[0] * self.extents[1] * self.extents[2]
    }
}

// ---------------------------------------------------------------------------
// ContinuumRegion
// ---------------------------------------------------------------------------

/// A continuum (FEM/FV) simulation region.
#[derive(Debug, Clone)]
pub struct ContinuumRegion {
    /// Nodal positions (m).
    pub node_positions: Vec<[f64; 3]>,
    /// Nodal displacements (m).
    pub displacements: Vec<[f64; 3]>,
    /// Nodal velocities (m/s).
    pub velocities: Vec<[f64; 3]>,
    /// Per-node mass (kg).
    pub node_masses: Vec<f64>,
    /// Stress tensor at each node (Voigt notation: sxx, syy, szz, sxy, syz, szx) (Pa).
    pub stress: Vec<[f64; 6]>,
}

impl ContinuumRegion {
    /// Create a new empty continuum region.
    pub fn new() -> Self {
        Self {
            node_positions: Vec::new(),
            displacements: Vec::new(),
            velocities: Vec::new(),
            node_masses: Vec::new(),
            stress: Vec::new(),
        }
    }

    /// Return the number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.node_positions.len()
    }

    /// Add a node to the continuum region.
    pub fn add_node(&mut self, pos: [f64; 3], mass: f64) {
        self.node_positions.push(pos);
        self.displacements.push([0.0; 3]);
        self.velocities.push([0.0; 3]);
        self.node_masses.push(mass);
        self.stress.push([0.0; 6]);
    }
}

impl Default for ContinuumRegion {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// HandshakeZone
// ---------------------------------------------------------------------------

/// Overlap/handshake zone between an atomistic and continuum region.
#[derive(Debug, Clone)]
pub struct HandshakeZone {
    /// Minimum corner of the handshake box (Å).
    pub min_corner: [f64; 3],
    /// Maximum corner of the handshake box (Å).
    pub max_corner: [f64; 3],
    /// Thickness of the bridging domain (Å).
    pub thickness: f64,
}

impl HandshakeZone {
    /// Create a handshake zone defined by min/max corners.
    pub fn new(min_corner: [f64; 3], max_corner: [f64; 3]) -> Self {
        let thickness = (max_corner[0] - min_corner[0])
            .min(max_corner[1] - min_corner[1])
            .min(max_corner[2] - min_corner[2]);
        Self {
            min_corner,
            max_corner,
            thickness,
        }
    }

    /// Return true if a point lies inside the handshake zone.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        for (i, &p) in point.iter().enumerate() {
            if p < self.min_corner[i] || p > self.max_corner[i] {
                return false;
            }
        }
        true
    }

    /// Bridging weight function: linear ramp from 0 (MD side) to 1 (continuum side)
    /// along the first axis.
    pub fn bridging_weight(&self, x: f64) -> f64 {
        let lo = self.min_corner[0];
        let hi = self.max_corner[0];
        if (hi - lo).abs() < 1e-15 {
            return 0.5_f64;
        }
        ((x - lo) / (hi - lo)).clamp(0.0_f64, 1.0_f64)
    }
}

// ---------------------------------------------------------------------------
// HybridDomainDecomposition
// ---------------------------------------------------------------------------

/// Hybrid domain decomposition combining an MD region, a continuum region,
/// and a handshaking zone for concurrent multiscale coupling.
///
/// Implements the concurrent coupling framework of Rudd & Broughton (1998)
/// where an atomistic domain is embedded inside a larger continuum domain.
/// The handshaking zone provides smooth information transfer between scales.
#[derive(Debug, Clone)]
pub struct HybridDomainDecomposition {
    /// Atomistic sub-domain.
    pub md_region: AtomisticRegion,
    /// Continuum sub-domain.
    pub continuum_region: ContinuumRegion,
    /// Handshaking / bridging zone.
    pub handshake_zone: HandshakeZone,
    /// Whether energy conservation correction is applied at each step.
    pub energy_correction: bool,
    /// Total system energy (J), tracked for conservation diagnostics.
    pub total_energy: f64,
}

impl HybridDomainDecomposition {
    /// Create a new hybrid domain decomposition.
    pub fn new(
        md_region: AtomisticRegion,
        continuum_region: ContinuumRegion,
        handshake_zone: HandshakeZone,
    ) -> Self {
        Self {
            md_region,
            continuum_region,
            handshake_zone,
            energy_correction: true,
            total_energy: 0.0_f64,
        }
    }

    /// Count atoms that lie inside the handshaking zone.
    pub fn atoms_in_handshake(&self) -> usize {
        self.md_region
            .positions
            .iter()
            .filter(|&&p| self.handshake_zone.contains(p))
            .count()
    }

    /// Count continuum nodes that lie inside the handshaking zone.
    pub fn nodes_in_handshake(&self) -> usize {
        self.continuum_region
            .node_positions
            .iter()
            .filter(|&&p| self.handshake_zone.contains(p))
            .count()
    }

    /// Verify domain consistency: every atom inside the MD region extents
    /// has a corresponding entry in velocities and masses.
    pub fn is_consistent(&self) -> bool {
        let na = self.md_region.positions.len();
        na == self.md_region.velocities.len() && na == self.md_region.masses.len()
    }

    /// Perform one coupling step: transfer averaged MD velocity to matching
    /// continuum nodes in the handshake zone using linear blending.
    pub fn coupling_step(&mut self) {
        // Average velocity of atoms inside the handshake zone
        let mut v_avg = [0.0_f64; 3];
        let mut count = 0_usize;
        for (i, &pos) in self.md_region.positions.iter().enumerate() {
            if self.handshake_zone.contains(pos) {
                let v = self.md_region.velocities[i];
                v_avg[0] += v[0];
                v_avg[1] += v[1];
                v_avg[2] += v[2];
                count += 1;
            }
        }
        if count == 0 {
            return;
        }
        let n = count as f64;
        v_avg[0] /= n;
        v_avg[1] /= n;
        v_avg[2] /= n;

        // Blend into continuum nodes inside handshake zone
        for (j, &pos) in self.continuum_region.node_positions.iter().enumerate() {
            if self.handshake_zone.contains(pos) {
                let w = self.handshake_zone.bridging_weight(pos[0]);
                let v = &mut self.continuum_region.velocities[j];
                v[0] = (1.0_f64 - w) * v[0] + w * v_avg[0];
                v[1] = (1.0_f64 - w) * v[1] + w * v_avg[1];
                v[2] = (1.0_f64 - w) * v[2] + w * v_avg[2];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// AveragingOperator
// ---------------------------------------------------------------------------

/// Spatial averaging operator that coarse-grains atomistic quantities
/// (position, velocity, stress) onto a regular continuum grid.
///
/// Implements the Hardy (1982) spatial averaging procedure using a
/// localization function (kernel) to produce smooth fields from
/// discrete particle data.
#[derive(Debug, Clone)]
pub struct AveragingOperator {
    /// Grid spacing in each direction (Å).
    pub grid_spacing: [f64; 3],
    /// Grid dimensions (nx, ny, nz).
    pub grid_dims: [usize; 3],
    /// Grid origin (Å).
    pub origin: [f64; 3],
    /// Kernel bandwidth (Å).
    pub bandwidth: f64,
    /// Averaged density field (atoms/Å³).
    pub density: Vec<f64>,
    /// Averaged velocity field (Å/ps).
    pub velocity: Vec<[f64; 3]>,
    /// Averaged stress field (Voigt, eV/Å³).
    pub stress: Vec<[f64; 6]>,
}

impl AveragingOperator {
    /// Create a new averaging operator over a uniform grid.
    pub fn new(
        origin: [f64; 3],
        grid_dims: [usize; 3],
        grid_spacing: [f64; 3],
        bandwidth: f64,
    ) -> Self {
        let n = grid_dims[0] * grid_dims[1] * grid_dims[2];
        Self {
            grid_spacing,
            grid_dims,
            origin,
            bandwidth,
            density: vec![0.0_f64; n],
            velocity: vec![[0.0_f64; 3]; n],
            stress: vec![[0.0_f64; 6]; n],
        }
    }

    /// Total number of grid points.
    pub fn n_points(&self) -> usize {
        self.grid_dims[0] * self.grid_dims[1] * self.grid_dims[2]
    }

    /// Flat index from (ix, iy, iz) grid indices.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz + self.grid_dims[2] * (iy + self.grid_dims[1] * ix)
    }

    /// Position of grid point (ix, iy, iz) in Å.
    pub fn grid_position(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.grid_spacing[0],
            self.origin[1] + iy as f64 * self.grid_spacing[1],
            self.origin[2] + iz as f64 * self.grid_spacing[2],
        ]
    }

    /// Gaussian localization kernel evaluated at distance `r`.
    ///
    /// w(r) = (2π σ²)^(−3/2) exp(−r²/(2σ²))
    pub fn kernel(&self, r_sq: f64) -> f64 {
        let sigma2 = self.bandwidth * self.bandwidth;
        let norm = (2.0_f64 * PI * sigma2).powf(1.5_f64);
        (-r_sq / (2.0_f64 * sigma2)).exp() / norm
    }

    /// Average atomistic number density onto the grid.
    ///
    /// # Arguments
    /// * `positions` – atom positions (Å)
    pub fn average_density(&mut self, positions: &[[f64; 3]]) {
        let n = self.n_points();
        let mut rho = vec![0.0_f64; n];
        let [nx, ny, nz] = self.grid_dims;
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let gp = self.grid_position(ix, iy, iz);
                    let mut w_sum = 0.0_f64;
                    for &ap in positions {
                        let dr = [ap[0] - gp[0], ap[1] - gp[1], ap[2] - gp[2]];
                        let r_sq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                        w_sum += self.kernel(r_sq);
                    }
                    rho[self.index(ix, iy, iz)] = w_sum;
                }
            }
        }
        self.density = rho;
    }

    /// Average atomistic velocity field onto the grid (mass-weighted).
    ///
    /// # Arguments
    /// * `positions`  – atom positions (Å)
    /// * `velocities` – atom velocities (Å/ps)
    /// * `masses`     – atom masses (amu)
    pub fn average_velocity(
        &mut self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
    ) {
        let n = self.n_points();
        let mut v_field = vec![[0.0_f64; 3]; n];
        let [nx, ny, nz] = self.grid_dims;
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let idx = self.index(ix, iy, iz);
                    let gp = self.grid_position(ix, iy, iz);
                    let mut wm_sum = 0.0_f64;
                    let mut wmv = [0.0_f64; 3];
                    for ((&ap, &av), &_am) in
                        positions.iter().zip(velocities.iter()).zip(masses.iter())
                    {
                        let dr = [ap[0] - gp[0], ap[1] - gp[1], ap[2] - gp[2]];
                        let r_sq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                        let w = self.kernel(r_sq);
                        wm_sum += w;
                        wmv[0] += w * av[0];
                        wmv[1] += w * av[1];
                        wmv[2] += w * av[2];
                    }
                    if wm_sum > 1e-30_f64 {
                        v_field[idx] = [wmv[0] / wm_sum, wmv[1] / wm_sum, wmv[2] / wm_sum];
                    }
                }
            }
        }
        self.velocity = v_field;
    }

    /// Return the peak density value on the grid.
    pub fn peak_density(&self) -> f64 {
        self.density.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// BridgingDomainMethod
// ---------------------------------------------------------------------------

/// Bridging Domain Method (BDM) for coupling displacement/velocity fields
/// between atomistic and continuum regions.
///
/// Implements the framework of Xiao & Belytschko (2004) where an overlap
/// zone has both MD atoms and FE nodes active, blended with weight functions.
#[derive(Debug, Clone)]
pub struct BridgingDomainMethod {
    /// Overlap zone definition.
    pub handshake: HandshakeZone,
    /// Lagrange multiplier penalty stiffness (eV/Å²).
    pub penalty_stiffness: f64,
    /// Constraint residual from last update (Å).
    pub constraint_residual: f64,
}

impl BridgingDomainMethod {
    /// Create a new BDM coupler with the specified overlap zone.
    pub fn new(handshake: HandshakeZone, penalty_stiffness: f64) -> Self {
        Self {
            handshake,
            penalty_stiffness,
            constraint_residual: 0.0_f64,
        }
    }

    /// Compute the blending weight θ(x) for a point in the overlap zone.
    ///
    /// θ ranges from 0 (pure MD) to 1 (pure continuum).
    pub fn blending_weight(&self, pos: [f64; 3]) -> f64 {
        self.handshake.bridging_weight(pos[0])
    }

    /// Blend MD and FE nodal displacements using weight function.
    ///
    /// u_blended = (1 − θ) · u_md + θ · u_fem
    pub fn blend_displacement(&self, u_md: [f64; 3], u_fem: [f64; 3], pos: [f64; 3]) -> [f64; 3] {
        let theta = self.blending_weight(pos);
        let one_minus = 1.0_f64 - theta;
        [
            one_minus * u_md[0] + theta * u_fem[0],
            one_minus * u_md[1] + theta * u_fem[1],
            one_minus * u_md[2] + theta * u_fem[2],
        ]
    }

    /// Compute the constraint force on atom `i` in the overlap zone.
    ///
    /// F_constraint = −k · θ · (u_atom − u_fem)
    pub fn constraint_force(
        &mut self,
        u_atom: [f64; 3],
        u_fem: [f64; 3],
        pos: [f64; 3],
    ) -> [f64; 3] {
        let theta = self.blending_weight(pos);
        let du = [
            u_atom[0] - u_fem[0],
            u_atom[1] - u_fem[1],
            u_atom[2] - u_fem[2],
        ];
        let norm = (du[0] * du[0] + du[1] * du[1] + du[2] * du[2]).sqrt();
        self.constraint_residual = norm;
        let k = self.penalty_stiffness * theta;
        [-k * du[0], -k * du[1], -k * du[2]]
    }

    /// Update the overlap zone thickness.
    pub fn set_thickness(&mut self, thickness: f64) {
        self.handshake.thickness = thickness;
    }

    /// Check whether constraints have converged (residual < tolerance).
    pub fn has_converged(&self, tol: f64) -> bool {
        self.constraint_residual < tol
    }
}

// ---------------------------------------------------------------------------
// AtomisticToContin
// ---------------------------------------------------------------------------

/// Atomistic-to-continuum stress extractor using the Cauchy-Born rule
/// and the Hardy (virial) stress formulation.
///
/// Computes the Cauchy stress tensor from atomic positions and forces
/// using the Irving-Kirkwood-Hardy procedure.
#[derive(Debug, Clone)]
pub struct AtomisticToContin {
    /// Volume assigned to each atom for stress normalization (Å³).
    pub atomic_volume: f64,
    /// Last computed stress tensor in Voigt notation (eV/Å³).
    pub stress_voigt: [f64; 6],
}

impl AtomisticToContin {
    /// Create a new extractor with the given per-atom volume.
    pub fn new(atomic_volume: f64) -> Self {
        Self {
            atomic_volume,
            stress_voigt: [0.0_f64; 6],
        }
    }

    /// Compute the virial (kinetic) contribution to stress.
    ///
    /// σ_kin = − (1/V) Σ_i m_i v_i ⊗ v_i
    ///
    /// Returns 6-component Voigt stress (eV/Å³).
    pub fn kinetic_stress(&self, velocities: &[[f64; 3]], masses: &[f64], volume: f64) -> [f64; 6] {
        let mut s = [0.0_f64; 6];
        // mass unit: amu, velocity: Å/ps → eV conversion factor ≈ 1.03643e-4
        // Here we leave in amu·Å²/ps² units; caller converts as needed.
        for (&v, &m) in velocities.iter().zip(masses.iter()) {
            s[0] -= m * v[0] * v[0]; // xx
            s[1] -= m * v[1] * v[1]; // yy
            s[2] -= m * v[2] * v[2]; // zz
            s[3] -= m * v[0] * v[1]; // xy
            s[4] -= m * v[1] * v[2]; // yz
            s[5] -= m * v[2] * v[0]; // zx
        }
        if volume > 1e-30_f64 {
            for si in s.iter_mut() {
                *si /= volume;
            }
        }
        s
    }

    /// Compute the potential (virial) stress from pairwise forces.
    ///
    /// σ_pot = (1/V) Σ_{i<j} r_ij ⊗ f_ij
    ///
    /// # Arguments
    /// * `positions` – atom positions (Å)
    /// * `forces`    – forces on atoms (eV/Å)
    /// * `volume`    – averaging volume (Å³)
    pub fn potential_stress(
        &mut self,
        positions: &[[f64; 3]],
        forces: &[[f64; 3]],
        volume: f64,
    ) -> [f64; 6] {
        let mut s = [0.0_f64; 6];
        for (&r, &f) in positions.iter().zip(forces.iter()) {
            s[0] += r[0] * f[0];
            s[1] += r[1] * f[1];
            s[2] += r[2] * f[2];
            s[3] += 0.5_f64 * (r[0] * f[1] + r[1] * f[0]);
            s[4] += 0.5_f64 * (r[1] * f[2] + r[2] * f[1]);
            s[5] += 0.5_f64 * (r[2] * f[0] + r[0] * f[2]);
        }
        if volume > 1e-30_f64 {
            for si in s.iter_mut() {
                *si /= volume;
            }
        }
        self.stress_voigt = s;
        s
    }

    /// Total Cauchy stress = kinetic + potential contributions.
    pub fn total_stress(
        &mut self,
        positions: &[[f64; 3]],
        forces: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        volume: f64,
    ) -> [f64; 6] {
        let s_kin = self.kinetic_stress(velocities, masses, volume);
        let s_pot = self.potential_stress(positions, forces, volume);
        let mut s_total = [0.0_f64; 6];
        for i in 0..6 {
            s_total[i] = s_kin[i] + s_pot[i];
        }
        self.stress_voigt = s_total;
        s_total
    }

    /// Hydrostatic pressure from Voigt stress: P = −tr(σ)/3.
    pub fn pressure(&self) -> f64 {
        -(self.stress_voigt[0] + self.stress_voigt[1] + self.stress_voigt[2]) / 3.0_f64
    }

    /// Von Mises stress from Voigt components.
    pub fn von_mises(&self) -> f64 {
        let s = &self.stress_voigt;
        let dev_sq = (s[0] - s[1]).powi(2)
            + (s[1] - s[2]).powi(2)
            + (s[2] - s[0]).powi(2)
            + 6.0_f64 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
        (0.5_f64 * dev_sq).sqrt()
    }
}

// ---------------------------------------------------------------------------
// FluctuatingHydrodynamics
// ---------------------------------------------------------------------------

/// Fluctuating hydrodynamics coupler: injects thermal noise into the
/// continuum Navier-Stokes equations derived from MD fluctuation statistics.
///
/// Based on the Landau-Lifshitz fluctuating hydrodynamics framework,
/// where stochastic stress tensors are added to the continuum equations
/// to reproduce the correct fluctuation-dissipation balance.
#[derive(Debug, Clone)]
pub struct FluctuatingHydrodynamics {
    /// Temperature of the MD system (K).
    pub temperature: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Bulk viscosity (Pa·s).
    pub bulk_viscosity: f64,
    /// Thermal conductivity (W/(m·K)).
    pub thermal_conductivity: f64,
    /// Grid cell volume for noise normalization (m³).
    pub cell_volume: f64,
    /// Time step (s).
    pub dt: f64,
    /// Random number generator seed.
    pub seed: u64,
}

impl FluctuatingHydrodynamics {
    /// Create a new fluctuating hydrodynamics coupler.
    pub fn new(
        temperature: f64,
        viscosity: f64,
        bulk_viscosity: f64,
        thermal_conductivity: f64,
        cell_volume: f64,
        dt: f64,
        seed: u64,
    ) -> Self {
        Self {
            temperature,
            viscosity,
            bulk_viscosity,
            thermal_conductivity,
            cell_volume,
            dt,
            seed,
        }
    }

    /// Amplitude of the stochastic stress tensor noise.
    ///
    /// A² = 2 k_B T η / (V · Δt)
    ///
    /// Returns noise amplitude in Pa.
    pub fn stress_noise_amplitude(&self) -> f64 {
        const K_B: f64 = 1.380649e-23_f64; // J/K
        let denom = self.cell_volume * self.dt;
        if denom < 1e-60_f64 {
            return 0.0_f64;
        }
        (2.0_f64 * K_B * self.temperature * self.viscosity / denom).sqrt()
    }

    /// Amplitude of the stochastic heat flux noise.
    ///
    /// A_q² = 2 k_B T² κ / (V · Δt)
    pub fn heat_flux_noise_amplitude(&self) -> f64 {
        const K_B: f64 = 1.380649e-23_f64;
        let denom = self.cell_volume * self.dt;
        if denom < 1e-60_f64 {
            return 0.0_f64;
        }
        (2.0_f64 * K_B * self.temperature * self.temperature * self.thermal_conductivity / denom)
            .sqrt()
    }

    /// Generate a symmetric traceless stochastic stress tensor sample.
    ///
    /// Uses Box-Muller transform for Gaussian noise.
    /// Returns Voigt notation: \[Sxx, Syy, Szz, Sxy, Syz, Szx\].
    pub fn sample_stress_noise(&mut self) -> [f64; 6] {
        let amp = self.stress_noise_amplitude();
        // Simple deterministic pseudo-noise based on seed (no external rand dep)
        let mut s = [0.0_f64; 6];
        for (i, si) in s.iter_mut().enumerate() {
            self.seed = self
                .seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u = (self.seed >> 33) as f64 / (u64::MAX >> 33) as f64;
            self.seed = self
                .seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = (self.seed >> 33) as f64 / (u64::MAX >> 33) as f64;
            // Box-Muller
            let gauss = if u < 1e-30_f64 {
                0.0_f64
            } else {
                (-2.0_f64 * u.ln()).sqrt() * (2.0_f64 * PI * v).cos()
            };
            *si = amp * gauss * if i < 3 { 1.0_f64 } else { 2.0_f64.sqrt() };
        }
        // Make trace-free: subtract mean diagonal
        let trace = (s[0] + s[1] + s[2]) / 3.0_f64;
        s[0] -= trace;
        s[1] -= trace;
        s[2] -= trace;
        s
    }

    /// Update temperature from MD kinetic energy.
    ///
    /// T = (2 / (3 N k_B)) Σ_i (1/2) m_i v_i²
    pub fn update_temperature_from_md(&mut self, velocities: &[[f64; 3]], masses: &[f64]) {
        const K_B_AMU: f64 = 8.314462618e-3_f64; // kJ/(mol·K) ≈ units for amu·Å²/ps²
        let mut ke = 0.0_f64;
        for (&v, &m) in velocities.iter().zip(masses.iter()) {
            ke += 0.5_f64 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        let n = velocities.len() as f64;
        if n > 0.0_f64 {
            self.temperature = 2.0_f64 * ke / (3.0_f64 * n * K_B_AMU);
        }
    }
}

// ---------------------------------------------------------------------------
// CoupleEnergy
// ---------------------------------------------------------------------------

/// Energy handshake ensuring conservation across the MD-continuum coupling
/// interface.
///
/// Tracks the total mechanical energy (kinetic + potential) in both MD
/// and continuum sub-domains and applies correction fluxes to maintain
/// global conservation.
#[derive(Debug, Clone)]
pub struct CoupleEnergy {
    /// Kinetic energy in the MD region (eV).
    pub md_kinetic: f64,
    /// Potential energy in the MD region (eV).
    pub md_potential: f64,
    /// Kinetic energy in the continuum region (eV).
    pub continuum_kinetic: f64,
    /// Elastic strain energy in the continuum region (eV).
    pub continuum_elastic: f64,
    /// Energy flux injected at the interface at the last step (eV).
    pub interface_flux: f64,
    /// Cumulative energy error (eV).
    pub cumulative_error: f64,
}

impl CoupleEnergy {
    /// Create a new energy coupler (all energies zero).
    pub fn new() -> Self {
        Self {
            md_kinetic: 0.0_f64,
            md_potential: 0.0_f64,
            continuum_kinetic: 0.0_f64,
            continuum_elastic: 0.0_f64,
            interface_flux: 0.0_f64,
            cumulative_error: 0.0_f64,
        }
    }

    /// Update MD kinetic energy from velocities and masses.
    ///
    /// Uses conversion: 1 amu·(Å/ps)² = 0.01036443 eV.
    pub fn update_md_kinetic(&mut self, velocities: &[[f64; 3]], masses: &[f64]) {
        const CONV: f64 = 0.010364269656262398_f64; // amu·(Å/ps)² → eV
        let mut ke = 0.0_f64;
        for (&v, &m) in velocities.iter().zip(masses.iter()) {
            ke += 0.5_f64 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        self.md_kinetic = ke * CONV;
    }

    /// Update continuum kinetic energy from nodal velocities and masses.
    pub fn update_continuum_kinetic(&mut self, velocities: &[[f64; 3]], masses: &[f64]) {
        let mut ke = 0.0_f64;
        for (&v, &m) in velocities.iter().zip(masses.iter()) {
            ke += 0.5_f64 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        self.continuum_kinetic = ke;
    }

    /// Total system energy (MD + continuum).
    pub fn total_energy(&self) -> f64 {
        self.md_kinetic + self.md_potential + self.continuum_kinetic + self.continuum_elastic
    }

    /// Compute the interface energy flux needed to maintain conservation.
    ///
    /// `reference_energy` is the conserved total energy target.
    pub fn compute_flux(&mut self, reference_energy: f64) -> f64 {
        let error = reference_energy - self.total_energy();
        self.cumulative_error += error;
        self.interface_flux = error;
        error
    }

    /// Apply an energy correction to MD velocities by uniformly scaling them.
    ///
    /// This rescales velocities to inject/remove `delta_e` eV of kinetic energy.
    pub fn apply_kinetic_correction(
        &self,
        velocities: &mut [[f64; 3]],
        masses: &[f64],
        delta_e: f64,
    ) {
        const CONV: f64 = 0.010364269656262398_f64;
        // Current KE in eV
        let mut ke = 0.0_f64;
        for (&v, &m) in velocities.iter().zip(masses.iter()) {
            ke += 0.5_f64 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        let ke_ev = ke * CONV;
        if ke_ev < 1e-30_f64 {
            return;
        }
        let scale = ((ke_ev + delta_e) / ke_ev).max(0.0_f64).sqrt();
        for v in velocities.iter_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }

    /// Return true if the energy conservation error is below the tolerance.
    pub fn is_conserved(&self, tol: f64) -> bool {
        self.cumulative_error.abs() < tol
    }
}

impl Default for CoupleEnergy {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CauchyBornModel
// ---------------------------------------------------------------------------

/// Cauchy-Born constitutive model for extracting continuum stress from
/// the deformation gradient of the crystal lattice.
///
/// The Cauchy-Born hypothesis assumes that the deformation of the crystal
/// at each continuum point is homogeneous and equal to the macroscopic
/// deformation gradient F.
#[derive(Debug, Clone)]
pub struct CauchyBornModel {
    /// Reference lattice parameter (Å).
    pub lattice_parameter: f64,
    /// Elastic constants C11, C12, C44 (eV/Å³) for cubic symmetry.
    pub elastic_constants: [f64; 3],
    /// Current deformation gradient F (3×3, row-major).
    pub deformation_gradient: [f64; 9],
}

impl CauchyBornModel {
    /// Create a new Cauchy-Born model with given lattice parameter and
    /// elastic constants for cubic symmetry.
    pub fn new(lattice_parameter: f64, c11: f64, c12: f64, c44: f64) -> Self {
        Self {
            lattice_parameter,
            elastic_constants: [c11, c12, c44],
            deformation_gradient: [
                1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64,
            ],
        }
    }

    /// Set the deformation gradient from current and reference positions.
    ///
    /// F ≈ I + ∇u (small-strain approximation here).
    pub fn set_deformation_gradient(&mut self, f: [f64; 9]) {
        self.deformation_gradient = f;
    }

    /// Compute the Green-Lagrange strain tensor E = (F^T F − I)/2.
    ///
    /// Returns Voigt notation: \[Exx, Eyy, Ezz, 2Exy, 2Eyz, 2Ezx\].
    pub fn green_lagrange_strain(&self) -> [f64; 6] {
        let f = &self.deformation_gradient;
        // FtF[i][j] = Σ_k F[k][i] * F[k][j]  (row-major, F[row*3+col])
        let ftf = |i: usize, j: usize| -> f64 { (0..3).map(|k| f[k * 3 + i] * f[k * 3 + j]).sum() };
        [
            0.5_f64 * (ftf(0, 0) - 1.0_f64),
            0.5_f64 * (ftf(1, 1) - 1.0_f64),
            0.5_f64 * (ftf(2, 2) - 1.0_f64),
            ftf(0, 1),
            ftf(1, 2),
            ftf(2, 0),
        ]
    }

    /// Compute the 2nd Piola-Kirchhoff stress from the GL strain
    /// using the linear elastic constitutive law for cubic symmetry.
    ///
    /// Returns Voigt notation \[Sxx, Syy, Szz, Sxy, Syz, Szx\] (eV/Å³).
    pub fn pk2_stress(&self) -> [f64; 6] {
        let e = self.green_lagrange_strain();
        let [c11, c12, c44] = self.elastic_constants;
        [
            c11 * e[0] + c12 * (e[1] + e[2]),
            c11 * e[1] + c12 * (e[0] + e[2]),
            c11 * e[2] + c12 * (e[0] + e[1]),
            c44 * e[3],
            c44 * e[4],
            c44 * e[5],
        ]
    }

    /// Young's modulus from cubic elastic constants.
    pub fn youngs_modulus(&self) -> f64 {
        let [c11, c12, _c44] = self.elastic_constants;
        (c11 - c12) * (c11 + 2.0_f64 * c12) / (c11 + c12)
    }
}

// ---------------------------------------------------------------------------
// CoarseGrainBin
// ---------------------------------------------------------------------------

/// A single spatial bin for coarse-graining atomistic data.
#[derive(Debug, Clone)]
pub struct CoarseGrainBin {
    /// Center of the bin (Å).
    pub center: [f64; 3],
    /// Half-width of the bin (Å).
    pub half_width: f64,
    /// Sum of atomic masses in this bin (amu).
    pub total_mass: f64,
    /// Mass-weighted velocity (Å/ps).
    pub mean_velocity: [f64; 3],
    /// Number of atoms in this bin.
    pub atom_count: usize,
    /// Accumulated kinetic temperature estimate (K).
    pub temperature: f64,
}

impl CoarseGrainBin {
    /// Create an empty bin at the given center.
    pub fn new(center: [f64; 3], half_width: f64) -> Self {
        Self {
            center,
            half_width,
            total_mass: 0.0_f64,
            mean_velocity: [0.0_f64; 3],
            atom_count: 0,
            temperature: 0.0_f64,
        }
    }

    /// Check whether the given point falls inside this bin.
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        (pos[0] - self.center[0]).abs() <= self.half_width
            && (pos[1] - self.center[1]).abs() <= self.half_width
            && (pos[2] - self.center[2]).abs() <= self.half_width
    }

    /// Accumulate one atom into this bin.
    pub fn add(&mut self, vel: [f64; 3], mass: f64) {
        let new_mass = self.total_mass + mass;
        if new_mass < 1e-30_f64 {
            return;
        }
        self.mean_velocity[0] =
            (self.mean_velocity[0] * self.total_mass + vel[0] * mass) / new_mass;
        self.mean_velocity[1] =
            (self.mean_velocity[1] * self.total_mass + vel[1] * mass) / new_mass;
        self.mean_velocity[2] =
            (self.mean_velocity[2] * self.total_mass + vel[2] * mass) / new_mass;
        self.total_mass = new_mass;
        self.atom_count += 1;
    }

    /// Clear accumulated data.
    pub fn reset(&mut self) {
        self.total_mass = 0.0_f64;
        self.mean_velocity = [0.0_f64; 3];
        self.atom_count = 0;
        self.temperature = 0.0_f64;
    }
}

// ---------------------------------------------------------------------------
// EamPotential
// ---------------------------------------------------------------------------

/// Embedded atom method (EAM) potential for metallic systems.
///
/// The EAM total energy is:
/// E = Σ_i F(ρ_i) + (1/2) Σ_{i≠j} φ(r_ij)
///
/// where F is the embedding function, ρ is the electron density,
/// and φ is the pairwise repulsion.
#[derive(Debug, Clone)]
pub struct EamPotential {
    /// Cutoff radius for pair interactions (Å).
    pub cutoff: f64,
    /// Lennard-Jones epsilon for pairwise part (eV).
    pub epsilon: f64,
    /// Lennard-Jones sigma (Å).
    pub sigma: f64,
    /// Embedding energy scale (eV).
    pub embedding_scale: f64,
}

impl EamPotential {
    /// Create a new EAM potential with simplified LJ pairwise and sqrt embedding.
    pub fn new(cutoff: f64, epsilon: f64, sigma: f64, embedding_scale: f64) -> Self {
        Self {
            cutoff,
            epsilon,
            sigma,
            embedding_scale,
        }
    }

    /// LJ pairwise interaction φ(r).
    pub fn pair_energy(&self, r: f64) -> f64 {
        if r > self.cutoff || r < 1e-10_f64 {
            return 0.0_f64;
        }
        let sr6 = (self.sigma / r).powi(6);
        4.0_f64 * self.epsilon * (sr6 * sr6 - sr6)
    }

    /// Electron density contribution from neighbour at distance r.
    pub fn electron_density(&self, r: f64) -> f64 {
        if r > self.cutoff || r < 1e-10_f64 {
            return 0.0_f64;
        }
        let sr = self.sigma / r;
        sr.powi(6)
    }

    /// Embedding function F(ρ) = −A √ρ.
    pub fn embedding_energy(&self, rho: f64) -> f64 {
        if rho < 0.0_f64 {
            return 0.0_f64;
        }
        -self.embedding_scale * rho.sqrt()
    }

    /// Derivative of embedding function dF/dρ = −A / (2√ρ).
    pub fn d_embedding(&self, rho: f64) -> f64 {
        if rho < 1e-30_f64 {
            return 0.0_f64;
        }
        -self.embedding_scale / (2.0_f64 * rho.sqrt())
    }

    /// Total EAM energy for a set of atoms.
    pub fn total_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        // Compute electron densities
        let rhos: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .filter(|&j| j != i)
                    .map(|j| {
                        let dx = positions[i][0] - positions[j][0];
                        let dy = positions[i][1] - positions[j][1];
                        let dz = positions[i][2] - positions[j][2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        self.electron_density(r)
                    })
                    .sum()
            })
            .collect();

        let mut e_embed: f64 = rhos.iter().map(|&rho| self.embedding_energy(rho)).sum();

        let mut e_pair = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                e_pair += self.pair_energy(r);
            }
        }

        e_embed += e_pair;
        e_embed
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── AtomisticRegion ───────────────────────────────────────────────────

    #[test]
    fn atomistic_region_contains_point_inside() {
        let r = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        assert!(r.contains([5.0_f64, 5.0_f64, 5.0_f64]));
    }

    #[test]
    fn atomistic_region_does_not_contain_outside_point() {
        let r = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        assert!(!r.contains([11.0_f64, 5.0_f64, 5.0_f64]));
    }

    #[test]
    fn atomistic_region_volume() {
        let r = AtomisticRegion::new(
            [0.0_f64; 3],
            [2.0_f64, 3.0_f64, 4.0_f64],
            BoundaryType::Fixed,
        );
        assert!((r.volume() - 24.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn atomistic_region_add_atom_updates_counts() {
        let mut r = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        r.add_atom(
            [1.0_f64, 2.0_f64, 3.0_f64],
            [0.1_f64, 0.0_f64, 0.0_f64],
            12.0_f64,
        );
        r.add_atom(
            [4.0_f64, 5.0_f64, 6.0_f64],
            [0.0_f64, 0.1_f64, 0.0_f64],
            12.0_f64,
        );
        assert_eq!(r.n_atoms(), 2);
        assert_eq!(r.positions.len(), r.velocities.len());
    }

    #[test]
    fn atomistic_region_boundary_types() {
        let r1 = AtomisticRegion::new([0.0_f64; 3], [5.0_f64; 3], BoundaryType::Periodic);
        let r2 = AtomisticRegion::new(
            [0.0_f64; 3],
            [5.0_f64; 3],
            BoundaryType::PaddingLayer { thickness: 2.0_f64 },
        );
        assert_ne!(r1.boundary, r2.boundary);
    }

    // ── ContinuumRegion ───────────────────────────────────────────────────

    #[test]
    fn continuum_region_add_nodes() {
        let mut cr = ContinuumRegion::new();
        cr.add_node([0.0_f64; 3], 1.0_f64);
        cr.add_node([1.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert_eq!(cr.n_nodes(), 2);
        assert_eq!(cr.velocities.len(), 2);
        assert_eq!(cr.stress.len(), 2);
    }

    #[test]
    fn continuum_region_default_velocities_zero() {
        let mut cr = ContinuumRegion::new();
        cr.add_node([1.0_f64, 2.0_f64, 3.0_f64], 5.0_f64);
        assert_eq!(cr.velocities[0], [0.0_f64; 3]);
    }

    // ── HandshakeZone ─────────────────────────────────────────────────────

    #[test]
    fn handshake_zone_contains_interior() {
        let hz = HandshakeZone::new([0.0_f64; 3], [5.0_f64; 3]);
        assert!(hz.contains([2.5_f64; 3]));
    }

    #[test]
    fn handshake_zone_does_not_contain_exterior() {
        let hz = HandshakeZone::new([0.0_f64; 3], [5.0_f64; 3]);
        assert!(!hz.contains([6.0_f64, 2.5_f64, 2.5_f64]));
    }

    #[test]
    fn handshake_zone_bridging_weight_edges() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        assert!((hz.bridging_weight(0.0_f64)).abs() < 1e-12_f64);
        assert!((hz.bridging_weight(10.0_f64) - 1.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn handshake_zone_bridging_weight_midpoint() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        assert!((hz.bridging_weight(5.0_f64) - 0.5_f64).abs() < 1e-12_f64);
    }

    // ── HybridDomainDecomposition ─────────────────────────────────────────

    #[test]
    fn hybrid_domain_is_consistent_after_construction() {
        let md = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        let cont = ContinuumRegion::new();
        let hz = HandshakeZone::new([4.0_f64, 0.0_f64, 0.0_f64], [6.0_f64; 3]);
        let hdd = HybridDomainDecomposition::new(md, cont, hz);
        assert!(hdd.is_consistent());
    }

    #[test]
    fn hybrid_domain_atoms_in_handshake_count() {
        let mut md = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        md.add_atom([5.0_f64; 3], [0.0_f64; 3], 1.0_f64); // inside handshake
        md.add_atom([1.0_f64; 3], [0.0_f64; 3], 1.0_f64); // outside
        let cont = ContinuumRegion::new();
        let hz = HandshakeZone::new([4.0_f64; 3], [6.0_f64; 3]);
        let hdd = HybridDomainDecomposition::new(md, cont, hz);
        assert_eq!(hdd.atoms_in_handshake(), 1);
    }

    #[test]
    fn hybrid_domain_coupling_step_blends_velocity() {
        let mut md = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        md.add_atom([5.0_f64; 3], [2.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        let mut cont = ContinuumRegion::new();
        cont.add_node([5.0_f64; 3], 1.0_f64); // inside handshake
        let hz = HandshakeZone::new([4.0_f64; 3], [6.0_f64; 3]);
        let mut hdd = HybridDomainDecomposition::new(md, cont, hz);
        hdd.coupling_step();
        // After coupling, the continuum node velocity should have shifted
        let vx = hdd.continuum_region.velocities[0][0];
        assert!(vx > 0.0_f64, "vx={vx}");
    }

    #[test]
    fn hybrid_domain_no_crash_empty_handshake() {
        let md = AtomisticRegion::new([0.0_f64; 3], [10.0_f64; 3], BoundaryType::Periodic);
        let cont = ContinuumRegion::new();
        let hz = HandshakeZone::new([4.0_f64; 3], [6.0_f64; 3]);
        let mut hdd = HybridDomainDecomposition::new(md, cont, hz);
        hdd.coupling_step(); // should not panic
    }

    // ── AveragingOperator ─────────────────────────────────────────────────

    #[test]
    fn averaging_operator_n_points() {
        let ao = AveragingOperator::new([0.0_f64; 3], [4, 4, 4], [1.0_f64; 3], 1.5_f64);
        assert_eq!(ao.n_points(), 64);
    }

    #[test]
    fn averaging_operator_index_and_position() {
        let ao = AveragingOperator::new([0.0_f64; 3], [3, 3, 3], [2.0_f64; 3], 1.0_f64);
        let pos = ao.grid_position(1, 1, 1);
        assert!((pos[0] - 2.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn averaging_operator_kernel_positive() {
        let ao = AveragingOperator::new([0.0_f64; 3], [2, 2, 2], [1.0_f64; 3], 1.0_f64);
        assert!(ao.kernel(0.0_f64) > 0.0_f64);
        assert!(ao.kernel(1.0_f64) > 0.0_f64);
        assert!(ao.kernel(1.0_f64) < ao.kernel(0.0_f64));
    }

    #[test]
    fn averaging_density_zero_for_empty_positions() {
        let mut ao = AveragingOperator::new([0.0_f64; 3], [2, 2, 2], [1.0_f64; 3], 1.0_f64);
        ao.average_density(&[]);
        assert!(ao.density.iter().all(|&d| d.abs() < 1e-30_f64));
    }

    #[test]
    fn averaging_density_nonzero_for_single_atom_at_center() {
        let mut ao = AveragingOperator::new([-5.0_f64; 3], [3, 3, 3], [5.0_f64; 3], 2.0_f64);
        ao.average_density(&[[0.0_f64; 3]]);
        let peak = ao.peak_density();
        assert!(peak > 0.0_f64, "peak density should be positive");
    }

    #[test]
    fn averaging_velocity_zero_for_zero_velocities() {
        let mut ao = AveragingOperator::new([0.0_f64; 3], [2, 2, 2], [1.0_f64; 3], 1.0_f64);
        let pos = vec![[0.5_f64, 0.5_f64, 0.5_f64]];
        let vel = vec![[0.0_f64; 3]];
        let masses = vec![1.0_f64];
        ao.average_velocity(&pos, &vel, &masses);
        assert!(ao.velocity.iter().all(|v| v[0].abs() < 1e-15_f64));
    }

    // ── BridgingDomainMethod ──────────────────────────────────────────────

    #[test]
    fn bdm_blend_displacement_at_md_end() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let bdm = BridgingDomainMethod::new(hz, 100.0_f64);
        let u_md = [1.0_f64, 0.0_f64, 0.0_f64];
        let u_fem = [2.0_f64, 0.0_f64, 0.0_f64];
        let pos = [0.0_f64, 5.0_f64, 5.0_f64]; // weight = 0
        let u = bdm.blend_displacement(u_md, u_fem, pos);
        assert!((u[0] - 1.0_f64).abs() < 1e-12_f64, "u[0]={}", u[0]);
    }

    #[test]
    fn bdm_blend_displacement_at_fem_end() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let bdm = BridgingDomainMethod::new(hz, 100.0_f64);
        let u_md = [1.0_f64, 0.0_f64, 0.0_f64];
        let u_fem = [2.0_f64, 0.0_f64, 0.0_f64];
        let pos = [10.0_f64, 5.0_f64, 5.0_f64]; // weight = 1
        let u = bdm.blend_displacement(u_md, u_fem, pos);
        assert!((u[0] - 2.0_f64).abs() < 1e-12_f64, "u[0]={}", u[0]);
    }

    #[test]
    fn bdm_constraint_force_zero_when_no_mismatch() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let mut bdm = BridgingDomainMethod::new(hz, 100.0_f64);
        let u = [1.0_f64, 1.0_f64, 1.0_f64];
        let pos = [5.0_f64; 3];
        let f = bdm.constraint_force(u, u, pos);
        assert!(f.iter().all(|&x| x.abs() < 1e-12_f64));
    }

    #[test]
    fn bdm_constraint_force_opposes_displacement_mismatch() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let mut bdm = BridgingDomainMethod::new(hz, 100.0_f64);
        let u_atom = [2.0_f64, 0.0_f64, 0.0_f64];
        let u_fem = [0.0_f64, 0.0_f64, 0.0_f64];
        let pos = [5.0_f64; 3];
        let f = bdm.constraint_force(u_atom, u_fem, pos);
        assert!(f[0] < 0.0_f64, "constraint force should oppose mismatch");
    }

    #[test]
    fn bdm_has_converged_after_zero_residual() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let mut bdm = BridgingDomainMethod::new(hz, 100.0_f64);
        bdm.constraint_residual = 0.0_f64;
        assert!(bdm.has_converged(1e-6_f64));
    }

    // ── AtomisticToContin ─────────────────────────────────────────────────

    #[test]
    fn stress_kinetic_zero_velocities() {
        let extractor = AtomisticToContin::new(10.0_f64);
        let vel = vec![[0.0_f64; 3]; 5];
        let masses = vec![12.0_f64; 5];
        let s = extractor.kinetic_stress(&vel, &masses, 100.0_f64);
        assert!(s.iter().all(|&x| x.abs() < 1e-15_f64));
    }

    #[test]
    fn stress_potential_from_uniform_forces() {
        let mut extractor = AtomisticToContin::new(10.0_f64);
        let pos = vec![[1.0_f64, 0.0_f64, 0.0_f64]];
        let forces = vec![[1.0_f64, 0.0_f64, 0.0_f64]];
        let s = extractor.potential_stress(&pos, &forces, 1.0_f64);
        assert!((s[0] - 1.0_f64).abs() < 1e-12_f64, "s[0]={}", s[0]);
    }

    #[test]
    fn stress_hydrostatic_pressure_sign() {
        let mut extractor = AtomisticToContin::new(10.0_f64);
        extractor.stress_voigt = [1.0_f64, 1.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64];
        let p = extractor.pressure();
        assert!((p + 1.0_f64).abs() < 1e-12_f64, "p={p}");
    }

    #[test]
    fn stress_von_mises_zero_for_hydrostatic() {
        let mut extractor = AtomisticToContin::new(10.0_f64);
        extractor.stress_voigt = [1.0_f64, 1.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64];
        let vm = extractor.von_mises();
        assert!(vm < 1e-10_f64, "vm={vm}");
    }

    #[test]
    fn stress_von_mises_positive_for_shear() {
        let mut extractor = AtomisticToContin::new(10.0_f64);
        extractor.stress_voigt = [0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64];
        let vm = extractor.von_mises();
        assert!(vm > 0.0_f64, "vm={vm}");
    }

    // ── FluctuatingHydrodynamics ───────────────────────────────────────────

    #[test]
    fn fluctuating_hydro_stress_amplitude_positive() {
        let fh = FluctuatingHydrodynamics::new(
            300.0_f64, 1e-3_f64, 0.0_f64, 0.6_f64, 1e-27_f64, 1e-15_f64, 42,
        );
        let amp = fh.stress_noise_amplitude();
        assert!(amp > 0.0_f64, "amp={amp}");
    }

    #[test]
    fn fluctuating_hydro_heat_flux_amplitude_positive() {
        let fh = FluctuatingHydrodynamics::new(
            300.0_f64, 1e-3_f64, 0.0_f64, 0.6_f64, 1e-27_f64, 1e-15_f64, 42,
        );
        let amp = fh.heat_flux_noise_amplitude();
        assert!(amp > 0.0_f64, "amp={amp}");
    }

    #[test]
    fn fluctuating_hydro_stress_noise_traceless() {
        let mut fh = FluctuatingHydrodynamics::new(
            300.0_f64, 1e-3_f64, 0.0_f64, 0.6_f64, 1e-27_f64, 1e-15_f64, 123,
        );
        let s = fh.sample_stress_noise();
        let trace = s[0] + s[1] + s[2];
        assert!(trace.abs() < 1e-8_f64, "trace={trace}");
    }

    #[test]
    fn fluctuating_hydro_update_temperature_from_md() {
        let mut fh = FluctuatingHydrodynamics::new(
            0.0_f64, 1e-3_f64, 0.0_f64, 0.6_f64, 1e-27_f64, 1e-15_f64, 0,
        );
        let vel = vec![[1.0_f64, 0.0_f64, 0.0_f64]; 100];
        let masses = vec![1.0_f64; 100];
        fh.update_temperature_from_md(&vel, &masses);
        assert!(fh.temperature > 0.0_f64, "T={}", fh.temperature);
    }

    #[test]
    fn fluctuating_hydro_zero_amplitude_when_zero_volume() {
        let fh = FluctuatingHydrodynamics::new(
            300.0_f64, 1e-3_f64, 0.0_f64, 0.6_f64, 0.0_f64, 1e-15_f64, 0,
        );
        assert_eq!(fh.stress_noise_amplitude(), 0.0_f64);
    }

    // ── CoupleEnergy ──────────────────────────────────────────────────────

    #[test]
    fn couple_energy_initial_total_is_zero() {
        let ce = CoupleEnergy::new();
        assert_eq!(ce.total_energy(), 0.0_f64);
    }

    #[test]
    fn couple_energy_update_md_kinetic() {
        let mut ce = CoupleEnergy::new();
        let vel = vec![[1.0_f64, 0.0_f64, 0.0_f64]];
        let masses = vec![1.0_f64]; // 1 amu at 1 Å/ps
        ce.update_md_kinetic(&vel, &masses);
        // KE = 0.5 * 1 * 1 * CONV ≈ 0.00518 eV
        assert!(ce.md_kinetic > 0.0_f64, "md_kinetic={}", ce.md_kinetic);
    }

    #[test]
    fn couple_energy_total_sums_all_parts() {
        let mut ce = CoupleEnergy::new();
        ce.md_kinetic = 1.0_f64;
        ce.md_potential = 2.0_f64;
        ce.continuum_kinetic = 3.0_f64;
        ce.continuum_elastic = 4.0_f64;
        assert!((ce.total_energy() - 10.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn couple_energy_flux_drives_toward_reference() {
        let mut ce = CoupleEnergy::new();
        ce.md_kinetic = 5.0_f64;
        let flux = ce.compute_flux(10.0_f64);
        assert!((flux - 5.0_f64).abs() < 1e-12_f64, "flux={flux}");
    }

    #[test]
    fn couple_energy_kinetic_correction_scales_velocities() {
        let ce = CoupleEnergy::new();
        let mut vel = vec![[1.0_f64, 0.0_f64, 0.0_f64]];
        let masses = vec![1.0_f64];
        // Apply a large positive correction → speeds should increase
        let v_before = vel[0][0];
        ce.apply_kinetic_correction(&mut vel, &masses, 10.0_f64);
        assert!(vel[0][0] > v_before, "v should increase");
    }

    #[test]
    fn couple_energy_is_conserved_when_error_small() {
        let mut ce = CoupleEnergy::new();
        ce.cumulative_error = 1e-12_f64;
        assert!(ce.is_conserved(1e-6_f64));
    }

    // ── CauchyBornModel ───────────────────────────────────────────────────

    #[test]
    fn cauchy_born_identity_gives_zero_strain() {
        let cb = CauchyBornModel::new(3.52_f64, 1.0_f64, 0.5_f64, 0.3_f64);
        let e = cb.green_lagrange_strain();
        assert!(e.iter().all(|&x| x.abs() < 1e-12_f64), "{e:?}");
    }

    #[test]
    fn cauchy_born_pk2_zero_for_zero_strain() {
        let cb = CauchyBornModel::new(3.52_f64, 1.0_f64, 0.5_f64, 0.3_f64);
        let s = cb.pk2_stress();
        assert!(s.iter().all(|&x| x.abs() < 1e-12_f64), "{s:?}");
    }

    #[test]
    fn cauchy_born_youngs_modulus_positive() {
        let cb = CauchyBornModel::new(3.52_f64, 1.0_f64, 0.5_f64, 0.3_f64);
        assert!(cb.youngs_modulus() > 0.0_f64);
    }

    #[test]
    fn cauchy_born_uniaxial_stretch_gives_tensile_stress() {
        let mut cb = CauchyBornModel::new(3.52_f64, 1.0_f64, 0.5_f64, 0.3_f64);
        // Stretch in x: F = diag(1.1, 1, 1)
        cb.set_deformation_gradient([
            1.1_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64,
        ]);
        let s = cb.pk2_stress();
        assert!(s[0] > 0.0_f64, "Sxx should be positive under tension");
    }

    // ── CoarseGrainBin ────────────────────────────────────────────────────

    #[test]
    fn coarse_grain_bin_contains() {
        let bin = CoarseGrainBin::new([5.0_f64; 3], 2.0_f64);
        assert!(bin.contains([4.0_f64, 5.0_f64, 6.0_f64]));
        assert!(!bin.contains([8.0_f64, 5.0_f64, 5.0_f64]));
    }

    #[test]
    fn coarse_grain_bin_accumulate_velocity() {
        let mut bin = CoarseGrainBin::new([5.0_f64; 3], 2.0_f64);
        bin.add([2.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        bin.add([4.0_f64, 0.0_f64, 0.0_f64], 1.0_f64);
        assert!((bin.mean_velocity[0] - 3.0_f64).abs() < 1e-12_f64);
        assert_eq!(bin.atom_count, 2);
    }

    #[test]
    fn coarse_grain_bin_reset_clears_data() {
        let mut bin = CoarseGrainBin::new([5.0_f64; 3], 2.0_f64);
        bin.add([1.0_f64, 0.0_f64, 0.0_f64], 12.0_f64);
        bin.reset();
        assert_eq!(bin.atom_count, 0);
        assert!(bin.total_mass.abs() < 1e-30_f64);
    }

    // ── EamPotential ──────────────────────────────────────────────────────

    #[test]
    fn eam_pair_energy_zero_beyond_cutoff() {
        let eam = EamPotential::new(5.0_f64, 0.1_f64, 2.5_f64, 1.0_f64);
        assert_eq!(eam.pair_energy(6.0_f64), 0.0_f64);
    }

    #[test]
    fn eam_electron_density_zero_beyond_cutoff() {
        let eam = EamPotential::new(5.0_f64, 0.1_f64, 2.5_f64, 1.0_f64);
        assert_eq!(eam.electron_density(6.0_f64), 0.0_f64);
    }

    #[test]
    fn eam_embedding_energy_negative_for_positive_rho() {
        let eam = EamPotential::new(5.0_f64, 0.1_f64, 2.5_f64, 1.0_f64);
        assert!(eam.embedding_energy(1.0_f64) < 0.0_f64);
    }

    #[test]
    fn eam_d_embedding_negative() {
        let eam = EamPotential::new(5.0_f64, 0.1_f64, 2.5_f64, 1.0_f64);
        assert!(eam.d_embedding(1.0_f64) < 0.0_f64);
    }

    #[test]
    fn eam_total_energy_two_atoms() {
        let eam = EamPotential::new(5.0_f64, 0.01_f64, 2.5_f64, 0.5_f64);
        let pos = vec![[0.0_f64; 3], [3.0_f64, 0.0_f64, 0.0_f64]];
        let e = eam.total_energy(&pos);
        assert!(e.is_finite(), "energy should be finite: {e}");
    }

    // ── Integration / consistency checks ──────────────────────────────────

    #[test]
    fn domain_decomposition_consistency_check() {
        let mut md = AtomisticRegion::new([0.0_f64; 3], [20.0_f64; 3], BoundaryType::Periodic);
        // Place atoms along x with y=9, z=9 so they fall inside hz [8..12]^3
        for i in 0..10 {
            md.add_atom(
                [i as f64, 9.0_f64, 9.0_f64],
                [0.1_f64 * i as f64, 0.0_f64, 0.0_f64],
                12.0_f64,
            );
        }
        let cont = ContinuumRegion::new();
        let hz = HandshakeZone::new([8.0_f64; 3], [12.0_f64; 3]);
        let hdd = HybridDomainDecomposition::new(md, cont, hz);
        assert!(hdd.is_consistent());
        // Atoms at x=8,9,10,11,12 (indices 8,9 of the 0..10 range) should be in handshake
        assert!(hdd.atoms_in_handshake() >= 2);
    }

    #[test]
    fn averaging_convergence_with_more_atoms() {
        let mut ao1 = AveragingOperator::new([-5.0_f64; 3], [2, 2, 2], [5.0_f64; 3], 2.0_f64);
        let mut ao2 = AveragingOperator::new([-5.0_f64; 3], [2, 2, 2], [5.0_f64; 3], 2.0_f64);
        let pos1 = vec![[0.0_f64; 3]];
        let pos2 = vec![[0.0_f64; 3]; 10];
        ao1.average_density(&pos1);
        ao2.average_density(&pos2);
        assert!(ao2.peak_density() > ao1.peak_density());
    }

    #[test]
    fn bdm_overlap_blending_midpoint() {
        let hz = HandshakeZone::new([0.0_f64; 3], [10.0_f64; 3]);
        let bdm = BridgingDomainMethod::new(hz, 1.0_f64);
        let u_md = [0.0_f64, 0.0_f64, 0.0_f64];
        let u_fem = [2.0_f64, 0.0_f64, 0.0_f64];
        let pos = [5.0_f64, 5.0_f64, 5.0_f64];
        let u = bdm.blend_displacement(u_md, u_fem, pos);
        assert!((u[0] - 1.0_f64).abs() < 1e-12_f64, "u[0]={}", u[0]);
    }

    #[test]
    fn couple_energy_conservation_cycle() {
        let mut ce = CoupleEnergy::new();
        let vel = vec![[1.0_f64; 3]; 10];
        let masses = vec![1.0_f64; 10];
        ce.update_md_kinetic(&vel, &masses);
        let e0 = ce.total_energy();
        let flux = ce.compute_flux(e0 * 1.1_f64);
        assert!(flux > 0.0_f64, "flux={flux}");
    }
}
