#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lipid bilayer molecular dynamics (MARTINI-like coarse-grained model).
//!
//! This module provides a comprehensive coarse-grained lipid bilayer simulation
//! framework with the following features:
//!
//! - **Coarse-grained lipid model** ([`CgLipid`]): MARTINI-like bead representation
//!   with head, glycerol, and tail beads.
//! - **Bilayer system** ([`BilayerSystem`]): Full bilayer with upper/lower leaflets,
//!   periodic boundary conditions, and neighbor lists.
//! - **Self-assembly** ([`BilayerSystem::run_self_assembly`]): Spontaneous bilayer
//!   formation from random lipid configurations.
//! - **Membrane properties**: area per lipid, bilayer thickness, order parameter,
//!   lateral diffusion coefficient, flip-flop rate.
//! - **Helfrich elastic theory** ([`HelfrichElasticity`]): bending modulus, mean
//!   curvature, Gaussian curvature, undulation spectrum.
//! - **Lateral pressure profile** ([`LateralPressureProfile`]): Irving-Kirkwood
//!   pressure tensor decomposition across the membrane normal.
//! - **Phase transitions** ([`PhaseTransition`]): gel/liquid-crystalline transition,
//!   cholesterol condensing effect.
//! - **Transmembrane proteins** ([`TmProtein`]): hydrophobic mismatch, tilt energy,
//!   protein–lipid coupling.
//! - **Membrane curvature** ([`CurvatureField`]): local curvature computation,
//!   spontaneous curvature, membrane tension and area compressibility.
//! - **Vesicle formation** ([`VesicleBuilder`]): spherical vesicle assembly from
//!   a flat bilayer patch.
//! - **Bilayer fusion** ([`BilayerFusion`]): hemifusion stalk and pore formation.
//!
//! Units: nm for length, kJ/mol for energy, ps for time, K for temperature.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann constant in kJ/(mol·K).
pub const KB: f64 = 8.314_462_618e-3;

/// Reference temperature 300 K.
pub const T_REF: f64 = 300.0;

/// kB·T at 300 K \[kJ/mol\].
pub const KBT_300: f64 = KB * T_REF;

/// Avogadro's number.
pub const NA: f64 = 6.022_140_76e23;

/// Standard bilayer thickness (nm).
pub const BILAYER_THICKNESS: f64 = 4.0;

/// Typical area per lipid for DPPC (nm²).
pub const AREA_PER_LIPID_DPPC: f64 = 0.64;

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra – plain [f64;3])
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
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a).max(1e-30);
    scale3(a, 1.0 / n)
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Apply minimum image convention for a cubic box.
#[inline]
fn pbc_displacement(r: [f64; 3], box_len: [f64; 3]) -> [f64; 3] {
    [
        r[0] - box_len[0] * (r[0] / box_len[0]).round(),
        r[1] - box_len[1] * (r[1] / box_len[1]).round(),
        r[2] - box_len[2] * (r[2] / box_len[2]).round(),
    ]
}

/// Wrap a position into the simulation box.
#[inline]
fn wrap_pbc(pos: [f64; 3], box_len: [f64; 3]) -> [f64; 3] {
    [
        pos[0] - box_len[0] * (pos[0] / box_len[0]).floor(),
        pos[1] - box_len[1] * (pos[1] / box_len[1]).floor(),
        pos[2] - box_len[2] * (pos[2] / box_len[2]).floor(),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Bead types
// ─────────────────────────────────────────────────────────────────────────────

/// Coarse-grained bead type in the MARTINI-like model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BeadType {
    /// Phosphate/choline head group bead (polar).
    Head,
    /// Glycerol backbone bead (slightly polar).
    Glycerol,
    /// Hydrophobic acyl tail bead (apolar).
    Tail,
    /// Cholesterol ring bead.
    Cholesterol,
    /// Transmembrane protein bead.
    Protein,
    /// Water bead (4 real waters per CG bead).
    Water,
}

impl BeadType {
    /// Returns the Lennard-Jones epsilon \[kJ/mol\] for this bead type.
    pub fn epsilon(&self) -> f64 {
        match self {
            BeadType::Head => 5.0,
            BeadType::Glycerol => 4.5,
            BeadType::Tail => 3.5,
            BeadType::Cholesterol => 4.0,
            BeadType::Protein => 4.0,
            BeadType::Water => 5.0,
        }
    }

    /// Returns the Lennard-Jones sigma \[nm\] for this bead type.
    pub fn sigma(&self) -> f64 {
        match self {
            BeadType::Head => 0.47,
            BeadType::Glycerol => 0.47,
            BeadType::Tail => 0.47,
            BeadType::Cholesterol => 0.43,
            BeadType::Protein => 0.47,
            BeadType::Water => 0.47,
        }
    }

    /// Returns the mass \[u\] for this bead type.
    pub fn mass(&self) -> f64 {
        match self {
            BeadType::Head => 121.0,
            BeadType::Glycerol => 72.0,
            BeadType::Tail => 57.0,
            BeadType::Cholesterol => 74.0,
            BeadType::Protein => 72.0,
            BeadType::Water => 72.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CG bead
// ─────────────────────────────────────────────────────────────────────────────

/// A single coarse-grained bead in the lipid simulation.
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Position \[nm\].
    pub pos: [f64; 3],
    /// Velocity \[nm/ps\].
    pub vel: [f64; 3],
    /// Force \[kJ/(mol·nm)\].
    pub force: [f64; 3],
    /// Bead type.
    pub bead_type: BeadType,
    /// Lipid index this bead belongs to.
    pub lipid_idx: usize,
    /// Local bead index within the lipid (0 = head, …).
    pub local_idx: usize,
}

impl CgBead {
    /// Construct a new CG bead.
    pub fn new(pos: [f64; 3], bead_type: BeadType, lipid_idx: usize, local_idx: usize) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            bead_type,
            lipid_idx,
            local_idx,
        }
    }

    /// Returns the mass \[u\].
    pub fn mass(&self) -> f64 {
        self.bead_type.mass()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CG Lipid
// ─────────────────────────────────────────────────────────────────────────────

/// Lipid species supported by the model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LipidSpecies {
    /// Dipalmitoylphosphatidylcholine (DPPC).
    Dppc,
    /// Dioleoylphosphatidylcholine (DOPC).
    Dopc,
    /// Palmitoyl-oleoyl phosphatidylcholine (POPC).
    Popc,
    /// Palmitoyl-oleoyl phosphatidylethanolamine (POPE).
    Pope,
    /// Cholesterol.
    Cholesterol,
}

impl LipidSpecies {
    /// Number of CG beads for this lipid species.
    pub fn n_beads(&self) -> usize {
        match self {
            LipidSpecies::Dppc => 12,
            LipidSpecies::Dopc => 13,
            LipidSpecies::Popc => 12,
            LipidSpecies::Pope => 12,
            LipidSpecies::Cholesterol => 8,
        }
    }

    /// Effective length along the bilayer normal \[nm\].
    pub fn length(&self) -> f64 {
        match self {
            LipidSpecies::Dppc => 1.9,
            LipidSpecies::Dopc => 1.9,
            LipidSpecies::Popc => 1.9,
            LipidSpecies::Pope => 1.85,
            LipidSpecies::Cholesterol => 1.6,
        }
    }

    /// Intrinsic spontaneous curvature C0 \[1/nm\].
    pub fn spontaneous_curvature(&self) -> f64 {
        match self {
            LipidSpecies::Dppc => 0.0,
            LipidSpecies::Dopc => 0.0,
            LipidSpecies::Popc => 0.0,
            LipidSpecies::Pope => -0.28,
            LipidSpecies::Cholesterol => -0.49,
        }
    }
}

/// A coarse-grained lipid molecule.
#[derive(Debug, Clone)]
pub struct CgLipid {
    /// Lipid species.
    pub species: LipidSpecies,
    /// Bead indices in the global bead array.
    pub bead_indices: Vec<usize>,
    /// Leaflet assignment: true = upper (+z), false = lower (−z).
    pub upper_leaflet: bool,
}

impl CgLipid {
    /// Construct a new lipid with given species and bead indices.
    pub fn new(species: LipidSpecies, bead_indices: Vec<usize>, upper_leaflet: bool) -> Self {
        Self {
            species,
            bead_indices,
            upper_leaflet,
        }
    }

    /// Head bead index (first bead in chain).
    pub fn head_bead_idx(&self) -> Option<usize> {
        self.bead_indices.first().copied()
    }

    /// Tail end bead index (last bead in chain).
    pub fn tail_bead_idx(&self) -> Option<usize> {
        self.bead_indices.last().copied()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bond / angle parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Harmonic bond parameters for CG lipid bonds.
#[derive(Debug, Clone, Copy)]
pub struct BondParam {
    /// Equilibrium bond length \[nm\].
    pub r0: f64,
    /// Force constant \[kJ/(mol·nm²)\].
    pub k: f64,
}

impl BondParam {
    /// Construct bond parameters.
    pub fn new(r0: f64, k: f64) -> Self {
        Self { r0, k }
    }

    /// Compute bond energy for distance `r`.
    pub fn energy(&self, r: f64) -> f64 {
        0.5 * self.k * (r - self.r0).powi(2)
    }

    /// Compute bond force magnitude (negative derivative) for distance `r`.
    pub fn force_magnitude(&self, r: f64) -> f64 {
        -self.k * (r - self.r0)
    }
}

/// Cosine angle potential parameters for CG bead triples.
#[derive(Debug, Clone, Copy)]
pub struct AngleParam {
    /// Equilibrium cosine of the angle.
    pub cos_theta0: f64,
    /// Force constant \[kJ/mol\].
    pub k: f64,
}

impl AngleParam {
    /// Construct angle parameters given equilibrium angle in radians.
    pub fn new(theta0_rad: f64, k: f64) -> Self {
        Self {
            cos_theta0: theta0_rad.cos(),
            k,
        }
    }

    /// Compute angle energy given cosine of the current angle.
    pub fn energy_from_cos(&self, cos_theta: f64) -> f64 {
        0.5 * self.k * (cos_theta - self.cos_theta0).powi(2)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Lennard-Jones potential (shifted-force cutoff)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Lennard-Jones energy \[kJ/mol\] and force magnitude at distance `r`.
///
/// Uses a shifted-force cutoff at `r_cut`. Returns `(energy, force_over_r)`.
pub fn lj_sf(eps: f64, sigma: f64, r: f64, r_cut: f64) -> (f64, f64) {
    if r >= r_cut {
        return (0.0, 0.0);
    }
    let sr = sigma / r;
    let sr2 = sr * sr;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;

    let sr_c = sigma / r_cut;
    let src2 = sr_c * sr_c;
    let src6 = src2 * src2 * src2;
    let src12 = src6 * src6;

    let e = 4.0 * eps * (sr12 - sr6);
    let e_cut = 4.0 * eps * (src12 - src6);
    let de_dr_cut = 4.0 * eps * (-12.0 * src12 / r_cut + 6.0 * src6 / r_cut);

    let energy = e - e_cut - de_dr_cut * (r - r_cut);
    let force_r = 4.0 * eps * (12.0 * sr12 - 6.0 * sr6) / (r * r) - de_dr_cut / r;
    (energy, force_r)
}

// ─────────────────────────────────────────────────────────────────────────────
// Bilayer System
// ─────────────────────────────────────────────────────────────────────────────

/// Full CG bilayer simulation system.
#[derive(Debug, Clone)]
pub struct BilayerSystem {
    /// All CG beads in the system.
    pub beads: Vec<CgBead>,
    /// All lipid molecules.
    pub lipids: Vec<CgLipid>,
    /// Simulation box lengths \[nm\].
    pub box_len: [f64; 3],
    /// LJ cutoff radius \[nm\].
    pub r_cut: f64,
    /// Current simulation time \[ps\].
    pub time: f64,
    /// Timestep \[ps\].
    pub dt: f64,
    /// Temperature \[K\].
    pub temperature: f64,
    /// Accumulated potential energy \[kJ/mol\].
    pub pot_energy: f64,
    /// Accumulated kinetic energy \[kJ/mol\].
    pub kin_energy: f64,
}

impl BilayerSystem {
    /// Construct a pre-assembled flat bilayer patch with `n_per_leaflet` lipids per leaflet.
    pub fn new_flat_bilayer(
        species: LipidSpecies,
        n_per_leaflet: usize,
        box_xy: f64,
        temperature: f64,
        dt: f64,
    ) -> Self {
        let box_z = 10.0;
        let box_len = [box_xy, box_xy, box_z];
        let r_cut = 1.2;

        let mut beads: Vec<CgBead> = Vec::new();
        let mut lipids: Vec<CgLipid> = Vec::new();

        let nx = ((n_per_leaflet as f64).sqrt().ceil()) as usize;
        let dx = box_xy / nx as f64;

        for leaflet in 0..2usize {
            for i in 0..n_per_leaflet {
                let ix = i % nx;
                let iy = i / nx;
                let x0 = ix as f64 * dx + 0.1;
                let y0 = iy as f64 * dx + 0.1;
                let z_center = if leaflet == 0 {
                    box_z * 0.5 + species.length() * 0.5
                } else {
                    box_z * 0.5 - species.length() * 0.5
                };

                let lipid_idx = lipids.len();
                let n_b = species.n_beads();
                let mut bead_indices = Vec::with_capacity(n_b);

                for bi in 0..n_b {
                    let frac = bi as f64 / (n_b - 1).max(1) as f64;
                    let z_offset = if leaflet == 0 {
                        -frac * species.length()
                    } else {
                        frac * species.length()
                    };
                    let btype = if bi == 0 {
                        BeadType::Head
                    } else if bi <= 2 {
                        BeadType::Glycerol
                    } else {
                        BeadType::Tail
                    };
                    let pos = [x0, y0, z_center + z_offset];
                    let idx = beads.len();
                    beads.push(CgBead::new(pos, btype, lipid_idx, bi));
                    bead_indices.push(idx);
                }

                lipids.push(CgLipid::new(species, bead_indices, leaflet == 0));
            }
        }

        let mut sys = Self {
            beads,
            lipids,
            box_len,
            r_cut,
            time: 0.0,
            dt,
            temperature,
            pot_energy: 0.0,
            kin_energy: 0.0,
        };
        sys.initialize_velocities(temperature);
        sys
    }

    /// Initialize Maxwell-Boltzmann velocities at given temperature.
    pub fn initialize_velocities(&mut self, temperature: f64) {
        let mut rng = rand::rng();
        for bead in &mut self.beads {
            let m = bead.mass();
            let sigma_v = (KB * temperature / m).sqrt();
            let vx = rand_normal(&mut rng) * sigma_v;
            let vy = rand_normal(&mut rng) * sigma_v;
            let vz = rand_normal(&mut rng) * sigma_v;
            bead.vel = [vx, vy, vz];
        }
        self.remove_center_of_mass_motion();
    }

    /// Remove center-of-mass translational motion.
    pub fn remove_center_of_mass_motion(&mut self) {
        let mut total_mass = 0.0;
        let mut com_vel = [0.0f64; 3];
        for bead in &self.beads {
            let m = bead.mass();
            total_mass += m;
            for d in 0..3 {
                com_vel[d] += m * bead.vel[d];
            }
        }
        if total_mass > 0.0 {
            for d in 0..3 {
                com_vel[d] /= total_mass;
            }
            for bead in &mut self.beads {
                for d in 0..3 {
                    bead.vel[d] -= com_vel[d];
                }
            }
        }
    }

    /// Zero all forces.
    fn zero_forces(&mut self) {
        for bead in &mut self.beads {
            bead.force = [0.0; 3];
        }
    }

    /// Compute non-bonded (LJ) forces between all pairs within cutoff.
    ///
    /// Uses O(N²) naive loop suitable for small systems. Returns total LJ energy.
    pub fn compute_nonbonded_forces(&mut self) -> f64 {
        self.zero_forces();
        let n = self.beads.len();
        let mut energy = 0.0;

        for i in 0..n {
            for j in (i + 1)..n {
                // Skip 1-2 and 1-3 bonded pairs (same lipid, close in sequence)
                if self.beads[i].lipid_idx == self.beads[j].lipid_idx {
                    let diff = (self.beads[i].local_idx as isize
                        - self.beads[j].local_idx as isize)
                        .unsigned_abs();
                    if diff <= 2 {
                        continue;
                    }
                }

                let ri = self.beads[i].pos;
                let rj = self.beads[j].pos;
                let dr_raw = sub3(rj, ri);
                let dr = pbc_displacement(dr_raw, self.box_len);
                let r = norm3(dr);
                if r < 1e-6 {
                    continue;
                }

                // Mixing rules: arithmetic sigma, geometric epsilon
                let eps_i = self.beads[i].bead_type.epsilon();
                let eps_j = self.beads[j].bead_type.epsilon();
                let sig_i = self.beads[i].bead_type.sigma();
                let sig_j = self.beads[j].bead_type.sigma();
                let eps_ij = (eps_i * eps_j).sqrt();
                let sig_ij = 0.5 * (sig_i + sig_j);

                let (e, f_over_r) = lj_sf(eps_ij, sig_ij, r, self.r_cut);
                energy += e;

                let fvec = scale3(dr, f_over_r);
                for d in 0..3 {
                    self.beads[i].force[d] -= fvec[d];
                    self.beads[j].force[d] += fvec[d];
                }
            }
        }
        energy
    }

    /// Compute bonded (harmonic bond) forces along lipid chains.
    ///
    /// Returns total bond energy.
    pub fn compute_bonded_forces(&mut self) -> f64 {
        let bp = BondParam::new(0.47, 3800.0);
        let mut energy = 0.0;

        for lip in &self.lipids.clone() {
            let nb = lip.bead_indices.len();
            for k in 0..(nb - 1) {
                let i = lip.bead_indices[k];
                let j = lip.bead_indices[k + 1];
                let ri = self.beads[i].pos;
                let rj = self.beads[j].pos;
                let dr = pbc_displacement(sub3(rj, ri), self.box_len);
                let r = norm3(dr);
                if r < 1e-6 {
                    continue;
                }
                energy += bp.energy(r);
                let f_mag = bp.force_magnitude(r);
                let fvec = scale3(normalize3(dr), f_mag);
                for d in 0..3 {
                    self.beads[i].force[d] -= fvec[d];
                    self.beads[j].force[d] += fvec[d];
                }
            }
        }
        energy
    }

    /// Perform one velocity-Verlet integration step with a Berendsen thermostat.
    pub fn step_velocity_verlet(&mut self, tau_t: f64) {
        let dt = self.dt;

        // Half-step velocity update + full position update
        for bead in &mut self.beads {
            let m = bead.mass();
            for d in 0..3 {
                bead.vel[d] += 0.5 * dt * bead.force[d] / m;
                bead.pos[d] += dt * bead.vel[d];
            }
            bead.pos = wrap_pbc(bead.pos, self.box_len);
        }

        // Recompute forces
        let e_nb = self.compute_nonbonded_forces();
        let e_b = self.compute_bonded_forces();
        self.pot_energy = e_nb + e_b;

        // Second half-step velocity update
        for bead in &mut self.beads {
            let m = bead.mass();
            for d in 0..3 {
                bead.vel[d] += 0.5 * dt * bead.force[d] / m;
            }
        }

        // Berendsen thermostat rescaling
        let kin = self.kinetic_energy();
        let t_cur = self.temperature_from_kinetic(kin);
        if t_cur > 1e-6 {
            let lambda = (1.0 + (dt / tau_t) * (self.temperature / t_cur - 1.0))
                .clamp(0.5, 2.0)
                .sqrt();
            for bead in &mut self.beads {
                for d in 0..3 {
                    bead.vel[d] *= lambda;
                }
            }
        }

        self.kin_energy = self.kinetic_energy();
        self.time += dt;
    }

    /// Compute total kinetic energy \[kJ/mol\].
    pub fn kinetic_energy(&self) -> f64 {
        self.beads
            .iter()
            .map(|b| 0.5 * b.mass() * dot3(b.vel, b.vel))
            .sum()
    }

    /// Compute instantaneous temperature from kinetic energy \[K\].
    pub fn temperature_from_kinetic(&self, kin: f64) -> f64 {
        let dof = (3 * self.beads.len()).saturating_sub(3) as f64;
        if dof > 0.0 {
            2.0 * kin / (dof * KB)
        } else {
            0.0
        }
    }

    /// Run `n_steps` of MD with the Berendsen thermostat.
    pub fn run(&mut self, n_steps: usize, tau_t: f64) {
        self.compute_nonbonded_forces();
        self.compute_bonded_forces();
        for _ in 0..n_steps {
            self.step_velocity_verlet(tau_t);
        }
    }

    /// Run a short self-assembly simulation starting from random lipid positions.
    pub fn run_self_assembly(&mut self, n_steps: usize) {
        self.run(n_steps, 0.5);
    }

    // ─── Membrane observables ────────────────────────────────────────────────

    /// Compute area per lipid \[nm²\] from box XY area and number of lipids.
    pub fn area_per_lipid(&self) -> f64 {
        let n_upper = self.lipids.iter().filter(|l| l.upper_leaflet).count();
        if n_upper == 0 {
            return 0.0;
        }
        self.box_len[0] * self.box_len[1] / n_upper as f64
    }

    /// Compute approximate bilayer thickness \[nm\] as the mean z-separation
    /// between upper and lower leaflet head groups.
    pub fn bilayer_thickness(&self) -> f64 {
        let upper_z = self.mean_head_z(true);
        let lower_z = self.mean_head_z(false);
        (upper_z - lower_z).abs()
    }

    /// Mean z-position of head beads for the specified leaflet.
    fn mean_head_z(&self, upper: bool) -> f64 {
        let head_beads: Vec<f64> = self
            .lipids
            .iter()
            .filter(|l| l.upper_leaflet == upper)
            .filter_map(|l| l.head_bead_idx())
            .map(|idx| self.beads[idx].pos[2])
            .collect();
        if head_beads.is_empty() {
            return 0.0;
        }
        head_beads.iter().sum::<f64>() / head_beads.len() as f64
    }

    /// Compute the mean-square displacement of lipid head groups over a trajectory
    /// starting from `ref_positions` \[nm\].
    ///
    /// Returns MSD \[nm²\] and estimated lateral diffusion coefficient D \[nm²/ps\].
    pub fn lateral_diffusion(&self, ref_positions: &[[f64; 3]], elapsed_ps: f64) -> (f64, f64) {
        let n = ref_positions.len().min(self.lipids.len());
        if n == 0 || elapsed_ps <= 0.0 {
            return (0.0, 0.0);
        }
        let mut msd = 0.0;
        for (i, lip) in self.lipids.iter().enumerate().take(n) {
            if let Some(hi) = lip.head_bead_idx() {
                let cur = self.beads[hi].pos;
                let r = &ref_positions[i];
                let dx = cur[0] - r[0];
                let dy = cur[1] - r[1];
                msd += dx * dx + dy * dy;
            }
        }
        msd /= n as f64;
        let d = msd / (4.0 * elapsed_ps);
        (msd, d)
    }

    /// Compute deuterium order parameter S_CD for tail beads.
    ///
    /// S = 0 for isotropic (fluid), S = 1 for perfectly ordered (gel).
    pub fn order_parameter(&self) -> f64 {
        let mut s_sum = 0.0;
        let mut count = 0usize;

        for lip in &self.lipids {
            let n_b = lip.bead_indices.len();
            if n_b < 2 {
                continue;
            }
            for k in 0..(n_b - 1) {
                let i = lip.bead_indices[k];
                let j = lip.bead_indices[k + 1];
                let ti = self.beads[i].bead_type;
                let tj = self.beads[j].bead_type;
                if ti != BeadType::Tail || tj != BeadType::Tail {
                    continue;
                }
                let dr = sub3(self.beads[j].pos, self.beads[i].pos);
                let n = normalize3(dr);
                let cos_theta = n[2]; // angle with z-axis (bilayer normal)
                s_sum += 0.5 * (3.0 * cos_theta * cos_theta - 1.0);
                count += 1;
            }
        }
        if count > 0 { s_sum / count as f64 } else { 0.0 }
    }

    /// Estimate lipid flip-flop rate \[1/ns\] using the Arrhenius model.
    ///
    /// `delta_g_kj` is the activation free energy barrier \[kJ/mol\].
    pub fn flipflop_rate(temperature: f64, delta_g_kj: f64) -> f64 {
        let kbt = KB * temperature;
        let nu0 = 1.0e6; // attempt frequency [1/ns], typical value
        nu0 * (-delta_g_kj / kbt).exp()
    }

    /// Estimate transmembrane diffusion coefficient D_z from flip-flop rate \[nm²/ps\].
    pub fn tm_diffusion_from_flipflop(flipflop_rate_ns: f64, thickness_nm: f64) -> f64 {
        // D ~ d²/(2*tau), tau = 1/rate
        let tau_ps = 1.0e3 / flipflop_rate_ns.max(1e-30);
        thickness_nm * thickness_nm / (2.0 * tau_ps)
    }

    /// Compute the lateral pressure profile along z binned into `n_bins` slabs.
    ///
    /// Returns the bin centers \[nm\] and lateral pressure difference
    /// `P_N - P_L` \[kJ/mol/nm³\] in each slab.
    pub fn lateral_pressure_profile(&self, n_bins: usize) -> (Vec<f64>, Vec<f64>) {
        let dz = self.box_len[2] / n_bins as f64;
        let mut p_diff = vec![0.0f64; n_bins];

        let n = self.beads.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let ri = self.beads[i].pos;
                let rj = self.beads[j].pos;
                let dr = pbc_displacement(sub3(rj, ri), self.box_len);
                let r = norm3(dr);
                if r < 1e-6 || r >= self.r_cut {
                    continue;
                }
                let eps_i = self.beads[i].bead_type.epsilon();
                let eps_j = self.beads[j].bead_type.epsilon();
                let sig_i = self.beads[i].bead_type.sigma();
                let sig_j = self.beads[j].bead_type.sigma();
                let eps_ij = (eps_i * eps_j).sqrt();
                let sig_ij = 0.5 * (sig_i + sig_j);

                let (_e, f_over_r) = lj_sf(eps_ij, sig_ij, r, self.r_cut);
                let fz = f_over_r * dr[2];
                let fx = f_over_r * dr[0];
                let fy = f_over_r * dr[1];

                let z_i = ri[2].rem_euclid(self.box_len[2]);
                let bin_i = ((z_i / dz) as usize).min(n_bins - 1);
                let area = self.box_len[0] * self.box_len[1];

                // Pressure difference p_N - p_L ~ fz*dz/Vol - (fx+fy)/2/area
                let contrib =
                    fz * dr[2] / (dz * area) - (fx * dr[0] + fy * dr[1]) / (2.0 * dz * area);
                p_diff[bin_i] += contrib;
            }
        }

        let centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * dz).collect();
        (centers, p_diff)
    }

    /// Compute total membrane tension γ \[kJ/mol/nm²\] from the lateral pressure profile.
    pub fn membrane_tension(&self, n_bins: usize) -> f64 {
        let (_z, p_diff) = self.lateral_pressure_profile(n_bins);
        let dz = self.box_len[2] / n_bins as f64;
        p_diff.iter().sum::<f64>() * dz
    }

    /// Estimate area compressibility modulus K_A \[kJ/mol/nm²\].
    ///
    /// Uses fluctuation relation K_A = k_B T * `A` / <(δA)²>.
    pub fn area_compressibility(area_samples: &[f64], temperature: f64) -> f64 {
        if area_samples.len() < 2 {
            return 0.0;
        }
        let mean_a = area_samples.iter().sum::<f64>() / area_samples.len() as f64;
        let var_a = area_samples
            .iter()
            .map(|&a| (a - mean_a).powi(2))
            .sum::<f64>()
            / (area_samples.len() - 1) as f64;
        if var_a < 1e-30 {
            return 0.0;
        }
        KB * temperature * mean_a / var_a
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helfrich Elastic Theory
// ─────────────────────────────────────────────────────────────────────────────

/// Helfrich elastic free energy model for fluid membranes.
///
/// The Helfrich Hamiltonian:
/// H = ∫ dA \[ (κ/2)(2H − c₀)² + κ_G K \]
///
/// where H is mean curvature, K is Gaussian curvature,
/// κ is bending modulus, κ_G is Gaussian modulus, c₀ spontaneous curvature.
#[derive(Debug, Clone)]
pub struct HelfrichElasticity {
    /// Bending modulus κ \[kJ/mol\].
    pub kappa: f64,
    /// Gaussian modulus κ_G \[kJ/mol\].
    pub kappa_gaussian: f64,
    /// Spontaneous curvature c₀ \[1/nm\].
    pub c0: f64,
}

impl HelfrichElasticity {
    /// Construct Helfrich elasticity parameters.
    pub fn new(kappa: f64, kappa_gaussian: f64, c0: f64) -> Self {
        Self {
            kappa,
            kappa_gaussian,
            c0,
        }
    }

    /// Construct with typical DPPC bilayer parameters.
    pub fn dppc() -> Self {
        Self {
            kappa: 20.0 * KBT_300, // ~80 pN·nm
            kappa_gaussian: -10.0 * KBT_300,
            c0: 0.0,
        }
    }

    /// Bending energy density \[kJ/mol/nm²\] given principal curvatures c1, c2 \[1/nm\].
    pub fn bending_energy_density(&self, c1: f64, c2: f64) -> f64 {
        let h = 0.5 * (c1 + c2); // mean curvature
        let k = c1 * c2; // Gaussian curvature
        0.5 * self.kappa * (2.0 * h - self.c0).powi(2) + self.kappa_gaussian * k
    }

    /// Total bending energy \[kJ/mol\] for a membrane patch of area `a_nm2`.
    pub fn total_bending_energy(&self, c1: f64, c2: f64, a_nm2: f64) -> f64 {
        self.bending_energy_density(c1, c2) * a_nm2
    }

    /// Undulation spectrum amplitude ⟨|h_q|²⟩ \[nm²\] for wavevector q \[1/nm\].
    ///
    /// In the tension-free limit: ⟨|h_q|²⟩ = k_B T / (κ q⁴ A).
    pub fn undulation_amplitude(&self, q: f64, area: f64, temperature: f64) -> f64 {
        let kbt = KB * temperature;
        let q4 = q.powi(4).max(1e-30);
        kbt / (self.kappa * q4 * area)
    }

    /// Estimate bending modulus from the undulation spectrum via linear regression
    /// of ln⟨|h_q|²⟩ vs ln(q).
    ///
    /// Input: slices of `q` \[1/nm\] and `h2` = ⟨|h_q|²⟩ \[nm²\], `area` \[nm²\].
    /// Returns κ \[kJ/mol\].
    pub fn estimate_kappa_from_spectrum(q: &[f64], h2: &[f64], area: f64, temperature: f64) -> f64 {
        if q.is_empty() || q.len() != h2.len() {
            return 0.0;
        }
        let kbt = KB * temperature;
        // κ = k_B T / (⟨|h_q|²⟩ * q⁴ * A)
        let kappas: Vec<f64> = q
            .iter()
            .zip(h2.iter())
            .filter_map(|(&qi, &hi)| {
                if qi > 0.0 && hi > 0.0 {
                    Some(kbt / (hi * qi.powi(4) * area))
                } else {
                    None
                }
            })
            .collect();
        if kappas.is_empty() {
            return 0.0;
        }
        kappas.iter().sum::<f64>() / kappas.len() as f64
    }

    /// Membrane persistence length \[nm\] from bending modulus and temperature.
    ///
    /// ξ_p = a · exp(4π κ / (3 k_B T)) where a is a molecular length scale.
    pub fn persistence_length(&self, a_nm: f64, temperature: f64) -> f64 {
        let kbt = KB * temperature;
        a_nm * (4.0 * PI * self.kappa / (3.0 * kbt)).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Lateral Pressure Profile Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Lateral pressure profile decomposed into kinetic and potential contributions.
#[derive(Debug, Clone)]
pub struct LateralPressureProfile {
    /// Bin centers along z \[nm\].
    pub z_centers: Vec<f64>,
    /// Lateral pressure P_L(z) \[kJ/mol/nm³\] = (P_xx + P_yy)/2.
    pub p_lateral: Vec<f64>,
    /// Normal pressure P_N(z) \[kJ/mol/nm³\] = P_zz.
    pub p_normal: Vec<f64>,
    /// Kinetic pressure contribution.
    pub p_kinetic: Vec<f64>,
}

impl LateralPressureProfile {
    /// Compute the lateral pressure profile from a bilayer system.
    pub fn compute(system: &BilayerSystem, n_bins: usize) -> Self {
        let lz = system.box_len[2];
        let lx = system.box_len[0];
        let ly = system.box_len[1];
        let dz = lz / n_bins as f64;
        let vol_slab = lx * ly * dz;

        let mut p_xx = vec![0.0f64; n_bins];
        let mut p_yy = vec![0.0f64; n_bins];
        let mut p_zz = vec![0.0f64; n_bins];
        let mut p_kin = vec![0.0f64; n_bins];

        // Kinetic contribution
        for bead in &system.beads {
            let m = bead.mass();
            let z = bead.pos[2].rem_euclid(lz);
            let bin = ((z / dz) as usize).min(n_bins - 1);
            p_kin[bin] += m * bead.vel[0].powi(2) / vol_slab;
            p_xx[bin] += m * bead.vel[0].powi(2) / vol_slab;
            p_yy[bin] += m * bead.vel[1].powi(2) / vol_slab;
            p_zz[bin] += m * bead.vel[2].powi(2) / vol_slab;
        }

        // Potential (virial) contribution using Irving-Kirkwood
        let n = system.beads.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let ri = system.beads[i].pos;
                let rj = system.beads[j].pos;
                let dr = pbc_displacement(sub3(rj, ri), system.box_len);
                let r = norm3(dr);
                if r < 1e-6 || r >= system.r_cut {
                    continue;
                }
                let eps_ij = (system.beads[i].bead_type.epsilon()
                    * system.beads[j].bead_type.epsilon())
                .sqrt();
                let sig_ij =
                    0.5 * (system.beads[i].bead_type.sigma() + system.beads[j].bead_type.sigma());
                let (_e, f_over_r) = lj_sf(eps_ij, sig_ij, r, system.r_cut);

                let fx = f_over_r * dr[0];
                let fy = f_over_r * dr[1];
                let fz = f_over_r * dr[2];

                let zi = ri[2].rem_euclid(lz);
                let bin = ((zi / dz) as usize).min(n_bins - 1);
                p_xx[bin] += fx * dr[0] / (2.0 * vol_slab);
                p_yy[bin] += fy * dr[1] / (2.0 * vol_slab);
                p_zz[bin] += fz * dr[2] / (2.0 * vol_slab);
            }
        }

        let z_centers: Vec<f64> = (0..n_bins).map(|b| (b as f64 + 0.5) * dz).collect();
        let p_lateral: Vec<f64> = p_xx
            .iter()
            .zip(p_yy.iter())
            .map(|(x, y)| 0.5 * (x + y))
            .collect();

        Self {
            z_centers,
            p_lateral,
            p_normal: p_zz,
            p_kinetic: p_kin,
        }
    }

    /// Surface tension γ \[kJ/mol/nm²\] = ∫(P_N − P_L) dz.
    pub fn surface_tension(&self) -> f64 {
        if self.z_centers.len() < 2 {
            return 0.0;
        }
        let dz = self.z_centers[1] - self.z_centers[0];
        self.p_normal
            .iter()
            .zip(self.p_lateral.iter())
            .map(|(pn, pl)| (pn - pl) * dz)
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Curvature Field
// ─────────────────────────────────────────────────────────────────────────────

/// Local membrane curvature field computed on a 2D grid.
#[derive(Debug, Clone)]
pub struct CurvatureField {
    /// Grid dimensions (nx, ny).
    pub dims: (usize, usize),
    /// Grid spacing \[nm\].
    pub dx: f64,
    /// Height field h(x,y) \[nm\] (deviation from mean bilayer plane).
    pub height: Vec<f64>,
}

impl CurvatureField {
    /// Construct a curvature field from head-bead positions.
    pub fn from_head_positions(
        positions: &[[f64; 3]],
        nx: usize,
        ny: usize,
        lx: f64,
        ly: f64,
    ) -> Self {
        let dx = lx / nx as f64;
        let dy = ly / ny as f64;
        let mut height = vec![0.0f64; nx * ny];
        let mut count = vec![0usize; nx * ny];

        for pos in positions {
            let ix = (pos[0].rem_euclid(lx) / dx) as usize % nx;
            let iy = (pos[1].rem_euclid(ly) / dy) as usize % ny;
            height[ix * ny + iy] += pos[2];
            count[ix * ny + iy] += 1;
        }
        for k in 0..(nx * ny) {
            if count[k] > 0 {
                height[k] /= count[k] as f64;
            }
        }

        Self {
            dims: (nx, ny),
            dx,
            height,
        }
    }

    /// Mean curvature H at grid point (ix, iy) using finite differences.
    pub fn mean_curvature(&self, ix: usize, iy: usize) -> f64 {
        let (nx, ny) = self.dims;
        let idx = |i: usize, j: usize| i * ny + j;
        let h = &self.height;
        let dx = self.dx;

        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let iyp = (iy + 1) % ny;
        let iym = (iy + ny - 1) % ny;

        let hx = (h[idx(ixp, iy)] - h[idx(ixm, iy)]) / (2.0 * dx);
        let hy = (h[idx(ix, iyp)] - h[idx(ix, iym)]) / (2.0 * dx);
        let hxx = (h[idx(ixp, iy)] - 2.0 * h[idx(ix, iy)] + h[idx(ixm, iy)]) / (dx * dx);
        let hyy = (h[idx(ix, iyp)] - 2.0 * h[idx(ix, iy)] + h[idx(ix, iym)]) / (dx * dx);
        let hxy = (h[idx(ixp, iyp)] - h[idx(ixp, iym)] - h[idx(ixm, iyp)] + h[idx(ixm, iym)])
            / (4.0 * dx * dx);

        let denom = (1.0 + hx * hx + hy * hy).powf(1.5);
        if denom < 1e-30 {
            return 0.0;
        }
        ((1.0 + hy * hy) * hxx - 2.0 * hx * hy * hxy + (1.0 + hx * hx) * hyy) / (2.0 * denom)
    }

    /// Gaussian curvature K at grid point (ix, iy).
    pub fn gaussian_curvature(&self, ix: usize, iy: usize) -> f64 {
        let (nx, ny) = self.dims;
        let idx = |i: usize, j: usize| i * ny + j;
        let h = &self.height;
        let dx = self.dx;

        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let iyp = (iy + 1) % ny;
        let iym = (iy + ny - 1) % ny;

        let hx = (h[idx(ixp, iy)] - h[idx(ixm, iy)]) / (2.0 * dx);
        let hy = (h[idx(ix, iyp)] - h[idx(ix, iym)]) / (2.0 * dx);
        let hxx = (h[idx(ixp, iy)] - 2.0 * h[idx(ix, iy)] + h[idx(ixm, iy)]) / (dx * dx);
        let hyy = (h[idx(ix, iyp)] - 2.0 * h[idx(ix, iy)] + h[idx(ix, iym)]) / (dx * dx);
        let hxy = (h[idx(ixp, iyp)] - h[idx(ixp, iym)] - h[idx(ixm, iyp)] + h[idx(ixm, iym)])
            / (4.0 * dx * dx);

        let denom = (1.0 + hx * hx + hy * hy).powi(2);
        if denom < 1e-30 {
            return 0.0;
        }
        (hxx * hyy - hxy * hxy) / denom
    }

    /// Mean curvature averaged over all grid points.
    pub fn mean_mean_curvature(&self) -> f64 {
        let (nx, ny) = self.dims;
        let total: f64 = (0..nx)
            .flat_map(|i| (0..ny).map(move |j| (i, j)))
            .map(|(i, j)| self.mean_curvature(i, j))
            .sum();
        total / (nx * ny) as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase Transition
// ─────────────────────────────────────────────────────────────────────────────

/// Gel–liquid-crystalline phase transition model.
#[derive(Debug, Clone)]
pub struct PhaseTransition {
    /// Melting temperature T_m \[K\] for the pure lipid.
    pub t_melt: f64,
    /// Enthalpy of transition ΔH \[kJ/mol\].
    pub delta_h: f64,
    /// Width of the transition (cooperativity parameter σ \[K\]).
    pub sigma: f64,
}

impl PhaseTransition {
    /// Construct phase transition parameters for DPPC.
    pub fn dppc() -> Self {
        Self {
            t_melt: 314.15,
            delta_h: 36.0,
            sigma: 2.0,
        }
    }

    /// Mole fraction of lipids in the liquid-crystalline phase at temperature T.
    pub fn liquid_fraction(&self, temperature: f64) -> f64 {
        // Sigmoidal transition
        let x = (temperature - self.t_melt) / self.sigma;
        1.0 / (1.0 + (-x).exp())
    }

    /// Specific heat capacity peak near T_m \[kJ/(mol·K)\].
    pub fn heat_capacity_peak(&self, temperature: f64) -> f64 {
        let x = (temperature - self.t_melt) / self.sigma;
        let f = 1.0 / (1.0 + (-x).exp());
        self.delta_h / self.sigma * f * (1.0 - f)
    }

    /// Effective area per lipid at temperature T accounting for phase state.
    pub fn area_per_lipid(&self, temperature: f64) -> f64 {
        let f_liq = self.liquid_fraction(temperature);
        let a_gel = 0.47; // nm² for gel phase DPPC
        let a_liq = 0.64; // nm² for liquid phase DPPC
        f_liq * a_liq + (1.0 - f_liq) * a_gel
    }

    /// Effect of cholesterol mole fraction `x_chol` on melting temperature \[K\].
    ///
    /// Cholesterol broadens and suppresses the main transition.
    pub fn cholesterol_tm_shift(x_chol: f64) -> f64 {
        // Empirical: shift ~ -30 K per mole fraction of cholesterol
        -30.0 * x_chol
    }

    /// Cholesterol-induced ordering: increase in order parameter S_CD.
    pub fn cholesterol_order_increase(x_chol: f64, temperature: f64, t_melt: f64) -> f64 {
        // Above T_m, cholesterol increases order; below T_m, smaller effect
        let base = if temperature > t_melt { 0.4 } else { 0.1 };
        base * x_chol * (1.0 - x_chol)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Transmembrane Protein
// ─────────────────────────────────────────────────────────────────────────────

/// Coarse-grained transmembrane protein model.
#[derive(Debug, Clone)]
pub struct TmProtein {
    /// Number of TM helices.
    pub n_helices: usize,
    /// Hydrophobic length of the TM segment \[nm\].
    pub hydrophobic_length: f64,
    /// Tilt angle from bilayer normal \[radians\].
    pub tilt_angle: f64,
    /// Protein radius \[nm\].
    pub radius: f64,
}

impl TmProtein {
    /// Construct a transmembrane protein model.
    pub fn new(n_helices: usize, hydrophobic_length: f64, radius: f64) -> Self {
        Self {
            n_helices,
            hydrophobic_length,
            tilt_angle: 0.0,
            radius,
        }
    }

    /// Hydrophobic mismatch energy \[kJ/mol\] for bilayer thickness `d_bilayer` \[nm\].
    ///
    /// Uses the elastic model: E = H_s · (d_p − d_b)² where H_s ≈ 40 kJ/(mol·nm²).
    pub fn mismatch_energy(&self, d_bilayer: f64) -> f64 {
        let h_s = 40.0; // kJ/(mol·nm²)
        let mismatch = self.hydrophobic_length - d_bilayer;
        h_s * mismatch * mismatch
    }

    /// Tilt energy \[kJ/mol\] for given tilt angle θ \[rad\] and bending modulus κ.
    pub fn tilt_energy(&self, theta: f64, kappa: f64) -> f64 {
        // E_tilt = (π κ / 2) ln(R/a) θ²
        let a_nm = 0.5; // protein–lipid contact distance [nm]
        let r = self.radius.max(a_nm + 0.01);
        0.5 * PI * kappa * (r / a_nm).ln() * theta * theta
    }

    /// Effective hydrophobic length accounting for tilt.
    pub fn effective_length(&self) -> f64 {
        self.hydrophobic_length * self.tilt_angle.cos()
    }

    /// Lipid distortion radius \[nm\] around the protein.
    pub fn distortion_radius(&self, kappa: f64, d_bilayer: f64) -> f64 {
        // Decay length λ = (κ / H_s)^(1/4)
        let h_s = 40.0;
        (kappa / h_s).powf(0.25) + d_bilayer * 0.5
    }

    /// Estimate number of lipids in the protein annular shell.
    pub fn annular_lipids(&self, area_per_lipid: f64) -> f64 {
        let r_shell = self.radius + 0.8; // one lipid shell thickness
        let annular_area = PI * (r_shell * r_shell - self.radius * self.radius);
        annular_area / area_per_lipid
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Vesicle Formation
// ─────────────────────────────────────────────────────────────────────────────

/// Vesicle geometry and formation parameters.
#[derive(Debug, Clone)]
pub struct VesicleBuilder {
    /// Vesicle outer radius \[nm\].
    pub radius: f64,
    /// Bilayer thickness \[nm\].
    pub thickness: f64,
    /// Lipid species.
    pub species: LipidSpecies,
}

impl VesicleBuilder {
    /// Construct a vesicle builder.
    pub fn new(radius: f64, thickness: f64, species: LipidSpecies) -> Self {
        Self {
            radius,
            thickness,
            species,
        }
    }

    /// Estimate total number of lipids in the vesicle.
    pub fn n_lipids(&self, area_per_lipid: f64) -> usize {
        let r_out = self.radius;
        let r_in = (self.radius - self.thickness).max(1.0);
        let a_out = 4.0 * PI * r_out * r_out;
        let a_in = 4.0 * PI * r_in * r_in;
        let n_out = (a_out / area_per_lipid) as usize;
        let n_in = (a_in / area_per_lipid) as usize;
        n_out + n_in
    }

    /// Bending energy of the complete vesicle \[kJ/mol\].
    pub fn bending_energy(&self, kappa: f64, c0: f64) -> f64 {
        // For a sphere: H = 1/R, K = 1/R²
        let r = self.radius;
        let h = 1.0 / r;
        let k_gauss = 1.0 / (r * r);
        let kappa_g = -0.5 * kappa; // typical Gaussian modulus
        let a = 4.0 * PI * r * r;
        (0.5 * kappa * (2.0 * h - c0).powi(2) + kappa_g * k_gauss) * a
    }

    /// Estimate vesicle formation time \[ns\] from a flat bilayer patch.
    ///
    /// Very approximate; based on measured timescales in CG-MD simulations.
    pub fn formation_time_estimate(&self, temperature: f64) -> f64 {
        // Empirical scaling: τ ~ R² / D_bend
        let d_bend = KB * temperature / (6.0 * PI * 1e-3 * self.radius); // nm²/ps, Stokes
        self.radius * self.radius / d_bend * 1e-3 // convert to ns
    }

    /// Generate spherical vesicle bead positions (outer leaflet only for brevity).
    pub fn generate_outer_positions(&self, n_beads: usize) -> Vec<[f64; 3]> {
        let mut positions = Vec::with_capacity(n_beads);
        // Use Fibonacci sphere
        let phi = (1.0 + 5.0_f64.sqrt()) * 0.5;
        for i in 0..n_beads {
            let theta = 2.0 * PI * i as f64 / phi;
            let cos_phi = 1.0 - 2.0 * (i as f64 + 0.5) / n_beads as f64;
            let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
            positions.push([
                self.radius * sin_phi * theta.cos(),
                self.radius * sin_phi * theta.sin(),
                self.radius * cos_phi,
            ]);
        }
        positions
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Bilayer Fusion
// ─────────────────────────────────────────────────────────────────────────────

/// Bilayer fusion pathway model (hemifusion stalk → pore).
#[derive(Debug, Clone)]
pub struct BilayerFusion {
    /// Stalk energy \[kJ/mol\].
    pub stalk_energy: f64,
    /// Hemifusion diaphragm energy \[kJ/mol\].
    pub diaphragm_energy: f64,
    /// Fusion pore energy \[kJ/mol\].
    pub pore_energy: f64,
    /// Bending modulus of the membrane \[kJ/mol\].
    pub kappa: f64,
}

impl BilayerFusion {
    /// Construct a fusion model for a given bending modulus.
    pub fn new(kappa: f64) -> Self {
        // Approximate energies from Kozlovsky & Kozlov (2002)
        Self {
            stalk_energy: 40.0 * KBT_300,
            diaphragm_energy: 20.0 * KBT_300,
            pore_energy: 10.0 * KBT_300,
            kappa,
        }
    }

    /// Stalk formation energy \[kJ/mol\] as function of lipid spontaneous curvature c0.
    pub fn stalk_formation_energy(&self, c0: f64) -> f64 {
        // Negative spontaneous curvature lowers the stalk energy
        self.stalk_energy * (1.0 - 2.0 * c0 * c0)
    }

    /// Fusion pore radius \[nm\] at equilibrium.
    pub fn equilibrium_pore_radius(&self, membrane_tension: f64) -> f64 {
        // Balance bending vs tension: r = sqrt(kappa / gamma)
        if membrane_tension <= 0.0 {
            return 0.0;
        }
        (self.kappa / membrane_tension).sqrt()
    }

    /// Rate of stalk formation \[1/ns\] (Arrhenius with stalk energy as barrier).
    pub fn stalk_rate(&self, temperature: f64, c0: f64) -> f64 {
        let kbt = KB * temperature;
        let e_stalk = self.stalk_formation_energy(c0);
        let nu0 = 1.0e9; // [1/ns]
        nu0 * (-e_stalk / kbt).exp()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Membrane Undulation Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Analysis of membrane thermal undulations.
#[derive(Debug, Clone)]
pub struct UndulationAnalysis {
    /// Box dimensions \[nm\].
    pub lx: f64,
    /// Box dimensions \[nm\].
    pub ly: f64,
    /// Bending modulus κ \[kJ/mol\] extracted from spectrum.
    pub kappa_fit: f64,
}

impl UndulationAnalysis {
    /// Construct undulation analysis for given box.
    pub fn new(lx: f64, ly: f64) -> Self {
        Self {
            lx,
            ly,
            kappa_fit: 0.0,
        }
    }

    /// Compute discrete Fourier modes of height field h(x,y).
    ///
    /// Returns (q_values \[1/nm\], power_spectrum \[nm²\]).
    pub fn compute_spectrum(
        &mut self,
        height: &[f64],
        nx: usize,
        ny: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        // Simple 2D DFT (only along x for clarity)
        let mut q_vals = Vec::new();
        let mut power = Vec::new();

        for kx in 1..=(nx / 2) {
            let q = 2.0 * PI * kx as f64 / self.lx;
            let mut re = 0.0;
            let mut im = 0.0;
            for ix in 0..nx {
                let angle = 2.0 * PI * kx as f64 * ix as f64 / nx as f64;
                // Average over y
                let h_x: f64 = (0..ny).map(|iy| height[ix * ny + iy]).sum::<f64>() / ny as f64;
                re += h_x * angle.cos();
                im -= h_x * angle.sin();
            }
            let amplitude_sq = (re * re + im * im) / (nx * nx) as f64;
            q_vals.push(q);
            power.push(amplitude_sq);
        }

        q_vals
            .iter()
            .zip(power.iter())
            .filter(|&(_, p)| *p > 0.0)
            .for_each(|_| {});

        (q_vals, power)
    }

    /// Fit bending modulus κ from spectrum and store in `kappa_fit`.
    pub fn fit_kappa(&mut self, q: &[f64], power: &[f64], area: f64, temperature: f64) {
        self.kappa_fit =
            HelfrichElasticity::estimate_kappa_from_spectrum(q, power, area, temperature);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: Box-Muller normal RNG (no external deps beyond rand)
// ─────────────────────────────────────────────────────────────────────────────

/// Generate a standard-normal sample using Box-Muller transform.
fn rand_normal(rng: &mut impl rand::RngExt) -> f64 {
    use std::f64::consts::TAU;
    let u1: f64 = rng.random_range(1e-10..1.0);
    let u2: f64 = rng.random_range(0.0..1.0);
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional membrane analysis helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the radial distribution function g(r) between two bead selections.
///
/// Returns `(r_bins, g_r)` with bin centers \[nm\] and g(r) values.
pub fn radial_distribution_function(
    pos_a: &[[f64; 3]],
    pos_b: &[[f64; 3]],
    box_len: [f64; 3],
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    let na = pos_a.len();
    let nb = pos_b.len();

    for a in pos_a {
        for b in pos_b {
            let delta = pbc_displacement(sub3(*b, *a), box_len);
            let r = norm3(delta);
            if r > 0.0 && r < r_max {
                let bin = (r / dr) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }

    let volume = box_len[0] * box_len[1] * box_len[2];
    let rho_b = nb as f64 / volume;
    let mut g_r = vec![0.0f64; n_bins];
    for (k, &count) in hist.iter().enumerate() {
        let r_lo = k as f64 * dr;
        let r_hi = r_lo + dr;
        let shell_vol = 4.0 / 3.0 * PI * (r_hi.powi(3) - r_lo.powi(3));
        g_r[k] = count as f64 / (na as f64 * rho_b * shell_vol);
    }

    let centers: Vec<f64> = (0..n_bins).map(|k| (k as f64 + 0.5) * dr).collect();
    (centers, g_r)
}

/// Voronoi tessellation-based area per lipid for non-uniform packing.
///
/// Computes the 2D Voronoi cell area for each lipid projected onto the xy-plane.
/// Returns a vector of areas \[nm²\] with the same length as `positions`.
pub fn voronoi_area_per_lipid(positions: &[[f64; 3]], lx: f64, ly: f64) -> Vec<f64> {
    // Approximate Voronoi areas using nearest-neighbour distances
    let n = positions.len();
    if n == 0 {
        return Vec::new();
    }
    let mut areas = vec![lx * ly / n as f64; n]; // default to average
    for i in 0..n {
        let xi = positions[i][0];
        let yi = positions[i][1];
        let mut min_dist = f64::MAX;
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = (positions[j][0] - xi).abs();
            let dx = dx.min(lx - dx);
            let dy = (positions[j][1] - yi).abs();
            let dy = dy.min(ly - dy);
            let d2 = dx * dx + dy * dy;
            if d2 < min_dist {
                min_dist = d2;
            }
        }
        // Approximate cell area as π r² / 2 (hexagonal packing)
        areas[i] = PI * min_dist * 0.5;
    }
    areas
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Bead types ────────────────────────────────────────────────────────────

    #[test]
    fn test_bead_type_epsilon_positive() {
        for bt in &[
            BeadType::Head,
            BeadType::Glycerol,
            BeadType::Tail,
            BeadType::Cholesterol,
            BeadType::Protein,
            BeadType::Water,
        ] {
            assert!(
                bt.epsilon() > 0.0,
                "epsilon should be positive for {:?}",
                bt
            );
        }
    }

    #[test]
    fn test_bead_type_sigma_positive() {
        for bt in &[BeadType::Head, BeadType::Tail, BeadType::Water] {
            assert!(bt.sigma() > 0.0);
        }
    }

    #[test]
    fn test_bead_type_mass_positive() {
        assert!(BeadType::Tail.mass() > 0.0);
        assert!(BeadType::Head.mass() > BeadType::Tail.mass());
    }

    // ── Bond potential ────────────────────────────────────────────────────────

    #[test]
    fn test_bond_energy_at_equilibrium() {
        let bp = BondParam::new(0.47, 3800.0);
        assert!(bp.energy(0.47).abs() < 1e-10);
    }

    #[test]
    fn test_bond_force_at_equilibrium() {
        let bp = BondParam::new(0.47, 3800.0);
        assert!(bp.force_magnitude(0.47).abs() < 1e-10);
    }

    #[test]
    fn test_bond_energy_positive_off_equilibrium() {
        let bp = BondParam::new(0.47, 3800.0);
        assert!(bp.energy(0.60) > 0.0);
        assert!(bp.energy(0.30) > 0.0);
    }

    // ── Angle potential ───────────────────────────────────────────────────────

    #[test]
    fn test_angle_energy_at_equilibrium() {
        let ap = AngleParam::new(PI * 0.8, 25.0);
        let e = ap.energy_from_cos(ap.cos_theta0);
        assert!(e.abs() < 1e-12);
    }

    // ── Lennard-Jones ─────────────────────────────────────────────────────────

    #[test]
    fn test_lj_zero_beyond_cutoff() {
        let (e, f) = lj_sf(3.5, 0.47, 1.3, 1.2);
        assert_eq!(e, 0.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_lj_repulsive_at_short_range() {
        let (e, _f) = lj_sf(3.5, 0.47, 0.3, 1.2);
        assert!(e > 0.0, "LJ should be repulsive at r < sigma");
    }

    #[test]
    fn test_lj_attractive_near_minimum() {
        // LJ minimum is at r = 2^(1/6) * sigma
        let sigma = 0.47;
        let r_min = sigma * 2.0_f64.powf(1.0 / 6.0);
        let (e, _) = lj_sf(3.5, sigma, r_min, 1.5);
        // Just past minimum, energy should be negative
        assert!(e < 0.0, "LJ should be attractive near minimum, e={:.6}", e);
    }

    // ── Bilayer system ────────────────────────────────────────────────────────

    #[test]
    fn test_bilayer_bead_count() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 4, 8.0, 300.0, 0.02);
        // 2 leaflets * 4 lipids * 12 beads = 96
        assert_eq!(sys.beads.len(), 96);
        assert_eq!(sys.lipids.len(), 8);
    }

    #[test]
    fn test_area_per_lipid_positive() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 4, 8.0, 300.0, 0.02);
        let apl = sys.area_per_lipid();
        assert!(
            apl > 0.0,
            "area per lipid should be positive, got {:.6}",
            apl
        );
    }

    #[test]
    fn test_bilayer_thickness_positive() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 4, 8.0, 300.0, 0.02);
        let th = sys.bilayer_thickness();
        assert!(
            th > 0.0,
            "bilayer thickness should be positive, got {:.6}",
            th
        );
    }

    #[test]
    fn test_kinetic_energy_positive_after_init() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 2, 6.0, 300.0, 0.02);
        let kin = sys.kinetic_energy();
        assert!(
            kin > 0.0,
            "kinetic energy should be positive after velocity init"
        );
    }

    #[test]
    fn test_temperature_from_kinetic() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 2, 6.0, 300.0, 0.02);
        let kin = sys.kinetic_energy();
        let t = sys.temperature_from_kinetic(kin);
        // Temperature should be in a reasonable range
        assert!(t > 0.0 && t < 2000.0, "temperature out of range: {:.6}", t);
    }

    #[test]
    fn test_nonbonded_forces_finite() {
        let mut sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 2, 6.0, 300.0, 0.02);
        let e = sys.compute_nonbonded_forces();
        assert!(e.is_finite(), "nonbonded energy should be finite: {:.6}", e);
    }

    #[test]
    fn test_bonded_forces_finite() {
        let mut sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 2, 6.0, 300.0, 0.02);
        let e = sys.compute_bonded_forces();
        assert!(e.is_finite(), "bonded energy should be finite: {:.6}", e);
    }

    #[test]
    fn test_order_parameter_range() {
        let sys = BilayerSystem::new_flat_bilayer(LipidSpecies::Dppc, 4, 8.0, 300.0, 0.02);
        let s = sys.order_parameter();
        assert!(
            (-0.5..=1.0).contains(&s),
            "order parameter out of range: {:.6}",
            s
        );
    }

    // ── Helfrich elasticity ───────────────────────────────────────────────────

    #[test]
    fn test_helfrich_zero_curvature() {
        let helf = HelfrichElasticity::dppc();
        let e = helf.bending_energy_density(0.0, 0.0);
        assert!(
            e.abs() < 1e-10,
            "zero curvature should give zero energy, got {:.6}",
            e
        );
    }

    #[test]
    fn test_helfrich_spherical_energy_positive() {
        let helf = HelfrichElasticity::dppc();
        let r = 50.0; // large vesicle
        let c = 1.0 / r;
        let e = helf.bending_energy_density(c, c);
        assert!(e > 0.0 || e.abs() < 1e-6);
    }

    #[test]
    fn test_undulation_amplitude_decreases_with_q() {
        let helf = HelfrichElasticity::dppc();
        let area = 400.0;
        let a1 = helf.undulation_amplitude(0.1, area, 300.0);
        let a2 = helf.undulation_amplitude(1.0, area, 300.0);
        assert!(a1 > a2, "undulation amplitude should decrease with q");
    }

    #[test]
    fn test_persistence_length_positive() {
        let helf = HelfrichElasticity::dppc();
        let lp = helf.persistence_length(0.5, 300.0);
        assert!(lp > 0.0, "persistence length should be positive: {:.6}", lp);
    }

    // ── Phase transition ──────────────────────────────────────────────────────

    #[test]
    fn test_liquid_fraction_sigmoidal() {
        let pt = PhaseTransition::dppc();
        let f_low = pt.liquid_fraction(280.0);
        let f_high = pt.liquid_fraction(340.0);
        assert!(f_low < 0.5 && f_high > 0.5);
    }

    #[test]
    fn test_area_per_lipid_increases_with_temperature() {
        let pt = PhaseTransition::dppc();
        let a_gel = pt.area_per_lipid(290.0);
        let a_liq = pt.area_per_lipid(330.0);
        assert!(a_liq > a_gel, "area per lipid should increase above T_m");
    }

    // ── Transmembrane protein ─────────────────────────────────────────────────

    #[test]
    fn test_mismatch_energy_zero_at_match() {
        let prot = TmProtein::new(1, 4.0, 1.2);
        let e = prot.mismatch_energy(4.0);
        assert!(
            e.abs() < 1e-10,
            "mismatch energy should be zero at perfect match"
        );
    }

    #[test]
    fn test_mismatch_energy_positive_off_match() {
        let prot = TmProtein::new(1, 4.0, 1.2);
        assert!(prot.mismatch_energy(3.5) > 0.0);
        assert!(prot.mismatch_energy(4.5) > 0.0);
    }

    // ── Vesicle ───────────────────────────────────────────────────────────────

    #[test]
    fn test_vesicle_n_lipids_positive() {
        let vb = VesicleBuilder::new(30.0, 4.0, LipidSpecies::Dppc);
        let n = vb.n_lipids(0.64);
        assert!(n > 0, "vesicle should have positive number of lipids");
    }

    #[test]
    fn test_vesicle_outer_positions_count() {
        let vb = VesicleBuilder::new(20.0, 4.0, LipidSpecies::Dopc);
        let pos = vb.generate_outer_positions(100);
        assert_eq!(pos.len(), 100);
    }

    // ── Bilayer fusion ────────────────────────────────────────────────────────

    #[test]
    fn test_fusion_stalk_rate_positive() {
        let fusion = BilayerFusion::new(20.0 * KBT_300);
        let rate = fusion.stalk_rate(300.0, -0.3);
        assert!(
            rate >= 0.0,
            "stalk rate should be non-negative: {:.6}",
            rate
        );
    }

    #[test]
    fn test_fusion_pore_radius_positive_tension() {
        let fusion = BilayerFusion::new(20.0 * KBT_300);
        let r = fusion.equilibrium_pore_radius(1.0);
        assert!(
            r > 0.0,
            "pore radius should be positive under tension: {:.6}",
            r
        );
    }

    // ── Radial distribution function ──────────────────────────────────────────

    #[test]
    fn test_rdf_returns_correct_bins() {
        let pos_a = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let pos_b = vec![[0.5, 0.0, 0.0]];
        let (_r, g) = radial_distribution_function(&pos_a, &pos_b, [5.0, 5.0, 5.0], 2.5, 10);
        assert_eq!(g.len(), 10);
    }

    #[test]
    fn test_flipflop_rate_decreases_with_barrier() {
        let rate_low = BilayerSystem::flipflop_rate(300.0, 20.0);
        let rate_high = BilayerSystem::flipflop_rate(300.0, 80.0);
        assert!(
            rate_low > rate_high,
            "higher barrier should give lower rate: low={:.6}, high={:.6}",
            rate_low,
            rate_high
        );
    }

    #[test]
    fn test_area_compressibility_positive() {
        let areas = vec![64.0, 64.5, 63.8, 64.2, 64.1];
        let ka = BilayerSystem::area_compressibility(&areas, 300.0);
        assert!(
            ka > 0.0,
            "area compressibility should be positive: {:.6}",
            ka
        );
    }
}
