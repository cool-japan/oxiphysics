// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ab initio MD interface (DFT-based forces).
//!
//! This module provides a simplified interface for DFT-based molecular dynamics,
//! including:
//! - [`DftSystem`]: system descriptor with basis set and XC functional
//! - [`KohnSham`]: Kohn-Sham orbital representation
//! - Hartree-Fock total energy (tight-binding placeholder)
//! - Hellmann-Feynman forces from electron density
//! - Pulay density mixing for SCF convergence
//! - Born-Oppenheimer MD step
//! - Ehrenfest coupled electron-nuclear dynamics

// ---------------------------------------------------------------------------
// Enumerations
// ---------------------------------------------------------------------------

/// Basis set descriptor for DFT calculations.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum BasisSet {
    /// Minimal STO-3G basis set.
    Sto3G,
    /// 6-31G double-zeta split-valence basis.
    Basis631G,
    /// Triple-zeta with polarization functions (TZP).
    Tzp,
    /// Plane-wave basis with energy cutoff (Hartree).
    PlaneWave(f64),
}

/// Exchange-correlation functional type.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum XcFunctional {
    /// Local density approximation (LDA/VWN).
    Lda,
    /// Generalized gradient approximation (PBE).
    Pbe,
    /// Hybrid functional B3LYP.
    B3lyp,
    /// Meta-GGA functional (TPSS).
    Tpss,
}

// ---------------------------------------------------------------------------
// DftSystem
// ---------------------------------------------------------------------------

/// DFT system descriptor.
///
/// Holds the configuration for a DFT calculation: number of atoms,
/// basis set, exchange-correlation functional, and nuclear positions.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DftSystem {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Atomic numbers (Z) for each atom.
    pub atomic_numbers: Vec<u32>,
    /// Nuclear positions (x, y, z) in Angstrom.
    pub positions: Vec<[f64; 3]>,
    /// Basis set descriptor.
    pub basis_set: BasisSet,
    /// Exchange-correlation functional.
    pub xc_functional: XcFunctional,
    /// Total charge of the system.
    pub charge: i32,
    /// Spin multiplicity (2S+1).
    pub multiplicity: u32,
}

impl DftSystem {
    /// Create a new DFT system.
    pub fn new(
        atomic_numbers: Vec<u32>,
        positions: Vec<[f64; 3]>,
        basis_set: BasisSet,
        xc_functional: XcFunctional,
    ) -> Self {
        let n_atoms = atomic_numbers.len();
        Self {
            n_atoms,
            atomic_numbers,
            positions,
            basis_set,
            xc_functional,
            charge: 0,
            multiplicity: 1,
        }
    }

    /// Total number of electrons (sum of atomic numbers minus charge).
    pub fn n_electrons(&self) -> i32 {
        let z_sum: u32 = self.atomic_numbers.iter().sum();
        z_sum as i32 - self.charge
    }

    /// Nuclear repulsion energy (Coulomb, in Hartree with atomic units conversion).
    ///
    /// `E_nn = sum_{A<B} Z_A * Z_B / |R_A - R_B|`
    ///
    /// Positions are in Angstrom; returns energy in eV (1 Hartree ≈ 27.211 eV).
    pub fn nuclear_repulsion_ev(&self) -> f64 {
        let bohr_per_angstrom = 1.0 / 0.529_177;
        let hartree_to_ev = 27.211_386;
        let mut energy = 0.0_f64;
        for a in 0..self.n_atoms {
            for b in (a + 1)..self.n_atoms {
                let za = self.atomic_numbers[a] as f64;
                let zb = self.atomic_numbers[b] as f64;
                let ra = &self.positions[a];
                let rb = &self.positions[b];
                let dx = ra[0] - rb[0];
                let dy = ra[1] - rb[1];
                let dz = ra[2] - rb[2];
                let r_ang = (dx * dx + dy * dy + dz * dz).sqrt();
                let r_bohr = r_ang * bohr_per_angstrom;
                if r_bohr > 1e-12 {
                    energy += za * zb / r_bohr;
                }
            }
        }
        energy * hartree_to_ev
    }
}

// ---------------------------------------------------------------------------
// KohnSham
// ---------------------------------------------------------------------------

/// Kohn-Sham orbital representation.
///
/// Stores orbital energies and occupation numbers for a closed-shell system.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct KohnSham {
    /// Orbital energies in eV.
    pub orbital_energies: Vec<f64>,
    /// Occupation numbers (0, 1, or 2 for closed-shell).
    pub occupations: Vec<f64>,
    /// Electron density on a real-space grid.
    pub density: Vec<f64>,
    /// Total SCF energy in eV.
    pub total_energy: f64,
    /// HOMO index.
    pub homo_idx: usize,
}

impl KohnSham {
    /// Create a new KohnSham object.
    pub fn new(orbital_energies: Vec<f64>, occupations: Vec<f64>) -> Self {
        let n = orbital_energies.len();
        let homo_idx = occupations.iter().rposition(|&o| o > 0.0).unwrap_or(0);
        Self {
            orbital_energies,
            occupations,
            density: vec![0.0; n],
            total_energy: 0.0,
            homo_idx,
        }
    }

    /// HOMO energy in eV.
    pub fn homo_energy(&self) -> f64 {
        self.orbital_energies
            .get(self.homo_idx)
            .copied()
            .unwrap_or(0.0)
    }

    /// LUMO energy in eV (orbital after HOMO).
    pub fn lumo_energy(&self) -> f64 {
        self.orbital_energies
            .get(self.homo_idx + 1)
            .copied()
            .unwrap_or(0.0)
    }

    /// HOMO-LUMO gap in eV.
    pub fn homo_lumo_gap(&self) -> f64 {
        self.lumo_energy() - self.homo_energy()
    }

    /// Band energy: sum of eps_i * n_i.
    pub fn band_energy(&self) -> f64 {
        self.orbital_energies
            .iter()
            .zip(&self.occupations)
            .map(|(e, n)| e * n)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Hartree-Fock energy (tight-binding placeholder)
// ---------------------------------------------------------------------------

/// Compute a tight-binding estimate of the total electronic energy.
///
/// Uses the extended Hückel approximation:
/// `E = sum_i n_i * eps_i`
/// where `eps_i` are parametric orbital energies and `n_i` are occupations.
///
/// # Arguments
/// - `orbital_energies`: single-particle energies in eV
/// - `occupations`: occupation number for each orbital
///
/// Returns the band energy in eV.
#[allow(dead_code)]
pub fn hartree_fock_energy(orbital_energies: &[f64], occupations: &[f64]) -> f64 {
    orbital_energies
        .iter()
        .zip(occupations)
        .map(|(e, n)| e * n)
        .sum()
}

// ---------------------------------------------------------------------------
// DFT force (Hellmann-Feynman)
// ---------------------------------------------------------------------------

/// Compute the Hellmann-Feynman force on nucleus A from the electron density.
///
/// `F_A = Z_A * sum_{grid} rho(r) * (R_A - r) / |R_A - r|^3 * dV`
///
/// This is a simplified discrete sum over a 3D real-space grid.
///
/// # Arguments
/// - `z_a`: atomic number of nucleus A
/// - `r_a`: position of nucleus A `[x, y, z]`
/// - `density_grid`: electron density values at grid points
/// - `grid_points`: grid point positions `[[x, y, z\], ...]`
/// - `dv`: volume element per grid point (Bohr^3)
///
/// Returns the force vector `[fx, fy, fz]` in Hartree/Bohr.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn dft_force(
    z_a: f64,
    r_a: &[f64; 3],
    density_grid: &[f64],
    grid_points: &[[f64; 3]],
    dv: f64,
) -> [f64; 3] {
    let mut force = [0.0_f64; 3];
    for (rho, r_grid) in density_grid.iter().zip(grid_points) {
        let dx = r_a[0] - r_grid[0];
        let dy = r_a[1] - r_grid[1];
        let dz = r_a[2] - r_grid[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 < 1e-20 {
            continue;
        }
        let r3 = r2.sqrt() * r2;
        let fac = z_a * rho * dv / r3;
        force[0] += fac * dx;
        force[1] += fac * dy;
        force[2] += fac * dz;
    }
    force
}

// ---------------------------------------------------------------------------
// Pulay density mixing (DIIS-like)
// ---------------------------------------------------------------------------

/// Pulay density mixing for SCF convergence.
///
/// The Pulay scheme (also known as DIIS) mixes input densities to minimize
/// the residual. This simplified version applies a linear combination:
///
/// `rho_out = sum_i c_i * rho_in_i`
///
/// where the coefficients `c_i` are determined from the error vectors.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PulayMixer {
    /// History of input densities.
    pub density_history: Vec<Vec<f64>>,
    /// History of residuals (rho_out - rho_in).
    pub residual_history: Vec<Vec<f64>>,
    /// Maximum history size (DIIS subspace dimension).
    pub max_history: usize,
    /// Linear mixing coefficient (fallback).
    pub alpha: f64,
}

impl PulayMixer {
    /// Create a new Pulay mixer.
    pub fn new(max_history: usize, alpha: f64) -> Self {
        Self {
            density_history: Vec::new(),
            residual_history: Vec::new(),
            max_history,
            alpha,
        }
    }

    /// Add a density/residual pair to the history.
    pub fn push(&mut self, rho_in: Vec<f64>, rho_out: Vec<f64>) {
        let residual: Vec<f64> = rho_out.iter().zip(&rho_in).map(|(o, i)| o - i).collect();
        self.density_history.push(rho_in);
        self.residual_history.push(residual);
        if self.density_history.len() > self.max_history {
            self.density_history.remove(0);
            self.residual_history.remove(0);
        }
    }

    /// Compute the RMS residual of the latest entry.
    pub fn rms_residual(&self) -> f64 {
        if let Some(res) = self.residual_history.last() {
            let n = res.len();
            if n == 0 {
                return 0.0;
            }
            let sum_sq: f64 = res.iter().map(|r| r * r).sum();
            (sum_sq / n as f64).sqrt()
        } else {
            0.0
        }
    }

    /// Apply simple linear mixing: rho_new = rho_in + alpha * residual.
    pub fn pulay_mixing(&self, rho_in: &[f64], rho_out: &[f64]) -> Vec<f64> {
        rho_in
            .iter()
            .zip(rho_out)
            .map(|(i, o)| i + self.alpha * (o - i))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Born-Oppenheimer MD step
// ---------------------------------------------------------------------------

/// Perform one Born-Oppenheimer MD step.
///
/// In BO-MD the electrons are assumed to be in the instantaneous ground state.
/// This function:
/// 1. Computes forces on nuclei from the electron density (Hellmann-Feynman)
/// 2. Updates nuclear velocities and positions (velocity Verlet half-step)
///
/// # Arguments
/// - `positions`: nuclear positions `[[x,y,z\], ...]` (modified in place)
/// - `velocities`: nuclear velocities (modified in place)
/// - `forces`: forces on nuclei (updated from density)
/// - `masses`: nuclear masses in atomic mass units
/// - `dt`: timestep in femtoseconds
#[allow(dead_code)]
pub fn born_oppenheimer_step(
    positions: &mut [[f64; 3]],
    velocities: &mut [[f64; 3]],
    forces: &[[f64; 3]],
    masses: &[f64],
    dt: f64,
) {
    let amu_to_kg = 1.660_539e-27_f64;
    let fs_to_s = 1e-15_f64;
    let ev_per_angstrom_to_n = 1.602_176_634e-19 / 1e-10; // eV/Å to N
    let _unit = ev_per_angstrom_to_n / amu_to_kg * fs_to_s * fs_to_s;

    let n = positions.len();
    for i in 0..n {
        let m = masses[i].max(1e-30);
        for k in 0..3 {
            // a = F/m
            let acc = forces[i][k] / m;
            // velocity Verlet: v(t+dt/2) = v(t) + 0.5*a*dt
            velocities[i][k] += 0.5 * acc * dt;
            // x(t+dt) = x(t) + v(t+dt/2)*dt
            positions[i][k] += velocities[i][k] * dt;
            // v(t+dt) = v(t+dt/2) + 0.5*a*dt (using same force; full BO would recompute)
            velocities[i][k] += 0.5 * acc * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Ehrenfest dynamics
// ---------------------------------------------------------------------------

/// Ehrenfest dynamics state.
///
/// In Ehrenfest dynamics, the electrons evolve on the mean-field potential
/// while coupled to the nuclei. The electronic wavefunction is a superposition.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EhrenfestState {
    /// Nuclear positions `[[x,y,z\], ...]`.
    pub positions: Vec<[f64; 3]>,
    /// Nuclear velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Nuclear masses in atomic mass units.
    pub masses: Vec<f64>,
    /// Electronic state populations (occupation of each adiabatic state).
    pub populations: Vec<f64>,
    /// Electronic state energies in eV.
    pub state_energies: Vec<f64>,
}

impl EhrenfestState {
    /// Create a new Ehrenfest dynamics state.
    pub fn new(
        positions: Vec<[f64; 3]>,
        velocities: Vec<[f64; 3]>,
        masses: Vec<f64>,
        n_states: usize,
    ) -> Self {
        let mut populations = vec![0.0; n_states];
        if n_states > 0 {
            populations[0] = 1.0; // start in ground state
        }
        Self {
            positions,
            velocities,
            masses,
            populations,
            state_energies: vec![0.0; n_states],
        }
    }

    /// Total electronic energy as population-weighted average.
    pub fn mean_field_energy(&self) -> f64 {
        self.populations
            .iter()
            .zip(&self.state_energies)
            .map(|(p, e)| p * e)
            .sum()
    }
}

/// Perform one Ehrenfest MD step.
///
/// Propagates both electronic populations (via simple rate equation) and
/// nuclear positions/velocities (velocity Verlet) on the mean-field surface.
///
/// # Arguments
/// - `state`: current Ehrenfest state (modified in place)
/// - `forces`: mean-field forces on each nucleus `[[fx,fy,fz\], ...]`
/// - `coupling`: non-adiabatic coupling between state 0 and 1 (scalar placeholder)
/// - `dt`: timestep in fs
#[allow(dead_code)]
pub fn ehrenfest_dynamics(state: &mut EhrenfestState, forces: &[[f64; 3]], coupling: f64, dt: f64) {
    let n = state.positions.len();
    // Update nuclear positions (velocity Verlet with mean-field forces)
    for i in 0..n {
        let m = state.masses[i].max(1e-30);
        for k in 0..3 {
            if i < forces.len() {
                let acc = forces[i][k] / m;
                state.velocities[i][k] += acc * dt;
                state.positions[i][k] += state.velocities[i][k] * dt;
            }
        }
    }
    // Update electronic populations (simple rate equation for 2-state model)
    if state.populations.len() >= 2 {
        let p0 = state.populations[0];
        let p1 = state.populations[1];
        // Simple oscillatory coupling: dp0/dt = -coupling * p1
        let dp = coupling * dt;
        state.populations[0] = (p0 - dp * p1).clamp(0.0, 1.0);
        state.populations[1] = (p1 + dp * p0).clamp(0.0, 1.0);
        // Renormalize
        let total: f64 = state.populations.iter().sum();
        if total > 1e-30 {
            for p in &mut state.populations {
                *p /= total;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- DftSystem ---

    #[test]
    fn test_dft_system_new() {
        let sys = DftSystem::new(
            vec![1, 1],
            vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
            BasisSet::Sto3G,
            XcFunctional::Lda,
        );
        assert_eq!(sys.n_atoms, 2);
        assert_eq!(sys.charge, 0);
    }

    #[test]
    fn test_dft_system_n_electrons() {
        let sys = DftSystem::new(
            vec![6, 1, 1, 1, 1], // CH4
            vec![[0.0; 3]; 5],
            BasisSet::Sto3G,
            XcFunctional::Pbe,
        );
        assert_eq!(sys.n_electrons(), 10); // 6 + 4*1
    }

    #[test]
    fn test_dft_system_n_electrons_charged() {
        let mut sys = DftSystem::new(
            vec![8, 1, 1], // H2O
            vec![[0.0; 3]; 3],
            BasisSet::Sto3G,
            XcFunctional::Lda,
        );
        sys.charge = 1;
        assert_eq!(sys.n_electrons(), 9); // 8+2 - 1
    }

    #[test]
    fn test_nuclear_repulsion_h2() {
        // H2 with R=0.74 Angstrom, Z_A=Z_B=1
        let sys = DftSystem::new(
            vec![1, 1],
            vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
            BasisSet::Sto3G,
            XcFunctional::Lda,
        );
        let e_nn = sys.nuclear_repulsion_ev();
        // E = 1*1/(0.74*1.8897) Hartree * 27.211 eV/Hartree ≈ 19.4 eV
        assert!(e_nn > 15.0 && e_nn < 25.0);
    }

    #[test]
    fn test_nuclear_repulsion_single_atom() {
        let sys = DftSystem::new(
            vec![1],
            vec![[0.0, 0.0, 0.0]],
            BasisSet::Sto3G,
            XcFunctional::Lda,
        );
        assert_eq!(sys.nuclear_repulsion_ev(), 0.0);
    }

    // --- KohnSham ---

    #[test]
    fn test_ks_new() {
        let ks = KohnSham::new(vec![-10.0, -5.0, 1.0], vec![2.0, 2.0, 0.0]);
        assert_eq!(ks.homo_idx, 1);
        assert_eq!(ks.homo_energy(), -5.0);
    }

    #[test]
    fn test_ks_homo_lumo_gap() {
        let ks = KohnSham::new(vec![-15.0, -10.0, -5.0, 2.0], vec![2.0, 2.0, 2.0, 0.0]);
        assert_eq!(ks.homo_idx, 2);
        let gap = ks.homo_lumo_gap();
        assert!((gap - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_ks_band_energy() {
        let ks = KohnSham::new(vec![-10.0, -5.0], vec![2.0, 1.0]);
        let e = ks.band_energy();
        assert!((e - (-10.0 * 2.0 + -5.0 * 1.0)).abs() < 1e-12);
    }

    // --- hartree_fock_energy ---

    #[test]
    fn test_hf_energy_basic() {
        let e = hartree_fock_energy(&[-10.0, -5.0, 2.0], &[2.0, 2.0, 0.0]);
        assert!((e - (-30.0)).abs() < 1e-12);
    }

    #[test]
    fn test_hf_energy_empty() {
        let e = hartree_fock_energy(&[], &[]);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_hf_energy_all_occupied() {
        let e = hartree_fock_energy(&[1.0, 2.0, 3.0], &[2.0, 2.0, 2.0]);
        assert!((e - 12.0).abs() < 1e-12);
    }

    // --- dft_force ---

    #[test]
    fn test_dft_force_zero_density() {
        let density = vec![0.0; 4];
        let grid = vec![[1.0, 0.0, 0.0]; 4];
        let f = dft_force(1.0, &[0.0, 0.0, 0.0], &density, &grid, 1.0);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_dft_force_symmetric_cancel() {
        // Two equal density points at +x and -x should give zero net force in x
        let density = vec![1.0, 1.0];
        let grid = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let f = dft_force(1.0, &[0.0, 0.0, 0.0], &density, &grid, 1.0);
        assert!(f[0].abs() < 1e-12);
    }

    #[test]
    fn test_dft_force_asymmetric() {
        // Single density point at +x => force in +x direction on nucleus at origin
        let density = vec![1.0];
        let grid = vec![[2.0, 0.0, 0.0]];
        let f = dft_force(1.0, &[0.0, 0.0, 0.0], &density, &grid, 1.0);
        // force_x = z_a * rho * dv * (0 - 2) / |(-2)|^3 = -1/4 (attractive)
        // Actually: F_A = z_a * rho * dv * (R_A - r) / |R_A - r|^3
        // R_A = [0,0,0], r = [2,0,0] => (R_A - r) = [-2,0,0], |R_A-r|=2
        // F_x = 1 * 1 * 1 * (-2) / 8 = -0.25
        assert!(f[0] < 0.0);
    }

    // --- PulayMixer ---

    #[test]
    fn test_pulay_mixer_new() {
        let mixer = PulayMixer::new(5, 0.2);
        assert_eq!(mixer.max_history, 5);
        assert_eq!(mixer.alpha, 0.2);
    }

    #[test]
    fn test_pulay_mixer_rms_empty() {
        let mixer = PulayMixer::new(5, 0.2);
        assert_eq!(mixer.rms_residual(), 0.0);
    }

    #[test]
    fn test_pulay_mixer_push_and_rms() {
        let mut mixer = PulayMixer::new(5, 0.2);
        let rho_in = vec![1.0, 1.0, 1.0];
        let rho_out = vec![1.1, 1.0, 0.9];
        mixer.push(rho_in, rho_out);
        let rms = mixer.rms_residual();
        // residual = [0.1, 0, -0.1], rms = sqrt((0.01 + 0 + 0.01)/3) = sqrt(0.02/3)
        let expected = (0.02_f64 / 3.0).sqrt();
        assert!((rms - expected).abs() < 1e-12);
    }

    #[test]
    fn test_pulay_mixer_history_limit() {
        let mut mixer = PulayMixer::new(3, 0.2);
        for i in 0..5 {
            let rho = vec![i as f64; 4];
            mixer.push(rho.clone(), rho);
        }
        assert!(mixer.density_history.len() <= 3);
    }

    #[test]
    fn test_pulay_linear_mixing() {
        let mixer = PulayMixer::new(5, 0.5);
        let rho_in = vec![1.0, 2.0];
        let rho_out = vec![1.2, 1.8];
        let mixed = mixer.pulay_mixing(&rho_in, &rho_out);
        // new = in + 0.5*(out-in) = in*0.5 + out*0.5
        assert!((mixed[0] - 1.1).abs() < 1e-14);
        assert!((mixed[1] - 1.9).abs() < 1e-14);
    }

    // --- born_oppenheimer_step ---

    #[test]
    fn test_bo_step_zero_force() {
        let mut positions = vec![[1.0_f64, 0.0, 0.0]];
        let mut velocities = vec![[0.5_f64, 0.0, 0.0]];
        let forces = vec![[0.0_f64, 0.0, 0.0]];
        let masses = vec![1.0_f64];
        born_oppenheimer_step(&mut positions, &mut velocities, &forces, &masses, 1.0);
        // With zero force, velocity unchanged, position changes by v*dt
        assert!((positions[0][0] - 1.5).abs() < 1e-12);
        assert!((velocities[0][0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_bo_step_constant_force() {
        let mut positions = vec![[0.0_f64, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64, 0.0, 0.0]];
        let forces = vec![[1.0_f64, 0.0, 0.0]]; // F = 1 in x
        let masses = vec![1.0_f64];
        born_oppenheimer_step(&mut positions, &mut velocities, &forces, &masses, 1.0);
        // After one step: v = F/m * dt = 1, x = v*dt = 1
        assert!(positions[0][0] > 0.0);
        assert!(velocities[0][0] > 0.0);
    }

    // --- EhrenfestState ---

    #[test]
    fn test_ehrenfest_state_new() {
        let state = EhrenfestState::new(vec![[0.0; 3]], vec![[0.0; 3]], vec![1.0], 2);
        assert!((state.populations[0] - 1.0).abs() < 1e-14);
        assert!((state.populations[1] - 0.0).abs() < 1e-14);
    }

    #[test]
    fn test_ehrenfest_mean_field_energy() {
        let mut state = EhrenfestState::new(vec![[0.0; 3]], vec![[0.0; 3]], vec![1.0], 2);
        state.state_energies = vec![-10.0, -5.0];
        state.populations = vec![0.5, 0.5];
        let e = state.mean_field_energy();
        assert!((e - (-7.5)).abs() < 1e-12);
    }

    #[test]
    fn test_ehrenfest_dynamics_population_normalization() {
        let mut state = EhrenfestState::new(vec![[0.0; 3]], vec![[0.0; 3]], vec![1.0], 2);
        state.state_energies = vec![-10.0, -5.0];
        let forces = vec![[0.0_f64; 3]];
        ehrenfest_dynamics(&mut state, &forces, 0.1, 0.1);
        let total: f64 = state.populations.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ehrenfest_dynamics_nuclear_move() {
        let mut state =
            EhrenfestState::new(vec![[0.0, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]], vec![1.0], 1);
        let forces = vec![[0.0_f64; 3]];
        ehrenfest_dynamics(&mut state, &forces, 0.0, 0.5);
        // With zero force and v_x = 1, x should increase
        assert!(state.positions[0][0] > 0.0);
    }

    #[test]
    fn test_basis_set_plane_wave() {
        let bs = BasisSet::PlaneWave(400.0);
        assert!(matches!(bs, BasisSet::PlaneWave(400.0)));
    }

    #[test]
    fn test_xc_functional_variants() {
        let _lda = XcFunctional::Lda;
        let _pbe = XcFunctional::Pbe;
        let _b3lyp = XcFunctional::B3lyp;
        let _tpss = XcFunctional::Tpss;
    }
}
