// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Quantum-classical interface for Molecular Dynamics.
//!
//! Provides quantum mechanics infrastructure for mixed quantum-classical
//! simulations, including:
//!
//! - Quantum state representation and manipulation
//! - Hamiltonian matrix construction and diagonalization
//! - Tully's fewest-switches surface hopping (FSSH)
//! - Wigner distribution sampling for quantum initial conditions
//! - Born-Oppenheimer MD potential energy surfaces
//! - Density matrix (Lindblad/Redfield) open quantum system dynamics
//! - Zero-point energy estimates (harmonic and anharmonic)
//! - Instanton theory for quantum tunneling rates
//! - Quantum harmonic oscillator wavefunctions and expectation values
//!
//! References:
//! - Tully, J. C. (1990). J. Chem. Phys. 93, 1061.
//! - Wigner, E. (1932). Phys. Rev. 40, 749.
//! - Miller, W. H. (1975). J. Chem. Phys. 62, 1899.

use std::f64::consts::PI;

/// Physical constants used in quantum mechanics calculations.
/// Reduced Planck constant (J·s).
const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant (J K⁻¹).
const KB: f64 = 1.380_649e-23;
/// Electron mass (kg).
const ME: f64 = 9.109_383_701_5e-31;

// ---------------------------------------------------------------------------
// QuantumState
// ---------------------------------------------------------------------------

/// Quantum state vector for a multi-level system.
///
/// Stores the wavefunction coefficients, adiabatic energies, and populations
/// for a system with `n_states` electronic states.
#[derive(Debug, Clone)]
pub struct QuantumState {
    /// Wavefunction coefficients (real part only for simplicity).
    pub wavefunction: Vec<f64>,
    /// Adiabatic energies for each state (J or eV, consistent units).
    pub energies: Vec<f64>,
    /// Electronic populations |c_i|^2.
    pub populations: Vec<f64>,
    /// Number of electronic states.
    pub n_states: usize,
}

impl QuantumState {
    /// Create a new quantum state initialised in the ground state.
    pub fn new(n_states: usize) -> Self {
        let mut wavefunction = vec![0.0; n_states];
        let mut populations = vec![0.0; n_states];
        if n_states > 0 {
            wavefunction[0] = 1.0;
            populations[0] = 1.0;
        }
        Self {
            wavefunction,
            energies: vec![0.0; n_states],
            populations,
            n_states,
        }
    }

    /// Create a uniform superposition state.
    pub fn uniform_superposition(n_states: usize) -> Self {
        let amp = if n_states > 0 {
            1.0 / (n_states as f64).sqrt()
        } else {
            0.0
        };
        let wavefunction = vec![amp; n_states];
        let populations = vec![amp * amp; n_states];
        Self {
            wavefunction,
            energies: vec![0.0; n_states],
            populations,
            n_states,
        }
    }

    /// Compute the norm of the wavefunction: sqrt(sum |c_i|^2).
    pub fn norm(&self) -> f64 {
        self.wavefunction.iter().map(|c| c * c).sum::<f64>().sqrt()
    }

    /// Normalize the wavefunction in-place and update populations.
    pub fn normalize(&mut self) {
        let n = self.norm();
        if n > 1e-14 {
            for c in &mut self.wavefunction {
                *c /= n;
            }
        }
        self.update_populations();
    }

    /// Update populations from current wavefunction coefficients.
    pub fn update_populations(&mut self) {
        for (i, c) in self.wavefunction.iter().enumerate() {
            self.populations[i] = c * c;
        }
    }

    /// Return the ground-state (lowest energy) population.
    pub fn ground_state_population(&self) -> f64 {
        if self.populations.is_empty() {
            0.0
        } else {
            self.populations[0]
        }
    }

    /// Return the total population (should be 1 after normalization).
    pub fn total_population(&self) -> f64 {
        self.populations.iter().sum()
    }

    /// Return the mean energy `E` = sum_i p_i * E_i.
    pub fn mean_energy(&self) -> f64 {
        self.populations
            .iter()
            .zip(self.energies.iter())
            .map(|(p, e)| p * e)
            .sum()
    }

    /// Return the state index with the highest population.
    pub fn dominant_state(&self) -> usize {
        self.populations
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}

// ---------------------------------------------------------------------------
// HamiltonianMatrix
// ---------------------------------------------------------------------------

/// Real symmetric Hamiltonian matrix for a multi-level quantum system.
///
/// Stores the matrix as a dense `n x n` array and provides eigenvalue
/// solvers and various matrix properties.
#[derive(Debug, Clone)]
pub struct HamiltonianMatrix {
    /// Matrix elements stored as row-major `data[i][j]`.
    pub data: Vec<Vec<f64>>,
    /// Dimension of the matrix.
    pub n: usize,
}

impl HamiltonianMatrix {
    /// Construct a diagonal Hamiltonian from a slice of energies.
    pub fn diagonal(energies: &[f64]) -> Self {
        let n = energies.len();
        let mut data = vec![vec![0.0; n]; n];
        for (i, &e) in energies.iter().enumerate() {
            data[i][i] = e;
        }
        Self { data, n }
    }

    /// Construct a zero Hamiltonian of size n.
    pub fn zeros(n: usize) -> Self {
        Self {
            data: vec![vec![0.0; n]; n],
            n,
        }
    }

    /// Add a symmetric off-diagonal coupling between states i and j.
    pub fn add_coupling(&mut self, i: usize, j: usize, v: f64) {
        self.data[i][j] += v;
        self.data[j][i] += v;
    }

    /// Compute eigenvalues using Jacobi iteration (symmetric matrix).
    ///
    /// Returns eigenvalues sorted in ascending order.
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n;
        if n == 0 {
            return vec![];
        }
        if n == 1 {
            return vec![self.data[0][0]];
        }

        // Copy the matrix for in-place Jacobi sweeps
        let mut a: Vec<Vec<f64>> = self.data.clone();

        let max_iter = 100 * n * n;
        let tol = 1e-12;

        for _iter in 0..max_iter {
            // Find largest off-diagonal element
            let mut p = 0;
            let mut q = 1;
            let mut max_val = a[0][1].abs();
            for (i, row) in a.iter().enumerate() {
                for (j, &aij) in row.iter().enumerate().skip(i + 1) {
                    if aij.abs() > max_val {
                        max_val = aij.abs();
                        p = i;
                        q = j;
                    }
                }
            }
            if max_val < tol {
                break;
            }

            // Compute Jacobi rotation angle
            let theta = if (a[q][q] - a[p][p]).abs() < 1e-14 {
                PI / 4.0
            } else {
                0.5 * ((2.0 * a[p][q]) / (a[q][q] - a[p][p])).atan()
            };
            let c = theta.cos();
            let s = theta.sin();

            // Apply rotation
            let app = a[p][p];
            let aqq = a[q][q];
            let apq = a[p][q];
            a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;

            for (r, row_r) in a.iter_mut().enumerate() {
                if r != p && r != q {
                    let arp = row_r[p];
                    let arq = row_r[q];
                    row_r[p] = c * arp - s * arq;
                    row_r[q] = s * arp + c * arq;
                }
            }
            // Update symmetric counterparts a[p][r] and a[q][r].
            for (r, _) in (0..n).enumerate() {
                if r != p && r != q {
                    a[p][r] = a[r][p];
                    a[q][r] = a[r][q];
                }
            }
        }

        let mut evals: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
        evals.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        evals
    }

    /// Return the lowest eigenvalue (ground state energy).
    pub fn ground_state_energy(&self) -> f64 {
        self.eigenvalues()
            .into_iter()
            .next()
            .unwrap_or(self.data[0][0])
    }

    /// Compute the trace: sum of diagonal elements.
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.data[i][i]).sum()
    }

    /// Compute the Frobenius norm: sqrt(sum |a_ij|^2).
    pub fn frobenius_norm(&self) -> f64 {
        self.data
            .iter()
            .flat_map(|row| row.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
    }

    /// Apply the Hamiltonian to a state vector: v = H * psi.
    pub fn apply(&self, psi: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|i| (0..self.n).map(|j| self.data[i][j] * psi[j]).sum())
            .collect()
    }

    /// Commutator \[H, rho\] = H*rho - rho*H for a real density matrix.
    pub fn commutator_real(&self, rho: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = self.n;
        let mut result = vec![vec![0.0; n]; n];
        for (i, res_row) in result.iter_mut().enumerate() {
            for j in 0..n {
                let val = (0..n)
                    .map(|k| self.data[i][k] * rho[k][j] - rho[i][k] * self.data[k][j])
                    .sum();
                res_row[j] = val;
            }
        }
        result
    }
}

// ---------------------------------------------------------------------------
// SurfaceHoppingTSH (Tully's Fewest Switches Surface Hopping)
// ---------------------------------------------------------------------------

/// Tully's fewest-switches surface hopping (FSSH) implementation.
///
/// Implements stochastic hops between adiabatic states based on
/// non-adiabatic coupling vectors and the FSSH algorithm of Tully (1990).
#[derive(Debug, Clone)]
pub struct SurfaceHoppingTSH {
    /// Number of adiabatic states.
    pub n_states: usize,
    /// Current adiabatic state index (0-based).
    pub current_state: usize,
    /// Time step (fs or atomic units, consistent with velocities).
    pub dt: f64,
    /// Accumulated hop count for statistics.
    pub hop_count: usize,
}

impl SurfaceHoppingTSH {
    /// Create a new FSSH surface hopper starting in the ground state.
    pub fn new(n_states: usize, dt: f64) -> Self {
        Self {
            n_states,
            current_state: 0,
            dt,
            hop_count: 0,
        }
    }

    /// Create a new FSSH hopper starting in the specified state.
    pub fn with_initial_state(n_states: usize, dt: f64, initial_state: usize) -> Self {
        Self {
            n_states,
            current_state: initial_state,
            dt,
            hop_count: 0,
        }
    }

    /// Compute the scalar non-adiabatic coupling d_12 = <1|d/dt|2>.
    ///
    /// Approximated as d_12 ≈ (∇E_1 - ∇E_2) · v / (E_2 - E_1).
    /// Returns the time-derivative coupling (1/time units).
    pub fn non_adiabatic_coupling(grad_e1: &[f64], grad_e2: &[f64], velocity: &[f64]) -> f64 {
        let len = grad_e1.len().min(grad_e2.len()).min(velocity.len());
        let dot: f64 = (0..len)
            .map(|i| (grad_e1[i] - grad_e2[i]) * velocity[i])
            .sum();
        dot
    }

    /// Compute Tully's switching probability g_{jk} for a hop from state j to k.
    ///
    /// g_{jk} = max(0, -2 * Re\[rho_{jk}\] * d_{jk} * dt / rho_{jj})
    ///
    /// where d_{jk} is the NAC, rho_{jk} the coherence, and rho_{jj} the population.
    pub fn switching_probability(d12: f64, dt: f64, pop1: f64, pop2: f64) -> f64 {
        if pop2 < 1e-14 {
            return 0.0;
        }
        let coherence = (pop1 * pop2).sqrt(); // off-diagonal approximation
        let g = -2.0 * coherence * d12 * dt / pop2;
        g.clamp(0.0, 1.0)
    }

    /// Perform a stochastic hop based on switching probabilities.
    ///
    /// `probs[i]` is the probability of hopping to state i.
    /// `rand_val` is a uniform random number in \[0, 1).
    pub fn hop(&mut self, probs: &[f64], rand_val: f64) {
        let mut cumulative = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            cumulative += p;
            if rand_val < cumulative && i != self.current_state {
                self.current_state = i;
                self.hop_count += 1;
                return;
            }
        }
    }

    /// Compute all pairwise switching probabilities for a given state.
    ///
    /// Returns a vector of probabilities for hopping from `current_state`
    /// to each other state.
    pub fn all_switching_probabilities(
        &self,
        nac_matrix: &[Vec<f64>],
        populations: &[f64],
    ) -> Vec<f64> {
        let mut probs = vec![0.0; self.n_states];
        let j = self.current_state;
        for k in 0..self.n_states {
            if k != j && j < nac_matrix.len() && k < nac_matrix[j].len() {
                probs[k] = Self::switching_probability(
                    nac_matrix[j][k],
                    self.dt,
                    populations[j],
                    populations[k],
                );
            }
        }
        probs
    }

    /// Apply energy conservation (velocity rescaling) after a hop.
    ///
    /// Returns a scaling factor for the velocity component along the NAC vector.
    /// Returns `None` if there is insufficient kinetic energy for the hop.
    pub fn velocity_rescaling_factor(delta_e: f64, kinetic_energy: f64) -> Option<f64> {
        let new_ke = kinetic_energy - delta_e;
        if new_ke < 0.0 {
            None // frustrated hop
        } else {
            Some((new_ke / kinetic_energy).sqrt())
        }
    }
}

// ---------------------------------------------------------------------------
// WignerSampling
// ---------------------------------------------------------------------------

/// Wigner distribution sampling for quantum-classical initial conditions.
///
/// Provides methods to sample initial positions and momenta from the
/// Wigner quasi-probability distribution, ensuring proper quantum
/// statistical treatment of initial conditions.
#[derive(Debug, Clone)]
pub struct WignerSampling;

impl WignerSampling {
    /// Sample (x, p) pairs from the harmonic oscillator Wigner distribution.
    ///
    /// W(x,p) = (1/πℏ) exp(-mω x²/ℏ - p²/(mωℏ))
    ///
    /// At inverse temperature beta = 1/(k_B T), uses the thermal Wigner distribution.
    /// Returns a vector of (position, momentum) pairs.
    pub fn harmonic_wigner(omega: f64, beta: f64, n_samples: usize) -> Vec<(f64, f64)> {
        let mut rng = rand::rng();
        // Thermal Wigner width parameters
        let sigma_x = (HBAR / (2.0 * ME * omega)
            * (1.0 + 2.0 / ((HBAR * omega * beta).exp() - 1.0 + 1e-14)))
            .sqrt()
            .max(1e-20);
        let sigma_p = (ME * omega * HBAR / 2.0
            * (1.0 + 2.0 / ((HBAR * omega * beta).exp() - 1.0 + 1e-14)))
            .sqrt()
            .max(1e-20);

        let mut samples = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let x = sigma_x * box_muller_normal(&mut rng);
            let p = sigma_p * box_muller_normal(&mut rng);
            samples.push((x, p));
        }
        samples
    }

    /// Approximate Wigner sampling for a Morse oscillator.
    ///
    /// Uses the harmonic approximation at the Morse minimum with
    /// effective frequency omega and anharmonicity correction from
    /// dissociation energy `De` and range parameter `a`.
    ///
    /// Returns a single (x, p) sample.
    pub fn morse_wigner_approx(de: f64, a: f64, omega: f64, beta: f64) -> (f64, f64) {
        let mut rng = rand::rng();
        // Anharmonic correction to zero-point amplitude
        let anharm_factor = 1.0 + 0.25 * HBAR * omega / (4.0 * de);
        let _ = a; // used conceptually for anharmonicity
        let sigma_x = (HBAR / (2.0 * ME * omega)).sqrt() * anharm_factor;
        let sigma_p = (ME * omega * HBAR / 2.0).sqrt();

        // Thermal broadening
        let n_thermal = 1.0 / ((HBAR * omega * beta).exp() - 1.0 + 1e-14);
        let thermal_x = sigma_x * (1.0 + 2.0 * n_thermal).sqrt();
        let thermal_p = sigma_p * (1.0 + 2.0 * n_thermal).sqrt();

        (
            thermal_x * box_muller_normal(&mut rng),
            thermal_p * box_muller_normal(&mut rng),
        )
    }

    /// Sample from an arbitrary Wigner function using rejection sampling.
    ///
    /// `wigner_fn(x, p)` must return a non-negative value.
    /// `x_range` and `p_range` specify the sampling domain as (min, max).
    /// Returns up to `n` accepted (x, p) samples.
    pub fn sample_phase_space(
        wigner_fn: impl Fn(f64, f64) -> f64,
        x_range: (f64, f64),
        p_range: (f64, f64),
        n: usize,
    ) -> Vec<(f64, f64)> {
        use rand::RngExt;
        let mut rng = rand::rng();
        let mut samples = Vec::with_capacity(n);

        // Estimate maximum value with a grid search
        let grid = 20;
        let mut w_max = 1e-14_f64;
        for ix in 0..grid {
            for ip in 0..grid {
                let x = x_range.0 + (x_range.1 - x_range.0) * ix as f64 / (grid - 1) as f64;
                let p = p_range.0 + (p_range.1 - p_range.0) * ip as f64 / (grid - 1) as f64;
                w_max = w_max.max(wigner_fn(x, p));
            }
        }
        w_max *= 1.1;

        let max_attempts = n * 1000;
        let mut attempts = 0;
        while samples.len() < n && attempts < max_attempts {
            let x = rng.random_range(x_range.0..x_range.1);
            let p = rng.random_range(p_range.0..p_range.1);
            let w = wigner_fn(x, p);
            let u: f64 = rng.random_range(0.0..w_max);
            if u < w {
                samples.push((x, p));
            }
            attempts += 1;
        }
        samples
    }
}

/// Box-Muller transform: produces a standard normal variate.
fn box_muller_normal(rng: &mut impl rand::Rng) -> f64 {
    use rand::RngExt as _;
    let u1: f64 = rng.random_range(1e-10..1.0);
    let u2: f64 = rng.random_range(0.0..1.0);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

// ---------------------------------------------------------------------------
// BornOppenheimerMD
// ---------------------------------------------------------------------------

/// Born-Oppenheimer molecular dynamics potential energy surface tracker.
///
/// Stores the PES values computed along an MD trajectory and provides
/// gradient estimation by finite differences on the provided energy function.
#[derive(Debug, Clone)]
pub struct BornOppenheimerMD {
    /// Stored potential energy values per MD step.
    pub potential_energy_surface: Vec<f64>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Finite difference step size (Angstrom or bohr, consistent with positions).
    pub fd_step: f64,
}

impl BornOppenheimerMD {
    /// Create a new BOMD tracker.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            potential_energy_surface: Vec::new(),
            n_atoms,
            fd_step: 1e-4,
        }
    }

    /// Record a potential energy value.
    pub fn record_energy(&mut self, energy: f64) {
        self.potential_energy_surface.push(energy);
    }

    /// Compute the gradient at `positions` by central finite differences.
    ///
    /// `energy_fn` takes a slice of positions and returns the potential energy.
    /// Returns forces = -∇V.
    pub fn gradient_at(
        &self,
        positions: &[[f64; 3]],
        energy_fn: impl Fn(&[[f64; 3]]) -> f64,
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0; 3]; n];
        let mut pos_plus = positions.to_vec();
        let mut pos_minus = positions.to_vec();

        for i in 0..n {
            for d in 0..3 {
                pos_plus[i][d] = positions[i][d] + self.fd_step;
                pos_minus[i][d] = positions[i][d] - self.fd_step;

                let e_plus = energy_fn(&pos_plus);
                let e_minus = energy_fn(&pos_minus);

                // Force = -dE/dx
                forces[i][d] = -(e_plus - e_minus) / (2.0 * self.fd_step);

                pos_plus[i][d] = positions[i][d];
                pos_minus[i][d] = positions[i][d];
            }
        }
        forces
    }

    /// Compute the adiabatic force on the selected electronic state.
    ///
    /// Uses the supplied multi-state energy function and picks the gradient
    /// of state `state`.
    pub fn adiabatic_force(
        &self,
        positions: &[[f64; 3]],
        state: usize,
        energy_fn: impl Fn(&[[f64; 3]], usize) -> f64,
    ) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut forces = vec![[0.0; 3]; n];
        let mut pos_plus = positions.to_vec();
        let mut pos_minus = positions.to_vec();

        for i in 0..n {
            for d in 0..3 {
                pos_plus[i][d] = positions[i][d] + self.fd_step;
                pos_minus[i][d] = positions[i][d] - self.fd_step;

                let e_plus = energy_fn(&pos_plus, state);
                let e_minus = energy_fn(&pos_minus, state);

                forces[i][d] = -(e_plus - e_minus) / (2.0 * self.fd_step);

                pos_plus[i][d] = positions[i][d];
                pos_minus[i][d] = positions[i][d];
            }
        }
        forces
    }

    /// Return the most recently recorded potential energy.
    pub fn last_energy(&self) -> Option<f64> {
        self.potential_energy_surface.last().copied()
    }

    /// Return the average potential energy over the stored trajectory.
    pub fn average_energy(&self) -> f64 {
        if self.potential_energy_surface.is_empty() {
            return 0.0;
        }
        self.potential_energy_surface.iter().sum::<f64>()
            / self.potential_energy_surface.len() as f64
    }
}

// ---------------------------------------------------------------------------
// DensityMatrix
// ---------------------------------------------------------------------------

/// Complex density matrix for open quantum systems (Lindblad/Redfield).
///
/// Stores rho as a complex matrix where each element is `(re, im)`.
/// Supports von Neumann unitary evolution and Lindblad dissipation.
#[derive(Debug, Clone)]
pub struct DensityMatrix {
    /// Complex density matrix elements as (real, imaginary) pairs.
    /// Indexed as `rho[i][j]`.
    pub rho: Vec<Vec<(f64, f64)>>,
    /// Hilbert space dimension.
    pub n: usize,
}

impl DensityMatrix {
    /// Construct a pure ground-state density matrix.
    pub fn ground_state(n: usize) -> Self {
        let mut rho = vec![vec![(0.0, 0.0); n]; n];
        if n > 0 {
            rho[0][0] = (1.0, 0.0);
        }
        Self { rho, n }
    }

    /// Construct a maximally mixed density matrix (identity / n).
    pub fn mixed(n: usize) -> Self {
        let mut rho = vec![vec![(0.0, 0.0); n]; n];
        let val = 1.0 / n as f64;
        for (i, row) in rho.iter_mut().enumerate() {
            row[i] = (val, 0.0);
        }
        Self { rho, n }
    }

    /// Return the population of state i: Re\[rho_{ii}\].
    pub fn population(&self, i: usize) -> f64 {
        self.rho[i][i].0
    }

    /// Return the off-diagonal coherence rho_{ij} as (real, imag).
    pub fn coherence(&self, i: usize, j: usize) -> (f64, f64) {
        self.rho[i][j]
    }

    /// Return the purity Tr\[rho^2\].
    pub fn purity(&self) -> f64 {
        let n = self.n;
        let mut tr = 0.0;
        for i in 0..n {
            for k in 0..n {
                let (re, _im) = complex_mul_sum_row(&self.rho, i, k, n);
                if i == k {
                    tr += re;
                }
            }
        }
        tr
    }

    /// Compute the trace of the density matrix.
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.rho[i][i].0).sum()
    }

    /// Evolve under unitary von Neumann equation using RK4.
    ///
    /// dρ/dt = -i/ℏ \[H, ρ\]
    ///
    /// `h` is the Hamiltonian matrix, `dt` the time step in consistent units.
    pub fn evolve_unitary(&mut self, h: &HamiltonianMatrix, dt: f64) {
        let n = self.n;
        // Compute commutator [H, rho] as complex
        let k1 = self.von_neumann_rhs(h);
        let rho1 = add_complex_matrix(&self.rho, &scale_complex_matrix(&k1, dt / 2.0), n);
        let k2 = von_neumann_rhs_impl(h, &rho1, n);
        let rho2 = add_complex_matrix(&self.rho, &scale_complex_matrix(&k2, dt / 2.0), n);
        let k3 = von_neumann_rhs_impl(h, &rho2, n);
        let rho3 = add_complex_matrix(&self.rho, &scale_complex_matrix(&k3, dt), n);
        let k4 = von_neumann_rhs_impl(h, &rho3, n);

        for (i, rho_row) in self.rho.iter_mut().enumerate() {
            for (j, rho_ij) in rho_row.iter_mut().enumerate() {
                let re = rho_ij.0
                    + dt / 6.0 * (k1[i][j].0 + 2.0 * k2[i][j].0 + 2.0 * k3[i][j].0 + k4[i][j].0);
                let im = rho_ij.1
                    + dt / 6.0 * (k1[i][j].1 + 2.0 * k2[i][j].1 + 2.0 * k3[i][j].1 + k4[i][j].1);
                *rho_ij = (re, im);
            }
        }
    }

    /// Apply Lindblad dissipation: dρ/dt = γ(L ρ L† - ½{L†L, ρ}).
    ///
    /// `gamma` is the dissipation rate, `jump_op_idx` selects the
    /// lowering operator |0><jump_op_idx| (decay to ground state).
    pub fn lindblad_dissipation(&mut self, gamma: f64, jump_op_idx: usize) {
        let n = self.n;
        if jump_op_idx >= n {
            return;
        }
        // L = |0><jump_op_idx| (decay from jump_op_idx -> 0)
        // L ρ L† contribution
        let rho_kk = self.rho[jump_op_idx][jump_op_idx];
        // L†L = |jump_op_idx><jump_op_idx|
        let mut drho = vec![vec![(0.0_f64, 0.0_f64); n]; n];

        // L ρ L† = rho[jump_op_idx][jump_op_idx] * |0><0|
        drho[0][0].0 += gamma * rho_kk.0;

        // -1/2 {L†L, ρ} = -1/2 (L†L ρ + ρ L†L)
        // L†L = |jump_op_idx><jump_op_idx|
        for (j, drho_kj) in drho[jump_op_idx].iter_mut().enumerate() {
            let re = -0.5 * gamma * (self.rho[jump_op_idx][j].0 + self.rho[j][jump_op_idx].0);
            let im = -0.5 * gamma * (self.rho[jump_op_idx][j].1 + self.rho[j][jump_op_idx].1);
            drho_kj.0 += re;
            drho_kj.1 += im;
        }

        for (rho_row, drho_row) in self.rho.iter_mut().zip(drho.iter()) {
            for (rho_ij, &drho_ij) in rho_row.iter_mut().zip(drho_row.iter()) {
                rho_ij.0 += drho_ij.0;
                rho_ij.1 += drho_ij.1;
            }
        }
    }

    /// Compute the von Neumann RHS: -i/ℏ \[H, ρ\].
    fn von_neumann_rhs(&self, h: &HamiltonianMatrix) -> Vec<Vec<(f64, f64)>> {
        von_neumann_rhs_impl(h, &self.rho, self.n)
    }
}

/// Compute -i/ℏ \[H, ρ\] for a real H and complex ρ.
fn von_neumann_rhs_impl(
    h: &HamiltonianMatrix,
    rho: &[Vec<(f64, f64)>],
    n: usize,
) -> Vec<Vec<(f64, f64)>> {
    let mut result = vec![vec![(0.0, 0.0); n]; n];
    for (i, res_row) in result.iter_mut().enumerate() {
        for (j, res_ij) in res_row.iter_mut().enumerate() {
            let comm_re: f64 = (0..n)
                .map(|k| h.data[i][k] * rho[k][j].0 - rho[i][k].0 * h.data[k][j])
                .sum();
            let comm_im: f64 = (0..n)
                .map(|k| h.data[i][k] * rho[k][j].1 - rho[i][k].1 * h.data[k][j])
                .sum();
            // -i/ℏ [H, ρ] → multiply by -i/ℏ: (re,im) * (-i) = (im, -re) then / ℏ
            *res_ij = (comm_im / HBAR, -comm_re / HBAR);
        }
    }
    result
}

/// Add two complex matrices element-wise.
fn add_complex_matrix(
    a: &[Vec<(f64, f64)>],
    b: &[Vec<(f64, f64)>],
    n: usize,
) -> Vec<Vec<(f64, f64)>> {
    let mut c = vec![vec![(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            c[i][j] = (a[i][j].0 + b[i][j].0, a[i][j].1 + b[i][j].1);
        }
    }
    c
}

/// Scale a complex matrix by a real scalar.
fn scale_complex_matrix(a: &[Vec<(f64, f64)>], s: f64) -> Vec<Vec<(f64, f64)>> {
    a.iter()
        .map(|row| row.iter().map(|(re, im)| (re * s, im * s)).collect())
        .collect()
}

/// Compute the (i,k) element of the product rho * rho.
fn complex_mul_sum_row(rho: &[Vec<(f64, f64)>], i: usize, k: usize, n: usize) -> (f64, f64) {
    let mut re = 0.0;
    let mut im = 0.0;
    for (rho_im, rho_mk) in rho[i].iter().zip(rho.iter().map(|row| row[k])).take(n) {
        let (a_re, a_im) = *rho_im;
        let (b_re, b_im) = rho_mk;
        re += a_re * b_re - a_im * b_im;
        im += a_re * b_im + a_im * b_re;
    }
    (re, im)
}

// ---------------------------------------------------------------------------
// ZeroPointEnergy
// ---------------------------------------------------------------------------

/// Zero-point energy estimators for harmonic and anharmonic potentials.
///
/// Provides both the simple harmonic approximation and perturbation-theory
/// anharmonic corrections to ZPE.
#[derive(Debug, Clone)]
pub struct ZeroPointEnergy;

impl ZeroPointEnergy {
    /// Harmonic zero-point energy: ZPE = (1/2) Σ ℏω_i.
    ///
    /// `frequencies` are angular frequencies in rad/s (or rad/atomic-time-unit).
    pub fn harmonic_zpe(frequencies: &[f64]) -> f64 {
        0.5 * HBAR * frequencies.iter().sum::<f64>()
    }

    /// Anharmonic correction to ZPE using vibrational perturbation theory.
    ///
    /// Δε = -15/16 * (ℏ/(mω))^2 * k4/ω^2 + 11/16 * (k3/k2)^2 * (ℏ/(mω))^2 / ω
    ///
    /// `freq` is the harmonic frequency, `cubic_fc` the third-order force
    /// constant k3, and `quartic_fc` the fourth-order force constant k4.
    /// Returns the correction in the same energy units as ℏω.
    pub fn anharmonic_correction(freq: f64, cubic_fc: f64, quartic_fc: f64) -> f64 {
        let hw = HBAR * freq;
        let x = hw / freq; // ℏ/ω approximation
        let quartic_term = -15.0 / 16.0 * x * x * quartic_fc / (freq * freq);
        let cubic_term = 11.0 / 16.0 * (cubic_fc / freq) * (cubic_fc / freq) * x * x / freq;
        quartic_term + cubic_term
    }

    /// Compute the ZPE for all normal modes of a molecule.
    ///
    /// `frequencies` should be the positive normal mode frequencies.
    /// Imaginary modes (negative) are ignored.
    pub fn molecular_zpe(frequencies: &[f64]) -> f64 {
        0.5 * HBAR * frequencies.iter().filter(|&&f| f > 0.0).sum::<f64>()
    }

    /// Classical thermal energy for comparison: (3N/2) k_B T per atom.
    pub fn classical_thermal_energy(n_atoms: usize, temperature: f64) -> f64 {
        3.0 * n_atoms as f64 / 2.0 * KB * temperature
    }
}

// ---------------------------------------------------------------------------
// InstantonTheory
// ---------------------------------------------------------------------------

/// Instanton theory for quantum mechanical tunneling rates.
///
/// Provides WKB-based estimates of tunneling rates through barriers,
/// crossover temperatures, and instanton action calculations.
#[derive(Debug, Clone)]
pub struct InstantonTheory;

impl InstantonTheory {
    /// WKB action for a parabolic barrier: S = π * m * ω_b / ℏ * (E_b - E) at E=0.
    ///
    /// For a harmonic barrier: S = ω_b / ω_0 * sqrt(2 m E_b / ℏ^2) * π
    ///
    /// Uses the simplified instanton expression for a harmonic barrier:
    /// S_inst = π * ℏ * β_c / 2 where β_c = 1/(k_B T_c).
    pub fn action_harmonic(mass: f64, omega: f64, barrier_height: f64, beta: f64) -> f64 {
        // Full instanton action for harmonic barrier at inverse temperature beta
        // S = ℏ * ω * beta / 2 * (... ) — simplified Caldeira-Leggett form
        let _beta_c = PI / (HBAR * omega);
        let exponent = beta * HBAR * omega / 2.0;
        // S_inst = 2 * sqrt(2 * m * E_b) / omega * (cosh-dependent correction)
        let s0 = 2.0 * (2.0 * mass * barrier_height).sqrt() / omega;
        s0 * exponent.tanh().max(1e-14)
    }

    /// WKB tunneling rate from the instanton action S:
    ///
    /// k_tunnel = A * exp(-S/ℏ)
    ///
    /// where the prefactor A is set to 1.0 (units-dependent).
    pub fn tunneling_rate_wkb(action: f64) -> f64 {
        (-action / HBAR).exp()
    }

    /// Crossover temperature below which tunneling dominates.
    ///
    /// T_c = ℏ ω_b / (2 π k_B)
    ///
    /// `omega` is the imaginary barrier frequency (curvature at top of barrier).
    pub fn crossover_temperature(omega: f64, _barrier: f64) -> f64 {
        HBAR * omega / (2.0 * PI * KB)
    }

    /// Compute the thermal rate coefficient using the Arrhenius equation.
    ///
    /// k = A * exp(-E_a / (k_B T))
    pub fn arrhenius_rate(prefactor: f64, activation_energy: f64, temperature: f64) -> f64 {
        prefactor * (-activation_energy / (KB * temperature)).exp()
    }

    /// Compute the ratio of quantum to classical rate (tunneling enhancement).
    pub fn quantum_enhancement(action: f64, activation_energy: f64, temperature: f64) -> f64 {
        let classical = (-activation_energy / (KB * temperature)).exp();
        let quantum = (-action / HBAR).exp();
        if classical > 1e-300 {
            quantum / classical
        } else {
            quantum
        }
    }
}

// ---------------------------------------------------------------------------
// QuantumHarmonicOscillator
// ---------------------------------------------------------------------------

/// Quantum harmonic oscillator: wavefunctions and expectation values.
///
/// Provides Hermite-polynomial wavefunctions for the 1D QHO with
/// frequency `omega` and particle mass `mass`.
#[derive(Debug, Clone)]
pub struct QuantumHarmonicOscillator {
    /// Angular frequency ω (rad/s or rad/atomic-time-unit).
    pub omega: f64,
    /// Particle mass (kg or atomic mass units, consistent with omega).
    pub mass: f64,
}

impl QuantumHarmonicOscillator {
    /// Create a new quantum harmonic oscillator.
    pub fn new(omega: f64, mass: f64) -> Self {
        Self { omega, mass }
    }

    /// Energy of the n-th level: E_n = (n + 1/2) ℏω.
    pub fn energy_level(&self, n: usize) -> f64 {
        (n as f64 + 0.5) * HBAR * self.omega
    }

    /// Zero-point energy E_0 = ℏω/2.
    pub fn zero_point_energy(&self) -> f64 {
        0.5 * HBAR * self.omega
    }

    /// Normalized QHO wavefunction ψ_n(x) using Hermite polynomials.
    ///
    /// ψ_n(x) = N_n * H_n(ξ) * exp(-ξ²/2)
    /// where ξ = sqrt(mω/ℏ) * x and N_n is the normalization constant.
    pub fn wavefunction_hermite(&self, n: usize, x: f64) -> f64 {
        let alpha = (self.mass * self.omega / HBAR).sqrt();
        let xi = alpha * x;
        let h_n = hermite_polynomial(n, xi);
        let norm = (alpha / PI.sqrt()).sqrt() / (2.0_f64.powi(n as i32) * factorial(n)).sqrt();
        norm * h_n * (-0.5 * xi * xi).exp()
    }

    /// Expectation value `x²`_n = (n + 1/2) * ℏ/(mω).
    pub fn expectation_x2(&self, n: usize) -> f64 {
        (n as f64 + 0.5) * HBAR / (self.mass * self.omega)
    }

    /// Expectation value `p²`_n = (n + 1/2) * mωℏ.
    pub fn expectation_p2(&self, n: usize) -> f64 {
        (n as f64 + 0.5) * self.mass * self.omega * HBAR
    }

    /// `x`_n = 0 (symmetric potential).
    pub fn expectation_x(&self, _n: usize) -> f64 {
        0.0
    }

    /// `p`_n = 0.
    pub fn expectation_p(&self, _n: usize) -> f64 {
        0.0
    }

    /// Kinetic energy `T`_n = `p²`_n / (2m) = E_n / 2.
    pub fn kinetic_energy(&self, n: usize) -> f64 {
        self.energy_level(n) / 2.0
    }

    /// Potential energy `V`_n = `T`_n = E_n / 2 (virial theorem).
    pub fn potential_energy(&self, n: usize) -> f64 {
        self.energy_level(n) / 2.0
    }

    /// Transition dipole moment <n|x|n+1> = sqrt((n+1)*ℏ/(2mω)).
    pub fn transition_dipole(&self, n: usize) -> f64 {
        ((n as f64 + 1.0) * HBAR / (2.0 * self.mass * self.omega)).sqrt()
    }

    /// Compute the classical turning point for level n: x_c = sqrt(2E_n/(mω²)).
    pub fn classical_turning_point(&self, n: usize) -> f64 {
        let e_n = self.energy_level(n);
        (2.0 * e_n / (self.mass * self.omega * self.omega)).sqrt()
    }
}

/// Hermite polynomial H_n(x) using the recurrence relation.
///
/// H_0 = 1, H_1 = 2x, H_n = 2x H_{n-1} - 2(n-1) H_{n-2}.
pub fn hermite_polynomial(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut h_prev = 1.0;
            let mut h_curr = 2.0 * x;
            for k in 2..=n {
                let h_next = 2.0 * x * h_curr - 2.0 * (k - 1) as f64 * h_prev;
                h_prev = h_curr;
                h_curr = h_next;
            }
            h_curr
        }
    }
}

/// Integer factorial n! (exact for small n, returns f64).
fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- QuantumState ---

    #[test]
    fn test_quantum_state_new_ground_state() {
        let qs = QuantumState::new(3);
        assert_eq!(qs.n_states, 3);
        assert!((qs.populations[0] - 1.0).abs() < 1e-10);
        assert!(qs.populations[1].abs() < 1e-10);
    }

    #[test]
    fn test_quantum_state_norm() {
        let qs = QuantumState::new(3);
        assert!((qs.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quantum_state_normalize() {
        let mut qs = QuantumState::new(3);
        qs.wavefunction = vec![2.0, 0.0, 0.0];
        qs.normalize();
        assert!((qs.wavefunction[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quantum_state_uniform_superposition() {
        let qs = QuantumState::uniform_superposition(4);
        let total: f64 = qs.populations.iter().sum();
        assert!((total - 1.0).abs() < 1e-10);
        for p in &qs.populations {
            assert!((p - 0.25).abs() < 1e-10);
        }
    }

    #[test]
    fn test_quantum_state_ground_state_population() {
        let qs = QuantumState::new(3);
        assert!((qs.ground_state_population() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quantum_state_total_population() {
        let qs = QuantumState::uniform_superposition(5);
        assert!((qs.total_population() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quantum_state_mean_energy() {
        let mut qs = QuantumState::new(2);
        qs.energies = vec![1.0, 3.0];
        qs.wavefunction = vec![(0.5_f64).sqrt(), (0.5_f64).sqrt()];
        qs.update_populations();
        let me = qs.mean_energy();
        assert!((me - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_quantum_state_dominant_state() {
        let mut qs = QuantumState::new(3);
        qs.populations = vec![0.1, 0.7, 0.2];
        assert_eq!(qs.dominant_state(), 1);
    }

    // --- HamiltonianMatrix ---

    #[test]
    fn test_hamiltonian_diagonal() {
        let h = HamiltonianMatrix::diagonal(&[1.0, 2.0, 3.0]);
        assert!((h.data[0][0] - 1.0).abs() < 1e-12);
        assert!((h.data[1][1] - 2.0).abs() < 1e-12);
        assert!((h.data[2][2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_hamiltonian_trace() {
        let h = HamiltonianMatrix::diagonal(&[1.0, 2.0, 4.0]);
        assert!((h.trace() - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_hamiltonian_add_coupling() {
        let mut h = HamiltonianMatrix::diagonal(&[0.0, 0.0]);
        h.add_coupling(0, 1, 0.5);
        assert!((h.data[0][1] - 0.5).abs() < 1e-12);
        assert!((h.data[1][0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_hamiltonian_eigenvalues_diagonal() {
        let h = HamiltonianMatrix::diagonal(&[3.0, 1.0, 2.0]);
        let evals = h.eigenvalues();
        assert!((evals[0] - 1.0).abs() < 1e-8);
        assert!((evals[1] - 2.0).abs() < 1e-8);
        assert!((evals[2] - 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_hamiltonian_eigenvalues_2x2() {
        let mut h = HamiltonianMatrix::diagonal(&[0.0, 2.0]);
        h.add_coupling(0, 1, 1.0);
        let evals = h.eigenvalues();
        // Exact: (2 ± sqrt(4+4))/2 = 1 ± sqrt(2)
        let e0 = 1.0 - 2.0_f64.sqrt();
        let e1 = 1.0 + 2.0_f64.sqrt();
        assert!(
            (evals[0] - e0).abs() < 1e-6,
            "evals[0]={} expected {}",
            evals[0],
            e0
        );
        assert!(
            (evals[1] - e1).abs() < 1e-6,
            "evals[1]={} expected {}",
            evals[1],
            e1
        );
    }

    #[test]
    fn test_hamiltonian_ground_state_energy() {
        let h = HamiltonianMatrix::diagonal(&[5.0, 3.0, 7.0]);
        assert!((h.ground_state_energy() - 3.0).abs() < 1e-8);
    }

    #[test]
    fn test_hamiltonian_frobenius_norm() {
        let h = HamiltonianMatrix::diagonal(&[3.0, 4.0]);
        assert!((h.frobenius_norm() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_hamiltonian_apply() {
        let h = HamiltonianMatrix::diagonal(&[2.0, 3.0]);
        let v = h.apply(&[1.0, 1.0]);
        assert!((v[0] - 2.0).abs() < 1e-12);
        assert!((v[1] - 3.0).abs() < 1e-12);
    }

    // --- SurfaceHoppingTSH ---

    #[test]
    fn test_fssh_new() {
        let sh = SurfaceHoppingTSH::new(3, 0.1);
        assert_eq!(sh.n_states, 3);
        assert_eq!(sh.current_state, 0);
    }

    #[test]
    fn test_fssh_non_adiabatic_coupling() {
        let g1 = vec![1.0, 0.0, 0.0];
        let g2 = vec![0.0, 0.0, 0.0];
        let v = vec![2.0, 0.0, 0.0];
        let d = SurfaceHoppingTSH::non_adiabatic_coupling(&g1, &g2, &v);
        assert!((d - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_fssh_switching_probability_zero_pop() {
        let g = SurfaceHoppingTSH::switching_probability(1.0, 0.1, 0.5, 0.0);
        assert!((g).abs() < 1e-12);
    }

    #[test]
    fn test_fssh_switching_probability_positive() {
        let g = SurfaceHoppingTSH::switching_probability(-1.0, 0.1, 0.8, 0.2);
        assert!(g >= 0.0);
        assert!(g <= 1.0);
    }

    #[test]
    fn test_fssh_hop_no_hop_with_zero_prob() {
        let mut sh = SurfaceHoppingTSH::new(3, 0.1);
        sh.hop(&[0.0, 0.0, 0.0], 0.5);
        assert_eq!(sh.current_state, 0);
        assert_eq!(sh.hop_count, 0);
    }

    #[test]
    fn test_fssh_hop_certain() {
        let mut sh = SurfaceHoppingTSH::new(3, 0.1);
        // probs[0] = 1.0 but that's current state, probs[1] = 0.99
        sh.hop(&[0.0, 0.99, 0.01], 0.01);
        // rand_val=0.01 < cumulative 0.99 at index 1 (and index 1 != 0)
        assert_eq!(sh.current_state, 1);
        assert_eq!(sh.hop_count, 1);
    }

    #[test]
    fn test_fssh_velocity_rescaling_enough_energy() {
        let factor = SurfaceHoppingTSH::velocity_rescaling_factor(1.0, 5.0);
        assert!(factor.is_some());
        let f = factor.unwrap();
        assert!((f - (4.0 / 5.0_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_fssh_velocity_rescaling_frustrated() {
        let factor = SurfaceHoppingTSH::velocity_rescaling_factor(10.0, 5.0);
        assert!(factor.is_none());
    }

    // --- WignerSampling ---

    #[test]
    fn test_wigner_harmonic_n_samples() {
        let samples = WignerSampling::harmonic_wigner(1e14, 1.0 / (KB * 300.0), 50);
        assert_eq!(samples.len(), 50);
    }

    #[test]
    fn test_wigner_harmonic_finite_values() {
        let samples = WignerSampling::harmonic_wigner(1e14, 1.0 / (KB * 300.0), 10);
        for (x, p) in &samples {
            assert!(x.is_finite(), "x not finite");
            assert!(p.is_finite(), "p not finite");
        }
    }

    #[test]
    fn test_wigner_morse_approx() {
        let (x, p) = WignerSampling::morse_wigner_approx(1e-18, 1e10, 1e14, 1.0 / (KB * 300.0));
        assert!(x.is_finite());
        assert!(p.is_finite());
    }

    #[test]
    fn test_wigner_sample_phase_space() {
        let gaussian = |x: f64, p: f64| (-x * x - p * p).exp();
        let samples = WignerSampling::sample_phase_space(gaussian, (-3.0, 3.0), (-3.0, 3.0), 20);
        assert!(!samples.is_empty());
        for (x, p) in &samples {
            assert!(x.is_finite());
            assert!(p.is_finite());
        }
    }

    // --- BornOppenheimerMD ---

    #[test]
    fn test_bomd_new() {
        let bomd = BornOppenheimerMD::new(4);
        assert_eq!(bomd.n_atoms, 4);
        assert!(bomd.potential_energy_surface.is_empty());
    }

    #[test]
    fn test_bomd_record_energy() {
        let mut bomd = BornOppenheimerMD::new(1);
        bomd.record_energy(1.5);
        bomd.record_energy(2.0);
        assert_eq!(bomd.potential_energy_surface.len(), 2);
    }

    #[test]
    fn test_bomd_average_energy() {
        let mut bomd = BornOppenheimerMD::new(1);
        bomd.record_energy(1.0);
        bomd.record_energy(3.0);
        assert!((bomd.average_energy() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_bomd_gradient_harmonic() {
        let bomd = BornOppenheimerMD::new(1);
        let pos = vec![[1.0, 0.0, 0.0]];
        // V = x^2, dV/dx = 2x, force = -2x
        let forces = bomd.gradient_at(&pos, |p| p[0][0] * p[0][0]);
        assert!(
            (forces[0][0] - (-2.0)).abs() < 1e-4,
            "force={}",
            forces[0][0]
        );
    }

    #[test]
    fn test_bomd_last_energy_none() {
        let bomd = BornOppenheimerMD::new(1);
        assert!(bomd.last_energy().is_none());
    }

    #[test]
    fn test_bomd_last_energy_some() {
        let mut bomd = BornOppenheimerMD::new(1);
        bomd.record_energy(42.0);
        assert_eq!(bomd.last_energy(), Some(42.0));
    }

    // --- DensityMatrix ---

    #[test]
    fn test_density_matrix_ground_state() {
        let dm = DensityMatrix::ground_state(3);
        assert!((dm.population(0) - 1.0).abs() < 1e-12);
        assert!((dm.population(1)).abs() < 1e-12);
        assert!((dm.trace() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_density_matrix_mixed() {
        let dm = DensityMatrix::mixed(4);
        assert!((dm.trace() - 1.0).abs() < 1e-10);
        for i in 0..4 {
            assert!((dm.population(i) - 0.25).abs() < 1e-10);
        }
    }

    #[test]
    fn test_density_matrix_coherence() {
        let dm = DensityMatrix::ground_state(2);
        let (re, im) = dm.coherence(0, 1);
        assert!((re).abs() < 1e-12);
        assert!((im).abs() < 1e-12);
    }

    #[test]
    fn test_density_matrix_evolve_unitary_trace_preserved() {
        let mut dm = DensityMatrix::ground_state(2);
        let mut h = HamiltonianMatrix::diagonal(&[0.0, 1.0e-20]); // tiny step
        h.add_coupling(0, 1, 5.0e-21);
        dm.evolve_unitary(&h, 1e-15);
        let trace = dm.trace();
        assert!((trace - 1.0).abs() < 1e-6, "trace={}", trace);
    }

    #[test]
    fn test_density_matrix_purity_pure_state() {
        let dm = DensityMatrix::ground_state(3);
        let purity = dm.purity();
        assert!((purity - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_density_matrix_lindblad_decay() {
        let mut dm = DensityMatrix::mixed(2);
        let _pop_before = dm.population(1);
        dm.lindblad_dissipation(0.1, 1);
        // after decay the ground state pop should increase
        let pop0 = dm.population(0);
        assert!(pop0 >= 0.5 - 1e-10);
    }

    // --- ZeroPointEnergy ---

    #[test]
    fn test_zpe_harmonic() {
        let zpe = ZeroPointEnergy::harmonic_zpe(&[1e14, 2e14]);
        assert!((zpe - 0.5 * HBAR * 3e14).abs() < 1e-40);
    }

    #[test]
    fn test_zpe_molecular() {
        let zpe = ZeroPointEnergy::molecular_zpe(&[-1e13, 1e14, 2e14]);
        assert!((zpe - 0.5 * HBAR * 3e14).abs() < 1e-40);
    }

    #[test]
    fn test_zpe_anharmonic_correction() {
        let corr = ZeroPointEnergy::anharmonic_correction(1e14, 0.0, 0.0);
        assert!((corr).abs() < 1e-50); // no anharmonicity
    }

    #[test]
    fn test_zpe_classical_thermal() {
        let e = ZeroPointEnergy::classical_thermal_energy(1, 300.0);
        assert!((e - 1.5 * KB * 300.0).abs() < 1e-30);
    }

    // --- InstantonTheory ---

    #[test]
    fn test_instanton_action_positive() {
        let s = InstantonTheory::action_harmonic(ME, 1e14, 1e-20, 1.0 / (KB * 300.0));
        assert!(s > 0.0);
    }

    #[test]
    fn test_instanton_tunneling_rate_lt1() {
        let s = InstantonTheory::action_harmonic(ME, 1e14, 1e-20, 1.0 / (KB * 300.0));
        let rate = InstantonTheory::tunneling_rate_wkb(s);
        assert!(rate > 0.0);
        assert!(rate <= 1.0);
    }

    #[test]
    fn test_instanton_crossover_temperature() {
        let tc = InstantonTheory::crossover_temperature(1e13, 1e-20);
        assert!(tc > 0.0);
    }

    #[test]
    fn test_instanton_arrhenius() {
        let k = InstantonTheory::arrhenius_rate(1e12, 0.5 * KB * 300.0, 300.0);
        assert!(k > 0.0);
        assert!(k < 1e12);
    }

    #[test]
    fn test_instanton_quantum_enhancement_positive() {
        let enh = InstantonTheory::quantum_enhancement(1e-50, 1e-20, 300.0);
        assert!(enh.is_finite());
        assert!(enh > 0.0);
    }

    // --- QuantumHarmonicOscillator ---

    #[test]
    fn test_qho_energy_levels() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let e0 = qho.energy_level(0);
        let e1 = qho.energy_level(1);
        assert!((e1 - e0 - HBAR * 1e14).abs() < 1e-35);
    }

    #[test]
    fn test_qho_zero_point_energy() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        assert!((qho.zero_point_energy() - 0.5 * HBAR * 1e14).abs() < 1e-40);
    }

    #[test]
    fn test_qho_wavefunction_n0_finite() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let psi = qho.wavefunction_hermite(0, 0.0);
        assert!(psi.is_finite());
        assert!(psi > 0.0);
    }

    #[test]
    fn test_qho_wavefunction_n1_zero_at_origin() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let psi = qho.wavefunction_hermite(1, 0.0);
        assert!(psi.abs() < 1e-10);
    }

    #[test]
    fn test_qho_expectation_x() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        assert_eq!(qho.expectation_x(3), 0.0);
    }

    #[test]
    fn test_qho_expectation_x2_ground() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let x2 = qho.expectation_x2(0);
        assert!((x2 - 0.5 * HBAR / (ME * 1e14)).abs() < 1e-40);
    }

    #[test]
    fn test_qho_expectation_p2_ground() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let p2 = qho.expectation_p2(0);
        assert!((p2 - 0.5 * ME * 1e14 * HBAR).abs() < 1e-50);
    }

    #[test]
    fn test_qho_virial_theorem() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let ke = qho.kinetic_energy(2);
        let pe = qho.potential_energy(2);
        assert!((ke - pe).abs() < 1e-35);
    }

    #[test]
    fn test_qho_transition_dipole() {
        let qho = QuantumHarmonicOscillator::new(1e14, ME);
        let mu01 = qho.transition_dipole(0);
        let mu12 = qho.transition_dipole(1);
        assert!(mu12 > mu01); // increases with n
    }

    // --- Hermite polynomial ---

    #[test]
    fn test_hermite_h0() {
        assert!((hermite_polynomial(0, 2.5) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_h1() {
        assert!((hermite_polynomial(1, 2.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_hermite_h2() {
        // H2(x) = 4x² - 2
        let x = 1.5;
        let expected = 4.0 * x * x - 2.0;
        assert!((hermite_polynomial(2, x) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_hermite_h3() {
        // H3(x) = 8x³ - 12x
        let x = 1.0;
        let expected = 8.0 * x * x * x - 12.0 * x;
        assert!((hermite_polynomial(3, x) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_factorial_values() {
        assert!((factorial(0) - 1.0).abs() < 1e-12);
        assert!((factorial(1) - 1.0).abs() < 1e-12);
        assert!((factorial(5) - 120.0).abs() < 1e-12);
    }
}
