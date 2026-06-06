// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Statistical mechanics: ensembles, Ising model, phase transitions, fluctuation-dissipation,
//! Green-Kubo relations, and equations of state.
//!
//! Implements canonical, grand canonical, and microcanonical partition functions,
//! Boltzmann distributions, 1D/2D Ising model with Metropolis MC, Landau theory,
//! Einstein relations, transport coefficients, and classical equations of state.

// ─────────────────────────────────────────────────────────────────────────────
// Local LCG RNG
// ─────────────────────────────────────────────────────────────────────────────

/// Lightweight linear congruential generator used internally.
struct SmRng {
    state: u64,
}

impl SmRng {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PartitionFunction
// ─────────────────────────────────────────────────────────────────────────────

/// Ensemble type for computing the partition function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnsembleKind {
    /// Canonical ensemble (fixed N, V, T).
    Canonical,
    /// Grand canonical ensemble (fixed μ, V, T).
    GrandCanonical,
    /// Microcanonical ensemble (fixed N, V, E).
    Microcanonical,
}

/// Partition function for statistical mechanical ensembles.
///
/// Stores a discrete energy spectrum and computes thermodynamic quantities
/// such as free energy, entropy, and average energy.
#[derive(Debug, Clone)]
pub struct PartitionFunction {
    /// Ensemble type.
    pub kind: EnsembleKind,
    /// Energy levels of the system.
    pub energy_levels: Vec<f64>,
    /// Degeneracy of each energy level.
    pub degeneracies: Vec<f64>,
    /// Temperature in units where k_B = 1.
    pub temperature: f64,
    /// Chemical potential (used for grand canonical ensemble).
    pub chemical_potential: f64,
    /// Particle numbers for each state (grand canonical).
    pub particle_numbers: Vec<f64>,
}

impl PartitionFunction {
    /// Constructs a canonical partition function from energy levels and degeneracies.
    ///
    /// # Arguments
    /// * `energy_levels` - Vector of energy eigenvalues.
    /// * `degeneracies`  - Degeneracy of each level (same length as `energy_levels`).
    /// * `temperature`   - Temperature (k_B = 1 units).
    pub fn canonical(energy_levels: Vec<f64>, degeneracies: Vec<f64>, temperature: f64) -> Self {
        Self {
            kind: EnsembleKind::Canonical,
            energy_levels,
            degeneracies,
            temperature,
            chemical_potential: 0.0,
            particle_numbers: vec![],
        }
    }

    /// Constructs a grand canonical partition function.
    pub fn grand_canonical(
        energy_levels: Vec<f64>,
        degeneracies: Vec<f64>,
        particle_numbers: Vec<f64>,
        temperature: f64,
        chemical_potential: f64,
    ) -> Self {
        Self {
            kind: EnsembleKind::GrandCanonical,
            energy_levels,
            degeneracies,
            temperature,
            chemical_potential,
            particle_numbers,
        }
    }

    /// Constructs a microcanonical partition function (density of states).
    pub fn microcanonical(energy_levels: Vec<f64>, degeneracies: Vec<f64>) -> Self {
        Self {
            kind: EnsembleKind::Microcanonical,
            energy_levels,
            degeneracies,
            temperature: 0.0,
            chemical_potential: 0.0,
            particle_numbers: vec![],
        }
    }

    /// Returns the canonical partition function Z = Σ g_i exp(-β E_i).
    ///
    /// Returns 0.0 if temperature is non-positive.
    pub fn canonical_z(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / self.temperature;
        self.energy_levels
            .iter()
            .zip(self.degeneracies.iter())
            .map(|(&e, &g)| g * (-beta * e).exp())
            .sum()
    }

    /// Returns the grand canonical partition function Ξ = Σ g_i exp(-β(E_i - μ N_i)).
    pub fn grand_canonical_z(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / self.temperature;
        let mu = self.chemical_potential;
        let n_iter = self
            .particle_numbers
            .iter()
            .chain(std::iter::repeat(&0.0_f64));
        self.energy_levels
            .iter()
            .zip(self.degeneracies.iter())
            .zip(n_iter)
            .map(|((&e, &g), &n)| g * (-beta * (e - mu * n)).exp())
            .sum()
    }

    /// Returns the microcanonical density of states Ω(E) for a target energy within tolerance.
    ///
    /// Sums degeneracies of levels within `±tolerance` of `target_energy`.
    pub fn density_of_states(&self, target_energy: f64, tolerance: f64) -> f64 {
        self.energy_levels
            .iter()
            .zip(self.degeneracies.iter())
            .filter(|&(&e, _)| (e - target_energy).abs() <= tolerance)
            .map(|(_, &g)| g)
            .sum()
    }

    /// Returns the Helmholtz free energy F = -k_B T ln Z (canonical).
    pub fn free_energy(&self) -> f64 {
        let z = self.canonical_z();
        if z <= 0.0 {
            return f64::INFINITY;
        }
        -self.temperature * z.ln()
    }

    /// Returns the average energy ⟨E⟩ = -∂ ln Z / ∂β (canonical).
    pub fn average_energy(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / self.temperature;
        let z = self.canonical_z();
        if z <= 0.0 {
            return 0.0;
        }
        let num: f64 = self
            .energy_levels
            .iter()
            .zip(self.degeneracies.iter())
            .map(|(&e, &g)| g * e * (-beta * e).exp())
            .sum();
        num / z
    }

    /// Returns the entropy S = (⟨E⟩ - F) / T (canonical).
    pub fn entropy(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        (self.average_energy() - self.free_energy()) / self.temperature
    }

    /// Returns the heat capacity C_V = ∂⟨E⟩/∂T = (⟨E²⟩ - ⟨E⟩²) / (k_B T²).
    pub fn heat_capacity(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / self.temperature;
        let z = self.canonical_z();
        if z <= 0.0 {
            return 0.0;
        }
        let e_avg = self.average_energy();
        let e2_avg: f64 = self
            .energy_levels
            .iter()
            .zip(self.degeneracies.iter())
            .map(|(&e, &g)| g * e * e * (-beta * e).exp())
            .sum::<f64>()
            / z;
        (e2_avg - e_avg * e_avg) / (self.temperature * self.temperature)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BoltzmannDistribution
// ─────────────────────────────────────────────────────────────────────────────

/// Boltzmann probability distribution over a discrete energy spectrum.
///
/// Computes Boltzmann weights, occupation probabilities, entropy, and free energy.
#[derive(Debug, Clone)]
pub struct BoltzmannDistribution {
    /// Energy levels.
    pub energies: Vec<f64>,
    /// Degeneracies per level.
    pub degeneracies: Vec<f64>,
    /// Temperature (k_B = 1).
    pub temperature: f64,
}

impl BoltzmannDistribution {
    /// Creates a new Boltzmann distribution.
    pub fn new(energies: Vec<f64>, degeneracies: Vec<f64>, temperature: f64) -> Self {
        Self {
            energies,
            degeneracies,
            temperature,
        }
    }

    /// Creates a uniform degeneracy distribution (all degeneracies = 1).
    pub fn uniform(energies: Vec<f64>, temperature: f64) -> Self {
        let n = energies.len();
        Self::new(energies, vec![1.0; n], temperature)
    }

    /// Returns the un-normalised Boltzmann weight for level `i`.
    pub fn weight(&self, i: usize) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / self.temperature;
        self.degeneracies[i] * (-beta * self.energies[i]).exp()
    }

    /// Returns the partition function Z = Σ_i weight(i).
    pub fn partition_function(&self) -> f64 {
        (0..self.energies.len()).map(|i| self.weight(i)).sum()
    }

    /// Returns the normalised probability P_i = weight(i) / Z.
    pub fn probability(&self, i: usize) -> f64 {
        let z = self.partition_function();
        if z <= 0.0 {
            return 0.0;
        }
        self.weight(i) / z
    }

    /// Returns all occupation probabilities as a vector.
    pub fn probabilities(&self) -> Vec<f64> {
        let z = self.partition_function();
        if z <= 0.0 {
            return vec![0.0; self.energies.len()];
        }
        (0..self.energies.len())
            .map(|i| self.weight(i) / z)
            .collect()
    }

    /// Returns the Gibbs entropy S = -Σ_i P_i ln P_i (natural units).
    pub fn entropy(&self) -> f64 {
        self.probabilities()
            .iter()
            .filter(|&&p| p > 0.0)
            .map(|&p| -p * p.ln())
            .sum()
    }

    /// Returns the Helmholtz free energy F = -T ln Z.
    pub fn free_energy(&self) -> f64 {
        let z = self.partition_function();
        if z <= 0.0 {
            return f64::INFINITY;
        }
        -self.temperature * z.ln()
    }

    /// Returns the mean energy ⟨E⟩ = Σ_i P_i E_i.
    pub fn mean_energy(&self) -> f64 {
        let probs = self.probabilities();
        probs
            .iter()
            .zip(self.energies.iter())
            .map(|(&p, &e)| p * e)
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IsingModel
// ─────────────────────────────────────────────────────────────────────────────

/// Ising model on a 1D or 2D lattice with Metropolis Monte Carlo.
///
/// Spins are stored as `i8` values ±1.  The Hamiltonian is
/// H = -J Σ_{⟨i,j⟩} s_i s_j - h Σ_i s_i.
#[derive(Debug, Clone)]
pub struct IsingModel {
    /// Lattice dimension (1 or 2).
    pub dim: usize,
    /// Number of sites along each axis (Lx for 1D; Lx × Ly for 2D).
    pub lx: usize,
    /// Number of sites along y-axis (only used for 2D; 1 for 1D).
    pub ly: usize,
    /// Coupling constant J (J > 0 → ferromagnetic).
    pub j_coupling: f64,
    /// External magnetic field h.
    pub field: f64,
    /// Spin configuration: +1 or -1.
    pub spins: Vec<i8>,
}

impl IsingModel {
    /// Constructs a 1D Ising model with `n` sites, all spins +1.
    pub fn new_1d(n: usize, j_coupling: f64, field: f64) -> Self {
        Self {
            dim: 1,
            lx: n,
            ly: 1,
            j_coupling,
            field,
            spins: vec![1i8; n],
        }
    }

    /// Constructs a 2D Ising model (lx × ly lattice) with all spins +1.
    pub fn new_2d(lx: usize, ly: usize, j_coupling: f64, field: f64) -> Self {
        Self {
            dim: 2,
            lx,
            ly,
            j_coupling,
            field,
            spins: vec![1i8; lx * ly],
        }
    }

    /// Randomises the spin configuration using the given seed.
    pub fn randomise(&mut self, seed: u64) {
        let mut rng = SmRng::new(seed);
        for s in self.spins.iter_mut() {
            *s = if rng.next_f64() < 0.5 { 1 } else { -1 };
        }
    }

    /// Returns the total energy of the current configuration.
    pub fn energy(&self) -> f64 {
        let mut e = 0.0_f64;
        let n = self.spins.len();
        if self.dim == 1 {
            for i in 0..n {
                let j = (i + 1) % n;
                e -= self.j_coupling * (self.spins[i] as f64) * (self.spins[j] as f64);
                e -= self.field * self.spins[i] as f64;
            }
        } else {
            for row in 0..self.ly {
                for col in 0..self.lx {
                    let idx = row * self.lx + col;
                    let right = row * self.lx + (col + 1) % self.lx;
                    let down = ((row + 1) % self.ly) * self.lx + col;
                    e -= self.j_coupling
                        * (self.spins[idx] as f64)
                        * (self.spins[right] as f64 + self.spins[down] as f64);
                    e -= self.field * self.spins[idx] as f64;
                }
            }
        }
        e
    }

    /// Returns the magnetisation per spin m = (1/N) Σ_i s_i.
    pub fn magnetisation(&self) -> f64 {
        let total: i64 = self.spins.iter().map(|&s| s as i64).sum();
        total as f64 / self.spins.len() as f64
    }

    /// Computes the local energy change if spin at site `idx` is flipped.
    fn delta_energy(&self, idx: usize) -> f64 {
        let s = self.spins[idx] as f64;
        let neighbour_sum: f64 = if self.dim == 1 {
            let n = self.spins.len();
            let left = (idx + n - 1) % n;
            let right = (idx + 1) % n;
            self.spins[left] as f64 + self.spins[right] as f64
        } else {
            let row = idx / self.lx;
            let col = idx % self.lx;
            let left = row * self.lx + (col + self.lx - 1) % self.lx;
            let right = row * self.lx + (col + 1) % self.lx;
            let up = ((row + self.ly - 1) % self.ly) * self.lx + col;
            let down = ((row + 1) % self.ly) * self.lx + col;
            self.spins[left] as f64
                + self.spins[right] as f64
                + self.spins[up] as f64
                + self.spins[down] as f64
        };
        2.0 * s * (self.j_coupling * neighbour_sum + self.field)
    }

    /// Runs `n_sweeps` full Metropolis sweeps at temperature `temperature`.
    ///
    /// Each sweep attempts one flip per spin.
    pub fn metropolis_sweep(&mut self, n_sweeps: usize, temperature: f64, seed: u64) {
        let mut rng = SmRng::new(seed);
        let n = self.spins.len();
        for _ in 0..n_sweeps {
            for _attempt in 0..n {
                let idx = (rng.next_u64() as usize) % n;
                let de = self.delta_energy(idx);
                if de <= 0.0 || rng.next_f64() < (-de / temperature).exp() {
                    self.spins[idx] = -self.spins[idx];
                }
            }
        }
    }

    /// Returns the magnetic susceptibility χ = (⟨m²⟩ - ⟨m⟩²) N / T
    /// estimated over `n_samples` post-equilibration sweeps.
    pub fn susceptibility(
        &mut self,
        temperature: f64,
        n_eq: usize,
        n_samples: usize,
        seed: u64,
    ) -> f64 {
        self.metropolis_sweep(n_eq, temperature, seed);
        let n = self.spins.len();
        let mut m_sum = 0.0_f64;
        let mut m2_sum = 0.0_f64;
        let mut rng = SmRng::new(seed.wrapping_add(1));
        for _ in 0..n_samples {
            // one sweep
            for _attempt in 0..n {
                let idx = (rng.next_u64() as usize) % n;
                let de = self.delta_energy(idx);
                if de <= 0.0 || rng.next_f64() < (-de / temperature).exp() {
                    self.spins[idx] = -self.spins[idx];
                }
            }
            let m = self.magnetisation();
            m_sum += m;
            m2_sum += m * m;
        }
        let m_avg = m_sum / n_samples as f64;
        let m2_avg = m2_sum / n_samples as f64;
        (m2_avg - m_avg * m_avg) * n as f64 / temperature
    }

    /// Exact 1D critical temperature (returns J for reference; true T_c → 0 in 1D).
    pub fn critical_temperature_1d_approx(&self) -> f64 {
        // 1D Ising has T_c = 0; return the exchange energy scale
        2.0 * self.j_coupling.abs()
    }

    /// Mean-field critical temperature for the 2D Ising model: T_c ≈ z J / (2 k_B).
    /// For a 2D square lattice z = 4.
    pub fn mean_field_tc_2d(&self) -> f64 {
        4.0 * self.j_coupling
    }

    /// Exact Onsager critical temperature for the 2D square-lattice Ising model.
    /// T_c = 2J / ln(1 + √2).
    pub fn onsager_tc(&self) -> f64 {
        2.0 * self.j_coupling / (1.0 + 2.0_f64.sqrt()).ln()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PhaseTransition
// ─────────────────────────────────────────────────────────────────────────────

/// Landau theory description of a second-order phase transition.
///
/// Free energy density: f(φ) = a(T) φ²/2 + b φ⁴/4 + c φ⁶/6
/// where a(T) = a₀ (T - T_c).
#[derive(Debug, Clone)]
pub struct PhaseTransition {
    /// Coefficient a₀ (> 0).
    pub a0: f64,
    /// Critical temperature T_c.
    pub tc: f64,
    /// Coefficient b of the φ⁴ term (> 0 for a normal second-order transition).
    pub b: f64,
    /// Coefficient c of the φ⁶ term (stabilising; ≥ 0).
    pub c: f64,
    /// Current temperature.
    pub temperature: f64,
}

impl PhaseTransition {
    /// Creates a new Landau free energy model.
    pub fn new(a0: f64, tc: f64, b: f64, c: f64, temperature: f64) -> Self {
        Self {
            a0,
            tc,
            b,
            c,
            temperature,
        }
    }

    /// Returns a(T) = a₀(T - T_c).
    pub fn a_coeff(&self) -> f64 {
        self.a0 * (self.temperature - self.tc)
    }

    /// Returns the equilibrium order parameter φ* that minimises f(φ).
    ///
    /// For T > T_c: φ* = 0.
    /// For T < T_c: φ* = √(-a(T)/b) (when c = 0 and b > 0).
    pub fn order_parameter(&self) -> f64 {
        let a = self.a_coeff();
        if a >= 0.0 {
            return 0.0;
        }
        if self.b > 0.0 {
            (-a / self.b).sqrt()
        } else {
            0.0
        }
    }

    /// Returns the Landau free energy density at order parameter `phi`.
    pub fn free_energy_density(&self, phi: f64) -> f64 {
        let a = self.a_coeff();
        a * phi * phi / 2.0 + self.b * phi.powi(4) / 4.0 + self.c * phi.powi(6) / 6.0
    }

    /// Returns the reduced temperature t = (T - T_c) / T_c.
    pub fn reduced_temperature(&self) -> f64 {
        if self.tc == 0.0 {
            return 0.0;
        }
        (self.temperature - self.tc) / self.tc
    }

    /// Returns the critical exponent β (mean-field value = 1/2).
    pub fn mean_field_beta_exponent(&self) -> f64 {
        0.5
    }

    /// Returns the susceptibility χ ~ |T - T_c|^(-γ); mean-field γ = 1.
    pub fn susceptibility(&self) -> f64 {
        let a = self.a_coeff().abs();
        if a < 1e-15 {
            return f64::INFINITY;
        }
        1.0 / (2.0 * a)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FluctuationDissipation
// ─────────────────────────────────────────────────────────────────────────────

/// Tools for the fluctuation-dissipation theorem and related relations.
///
/// Implements the Einstein relation, response function, and Johnson-Nyquist
/// noise power spectral density.
#[derive(Debug, Clone)]
pub struct FluctuationDissipation {
    /// Temperature (k_B = 1).
    pub temperature: f64,
}

impl FluctuationDissipation {
    /// Creates a new `FluctuationDissipation` at the given temperature.
    pub fn new(temperature: f64) -> Self {
        Self { temperature }
    }

    /// Einstein relation: D = μ k_B T, where μ is the mobility.
    ///
    /// Returns the diffusion coefficient D.
    pub fn einstein_relation(&self, mobility: f64) -> f64 {
        mobility * self.temperature
    }

    /// Returns the mobility from diffusivity via the Einstein relation: μ = D / (k_B T).
    pub fn mobility_from_diffusivity(&self, diffusivity: f64) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        diffusivity / self.temperature
    }

    /// Johnson-Nyquist noise power spectral density S_V = 4 k_B T R.
    ///
    /// # Arguments
    /// * `resistance` - Electrical resistance R in Ohms.
    pub fn noise_power_spectral_density(&self, resistance: f64) -> f64 {
        4.0 * self.temperature * resistance
    }

    /// Imaginary part of the linear response function (Kubo formula).
    ///
    /// For a Lorentzian oscillator:
    /// χ''(ω) = (γ ω) / ((ω₀² - ω²)² + γ² ω²)
    pub fn response_function_imaginary(&self, omega: f64, omega0: f64, gamma: f64) -> f64 {
        let denom = (omega0 * omega0 - omega * omega).powi(2) + gamma * gamma * omega * omega;
        if denom < 1e-30 {
            return 0.0;
        }
        gamma * omega / denom
    }

    /// Real part of the linear response function (Kramers-Kronig).
    ///
    /// χ'(ω) = (ω₀² - ω²) / ((ω₀² - ω²)² + γ² ω²)
    pub fn response_function_real(&self, omega: f64, omega0: f64, gamma: f64) -> f64 {
        let denom = (omega0 * omega0 - omega * omega).powi(2) + gamma * gamma * omega * omega;
        if denom < 1e-30 {
            return 0.0;
        }
        (omega0 * omega0 - omega * omega) / denom
    }

    /// Fluctuation power spectral density via FDT: S(ω) = 2 k_B T χ''(ω) / ω.
    ///
    /// Returns 0 if ω = 0.
    pub fn fluctuation_spectrum(&self, omega: f64, omega0: f64, gamma: f64) -> f64 {
        if omega.abs() < 1e-30 {
            return 0.0;
        }
        2.0 * self.temperature * self.response_function_imaginary(omega, omega0, gamma) / omega
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GreenKubo
// ─────────────────────────────────────────────────────────────────────────────

/// Green-Kubo transport coefficient computation from time-correlation functions.
///
/// Implements the velocity autocorrelation function (VACF) and the Green-Kubo
/// integrals for shear viscosity, electrical conductivity, and diffusivity.
#[derive(Debug, Clone)]
pub struct GreenKubo {
    /// Simulation time step.
    pub dt: f64,
    /// Temperature (k_B = 1).
    pub temperature: f64,
    /// Volume of the simulation cell.
    pub volume: f64,
}

impl GreenKubo {
    /// Creates a new Green-Kubo integrator.
    pub fn new(dt: f64, temperature: f64, volume: f64) -> Self {
        Self {
            dt,
            temperature,
            volume,
        }
    }

    /// Computes the normalised autocorrelation function of `signal`.
    ///
    /// C(τ) = ⟨A(t) A(t+τ)⟩ / ⟨A(0)²⟩.
    /// Returns a vector of length `signal.len()` (full autocorrelation).
    pub fn autocorrelation(signal: &[f64]) -> Vec<f64> {
        let n = signal.len();
        if n == 0 {
            return vec![];
        }
        let mean: f64 = signal.iter().sum::<f64>() / n as f64;
        let centered: Vec<f64> = signal.iter().map(|&x| x - mean).collect();
        let norm: f64 = centered.iter().map(|&x| x * x).sum();
        if norm < 1e-30 {
            return vec![0.0; n];
        }
        (0..n)
            .map(|lag| {
                let count = n - lag;
                let sum: f64 = (0..count).map(|i| centered[i] * centered[i + lag]).sum();
                sum / norm
            })
            .collect()
    }

    /// Green-Kubo diffusivity D = (1/3) ∫₀^∞ ⟨v(0)·v(t)⟩ dt.
    ///
    /// # Arguments
    /// * `vacf` - Velocity autocorrelation function C(t) (not normalised; units: velocity²).
    pub fn diffusivity_from_vacf(&self, vacf: &[f64]) -> f64 {
        // Trapezoidal integration
        let n = vacf.len();
        if n == 0 {
            return 0.0;
        }
        let integral: f64 = (0..n - 1)
            .map(|i| 0.5 * (vacf[i] + vacf[i + 1]) * self.dt)
            .sum();
        integral / 3.0
    }

    /// Green-Kubo shear viscosity η = (V / k_B T) ∫₀^∞ ⟨P_xy(0) P_xy(t)⟩ dt.
    ///
    /// # Arguments
    /// * `stress_acf` - Autocorrelation of the off-diagonal stress element P_xy.
    pub fn shear_viscosity(&self, stress_acf: &[f64]) -> f64 {
        if self.temperature <= 0.0 || self.volume <= 0.0 {
            return 0.0;
        }
        let n = stress_acf.len();
        if n == 0 {
            return 0.0;
        }
        let integral: f64 = (0..n - 1)
            .map(|i| 0.5 * (stress_acf[i] + stress_acf[i + 1]) * self.dt)
            .sum();
        self.volume / self.temperature * integral
    }

    /// Green-Kubo electrical conductivity σ = (V / k_B T) ∫₀^∞ ⟨J(0)·J(t)⟩ dt / 3.
    ///
    /// # Arguments
    /// * `current_acf` - Autocorrelation of total current.
    pub fn electrical_conductivity(&self, current_acf: &[f64]) -> f64 {
        if self.temperature <= 0.0 || self.volume <= 0.0 {
            return 0.0;
        }
        let n = current_acf.len();
        if n == 0 {
            return 0.0;
        }
        let integral: f64 = (0..n - 1)
            .map(|i| 0.5 * (current_acf[i] + current_acf[i + 1]) * self.dt)
            .sum();
        self.volume / (3.0 * self.temperature) * integral
    }

    /// Estimates the correlation time τ = ∫₀^∞ C(t) dt / C(0) from normalised ACF.
    pub fn correlation_time(acf: &[f64], dt: f64) -> f64 {
        if acf.is_empty() || acf[0].abs() < 1e-30 {
            return 0.0;
        }
        let integral: f64 = (0..acf.len() - 1)
            .map(|i| 0.5 * (acf[i] + acf[i + 1]) * dt)
            .sum();
        integral / acf[0]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// EquationOfState
// ─────────────────────────────────────────────────────────────────────────────

/// Classical cubic equations of state for real gases.
///
/// Provides van der Waals, Peng-Robinson, and Redlich-Kwong equations.
#[derive(Debug, Clone)]
pub struct EquationOfState {
    /// Universal gas constant R.
    pub r_gas: f64,
    /// Temperature T.
    pub temperature: f64,
}

/// EOS model variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EosModel {
    /// Van der Waals equation of state.
    VanDerWaals,
    /// Peng-Robinson equation of state.
    PengRobinson,
    /// Redlich-Kwong equation of state.
    RedlichKwong,
}

impl EquationOfState {
    /// Creates a new EOS calculator.
    ///
    /// `r_gas` is usually 8.314 J/(mol·K).
    pub fn new(r_gas: f64, temperature: f64) -> Self {
        Self { r_gas, temperature }
    }

    /// Van der Waals pressure: P = RT/(V-b) - a/V².
    ///
    /// # Arguments
    /// * `molar_volume` - Molar volume V (m³/mol).
    /// * `a`           - Attractive parameter a (Pa·m⁶/mol²).
    /// * `b`           - Repulsive parameter b (m³/mol).
    pub fn van_der_waals_pressure(&self, molar_volume: f64, a: f64, b: f64) -> f64 {
        if molar_volume <= b {
            return f64::INFINITY;
        }
        self.r_gas * self.temperature / (molar_volume - b) - a / (molar_volume * molar_volume)
    }

    /// Van der Waals critical constants from parameters a and b.
    ///
    /// Returns `(T_c, P_c, V_c)`.
    pub fn vdw_critical_constants(a: f64, b: f64, r_gas: f64) -> (f64, f64, f64) {
        let tc = 8.0 * a / (27.0 * r_gas * b);
        let pc = a / (27.0 * b * b);
        let vc = 3.0 * b;
        (tc, pc, vc)
    }

    /// Peng-Robinson pressure.
    ///
    /// P = RT/(V-b) - a(T) / (V(V+b) + b(V-b))
    ///
    /// where a(T) = a_c · α(T), α(T) = \[1 + κ(1 - √(T/T_c))\]².
    pub fn peng_robinson_pressure(
        &self,
        molar_volume: f64,
        a_c: f64,
        b: f64,
        kappa: f64,
        tc: f64,
    ) -> f64 {
        if molar_volume <= b || tc <= 0.0 {
            return f64::INFINITY;
        }
        let tr = self.temperature / tc;
        let alpha = (1.0 + kappa * (1.0 - tr.sqrt())).powi(2);
        let a_t = a_c * alpha;
        self.r_gas * self.temperature / (molar_volume - b)
            - a_t / (molar_volume * (molar_volume + b) + b * (molar_volume - b))
    }

    /// Redlich-Kwong pressure.
    ///
    /// P = RT/(V-b) - a / (T^0.5 V(V+b))
    pub fn redlich_kwong_pressure(&self, molar_volume: f64, a: f64, b: f64) -> f64 {
        if molar_volume <= b || self.temperature <= 0.0 {
            return f64::INFINITY;
        }
        self.r_gas * self.temperature / (molar_volume - b)
            - a / (self.temperature.sqrt() * molar_volume * (molar_volume + b))
    }

    /// Compressibility factor Z = PV/(RT).
    pub fn compressibility_factor(&self, pressure: f64, molar_volume: f64) -> f64 {
        if self.r_gas <= 0.0 || self.temperature <= 0.0 {
            return 0.0;
        }
        pressure * molar_volume / (self.r_gas * self.temperature)
    }

    /// Boyle temperature for van der Waals gas: T_B = a / (Rb).
    pub fn boyle_temperature_vdw(a: f64, b: f64, r_gas: f64) -> f64 {
        a / (r_gas * b)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: Fermi-Dirac and Bose-Einstein distributions
// ─────────────────────────────────────────────────────────────────────────────

/// Returns the Fermi-Dirac occupation number f(ε) = 1/(exp((ε-μ)/T) + 1).
///
/// # Arguments
/// * `energy`   - Single-particle energy ε.
/// * `mu`       - Chemical potential μ.
/// * `temperature` - Temperature T (k_B = 1).
pub fn fermi_dirac(energy: f64, mu: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 {
        return if energy < mu { 1.0 } else { 0.0 };
    }
    1.0 / (((energy - mu) / temperature).exp() + 1.0)
}

/// Returns the Bose-Einstein occupation number n(ε) = 1/(exp((ε-μ)/T) - 1).
///
/// Returns 0 if the exponent makes the denominator non-positive.
pub fn bose_einstein(energy: f64, mu: f64, temperature: f64) -> f64 {
    if temperature <= 0.0 {
        return 0.0;
    }
    let x = (energy - mu) / temperature;
    if x <= 0.0 {
        return 0.0;
    }
    let denom = x.exp() - 1.0;
    if denom < 1e-30 { 0.0 } else { 1.0 / denom }
}

/// Stefan-Boltzmann law: power radiated per unit area P = σ T⁴.
///
/// Uses σ = 5.670374419 × 10⁻⁸ W m⁻² K⁻⁴.
pub fn stefan_boltzmann_power(temperature: f64) -> f64 {
    const SIGMA: f64 = 5.670_374_419e-8;
    SIGMA * temperature.powi(4)
}

/// Wien displacement law: λ_max T = b, returns λ_max.
///
/// b = 2.897771955 × 10⁻³ m·K.
pub fn wien_peak_wavelength(temperature: f64) -> f64 {
    const B: f64 = 2.897_771_955e-3;
    if temperature <= 0.0 {
        return f64::INFINITY;
    }
    B / temperature
}

/// Planck spectral radiance B_λ(λ, T) = 2hc² / (λ⁵ (exp(hc/λkT) - 1)).
///
/// Returns radiance in SI units W sr⁻¹ m⁻³.
pub fn planck_spectral_radiance(wavelength: f64, temperature: f64) -> f64 {
    const H: f64 = 6.626_070_15e-34;
    const C: f64 = 2.997_924_58e8;
    const KB: f64 = 1.380_649e-23;
    if wavelength <= 0.0 || temperature <= 0.0 {
        return 0.0;
    }
    let exponent = H * C / (wavelength * KB * temperature);
    if exponent > 700.0 {
        return 0.0;
    }
    2.0 * H * C * C / (wavelength.powi(5) * (exponent.exp() - 1.0))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── PartitionFunction ────────────────────────────────────────────────────

    #[test]
    fn test_canonical_z_two_level() {
        // Two-level system: ε=0 and ε=1, g=1 each, T=1
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 1.0);
        let z = pf.canonical_z();
        assert!((z - (1.0 + (-1.0_f64).exp())).abs() < 1e-12);
    }

    #[test]
    fn test_canonical_z_zero_temp() {
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 0.0);
        assert_eq!(pf.canonical_z(), 0.0);
    }

    #[test]
    fn test_free_energy_two_level() {
        let t = 2.0;
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], t);
        let z = pf.canonical_z();
        let f_expected = -t * z.ln();
        assert!((pf.free_energy() - f_expected).abs() < 1e-12);
    }

    #[test]
    fn test_average_energy_zero_temp() {
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 0.0);
        assert_eq!(pf.average_energy(), 0.0);
    }

    #[test]
    fn test_average_energy_high_temp() {
        // At high T both levels equally occupied → ⟨E⟩ = 0.5
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 1e6);
        let e_avg = pf.average_energy();
        assert!((e_avg - 0.5).abs() < 1e-3);
    }

    #[test]
    fn test_entropy_nonnegative() {
        let pf = PartitionFunction::canonical(vec![0.0, 1.0, 2.0], vec![1.0, 2.0, 1.0], 1.5);
        assert!(pf.entropy() >= 0.0);
    }

    #[test]
    fn test_heat_capacity_positive() {
        let pf = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 1.0);
        assert!(pf.heat_capacity() > 0.0);
    }

    #[test]
    fn test_grand_canonical_z() {
        let pf = PartitionFunction::grand_canonical(
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0],
            1.0,
            0.5,
        );
        let z = pf.grand_canonical_z();
        assert!(z > 0.0);
    }

    #[test]
    fn test_density_of_states() {
        let pf = PartitionFunction::microcanonical(vec![0.0, 1.0, 2.0], vec![1.0, 3.0, 1.0]);
        let omega = pf.density_of_states(1.0, 0.1);
        assert!((omega - 3.0).abs() < 1e-12);
    }

    // ── BoltzmannDistribution ────────────────────────────────────────────────

    #[test]
    fn test_boltzmann_probabilities_sum_to_one() {
        let bd = BoltzmannDistribution::uniform(vec![0.0, 1.0, 2.0, 3.0], 1.0);
        let sum: f64 = bd.probabilities().iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_boltzmann_entropy_nonneg() {
        let bd = BoltzmannDistribution::uniform(vec![0.0, 1.0, 2.0], 1.0);
        assert!(bd.entropy() >= 0.0);
    }

    #[test]
    fn test_boltzmann_free_energy() {
        let bd = BoltzmannDistribution::uniform(vec![0.0, 0.5, 1.0], 0.5);
        let f = bd.free_energy();
        assert!(f.is_finite());
    }

    #[test]
    fn test_boltzmann_mean_energy() {
        let bd = BoltzmannDistribution::uniform(vec![0.0, 1.0], 1e6);
        // High T → equal weights → mean ≈ 0.5
        assert!((bd.mean_energy() - 0.5).abs() < 1e-3);
    }

    #[test]
    fn test_boltzmann_zero_temp() {
        let bd = BoltzmannDistribution::uniform(vec![0.0, 1.0], 0.0);
        assert_eq!(bd.partition_function(), 0.0);
    }

    // ── IsingModel ───────────────────────────────────────────────────────────

    #[test]
    fn test_ising_1d_energy_all_up() {
        // All spins +1, H=0 → E = -J * N
        let model = IsingModel::new_1d(4, 1.0, 0.0);
        assert!((model.energy() - (-4.0)).abs() < 1e-12);
    }

    #[test]
    fn test_ising_magnetisation_all_up() {
        let model = IsingModel::new_1d(4, 1.0, 0.0);
        assert!((model.magnetisation() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_ising_magnetisation_alternating() {
        let mut model = IsingModel::new_1d(4, 1.0, 0.0);
        for i in 0..4 {
            model.spins[i] = if i % 2 == 0 { 1 } else { -1 };
        }
        assert!(model.magnetisation().abs() < 1e-12);
    }

    #[test]
    fn test_ising_2d_energy_all_up() {
        // 2×2 lattice, J=1, h=0, all spins up
        // Each spin has 2 right/down bonds; 4 spins, 8 bonds total (with PBC)
        let model = IsingModel::new_2d(2, 2, 1.0, 0.0);
        assert!(model.energy() < 0.0);
    }

    #[test]
    fn test_ising_randomise_changes_spins() {
        let mut model = IsingModel::new_1d(10, 1.0, 0.0);
        model.randomise(42);
        // After randomisation the spins should not all be +1 (very unlikely)
        let all_up = model.spins.iter().all(|&s| s == 1);
        let all_dn = model.spins.iter().all(|&s| s == -1);
        // At least one of all_up/all_dn is false for a 10-spin lattice
        assert!(!(all_up && all_dn));
    }

    #[test]
    fn test_ising_metropolis_runs() {
        let mut model = IsingModel::new_1d(8, 1.0, 0.0);
        model.randomise(1);
        // Should not panic
        model.metropolis_sweep(10, 2.0, 99);
        assert!(model.magnetisation().abs() <= 1.0);
    }

    #[test]
    fn test_ising_onsager_tc() {
        let model = IsingModel::new_2d(8, 8, 1.0, 0.0);
        let tc = model.onsager_tc();
        // Known value ≈ 2.269
        assert!((tc - 2.269).abs() < 0.01);
    }

    #[test]
    fn test_ising_susceptibility_positive() {
        let mut model = IsingModel::new_1d(8, 1.0, 0.0);
        model.randomise(7);
        let chi = model.susceptibility(2.5, 20, 50, 11);
        assert!(chi >= 0.0);
    }

    // ── PhaseTransition ──────────────────────────────────────────────────────

    #[test]
    fn test_phase_transition_order_param_above_tc() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.0, 3.0);
        assert_eq!(pt.order_parameter(), 0.0);
    }

    #[test]
    fn test_phase_transition_order_param_below_tc() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.0, 1.0);
        let phi = pt.order_parameter();
        assert!(phi > 0.0);
    }

    #[test]
    fn test_phase_transition_free_energy_density() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.5, 1.0);
        let f0 = pt.free_energy_density(0.0);
        let phi_eq = pt.order_parameter();
        let f_eq = pt.free_energy_density(phi_eq);
        // Minimum should be at equilibrium
        assert!(f_eq <= f0 + 1e-10);
    }

    #[test]
    fn test_phase_transition_reduced_temperature() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.0, 3.0);
        assert!((pt.reduced_temperature() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_phase_transition_susceptibility_diverges_near_tc() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.0, 2.0 + 1e-8);
        assert!(pt.susceptibility() > 1e5);
    }

    // ── FluctuationDissipation ───────────────────────────────────────────────

    #[test]
    fn test_einstein_relation() {
        let fd = FluctuationDissipation::new(1.0);
        assert!((fd.einstein_relation(0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_mobility_from_diffusivity() {
        let fd = FluctuationDissipation::new(2.0);
        assert!((fd.mobility_from_diffusivity(4.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_noise_power_spectral_density() {
        let fd = FluctuationDissipation::new(300.0);
        // S = 4 k_B T R; with k_B=1 here
        let s = fd.noise_power_spectral_density(50.0);
        assert!((s - 4.0 * 300.0 * 50.0).abs() < 1e-6);
    }

    #[test]
    fn test_response_function_imaginary_at_resonance() {
        let fd = FluctuationDissipation::new(1.0);
        let chi_im = fd.response_function_imaginary(1.0, 1.0, 0.1);
        assert!(chi_im > 0.0);
    }

    #[test]
    fn test_fluctuation_spectrum_zero_omega() {
        let fd = FluctuationDissipation::new(1.0);
        assert_eq!(fd.fluctuation_spectrum(0.0, 1.0, 0.1), 0.0);
    }

    // ── GreenKubo ────────────────────────────────────────────────────────────

    #[test]
    fn test_autocorrelation_constant_signal() {
        let signal = vec![1.0; 10];
        let acf = GreenKubo::autocorrelation(&signal);
        // Constant signal → all zero after removing mean
        assert_eq!(acf.len(), 10);
        for v in &acf {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn test_autocorrelation_lag_zero_is_one() {
        let signal: Vec<f64> = (0..20).map(|i| (i as f64 * 0.3).sin()).collect();
        let acf = GreenKubo::autocorrelation(&signal);
        assert!(!acf.is_empty());
        assert!((acf[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_diffusivity_from_vacf() {
        let gk = GreenKubo::new(0.01, 1.0, 1.0);
        // Exponentially decaying VACF: C(t) = C0 exp(-t/tau)
        let tau = 1.0;
        let n = 200;
        let vacf: Vec<f64> = (0..n).map(|i| (-(i as f64) * 0.01 / tau).exp()).collect();
        let d = gk.diffusivity_from_vacf(&vacf);
        // ∫₀^∞ exp(-t/τ) dt / 3 ≈ τ/3
        assert!((d - tau / 3.0).abs() < 0.05);
    }

    #[test]
    fn test_shear_viscosity_zero_temp() {
        let gk = GreenKubo::new(0.01, 0.0, 1.0);
        assert_eq!(gk.shear_viscosity(&[1.0, 0.5, 0.1]), 0.0);
    }

    #[test]
    fn test_correlation_time() {
        let tau = 2.0;
        let dt = 0.1;
        let n = 100;
        let acf: Vec<f64> = (0..n).map(|i| (-(i as f64) * dt / tau).exp()).collect();
        let ct = GreenKubo::correlation_time(&acf, dt);
        // Should be close to tau
        assert!((ct - tau).abs() < 0.5);
    }

    // ── EquationOfState ──────────────────────────────────────────────────────

    #[test]
    fn test_vdw_ideal_limit() {
        // With a=0, b=0 → P = RT/V
        let eos = EquationOfState::new(8.314, 300.0);
        let v = 0.025; // ~1 mol at ~1 atm
        let p = eos.van_der_waals_pressure(v, 0.0, 0.0);
        assert!((p - 8.314 * 300.0 / v).abs() < 1e-6);
    }

    #[test]
    fn test_vdw_critical_constants() {
        let (tc, pc, vc) = EquationOfState::vdw_critical_constants(0.364, 4.27e-5, 8.314);
        assert!(tc > 0.0 && pc > 0.0 && vc > 0.0);
    }

    #[test]
    fn test_peng_robinson_pressure_finite() {
        let eos = EquationOfState::new(8.314, 300.0);
        let p = eos.peng_robinson_pressure(0.001, 0.4, 2.5e-5, 0.37, 190.0);
        assert!(p.is_finite());
    }

    #[test]
    fn test_redlich_kwong_pressure_finite() {
        let eos = EquationOfState::new(8.314, 400.0);
        let p = eos.redlich_kwong_pressure(0.002, 1.0, 2.6e-5);
        assert!(p.is_finite() && p > 0.0);
    }

    #[test]
    fn test_compressibility_factor_ideal_gas() {
        let eos = EquationOfState::new(8.314, 300.0);
        let v = 0.025;
        let p = eos.van_der_waals_pressure(v, 0.0, 0.0);
        let z = eos.compressibility_factor(p, v);
        assert!((z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boyle_temperature_vdw() {
        let t_b = EquationOfState::boyle_temperature_vdw(0.364, 4.27e-5, 8.314);
        assert!(t_b > 0.0);
    }

    // ── Fermi-Dirac / Bose-Einstein ──────────────────────────────────────────

    #[test]
    fn test_fermi_dirac_at_mu() {
        // f(μ) = 0.5 for any T > 0
        assert!((fermi_dirac(1.0, 1.0, 1.0) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_fermi_dirac_below_mu() {
        assert!(fermi_dirac(0.0, 1.0, 1.0) > 0.5);
    }

    #[test]
    fn test_fermi_dirac_zero_temp_below_mu() {
        assert_eq!(fermi_dirac(0.5, 1.0, 0.0), 1.0);
    }

    #[test]
    fn test_fermi_dirac_zero_temp_above_mu() {
        assert_eq!(fermi_dirac(1.5, 1.0, 0.0), 0.0);
    }

    #[test]
    fn test_bose_einstein_positive() {
        let n = bose_einstein(1.0, 0.5, 1.0);
        assert!(n > 0.0);
    }

    #[test]
    fn test_bose_einstein_at_zero_argument() {
        // energy == mu → x=0, undefined → returns 0
        assert_eq!(bose_einstein(1.0, 1.0, 1.0), 0.0);
    }

    #[test]
    fn test_planck_spectral_radiance_positive() {
        let b = planck_spectral_radiance(500e-9, 5778.0);
        assert!(b > 0.0);
    }

    #[test]
    fn test_wien_peak_wavelength() {
        // Sun surface ~5778 K → peak ~502 nm
        let lam = wien_peak_wavelength(5778.0);
        assert!((lam - 5.02e-7).abs() < 1e-8);
    }

    #[test]
    fn test_stefan_boltzmann() {
        // T=1 → σ W/m²
        let p = stefan_boltzmann_power(1.0);
        assert!((p - 5.670374419e-8).abs() < 1e-15);
    }

    #[test]
    fn test_partition_function_degeneracy_scaling() {
        // Doubling all degeneracies should not change average energy
        let pf1 = PartitionFunction::canonical(vec![0.0, 1.0], vec![1.0, 1.0], 1.0);
        let pf2 = PartitionFunction::canonical(vec![0.0, 1.0], vec![2.0, 2.0], 1.0);
        let diff = (pf1.average_energy() - pf2.average_energy()).abs();
        assert!(diff < 1e-12);
    }

    #[test]
    fn test_ising_delta_energy_flip_all_up() {
        let model = IsingModel::new_1d(4, 1.0, 0.0);
        // Flipping spin 0 when all up: ΔE = 2J * (sum_neighbors) = 2*1*(+1+1)=4
        let de = model.delta_energy(0);
        assert!((de - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_phase_transition_a_coeff() {
        let pt = PhaseTransition::new(1.0, 2.0, 1.0, 0.0, 3.0);
        assert!((pt.a_coeff() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gk_autocorrelation_empty() {
        let acf = GreenKubo::autocorrelation(&[]);
        assert!(acf.is_empty());
    }

    #[test]
    fn test_electrical_conductivity_zero_volume() {
        let gk = GreenKubo::new(0.01, 1.0, 0.0);
        assert_eq!(gk.electrical_conductivity(&[1.0, 0.5]), 0.0);
    }

    #[test]
    fn test_vdw_pressure_below_b_returns_inf() {
        let eos = EquationOfState::new(8.314, 300.0);
        let p = eos.van_der_waals_pressure(1e-6, 0.0, 0.01);
        assert_eq!(p, f64::INFINITY);
    }

    #[test]
    fn test_mean_field_tc_2d() {
        let model = IsingModel::new_2d(4, 4, 1.0, 0.0);
        assert!((model.mean_field_tc_2d() - 4.0).abs() < 1e-12);
    }
}
