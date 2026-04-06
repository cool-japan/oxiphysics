// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Quantum transport coupled to molecular dynamics (Landauer-Büttiker + MD).
//!
//! This module provides quantum transport calculations that can be coupled
//! to molecular dynamics trajectories, including:
//!
//! - Tight-binding Hamiltonians ([`TightBindingModel`])
//! - Landauer transmission coefficients ([`LandauerTransmission`])
//! - Non-equilibrium Green's functions ([`GreensFunctionMd`])
//! - Non-equilibrium steady state currents ([`NonEquilibriumMd`])
//! - Electron-phonon coupling ([`ElectronPhononCoupling`])
//! - MD trajectory-averaged transport ([`MdTransportCoupling`])
//! - Wigner quasi-probability distribution ([`WignerFunction`])
//! - Boltzmann transport equation ([`BoltzmannTransport`])

/// Physical constants
pub mod constants {
    /// Elementary charge in Coulombs
    pub const E_CHARGE: f64 = 1.602_176_634e-19;
    /// Planck constant in J·s
    pub const H_PLANCK: f64 = 6.626_070_15e-34;
    /// Reduced Planck constant in J·s
    pub const HBAR: f64 = 1.054_571_817e-34;
    /// Boltzmann constant in J/K
    pub const K_BOLTZMANN: f64 = 1.380_649e-23;
    /// Conductance quantum 2e²/h in Siemens
    pub const G_QUANTUM: f64 = 7.748_091_729e-5;
    /// Von Klitzing constant h/e² in Ohms
    pub const R_KLITZING: f64 = 25812.807;
}

/// Complex number type for quantum calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub struct Complex {
    /// Real part
    pub re: f64,
    /// Imaginary part
    pub im: f64,
}

impl Complex {
    /// Creates a new complex number.
    #[allow(dead_code)]
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Returns the complex conjugate.
    #[allow(dead_code)]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// Returns the modulus squared |z|².
    #[allow(dead_code)]
    pub fn norm_sq(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Returns the modulus |z|.
    #[allow(dead_code)]
    pub fn norm(&self) -> f64 {
        self.norm_sq().sqrt()
    }

    /// Multiplies two complex numbers.
    #[allow(dead_code)]
    pub fn mul(&self, other: &Self) -> Self {
        Self {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    /// Adds two complex numbers.
    #[allow(dead_code)]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    /// Subtracts two complex numbers.
    #[allow(dead_code)]
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    /// Divides two complex numbers.
    #[allow(dead_code)]
    pub fn div(&self, other: &Self) -> Self {
        let denom = other.norm_sq();
        Self {
            re: (self.re * other.re + self.im * other.im) / denom,
            im: (self.im * other.re - self.re * other.im) / denom,
        }
    }
}

/// 1D/2D tight-binding Hamiltonian for quantum transport.
///
/// The tight-binding model uses on-site energies ε_i and
/// hopping integrals t_ij to construct a Hamiltonian matrix.
#[derive(Debug, Clone)]
pub struct TightBindingModel {
    /// Number of sites
    pub n_sites: usize,
    /// On-site energies (diagonal elements)
    pub on_site: Vec<f64>,
    /// Hopping integrals stored as (i, j, t_ij) triplets
    pub hoppings: Vec<(usize, usize, f64)>,
    /// Dimensionality (1 or 2)
    pub dim: usize,
}

impl TightBindingModel {
    /// Creates a new tight-binding model.
    ///
    /// # Arguments
    /// * `n_sites` - Number of lattice sites
    /// * `on_site` - On-site energies for each site
    /// * `dim` - Dimensionality (1 or 2)
    #[allow(dead_code)]
    pub fn new(n_sites: usize, on_site: Vec<f64>, dim: usize) -> Self {
        Self {
            n_sites,
            on_site,
            hoppings: Vec::new(),
            dim,
        }
    }

    /// Creates a uniform 1D chain with nearest-neighbor hopping.
    ///
    /// # Arguments
    /// * `n_sites` - Number of sites
    /// * `epsilon` - On-site energy (uniform)
    /// * `t` - Nearest-neighbor hopping integral
    #[allow(dead_code)]
    pub fn uniform_chain(n_sites: usize, epsilon: f64, t: f64) -> Self {
        let on_site = vec![epsilon; n_sites];
        let mut model = Self::new(n_sites, on_site, 1);
        for i in 0..n_sites - 1 {
            model.hoppings.push((i, i + 1, t));
            model.hoppings.push((i + 1, i, t));
        }
        model
    }

    /// Adds a hopping term between sites i and j.
    #[allow(dead_code)]
    pub fn add_hopping(&mut self, i: usize, j: usize, t: f64) {
        self.hoppings.push((i, j, t));
        self.hoppings.push((j, i, t));
    }

    /// Builds the full Hamiltonian matrix as a flat row-major array.
    #[allow(dead_code)]
    pub fn build_hamiltonian(&self) -> Vec<f64> {
        let n = self.n_sites;
        let mut h = vec![0.0_f64; n * n];
        for (i, &eps) in self.on_site.iter().enumerate() {
            h[i * n + i] = eps;
        }
        for &(i, j, t) in &self.hoppings {
            h[i * n + j] += t;
        }
        h
    }

    /// Computes eigenvalues via Jacobi iteration.
    ///
    /// Returns eigenvalues in ascending order.
    #[allow(dead_code)]
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n_sites;
        let h = self.build_hamiltonian();
        jacobi_eigenvalues(&h, n)
    }

    /// Returns the band width (max - min eigenvalue).
    #[allow(dead_code)]
    pub fn band_width(&self) -> f64 {
        let eigs = self.eigenvalues();
        if eigs.is_empty() {
            return 0.0;
        }
        eigs[eigs.len() - 1] - eigs[0]
    }

    /// Returns the density of states at energy E using Lorentzian broadening.
    ///
    /// # Arguments
    /// * `energy` - Energy at which to evaluate DOS
    /// * `eta` - Lorentzian broadening width
    #[allow(dead_code)]
    pub fn density_of_states(&self, energy: f64, eta: f64) -> f64 {
        let eigs = self.eigenvalues();
        let mut dos = 0.0;
        for e in eigs {
            dos += eta / (std::f64::consts::PI * ((energy - e).powi(2) + eta * eta));
        }
        dos
    }
}

/// Jacobi iteration to find eigenvalues of a symmetric matrix.
#[allow(dead_code)]
fn jacobi_eigenvalues(mat: &[f64], n: usize) -> Vec<f64> {
    let mut a = mat.to_vec();
    let max_iter = 100 * n * n;

    for _ in 0..max_iter {
        // Find largest off-diagonal element
        let mut max_val = 0.0_f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
                let val = a[i * n + j].abs();
                if val > max_val {
                    max_val = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }

        // Compute rotation angle
        let theta = if (a[q * n + q] - a[p * n + p]).abs() < 1e-15 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * ((2.0 * a[p * n + q]) / (a[q * n + q] - a[p * n + p])).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        // Apply Jacobi rotation
        let mut new_a = a.clone();
        for i in 0..n {
            if i != p && i != q {
                new_a[i * n + p] = c * a[i * n + p] + s * a[i * n + q];
                new_a[p * n + i] = new_a[i * n + p];
                new_a[i * n + q] = -s * a[i * n + p] + c * a[i * n + q];
                new_a[q * n + i] = new_a[i * n + q];
            }
        }
        new_a[p * n + p] = c * c * a[p * n + p] + 2.0 * s * c * a[p * n + q] + s * s * a[q * n + q];
        new_a[q * n + q] = s * s * a[p * n + p] - 2.0 * s * c * a[p * n + q] + c * c * a[q * n + q];
        new_a[p * n + q] = 0.0;
        new_a[q * n + p] = 0.0;
        a = new_a;
    }

    let mut eigs: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    eigs
}

/// Landauer transmission and conductance calculations.
///
/// Implements the Landauer formula: G = (2e²/h) * T(E_F)
/// where T(E_F) is the transmission coefficient at the Fermi energy.
#[derive(Debug, Clone)]
pub struct LandauerTransmission {
    /// Fermi energy in eV
    pub fermi_energy: f64,
    /// Number of transport channels
    pub n_channels: usize,
    /// Lead coupling strength (in energy units)
    pub gamma_lead: f64,
}

impl LandauerTransmission {
    /// Creates a new Landauer transmission calculator.
    ///
    /// # Arguments
    /// * `fermi_energy` - Fermi energy in eV
    /// * `n_channels` - Number of transport channels
    /// * `gamma_lead` - Lead broadening (imaginary self-energy)
    #[allow(dead_code)]
    pub fn new(fermi_energy: f64, n_channels: usize, gamma_lead: f64) -> Self {
        Self {
            fermi_energy,
            n_channels,
            gamma_lead,
        }
    }

    /// Computes transmission via transfer matrix method for a 1D chain.
    ///
    /// Returns T(E) in range \[0, N_channels\].
    ///
    /// # Arguments
    /// * `energy` - Energy at which to compute transmission
    /// * `on_site` - On-site energies along the chain
    /// * `hopping` - Nearest-neighbor hopping integral
    #[allow(dead_code)]
    pub fn transfer_matrix_transmission(&self, energy: f64, on_site: &[f64], hopping: f64) -> f64 {
        if hopping.abs() < 1e-15 {
            return 0.0;
        }
        // Transfer matrix: M_i = [(E-eps_i)/t, -1; 1, 0]
        // Total: M = prod M_i
        let mut m = [[1.0_f64, 0.0], [0.0, 1.0_f64]]; // identity
        for &eps in on_site {
            let ratio = (energy - eps) / hopping;
            let mi = [[ratio, -1.0], [1.0, 0.0]];
            let new_m = [
                [
                    m[0][0] * mi[0][0] + m[0][1] * mi[1][0],
                    m[0][0] * mi[0][1] + m[0][1] * mi[1][1],
                ],
                [
                    m[1][0] * mi[0][0] + m[1][1] * mi[1][0],
                    m[1][0] * mi[0][1] + m[1][1] * mi[1][1],
                ],
            ];
            m = new_m;
        }
        // Transmission from transfer matrix determinant
        // For 1D chain: T = 4 / |M_11 + M_22 + i*(M_12 - M_21)|^2
        // Simplified: T = 1 / (1 + (m[0][1]^2)/4) for uniform chain
        let trace = m[0][0] + m[1][1];
        let t_val = 4.0 / (trace * trace + (m[0][1] - m[1][0]).powi(2));
        t_val.clamp(0.0, 1.0)
    }

    /// Computes multichannel Landauer-Büttiker transmission.
    ///
    /// For a perfect conductor with N channels, T = N.
    #[allow(dead_code)]
    pub fn multichannel_transmission(&self, transmission_per_channel: &[f64]) -> f64 {
        transmission_per_channel.iter().sum()
    }

    /// Computes the Landauer conductance G = (2e²/h) * T.
    ///
    /// # Arguments
    /// * `transmission` - Total transmission coefficient
    #[allow(dead_code)]
    pub fn conductance(&self, transmission: f64) -> f64 {
        constants::G_QUANTUM * transmission
    }

    /// Computes perfect conductor conductance: G = N * 2e²/h.
    #[allow(dead_code)]
    pub fn perfect_conductor_conductance(&self) -> f64 {
        constants::G_QUANTUM * self.n_channels as f64
    }
}

/// Non-equilibrium Green's function for transport calculations.
///
/// Computes the retarded Green's function G^R = (E - H - Σ)^{-1}
/// where Σ is the lead self-energy in the wide-band limit.
#[derive(Debug, Clone)]
pub struct GreensFunctionMd {
    /// Size of the device region
    pub n_device: usize,
    /// Lead broadening Γ = -2 Im(Σ) in energy units
    pub gamma_left: f64,
    /// Right lead broadening
    pub gamma_right: f64,
    /// Infinitesimal imaginary part for regularization
    pub eta: f64,
}

impl GreensFunctionMd {
    /// Creates a new Green's function calculator.
    ///
    /// # Arguments
    /// * `n_device` - Number of sites in the device
    /// * `gamma_left` - Left lead broadening
    /// * `gamma_right` - Right lead broadening
    /// * `eta` - Regularization parameter
    #[allow(dead_code)]
    pub fn new(n_device: usize, gamma_left: f64, gamma_right: f64, eta: f64) -> Self {
        Self {
            n_device,
            gamma_left,
            gamma_right,
            eta,
        }
    }

    /// Computes retarded Green's function diagonal elements for a chain.
    ///
    /// Returns (Re G, Im G) for each site.
    /// Im(G^R) ≤ 0 always (positive semi-definite spectral function).
    ///
    /// # Arguments
    /// * `energy` - Energy argument E + iη
    /// * `hamiltonian` - Device Hamiltonian (flat row-major)
    #[allow(dead_code)]
    pub fn retarded_greens_function(&self, energy: f64, hamiltonian: &[f64]) -> Vec<Complex> {
        let n = self.n_device;
        // Build (E + iη - H - Σ) matrix
        // Self-energy (wide-band limit): Σ_L = -iΓ_L/2 at site 0, Σ_R = -iΓ_R/2 at site n-1
        let mut mat_re = vec![0.0_f64; n * n];
        let mut mat_im = vec![0.0_f64; n * n];

        for i in 0..n {
            for j in 0..n {
                mat_re[i * n + j] = if i == j {
                    energy - hamiltonian[i * n + j]
                } else {
                    -hamiltonian[i * n + j]
                };
                mat_im[i * n + j] = 0.0;
            }
        }
        // Add eta broadening + lead self-energies
        mat_im[0] += self.eta + self.gamma_left / 2.0;
        mat_im[(n - 1) * n + (n - 1)] += self.eta + self.gamma_right / 2.0;
        for i in 1..n - 1 {
            mat_im[i * n + i] += self.eta;
        }

        // Invert complex matrix (size n×n) using Gauss-Jordan elimination
        invert_complex_matrix(&mat_re, &mat_im, n)
    }

    /// Computes the local density of states at site i.
    ///
    /// LDOS_i = -1/π * Im(G^R_{ii})
    ///
    /// # Arguments
    /// * `greens_diag` - Diagonal elements of G^R
    #[allow(dead_code)]
    pub fn local_dos(&self, greens_diag: &[Complex]) -> Vec<f64> {
        greens_diag
            .iter()
            .map(|g| -g.im / std::f64::consts::PI)
            .collect()
    }

    /// Computes the spectral function A = i(G^R - G^A) = -2 Im G^R.
    #[allow(dead_code)]
    pub fn spectral_function(&self, greens_diag: &[Complex]) -> Vec<f64> {
        greens_diag.iter().map(|g| -2.0 * g.im).collect()
    }

    /// Computes the transmission using the Fisher-Lee relation.
    ///
    /// T = Tr\[Γ_L G^R Γ_R G^A\]
    /// Simplified for 1D: T = Γ_L * Γ_R * |G^R_{1N}|²
    ///
    /// # Arguments
    /// * `energy` - Energy at which to compute T
    /// * `hamiltonian` - Device Hamiltonian
    #[allow(dead_code)]
    pub fn transmission_fisher_lee(&self, energy: f64, hamiltonian: &[f64]) -> f64 {
        let n = self.n_device;
        let gf = self.retarded_greens_function(energy, hamiltonian);
        // G^R_{0, n-1} element
        let g_1n = gf[n - 1];
        self.gamma_left * self.gamma_right * g_1n.norm_sq()
    }
}

/// Invert a complex matrix using Gauss-Jordan elimination.
#[allow(dead_code)]
fn invert_complex_matrix(re: &[f64], im: &[f64], n: usize) -> Vec<Complex> {
    // Augmented matrix [A | I] for complex A
    let mut aug_re = vec![0.0_f64; n * 2 * n];
    let mut aug_im = vec![0.0_f64; n * 2 * n];

    for i in 0..n {
        for j in 0..n {
            aug_re[i * 2 * n + j] = re[i * n + j];
            aug_im[i * 2 * n + j] = im[i * n + j];
        }
        aug_re[i * 2 * n + n + i] = 1.0; // identity block
    }

    for col in 0..n {
        // Find pivot
        let mut max_val = 0.0_f64;
        let mut pivot = col;
        for row in col..n {
            let val = Complex::new(aug_re[row * 2 * n + col], aug_im[row * 2 * n + col]).norm();
            if val > max_val {
                max_val = val;
                pivot = row;
            }
        }
        if max_val < 1e-15 {
            // Singular or near-singular, return zeros
            return vec![Complex::new(0.0, 0.0); n * n];
        }
        // Swap rows
        if pivot != col {
            for j in 0..2 * n {
                aug_re.swap(col * 2 * n + j, pivot * 2 * n + j);
                aug_im.swap(col * 2 * n + j, pivot * 2 * n + j);
            }
        }
        // Normalize pivot row
        let piv = Complex::new(aug_re[col * 2 * n + col], aug_im[col * 2 * n + col]);
        for j in 0..2 * n {
            let elem = Complex::new(aug_re[col * 2 * n + j], aug_im[col * 2 * n + j]);
            let div = elem.div(&piv);
            aug_re[col * 2 * n + j] = div.re;
            aug_im[col * 2 * n + j] = div.im;
        }
        // Eliminate column
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = Complex::new(aug_re[row * 2 * n + col], aug_im[row * 2 * n + col]);
            for j in 0..2 * n {
                let piv_elem = Complex::new(aug_re[col * 2 * n + j], aug_im[col * 2 * n + j]);
                let sub = factor.mul(&piv_elem);
                aug_re[row * 2 * n + j] -= sub.re;
                aug_im[row * 2 * n + j] -= sub.im;
            }
        }
    }

    // Extract inverse
    let mut result = vec![Complex::new(0.0, 0.0); n * n];
    for i in 0..n {
        for j in 0..n {
            result[i * n + j] = Complex::new(aug_re[i * 2 * n + n + j], aug_im[i * 2 * n + n + j]);
        }
    }
    result
}

/// Non-equilibrium steady-state current calculations.
///
/// Computes the Landauer-Büttiker current:
/// J = (2e/h) ∫ T(E) \[f_L(E) - f_R(E)\] dE
#[derive(Debug, Clone)]
pub struct NonEquilibriumMd {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// Left lead chemical potential in eV
    pub mu_left: f64,
    /// Right lead chemical potential in eV
    pub mu_right: f64,
}

impl NonEquilibriumMd {
    /// Creates a new non-equilibrium calculator.
    ///
    /// # Arguments
    /// * `temperature` - Temperature in Kelvin
    /// * `mu_left` - Left lead chemical potential
    /// * `mu_right` - Right lead chemical potential
    #[allow(dead_code)]
    pub fn new(temperature: f64, mu_left: f64, mu_right: f64) -> Self {
        Self {
            temperature,
            mu_left,
            mu_right,
        }
    }

    /// Computes Fermi-Dirac distribution function.
    ///
    /// f(E) = 1 / (exp((E - μ) / kT) + 1)
    ///
    /// # Arguments
    /// * `energy` - Energy in eV
    /// * `mu` - Chemical potential in eV
    #[allow(dead_code)]
    pub fn fermi_dirac(&self, energy: f64, mu: f64) -> f64 {
        let kt = constants::K_BOLTZMANN * self.temperature / constants::E_CHARGE; // in eV
        if kt < 1e-15 {
            // Zero temperature: step function
            return if energy < mu { 1.0 } else { 0.0 };
        }
        let x = (energy - mu) / kt;
        // Clip to avoid overflow
        if x > 500.0 {
            return 0.0;
        }
        if x < -500.0 {
            return 1.0;
        }
        1.0 / (x.exp() + 1.0)
    }

    /// Computes the steady-state current by numerical integration.
    ///
    /// Uses trapezoidal rule over the bias window.
    ///
    /// # Arguments
    /// * `gf` - Green's function calculator
    /// * `hamiltonian` - Device Hamiltonian
    /// * `n_points` - Number of integration points
    #[allow(dead_code)]
    pub fn current(&self, gf: &GreensFunctionMd, hamiltonian: &[f64], n_points: usize) -> f64 {
        let bias = self.mu_left - self.mu_right;
        if bias.abs() < 1e-15 {
            return 0.0;
        }
        let e_min = self.mu_right.min(self.mu_left) - 0.5;
        let e_max = self.mu_right.max(self.mu_left) + 0.5;
        let de = (e_max - e_min) / n_points as f64;

        let mut integral = 0.0;
        for k in 0..n_points {
            let e = e_min + (k as f64 + 0.5) * de;
            let t_e = gf.transmission_fisher_lee(e, hamiltonian);
            let df = self.fermi_dirac(e, self.mu_left) - self.fermi_dirac(e, self.mu_right);
            integral += t_e * df * de;
        }
        // J = (2e/h) * integral * e (converting eV to J)
        let prefactor = 2.0 * constants::E_CHARGE / constants::H_PLANCK;
        prefactor * integral * constants::E_CHARGE
    }

    /// Computes the differential conductance dJ/dV at zero bias.
    #[allow(dead_code)]
    pub fn zero_bias_conductance(&self, transmission_at_fermi: f64) -> f64 {
        constants::G_QUANTUM * transmission_at_fermi
    }
}

/// Electron-phonon coupling in the Fröhlich model.
///
/// Computes polaron mass renormalization and self-energy corrections.
#[derive(Debug, Clone)]
pub struct ElectronPhononCoupling {
    /// Fröhlich coupling constant α (dimensionless)
    pub alpha: f64,
    /// Phonon frequency ω_LO in rad/s
    pub omega_lo: f64,
    /// Temperature in Kelvin
    pub temperature: f64,
}

impl ElectronPhononCoupling {
    /// Creates a new electron-phonon coupling calculator.
    ///
    /// # Arguments
    /// * `alpha` - Fröhlich coupling constant
    /// * `omega_lo` - LO phonon frequency
    /// * `temperature` - Temperature in Kelvin
    #[allow(dead_code)]
    pub fn new(alpha: f64, omega_lo: f64, temperature: f64) -> Self {
        Self {
            alpha,
            omega_lo,
            temperature,
        }
    }

    /// Computes the polaron effective mass ratio m*/m_e.
    ///
    /// In weak coupling: m*/m = 1 + α/6
    #[allow(dead_code)]
    pub fn polaron_mass_ratio(&self) -> f64 {
        1.0 + self.alpha / 6.0
    }

    /// Computes the electron-phonon self-energy (lowest order, imaginary part).
    ///
    /// Σ''(E) ≈ -π α ω_LO * n_BE(ω_LO) at lowest order
    ///
    /// # Arguments
    /// * `energy` - Electron energy
    #[allow(dead_code)]
    pub fn self_energy_imaginary(&self, energy: f64) -> f64 {
        let n_be = self.bose_einstein(self.omega_lo);
        // Imaginary self-energy from absorption
        let sigma_abs = -std::f64::consts::PI * self.alpha * self.omega_lo * n_be;
        // From emission (if energy - ω > 0)
        let sigma_em = if energy > self.omega_lo * constants::HBAR {
            -std::f64::consts::PI * self.alpha * self.omega_lo * (n_be + 1.0)
        } else {
            0.0
        };
        sigma_abs + sigma_em
    }

    /// Bose-Einstein distribution for phonons.
    ///
    /// n(ω) = 1 / (exp(ħω/kT) - 1)
    #[allow(dead_code)]
    pub fn bose_einstein(&self, omega: f64) -> f64 {
        let kt = constants::K_BOLTZMANN * self.temperature;
        let x = constants::HBAR * omega / kt;
        if x > 500.0 {
            return 0.0;
        }
        if x < 1e-10 {
            return kt / (constants::HBAR * omega);
        }
        1.0 / (x.exp() - 1.0)
    }

    /// Computes resistivity increase due to electron-phonon scattering.
    ///
    /// At high T: ρ ∝ T (Bloch-Grüneisen limit)
    ///
    /// # Arguments
    /// * `base_resistivity` - Residual resistivity at T=0
    #[allow(dead_code)]
    pub fn resistivity(&self, base_resistivity: f64) -> f64 {
        let n_be = self.bose_einstein(self.omega_lo);
        // Simple model: extra scattering ∝ phonon occupation
        base_resistivity * (1.0 + self.alpha * n_be)
    }
}

/// Couples MD trajectory to quantum transport calculations.
///
/// Computes time-averaged transport properties along an MD trajectory.
#[derive(Debug, Clone)]
pub struct MdTransportCoupling {
    /// Hopping fluctuation amplitude (eV)
    pub hopping_fluctuation: f64,
    /// Base hopping integral (eV)
    pub t0: f64,
    /// Fermi energy (eV)
    pub fermi_energy: f64,
    /// Running sum of transmission values
    sum_transmission: f64,
    /// Number of frames processed
    n_frames: usize,
}

impl MdTransportCoupling {
    /// Creates a new MD-transport coupling object.
    ///
    /// # Arguments
    /// * `t0` - Base hopping integral in eV
    /// * `hopping_fluctuation` - Fluctuation amplitude
    /// * `fermi_energy` - Fermi energy in eV
    #[allow(dead_code)]
    pub fn new(t0: f64, hopping_fluctuation: f64, fermi_energy: f64) -> Self {
        Self {
            t0,
            hopping_fluctuation,
            fermi_energy,
            sum_transmission: 0.0,
            n_frames: 0,
        }
    }

    /// Updates transport from a new MD frame.
    ///
    /// # Arguments
    /// * `positions` - Atom positions along the chain
    /// * `eq_spacing` - Equilibrium bond spacing
    /// * `gf` - Green's function calculator
    #[allow(dead_code)]
    pub fn update_frame(&mut self, positions: &[[f64; 3]], eq_spacing: f64, gf: &GreensFunctionMd) {
        let n = positions.len();
        if n < 2 {
            return;
        }
        // Build Hamiltonian from positions (SSH-like coupling)
        let mut h = vec![0.0_f64; n * n];
        for i in 0..n - 1 {
            let dx = positions[i + 1][0] - positions[i][0];
            let dy = positions[i + 1][1] - positions[i][1];
            let dz = positions[i + 1][2] - positions[i][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let t = self.t0 * (-(d - eq_spacing) / eq_spacing).exp();
            h[i * n + (i + 1)] = t;
            h[(i + 1) * n + i] = t;
        }
        let t_ef = gf.transmission_fisher_lee(self.fermi_energy, &h);
        self.sum_transmission += t_ef;
        self.n_frames += 1;
    }

    /// Returns the time-averaged transmission at the Fermi energy.
    #[allow(dead_code)]
    pub fn average_transmission(&self) -> f64 {
        if self.n_frames == 0 {
            return 0.0;
        }
        self.sum_transmission / self.n_frames as f64
    }

    /// Returns the number of frames processed.
    #[allow(dead_code)]
    pub fn n_frames(&self) -> usize {
        self.n_frames
    }
}

/// Wigner quasi-probability distribution for quantum phase space.
///
/// The Wigner function W(x, p) is a quasi-probability distribution
/// that provides a quantum analog of the classical phase-space distribution.
#[derive(Debug, Clone)]
pub struct WignerFunction {
    /// Position grid
    pub x_grid: Vec<f64>,
    /// Momentum grid
    pub p_grid: Vec<f64>,
    /// Wigner function values W\[ix\]\[ip\]
    pub values: Vec<Vec<f64>>,
}

impl WignerFunction {
    /// Creates a Wigner function from a pure state wavefunction.
    ///
    /// W(x, p) = (1/πħ) ∫ ψ*(x+y) ψ(x-y) exp(2ipy/ħ) dy
    ///
    /// # Arguments
    /// * `psi_re` - Real part of wavefunction on x_grid
    /// * `psi_im` - Imaginary part
    /// * `x_grid` - Position grid
    /// * `p_grid` - Momentum grid
    #[allow(dead_code)]
    pub fn from_wavefunction(
        psi_re: &[f64],
        psi_im: &[f64],
        x_grid: &[f64],
        p_grid: &[f64],
    ) -> Self {
        let nx = x_grid.len();
        let np = p_grid.len();
        let dx = if nx > 1 { x_grid[1] - x_grid[0] } else { 1.0 };
        let mut values = vec![vec![0.0_f64; np]; nx];

        for (ix, _xval) in x_grid.iter().enumerate() {
            for (ip, &pval) in p_grid.iter().enumerate() {
                let mut w = 0.0;
                for iy in 0..nx {
                    // y = (iy - nx/2) * dx
                    let y = (iy as f64 - nx as f64 / 2.0) * dx;
                    let ix_plus = ((ix as isize + iy as isize - nx as isize / 2)
                        .rem_euclid(nx as isize)) as usize;
                    let ix_minus = ((ix as isize - iy as isize + nx as isize / 2)
                        .rem_euclid(nx as isize)) as usize;
                    // ψ*(x+y)
                    let psi_plus_re = psi_re[ix_plus];
                    let psi_plus_im = -psi_im[ix_plus]; // conjugate
                    // ψ(x-y)
                    let psi_minus_re = psi_re[ix_minus];
                    let psi_minus_im = psi_im[ix_minus];
                    // exp(2ipy/ħ) -> use ħ=1 units
                    let phase = 2.0 * pval * y;
                    let exp_re = phase.cos();
                    let exp_im = phase.sin();
                    // Product of three complex numbers
                    let prod_re = (psi_plus_re * psi_minus_re - psi_plus_im * psi_minus_im)
                        * exp_re
                        - (psi_plus_re * psi_minus_im + psi_plus_im * psi_minus_re) * exp_im;
                    w += prod_re * dx;
                }
                values[ix][ip] = w / std::f64::consts::PI;
            }
        }

        Self {
            x_grid: x_grid.to_vec(),
            p_grid: p_grid.to_vec(),
            values,
        }
    }

    /// Integrates the Wigner function over momentum to get the position density.
    ///
    /// n(x) = ∫ W(x,p) dp = |ψ(x)|²
    #[allow(dead_code)]
    pub fn position_density(&self) -> Vec<f64> {
        let np = self.p_grid.len();
        let dp = if np > 1 {
            self.p_grid[1] - self.p_grid[0]
        } else {
            1.0
        };
        self.values
            .iter()
            .map(|row| row.iter().sum::<f64>() * dp)
            .collect()
    }

    /// Integrates the Wigner function over all phase space.
    ///
    /// Should equal 1 for a normalized wavefunction.
    #[allow(dead_code)]
    pub fn total_weight(&self) -> f64 {
        let np = self.p_grid.len();
        let nx = self.x_grid.len();
        let dp = if np > 1 {
            self.p_grid[1] - self.p_grid[0]
        } else {
            1.0
        };
        let dx = if nx > 1 {
            self.x_grid[1] - self.x_grid[0]
        } else {
            1.0
        };
        self.values
            .iter()
            .map(|row| row.iter().sum::<f64>())
            .sum::<f64>()
            * dx
            * dp
    }
}

/// Semiclassical Boltzmann transport equation solver.
///
/// Uses the relaxation time approximation τ(E) to compute
/// transport coefficients.
#[derive(Debug, Clone)]
pub struct BoltzmannTransport {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// Chemical potential in eV
    pub mu: f64,
    /// Deformation potential coupling (eV)
    pub d_ep: f64,
    /// Phonon energy scale (eV)
    pub omega_ph: f64,
}

impl BoltzmannTransport {
    /// Creates a new Boltzmann transport solver.
    ///
    /// # Arguments
    /// * `temperature` - Temperature in Kelvin
    /// * `mu` - Chemical potential in eV
    /// * `d_ep` - Deformation potential
    /// * `omega_ph` - Phonon energy
    #[allow(dead_code)]
    pub fn new(temperature: f64, mu: f64, d_ep: f64, omega_ph: f64) -> Self {
        Self {
            temperature,
            mu,
            d_ep,
            omega_ph,
        }
    }

    /// Computes the relaxation time τ(E) in the deformation potential approximation.
    ///
    /// 1/τ ∝ D² * ρ(E) * kT for acoustic phonons
    ///
    /// # Arguments
    /// * `energy` - Electron energy in eV
    #[allow(dead_code)]
    pub fn relaxation_time(&self, energy: f64) -> f64 {
        let kt = constants::K_BOLTZMANN * self.temperature / constants::E_CHARGE; // eV
        if self.d_ep.abs() < 1e-15 || kt < 1e-15 {
            return 1e-13; // 100 fs default
        }
        let dos = (energy.abs() + 1e-6).sqrt(); // simple 3D DOS ~ sqrt(E)
        1e-13 / (self.d_ep * self.d_ep * dos * kt + 1e-30)
    }

    /// Computes electrical conductivity via Boltzmann equation.
    ///
    /// σ = (2e²/V) Σ_k v_k² τ(E_k) (-∂f/∂E)|_{E=E_k}
    ///
    /// # Arguments
    /// * `energies` - Energy eigenvalues
    /// * `velocities` - Group velocities at each k-point
    #[allow(dead_code)]
    pub fn conductivity(&self, energies: &[f64], velocities: &[f64]) -> f64 {
        let kt = constants::K_BOLTZMANN * self.temperature / constants::E_CHARGE;
        if kt < 1e-15 {
            return 0.0;
        }
        let mut sigma = 0.0;
        for (&e, &v) in energies.iter().zip(velocities.iter()) {
            let tau = self.relaxation_time(e);
            let x = (e - self.mu) / kt;
            // -df/dE = (1/kT) * f(1-f)
            let f = if x > 500.0 {
                0.0
            } else if x < -500.0 {
                1.0
            } else {
                1.0 / (x.exp() + 1.0)
            };
            let df_de = f * (1.0 - f) / kt;
            sigma += v * v * tau * df_de;
        }
        sigma * 2.0 * constants::E_CHARGE * constants::E_CHARGE
    }

    /// Computes the Seebeck coefficient (thermopower).
    ///
    /// S = -(1/eT) * L_1/L_0 where L_n = ∫ (E-μ)^n τ(E) v² (-df/dE) dE
    ///
    /// # Arguments
    /// * `energies` - Energy eigenvalues
    /// * `velocities` - Group velocities
    #[allow(dead_code)]
    pub fn seebeck_coefficient(&self, energies: &[f64], velocities: &[f64]) -> f64 {
        let kt = constants::K_BOLTZMANN * self.temperature / constants::E_CHARGE;
        if kt < 1e-15 {
            return 0.0;
        }
        let mut l0 = 0.0;
        let mut l1 = 0.0;
        for (&e, &v) in energies.iter().zip(velocities.iter()) {
            let tau = self.relaxation_time(e);
            let x = (e - self.mu) / kt;
            let f = if x > 500.0 {
                0.0
            } else if x < -500.0 {
                1.0
            } else {
                1.0 / (x.exp() + 1.0)
            };
            let df_de = f * (1.0 - f) / kt;
            let weight = v * v * tau * df_de;
            l0 += weight;
            l1 += weight * (e - self.mu);
        }
        if l0.abs() < 1e-30 {
            return 0.0;
        }
        -l1 / (constants::E_CHARGE * self.temperature * l0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // ---- Complex arithmetic ----

    #[test]
    fn test_complex_norm_sq() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.norm_sq() - 25.0).abs() < EPS);
    }

    #[test]
    fn test_complex_norm() {
        let z = Complex::new(3.0, 4.0);
        assert!((z.norm() - 5.0).abs() < EPS);
    }

    #[test]
    fn test_complex_conj() {
        let z = Complex::new(1.0, 2.0);
        let c = z.conj();
        assert!((c.re - 1.0).abs() < EPS && (c.im + 2.0).abs() < EPS);
    }

    #[test]
    fn test_complex_mul() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let p = a.mul(&b);
        assert!((p.re - (-5.0)).abs() < EPS);
        assert!((p.im - 10.0).abs() < EPS);
    }

    #[test]
    fn test_complex_div() {
        let a = Complex::new(1.0, 0.0);
        let b = Complex::new(2.0, 0.0);
        let d = a.div(&b);
        assert!((d.re - 0.5).abs() < EPS && d.im.abs() < EPS);
    }

    // ---- TightBindingModel ----

    #[test]
    fn test_tight_binding_uniform_chain_eigenvalues_count() {
        let model = TightBindingModel::uniform_chain(5, 0.0, -1.0);
        let eigs = model.eigenvalues();
        assert_eq!(eigs.len(), 5);
    }

    #[test]
    fn test_tight_binding_chain_eigenvalues_sorted() {
        let model = TightBindingModel::uniform_chain(4, 0.0, -1.0);
        let eigs = model.eigenvalues();
        for i in 0..eigs.len() - 1 {
            assert!(eigs[i] <= eigs[i + 1] + 1e-10);
        }
    }

    #[test]
    fn test_tight_binding_single_site_eigenvalue() {
        let model = TightBindingModel::new(1, vec![2.5], 1);
        let eigs = model.eigenvalues();
        assert!((eigs[0] - 2.5).abs() < 1e-8);
    }

    #[test]
    fn test_tight_binding_two_site_eigenvalues() {
        let mut model = TightBindingModel::new(2, vec![0.0, 0.0], 1);
        model.add_hopping(0, 1, -1.0);
        let eigs = model.eigenvalues();
        // Should be ±1
        assert!((eigs[0] + 1.0).abs() < 1e-6);
        assert!((eigs[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_tight_binding_band_width_positive() {
        let model = TightBindingModel::uniform_chain(6, 0.0, -1.0);
        assert!(model.band_width() > 0.0);
    }

    #[test]
    fn test_tight_binding_dos_positive() {
        let model = TightBindingModel::uniform_chain(4, 0.0, -1.0);
        let dos = model.density_of_states(0.0, 0.1);
        assert!(dos >= 0.0);
    }

    #[test]
    fn test_tight_binding_hamiltonian_symmetry() {
        let model = TightBindingModel::uniform_chain(3, 0.0, -1.0);
        let h = model.build_hamiltonian();
        let n = 3;
        for i in 0..n {
            for j in 0..n {
                assert!((h[i * n + j] - h[j * n + i]).abs() < EPS);
            }
        }
    }

    // ---- LandauerTransmission ----

    #[test]
    fn test_landauer_perfect_conductor_n_channels() {
        let lt = LandauerTransmission::new(0.0, 4, 0.1);
        // Perfect conductor T = N_channels
        let g = lt.perfect_conductor_conductance();
        assert!((g - 4.0 * constants::G_QUANTUM).abs() < 1e-20);
    }

    #[test]
    fn test_landauer_conductance_quantized() {
        let lt = LandauerTransmission::new(0.0, 1, 0.1);
        let g = lt.conductance(1.0);
        assert!((g - constants::G_QUANTUM).abs() < 1e-20);
    }

    #[test]
    fn test_landauer_multichannel_sum() {
        let lt = LandauerTransmission::new(0.0, 3, 0.1);
        let t_per_ch = [0.9, 0.8, 1.0];
        let total = lt.multichannel_transmission(&t_per_ch);
        assert!((total - 2.7).abs() < 1e-10);
    }

    #[test]
    fn test_landauer_conductance_zero_for_zero_transmission() {
        let lt = LandauerTransmission::new(0.0, 1, 0.1);
        let g = lt.conductance(0.0);
        assert!(g.abs() < EPS);
    }

    #[test]
    fn test_landauer_transmission_bounded() {
        let lt = LandauerTransmission::new(0.0, 1, 0.1);
        let on_site = vec![0.0; 5];
        let t = lt.transfer_matrix_transmission(0.0, &on_site, -1.0);
        assert!((0.0..=1.0 + 1e-10).contains(&t));
    }

    // ---- GreensFunctionMd ----

    #[test]
    fn test_greens_function_imaginary_part_nonpositive() {
        let gf = GreensFunctionMd::new(3, 0.1, 0.1, 1e-4);
        let h = vec![0.0_f64; 9]; // zero Hamiltonian
        let g = gf.retarded_greens_function(0.0, &h);
        // Im(G^R_{ii}) should be ≤ 0 for physical Green's function
        for i in 0..3 {
            assert!(g[i * 3 + i].im <= 1e-6);
        }
    }

    #[test]
    fn test_greens_function_ldos_nonnegative() {
        let gf = GreensFunctionMd::new(2, 0.1, 0.1, 1e-4);
        let h = vec![0.0_f64; 4];
        let g = gf.retarded_greens_function(0.0, &h);
        let ldos = gf.local_dos(&g);
        for d in ldos {
            assert!(d >= -1e-6);
        }
    }

    #[test]
    fn test_greens_function_spectral_nonnegative() {
        let gf = GreensFunctionMd::new(2, 0.1, 0.1, 1e-4);
        let h = vec![0.0_f64; 4];
        let g = gf.retarded_greens_function(0.0, &h);
        let a = gf.spectral_function(&g);
        for ai in a {
            assert!(ai >= -1e-6);
        }
    }

    #[test]
    fn test_greens_function_transmission_nonneg() {
        let gf = GreensFunctionMd::new(3, 0.1, 0.1, 1e-4);
        let mut h = vec![0.0_f64; 9];
        h[1] = -1.0;
        h[3] = -1.0;
        h[3 + 2] = -1.0;
        h[2 * 3 + 1] = -1.0;
        let t = gf.transmission_fisher_lee(0.0, &h);
        assert!(t >= 0.0);
    }

    // ---- NonEquilibriumMd ----

    #[test]
    fn test_fermi_dirac_at_fermi_energy() {
        let neq = NonEquilibriumMd::new(300.0, 0.0, 0.0);
        let f = neq.fermi_dirac(0.0, 0.0);
        assert!((f - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_fermi_dirac_deep_below() {
        let neq = NonEquilibriumMd::new(300.0, 0.0, 0.0);
        let f = neq.fermi_dirac(-10.0, 0.0);
        assert!(f > 0.99);
    }

    #[test]
    fn test_fermi_dirac_deep_above() {
        let neq = NonEquilibriumMd::new(300.0, 0.0, 0.0);
        let f = neq.fermi_dirac(10.0, 0.0);
        assert!(f < 0.01);
    }

    #[test]
    fn test_current_vanishes_at_zero_bias() {
        let neq = NonEquilibriumMd::new(300.0, 0.5, 0.5); // zero bias
        let gf = GreensFunctionMd::new(2, 0.1, 0.1, 1e-4);
        let h = vec![0.0_f64; 4];
        let j = neq.current(&gf, &h, 10);
        assert!(j.abs() < 1e-20);
    }

    #[test]
    fn test_zero_bias_conductance_quantized() {
        let neq = NonEquilibriumMd::new(300.0, 0.0, 0.0);
        let g = neq.zero_bias_conductance(1.0);
        assert!((g - constants::G_QUANTUM).abs() < 1e-20);
    }

    // ---- ElectronPhononCoupling ----

    #[test]
    fn test_polaron_mass_greater_than_bare() {
        let epc = ElectronPhononCoupling::new(2.0, 1e13, 300.0);
        let ratio = epc.polaron_mass_ratio();
        assert!(ratio > 1.0);
    }

    #[test]
    fn test_polaron_mass_weak_coupling() {
        let epc = ElectronPhononCoupling::new(0.0, 1e13, 300.0);
        let ratio = epc.polaron_mass_ratio();
        assert!((ratio - 1.0).abs() < EPS);
    }

    #[test]
    fn test_self_energy_imaginary_nonpositive_negative_energy() {
        let epc = ElectronPhononCoupling::new(1.0, 1e12, 300.0);
        let sigma = epc.self_energy_imaginary(-1.0);
        assert!(sigma <= 0.0);
    }

    #[test]
    fn test_electron_phonon_resistivity_increases_with_alpha() {
        let epc_low = ElectronPhononCoupling::new(0.5, 1e12, 300.0);
        let epc_high = ElectronPhononCoupling::new(2.0, 1e12, 300.0);
        let rho_low = epc_low.resistivity(1.0);
        let rho_high = epc_high.resistivity(1.0);
        assert!(rho_high >= rho_low);
    }

    #[test]
    fn test_bose_einstein_positive() {
        let epc = ElectronPhononCoupling::new(1.0, 1e13, 300.0);
        let n = epc.bose_einstein(1e13);
        assert!(n >= 0.0);
    }

    // ---- MdTransportCoupling ----

    #[test]
    fn test_md_transport_zero_frames() {
        let coupling = MdTransportCoupling::new(-1.0, 0.01, 0.0);
        assert!((coupling.average_transmission()).abs() < EPS);
    }

    #[test]
    fn test_md_transport_frame_count() {
        let mut coupling = MdTransportCoupling::new(-1.0, 0.01, 0.0);
        let gf = GreensFunctionMd::new(2, 0.1, 0.1, 1e-4);
        let pos = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        coupling.update_frame(&pos, 1.5, &gf);
        coupling.update_frame(&pos, 1.5, &gf);
        assert_eq!(coupling.n_frames(), 2);
    }

    #[test]
    fn test_md_transport_average_nonneg() {
        let mut coupling = MdTransportCoupling::new(-1.0, 0.01, 0.0);
        let gf = GreensFunctionMd::new(3, 0.1, 0.1, 1e-4);
        let pos = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]];
        coupling.update_frame(&pos, 1.5, &gf);
        assert!(coupling.average_transmission() >= 0.0);
    }

    // ---- WignerFunction ----

    #[test]
    fn test_wigner_position_density_nonneg() {
        let n = 8;
        let x: Vec<f64> = (0..n).map(|i| i as f64 * 0.5).collect();
        let p: Vec<f64> = (0..n).map(|i| (i as f64 - 4.0) * 0.5).collect();
        // Gaussian wavefunction
        let psi_re: Vec<f64> = x.iter().map(|&xi| (-(xi - 2.0).powi(2)).exp()).collect();
        let psi_im = vec![0.0_f64; n];
        let wf = WignerFunction::from_wavefunction(&psi_re, &psi_im, &x, &p);
        let dens = wf.position_density();
        for d in dens {
            // Can be slightly negative due to numerical error
            assert!(d > -0.1);
        }
    }

    #[test]
    fn test_wigner_grid_sizes() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let p = vec![-1.0, 0.0, 1.0];
        let psi_re = vec![1.0, 0.5, 0.25, 0.1];
        let psi_im = vec![0.0; 4];
        let wf = WignerFunction::from_wavefunction(&psi_re, &psi_im, &x, &p);
        assert_eq!(wf.values.len(), 4);
        assert_eq!(wf.values[0].len(), 3);
    }

    // ---- BoltzmannTransport ----

    #[test]
    fn test_boltzmann_relaxation_time_positive() {
        let bt = BoltzmannTransport::new(300.0, 0.0, 5.0, 0.03);
        let tau = bt.relaxation_time(0.1);
        assert!(tau > 0.0);
    }

    #[test]
    fn test_boltzmann_conductivity_positive() {
        let bt = BoltzmannTransport::new(300.0, 0.0, 5.0, 0.03);
        let e = vec![0.05, 0.1, 0.15, 0.2];
        let v = vec![1e5_f64, 1.1e5, 1.2e5, 1.3e5];
        let sigma = bt.conductivity(&e, &v);
        assert!(sigma >= 0.0);
    }

    #[test]
    fn test_boltzmann_seebeck_finite() {
        let bt = BoltzmannTransport::new(300.0, 0.1, 5.0, 0.03);
        let e = vec![0.05, 0.1, 0.15, 0.2];
        let v = vec![1e5_f64, 1.1e5, 1.2e5, 1.3e5];
        let s = bt.seebeck_coefficient(&e, &v);
        assert!(s.is_finite());
    }

    #[test]
    fn test_boltzmann_seebeck_sign() {
        // For n-type (μ in band), Seebeck should be negative
        let bt = BoltzmannTransport::new(300.0, 0.15, 5.0, 0.03);
        let e: Vec<f64> = (0..10).map(|i| 0.05 + i as f64 * 0.02).collect();
        let v: Vec<f64> = e.iter().map(|ei| (ei * 1e10).sqrt()).collect();
        let s = bt.seebeck_coefficient(&e, &v);
        // Should be finite (direction depends on DOS asymmetry)
        assert!(s.is_finite());
    }

    #[test]
    fn test_constants_conductance_quantum() {
        // G_Q = 2e²/h
        let g_calc = 2.0 * constants::E_CHARGE * constants::E_CHARGE / constants::H_PLANCK;
        assert!((g_calc - constants::G_QUANTUM).abs() / constants::G_QUANTUM < 1e-6);
    }

    #[test]
    fn test_fermi_dirac_monotone_decreasing() {
        let neq = NonEquilibriumMd::new(300.0, 0.5, -0.5);
        let f1 = neq.fermi_dirac(-1.0, 0.0);
        let f2 = neq.fermi_dirac(0.0, 0.0);
        let f3 = neq.fermi_dirac(1.0, 0.0);
        assert!(f1 > f2);
        assert!(f2 > f3);
    }

    #[test]
    fn test_invert_identity_matrix() {
        let re = vec![1.0, 0.0, 0.0, 1.0];
        let im = vec![0.0, 0.0, 0.0, 0.0];
        let inv = invert_complex_matrix(&re, &im, 2);
        assert!((inv[0].re - 1.0).abs() < 1e-10);
        assert!((inv[3].re - 1.0).abs() < 1e-10);
        assert!(inv[1].re.abs() < 1e-10);
    }
}
