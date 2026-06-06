// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Surface science models for MD simulations.
//!
//! Provides:
//! - Adsorption isotherms (Langmuir, Freundlich, BET)
//! - [`SurfaceAdsorption`] — Metropolis Monte Carlo site occupancy
//! - Surface energy / wetting (Young, Dupré, spreading coefficient)
//! - [`SurfaceDiffusion`] — 1D surface hopping via Arrhenius rates
//! - Debye–Waller factor and surface phonon frequency
//! - [`ThinFilmGrowth`] — ballistic deposition model
//! - Hertz contact mechanics and DMT pull-off force
//!
//! References:
//! - Adamson, A. W. & Gast, A. P. (1997). *Physical Chemistry of Surfaces*.
//! - Masel, R. I. (1996). *Principles of Adsorption and Reaction on Solid Surfaces*.
//! - Johnson, K. L. (1985). *Contact Mechanics*. Cambridge.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// AdsorptionSite
// ---------------------------------------------------------------------------

/// A single discrete adsorption site on a surface.
#[derive(Debug, Clone)]
pub struct AdsorptionSite {
    /// Position of the site in 3D space (m).
    pub pos: [f64; 3],
    /// Adsorption energy (negative = attractive, J).
    pub energy: f64,
    /// Whether the site is currently occupied by an adsorbate.
    pub occupied: bool,
}

impl AdsorptionSite {
    /// Create a new unoccupied adsorption site.
    pub fn new(pos: [f64; 3], energy: f64) -> Self {
        Self {
            pos,
            energy,
            occupied: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Isotherms
// ---------------------------------------------------------------------------

/// Langmuir adsorption isotherm.
///
/// Fractional surface coverage:
///
/// ```text
/// θ = K·p / (1 + K·p)
/// ```
///
/// where `K` is the Langmuir constant (Pa⁻¹) and `p` is the gas pressure
/// (Pa).  Returns a value in \[0, 1).
///
/// # Arguments
/// * `pressure`     — gas pressure (Pa)
/// * `k_langmuir`   — Langmuir equilibrium constant (Pa⁻¹)
pub fn langmuir_isotherm(pressure: f64, k_langmuir: f64) -> f64 {
    let kp = k_langmuir * pressure;
    kp / (1.0 + kp)
}

/// Freundlich adsorption isotherm.
///
/// Empirical multi-layer isotherm:
///
/// ```text
/// q = K_f · p^(1/n)
/// ```
///
/// where `K_f` is the Freundlich capacity factor and `n` is the
/// heterogeneity factor.  Common values: `n > 1` (favourable adsorption).
///
/// # Arguments
/// * `pressure` — partial pressure (Pa)
/// * `k_f`      — Freundlich capacity factor
/// * `n`        — Freundlich intensity parameter
pub fn freundlich_isotherm(pressure: f64, k_f: f64, n: f64) -> f64 {
    if pressure <= 0.0 || n == 0.0 {
        return 0.0;
    }
    k_f * pressure.powf(1.0 / n)
}

/// BET (Brunauer–Emmett–Teller) multilayer adsorption isotherm.
///
/// Amount adsorbed:
///
/// ```text
/// V / Vm = C·x / [(1−x)(1−x+C·x)]
/// ```
///
/// where `x = p / p_sat` and `C` is the BET constant.
///
/// # Arguments
/// * `pressure` — gas pressure (Pa)
/// * `p_sat`    — saturation vapour pressure (Pa)
/// * `c_bet`    — BET constant (dimensionless)
/// * `vm`       — monolayer capacity (mol/m² or cm³/g STP)
pub fn bet_isotherm(pressure: f64, p_sat: f64, c_bet: f64, vm: f64) -> f64 {
    if p_sat <= 0.0 || pressure <= 0.0 || pressure >= p_sat {
        return 0.0;
    }
    let x = pressure / p_sat;
    let numerator = c_bet * x;
    let denominator = (1.0 - x) * (1.0 - x + c_bet * x);
    if denominator.abs() < 1e-30 {
        return 0.0;
    }
    vm * numerator / denominator
}

// ---------------------------------------------------------------------------
// SurfaceAdsorption
// ---------------------------------------------------------------------------

/// Monte Carlo surface adsorption model using the Metropolis algorithm.
///
/// Sites are toggled between occupied/unoccupied states based on the
/// chemical potential difference between the gas phase and the adsorption
/// energy.
pub struct SurfaceAdsorption {
    /// Collection of discrete adsorption sites.
    pub sites: Vec<AdsorptionSite>,
    /// Temperature (K).
    pub temp: f64,
    /// Thermal energy k_B · T (J).
    pub kbt: f64,
}

impl SurfaceAdsorption {
    /// Create a new surface adsorption model.
    ///
    /// # Arguments
    /// * `sites` — adsorption sites
    /// * `temp`  — temperature (K)
    /// * `kbt`   — k_B · T (J)
    pub fn new(sites: Vec<AdsorptionSite>, temp: f64, kbt: f64) -> Self {
        Self { sites, temp, kbt }
    }

    /// Perform one Metropolis sweep over all adsorption sites.
    ///
    /// For each site, attempt to flip its occupancy.  The energy change is:
    ///
    /// * Adsorption (empty → occupied): `ΔE = E_site − μ_gas`
    /// * Desorption (occupied → empty): `ΔE = −E_site + μ_gas`
    ///
    /// The move is accepted if `ΔE < 0` or with probability `exp(−ΔE / kBT)`.
    ///
    /// # Arguments
    /// * `pressure`  — gas pressure (for ideal-gas μ reference, currently unused)
    /// * `mu_gas`    — gas chemical potential (J)
    /// * `rng_seed`  — seed for the pseudo-random generator
    pub fn metropolis_step(&mut self, pressure: f64, mu_gas: f64, rng_seed: u64) {
        let _ = pressure; // retained for API completeness
        // Simple LCG pseudo-random to avoid rand dependency duplication
        let mut state = rng_seed.wrapping_add(0x9e3779b97f4a7c15);
        let mut next_f64 = || -> f64 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / (u32::MAX as f64)
        };

        for site in self.sites.iter_mut() {
            let delta_e = if site.occupied {
                // desorption: release binding energy, lose μ_gas
                -site.energy + mu_gas
            } else {
                // adsorption: gain binding energy (energy is negative for attraction)
                site.energy - mu_gas
            };

            let accept = if delta_e <= 0.0 {
                true
            } else {
                let prob = (-delta_e / self.kbt).exp();
                next_f64() < prob
            };

            if accept {
                site.occupied = !site.occupied;
            }
        }
    }

    /// Compute fractional surface coverage θ = (occupied sites) / (total sites).
    pub fn coverage(&self) -> f64 {
        if self.sites.is_empty() {
            return 0.0;
        }
        let n_occ = self.sites.iter().filter(|s| s.occupied).count();
        n_occ as f64 / self.sites.len() as f64
    }

    /// Compute the mean binding energy of currently occupied sites (J).
    ///
    /// Returns 0 if no site is occupied.
    pub fn heat_of_adsorption(&self) -> f64 {
        let occupied: Vec<&AdsorptionSite> = self.sites.iter().filter(|s| s.occupied).collect();
        if occupied.is_empty() {
            return 0.0;
        }
        let sum: f64 = occupied.iter().map(|s| s.energy.abs()).sum();
        sum / occupied.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Surface energy / wetting
// ---------------------------------------------------------------------------

/// Compute the cosine of the contact angle from Young's equation.
///
/// ```text
/// γ_sv = γ_sl + γ_lv · cos θ
/// cos θ = (γ_sv − γ_sl) / γ_lv
/// ```
///
/// Returns a value clamped to \[−1, 1\] corresponding to θ ∈ \[0°, 180°\].
///
/// # Arguments
/// * `gamma_sv` — solid–vapour interfacial energy (J/m²)
/// * `gamma_sl` — solid–liquid interfacial energy (J/m²)
/// * `gamma_lv` — liquid–vapour interfacial energy (J/m²)
pub fn surface_energy(gamma_sv: f64, gamma_sl: f64, gamma_lv: f64) -> f64 {
    if gamma_lv == 0.0 {
        return 0.0;
    }
    ((gamma_sv - gamma_sl) / gamma_lv).clamp(-1.0, 1.0)
}

/// Compute the work of adhesion (Dupré work) W_ad (J/m²).
///
/// ```text
/// W_ad = γ_sv + γ_lv − γ_sl = γ_lv (1 + cos θ)
/// ```
///
/// # Arguments
/// * `gamma_sv`      — solid–vapour interfacial energy (J/m²)
/// * `gamma_lv`      — liquid–vapour interfacial energy (J/m²)
/// * `contact_angle` — contact angle θ (radians)
pub fn work_of_adhesion(gamma_sv: f64, gamma_lv: f64, contact_angle: f64) -> f64 {
    let _ = gamma_sv; // included for physical completeness; Dupré form uses γ_lv + cos θ
    gamma_lv * (1.0 + contact_angle.cos())
}

/// Compute the spreading coefficient S (J/m²).
///
/// ```text
/// S = γ_sv − γ_sl − γ_lv
/// ```
///
/// `S ≥ 0` means spontaneous spreading; `S < 0` means partial wetting.
///
/// # Arguments
/// * `gamma_sv` — solid–vapour interfacial energy (J/m²)
/// * `gamma_sl` — solid–liquid interfacial energy (J/m²)
/// * `gamma_lv` — liquid–vapour interfacial energy (J/m²)
pub fn spreading_coefficient(gamma_sv: f64, gamma_sl: f64, gamma_lv: f64) -> f64 {
    gamma_sv - gamma_sl - gamma_lv
}

// ---------------------------------------------------------------------------
// SurfaceDiffusion
// ---------------------------------------------------------------------------

/// 1D surface diffusion via thermally activated site hopping.
///
/// Adatoms reside on discrete sites.  Hops occur with rate `k = ν₀ exp(−E_b/kBT)`.
pub struct SurfaceDiffusion {
    /// Positions of adatoms (lattice site indices, stored as f64 for MSD).
    pub positions: Vec<f64>,
    /// Activation barriers at each lattice site (J).
    pub barriers: Vec<f64>,
    /// Temperature (K).
    pub temp: f64,
    /// Thermal energy k_B · T (J).
    kbt: f64,
    /// Attempt frequency ν₀ (Hz).
    nu0: f64,
}

impl SurfaceDiffusion {
    /// Create a new surface diffusion model.
    ///
    /// # Arguments
    /// * `positions` — initial adatom positions (lattice units)
    /// * `barriers`  — energy barrier at each lattice site (J)
    /// * `temp`      — temperature (K)
    /// * `kbt`       — k_B · T (J)
    /// * `nu0`       — attempt frequency (Hz)
    pub fn new(positions: Vec<f64>, barriers: Vec<f64>, temp: f64, kbt: f64, nu0: f64) -> Self {
        Self {
            positions,
            barriers,
            temp,
            kbt,
            nu0,
        }
    }

    /// Compute the Arrhenius hop rate for a given barrier.
    ///
    /// ```text
    /// k = ν₀ · exp(−E_b / k_BT)
    /// ```
    pub fn hop_rate(&self, barrier: f64) -> f64 {
        let _ = self.temp; // temperature enters only through kbt
        self.nu0 * (-barrier / self.kbt).exp()
    }

    /// Perform one Monte Carlo step for all adatoms.
    ///
    /// Each adatom attempts to hop left or right by one lattice unit with
    /// probability `k · dt`.  The number of accepted hops is returned.
    ///
    /// # Arguments
    /// * `dt`   — time step (s)
    /// * `seed` — random seed (LCG)
    pub fn monte_carlo_step(&mut self, dt: f64, seed: u64) -> usize {
        let nb = self.barriers.len();
        if nb == 0 {
            return 0;
        }
        let mut state = seed.wrapping_add(0xdeadbeefcafe1234);
        let mut next_f64 = || -> f64 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / (u32::MAX as f64)
        };

        let nu0 = self.nu0;
        let kbt = self.kbt;
        let mut hops = 0usize;
        for pos in self.positions.iter_mut() {
            let site = (pos.round() as usize).clamp(0, nb - 1);
            let barrier = self.barriers[site];
            let rate = nu0 * (-barrier / kbt).exp();
            let prob = (rate * dt).clamp(0.0, 1.0);
            if next_f64() < prob {
                // decide direction
                let dir = if next_f64() < 0.5 { 1.0 } else { -1.0 };
                *pos += dir;
                hops += 1;
            }
        }
        hops
    }

    /// Compute the mean square displacement relative to initial positions.
    ///
    /// `MSD = <(x − x₀)²>`
    pub fn mean_square_displacement(&self, initial: &[f64]) -> f64 {
        if self.positions.is_empty() || initial.is_empty() {
            return 0.0;
        }
        let n = self.positions.len().min(initial.len());
        let s: f64 = self.positions[..n]
            .iter()
            .zip(initial[..n].iter())
            .map(|(x, x0)| (x - x0).powi(2))
            .sum();
        s / n as f64
    }
}

// ---------------------------------------------------------------------------
// Debye-Waller factor
// ---------------------------------------------------------------------------

/// Compute the Debye–Waller factor.
///
/// The intensity of a diffraction peak is reduced by thermal vibrations:
///
/// ```text
/// I / I₀ = exp(−q² `u²` / 3)
/// ```
///
/// where `q` is the scattering vector magnitude (m⁻¹) and `<u²>` is the
/// mean-square displacement amplitude (m²).
///
/// # Arguments
/// * `q`     — scattering vector magnitude (m⁻¹)
/// * `u_rms` — root-mean-square atomic displacement (m)
pub fn debye_waller_factor(q: f64, u_rms: f64) -> f64 {
    (-q * q * u_rms * u_rms / 3.0).exp()
}

// ---------------------------------------------------------------------------
// Surface phonon frequency
// ---------------------------------------------------------------------------

/// Compute the harmonic surface phonon frequency.
///
/// ```text
/// ω = √(k / m)
/// ```
///
/// where `k` is the effective surface force constant (N/m) and `m` the
/// atomic mass (kg).  Returns angular frequency (rad/s).
///
/// # Arguments
/// * `force_const` — restoring force constant (N/m)
/// * `mass`        — atomic mass (kg)
pub fn surface_phonon_freq(force_const: f64, mass: f64) -> f64 {
    if mass <= 0.0 || force_const < 0.0 {
        return 0.0;
    }
    (force_const / mass).sqrt()
}

// ---------------------------------------------------------------------------
// ThinFilmGrowth
// ---------------------------------------------------------------------------

/// 2D ballistic deposition model for thin film growth.
///
/// Each deposited particle sticks to the first surface atom it encounters
/// (nearest-neighbour sticking) or sits on top of the column.
pub struct ThinFilmGrowth {
    /// Height field `height[iy * nx + ix]` — integer column heights.
    pub height: Vec<f64>,
    /// Number of sites in x.
    pub nx: usize,
    /// Number of sites in y.
    pub ny: usize,
}

impl ThinFilmGrowth {
    /// Create a new thin film growth model with all heights at zero.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            height: vec![0.0_f64; nx * ny],
            nx,
            ny,
        }
    }

    /// Deposit one particle at column `(x, y)`.
    ///
    /// The particle sticks if any nearest-neighbour column is taller (shadow
    /// effect); otherwise it lands on top of the current column.
    ///
    /// # Arguments
    /// * `x` — column index in x (will wrap modulo nx)
    /// * `y` — column index in y (will wrap modulo ny)
    pub fn deposit_particle(&mut self, x: usize, y: usize) {
        let nx = self.nx;
        let ny = self.ny;
        let ix = x % nx;
        let iy = y % ny;
        let idx = iy * nx + ix;

        // Check four nearest neighbours for shadow sticking
        let neighbours = [
            ((ix + 1) % nx, iy),
            ((ix + nx - 1) % nx, iy),
            (ix, (iy + 1) % ny),
            (ix, (iy + ny - 1) % ny),
        ];
        let max_neighbour = neighbours
            .iter()
            .map(|(nx_, ny_)| self.height[ny_ * nx + nx_])
            .fold(f64::NEG_INFINITY, f64::max);

        // Stick at the higher of current+1 or max neighbour
        self.height[idx] = (self.height[idx] + 1.0).max(max_neighbour);
    }

    /// Compute the RMS surface roughness σ.
    ///
    /// ```text
    /// σ = √( `h²` − `h`² )
    /// ```
    pub fn roughness(&self) -> f64 {
        let n = self.height.len() as f64;
        if n == 0.0 {
            return 0.0;
        }
        let mean = self.height.iter().sum::<f64>() / n;
        let var = self.height.iter().map(|h| (h - mean).powi(2)).sum::<f64>() / n;
        var.sqrt()
    }

    /// Count the number of connected islands (connected components of elevated
    /// sites, using 4-connectivity flood fill).
    ///
    /// A site is considered "elevated" if its height > 0.
    pub fn island_count(&self) -> usize {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut visited = vec![false; n];
        let mut count = 0usize;

        for start in 0..n {
            if self.height[start] > 0.0 && !visited[start] {
                // BFS
                let mut queue = vec![start];
                visited[start] = true;
                count += 1;
                while let Some(idx) = queue.pop() {
                    let ix = idx % nx;
                    let iy = idx / nx;
                    let neighbours = [
                        (ix + 1, iy),
                        (ix.wrapping_sub(1), iy),
                        (ix, iy + 1),
                        (ix, iy.wrapping_sub(1)),
                    ];
                    for (nx_, ny_) in neighbours {
                        if nx_ < nx && ny_ < ny {
                            let nidx = ny_ * nx + nx_;
                            if self.height[nidx] > 0.0 && !visited[nidx] {
                                visited[nidx] = true;
                                queue.push(nidx);
                            }
                        }
                    }
                }
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// Contact mechanics
// ---------------------------------------------------------------------------

/// Compute the Hertz contact radius for a sphere on a flat surface.
///
/// ```text
/// a = (3FR / 4E*)^(1/3)
/// ```
///
/// where `F` is the normal force (N), `E*` is the reduced elastic modulus
/// (Pa), and `R` is the sphere radius (m).
///
/// # Arguments
/// * `force`   — applied normal force (N)
/// * `e_star`  — reduced modulus E* (Pa)
/// * `r`       — sphere radius (m)
pub fn hertz_contact_radius(force: f64, e_star: f64, r: f64) -> f64 {
    if e_star <= 0.0 || r <= 0.0 || force <= 0.0 {
        return 0.0;
    }
    (3.0 * force * r / (4.0 * e_star)).powf(1.0 / 3.0)
}

/// Compute the DMT (Derjaguin–Muller–Toporov) pull-off adhesion force.
///
/// ```text
/// F_pull-off = 4π γ R
/// ```
///
/// where `γ` is the surface energy per unit area (J/m²) and `R` is the
/// sphere radius (m).
///
/// # Arguments
/// * `e_star` — reduced modulus (Pa; unused in basic DMT but kept for API symmetry)
/// * `r`      — sphere radius (m)
/// * `gamma`  — surface energy (J/m²)
pub fn derjaguin_muller_toporov_adhesion(e_star: f64, r: f64, gamma: f64) -> f64 {
    let _ = e_star; // not needed in the basic DMT pull-off expression
    4.0 * PI * gamma * r
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- AdsorptionSite ---

    #[test]
    fn test_site_new_unoccupied() {
        let s = AdsorptionSite::new([0.0, 0.0, 0.0], -1.0e-20);
        assert!(!s.occupied);
    }

    #[test]
    fn test_site_energy_stored() {
        let s = AdsorptionSite::new([1.0, 2.0, 3.0], -5.0e-21);
        assert!((s.energy - (-5.0e-21)).abs() < 1e-35);
    }

    // --- Langmuir isotherm ---

    #[test]
    fn test_langmuir_zero_pressure() {
        assert_eq!(langmuir_isotherm(0.0, 1.0), 0.0);
    }

    #[test]
    fn test_langmuir_high_pressure() {
        let theta = langmuir_isotherm(1e10, 1.0);
        assert!(theta > 0.99, "theta={theta}");
    }

    #[test]
    fn test_langmuir_half_coverage() {
        // θ = 0.5 when K*p = 1, i.e. p = 1/K
        let k = 2.0;
        let p = 1.0 / k;
        let theta = langmuir_isotherm(p, k);
        assert!((theta - 0.5).abs() < 1e-10, "theta={theta}");
    }

    #[test]
    fn test_langmuir_monotone() {
        let k = 1.0;
        let t1 = langmuir_isotherm(1.0, k);
        let t2 = langmuir_isotherm(2.0, k);
        assert!(t2 > t1);
    }

    // --- Freundlich isotherm ---

    #[test]
    fn test_freundlich_zero_pressure() {
        assert_eq!(freundlich_isotherm(0.0, 1.0, 2.0), 0.0);
    }

    #[test]
    fn test_freundlich_n_equals_1() {
        // When n=1 it reduces to q = K_f * p (linear)
        let q = freundlich_isotherm(3.0, 2.0, 1.0);
        assert!((q - 6.0).abs() < 1e-10, "q={q}");
    }

    #[test]
    fn test_freundlich_n_greater_1() {
        // n > 1: sublinear increase
        let q1 = freundlich_isotherm(1.0, 1.0, 2.0);
        let q8 = freundlich_isotherm(8.0, 1.0, 2.0);
        assert!(q8 < 8.0 * q1, "should be sublinear in pressure");
    }

    #[test]
    fn test_freundlich_positive_result() {
        let q = freundlich_isotherm(10.0, 1.5, 3.0);
        assert!(q > 0.0);
    }

    // --- BET isotherm ---

    #[test]
    fn test_bet_zero_pressure() {
        assert_eq!(bet_isotherm(0.0, 1.0e5, 100.0, 1.0), 0.0);
    }

    #[test]
    fn test_bet_at_saturation() {
        // p >= p_sat returns 0
        assert_eq!(bet_isotherm(1.0, 1.0, 100.0, 1.0), 0.0);
    }

    #[test]
    fn test_bet_positive_below_saturation() {
        let v = bet_isotherm(0.5e5, 1.0e5, 100.0, 1.0);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_bet_increases_with_pressure() {
        let p_sat = 1.0e5;
        let v1 = bet_isotherm(0.1e5, p_sat, 100.0, 1.0);
        let v2 = bet_isotherm(0.5e5, p_sat, 100.0, 1.0);
        assert!(v2 > v1, "v1={v1} v2={v2}");
    }

    // --- SurfaceAdsorption ---

    #[test]
    fn test_coverage_all_empty() {
        let sites: Vec<AdsorptionSite> = (0..10)
            .map(|i| AdsorptionSite::new([i as f64, 0.0, 0.0], -1.0e-20))
            .collect();
        let sa = SurfaceAdsorption::new(sites, 300.0, 4.14e-21);
        assert_eq!(sa.coverage(), 0.0);
    }

    #[test]
    fn test_coverage_all_occupied() {
        let sites: Vec<AdsorptionSite> = (0..10)
            .map(|i| {
                let mut s = AdsorptionSite::new([i as f64, 0.0, 0.0], -1.0e-20);
                s.occupied = true;
                s
            })
            .collect();
        let sa = SurfaceAdsorption::new(sites, 300.0, 4.14e-21);
        assert_eq!(sa.coverage(), 1.0);
    }

    #[test]
    fn test_heat_of_adsorption_empty() {
        let sites = vec![AdsorptionSite::new([0.0, 0.0, 0.0], -1.0e-20)];
        let sa = SurfaceAdsorption::new(sites, 300.0, 4.14e-21);
        assert_eq!(sa.heat_of_adsorption(), 0.0);
    }

    #[test]
    fn test_heat_of_adsorption_occupied() {
        let mut s = AdsorptionSite::new([0.0, 0.0, 0.0], -2.0e-20);
        s.occupied = true;
        let sa = SurfaceAdsorption::new(vec![s], 300.0, 4.14e-21);
        assert!((sa.heat_of_adsorption() - 2.0e-20).abs() < 1e-33);
    }

    #[test]
    fn test_metropolis_step_runs() {
        let sites: Vec<AdsorptionSite> = (0..50)
            .map(|i| AdsorptionSite::new([i as f64, 0.0, 0.0], -1.0e-20))
            .collect();
        let mut sa = SurfaceAdsorption::new(sites, 300.0, 4.14e-21);
        sa.metropolis_step(1.0e5, -5.0e-21, 42);
        // coverage should be in [0,1]
        let cov = sa.coverage();
        assert!((0.0..=1.0).contains(&cov));
    }

    #[test]
    fn test_coverage_range() {
        let sites: Vec<AdsorptionSite> = (0..20)
            .map(|i| AdsorptionSite::new([i as f64, 0.0, 0.0], -3.0e-20))
            .collect();
        let mut sa = SurfaceAdsorption::new(sites, 300.0, 4.14e-21);
        for seed in 0..10u64 {
            sa.metropolis_step(1e5, -1e-20, seed * 12345 + 1);
        }
        let cov = sa.coverage();
        assert!((0.0..=1.0).contains(&cov), "cov={cov}");
    }

    // --- Surface energy / wetting ---

    #[test]
    fn test_young_complete_wetting() {
        // cos θ = 1 when γ_sv - γ_sl = γ_lv
        let cos_theta = surface_energy(0.1, 0.05, 0.05);
        assert!((cos_theta - 1.0).abs() < 1e-10, "cos_theta={cos_theta}");
    }

    #[test]
    fn test_young_clamped() {
        // Even if γ_sv >> γ_sl, cos θ ≤ 1
        let cos_theta = surface_energy(1.0, 0.0, 0.01);
        assert!(cos_theta <= 1.0);
    }

    #[test]
    fn test_young_zero_gamma_lv() {
        assert_eq!(surface_energy(0.1, 0.05, 0.0), 0.0);
    }

    #[test]
    fn test_work_of_adhesion_theta_zero() {
        // cos 0 = 1 → W = 2γ_lv
        let w = work_of_adhesion(0.05, 0.05, 0.0);
        assert!((w - 2.0 * 0.05).abs() < 1e-12, "w={w}");
    }

    #[test]
    fn test_work_of_adhesion_theta_90() {
        // cos(π/2) = 0 → W = γ_lv
        let w = work_of_adhesion(0.05, 0.05, PI / 2.0);
        assert!((w - 0.05).abs() < 1e-12, "w={w}");
    }

    #[test]
    fn test_spreading_positive() {
        // S > 0: spontaneous spreading
        let s = spreading_coefficient(0.3, 0.1, 0.1);
        assert!(s > 0.0, "s={s}");
    }

    #[test]
    fn test_spreading_negative() {
        // S < 0: partial wetting
        let s = spreading_coefficient(0.1, 0.3, 0.3);
        assert!(s < 0.0, "s={s}");
    }

    // --- SurfaceDiffusion ---

    #[test]
    fn test_hop_rate_positive() {
        let sd = SurfaceDiffusion::new(vec![0.0], vec![1.0e-20], 300.0, 4.14e-21, 1e13);
        let rate = sd.hop_rate(1.0e-20);
        assert!(rate > 0.0, "rate={rate}");
    }

    #[test]
    fn test_hop_rate_decreases_with_barrier() {
        let sd = SurfaceDiffusion::new(vec![0.0], vec![0.1], 300.0, 4.14e-21, 1e13);
        let r1 = sd.hop_rate(0.01e-19);
        let r2 = sd.hop_rate(1.0e-19);
        assert!(r1 > r2, "r1={r1} r2={r2}");
    }

    #[test]
    fn test_msd_zero_at_start() {
        let pos = vec![0.0, 1.0, 2.0];
        let initial = pos.clone();
        let sd = SurfaceDiffusion::new(pos, vec![0.1; 3], 300.0, 4.14e-21, 1e13);
        assert_eq!(sd.mean_square_displacement(&initial), 0.0);
    }

    #[test]
    fn test_msd_after_displacement() {
        let pos = vec![1.0, 2.0, 3.0];
        let initial = vec![0.0, 0.0, 0.0];
        let sd = SurfaceDiffusion::new(pos, vec![0.01; 3], 300.0, 4.14e-21, 1e13);
        let msd = sd.mean_square_displacement(&initial);
        // (1+4+9)/3 = 14/3
        assert!((msd - 14.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_mc_step_returns_count() {
        let barriers = vec![0.01e-20; 10]; // very small barrier → many hops
        let mut sd = SurfaceDiffusion::new(vec![5.0; 10], barriers, 300.0, 4.14e-21, 1e13);
        let h = sd.monte_carlo_step(1e-12, 99);
        // count should be between 0 and 10
        assert!(h <= 10);
    }

    #[test]
    fn test_mc_step_empty() {
        let mut sd = SurfaceDiffusion::new(vec![], vec![], 300.0, 4.14e-21, 1e13);
        assert_eq!(sd.monte_carlo_step(1e-12, 1), 0);
    }

    // --- Debye-Waller ---

    #[test]
    fn test_debye_waller_zero_q() {
        assert_eq!(debye_waller_factor(0.0, 1.0e-11), 1.0);
    }

    #[test]
    fn test_debye_waller_zero_u() {
        assert_eq!(debye_waller_factor(1.0e10, 0.0), 1.0);
    }

    #[test]
    fn test_debye_waller_decreases_with_q() {
        let u = 1.0e-11;
        let dw1 = debye_waller_factor(1.0e9, u);
        let dw2 = debye_waller_factor(1.0e10, u);
        assert!(dw2 < dw1, "dw1={dw1} dw2={dw2}");
    }

    #[test]
    fn test_debye_waller_range() {
        let dw = debye_waller_factor(1.0e10, 5.0e-11);
        assert!((0.0..=1.0).contains(&dw), "dw={dw}");
    }

    // --- Surface phonon ---

    #[test]
    fn test_phonon_freq_positive() {
        let f = surface_phonon_freq(10.0, 1.0e-26);
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn test_phonon_freq_zero_mass() {
        assert_eq!(surface_phonon_freq(10.0, 0.0), 0.0);
    }

    #[test]
    fn test_phonon_freq_negative_force_const() {
        assert_eq!(surface_phonon_freq(-1.0, 1.0e-26), 0.0);
    }

    #[test]
    fn test_phonon_freq_sqrt_relation() {
        let k: f64 = 5.0;
        let m: f64 = 2.0e-26;
        let expected = (k / m).sqrt();
        let got = surface_phonon_freq(k, m);
        assert!((got - expected).abs() < 1e-10 * expected);
    }

    // --- ThinFilmGrowth ---

    #[test]
    fn test_film_new_zero_height() {
        let f = ThinFilmGrowth::new(8, 8);
        assert!(f.height.iter().all(|h| *h == 0.0));
    }

    #[test]
    fn test_film_deposit_increases_height() {
        let mut f = ThinFilmGrowth::new(8, 8);
        f.deposit_particle(3, 3);
        let idx = 3 * 8 + 3;
        assert!(f.height[idx] >= 1.0);
    }

    #[test]
    fn test_film_roughness_zero() {
        let f = ThinFilmGrowth::new(4, 4);
        assert_eq!(f.roughness(), 0.0);
    }

    #[test]
    fn test_film_roughness_nonzero() {
        let mut f = ThinFilmGrowth::new(8, 8);
        // deposit on one column only
        for _ in 0..5 {
            f.deposit_particle(0, 0);
        }
        assert!(f.roughness() > 0.0);
    }

    #[test]
    fn test_film_island_count_empty() {
        let f = ThinFilmGrowth::new(8, 8);
        assert_eq!(f.island_count(), 0);
    }

    #[test]
    fn test_film_island_count_one_island() {
        let mut f = ThinFilmGrowth::new(8, 8);
        f.deposit_particle(3, 3);
        // The deposit_particle with shadow sticking may grow multiple sites,
        // but they form one connected island
        let count = f.island_count();
        assert!(count >= 1, "count={count}");
    }

    #[test]
    fn test_film_island_count_two_separate() {
        let mut f = ThinFilmGrowth::new(16, 16);
        // Two isolated deposits far apart
        let idx0 = 0;
        let idx1 = 10 * 16 + 10;
        f.height[idx0] = 1.0;
        f.height[idx1] = 1.0;
        assert_eq!(f.island_count(), 2);
    }

    #[test]
    fn test_film_deposit_wraps() {
        let mut f = ThinFilmGrowth::new(4, 4);
        f.deposit_particle(5, 6); // 5%4=1, 6%4=2
        let idx = 2 * 4 + 1;
        assert!(f.height[idx] >= 1.0);
    }

    // --- Hertz contact ---

    #[test]
    fn test_hertz_positive() {
        let a = hertz_contact_radius(1.0, 1.0e9, 1.0e-3);
        assert!(a > 0.0, "a={a}");
    }

    #[test]
    fn test_hertz_zero_force() {
        assert_eq!(hertz_contact_radius(0.0, 1.0e9, 1.0e-3), 0.0);
    }

    #[test]
    fn test_hertz_increases_with_force() {
        let a1 = hertz_contact_radius(1.0, 1.0e9, 1.0e-3);
        let a2 = hertz_contact_radius(8.0, 1.0e9, 1.0e-3);
        // a ∝ F^(1/3) → a(8F) = 2*a(F)
        assert!((a2 / a1 - 2.0).abs() < 1e-8, "ratio={}", a2 / a1);
    }

    // --- DMT adhesion ---

    #[test]
    fn test_dmt_positive() {
        let f = derjaguin_muller_toporov_adhesion(1.0e9, 1.0e-6, 0.05);
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn test_dmt_linear_in_r() {
        let f1 = derjaguin_muller_toporov_adhesion(1.0e9, 1.0e-6, 0.05);
        let f2 = derjaguin_muller_toporov_adhesion(1.0e9, 2.0e-6, 0.05);
        assert!((f2 / f1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_dmt_formula() {
        let e_star = 1.0e9;
        let r = 1.0e-6;
        let gamma = 0.05;
        let expected = 4.0 * PI * gamma * r;
        let got = derjaguin_muller_toporov_adhesion(e_star, r, gamma);
        assert!((got - expected).abs() < 1e-20 * expected.abs().max(1.0));
    }
}
