// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Glass-forming molecular dynamics module.
//!
//! Covers:
//! - Vogel-Fulcher-Tammann (VFT) and Williams-Landel-Ferry (WLF) viscosity equations
//! - Kauzmann temperature estimation
//! - Kohlrausch-Williams-Watts (KWW) stretched-exponential relaxation
//! - Glass transition detection from volume–temperature data
//! - Fragility classification
//! - BKS (van Beest-Kramer-van Santen) pair potential for silica
//! - Radial distribution function, coordination number, and structure factor
//! - Mean-squared displacement and non-Gaussian parameter
//! - `GlassSimulation` driver with anneal and quench protocols

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Permittivity of free space (F m⁻¹).
const EPS_0: f64 = 8.854_187_817e-12;

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors (a − b).
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ---------------------------------------------------------------------------
// GlassFormingParams
// ---------------------------------------------------------------------------

/// Thermophysical parameters describing a glass-forming material.
#[derive(Debug, Clone)]
pub struct GlassFormingParams {
    /// Fragility index m (dimensionless; strong ≈ 16, fragile ≈ 200).
    pub fragility: f64,
    /// Glass-transition temperature (K).
    pub glass_transition_temp: f64,
    /// Melting (liquidus) temperature (K).
    pub melting_temp: f64,
    /// Thermal cooling rate (K s⁻¹).
    pub cooling_rate: f64,
    /// Mass density at Tg (kg m⁻³).
    pub density: f64,
}

impl Default for GlassFormingParams {
    fn default() -> Self {
        Self {
            fragility: 40.0,
            glass_transition_temp: 1473.0, // SiO₂-like
            melting_temp: 1973.0,
            cooling_rate: 10.0,
            density: 2200.0,
        }
    }
}

// ---------------------------------------------------------------------------
// GlassType
// ---------------------------------------------------------------------------

/// Classification of glass-forming material type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GlassType {
    /// Amorphous SiO₂ (vitreous silica).
    Silica,
    /// Metallic glass (e.g. Zr-Cu alloys).
    Metallic,
    /// Amorphous polymer (e.g. PMMA, polystyrene).
    Polymer,
    /// Chalcogenide glass (e.g. As₂S₃, GeSe₂).
    Chalcogenide,
    /// Oxide glass other than silica (e.g. GeO₂, B₂O₃).
    Oxide,
    /// Borosilicate glass (e.g. Pyrex).
    Borosilicate,
}

// ---------------------------------------------------------------------------
// Vogel-Fulcher-Tammann viscosity equation
// ---------------------------------------------------------------------------

/// Vogel-Fulcher-Tammann (VFT) viscosity model.
pub struct VogelFulcherTammann;

impl VogelFulcherTammann {
    /// Compute the VFT viscosity at temperature `t`.
    ///
    /// η(T) = exp(A + B / (T − T₀))
    ///
    /// # Arguments
    /// * `t`  – temperature (K).
    /// * `a`  – prefactor (dimensionless log-viscosity offset).
    /// * `b`  – activation-like parameter (K).
    /// * `t0` – Vogel temperature (K).
    pub fn viscosity(t: f64, a: f64, b: f64, t0: f64) -> f64 {
        let denom = t - t0;
        if denom.abs() < 1.0e-9 {
            return f64::INFINITY;
        }
        (a + b / denom).exp()
    }

    /// Fit VFT parameters (A, B, T₀) to measured temperature–viscosity data
    /// using a simple iterative least-squares approach (Gauss–Newton, 50 iterations).
    ///
    /// # Arguments
    /// * `temps`       – measured temperatures (K), length ≥ 3.
    /// * `viscosities` – measured viscosities (Pa s), same length.
    ///
    /// # Returns
    /// `(A, B, T₀)` parameters.
    pub fn fit_vft(temps: &[f64], viscosities: &[f64]) -> (f64, f64, f64) {
        assert_eq!(temps.len(), viscosities.len());
        assert!(temps.len() >= 3);

        // Working in log-space: log(η) = A + B/(T - T0)
        let log_vis: Vec<f64> = viscosities.iter().map(|v| v.ln()).collect();

        // Initial guesses
        let mut a_val = 0.0_f64;
        let mut b_val = 4000.0_f64;
        let mut t0 = temps.iter().cloned().fold(f64::INFINITY, f64::min) * 0.8;

        let n = temps.len() as f64;
        let lr = 1.0e-4;

        for _ in 0..100 {
            let mut grad_a = 0.0_f64;
            let mut grad_b = 0.0_f64;
            let mut grad_t0 = 0.0_f64;

            for (i, &ti) in temps.iter().enumerate() {
                let dt = ti - t0;
                if dt.abs() < 1.0e-9 {
                    continue;
                }
                let pred = a_val + b_val / dt;
                let err = pred - log_vis[i];
                grad_a += 2.0 * err / n;
                grad_b += 2.0 * err / (n * dt);
                grad_t0 += 2.0 * err * b_val / (n * dt * dt);
            }

            a_val -= lr * grad_a;
            b_val -= lr * grad_b;
            t0 -= lr * grad_t0;
        }

        (a_val, b_val, t0)
    }
}

// ---------------------------------------------------------------------------
// Williams-Landel-Ferry (WLF) equation
// ---------------------------------------------------------------------------

/// Williams-Landel-Ferry (WLF) shift factor model for time-temperature superposition.
pub struct WilliamsLandelFerry;

impl WilliamsLandelFerry {
    /// Compute the WLF log₁₀ shift factor log₁₀(aT).
    ///
    /// log₁₀(aT) = −C₁·(T − Tref) / (C₂ + T − Tref)
    ///
    /// # Arguments
    /// * `t`     – temperature of interest (K).
    /// * `t_ref` – reference temperature (K).
    /// * `c1`    – WLF constant C₁.
    /// * `c2`    – WLF constant C₂ (K).
    ///
    /// # Returns
    /// log₁₀(aT).
    pub fn shift_factor(t: f64, t_ref: f64, c1: f64, c2: f64) -> f64 {
        let delta = t - t_ref;
        let denom = c2 + delta;
        if denom.abs() < 1.0e-9 {
            return f64::NEG_INFINITY;
        }
        -c1 * delta / denom
    }
}

// ---------------------------------------------------------------------------
// Kauzmann temperature
// ---------------------------------------------------------------------------

/// Estimate the Kauzmann temperature Tk from the glass transition temperature
/// and the fragility index m using the empirical relation
///
/// Tk ≈ Tg · (1 − 1/m).
///
/// # Arguments
/// * `tg` – glass-transition temperature (K).
/// * `m`  – kinetic fragility index (dimensionless, typically 16–200).
pub fn kauzmann_temp(tg: f64, m: f64) -> f64 {
    if m.abs() < 1.0e-9 {
        return tg;
    }
    tg * (1.0 - 1.0 / m)
}

// ---------------------------------------------------------------------------
// Structural relaxation
// ---------------------------------------------------------------------------

/// Structural relaxation models.
pub struct StructuralRelaxation;

impl StructuralRelaxation {
    /// Kohlrausch-Williams-Watts (KWW) stretched-exponential decay function.
    ///
    /// Φ(t) = exp(−(t/τ)^β)
    ///
    /// # Arguments
    /// * `t`    – time (s).
    /// * `tau`  – characteristic relaxation time (s).
    /// * `beta` – stretching exponent (0 < β ≤ 1).
    pub fn kww_relaxation(t: f64, tau: f64, beta: f64) -> f64 {
        if tau <= 0.0 || t < 0.0 {
            return 1.0;
        }
        (-(t / tau).powf(beta)).exp()
    }

    /// Average relaxation time for a KWW process, computed as
    ///
    /// <τ> = (τ/β) · Γ(1/β)
    ///
    /// where Γ is the Gamma function approximated via Lanczos.
    ///
    /// # Arguments
    /// * `tau`  – KWW time parameter (s).
    /// * `beta` – KWW stretching exponent.
    pub fn average_relaxation_time(tau: f64, beta: f64) -> f64 {
        if beta <= 0.0 || tau <= 0.0 {
            return 0.0;
        }
        (tau / beta) * gamma_approx(1.0 / beta)
    }
}

/// Gamma function approximation via the Lanczos series.
fn gamma_approx(z: f64) -> f64 {
    // Lanczos approximation coefficients (g=7)
    let g = 7.0_f64;
    let c: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_9,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_571e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        PI / ((PI * z).sin() * gamma_approx(1.0 - z))
    } else {
        let z = z - 1.0;
        let mut x = c[0];
        for (i, ci) in c.iter().enumerate().skip(1) {
            x += ci / (z + i as f64);
        }
        let t = z + g + 0.5;
        (2.0 * PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
    }
}

// ---------------------------------------------------------------------------
// Glass transition detection
// ---------------------------------------------------------------------------

/// Detect the glass transition temperature from volume–temperature data.
///
/// Searches for the largest discontinuity in the slope dV/dT, which corresponds
/// to the kink (change in thermal expansion coefficient) at Tg.
///
/// # Arguments
/// * `temps`   – temperature array (K), must be monotonically decreasing (cooling scan).
/// * `volumes` – specific volume array (m³ kg⁻¹), same length as `temps`.
///
/// # Returns
/// `Some(Tg)` if a kink is detected, `None` if the array is too short.
pub fn detect_glass_transition(temps: &[f64], volumes: &[f64]) -> Option<f64> {
    let n = temps.len();
    if n < 4 || volumes.len() != n {
        return None;
    }

    // Compute finite-difference slopes dV/dT
    let mut slopes = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let dt = temps[i + 1] - temps[i];
        let dv = volumes[i + 1] - volumes[i];
        if dt.abs() < 1.0e-30 {
            slopes.push(0.0);
        } else {
            slopes.push(dv / dt);
        }
    }

    // Find maximum change in slope (second derivative peak)
    let mut max_change = 0.0_f64;
    let mut max_idx = 1_usize;
    for i in 1..slopes.len() {
        let change = (slopes[i] - slopes[i - 1]).abs();
        if change > max_change {
            max_change = change;
            max_idx = i;
        }
    }

    Some((temps[max_idx] + temps[max_idx + 1]) / 2.0)
}

// ---------------------------------------------------------------------------
// Fragility classification
// ---------------------------------------------------------------------------

/// Return a qualitative fragility label based on the kinetic fragility index m.
///
/// * m < 30     → "strong"
/// * 30 ≤ m < 80 → "intermediate"
/// * m ≥ 80     → "fragile"
pub fn fragility_index(m: f64) -> &'static str {
    if m < 30.0 {
        "strong"
    } else if m < 80.0 {
        "intermediate"
    } else {
        "fragile"
    }
}

// ---------------------------------------------------------------------------
// BKS potential for SiO₂
// ---------------------------------------------------------------------------

/// BKS (van Beest-Kramer-van Santen) pair potential used for classical silica MD.
pub struct BksPotential;

impl BksPotential {
    /// Compute the BKS pair energy between atoms of species i and j.
    ///
    /// V(r) = q₁·q₂·e² / (4πε₀·r) + A·exp(−B·r) − C/r⁶
    ///
    /// # Arguments
    /// * `r`       – interatomic distance (m).
    /// * `a`       – repulsive pre-factor (J).
    /// * `b_param` – repulsive exponent (m⁻¹).
    /// * `c`       – dispersion coefficient (J m⁶).
    /// * `q1`      – charge on atom 1 (in units of e).
    /// * `q2`      – charge on atom 2 (in units of e).
    pub fn bks_energy(r: f64, a: f64, b_param: f64, c: f64, q1: f64, q2: f64) -> f64 {
        if r < 1.0e-15 {
            return 0.0;
        }
        let coulomb = q1 * q2 * E_CHARGE * E_CHARGE / (4.0 * PI * EPS_0 * r);
        let repulsive = a * (-b_param * r).exp();
        let dispersion = c / r.powi(6);
        coulomb + repulsive - dispersion
    }
}

// ---------------------------------------------------------------------------
// GlassStructure analysis
// ---------------------------------------------------------------------------

/// Structural analysis utilities for glass simulations.
pub struct GlassStructure;

impl GlassStructure {
    /// Compute the pair radial distribution function g(r) from particle positions.
    ///
    /// # Arguments
    /// * `positions` – atomic positions (m).
    /// * `r_max`     – maximum radial distance (m).
    /// * `n_bins`    – number of histogram bins.
    ///
    /// # Returns
    /// Normalised g(r) histogram (length = `n_bins`).
    pub fn radial_distribution_function(
        positions: &[[f64; 3]],
        r_max: f64,
        n_bins: usize,
    ) -> Vec<f64> {
        let n = positions.len();
        if n < 2 || n_bins == 0 {
            return vec![0.0; n_bins];
        }

        let dr = r_max / n_bins as f64;
        let mut hist = vec![0_u64; n_bins];

        for i in 0..n {
            for j in (i + 1)..n {
                let r = norm3(sub3(positions[i], positions[j]));
                if r < r_max {
                    let bin = (r / dr) as usize;
                    if bin < n_bins {
                        hist[bin] += 2; // count pair twice (i→j and j→i)
                    }
                }
            }
        }

        // Normalize by ideal gas shell volume
        let rho = n as f64 / (4.0 / 3.0 * PI * r_max.powi(3));
        let n_f = n as f64;
        let mut rdf = vec![0.0_f64; n_bins];
        for (k, &count) in hist.iter().enumerate() {
            let r_inner = k as f64 * dr;
            let r_outer = r_inner + dr;
            let shell_vol = 4.0 / 3.0 * PI * (r_outer.powi(3) - r_inner.powi(3));
            let ideal = rho * n_f * shell_vol;
            rdf[k] = if ideal > 0.0 {
                count as f64 / ideal
            } else {
                0.0
            };
        }
        rdf
    }

    /// Compute the coordination number by integrating g(r) up to the first minimum.
    ///
    /// N_coord = 4πρ ∫_{r_min}^{r_max} g(r)·r²·dr
    ///
    /// # Arguments
    /// * `rdf`   – g(r) histogram values.
    /// * `r_min` – integration lower bound (m).
    /// * `r_max` – integration upper bound (m).
    /// * `dr`    – bin width (m).
    /// * `rho`   – number density (m⁻³).
    pub fn coordination_number(rdf: &[f64], r_min: f64, r_max: f64, dr: f64, rho: f64) -> f64 {
        let mut integral = 0.0_f64;
        for (k, &gk) in rdf.iter().enumerate() {
            let r = (k as f64 + 0.5) * dr;
            if r >= r_min && r <= r_max {
                integral += gk * r * r * dr;
            }
        }
        4.0 * PI * rho * integral
    }

    /// Compute the static structure factor S(q) at a single wave vector magnitude q.
    ///
    /// S(q) = (1/N) |Σ_j exp(i·q·r_j)|²  (isotropic average along x-axis)
    ///
    /// # Arguments
    /// * `positions` – atomic positions (m).
    /// * `q`         – wave vector magnitude (m⁻¹).
    pub fn structure_factor_q(positions: &[[f64; 3]], q: f64) -> f64 {
        let n = positions.len();
        if n == 0 {
            return 0.0;
        }
        let mut re_sum = 0.0_f64;
        let mut im_sum = 0.0_f64;
        for pos in positions {
            // Average over three Cartesian directions
            for &ri in pos.iter() {
                let phase = q * ri;
                re_sum += phase.cos();
                im_sum += phase.sin();
            }
        }
        let norm = 3.0 * n as f64;
        (re_sum * re_sum + im_sum * im_sum) / norm
    }
}

// ---------------------------------------------------------------------------
// Mean-squared displacement
// ---------------------------------------------------------------------------

/// Compute the mean-squared displacement (MSD) between two configurations.
///
/// MSD = (1/N) Σ_i |r_i(t) − r_i(0)|²
///
/// # Arguments
/// * `positions_t0` – initial positions (m).
/// * `positions_t`  – positions at time t (m).
pub fn msd_glass(positions_t0: &[[f64; 3]], positions_t: &[[f64; 3]]) -> f64 {
    let n = positions_t0.len().min(positions_t.len());
    if n == 0 {
        return 0.0;
    }
    let sum: f64 = positions_t0[..n]
        .iter()
        .zip(positions_t[..n].iter())
        .map(|(r0, r)| {
            let dr = sub3(*r, *r0);
            dot3(dr, dr)
        })
        .sum();
    sum / n as f64
}

// ---------------------------------------------------------------------------
// Non-Gaussian parameter
// ---------------------------------------------------------------------------

/// Compute the non-Gaussian parameter α₂(t) characterising heterogeneous dynamics.
///
/// α₂ = (3/5) · <r⁴> / `r²`² − 1
///
/// # Arguments
/// * `positions_t0` – initial positions (m).
/// * `positions_t`  – positions at time t (m).
pub fn non_gaussian_param(positions_t0: &[[f64; 3]], positions_t: &[[f64; 3]]) -> f64 {
    let n = positions_t0.len().min(positions_t.len());
    if n == 0 {
        return 0.0;
    }
    let mut r2_sum = 0.0_f64;
    let mut r4_sum = 0.0_f64;
    for (r0, r) in positions_t0[..n].iter().zip(positions_t[..n].iter()) {
        let dr = sub3(*r, *r0);
        let r2 = dot3(dr, dr);
        r2_sum += r2;
        r4_sum += r2 * r2;
    }
    let mean_r2 = r2_sum / n as f64;
    let mean_r4 = r4_sum / n as f64;
    if mean_r2.abs() < 1.0e-30 {
        return 0.0;
    }
    3.0 / 5.0 * mean_r4 / (mean_r2 * mean_r2) - 1.0
}

// ---------------------------------------------------------------------------
// Vitrification criterion
// ---------------------------------------------------------------------------

/// Return `true` if the material vitrifies (forms a glass) at the given
/// cooling rate.
///
/// A material vitrifies when the cooling rate exceeds the critical cooling
/// rate required for crystallisation to be bypassed.
///
/// # Arguments
/// * `cooling_rate`   – actual cooling rate (K s⁻¹).
/// * `critical_rate`  – minimum cooling rate for glass formation (K s⁻¹).
pub fn is_glass(cooling_rate: f64, critical_rate: f64) -> bool {
    cooling_rate >= critical_rate
}

// ---------------------------------------------------------------------------
// GlassSimulation
// ---------------------------------------------------------------------------

/// A simple MD driver for glass-forming simulations.
///
/// Uses velocity-Verlet integration with a Lennard-Jones pair potential for
/// demonstration; real-world usage should swap in BKS or another potential.
#[derive(Debug, Clone)]
pub struct GlassSimulation {
    /// Atomic positions (m).
    pub positions: Vec<[f64; 3]>,
    /// Atomic velocities (m s⁻¹).
    pub velocities: Vec<[f64; 3]>,
    /// Atomic masses (kg).
    pub masses: Vec<f64>,
    /// Current temperature (K).
    pub temperature: f64,
    /// Glass-forming parameters.
    pub params: GlassFormingParams,
    /// Glass type.
    pub glass_type: GlassType,
    /// Time step (s).
    pub dt: f64,
}

impl GlassSimulation {
    /// Construct a new `GlassSimulation`.
    pub fn new(
        positions: Vec<[f64; 3]>,
        velocities: Vec<[f64; 3]>,
        masses: Vec<f64>,
        temperature: f64,
        params: GlassFormingParams,
        glass_type: GlassType,
        dt: f64,
    ) -> Self {
        Self {
            positions,
            velocities,
            masses,
            temperature,
            params,
            glass_type,
            dt,
        }
    }

    /// Compute simple Lennard-Jones forces for all pairs.
    fn compute_forces_lj(&self, eps: f64, sigma: f64) -> Vec<[f64; 3]> {
        let n = self.positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let dr = sub3(self.positions[j], self.positions[i]);
                let r2 = dot3(dr, dr).max(1.0e-30);
                let r = r2.sqrt();
                let sr = sigma / r;
                let sr6 = sr.powi(6);
                let sr12 = sr6 * sr6;
                // F = -dU/dr hat{r}; magnitude = 4ε(12σ¹²/r¹³ - 6σ⁶/r⁷)
                let f_mag = 24.0 * eps * (2.0 * sr12 - sr6) / (r2);
                let fvec = scale3(dr, f_mag);
                forces[i] = add3(forces[i], scale3(fvec, -1.0));
                forces[j] = add3(forces[j], fvec);
            }
        }
        forces
    }

    /// Perform a single velocity-Verlet MD step.
    pub fn step(&mut self) {
        let forces = self.compute_forces_lj(1.0e-21, 3.4e-10);
        let n = self.positions.len();
        let dt = self.dt;
        // Velocity Verlet
        for (i, (vel, pos)) in self
            .velocities
            .iter_mut()
            .zip(self.positions.iter_mut())
            .enumerate()
            .take(n)
        {
            let m = self.masses[i];
            let acc = scale3(forces[i], 1.0 / m);
            // half-kick velocity
            *vel = add3(*vel, scale3(acc, 0.5 * dt));
            // update position
            *pos = add3(*pos, scale3(*vel, dt));
        }
        let forces2 = self.compute_forces_lj(1.0e-21, 3.4e-10);
        for (i, vel) in self.velocities.iter_mut().enumerate().take(n) {
            let m = self.masses[i];
            let acc = scale3(forces2[i], 1.0 / m);
            *vel = add3(*vel, scale3(acc, 0.5 * dt));
        }

        // Rescale to current temperature (simple velocity rescaling thermostat)
        self.rescale_velocities();
    }

    /// Rescale velocities to match `self.temperature`.
    fn rescale_velocities(&mut self) {
        let n = self.positions.len();
        if n == 0 {
            return;
        }
        let ke: f64 = self
            .velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * dot3(*v, *v))
            .sum();
        let t_current = 2.0 * ke / (3.0 * n as f64 * K_B);
        if t_current < 1.0e-30 {
            return;
        }
        let scale = (self.temperature / t_current).sqrt();
        for v in self.velocities.iter_mut() {
            *v = scale3(*v, scale);
        }
    }

    /// Anneal the simulation: linearly ramp temperature from current value to
    /// `t_final` over `n_steps` MD steps.
    pub fn anneal(&mut self, t_final: f64, n_steps: usize) {
        let t_start = self.temperature;
        for step_idx in 0..n_steps {
            let frac = step_idx as f64 / n_steps.max(1) as f64;
            self.temperature = t_start + (t_final - t_start) * frac;
            self.step();
        }
        self.temperature = t_final;
    }

    /// Quench the simulation to the target temperature at the specified cooling
    /// rate, applying one MD step per degree.
    pub fn quench(&mut self, t_final: f64) {
        let cooling_rate = self.params.cooling_rate;
        let dt_therm = 1.0 / cooling_rate.max(1.0e-30);
        let n_steps = ((self.temperature - t_final) / (cooling_rate * self.dt + 1.0e-30)) as usize;
        let n_steps = n_steps.min(10_000);
        let t_start = self.temperature;
        for step_idx in 0..n_steps {
            let frac = step_idx as f64 / n_steps.max(1) as f64;
            self.temperature = (t_start - frac * (t_start - t_final)).max(t_final);
            let _ = dt_therm;
            self.step();
        }
        self.temperature = t_final;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- vector helpers ----

    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_norm3() {
        assert!((norm3([3.0, 4.0, 0.0]) - 5.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_add3_sub3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_eq!(add3(a, b), [5.0, 7.0, 9.0]);
        assert_eq!(sub3(b, a), [3.0, 3.0, 3.0]);
    }

    // ---- GlassFormingParams ----

    #[test]
    fn test_params_default() {
        let p = GlassFormingParams::default();
        assert!(p.glass_transition_temp > 0.0);
        assert!(p.melting_temp > p.glass_transition_temp);
    }

    // ---- VFT viscosity ----

    #[test]
    fn test_vft_viscosity_finite() {
        let eta = VogelFulcherTammann::viscosity(1800.0, -3.5, 4000.0, 500.0);
        assert!(eta.is_finite() && eta > 0.0);
    }

    #[test]
    fn test_vft_viscosity_at_t0_is_inf() {
        let eta = VogelFulcherTammann::viscosity(500.0, -3.5, 4000.0, 500.0);
        assert!(eta.is_infinite());
    }

    #[test]
    fn test_vft_viscosity_decreases_with_temperature() {
        let eta_low = VogelFulcherTammann::viscosity(1500.0, -3.5, 4000.0, 500.0);
        let eta_high = VogelFulcherTammann::viscosity(2000.0, -3.5, 4000.0, 500.0);
        assert!(eta_high < eta_low);
    }

    #[test]
    fn test_vft_fit_returns_triple() {
        let temps = [1500.0, 1700.0, 2000.0, 2500.0];
        let vis = [1.0e6, 1.0e4, 100.0, 1.0];
        let (a, b, t0) = VogelFulcherTammann::fit_vft(&temps, &vis);
        assert!(a.is_finite());
        assert!(b.is_finite());
        assert!(t0.is_finite());
    }

    // ---- WLF ----

    #[test]
    fn test_wlf_at_reference_is_zero() {
        let log_at = WilliamsLandelFerry::shift_factor(400.0, 400.0, 8.86, 101.6);
        assert!(log_at.abs() < 1.0e-12);
    }

    #[test]
    fn test_wlf_above_reference_negative() {
        let log_at = WilliamsLandelFerry::shift_factor(450.0, 400.0, 8.86, 101.6);
        assert!(log_at < 0.0);
    }

    #[test]
    fn test_wlf_below_reference_positive() {
        let log_at = WilliamsLandelFerry::shift_factor(350.0, 400.0, 8.86, 101.6);
        assert!(log_at > 0.0);
    }

    // ---- Kauzmann temperature ----

    #[test]
    fn test_kauzmann_temp() {
        let tk = kauzmann_temp(1200.0, 40.0);
        assert!(tk < 1200.0 && tk > 0.0);
    }

    #[test]
    fn test_kauzmann_temp_zero_fragility() {
        let tk = kauzmann_temp(1200.0, 0.0);
        assert_eq!(tk, 1200.0);
    }

    // ---- KWW relaxation ----

    #[test]
    fn test_kww_at_zero_is_one() {
        assert!((StructuralRelaxation::kww_relaxation(0.0, 100.0, 0.7) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_kww_decays_over_time() {
        let phi0 = StructuralRelaxation::kww_relaxation(0.0, 100.0, 0.7);
        let phi1 = StructuralRelaxation::kww_relaxation(100.0, 100.0, 0.7);
        let phi2 = StructuralRelaxation::kww_relaxation(500.0, 100.0, 0.7);
        assert!(phi0 > phi1 && phi1 > phi2);
    }

    #[test]
    fn test_average_relaxation_time_positive() {
        let t_avg = StructuralRelaxation::average_relaxation_time(100.0, 0.7);
        assert!(t_avg > 0.0 && t_avg.is_finite());
    }

    // ---- Glass transition detection ----

    #[test]
    fn test_detect_glass_transition_finds_kink() {
        // Simulate a kink in V(T) at T=800 K
        let mut temps: Vec<f64> = (0..40).map(|i| 1200.0 - i as f64 * 10.0).collect();
        let mut volumes: Vec<f64> = temps
            .iter()
            .map(|&t| {
                if t > 800.0 {
                    0.5e-3 + t * 2.0e-7 // liquid: larger thermal expansion
                } else {
                    0.5e-3 + t * 1.0e-7 // glass: smaller thermal expansion
                }
            })
            .collect();
        // Reverse so temps are decreasing (cooling)
        temps.reverse();
        volumes.reverse();
        let tg = detect_glass_transition(&temps, &volumes);
        assert!(tg.is_some());
        let tg = tg.unwrap();
        assert!(tg > 600.0 && tg < 1000.0, "Tg={tg:.1} not near 800 K");
    }

    #[test]
    fn test_detect_glass_transition_short_array() {
        assert!(detect_glass_transition(&[1.0, 2.0], &[1.0, 2.0]).is_none());
    }

    // ---- Fragility index ----

    #[test]
    fn test_fragility_strong() {
        assert_eq!(fragility_index(20.0), "strong");
    }

    #[test]
    fn test_fragility_intermediate() {
        assert_eq!(fragility_index(50.0), "intermediate");
    }

    #[test]
    fn test_fragility_fragile() {
        assert_eq!(fragility_index(120.0), "fragile");
    }

    // ---- BKS potential ----

    #[test]
    fn test_bks_energy_finite() {
        let e = BksPotential::bks_energy(2.0e-10, 1.8e3, 2.9e10, 1.3e-77, 2.4, -1.2);
        assert!(e.is_finite());
    }

    #[test]
    fn test_bks_energy_zero_distance() {
        let e = BksPotential::bks_energy(0.0, 1.8e3, 2.9e10, 1.3e-77, 2.4, -1.2);
        assert_eq!(e, 0.0);
    }

    // ---- RDF ----

    #[test]
    fn test_rdf_length() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0e-10, 0.0, 0.0], [2.0e-10, 0.0, 0.0]];
        let rdf = GlassStructure::radial_distribution_function(&positions, 5.0e-10, 50);
        assert_eq!(rdf.len(), 50);
    }

    #[test]
    fn test_rdf_empty() {
        let rdf = GlassStructure::radial_distribution_function(&[], 5.0e-10, 10);
        assert!(rdf.iter().all(|&v| v == 0.0));
    }

    // ---- Coordination number ----

    #[test]
    fn test_coordination_number_nonnegative() {
        let rdf = vec![1.0_f64; 100];
        let cn = GlassStructure::coordination_number(&rdf, 0.5e-10, 3.0e-10, 0.1e-10, 1.0e28);
        assert!(cn >= 0.0);
    }

    // ---- Structure factor ----

    #[test]
    fn test_structure_factor_nonnegative() {
        let positions = vec![[0.0, 0.0, 0.0], [3.4e-10, 0.0, 0.0]];
        let sq = GlassStructure::structure_factor_q(&positions, 1.8e10);
        assert!(sq >= 0.0);
    }

    #[test]
    fn test_structure_factor_empty() {
        assert_eq!(GlassStructure::structure_factor_q(&[], 1.0e10), 0.0);
    }

    // ---- MSD ----

    #[test]
    fn test_msd_zero_displacement() {
        let pos = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!((msd_glass(&pos, &pos)).abs() < 1.0e-15);
    }

    #[test]
    fn test_msd_known_value() {
        let p0 = vec![[0.0, 0.0, 0.0]];
        let p1 = vec![[1.0, 0.0, 0.0]];
        assert!((msd_glass(&p0, &p1) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_msd_empty() {
        assert_eq!(msd_glass(&[], &[]), 0.0);
    }

    // ---- Non-Gaussian parameter ----

    #[test]
    fn test_non_gaussian_uniform_motion() {
        // Uniform displacement → all r^4/<r^2>^2 = 1 → α₂ = 3/5 - 1 = -0.4.
        let p0 = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let p1 = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let ng = non_gaussian_param(&p0, &p1);
        assert!((ng - (-0.4)).abs() < 1.0e-10);
    }

    #[test]
    fn test_non_gaussian_empty() {
        assert_eq!(non_gaussian_param(&[], &[]), 0.0);
    }

    // ---- is_glass ----

    #[test]
    fn test_is_glass_true() {
        assert!(is_glass(100.0, 10.0));
    }

    #[test]
    fn test_is_glass_false() {
        assert!(!is_glass(1.0, 10.0));
    }

    #[test]
    fn test_is_glass_equal() {
        assert!(is_glass(10.0, 10.0));
    }

    // ---- GlassSimulation ----

    #[test]
    fn test_simulation_step_changes_positions() {
        let pos = vec![[0.0, 0.0, 0.0], [3.4e-10, 0.0, 0.0]];
        let vel = vec![[100.0, 0.0, 0.0], [-100.0, 0.0, 0.0]];
        let masses = vec![2.8e-26_f64; 2];
        let params = GlassFormingParams::default();
        let mut sim = GlassSimulation::new(
            pos.clone(),
            vel,
            masses,
            2000.0,
            params,
            GlassType::Silica,
            1.0e-15,
        );
        sim.step();
        // After a step, at least one position should differ.
        let changed = sim
            .positions
            .iter()
            .zip(pos.iter())
            .any(|(a, b)| (a[0] - b[0]).abs() > 1.0e-20);
        assert!(changed);
    }

    #[test]
    fn test_glass_type_variants() {
        assert_eq!(GlassType::Silica, GlassType::Silica);
        assert_ne!(GlassType::Metallic, GlassType::Polymer);
    }

    #[test]
    fn test_simulation_anneal_changes_temperature() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let vel = vec![[0.0, 0.0, 0.0]];
        let masses = vec![2.8e-26_f64];
        let params = GlassFormingParams::default();
        let mut sim =
            GlassSimulation::new(pos, vel, masses, 2000.0, params, GlassType::Silica, 1.0e-15);
        sim.anneal(1500.0, 5);
        assert!((sim.temperature - 1500.0).abs() < 1.0e-9);
    }
}
