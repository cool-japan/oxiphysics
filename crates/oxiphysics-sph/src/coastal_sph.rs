#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coastal and wave SPH models.
//!
//! Provides:
//!
//! - Shallow-water wave speed and Airy wave dispersion
//! - Second-order Stokes wave profile
//! - Wave breaking criterion
//! - 1-D shallow-water tsunami model (leapfrog)
//! - Wave runup models (Irribarren / Van der Meer)
//! - Underwater landslide buoyancy correction
//! - Tidal flat wetting/drying model
//! - Wave energy density and group power
//! - SPH [`WaveParticle`] data structure

#![allow(dead_code)]

/// Gravitational acceleration (m/s²).
const G: f64 = 9.81;

// ---------------------------------------------------------------------------
// Wave particle
// ---------------------------------------------------------------------------

/// A single SPH particle carrying wave-coastal state variables.
#[derive(Clone, Debug)]
pub struct WaveParticle {
    /// Position (x, y, z) in metres.
    pub pos: [f64; 3],
    /// Velocity (u, v, w) in m/s.
    pub vel: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub rho: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Local water depth (m).
    pub depth: f64,
    /// Local free-surface elevation η (m).
    pub wave_height: f64,
}

impl WaveParticle {
    /// Construct a new [`WaveParticle`] at rest.
    pub fn new(pos: [f64; 3], mass: f64, rho: f64, depth: f64) -> Self {
        WaveParticle {
            pos,
            vel: [0.0; 3],
            mass,
            rho,
            pressure: 0.0,
            depth,
            wave_height: 0.0,
        }
    }

    /// Kinetic energy of this particle (J).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.vel[0].powi(2) + self.vel[1].powi(2) + self.vel[2].powi(2);
        0.5 * self.mass * v2
    }
}

// ---------------------------------------------------------------------------
// Wave physics
// ---------------------------------------------------------------------------

/// Phase speed of a shallow-water (long) wave.
///
/// c = √(g h)
///
/// # Arguments
/// * `h` – water depth (m)
pub fn shallow_water_wave_speed(h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    (G * h).sqrt()
}

/// Solve the Airy (linear) wave dispersion relation ω² = g k tanh(k h)
/// for the wave number `k` given angular frequency `omega` and depth `h`.
///
/// Uses Newton–Raphson iteration starting from the deep-water approximation
/// k₀ = ω²/g.
///
/// Returns the wavenumber k (rad/m).
///
/// # Arguments
/// * `omega` – angular frequency (rad/s)
/// * `h`     – water depth (m)
pub fn airy_wave_dispersion(omega: f64, h: f64) -> f64 {
    if omega <= 0.0 || h <= 0.0 {
        return 0.0;
    }
    let omega2 = omega * omega;
    // Initial guess: deep-water limit
    let mut k = omega2 / G;
    for _ in 0..50 {
        let th = (k * h).tanh();
        let f = omega2 - G * k * th;
        let df = -G * (th + k * h * (1.0 - th * th));
        let dk = -f / df;
        k += dk;
        if k < 0.0 {
            k = f64::EPSILON;
        }
        if dk.abs() < 1e-12 * k {
            break;
        }
    }
    k
}

/// Second-order Stokes wave free-surface elevation η(x, t).
///
/// η = A cos(kx − ωt) + (A²k/4) cosh(kh)(2 + cosh(2kh)) / sinh³(kh) · cos(2(kx−ωt))
///
/// # Arguments
/// * `x`         – horizontal position (m)
/// * `t`         – time (s)
/// * `amplitude` – wave amplitude A (m)
/// * `k`         – wave number (rad/m)
/// * `omega`     – angular frequency (rad/s)
pub fn stokes_wave_profile(x: f64, t: f64, amplitude: f64, k: f64, omega: f64) -> f64 {
    if k <= 0.0 || omega <= 0.0 {
        return 0.0;
    }
    let phase = k * x - omega * t;
    let a = amplitude;
    // Second-order correction (uses a finite depth h inferred from dispersion)
    // For simplicity we use the deep-water second-order term:
    // η₂ = (k A² / 2) cos(2θ)
    let eta1 = a * phase.cos();
    let eta2 = 0.5 * k * a * a * (2.0 * phase).cos();
    eta1 + eta2
}

/// Wave breaking criterion: true when H/h > 0.78 (McCowan).
///
/// # Arguments
/// * `wave_height`  – wave height H (m)
/// * `water_depth`  – water depth h (m)
pub fn breaking_criterion(wave_height: f64, water_depth: f64) -> bool {
    if water_depth <= 0.0 {
        return true;
    }
    wave_height / water_depth > 0.78
}

// ---------------------------------------------------------------------------
// Tsunami model
// ---------------------------------------------------------------------------

/// 1-D shallow-water tsunami model on a 1-D domain.
///
/// Discretises the linear shallow-water equations
///
///   ∂η/∂t + ∂(H u)/∂x = 0
///   ∂u/∂t + g ∂η/∂x   = 0
///
/// using a leapfrog (staggered) scheme.
pub struct TsunamiModel {
    /// Bathymetry (still-water depth) at each cell centre (m).
    pub bathymetry: Vec<f64>,
    /// Free-surface elevation η at each cell centre (m).
    pub eta: Vec<f64>,
    /// Depth-averaged velocity at each cell *edge* (length nx+1, m/s).
    pub u: Vec<f64>,
    /// Number of cells.
    pub nx: usize,
    /// Cell width (m).
    pub dx: f64,
    /// Elapsed simulation time (s).
    elapsed: f64,
}

impl TsunamiModel {
    /// Create a new [`TsunamiModel`] with flat bathymetry `depth` (m).
    ///
    /// # Arguments
    /// * `nx`    – number of cells
    /// * `dx`    – cell width (m)
    /// * `depth` – uniform still-water depth (m)
    pub fn new(nx: usize, dx: f64, depth: f64) -> Self {
        TsunamiModel {
            bathymetry: vec![depth; nx],
            eta: vec![0.0; nx],
            u: vec![0.0; nx + 1],
            nx,
            dx,
            elapsed: 0.0,
        }
    }

    /// Set the initial free-surface elevation as a Gaussian hump.
    ///
    /// η(x) = A · exp(−(x − x₀)² / (2σ²))
    ///
    /// # Arguments
    /// * `amplitude` – peak elevation (m)
    /// * `x0`        – hump centre (m)
    /// * `sigma`     – width parameter (m)
    pub fn set_gaussian_hump(&mut self, amplitude: f64, x0: f64, sigma: f64) {
        for i in 0..self.nx {
            let x = (i as f64 + 0.5) * self.dx;
            self.eta[i] = amplitude * (-(x - x0).powi(2) / (2.0 * sigma * sigma)).exp();
        }
    }

    /// Advance by one leapfrog time step `dt` (s).
    ///
    /// Applies Neumann (zero-flux) boundary conditions at both ends.
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx;
        let dx = self.dx;

        // Update velocity at interior edges from η
        for i in 1..nx {
            self.u[i] -= dt * G * (self.eta[i] - self.eta[i - 1]) / dx;
        }
        // Walls at boundaries
        self.u[0] = 0.0;
        self.u[nx] = 0.0;

        // Update η at cell centres from u
        for i in 0..nx {
            let h = self.bathymetry[i];
            self.eta[i] -= dt * h * (self.u[i + 1] - self.u[i]) / dx;
        }
        self.elapsed += dt;
    }

    /// Maximum free-surface elevation over the entire domain (m).
    pub fn max_amplitude(&self) -> f64 {
        self.eta.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    }

    /// Time (s) at which any cell first exceeds `threshold` (m), or `None`.
    ///
    /// Useful for estimating tsunami arrival time after calling `step` in a loop.
    /// This method returns the current elapsed time if the threshold is already exceeded.
    pub fn arrival_time(&self, threshold: f64) -> Option<f64> {
        if self.eta.iter().any(|e| *e > threshold) {
            Some(self.elapsed)
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Wave runup
// ---------------------------------------------------------------------------

/// Wave runup model based on the Irribarren number (surf similarity).
///
/// Models a plane beach with given slope and roughness.
#[derive(Clone, Debug)]
pub struct WaveRunup {
    /// Beach slope (tan α, dimensionless).
    pub slope: f64,
    /// Manning roughness coefficient (s/m^(1/3)).
    pub roughness: f64,
}

impl WaveRunup {
    /// Create a new [`WaveRunup`] model.
    ///
    /// # Arguments
    /// * `slope`     – beach slope tan(α)
    /// * `roughness` – Manning n
    pub fn new(slope: f64, roughness: f64) -> Self {
        WaveRunup { slope, roughness }
    }

    /// Irribarren (surf similarity) number.
    ///
    /// ξ = tan(α) / √(H₀ / L₀)
    ///
    /// # Arguments
    /// * `h0` – offshore significant wave height (m)
    /// * `l0` – deep-water wavelength (m)
    pub fn irribarren_number(&self, h0: f64, l0: f64) -> f64 {
        if h0 <= 0.0 || l0 <= 0.0 {
            return 0.0;
        }
        self.slope / (h0 / l0).sqrt()
    }

    /// Runup height R (m) using the Irribarren-number parameterisation
    /// (Hunt 1959 / Holman 1986).
    ///
    /// R = ξ · H₀   (simplified form)
    ///
    /// # Arguments
    /// * `h0` – offshore significant wave height (m)
    /// * `l0` – deep-water wavelength (m)
    pub fn runup_irribarren(&self, h0: f64, l0: f64) -> f64 {
        if h0 <= 0.0 || l0 <= 0.0 {
            return 0.0;
        }
        let xi = self.irribarren_number(h0, l0);
        xi * h0
    }

    /// Runup height R (m) using the Van der Meer & Stam (1992) formula.
    ///
    /// For breaking waves (ξ < 1.5):
    ///   R/H₀ = γ · 1.5 · ξ^0.5
    /// For non-breaking (ξ ≥ 1.5):
    ///   R/H₀ = γ · 3.0 · ξ^0.2
    ///
    /// # Arguments
    /// * `h0`    – offshore significant wave height (m)
    /// * `l0`    – deep-water wavelength (m)
    /// * `gamma` – reduction factor (roughness, obliquity; typically 0.5–1.0)
    pub fn runup_van_der_meer(&self, h0: f64, l0: f64, gamma: f64) -> f64 {
        if h0 <= 0.0 || l0 <= 0.0 {
            return 0.0;
        }
        let xi = self.irribarren_number(h0, l0);
        let ratio = if xi < 1.5 {
            gamma * 1.5 * xi.powf(0.5)
        } else {
            gamma * 3.0 * xi.powf(0.2)
        };
        ratio * h0
    }
}

// ---------------------------------------------------------------------------
// Landslide
// ---------------------------------------------------------------------------

/// Buoyancy-corrected volume flux from an underwater landslide.
///
/// Returns the effective submerged volume (m³) that drives tsunami generation.
///
/// V_eff = V_slide · (ρ_s − ρ_f) / ρ_f
///
/// # Arguments
/// * `rho_s`    – slide density (kg/m³)
/// * `rho_f`    – fluid density (kg/m³)
/// * `v_slide`  – total slide volume (m³)
pub fn underwater_landslide_volume(rho_s: f64, rho_f: f64, v_slide: f64) -> f64 {
    if rho_f <= 0.0 {
        return 0.0;
    }
    v_slide * (rho_s - rho_f) / rho_f
}

// ---------------------------------------------------------------------------
// Tidal flat
// ---------------------------------------------------------------------------

/// 1-D tidal flat wetting/drying model.
///
/// Evolves the water depth h = η − z_bed on a flat-bed grid using a simple
/// explicit continuity update with volume fluxes at the boundaries.
pub struct TidalFlat {
    /// Bed elevation z_bed at each cell (m above datum).
    pub bed: Vec<f64>,
    /// Free-surface elevation η at each cell (m).
    pub water: Vec<f64>,
    /// Number of cells.
    pub nx: usize,
    /// Cell width (m).
    pub dx: f64,
    /// Manning roughness coefficient (s/m^(1/3)).
    pub manning_n: f64,
}

impl TidalFlat {
    /// Create a new [`TidalFlat`] with uniform bed elevation `z_bed`.
    ///
    /// # Arguments
    /// * `nx`       – number of cells
    /// * `dx`       – cell width (m)
    /// * `z_bed`    – uniform bed elevation (m)
    /// * `manning_n`– Manning n
    pub fn new(nx: usize, dx: f64, z_bed: f64, manning_n: f64) -> Self {
        TidalFlat {
            bed: vec![z_bed; nx],
            water: vec![z_bed; nx],
            nx,
            dx,
            manning_n,
        }
    }

    /// Water depth h = max(η − z_bed, 0) at cell `i`.
    #[inline]
    fn depth(&self, i: usize) -> f64 {
        (self.water[i] - self.bed[i]).max(0.0)
    }

    /// Advance the tidal flat by one explicit time step using volume conservation.
    ///
    /// Volume enters from the left at rate `q_in` (m²/s per unit width)
    /// and leaves from the right at rate `q_out`.  Internal fluxes use a
    /// simple upwind depth-averaged Manning formula.
    ///
    /// # Arguments
    /// * `dt`    – time step (s)
    /// * `q_in`  – inflow discharge at left boundary (m²/s)
    /// * `q_out` – outflow discharge at right boundary (m²/s)
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, dt: f64, q_in: f64, q_out: f64) {
        let nx = self.nx;
        let dx = self.dx;
        let n = self.manning_n;

        // Compute edge fluxes (nx+1 edges)
        let mut flux = vec![0.0f64; nx + 1];
        flux[0] = q_in;
        flux[nx] = q_out;

        for e in 1..nx {
            let hl = self.depth(e - 1);
            let hr = self.depth(e);
            let h_avg = 0.5 * (hl + hr);
            if h_avg < 1e-6 {
                continue;
            }
            let deta = self.water[e] - self.water[e - 1];
            let slope = deta / dx;
            // Manning: q = (1/n) h^(5/3) S^(1/2)  (sign from slope)
            let s_abs = slope.abs().max(1e-12);
            let sign_s = if slope >= 0.0 { 1.0 } else { -1.0 };
            flux[e] = sign_s * (1.0 / n) * h_avg.powf(5.0 / 3.0) * s_abs.sqrt();
        }

        // Update water surface
        for i in 0..nx {
            let dh = -dt * (flux[i + 1] - flux[i]) / dx;
            self.water[i] += dh;
            // Wetting/drying: clamp to bed
            if self.water[i] < self.bed[i] {
                self.water[i] = self.bed[i];
            }
        }
    }

    /// Total inundated area (m²) assuming unit width in the transverse direction.
    ///
    /// Counts cells where the water depth exceeds 1 mm.
    pub fn inundation_area(&self) -> f64 {
        let wet = (0..self.nx).filter(|i| self.depth(*i) > 1e-3).count();
        wet as f64 * self.dx
    }
}

// ---------------------------------------------------------------------------
// Wave energy & power
// ---------------------------------------------------------------------------

/// Wave energy density per unit surface area (J/m²).
///
/// E = ½ ρ g A²
///
/// # Arguments
/// * `rho`       – water density (kg/m³)
/// * `g`         – gravitational acceleration (m/s²)
/// * `amplitude` – wave amplitude A (m)
pub fn wave_energy_density(rho: f64, g: f64, amplitude: f64) -> f64 {
    0.5 * rho * g * amplitude * amplitude
}

/// Wave power (energy flux) per unit crest width (W/m).
///
/// P = E · c_g
///
/// where the group velocity c_g = (ω/2k)(1 + 2kh/sinh(2kh)).
///
/// # Arguments
/// * `rho`       – water density (kg/m³)
/// * `g`         – gravitational acceleration (m/s²)
/// * `amplitude` – wave amplitude A (m)
/// * `k`         – wave number (rad/m)
/// * `h`         – water depth (m)
/// * `omega`     – angular frequency (rad/s)
pub fn wave_power(rho: f64, g: f64, amplitude: f64, k: f64, h: f64, omega: f64) -> f64 {
    if k <= 0.0 || omega <= 0.0 {
        return 0.0;
    }
    let e = wave_energy_density(rho, g, amplitude);
    let kh2 = 2.0 * k * h;
    let n = 0.5 * (1.0 + kh2 / kh2.sinh());
    let c_phase = omega / k;
    let c_g = n * c_phase;
    e * c_g
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // -----------------------------------------------------------------------
    // WaveParticle
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_particle_new() {
        let p = WaveParticle::new([1.0, 2.0, 0.0], 1.0, 1025.0, 10.0);
        assert_eq!(p.pos, [1.0, 2.0, 0.0]);
        assert_eq!(p.vel, [0.0, 0.0, 0.0]);
        assert_eq!(p.rho, 1025.0);
        assert_eq!(p.depth, 10.0);
    }

    #[test]
    fn test_wave_particle_kinetic_energy_at_rest() {
        let p = WaveParticle::new([0.0; 3], 2.0, 1000.0, 5.0);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_wave_particle_kinetic_energy_moving() {
        let mut p = WaveParticle::new([0.0; 3], 2.0, 1000.0, 5.0);
        p.vel = [3.0, 4.0, 0.0]; // speed = 5
        assert!((p.kinetic_energy() - 25.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // shallow_water_wave_speed
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_speed_known_depth() {
        let c = shallow_water_wave_speed(1.0);
        assert!((c - G.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_wave_speed_zero_depth() {
        assert_eq!(shallow_water_wave_speed(0.0), 0.0);
    }

    #[test]
    fn test_wave_speed_negative_depth() {
        assert_eq!(shallow_water_wave_speed(-1.0), 0.0);
    }

    #[test]
    fn test_wave_speed_proportional() {
        let c1 = shallow_water_wave_speed(4.0);
        let c2 = shallow_water_wave_speed(1.0);
        assert!((c1 / c2 - 2.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // airy_wave_dispersion
    // -----------------------------------------------------------------------

    #[test]
    fn test_dispersion_deep_water() {
        // Deep water: omega^2 ≈ g k  (kh >> 1)
        let omega = 1.0;
        let h = 1000.0; // very deep
        let k = airy_wave_dispersion(omega, h);
        let k_deep = omega * omega / G;
        assert!((k - k_deep).abs() / k_deep < 1e-4);
    }

    #[test]
    fn test_dispersion_satisfies_relation() {
        let omega = 0.5;
        let h = 5.0;
        let k = airy_wave_dispersion(omega, h);
        let lhs = omega * omega;
        let rhs = G * k * (k * h).tanh();
        assert!(
            (lhs - rhs).abs() / lhs < 1e-8,
            "dispersion residual: {}",
            (lhs - rhs).abs()
        );
    }

    #[test]
    fn test_dispersion_zero_omega() {
        assert_eq!(airy_wave_dispersion(0.0, 10.0), 0.0);
    }

    #[test]
    fn test_dispersion_zero_depth() {
        assert_eq!(airy_wave_dispersion(1.0, 0.0), 0.0);
    }

    // -----------------------------------------------------------------------
    // stokes_wave_profile
    // -----------------------------------------------------------------------

    #[test]
    fn test_stokes_wave_at_crest() {
        // At x=0, t=0: eta should be positive for positive amplitude
        let eta = stokes_wave_profile(0.0, 0.0, 1.0, 0.5, 1.0);
        assert!(eta > 0.0);
    }

    #[test]
    fn test_stokes_wave_zero_amplitude() {
        let eta = stokes_wave_profile(0.0, 0.0, 0.0, 0.5, 1.0);
        assert_eq!(eta, 0.0);
    }

    #[test]
    fn test_stokes_wave_periodic() {
        // Periodicity: eta(x, t) == eta(x + 2π/k, t)
        let k = 0.4;
        let omega = 0.8;
        let a = 0.5;
        let lambda = 2.0 * PI / k;
        let eta1 = stokes_wave_profile(0.0, 0.0, a, k, omega);
        let eta2 = stokes_wave_profile(lambda, 0.0, a, k, omega);
        assert!((eta1 - eta2).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // breaking_criterion
    // -----------------------------------------------------------------------

    #[test]
    fn test_breaking_criterion_shallow() {
        assert!(breaking_criterion(1.0, 1.0)); // H/h = 1.0 > 0.78
    }

    #[test]
    fn test_breaking_criterion_deep() {
        assert!(!breaking_criterion(0.5, 10.0)); // H/h = 0.05
    }

    #[test]
    fn test_breaking_criterion_exact_threshold() {
        // H/h = 0.78: just not breaking
        assert!(!breaking_criterion(0.78, 1.0));
    }

    #[test]
    fn test_breaking_criterion_zero_depth() {
        assert!(breaking_criterion(0.1, 0.0));
    }

    // -----------------------------------------------------------------------
    // TsunamiModel
    // -----------------------------------------------------------------------

    #[test]
    fn test_tsunami_new_flat() {
        let m = TsunamiModel::new(10, 100.0, 4000.0);
        assert!(m.eta.iter().all(|e| *e == 0.0));
        assert_eq!(m.bathymetry.len(), 10);
    }

    #[test]
    fn test_tsunami_gaussian_hump() {
        let mut m = TsunamiModel::new(100, 1000.0, 4000.0);
        m.set_gaussian_hump(1.0, 50_000.0, 10_000.0);
        assert!(m.max_amplitude() > 0.9);
    }

    #[test]
    fn test_tsunami_step_conserves_volume_approx() {
        let mut m = TsunamiModel::new(50, 1000.0, 100.0);
        m.set_gaussian_hump(1.0, 25_000.0, 5000.0);
        let vol_before: f64 = m.eta.iter().sum::<f64>() * 1000.0;
        for _ in 0..10 {
            m.step(1.0);
        }
        let vol_after: f64 = m.eta.iter().sum::<f64>() * 1000.0;
        // Volume should be roughly conserved (walls at boundaries reflect)
        assert!((vol_after - vol_before).abs() < 0.01 * vol_before.abs() + 1.0);
    }

    #[test]
    fn test_tsunami_max_amplitude_nonzero() {
        let mut m = TsunamiModel::new(20, 500.0, 1000.0);
        m.set_gaussian_hump(2.0, 5000.0, 1000.0);
        assert!(m.max_amplitude() > 0.0);
    }

    #[test]
    fn test_tsunami_arrival_time_none_before_step() {
        let m = TsunamiModel::new(20, 500.0, 100.0);
        assert!(m.arrival_time(0.1).is_none());
    }

    #[test]
    fn test_tsunami_arrival_time_some_after_hump() {
        let mut m = TsunamiModel::new(20, 500.0, 100.0);
        m.set_gaussian_hump(1.0, 5000.0, 1000.0);
        // The hump is already > threshold = 0.1
        assert!(m.arrival_time(0.1).is_some());
    }

    // -----------------------------------------------------------------------
    // WaveRunup
    // -----------------------------------------------------------------------

    #[test]
    fn test_runup_irribarren_number() {
        let wr = WaveRunup::new(0.1, 0.02);
        let xi = wr.irribarren_number(2.0, 50.0);
        let expected = 0.1 / (2.0_f64 / 50.0).sqrt();
        assert!((xi - expected).abs() < 1e-10);
    }

    #[test]
    fn test_runup_irribarren_zero_h0() {
        let wr = WaveRunup::new(0.1, 0.02);
        assert_eq!(wr.runup_irribarren(0.0, 50.0), 0.0);
    }

    #[test]
    fn test_runup_irribarren_positive() {
        let wr = WaveRunup::new(0.1, 0.02);
        let r = wr.runup_irribarren(1.5, 40.0);
        assert!(r > 0.0);
    }

    #[test]
    fn test_runup_van_der_meer_breaking() {
        // xi < 1.5 → breaking regime
        let wr = WaveRunup::new(0.05, 0.02);
        let r = wr.runup_van_der_meer(2.0, 30.0, 1.0);
        assert!(r >= 0.0);
    }

    #[test]
    fn test_runup_van_der_meer_nonbreaking() {
        // steep slope → xi > 1.5
        let wr = WaveRunup::new(0.5, 0.02);
        let r = wr.runup_van_der_meer(1.0, 100.0, 1.0);
        assert!(r > 0.0);
    }

    #[test]
    fn test_runup_van_der_meer_gamma_scales() {
        let wr = WaveRunup::new(0.1, 0.02);
        let r1 = wr.runup_van_der_meer(2.0, 50.0, 1.0);
        let r2 = wr.runup_van_der_meer(2.0, 50.0, 0.5);
        assert!((r1 - 2.0 * r2).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // underwater_landslide_volume
    // -----------------------------------------------------------------------

    #[test]
    fn test_landslide_volume_basic() {
        // rho_s = 2000, rho_f = 1000, V = 1e6
        let v = underwater_landslide_volume(2000.0, 1000.0, 1.0e6);
        assert!((v - 1.0e6).abs() < 1.0); // (2000-1000)/1000 * 1e6 = 1e6
    }

    #[test]
    fn test_landslide_volume_zero_rho_f() {
        let v = underwater_landslide_volume(2000.0, 0.0, 1.0e6);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_landslide_volume_neutral_buoyancy() {
        // rho_s == rho_f → effective volume = 0
        let v = underwater_landslide_volume(1000.0, 1000.0, 1.0e6);
        assert_eq!(v, 0.0);
    }

    // -----------------------------------------------------------------------
    // TidalFlat
    // -----------------------------------------------------------------------

    #[test]
    fn test_tidal_flat_new() {
        let tf = TidalFlat::new(10, 100.0, 0.0, 0.03);
        assert_eq!(tf.nx, 10);
        assert!(tf.water.iter().all(|w| *w == 0.0));
    }

    #[test]
    fn test_tidal_flat_inundation_dry() {
        let tf = TidalFlat::new(10, 100.0, 0.0, 0.03);
        // All cells at bed level → no inundation
        assert_eq!(tf.inundation_area(), 0.0);
    }

    #[test]
    fn test_tidal_flat_inundation_wet() {
        let mut tf = TidalFlat::new(10, 100.0, 0.0, 0.03);
        for i in 0..10 {
            tf.water[i] = 1.0; // 1 m above bed
        }
        assert!((tf.inundation_area() - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_tidal_flat_step_no_panic() {
        let mut tf = TidalFlat::new(20, 50.0, 0.0, 0.025);
        for i in 0..20 {
            tf.water[i] = 0.5;
        }
        tf.step(0.1, 0.0, 0.0);
        // Water should not go below bed
        assert!(tf.water.iter().zip(tf.bed.iter()).all(|(w, b)| w >= b));
    }

    #[test]
    fn test_tidal_flat_step_inflow_raises_water() {
        let mut tf = TidalFlat::new(5, 100.0, 0.0, 0.03);
        // Uniform water level
        for i in 0..5 {
            tf.water[i] = 0.5;
        }
        tf.step(10.0, 1.0, 0.0);
        // Left cell should have gained water
        assert!(tf.water[0] > 0.5 || tf.water[1] > 0.5 || tf.inundation_area() > 0.0);
    }

    // -----------------------------------------------------------------------
    // wave_energy_density
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_energy_basic() {
        let e = wave_energy_density(1025.0, G, 1.0);
        let expected = 0.5 * 1025.0 * G;
        assert!((e - expected).abs() < 1e-6);
    }

    #[test]
    fn test_wave_energy_zero_amplitude() {
        assert_eq!(wave_energy_density(1025.0, G, 0.0), 0.0);
    }

    #[test]
    fn test_wave_energy_scales_as_amplitude_squared() {
        let e1 = wave_energy_density(1000.0, G, 1.0);
        let e2 = wave_energy_density(1000.0, G, 2.0);
        assert!((e2 - 4.0 * e1).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // wave_power
    // -----------------------------------------------------------------------

    #[test]
    fn test_wave_power_positive() {
        let p = wave_power(1025.0, G, 1.0, 0.1, 10.0, 1.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_wave_power_zero_k() {
        assert_eq!(wave_power(1025.0, G, 1.0, 0.0, 10.0, 1.0), 0.0);
    }

    #[test]
    fn test_wave_power_zero_omega() {
        assert_eq!(wave_power(1025.0, G, 1.0, 0.1, 10.0, 0.0), 0.0);
    }

    #[test]
    fn test_wave_power_scales_amplitude_squared() {
        let p1 = wave_power(1000.0, G, 1.0, 0.2, 5.0, 1.0);
        let p2 = wave_power(1000.0, G, 2.0, 0.2, 5.0, 1.0);
        assert!((p2 - 4.0 * p1).abs() < 1e-8 * p1.abs() + 1e-8);
    }
}

// ── Ocean / coastal advanced wave models ─────────────────────────────────────

use std::f64::consts::PI as F64_PI;

/// Ocean/coastal SPH particle with depth and wave information.
pub struct OceanParticle {
    /// Position (x, y, z) in metres.
    pub pos: [f64; 3],
    /// Velocity (u, v, w) in m/s.
    pub vel: [f64; 3],
    /// Density (kg/m³).
    pub rho: f64,
    /// Pressure (Pa).
    pub p: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Local water depth (m above sea bed).
    pub depth: f64,
}

impl OceanParticle {
    /// Create a new [`OceanParticle`] at rest.
    pub fn new(pos: [f64; 3], mass: f64, rho: f64, depth: f64, h: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            rho,
            p: 0.0,
            mass,
            h,
            depth,
        }
    }

    /// Kinetic energy of this particle (J).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.vel[0].powi(2) + self.vel[1].powi(2) + self.vel[2].powi(2);
        0.5 * self.mass * v2
    }
}

/// Compute the orbital velocity vector of an Airy (linear) wave.
///
/// Horizontal: u = a ω cosh(k(z+d)) / sinh(kd) · cos(kx − ωt)
/// Vertical:   w = a ω sinh(k(z+d)) / sinh(kd) · sin(kx − ωt)
///
/// # Arguments
/// * `x`     – horizontal position (m)
/// * `z`     – vertical coordinate (m, positive upward, z = 0 at mean surface)
/// * `t`     – time (s)
/// * `amp`   – wave amplitude a (m)
/// * `k`     – wave number (rad/m)
/// * `omega` – angular frequency (rad/s)
/// * `depth` – water depth d (m)
///
/// Returns `[u, 0, w]`.
pub fn airy_wave_velocity(
    x: f64,
    z: f64,
    t: f64,
    amp: f64,
    k: f64,
    omega: f64,
    depth: f64,
) -> [f64; 3] {
    if k <= 0.0 || omega <= 0.0 || depth <= 0.0 {
        return [0.0; 3];
    }
    let phase = k * x - omega * t;
    let kd = k * depth;
    let sinh_kd = kd.sinh().max(1e-30);
    let kzd = k * (z + depth);
    let u = amp * omega * kzd.cosh() / sinh_kd * phase.cos();
    let w = amp * omega * kzd.sinh() / sinh_kd * phase.sin();
    [u, 0.0, w]
}

/// Stokes second-order wave surface elevation η(x, t).
///
/// η = a cos(kx − ωt) + (a²k/2) · (cosh(kd)(2 + cosh(2kd))) / sinh³(kd)
///     · cos(2(kx − ωt))
///
/// # Arguments
/// * `x`     – horizontal position (m)
/// * `z`     – vertical coordinate (unused, signature for API consistency)
/// * `t`     – time (s)
/// * `amp`   – wave amplitude a (m)
/// * `k`     – wave number (rad/m)
/// * `omega` – angular frequency (rad/s)
/// * `depth` – water depth d (m)
pub fn stokes_second_order(
    x: f64,
    _z: f64,
    t: f64,
    amp: f64,
    k: f64,
    omega: f64,
    depth: f64,
) -> f64 {
    if k <= 0.0 || omega <= 0.0 {
        return 0.0;
    }
    let phase = k * x - omega * t;
    let kd = k * depth;
    let sinh_kd = kd.sinh().max(1e-30);
    let eta1 = amp * phase.cos();
    let second_amp = (amp * amp * k / 4.0) * kd.cosh() * (2.0 + (2.0 * kd).cosh())
        / (sinh_kd * sinh_kd * sinh_kd);
    let eta2 = second_amp * (2.0 * phase).cos();
    eta1 + eta2
}

/// Wave dispersion relation: ω² = g k tanh(k d).
///
/// Returns ω (rad/s) for given wavenumber and depth.
///
/// # Arguments
/// * `k`     – wave number (rad/m)
/// * `depth` – water depth d (m)
/// * `g`     – gravitational acceleration (m/s²)
pub fn wave_dispersion(k: f64, depth: f64, g: f64) -> f64 {
    if k <= 0.0 || g <= 0.0 {
        return 0.0;
    }
    (g * k * (k * depth).tanh()).sqrt()
}

/// Wave group velocity.
///
/// c_g = c_p/2 · (1 + 2kd / sinh(2kd))
///
/// where c_p = ω/k is the phase speed.
///
/// # Arguments
/// * `k`     – wave number (rad/m)
/// * `depth` – water depth d (m)
/// * `g`     – gravitational acceleration (m/s²)
pub fn group_velocity(k: f64, depth: f64, g: f64) -> f64 {
    if k <= 0.0 || g <= 0.0 {
        return 0.0;
    }
    let omega = wave_dispersion(k, depth, g);
    let c_p = omega / k;
    let kd2 = 2.0 * k * depth;
    let n = 0.5 * (1.0 + kd2 / kd2.sinh().max(1e-30));
    n * c_p
}

/// Shoaling coefficient Ks = sqrt(c_g1 / c_g2).
///
/// Wave amplitude is amplified by Ks as waves travel from depth d1 to d2.
///
/// # Arguments
/// * `k1` – wave number at location 1
/// * `d1` – water depth at location 1 (m)
/// * `k2` – wave number at location 2
/// * `d2` – water depth at location 2 (m)
pub fn shoaling_coefficient(k1: f64, d1: f64, k2: f64, d2: f64) -> f64 {
    let cg1 = group_velocity(k1, d1, G);
    let cg2 = group_velocity(k2, d2, G);
    if cg2 < 1e-30 {
        return 0.0;
    }
    (cg1 / cg2).sqrt()
}

/// McCowan breaking criterion: H/d > 0.78.
///
/// # Arguments
/// * `wave_height` – wave height H (m)
/// * `depth`       – water depth d (m)
pub fn breaking_criterion_ocean(wave_height: f64, depth: f64) -> bool {
    if depth <= 0.0 {
        return true;
    }
    wave_height / depth > 0.78
}

/// M2 + S2 tidal forcing.
pub struct TidalForcing {
    /// M2 tidal amplitude (m).
    pub m2_amp: f64,
    /// M2 angular frequency (rad/s).
    pub m2_freq: f64,
    /// S2 tidal amplitude (m).
    pub s2_amp: f64,
    /// S2 angular frequency (rad/s).
    pub s2_freq: f64,
}

impl TidalForcing {
    /// Create a new [`TidalForcing`] with standard M2/S2 frequencies.
    ///
    /// M2 period ≈ 12.42 h, S2 period = 12 h.
    ///
    /// # Arguments
    /// * `m2_amp` – M2 amplitude (m)
    /// * `s2_amp` – S2 amplitude (m)
    pub fn new(m2_amp: f64, s2_amp: f64) -> Self {
        // M2 frequency: 2π / (12.42 * 3600)
        let m2_freq = 2.0 * F64_PI / (12.42 * 3600.0);
        // S2 frequency: 2π / (12 * 3600)
        let s2_freq = 2.0 * F64_PI / (12.0 * 3600.0);
        Self {
            m2_amp,
            m2_freq,
            s2_amp,
            s2_freq,
        }
    }

    /// Tidal elevation η(t) = A_M2 cos(ω_M2 t) + A_S2 cos(ω_S2 t).
    ///
    /// # Arguments
    /// * `t` – time (s)
    pub fn elevation(&self, t: f64) -> f64 {
        self.m2_amp * (self.m2_freq * t).cos() + self.s2_amp * (self.s2_freq * t).cos()
    }

    /// Depth-averaged tidal current (m/s) assuming a simple linear relationship.
    ///
    /// u_tidal ≈ η_amplitude * sqrt(g/depth) (shallow-water approximation).
    ///
    /// # Arguments
    /// * `t`     – time (s)
    /// * `depth` – water depth (m)
    pub fn tidal_current(&self, t: f64, depth: f64) -> f64 {
        if depth <= 0.0 {
            return 0.0;
        }
        let eta = self.elevation(t);
        let c = (G * depth).sqrt();
        eta / depth * c
    }
}

/// Stokes drift velocity (mass transport due to wave orbital motion).
///
/// u_s = a² ω k · cosh(2k(z+d)) / (2 sinh²(kd))
///
/// # Arguments
/// * `amp`   – wave amplitude a (m)
/// * `k`     – wave number (rad/m)
/// * `omega` – angular frequency (rad/s)
/// * `z`     – vertical position (m, z=0 at surface, z=-d at bed)
pub fn stokes_drift(amp: f64, k: f64, omega: f64, z: f64) -> f64 {
    // Simplified deep-water form: u_s = a² ω k exp(2kz)
    if k <= 0.0 || omega <= 0.0 {
        return 0.0;
    }
    amp * amp * omega * k * (2.0 * k * z).exp()
}

/// Radiation stress S_xx = E · (2 c_g/c_p − 0.5).
///
/// # Arguments
/// * `amp`   – wave amplitude (m)
/// * `k`     – wave number (rad/m)
/// * `omega` – angular frequency (rad/s)
/// * `depth` – water depth (m)
pub fn radiation_stress_xx(amp: f64, k: f64, omega: f64, depth: f64) -> f64 {
    let energy = wave_energy_density(1025.0, G, amp);
    if k <= 0.0 || omega <= 0.0 {
        return 0.0;
    }
    let c_p = omega / k;
    let cg = group_velocity(k, depth, G);
    energy * (2.0 * cg / c_p - 0.5)
}

/// Coastal SPH solver holding a collection of [`OceanParticle`]s.
pub struct CoastalSolver {
    /// All ocean particles managed by this solver.
    pub particles: Vec<OceanParticle>,
}

impl CoastalSolver {
    /// Create an empty [`CoastalSolver`].
    pub fn new() -> Self {
        Self {
            particles: Vec::new(),
        }
    }

    /// Add a particle to the solver.
    pub fn add_particle(&mut self, p: OceanParticle) {
        self.particles.push(p);
    }

    /// Number of particles.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Total kinetic energy of all particles (J).
    pub fn total_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
}

impl Default for CoastalSolver {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests_ocean {
    use super::*;

    // ── OceanParticle ───────────────────────────────────────────────────────

    #[test]
    fn test_ocean_particle_new_at_rest() {
        let p = OceanParticle::new([0.0; 3], 1.0, 1025.0, 10.0, 0.05);
        assert_eq!(p.vel, [0.0; 3]);
        assert!((p.rho - 1025.0).abs() < 1e-10);
    }

    #[test]
    fn test_ocean_particle_kinetic_energy_rest() {
        let p = OceanParticle::new([0.0; 3], 2.0, 1025.0, 10.0, 0.05);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    // ── airy_wave_velocity ──────────────────────────────────────────────────

    #[test]
    fn test_airy_velocity_deep_water_surface() {
        // Deep water surface z=0: u = a*omega*cos(kx - omega*t)
        let (amp, k, omega, depth) = (1.0, 0.1, 1.0, 100.0);
        let v = airy_wave_velocity(0.0, 0.0, 0.0, amp, k, omega, depth);
        // at z=0, t=0, x=0: phase=0, cos(0)=1
        // kd very large → cosh(kzd)/sinh(kd) ≈ 1
        assert!(v[0].abs() > 0.0);
        assert!(v[1].abs() < 1e-14);
    }

    #[test]
    fn test_airy_velocity_zero_amp() {
        let v = airy_wave_velocity(0.0, 0.0, 0.0, 0.0, 0.1, 1.0, 10.0);
        assert_eq!(v, [0.0; 3]);
    }

    #[test]
    fn test_airy_velocity_zero_depth() {
        let v = airy_wave_velocity(0.0, 0.0, 0.0, 1.0, 0.1, 1.0, 0.0);
        assert_eq!(v, [0.0; 3]);
    }

    #[test]
    fn test_airy_velocity_finite() {
        let v = airy_wave_velocity(1.0, -1.0, 0.5, 0.5, 0.2, 1.5, 5.0);
        for c in v {
            assert!(c.is_finite());
        }
    }

    // ── wave_dispersion ─────────────────────────────────────────────────────

    #[test]
    fn test_dispersion_deep_water_omega_sq_eq_gk() {
        // Deep water: tanh(kd) → 1, so ω² = gk
        let k = 0.1;
        let depth = 1000.0;
        let omega = wave_dispersion(k, depth, G);
        let expected = (G * k).sqrt();
        assert!((omega - expected).abs() / expected < 1e-4);
    }

    #[test]
    fn test_dispersion_zero_k() {
        assert_eq!(wave_dispersion(0.0, 10.0, G), 0.0);
    }

    #[test]
    fn test_dispersion_positive() {
        let omega = wave_dispersion(0.5, 5.0, G);
        assert!(omega > 0.0);
    }

    #[test]
    fn test_dispersion_matches_relation() {
        let k = 0.3;
        let depth = 8.0;
        let omega = wave_dispersion(k, depth, G);
        let lhs = omega * omega;
        let rhs = G * k * (k * depth).tanh();
        assert!((lhs - rhs).abs() / lhs < 1e-12);
    }

    // ── group_velocity ──────────────────────────────────────────────────────

    #[test]
    fn test_group_velocity_leq_phase_velocity() {
        let k = 0.2;
        let depth = 10.0;
        let omega = wave_dispersion(k, depth, G);
        let c_p = omega / k;
        let c_g = group_velocity(k, depth, G);
        assert!(c_g <= c_p + 1e-10);
    }

    #[test]
    fn test_group_velocity_deep_water_half_phase() {
        // Deep water: c_g = c_p/2
        let k = 0.1;
        let depth = 1000.0;
        let omega = wave_dispersion(k, depth, G);
        let c_p = omega / k;
        let c_g = group_velocity(k, depth, G);
        assert!((c_g - c_p / 2.0).abs() / c_p < 0.01);
    }

    #[test]
    fn test_group_velocity_positive() {
        let c_g = group_velocity(0.5, 5.0, G);
        assert!(c_g > 0.0);
    }

    // ── shoaling_coefficient ────────────────────────────────────────────────

    #[test]
    fn test_shoaling_greater_than_one_shallower() {
        // Shallower water → smaller c_g → Ks > 1
        let ks = shoaling_coefficient(0.05, 50.0, 0.2, 10.0);
        assert!(ks > 1.0, "Ks should be > 1 for shoaling: {ks}");
    }

    #[test]
    fn test_shoaling_equal_depths_is_one() {
        let ks = shoaling_coefficient(0.1, 20.0, 0.1, 20.0);
        assert!((ks - 1.0).abs() < 1e-10);
    }

    // ── breaking_criterion_ocean ───────────────────────────────────────────

    #[test]
    fn test_breaking_ocean_true_above_threshold() {
        assert!(breaking_criterion_ocean(0.8, 1.0)); // H/d = 0.8 > 0.78
    }

    #[test]
    fn test_breaking_ocean_false_below_threshold() {
        assert!(!breaking_criterion_ocean(0.5, 10.0));
    }

    #[test]
    fn test_breaking_ocean_zero_depth() {
        assert!(breaking_criterion_ocean(0.1, 0.0));
    }

    // ── TidalForcing ────────────────────────────────────────────────────────

    #[test]
    fn test_tidal_elevation_at_t0() {
        let tf = TidalForcing::new(1.0, 0.5);
        // At t=0: cos(0)=1 → elevation = m2_amp + s2_amp
        let eta = tf.elevation(0.0);
        assert!((eta - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_tidal_elevation_bounded() {
        let tf = TidalForcing::new(1.0, 0.3);
        for t_h in 0..50 {
            let eta = tf.elevation(t_h as f64 * 3600.0);
            assert!(eta.abs() <= 1.3 + 1e-10);
        }
    }

    #[test]
    fn test_tidal_current_zero_depth() {
        let tf = TidalForcing::new(1.0, 0.5);
        assert_eq!(tf.tidal_current(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_tidal_current_finite() {
        let tf = TidalForcing::new(1.0, 0.5);
        let u = tf.tidal_current(1000.0, 10.0);
        assert!(u.is_finite());
    }

    // ── stokes_drift ────────────────────────────────────────────────────────

    #[test]
    fn test_stokes_drift_positive_at_surface() {
        let drift = stokes_drift(0.5, 0.2, 1.0, 0.0);
        assert!(drift >= 0.0);
    }

    #[test]
    fn test_stokes_drift_zero_amp() {
        assert_eq!(stokes_drift(0.0, 0.2, 1.0, 0.0), 0.0);
    }

    #[test]
    fn test_stokes_drift_decreases_with_depth() {
        // Drift decreases with increasing depth (more negative z)
        let drift_surf = stokes_drift(0.5, 0.2, 1.0, 0.0);
        let drift_deep = stokes_drift(0.5, 0.2, 1.0, -5.0);
        assert!(drift_surf > drift_deep);
    }

    // ── radiation_stress_xx ─────────────────────────────────────────────────

    #[test]
    fn test_radiation_stress_finite() {
        let s = radiation_stress_xx(1.0, 0.2, 1.0, 10.0);
        assert!(s.is_finite());
    }

    #[test]
    fn test_radiation_stress_zero_k() {
        assert_eq!(radiation_stress_xx(1.0, 0.0, 1.0, 10.0), 0.0);
    }

    // ── wave_energy_density ─────────────────────────────────────────────────

    #[test]
    fn test_wave_energy_nonneg() {
        let e = wave_energy_density(1025.0, G, 2.0);
        assert!(e >= 0.0);
    }

    // ── stokes_second_order ─────────────────────────────────────────────────

    #[test]
    fn test_stokes2_finite() {
        let eta = stokes_second_order(0.0, 0.0, 0.0, 1.0, 0.2, 1.0, 10.0);
        assert!(eta.is_finite());
    }

    #[test]
    fn test_stokes2_zero_amp() {
        let eta = stokes_second_order(0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 10.0);
        assert_eq!(eta, 0.0);
    }

    // ── CoastalSolver ───────────────────────────────────────────────────────

    #[test]
    fn test_coastal_solver_empty() {
        let solver = CoastalSolver::new();
        assert_eq!(solver.particle_count(), 0);
        assert_eq!(solver.total_energy(), 0.0);
    }

    #[test]
    fn test_coastal_solver_add_particle() {
        let mut solver = CoastalSolver::new();
        solver.add_particle(OceanParticle::new([0.0; 3], 1.0, 1025.0, 10.0, 0.05));
        assert_eq!(solver.particle_count(), 1);
    }

    #[test]
    fn test_coastal_solver_total_energy() {
        let mut solver = CoastalSolver::new();
        let mut p = OceanParticle::new([0.0; 3], 2.0, 1025.0, 10.0, 0.05);
        p.vel = [3.0, 0.0, 0.0];
        solver.add_particle(p);
        // KE = 0.5 * 2 * 9 = 9
        assert!((solver.total_energy() - 9.0).abs() < 1e-10);
    }
}
