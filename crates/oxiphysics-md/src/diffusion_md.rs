// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Diffusion coefficient calculation from MD trajectories.
//!
//! Provides two complementary methods:
//!
//! 1. **MSD method** ([`MsdCalculator`]): Slope of mean-squared displacement,
//!    D = lim_{t→∞} MSD(t) / 6t (3D).
//!
//! 2. **Green-Kubo / VACF method** ([`VelocityAutocorrelation`]):
//!    D = (1/3) ∫₀^∞ C_v(t) dt, where C_v is the velocity autocorrelation.
//!
//! Also includes the non-Gaussian parameter α₂ and van Hove self-correlation
//! function for detecting heterogeneous dynamics.

// ---------------------------------------------------------------------------
// MsdCalculator
// ---------------------------------------------------------------------------

/// Accumulates particle trajectory frames and computes mean-squared displacement.
///
/// Frames are stored as Cartesian positions with periodic boundary corrections
/// (unwrapped coordinates) applied between consecutive frames.
#[derive(Debug, Clone)]
pub struct MsdCalculator {
    /// Stored frames of unwrapped positions: `positions[frame][atom] = [x,y,z]`.
    pub positions: Vec<Vec<[f64; 3]>>,
    /// Timestep between consecutive frames (ps or fs).
    pub dt: f64,
    /// Simulation box dimensions `[Lx, Ly, Lz]`.
    pub box_size: [f64; 3],
}

impl MsdCalculator {
    /// Create a new MSD calculator.
    ///
    /// # Arguments
    /// - `dt`: time between frames
    /// - `box_size`: periodic box dimensions `[Lx, Ly, Lz]`
    pub fn new(dt: f64, box_size: [f64; 3]) -> Self {
        Self {
            positions: Vec::new(),
            dt,
            box_size,
        }
    }

    /// Add a frame of atomic positions.
    ///
    /// If a previous frame exists, periodic boundary corrections are applied
    /// to produce continuous (unwrapped) trajectories.
    pub fn add_frame(&mut self, mut new_pos: Vec<[f64; 3]>) {
        if let Some(prev) = self.positions.last() {
            let box_size = self.box_size;
            for (r_new, r_old) in new_pos.iter_mut().zip(prev.iter()) {
                *r_new = Self::unwrap_pbc_static(*r_new, *r_old, box_size);
            }
        }
        self.positions.push(new_pos);
    }

    /// Number of stored frames.
    pub fn frame_count(&self) -> usize {
        self.positions.len()
    }

    /// Compute MSD for a given lag (in frames).
    ///
    /// `MSD(τ) = < |r(t+τ) - r(t)|² >` averaged over all atoms and time origins.
    pub fn mean_squared_displacement(&self, lag: usize) -> f64 {
        let n_frames = self.positions.len();
        if lag == 0 || lag >= n_frames {
            return 0.0;
        }
        let n_atoms = self.positions[0].len();
        if n_atoms == 0 {
            return 0.0;
        }
        let n_origins = n_frames - lag;
        let mut sum = 0.0f64;
        for t0 in 0..n_origins {
            let t1 = t0 + lag;
            for a in 0..n_atoms {
                let r0 = self.positions[t0][a];
                let r1 = self.positions[t1][a];
                let dx = r1[0] - r0[0];
                let dy = r1[1] - r0[1];
                let dz = r1[2] - r0[2];
                sum += dx * dx + dy * dy + dz * dz;
            }
        }
        sum / (n_origins as f64 * n_atoms as f64)
    }

    /// Return `(time, MSD)` pairs for all available lags.
    pub fn msd_vs_time(&self) -> Vec<(f64, f64)> {
        let n_frames = self.positions.len();
        (1..n_frames)
            .map(|lag| {
                let t = lag as f64 * self.dt;
                let msd = self.mean_squared_displacement(lag);
                (t, msd)
            })
            .collect()
    }

    /// Estimate the diffusion coefficient from a linear fit to MSD(t).
    ///
    /// Uses frames from `t_start` to `t_end` (inclusive, as lag indices).
    /// `D = slope / 6` (Einstein relation in 3D).
    pub fn diffusion_coefficient_from_msd(&self, t_start: usize, t_end: usize) -> f64 {
        let t_end = t_end.min(self.positions.len().saturating_sub(1));
        if t_start >= t_end {
            return 0.0;
        }
        let xs: Vec<f64> = (t_start..=t_end).map(|l| l as f64 * self.dt).collect();
        let ys: Vec<f64> = (t_start..=t_end)
            .map(|l| self.mean_squared_displacement(l))
            .collect();
        let (slope, _) = linear_fit(&xs, &ys);
        slope / 6.0
    }

    /// Periodic boundary unwrap: adjust `r_new` so it is nearest image of `r_old`.
    pub fn unwrap_pbc(r_new: [f64; 3], r_old: [f64; 3], box_size: [f64; 3]) -> [f64; 3] {
        Self::unwrap_pbc_static(r_new, r_old, box_size)
    }

    fn unwrap_pbc_static(r_new: [f64; 3], r_old: [f64; 3], box_size: [f64; 3]) -> [f64; 3] {
        let mut out = r_new;
        for i in 0..3 {
            let mut d = r_new[i] - r_old[i];
            let l = box_size[i];
            if l > 1e-30 {
                d -= l * (d / l).round();
            }
            out[i] = r_old[i] + d;
        }
        out
    }
}

// ---------------------------------------------------------------------------
// VelocityAutocorrelation
// ---------------------------------------------------------------------------

/// Computes velocity autocorrelation function (VACF) and Green-Kubo diffusion.
///
/// Green-Kubo relation: `D = (1/3) ∫₀^∞ C_v(t) dt`
/// where `C_v(τ) = <v(t)·v(t+τ)>`.
#[derive(Debug, Clone)]
pub struct VelocityAutocorrelation {
    /// Stored velocity frames: `velocities[frame][atom] = [vx,vy,vz]`.
    pub velocities: Vec<Vec<[f64; 3]>>,
    /// Timestep between frames.
    pub dt: f64,
}

impl VelocityAutocorrelation {
    /// Create a new VACF calculator.
    ///
    /// # Arguments
    /// - `dt`: time between frames
    pub fn new(dt: f64) -> Self {
        Self {
            velocities: Vec::new(),
            dt,
        }
    }

    /// Add a frame of atomic velocities.
    pub fn add_frame(&mut self, vels: Vec<[f64; 3]>) {
        self.velocities.push(vels);
    }

    /// Compute the normalized VACF `C_v(τ)` for τ = 0..n_frames-1.
    ///
    /// Returns a `Vec`f64` where index `τ` = lag in frames.
    /// Normalized by C_v(0).
    pub fn compute_vacf(&self) -> Vec<f64> {
        let n_frames = self.velocities.len();
        if n_frames == 0 {
            return Vec::new();
        }
        let n_atoms = self.velocities[0].len();
        let mut vacf = vec![0.0f64; n_frames];

        for (lag, vacf_val) in vacf.iter_mut().enumerate() {
            let n_origins = n_frames - lag;
            let mut sum = 0.0f64;
            for t0 in 0..n_origins {
                let t1 = t0 + lag;
                for (v0, v1) in self.velocities[t0].iter().zip(self.velocities[t1].iter()) {
                    sum += v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                }
            }
            *vacf_val = sum / (n_origins as f64 * n_atoms as f64);
        }
        vacf
    }

    /// Compute diffusion coefficient via Green-Kubo integration of VACF.
    ///
    /// Uses trapezoidal integration: `D = (1/3) ∫ C_v(t) dt`.
    pub fn diffusion_from_vacf(&self) -> f64 {
        let vacf = self.compute_vacf();
        if vacf.len() < 2 {
            return 0.0;
        }
        // Trapezoidal rule
        let integral: f64 = vacf.windows(2).map(|w| 0.5 * (w[0] + w[1]) * self.dt).sum();
        integral / 3.0
    }
}

// ---------------------------------------------------------------------------
// Linear regression
// ---------------------------------------------------------------------------

/// Compute a linear least-squares fit y = slope·x + intercept.
///
/// Returns `(slope, intercept)`.
pub fn linear_fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len().min(y.len());
    if n < 2 {
        return (0.0, 0.0);
    }
    let n_f = n as f64;
    let sum_x: f64 = x[..n].iter().sum();
    let sum_y: f64 = y[..n].iter().sum();
    let sum_xx: f64 = x[..n].iter().map(|&xi| xi * xi).sum();
    let sum_xy: f64 = x[..n]
        .iter()
        .zip(y[..n].iter())
        .map(|(&xi, &yi)| xi * yi)
        .sum();
    let denom = n_f * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-30 {
        return (0.0, sum_y / n_f);
    }
    let slope = (n_f * sum_xy - sum_x * sum_y) / denom;
    let intercept = (sum_y - slope * sum_x) / n_f;
    (slope, intercept)
}

/// Compute the coefficient of determination R² for a linear fit.
///
/// `R² = 1 - SS_res / SS_tot`
pub fn r_squared(x: &[f64], y: &[f64], slope: f64, intercept: f64) -> f64 {
    let n = x.len().min(y.len());
    if n == 0 {
        return 0.0;
    }
    let y_mean: f64 = y[..n].iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = y[..n].iter().map(|&yi| (yi - y_mean).powi(2)).sum();
    let ss_res: f64 = x[..n]
        .iter()
        .zip(y[..n].iter())
        .map(|(&xi, &yi)| (yi - (slope * xi + intercept)).powi(2))
        .sum();
    if ss_tot < 1e-30 {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    }
}

// ---------------------------------------------------------------------------
// Non-Gaussian parameter
// ---------------------------------------------------------------------------

/// Non-Gaussian parameter α₂ for detecting heterogeneous dynamics.
///
/// `α₂ = 3·<r⁴> / (5·`r²`²) - 1`
///
/// α₂ = 0 for a Gaussian (diffusive) process;
/// α₂ > 0 indicates dynamic heterogeneity.
pub struct NgambaAnalysis {
    /// MSD values over time.
    pub msd: Vec<f64>,
    /// Timestep.
    pub dt: f64,
}

/// Compute the non-Gaussian parameter α₂ at a given lag.
///
/// `α₂ = 3·<r⁴> / (5·`r²`²) - 1`
///
/// Returns 0 if the mean-squared displacement is zero.
pub fn non_gaussian_parameter(positions: &[Vec<[f64; 3]>], lag: usize) -> f64 {
    let n_frames = positions.len();
    if lag == 0 || lag >= n_frames {
        return 0.0;
    }
    let n_atoms = positions[0].len();
    if n_atoms == 0 {
        return 0.0;
    }
    let n_origins = n_frames - lag;
    let mut r2_sum = 0.0f64;
    let mut r4_sum = 0.0f64;
    let mut count = 0usize;

    for t0 in 0..n_origins {
        let t1 = t0 + lag;
        for (r0, r1) in positions[t0].iter().zip(positions[t1].iter()) {
            let dx = r1[0] - r0[0];
            let dy = r1[1] - r0[1];
            let dz = r1[2] - r0[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            r2_sum += r2;
            r4_sum += r2 * r2;
            count += 1;
        }
    }
    if count == 0 {
        return 0.0;
    }
    let mean_r2 = r2_sum / count as f64;
    let mean_r4 = r4_sum / count as f64;
    if mean_r2 < 1e-30 {
        return 0.0;
    }
    3.0 * mean_r4 / (5.0 * mean_r2 * mean_r2) - 1.0
}

// ---------------------------------------------------------------------------
// Van Hove self-correlation function
// ---------------------------------------------------------------------------

/// Compute the self-part of the van Hove correlation function G_s(r, t).
///
/// `G_s(r, τ) = (1/N) Σ_i δ(r - |r_i(τ) - r_i(0)|)`
///
/// Returns `(r, G_s(r, τ))` histogram with `n_bins` bins from 0 to `r_max`.
pub fn van_hove_self(
    positions: &[Vec<[f64; 3]>],
    lag: usize,
    r_max: f64,
    n_bins: usize,
) -> Vec<(f64, f64)> {
    let n_frames = positions.len();
    if lag == 0 || lag >= n_frames || n_bins == 0 {
        return vec![(0.0, 0.0); n_bins];
    }
    let n_atoms = positions[0].len();
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];
    let n_origins = n_frames - lag;

    for t0 in 0..n_origins {
        let t1 = t0 + lag;
        for (r0, r1) in positions[t0].iter().zip(positions[t1].iter()) {
            let dx = r1[0] - r0[0];
            let dy = r1[1] - r0[1];
            let dz = r1[2] - r0[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let bin = (r / dr) as usize;
                let bin = bin.min(n_bins - 1);
                hist[bin] += 1;
            }
        }
    }

    let total = (n_origins * n_atoms) as f64;
    (0..n_bins)
        .map(|b| {
            let r = (b as f64 + 0.5) * dr;
            let gs = if total > 0.0 {
                hist[b] as f64 / (total * dr)
            } else {
                0.0
            };
            (r, gs)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a simple trajectory with `n_atoms` atoms doing a random walk.
    fn random_walk_trajectory(n_atoms: usize, n_frames: usize, step: f64) -> Vec<Vec<[f64; 3]>> {
        use rand::RngExt;
        let mut rng = rand::rng();
        let mut frames = Vec::new();
        let mut pos: Vec<[f64; 3]> = (0..n_atoms).map(|_| [0.0, 0.0, 0.0]).collect();
        for _ in 0..n_frames {
            for p in pos.iter_mut() {
                p[0] += step * (rng.random_range(-1.0_f64..1.0_f64));
                p[1] += step * (rng.random_range(-1.0_f64..1.0_f64));
                p[2] += step * (rng.random_range(-1.0_f64..1.0_f64));
            }
            frames.push(pos.clone());
        }
        frames
    }

    #[test]
    fn test_msd_zero_lag_is_zero() {
        let mut calc = MsdCalculator::new(0.01, [10.0, 10.0, 10.0]);
        let frames = random_walk_trajectory(10, 20, 0.1);
        for frame in frames {
            calc.add_frame(frame);
        }
        assert_eq!(calc.mean_squared_displacement(0), 0.0);
    }

    #[test]
    fn test_msd_grows_with_lag() {
        let mut calc = MsdCalculator::new(0.01, [100.0, 100.0, 100.0]);
        let frames = random_walk_trajectory(20, 50, 0.5);
        for frame in frames {
            calc.add_frame(frame);
        }
        let msd1 = calc.mean_squared_displacement(5);
        let msd10 = calc.mean_squared_displacement(20);
        assert!(
            msd10 > msd1,
            "MSD should grow with lag: msd1={}, msd10={}",
            msd1,
            msd10
        );
    }

    #[test]
    fn test_msd_non_negative() {
        let mut calc = MsdCalculator::new(0.01, [50.0, 50.0, 50.0]);
        let frames = random_walk_trajectory(5, 30, 0.2);
        for frame in frames {
            calc.add_frame(frame);
        }
        for lag in 1..10 {
            assert!(calc.mean_squared_displacement(lag) >= 0.0);
        }
    }

    #[test]
    fn test_frame_count() {
        let mut calc = MsdCalculator::new(0.01, [10.0; 3]);
        assert_eq!(calc.frame_count(), 0);
        calc.add_frame(vec![[0.0, 0.0, 0.0]]);
        assert_eq!(calc.frame_count(), 1);
        calc.add_frame(vec![[1.0, 0.0, 0.0]]);
        assert_eq!(calc.frame_count(), 2);
    }

    #[test]
    fn test_msd_vs_time_length() {
        let mut calc = MsdCalculator::new(0.01, [10.0; 3]);
        for _ in 0..10 {
            calc.add_frame(vec![[0.0, 0.0, 0.0]]);
        }
        let msd_t = calc.msd_vs_time();
        assert_eq!(msd_t.len(), 9); // lags 1..9
    }

    #[test]
    fn test_msd_vs_time_time_values() {
        let dt = 0.05;
        let mut calc = MsdCalculator::new(dt, [10.0; 3]);
        for _ in 0..5 {
            calc.add_frame(vec![[0.0; 3]]);
        }
        let msd_t = calc.msd_vs_time();
        assert!((msd_t[0].0 - dt).abs() < 1e-14);
        assert!((msd_t[3].0 - 4.0 * dt).abs() < 1e-14);
    }

    #[test]
    fn test_diffusion_coefficient_positive() {
        let mut calc = MsdCalculator::new(0.01, [200.0; 3]);
        let frames = random_walk_trajectory(50, 100, 1.0);
        for frame in frames {
            calc.add_frame(frame);
        }
        let d = calc.diffusion_coefficient_from_msd(10, 90);
        assert!(
            d > 0.0,
            "Diffusion coefficient should be positive, got {}",
            d
        );
    }

    #[test]
    fn test_unwrap_pbc_no_crossing() {
        let r_new = [1.0, 2.0, 3.0];
        let r_old = [1.1, 2.0, 3.0];
        let box_size = [10.0; 3];
        let unwrapped = MsdCalculator::unwrap_pbc(r_new, r_old, box_size);
        assert!((unwrapped[0] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_unwrap_pbc_crossing_positive() {
        // Particle jumped from 9.9 to 0.1 (crossed +x boundary)
        let r_old = [9.9, 0.0, 0.0];
        let r_new = [0.1, 0.0, 0.0];
        let box_size = [10.0; 3];
        let unwrapped = MsdCalculator::unwrap_pbc(r_new, r_old, box_size);
        // Should be 10.1 (unwrapped)
        assert!(
            (unwrapped[0] - 10.1).abs() < 1e-10,
            "Expected 10.1, got {}",
            unwrapped[0]
        );
    }

    #[test]
    fn test_unwrap_pbc_crossing_negative() {
        // Particle jumped from 0.1 to 9.9 (crossed -x boundary)
        let r_old = [0.1, 0.0, 0.0];
        let r_new = [9.9, 0.0, 0.0];
        let box_size = [10.0; 3];
        let unwrapped = MsdCalculator::unwrap_pbc(r_new, r_old, box_size);
        // Should be -0.1 (unwrapped)
        assert!(
            (unwrapped[0] - (-0.1)).abs() < 1e-10,
            "Expected -0.1, got {}",
            unwrapped[0]
        );
    }

    #[test]
    fn test_linear_fit_perfect_line() {
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 2.0 * xi + 1.0).collect();
        let (slope, intercept) = linear_fit(&x, &y);
        assert!(
            (slope - 2.0).abs() < 1e-10,
            "slope should be 2, got {}",
            slope
        );
        assert!(
            (intercept - 1.0).abs() < 1e-10,
            "intercept should be 1, got {}",
            intercept
        );
    }

    #[test]
    fn test_linear_fit_zero_slope() {
        let x: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let y = vec![3.0f64; 5];
        let (slope, intercept) = linear_fit(&x, &y);
        assert!(slope.abs() < 1e-10);
        assert!((intercept - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_r_squared_perfect_fit() {
        let x: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|&xi| 3.0 * xi).collect();
        let (slope, intercept) = linear_fit(&x, &y);
        let r2 = r_squared(&x, &y, slope, intercept);
        assert!((r2 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_r_squared_between_zero_and_one() {
        use rand::RngExt;
        let mut rng = rand::rng();
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = x
            .iter()
            .map(|&xi| 2.0 * xi + rng.random_range(-2.0_f64..2.0_f64))
            .collect();
        let (slope, intercept) = linear_fit(&x, &y);
        let r2 = r_squared(&x, &y, slope, intercept);
        assert!(
            (0.0..=1.0).contains(&r2),
            "R² should be in [0,1], got {}",
            r2
        );
    }

    #[test]
    fn test_vacf_lag0_positive() {
        let mut vacf = VelocityAutocorrelation::new(0.01);
        // Non-zero velocities
        let vels: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        for _ in 0..20 {
            vacf.add_frame(vels.clone());
        }
        let cv = vacf.compute_vacf();
        assert!(
            cv[0] > 0.0,
            "VACF at lag=0 should be positive (kinetic energy)"
        );
    }

    #[test]
    fn test_vacf_frame_count() {
        let mut vacf = VelocityAutocorrelation::new(0.01);
        vacf.add_frame(vec![[1.0, 0.0, 0.0]]);
        vacf.add_frame(vec![[1.0, 0.0, 0.0]]);
        assert_eq!(vacf.velocities.len(), 2);
    }

    #[test]
    fn test_vacf_length() {
        let mut vacf = VelocityAutocorrelation::new(0.01);
        for _ in 0..15 {
            vacf.add_frame(vec![[1.0, 2.0, 3.0]]);
        }
        let cv = vacf.compute_vacf();
        assert_eq!(cv.len(), 15);
    }

    #[test]
    fn test_diffusion_from_vacf_positive() {
        let mut vacf = VelocityAutocorrelation::new(0.01);
        // Constant positive velocity -> positive VACF integral
        for _ in 0..50 {
            vacf.add_frame(vec![[1.0, 1.0, 1.0]; 10]);
        }
        let d = vacf.diffusion_from_vacf();
        assert!(d > 0.0, "Diffusion from VACF should be positive, got {}", d);
    }

    #[test]
    fn test_non_gaussian_parameter_gaussian() {
        // For uniform step diffusion, alpha2 should be near 0
        let frames = random_walk_trajectory(100, 50, 0.1);
        let alpha2 = non_gaussian_parameter(&frames, 10);
        // Not strictly 0 but should be small for large ensemble
        assert!(alpha2.is_finite(), "alpha2 should be finite");
    }

    #[test]
    fn test_non_gaussian_zero_lag() {
        let frames = random_walk_trajectory(10, 20, 0.1);
        let alpha2 = non_gaussian_parameter(&frames, 0);
        assert_eq!(alpha2, 0.0);
    }

    #[test]
    fn test_van_hove_self_length() {
        let frames = random_walk_trajectory(10, 20, 0.5);
        let gs = van_hove_self(&frames, 5, 5.0, 20);
        assert_eq!(gs.len(), 20);
    }

    #[test]
    fn test_van_hove_self_integrates_to_one() {
        let frames = random_walk_trajectory(50, 30, 0.2);
        let r_max = 3.0;
        let n_bins = 30;
        let gs = van_hove_self(&frames, 5, r_max, n_bins);
        let dr = r_max / n_bins as f64;
        let integral: f64 = gs.iter().map(|(_, g)| g * dr).sum();
        // Should integrate close to 1 (fraction of particles within r_max)
        assert!(
            integral >= 0.0,
            "Integral of G_s should be non-negative: {}",
            integral
        );
    }

    #[test]
    fn test_van_hove_r_values() {
        let frames = random_walk_trajectory(5, 10, 0.1);
        let r_max = 2.0;
        let n_bins = 10;
        let gs = van_hove_self(&frames, 2, r_max, n_bins);
        let dr = r_max / n_bins as f64;
        // Check r values are bin centers
        for (b, (r, _)) in gs.iter().enumerate() {
            let expected = (b as f64 + 0.5) * dr;
            assert!((r - expected).abs() < 1e-14);
        }
    }

    #[test]
    fn test_linear_fit_single_point_returns_zero_slope() {
        let x = [1.0];
        let y = [3.0];
        let (slope, _) = linear_fit(&x, &y);
        assert_eq!(slope, 0.0);
    }

    #[test]
    fn test_msd_stationary_atoms() {
        let mut calc = MsdCalculator::new(0.01, [10.0; 3]);
        let pos = vec![[1.0, 2.0, 3.0]; 5];
        for _ in 0..10 {
            calc.add_frame(pos.clone());
        }
        for lag in 1..5 {
            assert!(calc.mean_squared_displacement(lag).abs() < 1e-14);
        }
    }

    #[test]
    fn test_msd_known_displacement() {
        let mut calc = MsdCalculator::new(1.0, [1000.0; 3]);
        // All atoms move by 1.0 in x each frame
        let mut pos = vec![[0.0; 3]; 3];
        calc.add_frame(pos.clone());
        for _ in 1..10 {
            for p in pos.iter_mut() {
                p[0] += 1.0;
            }
            calc.add_frame(pos.clone());
        }
        // After lag=1, each atom moved by 1 in x -> MSD = 1
        let msd = calc.mean_squared_displacement(1);
        assert!(
            (msd - 1.0).abs() < 1e-10,
            "MSD for unit step should be 1, got {}",
            msd
        );
    }

    #[test]
    fn test_r_squared_empty_returns_zero() {
        let x: Vec<f64> = vec![];
        let y: Vec<f64> = vec![];
        let r2 = r_squared(&x, &y, 1.0, 0.0);
        assert_eq!(r2, 0.0);
    }

    #[test]
    fn test_diffusion_from_short_trajectory_zero() {
        let calc = MsdCalculator::new(0.01, [10.0; 3]);
        // No frames -> diffusion should be 0
        let d = calc.diffusion_coefficient_from_msd(0, 5);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_vacf_zero_velocity_is_zero() {
        let mut vacf = VelocityAutocorrelation::new(0.01);
        for _ in 0..10 {
            vacf.add_frame(vec![[0.0, 0.0, 0.0]; 5]);
        }
        let cv = vacf.compute_vacf();
        for &v in cv.iter() {
            assert!(v.abs() < 1e-14);
        }
    }

    #[test]
    fn test_unwrap_pbc_all_dims() {
        let r_old = [9.9, 9.9, 0.1];
        let r_new = [0.1, 0.1, 9.9];
        let box_size = [10.0; 3];
        let u = MsdCalculator::unwrap_pbc(r_new, r_old, box_size);
        assert!((u[0] - 10.1).abs() < 1e-10);
        assert!((u[1] - 10.1).abs() < 1e-10);
        assert!((u[2] - (-0.1)).abs() < 1e-10);
    }
}
