//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::KB;

/// Radial distribution function g(r) for a set of positions in a periodic box.
#[derive(Debug, Clone)]
pub struct RadialDistributionFunction {
    /// Bin edges (left edge of each bin), length == n_bins.
    pub bins: Vec<f64>,
    /// g(r) values, one per bin.
    pub g_r: Vec<f64>,
}
impl RadialDistributionFunction {
    /// Number of bins.
    pub fn n_bins(&self) -> usize {
        self.bins.len()
    }
    /// Bin width (uniform).
    pub fn bin_width(&self) -> f64 {
        if self.bins.len() < 2 {
            return 0.0;
        }
        self.bins[1] - self.bins[0]
    }
    /// Returns the bin centre for bin index `i`.
    pub fn bin_center(&self, i: usize) -> f64 {
        self.bins[i] + self.bin_width() * 0.5
    }
    /// Returns the (r, g_r) pair for each bin as an iterator.
    pub fn iter(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        let bw = self.bin_width();
        self.bins
            .iter()
            .zip(self.g_r.iter())
            .map(move |(&b, &g)| (b + bw * 0.5, g))
    }
}
/// Mean squared displacement MSD(t) computed from a set of particle trajectories.
#[derive(Debug, Clone)]
pub struct MeanSquaredDisplacement {
    /// Time values (same units as `dt` passed to `compute_msd`).
    pub times: Vec<f64>,
    /// MSD values ⟨|r(t) − r(0)|²⟩ in units of length².
    pub msd: Vec<f64>,
}
impl MeanSquaredDisplacement {
    /// Number of time points.
    pub fn len(&self) -> usize {
        self.times.len()
    }
    /// True if there are no time points.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
    /// Estimate the diffusion coefficient from the long-time slope of the MSD.
    ///
    /// Uses the Einstein relation D = MSD(t) / (6 t) in 3D.
    /// Performs a simple linear regression over the second half of the data
    /// (avoids early-time ballistic regime).
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.times.len();
        if n < 4 {
            return 0.0;
        }
        let start = n / 2;
        let t_slice = &self.times[start..];
        let m_slice = &self.msd[start..];
        let n_pts = t_slice.len() as f64;
        let sum_t: f64 = t_slice.iter().sum();
        let sum_m: f64 = m_slice.iter().sum();
        let sum_t2: f64 = t_slice.iter().map(|&x| x * x).sum();
        let sum_tm: f64 = t_slice
            .iter()
            .zip(m_slice.iter())
            .map(|(&t, &m)| t * m)
            .sum();
        let denom = n_pts * sum_t2 - sum_t * sum_t;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        let slope = (n_pts * sum_tm - sum_t * sum_m) / denom;
        slope / 6.0
    }
}
/// Heat flux autocorrelation function J(t) = ⟨J(t)·J(0)⟩ / ⟨J(0)²⟩.
#[derive(Debug, Clone)]
pub struct HeatFluxAcf {
    /// Lag times.
    pub times: Vec<f64>,
    /// Normalised autocorrelation of the heat flux magnitude.
    pub acf: Vec<f64>,
    /// Unnormalised ⟨J(0)²⟩ value used for normalisation.
    pub j0_sq: f64,
}
impl HeatFluxAcf {
    /// Number of lag-time points.
    pub fn len(&self) -> usize {
        self.times.len()
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
    /// Estimate the thermal conductivity via the Green-Kubo relation:
    ///
    /// κ = V / (3 kB T²) ∫₀^∞ ⟨J(t)·J(0)⟩ dt
    ///
    /// This integrates the raw (un-normalised) ACF using the stored `j0_sq`.
    /// The `volume` must be in the same units as the heat flux.
    pub fn thermal_conductivity(&self, temperature: f64, volume: f64) -> f64 {
        if self.times.len() < 2 || temperature < 1e-10 || volume < 1e-30 {
            return 0.0;
        }
        let dt = self.times[1] - self.times[0];
        let integral: f64 = self
            .acf
            .windows(2)
            .map(|w| (w[0] + w[1]) * 0.5 * dt * self.j0_sq)
            .sum();
        volume / (3.0 * KB * temperature * temperature) * integral
    }
}
/// Angular momentum and rotational dynamics analysis.
pub struct AngularMomentumAnalysis;
impl AngularMomentumAnalysis {
    /// Compute the rotation correlation function C_rot(t) from angular velocity
    /// trajectories.
    ///
    /// ```text
    /// C_rot(lag) = <ω(t₀) · ω(t₀+lag)> / <ω(t₀)²>
    /// ```
    ///
    /// Averaged over all molecules and time origins.
    ///
    /// # Arguments
    /// * `angular_velocities` – `angular_velocities[molecule][frame]` = `[ωx, ωy, ωz]`.
    /// * `dt`                 – time step between frames.
    ///
    /// # Returns
    /// `(times, C_rot)` normalised so `C_rot(0) = 1`.
    pub fn compute_rotation_correlation(
        angular_velocities: &[Vec<[f64; 3]>],
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        if angular_velocities.is_empty() {
            return (vec![], vec![]);
        }
        let n_frames = angular_velocities[0].len();
        if n_frames == 0 {
            return (vec![], vec![]);
        }
        let mut c_raw = vec![0.0_f64; n_frames];
        let mut counts = vec![0u64; n_frames];
        for omega_traj in angular_velocities {
            if omega_traj.len() < n_frames {
                continue;
            }
            for lag in 0..n_frames {
                let n_origins = n_frames - lag;
                for t0 in 0..n_origins {
                    let o0 = omega_traj[t0];
                    let o1 = omega_traj[t0 + lag];
                    c_raw[lag] += o0[0] * o1[0] + o0[1] * o1[1] + o0[2] * o1[2];
                    counts[lag] += 1;
                }
            }
        }
        for (c, cnt) in c_raw.iter_mut().zip(counts.iter()) {
            if *cnt > 0 {
                *c /= *cnt as f64;
            }
        }
        let c0 = c_raw[0];
        let c_rot: Vec<f64> = if c0.abs() > 1e-30 {
            c_raw.iter().map(|&c| c / c0).collect()
        } else {
            vec![0.0; n_frames]
        };
        let times: Vec<f64> = (0..n_frames).map(|i| i as f64 * dt).collect();
        (times, c_rot)
    }
    /// Rotational diffusion coefficient from the long-time decay of C_rot(t).
    ///
    /// Fits an exponential exp(-6 D_rot t) to C_rot(t) using the ratio
    /// between adjacent points:
    ///
    /// ```text
    /// D_rot = -ln(C_rot(τ) / C_rot(0)) / (6 τ)
    /// ```
    ///
    /// Uses the midpoint τ = T/2 as the fit point.
    pub fn rotational_diffusion_coefficient(times: &[f64], c_rot: &[f64]) -> f64 {
        let n = times.len().min(c_rot.len());
        if n < 2 || c_rot[0].abs() < 1e-30 {
            return 0.0;
        }
        let mid = n / 2;
        let tau = times[mid];
        if tau < 1e-30 {
            return 0.0;
        }
        let ratio = c_rot[mid] / c_rot[0];
        if ratio <= 0.0 {
            return 0.0;
        }
        -ratio.ln() / (6.0 * tau)
    }
}
/// Velocity autocorrelation function VACF(t).
#[derive(Debug, Clone)]
pub struct VelocityAutocorrelation {
    /// Lag times (same units as `dt`).
    pub times: Vec<f64>,
    /// VACF values C(t) = ⟨v(t) · v(0)⟩ / ⟨v(0)²⟩ (normalised).
    pub vacf: Vec<f64>,
}
impl VelocityAutocorrelation {
    /// Number of lag-time points.
    pub fn len(&self) -> usize {
        self.times.len()
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
    /// Estimate the diffusion coefficient via the Green-Kubo relation.
    ///
    /// D = (1/3) ∫₀^∞ ⟨v(t) · v(0)⟩ dt
    ///
    /// Integrates using the trapezoidal rule over the stored VACF.
    /// The VACF must *not* be normalised (use the raw version) for this to
    /// return the correct D in physical units.  Here we return the normalised
    /// integral as a proxy (dimensionless).
    pub fn green_kubo_integral(&self) -> f64 {
        if self.times.len() < 2 {
            return 0.0;
        }
        let dt = self.times[1] - self.times[0];
        let sum: f64 = self.vacf.windows(2).map(|w| (w[0] + w[1]) * 0.5 * dt).sum();
        sum / 3.0
    }
}
impl VelocityAutocorrelation {
    /// Compute the vibrational density of states (VDOS) via the cosine transform
    /// of the VACF.
    ///
    /// The VDOS is defined as:
    ///
    /// ```text
    /// g(ω) = (2/π) ∫₀^∞ C(t) cos(ω t) dt
    /// ```
    ///
    /// This implementation uses the discrete cosine transform (type II) via
    /// direct summation for correctness and simplicity.
    ///
    /// # Arguments
    /// * `n_freq` – number of frequency points to evaluate.
    /// * `omega_max` – maximum angular frequency (rad/time unit of `dt`).
    ///
    /// # Returns
    /// `(frequencies, vdos)` each of length `n_freq`.
    pub fn compute_vibrational_density_of_states(
        &self,
        n_freq: usize,
        omega_max: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        if self.vacf.is_empty() || n_freq == 0 || omega_max <= 0.0 {
            return (vec![], vec![]);
        }
        let n_t = self.vacf.len();
        let dt = if self.times.len() >= 2 {
            self.times[1] - self.times[0]
        } else {
            1.0
        };
        let d_omega = omega_max / n_freq as f64;
        let mut frequencies = Vec::with_capacity(n_freq);
        let mut vdos = Vec::with_capacity(n_freq);
        for k in 0..n_freq {
            let omega = (k as f64 + 0.5) * d_omega;
            let mut integral = 0.0_f64;
            for t in 0..n_t {
                let time = self.times[t];
                let cos_term = (omega * time).cos();
                let weight = if t == 0 || t == n_t - 1 { 0.5 } else { 1.0 };
                integral += weight * self.vacf[t] * cos_term * dt;
            }
            frequencies.push(omega);
            vdos.push((2.0 / PI) * integral);
        }
        (frequencies, vdos)
    }
}
/// Structure analysis utilities.
pub struct StructureAnalysis;
impl StructureAnalysis {
    /// 3D pair correlation function (directional g(r)).
    ///
    /// Returns three histograms gx(r), gy(r), gz(r) – the one-dimensional
    /// pair correlations projected onto the x, y, and z axes respectively.
    /// Each histogram is normalised to g(r) → 1 at large r.
    ///
    /// # Arguments
    /// * `positions` – atomic positions (Å or nm, consistent with `box_size`).
    /// * `box_size`  – side length of a cubic periodic box.
    /// * `n_bins`    – number of histogram bins per axis.
    /// * `cutoff`    – maximum |Δx|, |Δy|, or |Δz| to include.
    ///
    /// # Returns
    /// `(bin_centers, gx, gy, gz)` each of length `n_bins`.
    pub fn compute_pair_correlation_3d(
        positions: &[[f64; 3]],
        box_size: f64,
        n_bins: usize,
        cutoff: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = positions.len();
        if n < 2 || n_bins == 0 || box_size <= 0.0 || cutoff <= 0.0 {
            let empty = vec![0.0; n_bins];
            let centers: Vec<f64> = (0..n_bins)
                .map(|i| (i as f64 + 0.5) * cutoff / n_bins as f64)
                .collect();
            return (centers, empty.clone(), empty.clone(), empty);
        }
        let bin_width = cutoff / n_bins as f64;
        let mut hx = vec![0u64; n_bins];
        let mut hy = vec![0u64; n_bins];
        let mut hz = vec![0u64; n_bins];
        for i in 0..n {
            for j in (i + 1)..n {
                let mut dx = positions[i][0] - positions[j][0];
                let mut dy = positions[i][1] - positions[j][1];
                let mut dz = positions[i][2] - positions[j][2];
                dx -= box_size * (dx / box_size).round();
                dy -= box_size * (dy / box_size).round();
                dz -= box_size * (dz / box_size).round();
                let ax = dx.abs();
                let ay = dy.abs();
                let az = dz.abs();
                if ax < cutoff {
                    hx[((ax / bin_width) as usize).min(n_bins - 1)] += 1;
                }
                if ay < cutoff {
                    hy[((ay / bin_width) as usize).min(n_bins - 1)] += 1;
                }
                if az < cutoff {
                    hz[((az / bin_width) as usize).min(n_bins - 1)] += 1;
                }
            }
        }
        let n_pairs = (n * (n - 1) / 2) as f64;
        let rho = n as f64 / (box_size * box_size * box_size);
        let norm_factor = n_pairs * rho * bin_width;
        let bin_centers: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * bin_width).collect();
        let normalize = |hist: &[u64]| -> Vec<f64> {
            hist.iter()
                .map(|&c| {
                    if norm_factor > 1e-30 {
                        c as f64 / norm_factor
                    } else {
                        0.0
                    }
                })
                .collect()
        };
        (bin_centers, normalize(&hx), normalize(&hy), normalize(&hz))
    }
}
