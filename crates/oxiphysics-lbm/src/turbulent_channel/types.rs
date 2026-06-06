//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    B_LL, BB, CX, CY, KAPPA_VK, W9, Y_PLUS_LOG_UPPER, Y_PLUS_TRANSITION, equilibrium_d2q9,
    moments_d2q9,
};

/// D2Q9 LBM solver for turbulent channel flow with wall model.
///
/// The flow is driven by a constant body force (`forcing`) in the x-direction.
/// No-slip walls are applied at j=0 and j=ny-1 via bounce-back.
#[derive(Debug, Clone)]
pub struct ChannelFlow {
    /// Number of lattice nodes in the streamwise direction.
    pub nx: usize,
    /// Number of lattice nodes in the wall-normal direction (including walls).
    pub ny: usize,
    /// Friction Reynolds number Re_τ = u_τ·h/ν.
    pub re_tau: f64,
    /// Friction velocity u_τ.
    pub u_tau: f64,
    /// Kinematic viscosity ν.
    pub nu: f64,
    /// Distribution functions (size nx*ny*9).
    pub f_dist: Vec<f64>,
    /// Velocity field \[(ux, uy)\] at each node.
    pub u: Vec<[f64; 2]>,
    /// Density field at each node.
    pub rho: Vec<f64>,
    /// Constant body force (pressure gradient) in x-direction.
    pub forcing: f64,
}
impl ChannelFlow {
    /// Create a new channel flow solver.
    ///
    /// # Arguments
    /// - `nx`: streamwise lattice size
    /// - `ny`: wall-normal lattice size (walls at j=0 and j=ny-1)
    /// - `re_tau`: target friction Reynolds number
    pub fn new(nx: usize, ny: usize, re_tau: f64) -> Self {
        let u_tau = 0.01_f64;
        let h = ny as f64 / 2.0;
        let nu = u_tau * h / re_tau;
        let forcing = 2.0 * u_tau * u_tau / (ny as f64);
        let n = nx * ny;
        let rho = vec![1.0f64; n];
        let u_field = vec![[0.0f64; 2]; n];
        let mut f_dist = vec![0.0f64; n * 9];
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                for q in 0..9 {
                    f_dist[idx * 9 + q] = W9[q] * rho[idx];
                }
            }
        }
        Self {
            nx,
            ny,
            re_tau,
            u_tau,
            nu,
            f_dist,
            u: u_field,
            rho,
            forcing,
        }
    }
    /// Flat node index for (i, j).
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }
    /// Flat distribution function index for node (i, j), direction q.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (i * self.ny + j) * 9 + q
    }
    /// Initialize velocity field with a parabolic Poiseuille profile.
    pub fn init_parabolic(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let u_max = 1.5 * self.bulk_velocity_from_forcing();
        for i in 0..nx {
            for j in 0..ny {
                let yc = 2.0 * j as f64 / (ny as f64 - 1.0) - 1.0;
                let ux = u_max * (1.0 - yc * yc);
                let idx = self.index(i, j);
                self.u[idx] = [ux, 0.0];
                let rho = self.rho[idx];
                for q in 0..9 {
                    let fi = self.fi(i, j, q);
                    self.f_dist[fi] = equilibrium_d2q9(rho, ux, 0.0, q);
                }
            }
        }
    }
    /// Estimate bulk velocity from forcing and viscosity (Poiseuille formula).
    fn bulk_velocity_from_forcing(&self) -> f64 {
        let h = (self.ny as f64) / 2.0;
        self.forcing * h * h / (3.0 * self.nu)
    }
    /// Compute wall shear stress τ_w = ρ u_τ².
    pub fn wall_shear_stress(&self) -> f64 {
        let u_tau = self.friction_velocity();
        let rho_mean: f64 = self.rho.iter().sum::<f64>() / self.rho.len() as f64;
        rho_mean * u_tau * u_tau
    }
    /// Viscous sublayer thickness δ_v = ν / u_τ.
    pub fn viscous_sublayer_thickness(&self) -> f64 {
        self.nu / self.u_tau.max(1e-30)
    }
    /// Log-law velocity u⁺ = (1/κ)·ln(y⁺) + B.
    pub fn log_law_u_plus(&self, y_plus: f64) -> f64 {
        if y_plus <= 0.0 {
            return 0.0;
        }
        if y_plus < Y_PLUS_TRANSITION {
            y_plus
        } else {
            (1.0 / KAPPA_VK) * y_plus.ln() + B_LL
        }
    }
    /// Volume-averaged (bulk) velocity in the x-direction.
    pub fn bulk_velocity(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for i in 0..nx {
            for j in 1..ny - 1 {
                sum += self.u[self.index(i, j)][0];
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
    /// Velocity at the channel centerline (j = ny/2).
    pub fn centerline_velocity(&self) -> f64 {
        let j_center = self.ny / 2;
        let mut sum = 0.0f64;
        for i in 0..self.nx {
            sum += self.u[self.index(i, j_center)][0];
        }
        sum / self.nx as f64
    }
    /// Return the wall-normal velocity profile as (y, u_x) pairs.
    pub fn velocity_profile(&self) -> Vec<(f64, f64)> {
        let nx = self.nx;
        let ny = self.ny;
        (0..ny)
            .map(|j| {
                let mean_u: f64 =
                    (0..nx).map(|i| self.u[self.index(i, j)][0]).sum::<f64>() / nx as f64;
                (j as f64, mean_u)
            })
            .collect()
    }
    /// Apply the constant body force to the x-momentum.
    pub fn apply_forcing(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 1..ny - 1 {
                let idx = self.index(i, j);
                self.u[idx][0] += self.forcing;
            }
        }
    }
    /// Apply no-slip (bounce-back) boundary conditions at j=0 and j=ny-1.
    pub fn no_slip_walls(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for &j in &[0, ny - 1] {
                for (q, &bb_q) in BB.iter().enumerate() {
                    let fi_src = self.fi(i, j, q);
                    let fi_bb = self.fi(i, j, bb_q);
                    let val = self.f_dist[fi_src];
                    self.f_dist[fi_bb] = val;
                }
            }
        }
    }
    /// D2Q9 equilibrium distribution.
    pub fn equilibrium(&self, rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
        equilibrium_d2q9(rho, ux, uy, q)
    }
    /// BGK collision, streaming, and wall forcing.
    pub fn collide_stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nu = self.nu;
        let tau = 3.0 * nu + 0.5;
        let omega = 1.0 / tau;
        let forcing = self.forcing;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let base = idx * 9;
                let (rho, mut ux, uy) = moments_d2q9(&self.f_dist, base);
                if j > 0 && j < ny - 1 {
                    ux += forcing;
                }
                self.rho[idx] = rho;
                self.u[idx] = [ux, uy];
                for q in 0..9 {
                    let fi = base + q;
                    let feq = equilibrium_d2q9(rho, ux, uy, q);
                    self.f_dist[fi] += omega * (feq - self.f_dist[fi]);
                }
            }
        }
        let mut f_new = vec![0.0f64; nx * ny * 9];
        for i in 0..nx {
            for j in 0..ny {
                for q in 0..9 {
                    let nj = j as i32 + CY[q];
                    if nj < 0 || nj >= ny as i32 {
                        let fi_src = self.fi(i, j, q);
                        let fi_dst = self.fi(i, j, BB[q]);
                        f_new[fi_dst] = self.f_dist[fi_src];
                    } else {
                        let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
                        let src = self.fi(i, j, q);
                        let dst = self.fi(ni, nj as usize, q);
                        f_new[dst] = self.f_dist[src];
                    }
                }
            }
        }
        self.f_dist = f_new;
    }
    /// Perform one complete channel flow timestep.
    pub fn step(&mut self) {
        self.collide_stream();
    }
    /// Estimate the Reynolds stress -<u'v'> at wall-normal row j.
    pub fn reynolds_stress(&self, j: usize) -> f64 {
        let nx = self.nx;
        let mean_u: f64 = (0..nx).map(|i| self.u[self.index(i, j)][0]).sum::<f64>() / nx as f64;
        let mean_v: f64 = (0..nx).map(|i| self.u[self.index(i, j)][1]).sum::<f64>() / nx as f64;
        let cov: f64 = (0..nx)
            .map(|i| {
                let idx = self.index(i, j);
                (self.u[idx][0] - mean_u) * (self.u[idx][1] - mean_v)
            })
            .sum::<f64>()
            / nx as f64;
        -cov
    }
    /// Estimate the friction velocity u_τ from the mean body force balance.
    pub fn friction_velocity(&self) -> f64 {
        let h = self.ny as f64 / 2.0;
        let tau_w = self.forcing * h;
        let rho_mean: f64 = self.rho.iter().sum::<f64>() / self.rho.len() as f64;
        (tau_w / rho_mean).sqrt()
    }
}
/// Turbulent kinetic energy budget terms at a wall-normal location.
///
/// Follows the standard decomposition: Production – Dissipation + Transport = 0.
#[derive(Debug, Clone, Copy, Default)]
pub struct TkeProductionBudget {
    /// TKE production P = -<u'v'> ∂`u`/∂y.
    pub production: f64,
    /// Viscous dissipation ε (modelled).
    pub dissipation: f64,
    /// Turbulent transport T (pressure + velocity diffusion).
    pub transport: f64,
    /// Viscous diffusion D.
    pub viscous_diffusion: f64,
}
impl TkeProductionBudget {
    /// Estimate budget from LBM fields at row j.
    ///
    /// Uses finite differences for gradients.
    ///
    /// # Arguments
    /// * `rs_j`   — Reynolds stress -<u'v'> at j
    /// * `du_dy`  — mean velocity gradient ∂`u`/∂y at j
    /// * `nu`     — kinematic viscosity
    /// * `k_j`    — TKE at j
    /// * `lm`     — mixing length at j (for dissipation model)
    pub fn estimate(rs_j: f64, du_dy: f64, nu: f64, k_j: f64, lm: f64) -> Self {
        let production = rs_j * du_dy;
        let c_mu = 0.09_f64;
        let dissipation = if lm > 1e-30 {
            c_mu * k_j.max(0.0).powf(1.5) / lm
        } else {
            0.0
        };
        let delta_nu = nu.sqrt().max(1e-30);
        let viscous_diffusion = nu * k_j / (delta_nu * delta_nu);
        let transport = dissipation - production - viscous_diffusion;
        Self {
            production,
            dissipation,
            transport,
            viscous_diffusion,
        }
    }
    /// Net residual (should be ≈ 0 in balanced budget).
    pub fn residual(&self) -> f64 {
        self.production - self.dissipation + self.transport + self.viscous_diffusion
    }
}
/// Collection of turbulent channel flow statistics.
///
/// Stores mean, rms, and Reynolds stress profiles obtained from time-averaging
/// or ensemble-averaging an LBM simulation.
#[derive(Debug, Clone)]
pub struct TurbulentStatistics {
    /// Number of wall-normal points.
    pub ny: usize,
    /// Accumulated mean ux at each j.
    pub sum_u: Vec<f64>,
    /// Accumulated mean ux² at each j (for rms).
    pub sum_u2: Vec<f64>,
    /// Accumulated mean uy² at each j (for v-rms).
    pub sum_v2: Vec<f64>,
    /// Accumulated -<u'v'> at each j (Reynolds stress).
    pub sum_uv: Vec<f64>,
    /// Number of samples accumulated.
    pub count: u64,
}
impl TurbulentStatistics {
    /// Create a new empty statistics collector for `ny` wall-normal levels.
    pub fn new(ny: usize) -> Self {
        Self {
            ny,
            sum_u: vec![0.0; ny],
            sum_u2: vec![0.0; ny],
            sum_v2: vec![0.0; ny],
            sum_uv: vec![0.0; ny],
            count: 0,
        }
    }
    /// Add one velocity snapshot to the running statistics.
    ///
    /// # Arguments
    /// * `u_field` — slice of `[ux, uy]` at each (i, j) node (row-major over i)
    /// * `nx`      — streamwise lattice size
    pub fn accumulate(&mut self, u_field: &[[f64; 2]], nx: usize) {
        let ny = self.ny;
        for j in 0..ny {
            let mean_u: f64 = (0..nx).map(|i| u_field[i * ny + j][0]).sum::<f64>() / nx as f64;
            let mean_v: f64 = (0..nx).map(|i| u_field[i * ny + j][1]).sum::<f64>() / nx as f64;
            let mean_u2: f64 =
                (0..nx).map(|i| u_field[i * ny + j][0].powi(2)).sum::<f64>() / nx as f64;
            let mean_v2: f64 =
                (0..nx).map(|i| u_field[i * ny + j][1].powi(2)).sum::<f64>() / nx as f64;
            let mean_uv: f64 = (0..nx)
                .map(|i| {
                    let ux = u_field[i * ny + j][0];
                    let uy = u_field[i * ny + j][1];
                    (ux - mean_u) * (uy - mean_v)
                })
                .sum::<f64>()
                / nx as f64;
            self.sum_u[j] += mean_u;
            self.sum_u2[j] += mean_u2;
            self.sum_v2[j] += mean_v2;
            self.sum_uv[j] += mean_uv;
        }
        self.count += 1;
    }
    /// Mean streamwise velocity profile `u`(j).
    pub fn mean_u_profile(&self) -> Vec<f64> {
        let c = self.count.max(1) as f64;
        self.sum_u.iter().map(|&s| s / c).collect()
    }
    /// Streamwise rms velocity u'_rms(j) = sqrt(`u²` - `u`²).
    pub fn u_rms_profile(&self) -> Vec<f64> {
        let c = self.count.max(1) as f64;
        (0..self.ny)
            .map(|j| {
                let mean = self.sum_u[j] / c;
                let mean2 = self.sum_u2[j] / c;
                (mean2 - mean * mean).max(0.0).sqrt()
            })
            .collect()
    }
    /// Wall-normal rms velocity v'_rms(j) = sqrt(`v²`).
    pub fn v_rms_profile(&self) -> Vec<f64> {
        let c = self.count.max(1) as f64;
        self.sum_v2
            .iter()
            .map(|&s| (s / c).max(0.0).sqrt())
            .collect()
    }
    /// Reynolds stress profile -<u'v'>(j).
    pub fn reynolds_stress_profile(&self) -> Vec<f64> {
        let c = self.count.max(1) as f64;
        self.sum_uv.iter().map(|&s| -(s / c)).collect()
    }
    /// Turbulent kinetic energy k = 0.5 * (<u'²> + <v'²>) at each j.
    pub fn tke_profile(&self) -> Vec<f64> {
        let u_rms = self.u_rms_profile();
        let v_rms = self.v_rms_profile();
        u_rms
            .iter()
            .zip(v_rms.iter())
            .map(|(&ur, &vr)| 0.5 * (ur * ur + vr * vr))
            .collect()
    }
    /// Maximum turbulent kinetic energy.
    pub fn max_tke(&self) -> f64 {
        self.tke_profile().iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Total Reynolds stress gradient d(-<u'v'>)/dy (finite difference).
    pub fn reynolds_stress_gradient(&self) -> Vec<f64> {
        let rs = self.reynolds_stress_profile();
        let ny = self.ny;
        let mut grad = vec![0.0; ny];
        for j in 1..ny - 1 {
            grad[j] = (rs[j + 1] - rs[j - 1]) * 0.5;
        }
        if ny >= 2 {
            grad[0] = rs[1] - rs[0];
            grad[ny - 1] = rs[ny - 1] - rs[ny - 2];
        }
        grad
    }
    /// Reset all accumulated statistics.
    pub fn reset(&mut self) {
        for v in self.sum_u.iter_mut() {
            *v = 0.0;
        }
        for v in self.sum_u2.iter_mut() {
            *v = 0.0;
        }
        for v in self.sum_v2.iter_mut() {
            *v = 0.0;
        }
        for v in self.sum_uv.iter_mut() {
            *v = 0.0;
        }
        self.count = 0;
    }
}
/// Dimensionless parameters for a turbulent channel flow.
///
/// These govern both the physical setup and the LBM grid resolution.
#[derive(Debug, Clone)]
pub struct TurbulentChannelParams {
    /// Friction Reynolds number Re_τ = u_τ h / ν.
    pub re_tau: f64,
    /// Friction velocity u_τ (lattice units).
    pub u_tau: f64,
    /// Channel half-width h (lattice units).
    pub half_width: f64,
    /// Bulk (mean) velocity U_b (lattice units).
    pub bulk_velocity: f64,
    /// Kinematic viscosity ν (lattice units).
    pub nu: f64,
}
impl TurbulentChannelParams {
    /// Construct parameter set from Re_τ and lattice half-width.
    ///
    /// The friction velocity is chosen to keep the Mach number low
    /// (`u_τ = 0.01`), and ν is derived from Re_τ.
    ///
    /// # Arguments
    /// * `re_tau`     — target friction Reynolds number
    /// * `half_width` — channel half-width in lattice nodes
    pub fn new(re_tau: f64, half_width: f64) -> Self {
        let u_tau = 0.01_f64;
        let nu = u_tau * half_width / re_tau;
        let bulk_velocity = u_tau * ((1.0 / KAPPA_VK) * re_tau.ln() + 3.0);
        Self {
            re_tau,
            u_tau,
            half_width,
            bulk_velocity,
            nu,
        }
    }
    /// Viscous length scale δ_ν = ν / u_τ.
    pub fn viscous_length(&self) -> f64 {
        self.nu / self.u_tau.max(1e-30)
    }
    /// Bulk Reynolds number Re_b = U_b * h / ν.
    pub fn re_bulk(&self) -> f64 {
        self.bulk_velocity * self.half_width / self.nu.max(1e-30)
    }
    /// Skin-friction coefficient C_f = 2 τ_w / (ρ U_b²) = 2 (u_τ/U_b)².
    pub fn friction_coefficient(&self) -> f64 {
        2.0 * (self.u_tau / self.bulk_velocity.max(1e-30)).powi(2)
    }
    /// y⁺ value at wall-normal position y (lattice units from wall).
    pub fn y_plus(&self, y: f64) -> f64 {
        y * self.u_tau / self.nu.max(1e-30)
    }
    /// Outer coordinate y/h at wall-normal position y.
    pub fn y_outer(&self, y: f64) -> f64 {
        y / self.half_width.max(1e-30)
    }
}
/// Analysis of a wall-normal velocity profile in inner and outer scaling.
///
/// Given a discrete (y, u) profile and the wall parameters, decomposes
/// the profile into viscous sublayer, buffer layer, log region, and wake.
#[derive(Debug, Clone)]
pub struct VelocityProfileAnalysis {
    /// y⁺ coordinate at each profile point.
    pub y_plus: Vec<f64>,
    /// u⁺ = u / u_τ at each profile point.
    pub u_plus: Vec<f64>,
    /// Outer coordinate y/h at each profile point.
    pub y_outer: Vec<f64>,
    /// Friction velocity u_τ used for scaling.
    pub u_tau: f64,
    /// Kinematic viscosity.
    pub nu: f64,
    /// Channel half-width h.
    pub half_width: f64,
}
impl VelocityProfileAnalysis {
    /// Build an analysis from a raw (y, u) profile and flow parameters.
    ///
    /// # Arguments
    /// * `profile`     — slice of (y_lattice, u_x) pairs
    /// * `u_tau`       — friction velocity
    /// * `nu`          — kinematic viscosity
    /// * `half_width`  — channel half-width
    pub fn from_profile(profile: &[(f64, f64)], u_tau: f64, nu: f64, half_width: f64) -> Self {
        let delta_nu = nu / u_tau.max(1e-30);
        let y_plus: Vec<f64> = profile.iter().map(|(y, _)| y / delta_nu).collect();
        let u_plus: Vec<f64> = profile.iter().map(|(_, u)| u / u_tau.max(1e-30)).collect();
        let y_outer: Vec<f64> = profile
            .iter()
            .map(|(y, _)| y / half_width.max(1e-30))
            .collect();
        Self {
            y_plus,
            u_plus,
            y_outer,
            u_tau,
            nu,
            half_width,
        }
    }
    /// Return u⁺ in the viscous sublayer region (y⁺ ≤ 5).
    pub fn viscous_sublayer_points(&self) -> Vec<(f64, f64)> {
        self.y_plus
            .iter()
            .zip(self.u_plus.iter())
            .filter(|&(yp, _)| *yp <= 5.0)
            .map(|(yp, up)| (*yp, *up))
            .collect()
    }
    /// Return profile points in the logarithmic region (30 ≤ y⁺ ≤ 300).
    pub fn log_region_points(&self) -> Vec<(f64, f64)> {
        self.y_plus
            .iter()
            .zip(self.u_plus.iter())
            .filter(|&(yp, _)| *yp >= 30.0 && *yp <= Y_PLUS_LOG_UPPER)
            .map(|(yp, up)| (*yp, *up))
            .collect()
    }
    /// Return profile points in the outer (wake) region (y/h ≥ 0.3).
    pub fn wake_region_points(&self) -> Vec<(f64, f64)> {
        self.y_outer
            .iter()
            .zip(self.u_plus.iter())
            .filter(|&(yo, _)| *yo >= 0.3)
            .map(|(yo, up)| (*yo, *up))
            .collect()
    }
    /// Log-law u⁺ = (1/κ)·ln(y⁺) + B for a given y⁺.
    pub fn log_law_u_plus(y_plus: f64) -> f64 {
        if y_plus <= 0.0 {
            0.0
        } else if y_plus < Y_PLUS_TRANSITION {
            y_plus
        } else {
            (1.0 / KAPPA_VK) * y_plus.ln() + B_LL
        }
    }
    /// Defect law: u_c⁺ - u⁺ where u_c is the centerline velocity.
    pub fn velocity_defect(&self) -> Vec<(f64, f64)> {
        if self.u_plus.is_empty() {
            return vec![];
        }
        let u_cl = self
            .u_plus
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        self.y_outer
            .iter()
            .zip(self.u_plus.iter())
            .map(|(&yo, &up)| (yo, u_cl - up))
            .collect()
    }
    /// Estimate log-law fit quality: mean squared deviation in log region.
    pub fn log_law_residual(&self) -> f64 {
        let log_pts = self.log_region_points();
        if log_pts.is_empty() {
            return 0.0;
        }
        let mse: f64 = log_pts
            .iter()
            .map(|&(yp, up)| {
                let u_ll = Self::log_law_u_plus(yp);
                (up - u_ll).powi(2)
            })
            .sum::<f64>()
            / log_pts.len() as f64;
        mse.sqrt()
    }
    /// Viscous sublayer linear fit slope (should be ≈ 1 for u⁺ = y⁺).
    pub fn sublayer_slope(&self) -> f64 {
        let pts = self.viscous_sublayer_points();
        if pts.len() < 2 {
            return 1.0;
        }
        let sum_yy: f64 = pts.iter().map(|(y, _)| y * y).sum();
        let sum_yu: f64 = pts.iter().map(|(y, u)| y * u).sum();
        if sum_yy < 1e-30 { 1.0 } else { sum_yu / sum_yy }
    }
}
/// DNS-resolution D2Q9 LBM channel flow solver.
///
/// Integrates a constant body force to drive streamwise flow, accumulates
/// running mean statistics, and provides friction-velocity diagnostics
/// consistent with the Kim–Moin–Moser benchmark.
#[derive(Debug, Clone)]
pub struct DnsDrivenLbm {
    /// Streamwise lattice size.
    pub nx: usize,
    /// Wall-normal lattice size (including solid wall nodes).
    pub ny: usize,
    /// Physical parameter set.
    pub params: TurbulentChannelParams,
    /// Distribution functions f_q at each (i, j) node, row-major, D2Q9.
    pub f_dist: Vec<f64>,
    /// Instantaneous velocity \[(ux, uy)\] at each node.
    pub u: Vec<[f64; 2]>,
    /// Density field at each node.
    pub rho: Vec<f64>,
    /// Accumulated mean streamwise velocity (for statistics).
    pub mean_u: Vec<f64>,
    /// Number of steps accumulated in mean_u.
    pub stat_count: u64,
    /// Constant body force (per unit mass, lattice units).
    pub body_force: f64,
}
impl DnsDrivenLbm {
    /// Create a new DNS-driven LBM channel flow solver.
    ///
    /// # Arguments
    /// * `nx`  — streamwise lattice size
    /// * `ny`  — wall-normal lattice size (walls at j=0 and j=ny-1)
    /// * `re_tau` — target friction Reynolds number
    pub fn new(nx: usize, ny: usize, re_tau: f64) -> Self {
        let half_width = ny as f64 / 2.0;
        let params = TurbulentChannelParams::new(re_tau, half_width);
        let body_force = params.u_tau * params.u_tau / half_width;
        let n = nx * ny;
        let mut f_dist = vec![0.0f64; n * 9];
        for idx in 0..n {
            for q in 0..9 {
                f_dist[idx * 9 + q] = W9[q];
            }
        }
        Self {
            nx,
            ny,
            params,
            f_dist,
            u: vec![[0.0; 2]; n],
            rho: vec![1.0; n],
            mean_u: vec![0.0; ny],
            stat_count: 0,
            body_force,
        }
    }
    /// Flat index for node (i, j).
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }
    /// Initialize with a turbulent-like superposition of parabola + perturbation.
    pub fn init_turbulent_profile(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let u_b = self.params.bulk_velocity;
        for i in 0..nx {
            for j in 0..ny {
                let yc = 2.0 * j as f64 / (ny as f64 - 1.0) - 1.0;
                let y_frac = (1.0 - yc.abs()).max(0.0);
                let ux = u_b * y_frac.powf(1.0 / 7.0);
                let idx = self.index(i, j);
                self.u[idx] = [ux, 0.0];
                for q in 0..9 {
                    self.f_dist[idx * 9 + q] = equilibrium_d2q9(1.0, ux, 0.0, q);
                }
            }
        }
    }
    /// Perform one BGK collision + streaming + bounce-back step.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nu = self.params.nu;
        let tau = 3.0 * nu + 0.5;
        let omega = 1.0 / tau;
        let bf = self.body_force;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let base = idx * 9;
                let (rho, mut ux, uy) = moments_d2q9(&self.f_dist, base);
                if j > 0 && j < ny - 1 {
                    ux += bf;
                }
                self.rho[idx] = rho;
                self.u[idx] = [ux, uy];
                for q in 0..9 {
                    let feq = equilibrium_d2q9(rho, ux, uy, q);
                    self.f_dist[base + q] += omega * (feq - self.f_dist[base + q]);
                }
            }
        }
        let mut f_new = vec![0.0f64; nx * ny * 9];
        for i in 0..nx {
            for j in 0..ny {
                let base = self.index(i, j) * 9;
                for q in 0..9 {
                    let nj = j as i32 + CY[q];
                    if nj < 0 || nj >= ny as i32 {
                        let bb_base = self.index(i, j) * 9;
                        f_new[bb_base + BB[q]] = self.f_dist[base + q];
                    } else {
                        let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
                        let dst = self.index(ni, nj as usize) * 9 + q;
                        f_new[dst] = self.f_dist[base + q];
                    }
                }
            }
        }
        self.f_dist = f_new;
    }
    /// Accumulate statistics for mean velocity profile.
    pub fn accumulate_stats(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for j in 0..ny {
            let mean: f64 = (0..nx).map(|i| self.u[self.index(i, j)][0]).sum::<f64>() / nx as f64;
            self.mean_u[j] += mean;
        }
        self.stat_count += 1;
    }
    /// Return time-averaged mean velocity profile (y, `u`).
    pub fn mean_velocity_profile(&self) -> Vec<(f64, f64)> {
        let ny = self.ny;
        let count = self.stat_count.max(1) as f64;
        (0..ny)
            .map(|j| (j as f64, self.mean_u[j] / count))
            .collect()
    }
    /// Reset accumulated statistics.
    pub fn reset_stats(&mut self) {
        for v in self.mean_u.iter_mut() {
            *v = 0.0;
        }
        self.stat_count = 0;
    }
    /// Return friction velocity from force balance.
    pub fn friction_velocity(&self) -> f64 {
        let h = self.params.half_width;
        let tau_w = self.body_force * h;
        let rho_mean: f64 = self.rho.iter().sum::<f64>() / self.rho.len() as f64;
        (tau_w / rho_mean.max(1e-30)).sqrt()
    }
    /// Bulk streamwise velocity (averaged over interior nodes).
    pub fn bulk_velocity(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut sum = 0.0;
        let mut n = 0usize;
        for i in 0..nx {
            for j in 1..ny - 1 {
                sum += self.u[self.index(i, j)][0];
                n += 1;
            }
        }
        if n == 0 { 0.0 } else { sum / n as f64 }
    }
    /// Friction Reynolds number Re_τ = u_τ h / ν from current state.
    pub fn computed_re_tau(&self) -> f64 {
        let u_tau = self.friction_velocity();
        u_tau * self.params.half_width / self.params.nu.max(1e-30)
    }
}
/// Full 2D Reynolds stress tensor <u_i' u_j'> at a single location.
///
/// Stores the symmetric tensor components in 2D.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReynoldsStressTensor {
    /// <u'u'> streamwise normal stress.
    pub uu: f64,
    /// <v'v'> wall-normal normal stress.
    pub vv: f64,
    /// <u'v'> shear stress (off-diagonal component).
    pub uv: f64,
}
impl ReynoldsStressTensor {
    /// Create from fluctuation samples.
    pub fn from_samples(u_flucts: &[f64], v_flucts: &[f64]) -> Self {
        let n = u_flucts.len().min(v_flucts.len());
        if n == 0 {
            return Self::default();
        }
        let nf = n as f64;
        let uu = u_flucts.iter().map(|&u| u * u).sum::<f64>() / nf;
        let vv = v_flucts.iter().map(|&v| v * v).sum::<f64>() / nf;
        let uv = u_flucts
            .iter()
            .zip(v_flucts.iter())
            .map(|(&u, &v)| u * v)
            .sum::<f64>()
            / nf;
        Self { uu, vv, uv }
    }
    /// Turbulent kinetic energy k = 0.5*(uu + vv).
    pub fn tke(&self) -> f64 {
        0.5 * (self.uu + self.vv)
    }
    /// Anisotropy: ratio of shear to TKE.
    pub fn anisotropy(&self) -> f64 {
        let k = self.tke();
        if k < 1e-30 { 0.0 } else { self.uv.abs() / k }
    }
}
/// Coarse-grid LBM channel flow with log-law wall function boundary.
///
/// The first interior node is treated as the "wall" node at y⁺ ≈ y_plus_first,
/// and the wall shear stress is computed via the log-law rather than resolved.
#[derive(Debug, Clone)]
pub struct WallFunctionLbm {
    /// Streamwise lattice size.
    pub nx: usize,
    /// Wall-normal lattice size (coarse grid).
    pub ny: usize,
    /// Physical parameters.
    pub params: TurbulentChannelParams,
    /// Distribution functions.
    pub f_dist: Vec<f64>,
    /// Velocity field.
    pub u: Vec<[f64; 2]>,
    /// Density field.
    pub rho: Vec<f64>,
    /// y⁺ at the first interior node (used by the wall function).
    pub y_plus_first: f64,
    /// Body force.
    pub body_force: f64,
}
impl WallFunctionLbm {
    /// Create a new wall-function LBM solver.
    ///
    /// # Arguments
    /// * `nx`     — streamwise lattice size
    /// * `ny`     — wall-normal lattice size (coarse)
    /// * `re_tau` — friction Reynolds number
    pub fn new(nx: usize, ny: usize, re_tau: f64) -> Self {
        let half_width = ny as f64 / 2.0;
        let params = TurbulentChannelParams::new(re_tau, half_width);
        let dy = (ny as f64 - 1.0).max(1.0);
        let y_first = half_width / dy;
        let y_plus_first = params.y_plus(y_first);
        let body_force = params.u_tau * params.u_tau / half_width;
        let n = nx * ny;
        let mut f_dist = vec![0.0f64; n * 9];
        for idx in 0..n {
            for q in 0..9 {
                f_dist[idx * 9 + q] = W9[q];
            }
        }
        Self {
            nx,
            ny,
            params,
            f_dist,
            u: vec![[0.0; 2]; n],
            rho: vec![1.0; n],
            y_plus_first,
            body_force,
        }
    }
    /// Compute y⁺ at a wall-normal node j.
    pub fn y_plus_at(&self, j: usize) -> f64 {
        let ny = self.ny;
        let y = j.min(ny - 1 - j) as f64;
        self.params.y_plus(y)
    }
    /// Wall shear stress from log-law: τ_w = ρ·ν·u / (δ_ν·u⁺(y⁺)).
    ///
    /// Given the tangential velocity `u_tan` at y⁺, return τ_w.
    pub fn wall_shear_from_log_law(&self, u_tan: f64, y_plus: f64) -> f64 {
        let u_plus = if y_plus < Y_PLUS_TRANSITION {
            y_plus
        } else {
            (1.0 / KAPPA_VK) * y_plus.ln() + B_LL
        };
        if u_plus < 1e-30 {
            return 0.0;
        }
        let u_tau = u_tan / u_plus;
        u_tau * u_tau
    }
    /// Apply wall-function correction at the near-wall nodes.
    ///
    /// Adjusts the velocity at j=1 and j=ny-2 to match the log-law.
    pub fn apply_wall_function(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let y_plus = self.y_plus_first;
        let u_plus = if y_plus < Y_PLUS_TRANSITION {
            y_plus
        } else {
            (1.0 / KAPPA_VK) * y_plus.ln() + B_LL
        };
        let u_tau = self.params.u_tau;
        let u_wall = u_plus * u_tau;
        for i in 0..nx {
            let idx_bot = i * ny + 1;
            let idx_top = i * ny + (ny - 2);
            self.u[idx_bot][0] = u_wall;
            self.u[idx_top][0] = u_wall;
        }
    }
    /// Perform one LBM step with wall function.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let nu = self.params.nu;
        let tau = 3.0 * nu + 0.5;
        let omega = 1.0 / tau;
        let bf = self.body_force;
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let base = idx * 9;
                let (rho, mut ux, uy) = moments_d2q9(&self.f_dist, base);
                if j > 0 && j < ny - 1 {
                    ux += bf;
                }
                self.rho[idx] = rho;
                self.u[idx] = [ux, uy];
                for q in 0..9 {
                    let feq = equilibrium_d2q9(rho, ux, uy, q);
                    self.f_dist[base + q] += omega * (feq - self.f_dist[base + q]);
                }
            }
        }
        let mut f_new = vec![0.0f64; nx * ny * 9];
        for i in 0..nx {
            for j in 0..ny {
                let base = (i * ny + j) * 9;
                for q in 0..9 {
                    let nj = j as i32 + CY[q];
                    if nj < 0 || nj >= ny as i32 {
                        f_new[(i * ny + j) * 9 + BB[q]] = self.f_dist[base + q];
                    } else {
                        let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
                        let dst = (ni * ny + nj as usize) * 9 + q;
                        f_new[dst] = self.f_dist[base + q];
                    }
                }
            }
        }
        self.f_dist = f_new;
        self.apply_wall_function();
    }
    /// Friction velocity computed from force balance.
    pub fn friction_velocity(&self) -> f64 {
        let h = self.params.half_width;
        let tau_w = self.body_force * h;
        (tau_w).sqrt()
    }
    /// Bulk velocity.
    pub fn bulk_velocity(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in 0..nx {
            for j in 1..ny - 1 {
                sum += self.u[i * ny + j][0];
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
}
/// Validation of LBM channel flow results against Kim–Moin–Moser DNS data.
///
/// Provides reference profiles interpolated from KMM (1987) at Re_τ = 180
/// and MKM (1999) at Re_τ = 590.
#[derive(Debug, Clone)]
pub struct ChannelFlowValidation {
    /// Friction Reynolds number of the reference DNS.
    pub re_tau_dns: f64,
    /// Reference y⁺ points from DNS.
    pub dns_y_plus: Vec<f64>,
    /// Reference mean u⁺ from DNS.
    pub dns_u_plus: Vec<f64>,
    /// Reference u_rms⁺ from DNS.
    pub dns_u_rms: Vec<f64>,
}
impl ChannelFlowValidation {
    /// Create a reference dataset mimicking KMM Re_τ = 180 behaviour.
    ///
    /// Uses the analytical approximation: log-law + parabolic wake.
    pub fn kmm_re180() -> Self {
        let n = 30usize;
        let y_plus: Vec<f64> = (1..=n).map(|k| (k as f64) * 6.0).collect();
        let u_plus: Vec<f64> = y_plus
            .iter()
            .map(|&yp| {
                if yp < Y_PLUS_TRANSITION {
                    yp
                } else {
                    (1.0 / KAPPA_VK) * yp.ln() + B_LL
                }
            })
            .collect();
        let u_rms: Vec<f64> = y_plus
            .iter()
            .map(|&yp| {
                let peak_yp = 15.0_f64;
                let sigma = 25.0_f64;
                2.7 * (-((yp - peak_yp) / sigma).powi(2)).exp()
            })
            .collect();
        Self {
            re_tau_dns: 180.0,
            dns_y_plus: y_plus,
            dns_u_plus: u_plus,
            dns_u_rms: u_rms,
        }
    }
    /// Create a reference dataset mimicking MKM Re_τ = 590.
    pub fn mkm_re590() -> Self {
        let n = 40usize;
        let y_plus: Vec<f64> = (1..=n).map(|k| (k as f64) * 15.0).collect();
        let u_plus: Vec<f64> = y_plus
            .iter()
            .map(|&yp| {
                if yp < Y_PLUS_TRANSITION {
                    yp
                } else {
                    (1.0 / KAPPA_VK) * yp.ln() + B_LL
                }
            })
            .collect();
        let u_rms: Vec<f64> = y_plus
            .iter()
            .map(|&yp| {
                let peak_yp = 15.0_f64;
                let sigma = 30.0_f64;
                2.9 * (-((yp - peak_yp) / sigma).powi(2)).exp()
            })
            .collect();
        Self {
            re_tau_dns: 590.0,
            dns_y_plus: y_plus,
            dns_u_plus: u_plus,
            dns_u_rms: u_rms,
        }
    }
    /// Interpolate reference u⁺ at a given y⁺ (linear interpolation).
    pub fn interpolate_u_plus(&self, y_plus: f64) -> f64 {
        let pts = &self.dns_y_plus;
        let us = &self.dns_u_plus;
        if pts.is_empty() {
            return 0.0;
        }
        if y_plus <= pts[0] {
            return us[0];
        }
        let n = pts.len();
        if y_plus >= pts[n - 1] {
            return us[n - 1];
        }
        for k in 0..n - 1 {
            if y_plus >= pts[k] && y_plus < pts[k + 1] {
                let t = (y_plus - pts[k]) / (pts[k + 1] - pts[k]);
                return us[k] + t * (us[k + 1] - us[k]);
            }
        }
        us[n - 1]
    }
    /// Compare an LBM profile against DNS: return L2 error in u⁺.
    ///
    /// # Arguments
    /// * `lbm_y_plus` — y⁺ values from LBM
    /// * `lbm_u_plus` — u⁺ values from LBM
    pub fn l2_error_u_plus(&self, lbm_y_plus: &[f64], lbm_u_plus: &[f64]) -> f64 {
        if lbm_y_plus.is_empty() {
            return 0.0;
        }
        let mse: f64 = lbm_y_plus
            .iter()
            .zip(lbm_u_plus.iter())
            .map(|(&yp, &up)| {
                let u_dns = self.interpolate_u_plus(yp);
                (up - u_dns).powi(2)
            })
            .sum::<f64>()
            / lbm_y_plus.len() as f64;
        mse.sqrt()
    }
    /// Friction coefficient C_f = 2 / (U_b⁺)² where U_b⁺ = U_b / u_τ.
    pub fn friction_coefficient(u_b_plus: f64) -> f64 {
        if u_b_plus < 1e-30 {
            return 0.0;
        }
        2.0 / (u_b_plus * u_b_plus)
    }
    /// Dean (1978) correlation: C_f ≈ 0.073 Re_b^{-0.25}.
    pub fn dean_friction_coefficient(re_bulk: f64) -> f64 {
        0.073 * re_bulk.powf(-0.25)
    }
}
/// RANS-LBM: Reynolds-Averaged Navier–Stokes coupled with LBM.
///
/// Uses a mixing-length eddy-viscosity model to close the turbulence.
/// The effective relaxation time τ_eff = 0.5 + 3(ν + ν_t) where ν_t
/// is the turbulent (eddy) viscosity from the mixing-length model.
#[derive(Debug, Clone)]
pub struct ReynoldsAveragedLbm {
    /// Streamwise lattice size.
    pub nx: usize,
    /// Wall-normal lattice size.
    pub ny: usize,
    /// Physical parameters.
    pub params: TurbulentChannelParams,
    /// Distribution functions.
    pub f_dist: Vec<f64>,
    /// Velocity field.
    pub u: Vec<[f64; 2]>,
    /// Density field.
    pub rho: Vec<f64>,
    /// Turbulent (eddy) viscosity field.
    pub nu_t: Vec<f64>,
    /// Body force.
    pub body_force: f64,
    /// Van Driest damping constant A⁺.
    pub a_plus: f64,
}
impl ReynoldsAveragedLbm {
    /// Create a new RANS-LBM solver with Van Driest damping.
    ///
    /// # Arguments
    /// * `nx`     — streamwise lattice size
    /// * `ny`     — wall-normal lattice size
    /// * `re_tau` — friction Reynolds number
    pub fn new(nx: usize, ny: usize, re_tau: f64) -> Self {
        let half_width = ny as f64 / 2.0;
        let params = TurbulentChannelParams::new(re_tau, half_width);
        let body_force = params.u_tau * params.u_tau / half_width;
        let n = nx * ny;
        let mut f_dist = vec![0.0f64; n * 9];
        for idx in 0..n {
            for q in 0..9 {
                f_dist[idx * 9 + q] = W9[q];
            }
        }
        Self {
            nx,
            ny,
            params,
            f_dist,
            u: vec![[0.0; 2]; n],
            rho: vec![1.0; n],
            nu_t: vec![0.0; n],
            body_force,
            a_plus: 26.0,
        }
    }
    /// Compute mixing-length at wall-normal position j.
    ///
    /// Uses Van Driest damping: l_m = κ·y⁺·(1 - exp(-y⁺/A⁺)).
    pub fn mixing_length(&self, j: usize) -> f64 {
        let y = j.min(self.ny - 1 - j) as f64;
        let y_plus = self.params.y_plus(y);
        let damp = 1.0 - (-y_plus / self.a_plus).exp();
        KAPPA_VK * y * damp
    }
    /// Update turbulent viscosity field from local velocity gradients.
    ///
    /// ν_t = l_m² |∂u/∂y|
    pub fn update_nu_t(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 1..ny - 1 {
                let u_up = self.u[i * ny + j + 1][0];
                let u_dn = self.u[i * ny + j - 1][0];
                let du_dy = (u_up - u_dn) * 0.5;
                let lm = self.mixing_length(j);
                self.nu_t[i * ny + j] = lm * lm * du_dy.abs();
            }
        }
    }
    /// Perform one RANS-LBM timestep with position-dependent relaxation.
    pub fn step(&mut self) {
        self.update_nu_t();
        let nx = self.nx;
        let ny = self.ny;
        let nu = self.params.nu;
        let bf = self.body_force;
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let base = idx * 9;
                let nu_eff = nu + self.nu_t[idx];
                let tau = 3.0 * nu_eff + 0.5;
                let omega = 1.0 / tau;
                let (rho, mut ux, uy) = moments_d2q9(&self.f_dist, base);
                if j > 0 && j < ny - 1 {
                    ux += bf;
                }
                self.rho[idx] = rho;
                self.u[idx] = [ux, uy];
                for q in 0..9 {
                    let feq = equilibrium_d2q9(rho, ux, uy, q);
                    self.f_dist[base + q] += omega * (feq - self.f_dist[base + q]);
                }
            }
        }
        let mut f_new = vec![0.0f64; nx * ny * 9];
        for i in 0..nx {
            for j in 0..ny {
                let base = (i * ny + j) * 9;
                for q in 0..9 {
                    let nj = j as i32 + CY[q];
                    if nj < 0 || nj >= ny as i32 {
                        f_new[(i * ny + j) * 9 + BB[q]] = self.f_dist[base + q];
                    } else {
                        let ni = (i as i32 + CX[q]).rem_euclid(nx as i32) as usize;
                        let dst = (ni * ny + nj as usize) * 9 + q;
                        f_new[dst] = self.f_dist[base + q];
                    }
                }
            }
        }
        self.f_dist = f_new;
    }
    /// Mean eddy viscosity at wall-normal row j.
    pub fn mean_nu_t(&self, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        (0..nx).map(|i| self.nu_t[i * ny + j]).sum::<f64>() / nx as f64
    }
    /// Eddy diffusivity at wall-normal row j (ν_t / Pr_t, Pr_t = 0.9).
    pub fn eddy_diffusivity(&self, j: usize) -> f64 {
        self.mean_nu_t(j) / 0.9
    }
    /// Bulk streamwise velocity.
    pub fn bulk_velocity(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in 0..nx {
            for j in 1..ny - 1 {
                sum += self.u[i * ny + j][0];
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }
}
