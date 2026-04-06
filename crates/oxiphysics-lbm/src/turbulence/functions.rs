//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

use crate::grid::{LbmGrid2D, equilibrium_2d};
use crate::lattice::CS2;

/// Compute the local effective omega at cell (x, y) using the Smagorinsky model.
///
/// The non-equilibrium stress tensor is computed from `f_neq = f - f_eq`.
/// The strain rate magnitude `|S|` is derived from the second-order moment
/// of the non-equilibrium part. The effective relaxation time is then:
///
/// `tau_eff = 0.5 * (tau + sqrt(tau² + 18 * cs_smag² * |S| / rho))`
///
/// Returns `omega_eff = 1 / tau_eff`.
pub fn smagorinsky_omega(
    grid: &LbmGrid2D,
    cs_smag: f64,
    base_omega: f64,
    x: usize,
    y: usize,
) -> f64 {
    let k = grid.idx(x, y);
    let q = grid.lattice.q();
    let rho = grid.rho[k];
    let ux = grid.ux[k];
    let uy = grid.uy[k];
    let mut pi_xx = 0.0;
    let mut pi_xy = 0.0;
    let mut pi_yy = 0.0;
    for i in 0..q {
        let w = grid.lattice.weight(i);
        let c = grid.lattice.velocity_2d(i);
        let cx = c[0] as f64;
        let cy = c[1] as f64;
        let feq = equilibrium_2d(w, rho, ux, uy, cx, cy);
        let fneq = grid.f[i][k] - feq;
        pi_xx += fneq * cx * cx;
        pi_xy += fneq * cx * cy;
        pi_yy += fneq * cy * cy;
    }
    let pi_mag = (pi_xx * pi_xx + 2.0 * pi_xy * pi_xy + pi_yy * pi_yy).sqrt();
    let tau = 1.0 / base_omega;
    let tau_eff =
        0.5 * (tau + (tau * tau + 18.0 * cs_smag * cs_smag * pi_mag / (rho * CS2 * CS2)).sqrt());
    1.0 / tau_eff
}
/// Perform a BGK collision step with Smagorinsky turbulence model on a 2D grid.
///
/// Uses a locally varying omega computed from the Smagorinsky model.
pub fn bgk_collide_smagorinsky_2d(grid: &mut LbmGrid2D, base_omega: f64, cs_smag: f64) {
    let nx = grid.nx;
    let ny = grid.ny;
    let q = grid.lattice.q();
    grid.compute_macroscopic();
    for y in 0..ny {
        for x in 0..nx {
            let omega = smagorinsky_omega(grid, cs_smag, base_omega, x, y);
            let k = grid.idx(x, y);
            let rho_k = grid.rho[k];
            let ux_k = grid.ux[k];
            let uy_k = grid.uy[k];
            for i in 0..q {
                let w = grid.lattice.weight(i);
                let c = grid.lattice.velocity_2d(i);
                let feq = equilibrium_2d(w, rho_k, ux_k, uy_k, c[0] as f64, c[1] as f64);
                grid.f[i][k] -= omega * (grid.f[i][k] - feq);
            }
        }
    }
}
/// Compute the effective strain rate magnitude from D3Q19 non-equilibrium
/// distributions using the second-order stress tensor approach.
///
/// `f_neq` must contain `f_i - f_eq_i` for all 19 directions.
/// Returns `|S| = sqrt(2 * S_ab * S_ab)` where
/// `S_ab = -1/(2 * rho * cs² * tau) * Pi_ab_neq`.
#[allow(dead_code)]
pub fn compute_strain_rate_d3q19(f_neq: &[f64; 19], rho: f64) -> f64 {
    use crate::lattice::D3Q19_VELOCITIES;
    let mut pi = [[0.0f64; 3]; 3];
    for i in 0..19 {
        let c = D3Q19_VELOCITIES[i];
        for a in 0..3 {
            for b in 0..3 {
                pi[a][b] += f_neq[i] * c[a] as f64 * c[b] as f64;
            }
        }
    }
    let mut norm_sq = 0.0f64;
    for a in 0..3 {
        for b in 0..3 {
            norm_sq += pi[a][b] * pi[a][b];
        }
    }
    let cs4 = CS2 * CS2;
    (norm_sq / (2.0 * rho * rho * cs4)).sqrt()
}
/// Compute the effective omega using the 3D Smagorinsky formula.
///
/// `tau_eff = 0.5 * (tau_base + sqrt(tau_base² + 18*(Cs*dx)²*|S|))`
///
/// `dx` is the lattice spacing (typically 1.0 in LBM units),
/// `dt` is the time step (typically 1.0).
#[allow(dead_code)]
pub fn compute_effective_omega(omega_base: f64, strain_rate: f64, dx: f64, _dt: f64) -> f64 {
    let tau = 1.0 / omega_base;
    let cs_dx = 0.1 * dx;
    let discriminant = tau * tau + 18.0 * cs_dx * cs_dx * strain_rate;
    let tau_eff = 0.5 * (tau + discriminant.max(0.0).sqrt());
    1.0 / tau_eff
}
/// Compute turbulent (eddy) viscosity from strain rate and mixing length.
///
/// `nu_t = l_mix² * |S|`
///
/// where `l_mix = Cs * Delta` is the Smagorinsky mixing length.
#[allow(dead_code)]
pub fn compute_turbulent_viscosity(strain_rate: f64, mixing_length: f64) -> f64 {
    mixing_length * mixing_length * strain_rate
}
/// Compute the WALE (Wall-Adaptive Large-Eddy Simulation) eddy viscosity.
///
/// Given a 3x3 velocity-gradient tensor `g_ij = du_i/dx_j` (row-major),
/// the WALE model computes:
///
/// `nu_t = (Cw * Delta)^2 * (S^d_ij S^d_ij)^(3/2) / (S_ij S_ij)^(5/2) + (S^d_ij S^d_ij)^(5/4)`
///
/// where `S_ij = (g_ij + g_ji)/2` is the symmetric rate-of-strain tensor and
/// `S^d_ij` is its traceless symmetric square.
///
/// The `velocity_gradient` argument is a flat 9-element array in row-major order
/// (`[du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx, dw/dy, dw/dz]`).
///
/// Returns `(Cw * Delta)^2` scaled numerator / denominator; the caller should
/// multiply by `(Cw * Delta)^2` to get the dimensional eddy viscosity.
#[allow(dead_code)]
pub fn wale_model(velocity_gradient: &[f64; 9]) -> f64 {
    let g = |i: usize, j: usize| velocity_gradient[i * 3 + j];
    let s = |i: usize, j: usize| 0.5 * (g(i, j) + g(j, i));
    let mut sd = [[0.0f64; 3]; 3];
    let mut trace_gsq = 0.0f64;
    let mut gsq = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut val = 0.0f64;
            for k in 0..3 {
                val += g(i, k) * g(k, j);
            }
            gsq[i][j] = val;
            if i == j {
                trace_gsq += val;
            }
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let sym_gsq = 0.5 * (gsq[i][j] + gsq[j][i]);
            let delta_ij = if i == j { 1.0 } else { 0.0 };
            sd[i][j] = sym_gsq - (1.0 / 3.0) * delta_ij * trace_gsq;
        }
    }
    let mut sd_sq = 0.0f64;
    let mut s_sq = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            sd_sq += sd[i][j] * sd[i][j];
            let sij = s(i, j);
            s_sq += sij * sij;
        }
    }
    let numerator = sd_sq.powf(1.5);
    let denominator = s_sq.powf(2.5) + sd_sq.powf(1.25) + 1e-30;
    numerator / denominator
}
/// Compute eddy viscosity using the Vreman sub-grid scale model.
///
/// The Vreman model is defined as:
/// `nu_t = Cv * sqrt(B_beta / (alpha_ij * alpha_ij))`
///
/// where `alpha_ij = du_i/dx_j` and `B_beta` involves the invariants of
/// `beta_ij = Delta_m^2 * alpha_mi * alpha_mj`.
///
/// `velocity_gradient` is row-major `[du/dx, du/dy, du/dz, dv/dx, ...]`.
/// `cv` is the Vreman constant (typically ~0.07).
#[allow(dead_code)]
pub fn vreman_model(velocity_gradient: &[f64; 9], cv: f64) -> f64 {
    let g = |i: usize, j: usize| velocity_gradient[i * 3 + j];
    let mut alpha_sq = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let aij = g(i, j);
            alpha_sq += aij * aij;
        }
    }
    if alpha_sq < 1e-30 {
        return 0.0;
    }
    let mut beta = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for m in 0..3 {
                beta[i][j] += g(m, i) * g(m, j);
            }
        }
    }
    let b_beta = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] + beta[0][0] * beta[2][2]
        - beta[0][2] * beta[0][2]
        + beta[1][1] * beta[2][2]
        - beta[1][2] * beta[1][2];
    if b_beta <= 0.0 {
        return 0.0;
    }
    cv * (b_beta / alpha_sq).sqrt()
}
/// Compute the friction velocity `u_tau` from the log-law wall function.
///
/// Uses Newton iteration to solve `u+ = (1/kappa) * ln(y+) + B`
/// where `u+ = u/u_tau` and `y+ = y * u_tau / nu`.
///
/// Returns the friction velocity `u_tau`.
#[allow(dead_code)]
pub fn wall_function_log_law(u_tangential: f64, y: f64, nu: f64) -> f64 {
    pub(super) const KAPPA: f64 = 0.41;
    pub(super) const B: f64 = 5.2;
    if u_tangential.abs() < 1e-30 || y < 1e-30 || nu < 1e-30 {
        return 0.0;
    }
    let mut u_tau = (nu * u_tangential.abs() / y).sqrt();
    for _ in 0..20 {
        if u_tau < 1e-30 {
            u_tau = 1e-10;
        }
        let y_plus = y * u_tau / nu;
        let u_plus = u_tangential.abs() / u_tau;
        let f_val = u_plus - (1.0 / KAPPA) * (y_plus.max(1e-10)).ln() - B;
        let df = -u_tangential.abs() / (u_tau * u_tau) - 1.0 / (KAPPA * u_tau);
        let delta = f_val / df;
        u_tau -= delta;
        u_tau = u_tau.max(1e-15);
        if delta.abs() < 1e-10 * u_tau {
            break;
        }
    }
    u_tau
}
/// Wall shear stress from friction velocity: `tau_w = rho * u_tau^2`.
#[allow(dead_code)]
pub fn wall_shear_stress(rho: f64, u_tau: f64) -> f64 {
    rho * u_tau * u_tau
}
/// Compute y+ (dimensionless wall distance).
#[allow(dead_code)]
pub fn y_plus(y: f64, u_tau: f64, nu: f64) -> f64 {
    y * u_tau / nu
}
/// Turbulent boundary layer velocity profile using Spalding's law.
///
/// `u+ = y+ + exp(-kappa*B) * [exp(kappa*u+) - 1 - kappa*u+ - (kappa*u+)^2/2 - (kappa*u+)^3/6]`
///
/// Given y+, find u+ iteratively.
#[allow(dead_code)]
pub fn spalding_profile(y_plus_val: f64) -> f64 {
    pub(super) const KAPPA: f64 = 0.41;
    pub(super) const B: f64 = 5.2;
    let mut u_plus = if y_plus_val < 11.0 {
        y_plus_val
    } else {
        (1.0 / KAPPA) * y_plus_val.ln() + B
    };
    let exp_kb = (-KAPPA * B).exp();
    for _ in 0..30 {
        let ku = KAPPA * u_plus;
        let rhs = u_plus + exp_kb * (ku.exp() - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0);
        let drhs = 1.0 + exp_kb * (KAPPA * ku.exp() - KAPPA - KAPPA * ku - KAPPA * ku * ku / 2.0);
        let delta = (rhs - y_plus_val) / drhs;
        u_plus -= delta;
        if delta.abs() < 1e-10 {
            break;
        }
    }
    u_plus
}
/// Compute turbulent kinetic energy from velocity fluctuation components.
///
/// `TKE = 0.5 * (u'^2 + v'^2 + w'^2)`
///
/// `fluctuations` is a 3-element array `[u', v', w']`.
#[allow(dead_code)]
pub fn turbulent_kinetic_energy(fluctuations: &[f64; 3]) -> f64 {
    0.5 * (fluctuations[0].powi(2) + fluctuations[1].powi(2) + fluctuations[2].powi(2))
}
/// Turbulence intensity: TI = u_rms / U_mean.
///
/// Returns `f64::INFINITY` when `u_mean` is zero (caller should guard).
#[allow(dead_code)]
pub fn turbulence_intensity(u_rms: f64, u_mean: f64) -> f64 {
    if u_mean.abs() < 1e-30 {
        f64::INFINITY
    } else {
        u_rms / u_mean
    }
}
/// Estimate the one-sided power spectrum of a signal using a simple DFT.
///
/// Returns `n/2` values corresponding to frequencies `k / n` for `k = 0..n/2`.
/// This is an O(n²) implementation suitable for moderate-length signals.
#[allow(dead_code)]
pub fn estimate_power_spectrum(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let half = n / 2;
    let mut spectrum = vec![0.0_f64; half];
    use std::f64::consts::PI;
    for k in 0..half {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for (j, &s) in signal.iter().enumerate() {
            let angle = 2.0 * PI * k as f64 * j as f64 / n as f64;
            re += s * angle.cos();
            im -= s * angle.sin();
        }
        spectrum[k] = (re * re + im * im) / (n as f64 * n as f64);
    }
    spectrum
}
/// Turbulent kinetic energy production: `P_k = nu_t * S^2`
///
/// where `S^2 = 2 * S_ij * S_ij` is twice the resolved strain rate invariant.
#[allow(dead_code)]
pub fn tke_production(nu_t: f64, s_sq: f64) -> f64 {
    nu_t * s_sq
}
/// Turbulent diffusion of TKE: `D_t = sigma_k * nu_t * d^2k/dx^2`
///
/// Simplified 1D version using second derivative estimate.
#[allow(dead_code)]
pub fn tke_turbulent_diffusion(nu_t: f64, k_laplacian: f64, sigma_k: f64) -> f64 {
    sigma_k * nu_t * k_laplacian
}
/// Taylor microscale: `lambda = sqrt(15 * nu * u_rms^2 / epsilon)`.
///
/// Represents the length scale at which viscous effects begin to affect the
/// turbulent velocity fluctuations.
#[allow(dead_code)]
pub fn taylor_microscale(u_rms: f64, nu: f64, epsilon: f64) -> f64 {
    if epsilon < 1e-30 {
        return f64::INFINITY;
    }
    (15.0 * nu * u_rms * u_rms / epsilon).sqrt()
}
/// Kolmogorov length scale: `eta = (nu^3 / epsilon)^0.25`.
///
/// The smallest scale in turbulent flow where viscosity dominates.
#[allow(dead_code)]
pub fn kolmogorov_length_scale(nu: f64, epsilon: f64) -> f64 {
    if epsilon < 1e-30 {
        return f64::INFINITY;
    }
    (nu * nu * nu / epsilon).powf(0.25)
}
/// Kolmogorov time scale: `tau_eta = sqrt(nu / epsilon)`.
#[allow(dead_code)]
pub fn kolmogorov_time_scale(nu: f64, epsilon: f64) -> f64 {
    if epsilon < 1e-30 {
        return f64::INFINITY;
    }
    (nu / epsilon).sqrt()
}
/// Kolmogorov velocity scale: `u_eta = (nu * epsilon)^0.25`.
#[allow(dead_code)]
pub fn kolmogorov_velocity_scale(nu: f64, epsilon: f64) -> f64 {
    (nu * epsilon).powf(0.25)
}
/// Integral length scale: `L = C_mu^{3/4} * k^{3/2} / epsilon`.
///
/// Uses the standard coefficient `C_mu^{3/4}` absorbed into the leading factor.
/// For simplicity this returns `k^{3/2} / epsilon`.
#[allow(dead_code)]
pub fn integral_length_scale(k: f64, epsilon: f64) -> f64 {
    if epsilon < 1e-30 {
        return f64::INFINITY;
    }
    k.powf(1.5) / epsilon
}
/// Turbulent Reynolds number: `Re_t = k^2 / (nu * epsilon)`.
#[allow(dead_code)]
pub fn turbulent_reynolds_number(k: f64, nu: f64, epsilon: f64) -> f64 {
    if nu < 1e-30 || epsilon < 1e-30 {
        return 0.0;
    }
    k * k / (nu * epsilon)
}
/// Kolmogorov inertial-subrange energy spectrum: `E(k) = C_k * epsilon^{2/3} * k^{-5/3}`.
///
/// `wavenumber` is the wavenumber magnitude, `epsilon` is the dissipation rate,
/// and `c_k` is the Kolmogorov constant (typically ~1.5).
#[allow(dead_code)]
pub fn kolmogorov_energy_spectrum(wavenumber: f64, epsilon: f64, c_k: f64) -> f64 {
    if wavenumber < 1e-30 {
        return 0.0;
    }
    c_k * epsilon.powf(2.0 / 3.0) * wavenumber.powf(-5.0 / 3.0)
}
/// Von Kármán energy spectrum (interpolation between energy-containing and inertial ranges).
///
/// `E(k) = alpha * epsilon^{2/3} * k^{-5/3} * (k * L)^4 / (1 + (k*L)^2)^{17/6}`
///
/// where `L` is the integral length scale and `alpha` ≈ 1.5.
#[allow(dead_code)]
pub fn von_karman_spectrum(wavenumber: f64, epsilon: f64, integral_scale: f64) -> f64 {
    pub(super) const ALPHA: f64 = 1.5;
    if wavenumber < 1e-30 || integral_scale < 1e-30 {
        return 0.0;
    }
    let kl = wavenumber * integral_scale;
    let numerator = kl.powi(4);
    let denominator = (1.0 + kl * kl).powf(17.0 / 6.0);
    ALPHA * epsilon.powf(2.0 / 3.0) * wavenumber.powf(-5.0 / 3.0) * numerator / denominator
}
/// Compute the strain-rate magnitude from a velocity-gradient tensor.
///
/// `|S| = sqrt(S_ij S_ij)` where `S_ij = (g_ij + g_ji) / 2`.
///
/// `velocity_gradient` is row-major `[du/dx, du/dy, du/dz, dv/dx, ...]`.
#[allow(dead_code)]
pub fn strain_rate_magnitude(velocity_gradient: &[f64; 9]) -> f64 {
    let g = |i: usize, j: usize| velocity_gradient[i * 3 + j];
    let mut s_sq = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let sij = 0.5 * (g(i, j) + g(j, i));
            s_sq += sij * sij;
        }
    }
    s_sq.sqrt()
}
/// Turbulent Prandtl number: `Pr_t = nu_t / alpha_t` where `alpha_t` is the
/// turbulent thermal diffusivity.
#[allow(dead_code)]
pub fn turbulent_prandtl_number(nu_t: f64, alpha_t: f64) -> f64 {
    if alpha_t.abs() < 1e-30 {
        return 0.0;
    }
    nu_t / alpha_t
}
/// Eddy viscosity from k-omega model: `nu_t = k / omega`.
#[allow(dead_code)]
pub fn k_omega_eddy_viscosity(k: f64, omega: f64) -> f64 {
    if omega.abs() < 1e-30 {
        return 0.0;
    }
    k / omega
}
/// k-omega specific dissipation rate from epsilon: `omega = epsilon / (C_mu * k)`.
#[allow(dead_code)]
pub fn k_omega_from_k_epsilon(k: f64, epsilon: f64) -> f64 {
    pub(super) const C_MU: f64 = 0.09;
    if k.abs() < 1e-30 {
        return 0.0;
    }
    epsilon / (C_MU * k)
}
/// k-omega production term: `P_omega = alpha * omega / k * P_k`.
#[allow(dead_code)]
pub fn k_omega_production(k: f64, omega: f64, p_k: f64) -> f64 {
    pub(super) const ALPHA: f64 = 5.0 / 9.0;
    if k.abs() < 1e-30 {
        return 0.0;
    }
    ALPHA * omega / k * p_k
}
/// k-omega destruction term: `D_omega = beta * omega^2`.
#[allow(dead_code)]
pub fn k_omega_destruction(omega: f64) -> f64 {
    pub(super) const BETA: f64 = 3.0 / 40.0;
    BETA * omega * omega
}
/// Estimate backscatter (reverse energy transfer to resolved scales) using
/// the spectral eddy viscosity concept.
///
/// Returns a non-dimensional backscatter indicator ≥ 0.
/// In practice, backscatter occurs when the local SGS stress does negative work.
/// Here we return the absolute difference between the symmetric and anti-symmetric
/// parts of the velocity gradient, normalized by the Frobenius norm.
#[allow(dead_code)]
pub fn backscatter_indicator(velocity_gradient: &[f64; 9]) -> f64 {
    let g = |i: usize, j: usize| velocity_gradient[i * 3 + j];
    let mut sym_sq = 0.0f64;
    let mut asym_sq = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let sij = 0.5 * (g(i, j) + g(j, i));
            let oij = 0.5 * (g(i, j) - g(j, i));
            sym_sq += sij * sij;
            asym_sq += oij * oij;
        }
    }
    (asym_sq - sym_sq).abs() / (sym_sq + asym_sq + 1e-30)
}
/// Compute the SST blending function F1 for the k-ω SST model.
///
/// F1 blends between the k-ω model (near the wall, F1 → 1) and the k-ε model
/// (in the free stream, F1 → 0).
///
/// `F1 = tanh(arg1^4)`
///
/// where:
/// `arg1 = min( max(sqrt(k)/(C_mu*omega*y), 500*nu/(omega*y^2)),
///              4*sigma_w2*k/(CD_kw*y^2) )`
///
/// and `CD_kw = max(2*sigma_w2*(dk/dx_j * domega/dx_j)/omega, 1e-10)`.
///
/// # Arguments
/// - `k`:      turbulent kinetic energy
/// - `omega`:  specific dissipation rate
/// - `y`:      wall-normal distance
/// - `nu`:     kinematic viscosity
/// - `dk_domega_dot`: dot product `grad(k) · grad(omega)` (cross-diffusion term)
#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn sst_blending_f1(k: f64, omega: f64, y: f64, nu: f64, dk_domega_dot: f64) -> f64 {
    pub(super) const C_MU: f64 = 0.09;
    pub(super) const SIGMA_W2: f64 = 0.856;
    pub(super) const NU_DISS: f64 = 500.0;
    if omega < 1e-30 || y < 1e-30 {
        return 1.0;
    }
    let sqrt_k = k.max(0.0).sqrt();
    let term1 = sqrt_k / (C_MU * omega * y);
    let term2 = NU_DISS * nu / (omega * y * y);
    let CD_kw = (2.0 * SIGMA_W2 * dk_domega_dot / omega.max(1e-30)).max(1e-10);
    let term3 = 4.0 * SIGMA_W2 * k / (CD_kw * y * y);
    let arg1 = term1.max(term2).min(term3);
    arg1.powi(4).tanh()
}
/// Compute the SST blending function F2.
///
/// `F2 = tanh(arg2^2)` where
/// `arg2 = max( 2*sqrt(k)/(C_mu*omega*y), 500*nu/(omega*y^2) )`
///
/// F2 is used to limit the eddy-viscosity near walls.
#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn sst_blending_f2(k: f64, omega: f64, y: f64, nu: f64) -> f64 {
    pub(super) const C_MU: f64 = 0.09;
    pub(super) const NU_DISS: f64 = 500.0;
    if omega < 1e-30 || y < 1e-30 {
        return 1.0;
    }
    let sqrt_k = k.max(0.0).sqrt();
    let term1 = 2.0 * sqrt_k / (C_MU * omega * y);
    let term2 = NU_DISS * nu / (omega * y * y);
    let arg2 = term1.max(term2);
    arg2.powi(2).tanh()
}
/// Blend a pair of k-ε and k-ω coefficients using the SST F1 function.
///
/// `phi_sst = F1 * phi_kw + (1 - F1) * phi_ke`
#[allow(dead_code)]
pub fn sst_blend_coefficient(phi_kw: f64, phi_ke: f64, f1: f64) -> f64 {
    f1 * phi_kw + (1.0 - f1) * phi_ke
}
/// Compute the k-ω SST eddy viscosity with the Bradshaw-limitation.
///
/// `nu_t = a1 * k / max(a1 * omega, S * F2)`
///
/// where `a1 = 0.31`, `S` is the strain-rate magnitude, and `F2` limits the
/// eddy viscosity in adverse pressure gradient regions.
#[allow(dead_code)]
pub fn sst_eddy_viscosity(k: f64, omega: f64, strain_rate_s: f64, f2: f64) -> f64 {
    pub(super) const A1: f64 = 0.31;
    let denom = (A1 * omega).max(strain_rate_s * f2);
    if denom < 1e-30 {
        return 0.0;
    }
    A1 * k / denom
}
/// Evaluate the log-law velocity profile `u+ = (1/kappa) * ln(y+) + B`.
///
/// Valid in the log-law region (approximately `y+ > 30`).
///
/// # Arguments
/// - `y_plus`: dimensionless wall distance
/// - `kappa`:  von Kármán constant (default 0.41)
/// - `B`:      log-law intercept (default 5.2)
///
/// Returns `u+`.
#[allow(dead_code)]
pub fn log_law_u_plus(y_plus: f64, kappa: f64, b: f64) -> f64 {
    if y_plus < 1e-30 {
        return 0.0;
    }
    (1.0 / kappa) * y_plus.ln() + b
}
/// Evaluate the composite wall law (viscous sub-layer + log-law).
///
/// - For `y+ < y_plus_match`: returns `y+` (linear sub-layer)
/// - For `y+ >= y_plus_match`: returns `(1/kappa)*ln(y+) + B`
///
/// `y_plus_match` is typically set to 11.3 (intersection of sub-layer and log-law).
#[allow(dead_code)]
pub fn composite_wall_law(y_plus: f64) -> f64 {
    pub(super) const KAPPA: f64 = 0.41;
    pub(super) const B: f64 = 5.2;
    pub(super) const Y_MATCH: f64 = 11.3;
    if y_plus < Y_MATCH {
        y_plus
    } else {
        log_law_u_plus(y_plus, KAPPA, B)
    }
}
/// Compute the friction Reynolds number `Re_tau = delta * u_tau / nu`.
#[allow(dead_code)]
pub fn friction_reynolds_number(delta: f64, u_tau: f64, nu: f64) -> f64 {
    if nu < 1e-30 {
        return 0.0;
    }
    delta * u_tau / nu
}
/// Apply a sharp spectral cutoff filter to a 1D signal in physical space.
///
/// Modes with wavenumber `|k| > k_cutoff` are zeroed out.
/// This is an O(N²) direct implementation for moderate N.
///
/// Returns the filtered signal.
#[allow(dead_code)]
pub fn spectral_cutoff_filter(signal: &[f64], k_cutoff: usize) -> Vec<f64> {
    let n = signal.len();
    let mut re = vec![0.0_f64; n];
    let mut im = vec![0.0_f64; n];
    for k in 0..n {
        for (j, &s) in signal.iter().enumerate() {
            let angle = 2.0 * PI * k as f64 * j as f64 / n as f64;
            re[k] += s * angle.cos();
            im[k] -= s * angle.sin();
        }
    }
    for k in 0..n {
        let k_eff = k.min(n - k);
        if k_eff > k_cutoff {
            re[k] = 0.0;
            im[k] = 0.0;
        }
    }
    let mut out = vec![0.0_f64; n];
    let n_f64 = n as f64;
    for j in 0..n {
        let mut val = 0.0_f64;
        for k in 0..n {
            let angle = 2.0 * PI * k as f64 * j as f64 / n_f64;
            val += re[k] * angle.cos() - im[k] * angle.sin();
        }
        out[j] = val / n_f64;
    }
    out
}
/// Apply a Gaussian LES filter to a 1D signal in physical space.
///
/// The Gaussian filter width `Delta` is the standard deviation of the kernel:
/// `G(x) = sqrt(6/(pi*Delta^2)) * exp(-6*x^2/Delta^2)`
///
/// Here we implement the discrete convolution version using a truncated kernel
/// of half-width `trunc` grid points.
#[allow(dead_code)]
pub fn gaussian_les_filter(signal: &[f64], delta: f64, dx: f64, trunc: usize) -> Vec<f64> {
    let n = signal.len();
    let sigma2 = delta * delta / 6.0;
    let half = trunc.min(n / 2);
    let mut kernel: Vec<f64> = (0..=half)
        .map(|i| {
            let x = i as f64 * dx;
            (-x * x / (2.0 * sigma2)).exp()
        })
        .collect();
    let norm: f64 = kernel[0] + 2.0 * kernel[1..].iter().sum::<f64>();
    for k in &mut kernel {
        *k /= norm;
    }
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let mut val = kernel[0] * signal[i];
        for (j, &kw) in kernel[1..].iter().enumerate() {
            let j = j + 1;
            let i_p = (i + j) % n;
            let i_m = (i + n - j) % n;
            val += kw * (signal[i_p] + signal[i_m]);
        }
        out[i] = val;
    }
    out
}
/// Kays-Crawford correlation for the turbulent Prandtl number.
///
/// `Pr_t = 1 / (0.5882 + 0.228 * nu_t/nu - 0.0441 * (nu_t/nu)^2 * (1 - exp(-5.165*nu/(nu_t))))`
///
/// Valid for moderate turbulent-to-molecular viscosity ratios.
/// Reduces to the constant `Pr_t ≈ 0.85` for high `nu_t/nu`.
///
/// Returns `Pr_t` (typically 0.7–0.9).
#[allow(dead_code)]
pub fn turbulent_prandtl_kays_crawford(nu_t: f64, nu: f64) -> f64 {
    if nu < 1e-30 {
        return 0.85;
    }
    let r = nu_t / nu;
    if r < 1e-30 {
        return 0.85;
    }
    let exp_term = (-5.165 / r).exp();
    let denom = 0.5882 + 0.228 * r - 0.0441 * r * r * (1.0 - exp_term);
    if denom.abs() < 1e-30 {
        return 0.85;
    }
    1.0 / denom
}
/// Check whether a computed energy spectrum follows the Kolmogorov -5/3 law
/// over a given inertial-subrange wavenumber band `[k_lo, k_hi]`.
///
/// Returns the least-squares slope of `log(E)` vs `log(k)` over the band.
/// A return value close to -5/3 ≈ -1.667 indicates a resolved inertial subrange.
///
/// `spectrum[i]` is the energy at wavenumber `i+1` (1-indexed to avoid log(0)).
#[allow(dead_code)]
pub fn les_spectrum_slope(spectrum: &[f64], k_lo: usize, k_hi: usize) -> f64 {
    if k_lo >= k_hi || k_hi > spectrum.len() {
        return 0.0;
    }
    let mut sum_x = 0.0_f64;
    let mut sum_y = 0.0_f64;
    let mut sum_xx = 0.0_f64;
    let mut sum_xy = 0.0_f64;
    let mut n = 0_usize;
    for k in k_lo..=k_hi {
        if k == 0 || k > spectrum.len() {
            continue;
        }
        let e = spectrum[k - 1];
        if e <= 0.0 {
            continue;
        }
        let lk = (k as f64).ln();
        let le = e.ln();
        sum_x += lk;
        sum_y += le;
        sum_xx += lk * lk;
        sum_xy += lk * le;
        n += 1;
    }
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    let denom = nf * sum_xx - sum_x * sum_x;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    (nf * sum_xy - sum_x * sum_y) / denom
}
/// Verify that the energy spectrum has the Kolmogorov -5/3 slope within a tolerance.
///
/// Returns `true` if `|slope + 5/3| < tolerance`.
#[allow(dead_code)]
pub fn spectrum_is_kolmogorov(spectrum: &[f64], k_lo: usize, k_hi: usize, tolerance: f64) -> bool {
    let slope = les_spectrum_slope(spectrum, k_lo, k_hi);
    (slope + 5.0 / 3.0).abs() < tolerance
}
/// Generate a synthetic Kolmogorov spectrum for testing purposes.
///
/// `E(k) = C_k * epsilon^(2/3) * k^(-5/3)` for `k = 1..n`.
/// Returns a `Vec`f64` of length `n`.
#[allow(dead_code)]
pub fn synthetic_kolmogorov_spectrum(n: usize, epsilon: f64, c_k: f64) -> Vec<f64> {
    (1..=n)
        .map(|k| c_k * epsilon.powf(2.0 / 3.0) * (k as f64).powf(-5.0 / 3.0))
        .collect()
}
/// Compute the Vreman `B_beta` invariant from a velocity-gradient tensor.
///
/// `B_beta = beta_11*beta_22 - beta_12^2 + beta_11*beta_33 - beta_13^2 + beta_22*beta_33 - beta_23^2`
///
/// where `beta_ij = sum_m alpha_mi * alpha_mj` and `alpha_ij = du_i/dx_j`.
///
/// Returns `B_beta` (non-negative for physical flows).
#[allow(dead_code)]
pub fn vreman_b_beta(velocity_gradient: &[f64; 9]) -> f64 {
    let g = |i: usize, j: usize| velocity_gradient[i * 3 + j];
    let mut beta = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for m in 0..3 {
                beta[i][j] += g(m, i) * g(m, j);
            }
        }
    }
    beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] + beta[0][0] * beta[2][2]
        - beta[0][2] * beta[0][2]
        + beta[1][1] * beta[2][2]
        - beta[1][2] * beta[1][2]
}
/// Compute `alpha_ij^2 = sum_ij (du_i/dx_j)^2` from the velocity-gradient tensor.
#[allow(dead_code)]
pub fn vreman_alpha_sq(velocity_gradient: &[f64; 9]) -> f64 {
    velocity_gradient.iter().map(|&v| v * v).sum()
}
/// Vreman eddy-viscosity from pre-computed `B_beta` and `alpha_sq`.
///
/// `nu_t = cv * sqrt(B_beta / alpha_sq)` when `B_beta > 0`.
#[allow(dead_code)]
pub fn vreman_from_invariants(b_beta: f64, alpha_sq: f64, cv: f64) -> f64 {
    if alpha_sq < 1e-30 || b_beta <= 0.0 {
        return 0.0;
    }
    cv * (b_beta / alpha_sq).sqrt()
}
