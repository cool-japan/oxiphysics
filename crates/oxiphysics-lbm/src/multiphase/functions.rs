//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{D3Q19_VELOCITIES, D3Q19_WEIGHTS};

use super::types::{FreeEnergyModel, PsiType, SpinodалDecompositionParams};

#[cfg(test)]
use crate::grid::LbmGrid2D;
#[cfg(test)]
use crate::multiphase::types::*;

/// Compute the Shan-Chen pseudo-potential `psi(rho)`.
///
/// # Arguments
/// - `rho`:      local density
/// - `rho_0`:    reference density (used by Exponential and Sukop-Thorne)
/// - `psi_type`: choice of functional form
pub fn compute_psi(rho: f64, rho_0: f64, psi_type: PsiType) -> f64 {
    match psi_type {
        PsiType::Linear => rho,
        PsiType::Exponential => rho_0 * (1.0 - (-rho / rho_0).exp()),
        PsiType::SukopThorne => (-rho_0 / rho.max(1e-30)).exp(),
        PsiType::YuanSchaefer => rho,
    }
}
/// Compute psi from a custom equation of state.
///
/// `psi = sqrt(2 * (p_eos(rho) - rho * cs^2) / (G * cs^2))`
///
/// Returns zero if the argument under the square root is non-positive.
pub fn compute_psi_from_eos(p_eos: f64, rho: f64, g: f64) -> f64 {
    let cs2 = 1.0 / 3.0;
    let arg = 2.0 * (p_eos - rho * cs2) / (g * cs2);
    if arg > 0.0 { arg.sqrt() } else { 0.0 }
}
/// Carnahan-Starling equation of state.
///
/// `p = rho * R * T * (1 + eta + eta^2 - eta^3) / (1 - eta)^3 - a * rho^2`
///
/// where `eta = b * rho / 4`.
pub fn carnahan_starling_eos(rho: f64, a: f64, b: f64, r_t: f64) -> f64 {
    let eta = b * rho / 4.0;
    let denom = (1.0 - eta).powi(3);
    if denom.abs() < 1e-30 {
        return rho * r_t;
    }
    rho * r_t * (1.0 + eta + eta * eta - eta * eta * eta) / denom - a * rho * rho
}
/// Peng-Robinson equation of state.
///
/// `p = rho * R * T / (1 - b*rho) - a * alpha(T) * rho^2 / (1 + 2*b*rho - b^2*rho^2)`
///
/// Simplified version with alpha = 1.
pub fn peng_robinson_eos(rho: f64, a: f64, b: f64, r_t: f64) -> f64 {
    let denom1 = 1.0 - b * rho;
    let denom2 = 1.0 + 2.0 * b * rho - b * b * rho * rho;
    if denom1.abs() < 1e-30 || denom2.abs() < 1e-30 {
        return rho * r_t;
    }
    rho * r_t / denom1 - a * rho * rho / denom2
}
/// Multi-range Shan-Chen interaction force.
///
/// Uses both nearest-neighbor and next-nearest-neighbor interactions
/// for improved isotropy of the interface:
///
/// `F = -psi_here * (G1 * sum_nn + G2 * sum_nnn)`
///
/// Returns (Fx, Fy).
pub fn multi_range_interaction_force(
    psi_field: &[f64],
    x: usize,
    y: usize,
    nx: usize,
    ny: usize,
    g1: f64,
    g2: f64,
) -> (f64, f64) {
    let idx = |xi: usize, yi: usize| yi * nx + xi;
    let psi_here = psi_field[idx(x, y)];
    let nn_offsets: [(i64, i64, f64); 4] = [
        (1, 0, 1.0 / 3.0),
        (-1, 0, 1.0 / 3.0),
        (0, 1, 1.0 / 3.0),
        (0, -1, 1.0 / 3.0),
    ];
    let diag_offsets: [(i64, i64, f64); 4] = [
        (1, 1, 1.0 / 12.0),
        (-1, 1, 1.0 / 12.0),
        (-1, -1, 1.0 / 12.0),
        (1, -1, 1.0 / 12.0),
    ];
    let nnn_offsets: [(i64, i64, f64); 4] = [
        (2, 0, 1.0 / 12.0),
        (-2, 0, 1.0 / 12.0),
        (0, 2, 1.0 / 12.0),
        (0, -2, 1.0 / 12.0),
    ];
    let mut fx = 0.0;
    let mut fy = 0.0;
    for &(dx, dy, w) in nn_offsets.iter().chain(diag_offsets.iter()) {
        let xi = ((x as i64 + dx).rem_euclid(nx as i64)) as usize;
        let yi = ((y as i64 + dy).rem_euclid(ny as i64)) as usize;
        let psi_nb = psi_field[idx(xi, yi)];
        fx += g1 * w * psi_nb * dx as f64;
        fy += g1 * w * psi_nb * dy as f64;
    }
    for &(dx, dy, w) in &nnn_offsets {
        let xi = ((x as i64 + dx).rem_euclid(nx as i64)) as usize;
        let yi = ((y as i64 + dy).rem_euclid(ny as i64)) as usize;
        let psi_nb = psi_field[idx(xi, yi)];
        fx += g2 * w * psi_nb * dx as f64;
        fy += g2 * w * psi_nb * dy as f64;
    }
    fx *= -psi_here;
    fy *= -psi_here;
    (fx, fy)
}
/// Compute the contact angle wall density for a wetting boundary.
///
/// The wetting condition is enforced by setting a virtual density on
/// solid nodes such that the resulting pseudo-potential gives the desired
/// contact angle theta:
///
/// `rho_wall = rho_vapor + (rho_liquid - rho_vapor) * (1 + cos(theta)) / 2`
///
/// where theta is the contact angle (0 = fully wetting, pi = non-wetting).
pub fn contact_angle_wall_density(rho_liquid: f64, rho_vapor: f64, theta: f64) -> f64 {
    rho_vapor + (rho_liquid - rho_vapor) * 0.5 * (1.0 + theta.cos())
}
/// Apply geometric wetting boundary condition to solid nodes.
///
/// Sets the pseudo-potential on solid boundary nodes to achieve a desired
/// contact angle. `is_solid` marks solid nodes, and the `psi_field` is
/// modified in-place for solid nodes.
pub fn apply_wetting_boundary(psi_field: &mut [f64], is_solid: &[bool], psi_wall: f64) {
    for (k, &solid) in is_solid.iter().enumerate() {
        if solid {
            psi_field[k] = psi_wall;
        }
    }
}
/// Compute the Shan-Chen interaction force at voxel `(x, y, z)` using D3Q19.
///
/// `psi_field` is a flat 3D array of pseudo-potential values indexed
/// `psi_field[z * ny * nx + y * nx + x]`.
///
/// Returns the force vector `[Fx, Fy, Fz]`.
pub fn compute_interaction_force_d3q19(
    psi_field: &[f64],
    x: usize,
    y: usize,
    z: usize,
    nx: usize,
    ny: usize,
    nz: usize,
    g: f64,
) -> [f64; 3] {
    let idx = |xi: usize, yi: usize, zi: usize| zi * ny * nx + yi * nx + xi;
    let psi_here = psi_field[idx(x, y, z)];
    let mut force = [0.0f64; 3];
    for (c, &w) in D3Q19_VELOCITIES[1..19].iter().zip(&D3Q19_WEIGHTS[1..19]) {
        let xi = ((x as i64 + c[0] as i64).rem_euclid(nx as i64)) as usize;
        let yi = ((y as i64 + c[1] as i64).rem_euclid(ny as i64)) as usize;
        let zi = ((z as i64 + c[2] as i64).rem_euclid(nz as i64)) as usize;
        let psi_nb = psi_field[idx(xi, yi, zi)];
        for (force_d, &cd) in force.iter_mut().zip(c.iter()) {
            *force_d += w * psi_nb * cd as f64;
        }
    }
    for force_d in &mut force {
        *force_d *= -g * psi_here;
    }
    force
}
/// Compute the free-energy thermodynamic pressure (standalone function).
///
/// See `FreeEnergyModel::compute_free_energy_pressure`.
pub fn compute_free_energy_pressure(rho: f64, grad_rho_sq: f64, kappa: f64, a: f64, b: f64) -> f64 {
    let model = FreeEnergyModel::new(kappa, a, b);
    model.compute_free_energy_pressure(rho, grad_rho_sq)
}
/// Compute the density ratio between liquid and vapor phases.
///
/// Uses the equilibrium densities from the free-energy model.
pub fn density_ratio(rho_liquid: f64, rho_vapor: f64) -> f64 {
    if rho_vapor.abs() < 1e-30 {
        return f64::INFINITY;
    }
    (rho_liquid / rho_vapor).abs()
}
/// Compute effective density using harmonic mean (for momentum exchange).
///
/// `rho_eff = 2 * rho_1 * rho_2 / (rho_1 + rho_2)`
pub fn harmonic_mean_density(rho_1: f64, rho_2: f64) -> f64 {
    let sum = rho_1 + rho_2;
    if sum.abs() < 1e-30 {
        return 0.0;
    }
    2.0 * rho_1 * rho_2 / sum
}
/// Clip density to physically meaningful range.
///
/// Prevents negative or excessively low densities that can cause
/// numerical instability in multiphase simulations.
pub fn clip_density(rho: f64, rho_min: f64, rho_max: f64) -> f64 {
    rho.clamp(rho_min, rho_max)
}
/// Initialize a flat interface density profile.
///
/// Sets up a 1D density distribution with liquid on one side and vapor
/// on the other, connected by a tanh interface.
///
/// Returns density as a function of position.
pub fn flat_interface_density(
    y: f64,
    y_interface: f64,
    rho_liquid: f64,
    rho_vapor: f64,
    interface_width: f64,
) -> f64 {
    let rho_mean = 0.5 * (rho_liquid + rho_vapor);
    let rho_diff = 0.5 * (rho_liquid - rho_vapor);
    rho_mean + rho_diff * ((y - y_interface) / interface_width).tanh()
}
/// Initialize a circular droplet density profile.
pub fn circular_droplet_density(
    x: f64,
    y: f64,
    cx: f64,
    cy: f64,
    radius: f64,
    rho_liquid: f64,
    rho_vapor: f64,
    interface_width: f64,
) -> f64 {
    let r = ((x - cx).powi(2) + (y - cy).powi(2)).sqrt();
    let rho_mean = 0.5 * (rho_liquid + rho_vapor);
    let rho_diff = 0.5 * (rho_liquid - rho_vapor);
    rho_mean - rho_diff * ((r - radius) / interface_width).tanh()
}
/// Compute the order parameter for phase separation analysis.
///
/// `m = (1/N) * sum_i (rho_i - rho_mean)^2`
///
/// This quantity increases as phase separation progresses.
pub fn phase_separation_order_parameter(rho: &[f64]) -> f64 {
    if rho.is_empty() {
        return 0.0;
    }
    let n = rho.len() as f64;
    let rho_mean = rho.iter().sum::<f64>() / n;
    rho.iter().map(|&r| (r - rho_mean).powi(2)).sum::<f64>() / n
}
/// Compute the spinodal decomposition growth rate.
///
/// For the Cahn-Hilliard equation, the growth rate of a perturbation
/// with wavenumber k is:
///
/// `sigma(k) = -M * k^2 * (f''(phi0) + 2*kappa*k^2)`
///
/// where `f''(phi0) = A + 3*B*phi0^2`.
///
/// Returns the growth rate (positive = unstable).
pub fn spinodal_growth_rate(
    k_wave: f64,
    mobility: f64,
    a: f64,
    b: f64,
    kappa: f64,
    phi0: f64,
) -> f64 {
    let f_double_prime = a + 3.0 * b * phi0 * phi0;
    -mobility * k_wave * k_wave * (f_double_prime + 2.0 * kappa * k_wave * k_wave)
}
/// Check if a state is inside the spinodal region.
///
/// The spinodal is defined by `f''(phi) = A + 3*B*phi^2 < 0`.
pub fn is_spinodal(phi: f64, a: f64, b: f64) -> bool {
    let f_double_prime = a + 3.0 * b * phi * phi;
    f_double_prime < 0.0
}
/// Compute the Maxwell construction (equal area rule) for coexistence densities.
///
/// For a symmetric double-well potential: `phi_eq = sqrt(-A / B)`.
/// Returns `(phi_liquid, phi_vapor)`.
pub fn maxwell_construction(a: f64, b: f64) -> (f64, f64) {
    if a < 0.0 && b > 0.0 {
        let phi_eq = (-a / b).sqrt();
        (phi_eq, -phi_eq)
    } else {
        (0.0, 0.0)
    }
}
/// Compute the interface profile from the tanh solution.
///
/// `phi(x) = phi_eq * tanh(x / (sqrt(2) * xi))`
///
/// where `xi = sqrt(-kappa / A)` is the interface width.
pub fn tanh_interface_profile(x: f64, a: f64, b: f64, kappa: f64) -> f64 {
    if a >= 0.0 || b <= 0.0 || kappa <= 0.0 {
        return 0.0;
    }
    let phi_eq = (-a / b).sqrt();
    let xi = (-kappa / a).sqrt();
    phi_eq * (x / (std::f64::consts::SQRT_2 * xi)).tanh()
}
/// Determine whether a cell with density `rho` is an interface cell.
///
/// A cell is considered an interface cell if its density lies strictly between
/// the liquid and vapor equilibrium densities.
pub fn is_interface_cell(rho: f64) -> bool {
    let rho_low = 0.4;
    let rho_high = 1.6;
    rho > rho_low && rho < rho_high
}
/// Determine interface cells with custom thresholds.
pub fn is_interface_cell_custom(rho: f64, rho_low: f64, rho_high: f64) -> bool {
    rho > rho_low && rho < rho_high
}
/// Count the number of cells in each phase.
///
/// Returns `(n_liquid, n_interface, n_vapor)`.
pub fn count_phases(rho: &[f64], rho_low: f64, rho_high: f64) -> (usize, usize, usize) {
    let mut n_liquid = 0;
    let mut n_interface = 0;
    let mut n_vapor = 0;
    for &r in rho {
        if r >= rho_high {
            n_liquid += 1;
        } else if r <= rho_low {
            n_vapor += 1;
        } else {
            n_interface += 1;
        }
    }
    (n_liquid, n_interface, n_vapor)
}
/// Compute the Laplacian of a 2D scalar field using 5-point stencil with periodic BCs.
///
/// Returns laplacian\[y * nx + x\].
pub fn laplacian_2d_periodic(field: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let mut lap = vec![0.0f64; nx * ny];
    let idx = |xi: usize, yi: usize| yi * nx + xi;
    for y in 0..ny {
        for x in 0..nx {
            let xp = (x + 1) % nx;
            let xm = (x + nx - 1) % nx;
            let yp = (y + 1) % ny;
            let ym = (y + ny - 1) % ny;
            let c = field[idx(x, y)];
            lap[idx(x, y)] =
                field[idx(xp, y)] + field[idx(xm, y)] + field[idx(x, yp)] + field[idx(x, ym)]
                    - 4.0 * c;
        }
    }
    lap
}
/// Perform one forward-Euler Cahn-Hilliard step on a 2D order-parameter field `phi`.
///
/// `phi_{n+1} = phi_n + dt * M * lap( A*phi + B*phi^3 - kappa * lap(phi) )`
///
/// Modifies `phi` in-place.
pub fn cahn_hilliard_step(
    phi: &mut [f64],
    nx: usize,
    ny: usize,
    dt: f64,
    params: &SpinodалDecompositionParams,
) {
    let lap_phi = laplacian_2d_periodic(phi, nx, ny);
    let mut mu = vec![0.0f64; nx * ny];
    for (k, mu_k) in mu.iter_mut().enumerate() {
        *mu_k = params.a * phi[k] + params.b * phi[k].powi(3) - params.kappa * lap_phi[k];
    }
    let lap_mu = laplacian_2d_periodic(&mu, nx, ny);
    for (phi_k, &lap_mu_k) in phi.iter_mut().zip(&lap_mu) {
        *phi_k += dt * params.mobility * lap_mu_k;
    }
}
/// Compute the total free energy for a 2D phi field.
///
/// `F = sum [ A/2 * phi^2 + B/4 * phi^4 + kappa/2 * |grad phi|^2 ]`
pub fn total_free_energy_2d(phi: &[f64], nx: usize, ny: usize, a: f64, b: f64, kappa: f64) -> f64 {
    let idx = |xi: usize, yi: usize| yi * nx + xi;
    let mut energy = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            let xp = (x + 1) % nx;
            let yp = (y + 1) % ny;
            let phi_c = phi[idx(x, y)];
            let dphi_dx = phi[idx(xp, y)] - phi_c;
            let dphi_dy = phi[idx(x, yp)] - phi_c;
            let grad_sq = dphi_dx * dphi_dx + dphi_dy * dphi_dy;
            energy += 0.5 * a * phi_c * phi_c + 0.25 * b * phi_c.powi(4) + 0.5 * kappa * grad_sq;
        }
    }
    energy
}
/// Compute the Laplace pressure across a 2D droplet.
///
/// Using the Young-Laplace equation: `ΔP = 2 * sigma / R` for a 2D circular droplet
/// (factors of 2 apply in 2D: single interface in each direction).
///
/// Returns `2 * sigma / radius`.
pub fn laplace_pressure_2d(sigma: f64, radius: f64) -> f64 {
    if radius < 1e-30 {
        f64::INFINITY
    } else {
        2.0 * sigma / radius
    }
}
/// Compute the Laplace pressure across a 3D spherical droplet.
///
/// `ΔP = 4 * sigma / D` (two interfaces for a bubble or `2*sigma/R` for a droplet).
pub fn laplace_pressure_3d_droplet(sigma: f64, radius: f64) -> f64 {
    if radius < 1e-30 {
        f64::INFINITY
    } else {
        2.0 * sigma / radius
    }
}
/// Rayleigh–Plesset simplified: estimate droplet oscillation frequency.
///
/// For small oscillations of a bubble of radius R0:
/// `omega_n = sqrt(3 * kappa_ad * P_inf / (rho_liq * R0^2))`
///
/// where `kappa_ad` is the polytropic index and `P_inf` the ambient pressure.
pub fn bubble_natural_frequency(r0: f64, rho_liquid: f64, p_inf: f64, kappa_ad: f64) -> f64 {
    if r0 < 1e-30 || rho_liquid < 1e-30 {
        return 0.0;
    }
    (3.0 * kappa_ad * p_inf / (rho_liquid * r0 * r0)).sqrt()
}
/// Estimate droplet terminal velocity using Hadamard–Rybczynski model.
///
/// For a spherical droplet of radius R rising in a fluid of viscosity mu_c:
/// `U_t = (2/3) * R^2 * (rho_c - rho_d) * g / (mu_c * (mu_d + mu_c) / (mu_d + 1.5*mu_c))`
pub fn droplet_terminal_velocity(
    radius: f64,
    rho_continuous: f64,
    rho_droplet: f64,
    mu_continuous: f64,
    mu_droplet: f64,
    gravity: f64,
) -> f64 {
    let delta_rho = rho_continuous - rho_droplet;
    let factor = (mu_droplet + mu_continuous) / (mu_droplet + 1.5 * mu_continuous);
    let denom = mu_continuous;
    if denom < 1e-30 {
        return 0.0;
    }
    (2.0 / 3.0) * radius * radius * delta_rho * gravity / denom * factor
}
/// Spherical droplet drag coefficient (Hadamard–Rybczynski).
///
/// `Cd = (8/Re) * (2 + 3 * lambda) / (3 * (1 + lambda))`
/// where `lambda = mu_droplet / mu_continuous`.
pub fn droplet_drag_coefficient(re: f64, mu_droplet: f64, mu_continuous: f64) -> f64 {
    if re < 1e-30 {
        return f64::INFINITY;
    }
    let lambda = mu_droplet / mu_continuous;
    (8.0 / re) * (2.0 + 3.0 * lambda) / (3.0 * (1.0 + lambda))
}
/// Leidenfrost number estimation: superheat for film boiling onset.
///
/// Approximation using the Berenson correlation for horizontal surfaces:
/// `T_Leidenfrost ≈ T_sat + C_L * (sigma^2 / (rho_v^2 * cp_v * lambda))^(1/3)`
///
/// Here we use a simplified non-dimensional form.
/// Returns the estimated superheat temperature difference `ΔT`.
pub fn leidenfrost_superheat(
    sigma: f64,
    rho_vapor: f64,
    cp_vapor: f64,
    thermal_conductivity: f64,
) -> f64 {
    if rho_vapor < 1e-30 || cp_vapor < 1e-30 {
        return 0.0;
    }
    (sigma / rho_vapor).sqrt() / cp_vapor * thermal_conductivity
}
/// Film boiling heat transfer coefficient (Bromley correlation simplified).
///
/// `h_film = C * [ k_v^3 * rho_v * (rho_l - rho_v) * g * h_fg ] / [ mu_v * delta_T * D ]^0.25`
///
/// This is the functional form; we expose a simplified version:
/// `h ~ [k^3 * rho^2 * g * h_fg / (mu * delta_T)]^0.25`
pub fn film_boiling_heat_transfer(
    k_vapor: f64,
    rho_vapor: f64,
    rho_liquid: f64,
    mu_vapor: f64,
    g: f64,
    h_fg: f64,
    delta_t: f64,
    char_length: f64,
) -> f64 {
    if mu_vapor < 1e-30 || delta_t < 1e-30 || char_length < 1e-30 {
        return 0.0;
    }
    let num = k_vapor.powi(3) * rho_vapor * (rho_liquid - rho_vapor).abs() * g * h_fg;
    let den = mu_vapor * delta_t * char_length;
    0.62 * (num / den).powf(0.25)
}
/// Compute the local phase change source term for a Lee-Liu-like model.
///
/// The mass source is driven by the temperature departure from saturation:
/// `S_mass = Gamma * (T - T_sat)`
///
/// where `Gamma` is the volumetric phase change rate coefficient.
pub fn phase_change_source(temperature: f64, t_sat: f64, gamma: f64) -> f64 {
    gamma * (temperature - t_sat)
}
/// Compute the energy source term due to latent heat release/absorption.
///
/// `Q_latent = L_fg * S_mass`
///
/// where `L_fg` is the latent heat of vaporization.
pub fn latent_heat_source(s_mass: f64, l_fg: f64) -> f64 {
    l_fg * s_mass
}
/// Apply phase change to a density field (condensation/evaporation).
///
/// Modifies the density field in-place by adding `dt * S_mass` at each node.
/// The rate `S_mass` is computed from the temperature field.
///
/// This implements the simplified pseudopotential phase change:
/// phase-change occurs when `rho * (T - T_sat) * gamma * dt` is non-negligible.
pub fn apply_phase_change(rho: &mut [f64], temperature: &[f64], t_sat: f64, gamma: f64, dt: f64) {
    for (r, &t) in rho.iter_mut().zip(temperature.iter()) {
        let s = phase_change_source(t, t_sat, gamma);
        *r += dt * s;
        if *r < 0.0 {
            *r = 0.0;
        }
    }
}
/// Compute the evaporation rate from a liquid-vapor interface.
///
/// Using the Hertz-Knudsen-Schrage correlation:
/// `J_evap = alpha_e * sqrt(M / (2 * pi * R * T_sat)) * (P_sat(T) - P_v)`
///
/// Simplified version with alpha_e = 1.
pub fn hertz_knudsen_evaporation_rate(
    t_surface: f64,
    t_sat: f64,
    p_vapor_partial: f64,
    molar_mass: f64,
    r_gas: f64,
) -> f64 {
    let p_sat_ref = 1.0;
    let p_sat = p_sat_ref * (molar_mass * (t_surface - t_sat) / (r_gas * t_sat * t_sat)).exp();
    let dp = p_sat - p_vapor_partial;
    if r_gas < 1e-30 || t_sat < 1e-30 {
        return 0.0;
    }
    (molar_mass / (std::f64::consts::TAU * r_gas * t_sat)).sqrt() * dp
}
/// Phase indicator (volume fraction): smoothly transitions from 0 (vapor) to 1 (liquid).
///
/// `C = 0.5 * (1 + tanh((rho - rho_mid) / delta))`
pub fn volume_fraction(rho: f64, rho_liquid: f64, rho_vapor: f64, interface_width: f64) -> f64 {
    let rho_mid = 0.5 * (rho_liquid + rho_vapor);
    let delta = interface_width.max(1e-30);
    0.5 * (1.0 + ((rho - rho_mid) / delta).tanh())
}
/// Compute the interface normal vector at a 2D node using central differences.
///
/// Returns the unit normal `[nx, ny]` pointing from vapor to liquid.
pub fn interface_normal_2d(
    phi: &[f64],
    x: usize,
    y: usize,
    grid_nx: usize,
    grid_ny: usize,
) -> [f64; 2] {
    let idx = |xi: usize, yi: usize| yi * grid_nx + xi;
    let xp = (x + 1) % grid_nx;
    let xm = (x + grid_nx - 1) % grid_nx;
    let yp = (y + 1) % grid_ny;
    let ym = (y + grid_ny - 1) % grid_ny;
    let gx = 0.5 * (phi[idx(xp, y)] - phi[idx(xm, y)]);
    let gy = 0.5 * (phi[idx(x, yp)] - phi[idx(x, ym)]);
    let mag = (gx * gx + gy * gy).sqrt();
    if mag > 1e-30 {
        [gx / mag, gy / mag]
    } else {
        [0.0, 0.0]
    }
}
/// Compute the interface curvature at a 2D node from the phase field gradient.
///
/// Uses the divergence of the unit normal: `kappa = div(n)`.
/// Returns a signed curvature (positive = concave toward vapor).
pub fn interface_curvature_2d(
    phi: &[f64],
    x: usize,
    y: usize,
    grid_nx: usize,
    grid_ny: usize,
) -> f64 {
    let xp = (x + 1) % grid_nx;
    let xm = (x + grid_nx - 1) % grid_nx;
    let yp = (y + 1) % grid_ny;
    let ym = (y + grid_ny - 1) % grid_ny;
    let n_xp = interface_normal_2d(phi, xp, y, grid_nx, grid_ny);
    let n_xm = interface_normal_2d(phi, xm, y, grid_nx, grid_ny);
    let n_yp = interface_normal_2d(phi, x, yp, grid_nx, grid_ny);
    let n_ym = interface_normal_2d(phi, x, ym, grid_nx, grid_ny);
    0.5 * (n_xp[0] - n_xm[0]) + 0.5 * (n_yp[1] - n_ym[1])
}
/// Compute the Sukop-Thorne pseudo-potential with explicit `psi0`:
///
/// `psi(rho) = psi0 * exp(-rho0 / rho)`
///
/// This form allows independent control of the saturation value `psi0`
/// and the density scale `rho0`.  Setting `psi0 = 1` recovers the
/// original Shan-Chen form used elsewhere in this module.
pub fn sukop_thorne_psi(rho: f64, rho0: f64, psi0: f64) -> f64 {
    psi0 * (-rho0 / rho.max(1e-30)).exp()
}
/// Estimate the critical coupling `G_c` for the Sukop-Thorne model.
///
/// Phase separation occurs when `|G| > G_c`.  For the Sukop-Thorne
/// pseudo-potential the spinodal condition gives (in lattice units):
///
/// `G_c = -1 / (6 * psi0^2 * exp(-2))`
///
/// where `exp(-2)` comes from evaluating `d(psi^2)/d(rho)` at the
/// inflection point `rho = rho0 / 2`.
pub fn sukop_thorne_critical_g(psi0: f64) -> f64 {
    -1.0 / (6.0 * psi0 * psi0 * (-2.0_f64).exp())
}
/// Detect if two circular droplets are coalescing.
///
/// Coalescence is triggered when the separation between droplet centres
/// drops below `r1 + r2 + tolerance` (centres closer than the sum of radii).
///
/// Returns `true` if the droplets overlap or are within `tolerance` of contact.
pub fn droplets_coalescing(
    cx1: f64,
    cy1: f64,
    r1: f64,
    cx2: f64,
    cy2: f64,
    r2: f64,
    tolerance: f64,
) -> bool {
    let dist = ((cx2 - cx1).powi(2) + (cy2 - cy1).powi(2)).sqrt();
    dist < r1 + r2 + tolerance
}
/// Estimate the coalescence time scale for two approaching droplets.
///
/// Uses the lubrication approximation result for film drainage:
///
/// `t_c ~ (3 * mu * R_eff) / (2 * sigma * h_0)`
///
/// where `h_0` is the initial gap, `R_eff = 2*R1*R2/(R1+R2)` and
/// `sigma` is the surface tension.
pub fn coalescence_time_scale(r1: f64, r2: f64, mu: f64, sigma: f64, h0: f64) -> f64 {
    if sigma < 1e-30 || h0 < 1e-30 {
        return f64::INFINITY;
    }
    let r_eff = 2.0 * r1 * r2 / (r1 + r2).max(1e-30);
    3.0 * mu * r_eff / (2.0 * sigma * h0)
}
/// Estimate the new droplet radius after coalescence (volume conservation).
///
/// For 2D: `r_new = sqrt(r1^2 + r2^2)`
/// For 3D: `r_new = (r1^3 + r2^3)^(1/3)`
pub fn coalesced_radius_2d(r1: f64, r2: f64) -> f64 {
    (r1 * r1 + r2 * r2).sqrt()
}
/// 3D coalesced radius.
pub fn coalesced_radius_3d(r1: f64, r2: f64) -> f64 {
    (r1.powi(3) + r2.powi(3)).powf(1.0 / 3.0)
}
/// Phase separation rate constant based on Cahn-Hilliard linear stability.
///
/// The rate at which phase separation grows from a small perturbation
/// of wavenumber `k` is:
///
/// `rate(k) = -M * k^2 * (f''(phi_0) + 2 * kappa * k^2)`
///
/// For the most unstable wavenumber `k* = sqrt(-f''/(2*kappa))`, the
/// maximum growth rate is:
///
/// `rate_max = M * f''^2 / (4 * kappa)`
///
/// This function returns the dimensional rate constant given the
/// mobility `M`, the curvature of the free energy `f_double_prime`, and
/// the gradient coefficient `kappa`.
///
/// A **positive** return value means the uniform phase is unstable.
pub fn phase_separation_rate_constant(mobility: f64, f_double_prime: f64, kappa: f64) -> f64 {
    if kappa < 1e-30 || f_double_prime >= 0.0 {
        return 0.0;
    }
    mobility * f_double_prime * f_double_prime / (4.0 * kappa)
}
/// Dimensionless spinodal decomposition number.
///
/// `Sp = |f''(phi_0)| / (kappa * k_min^2)`
///
/// A value > 1 indicates that the system is in the spinodal region and
/// will phase-separate spontaneously.
pub fn spinodal_number(f_double_prime: f64, kappa: f64, k_min: f64) -> f64 {
    if kappa < 1e-30 || k_min < 1e-30 {
        return 0.0;
    }
    f_double_prime.abs() / (kappa * k_min * k_min)
}
/// Extract surface tension from a Young-Laplace test simulation.
///
/// Given a 2D circular droplet in mechanical equilibrium, the Young-Laplace
/// equation gives:
///
/// `sigma = (P_inside - P_outside) * R / (D - 1)`
///
/// For 2D (D=2): `sigma = ΔP * R`
/// For 3D (D=3): `sigma = ΔP * R / 2`
///
/// This function takes the measured interior and exterior pressures and the
/// droplet radius, and returns the surface tension.
///
/// `dimension` must be 2 or 3.
pub fn surface_tension_from_laplace(
    p_inside: f64,
    p_outside: f64,
    radius: f64,
    dimension: usize,
) -> f64 {
    let delta_p = p_inside - p_outside;
    match dimension {
        2 => delta_p * radius,
        3 => delta_p * radius / 2.0,
        _ => 0.0,
    }
}
/// Estimate droplet equilibrium radius from Laplace pressure balance.
///
/// Given surface tension `sigma`, external pressure `p_ext`, and gas
/// pressure `p_gas`, the equilibrium radius satisfies:
///
/// `R_eq = 2 * sigma / (p_gas - p_ext)` (3D spherical)
pub fn equilibrium_droplet_radius(sigma: f64, p_gas: f64, p_ext: f64) -> f64 {
    let dp = p_gas - p_ext;
    if dp < 1e-30 {
        return f64::INFINITY;
    }
    2.0 * sigma / dp
}
/// Compute the diffuse interface width from Cahn-Hilliard parameters.
///
/// The tanh interface profile gives a width (10–90% transition) of:
///
/// `W_90 = 2 * sqrt(2) * xi` where `xi = sqrt(kappa / |A|)`.
///
/// This width is the relevant parameter for mesh resolution requirements:
/// at least 4–5 grid points should span `W_90`.
pub fn diffuse_interface_width_90(kappa: f64, a: f64) -> f64 {
    if a >= 0.0 || kappa <= 0.0 {
        return 0.0;
    }
    2.0 * std::f64::consts::SQRT_2 * (kappa / (-a)).sqrt()
}
/// Check whether the grid spacing `dx` resolves the diffuse interface.
///
/// Returns `true` if at least `n_min` grid points span the 10–90% width.
pub fn interface_is_resolved(kappa: f64, a: f64, dx: f64, n_min: usize) -> bool {
    let w90 = diffuse_interface_width_90(kappa, a);
    w90 / dx >= n_min as f64
}
/// Locate the tie-line endpoints from the common-tangent construction.
///
/// For a symmetric double-well `f(phi) = A/2 * phi^2 + B/4 * phi^4`,
/// the coexistence densities satisfy:
///
/// - `df/dphi = 0` → `phi_coex = ±sqrt(-A/B)` (bulk phases)
/// - Common-tangent: `mu_1 = mu_2` and `f_1 - mu*phi_1 = f_2 - mu*phi_2`
///
/// For the symmetric potential the common tangent passes through `mu = 0`
/// and touches both minima, giving `phi_L = +sqrt(-A/B)`, `phi_V = -sqrt(-A/B)`.
///
/// Returns `(phi_vapor, phi_liquid)` sorted with vapor < liquid.
pub fn tie_line_endpoints(a: f64, b: f64) -> (f64, f64) {
    if a < 0.0 && b > 0.0 {
        let phi_eq = (-a / b).sqrt();
        (-phi_eq, phi_eq)
    } else {
        (0.0, 0.0)
    }
}
/// Compute the lever rule volume fraction of the liquid phase.
///
/// Given the overall composition `phi_0` and the tie-line endpoints
/// `(phi_v, phi_l)`, the liquid volume fraction is:
///
/// `x_l = (phi_0 - phi_v) / (phi_l - phi_v)`
pub fn lever_rule_liquid_fraction(phi_0: f64, phi_v: f64, phi_l: f64) -> f64 {
    let denom = phi_l - phi_v;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    ((phi_0 - phi_v) / denom).clamp(0.0, 1.0)
}
/// Add a small Gaussian perturbation to a scalar field at a given centre.
///
/// `field[k] += amplitude * exp(-((x - cx)^2 + (y - cy)^2) / (2 * sigma^2))`
///
/// Used to seed nucleation events in phase-separation simulations.
pub fn add_nucleation_seed(
    field: &mut [f64],
    nx: usize,
    ny: usize,
    cx: f64,
    cy: f64,
    amplitude: f64,
    sigma: f64,
) {
    let sigma2 = (sigma * sigma).max(1e-30);
    for y in 0..ny {
        for x in 0..nx {
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let r2 = dx * dx + dy * dy;
            field[y * nx + x] += amplitude * (-r2 / (2.0 * sigma2)).exp();
        }
    }
}
/// Add multiple random nucleation seeds to a field.
///
/// Seeds are placed at positions `(x_k, y_k)` drawn from the `centres` slice.
/// Each seed has the same amplitude and width.
pub fn add_nucleation_seeds(
    field: &mut [f64],
    nx: usize,
    ny: usize,
    centres: &[(f64, f64)],
    amplitude: f64,
    sigma: f64,
) {
    for &(cx, cy) in centres {
        add_nucleation_seed(field, nx, ny, cx, cy, amplitude, sigma);
    }
}
/// Wetting boundary condition via a modified pseudo-potential on solid nodes.
///
/// This extends `apply_wetting_boundary` to compute `psi_wall` automatically
/// from the desired contact angle `theta`, given the bulk liquid/vapor
/// pseudo-potential values `psi_liq` and `psi_vap`.
///
/// `psi_wall = psi_vap + (psi_liq - psi_vap) * (1 + cos(theta)) / 2`
///
/// This sets a virtual psi on the solid that interpolates between the two
/// bulk phases based on the wettability.
pub fn wetting_psi_from_angle(psi_liq: f64, psi_vap: f64, theta: f64) -> f64 {
    psi_vap + (psi_liq - psi_vap) * 0.5 * (1.0 + theta.cos())
}
/// Apply a contact-angle wetting boundary using the psi field.
///
/// For each solid node, sets `psi_field[k]` to `psi_wall`.
/// Non-solid nodes are left unchanged.
pub fn apply_contact_angle_bc(
    psi_field: &mut [f64],
    is_solid: &[bool],
    psi_liq: f64,
    psi_vap: f64,
    theta: f64,
) {
    let psi_wall = wetting_psi_from_angle(psi_liq, psi_vap, theta);
    for (k, &solid) in is_solid.iter().enumerate() {
        if solid {
            psi_field[k] = psi_wall;
        }
    }
}
/// Compute the apparent contact angle from the psi values at a wall.
///
/// Inverts the wetting formula:
/// `theta = acos(2 * (psi_wall - psi_vap) / (psi_liq - psi_vap) - 1)`
pub fn contact_angle_from_psi(psi_wall: f64, psi_liq: f64, psi_vap: f64) -> f64 {
    let denom = psi_liq - psi_vap;
    if denom.abs() < 1e-30 {
        return std::f64::consts::FRAC_PI_2;
    }
    let cos_theta = 2.0 * (psi_wall - psi_vap) / denom - 1.0;
    cos_theta.clamp(-1.0, 1.0).acos()
}
/// Compute the coupling constant `G_cc` required for a target density ratio.
///
/// For the symmetric Shan-Chen model in D2Q9, the density ratio is
/// approximately:
///
/// `D = exp(2 * |G_cc| * sum_i w_i * c_i^2 * psi_0^2)`
///
/// where the sum over D2Q9 gives `sum = 1` (nearest neighbours).
/// Solving for `G_cc`:
///
/// `G_cc = -ln(D) / (2 * psi_0^2)`
///
/// This estimate holds for moderate density ratios (up to ~10).
pub fn g_cc_for_density_ratio(density_ratio_target: f64, psi_0: f64) -> f64 {
    if density_ratio_target <= 1.0 + 1e-10 || psi_0 < 1e-30 {
        return 0.0;
    }
    -density_ratio_target.ln() / (2.0 * psi_0 * psi_0)
}
/// Validate that the chosen `G_cc` and `psi_0` are below the stability limit.
///
/// For D2Q9, the BGK scheme is stable when
/// `|G_cc| * psi_0^2 < 1/6` (heuristic upper bound).
pub fn g_cc_is_stable(g_cc: f64, psi_0: f64) -> bool {
    g_cc.abs() * psi_0 * psi_0 < 1.0 / 6.0
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::LatticeType;
    #[test]
    fn test_shan_chen_uniform_zero_force() {
        let nx = 8;
        let ny = 8;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        grid.compute_macroscopic();
        let sc = ShanChenModel::new(-1.0, 1.0);
        let (fx, fy) = sc.interaction_force(&grid, 4, 4);
        assert!(
            fx.abs() < 1e-14 && fy.abs() < 1e-14,
            "Uniform density: force should be zero, got ({fx}, {fy})"
        );
    }
    #[test]
    fn test_shan_chen_gradient_nonzero_force() {
        let nx = 10;
        let ny = 10;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        for y in 0..ny {
            for x in 0..nx {
                let rho = 1.0 + 0.1 * x as f64;
                grid.set_equilibrium(x, y, rho, 0.0, 0.0);
            }
        }
        grid.compute_macroscopic();
        let sc = ShanChenModel::new(-1.0, 1.0);
        let (fx, _) = sc.interaction_force(&grid, 5, 5);
        assert!(
            fx.abs() > 1e-10,
            "Force should be non-zero for density gradient"
        );
    }
    #[test]
    fn test_psi_monotonicity() {
        let rho_0 = 1.0;
        let rho_vals: [f64; 5] = [0.1, 0.5, 1.0, 2.0, 5.0];
        let psi_vals: Vec<f64> = rho_vals
            .iter()
            .map(|&r| compute_psi(r, rho_0, PsiType::Exponential))
            .collect();
        for i in 1..psi_vals.len() {
            assert!(
                psi_vals[i] > psi_vals[i - 1],
                "psi not monotone: psi[{i}]={} <= psi[{}]={}",
                psi_vals[i],
                i - 1,
                psi_vals[i - 1]
            );
        }
    }
    #[test]
    fn test_psi_linear() {
        for rho in [0.1, 0.5, 1.0, 2.5] {
            let psi = compute_psi(rho, 1.0, PsiType::Linear);
            assert!(
                (psi - rho).abs() < 1e-14,
                "Linear psi != rho: psi={psi}, rho={rho}"
            );
        }
    }
    #[test]
    fn test_free_energy_pressure_gradient_correction() {
        let model = FreeEnergyModel::new(0.1, -0.1, 0.1);
        let rho = 0.5;
        let p_bulk = model.compute_free_energy_pressure(rho, 0.0);
        let p_grad = model.compute_free_energy_pressure(rho, 1.0);
        assert!(
            p_bulk > p_grad,
            "Gradient correction should reduce pressure: bulk={p_bulk}, grad={p_grad}"
        );
    }
    #[test]
    fn test_free_energy_surface_tension_positive() {
        let model = FreeEnergyModel::new(0.01, -0.1, 0.1);
        let sigma = model.surface_tension();
        assert!(sigma > 0.0, "Surface tension should be positive: {sigma}");
    }
    #[test]
    fn test_is_interface_cell() {
        assert!(!is_interface_cell(2.0), "rho=2.0 should be pure liquid");
        assert!(!is_interface_cell(0.05), "rho=0.05 should be pure vapor");
        assert!(is_interface_cell(1.0), "rho=1.0 should be interface");
        assert!(is_interface_cell(0.5), "rho=0.5 should be interface");
    }
    #[test]
    fn test_3d_force_uniform_zero() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let psi_field = vec![1.0f64; nx * ny * nz];
        let force = compute_interaction_force_d3q19(&psi_field, 2, 2, 2, nx, ny, nz, -1.0);
        for (d, &force_d) in force.iter().enumerate() {
            assert!(
                force_d.abs() < 1e-14,
                "3D uniform force should be zero: F[{d}] = {}",
                force_d
            );
        }
    }
    #[test]
    fn test_free_energy_equilibrium_densities() {
        let a = -0.2;
        let b = 0.1;
        let model = FreeEnergyModel::new(0.01, a, b);
        let rho_liq = model.rho_liquid();
        let expected_sq = -a / b;
        assert!(
            (rho_liq * rho_liq - expected_sq).abs() < 1e-14,
            "rho_liquid^2 should equal -A/B: got {}, expected {}",
            rho_liq * rho_liq,
            expected_sq
        );
    }
    #[test]
    fn test_carnahan_starling_ideal_limit() {
        let rho = 1.0;
        let r_t = 1.0;
        let p = carnahan_starling_eos(rho, 0.0, 1e-10, r_t);
        assert!(
            (p - rho * r_t).abs() < 0.01,
            "CS EOS should approach ideal gas: p={p}"
        );
    }
    #[test]
    fn test_peng_robinson_ideal_limit() {
        let rho = 0.5;
        let r_t = 1.0;
        let p = peng_robinson_eos(rho, 0.0, 0.0, r_t);
        assert!(
            (p - rho * r_t).abs() < 1e-10,
            "PR EOS should be ideal gas: p={p}"
        );
    }
    #[test]
    fn test_contact_angle_wall_density() {
        let rho_l = 2.0;
        let rho_v = 0.1;
        let rho_wall_0 = contact_angle_wall_density(rho_l, rho_v, 0.0);
        assert!(
            (rho_wall_0 - rho_l).abs() < 1e-14,
            "theta=0: rho_wall should be rho_l: {rho_wall_0}"
        );
        let rho_wall_pi = contact_angle_wall_density(rho_l, rho_v, std::f64::consts::PI);
        assert!(
            (rho_wall_pi - rho_v).abs() < 1e-14,
            "theta=pi: rho_wall should be rho_v: {rho_wall_pi}"
        );
        let rho_wall_90 = contact_angle_wall_density(rho_l, rho_v, std::f64::consts::FRAC_PI_2);
        let expected = 0.5 * (rho_l + rho_v);
        assert!(
            (rho_wall_90 - expected).abs() < 1e-12,
            "theta=90: {rho_wall_90} vs {expected}"
        );
    }
    #[test]
    fn test_density_ratio() {
        let ratio = density_ratio(2.0, 0.1);
        assert!((ratio - 20.0).abs() < 1e-14, "ratio should be 20: {ratio}");
    }
    #[test]
    fn test_harmonic_mean_density() {
        let rho_eff = harmonic_mean_density(1.0, 1.0);
        assert!(
            (rho_eff - 1.0).abs() < 1e-14,
            "Equal densities: harmonic mean = same"
        );
        let rho_eff2 = harmonic_mean_density(2.0, 0.0);
        assert!(
            rho_eff2.abs() < 1e-14,
            "One zero density: harmonic mean = 0"
        );
    }
    #[test]
    fn test_order_parameter_uniform() {
        let rho = vec![1.0; 100];
        let op = phase_separation_order_parameter(&rho);
        assert!(
            op.abs() < 1e-14,
            "Uniform density → zero order parameter: {op}"
        );
    }
    #[test]
    fn test_order_parameter_separated() {
        let mut rho = vec![0.1; 50];
        rho.extend(vec![2.0; 50]);
        let op = phase_separation_order_parameter(&rho);
        assert!(
            op > 0.0,
            "Separated phases → positive order parameter: {op}"
        );
    }
    #[test]
    fn test_spinodal_region() {
        let a = -0.1;
        let b = 0.1;
        assert!(
            is_spinodal(0.0, a, b),
            "phi=0 should be in spinodal for A<0"
        );
        let phi_eq = (-a / b).sqrt();
        assert!(
            !is_spinodal(2.0 * phi_eq, a, b),
            "phi>>phi_eq should be outside spinodal"
        );
    }
    #[test]
    fn test_maxwell_construction() {
        let a = -0.2;
        let b = 0.1;
        let (phi_l, phi_v) = maxwell_construction(a, b);
        assert!(
            phi_l > 0.0 && phi_v < 0.0,
            "phi_l should be positive, phi_v negative"
        );
        assert!(
            (phi_l + phi_v).abs() < 1e-14,
            "Symmetric: phi_l + phi_v = 0"
        );
        assert!((phi_l * phi_l - (-a / b)).abs() < 1e-14, "phi_l^2 = -A/B");
    }
    #[test]
    fn test_flat_interface_density_limits() {
        let rho_l = 2.0;
        let rho_v = 0.1;
        let w = 2.0;
        let y_int = 50.0;
        let rho_above = flat_interface_density(100.0, y_int, rho_l, rho_v, w);
        assert!(
            (rho_above - rho_l).abs() < 0.01,
            "Far above: {rho_above} should be ~{rho_l}"
        );
        let rho_below = flat_interface_density(0.0, y_int, rho_l, rho_v, w);
        assert!(
            (rho_below - rho_v).abs() < 0.01,
            "Far below: {rho_below} should be ~{rho_v}"
        );
    }
    #[test]
    fn test_circular_droplet_density() {
        let rho_l = 2.0;
        let rho_v = 0.1;
        let rho_center = circular_droplet_density(50.0, 50.0, 50.0, 50.0, 10.0, rho_l, rho_v, 2.0);
        assert!(
            (rho_center - rho_l).abs() < 0.1,
            "Center should be ~liquid: {rho_center}"
        );
        let rho_far = circular_droplet_density(0.0, 0.0, 50.0, 50.0, 10.0, rho_l, rho_v, 2.0);
        assert!(
            (rho_far - rho_v).abs() < 0.1,
            "Far away should be ~vapor: {rho_far}"
        );
    }
    #[test]
    fn test_tanh_interface_profile() {
        let a = -0.1;
        let b = 0.1;
        let kappa = 0.01;
        let phi_0 = tanh_interface_profile(0.0, a, b, kappa);
        assert!(phi_0.abs() < 1e-14, "At x=0, tanh(0)=0: {phi_0}");
        let phi_far = tanh_interface_profile(100.0, a, b, kappa);
        let phi_eq = (-a / b).sqrt();
        assert!(
            (phi_far - phi_eq).abs() < 0.01,
            "Far from interface: {phi_far} should approach {phi_eq}"
        );
    }
    #[test]
    fn test_multi_range_uniform_zero() {
        let nx = 8;
        let ny = 8;
        let psi = vec![1.0; nx * ny];
        let (fx, fy) = multi_range_interaction_force(&psi, 4, 4, nx, ny, -1.0, -0.5);
        assert!(fx.abs() < 1e-12, "Uniform: fx should be zero: {fx}");
        assert!(fy.abs() < 1e-12, "Uniform: fy should be zero: {fy}");
    }
    #[test]
    fn test_count_phases() {
        let rho = vec![0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5];
        let (n_l, n_i, n_v) = count_phases(&rho, 0.4, 1.6);
        assert_eq!(n_v, 2, "vapor count");
        assert_eq!(n_i, 3, "interface count");
        assert_eq!(n_l, 2, "liquid count");
    }
    #[test]
    fn test_bulk_free_energy_at_equilibrium() {
        let a = -0.2;
        let b = 0.1;
        let model = FreeEnergyModel::new(0.01, a, b);
        let phi_eq = model.rho_liquid();
        let f_eq = model.bulk_free_energy_density(phi_eq);
        let f_0 = model.bulk_free_energy_density(0.0);
        assert!(
            f_eq < f_0,
            "Free energy at equilibrium ({f_eq}) should be lower than at phi=0 ({f_0})"
        );
    }
    #[test]
    fn test_chemical_potential_at_equilibrium() {
        let a = -0.2;
        let b = 0.1;
        let model = FreeEnergyModel::new(0.01, a, b);
        let phi_eq = model.rho_liquid();
        let mu = model.chemical_potential(phi_eq, 0.0);
        assert!(
            mu.abs() < 1e-12,
            "Chemical potential at equilibrium should be zero: {mu}"
        );
    }
    #[test]
    fn test_clip_density() {
        assert!((clip_density(0.5, 0.1, 2.0) - 0.5).abs() < 1e-14);
        assert!((clip_density(-0.1, 0.01, 5.0) - 0.01).abs() < 1e-14);
        assert!((clip_density(10.0, 0.01, 5.0) - 5.0).abs() < 1e-14);
    }
    #[test]
    fn test_spinodal_params_dominant_wavenumber() {
        let p = SpinodалDecompositionParams::new(-0.1, 0.1, 0.01, 1.0);
        let k = p.dominant_wavenumber();
        let expected = (0.1_f64 / 0.02).sqrt();
        assert!(
            (k - expected).abs() < 1e-10,
            "k* = {k}, expected {expected}"
        );
    }
    #[test]
    fn test_spinodal_params_dominant_wavelength() {
        let p = SpinodалDecompositionParams::new(-0.1, 0.1, 0.01, 1.0);
        let lambda = p.dominant_wavelength();
        let k = p.dominant_wavenumber();
        let expected = std::f64::consts::TAU / k;
        assert!(
            (lambda - expected).abs() < 1e-8,
            "lambda = {lambda}, expected {expected}"
        );
    }
    #[test]
    fn test_spinodal_params_no_unstable_state() {
        let p = SpinodалDecompositionParams::new(0.1, 0.1, 0.01, 1.0);
        assert!((p.dominant_wavenumber()).abs() < 1e-14);
        assert_eq!(p.dominant_wavelength(), f64::INFINITY);
    }
    #[test]
    fn test_laplacian_2d_constant_field() {
        let nx = 8;
        let ny = 8;
        let phi = vec![2.0f64; nx * ny];
        let lap = laplacian_2d_periodic(&phi, nx, ny);
        for (k, &l) in lap.iter().enumerate() {
            assert!(l.abs() < 1e-14, "Laplacian of constant field at {k}: {l}");
        }
    }
    #[test]
    fn test_cahn_hilliard_step_conserves_mass() {
        let nx = 16;
        let ny = 16;
        let mut phi: Vec<f64> = (0..nx * ny)
            .map(|k| {
                let x = (k % nx) as f64;
                0.01 * (x * 0.7).sin()
            })
            .collect();
        let params = SpinodалDecompositionParams::new(-0.1, 0.1, 0.01, 0.1);
        let mass_before: f64 = phi.iter().sum();
        cahn_hilliard_step(&mut phi, nx, ny, 0.001, &params);
        let mass_after: f64 = phi.iter().sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "Cahn-Hilliard must conserve mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_total_free_energy_decreases() {
        let nx = 16;
        let ny = 16;
        let a = -0.1;
        let b = 0.1;
        let kappa = 0.01;
        let mut phi: Vec<f64> = (0..nx * ny)
            .map(|k| {
                let x = (k % nx) as f64 / nx as f64;
                0.05 * (x * std::f64::consts::TAU).sin()
            })
            .collect();
        let params = SpinodалDecompositionParams::new(a, b, kappa, 0.1);
        let f0 = total_free_energy_2d(&phi, nx, ny, a, b, kappa);
        for _ in 0..10 {
            cahn_hilliard_step(&mut phi, nx, ny, 1e-4, &params);
        }
        let f1 = total_free_energy_2d(&phi, nx, ny, a, b, kappa);
        assert!(
            f1 <= f0 + 1e-10,
            "Free energy should not increase: before={f0}, after={f1}"
        );
    }
    #[test]
    fn test_laplace_pressure_2d() {
        let sigma = 0.072;
        let r = 1e-3;
        let dp = laplace_pressure_2d(sigma, r);
        assert!(
            (dp - 144.0).abs() < 1e-6,
            "Laplace pressure = {dp}, expected 144"
        );
    }
    #[test]
    fn test_laplace_pressure_3d_droplet() {
        let sigma = 0.072;
        let r = 0.001;
        let dp = laplace_pressure_3d_droplet(sigma, r);
        assert!((dp - 144.0).abs() < 1e-6, "3D Laplace pressure = {dp}");
    }
    #[test]
    fn test_bubble_natural_frequency_positive() {
        let f = bubble_natural_frequency(1e-4, 1000.0, 1e5, 1.4);
        assert!(f > 0.0, "Bubble frequency should be positive: {f}");
    }
    #[test]
    fn test_droplet_terminal_velocity_sign() {
        let u = droplet_terminal_velocity(1e-3, 1000.0, 800.0, 1e-3, 5e-3, 9.81);
        assert!(
            u > 0.0,
            "Terminal velocity should be positive (rising): {u}"
        );
    }
    #[test]
    fn test_droplet_drag_coefficient_stokes() {
        let cd = droplet_drag_coefficient(1.0, 1e10, 1.0);
        assert!(cd.is_finite() && cd > 0.0, "Cd = {cd}");
    }
    #[test]
    fn test_bubble_state_equilibrium() {
        let b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 1e-3);
        let p_gas = b.gas_pressure();
        let p_expected = b.p_inf + 2.0 * b.sigma / b.r_eq;
        assert!(
            (p_gas - p_expected).abs() < 1e-6,
            "Gas pressure at eq: {p_gas} vs {p_expected}"
        );
    }
    #[test]
    fn test_bubble_step_finite() {
        let mut b = BubbleState::new(1e-3, 1e5, 1000.0, 0.072, 1e-3);
        b.radius = 1.1e-3;
        for _ in 0..10 {
            b.step(1e-9);
        }
        assert!(
            b.radius.is_finite() && b.radius > 0.0,
            "Bubble radius = {}",
            b.radius
        );
    }
    #[test]
    fn test_collapse_time_estimate() {
        let r_eq = 1e-3;
        let sigma = 0.0;
        let p_inf = 1e5;
        let rho_liq = 1000.0;
        let mut b = BubbleState::new(r_eq, p_inf, rho_liq, sigma, 1e-3);
        b.radius = 1.5 * r_eq;
        let p_gas = b.gas_pressure();
        let p_drive = p_inf - p_gas;
        assert!(p_drive > 0.0, "p_drive = {p_drive}");
        let t = b.collapse_time_estimate();
        assert!(t > 0.0 && t.is_finite(), "Collapse time = {t}");
        let mut b_big = BubbleState::new(2.0 * r_eq, p_inf, rho_liq, sigma, 1e-3);
        b_big.radius = 1.5 * (2.0 * r_eq);
        let t_big = b_big.collapse_time_estimate();
        assert!(
            t_big > t,
            "Larger bubble collapses slower: t={t}, t_big={t_big}"
        );
    }
    #[test]
    fn test_phase_change_source_above_sat() {
        let s = phase_change_source(105.0, 100.0, 1.0);
        assert!(s > 0.0, "Above T_sat: evaporation (positive source): {s}");
    }
    #[test]
    fn test_phase_change_source_below_sat() {
        let s = phase_change_source(95.0, 100.0, 1.0);
        assert!(s < 0.0, "Below T_sat: condensation (negative source): {s}");
    }
    #[test]
    fn test_latent_heat_source() {
        let q = latent_heat_source(2.0, 2257.0);
        assert!((q - 4514.0).abs() < 1e-6, "Q = {q}, expected 4514");
    }
    #[test]
    fn test_apply_phase_change_evaporates() {
        let mut rho = vec![1.0; 5];
        let temp = vec![105.0; 5];
        apply_phase_change(&mut rho, &temp, 100.0, 0.01, 1.0);
        for &r in &rho {
            assert!(r > 1.0, "Evaporation should increase vapor density");
        }
    }
    #[test]
    fn test_apply_phase_change_no_negative() {
        let mut rho = vec![0.001; 5];
        let temp = vec![50.0; 5];
        apply_phase_change(&mut rho, &temp, 100.0, 1000.0, 1.0);
        for &r in &rho {
            assert!(r >= 0.0, "Density should never go negative: {r}");
        }
    }
    #[test]
    fn test_volume_fraction_liquid() {
        let vf = volume_fraction(2.0, 2.0, 0.1, 0.1);
        assert!(
            vf > 0.9,
            "At rho_liquid, volume fraction should be ~1: {vf}"
        );
    }
    #[test]
    fn test_volume_fraction_vapor() {
        let vf = volume_fraction(0.1, 2.0, 0.1, 0.1);
        assert!(vf < 0.1, "At rho_vapor, volume fraction should be ~0: {vf}");
    }
    #[test]
    fn test_interface_normal_2d_uniform() {
        let nx = 4;
        let ny = 4;
        let phi = vec![1.0f64; nx * ny];
        let n = interface_normal_2d(&phi, 2, 2, nx, ny);
        assert!(
            n[0].abs() < 1e-14 && n[1].abs() < 1e-14,
            "Uniform field: normal = 0"
        );
    }
    #[test]
    fn test_interface_curvature_2d_uniform() {
        let nx = 8;
        let ny = 8;
        let phi = vec![1.0f64; nx * ny];
        let kap = interface_curvature_2d(&phi, 2, 2, nx, ny);
        assert!(kap.abs() < 1e-14, "Uniform field: curvature = 0, got {kap}");
    }
    #[test]
    fn test_film_boiling_heat_transfer_positive() {
        let h = film_boiling_heat_transfer(0.025, 0.6, 1000.0, 1.2e-5, 9.81, 2257e3, 20.0, 0.01);
        assert!(h > 0.0, "Film boiling h should be positive: {h}");
    }
    #[test]
    fn test_hertz_knudsen_evaporation_rate_above_sat() {
        let j = hertz_knudsen_evaporation_rate(105.0, 100.0, 0.5, 0.018, 8.314);
        assert!(j.is_finite(), "Evaporation rate should be finite");
    }
    #[test]
    fn test_leidenfrost_superheat_positive() {
        let dt = leidenfrost_superheat(0.058, 0.6, 2000.0, 0.025);
        assert!(dt.is_finite() && dt >= 0.0, "Leidenfrost superheat = {dt}");
    }
    #[test]
    fn test_sukop_thorne_psi_range() {
        let rho0 = 1.0;
        let psi0 = 4.0;
        let psi_low = sukop_thorne_psi(0.1, rho0, psi0);
        let psi_mid = sukop_thorne_psi(1.0, rho0, psi0);
        let psi_high = sukop_thorne_psi(10.0, rho0, psi0);
        assert!(
            psi_low < psi_mid && psi_mid < psi_high,
            "ST psi not monotone"
        );
        assert!(psi_high < psi0, "ST psi approaches psi0 from below");
    }
    #[test]
    fn test_sukop_thorne_critical_g_negative() {
        let g_c = sukop_thorne_critical_g(1.0);
        assert!(g_c < 0.0, "G_c should be negative: {g_c}");
    }
    #[test]
    fn test_multicomponent_sc_psi() {
        let mc = MultiComponentSC::new(-0.5, -0.5, -1.5, 1.0);
        let psi_small = mc.psi(0.01);
        let psi_large = mc.psi(100.0);
        assert!(psi_small < psi_large, "psi should increase with density");
    }
    #[test]
    fn test_multicomponent_sc_density_ratio() {
        let mc = MultiComponentSC::new(0.0, 0.0, -1.5, 1.0);
        let dr = mc.density_ratio_estimate();
        assert!(dr > 1.0, "Density ratio should exceed 1 for g_cc < 0: {dr}");
    }
    #[test]
    fn test_multicomponent_sc_uniform_zero_force() {
        let nx = 6;
        let ny = 6;
        let mc = MultiComponentSC::new(0.0, 0.0, -1.5, 1.0);
        let psi_a = vec![1.0_f64; nx * ny];
        let psi_b = vec![1.0_f64; nx * ny];
        let (fx, fy) = mc.interaction_force(&psi_a, &psi_b, 3, 3, nx, ny);
        assert!(fx.abs() < 1e-13, "Uniform: fx should be ~0: {fx}");
        assert!(fy.abs() < 1e-13, "Uniform: fy should be ~0: {fy}");
    }
    #[test]
    fn test_droplets_coalescing_overlapping() {
        assert!(
            droplets_coalescing(0.0, 0.0, 0.8, 1.0, 0.0, 0.8, 0.0),
            "Overlapping droplets should be detected as coalescing"
        );
    }
    #[test]
    fn test_droplets_coalescing_separated() {
        assert!(
            !droplets_coalescing(0.0, 0.0, 1.0, 10.0, 0.0, 1.0, 0.0),
            "Separated droplets should not be coalescing"
        );
    }
    #[test]
    fn test_coalesced_radius_2d_conservation() {
        let r_new = coalesced_radius_2d(3.0, 4.0);
        assert!(
            (r_new - 5.0).abs() < 1e-12,
            "coalesced radius = {r_new}, expected 5"
        );
    }
    #[test]
    fn test_coalesced_radius_3d_conservation() {
        let r_new = coalesced_radius_3d(1.0, 1.0);
        let expected = 2.0_f64.powf(1.0 / 3.0);
        assert!(
            (r_new - expected).abs() < 1e-12,
            "r_new = {r_new}, expected {expected}"
        );
    }
    #[test]
    fn test_coalescence_time_scale_positive() {
        let t = coalescence_time_scale(1e-3, 1e-3, 1e-3, 0.072, 1e-6);
        assert!(
            t > 0.0 && t.is_finite(),
            "Coalescence time should be positive: {t}"
        );
    }
    #[test]
    fn test_phase_separation_rate_constant_positive() {
        let rate = phase_separation_rate_constant(1.0, -0.1, 0.01);
        assert!(
            rate > 0.0,
            "Rate constant should be positive inside spinodal: {rate}"
        );
    }
    #[test]
    fn test_phase_separation_rate_constant_stable() {
        let rate = phase_separation_rate_constant(1.0, 0.1, 0.01);
        assert!(
            rate.abs() < 1e-14,
            "Rate constant should be zero outside spinodal: {rate}"
        );
    }
    #[test]
    fn test_spinodal_number_gt_one_inside() {
        let sp = spinodal_number(-1.0, 0.01, 0.1);
        assert!(
            sp > 1.0,
            "Spinodal number should be > 1 inside spinodal: {sp}"
        );
    }
    #[test]
    fn test_surface_tension_from_laplace_2d() {
        let sigma = surface_tension_from_laplace(1.072, 1.0, 1e-3, 2);
        assert!((sigma - 0.072e-3).abs() < 1e-12, "sigma = {sigma}");
    }
    #[test]
    fn test_equilibrium_droplet_radius() {
        let r_eq = equilibrium_droplet_radius(0.072, 1000.144, 1000.0);
        assert!((r_eq - 1.0).abs() < 1e-10, "R_eq = {r_eq}, expected ~1.0");
    }
    #[test]
    fn test_diffuse_interface_width_90_positive() {
        let w = diffuse_interface_width_90(0.01, -0.1);
        assert!(w > 0.0, "Interface width should be positive: {w}");
    }
    #[test]
    fn test_diffuse_interface_width_90_stable() {
        let w = diffuse_interface_width_90(0.01, 0.1);
        assert!(w.abs() < 1e-14, "No interface for A > 0: {w}");
    }
    #[test]
    fn test_interface_is_resolved() {
        assert!(
            interface_is_resolved(0.01, -0.1, 0.1, 4),
            "Interface should be resolved"
        );
        assert!(
            !interface_is_resolved(0.01, -0.1, 0.5, 4),
            "Interface should NOT be resolved"
        );
    }
    #[test]
    fn test_tie_line_symmetric() {
        let (phi_v, phi_l) = tie_line_endpoints(-0.2, 0.1);
        assert!(phi_l > 0.0 && phi_v < 0.0, "phi_l > 0, phi_v < 0");
        assert!((phi_l + phi_v).abs() < 1e-14, "Symmetric: phi_l = -phi_v");
    }
    #[test]
    fn test_lever_rule_midpoint() {
        let (phi_v, phi_l) = tie_line_endpoints(-0.2, 0.1);
        let phi_mid = 0.5 * (phi_v + phi_l);
        let frac = lever_rule_liquid_fraction(phi_mid, phi_v, phi_l);
        assert!((frac - 0.5).abs() < 1e-12, "Midpoint → 50/50 split: {frac}");
    }
    #[test]
    fn test_lever_rule_clamp() {
        let (phi_v, phi_l) = tie_line_endpoints(-0.2, 0.1);
        let frac = lever_rule_liquid_fraction(phi_l + 1.0, phi_v, phi_l);
        assert!((frac - 1.0).abs() < 1e-12, "Clamped to 1.0: {frac}");
    }
    #[test]
    fn test_nucleation_seed_peak() {
        let nx = 16;
        let ny = 16;
        let mut field = vec![0.0_f64; nx * ny];
        add_nucleation_seed(&mut field, nx, ny, 8.0, 8.0, 1.0, 1.0);
        let peak = field[8 * nx + 8];
        assert!(
            (peak - 1.0).abs() < 1e-12,
            "Peak should be at amplitude=1.0: {peak}"
        );
    }
    #[test]
    fn test_nucleation_seed_decay() {
        let nx = 32;
        let ny = 32;
        let mut field = vec![0.0_f64; nx * ny];
        add_nucleation_seed(&mut field, nx, ny, 16.0, 16.0, 1.0, 2.0);
        let peak = field[16 * nx + 16];
        let off = field[16 * nx + 20];
        assert!(
            peak > off,
            "Seed should decay away from centre: peak={peak}, off={off}"
        );
    }
    #[test]
    fn test_nucleation_seeds_multiple() {
        let nx = 16;
        let ny = 16;
        let mut field = vec![0.0_f64; nx * ny];
        let centres = [(4.0_f64, 4.0_f64), (12.0, 12.0)];
        add_nucleation_seeds(&mut field, nx, ny, &centres, 0.5, 1.0);
        let p1 = field[4 * nx + 4];
        let p2 = field[12 * nx + 12];
        assert!((p1 - 0.5).abs() < 1e-12, "Seed 1 peak: {p1}");
        assert!((p2 - 0.5).abs() < 1e-12, "Seed 2 peak: {p2}");
    }
    #[test]
    fn test_wetting_psi_from_angle_fully_wetting() {
        let psi_w = wetting_psi_from_angle(2.0, 0.1, 0.0);
        assert!(
            (psi_w - 2.0).abs() < 1e-12,
            "theta=0: psi_wall = psi_liq: {psi_w}"
        );
    }
    #[test]
    fn test_wetting_psi_from_angle_nonwetting() {
        let psi_w = wetting_psi_from_angle(2.0, 0.1, std::f64::consts::PI);
        assert!(
            (psi_w - 0.1).abs() < 1e-12,
            "theta=pi: psi_wall = psi_vap: {psi_w}"
        );
    }
    #[test]
    fn test_contact_angle_round_trip() {
        let theta = std::f64::consts::FRAC_PI_4;
        let psi_w = wetting_psi_from_angle(2.0, 0.1, theta);
        let theta2 = contact_angle_from_psi(psi_w, 2.0, 0.1);
        assert!(
            (theta - theta2).abs() < 1e-12,
            "Round-trip: theta={theta}, theta2={theta2}"
        );
    }
    #[test]
    fn test_apply_contact_angle_bc() {
        let mut psi = vec![0.5_f64; 4];
        let is_solid = vec![false, true, false, true];
        apply_contact_angle_bc(&mut psi, &is_solid, 2.0, 0.1, 0.0);
        assert!(
            (psi[1] - 2.0).abs() < 1e-12,
            "psi[1] should be psi_liq: {}",
            psi[1]
        );
        assert!(
            (psi[3] - 2.0).abs() < 1e-12,
            "psi[3] should be psi_liq: {}",
            psi[3]
        );
        assert!((psi[0] - 0.5).abs() < 1e-12, "psi[0] unchanged: {}", psi[0]);
    }
    #[test]
    fn test_g_cc_for_density_ratio_negative() {
        let g_cc = g_cc_for_density_ratio(5.0, 1.0);
        assert!(
            g_cc < 0.0,
            "G_cc for immiscibility should be negative: {g_cc}"
        );
    }
    #[test]
    fn test_g_cc_identity_for_dr_one() {
        let g_cc = g_cc_for_density_ratio(1.0, 1.0);
        assert!(g_cc.abs() < 1e-10, "G_cc for DR=1 should be ~0: {g_cc}");
    }
    #[test]
    fn test_g_cc_is_stable_small() {
        assert!(g_cc_is_stable(-0.1, 1.0), "Small G_cc should be stable");
    }
    #[test]
    fn test_g_cc_is_stable_large() {
        assert!(!g_cc_is_stable(-1.0, 1.0), "Large G_cc should be unstable");
    }
}
