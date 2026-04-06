//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use std::f64::consts::PI;

/// Constant mobility model.
pub fn constant_mobility(m: f64, _phi: f64) -> f64 {
    m
}
/// Degenerate mobility: M(phi) = M_0 * phi * (1 - phi).
///
/// This ensures zero mobility in the bulk phases (phi=0 and phi=1),
/// confining diffusion to the interface region.
pub fn degenerate_mobility(m0: f64, phi: f64) -> f64 {
    m0 * phi * (1.0 - phi)
}
/// Concentration-dependent mobility: M(phi) = M_0 * (1 + beta * phi).
pub fn concentration_dependent_mobility(m0: f64, beta: f64, phi: f64) -> f64 {
    m0 * (1.0 + beta * phi)
}
/// Compute the body force on the fluid due to the phase field
/// (Korteweg stress / capillary force).
///
/// ```text
/// F_x = -phi * d(mu)/dx
/// F_y = -phi * d(mu)/dy
/// ```
///
/// This couples the phase field evolution to the Navier-Stokes momentum equation.
pub fn phase_field_body_force(phi: &[f64], mu: &[f64], nx: usize, ny: usize) -> Vec<[f64; 2]> {
    let n = nx * ny;
    let mut force = vec![[0.0_f64; 2]; n];
    let idx = |x: usize, y: usize| y * nx + x;
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let xm = if x > 0 { x - 1 } else { x };
            let xp = if x < nx - 1 { x + 1 } else { x };
            let ym = if y > 0 { y - 1 } else { y };
            let yp = if y < ny - 1 { y + 1 } else { y };
            let dx_span = (xp - xm).max(1) as f64;
            let dy_span = (yp - ym).max(1) as f64;
            let dmu_dx = (mu[idx(xp, y)] - mu[idx(xm, y)]) / dx_span;
            let dmu_dy = (mu[idx(x, yp)] - mu[idx(x, ym)]) / dy_span;
            force[k][0] = -phi[k] * dmu_dx;
            force[k][1] = -phi[k] * dmu_dy;
        }
    }
    force
}
/// Advect the phase field by a velocity field using upwind differencing.
///
/// ```text
/// d(phi)/dt + u . grad(phi) = 0
/// ```
///
/// Returns the advection contribution (to be subtracted from the RHS).
pub fn advect_phase_field(phi: &[f64], ux: &[f64], uy: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    let n = nx * ny;
    let mut advection = vec![0.0_f64; n];
    let idx = |x: usize, y: usize| y * nx + x;
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let xm = if x > 0 { x - 1 } else { x };
            let xp = if x < nx - 1 { x + 1 } else { x };
            let ym = if y > 0 { y - 1 } else { y };
            let yp = if y < ny - 1 { y + 1 } else { y };
            let dphi_dx = if ux[k] > 0.0 {
                phi[k] - phi[idx(xm, y)]
            } else {
                phi[idx(xp, y)] - phi[k]
            };
            let dphi_dy = if uy[k] > 0.0 {
                phi[k] - phi[idx(x, ym)]
            } else {
                phi[idx(x, yp)] - phi[k]
            };
            advection[k] = ux[k] * dphi_dx + uy[k] * dphi_dy;
        }
    }
    advection
}
/// Compute the density field from the phase field using linear interpolation.
///
/// rho(phi) = rho_1 + (rho_2 - rho_1) * phi
pub fn density_from_phase_field(phi: &[f64], rho_1: f64, rho_2: f64) -> Vec<f64> {
    phi.iter().map(|&p| rho_1 + (rho_2 - rho_1) * p).collect()
}
/// Compute the dynamic viscosity field from the phase field.
///
/// mu(phi) = mu_1 + (mu_2 - mu_1) * phi
pub fn viscosity_from_phase_field(phi: &[f64], mu_1: f64, mu_2: f64) -> Vec<f64> {
    phi.iter().map(|&p| mu_1 + (mu_2 - mu_1) * p).collect()
}
/// Identify interface cells where phi crosses a threshold.
///
/// A cell is an interface cell if any of its neighbors has a phi value
/// on the opposite side of the threshold.
pub fn find_interface_cells(phi: &[f64], nx: usize, ny: usize, threshold: f64) -> Vec<usize> {
    let idx = |x: usize, y: usize| y * nx + x;
    let mut interface = Vec::new();
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let above = phi[k] > threshold;
            let mut is_interface = false;
            if x > 0 && (phi[idx(x - 1, y)] > threshold) != above {
                is_interface = true;
            }
            if x < nx - 1 && (phi[idx(x + 1, y)] > threshold) != above {
                is_interface = true;
            }
            if y > 0 && (phi[idx(x, y - 1)] > threshold) != above {
                is_interface = true;
            }
            if y < ny - 1 && (phi[idx(x, y + 1)] > threshold) != above {
                is_interface = true;
            }
            if is_interface {
                interface.push(k);
            }
        }
    }
    interface
}
/// Compute the phase field interface width from a 1D profile.
///
/// Measures the distance over which phi goes from 0.1 to 0.9 along a column.
pub fn measure_interface_width(phi: &[f64], nx: usize, ny: usize, x_col: usize) -> f64 {
    let idx = |y: usize| y * nx + x_col;
    let mut y_low = None;
    let mut y_high = None;
    for y in 0..ny {
        let p = phi[idx(y)];
        if y_low.is_none() && p > 0.1 {
            y_low = Some(y as f64);
        }
        if y_high.is_none() && p > 0.9 {
            y_high = Some(y as f64);
        }
    }
    match (y_low, y_high) {
        (Some(l), Some(h)) => (h - l).abs(),
        _ => 0.0,
    }
}
/// Advance the order-parameter field `phi` by one step under the
/// convective Cahn-Hilliard equation:
///
/// `d(phi)/dt + u . grad(phi) = M * lap(mu)`
///
/// where `mu = A*phi + B*phi^3 - kappa * lap(phi)`.
///
/// Uses a central-difference spatial discretisation with periodic BCs.
/// Returns the updated phi field.
#[allow(dead_code)]
pub fn convective_cahn_hilliard_step(
    phi: &[f64],
    ux: &[f64],
    uy: &[f64],
    nx: usize,
    ny: usize,
    dt: f64,
    a: f64,
    b: f64,
    kappa: f64,
    mobility: f64,
) -> Vec<f64> {
    let n = nx * ny;
    let idx = |x: usize, y: usize| y * nx + x;
    let lap = |field: &[f64], x: usize, y: usize| -> f64 {
        let xp = (x + 1) % nx;
        let xm = (x + nx - 1) % nx;
        let yp = (y + 1) % ny;
        let ym = (y + ny - 1) % ny;
        field[idx(xp, y)] + field[idx(xm, y)] + field[idx(x, yp)] + field[idx(x, ym)]
            - 4.0 * field[idx(x, y)]
    };
    let mut mu = vec![0.0_f64; n];
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let p = phi[k];
            mu[k] = a * p + b * p * p * p - kappa * lap(phi, x, y);
        }
    }
    let mut phi_new = phi.to_vec();
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let diff = mobility * lap(&mu, x, y);
            let xp = (x + 1) % nx;
            let xm = (x + nx - 1) % nx;
            let yp = (y + 1) % ny;
            let ym = (y + ny - 1) % ny;
            let dphi_dx = 0.5 * (phi[idx(xp, y)] - phi[idx(xm, y)]);
            let dphi_dy = 0.5 * (phi[idx(x, yp)] - phi[idx(x, ym)]);
            let adv = ux[k] * dphi_dx + uy[k] * dphi_dy;
            phi_new[k] += dt * (diff - adv);
        }
    }
    phi_new
}
/// Compute the Cahn number: `Cn = epsilon / L` where `L` is the macroscopic
/// length scale.
///
/// The Cahn number controls the ratio of interface thickness to domain size.
/// For convergence to the sharp-interface limit, `Cn → 0` is required.
/// Typical values: `Cn = 0.01 – 0.1`.
#[allow(dead_code)]
pub fn cahn_number(epsilon: f64, domain_length: f64) -> f64 {
    if domain_length < 1e-30 {
        return f64::INFINITY;
    }
    epsilon / domain_length
}
/// Compute the Peclet number for phase-field convection.
///
/// `Pe = U * L / (M * |A|)` where `U` is the characteristic velocity,
/// `L` is the domain length, and `M * |A|` is the effective diffusivity.
#[allow(dead_code)]
pub fn phase_field_peclet(u_ref: f64, domain_length: f64, mobility: f64, a: f64) -> f64 {
    let diff = mobility * a.abs();
    if diff < 1e-30 {
        return f64::INFINITY;
    }
    u_ref * domain_length / diff
}
/// Compute the required grid resolution to resolve the diffuse interface.
///
/// A minimum of `n_points_per_xi` grid points should span the interface
/// width `xi = sqrt(kappa / |A|)`.
///
/// Returns the maximum allowed grid spacing `dx_max`.
#[allow(dead_code)]
pub fn max_grid_spacing(kappa: f64, a: f64, n_points_per_xi: usize) -> f64 {
    if a >= 0.0 || kappa <= 0.0 || n_points_per_xi == 0 {
        return f64::INFINITY;
    }
    let xi = (kappa / (-a)).sqrt();
    xi / n_points_per_xi as f64
}
/// Compute the free energy density for the double-well potential.
///
/// `f(phi) = A/2 * phi^2 + B/4 * phi^4`
#[allow(dead_code)]
pub fn double_well_free_energy(phi: f64, a: f64, b: f64) -> f64 {
    0.5 * a * phi * phi + 0.25 * b * phi.powi(4)
}
/// Compute the nucleation energy barrier for a spherical nucleus of radius R.
///
/// Using classical nucleation theory:
/// `DeltaG = -4/3 * pi * R^3 * Delta f_v + 4 * pi * R^2 * sigma`
///
/// where `Delta f_v` is the bulk free energy driving force and `sigma` the
/// surface tension.
///
/// Returns the activation energy `Delta G`.
#[allow(dead_code)]
pub fn nucleation_energy_barrier(radius: f64, delta_f_v: f64, sigma: f64) -> f64 {
    use std::f64::consts::PI;
    -4.0 / 3.0 * PI * radius.powi(3) * delta_f_v + 4.0 * PI * radius * radius * sigma
}
/// Critical nucleus radius from classical nucleation theory.
///
/// `R_crit = 2 * sigma / Delta_f_v`
///
/// Nuclei smaller than `R_crit` dissolve; larger nuclei grow.
#[allow(dead_code)]
pub fn critical_nucleus_radius(sigma: f64, delta_f_v: f64) -> f64 {
    if delta_f_v < 1e-30 {
        return f64::INFINITY;
    }
    2.0 * sigma / delta_f_v
}
/// Maximum nucleation energy barrier at the critical radius.
///
/// `Delta G_crit = (16 * pi * sigma^3) / (3 * Delta_f_v^2)`
#[allow(dead_code)]
pub fn critical_nucleation_barrier(sigma: f64, delta_f_v: f64) -> f64 {
    if delta_f_v < 1e-30 {
        return f64::INFINITY;
    }
    16.0 * PI * sigma.powi(3) / (3.0 * delta_f_v * delta_f_v)
}
/// Compute the Tolman length correction to surface tension for a curved interface.
///
/// `sigma(R) = sigma_flat / (1 + 2 * delta / R)`
///
/// where `delta` is the Tolman length (typically small).
#[allow(dead_code)]
pub fn tolman_surface_tension(sigma_flat: f64, radius: f64, tolman_length: f64) -> f64 {
    let corr = 1.0 + 2.0 * tolman_length / radius.max(1e-30);
    sigma_flat / corr
}
/// Volume-fraction-weighted mixture density.
///
/// `rho_mix = sum_k (phi_k * rho_k)`
///
/// where `phi_k` are the volume fractions of each component.
#[allow(dead_code)]
pub fn mixture_density(phi_components: &[f64], rho_components: &[f64]) -> f64 {
    phi_components
        .iter()
        .zip(rho_components.iter())
        .map(|(&p, &r)| p * r)
        .sum()
}
/// Volume-fraction-weighted mixture viscosity (linear blending).
#[allow(dead_code)]
pub fn mixture_viscosity(phi_components: &[f64], mu_components: &[f64]) -> f64 {
    phi_components
        .iter()
        .zip(mu_components.iter())
        .map(|(&p, &m)| p * m)
        .sum()
}
/// Check that the phase-field order parameters sum to unity (partition of unity).
///
/// Returns the maximum deviation from 1 across all grid cells.
#[allow(dead_code)]
pub fn partition_of_unity_error(phi_fields: &[&[f64]], n_cells: usize) -> f64 {
    let mut max_err = 0.0_f64;
    for k in 0..n_cells {
        let sum: f64 = phi_fields.iter().map(|phi| phi[k]).sum();
        max_err = max_err.max((sum - 1.0).abs());
    }
    max_err
}
/// Enforce the partition-of-unity constraint by renormalising the phase fields.
///
/// After each step, `phi_k` may drift; this function redistributes the
/// residual equally among all components.
#[allow(dead_code)]
pub fn enforce_partition_of_unity(phi_fields: &mut [Vec<f64>], n_cells: usize) {
    let n_comp = phi_fields.len();
    if n_comp == 0 {
        return;
    }
    for k in 0..n_cells {
        let sum: f64 = phi_fields.iter().map(|phi| phi[k]).sum();
        if sum.abs() > 1e-30 {
            for phi in phi_fields.iter_mut() {
                phi[k] /= sum;
            }
        }
    }
}
/// Cahn-Hilliard free energy functional (discrete, 2D).
///
/// `F[phi] = sum_{x,y} [ f_bulk(phi) + kappa/2 * |grad phi|^2 ]`
///
/// Uses forward-difference gradients with periodic BCs.
#[allow(dead_code)]
pub fn cahn_hilliard_free_energy(
    phi: &[f64],
    nx: usize,
    ny: usize,
    a: f64,
    b: f64,
    kappa: f64,
) -> f64 {
    let idx = |x: usize, y: usize| y * nx + x;
    let mut energy = 0.0_f64;
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let xp = (x + 1) % nx;
            let yp = (y + 1) % ny;
            let p = phi[k];
            let dphi_dx = phi[idx(xp, y)] - p;
            let dphi_dy = phi[idx(x, yp)] - p;
            energy += double_well_free_energy(p, a, b)
                + 0.5 * kappa * (dphi_dx * dphi_dx + dphi_dy * dphi_dy);
        }
    }
    energy
}
/// Compute the variational derivative `delta F / delta phi = mu`.
///
/// `mu(x) = A*phi + B*phi^3 - kappa * lap(phi)`
///
/// Returns the chemical potential field.
#[allow(dead_code)]
pub fn cahn_hilliard_chemical_potential(
    phi: &[f64],
    nx: usize,
    ny: usize,
    a: f64,
    b: f64,
    kappa: f64,
) -> Vec<f64> {
    let n = nx * ny;
    let idx = |x: usize, y: usize| y * nx + x;
    let mut mu = vec![0.0_f64; n];
    for y in 0..ny {
        for x in 0..nx {
            let k = idx(x, y);
            let p = phi[k];
            let xp = (x + 1) % nx;
            let xm = (x + nx - 1) % nx;
            let yp = (y + 1) % ny;
            let ym = (y + ny - 1) % ny;
            let lap =
                phi[idx(xp, y)] + phi[idx(xm, y)] + phi[idx(x, yp)] + phi[idx(x, ym)] - 4.0 * p;
            mu[k] = a * p + b * p * p * p - kappa * lap;
        }
    }
    mu
}
