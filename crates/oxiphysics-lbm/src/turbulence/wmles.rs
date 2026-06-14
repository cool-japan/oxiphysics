// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wall-Modeled LES (WMLES) boundary condition coupling for LBM channel flow.
//!
//! Provides:
//! - [`WmlesChannelCoupling`]: Coupling struct with WMLES model and Smagorinsky constant.
//! - [`collide_stream_wmles`]: Combined collision and streaming with WMLES closure.
//! - [`wmles_eddy_viscosity_profile`]: Wall-normal eddy viscosity profile.

use crate::turbulent_channel::ChannelFlow;
use crate::wall_model::WallModeledLes;

// D2Q9 constants (local copies since turbulent_channel's are pub(super))
const W9_LOCAL: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
const CX_LOCAL: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const CY_LOCAL: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const BB_LOCAL: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// WMLES coupling parameters for an LBM channel flow simulation.
///
/// Combines a [`WallModeledLes`] model for near-wall cells with a Smagorinsky
/// SGS closure for the interior of the domain.
pub struct WmlesChannelCoupling {
    /// The wall-modeled LES instance providing wall stress and eddy viscosity.
    pub wall_model: WallModeledLes,
    /// Smagorinsky constant Cs for the SGS model in interior cells.
    pub cs_smag: f64,
    /// Molecular kinematic viscosity.
    pub nu: f64,
}

impl WmlesChannelCoupling {
    /// Create a new WMLES channel coupling.
    ///
    /// # Arguments
    /// * `nu`      - Molecular kinematic viscosity (lattice units)
    /// * `cs_smag` - Smagorinsky constant (typically 0.1–0.2)
    pub fn new(nu: f64, cs_smag: f64) -> Self {
        Self {
            wall_model: WallModeledLes::new(nu),
            cs_smag,
            nu,
        }
    }
}

/// D2Q9 equilibrium distribution.
#[inline]
fn equilibrium_local(rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
    let cx = CX_LOCAL[q] as f64;
    let cy = CY_LOCAL[q] as f64;
    let eu = cx * ux + cy * uy;
    let u2 = ux * ux + uy * uy;
    W9_LOCAL[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2)
}

/// Compute macroscopic moments from a 9-element slice.
#[inline]
fn moments_local(f: &[f64], base: usize) -> (f64, f64, f64) {
    let mut rho = 0.0_f64;
    let mut mx = 0.0_f64;
    let mut my = 0.0_f64;
    for q in 0..9 {
        let fq = f[base + q];
        rho += fq;
        mx += fq * CX_LOCAL[q] as f64;
        my += fq * CY_LOCAL[q] as f64;
    }
    let ux = if rho > 1e-30 { mx / rho } else { 0.0 };
    let uy = if rho > 1e-30 { my / rho } else { 0.0 };
    (rho, ux, uy)
}

/// Compute the non-equilibrium stress magnitude from f and feq.
///
/// Returns |Pi_neq| = sqrt(Pi_xx^2 + 2*Pi_xy^2 + Pi_yy^2)
#[inline]
fn pi_neq_magnitude(f: &[f64], base: usize, rho: f64, ux: f64, uy: f64) -> f64 {
    let mut pi_xx = 0.0_f64;
    let mut pi_xy = 0.0_f64;
    let mut pi_yy = 0.0_f64;
    for q in 0..9 {
        let cx = CX_LOCAL[q] as f64;
        let cy = CY_LOCAL[q] as f64;
        let feq = equilibrium_local(rho, ux, uy, q);
        let fneq = f[base + q] - feq;
        pi_xx += fneq * cx * cx;
        pi_xy += fneq * cx * cy;
        pi_yy += fneq * cy * cy;
    }
    (pi_xx * pi_xx + 2.0 * pi_xy * pi_xy + pi_yy * pi_yy).sqrt()
}

/// Combined collision and streaming with WMLES closure for a channel flow.
///
/// - Skips wall nodes (j == 0 or j == ny-1): they are treated as solid.
/// - Near-wall nodes (j == 1 or j == ny-2): eddy viscosity from [`WallModeledLes`].
/// - Interior nodes: Smagorinsky SGS closure with `coupling.cs_smag`.
/// - Streaming: bounce-back at walls, periodic in x.
///
/// # Arguments
/// * `channel`  - Mutable reference to the channel flow state.
/// * `coupling` - WMLES coupling parameters.
pub fn collide_stream_wmles(channel: &mut ChannelFlow, coupling: &WmlesChannelCoupling) {
    let nx = channel.nx;
    let ny = channel.ny;
    let nu = coupling.nu;
    let cs = coupling.cs_smag;
    let forcing = channel.forcing;

    // cs^2 / (2*cs4): Smagorinsky coefficient for tau computation
    // tau_eff = 0.5*(tau + sqrt(tau^2 + 18*cs^2*|Pi_neq|/(rho*cs4)))
    let cs4 = 1.0 / 9.0; // cs^2 in LBM = 1/3; cs4 = (1/3)^2 = 1/9

    // Collision step
    for i in 0..nx {
        for j in 0..ny {
            let idx = i * ny + j;
            let base = idx * 9;
            let (rho, mut ux, uy) = moments_local(&channel.f_dist, base);
            if j > 0 && j < ny - 1 {
                ux += forcing;
            }
            channel.rho[idx] = rho;
            channel.u[idx] = [ux, uy];

            // Skip wall nodes
            if j == 0 || j == ny - 1 {
                continue;
            }

            let omega = if j == 1 || j == ny - 2 {
                // Near-wall: use WMLES eddy viscosity
                let y_match = j.min(ny - 1 - j) as f64;
                let u_match = ux.abs();
                let nu_t = coupling
                    .wall_model
                    .eddy_viscosity(y_match.max(1e-30), u_match);
                let nu_eff = nu + nu_t;
                let tau_wall = 3.0 * nu_eff + 0.5;
                1.0 / tau_wall.max(0.5 + 1e-10)
            } else {
                // Interior: Smagorinsky SGS
                let tau_base = 3.0 * nu + 0.5;
                let pi_mag = pi_neq_magnitude(&channel.f_dist, base, rho, ux, uy);
                let tau_eff = 0.5
                    * (tau_base
                        + (tau_base * tau_base + 18.0 * cs * cs * pi_mag / (rho * cs4)).sqrt());
                1.0 / tau_eff.max(0.5 + 1e-10)
            };

            for q in 0..9 {
                let feq = equilibrium_local(rho, ux, uy, q);
                channel.f_dist[base + q] += omega * (feq - channel.f_dist[base + q]);
            }
        }
    }

    // Streaming step
    let mut f_new = vec![0.0_f64; nx * ny * 9];
    for i in 0..nx {
        for j in 0..ny {
            let base = (i * ny + j) * 9;
            for q in 0..9 {
                let nj = j as i32 + CY_LOCAL[q];
                if nj < 0 || nj >= ny as i32 {
                    // Bounce-back at walls
                    f_new[(i * ny + j) * 9 + BB_LOCAL[q]] = channel.f_dist[base + q];
                } else {
                    let ni = (i as i32 + CX_LOCAL[q]).rem_euclid(nx as i32) as usize;
                    let dst = (ni * ny + nj as usize) * 9 + q;
                    f_new[dst] = channel.f_dist[base + q];
                }
            }
        }
    }
    channel.f_dist = f_new;
}

/// Compute the wall-normal eddy viscosity profile using the WMLES model.
///
/// For each wall-normal index j, averages the streamwise velocity over i,
/// then queries the WMLES model for the eddy viscosity.
///
/// Returns a vector of `(y, nu_t)` pairs where y is the lattice index.
///
/// # Arguments
/// * `channel`  - Reference to the channel flow state.
/// * `coupling` - WMLES coupling parameters.
pub fn wmles_eddy_viscosity_profile(
    channel: &ChannelFlow,
    coupling: &WmlesChannelCoupling,
) -> Vec<(f64, f64)> {
    let nx = channel.nx;
    let ny = channel.ny;
    (0..ny)
        .map(|j| {
            let y_match = j.min(ny - 1 - j) as f64;
            let mean_u: f64 = (0..nx).map(|i| channel.u[i * ny + j][0]).sum::<f64>() / nx as f64;
            let nu_t = if j == 0 || j == ny - 1 {
                0.0
            } else {
                coupling
                    .wall_model
                    .eddy_viscosity(y_match.max(1e-30), mean_u.abs())
            };
            (j as f64, nu_t)
        })
        .collect()
}
