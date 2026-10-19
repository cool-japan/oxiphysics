// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Guo body-force scheme shared by the D2Q9 and D3Q19 simulations.
//!
//! A uniform body acceleration `a` (force per unit mass, lattice units, time
//! step 1) acts on the fluid with the force density `F = rho a`. BGK collision
//! alone cannot inject that force: shifting the equilibrium velocity by
//! `a / (2 omega)` adds only `F / 2` of momentum per step. The scheme of Guo,
//! Zheng and Shi (Phys. Rev. E 65, 046308, 2002) adds exactly `F` per step:
//!
//! * the physical velocity of a cell is `u = (sum_q e_q f_q + F / 2) / rho`,
//!   that is `sum_q e_q f_q / rho + a / 2`; the equilibrium is evaluated at
//!   this `u`, and so is every velocity the API reports;
//! * collision relaxes towards that equilibrium and adds the source term
//!   `S_q = (1 - omega / 2) w_q [(e_q - u) / c_s^2 + (e_q . u) e_q / c_s^4] . F`
//!   with `c_s^2 = 1/3`.
//!
//! The source term carries no mass (`sum_q S_q = 0`) and the momentum
//! `sum_q e_q S_q = (1 - omega / 2) F`, which together with the shift of `u`
//! in the equilibrium gives a momentum gain of `F` per step for every `omega`.
//! The populations of a cell therefore hold the momentum `rho (u - a / 2)`.

/// Densities at or below this value mark a cell as vacuum.
const VACUUM_DENSITY: f64 = 1e-15;

/// Physical velocity `u = j / rho + a / 2` of a cell.
///
/// `momentum` is `j = sum_q e_q f_q` of the cell's populations and `accel` the
/// body acceleration `a`. A vacuum cell (`rho <= 1e-15`) has no velocity and
/// yields zeros.
#[inline]
pub(super) fn forced_velocity<const D: usize>(
    momentum: [f64; D],
    rho: f64,
    accel: [f64; D],
) -> [f64; D] {
    if rho <= VACUUM_DENSITY {
        return [0.0; D];
    }
    std::array::from_fn(|k| momentum[k] / rho + 0.5 * accel[k])
}

/// Velocity `u - a / 2` of the equilibrium that stores the physical velocity
/// `u` in the populations, so that [`forced_velocity`] reads `u` back.
#[inline]
pub(super) fn stored_velocity<const D: usize>(velocity: [f64; D], accel: [f64; D]) -> [f64; D] {
    std::array::from_fn(|k| velocity[k] - 0.5 * accel[k])
}

/// Guo source term `S_q` of lattice direction `e_q` with weight `w_q`.
///
/// `u` is the physical velocity of the cell (see [`forced_velocity`]) and
/// `force` the force density `F = rho a`; `omega` is the BGK relaxation rate.
#[inline]
pub(super) fn guo_source<const D: usize>(
    omega: f64,
    weight: f64,
    e_q: [f64; D],
    u: [f64; D],
    force: [f64; D],
) -> f64 {
    let e_dot_u: f64 = e_q.iter().zip(&u).map(|(e, v)| e * v).sum();
    let projected: f64 = e_q
        .iter()
        .zip(&u)
        .zip(&force)
        .map(|((e, v), f)| (3.0 * (e - v) + 9.0 * e_dot_u * e) * f)
        .sum();
    (1.0 - 0.5 * omega) * weight * projected
}
