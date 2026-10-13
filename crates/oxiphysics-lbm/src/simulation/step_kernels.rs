//! Shared per-cell kernels of the high-level [`LbmSimulation`](super::LbmSimulation)
//! driver.
//!
//! The driver supports D2Q9, D3Q19 and D3Q27 lattices.  All three store their
//! populations as `f[direction][cell]`, so the kernels in this module work on
//! that layout directly and take the lattice description (weights, velocity
//! vectors, opposite directions) from a [`LatticeTable`].  Two-dimensional
//! lattices are embedded in three dimensions with `c_z = 0`, which lets one
//! implementation serve every lattice.
//!
//! # Contents
//!
//! * [`collide_field`] — fused BGK collision with the Guo forcing scheme and
//!   an optional per-cell Smagorinsky relaxation time.
//! * [`smagorinsky_tau`] — local relaxation time of the Smagorinsky model.
//! * [`rescale_field`] — re-expresses the populations in the lattice units of
//!   a different time step (used when the caller changes `dt` between steps).
//! * [`build_wall_links_3d`], [`apply_wall_links`], [`reset_solid_cells`] —
//!   half-way bounce-back walls for 3-D simulations.

use super::types::GuoForceParams;
use crate::lattice::{CS2, Lattice};

/// Largest number of discrete velocities supported (D3Q27).
pub(super) const MAX_Q: usize = 27;

/// Cells whose density does not exceed this value are treated as vacuum and
/// are left untouched by the collision kernel.
const RHO_FLOOR: f64 = 1e-15;

/// Weights, velocity vectors and opposite directions of a lattice, with the
/// velocities stored as floating-point triples (`c_z = 0` for D2Q9).
#[derive(Debug, Clone)]
pub(super) struct LatticeTable {
    /// Number of discrete velocities.
    pub(super) q: usize,
    /// Quadrature weights `w_i`.
    pub(super) weights: Vec<f64>,
    /// Discrete velocities `c_i` as `[cx, cy, cz]`.
    pub(super) velocities: Vec<[f64; 3]>,
    /// Index of the direction opposite to each direction.
    pub(super) opposite: Vec<usize>,
}

impl LatticeTable {
    /// Build the table from a lattice descriptor.
    ///
    /// The opposite-direction table is derived from the velocity vectors
    /// themselves (`c_opposite(i) = -c_i`), so the bounce-back walls cannot
    /// disagree with the velocities that streaming uses.
    pub(super) fn from_lattice(lattice: &Lattice) -> Self {
        let q = lattice.q();
        let dim = lattice.dim();
        let mut weights = Vec::with_capacity(q);
        let mut velocities: Vec<[f64; 3]> = Vec::with_capacity(q);
        for i in 0..q {
            weights.push(lattice.weight(i));
            let c = if dim == 2 {
                let v = lattice.velocity_2d(i);
                [f64::from(v[0]), f64::from(v[1]), 0.0]
            } else {
                let v = lattice.velocity_3d(i);
                [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])]
            };
            velocities.push(c);
        }
        let opposite: Vec<usize> = velocities
            .iter()
            .enumerate()
            .map(|(i, c)| {
                velocities
                    .iter()
                    .position(|d| d[0] == -c[0] && d[1] == -c[1] && d[2] == -c[2])
                    .unwrap_or(i)
            })
            .collect();
        Self {
            q,
            weights,
            velocities,
            opposite,
        }
    }
}

/// Inputs of one collision sweep.
#[derive(Debug, Clone, Copy)]
pub(super) struct CollisionSetup {
    /// Base BGK relaxation frequency `omega = 1 / tau_0`.
    pub(super) omega: f64,
    /// Smagorinsky constant; the model is active when this is `Some(c)` with
    /// `c > 0`.
    pub(super) smagorinsky_cs: Option<f64>,
    /// Body force per unit volume acting during this step, in lattice units
    /// (all zeros when no force is applied).
    pub(super) force: [f64; 3],
}

/// Second-order BGK equilibrium of one direction.
#[inline]
fn equilibrium(w: f64, rho: f64, c: [f64; 3], u: [f64; 3], u_sq: f64) -> f64 {
    let cu = c[0] * u[0] + c[1] * u[1] + c[2] * u[2];
    w * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2))
}

/// Zeroth and first moments `(rho, rho * u)` of one cell's populations.
#[inline]
fn moments(pop: &[f64], velocities: &[[f64; 3]]) -> (f64, [f64; 3]) {
    let mut rho = 0.0;
    let mut mom = [0.0f64; 3];
    for (fi, c) in pop.iter().zip(velocities) {
        rho += fi;
        mom[0] += fi * c[0];
        mom[1] += fi * c[1];
        mom[2] += fi * c[2];
    }
    (rho, mom)
}

/// Relaxation time of the Smagorinsky sub-grid-scale model in lattice units.
///
/// With the filter width equal to the lattice spacing, the eddy viscosity is
/// `nu_t = (C_s)^2 |S|` and the non-equilibrium momentum flux satisfies
/// `Pi_neq = -2 rho c_s^2 tau_eff S`, so that `|S| = ||Pi_neq|| /
/// (sqrt(2) rho c_s^2 tau_eff)` where `||Pi_neq||` is the Frobenius norm.
/// Requiring `tau_eff = tau_0 + C_s^2 |S| / c_s^2` gives the quadratic
/// `tau_eff (tau_eff - tau_0) = C_s^2 ||Pi_neq|| / (sqrt(2) rho c_s^4)` whose
/// positive root is returned.
///
/// # Arguments
/// * `tau0`       – molecular relaxation time `1 / omega`
/// * `cs_smag`    – Smagorinsky constant `C_s`
/// * `pi_norm`    – Frobenius norm of the non-equilibrium momentum flux
/// * `rho`        – local density
pub(super) fn smagorinsky_tau(tau0: f64, cs_smag: f64, pi_norm: f64, rho: f64) -> f64 {
    let discriminant = tau0 * tau0
        + 4.0 * cs_smag * cs_smag * pi_norm / (std::f64::consts::SQRT_2 * rho * CS2 * CS2);
    0.5 * (tau0 + discriminant.sqrt())
}

/// Local relaxation frequency of a cell.
///
/// `feq` must be the equilibrium evaluated with the same (force-shifted)
/// velocity that the collision uses.
fn local_omega(
    pop: &[f64],
    feq: &[f64],
    velocities: &[[f64; 3]],
    setup: &CollisionSetup,
    rho: f64,
) -> f64 {
    let Some(cs) = setup.smagorinsky_cs.filter(|&c| c > 0.0) else {
        return setup.omega;
    };
    let (mut xx, mut yy, mut zz) = (0.0, 0.0, 0.0);
    let (mut xy, mut xz, mut yz) = (0.0, 0.0, 0.0);
    for ((fi, fe), c) in pop.iter().zip(feq).zip(velocities) {
        let fneq = fi - fe;
        xx += fneq * c[0] * c[0];
        yy += fneq * c[1] * c[1];
        zz += fneq * c[2] * c[2];
        xy += fneq * c[0] * c[1];
        xz += fneq * c[0] * c[2];
        yz += fneq * c[1] * c[2];
    }
    let pi_norm = (xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + xz * xz + yz * yz)).sqrt();
    1.0 / smagorinsky_tau(1.0 / setup.omega, cs, pi_norm, rho)
}

/// BGK collision of one cell with Guo forcing.
///
/// The density is `rho = sum f_i`; the fluid velocity is the Guo velocity
/// `u = (sum f_i c_i + F/2) / rho`, which is used both for the equilibrium and
/// for the forcing term.  The post-collision populations are
/// `f_i + omega (feq_i - f_i) + F_i` with
/// `F_i = w_i (1 - omega/2) [(c_i - u)/c_s^2 + (c_i . u) c_i / c_s^4] . F`,
/// so the collision conserves mass and adds exactly `F` to the momentum.
/// When Smagorinsky is active, `omega` is the local frequency of the cell.
pub(super) fn collide_cell(pop: &mut [f64], table: &LatticeTable, setup: &CollisionSetup) {
    let q = table.q;
    let pop = &mut pop[..q];
    let (rho, mom) = moments(pop, &table.velocities);
    if rho.is_nan() || rho <= RHO_FLOOR {
        return;
    }
    let force = setup.force;
    let inv_rho = 1.0 / rho;
    let u = [
        (mom[0] + 0.5 * force[0]) * inv_rho,
        (mom[1] + 0.5 * force[1]) * inv_rho,
        (mom[2] + 0.5 * force[2]) * inv_rho,
    ];
    let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    let mut feq = [0.0f64; MAX_Q];
    for ((e, &w), &c) in feq.iter_mut().zip(&table.weights).zip(&table.velocities) {
        *e = equilibrium(w, rho, c, u, u_sq);
    }
    let omega = local_omega(pop, &feq[..q], &table.velocities, setup, rho);
    let guo = GuoForceParams::new(force, omega);
    for (((fi, &fe), &w), &c) in pop
        .iter_mut()
        .zip(&feq)
        .zip(&table.weights)
        .zip(&table.velocities)
    {
        *fi += omega * (fe - *fi) + guo.force_term(w, c, u);
    }
}

/// Collide every non-solid cell of a field stored as `f[direction][cell]`.
pub(super) fn collide_field(
    f: &mut [Vec<f64>],
    n_cells: usize,
    table: &LatticeTable,
    setup: &CollisionSetup,
    solid: Option<&[bool]>,
) {
    let q = table.q;
    let mut pop = [0.0f64; MAX_Q];
    for k in 0..n_cells {
        if is_solid(solid, k) {
            continue;
        }
        for (p, plane) in pop.iter_mut().zip(f.iter()) {
            *p = plane[k];
        }
        collide_cell(&mut pop[..q], table, setup);
        for (p, plane) in pop.iter().zip(f.iter_mut()) {
            plane[k] = *p;
        }
    }
}

/// Whether cell `k` is masked as solid.
#[inline]
pub(super) fn is_solid(solid: Option<&[bool]>, k: usize) -> bool {
    solid.is_some_and(|mask| mask.get(k).copied().unwrap_or(false))
}

/// Add the second half of the force impulse to a macroscopic velocity field.
///
/// The Guo scheme defines the physical velocity as
/// `u = (sum f_i c_i + F/2) / rho`; `compute_macroscopic` only evaluates the
/// plain momentum moment, so the fields reported by the driver are corrected
/// with this routine.  Solid cells are skipped.
pub(super) fn add_half_force_velocity(
    rho: &[f64],
    velocity: &mut [&mut [f64]],
    force: [f64; 3],
    solid: Option<&[bool]>,
) {
    if force.iter().all(|&component| component == 0.0) {
        return;
    }
    for (k, &density) in rho.iter().enumerate() {
        if is_solid(solid, k) || density.abs() <= RHO_FLOOR {
            continue;
        }
        let inv_rho = 0.5 / density;
        for (component, field) in velocity.iter_mut().enumerate() {
            field[k] += force[component] * inv_rho;
        }
    }
}

/// Re-express a field of pre-collision populations in the lattice units of a
/// different time step.
///
/// When the physical time covered by one lattice step changes by the factor
/// `ratio = dt_new / dt_old` while the lattice spacing is fixed, the same
/// physical flow has the following representation:
///
/// * lattice velocities scale with `ratio`,
/// * lattice viscosities scale with `ratio` (diffusive scaling),
/// * lattice forces scale with `ratio^2`,
/// * the density deviation from its mean (the pressure) scales with
///   `ratio^2`, being of order `Ma^2`,
/// * the non-equilibrium part of the populations, which is proportional to
///   `tau * du/dx`, scales with `neq_factor = (tau_new / tau_old) * ratio`.
///
/// `force_prev` is the lattice body force that acted during the previous step;
/// it enters through the Guo half-force shift of the velocity.
pub(super) fn rescale_field(
    f: &mut [Vec<f64>],
    n_cells: usize,
    table: &LatticeTable,
    ratio: f64,
    neq_factor: f64,
    force_prev: [f64; 3],
    solid: Option<&[bool]>,
) {
    let q = table.q;
    let mut pop = [0.0f64; MAX_Q];

    // Mean density of the fluid cells: the reference of the pressure scaling.
    let mut mass = 0.0;
    let mut fluid_cells = 0usize;
    for k in 0..n_cells {
        if is_solid(solid, k) {
            continue;
        }
        mass += f.iter().map(|plane| plane[k]).sum::<f64>();
        fluid_cells += 1;
    }
    if fluid_cells == 0 {
        return;
    }
    let rho_mean = mass / fluid_cells as f64;

    let force_new_scale = ratio * ratio;
    for k in 0..n_cells {
        if is_solid(solid, k) {
            continue;
        }
        for (p, plane) in pop.iter_mut().zip(f.iter()) {
            *p = plane[k];
        }
        let (rho, mom) = moments(&pop[..q], &table.velocities);
        if rho.is_nan() || rho <= RHO_FLOOR {
            continue;
        }
        let rho_new = rho_mean + ratio * ratio * (rho - rho_mean);
        // Physical velocity before the change, scaled, minus the new half-force.
        let mut mom_new = [0.0f64; 3];
        for d in 0..3 {
            let u_phys_old = (mom[d] + 0.5 * force_prev[d]) / rho;
            mom_new[d] = rho_new * ratio * u_phys_old - 0.5 * force_new_scale * force_prev[d];
        }
        let u_old = [mom[0] / rho, mom[1] / rho, mom[2] / rho];
        let u_new = [
            mom_new[0] / rho_new,
            mom_new[1] / rho_new,
            mom_new[2] / rho_new,
        ];
        let u_old_sq = u_old[0] * u_old[0] + u_old[1] * u_old[1] + u_old[2] * u_old[2];
        let u_new_sq = u_new[0] * u_new[0] + u_new[1] * u_new[1] + u_new[2] * u_new[2];
        for ((p, &w), &c) in pop.iter_mut().zip(&table.weights).zip(&table.velocities) {
            let fneq = *p - equilibrium(w, rho, c, u_old, u_old_sq);
            *p = equilibrium(w, rho_new, c, u_new, u_new_sq) + neq_factor * fneq;
        }
        for (p, plane) in pop.iter().zip(f.iter_mut()) {
            plane[k] = *p;
        }
    }
}

/// One link between a fluid cell and a solid neighbour of a 3-D simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct WallLink {
    /// Flat index of the fluid cell.
    pub(super) fluid: usize,
    /// Flat index of the solid neighbour, `fluid + c_dir` (periodic).
    pub(super) solid: usize,
    /// Direction `dir` that points from the fluid cell into the solid cell.
    pub(super) dir: usize,
}

/// Collect every fluid-to-solid link of a periodic `nx * ny * nz` domain.
pub(super) fn build_wall_links_3d(
    solid: &[bool],
    dims: [usize; 3],
    table: &LatticeTable,
) -> Vec<WallLink> {
    let [nx, ny, nz] = dims;
    let mut links = Vec::new();
    if nx == 0 || ny == 0 || nz == 0 {
        return links;
    }
    let wrap = |coordinate: usize, step: f64, extent: usize| -> usize {
        // Lattice velocities are -1, 0 or +1.
        let shifted = coordinate as isize + step as isize;
        shifted.rem_euclid(extent as isize) as usize
    };
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let fluid = z * ny * nx + y * nx + x;
                if solid.get(fluid).copied().unwrap_or(false) {
                    continue;
                }
                for (dir, c) in table.velocities.iter().enumerate() {
                    if c[0] == 0.0 && c[1] == 0.0 && c[2] == 0.0 {
                        continue;
                    }
                    let xn = wrap(x, c[0], nx);
                    let yn = wrap(y, c[1], ny);
                    let zn = wrap(z, c[2], nz);
                    let neighbour = zn * ny * nx + yn * nx + xn;
                    if solid.get(neighbour).copied().unwrap_or(false) {
                        links.push(WallLink {
                            fluid,
                            solid: neighbour,
                            dir,
                        });
                    }
                }
            }
        }
    }
    links
}

/// Half-way bounce-back, to be applied right after streaming.
///
/// Streaming has moved the post-collision population of direction `dir` of a
/// fluid cell into its solid neighbour.  The wall reflects it: it is read back
/// from the solid cell and becomes the population of the opposite direction of
/// the fluid cell, replacing the meaningless value that streaming brought in
/// from the solid side.  The wall therefore sits half-way between the two
/// cells and the fluid mass is conserved exactly.
pub(super) fn apply_wall_links(f: &mut [Vec<f64>], links: &[WallLink], opposite: &[usize]) {
    for link in links {
        let reflected = f[link.dir][link.solid];
        f[opposite[link.dir]][link.fluid] = reflected;
    }
}

/// Put every solid cell back to the rest state `f_i = w_i` (density 1,
/// velocity 0) so that the values that streaming moved through the walls never
/// accumulate and the macroscopic fields of solid cells stay well defined.
pub(super) fn reset_solid_cells(f: &mut [Vec<f64>], solid: &[bool], table: &LatticeTable) {
    for (k, &masked) in solid.iter().enumerate() {
        if !masked {
            continue;
        }
        for (plane, &w) in f.iter_mut().zip(&table.weights) {
            plane[k] = w;
        }
    }
}
