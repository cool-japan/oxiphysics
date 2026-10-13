// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Weakly compressible SPH force evaluation for
//! [`WasmSphSimulation`](super::WasmSphSimulation).
//!
//! All sums run over the particles inside the compact support radius `h`
//! (`WasmSphConfig::smoothing_length`) of the cubic-spline kernel
//! [`cubic_spline_compact`], found with a uniform cell list ([`NeighborGrid`]),
//! so a force pass costs `O(n * k)` for `k` neighbours instead of `O(n^2)`.
//!
//! # Formulation
//!
//! With densities `rho` from SPH summation and the Tait equation of state
//! `P = B max(0, (rho / rho0)^gamma - 1)` (pressure is clamped at zero: a
//! sparsely populated free surface or wall layer is under-dense, and a
//! negative pressure would pull neighbours together, the tensile instability):
//!
//! * pressure: `a_i = -sum_j m_j (P_i / rho_i^2 + P_j / rho_j^2) grad_i W_ij`,
//!   which is antisymmetric in `i <-> j` and therefore conserves momentum;
//! * viscosity (Morris et al. 1997):
//!   `a_i = sum_j 2 mu m_j (r_ij . grad_i W_ij) / (rho_i rho_j (r_ij^2 + 0.01 h^2)) (v_i - v_j)`;
//! * surface tension (continuum surface force, Brackbill et al. 1992):
//!   with the colour field `c_i = sum_j (m_j / rho_j) W_ij` and the interface
//!   normal `n_i = grad c_i = sum_j (m_j / rho_j) grad_i W_ij` (it points into
//!   the fluid), particles with `|n_i|` above a small threshold lie on an
//!   interface.  Their curvature is `kappa_i = -div(n / |n|)_i
//!   = -sum_j (m_j / rho_j) (nhat_j - nhat_i) . grad_i W_ij` over the
//!   interface neighbours `j` (it vanishes on a flat interface) and the
//!   acceleration is `a_i = gamma kappa_i n_i / rho_i`;
//! * XSPH smoothing: `dv_i = eps sum_j (m_j / rho_ij) (v_j - v_i) W_ij` with
//!   `rho_ij = (rho_i + rho_j) / 2`, evaluated from a velocity snapshot;
//! * boundary contact: a fluid particle that comes closer than
//!   `particle_radius` to a boundary particle is pushed back out and the
//!   velocity component into the wall is reflected with
//!   `boundary_restitution`.
//!
//! Boundary particles take part in every sum (they fill the density and
//! pressure fields near walls) but are never moved.

use std::collections::HashMap;
use std::f64::consts::PI;

use super::{WasmSphConfig, WasmSphParticle, cubic_spline_compact};

/// Plain 3-vector.
type Vec3 = [f64; 3];

/// Interface-normal magnitude `|grad c|`, in units of `1 / h`, above which a
/// particle counts as lying on an interface for the surface-tension force.
/// On a regularly packed fluid the outermost particle layer gives roughly
/// `1 / h` to `1.7 / h`, the layer below it up to about `0.3 / h`, and the
/// bulk gives zero, so `0.2 / h` selects the two outer layers.
const SURFACE_NORMAL_THRESHOLD: f64 = 0.2;

#[inline]
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `acc += s * v`.
#[inline]
pub(super) fn axpy(acc: &mut Vec3, s: f64, v: Vec3) {
    acc[0] += s * v[0];
    acc[1] += s * v[1];
    acc[2] += s * v[2];
}

/// `true` when the squared distance `dist2` lies strictly inside `limit2`;
/// `false` for NaN, so degenerate coordinates are treated as "not a neighbour".
#[inline]
fn is_within(dist2: f64, limit2: f64) -> bool {
    dist2 < limit2
}

/// `1 / density`, or `0` for a particle without a usable density (zero mass
/// or non-finite values) so that it neither exerts nor feels SPH forces
/// instead of poisoning the simulation with NaN.
#[inline]
fn inverse_density(density: f64) -> f64 {
    if density.is_finite() && density > 0.0 {
        1.0 / density
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Kernel derivatives
// ---------------------------------------------------------------------------

/// Gradient `grad_i W(x_i - x_j, h)` of [`cubic_spline_compact`] with respect
/// to the position of particle `i`.
///
/// `d = x_i - x_j` and `r = |d|`.  The gradient vanishes at `r = 0` (by
/// symmetry) and outside the support `r >= h`.
pub(super) fn cubic_spline_gradient(d: Vec3, r: f64, h: f64) -> Vec3 {
    if !(r > 0.0 && r < h) {
        return [0.0; 3];
    }
    let sigma = 8.0 / (PI * h * h * h);
    let q = r / h;
    let dw_dq = if q <= 0.5 {
        sigma * (18.0 * q * q - 12.0 * q)
    } else {
        let t = 1.0 - q;
        -6.0 * sigma * t * t
    };
    let s = dw_dq / (h * r);
    [d[0] * s, d[1] * s, d[2] * s]
}

// ---------------------------------------------------------------------------
// Equation of state
// ---------------------------------------------------------------------------

/// `true` when `config` describes a usable kernel and equation of state
/// (`smoothing_length` and `rest_density` finite and positive).  Otherwise
/// the particles interact with nothing but gravity.
pub(super) fn interaction_enabled(config: &WasmSphConfig) -> bool {
    config.smoothing_length.is_finite()
        && config.smoothing_length > 0.0
        && config.rest_density.is_finite()
        && config.rest_density > 0.0
}

/// Tait equation of state `P = B ((rho / rho0)^gamma - 1)`, clamped at zero.
pub(super) fn tait_pressure(density: f64, config: &WasmSphConfig) -> f64 {
    let ratio = density / config.rest_density;
    (config.stiffness * (ratio.powf(config.tait_exponent) - 1.0)).max(0.0)
}

// ---------------------------------------------------------------------------
// Neighbour search
// ---------------------------------------------------------------------------

/// Uniform cell list over particle positions.
///
/// The cell edge equals the kernel support radius, so every particle within
/// that radius of a point lies in the point's own cell or one of the 26
/// adjacent cells.  Candidates are visited in a fixed order (cell offsets in
/// order, particles by ascending index), which keeps every floating-point
/// sum deterministic.
pub(super) struct NeighborGrid {
    inv_cell: f64,
    cells: HashMap<[i64; 3], Vec<usize>>,
}

impl NeighborGrid {
    /// Cell list of all particles for support radius `support`.
    pub(super) fn build(particles: &[WasmSphParticle], support: f64) -> Self {
        Self::build_filtered(particles, support, |_| true)
    }

    /// Cell list of the particles accepted by `keep`.
    pub(super) fn build_filtered(
        particles: &[WasmSphParticle],
        support: f64,
        keep: impl Fn(&WasmSphParticle) -> bool,
    ) -> Self {
        let inv_cell = 1.0 / support;
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, p) in particles.iter().enumerate() {
            if keep(p) {
                cells
                    .entry(Self::cell_of(p.position, inv_cell))
                    .or_default()
                    .push(i);
            }
        }
        NeighborGrid { inv_cell, cells }
    }

    /// Integer cell coordinates of `pos`.  The float-to-integer cast
    /// saturates and maps NaN to zero, so absurd positions cannot panic.
    fn cell_of(pos: Vec3, inv_cell: f64) -> [i64; 3] {
        [
            (pos[0] * inv_cell).floor() as i64,
            (pos[1] * inv_cell).floor() as i64,
            (pos[2] * inv_cell).floor() as i64,
        ]
    }

    /// `true` when no particle was inserted.
    pub(super) fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    /// Call `visit(j)` for every stored particle in the 27 cells around
    /// `pos`.  This includes particles farther than the support radius; the
    /// caller applies the exact distance cutoff.
    pub(super) fn for_each_candidate(&self, pos: Vec3, mut visit: impl FnMut(usize)) {
        let c = Self::cell_of(pos, self.inv_cell);
        for dz in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dx in -1..=1_i64 {
                    let key = [
                        c[0].saturating_add(dx),
                        c[1].saturating_add(dy),
                        c[2].saturating_add(dz),
                    ];
                    if let Some(list) = self.cells.get(&key) {
                        for &j in list {
                            visit(j);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Densities and forces
// ---------------------------------------------------------------------------

/// SPH density `rho_i = sum_j m_j W(|x_i - x_j|, h)` of every particle, fluid
/// and boundary, summed over the neighbours found by `grid`.
///
/// `grid` must index the current positions of `particles`.  Particles with
/// non-finite coordinates are never neighbours of anything (their distance
/// comparison fails), so they cannot poison the densities of the others.
pub(super) fn compute_densities(
    particles: &[WasmSphParticle],
    h: f64,
    grid: &NeighborGrid,
) -> Vec<f64> {
    let h2 = h * h;
    particles
        .iter()
        .map(|pi| {
            let mut rho = 0.0_f64;
            grid.for_each_candidate(pi.position, |j| {
                let pj = &particles[j];
                let d = sub(pi.position, pj.position);
                let r2 = dot(d, d);
                if is_within(r2, h2) {
                    rho += pj.mass * cubic_spline_compact(r2.sqrt(), h);
                }
            });
            rho
        })
        .collect()
}

/// Result of one force evaluation.
pub(super) struct FluidForces {
    /// Total acceleration of every fluid particle (gravity, pressure,
    /// viscosity and surface tension); zero for boundary particles.
    pub(super) acceleration: Vec<Vec3>,
    /// Smoothed colour field `c_i = sum_j (m_j / rho_j) W_ij` of every
    /// particle (about `1` inside the fluid; its gradient marks the surface).
    pub(super) colour: Vec<f64>,
}

/// Colour field and interface normals of all particles.
struct InterfaceFields {
    /// Colour field `c_i`.
    colour: Vec<f64>,
    /// `grad c_i`, pointing into the fluid (empty unless surface tension is on).
    gradient: Vec<Vec3>,
    /// Unit normal of the particles that lie on an interface (`|grad c|`
    /// above [`SURFACE_NORMAL_THRESHOLD`]); `None` for all others (empty
    /// unless surface tension is on).
    unit_normal: Vec<Option<Vec3>>,
}

/// Colour field of every particle and, when `with_normals` is set, the
/// interface normals needed by the surface-tension force.
fn interface_fields(
    particles: &[WasmSphParticle],
    config: &WasmSphConfig,
    grid: &NeighborGrid,
    with_normals: bool,
) -> InterfaceFields {
    let h = config.smoothing_length;
    let h2 = h * h;
    let threshold = SURFACE_NORMAL_THRESHOLD / h;
    let mut colour = Vec::with_capacity(particles.len());
    let mut gradient = Vec::new();
    let mut unit_normal = Vec::new();
    for pi in particles {
        let mut colour_i = 0.0_f64;
        let mut normal = [0.0_f64; 3];
        grid.for_each_candidate(pi.position, |j| {
            let pj = &particles[j];
            let d = sub(pi.position, pj.position);
            let r2 = dot(d, d);
            if !is_within(r2, h2) {
                return;
            }
            let r = r2.sqrt();
            let volume_j = pj.mass * inverse_density(pj.density);
            colour_i += volume_j * cubic_spline_compact(r, h);
            if with_normals {
                axpy(&mut normal, volume_j, cubic_spline_gradient(d, r, h));
            }
        });
        colour.push(colour_i);
        if with_normals {
            let length = dot(normal, normal).sqrt();
            unit_normal.push(
                (length > threshold)
                    .then(|| [normal[0] / length, normal[1] / length, normal[2] / length]),
            );
            gradient.push(normal);
        }
    }
    InterfaceFields {
        colour,
        gradient,
        unit_normal,
    }
}

/// Evaluate the accelerations of all fluid particles.
///
/// The particles' `density` and `pressure` fields must be up to date, and
/// `grid` must index their current positions.
pub(super) fn compute_forces(
    particles: &[WasmSphParticle],
    config: &WasmSphConfig,
    grid: &NeighborGrid,
) -> FluidForces {
    let h = config.smoothing_length;
    let h2 = h * h;
    let eta2 = 0.01 * h2;
    let two_mu = 2.0 * config.viscosity;
    // A non-finite coefficient disables the force instead of poisoning the run.
    let gamma = if config.surface_tension.is_finite() {
        config.surface_tension
    } else {
        0.0
    };
    let with_surface_tension = gamma != 0.0;
    let interface = interface_fields(particles, config, grid, with_surface_tension);
    let mut acceleration = vec![[0.0_f64; 3]; particles.len()];

    for (i, pi) in particles.iter().enumerate() {
        if pi.is_boundary {
            continue;
        }
        let inv_rho_i = inverse_density(pi.density);
        let p_over_rho2_i = pi.pressure * inv_rho_i * inv_rho_i;
        let unit_normal_i = interface.unit_normal.get(i).copied().flatten();
        let mut a_pressure = [0.0_f64; 3];
        let mut a_viscous = [0.0_f64; 3];
        let mut normal_divergence = 0.0_f64;

        grid.for_each_candidate(pi.position, |j| {
            let pj = &particles[j];
            let d = sub(pi.position, pj.position);
            let r2 = dot(d, d);
            if !is_within(r2, h2) {
                return;
            }
            let r = r2.sqrt();
            let inv_rho_j = inverse_density(pj.density);
            let grad = cubic_spline_gradient(d, r, h);

            // Symmetric pressure force.
            let p_over_rho2_j = pj.pressure * inv_rho_j * inv_rho_j;
            axpy(
                &mut a_pressure,
                -pj.mass * (p_over_rho2_i + p_over_rho2_j),
                grad,
            );

            // Morris viscosity.
            let weight = two_mu * pj.mass * inv_rho_i * inv_rho_j * dot(d, grad) / (r2 + eta2);
            axpy(&mut a_viscous, weight, sub(pi.velocity, pj.velocity));

            // Divergence of the unit interface normal over interface neighbours.
            if let Some(n_i) = unit_normal_i
                && let Some(Some(n_j)) = interface.unit_normal.get(j)
            {
                normal_divergence += pj.mass * inv_rho_j * dot(sub(*n_j, n_i), grad);
            }
        });

        let mut total = config.gravity;
        axpy(&mut total, 1.0, a_pressure);
        axpy(&mut total, 1.0, a_viscous);
        if unit_normal_i.is_some() {
            // Curvature kappa = -div(nhat); force density gamma kappa grad c.
            axpy(
                &mut total,
                -gamma * normal_divergence * inv_rho_i,
                interface.gradient[i],
            );
        }
        acceleration[i] = total;
    }

    FluidForces {
        acceleration,
        colour: interface.colour,
    }
}

/// XSPH velocity corrections `eps sum_j (m_j / rho_ij) (v_j - v_i) W_ij` of
/// every fluid particle (zero for boundary particles).
///
/// They are computed from the velocities currently stored in `particles`
/// (a snapshot: nothing is updated while the sums run), and `grid` must index
/// the current positions.
pub(super) fn xsph_corrections(
    particles: &[WasmSphParticle],
    config: &WasmSphConfig,
    grid: &NeighborGrid,
) -> Vec<Vec3> {
    let h = config.smoothing_length;
    let h2 = h * h;
    let epsilon = config.xsph_epsilon;
    particles
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            let mut correction = [0.0_f64; 3];
            if pi.is_boundary {
                return correction;
            }
            grid.for_each_candidate(pi.position, |j| {
                if j == i {
                    return;
                }
                let pj = &particles[j];
                let d = sub(pi.position, pj.position);
                let r2 = dot(d, d);
                if !is_within(r2, h2) {
                    return;
                }
                let mean_density = 0.5 * (pi.density + pj.density);
                if !(mean_density.is_finite() && mean_density > 0.0) {
                    return;
                }
                let weight = epsilon * pj.mass / mean_density * cubic_spline_compact(r2.sqrt(), h);
                axpy(&mut correction, weight, sub(pj.velocity, pi.velocity));
            });
            correction
        })
        .collect()
}

/// Push fluid particles that penetrated the boundary layer back out.
///
/// A fluid particle closer than `particle_radius` to a boundary particle
/// is moved back to that distance along the line of centres, and the
/// velocity component into the boundary particle is reflected and scaled by
/// `boundary_restitution` (clamped to `[0, 1]`).  Boundary particles do not
/// move.
pub(super) fn resolve_boundary_contacts(particles: &mut [WasmSphParticle], config: &WasmSphConfig) {
    let contact = config.particle_radius.min(config.smoothing_length);
    if !(contact.is_finite() && contact > 0.0) {
        return;
    }
    let restitution = if config.boundary_restitution.is_finite() {
        config.boundary_restitution.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let grid = NeighborGrid::build_filtered(particles, config.smoothing_length, |p| p.is_boundary);
    if grid.is_empty() {
        return;
    }
    let contact2 = contact * contact;
    let mut candidates: Vec<usize> = Vec::new();
    for i in 0..particles.len() {
        if particles[i].is_boundary {
            continue;
        }
        candidates.clear();
        grid.for_each_candidate(particles[i].position, |j| candidates.push(j));
        let mut position = particles[i].position;
        let mut velocity = particles[i].velocity;
        for &j in &candidates {
            let d = sub(position, particles[j].position);
            let r2 = dot(d, d);
            if !is_within(r2, contact2) {
                continue;
            }
            let r = r2.sqrt();
            // Contact normal: from the boundary particle towards the fluid
            // particle; for coincident centres use the reversed velocity (or
            // "up" for a particle at rest).
            let normal = if r > 0.0 {
                [d[0] / r, d[1] / r, d[2] / r]
            } else {
                let speed = dot(velocity, velocity).sqrt();
                if speed > 0.0 {
                    [
                        -velocity[0] / speed,
                        -velocity[1] / speed,
                        -velocity[2] / speed,
                    ]
                } else {
                    [0.0, 1.0, 0.0]
                }
            };
            position = particles[j].position;
            axpy(&mut position, contact, normal);
            let into_wall = dot(velocity, normal);
            if into_wall < 0.0 {
                axpy(&mut velocity, -(1.0 + restitution) * into_wall, normal);
            }
        }
        particles[i].position = position;
        particles[i].velocity = velocity;
    }
}
