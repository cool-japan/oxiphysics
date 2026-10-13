//! Weakly-compressible colour-function SPH operators for the 2-D
//! [`ColorParticle`] system advanced by
//! [`MultiphaseSphSimulation`](super::types::MultiphaseSphSimulation).
//!
//! Every operator is a neighbour-list sum over the cubic-spline support radius
//! `2h` (see [`support_neighbors`]), so a full step costs `O(N · k)` for `k`
//! neighbours per particle.  The formulation follows the usual air–water SPH
//! recipe:
//!
//! * **Density** by the particle-number summation `ρᵢ = mᵢ Σⱼ W(|xᵢ − xⱼ|, h)`
//!   (Hu & Adams 2006).  Using only the particle's *own* mass keeps the
//!   estimate free of the spurious density jump a mass-weighted summation
//!   suffers across a 1000 : 1 liquid–gas interface.
//! * **Pressure** from the Tait equation of state with the colour-interpolated
//!   rest density `ρ₀` of the particle ([`ColorParticle::interpolated_density`]),
//!   clamped to be non-negative (no tension).
//! * **Pressure acceleration** in the momentum-conserving two-phase form
//!   `aᵢ = −Σⱼ mⱼ (Pᵢ + Pⱼ)/(ρᵢ ρⱼ) ∇ᵢWᵢⱼ` (Colagrossi & Landrini 2003), which
//!   keeps the interface pressure continuous for large density ratios.
//! * **Surface tension** as the continuum surface force of Brackbill et al.
//!   (1992) in the SPH form of Morris (2000): colour gradient `∇Cᵢ`, unit normal
//!   `n̂ᵢ = ∇Cᵢ/|∇Cᵢ|` where `|∇Cᵢ| > 0.01/h`, curvature `κᵢ = −∇·n̂ᵢ` and
//!   acceleration `aᵢ = σ κᵢ ∇Cᵢ / ρᵢ`.

use super::functions::{cubic_spline_gradient, cubic_spline_kernel};
use super::types::ColorParticle;
use crate::neighbor::SpatialHash3D;

/// Relative threshold on `|∇C|` (in units of `1/h`) below which a particle is
/// regarded as lying away from every interface and receives the zero normal.
const NORMAL_THRESHOLD: f64 = 0.01;

/// Exponent `γ` of the Tait equation of state.
pub(super) const TAIT_GAMMA: f64 = 7.0;

/// Neighbour lists of `particles` within the cubic-spline support radius `2h`
/// (self excluded, indices sorted).
///
/// A non-positive or non-finite `h` yields empty lists.
pub(super) fn support_neighbors(particles: &[ColorParticle], h: f64) -> Vec<Vec<usize>> {
    let radius = 2.0 * h;
    if !(radius > 0.0 && radius.is_finite()) {
        return vec![Vec::new(); particles.len()];
    }
    let positions: Vec<[f64; 3]> = particles
        .iter()
        .map(|p| [p.pos[0], p.pos[1], 0.0])
        .collect();
    SpatialHash3D::all_neighbors(&positions, radius)
}

/// `true` for strictly positive (hence non-NaN) values.
#[inline]
fn is_positive(x: f64) -> bool {
    x > 0.0
}

/// Separation `x_i − x_j` and its length for two particles.
#[inline]
fn separation(pi: &ColorParticle, pj: &ColorParticle) -> (f64, f64, f64) {
    let dx = pi.pos[0] - pj.pos[0];
    let dy = pi.pos[1] - pj.pos[1];
    (dx, dy, (dx * dx + dy * dy).sqrt())
}

/// Particle-number density summation `ρᵢ = mᵢ Σⱼ W(|xᵢ − xⱼ|, h)` (the particle
/// itself included), stored in [`ColorParticle::density`].
pub(super) fn number_density_summation(
    particles: &mut [ColorParticle],
    h: f64,
    neighbors: &[Vec<usize>],
) {
    let self_weight = cubic_spline_kernel(0.0, h);
    let number_densities: Vec<f64> = particles
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            let mut sigma = self_weight;
            if let Some(list) = neighbors.get(i) {
                for &j in list {
                    if j == i {
                        continue;
                    }
                    let Some(pj) = particles.get(j) else {
                        continue;
                    };
                    let (_, _, r) = separation(pi, pj);
                    sigma += cubic_spline_kernel(r, h);
                }
            }
            sigma
        })
        .collect();
    for (p, sigma) in particles.iter_mut().zip(number_densities) {
        p.density = p.mass * sigma;
    }
}

/// Tait pressure `P = B ((ρ/ρ₀)^γ − 1)`, `B = ρ₀ c₀²/γ`, with the
/// colour-interpolated rest density `ρ₀` of each particle, clamped to `P ≥ 0`.
///
/// Particles with a non-positive rest density, density or sound speed get
/// `P = 0`.
pub(super) fn tait_pressures(
    particles: &mut [ColorParticle],
    rho_l: f64,
    rho_g: f64,
    sound_speed: f64,
) {
    for p in particles.iter_mut() {
        let rho0 = p.interpolated_density(rho_l, rho_g);
        p.pressure = if rho0 > 0.0 && p.density > 0.0 && sound_speed > 0.0 {
            let stiffness = rho0 * sound_speed * sound_speed / TAIT_GAMMA;
            (stiffness * ((p.density / rho0).powf(TAIT_GAMMA) - 1.0)).max(0.0)
        } else {
            0.0
        };
    }
}

/// Curvature `κᵢ = −∇·n̂ᵢ` of every particle from the stored colour gradients
/// ([`ColorParticle::color_grad`], computed by `compute_color_gradients`),
/// stored in [`ColorParticle::curvature`].
///
/// Following Morris (2000), the unit normal `n̂ = ∇C/|∇C|` exists only where
/// `|∇C| > 0.01/h`, and the divergence is taken over the neighbours that have a
/// normal, normalised by the SPH estimate of unity over that same set:
///
/// `κᵢ = −[Σⱼ Vⱼ (n̂ⱼ − n̂ᵢ)·∇ᵢWᵢⱼ] / [Vᵢ W(0) + Σⱼ Vⱼ Wᵢⱼ]`, `Vⱼ = mⱼ/ρⱼ`.
///
/// Letting the normal-free particles outside the interface band count as
/// "zero normals" would instead produce a large spurious curvature at the edge
/// of the band, even for a flat interface.  Particles without a normal get
/// `κ = 0`.  `κ` is positive where the `C = 1` phase is convex (`1/R` for a disc
/// of radius `R`) and vanishes for a straight interface.
pub(super) fn curvatures_from_gradients(
    particles: &mut [ColorParticle],
    h: f64,
    neighbors: &[Vec<usize>],
) {
    let threshold = NORMAL_THRESHOLD / h;
    let normals: Vec<Option<[f64; 2]>> = particles
        .iter()
        .map(|p| {
            let magnitude = p.color_grad[0].hypot(p.color_grad[1]);
            (magnitude > threshold)
                .then(|| [p.color_grad[0] / magnitude, p.color_grad[1] / magnitude])
        })
        .collect();
    let self_weight = cubic_spline_kernel(0.0, h);
    let curvatures: Vec<f64> = particles
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            let Some(n_i) = normals[i] else {
                return 0.0;
            };
            let mut divergence = 0.0_f64;
            let mut unity = pi.mass / pi.density.max(1e-30) * self_weight;
            if let Some(list) = neighbors.get(i) {
                for &j in list {
                    if j == i {
                        continue;
                    }
                    let (Some(pj), Some(n_j)) =
                        (particles.get(j), normals.get(j).copied().flatten())
                    else {
                        continue;
                    };
                    let (dx, dy, r) = separation(pi, pj);
                    if r < 1e-14 || r > 2.0 * h {
                        continue;
                    }
                    let (gx, gy) = cubic_spline_gradient(dx, dy, r, h);
                    let volume = pj.mass / pj.density.max(1e-30);
                    divergence += volume * ((n_j[0] - n_i[0]) * gx + (n_j[1] - n_i[1]) * gy);
                    unity += volume * cubic_spline_kernel(r, h);
                }
            }
            if unity > 0.0 {
                -divergence / unity
            } else {
                0.0
            }
        })
        .collect();
    for (p, kappa) in particles.iter_mut().zip(curvatures) {
        p.curvature = kappa;
    }
}

/// Momentum-conserving two-phase pressure acceleration
/// `aᵢ = −Σⱼ mⱼ (Pᵢ + Pⱼ)/(ρᵢ ρⱼ) ∇ᵢWᵢⱼ` of every particle.
///
/// For equal densities this is the familiar `Pᵢ/ρᵢ² + Pⱼ/ρⱼ²` form; for a
/// liquid–gas pair it keeps the force on each side proportional to the common
/// interface pressure instead of to `P/ρ²` of the lighter phase.  Each pair
/// contributes equal and opposite momentum `mᵢ aᵢ = −mⱼ aⱼ`.
pub(super) fn pressure_accelerations(
    particles: &[ColorParticle],
    h: f64,
    neighbors: &[Vec<usize>],
) -> Vec<[f64; 2]> {
    particles
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            let mut acc = [0.0_f64; 2];
            if !is_positive(pi.density) {
                return acc;
            }
            if let Some(list) = neighbors.get(i) {
                for &j in list {
                    if j == i {
                        continue;
                    }
                    let Some(pj) = particles.get(j) else {
                        continue;
                    };
                    if !is_positive(pj.density) {
                        continue;
                    }
                    let (dx, dy, r) = separation(pi, pj);
                    if r < 1e-14 || r > 2.0 * h {
                        continue;
                    }
                    let (gx, gy) = cubic_spline_gradient(dx, dy, r, h);
                    let coeff = -pj.mass * (pi.pressure + pj.pressure) / (pi.density * pj.density);
                    acc[0] += coeff * gx;
                    acc[1] += coeff * gy;
                }
            }
            acc
        })
        .collect()
}

/// Continuum-surface-force acceleration `aᵢ = σ κᵢ ∇Cᵢ / ρᵢ` of every particle
/// from the stored colour gradients and curvatures.
///
/// `∇C` points into the `C = 1` phase and `κ > 0` where that phase is convex,
/// so the acceleration pulls a convex phase inwards.
pub(super) fn csf_accelerations(particles: &[ColorParticle], sigma: f64) -> Vec<[f64; 2]> {
    particles
        .iter()
        .map(|p| {
            if p.density > 0.0 {
                let factor = sigma * p.curvature / p.density;
                [factor * p.color_grad[0], factor * p.color_grad[1]]
            } else {
                [0.0; 2]
            }
        })
        .collect()
}
