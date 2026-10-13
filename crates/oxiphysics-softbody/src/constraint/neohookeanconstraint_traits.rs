//! # NeoHookeanConstraint - Trait Implementations
//!
//! This module contains trait implementations for `NeoHookeanConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
use super::types::NeoHookeanConstraint;

/// Frobenius norm of the deformation gradient below which the direction of the
/// deviatoric gradient is undefined; the deviatoric constraint is skipped there
/// (its value is zero as well).
const MIN_FROBENIUS_NORM: Real = 1.0e-12;

/// Smallest weighted-gradient-plus-compliance denominator that is divided by.
const MIN_DENOMINATOR: Real = 1.0e-14;

/// Relative size of the coupled 2 x 2 determinant below which the two
/// constraints are solved independently instead of jointly.
const MIN_RELATIVE_DETERMINANT: Real = 1.0e-12;

/// Armijo constant of the step-length control: a step is accepted when it
/// achieves this fraction of the decrease that the linear model predicts.
const SUFFICIENT_DECREASE: Real = 1.0e-4;

/// Largest number of step halvings; the smallest fraction tried is `2^-30`.
const MAX_STEP_HALVINGS: usize = 30;

/// Predicted decrease, relative to the potential itself, below which the
/// element counts as converged and the full step is taken unchecked: the
/// decrease is then at the round-off level of the potential.
const NEGLIGIBLE_DECREASE: Real = 1.0e-12;

/// Cofactor matrix of a column-major 3 x 3 matrix.
///
/// For `F = [f1 f2 f3]` the columns of `cof(F) = dJ/dF` are
/// `[f2 x f3, f3 x f1, f1 x f2]`.
fn cofactor(f: [Real; 9]) -> [Real; 9] {
    let f1 = Vec3::new(f[0], f[1], f[2]);
    let f2 = Vec3::new(f[3], f[4], f[5]);
    let f3 = Vec3::new(f[6], f[7], f[8]);
    let c1 = f2.cross(&f3);
    let c2 = f3.cross(&f1);
    let c3 = f1.cross(&f2);
    [c1.x, c1.y, c1.z, c2.x, c2.y, c2.z, c3.x, c3.y, c3.z]
}

/// `numerator / denominator`, or `0` when the denominator is not safely positive.
fn safe_ratio(numerator: Real, denominator: Real) -> Real {
    if denominator > MIN_DENOMINATOR {
        numerator / denominator
    } else {
        0.0
    }
}

/// Values and vertex gradients of the two constraints of one tetrahedron.
#[derive(Debug, Clone, Copy)]
pub(super) struct NeoHookeanTerms {
    /// Deviatoric constraint `C_D = sqrt(tr(F^T F))`.
    pub(super) c_dev: Real,
    /// Hydrostatic constraint `C_H = det(F) - gamma`.
    pub(super) c_hyd: Real,
    /// Gradient of `C_D` with respect to the four vertex positions.
    pub(super) grad_dev: [Vec3; 4],
    /// Gradient of `C_H` with respect to the four vertex positions.
    pub(super) grad_hyd: [Vec3; 4],
    /// `false` when `F` is numerically zero and the gradient of `C_D` is undefined.
    pub(super) dev_defined: bool,
}

impl NeoHookeanConstraint {
    /// Stiffness `K = lambda + mu` of the hydrostatic constraint.
    ///
    /// With `K` the first Piola stress `P = mu F + K (J - gamma) cof(F)` has the
    /// linear-elastic small-strain limit of the Lamé parameters `(mu, lambda)` and
    /// the rest state is stable for every Poisson ratio below 1/2.  Using `lambda`
    /// instead would give the bulk modulus `lambda - mu/3` and an unstable rest
    /// state for Poisson ratios below 1/8.  A non-positive result disables the
    /// hydrostatic constraint.
    pub fn hydrostatic_modulus(&self) -> Real {
        self.lambda_lame + self.mu.max(0.0)
    }

    /// Value `gamma = 1 + mu / K` that `det(F)` is driven towards.
    ///
    /// It makes the rest state `F = I` stress free.
    pub fn hydrostatic_target(&self) -> Real {
        let modulus = self.hydrostatic_modulus();
        if modulus > 0.0 {
            1.0 + self.mu.max(0.0) / modulus
        } else {
            1.0
        }
    }

    /// First Piola-Kirchhoff stress `P = dPsi/dF` (column-major) of the current
    /// configuration, with `Psi = mu/2 (I1 - 3) + K/2 (J - gamma)^2 + const`.
    ///
    /// It vanishes at rest, and for small strains it equals Hooke's law with the
    /// Lamé parameters `(mu, lambda)`.
    pub fn first_piola_stress(&self, particles: &[SoftParticle]) -> [Real; 9] {
        let f = self.deformation_gradient(particles);
        let shear = self.mu.max(0.0);
        let modulus = self.hydrostatic_modulus();
        let pressure = if modulus > 0.0 {
            modulus * (Self::compute_j(f) - self.hydrostatic_target())
        } else {
            0.0
        };
        let cof = cofactor(f);
        std::array::from_fn(|n| shear * f[n] + pressure * cof[n])
    }

    /// Elastic energy `V0 * Psi(F)` of the element in joules, zero at rest.
    pub fn elastic_energy(&self, particles: &[SoftParticle]) -> Real {
        let f = self.deformation_gradient(particles);
        let shear = self.mu.max(0.0);
        let modulus = self.hydrostatic_modulus();
        let mut density = 0.5 * shear * (Self::compute_i1(f) - 3.0);
        if modulus > 0.0 {
            let gamma = self.hydrostatic_target();
            let volumetric = Self::compute_j(f) - gamma;
            density += 0.5 * modulus * (volumetric * volumetric - (1.0 - gamma) * (1.0 - gamma));
        }
        self.rest_volume * density
    }

    /// Vertex gradients of a constraint with `P = dC/dF` (column-major).
    ///
    /// Vertices 1..3 receive the columns of `P Dm^-T`; vertex 0 gets the
    /// negated sum so that the gradient is translation invariant.
    fn vertex_gradients(&self, p: [Real; 9]) -> [Vec3; 4] {
        let b = &self.dm_inv;
        let mut gradients = [Vec3::zeros(); 4];
        for m in 0..3 {
            let mut gradient = Vec3::zeros();
            for j in 0..3 {
                let column = Vec3::new(p[3 * j], p[3 * j + 1], p[3 * j + 2]);
                gradient += column * b[m + 3 * j];
            }
            gradients[m + 1] = gradient;
        }
        gradients[0] = -(gradients[1] + gradients[2] + gradients[3]);
        gradients
    }

    /// Evaluate both constraints and their vertex gradients.
    ///
    /// The tet indices must be valid for `particles`.
    pub(super) fn evaluate_terms(&self, particles: &[SoftParticle]) -> NeoHookeanTerms {
        let f = self.deformation_gradient(particles);
        let frobenius = Self::compute_i1(f).sqrt();
        let dev_defined = frobenius > MIN_FROBENIUS_NORM;
        let dev_direction = if dev_defined {
            f.map(|entry| entry / frobenius)
        } else {
            [0.0; 9]
        };
        NeoHookeanTerms {
            c_dev: frobenius,
            c_hyd: Self::compute_j(f) - self.hydrostatic_target(),
            grad_dev: self.vertex_gradients(dev_direction),
            grad_hyd: self.vertex_gradients(cofactor(f)),
            dev_defined,
        }
    }

    /// Values `(C_D, C_H)` of the two constraints for an element whose corners
    /// are at `positions`, given in the order of `indices`.
    fn constraint_values_at(&self, positions: &[Vec3; 4]) -> (Real, Real) {
        let f = self.deformation_gradient_of(positions);
        (
            Self::compute_i1(f).sqrt(),
            Self::compute_j(f) - self.hydrostatic_target(),
        )
    }
}

/// The element-local implicit-Euler problem that one coupled projection
/// linearises.
///
/// XPBD reads `x = x~ + M^-1 sum_i grad C_i lambda_i` with `lambda_i = -C_i / a_i`
/// (`a_i` the compliance divided by `dt^2`) as the stationarity condition of
///
/// ```text
/// G(z) = 1/2 |z - x~|_M^2 + sum_i C_i(z)^2 / (2 a_i),
/// ```
///
/// where `x~` is the inertial prediction.  The constraint only knows the
/// multipliers accumulated so far, so it reconstructs the prediction as
/// `x~ = x - M^-1 sum_i grad C_i lambda_i`; at the first projection of a sub-step
/// the multipliers are zero and `x~` is simply the current position.  The XPBD
/// update is the Gauss-Newton step of `G`, and `G` is the merit function for its
/// step-length control.
struct LocalPotential<'a> {
    /// The constraint whose element is being projected.
    constraint: &'a NeoHookeanConstraint,
    /// Current positions of the four vertices.
    start: [Vec3; 4],
    /// Inverse masses of the four vertices (zero for pinned vertices).
    weights: [Real; 4],
    /// Displacement `M^-1 sum_i grad C_i lambda_i` of every vertex from the
    /// reconstructed prediction.
    shift: [Vec3; 4],
    /// Compliance of the deviatoric constraint divided by `dt^2`, `None` when
    /// the constraint is inactive.
    alpha_dev: Option<Real>,
    /// Compliance of the hydrostatic constraint divided by `dt^2`, `None` when
    /// the constraint is inactive.
    alpha_hyd: Option<Real>,
}

impl LocalPotential<'_> {
    /// `G` for vertices that are `fraction * step` away from `start` and have the
    /// constraint values `c_dev` and `c_hyd` there.
    fn assemble(&self, step: &[Vec3; 4], fraction: Real, c_dev: Real, c_hyd: Real) -> Real {
        let mut value = 0.0;
        for ((&weight, shift), step) in self.weights.iter().zip(&self.shift).zip(step) {
            if weight > 0.0 {
                let offset = shift + step * fraction;
                value += 0.5 * offset.norm_squared() / weight;
            }
        }
        if let Some(alpha) = self.alpha_dev {
            value += 0.5 * c_dev * c_dev / alpha;
        }
        if let Some(alpha) = self.alpha_hyd {
            value += 0.5 * c_hyd * c_hyd / alpha;
        }
        value
    }

    /// `G` at the current positions, given the constraint values there.
    fn value_at_start(&self, step: &[Vec3; 4], c_dev: Real, c_hyd: Real) -> Real {
        self.assemble(step, 0.0, c_dev, c_hyd)
    }

    /// `G` after moving the vertices by `fraction * step`.
    fn value_after(&self, step: &[Vec3; 4], fraction: Real) -> Real {
        let moved: [Vec3; 4] = std::array::from_fn(|k| self.start[k] + step[k] * fraction);
        let (c_dev, c_hyd) = self.constraint.constraint_values_at(&moved);
        self.assemble(step, fraction, c_dev, c_hyd)
    }

    /// The largest fraction `2^-n` of `step` that satisfies the Armijo
    /// condition `G(fraction) <= current + c1 * fraction * slope`, or `None`
    /// when even the smallest tried fraction does not (or `G` is not a number).
    fn accepted_fraction(&self, step: &[Vec3; 4], current: Real, slope: Real) -> Option<Real> {
        let mut fraction = 1.0;
        for _ in 0..=MAX_STEP_HALVINGS {
            if self.value_after(step, fraction) <= current + SUFFICIENT_DECREASE * fraction * slope
            {
                return Some(fraction);
            }
            fraction *= 0.5;
        }
        None
    }
}

impl SoftConstraint for NeoHookeanConstraint {
    /// One coupled XPBD update of the deviatoric and hydrostatic constraints.
    ///
    /// Both constraints are projected together by solving the 2 x 2 system
    ///
    /// ```text
    /// [ S_DD  S_DH ] [ dl_D ]   [ -C_D - a_D l_D ]
    /// [ S_DH  S_HH ] [ dl_H ] = [ -C_H - a_H l_H ]
    /// ```
    ///
    /// with `S_ij = sum_k w_k grad_k C_i . grad_k C_j + a_i`, the compliances
    /// `a_D = 1 / (mu V0 dt^2)` and `a_H = (1 / (K V0) + compliance) / dt^2`, and
    /// then moving every vertex by `w_k (grad_k C_D dl_D + grad_k C_H dl_H)`.
    /// At rest the two constraint forces cancel, so the update is zero.
    ///
    /// The step is the Gauss-Newton step of the element-local implicit-Euler
    /// problem and is applied in full whenever it sufficiently decreases that
    /// problem's potential (Armijo condition).  A step that would overshoot,
    /// as a stiff element does over a long sub-step, is halved until it does,
    /// and the multipliers advance by the same fraction.  Consequently, up to
    /// round-off, a projection never raises the elastic energy of an element
    /// that starts without accumulated multipliers (after `reset_lambda`).
    ///
    /// The constraint does nothing for a non-positive `dt_sub`, an invalid tet
    /// index or a degenerate rest shape.
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let usable = dt_sub.is_finite() && dt_sub > 0.0 && self.rest_volume > 0.0;
        if !usable {
            return;
        }
        if self.indices.iter().any(|&index| index >= particles.len()) {
            return;
        }
        let terms = self.evaluate_terms(particles);
        let weights = self.indices.map(|index| particles[index].inverse_mass);
        let dt_sq = dt_sub * dt_sub;
        let shear = self.mu.max(0.0);
        let modulus = self.hydrostatic_modulus();
        let dev_active = shear > 0.0 && terms.dev_defined;
        let hyd_active = modulus > 0.0;
        if !dev_active && !hyd_active {
            return;
        }
        let alpha_dev = if dev_active {
            1.0 / (shear * self.rest_volume * dt_sq)
        } else {
            0.0
        };
        let alpha_hyd = if hyd_active {
            (1.0 / (modulus * self.rest_volume) + self.compliance.max(0.0)) / dt_sq
        } else {
            0.0
        };

        let gram = |a: &[Vec3; 4], b: &[Vec3; 4]| -> Real {
            (0..4).map(|k| weights[k] * a[k].dot(&b[k])).sum()
        };
        let s_dd = gram(&terms.grad_dev, &terms.grad_dev) + alpha_dev;
        let s_hh = gram(&terms.grad_hyd, &terms.grad_hyd) + alpha_hyd;
        let s_dh = gram(&terms.grad_dev, &terms.grad_hyd);
        let rhs_dev = -terms.c_dev - alpha_dev * self.lambda_dev;
        let rhs_hyd = -terms.c_hyd - alpha_hyd * self.lambda_vol;

        let (delta_dev, delta_hyd) = match (dev_active, hyd_active) {
            (true, true) => {
                let determinant = s_dd * s_hh - s_dh * s_dh;
                if determinant > MIN_RELATIVE_DETERMINANT * s_dd * s_hh {
                    (
                        (rhs_dev * s_hh - rhs_hyd * s_dh) / determinant,
                        (s_dd * rhs_hyd - s_dh * rhs_dev) / determinant,
                    )
                } else {
                    (safe_ratio(rhs_dev, s_dd), safe_ratio(rhs_hyd, s_hh))
                }
            }
            (true, false) => (safe_ratio(rhs_dev, s_dd), 0.0),
            (false, true) => (0.0, safe_ratio(rhs_hyd, s_hh)),
            (false, false) => return,
        };
        if !(delta_dev.is_finite() && delta_hyd.is_finite()) {
            return;
        }

        // Full Gauss-Newton displacement of every vertex.
        let displacement: [Vec3; 4] = std::array::from_fn(|k| {
            (terms.grad_dev[k] * delta_dev + terms.grad_hyd[k] * delta_hyd) * weights[k]
        });
        // Step-length control on the element-local potential.  Its gradient is
        // `grad C_D (l_D + C_D / a_D) + grad C_H (l_H + C_H / a_H)` per vertex.
        let fraction = {
            let potential = LocalPotential {
                constraint: self,
                start: self.indices.map(|index| particles[index].position),
                weights,
                shift: std::array::from_fn(|k| {
                    (terms.grad_dev[k] * self.lambda_dev + terms.grad_hyd[k] * self.lambda_vol)
                        * weights[k]
                }),
                alpha_dev: dev_active.then_some(alpha_dev),
                alpha_hyd: hyd_active.then_some(alpha_hyd),
            };
            let pull_dev = if dev_active {
                -rhs_dev / alpha_dev
            } else {
                0.0
            };
            let pull_hyd = if hyd_active {
                -rhs_hyd / alpha_hyd
            } else {
                0.0
            };
            let slope: Real = (0..4)
                .map(|k| {
                    (terms.grad_dev[k] * pull_dev + terms.grad_hyd[k] * pull_hyd)
                        .dot(&displacement[k])
                })
                .sum();
            let current = potential.value_at_start(&displacement, terms.c_dev, terms.c_hyd);
            if -slope <= NEGLIGIBLE_DECREASE * current.abs() {
                1.0
            } else {
                match potential.accepted_fraction(&displacement, current, slope) {
                    Some(fraction) => fraction,
                    None => return,
                }
            }
        };

        self.lambda_dev += fraction * delta_dev;
        self.lambda_vol += fraction * delta_hyd;
        for (k, &index) in self.indices.iter().enumerate() {
            particles[index].position += displacement[k] * fraction;
        }
    }

    fn reset_lambda(&mut self) {
        Self::reset_lambda(self);
    }
}
