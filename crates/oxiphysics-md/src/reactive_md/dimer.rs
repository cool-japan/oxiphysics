// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dimer method for saddle-point (transition-state) search.
//!
//! The types are re-exported from the parent module; the image and curvature
//! conventions are documented on [`TransitionStateSearch`].

use crate::MdError;
use std::f64::consts::FRAC_PI_2;

/// Outcome of [`TransitionStateSearch::rotate`].
#[derive(Clone, Debug, PartialEq)]
pub struct DimerRotation {
    /// Curvature `C = Nᵀ H N` along the final dimer axis (energy / length²).
    pub curvature: f64,
    /// Number of rotation iterations performed (each rotates the axis once).
    pub steps: usize,
    /// `true` if the rotation converged: the rotation angle that would still
    /// lower the curvature fell below `rotation_tol`, or the rotational force
    /// vanished.  `false` if `max_rotation_steps` was exhausted.
    pub converged: bool,
    /// Number of calls of the force routine made by the rotation.
    pub force_evaluations: usize,
}

/// Outcome of [`TransitionStateSearch::find_saddle`].
#[derive(Clone, Debug, PartialEq)]
pub struct DimerSearchResult {
    /// Final midpoint of the dimer (the saddle-point estimate).
    pub midpoint: Vec<[f64; 3]>,
    /// Final unit dimer axis, aligned with the lowest-curvature mode.
    pub axis: Vec<[f64; 3]>,
    /// Curvature along the final axis; negative at a first-order saddle point.
    pub curvature: f64,
    /// Largest absolute force component at the final midpoint.
    pub max_force: f64,
    /// Number of translation iterations performed.
    pub iterations: usize,
    /// `true` if `max_force` fell below the tolerance with negative curvature.
    pub converged: bool,
    /// Total number of calls of the force routine.
    pub force_evaluations: usize,
}

/// Dimer method for saddle point search (transition state finding).
///
/// Reference: Henkelman & Jonsson, J. Chem. Phys. 111, 7010 (1999).
///
/// # Conventions
///
/// The dimer consists of two images placed symmetrically about the midpoint
/// `R₀` along the unit axis `N` (a vector with one 3-vector per atom):
///
/// ```text
/// R₁ = R₀ + L·N,     R₂ = R₀ − L·N,     L = dimer_length (half the image distance)
/// ```
///
/// With the forces `F = −∇V` at the images, the curvature of the potential along
/// the axis is `C = (F₂ − F₁)·N / (2L) ≈ Nᵀ H N` (`H` = Hessian), so that a
/// harmonic well `V = ½ k |x|²` has `C = +k` and the downhill direction of a
/// saddle has `C < 0`.  The rotational force `F⊥ = (F₁ − F₂) − ((F₁ − F₂)·N) N`
/// points towards smaller curvature; rotating `N` along it until `C` is minimal
/// ([`TransitionStateSearch::rotate`]) aligns the dimer with the lowest-curvature
/// mode of the Hessian.  The midpoint then moves with the effective force
/// ([`TransitionStateSearch::translate_midpoint`])
///
/// ```text
/// F† = −F∥            if C > 0   (climb along the mode out of the convex region)
/// F† = F − 2 F∥       if C < 0   (climb along the mode, relax in all other directions)
/// ```
///
/// where `F∥ = (F·N) N`, which converges to a first-order saddle point
/// ([`TransitionStateSearch::find_saddle`]).
#[derive(Clone, Debug)]
pub struct TransitionStateSearch {
    /// Half-distance `L` between the dimer midpoint and each image (Angstrom).
    pub dimer_length: f64,
    /// Trial rotation angle `Δθ` (radians, `0 < Δθ < π/2`) used to probe the
    /// curvature along a rotated axis in the modified-Newton rotation.
    pub rotation_angle: f64,
    /// Translation step size: the midpoint moves by `translation_step · F†`.
    pub translation_step: f64,
    /// Maximum rotation iterations per call of [`TransitionStateSearch::rotate`].
    pub max_rotation_steps: usize,
    /// Rotation convergence criterion: the rotation is converged once the angle
    /// to the curvature minimum is below this value (radians).
    pub rotation_tol: f64,
}

/// Flatten per-atom 3-vectors into one `3N` vector.
fn flatten(vectors: &[[f64; 3]]) -> Vec<f64> {
    vectors.iter().flat_map(|v| v.iter().copied()).collect()
}

/// Dot product of two flat vectors.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

impl TransitionStateSearch {
    /// Create dimer method with default parameters.
    pub fn new() -> Self {
        Self {
            dimer_length: 0.01,
            rotation_angle: 0.01,
            translation_step: 0.1,
            max_rotation_steps: 10,
            rotation_tol: 1e-4,
        }
    }

    /// Compute the effective curvature along the dimer axis.
    ///
    /// C = (F2 - F1) · N / (2 * dimer_length)
    ///
    /// where `N` is the dimer unit vector, `F1` the force at image 1
    /// (`R0 + L·N`, see [`TransitionStateSearch::build_image1`]) and `F2` the
    /// force at image 2 (`R0 − L·N`, see [`TransitionStateSearch::build_image2`]).
    /// For `V = ½ k |x|²` this returns `+k`.
    pub fn curvature(
        &self,
        force_r1: &[[f64; 3]],
        force_r2: &[[f64; 3]],
        dimer_axis: &[[f64; 3]],
    ) -> f64 {
        let delta_f: f64 = force_r1
            .iter()
            .zip(force_r2.iter())
            .zip(dimer_axis.iter())
            .map(|((f1, f2), n)| {
                (f2[0] - f1[0]) * n[0] + (f2[1] - f1[1]) * n[1] + (f2[2] - f1[2]) * n[2]
            })
            .sum();
        delta_f / (2.0 * self.dimer_length)
    }

    /// Construct dimer image 1, `R0 + L·N`, from the midpoint and axis direction.
    pub fn build_image1(
        midpoint: &[[f64; 3]],
        dimer_axis: &[[f64; 3]],
        dimer_length: f64,
    ) -> Vec<[f64; 3]> {
        midpoint
            .iter()
            .zip(dimer_axis.iter())
            .map(|(r, n)| {
                [
                    r[0] + dimer_length * n[0],
                    r[1] + dimer_length * n[1],
                    r[2] + dimer_length * n[2],
                ]
            })
            .collect()
    }

    /// Construct dimer image 2, `R0 − L·N`, from the midpoint and axis direction.
    pub fn build_image2(
        midpoint: &[[f64; 3]],
        dimer_axis: &[[f64; 3]],
        dimer_length: f64,
    ) -> Vec<[f64; 3]> {
        midpoint
            .iter()
            .zip(dimer_axis.iter())
            .map(|(r, n)| {
                [
                    r[0] - dimer_length * n[0],
                    r[1] - dimer_length * n[1],
                    r[2] - dimer_length * n[2],
                ]
            })
            .collect()
    }

    /// Translate the dimer midpoint with the Henkelman-Jónsson effective force.
    ///
    /// With `F∥ = (F·N) N` the projection of `force` onto the dimer axis, the
    /// midpoint moves by `translation_step · F†` where
    ///
    /// * `F† = F − 2 F∥ = F⊥ − F∥` for negative curvature (climb along the mode,
    ///   relax perpendicular to it), and
    /// * `F† = −F∥` for non-negative curvature (the perpendicular force is
    ///   dropped and the dimer climbs out of the convex region along the mode;
    ///   following the full force there would only slide back into the minimum).
    ///
    /// `dimer_axis` need not be normalised; only its direction matters.
    pub fn translate_midpoint(
        &self,
        midpoint: &mut [[f64; 3]],
        force: &[[f64; 3]],
        dimer_axis: &[[f64; 3]],
        curvature: f64,
    ) {
        let axis_norm_sq: f64 = dimer_axis
            .iter()
            .map(|n| n[0] * n[0] + n[1] * n[1] + n[2] * n[2])
            .sum();
        let f_parallel: f64 = if axis_norm_sq > 0.0 {
            force
                .iter()
                .zip(dimer_axis.iter())
                .map(|(f, n)| f[0] * n[0] + f[1] * n[1] + f[2] * n[2])
                .sum::<f64>()
                / axis_norm_sq
        } else {
            0.0
        };

        for ((r, f), n) in midpoint.iter_mut().zip(force.iter()).zip(dimer_axis.iter()) {
            for d in 0..3 {
                let f_para = f_parallel * n[d];
                let effective = if curvature < 0.0 {
                    f[d] - 2.0 * f_para
                } else {
                    -f_para
                };
                r[d] += self.translation_step * effective;
            }
        }
    }

    /// Check that the dimer parameters are usable.
    fn validate_parameters(&self) -> Result<(), MdError> {
        if !(self.dimer_length.is_finite() && self.dimer_length > 0.0) {
            return Err(MdError::General(format!(
                "dimer_length must be positive and finite, got {}",
                self.dimer_length
            )));
        }
        if !(self.rotation_angle > 0.0 && self.rotation_angle < FRAC_PI_2) {
            return Err(MdError::General(format!(
                "rotation_angle must lie in (0, pi/2) radians, got {}",
                self.rotation_angle
            )));
        }
        Ok(())
    }

    /// Rotate the dimer axis towards the lowest-curvature mode.
    ///
    /// `midpoint` is the dimer midpoint `R₀`, `force_midpoint` the force `F₀`
    /// there (which the translation needs anyway) and `axis` the dimer axis,
    /// which is normalised and rotated **in place**.  `force_fn` returns the
    /// forces `−∇V` for a configuration and is called once for image 1 at the
    /// start and then twice per iteration (one trial rotation, one evaluation at
    /// the rotated axis); the forces of image 2 are obtained from
    /// `F₂ ≈ 2 F₀ − F₁`, which is exact for a force that is linear along the axis
    /// (a quadratic potential) and accurate to `O(L²)` otherwise.
    ///
    /// Every iteration
    ///
    /// 1. forms `G = (F₀ − F₁)/L ≈ H N`, the curvature `C₀ = N·G` and the
    ///    rotational force `F⊥ ∝ −(G − C₀ N)` that defines the unit direction `Θ`
    ///    perpendicular to `N`,
    /// 2. measures the curvature `C_Δθ` along the trial axis
    ///    `N cos Δθ + Θ sin Δθ` (`Δθ = rotation_angle`),
    /// 3. fits `C(θ) = a₀/2 + a₁ cos 2θ + b₁ sin 2θ` to `C₀`, `C'(0) = 2 b₁ = −2|G⊥|`
    ///    and `C_Δθ` (Henkelman-Jónsson modified Newton),
    ///    `a₁ = (C₀ − C_Δθ + b₁ sin 2Δθ) / (1 − cos 2Δθ)`, and jumps to its minimum
    ///    `θ_min = ½ atan2(b₁, a₁) + π/2 ∈ (0, π/2)`.
    ///
    /// The rotation stops when `θ_min < rotation_tol` (or the rotational force
    /// vanishes), or after `max_rotation_steps` iterations.  The curvature never
    /// increases from iteration to iteration for a quadratic potential, where the
    /// fit is exact within the rotation plane.
    ///
    /// # Errors
    /// [`MdError::General`] for a non-positive `dimer_length`, a `rotation_angle`
    /// outside `(0, π/2)`, mismatched lengths of `midpoint`, `axis` and
    /// `force_midpoint`, a zero or non-finite axis, a `force_fn` result with the
    /// wrong number of atoms, or non-finite forces.
    pub fn rotate<F>(
        &self,
        midpoint: &[[f64; 3]],
        axis: &mut [[f64; 3]],
        force_midpoint: &[[f64; 3]],
        mut force_fn: F,
    ) -> Result<DimerRotation, MdError>
    where
        F: FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
    {
        self.validate_parameters()?;
        let n_atoms = midpoint.len();
        if axis.len() != n_atoms || force_midpoint.len() != n_atoms {
            return Err(MdError::General(format!(
                "dimer rotation needs matching lengths: {n_atoms} midpoint atoms, {} axis \
                 vectors, {} midpoint forces",
                axis.len(),
                force_midpoint.len()
            )));
        }
        let r0 = flatten(midpoint);
        let f0 = flatten(force_midpoint);
        let mut n = flatten(axis);
        let axis_norm = dot(&n, &n).sqrt();
        if !(axis_norm.is_finite() && axis_norm > 0.0) {
            return Err(MdError::General(
                "the dimer axis must be a non-zero finite vector".into(),
            ));
        }
        for component in n.iter_mut() {
            *component /= axis_norm;
        }

        let length = self.dimer_length;
        let trial = self.rotation_angle;
        let (sin_trial, cos_trial) = trial.sin_cos();
        let (sin_2trial, cos_2trial) = (2.0 * trial).sin_cos();
        let mut evaluations = 0_usize;

        let mut f1 = Self::image_force(&mut force_fn, &r0, &n, length, &mut evaluations)?;
        let mut steps = 0_usize;
        let mut converged_curvature: Option<f64> = None;

        while steps < self.max_rotation_steps {
            let g = Self::curvature_vector(&f0, &f1, length);
            let c0 = dot(&n, &g);
            if !c0.is_finite() {
                return Err(MdError::General(
                    "non-finite dimer curvature (check the forces)".into(),
                ));
            }
            let g_perp: Vec<f64> = g
                .iter()
                .zip(n.iter())
                .map(|(gi, ni)| gi - c0 * ni)
                .collect();
            let g_perp_norm = dot(&g_perp, &g_perp).sqrt();
            if g_perp_norm <= 1e-12 * (1.0 + c0.abs()) {
                // N is an eigenvector: the rotational force vanishes.
                converged_curvature = Some(c0);
                break;
            }
            let theta_dir: Vec<f64> = g_perp.iter().map(|x| -x / g_perp_norm).collect();

            // Curvature along the trial axis N cos Δθ + Θ sin Δθ.
            let trial_axis: Vec<f64> = n
                .iter()
                .zip(theta_dir.iter())
                .map(|(ni, ti)| ni * cos_trial + ti * sin_trial)
                .collect();
            let f1_trial =
                Self::image_force(&mut force_fn, &r0, &trial_axis, length, &mut evaluations)?;
            let g_trial = Self::curvature_vector(&f0, &f1_trial, length);
            let c_trial = dot(&trial_axis, &g_trial);

            // C(θ) = a0/2 + a1 cos 2θ + b1 sin 2θ with C'(0) = 2 b1 = −2 |G⊥|.
            let b1 = -g_perp_norm;
            let a1 = (c0 - c_trial + b1 * sin_2trial) / (1.0 - cos_2trial);
            let theta = 0.5 * b1.atan2(a1) + FRAC_PI_2;
            let c_min = (c0 - a1) + a1 * (2.0 * theta).cos() + b1 * (2.0 * theta).sin();
            if !(theta.is_finite() && c_min.is_finite()) {
                return Err(MdError::General(
                    "non-finite dimer rotation angle (check the forces)".into(),
                ));
            }

            let (sin_theta, cos_theta) = theta.sin_cos();
            let mut rotated: Vec<f64> = n
                .iter()
                .zip(theta_dir.iter())
                .map(|(ni, ti)| ni * cos_theta + ti * sin_theta)
                .collect();
            let rotated_norm = dot(&rotated, &rotated).sqrt();
            for component in rotated.iter_mut() {
                *component /= rotated_norm;
            }
            n = rotated;
            steps += 1;

            if theta < self.rotation_tol {
                converged_curvature = Some(c_min);
                break;
            }
            f1 = Self::image_force(&mut force_fn, &r0, &n, length, &mut evaluations)?;
        }

        let (curvature, converged) = match converged_curvature {
            Some(curvature) => (curvature, true),
            None => {
                // Out of iterations: report the curvature along the final axis.
                let g = Self::curvature_vector(&f0, &f1, length);
                let curvature = dot(&n, &g);
                if !curvature.is_finite() {
                    return Err(MdError::General(
                        "non-finite dimer curvature (check the forces)".into(),
                    ));
                }
                (curvature, false)
            }
        };

        for (slot, chunk) in axis.iter_mut().zip(n.chunks_exact(3)) {
            *slot = [chunk[0], chunk[1], chunk[2]];
        }
        Ok(DimerRotation {
            curvature,
            steps,
            converged,
            force_evaluations: evaluations,
        })
    }

    /// `G = (F₀ − F₁)/L ≈ H N` for the extrapolated image forces `F₂ = 2F₀ − F₁`.
    fn curvature_vector(f0: &[f64], f1: &[f64], length: f64) -> Vec<f64> {
        f0.iter()
            .zip(f1.iter())
            .map(|(a, b)| (a - b) / length)
            .collect()
    }

    /// Force at image 1 (`R₀ + L·N`), flattened, with a length check.
    fn image_force<F>(
        force_fn: &mut F,
        r0: &[f64],
        axis: &[f64],
        length: f64,
        evaluations: &mut usize,
    ) -> Result<Vec<f64>, MdError>
    where
        F: FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
    {
        let image: Vec<[f64; 3]> = r0
            .chunks_exact(3)
            .zip(axis.chunks_exact(3))
            .map(|(r, n)| {
                [
                    r[0] + length * n[0],
                    r[1] + length * n[1],
                    r[2] + length * n[2],
                ]
            })
            .collect();
        let forces = force_fn(&image);
        *evaluations += 1;
        if forces.len() != image.len() {
            return Err(MdError::General(format!(
                "the force routine returned {} atoms for a dimer image of {} atoms",
                forces.len(),
                image.len()
            )));
        }
        Ok(flatten(&forces))
    }

    /// Search for a first-order saddle point with the dimer method.
    ///
    /// Starting from `start` with the initial axis `initial_axis` (any non-zero
    /// vector; it is normalised), each iteration evaluates `F₀ = force_fn(R₀)`,
    /// rotates the dimer into the lowest-curvature mode
    /// ([`TransitionStateSearch::rotate`]) and translates the midpoint with the
    /// effective force ([`TransitionStateSearch::translate_midpoint`]).  The
    /// search stops when the largest force component at the midpoint is at most
    /// `force_tol` **and** the curvature along the axis is negative, or after
    /// `max_iter` iterations.
    ///
    /// The dimer follows the lowest-curvature mode, so it needs a start from which
    /// that mode leads to the saddle; the search cannot leave a minimum whose
    /// force vanishes exactly.
    ///
    /// # Errors
    /// See [`TransitionStateSearch::rotate`].
    pub fn find_saddle<F>(
        &self,
        start: &[[f64; 3]],
        initial_axis: &[[f64; 3]],
        force_tol: f64,
        max_iter: usize,
        mut force_fn: F,
    ) -> Result<DimerSearchResult, MdError>
    where
        F: FnMut(&[[f64; 3]]) -> Vec<[f64; 3]>,
    {
        let mut midpoint = start.to_vec();
        let mut axis = initial_axis.to_vec();
        let mut force_evaluations = 0_usize;
        let mut curvature = f64::NAN;
        let mut max_force = f64::INFINITY;
        let mut iterations = 0_usize;
        let mut converged = false;

        while iterations <= max_iter {
            let f0 = force_fn(&midpoint);
            force_evaluations += 1;
            if f0.len() != midpoint.len() {
                return Err(MdError::General(format!(
                    "the force routine returned {} atoms for a system of {} atoms",
                    f0.len(),
                    midpoint.len()
                )));
            }
            let rotation = self.rotate(&midpoint, &mut axis, &f0, &mut force_fn)?;
            force_evaluations += rotation.force_evaluations;
            curvature = rotation.curvature;
            max_force = f0
                .iter()
                .flat_map(|f| f.iter())
                .fold(0.0_f64, |acc, c| acc.max(c.abs()));
            if max_force <= force_tol && curvature < 0.0 {
                converged = true;
                break;
            }
            if iterations == max_iter {
                break;
            }
            self.translate_midpoint(&mut midpoint, &f0, &axis, curvature);
            iterations += 1;
        }

        Ok(DimerSearchResult {
            midpoint,
            axis,
            curvature,
            max_force,
            iterations,
            converged,
            force_evaluations,
        })
    }
}

impl Default for TransitionStateSearch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests;
