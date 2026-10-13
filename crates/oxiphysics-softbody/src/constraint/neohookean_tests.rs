// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the stable Neo-Hookean XPBD constraint.
//!
//! The old projection shifted the `y` coordinate of every vertex (a rigid
//! translation that cannot change `det(F)`), pushed vertices radially away from
//! the centroid, treated the stiffness as a compliance and never read
//! `compliance`.  The expected values below come from the continuum
//! mechanics of the model itself: finite-difference gradients, Hooke's law in
//! the small-strain limit, and the exact rest shape that an energy-conserving
//! relaxation must reach.

use super::functions::SoftConstraint;
use super::neohookeanconstraint_traits::NeoHookeanTerms;
use super::types::NeoHookeanConstraint;
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Mat3, Real, Vec3};

/// Lamé parameters `(mu, lambda)` from Young's modulus and Poisson's ratio.
fn lame(youngs: Real, poisson: Real) -> (Real, Real) {
    let mu = youngs / (2.0 * (1.0 + poisson));
    let lambda = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    (mu, lambda)
}

fn tet(points: [[Real; 3]; 4]) -> Vec<SoftParticle> {
    points
        .iter()
        .map(|p| SoftParticle::new(Vec3::new(p[0], p[1], p[2]), 1.0))
        .collect()
}

/// Right tetrahedron with unit legs, volume 1/6.
fn unit_tet() -> Vec<SoftParticle> {
    tet([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
}

/// A tetrahedron that is not aligned with the axes.
fn general_tet() -> Vec<SoftParticle> {
    tet([
        [0.3, 0.1, -0.2],
        [1.4, 0.2, 0.1],
        [0.2, 1.1, 0.3],
        [-0.1, 0.2, 1.3],
    ])
}

/// Signed volume of the tetrahedron `particles[0..4]`.
fn signed_volume(particles: &[SoftParticle]) -> Real {
    let e1 = particles[1].position - particles[0].position;
    let e2 = particles[2].position - particles[0].position;
    let e3 = particles[3].position - particles[0].position;
    e1.dot(&e2.cross(&e3)) / 6.0
}

fn mass_weighted_centroid(particles: &[SoftParticle]) -> Vec3 {
    let mut total_mass = 0.0;
    let mut sum = Vec3::zeros();
    for p in particles {
        let mass = 1.0 / p.inverse_mass;
        total_mass += mass;
        sum += p.position * mass;
    }
    sum / total_mass
}

/// Deform `rest` by `x -> map * x + shift`, which gives `F = map`.
fn deform(rest: &[SoftParticle], map: &Mat3, shift: Vec3) -> Vec<SoftParticle> {
    rest.iter()
        .map(|p| SoftParticle::new(map * p.position + shift, 1.0 / p.inverse_mass))
        .collect()
}

fn constraint_for(
    particles: &[SoftParticle],
    youngs: Real,
    poisson: Real,
    compliance: Real,
) -> NeoHookeanConstraint {
    let (mu, lambda) = lame(youngs, poisson);
    NeoHookeanConstraint::from_particles([0, 1, 2, 3], particles, mu, lambda, compliance)
}

/// Mimic the XPBD sub-step loop: predict, reset the multipliers, project every
/// constraint `iterations` times (Gauss-Seidel), derive the velocities and damp
/// them.  There is no external force.
fn relax_mesh_iterated(
    constraints: &mut [NeoHookeanConstraint],
    particles: &mut [SoftParticle],
    frame_dt: Real,
    substeps: usize,
    frames: usize,
    damping: Real,
    iterations: usize,
) {
    let dt = frame_dt / substeps as Real;
    for _ in 0..frames {
        for _ in 0..substeps {
            for p in particles.iter_mut() {
                p.prev_position = p.position;
                if !p.is_static() {
                    p.position += p.velocity * dt;
                }
            }
            for constraint in constraints.iter_mut() {
                constraint.reset_lambda();
            }
            for _ in 0..iterations {
                for constraint in constraints.iter_mut() {
                    constraint.project(particles, dt);
                }
            }
            for p in particles.iter_mut() {
                if !p.is_static() {
                    p.velocity = (p.position - p.prev_position) / dt * (1.0 - damping);
                }
            }
        }
    }
}

/// [`relax_mesh_iterated`] with one iteration per sub-step.
fn relax_mesh(
    constraints: &mut [NeoHookeanConstraint],
    particles: &mut [SoftParticle],
    frame_dt: Real,
    substeps: usize,
    frames: usize,
    damping: Real,
) {
    relax_mesh_iterated(
        constraints,
        particles,
        frame_dt,
        substeps,
        frames,
        damping,
        1,
    );
}

/// [`relax_mesh`] for a single element.
fn relax(
    constraint: &mut NeoHookeanConstraint,
    particles: &mut [SoftParticle],
    frame_dt: Real,
    substeps: usize,
    frames: usize,
    damping: Real,
) {
    relax_mesh(
        std::slice::from_mut(constraint),
        particles,
        frame_dt,
        substeps,
        frames,
        damping,
    );
}

/// `max |F^T F - I|` entry of the deformation gradient: zero for a pure rotation.
fn deviation_from_rotation(constraint: &NeoHookeanConstraint, particles: &[SoftParticle]) -> Real {
    let f = constraint.deformation_gradient(particles);
    let mut worst: Real = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let dot: Real = (0..3).map(|k| f[k + 3 * i] * f[k + 3 * j]).sum();
            let identity = if i == j { 1.0 } else { 0.0 };
            worst = worst.max((dot - identity).abs());
        }
    }
    worst
}

/// Central-difference gradient of a scalar function of the four vertex positions.
fn numeric_gradient(
    particles: &[SoftParticle],
    value: impl Fn(&[SoftParticle]) -> Real,
) -> [Vec3; 4] {
    let h = 1e-6;
    let mut work = particles.to_vec();
    let mut gradient = [Vec3::zeros(); 4];
    for (k, grad) in gradient.iter_mut().enumerate() {
        for axis in 0..3 {
            let original = work[k].position[axis];
            work[k].position[axis] = original + h;
            let plus = value(&work);
            work[k].position[axis] = original - h;
            let minus = value(&work);
            work[k].position[axis] = original;
            grad[axis] = (plus - minus) / (2.0 * h);
        }
    }
    gradient
}

fn deformed_general_tet() -> Vec<SoftParticle> {
    let map = Mat3::new(1.2, 0.3, -0.1, 0.05, 0.8, 0.2, -0.15, 0.1, 1.1);
    deform(&general_tet(), &map, Vec3::new(2.0, -1.0, 0.5))
}

#[test]
fn deformation_gradient_reproduces_an_applied_linear_map() {
    let rest = general_tet();
    let constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    let map = Mat3::new(1.2, 0.3, -0.1, 0.05, 0.8, 0.2, -0.15, 0.1, 1.1);
    let deformed = deform(&rest, &map, Vec3::new(2.0, -1.0, 0.5));
    let f = constraint.deformation_gradient(&deformed);
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (f[i + 3 * j] - map[(i, j)]).abs() < 1e-12,
                "F[{i},{j}] = {} but the applied map has {}",
                f[i + 3 * j],
                map[(i, j)]
            );
        }
    }
}

#[test]
fn constraint_values_at_rest() {
    let rest = general_tet();
    let (mu, lambda) = lame(1e5, 0.3);
    let constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    let terms = constraint.evaluate_terms(&rest);
    assert!(
        (terms.c_dev - 3.0_f64.sqrt()).abs() < 1e-12,
        "C_D = sqrt(tr(F^T F)) = sqrt(3) at rest, got {}",
        terms.c_dev
    );
    let gamma = 1.0 + mu / (lambda + mu);
    assert!((constraint.hydrostatic_target() - gamma).abs() < 1e-12);
    assert!(
        (terms.c_hyd - (1.0 - gamma)).abs() < 1e-12,
        "C_H = det(F) - gamma = -mu/K at rest, got {}",
        terms.c_hyd
    );
    assert!(terms.dev_defined);
    assert!((constraint.rest_volume - signed_volume(&rest)).abs() < 1e-12);
    // The rest state is stress free and carries no energy.
    for entry in constraint.first_piola_stress(&rest) {
        assert!(entry.abs() < 1e-8 * mu, "rest stress {entry}");
    }
    assert!(constraint.elastic_energy(&rest).abs() < 1e-9);
}

#[test]
fn vertex_gradients_match_finite_differences() {
    let deformed = deformed_general_tet();
    let constraint = constraint_for(&general_tet(), 1e5, 0.3, 0.0);
    let terms = constraint.evaluate_terms(&deformed);
    let numeric_dev = numeric_gradient(&deformed, |p| constraint.evaluate_terms(p).c_dev);
    let numeric_hyd = numeric_gradient(&deformed, |p| constraint.evaluate_terms(p).c_hyd);
    for (k, (dev, hyd)) in numeric_dev.iter().zip(&numeric_hyd).enumerate() {
        assert!(
            (terms.grad_dev[k] - dev).norm() < 1e-7,
            "deviatoric gradient of vertex {k}: {:?} vs {dev:?}",
            terms.grad_dev[k]
        );
        assert!(
            (terms.grad_hyd[k] - hyd).norm() < 1e-7,
            "hydrostatic gradient of vertex {k}: {:?} vs {hyd:?}",
            terms.grad_hyd[k]
        );
    }
}

#[test]
fn gradients_are_translation_and_rotation_invariant() {
    let deformed = deformed_general_tet();
    let constraint = constraint_for(&general_tet(), 1e5, 0.3, 0.0);
    let terms = constraint.evaluate_terms(&deformed);
    for (name, gradient) in [
        ("deviatoric", &terms.grad_dev),
        ("hydrostatic", &terms.grad_hyd),
    ] {
        let force = gradient.iter().fold(Vec3::zeros(), |acc, g| acc + g);
        assert!(
            force.norm() < 1e-12,
            "{name} gradient is not translation invariant: {force:?}"
        );
        let torque = gradient
            .iter()
            .zip(&deformed)
            .fold(Vec3::zeros(), |acc, (g, p)| acc + p.position.cross(g));
        assert!(
            torque.norm() < 1e-12,
            "{name} gradient is not rotation invariant: {torque:?}"
        );
    }
}

#[test]
fn elastic_forces_equal_the_sum_of_constraint_forces() {
    // dW/dx_k = C_D / a_D * grad_k C_D + C_H / a_H * grad_k C_H with the
    // compliances a_D = 1 / (mu V0) and a_H = 1 / (K V0).
    let deformed = deformed_general_tet();
    let constraint = constraint_for(&general_tet(), 1e5, 0.3, 0.0);
    let terms = constraint.evaluate_terms(&deformed);
    let v0 = constraint.rest_volume;
    let shear_weight = constraint.mu * v0;
    let bulk_weight = constraint.hydrostatic_modulus() * v0;
    let numeric = numeric_gradient(&deformed, |p| constraint.elastic_energy(p));
    for (k, expected) in numeric.iter().enumerate() {
        let analytic = terms.grad_dev[k] * (shear_weight * terms.c_dev)
            + terms.grad_hyd[k] * (bulk_weight * terms.c_hyd);
        assert!(
            (analytic - expected).norm() < 1e-5 * analytic.norm().max(1.0),
            "vertex {k}: constraint forces {analytic:?} vs energy gradient {expected:?}"
        );
    }
}

#[test]
fn small_strain_stress_is_hookes_law() {
    let rest = unit_tet();
    let (mu, lambda) = lame(2e5, 0.3);
    let constraint = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, mu, lambda, 0.0);
    let eps = 1e-6;
    // Uniaxial strain: sigma_xx = (lambda + 2 mu) eps and sigma_yy = sigma_zz = lambda eps.
    let stretched = deform(
        &rest,
        &Mat3::new(1.0 + eps, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        Vec3::zeros(),
    );
    let p = constraint.first_piola_stress(&stretched);
    let scale = (lambda + 2.0 * mu) * eps;
    assert!(
        (p[0] - (lambda + 2.0 * mu) * eps).abs() < 1e-4 * scale,
        "P_xx = {}",
        p[0]
    );
    assert!(
        (p[4] - lambda * eps).abs() < 1e-4 * scale,
        "P_yy = {}",
        p[4]
    );
    assert!(
        (p[8] - lambda * eps).abs() < 1e-4 * scale,
        "P_zz = {}",
        p[8]
    );
    // Simple shear: sigma_xy = sigma_yx = mu gamma.
    let gamma = 1e-6;
    let sheared = deform(
        &rest,
        &Mat3::new(1.0, gamma, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        Vec3::zeros(),
    );
    let p = constraint.first_piola_stress(&sheared);
    assert!(
        (p[3] - mu * gamma).abs() < 1e-4 * mu * gamma,
        "P_xy = {}",
        p[3]
    );
    assert!(
        (p[1] - mu * gamma).abs() < 1e-4 * mu * gamma,
        "P_yx = {}",
        p[1]
    );
    // The energy of a uniaxial strain is exactly V0 (lambda + 2 mu) eps^2 / 2.
    let strain = 1e-4;
    let stretched = deform(
        &rest,
        &Mat3::new(1.0 + strain, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        Vec3::zeros(),
    );
    let energy = constraint.elastic_energy(&stretched);
    let expected = 0.5 * constraint.rest_volume * (lambda + 2.0 * mu) * strain * strain;
    assert!(
        (energy - expected).abs() < 1e-6 * expected,
        "W = {energy} vs {expected}"
    );
}

#[test]
fn rest_state_is_a_fixed_point_for_every_poisson_ratio_and_time_step() {
    for &(youngs, poisson) in &[
        (1e3, 0.3),
        (1e5, 0.45),
        (1e7, 0.0),
        (1e9, 0.49),
        (1e5, -0.4),
        (1e6, 0.1),
    ] {
        for &dt in &[1e-4, 1e-2, 1.0 / 60.0] {
            let rest = general_tet();
            let mut particles = rest.clone();
            let mut constraint = constraint_for(&rest, youngs, poisson, 0.0);
            for _ in 0..30 {
                constraint.project(&mut particles, dt);
            }
            for (p, q) in particles.iter().zip(&rest) {
                assert!(
                    (p.position - q.position).norm() < 1e-8,
                    "E = {youngs}, nu = {poisson}, dt = {dt}: rest state drifted by {}",
                    (p.position - q.position).norm()
                );
            }
        }
    }
}

#[test]
fn x_compressed_tet_recovers_rest_volume_with_unchanged_centroid() {
    let rest = unit_tet();
    let mut particles = rest.clone();
    for p in &mut particles {
        p.position.x *= 0.5;
    }
    let rest_volume = signed_volume(&rest);
    assert!((signed_volume(&particles) / rest_volume - 0.5).abs() < 1e-12);
    let centroid = mass_weighted_centroid(&particles);
    let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    let initial_energy = constraint.elastic_energy(&particles);
    assert!(initial_energy > 0.0);
    relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    let ratio = signed_volume(&particles) / rest_volume;
    assert!(
        (ratio - 1.0).abs() < 0.01,
        "V / V0 = {ratio} should return to 1"
    );
    assert!(
        (mass_weighted_centroid(&particles) - centroid).norm() < 1e-9,
        "the centroid must not move"
    );
    assert!(
        deviation_from_rotation(&constraint, &particles) < 1e-6,
        "the shape must return to a rotated copy of the rest shape"
    );
    assert!(constraint.elastic_energy(&particles) < 1e-9 * initial_energy);
}

#[test]
fn stretched_sheared_and_general_tets_recover_the_rest_shape() {
    let cases = [
        (
            "stretched",
            unit_tet(),
            Mat3::new(1.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        ),
        (
            "sheared",
            unit_tet(),
            Mat3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.7, 0.0, 0.0, 1.0),
        ),
        (
            "general tet, anisotropic",
            general_tet(),
            Mat3::new(0.6, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.2),
        ),
    ];
    for (name, rest, map) in cases {
        let mut particles = deform(&rest, &map, Vec3::zeros());
        let centroid = mass_weighted_centroid(&particles);
        let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
        relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
        let ratio = signed_volume(&particles) / signed_volume(&rest);
        assert!((ratio - 1.0).abs() < 0.01, "{name}: V / V0 = {ratio}");
        assert!(
            deviation_from_rotation(&constraint, &particles) < 1e-6,
            "{name}"
        );
        assert!(
            (mass_weighted_centroid(&particles) - centroid).norm() < 1e-9,
            "{name}"
        );
    }
}

#[test]
fn inverted_tet_is_untangled() {
    let rest = unit_tet();
    let mut particles = rest.clone();
    particles[3].position.z = -1.0;
    assert!(
        signed_volume(&particles) < 0.0,
        "the test tet must start inverted"
    );
    let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    let ratio = signed_volume(&particles) / signed_volume(&rest);
    assert!(
        (ratio - 1.0).abs() < 0.01,
        "an inverted tet must recover V0, got V / V0 = {ratio}"
    );
    assert!(deviation_from_rotation(&constraint, &particles) < 1e-6);
}

#[test]
fn low_poisson_ratio_materials_stay_stable() {
    // lambda = 0 (nu = 0) and auxetic lambda < 0 must neither blow up nor collapse.
    for poisson in [0.0, -0.3] {
        let rest = unit_tet();
        let mut particles = rest.clone();
        for p in &mut particles {
            p.position.x *= 0.9;
        }
        let mut constraint = constraint_for(&rest, 1e5, poisson, 0.0);
        relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
        let ratio = signed_volume(&particles) / signed_volume(&rest);
        assert!(
            (ratio - 1.0).abs() < 0.01,
            "nu = {poisson}: V / V0 = {ratio}"
        );
    }
}

#[test]
fn a_single_projection_conserves_momentum() {
    let rest = general_tet();
    let masses = [1.0, 2.0, 3.0, 4.0];
    let mut particles: Vec<SoftParticle> = deformed_general_tet()
        .iter()
        .zip(masses)
        .map(|(p, m)| SoftParticle::new(p.position, m))
        .collect();
    let before: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
    let mut constraint = constraint_for(&rest, 1e4, 0.3, 0.0);
    constraint.project(&mut particles, 1.0 / 600.0);
    let momentum = particles
        .iter()
        .zip(&before)
        .fold(Vec3::zeros(), |acc, (p, b)| {
            acc + (p.position - b) / p.inverse_mass
        });
    assert!(
        momentum.norm() < 1e-9,
        "mass-weighted displacement {momentum:?}"
    );
    let moved: Real = particles
        .iter()
        .zip(&before)
        .map(|(p, b)| (p.position - b).norm())
        .sum();
    assert!(moved > 1e-6, "a deformed tet must be corrected");
}

#[test]
fn pinned_vertex_stays_put_while_the_rest_recovers() {
    let rest = unit_tet();
    let mut particles = rest.clone();
    particles[0] = SoftParticle::new_static(rest[0].position);
    for p in &mut particles {
        p.position.x *= 0.5;
    }
    let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    assert_eq!(
        particles[0].position, rest[0].position,
        "the pinned vertex moved"
    );
    let ratio = signed_volume(&particles) / signed_volume(&rest);
    assert!((ratio - 1.0).abs() < 0.01, "V / V0 = {ratio}");
}

#[test]
fn volume_only_material_preserves_volume() {
    // mu = 0 leaves just the hydrostatic constraint with gamma = 1.
    let rest = unit_tet();
    let mut particles = rest.clone();
    for p in &mut particles {
        p.position.x *= 0.5;
    }
    let mut constraint = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, 0.0, 1e5, 0.0);
    assert!((constraint.hydrostatic_target() - 1.0).abs() < 1e-15);
    let centroid = mass_weighted_centroid(&particles);
    relax(&mut constraint, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    let ratio = signed_volume(&particles) / signed_volume(&rest);
    assert!((ratio - 1.0).abs() < 0.01, "V / V0 = {ratio}");
    assert!((mass_weighted_centroid(&particles) - centroid).norm() < 1e-9);
}

#[test]
fn compliance_softens_the_volume_response() {
    // Pure volume material (mu = 0, lambda = 1e7): one projection of a 20 %
    // compressed tet.  The extra compliance is added in series to 1 / (K V0).
    let rest = unit_tet();
    let rest_volume = signed_volume(&rest);
    let volume_error_after_one_step = |compliance: Real| -> Real {
        let mut particles = rest.clone();
        for p in &mut particles {
            p.position.x *= 0.8;
        }
        let mut constraint =
            NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, 0.0, 1e7, compliance);
        constraint.project(&mut particles, 1.0 / 600.0);
        (signed_volume(&particles) / rest_volume - 1.0).abs()
    };
    let stiff = volume_error_after_one_step(0.0);
    let medium = volume_error_after_one_step(1e-5);
    let soft = volume_error_after_one_step(1e-3);
    assert!(
        stiff < 0.01,
        "the stiff constraint must nearly restore the volume, error {stiff}"
    );
    assert!(
        medium > 5.0 * stiff && medium < 0.15,
        "a medium compliance corrects only part of the error, error {medium}"
    );
    assert!(
        soft > 0.15 && soft < 0.2 - 1e-6,
        "a large compliance barely moves but still corrects towards V0, error {soft}"
    );
    assert!(
        stiff < medium && medium < soft,
        "{stiff} < {medium} < {soft}"
    );
}

#[test]
fn multipliers_accumulate_and_reset() {
    let mut particles = deformed_general_tet();
    let mut constraint = constraint_for(&general_tet(), 1e5, 0.3, 0.0);
    assert_eq!(constraint.lambda_dev, 0.0);
    assert_eq!(constraint.lambda_vol, 0.0);
    constraint.project(&mut particles, 1.0 / 600.0);
    assert!(constraint.lambda_dev != 0.0 && constraint.lambda_vol != 0.0);
    let boxed: &mut dyn SoftConstraint = &mut constraint;
    boxed.reset_lambda();
    assert_eq!(constraint.lambda_dev, 0.0);
    assert_eq!(constraint.lambda_vol, 0.0);
}

#[test]
fn invalid_input_leaves_the_particles_untouched() {
    let rest = general_tet();
    let deformed = deformed_general_tet();
    let snapshot: Vec<Vec3> = deformed.iter().map(|p| p.position).collect();
    let unchanged = |particles: &[SoftParticle]| {
        particles
            .iter()
            .zip(&snapshot)
            .all(|(p, s)| p.position == *s)
    };

    let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    for bad_dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut particles = deformed.clone();
        constraint.project(&mut particles, bad_dt);
        assert!(unchanged(&particles), "dt = {bad_dt}");
    }

    // Index outside the particle array.
    let mut out_of_range = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, 1e4, 1e4, 0.0);
    out_of_range.indices = [0, 1, 2, 99];
    let mut particles = deformed.clone();
    out_of_range.project(&mut particles, 1.0 / 600.0);
    assert!(unchanged(&particles));

    // Degenerate (flat) rest shape has no volume to preserve.
    let flat = tet([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ]);
    let mut degenerate = constraint_for(&flat, 1e5, 0.3, 0.0);
    assert_eq!(degenerate.rest_volume, 0.0);
    let mut particles = deformed.clone();
    degenerate.project(&mut particles, 1.0 / 600.0);
    assert!(unchanged(&particles));

    // No stiffness at all: nothing to enforce.
    let mut soft = NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, 0.0, 0.0, 0.0);
    let mut particles = deformed.clone();
    soft.project(&mut particles, 1.0 / 600.0);
    assert!(unchanged(&particles));
}

#[test]
fn fully_collapsed_element_stays_finite() {
    let rest = unit_tet();
    let mut constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    let mut particles = tet([[0.4, 0.4, 0.4]; 4]);
    for _ in 0..10 {
        constraint.project(&mut particles, 1.0 / 600.0);
    }
    for p in &particles {
        assert!(p.position.iter().all(|c| c.is_finite()), "{:?}", p.position);
    }
}

#[test]
fn unusual_parameter_combinations_stay_finite() {
    let rest = unit_tet();
    // K = lambda + mu = 0 disables the hydrostatic constraint; mu = 0 disables the other.
    for (mu, lambda) in [(1e4, -1e4), (1e4, -2e4), (0.0, -5.0), (1e4, 0.0)] {
        let mut constraint =
            NeoHookeanConstraint::from_particles([0, 1, 2, 3], &rest, mu, lambda, 0.0);
        let mut particles = deformed_general_tet();
        for _ in 0..5 {
            constraint.project(&mut particles, 1.0 / 600.0);
        }
        for p in &particles {
            assert!(
                p.position.iter().all(|c| c.is_finite()),
                "mu = {mu}, lambda = {lambda}"
            );
        }
    }
}

#[test]
fn terms_report_an_undefined_deviatoric_gradient_for_a_collapsed_element() {
    let rest = unit_tet();
    let constraint = constraint_for(&rest, 1e5, 0.3, 0.0);
    let collapsed = tet([[1.0, 2.0, 3.0]; 4]);
    let terms: NeoHookeanTerms = constraint.evaluate_terms(&collapsed);
    assert!(!terms.dev_defined);
    assert_eq!(terms.c_dev, 0.0);
    assert!(terms.grad_dev.iter().all(|g| g.norm() == 0.0));
}

/// The unit cube split into five tetrahedra that share their corner vertices:
/// the central element `(1, 3, 4, 6)` and the four corner elements.  The
/// elements do not all have the same orientation, so volumes are compared as
/// ratios of signed volumes.
fn five_tet_cube() -> (Vec<SoftParticle>, Vec<[usize; 4]>) {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let particles = corners
        .iter()
        .map(|c| SoftParticle::new(Vec3::new(c[0], c[1], c[2]), 1.0))
        .collect();
    let elements = vec![
        [1, 3, 4, 6],
        [0, 1, 3, 4],
        [2, 1, 3, 6],
        [5, 1, 4, 6],
        [7, 3, 4, 6],
    ];
    (particles, elements)
}

fn cube_constraints(
    rest: &[SoftParticle],
    elements: &[[usize; 4]],
    youngs: Real,
    poisson: Real,
) -> Vec<NeoHookeanConstraint> {
    let (mu, lambda) = lame(youngs, poisson);
    elements
        .iter()
        .map(|&element| NeoHookeanConstraint::from_particles(element, rest, mu, lambda, 0.0))
        .collect()
}

/// Signed volume of the tetrahedron `element` of `particles`.
fn element_volume(particles: &[SoftParticle], element: &[usize; 4]) -> Real {
    let corner = |k: usize| particles[element[k]].position;
    let e1 = corner(1) - corner(0);
    let e2 = corner(2) - corner(0);
    let e3 = corner(3) - corner(0);
    e1.dot(&e2.cross(&e3)) / 6.0
}

fn total_energy(constraints: &[NeoHookeanConstraint], particles: &[SoftParticle]) -> Real {
    constraints
        .iter()
        .map(|c| c.elastic_energy(particles))
        .sum()
}

#[test]
fn compressed_five_tet_cube_recovers_its_shape_across_shared_vertices() {
    let (rest, elements) = five_tet_cube();
    let mut constraints = cube_constraints(&rest, &elements, 1e5, 0.3);
    let mut particles = rest.clone();
    for p in &mut particles {
        p.position.x *= 0.6;
    }
    let centroid = mass_weighted_centroid(&particles);
    let initial_energy = total_energy(&constraints, &particles);
    assert!(initial_energy > 0.0);
    relax_mesh(&mut constraints, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    for (k, element) in elements.iter().enumerate() {
        let ratio = element_volume(&particles, element) / element_volume(&rest, element);
        assert!(
            (ratio - 1.0).abs() < 0.01,
            "element {k}: V / V0 = {ratio} should return to 1"
        );
        assert!(
            deviation_from_rotation(&constraints[k], &particles) < 1e-6,
            "element {k} must be a rotated copy of its rest shape"
        );
    }
    assert!(
        (mass_weighted_centroid(&particles) - centroid).norm() < 1e-9,
        "internal forces must not move the centroid"
    );
    assert!(total_energy(&constraints, &particles) < 1e-9 * initial_energy);
}

#[test]
fn stretched_and_sheared_five_tet_cubes_recover_their_shape() {
    let cases = [
        (
            "stretched",
            Mat3::new(1.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        ),
        (
            "sheared",
            Mat3::new(1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        ),
        (
            "compressed",
            Mat3::new(0.7, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.7),
        ),
    ];
    for (name, map) in cases {
        let (rest, elements) = five_tet_cube();
        let mut constraints = cube_constraints(&rest, &elements, 1e5, 0.3);
        let mut particles = deform(&rest, &map, Vec3::zeros());
        relax_mesh(&mut constraints, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
        for (k, element) in elements.iter().enumerate() {
            let ratio = element_volume(&particles, element) / element_volume(&rest, element);
            assert!(
                (ratio - 1.0).abs() < 0.01,
                "{name}, element {k}: V / V0 = {ratio}"
            );
            assert!(
                deviation_from_rotation(&constraints[k], &particles) < 1e-6,
                "{name}, element {k}"
            );
        }
    }
}

#[test]
fn five_tet_cube_with_a_pinned_face_springs_back_to_the_rest_shape() {
    let (rest, elements) = five_tet_cube();
    let mut constraints = cube_constraints(&rest, &elements, 1e5, 0.3);
    let mut particles = rest.clone();
    // The x = 0 face is clamped; the other vertices are pushed towards it.
    let clamped = [0, 3, 4, 7];
    for (k, p) in particles.iter_mut().enumerate() {
        if clamped.contains(&k) {
            *p = SoftParticle::new_static(rest[k].position);
        } else {
            p.position.x *= 0.6;
        }
    }
    relax_mesh(&mut constraints, &mut particles, 1.0 / 60.0, 20, 300, 0.05);
    for (k, (p, q)) in particles.iter().zip(&rest).enumerate() {
        assert!(
            (p.position - q.position).norm() < 1e-3,
            "vertex {k} ended {} away from the rest position",
            (p.position - q.position).norm()
        );
    }
    for &k in &clamped {
        assert_eq!(
            particles[k].position, rest[k].position,
            "clamped vertex {k}"
        );
    }
}

#[test]
fn rigidly_moved_five_tet_cube_is_a_fixed_point() {
    // A rotated and translated copy of the rest shape is stress free as well.
    let (rest, elements) = five_tet_cube();
    let rotation = yaw_pitch_rotation(0.7, 0.4);
    let mut particles = deform(&rest, &rotation, Vec3::new(3.0, -2.0, 5.0));
    let snapshot: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
    let mut constraints = cube_constraints(&rest, &elements, 1e6, 0.4);
    for _ in 0..30 {
        for constraint in constraints.iter_mut() {
            constraint.reset_lambda();
        }
        for constraint in constraints.iter_mut() {
            constraint.project(&mut particles, 1.0 / 600.0);
        }
    }
    for (p, s) in particles.iter().zip(&snapshot) {
        assert!(
            (p.position - s).norm() < 1e-9,
            "a rigid motion must not be corrected, moved by {}",
            (p.position - s).norm()
        );
    }
    assert!(total_energy(&constraints, &particles) < 1e-9);
}

/// Rotation by `yaw` about `z` followed by `pitch` about `x`.
fn yaw_pitch_rotation(yaw: Real, pitch: Real) -> Mat3 {
    let (sy, cy) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    let about_z = Mat3::new(cy, -sy, 0.0, sy, cy, 0.0, 0.0, 0.0, 1.0);
    let about_x = Mat3::new(1.0, 0.0, 0.0, 0.0, cp, -sp, 0.0, sp, cp);
    about_z * about_x
}

/// Deterministic 64-bit linear congruential generator for randomised tests.
struct Lcg(u64);

impl Lcg {
    /// Uniform sample in `[0, 1)`.
    fn next_unit(&mut self) -> Real {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as Real / (1_u64 << 53) as Real
    }
}

/// Energy scale `V0 (mu + K)` against which energy changes are measured.
fn energy_scale(constraint: &NeoHookeanConstraint) -> Real {
    constraint.rest_volume * (constraint.mu.max(0.0) + constraint.hydrostatic_modulus().max(0.0))
}

#[test]
fn stiff_elements_relax_monotonically_over_long_time_steps() {
    // One full Gauss-Newton step of a stiff element over a long sub-step
    // overshoots: without step-length control a single projection of a tet
    // compressed to half its volume takes it down to a few percent of its
    // volume and further projections oscillate wildly.  Every projection must
    // instead lower the elastic energy, and the element must reach its rest
    // shape.
    let cases = [
        (1e5, 0.3, 1.0 / 60.0, 300),
        (1e6, 0.45, 1.0 / 60.0, 1500),
        (1e6, 0.3, 1.0 / 240.0, 300),
        (1e3, 0.3, 1.0 / 60.0, 300),
    ];
    for (youngs, poisson, dt, iterations) in cases {
        let rest = unit_tet();
        let mut particles = rest.clone();
        for p in &mut particles {
            p.position.x *= 0.5;
        }
        let centroid = mass_weighted_centroid(&particles);
        let mut constraint = constraint_for(&rest, youngs, poisson, 0.0);
        let tolerance = 1e-9 * energy_scale(&constraint);
        let mut previous = constraint.elastic_energy(&particles);
        for step in 0..iterations {
            constraint.reset_lambda();
            constraint.project(&mut particles, dt);
            let energy = constraint.elastic_energy(&particles);
            assert!(
                energy <= previous + tolerance,
                "E = {youngs}, nu = {poisson}, dt = {dt}, step {step}: energy rose from {previous} to {energy}"
            );
            previous = energy;
        }
        let ratio = signed_volume(&particles) / signed_volume(&rest);
        assert!(
            (ratio - 1.0).abs() < 0.01,
            "E = {youngs}, nu = {poisson}, dt = {dt}: V / V0 = {ratio}"
        );
        assert!(
            deviation_from_rotation(&constraint, &particles) < 1e-3,
            "E = {youngs}, nu = {poisson}, dt = {dt}: the shape did not recover"
        );
        assert!(
            (mass_weighted_centroid(&particles) - centroid).norm() < 1e-9,
            "E = {youngs}, nu = {poisson}, dt = {dt}: the centroid moved"
        );
    }
}

#[test]
fn projection_never_raises_the_elastic_energy_for_random_materials_and_deformations() {
    // Random stiffness (100 Pa to 100 MPa), Poisson ratio (-0.45 to 0.49),
    // sub-step (0.1 ms to 0.1 s) and deformation, including inverted tets.
    let rest = general_tet();
    let mut rng = Lcg(0x5EED_CAFE);
    for trial in 0..200 {
        let youngs = 10.0_f64.powf(2.0 + 6.0 * rng.next_unit());
        let poisson = -0.45 + 0.94 * rng.next_unit();
        let dt = 10.0_f64.powf(-4.0 + 3.0 * rng.next_unit());
        let map = Mat3::from_fn(|i, j| {
            let identity = if i == j { 1.0 } else { 0.0 };
            identity + 1.8 * (rng.next_unit() - 0.5)
        });
        let mut particles = deform(&rest, &map, Vec3::zeros());
        let mut constraint = constraint_for(&rest, youngs, poisson, 0.0);
        let tolerance = 1e-9 * energy_scale(&constraint);
        let mut previous = constraint.elastic_energy(&particles);
        for step in 0..60 {
            constraint.reset_lambda();
            constraint.project(&mut particles, dt);
            assert!(
                particles
                    .iter()
                    .all(|p| p.position.iter().all(|c| c.is_finite())),
                "trial {trial}, step {step}: non-finite position"
            );
            let energy = constraint.elastic_energy(&particles);
            assert!(
                energy <= previous + tolerance,
                "trial {trial} (E = {youngs:e}, nu = {poisson}, dt = {dt:e}), step {step}: \
                 energy rose from {previous} to {energy}"
            );
            previous = energy;
        }
    }
}

#[test]
fn well_conditioned_steps_are_applied_in_full() {
    // For a soft enough element over a short sub-step the Gauss-Newton step
    // already decreases the potential, so the step-length control must leave it
    // untouched: the multipliers are the exact solution of the coupled 2 x 2
    // system and the vertices move by `w grad C . dl`.
    let rest = general_tet();
    let deformed = deformed_general_tet();
    let dt = 1.0 / 600.0;
    let mut constraint = constraint_for(&rest, 1e4, 0.3, 0.0);
    let terms = constraint.evaluate_terms(&deformed);
    let alpha_dev = 1.0 / (constraint.mu * constraint.rest_volume * dt * dt);
    let alpha_hyd = 1.0 / (constraint.hydrostatic_modulus() * constraint.rest_volume * dt * dt);
    let gram =
        |a: &[Vec3; 4], b: &[Vec3; 4]| -> Real { a.iter().zip(b).map(|(x, y)| x.dot(y)).sum() };
    let s_dd = gram(&terms.grad_dev, &terms.grad_dev) + alpha_dev;
    let s_hh = gram(&terms.grad_hyd, &terms.grad_hyd) + alpha_hyd;
    let s_dh = gram(&terms.grad_dev, &terms.grad_hyd);
    let (rhs_dev, rhs_hyd) = (-terms.c_dev, -terms.c_hyd);
    let determinant = s_dd * s_hh - s_dh * s_dh;
    let delta_dev = (rhs_dev * s_hh - rhs_hyd * s_dh) / determinant;
    let delta_hyd = (s_dd * rhs_hyd - s_dh * rhs_dev) / determinant;

    let mut particles = deformed.clone();
    constraint.project(&mut particles, dt);
    assert!(
        (constraint.lambda_dev - delta_dev).abs() <= 1e-12 * delta_dev.abs(),
        "deviatoric multiplier {} vs {delta_dev}",
        constraint.lambda_dev
    );
    assert!(
        (constraint.lambda_vol - delta_hyd).abs() <= 1e-12 * delta_hyd.abs(),
        "hydrostatic multiplier {} vs {delta_hyd}",
        constraint.lambda_vol
    );
    for (k, ((moved, start), (grad_dev, grad_hyd))) in particles
        .iter()
        .zip(&deformed)
        .zip(terms.grad_dev.iter().zip(&terms.grad_hyd))
        .enumerate()
    {
        let expected = start.position + grad_dev * delta_dev + grad_hyd * delta_hyd;
        assert!(
            (moved.position - expected).norm() < 1e-12,
            "vertex {k} was not moved by the full step"
        );
    }
}

#[test]
fn several_iterations_make_a_stiff_one_substep_scheme_converge() {
    // E = 1 MPa with a single 1/60 s sub-step is far stiffer than one
    // Gauss-Newton iteration per sub-step resolves; with the step-length
    // control more iterations per sub-step converge to the implicit-Euler
    // solution, so the compressed tet still springs back to its rest shape.
    let rest = unit_tet();
    let mut particles = rest.clone();
    for p in &mut particles {
        p.position.x *= 0.5;
    }
    let centroid = mass_weighted_centroid(&particles);
    let mut constraints = vec![constraint_for(&rest, 1e6, 0.3, 0.0)];
    relax_mesh_iterated(
        &mut constraints,
        &mut particles,
        1.0 / 60.0,
        1,
        150,
        0.05,
        20,
    );
    let ratio = signed_volume(&particles) / signed_volume(&rest);
    assert!((ratio - 1.0).abs() < 0.01, "V / V0 = {ratio}");
    assert!(deviation_from_rotation(&constraints[0], &particles) < 1e-3);
    assert!((mass_weighted_centroid(&particles) - centroid).norm() < 1e-9);
}
