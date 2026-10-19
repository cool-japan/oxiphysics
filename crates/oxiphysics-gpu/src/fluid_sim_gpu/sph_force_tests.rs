// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the SPH force densities.
//!
//! `sph_compute_forces` used to add the pressure term as an *acceleration*
//! (`-Σ m_j (p_i/ρ_i² + p_j/ρ_j²) ∇W`) to the viscosity and gravity *force
//! densities*, while `sph_integrate` divides the sum by `ρ_i`: the pressure
//! was a factor `ρ_i` too weak.  The surface tension promised by the doc
//! comment and by `SphConfig::surface_tension` was not implemented at all.
//! The tests below pin the unit of every term against closed-form pair
//! interactions and check the Akinci et al. (2013) surface tension.

use super::*;
use std::f64::consts::PI;

/// Reference configuration with every force switched off except what a test
/// enables explicitly.
fn quiet_config() -> SphConfig {
    SphConfig {
        h: 0.1,
        rest_density: 1000.0,
        pressure_k: 100.0,
        viscosity: 0.0,
        surface_tension: 0.0,
        gravity: [0.0; 3],
        dt: 1.0e-3,
    }
}

/// A particle with explicitly chosen state (density and pressure are *not*
/// recomputed, so the analytic expectations below hold exactly).
fn particle(
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
    density: f64,
    pressure: f64,
) -> SphParticle {
    SphParticle {
        position,
        velocity,
        mass,
        density,
        pressure,
        force: [0.0; 3],
    }
}

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: {actual} vs {expected} (tol {tol})"
    );
}

/// `|∇W_spiky(r)| = 45 / (π h⁶) (h - r)²`.
fn spiky_gradient_magnitude(r: f64, h: f64) -> f64 {
    45.0 / (PI * h.powi(6)) * (h - r).powi(2)
}

/// `∇²W_viscosity(r) = 45 / (π h⁶) (h - r)`.
fn viscosity_laplacian_value(r: f64, h: f64) -> f64 {
    45.0 / (PI * h.powi(6)) * (h - r)
}

/// `|W'_poly6(r)| = 945 / (32 π h⁹) (h² - r²)² r`.
fn poly6_slope_magnitude(r: f64, h: f64) -> f64 {
    945.0 / (32.0 * PI * h.powi(9)) * (h * h - r * r).powi(2) * r
}

/// Akinci et al. cohesion spline written out in units of `h`
/// (`C = ĉ(r/h) / h³`).
fn cohesion_reference(r: f64, h: f64) -> f64 {
    let q = r / h;
    let shape = if q > 1.0 || q <= 0.0 {
        0.0
    } else if q > 0.5 {
        (1.0 - q).powi(3) * q.powi(3)
    } else {
        2.0 * (1.0 - q).powi(3) * q.powi(3) - 1.0 / 64.0
    };
    32.0 / PI * shape / h.powi(3)
}

// ── pressure ─────────────────────────────────────────────────────────────────

/// Two particles at `r = h/2` with positive pressure repel each other with the
/// analytic SPH pressure acceleration: the force is `ρ a` and, after
/// `sph_integrate`, the velocity change is exactly `dt a`.
#[test]
fn pressure_pair_velocity_change_matches_the_analytic_acceleration() {
    let config = quiet_config();
    let (h, dt) = (config.h, config.dt);
    let r = 0.5 * h;
    let (mass, density, pressure) = (0.02, 900.0, 150.0);
    let mut particles = vec![
        particle([0.0; 3], [0.0; 3], mass, density, pressure),
        particle([r, 0.0, 0.0], [0.0; 3], mass, density, pressure),
    ];

    sph_compute_forces(&mut particles, &config);
    let accel = mass * (2.0 * pressure / (density * density)) * spiky_gradient_magnitude(r, h);
    assert!(accel > 0.1, "test setup: acceleration {accel}");
    // Force density, not acceleration.
    assert_close(
        particles[0].force[0],
        -density * accel,
        1e-12 * density * accel,
        "f0",
    );
    assert_close(
        particles[1].force[0],
        density * accel,
        1e-12 * density * accel,
        "f1",
    );

    sph_integrate(&mut particles, &config);
    assert_close(
        particles[0].velocity[0],
        -dt * accel,
        1e-12 * dt * accel,
        "v0",
    );
    assert_close(
        particles[1].velocity[0],
        dt * accel,
        1e-12 * dt * accel,
        "v1",
    );
    for p in &particles {
        assert_eq!(p.velocity[1], 0.0);
        assert_eq!(p.velocity[2], 0.0);
    }
    // Semi-implicit Euler: the position uses the updated velocity.
    assert_close(
        particles[0].position[0],
        -dt * dt * accel,
        1e-12 * dt * dt * accel,
        "x0",
    );
    assert_close(particles[1].position[0], r + dt * dt * accel, 1e-12, "x1");
}

/// The pressure acceleration is independent of the particle's own density once
/// `p/ρ²` is fixed: the force density scales with `ρ_i` and `sph_integrate`
/// divides it out again.
#[test]
fn pressure_acceleration_does_not_depend_on_the_force_density_scale() {
    let config = quiet_config();
    let (h, dt) = (config.h, config.dt);
    let r = 0.4 * h;
    let mass = 0.03;
    let mut expected = None;
    for (density, pressure) in [(500.0, 40.0), (1000.0, 160.0), (2000.0, 640.0)] {
        // p / ρ² = 1.6e-4 in all three cases.
        let mut particles = vec![
            particle([0.0; 3], [0.0; 3], mass, density, pressure),
            particle([r, 0.0, 0.0], [0.0; 3], mass, density, pressure),
        ];
        sph_compute_forces(&mut particles, &config);
        sph_integrate(&mut particles, &config);
        let dv = particles[0].velocity[0];
        assert!(dv < 0.0);
        match expected {
            None => expected = Some(dv),
            Some(first) => assert_close(dv, first, 1e-12 * dt, "velocity change"),
        }
    }
}

/// Pressure, viscosity and gravity of an asymmetric pair (different masses,
/// densities, pressures and velocities): each term must enter `sph_integrate`
/// as a force density.
#[test]
fn pressure_viscosity_and_gravity_share_one_unit() {
    let config = SphConfig {
        viscosity: 0.8,
        gravity: [0.0, -9.81, 0.3],
        ..quiet_config()
    };
    let (h, dt) = (config.h, config.dt);
    let r = 0.6 * h;
    let (m0, m1) = (0.02, 0.03);
    let (rho0, rho1) = (800.0, 1100.0);
    let (p0, p1) = (120.0, 200.0);
    let (v0, v1) = ([0.0, 0.0, 0.0], [0.4, 0.1, -0.2]);
    let mut particles = vec![
        particle([0.0; 3], v0, m0, rho0, p0),
        particle([r, 0.0, 0.0], v1, m1, rho1, p1),
    ];

    sph_compute_forces(&mut particles, &config);
    sph_integrate(&mut particles, &config);

    let pressure_factor = p0 / (rho0 * rho0) + p1 / (rho1 * rho1);
    let grad = spiky_gradient_magnitude(r, h);
    let lap = viscosity_laplacian_value(r, h);
    // Particle 0 is pushed towards -x by the pressure, particle 1 towards +x.
    let mut expected0 = [-m1 * pressure_factor * grad, 0.0, 0.0];
    let mut expected1 = [m0 * pressure_factor * grad, 0.0, 0.0];
    for c in 0..3 {
        // Viscosity force density mu m_j (v_j - v_i) lap / rho_j, divided by rho_i.
        expected0[c] += config.viscosity * m1 * (v1[c] - v0[c]) * lap / rho1 / rho0;
        expected1[c] += config.viscosity * m0 * (v0[c] - v1[c]) * lap / rho0 / rho1;
        expected0[c] += config.gravity[c];
        expected1[c] += config.gravity[c];
    }
    for c in 0..3 {
        assert_close(
            particles[0].velocity[c],
            v0[c] + dt * expected0[c],
            1e-13,
            &format!("particle 0 component {c}"),
        );
        assert_close(
            particles[1].velocity[c],
            v1[c] + dt * expected1[c],
            1e-13,
            &format!("particle 1 component {c}"),
        );
    }
}

#[test]
fn particles_beyond_the_support_radius_do_not_interact() {
    let config = SphConfig {
        viscosity: 1.0,
        surface_tension: 1.0,
        ..quiet_config()
    };
    let h = config.h;
    for separation in [1.0 * h, 1.2 * h, 5.0 * h] {
        let mut particles = vec![
            particle([0.0; 3], [1.0, 0.0, 0.0], 0.05, 1000.0, 500.0),
            particle(
                [separation, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                0.05,
                1000.0,
                500.0,
            ),
        ];
        sph_compute_forces(&mut particles, &config);
        for p in &particles {
            assert_eq!(p.force, [0.0; 3], "separation {separation}");
        }
    }
}

#[test]
fn isolated_particle_feels_only_gravity() {
    let config = SphConfig {
        gravity: [0.0, -9.81, 0.0],
        surface_tension: 0.5,
        viscosity: 0.3,
        ..quiet_config()
    };
    let mut particles = vec![particle(
        [0.3, 0.2, 0.1],
        [0.5, 0.0, 0.0],
        0.04,
        640.0,
        90.0,
    )];
    sph_compute_forces(&mut particles, &config);
    assert_eq!(particles[0].force, [0.0, -9.81 * 640.0, 0.0]);
    assert_eq!(sph_surface_normals(&particles, &config), vec![[0.0; 3]]);
}

// ── cohesion spline ──────────────────────────────────────────────────────────

#[test]
fn cohesion_kernel_matches_the_published_spline() {
    let h = 0.1;
    for q in [0.05, 0.1, 0.2, 0.3, 0.45, 0.5, 0.55, 0.75, 0.9, 1.0] {
        let r = q * h;
        let got = SphKernels::cohesion(r, h);
        let expected = cohesion_reference(r, h);
        assert_close(
            got,
            expected,
            1e-12 * expected.abs().max(1.0),
            &format!("C({q} h)"),
        );
    }
    // Hand-evaluated values: C(0.75 h) = 32/π (0.25)³ (0.75)³ / h³ and
    // C(0.25 h) = 32/π (2 (0.75)³ (0.25)³ - 1/64) / h³.
    let c_outer = 32.0 / PI * 0.015625 * 0.421875 / h.powi(3);
    assert_close(
        SphKernels::cohesion(0.75 * h, h),
        c_outer,
        1e-9 * c_outer,
        "outer branch",
    );
    let c_inner = 32.0 / PI * (2.0 * 0.421875 * 0.015625 - 1.0 / 64.0) / h.powi(3);
    assert_close(
        SphKernels::cohesion(0.25 * h, h),
        c_inner,
        1e-9 * c_inner.abs(),
        "inner branch",
    );
}

#[test]
fn cohesion_kernel_is_continuous_attractive_in_range_and_repulsive_at_the_core() {
    let h = 0.2;
    // Continuous across the branch change at h/2 and compactly supported.
    let below = SphKernels::cohesion(0.5 * h - 1e-9, h);
    let above = SphKernels::cohesion(0.5 * h + 1e-9, h);
    assert_close(below, above, 1e-5 * above, "continuity at h/2");
    assert_eq!(SphKernels::cohesion(h, h), 0.0);
    assert_eq!(SphKernels::cohesion(1.5 * h, h), 0.0);
    assert_eq!(SphKernels::cohesion(0.0, h), 0.0);
    assert_eq!(SphKernels::cohesion(-0.1, h), 0.0);
    // Repulsive core below ~0.27 h, attractive from there to the support edge.
    assert!(SphKernels::cohesion(0.1 * h, h) < 0.0);
    assert!(SphKernels::cohesion(0.2 * h, h) < 0.0);
    for q in [0.3, 0.4, 0.5, 0.6, 0.8, 0.95] {
        assert!(
            SphKernels::cohesion(q * h, h) > 0.0,
            "C({q} h) must attract"
        );
    }
}

// ── surface tension ──────────────────────────────────────────────────────────

/// Two particles of equal mass but different densities: closed form of the
/// cohesion plus curvature acceleration,
/// `a_i = γ K [ m C(r) + h |W'(r)| (m/ρ_i + m/ρ_j) ]` towards the partner with
/// `K = 2 ρ_0 / (ρ_i + ρ_j)`, equal and opposite for the two particles.
#[test]
fn surface_tension_of_a_pair_matches_the_closed_form() {
    let config = SphConfig {
        surface_tension: 0.5,
        ..quiet_config()
    };
    let h = config.h;
    let mass = 0.05;
    let (rho_i, rho_j) = (700.0, 1300.0);
    for q in [0.2, 0.3, 0.45, 0.6, 0.9] {
        let r = q * h;
        let mut particles = vec![
            particle([0.0; 3], [0.0; 3], mass, rho_i, 0.0),
            particle([r, 0.0, 0.0], [0.0; 3], mass, rho_j, 0.0),
        ];
        sph_compute_forces(&mut particles, &config);

        let k_ij = 2.0 * config.rest_density / (rho_i + rho_j);
        let magnitude = config.surface_tension
            * k_ij
            * (mass * cohesion_reference(r, h)
                + h * poly6_slope_magnitude(r, h) * (mass / rho_i + mass / rho_j));
        let a_i = particles[0].force[0] / rho_i;
        let a_j = particles[1].force[0] / rho_j;
        let tol = 1e-12 * magnitude.abs().max(1e-12);
        assert_close(a_i, magnitude, tol, &format!("a_i at {q} h"));
        assert_close(a_j, -magnitude, tol, &format!("a_j at {q} h"));
        // Equal masses: the pair exchanges equal and opposite momentum.
        assert_close(
            mass * (a_i + a_j),
            0.0,
            1e-12 * mass * magnitude.abs(),
            "momentum",
        );
        for p in &particles {
            assert_eq!(p.force[1], 0.0);
            assert_eq!(p.force[2], 0.0);
        }
    }
}

/// Beyond the cohesion sign change the pull is attractive: after one step the
/// pair moves towards each other; inside the repulsive core the cohesion
/// pushes them apart while the curvature term still pulls, so the sign of the
/// total follows the closed form (checked above) rather than a fixed rule.
#[test]
fn surface_tension_pulls_a_pair_in_the_attractive_range_together() {
    let config = SphConfig {
        surface_tension: 0.8,
        ..quiet_config()
    };
    let h = config.h;
    let mut particles = vec![
        particle([0.0; 3], [0.0; 3], 0.05, 1000.0, 0.0),
        particle([0.7 * h, 0.0, 0.0], [0.0; 3], 0.05, 1000.0, 0.0),
    ];
    sph_step_forces_only(&mut particles, &config);
    assert!(
        particles[0].velocity[0] > 0.0,
        "particle 0 must move towards +x"
    );
    assert!(
        particles[1].velocity[0] < 0.0,
        "particle 1 must move towards -x"
    );
    let gap = particles[1].position[0] - particles[0].position[0];
    assert!(gap < 0.7 * h, "gap {gap} should shrink");
}

/// Force and integration without recomputing density and pressure.
fn sph_step_forces_only(particles: &mut [SphParticle], config: &SphConfig) {
    sph_compute_forces(particles, config);
    sph_integrate(particles, config);
}

#[test]
fn surface_tension_is_off_for_zero_and_negative_coefficients() {
    let h = quiet_config().h;
    for gamma in [0.0, -0.5] {
        let config = SphConfig {
            surface_tension: gamma,
            ..quiet_config()
        };
        let mut particles = vec![
            particle([0.0; 3], [0.0; 3], 0.05, 1000.0, 0.0),
            particle([0.7 * h, 0.0, 0.0], [0.0; 3], 0.05, 1000.0, 0.0),
        ];
        sph_compute_forces(&mut particles, &config);
        for p in &particles {
            assert_eq!(p.force, [0.0; 3], "gamma = {gamma}");
        }
    }
}

/// Lattice of `3 x 3 x 3` particles centred on the origin with spacing
/// `spacing`, optionally jittered deterministically; densities from the SPH
/// summation.
fn cluster(spacing: f64, jitter: f64, config: &SphConfig) -> Vec<SphParticle> {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    let mut particles = Vec::new();
    for k in -1_i32..=1 {
        for j in -1_i32..=1 {
            for i in -1_i32..=1 {
                let position = [
                    f64::from(i) * spacing + jitter * next(),
                    f64::from(j) * spacing + jitter * next(),
                    f64::from(k) * spacing + jitter * next(),
                ];
                particles.push(SphParticle::new(position, 0.05));
            }
        }
    }
    sph_compute_density(&mut particles, config);
    particles
}

/// The normals point out of the cluster: parallel to the outward direction of
/// every non-central lattice particle (each coordinate has the sign of the
/// particle's own coordinate) and vanish at the centre by symmetry.
#[test]
fn surface_normals_point_out_of_a_cluster() {
    let config = SphConfig {
        surface_tension: 0.5,
        ..quiet_config()
    };
    let h = config.h;
    let particles = cluster(0.4 * h, 0.0, &config);
    let normals = sph_surface_normals(&particles, &config);
    assert_eq!(normals.len(), particles.len());
    let corner_size = normals.iter().map(|n| length3(*n)).fold(0.0_f64, f64::max);
    assert!(
        corner_size > 1e-3,
        "test setup: normals vanish ({corner_size})"
    );
    for (p, n) in particles.iter().zip(normals.iter()) {
        let is_centre = p.position.iter().all(|c| c.abs() < 1e-12);
        if is_centre {
            assert!(length3(*n) < 1e-12 * corner_size, "centre normal {n:?}");
            continue;
        }
        let outward = dot3(*n, p.position);
        assert!(
            outward > 0.0,
            "normal {n:?} of the particle at {:?} points inwards",
            p.position
        );
        for c in 0..3 {
            // Each normal component has the sign of the coordinate (0 on the
            // symmetry planes).
            if p.position[c].abs() < 1e-12 {
                assert!(n[c].abs() < 1e-12 * corner_size, "component {c} of {n:?}");
            } else {
                assert!(n[c] * p.position[c] > 0.0, "component {c} of {n:?}");
            }
        }
    }
}

/// Surface tension compresses a cluster: every non-central particle of the
/// symmetric lattice is accelerated towards the centre, and the central
/// particle, which is surrounded symmetrically, feels no net surface tension.
#[test]
fn surface_tension_pulls_the_particles_of_a_cluster_inwards() {
    let config = SphConfig {
        surface_tension: 0.5,
        ..quiet_config()
    };
    let h = config.h;
    let mut particles = cluster(0.4 * h, 0.0, &config);
    sph_compute_forces(&mut particles, &config);
    let accelerations: Vec<[f64; 3]> = particles
        .iter()
        .map(|p| scale3(p.force, 1.0 / p.density))
        .collect();
    let scale = accelerations
        .iter()
        .map(|a| length3(*a))
        .fold(0.0_f64, f64::max);
    assert!(scale > 1.0, "test setup: no surface tension acts ({scale})");
    for (p, accel) in particles.iter().zip(accelerations.iter()) {
        let is_centre = p.position.iter().all(|c| c.abs() < 1e-12);
        if is_centre {
            assert!(
                length3(*accel) < 1e-12 * scale,
                "centre acceleration {accel:?}"
            );
        } else {
            assert!(
                dot3(*accel, p.position) < 0.0,
                "particle at {:?} is pushed outwards: {accel:?}",
                p.position
            );
        }
    }
}

/// Pairwise antisymmetry (equal masses): the surface tension of a jittered
/// cluster with varying densities exchanges momentum between particles only,
/// so the total `Σ m a` vanishes.
#[test]
fn surface_tension_conserves_momentum_of_a_cluster() {
    let config = SphConfig {
        surface_tension: 0.7,
        ..quiet_config()
    };
    let h = config.h;
    let mut particles = cluster(0.35 * h, 0.2 * h, &config);
    let densities: Vec<f64> = particles.iter().map(|p| p.density).collect();
    assert!(
        densities
            .iter()
            .any(|d| (d - densities[0]).abs() > 1e-6 * densities[0]),
        "test setup: densities must differ"
    );
    sph_compute_forces(&mut particles, &config);
    let mut total = [0.0_f64; 3];
    let mut scale = 0.0_f64;
    for p in &particles {
        let accel = scale3(p.force, 1.0 / p.density);
        total = add3(total, scale3(accel, p.mass));
        scale = scale.max(p.mass * length3(accel));
    }
    assert!(
        scale > 1e-4,
        "test setup: no surface tension acts ({scale})"
    );
    assert!(
        length3(total) < 1e-12 * scale,
        "net momentum rate {total:?} vs scale {scale}"
    );
}

/// The surface tension term only adds to the other forces: with pressure,
/// viscosity and gravity present the difference to the `γ = 0` forces equals
/// the surface-tension-only forces.
#[test]
fn surface_tension_is_additive_to_the_other_force_densities() {
    let base = SphConfig {
        viscosity: 0.4,
        gravity: [0.0, -9.81, 0.0],
        ..quiet_config()
    };
    let with_tension = SphConfig {
        surface_tension: 0.6,
        ..base.clone()
    };
    let tension_only = SphConfig {
        viscosity: 0.0,
        gravity: [0.0; 3],
        ..with_tension.clone()
    };
    let h = base.h;
    let mut reference = cluster(0.4 * h, 0.1 * h, &base);
    sph_compute_pressure(&mut reference, &base);
    for (n, p) in reference.iter_mut().enumerate() {
        p.velocity = [0.01 * n as f64, -0.02 * n as f64, 0.005];
    }
    let mut a = reference.clone();
    let mut b = reference.clone();
    let mut c = reference;
    sph_compute_forces(&mut a, &base);
    sph_compute_forces(&mut b, &with_tension);
    // Pressure is part of `a` and `b`; subtract it by running the tension-only
    // configuration on zero-pressure particles.
    for p in c.iter_mut() {
        p.pressure = 0.0;
    }
    sph_compute_forces(&mut c, &tension_only);
    for ((pa, pb), pc) in a.iter().zip(b.iter()).zip(c.iter()) {
        for axis in 0..3 {
            let difference = pb.force[axis] - pa.force[axis];
            assert_close(
                difference,
                pc.force[axis],
                1e-9 * (1.0 + pc.force[axis].abs()),
                "tension contribution",
            );
        }
    }
}
