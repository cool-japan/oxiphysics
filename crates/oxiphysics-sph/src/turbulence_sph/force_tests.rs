//! Regression tests for the turbulent viscosity force and the SPH velocity
//! gradient: direction (viscous, not anti-viscous), magnitude (`ν_t ∇²v`),
//! momentum conservation and the use of the real particle mass.

use super::*;

const RHO: f64 = 1000.0;

fn particle(pos: [f64; 3], vel: [f64; 3], mass: f64, nut: f64) -> TurbulentParticle {
    TurbulentParticle {
        pos,
        vel,
        mass,
        density: RHO,
        pressure: 0.0,
        turbulent_ke: 0.0,
        turbulent_dissipation: 0.0,
        eddy_viscosity: nut,
        smoothing_length: 0.1,
    }
}

/// `n³` cubic lattice centred on the origin with spacing `dx` and mass `ρ dx³`.
/// Returns the particles and the index of the central one.
fn lattice(
    n: usize,
    dx: f64,
    nut: f64,
    velocity: impl Fn([f64; 3]) -> [f64; 3],
) -> (Vec<TurbulentParticle>, usize) {
    let half = (n as f64 - 1.0) / 2.0;
    let mut particles = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [
                    (ix as f64 - half) * dx,
                    (iy as f64 - half) * dx,
                    (iz as f64 - half) * dx,
                ];
                particles.push(particle(pos, velocity(pos), RHO * dx * dx * dx, nut));
            }
        }
    }
    let centre = (n / 2) * n * n + (n / 2) * n + n / 2;
    (particles, centre)
}

/// Deterministic pseudo-random numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_cloud(count: usize, seed: u64) -> Vec<TurbulentParticle> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|_| {
            let pos = [
                0.3 * rng.next_f64(),
                0.3 * rng.next_f64(),
                0.3 * rng.next_f64(),
            ];
            let vel = [
                rng.next_f64() - 0.5,
                rng.next_f64() - 0.5,
                rng.next_f64() - 0.5,
            ];
            let mass = 0.5 + rng.next_f64();
            let nut = 1.0e-3 * (0.5 + rng.next_f64());
            particle(pos, vel, mass, nut)
        })
        .collect()
}

fn brute_force_neighbors(particles: &[TurbulentParticle], h: f64) -> Vec<Vec<usize>> {
    (0..particles.len())
        .map(|i| {
            (0..particles.len())
                .filter(|&j| j != i && length(sub(particles[i].pos, particles[j].pos)) <= 2.0 * h)
                .collect()
        })
        .collect()
}

fn momentum(particles: &[TurbulentParticle]) -> [f64; 3] {
    particles
        .iter()
        .fold([0.0; 3], |acc, p| add(acc, scale(p.vel, p.mass)))
}

fn kinetic_energy(particles: &[TurbulentParticle]) -> f64 {
    particles
        .iter()
        .map(|p| 0.5 * p.mass * dot(p.vel, p.vel))
        .sum()
}

#[test]
fn pair_force_reduces_velocity_difference_and_conserves_momentum() {
    let h = 0.1;
    let dt = 0.01;
    let nut = 1.0e-3;
    let (m0, m1) = (0.2, 0.1);
    let mut particles = [
        particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], m0, nut),
        particle([0.05, 0.0, 0.0], [-1.0, 0.0, 0.0], m1, nut),
    ];
    let momentum_before = momentum(&particles);
    let diff_before = particles[0].vel[0] - particles[1].vel[0];
    TurbulentViscosityForce::apply(&mut particles, h, dt);
    let diff_after = particles[0].vel[0] - particles[1].vel[0];

    assert!(
        diff_after < diff_before,
        "viscosity must damp the velocity difference: {diff_before} -> {diff_after}"
    );
    assert!(
        diff_after > 0.0,
        "an explicit step of this size must not overshoot: {diff_after}"
    );
    assert!(particles[0].vel[0] < 1.0 && particles[1].vel[0] > -1.0);

    let momentum_after = momentum(&particles);
    for axis in 0..3 {
        assert!(
            (momentum_after[axis] - momentum_before[axis]).abs() < 1e-15,
            "momentum[{axis}] changed: {} -> {}",
            momentum_before[axis],
            momentum_after[axis]
        );
    }

    // Closed form of the Morris pair term with the real masses.
    let r = 0.05;
    let r_ij = [-r, 0.0, 0.0];
    let r_dot_grad_w = dot(r_ij, kernel_gradient(r_ij, h));
    assert!(r_dot_grad_w < 0.0, "r_ij·∇_iW_ij must be negative");
    let k = 2.0 * nut * r_dot_grad_w / (RHO * (r * r + 0.01 * h * h));
    let v0_expected = 1.0 + dt * m1 * k * 2.0;
    let v1_expected = -1.0 + dt * m0 * k * (-2.0);
    assert!(
        (particles[0].vel[0] - v0_expected).abs() < 1e-14,
        "v0 = {} but the Morris pair term gives {v0_expected}",
        particles[0].vel[0]
    );
    assert!(
        (particles[1].vel[0] - v1_expected).abs() < 1e-14,
        "v1 = {} but the Morris pair term gives {v1_expected}",
        particles[1].vel[0]
    );
    for p in &particles {
        assert_eq!(p.vel[1], 0.0);
        assert_eq!(p.vel[2], 0.0);
    }
}

#[test]
fn pair_force_is_antisymmetric_for_unequal_viscosities_and_masses() {
    let h = 0.1;
    let (m0, m1) = (0.7, 0.3);
    let p0 = particle([0.0, 0.0, 0.0], [0.4, -0.2, 0.1], m0, 2.0e-3);
    let p1 = particle([0.03, 0.04, 0.0], [-0.3, 0.5, 0.2], m1, 5.0e-4);
    let s0 = TurbParticleSnapshot::from(&p0);
    let s1 = TurbParticleSnapshot::from(&p1);
    let a0 = TurbulentViscosityForce::single_particle_force(
        p0.pos,
        p0.vel,
        RHO,
        p0.eddy_viscosity,
        &[s1],
        h,
    );
    let a1 = TurbulentViscosityForce::single_particle_force(
        p1.pos,
        p1.vel,
        RHO,
        p1.eddy_viscosity,
        &[s0],
        h,
    );
    for axis in 0..3 {
        let f0 = m0 * a0[axis];
        let f1 = m1 * a1[axis];
        assert!(
            (f0 + f1).abs() < 1e-14 * (f0.abs() + f1.abs()).max(1e-30),
            "axis {axis}: m0 a0 = {f0} must cancel m1 a1 = {f1}"
        );
    }
    // The acceleration of particle 0 points along (v1 - v0): it is dragged towards v1.
    let dv = sub(p1.vel, p0.vel);
    assert!(dot(a0, dv) > 0.0, "a0 = {a0:?}, v1 - v0 = {dv:?}");
}

#[test]
fn particles_beyond_the_kernel_support_do_not_interact() {
    let h = 0.1;
    let mut particles = [
        particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0e-2),
        particle([0.25, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.1, 1.0e-2),
    ];
    TurbulentViscosityForce::apply(&mut particles, h, 0.1);
    assert_eq!(particles[0].vel, [1.0, 0.0, 0.0]);
    assert_eq!(particles[1].vel, [-1.0, 0.0, 0.0]);
}

#[test]
fn nearly_coincident_particles_receive_a_bounded_regularised_force() {
    let h = 0.1;
    let mut particles = [
        particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, 1.0e-3),
        particle([1.0e-6, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.1, 1.0e-3),
    ];
    TurbulentViscosityForce::apply(&mut particles, h, 0.01);
    for p in &particles {
        assert!(p.vel[0].is_finite());
        assert!(
            p.vel[0].abs() <= 1.0,
            "no spurious acceleration: {}",
            p.vel[0]
        );
    }
    assert!(particles[0].vel[0] < 1.0 && particles[1].vel[0] > -1.0);
}

#[test]
fn quadratic_velocity_field_gives_two_nu_t_acceleration() {
    // u_x = y²  =>  ∇²u_x = 2  =>  a_x = ν_t ∇²u_x = 2 ν_t at the lattice centre.
    let dx = 0.05;
    let h = 1.25 * dx;
    let nut = 1.0e-3;
    let (particles, centre) = lattice(7, dx, nut, |x| [x[1] * x[1], 0.0, 0.0]);
    let snapshot: Vec<TurbParticleSnapshot> =
        particles.iter().map(TurbParticleSnapshot::from).collect();
    let c = &particles[centre];
    let acc =
        TurbulentViscosityForce::single_particle_force(c.pos, c.vel, c.density, nut, &snapshot, h);
    let expected = 2.0 * nut;
    assert!(
        (acc[0] - expected).abs() < 0.05 * expected,
        "a_x = {} but ν_t ∇²u_x = {expected}",
        acc[0]
    );
    assert!(
        acc[1].abs() < 1e-9 * expected && acc[2].abs() < 1e-9 * expected,
        "{acc:?}"
    );

    // `apply` produces the same increment: Δv = dt a.
    let dt = 1.0e-3;
    let mut stepped = particles;
    let vx_before = stepped[centre].vel[0];
    TurbulentViscosityForce::apply(&mut stepped, h, dt);
    let dv = stepped[centre].vel[0] - vx_before;
    assert!(
        (dv - dt * acc[0]).abs() < 1e-12 * dt * acc[0].abs(),
        "apply: Δv = {dv}, dt·a = {}",
        dt * acc[0]
    );
}

#[test]
fn uniform_and_rigidly_translating_flow_feels_no_force() {
    let dx = 0.05;
    let h = 1.25 * dx;
    let (mut particles, _) = lattice(5, dx, 1.0e-3, |_| [0.7, -0.2, 0.4]);
    TurbulentViscosityForce::apply(&mut particles, h, 0.01);
    for p in &particles {
        assert!((p.vel[0] - 0.7).abs() < 1e-15);
        assert!((p.vel[1] + 0.2).abs() < 1e-15);
        assert!((p.vel[2] - 0.4).abs() < 1e-15);
    }
}

#[test]
fn neighbour_list_version_matches_brute_force_version() {
    let h = 0.1;
    let dt = 5.0e-3;
    let mut brute = random_cloud(60, 7);
    let mut listed = random_cloud(60, 7);
    let neighbors = brute_force_neighbors(&listed, h);
    assert!(
        neighbors.iter().any(|list| list.len() > 3),
        "the cloud must be dense enough to produce interactions"
    );
    TurbulentViscosityForce::apply(&mut brute, h, dt);
    TurbulentViscosityForce::apply_with_neighbors(&mut listed, &neighbors, h, dt);
    for (a, b) in brute.iter().zip(listed.iter()) {
        for axis in 0..3 {
            assert!(
                (a.vel[axis] - b.vel[axis]).abs() < 1e-14,
                "apply {:?} vs apply_with_neighbors {:?}",
                a.vel,
                b.vel
            );
        }
    }
}

#[test]
fn neighbour_list_version_tolerates_malformed_lists() {
    let h = 0.1;
    let mut particles = random_cloud(5, 3);
    let reference: Vec<[f64; 3]> = particles.iter().map(|p| p.vel).collect();
    // Too few lists, self references and out-of-range indices must be ignored.
    let neighbors = [vec![0, 99], vec![1, 42, 77]];
    TurbulentViscosityForce::apply_with_neighbors(&mut particles, &neighbors, h, 1.0e-3);
    for (p, v) in particles.iter().zip(reference.iter()) {
        assert_eq!(&p.vel, v, "no valid neighbours: velocities must not change");
    }
}

#[test]
fn random_cloud_conserves_momentum_and_dissipates_kinetic_energy() {
    let h = 0.1;
    let dt = 1.0e-4;
    let mut particles = random_cloud(80, 11);
    let neighbors = brute_force_neighbors(&particles, h);
    let momentum_before = momentum(&particles);
    let scale_p: f64 = particles.iter().map(|p| p.mass * length(p.vel)).sum();
    let energy_before = kinetic_energy(&particles);
    TurbulentViscosityForce::apply_with_neighbors(&mut particles, &neighbors, h, dt);
    let momentum_after = momentum(&particles);
    for axis in 0..3 {
        assert!(
            (momentum_after[axis] - momentum_before[axis]).abs() < 1e-13 * scale_p,
            "momentum[{axis}]: {} -> {}",
            momentum_before[axis],
            momentum_after[axis]
        );
    }
    let energy_after = kinetic_energy(&particles);
    assert!(
        energy_after < energy_before,
        "turbulent viscosity must dissipate kinetic energy: {energy_before} -> {energy_after}"
    );
}

#[test]
fn velocity_gradient_recovers_a_linear_field_with_the_real_mass() {
    // v = A x  =>  ∇v = A.  With the former fabricated mass ρh³ (≈ 1.95 ρ dx³
    // at h = 1.25 dx) the result was about twice too large.
    let a = [[0.3, 0.2, -0.1], [0.0, -0.5, 0.4], [0.7, 0.1, 0.2]];
    let dx = 0.05;
    let h = 1.25 * dx;
    let (particles, centre) = lattice(7, dx, 1.0e-3, |x| {
        [
            a[0][0] * x[0] + a[0][1] * x[1] + a[0][2] * x[2],
            a[1][0] * x[0] + a[1][1] * x[1] + a[1][2] * x[2],
            a[2][0] * x[0] + a[2][1] * x[1] + a[2][2] * x[2],
        ]
    });
    let c = &particles[centre];
    let neighbors: Vec<(&TurbulentParticle, [f64; 3])> = particles
        .iter()
        .enumerate()
        .filter(|(j, p)| *j != centre && length(sub(c.pos, p.pos)) <= 2.0 * h)
        .map(|(_, p)| (p, kernel_gradient(sub(c.pos, p.pos), h)))
        .collect();
    let l = velocity_gradient_sph(c.pos, c.vel, &neighbors, h);
    // SPH consistency error on the lattice: below 2 % of the largest entry of A.
    let tolerance = 0.03 * 0.7;
    for row in 0..3 {
        for col in 0..3 {
            assert!(
                (l[row][col] - a[row][col]).abs() < tolerance,
                "∇v[{row}][{col}] = {} but A[{row}][{col}] = {}",
                l[row][col],
                a[row][col]
            );
        }
    }
}

#[test]
fn velocity_gradient_scales_linearly_with_particle_mass() {
    let h = 0.1;
    let make = |mass: f64| {
        let p = particle([0.05, 0.0, 0.0], [0.0, 1.0, 0.0], mass, 1.0e-3);
        let grad_w = kernel_gradient([-0.05, 0.0, 0.0], h);
        (p, grad_w)
    };
    let (light, g_light) = make(0.1);
    let (heavy, g_heavy) = make(0.3);
    let l_light = velocity_gradient_sph([0.0; 3], [0.0; 3], &[(&light, g_light)], h);
    let l_heavy = velocity_gradient_sph([0.0; 3], [0.0; 3], &[(&heavy, g_heavy)], h);
    // L_yx = (m/ρ) Δv_y ∂_xW with Δv_y = 1.
    let expected_light = 0.1 / RHO * g_light[0];
    assert!((l_light[1][0] - expected_light).abs() < 1e-15 * expected_light.abs());
    assert!((l_heavy[1][0] - 3.0 * l_light[1][0]).abs() < 1e-14 * l_heavy[1][0].abs());
}
