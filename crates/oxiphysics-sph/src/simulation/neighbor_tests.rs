//! Tests of the neighbour-list based pair sums of [`SphSim`] against all-pairs
//! reference implementations of the same formulas.

use super::sph_sim_neighbors::SphSimNeighbors;
use super::types::{SphSim, SphSimConfig};

type Vec3 = [f64; 3];

/// Deterministic uniform pseudo-random numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// A jittered lattice of `n` particles with random velocities and, optionally,
/// random masses.  The jitter keeps the density field non-uniform so that every
/// pair sum (pressure, viscosity, surface tension) is non-trivial.
fn jittered_sim(n: usize, seed: u64, random_masses: bool) -> SphSim {
    let mut sim = SphSim::new(SphSimConfig::default_water(), n);
    let mut rng = Lcg(seed);
    let spacing = (1.0 / (n as f64).cbrt().ceil()).min(1.0);
    for p in sim.positions.iter_mut() {
        for c in p.iter_mut() {
            *c += rng.range(-0.15, 0.15) * spacing;
        }
    }
    for v in sim.velocities.iter_mut() {
        *v = [
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        ];
    }
    if random_masses {
        for m in sim.masses.iter_mut() {
            *m *= rng.range(0.5, 1.5);
        }
    }
    sim
}

// ---------------------------------------------------------------------------
// All-pairs reference implementations (the pre-neighbour-list formulas, with
// the reaction force of a pair equal and opposite to the action)
// ---------------------------------------------------------------------------

fn reference_density(sim: &SphSim) -> Vec<f64> {
    let n = sim.len();
    let h = sim.config.kernel_radius;
    let h2 = h * h;
    (0..n)
        .map(|i| {
            let mut rho = sim.masses[i] * SphSim::cubic_w(0.0, h);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let d2 = SphSim::dist2(&sim.positions[i], &sim.positions[j]);
                if d2 > 4.0 * h2 {
                    continue;
                }
                rho += sim.masses[j] * SphSim::cubic_w(d2.sqrt(), h);
            }
            rho
        })
        .collect()
}

/// Tait pressure as computed by `SphSim::step`.
fn tait_pressures(sim: &SphSim) -> Vec<f64> {
    let rho0 = sim.config.rest_density;
    let c0 = sim.config.estimated_sound_speed().max(1.0);
    let gamma = 7.0_f64;
    let b = rho0 * c0 * c0 / gamma;
    sim.densities
        .iter()
        .map(|rho| b * ((rho / rho0).powf(gamma) - 1.0))
        .collect()
}

fn reference_pressure_forces(sim: &SphSim) -> Vec<Vec3> {
    let n = sim.len();
    let h = sim.config.kernel_radius;
    let h2 = h * h;
    let (x, rho, pressure, mass) = (&sim.positions, &sim.densities, &sim.pressures, &sim.masses);
    let mut forces = vec![[0.0; 3]; n];
    for i in 0..n {
        let rhoi = rho[i].max(1e-14);
        let pi = pressure[i];
        for j in (i + 1)..n {
            let d2 = SphSim::dist2(&x[i], &x[j]);
            if d2 > 4.0 * h2 {
                continue;
            }
            let r = d2.sqrt();
            if r < 1e-14 {
                continue;
            }
            let rhoj = rho[j].max(1e-14);
            let pj = pressure[j];
            let grad = SphSim::cubic_grad_w(r, h);
            let factor = -mass[j] * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * grad;
            for k in 0..3 {
                let rhat_k = (x[i][k] - x[j][k]) / r;
                let fk = factor * rhat_k;
                forces[i][k] += fk * mass[i];
                forces[j][k] -= fk * mass[i];
            }
        }
    }
    forces
}

fn reference_viscosity_forces(sim: &SphSim) -> Vec<Vec3> {
    let n = sim.len();
    let h = sim.config.kernel_radius;
    let h2 = h * h;
    let mu = sim.config.viscosity;
    let (x, v, rho, mass) = (&sim.positions, &sim.velocities, &sim.densities, &sim.masses);
    let mut forces = vec![[0.0; 3]; n];
    for i in 0..n {
        let rhoi = rho[i].max(1e-14);
        for j in (i + 1)..n {
            let d2 = SphSim::dist2(&x[i], &x[j]);
            if d2 > 4.0 * h2 || d2 < 1e-28 {
                continue;
            }
            let r = d2.sqrt();
            let rhoj = rho[j].max(1e-14);
            let rho_avg = 0.5 * (rhoi + rhoj);
            let mut vr_dot = 0.0_f64;
            for k in 0..3 {
                vr_dot += (v[i][k] - v[j][k]) * (x[i][k] - x[j][k]);
            }
            if vr_dot >= 0.0 {
                continue;
            }
            let nu = h * vr_dot / (d2 + 0.01 * h2);
            let pi_visc = -mu * nu / rho_avg;
            let grad = SphSim::cubic_grad_w(r, h);
            for k in 0..3 {
                let rhat_k = (x[i][k] - x[j][k]) / r;
                let fk = -mass[j] * pi_visc * grad * rhat_k;
                forces[i][k] += fk * mass[i];
                forces[j][k] -= fk * mass[i];
            }
        }
    }
    forces
}

fn reference_surface_tension_forces(sim: &SphSim) -> Vec<Vec3> {
    let n = sim.len();
    let h = sim.config.kernel_radius;
    let h2 = h * h;
    let sigma = sim.config.surface_tension_coeff;
    let (x, rho, mass) = (&sim.positions, &sim.densities, &sim.masses);
    let mut forces = vec![[0.0; 3]; n];
    let mut grad_c = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        let rhoi = rho[i].max(1e-14);
        for j in 0..n {
            if j == i {
                continue;
            }
            let d2 = SphSim::dist2(&x[i], &x[j]);
            if d2 > 4.0 * h2 || d2 < 1e-28 {
                continue;
            }
            let r = d2.sqrt();
            let rhoj = rho[j].max(1e-14);
            let ci = rhoi / sim.config.rest_density;
            let cj = rhoj / sim.config.rest_density;
            let grad_w = SphSim::cubic_grad_w(r, h);
            let factor = mass[j] / rhoj * (cj - ci) * grad_w / r;
            for k in 0..3 {
                grad_c[i][k] += factor * (x[j][k] - x[i][k]);
            }
        }
    }
    for i in 0..n {
        let gc_mag2 = grad_c[i].iter().map(|c| c * c).sum::<f64>();
        if gc_mag2 < 1e-20 {
            continue;
        }
        let gc_mag = gc_mag2.sqrt();
        let n_hat = [
            grad_c[i][0] / gc_mag,
            grad_c[i][1] / gc_mag,
            grad_c[i][2] / gc_mag,
        ];
        let rhoi = rho[i].max(1e-14);
        let mut kappa = 0.0_f64;
        for j in 0..n {
            if j == i {
                continue;
            }
            let d2 = SphSim::dist2(&x[i], &x[j]);
            if d2 > 4.0 * h2 || d2 < 1e-28 {
                continue;
            }
            let r = d2.sqrt();
            let rhoj = rho[j].max(1e-14);
            let gc_j_mag2 = grad_c[j].iter().map(|c| c * c).sum::<f64>();
            if gc_j_mag2 < 1e-20 {
                continue;
            }
            let gc_j_mag = gc_j_mag2.sqrt();
            let n_hat_j = [
                grad_c[j][0] / gc_j_mag,
                grad_c[j][1] / gc_j_mag,
                grad_c[j][2] / gc_j_mag,
            ];
            let grad_w = SphSim::cubic_grad_w(r, h);
            let mut dn_dot_r = 0.0_f64;
            for k in 0..3 {
                dn_dot_r += (n_hat_j[k] - n_hat[k]) * (x[j][k] - x[i][k]);
            }
            kappa -= mass[j] / rhoj * dn_dot_r * grad_w / r;
        }
        let force_scale = sigma * kappa * gc_mag * mass[i] / rhoi;
        for (k, nk) in n_hat.iter().enumerate() {
            forces[i][k] += force_scale * nk;
        }
    }
    forces
}

/// All-pairs reference of `SphSim::step` (density → Tait → gravity, pressure and
/// viscosity forces → Euler).
fn reference_step(sim: &mut SphSim) {
    sim.densities = reference_density(sim);
    sim.pressures = tait_pressures(sim);
    let pressure = reference_pressure_forces(sim);
    let viscosity = reference_viscosity_forces(sim);
    let dt = sim.config.dt;
    let g = sim.config.gravity;
    for i in 0..sim.len() {
        let m = sim.masses[i];
        for (k, gk) in g.iter().enumerate() {
            let force = m * gk + pressure[i][k] + viscosity[i][k];
            sim.velocities[i][k] += force / m.max(1e-30) * dt;
            sim.positions[i][k] += sim.velocities[i][k] * dt;
        }
    }
    sim.state.advance(dt);
}

fn max_abs(values: &[Vec3]) -> f64 {
    values
        .iter()
        .flat_map(|v| v.iter())
        .fold(0.0_f64, |m, x| m.max(x.abs()))
}

fn assert_vec_fields_close(actual: &[Vec3], expected: &[Vec3], rel: f64, what: &str) {
    assert_eq!(actual.len(), expected.len());
    let scale = max_abs(expected).max(f64::MIN_POSITIVE);
    for (i, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        for k in 0..3 {
            assert!(
                (a[k] - e[k]).abs() <= rel * scale,
                "{what}: particle {i} axis {k}: {} vs {} (scale {scale})",
                a[k],
                e[k]
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Equality with the all-pairs sums
// ---------------------------------------------------------------------------

#[test]
fn density_sum_equals_the_all_pairs_sum_on_1000_particles() {
    let mut sim = jittered_sim(1000, 1, true);
    assert_eq!(sim.len(), 1000);
    sim.compute_density_sum();
    let expected = reference_density(&sim);
    let scale = expected.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    assert!(scale > 0.0);
    for (i, (&rho, &reference)) in sim.densities.iter().zip(expected.iter()).enumerate() {
        assert!(
            (rho - reference).abs() <= 1e-12 * scale,
            "particle {i}: {rho} vs {reference}"
        );
    }
}

#[test]
fn pressure_forces_equal_the_all_pairs_sum_on_1000_particles() {
    let mut sim = jittered_sim(1000, 2, true);
    sim.compute_density_sum();
    sim.pressures = tait_pressures(&sim);
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_pressure_forces();
    let expected = reference_pressure_forces(&sim);
    assert!(max_abs(&expected) > 0.0);
    assert_vec_fields_close(&sim.forces, &expected, 1e-12, "pressure force");
}

#[test]
fn viscosity_forces_equal_the_all_pairs_sum_on_1000_particles() {
    let mut sim = jittered_sim(1000, 3, true);
    sim.compute_density_sum();
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_viscosity_forces();
    let expected = reference_viscosity_forces(&sim);
    assert!(max_abs(&expected) > 0.0);
    assert_vec_fields_close(&sim.forces, &expected, 1e-12, "viscosity force");
}

#[test]
fn surface_tension_forces_equal_the_all_pairs_sum_on_1000_particles() {
    let mut sim = jittered_sim(1000, 4, true);
    sim.compute_density_sum();
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_surface_tension_forces();
    let expected = reference_surface_tension_forces(&sim);
    assert!(max_abs(&expected) > 0.0);
    assert_vec_fields_close(&sim.forces, &expected, 1e-12, "surface tension force");
}

#[test]
fn surface_tension_is_a_no_op_for_a_vanishing_coefficient() {
    let mut sim = jittered_sim(125, 5, false);
    sim.config.surface_tension_coeff = 0.0;
    sim.compute_density_sum();
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_surface_tension_forces();
    assert_eq!(max_abs(&sim.forces), 0.0);
}

#[test]
fn full_steps_follow_the_all_pairs_reference() {
    let mut fast = jittered_sim(343, 6, true);
    let mut slow = fast.clone();
    for step in 0..3 {
        fast.step();
        reference_step(&mut slow);
        assert_eq!(fast.state.step, slow.state.step);
        assert_vec_fields_close(&fast.positions, &slow.positions, 1e-12, "position");
        assert_vec_fields_close(&fast.velocities, &slow.velocities, 1e-12, "velocity");
        let scale = slow.densities.iter().fold(0.0_f64, |m, r| m.max(*r));
        for (i, (a, b)) in fast.densities.iter().zip(slow.densities.iter()).enumerate() {
            assert!((a - b).abs() <= 1e-12 * scale, "step {step} particle {i}");
        }
    }
}

// ---------------------------------------------------------------------------
// Newton's third law (the reaction used mass j instead of mass i)
// ---------------------------------------------------------------------------

fn total_force(forces: &[Vec3]) -> (Vec3, f64) {
    let mut total = [0.0; 3];
    let mut scale = 0.0;
    for f in forces {
        for k in 0..3 {
            total[k] += f[k];
            scale += f[k].abs();
        }
    }
    (total, scale)
}

#[test]
fn a_pressure_pair_exerts_equal_and_opposite_forces_for_unequal_masses() {
    let mut sim = SphSim::new(SphSimConfig::default_water(), 2);
    sim.positions = vec![[0.0; 3], [0.1, 0.0, 0.0]];
    sim.masses = vec![1.0, 4.0];
    sim.densities = vec![1000.0, 1000.0];
    sim.pressures = vec![5000.0, 1000.0];
    sim.forces = vec![[0.0; 3]; 2];
    sim.compute_pressure_forces();
    // F_0 = -m_0 m_1 (P_0/rho_0^2 + P_1/rho_1^2) W'(r) (x_0 - x_1)/r.
    let h = sim.config.kernel_radius;
    let w_prime = SphSim::cubic_grad_w(0.1, h);
    let magnitude = 4.0 * (5000.0 / 1.0e6 + 1000.0 / 1.0e6) * w_prime;
    // (x_0 - x_1)/r = -x and W' < 0, so F_0 points to -x: away from the neighbour.
    assert!((sim.forces[0][0] - magnitude).abs() <= 1e-12 * magnitude.abs());
    assert!(sim.forces[0][0] < 0.0);
    assert_eq!(sim.forces[1][0], -sim.forces[0][0]);
    for k in 1..3 {
        assert_eq!(sim.forces[0][k], 0.0);
        assert_eq!(sim.forces[1][k], 0.0);
    }
}

#[test]
fn a_viscous_pair_exerts_equal_and_opposite_forces_for_unequal_masses() {
    let mut sim = SphSim::new(SphSimConfig::default_water(), 2);
    sim.positions = vec![[0.0; 3], [0.1, 0.0, 0.0]];
    sim.velocities = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    sim.masses = vec![1.0, 4.0];
    sim.densities = vec![1000.0, 1000.0];
    sim.forces = vec![[0.0; 3]; 2];
    sim.compute_viscosity_forces();
    assert!(
        sim.forces[0][0].abs() > 0.0,
        "approaching pair must interact"
    );
    // Approaching along x: the viscous force decelerates the approach, i.e.
    // it pushes particle 0 to -x.
    assert!(sim.forces[0][0] < 0.0);
    assert_eq!(sim.forces[1][0], -sim.forces[0][0]);
}

#[test]
fn pressure_and_viscosity_forces_conserve_momentum_for_random_masses() {
    let mut sim = jittered_sim(512, 7, true);
    sim.compute_density_sum();
    sim.pressures = tait_pressures(&sim);
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_pressure_forces();
    let (total, scale) = total_force(&sim.forces);
    assert!(scale > 0.0);
    for component in total {
        assert!(
            component.abs() <= 1e-10 * scale,
            "net pressure force {total:?} (scale {scale})"
        );
    }
    sim.forces = vec![[0.0; 3]; sim.len()];
    sim.compute_viscosity_forces();
    let (total, scale) = total_force(&sim.forces);
    assert!(scale > 0.0);
    for component in total {
        assert!(
            component.abs() <= 1e-10 * scale,
            "net viscous force {total:?} (scale {scale})"
        );
    }
}

// ---------------------------------------------------------------------------
// Cached neighbour lists
// ---------------------------------------------------------------------------

#[test]
fn neighbour_lists_follow_changes_of_the_public_position_field() {
    let mut sim = jittered_sim(216, 8, false);
    sim.compute_density_sum();
    let before = sim.densities.clone();
    // Move one particle far away through the public field.
    sim.positions[10] = [5.0, 5.0, 5.0];
    sim.compute_density_sum();
    let expected = reference_density(&sim);
    let scale = expected.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    for (i, (&rho, &reference)) in sim.densities.iter().zip(expected.iter()).enumerate() {
        assert!((rho - reference).abs() <= 1e-12 * scale, "particle {i}");
    }
    assert!(
        sim.densities
            .iter()
            .zip(before.iter())
            .any(|(a, b)| (a - b).abs() > 1e-9 * scale),
        "the move must change some density"
    );
}

#[test]
fn neighbour_lists_follow_changes_of_the_kernel_radius() {
    let mut sim = jittered_sim(216, 9, false);
    sim.compute_density_sum();
    let entries_small = sim.refresh_neighbors().entry_count();
    sim.config.kernel_radius *= 1.5;
    sim.compute_density_sum();
    let entries_large = sim.refresh_neighbors().entry_count();
    assert!(
        entries_large > entries_small,
        "{entries_large} <= {entries_small}"
    );
    let expected = reference_density(&sim);
    let scale = expected.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    for (i, (&rho, &reference)) in sim.densities.iter().zip(expected.iter()).enumerate() {
        assert!((rho - reference).abs() <= 1e-12 * scale, "particle {i}");
    }
}

#[test]
fn refresh_neighbors_reuses_lists_until_something_changes() {
    let mut sim = jittered_sim(125, 10, false);
    let h = sim.config.kernel_radius;
    assert!(!sim.refresh_neighbors().is_empty());
    let snapshot = sim.positions.clone();
    assert!(sim.refresh_neighbors().is_current(&snapshot, h));
    assert!(!sim.refresh_neighbors().is_current(&snapshot, 2.0 * h));
    let mut moved = snapshot.clone();
    moved[0][0] += 1e-9;
    assert!(!sim.refresh_neighbors().is_current(&moved, h));
    sim.step();
    // The step moved the particles, so the lists were rebuilt for the new positions.
    let positions = sim.positions.clone();
    assert!(sim.refresh_neighbors().is_current(&positions, h));
}

#[test]
fn neighbour_lists_are_symmetric_sorted_and_free_of_self_entries() {
    let sim = jittered_sim(343, 11, false);
    let h = sim.config.kernel_radius;
    let lists = SphSimNeighbors::build(&sim.positions, h);
    assert_eq!(lists.len(), sim.len());
    let radius_sq = 4.0 * h * h;
    for i in 0..sim.len() {
        let row = lists.neighbors_of(i);
        assert!(
            row.windows(2).all(|w| w[0] < w[1]),
            "row {i} must be strictly ascending"
        );
        assert!(!row.contains(&i));
        for &j in row {
            assert!(
                lists.neighbors_of(j).contains(&i),
                "pair ({i}, {j}) not symmetric"
            );
        }
        let expected: Vec<usize> = (0..sim.len())
            .filter(|&j| j != i && SphSim::dist2(&sim.positions[i], &sim.positions[j]) <= radius_sq)
            .collect();
        assert_eq!(row, expected.as_slice(), "row {i}");
    }
    assert!(lists.neighbors_of(sim.len()).is_empty());
    assert!(lists.neighbors_of(usize::MAX).is_empty());
}

#[test]
fn neighbour_lists_scale_linearly_with_the_particle_count() {
    // Same spacing and kernel radius, 8x as many particles: the number of
    // neighbours per particle stays bounded instead of growing with N.
    let h = 0.1;
    let per_particle = |n_side: usize| {
        let spacing = 0.1;
        let mut positions = Vec::new();
        for i in 0..n_side {
            for j in 0..n_side {
                for k in 0..n_side {
                    positions.push([
                        (i as f64 + 0.5) * spacing,
                        (j as f64 + 0.5) * spacing,
                        (k as f64 + 0.5) * spacing,
                    ]);
                }
            }
        }
        let lists = SphSimNeighbors::build(&positions, h);
        lists.entry_count() as f64 / positions.len() as f64
    };
    let small = per_particle(10);
    let large = per_particle(20);
    assert!(small > 10.0 && large < 40.0, "{small}, {large}");
    assert!(
        large / small < 2.0,
        "neighbours per particle grew from {small} to {large}"
    );
}

#[test]
fn degenerate_inputs_give_empty_lists_and_do_not_panic() {
    let mut empty = SphSim::new(SphSimConfig::default_water(), 0);
    empty.step();
    empty.step_with_surface_tension();
    assert!(empty.is_empty());
    assert_eq!(empty.state.step, 2);

    for radius in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let lists = SphSimNeighbors::build(&[[0.0; 3], [0.01, 0.0, 0.0]], radius);
        assert_eq!(lists.len(), 2);
        assert!(lists.neighbors_of(0).is_empty());
        assert!(lists.neighbors_of(1).is_empty());
    }
    assert!(SphSimNeighbors::empty().is_empty());
    assert!(!SphSimNeighbors::empty().is_current(&[], 0.1));
}

#[test]
fn the_floor_reflection_is_seen_by_the_next_pair_sum() {
    let mut sim = jittered_sim(125, 12, false);
    sim.compute_density_sum();
    for p in sim.positions.iter_mut() {
        p[1] -= 0.4;
    }
    sim.apply_floor_reflection(0.0, 0.5);
    sim.compute_density_sum();
    let expected = reference_density(&sim);
    let scale = expected.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    for (i, (&rho, &reference)) in sim.densities.iter().zip(expected.iter()).enumerate() {
        assert!((rho - reference).abs() <= 1e-12 * scale, "particle {i}");
    }
}
