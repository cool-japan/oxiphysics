//! Tests of the Akinci et al. (2013) pairwise surface tension model: cohesion plus
//! curvature minimisation with the normal field `n_i` and the density correction
//! `K_ij = 2 ρ₀/(ρ_i + ρ_j)`.

use super::types::PairwiseSurfaceTension;
use crate::kernel::{CubicSplineKernel, SphKernel};
use crate::neighbor::SpatialHash;
use oxiphysics_core::math::Vec3;

/// Particle spacing (m).
const DX: f64 = 0.1;
/// Droplet radius (m).
const RADIUS: f64 = 0.5;
/// Support radius of the cohesion and normal kernels: `h = 2 dx` (Akinci: `h = 4 r`).
const SUPPORT: f64 = 2.0 * DX;
/// Nominal water density (kg/m³).
const RHO_WATER: f64 = 1000.0;

/// Spherical droplet cut out of a cubic lattice.
struct Droplet {
    positions: Vec<Vec3>,
    masses: Vec<f64>,
    densities: Vec<f64>,
    neighbors: Vec<Vec<usize>>,
    /// Index of the particle at the centre of the sphere.
    centre: usize,
}

fn droplet() -> Droplet {
    let n = (RADIUS / DX).ceil() as i32 + 1;
    let mut positions = Vec::new();
    for ix in -n..=n {
        for iy in -n..=n {
            for iz in -n..=n {
                let p = Vec3::new(ix as f64 * DX, iy as f64 * DX, iz as f64 * DX);
                if p.norm() <= RADIUS + 1.0e-9 {
                    positions.push(p);
                }
            }
        }
    }
    let count = positions.len();
    let mass = RHO_WATER * DX * DX * DX;
    let neighbors = SpatialHash::find_all_neighbors(&positions, SUPPORT);
    // SPH summation density with the cubic spline of support `SUPPORT` (smoothing DX).
    let kernel = CubicSplineKernel;
    let densities: Vec<f64> = (0..count)
        .map(|i| {
            let mut rho = mass * kernel.w(0.0, 0.5 * SUPPORT);
            for &j in &neighbors[i] {
                rho += mass * kernel.w((positions[i] - positions[j]).norm(), 0.5 * SUPPORT);
            }
            rho
        })
        .collect();
    let centre = (0..count)
        .min_by(|&a, &b| positions[a].norm().total_cmp(&positions[b].norm()))
        .unwrap_or(0);
    Droplet {
        positions,
        masses: vec![mass; count],
        densities,
        neighbors,
        centre,
    }
}

impl Droplet {
    /// Model whose rest density is the bulk density measured at the droplet centre.
    fn model(&self, gamma: f64) -> PairwiseSurfaceTension {
        PairwiseSurfaceTension::new(gamma, SUPPORT).with_rest_density(self.densities[self.centre])
    }

    /// Indices of the outermost lattice layer.
    fn surface(&self) -> Vec<usize> {
        (0..self.positions.len())
            .filter(|&i| self.positions[i].norm() > RADIUS - 0.6 * DX)
            .collect()
    }

    /// Radial component of `vectors[i]` for particle `i`.
    fn radial(&self, vectors: &[Vec3], i: usize) -> f64 {
        vectors[i].dot(&(self.positions[i] / self.positions[i].norm()))
    }
}

#[test]
fn droplet_fixture_is_a_sphere_with_bulk_density_near_rest_density() {
    let d = droplet();
    assert!(d.positions.len() > 400, "{} particles", d.positions.len());
    assert!(d.surface().len() > 100);
    assert!(d.positions[d.centre].norm() < 1.0e-12);
    assert!((d.densities[d.centre] - RHO_WATER).abs() < 0.01 * RHO_WATER);
    // Free-surface particles are under-dense (missing neighbours).
    for &i in d.surface().iter() {
        assert!(d.densities[i] < 0.9 * d.densities[d.centre]);
    }
}

#[test]
fn normals_point_out_of_the_fluid_and_vanish_in_the_bulk() {
    let d = droplet();
    let model = d.model(1.0);
    let normals = model.compute_normals(&d.positions, &d.masses, &d.densities, &d.neighbors);
    for &i in d.surface().iter() {
        let radial = d.radial(&normals, i);
        assert!(
            radial > 0.8 && radial < 1.8,
            "surface normal must be outward and O(1): {radial}"
        );
        // The normal of a spherical surface is radial.
        let tangential = (normals[i] - radial * d.positions[i] / d.positions[i].norm()).norm();
        assert!(tangential < 0.7 * radial, "tangential part {tangential}");
    }
    for (position, normal) in d.positions.iter().zip(normals.iter()) {
        if position.norm() <= 2.0 * DX + 1.0e-9 {
            assert!(
                normal.norm() < 1.0e-3,
                "bulk normal must vanish: {}",
                normal.norm()
            );
        }
    }
}

#[test]
fn surface_particles_of_a_droplet_are_pulled_inward() {
    let d = droplet();
    let model = d.model(1.0);
    let forces = model.compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    let surface = d.surface();
    let mut mean = 0.0;
    for &i in surface.iter() {
        let radial = d.radial(&forces, i);
        assert!(
            radial < 0.0,
            "surface particle {i} is pushed outward: {radial}"
        );
        mean += radial / surface.len() as f64;
    }
    assert!(mean < -50.0, "mean inward force {mean}");
}

#[test]
fn curvature_minimisation_alone_contracts_the_droplet() {
    // The curvature term must act like surface tension, not against it: an inward
    // pointing normal field would turn it into a force that expands the droplet.
    let d = droplet();
    let model = d.model(1.0);
    let curvature = model.pair_forces(
        &d.positions,
        &d.masses,
        &d.densities,
        &d.neighbors,
        false,
        true,
    );
    let cohesion = model.pair_forces(
        &d.positions,
        &d.masses,
        &d.densities,
        &d.neighbors,
        true,
        false,
    );
    let total = model.compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    for &i in d.surface().iter() {
        assert!(
            d.radial(&curvature, i) < 0.0,
            "curvature force on surface particle {i} must point inward: {}",
            d.radial(&curvature, i)
        );
        assert!(d.radial(&cohesion, i) < 0.0);
        let sum = curvature[i] + cohesion[i];
        assert!((sum - total[i]).norm() <= 1.0e-12 * total[i].norm().max(1.0e-30));
    }
}

#[test]
fn forces_sum_to_zero_for_equal_masses() {
    let d = droplet();
    let model = d.model(0.72);
    let forces = model.compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    let net = forces.iter().fold(Vec3::zeros(), |acc, f| acc + f);
    let scale = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
    assert!(scale > 10.0, "the forces are not trivially zero: {scale}");
    assert!(
        net.norm() < 1.0e-9 * scale,
        "net force {} vs {scale}",
        net.norm()
    );
}

#[test]
fn centre_particle_feels_no_force_by_symmetry() {
    let d = droplet();
    let model = d.model(1.0);
    let forces = model.compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    let scale = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
    assert!(forces[d.centre].norm() < 1.0e-9 * scale);
}

#[test]
fn forces_are_linear_in_gamma_and_in_rest_density() {
    let d = droplet();
    let base = d.model(1.0);
    let f1 = base.compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    let f3 = d
        .model(3.0)
        .compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    // K_ij = 2 ρ₀/(ρ_i+ρ_j) is proportional to ρ₀: doubling ρ₀ doubles every force.
    let heavy = base
        .clone()
        .with_rest_density(2.0 * base.rest_density)
        .compute_forces(&d.positions, &d.masses, &d.densities, &d.neighbors);
    let mut checked = 0;
    for i in 0..f1.len() {
        if f1[i].norm() > 1.0e-6 {
            assert!((f3[i] - 3.0 * f1[i]).norm() < 1.0e-9 * f1[i].norm());
            assert!((heavy[i] - 2.0 * f1[i]).norm() < 1.0e-9 * f1[i].norm());
            checked += 1;
        }
    }
    assert!(checked > 100);
}

#[test]
fn pair_force_matches_the_akinci_formula() {
    // Two particles on the x axis, r = 0.6 h apart, with different densities.
    let h = 0.2;
    let r = 0.6 * h;
    let (m, rho_a, rho_b, rho0, gamma) = (0.4, 900.0, 1100.0, 1200.0, 0.7);
    let model = PairwiseSurfaceTension::new(gamma, h).with_rest_density(rho0);
    let positions = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(r, 0.0, 0.0)];
    let neighbors = vec![vec![1_usize], vec![0_usize]];
    let forces = model.compute_forces(&positions, &[m, m], &[rho_a, rho_b], &neighbors);

    let kernel = CubicSplineKernel;
    let dw_dr = kernel.grad_w(r, 0.5 * h);
    assert!(dw_dr < 0.0);
    // Outward normals: a points to −x (away from b), b to +x (away from a).
    let n_a = Vec3::new(h * (m / rho_b) * dw_dr, 0.0, 0.0);
    let n_b = Vec3::new(-h * (m / rho_a) * dw_dr, 0.0, 0.0);
    assert!(n_a.x < 0.0 && n_b.x > 0.0);
    let k_ij = 2.0 * rho0 / (rho_a + rho_b);
    assert!((k_ij - 1.2).abs() < 1.0e-12);
    let c = model.cohesion_kernel(r);
    assert!(c > 0.0, "r = 0.6 h lies in the attractive range of C");
    // F_a←b = K [−γ m² C x̂_ab − γ m (n_a − n_b)],  x̂_ab = (x_a − x_b)/r = −x̂.
    let expected_a = k_ij * (gamma * m * m * c * Vec3::x() - gamma * m * (n_a - n_b));
    let expected_b = k_ij * (-gamma * m * m * c * Vec3::x() - gamma * m * (n_b - n_a));
    assert!(
        (forces[0] - expected_a).norm() < 1.0e-12 * expected_a.norm(),
        "{:?} vs {:?}",
        forces[0],
        expected_a
    );
    assert!((forces[1] - expected_b).norm() < 1.0e-12 * expected_b.norm());
    // Attraction towards the partner, equal and opposite for equal masses.
    assert!(forces[0].x > 0.0 && forces[1].x < 0.0);
    assert!((forces[0] + forces[1]).norm() < 1.0e-12 * forces[0].norm());
}

#[test]
fn pair_beyond_the_smoothing_length_does_not_interact() {
    let h = 0.2;
    let model = PairwiseSurfaceTension::new(1.0, h);
    let positions = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5 * h, 0.0, 0.0)];
    let neighbors = vec![vec![1_usize], vec![0_usize]];
    let forces = model.compute_forces(&positions, &[0.4, 0.4], &[1000.0, 1000.0], &neighbors);
    assert_eq!(forces[0], Vec3::zeros());
    assert_eq!(forces[1], Vec3::zeros());
}

#[test]
fn isolated_particle_and_empty_input_are_handled() {
    let model = PairwiseSurfaceTension::new(1.0, 0.2);
    let single = model.compute_forces(&[Vec3::zeros()], &[1.0], &[1000.0], &[Vec::new()]);
    assert_eq!(single, vec![Vec3::zeros()]);
    assert!(model.compute_forces(&[], &[], &[], &[]).is_empty());
}

#[test]
fn default_rest_density_is_water_and_the_builder_overrides_it() {
    let model = PairwiseSurfaceTension::new(0.072, 0.1);
    assert_eq!(model.rest_density, 1000.0);
    assert_eq!(model.with_rest_density(850.0).rest_density, 850.0);
}

#[test]
#[should_panic(expected = "rest density")]
fn rest_density_must_be_positive() {
    let _ = PairwiseSurfaceTension::new(1.0, 0.1).with_rest_density(0.0);
}

#[test]
#[should_panic(expected = "rest density")]
fn rest_density_must_be_finite() {
    let _ = PairwiseSurfaceTension::new(1.0, 0.1).with_rest_density(f64::NAN);
}
