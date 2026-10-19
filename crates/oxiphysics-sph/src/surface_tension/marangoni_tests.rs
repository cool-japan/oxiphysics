//! Tests of the Marangoni force: a tangential interface stress `(dσ/dT) ∇_s T δ_s`,
//! not a bulk force along `∇T`.

use super::functions::{compute_multiphase_color_gradient, marangoni_forces};
use crate::kernel::CubicSplineKernel;
use crate::neighbor::SpatialHash;
use oxiphysics_core::math::Vec3;

const DX: f64 = 0.1;
const H: f64 = 1.3 * DX;
const NX: usize = 9;
const NY: usize = 8;
const NZ: usize = 5;
/// Water-like temperature coefficient of the surface tension (N/(m K)).
const DSIGMA_DT: f64 = -1.5e-4;

/// Two stacked fluid phases (`y < NY/2` is phase 0) in a regular box; the flat
/// interface between them has the normal `±ŷ`.
struct Layers {
    positions: Vec<Vec3>,
    masses: Vec<f64>,
    densities: Vec<f64>,
    phase: Vec<u32>,
    neighbors: Vec<Vec<usize>>,
    color_grad: Vec<Vec3>,
}

fn layers() -> Layers {
    let mut positions = Vec::new();
    let mut phase = Vec::new();
    for ix in 0..NX {
        for iy in 0..NY {
            for iz in 0..NZ {
                positions.push(Vec3::new(ix as f64 * DX, iy as f64 * DX, iz as f64 * DX));
                phase.push(if iy < NY / 2 { 0 } else { 1 });
            }
        }
    }
    let n = positions.len();
    let masses = vec![1000.0 * DX * DX * DX; n];
    let densities = vec![1000.0; n];
    let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * H);
    let color_grad = compute_multiphase_color_gradient(
        &positions,
        &masses,
        &densities,
        &phase,
        &neighbors,
        &CubicSplineKernel,
        H,
    );
    Layers {
        positions,
        masses,
        densities,
        phase,
        neighbors,
        color_grad,
    }
}

impl Layers {
    fn forces(&self, temperatures: &[f64], dsigma_dt: f64) -> Vec<Vec3> {
        marangoni_forces(
            &self.positions,
            &self.masses,
            &self.densities,
            temperatures,
            &self.color_grad,
            &self.neighbors,
            H,
            dsigma_dt,
        )
    }

    fn temperatures(&self, field: impl Fn(&Vec3) -> f64) -> Vec<f64> {
        self.positions.iter().map(field).collect()
    }

    /// Particles of the two layers next to the interface, away from the box edges in x and z.
    fn interface_core(&self) -> Vec<usize> {
        (0..self.positions.len())
            .filter(|&i| {
                let p = self.positions[i];
                let iy = (p.y / DX).round() as usize;
                let ix = (p.x / DX).round() as usize;
                let iz = (p.z / DX).round() as usize;
                (iy == NY / 2 - 1 || iy == NY / 2)
                    && (3..=NX - 4).contains(&ix)
                    && (2..=NZ - 3).contains(&iz)
            })
            .collect()
    }
}

#[test]
fn fixture_has_interface_particles_with_unit_normals_along_y() {
    let l = layers();
    let core = l.interface_core();
    assert!(core.len() >= 6);
    for &i in &core {
        let g = l.color_grad[i];
        assert!(
            g.norm() > 1.0,
            "interface particle {i} has |∇c| = {}",
            g.norm()
        );
        assert!(
            g.y.abs() > 0.999 * g.norm(),
            "flat interface normal is ±y: {g:?}"
        );
    }
    // Deep in either phase the colour gradient vanishes.
    for i in 0..l.positions.len() {
        let iy = (l.positions[i].y / DX).round() as usize;
        if iy <= NY / 2 - 3 || iy >= NY / 2 + 2 {
            assert!(
                l.color_grad[i].norm() < 1.0e-12,
                "bulk |∇c| = {}",
                l.color_grad[i].norm()
            );
        }
    }
    assert!(l.phase.contains(&0) && l.phase.contains(&1));
}

#[test]
fn force_points_towards_the_cold_side_for_negative_dsigma_dt() {
    // T grows with +x; σ falls with T, so the interface is pulled towards −x (higher σ).
    let l = layers();
    let gradient = 100.0; // K/m
    let temperatures = l.temperatures(|p| 300.0 + gradient * p.x);
    let forces = l.forces(&temperatures, DSIGMA_DT);
    for &i in &l.interface_core() {
        assert!(
            forces[i].x < 0.0,
            "interface particle {i}: f = {:?}",
            forces[i]
        );
        // f = (dσ/dT) (∂T/∂x) |∇c|; the SPH gradient of a linear field on this lattice
        // is within 1 % of the exact slope.
        let expected = DSIGMA_DT * gradient * l.color_grad[i].norm();
        assert!(
            (forces[i].x - expected).abs() < 0.02 * expected.abs(),
            "particle {i}: f_x = {} vs {expected}",
            forces[i].x
        );
        assert!(forces[i].y.abs() < 1.0e-9 * forces[i].x.abs());
        assert!(forces[i].z.abs() < 1.0e-9 * forces[i].x.abs());
    }
}

#[test]
fn force_flips_with_the_sign_of_dsigma_dt_and_scales_linearly() {
    let l = layers();
    let temperatures = l.temperatures(|p| 300.0 + 40.0 * p.x - 25.0 * p.z);
    let neg = l.forces(&temperatures, -2.0e-4);
    let pos = l.forces(&temperatures, 2.0e-4);
    let double = l.forces(&temperatures, -4.0e-4);
    for &i in &l.interface_core() {
        assert!((neg[i] + pos[i]).norm() < 1.0e-12 * neg[i].norm());
        assert!((double[i] - 2.0 * neg[i]).norm() < 1.0e-12 * neg[i].norm());
        assert!(neg[i].norm() > 0.0);
    }
}

#[test]
fn particles_away_from_the_interface_feel_no_marangoni_force() {
    // A temperature gradient exists everywhere, but only the interface carries a stress.
    let l = layers();
    let temperatures = l.temperatures(|p| 300.0 + 100.0 * p.x + 30.0 * p.y);
    let forces = l.forces(&temperatures, DSIGMA_DT);
    let mut bulk = 0;
    for (i, force) in forces.iter().enumerate() {
        if l.color_grad[i].norm() < 1.0e-12 {
            assert_eq!(*force, Vec3::zeros(), "bulk particle {i}");
            bulk += 1;
        }
    }
    assert!(bulk > l.positions.len() / 3, "{bulk} bulk particles");
    let on_interface = l
        .interface_core()
        .iter()
        .map(|&i| forces[i].norm())
        .fold(0.0_f64, f64::max);
    assert!(on_interface > 1.0e-3, "interface force {on_interface}");
}

#[test]
fn force_is_tangential_to_the_interface() {
    // Only the surface gradient ∇_s T drives the flow: a temperature gradient along
    // the interface normal produces no force at all, and a mixed gradient only the
    // tangential part.
    let l = layers();
    let normal_only = l.temperatures(|p| 300.0 + 80.0 * p.y);
    let forces = l.forces(&normal_only, DSIGMA_DT);
    for &i in &l.interface_core() {
        assert!(
            forces[i].norm() < 1.0e-9 * DSIGMA_DT.abs() * 80.0 * l.color_grad[i].norm(),
            "normal gradient must not push the interface: {:?}",
            forces[i]
        );
    }
    let mixed = l.temperatures(|p| 300.0 + 80.0 * p.y + 100.0 * p.x + 60.0 * p.z);
    let forces = l.forces(&mixed, DSIGMA_DT);
    for &i in &l.interface_core() {
        let n = l.color_grad[i] / l.color_grad[i].norm();
        assert!(
            forces[i].dot(&n).abs() < 1.0e-9 * forces[i].norm(),
            "f·n̂ = {}",
            forces[i].dot(&n)
        );
        let expected_x = DSIGMA_DT * 100.0 * l.color_grad[i].norm();
        let expected_z = DSIGMA_DT * 60.0 * l.color_grad[i].norm();
        assert!((forces[i].x - expected_x).abs() < 0.02 * expected_x.abs());
        assert!((forces[i].z - expected_z).abs() < 0.02 * expected_z.abs());
    }
}

#[test]
fn uniform_temperature_gives_no_force_even_on_the_interface() {
    let l = layers();
    let forces = l.forces(&vec![310.0; l.positions.len()], DSIGMA_DT);
    for f in forces {
        assert!(f.norm() < 1.0e-12, "{f:?}");
    }
}

#[test]
fn faint_colour_gradients_below_the_interface_threshold_are_ignored() {
    let l = layers();
    let temperatures = l.temperatures(|p| 300.0 + 100.0 * p.x);
    let faint = vec![Vec3::new(0.0, 0.5e-3 / H, 0.0); l.positions.len()];
    let forces = marangoni_forces(
        &l.positions,
        &l.masses,
        &l.densities,
        &temperatures,
        &faint,
        &l.neighbors,
        H,
        DSIGMA_DT,
    );
    assert!(forces.iter().all(|f| *f == Vec3::zeros()));
}

#[test]
#[should_panic(expected = "colour gradient per particle")]
fn colour_gradient_length_must_match_the_particle_count() {
    let l = layers();
    let temperatures = vec![300.0; l.positions.len()];
    let _ = marangoni_forces(
        &l.positions,
        &l.masses,
        &l.densities,
        &temperatures,
        &[Vec3::zeros()],
        &l.neighbors,
        H,
        DSIGMA_DT,
    );
}
