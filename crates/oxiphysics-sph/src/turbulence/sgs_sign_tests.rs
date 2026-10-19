//! Sign, scale and consistency tests for the sub-particle-scale (SPS) stress and the
//! Smagorinsky SGS accelerations.
//!
//! Convention used throughout: `∇_iW_ij = (dW/dr)(x_i − x_j)/|x_i − x_j|` with the
//! cubic-spline kernel, so the gradient points from `i` towards `j`.  The SPS stress
//! `τ = ρ[2ν_t S_dev − ⅔ C_I Δ²|S|² I]` is *added* to the momentum equation, while the
//! Smagorinsky SGS stress `τ = −2ρν_t S` is *subtracted* from it; both therefore
//! produce a viscous acceleration that opposes the relative velocity of the particles.

use super::types::{SmagorinskyLes, SpsModel, TurbulentDiffusionSph};
use crate::kernel::{CubicSplineKernel, SphKernel};
use crate::neighbor::SpatialHash;
use oxiphysics_core::math::Vec3;

type Mat3 = [[f64; 3]; 3];

const RHO: f64 = 1000.0;

fn zero_mat() -> Mat3 {
    [[0.0; 3]; 3]
}

fn trace(a: &Mat3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

/// `|S| = sqrt(2 S_ij S_ij)`.
fn strain_magnitude(s: &Mat3) -> f64 {
    let mut sum = 0.0;
    for row in s {
        for v in row {
            sum += v * v;
        }
    }
    (2.0 * sum).sqrt()
}

fn rel_err(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs()
}

/// `∇_i W_ij` of the cubic-spline kernel with smoothing length `h`.
fn kernel_grad(xi: [f64; 3], xj: [f64; 3], h: f64) -> [f64; 3] {
    let d = [xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]];
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if r < 1.0e-14 {
        return [0.0; 3];
    }
    let dw_dr = CubicSplineKernel.grad_w(r, h);
    [dw_dr * d[0] / r, dw_dr * d[1] / r, dw_dr * d[2] / r]
}

/// Symmetric strain rate of particle `i` from neighbour contributions
/// `(V_j, v_j − v_i, ∇_iW_ij)`.
fn strain_from_pairs(pairs: &[(f64, [f64; 3], [f64; 3])]) -> Mat3 {
    let mut dv = zero_mat();
    for &(vol, dvel, grad) in pairs {
        for (row, &dv_a) in dv.iter_mut().zip(dvel.iter()) {
            for (entry, &g_b) in row.iter_mut().zip(grad.iter()) {
                *entry += vol * dv_a * g_b;
            }
        }
    }
    let mut s = zero_mat();
    for a in 0..3 {
        for b in 0..3 {
            s[a][b] = 0.5 * (dv[a][b] + dv[b][a]);
        }
    }
    s
}

/// `(A · g)_α = Σ_β A_αβ g_β`.
fn mat_vec(a: &Mat3, g: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for (component, row) in out.iter_mut().zip(a.iter()) {
        *component = dot(*row, g);
    }
    out
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// SPS stress (F4703)
// ---------------------------------------------------------------------------

fn simple_shear_strain(a: f64) -> Mat3 {
    [[0.0, a, 0.0], [a, 0.0, 0.0], [0.0, 0.0, 0.0]]
}

#[test]
fn sps_stress_for_simple_shear_matches_dalrymple_rogers() {
    // S_xy = 5 s⁻¹  ⇒  |S| = 10 s⁻¹, ν_t = (0.12·0.1)²·10 = 1.44e-3 m²/s.
    let model = SpsModel::new(0.12, 0.0066, 0.1);
    let s = simple_shear_strain(5.0);
    let tau = model.compute_sps_stress(&s, RHO);

    let nu_t = (0.12_f64 * 0.1).powi(2) * 10.0;
    assert!(rel_err(nu_t, 1.44e-3) < 1.0e-12);
    // τ_xy = 2 ρ ν_t S_xy = 14.4 Pa.
    assert!(rel_err(tau[0][1], 2.0 * RHO * nu_t * 5.0) < 1.0e-12);
    assert!(rel_err(tau[0][1], 14.4) < 1.0e-12, "τ_xy = {}", tau[0][1]);
    assert_eq!(tau[0][1], tau[1][0]);
    // Isotropic part −⅔ ρ C_I Δ² |S|² = −⅔ · 1000 · 0.0066 · 0.01 · 100 = −4.4 Pa.
    for (k, row) in tau.iter().enumerate() {
        assert!(rel_err(row[k], -4.4) < 1.0e-12, "τ_{k}{k} = {}", row[k]);
    }
    assert_eq!(tau[0][2], 0.0);
    assert_eq!(tau[1][2], 0.0);
    // trace(τ) = −2 ρ C_I Δ² |S|².
    assert!(rel_err(trace(&tau), -2.0 * RHO * 0.0066 * 0.01 * 100.0) < 1.0e-12);
}

#[test]
fn sps_stress_is_proportional_to_density() {
    let model = SpsModel::default();
    let s: Mat3 = [[0.4, 1.2, -0.3], [1.2, -0.9, 0.5], [-0.3, 0.5, 0.2]];
    let light = model.compute_sps_stress(&s, 1000.0);
    let heavy = model.compute_sps_stress(&s, 3000.0);
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (heavy[a][b] - 3.0 * light[a][b]).abs() < 1.0e-12 * light[a][b].abs().max(1.0),
                "τ[{a}][{b}] must scale with ρ: {} vs {}",
                heavy[a][b],
                light[a][b]
            );
        }
    }
    assert!(
        light[0][1].abs() > 1.0e-6,
        "the reference stress is non-trivial"
    );
}

#[test]
fn sps_isotropic_term_is_linear_in_the_yoshizawa_constant() {
    let s = simple_shear_strain(3.0);
    let single = SpsModel::new(0.12, 0.0066, 0.1).compute_sps_stress(&s, RHO);
    let doubled = SpsModel::new(0.12, 0.0132, 0.1).compute_sps_stress(&s, RHO);
    // C_I enters linearly (not squared): doubling it doubles the isotropic stress.
    assert!(rel_err(trace(&doubled), 2.0 * trace(&single)) < 1.0e-12);
    // The viscous (deviatoric) part does not depend on C_I.
    assert!((doubled[0][1] - single[0][1]).abs() < 1.0e-12);
}

#[test]
fn sps_stress_trace_does_not_depend_on_the_normal_strain() {
    // S_kk ≠ 0: the viscous part is trace-free, so trace(τ) = −2 ρ C_I Δ² |S|² and the
    // deviatoric normal stresses follow 2 ρ ν_t (S_ii − S_jj).
    let model = SpsModel::new(0.12, 0.0066, 0.1);
    let s: Mat3 = [[2.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, 1.0]];
    let tau = model.compute_sps_stress(&s, RHO);
    let s_mag = strain_magnitude(&s);
    let nu_t = (0.12_f64 * 0.1).powi(2) * s_mag;
    assert!(rel_err(trace(&tau), -2.0 * RHO * 0.0066 * 0.01 * s_mag * s_mag) < 1.0e-12);
    assert!(rel_err(tau[0][0] - tau[1][1], 2.0 * RHO * nu_t * 2.5) < 1.0e-12);
    assert!(rel_err(tau[0][0] - tau[2][2], 2.0 * RHO * nu_t * 1.0) < 1.0e-12);
}

#[test]
fn sps_stress_vanishes_without_strain() {
    let tau = SpsModel::default().compute_sps_stress(&zero_mat(), RHO);
    assert_eq!(tau, zero_mat());
}

#[test]
fn sps_force_of_a_pair_matches_the_hand_computed_stress_divergence() {
    // Particle 0 at rest, particle 1 moving with (U, 0, 0), separated along y.
    let model = SpsModel::new(0.12, 0.0066, 0.05);
    let h = 0.065;
    let (m, vol) = (RHO * 0.05_f64.powi(3), 0.05_f64.powi(3));
    let positions = [[0.0, 0.0, 0.0], [0.0, 0.05, 0.0]];
    let velocities = [[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]];
    let g01 = kernel_grad(positions[0], positions[1], h);
    let g10 = kernel_grad(positions[1], positions[0], h);
    let neighbors = vec![vec![1_usize], vec![0_usize]];
    let grads = vec![vec![g01], vec![g10]];
    let forces = model.apply_sps_forces(
        &positions,
        &velocities,
        &[m, m],
        &[RHO, RHO],
        &neighbors,
        &grads,
    );

    let s0 = strain_from_pairs(&[(vol, [0.8, 0.0, 0.0], g01)]);
    let s1 = strain_from_pairs(&[(vol, [-0.8, 0.0, 0.0], g10)]);
    let tau0 = model.compute_sps_stress(&s0, RHO);
    let tau1 = model.compute_sps_stress(&s1, RHO);
    for (i, force) in forces.iter().enumerate() {
        let (g, own, other) = if i == 0 {
            (g01, &tau0, &tau1)
        } else {
            (g10, &tau1, &tau0)
        };
        let a = mat_vec(own, g);
        let b = mat_vec(other, g);
        for c in 0..3 {
            let expected = m * m * (a[c] + b[c]) / (RHO * RHO);
            assert!(
                (force[c] - expected).abs() <= 1.0e-12 * expected.abs().max(1.0e-30),
                "force[{i}][{c}] = {} vs {expected}",
                force[c]
            );
        }
    }
    // Internal forces: momentum is conserved and the force on the fast particle
    // decelerates it (viscous), the force on the slow one accelerates it.
    for (f0, f1) in forces[0].iter().zip(forces[1].iter()) {
        assert!((f0 + f1).abs() <= 1.0e-12 * f0.abs().max(1.0e-30));
    }
    assert!(
        forces[1][0] < 0.0,
        "fast particle must be slowed: {}",
        forces[1][0]
    );
    assert!(
        forces[0][0] > 0.0,
        "slow particle must be dragged along: {}",
        forces[0][0]
    );
}

// ---------------------------------------------------------------------------
// Sheared lattice: u_x = U sin(k y)
// ---------------------------------------------------------------------------

/// Cubic lattice carrying the velocity field `u_x = U sin(k y)` (a periodic shear
/// layer).  The core is one x/z cell and one wavelength in y; a halo of `HALO` layers
/// surrounds it so that the second neighbourhood of every core particle is complete.
struct ShearedLattice {
    positions: Vec<[f64; 3]>,
    velocities: Vec<[f64; 3]>,
    masses: Vec<f64>,
    densities: Vec<f64>,
    neighbors: Vec<Vec<usize>>,
    grads: Vec<Vec<[f64; 3]>>,
    core: Vec<usize>,
    dx: f64,
    h: f64,
    k: f64,
    amplitude: f64,
}

const HALO: i32 = 6;

fn sheared_lattice(wavelength_cells: i32, amplitude: f64) -> ShearedLattice {
    let dx = 0.1;
    let h = 1.3 * dx;
    let k = 2.0 * std::f64::consts::PI / (wavelength_cells as f64 * dx);
    let mut positions = Vec::new();
    let mut core = Vec::new();
    for ix in -HALO..(1 + HALO) {
        for iy in -HALO..(wavelength_cells + HALO) {
            for iz in -HALO..(1 + HALO) {
                if ix == 0 && iz == 0 && (0..wavelength_cells).contains(&iy) {
                    core.push(positions.len());
                }
                positions.push([ix as f64 * dx, iy as f64 * dx, iz as f64 * dx]);
            }
        }
    }
    let velocities: Vec<[f64; 3]> = positions
        .iter()
        .map(|p| [amplitude * (k * p[1]).sin(), 0.0, 0.0])
        .collect();
    let n = positions.len();
    let hash_positions: Vec<Vec3> = positions
        .iter()
        .map(|p| Vec3::new(p[0], p[1], p[2]))
        .collect();
    let neighbors = SpatialHash::find_all_neighbors(&hash_positions, 2.0 * h);
    let grads = (0..n)
        .map(|i| {
            neighbors[i]
                .iter()
                .map(|&j| kernel_grad(positions[i], positions[j], h))
                .collect()
        })
        .collect();
    ShearedLattice {
        positions,
        velocities,
        masses: vec![RHO * dx * dx * dx; n],
        densities: vec![RHO; n],
        neighbors,
        grads,
        core,
        dx,
        h,
        k,
        amplitude,
    }
}

impl ShearedLattice {
    /// Resolved strain rate of every particle from the SPH velocity gradient.
    fn strain_rates(&self) -> Vec<Mat3> {
        (0..self.positions.len())
            .map(|i| {
                let pairs: Vec<(f64, [f64; 3], [f64; 3])> = self.neighbors[i]
                    .iter()
                    .zip(self.grads[i].iter())
                    .map(|(&j, &g)| {
                        let dv = [
                            self.velocities[j][0] - self.velocities[i][0],
                            self.velocities[j][1] - self.velocities[i][1],
                            self.velocities[j][2] - self.velocities[i][2],
                        ];
                        (self.masses[j] / self.densities[j], dv, g)
                    })
                    .collect();
                strain_from_pairs(&pairs)
            })
            .collect()
    }

    fn y(&self, i: usize) -> f64 {
        self.positions[i][1]
    }
}

/// Analytic `∂_y(2 ν_t S_xy)` for `u_x = U sin(k y)`, `ν_t = (C Δ)² |S|`.
fn analytic_viscous_ax(c_delta_sq: f64, amplitude: f64, k: f64, y: f64) -> f64 {
    -2.0 * c_delta_sq * amplitude * amplitude * k.powi(3) * (k * y).cos().abs() * (k * y).sin()
}

#[test]
fn sheared_lattice_has_the_expected_strain_rate() {
    // Guards the fixture: S_xy = ½ U k cos(k y) at the core particles (to the accuracy
    // of the SPH gradient).
    let lattice = sheared_lattice(20, 1.0);
    let strain = lattice.strain_rates();
    let i = lattice.core[2]; // y = 2 dx
    let expected = 0.5 * lattice.amplitude * lattice.k * (lattice.k * lattice.y(i)).cos();
    assert!(
        rel_err(strain[i][0][1], expected) < 0.05,
        "S_xy = {} vs {expected}",
        strain[i][0][1]
    );
    assert!(strain[i][0][0].abs() < 1.0e-9 && strain[i][1][1].abs() < 1.0e-9);
}

#[test]
fn sps_force_on_sheared_lattice_matches_the_stress_divergence() {
    // a = (1/ρ) ∇·τ with τ_xy = 2ρν_t S_xy and τ_yy = −⅔ρ C_I Δ²|S|², so
    //   a_x = −2 (C_sΔ)² U² k³ |cos ky| sin ky
    //   a_y = +4/3 C_I Δ² U² k³ sin ky cos ky.
    let lattice = sheared_lattice(20, 1.0);
    let model = SpsModel::new(0.12, 0.0066, lattice.dx);
    let forces = model.apply_sps_forces(
        &lattice.positions,
        &lattice.velocities,
        &lattice.masses,
        &lattice.densities,
        &lattice.neighbors,
        &lattice.grads,
    );
    let (u, k, dx) = (lattice.amplitude, lattice.k, lattice.dx);
    let c1 = (0.12 * dx) * (0.12 * dx);
    let mut compared = 0;
    for &i in &lattice.core {
        let y = lattice.y(i);
        let accel = [
            forces[i][0] / lattice.masses[i],
            forces[i][1] / lattice.masses[i],
            forces[i][2] / lattice.masses[i],
        ];
        let ax = analytic_viscous_ax(c1, u, k, y);
        let ay = 4.0 / 3.0 * 0.0066 * dx * dx * u * u * k.powi(3) * (k * y).sin() * (k * y).cos();
        // Skip the nodes of the analytic profile.
        if ax.abs() > 1.0e-3 {
            let ratio = accel[0] / ax;
            assert!(
                ratio > 0.7 && ratio < 1.05,
                "y = {y}: a_x = {} vs analytic {ax} (ratio {ratio})",
                accel[0]
            );
            compared += 1;
        }
        if ay.abs() > 5.0e-4 {
            let ratio = accel[1] / ay;
            assert!(
                ratio > 0.7 && ratio < 1.05,
                "y = {y}: a_y = {} vs analytic {ay} (ratio {ratio})",
                accel[1]
            );
        }
        assert!(accel[2].abs() < 1.0e-9, "no spanwise force: {}", accel[2]);
    }
    assert!(compared >= 10, "only {compared} particles were compared");
}

#[test]
fn sps_force_is_dissipative_for_the_shear_layer() {
    let lattice = sheared_lattice(20, 1.0);
    let model = SpsModel::new(0.12, 0.0066, lattice.dx);
    let forces = model.apply_sps_forces(
        &lattice.positions,
        &lattice.velocities,
        &lattice.masses,
        &lattice.densities,
        &lattice.neighbors,
        &lattice.grads,
    );
    let power: f64 = lattice
        .core
        .iter()
        .map(|&i| dot(lattice.velocities[i], forces[i]))
        .sum();
    assert!(
        power < 0.0,
        "SPS work on the shear flow must be negative: {power}"
    );
    // The force decelerates the fast side and accelerates the slow one.
    let fast = lattice.core[2]; // y = 2 dx: u_x > 0, force along −x
    let slow = lattice.core[12]; // y = 12 dx: u_x < 0, force along +x
    assert!(lattice.velocities[fast][0] > 0.0 && forces[fast][0] < 0.0);
    assert!(lattice.velocities[slow][0] < 0.0 && forces[slow][0] > 0.0);
}

#[test]
fn sps_acceleration_does_not_depend_on_the_density_scale() {
    // The stress is ∝ ρ, so the acceleration (1/ρ)∇·τ is independent of the density as
    // long as the particle volumes V = m/ρ (and hence the strain rates) are unchanged: a
    // water-like fluid (ρ = 1000, mass m) and a 1000× lighter one (ρ = 1, mass m/1000)
    // must accelerate identically.  A stress that ignored ρ would make the light fluid
    // accelerate 1000× faster.
    let lattice = sheared_lattice(20, 1.0);
    let model = SpsModel::new(0.12, 0.0066, lattice.dx);
    let n = lattice.positions.len();
    let heavy_masses = lattice.masses.clone();
    let light_masses: Vec<f64> = lattice.masses.iter().map(|m| m * 1.0e-3).collect();
    let heavy = model.apply_sps_forces(
        &lattice.positions,
        &lattice.velocities,
        &heavy_masses,
        &vec![1000.0; n],
        &lattice.neighbors,
        &lattice.grads,
    );
    let light = model.apply_sps_forces(
        &lattice.positions,
        &lattice.velocities,
        &light_masses,
        &vec![1.0; n],
        &lattice.neighbors,
        &lattice.grads,
    );
    let mut compared = 0;
    for &i in &lattice.core {
        for c in 0..2 {
            let a_heavy = heavy[i][c] / heavy_masses[i];
            let a_light = light[i][c] / light_masses[i];
            if a_heavy.abs() > 1.0e-4 {
                assert!(
                    rel_err(a_light, a_heavy) < 1.0e-9,
                    "particle {i}, component {c}: {a_light} vs {a_heavy}"
                );
                compared += 1;
            }
        }
    }
    assert!(
        compared >= 20,
        "only {compared} accelerations were compared"
    );
}

// ---------------------------------------------------------------------------
// Smagorinsky SGS acceleration (F4704)
// ---------------------------------------------------------------------------

/// Two particles 0.05 m apart (along y); particle 1 moves with `(U, 0, 0)`.
struct ShearPair {
    velocities: [[f64; 3]; 2],
    g01: [f64; 3],
    g10: [f64; 3],
    mass: f64,
    strain: [Mat3; 2],
}

fn shear_pair(speed: f64) -> ShearPair {
    let dx = 0.05;
    let h = 0.065;
    let mass = RHO * dx * dx * dx;
    let vol = dx * dx * dx;
    let positions = [[0.0, 0.0, 0.0], [0.0, dx, 0.0]];
    let velocities = [[0.0; 3], [speed, 0.0, 0.0]];
    let g01 = kernel_grad(positions[0], positions[1], h);
    let g10 = kernel_grad(positions[1], positions[0], h);
    let strain = [
        strain_from_pairs(&[(vol, [speed, 0.0, 0.0], g01)]),
        strain_from_pairs(&[(vol, [-speed, 0.0, 0.0], g10)]),
    ];
    ShearPair {
        velocities,
        g01,
        g10,
        mass,
        strain,
    }
}

#[test]
fn smagorinsky_acceleration_opposes_the_relative_velocity_of_a_pair() {
    let pair = shear_pair(0.8);
    let model = SmagorinskyLes::new(0.17, 0.05);
    assert!(pair.strain[0][0][1].abs() > 1.0e-3, "the pair is sheared");

    let tau0 = model.sgs_stress(&pair.strain[0], RHO);
    let tau1 = model.sgs_stress(&pair.strain[1], RHO);
    let a0 = model.sgs_acceleration(&tau0, RHO, &[tau1], &[RHO], &[pair.mass], &[pair.g01]);
    let a1 = model.sgs_acceleration(&tau1, RHO, &[tau0], &[RHO], &[pair.mass], &[pair.g10]);

    // Each particle is pulled towards the velocity of the other one.
    let dv01 = [0.8, 0.0, 0.0]; // v_1 − v_0
    assert!(dot(a0, dv01) > 0.0, "a_0 = {a0:?} must follow v_1 − v_0");
    assert!(dot(a1, dv01) < 0.0, "a_1 = {a1:?} must follow v_0 − v_1");
    // Equal masses: the internal force pair is equal and opposite.
    for c in 0..3 {
        assert!((a0[c] + a1[c]).abs() <= 1.0e-12 * a0[c].abs().max(1.0e-30));
    }
    // The momentum-conserving pair exchange removes resolved kinetic energy.
    let dt = 1.0e-3;
    let kinetic = |v: [[f64; 3]; 2]| 0.5 * pair.mass * (dot(v[0], v[0]) + dot(v[1], v[1]));
    let before = kinetic(pair.velocities);
    let mut after_v = pair.velocities;
    for c in 0..3 {
        after_v[0][c] += a0[c] * dt;
        after_v[1][c] += a1[c] * dt;
    }
    let after = kinetic(after_v);
    assert!(
        after < before,
        "kinetic energy must drop: {before} -> {after}"
    );
}

#[test]
fn smagorinsky_acceleration_matches_the_hand_computed_stress_divergence() {
    let pair = shear_pair(0.8);
    let model = SmagorinskyLes::new(0.17, 0.05);
    let tau0 = model.sgs_stress(&pair.strain[0], RHO);
    let tau1 = model.sgs_stress(&pair.strain[1], RHO);
    let a0 = model.sgs_acceleration(&tau0, RHO, &[tau1], &[RHO], &[pair.mass], &[pair.g01]);
    // a_0 = −m (τ_0/ρ² + τ_1/ρ²)·∇_0W_01
    let own = mat_vec(&tau0, pair.g01);
    let other = mat_vec(&tau1, pair.g01);
    for c in 0..3 {
        let expected = -pair.mass * (own[c] + other[c]) / (RHO * RHO);
        assert!(
            (a0[c] - expected).abs() <= 1.0e-12 * expected.abs().max(1.0e-30),
            "a[{c}] = {} vs {expected}",
            a0[c]
        );
    }
    // τ = −2ρν_t S, so a = +2 m (ν_t,0 S_0 + ν_t,1 S_1)/ρ · ∇W: component along +x positive.
    assert!(a0[0] > 0.0);
}

#[test]
fn smagorinsky_acceleration_cancels_for_a_uniform_stress() {
    let model = SmagorinskyLes::new(0.17, 0.1);
    let tau: Mat3 = [[1.0, 0.3, 0.0], [0.3, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let nbr_tau = vec![tau; 2];
    let nbr_rho = vec![RHO; 2];
    let nbr_mass = vec![0.001_f64; 2];
    let grads = vec![[1.0, 0.5, 0.0_f64], [-1.0, -0.5, 0.0]];
    let acc = model.sgs_acceleration(&tau, RHO, &nbr_tau, &nbr_rho, &nbr_mass, &grads);
    for a in acc {
        assert!(a.abs() < 1.0e-18, "opposite kernel gradients cancel: {a}");
    }
}

#[test]
fn smagorinsky_acceleration_on_sheared_lattice_matches_the_stress_divergence() {
    // a = −(1/ρ) ∇·τ = ∂_y(2 ν_t S_xy) = −2 (C_sΔ)² U² k³ |cos ky| sin ky.
    let lattice = sheared_lattice(20, 1.0);
    let model = SmagorinskyLes::new(0.17, lattice.dx);
    let strain = lattice.strain_rates();
    let stress: Vec<Mat3> = strain
        .iter()
        .zip(lattice.densities.iter())
        .map(|(s, &rho)| model.sgs_stress(s, rho))
        .collect();
    let c1 = (0.17 * lattice.dx) * (0.17 * lattice.dx);
    let mut compared = 0;
    let mut power = 0.0;
    for &i in &lattice.core {
        let nbr_tau: Vec<Mat3> = lattice.neighbors[i].iter().map(|&j| stress[j]).collect();
        let nbr_rho: Vec<f64> = lattice.neighbors[i]
            .iter()
            .map(|&j| lattice.densities[j])
            .collect();
        let nbr_mass: Vec<f64> = lattice.neighbors[i]
            .iter()
            .map(|&j| lattice.masses[j])
            .collect();
        let acc = model.sgs_acceleration(
            &stress[i],
            lattice.densities[i],
            &nbr_tau,
            &nbr_rho,
            &nbr_mass,
            &lattice.grads[i],
        );
        let ax = analytic_viscous_ax(c1, lattice.amplitude, lattice.k, lattice.y(i));
        if ax.abs() > 1.0e-3 {
            let ratio = acc[0] / ax;
            assert!(
                ratio > 0.7 && ratio < 1.05,
                "y = {}: a_x = {} vs analytic {ax} (ratio {ratio})",
                lattice.y(i),
                acc[0]
            );
            compared += 1;
        }
        power += dot(lattice.velocities[i], acc);
    }
    assert!(compared >= 10, "only {compared} particles were compared");
    assert!(
        power < 0.0,
        "SGS work on the shear layer must be negative: {power}"
    );
}

#[test]
fn turbulent_diffusion_acceleration_opposes_the_relative_velocity_of_a_pair() {
    let pair = shear_pair(0.8);
    let model = TurbulentDiffusionSph::new(0.17, 0.05);
    let a0 = model.turbulent_acceleration(
        pair.strain[0],
        pair.strain[1],
        RHO,
        RHO,
        pair.mass,
        pair.g01,
    );
    let a1 = model.turbulent_acceleration(
        pair.strain[1],
        pair.strain[0],
        RHO,
        RHO,
        pair.mass,
        pair.g10,
    );
    assert!(a0[0] > 0.0, "a_0 must follow the faster neighbour: {a0:?}");
    assert!(a1[0] < 0.0, "a_1 must follow the slower neighbour: {a1:?}");
    let dt = 1.0e-3;
    let mut v = pair.velocities;
    let before = 0.5 * pair.mass * (dot(v[0], v[0]) + dot(v[1], v[1]));
    for c in 0..3 {
        v[0][c] += a0[c] * dt;
        v[1][c] += a1[c] * dt;
    }
    let after = 0.5 * pair.mass * (dot(v[0], v[0]) + dot(v[1], v[1]));
    assert!(
        after < before,
        "kinetic energy must drop: {before} -> {after}"
    );
}

#[test]
fn turbulent_diffusion_agrees_with_the_smagorinsky_acceleration() {
    // Summing the pair contributions of TurbulentDiffusionSph gives exactly the
    // SmagorinskyLes acceleration (same stress closure, same sign convention).
    let lattice = sheared_lattice(20, 1.0);
    let delta = lattice.dx;
    let les = SmagorinskyLes::new(0.17, delta);
    let tds = TurbulentDiffusionSph::new(0.17, delta);
    let strain = lattice.strain_rates();
    let stress: Vec<Mat3> = strain.iter().map(|s| les.sgs_stress(s, RHO)).collect();
    for &i in lattice.core.iter().take(8) {
        let nbr_tau: Vec<Mat3> = lattice.neighbors[i].iter().map(|&j| stress[j]).collect();
        let nbr_rho = vec![RHO; nbr_tau.len()];
        let nbr_mass: Vec<f64> = lattice.neighbors[i]
            .iter()
            .map(|&j| lattice.masses[j])
            .collect();
        let reference = les.sgs_acceleration(
            &stress[i],
            RHO,
            &nbr_tau,
            &nbr_rho,
            &nbr_mass,
            &lattice.grads[i],
        );
        let mut summed = [0.0; 3];
        for (k, &j) in lattice.neighbors[i].iter().enumerate() {
            let pair = tds.turbulent_acceleration(
                strain[i],
                strain[j],
                RHO,
                RHO,
                lattice.masses[j],
                lattice.grads[i][k],
            );
            for c in 0..3 {
                summed[c] += pair[c];
            }
        }
        for c in 0..3 {
            // Relative agreement, with an absolute floor for the nodes of the profile where
            // the acceleration cancels to rounding noise (the typical value is 1e-3).
            assert!(
                (summed[c] - reference[c]).abs() <= 1.0e-10 * reference[c].abs() + 1.0e-14,
                "particle {i}, component {c}: {} vs {}",
                summed[c],
                reference[c]
            );
        }
    }
    // Sanity: the lattice kernel radius and spacing are the documented ones.
    assert!((lattice.h - 1.3 * lattice.dx).abs() < 1.0e-15);
}
