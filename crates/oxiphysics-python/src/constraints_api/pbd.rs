// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Position-based dynamics (PBD/XPBD) solver of the constraint API.
//!
//! [`PyPbdSolver`] integrates particles with gravity and projects distance
//! (stretch) and tetrahedron volume constraints with XPBD (compliance and
//! accumulated Lagrange multipliers per substep).

use super::cross;
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PyPbdSolver
// ---------------------------------------------------------------------------

/// PBD/XPBD constraint type.
///
/// Every kind listed here has a constructor on [`PyPbdSolver`] and is solved
/// by [`PyPbdSolver::substep`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PbdConstraintType {
    /// Distance / stretch constraint between two particles.
    Stretch { rest_length: f64, compliance: f64 },
    /// Volume conservation constraint of a tetrahedron (four particles).
    ///
    /// `rest_volume` is the *signed* volume
    /// `V0 = (x1 - x0) . ((x2 - x0) x (x3 - x0)) / 6` of the tetrahedron when
    /// the constraint was created, so the orientation is preserved.
    Volume { rest_volume: f64, compliance: f64 },
}

/// A PBD/XPBD constraint connecting particle indices.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PbdConstraint {
    /// Particle indices involved (two for a stretch, four for a volume constraint).
    pub particles: Vec<usize>,
    /// Constraint type and parameters.
    pub kind: PbdConstraintType,
    /// Accumulated Lagrange multiplier (XPBD).
    pub lambda: f64,
}

/// Constraint definitions rejected by [`PyPbdSolver`].
///
/// A bad index surfaces in Python as an `IndexError`, a repeated particle as
/// a `ValueError`.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PbdConstraintError {
    /// A particle index is not smaller than the number of particles.
    #[error("particle index {index} is out of range for a solver with {particle_count} particles")]
    ParticleOutOfRange {
        /// The offending index.
        index: usize,
        /// Number of particles the solver holds.
        particle_count: usize,
    },
    /// The same particle appears twice in one constraint.
    #[error("constraint particles must be distinct, but particle {0} appears twice")]
    DuplicateParticle(usize),
}

impl From<PbdConstraintError> for PyErr {
    fn from(err: PbdConstraintError) -> PyErr {
        let message = err.to_string();
        match err {
            PbdConstraintError::ParticleOutOfRange { .. } => {
                pyo3::exceptions::PyIndexError::new_err(message)
            }
            PbdConstraintError::DuplicateParticle(_) => {
                pyo3::exceptions::PyValueError::new_err(message)
            }
        }
    }
}

/// Range of `positions` holding particle `i`, or `None` on index overflow.
fn particle_range(i: usize) -> Option<std::ops::Range<usize>> {
    let start = i.checked_mul(3)?;
    Some(start..start.checked_add(3)?)
}

/// Position of particle `i`, or `None` if `positions` does not hold it.
fn particle_position(positions: &[f64], i: usize) -> Option<[f64; 3]> {
    match positions.get(particle_range(i)?) {
        Some(&[x, y, z]) => Some([x, y, z]),
        _ => None,
    }
}

/// Add `delta` to the position of particle `i` (no-op if it does not exist).
fn displace_particle(positions: &mut [f64], i: usize, delta: [f64; 3]) {
    if let Some(range) = particle_range(i)
        && let Some(slot) = positions.get_mut(range)
    {
        for (coordinate, d) in slot.iter_mut().zip(delta) {
            *coordinate += d;
        }
    }
}

fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Signed volume `(x1 - x0) . ((x2 - x0) x (x3 - x0)) / 6` of a tetrahedron.
fn signed_tet_volume(x: &[[f64; 3]; 4]) -> f64 {
    let e1 = vec_sub(x[1], x[0]);
    let e2 = vec_sub(x[2], x[0]);
    let e3 = vec_sub(x[3], x[0]);
    vec_dot(e1, cross(e2, e3)) / 6.0
}

/// One XPBD iteration of a distance constraint between particles `[i, j]`.
fn solve_stretch(
    positions: &mut [f64],
    inv_masses: &[f64],
    particles: &[usize],
    rest_length: f64,
    compliance: f64,
    h: f64,
    lambda: &mut f64,
) {
    let (Some(&i), Some(&j)) = (particles.first(), particles.get(1)) else {
        return;
    };
    let (Some(xi), Some(xj), Some(&wi), Some(&wj)) = (
        particle_position(positions, i),
        particle_position(positions, j),
        inv_masses.get(i),
        inv_masses.get(j),
    ) else {
        return;
    };
    let w_sum = wi + wj;
    if w_sum < 1e-15 {
        return;
    }
    let d = vec_sub(xi, xj);
    let len = vec_dot(d, d).sqrt();
    if len < 1e-12 {
        return;
    }
    let alpha = compliance / (h * h);
    let c_val = len - rest_length;
    let d_lambda = (-c_val - alpha * *lambda) / (w_sum + alpha);
    *lambda += d_lambda;
    let n = [d[0] / len, d[1] / len, d[2] / len];
    displace_particle(
        positions,
        i,
        [
            wi * d_lambda * n[0],
            wi * d_lambda * n[1],
            wi * d_lambda * n[2],
        ],
    );
    displace_particle(
        positions,
        j,
        [
            -wj * d_lambda * n[0],
            -wj * d_lambda * n[1],
            -wj * d_lambda * n[2],
        ],
    );
}

/// One XPBD iteration of a tetrahedron volume constraint on `[p0, p1, p2, p3]`.
///
/// With `V = (x1 - x0) . ((x2 - x0) x (x3 - x0)) / 6` and `C = V - V0` the
/// gradients are `grad_1 = (x2 - x0) x (x3 - x0) / 6`,
/// `grad_2 = (x3 - x0) x (x1 - x0) / 6`, `grad_3 = (x1 - x0) x (x2 - x0) / 6`
/// and `grad_0 = -(grad_1 + grad_2 + grad_3)`.
fn solve_volume(
    positions: &mut [f64],
    inv_masses: &[f64],
    particles: &[usize],
    rest_volume: f64,
    compliance: f64,
    h: f64,
    lambda: &mut f64,
) {
    let &[p0, p1, p2, p3, ..] = particles else {
        return;
    };
    let ids = [p0, p1, p2, p3];
    let mut x = [[0.0f64; 3]; 4];
    let mut w = [0.0f64; 4];
    for (slot, &id) in ids.iter().enumerate() {
        let (Some(position), Some(&inv_mass)) =
            (particle_position(positions, id), inv_masses.get(id))
        else {
            return;
        };
        x[slot] = position;
        w[slot] = inv_mass;
    }
    let e1 = vec_sub(x[1], x[0]);
    let e2 = vec_sub(x[2], x[0]);
    let e3 = vec_sub(x[3], x[0]);
    let scale = 1.0 / 6.0;
    let g1 = cross(e2, e3).map(|v| v * scale);
    let g2 = cross(e3, e1).map(|v| v * scale);
    let g3 = cross(e1, e2).map(|v| v * scale);
    let g0 = [
        -(g1[0] + g2[0] + g3[0]),
        -(g1[1] + g2[1] + g3[1]),
        -(g1[2] + g2[2] + g3[2]),
    ];
    let grads = [g0, g1, g2, g3];
    let weighted_norm_sq: f64 = grads
        .iter()
        .zip(w)
        .map(|(g, inv_mass)| inv_mass * vec_dot(*g, *g))
        .sum();
    let alpha = compliance / (h * h);
    let denominator = weighted_norm_sq + alpha;
    if denominator < 1e-30 {
        return;
    }
    let c_val = signed_tet_volume(&x) - rest_volume;
    let d_lambda = (-c_val - alpha * *lambda) / denominator;
    *lambda += d_lambda;
    for ((id, g), inv_mass) in ids.iter().zip(grads).zip(w) {
        let k = inv_mass * d_lambda;
        displace_particle(positions, *id, [k * g[0], k * g[1], k * g[2]]);
    }
}

/// Position-based dynamics (XPBD) solver.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPbdSolver {
    /// Current particle positions (flat \[x0, y0, z0, x1, ...\]).
    pub positions: Vec<f64>,
    /// Previous particle positions (for XPBD velocity estimation).
    pub prev_positions: Vec<f64>,
    /// Per-particle inverse mass (0 = kinematic/fixed).
    pub inv_masses: Vec<f64>,
    /// Constraints.
    pub constraints: Vec<PbdConstraint>,
    /// Number of substeps per time step.
    pub substeps: u32,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

#[pymethods]
impl PyPbdSolver {
    /// Create a new PBD solver.
    #[new]
    pub fn new(positions: Vec<f64>, inv_masses: Vec<f64>) -> Self {
        let prev = positions.clone();
        Self {
            positions,
            prev_positions: prev,
            inv_masses,
            constraints: Vec::new(),
            substeps: 8,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add a stretch constraint between two particles.
    pub fn add_stretch(&mut self, i: usize, j: usize, compliance: f64) {
        let dx = self.positions[3 * i] - self.positions[3 * j];
        let dy = self.positions[3 * i + 1] - self.positions[3 * j + 1];
        let dz = self.positions[3 * i + 2] - self.positions[3 * j + 2];
        let rest = (dx * dx + dy * dy + dz * dz).sqrt();
        self.constraints.push(PbdConstraint {
            particles: vec![i, j],
            kind: PbdConstraintType::Stretch {
                rest_length: rest,
                compliance,
            },
            lambda: 0.0,
        });
    }

    /// Add a volume conservation constraint for the tetrahedron `(i, j, k, l)`.
    ///
    /// The rest volume is the signed volume of the tetrahedron in the
    /// *current* particle positions, so the constraint preserves the volume
    /// (and orientation) the tetrahedron has when it is added.
    ///
    /// Raises `IndexError` when an index is not a valid particle and
    /// `ValueError` when the same particle is listed twice.
    pub fn add_volume(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        compliance: f64,
    ) -> PyResult<()> {
        Ok(self.try_add_volume(i, j, k, l, compliance)?)
    }

    /// Perform one substep of XPBD simulation.
    pub fn substep(&mut self, dt: f64) {
        let h = dt / self.substeps as f64;
        let np = self.positions.len() / 3;
        // Apply gravity and predict positions
        for i in 0..np {
            if self.inv_masses[i] > 0.0 {
                let vx = self.positions[3 * i] - self.prev_positions[3 * i];
                let vy = self.positions[3 * i + 1] - self.prev_positions[3 * i + 1];
                let vz = self.positions[3 * i + 2] - self.prev_positions[3 * i + 2];
                self.prev_positions[3 * i] = self.positions[3 * i];
                self.prev_positions[3 * i + 1] = self.positions[3 * i + 1];
                self.prev_positions[3 * i + 2] = self.positions[3 * i + 2];
                self.positions[3 * i] += vx + self.gravity[0] * h * h;
                self.positions[3 * i + 1] += vy + self.gravity[1] * h * h;
                self.positions[3 * i + 2] += vz + self.gravity[2] * h * h;
            }
        }
        // Reset lambdas, then solve every constraint once (Gauss-Seidel).
        let Self {
            positions,
            inv_masses,
            constraints,
            ..
        } = self;
        for constraint in constraints.iter_mut() {
            constraint.lambda = 0.0;
            match &constraint.kind {
                PbdConstraintType::Stretch {
                    rest_length,
                    compliance,
                } => solve_stretch(
                    positions,
                    inv_masses,
                    &constraint.particles,
                    *rest_length,
                    *compliance,
                    h,
                    &mut constraint.lambda,
                ),
                PbdConstraintType::Volume {
                    rest_volume,
                    compliance,
                } => solve_volume(
                    positions,
                    inv_masses,
                    &constraint.particles,
                    *rest_volume,
                    *compliance,
                    h,
                    &mut constraint.lambda,
                ),
            }
        }
    }

    /// Step the simulation for one full time step (using substeps).
    pub fn step(&mut self, dt: f64) {
        for _ in 0..self.substeps {
            self.substep(dt);
        }
    }

    /// Return particle count.
    pub fn particle_count(&self) -> usize {
        self.positions.len() / 3
    }

    /// Get positions.
    #[getter]
    pub fn get_positions(&self) -> Vec<f64> {
        self.positions.clone()
    }

    /// Set positions.
    #[setter]
    pub fn set_positions(&mut self, v: Vec<f64>) {
        self.positions = v;
    }

    /// Get prev_positions.
    #[getter]
    pub fn get_prev_positions(&self) -> Vec<f64> {
        self.prev_positions.clone()
    }

    /// Get inv_masses.
    #[getter]
    pub fn get_inv_masses(&self) -> Vec<f64> {
        self.inv_masses.clone()
    }

    /// Set inv_masses.
    #[setter]
    pub fn set_inv_masses(&mut self, v: Vec<f64>) {
        self.inv_masses = v;
    }

    /// Get substeps.
    #[getter]
    pub fn get_substeps(&self) -> u32 {
        self.substeps
    }

    /// Set substeps.
    #[setter]
    pub fn set_substeps(&mut self, v: u32) {
        self.substeps = v;
    }

    /// Get gravity.
    #[getter]
    pub fn get_gravity(&self) -> [f64; 3] {
        self.gravity
    }

    /// Set gravity.
    #[setter]
    pub fn set_gravity(&mut self, v: [f64; 3]) {
        self.gravity = v;
    }
}

impl PyPbdSolver {
    /// Rust-level core of [`Self::add_volume`] with a typed error.
    pub fn try_add_volume(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        compliance: f64,
    ) -> Result<(), PbdConstraintError> {
        let ids = [i, j, k, l];
        let particle_count = (self.positions.len() / 3).min(self.inv_masses.len());
        let mut x = [[0.0f64; 3]; 4];
        for (slot, &id) in ids.iter().enumerate() {
            let position = particle_position(&self.positions, id)
                .filter(|_| id < particle_count)
                .ok_or(PbdConstraintError::ParticleOutOfRange {
                    index: id,
                    particle_count,
                })?;
            if ids[..slot].contains(&id) {
                return Err(PbdConstraintError::DuplicateParticle(id));
            }
            x[slot] = position;
        }
        self.constraints.push(PbdConstraint {
            particles: ids.to_vec(),
            kind: PbdConstraintType::Volume {
                rest_volume: signed_tet_volume(&x),
                compliance,
            },
            lambda: 0.0,
        });
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pbd_solver_stretch() {
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let inv_masses = vec![1.0, 1.0];
        let mut solver = PyPbdSolver::new(positions, inv_masses);
        solver.add_stretch(0, 1, 1e-4);
        solver.step(0.016);
        // Particles should have moved due to gravity
        assert!(solver.positions[1] < 0.0); // y should be negative
    }

    #[test]
    fn test_pbd_particle_count() {
        let pos = vec![0.0; 9]; // 3 particles
        let inv = vec![1.0; 3];
        let s = PyPbdSolver::new(pos, inv);
        assert_eq!(s.particle_count(), 3);
    }

    /// Unit right tetrahedron (volume 1/6), unit masses, no gravity.
    fn tetra_solver() -> PyPbdSolver {
        let positions = vec![
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ];
        let mut solver = PyPbdSolver::new(positions, vec![1.0; 4]);
        solver.gravity = [0.0; 3];
        solver
    }

    fn tet_corners(solver: &PyPbdSolver) -> [[f64; 3]; 4] {
        let mut x = [[0.0; 3]; 4];
        for (i, corner) in x.iter_mut().enumerate() {
            *corner = particle_position(&solver.positions, i).expect("tetrahedron particle");
        }
        x
    }

    fn tet_volume(solver: &PyPbdSolver) -> f64 {
        signed_tet_volume(&tet_corners(solver))
    }

    fn tet_centroid(solver: &PyPbdSolver) -> [f64; 3] {
        let x = tet_corners(solver);
        let mut c = [0.0; 3];
        for corner in &x {
            for k in 0..3 {
                c[k] += corner[k] / 4.0;
            }
        }
        c
    }

    /// Scale the tetrahedron about its centroid and bring it to rest.
    fn compress_about_centroid(solver: &mut PyPbdSolver, factor: f64) {
        let c = tet_centroid(solver);
        for particle in solver.positions.chunks_exact_mut(3) {
            for k in 0..3 {
                particle[k] = c[k] + factor * (particle[k] - c[k]);
            }
        }
        solver.prev_positions = solver.positions.clone();
    }

    fn rest_volume_of(solver: &PyPbdSolver, constraint: usize) -> f64 {
        match solver.constraints[constraint].kind {
            PbdConstraintType::Volume { rest_volume, .. } => rest_volume,
            PbdConstraintType::Stretch { .. } => panic!("expected a volume constraint"),
        }
    }

    #[test]
    fn test_pbd_add_volume_records_the_real_rest_volume() {
        let mut solver = tetra_solver();
        solver
            .add_volume(0, 1, 2, 3, 0.0)
            .expect("valid tetrahedron");
        assert!((rest_volume_of(&solver, 0) - 1.0 / 6.0).abs() < 1e-15);
        // Twice the edge length: eight times the volume, not a hard-coded 1.0.
        let mut big = tetra_solver();
        big.positions.iter_mut().for_each(|v| *v *= 2.0);
        big.add_volume(0, 1, 2, 3, 0.0).expect("valid tetrahedron");
        assert!((rest_volume_of(&big, 0) - 8.0 / 6.0).abs() < 1e-14);
        // Swapping two vertices flips the orientation, so the signed volume.
        let mut flipped = tetra_solver();
        flipped
            .add_volume(0, 2, 1, 3, 0.0)
            .expect("valid tetrahedron");
        assert!((rest_volume_of(&flipped, 0) + 1.0 / 6.0).abs() < 1e-15);
    }

    #[test]
    fn test_pbd_volume_constraint_restores_a_compressed_tetrahedron() {
        let mut solver = tetra_solver();
        solver
            .add_volume(0, 1, 2, 3, 0.0)
            .expect("valid tetrahedron");
        let rest = tet_volume(&solver);
        let centroid_before = tet_centroid(&solver);
        // Compress by 20 % of the volume (linear factor 0.8^(1/3)).
        compress_about_centroid(&mut solver, 0.8_f64.cbrt());
        let squeezed = tet_volume(&solver);
        assert!(
            ((rest - squeezed) / rest - 0.2).abs() < 1e-12,
            "setup must remove 20 % of the volume"
        );
        for _ in 0..20 {
            solver.step(1.0 / 60.0);
        }
        let recovered = tet_volume(&solver);
        assert!(
            ((recovered - rest) / rest).abs() < 0.01,
            "volume {recovered} should be back within 1 % of {rest}"
        );
        let centroid_after = tet_centroid(&solver);
        for k in 0..3 {
            assert!(
                (centroid_after[k] - centroid_before[k]).abs() < 1e-9,
                "internal constraint forces must not move the centroid"
            );
        }
    }

    #[test]
    fn test_pbd_volume_constraint_restores_an_inverted_tetrahedron() {
        let mut solver = tetra_solver();
        solver
            .add_volume(0, 2, 1, 3, 0.0)
            .expect("valid tetrahedron");
        let rest = rest_volume_of(&solver, 0);
        assert!(rest < 0.0);
        compress_about_centroid(&mut solver, 0.8_f64.cbrt());
        // Same particles, same orientation as the constraint: C = V - V0 > 0.
        for _ in 0..20 {
            solver.step(1.0 / 60.0);
        }
        // The constraint orders the particles (0, 2, 1, 3).
        let x = tet_corners(&solver);
        let recovered = signed_tet_volume(&[x[0], x[2], x[1], x[3]]);
        assert!(
            ((recovered - rest) / rest).abs() < 0.01,
            "signed volume {recovered} should be back within 1 % of {rest}"
        );
    }

    #[test]
    fn test_pbd_volume_constraint_respects_fixed_particles() {
        let mut solver = tetra_solver();
        solver.inv_masses[0] = 0.0;
        solver
            .add_volume(0, 1, 2, 3, 0.0)
            .expect("valid tetrahedron");
        let rest = tet_volume(&solver);
        compress_about_centroid(&mut solver, 0.8_f64.cbrt());
        let pinned = [
            solver.positions[0],
            solver.positions[1],
            solver.positions[2],
        ];
        for _ in 0..40 {
            solver.step(1.0 / 60.0);
        }
        assert_eq!(
            [
                solver.positions[0],
                solver.positions[1],
                solver.positions[2]
            ],
            pinned,
            "a particle with zero inverse mass must not move"
        );
        assert!(((tet_volume(&solver) - rest) / rest).abs() < 0.01);
    }

    #[test]
    fn test_pbd_soft_volume_constraint_recovers_less_than_a_stiff_one() {
        let run = |compliance: f64| {
            let mut solver = tetra_solver();
            solver
                .add_volume(0, 1, 2, 3, compliance)
                .expect("valid tetrahedron");
            let rest = tet_volume(&solver);
            compress_about_centroid(&mut solver, 0.8_f64.cbrt());
            solver.step(1.0 / 60.0);
            ((tet_volume(&solver) - rest) / rest).abs()
        };
        let stiff_error = run(0.0);
        let soft_error = run(1.0);
        assert!(stiff_error < 0.05, "stiff error {stiff_error}");
        assert!(
            soft_error > 0.15,
            "a very compliant constraint barely reacts, error {soft_error}"
        );
    }

    #[test]
    fn test_pbd_volume_constraint_holds_the_volume_against_gravity() {
        // The apex is pinned and the three base particles fall. Without a
        // constraint the volume (1 - d) / 6 shrinks with the fall distance d;
        // the volume constraint keeps it at its rest value.
        let run = |constrained: bool| {
            let mut solver = tetra_solver();
            solver.gravity = [0.0, -9.81, 0.0];
            solver.inv_masses[0] = 0.0;
            if constrained {
                solver
                    .add_volume(0, 1, 2, 3, 0.0)
                    .expect("valid tetrahedron");
            }
            for _ in 0..20 {
                solver.step(1.0 / 60.0);
            }
            tet_volume(&solver)
        };
        let rest = 1.0 / 6.0;
        let free = run(false);
        let held = run(true);
        assert!(
            ((free - rest) / rest).abs() > 0.3,
            "unconstrained volume {free} must change visibly"
        );
        assert!(
            ((held - rest) / rest).abs() < 0.05,
            "constrained volume {held} must stay near {rest}"
        );
    }

    #[test]
    fn test_pbd_add_volume_rejects_bad_particles() {
        let mut solver = tetra_solver();
        assert_eq!(
            solver.try_add_volume(0, 1, 2, 4, 0.0),
            Err(PbdConstraintError::ParticleOutOfRange {
                index: 4,
                particle_count: 4
            })
        );
        assert_eq!(
            solver.try_add_volume(0, 1, 1, 3, 0.0),
            Err(PbdConstraintError::DuplicateParticle(1))
        );
        assert!(solver.add_volume(0, 1, 2, usize::MAX, 0.0).is_err());
        assert!(solver.add_volume(3, 3, 3, 3, 0.0).is_err());
        assert!(
            solver.constraints.is_empty(),
            "rejected constraints must not be stored"
        );
    }

    #[test]
    fn test_pbd_volume_constraint_with_missing_particles_does_not_panic() {
        let mut solver = tetra_solver();
        solver
            .add_volume(0, 1, 2, 3, 0.0)
            .expect("valid tetrahedron");
        // Shrink the particle arrays behind the constraint's back.
        solver.positions.truncate(6);
        solver.prev_positions.truncate(6);
        solver.inv_masses.truncate(2);
        solver.step(1.0 / 60.0);
        assert_eq!(solver.particle_count(), 2);
    }
}
