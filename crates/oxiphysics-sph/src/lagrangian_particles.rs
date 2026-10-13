// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lagrangian particle tracking in SPH flow fields.
//!
//! Provides tracer particles with RK4 integration, Finite-Time Lyapunov
//! Exponents (FTLE), Lagrangian Coherent Structures (LCS) ridge detection,
//! deformation gradient computation, and mixing diagnostics.
//!
//! Every [`TracerParticle`] records its creation position (the Lagrangian
//! label `X`).  The diagnostics measure displacements against that recorded
//! reference configuration: the pair FTLE compares the current and initial
//! tracer separation, the deformation gradient is `F = ∂x/∂X` fitted to four
//! tracers, and the dispersion is the mean `|x − X|²`.

// ══════════════════════════════════════════════════════════════════════════════
// § 1  VECTOR HELPERS
// ══════════════════════════════════════════════════════════════════════════════

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

/// 3×3 identity matrix in row-major layout.
const MAT3_IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Row-major 3×3 matrix whose columns are the vectors `c0`, `c1`, `c2`.
#[inline]
fn mat3_from_columns(c0: [f64; 3], c1: [f64; 3], c2: [f64; 3]) -> [f64; 9] {
    [
        c0[0], c1[0], c2[0], //
        c0[1], c1[1], c2[1], //
        c0[2], c1[2], c2[2],
    ]
}

/// Determinant of a row-major 3×3 matrix.
#[inline]
fn mat3_det(m: &[f64; 9]) -> f64 {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}

/// Inverse of a row-major 3×3 matrix through its adjugate; `None` when the
/// determinant is not finite or vanishes.
fn mat3_inverse(m: &[f64; 9]) -> Option<[f64; 9]> {
    let det = mat3_det(m);
    if !det.is_finite() || det == 0.0 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ])
}

/// Product `a · b` of two row-major 3×3 matrices.
fn mat3_mul(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    let mut out = [0.0_f64; 9];
    for r in 0..3 {
        for c in 0..3 {
            out[r * 3 + c] = (0..3).map(|k| a[r * 3 + k] * b[k * 3 + c]).sum();
        }
    }
    out
}

/// Transpose of a row-major 3×3 matrix.
fn mat3_transpose(m: &[f64; 9]) -> [f64; 9] {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

/// Largest eigenvalue of a symmetric row-major 3×3 matrix (closed-form
/// trigonometric solution of the characteristic cubic).
fn sym3_max_eigenvalue(a: &[f64; 9]) -> f64 {
    let off_diagonal = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
    if off_diagonal == 0.0 {
        return a[0].max(a[4]).max(a[8]);
    }
    let mean = (a[0] + a[4] + a[8]) / 3.0;
    let spread_sq =
        (a[0] - mean).powi(2) + (a[4] - mean).powi(2) + (a[8] - mean).powi(2) + 2.0 * off_diagonal;
    let spread = (spread_sq / 6.0).sqrt();
    let shifted = [
        (a[0] - mean) / spread,
        a[1] / spread,
        a[2] / spread,
        a[3] / spread,
        (a[4] - mean) / spread,
        a[5] / spread,
        a[6] / spread,
        a[7] / spread,
        (a[8] - mean) / spread,
    ];
    let cos_arg = (mat3_det(&shifted) / 2.0).clamp(-1.0, 1.0);
    let angle = cos_arg.acos() / 3.0;
    mean + 2.0 * spread * angle.cos()
}

// ══════════════════════════════════════════════════════════════════════════════
// § 2  TRACER PARTICLE
// ══════════════════════════════════════════════════════════════════════════════

/// A massless tracer particle advected by a flow field.
#[derive(Debug, Clone)]
pub struct TracerParticle {
    /// Current position \[m\].
    pub position: [f64; 3],
    /// Current velocity \[m/s\] (updated after each integration step).
    pub velocity: [f64; 3],
    /// Unique particle identifier.
    pub id: usize,
    /// Accumulated time the particle has been tracked \[s\].
    pub age: f64,
    /// Position at creation \[m\] (the Lagrangian label `X` of the particle).
    ///
    /// Set by [`TracerParticle::new`] and never modified by
    /// [`TracerParticle::integrate_rk4`]; FTLE, deformation gradient and
    /// dispersion diagnostics use it as the reference configuration.
    pub initial_position: [f64; 3],
}

impl TracerParticle {
    /// Create a new tracer at `position` with the given `id`.
    ///
    /// `position` is also recorded as the tracer's [`initial_position`](Self::initial_position).
    pub fn new(position: [f64; 3], id: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            id,
            age: 0.0,
            initial_position: position,
        }
    }

    /// Advance the tracer by one RK4 step of size `dt` seconds at current
    /// time `t` using the provided velocity field closure.
    ///
    /// The velocity field is `vel_field(pos, t) -> velocity`.
    pub fn integrate_rk4(
        &mut self,
        vel_field: &dyn Fn([f64; 3], f64) -> [f64; 3],
        dt: f64,
        t: f64,
    ) {
        let p = self.position;

        let k1 = vel_field(p, t);
        let k2 = vel_field(vec3_add(p, vec3_scale(k1, 0.5 * dt)), t + 0.5 * dt);
        let k3 = vel_field(vec3_add(p, vec3_scale(k2, 0.5 * dt)), t + 0.5 * dt);
        let k4 = vel_field(vec3_add(p, vec3_scale(k3, dt)), t + dt);

        // Weighted average
        let dp = vec3_scale(
            vec3_add(
                vec3_add(k1, vec3_scale(k2, 2.0)),
                vec3_add(vec3_scale(k3, 2.0), k4),
            ),
            dt / 6.0,
        );

        self.position = vec3_add(p, dp);
        self.velocity = vel_field(self.position, t + dt);
        self.age += dt;
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 3  PARTICLE TRACKER
// ══════════════════════════════════════════════════════════════════════════════

/// Tracks a collection of `TracerParticle`s through a time-varying flow field.
#[derive(Debug, Clone)]
pub struct ParticleTracker {
    /// All tracer particles managed by this tracker.
    pub tracers: Vec<TracerParticle>,
}

impl ParticleTracker {
    /// Create a tracker from a list of initial particle positions.
    ///
    /// Each particle receives an id equal to its index in the list.
    pub fn new(positions: Vec<[f64; 3]>) -> Self {
        let tracers = positions
            .into_iter()
            .enumerate()
            .map(|(id, pos)| TracerParticle::new(pos, id))
            .collect();
        Self { tracers }
    }

    /// Advance all tracers by one RK4 step at current time `t` with step `dt`.
    pub fn step(&mut self, vel_field: &dyn Fn([f64; 3], f64) -> [f64; 3], dt: f64, t: f64) {
        for p in &mut self.tracers {
            p.integrate_rk4(vel_field, dt, t);
        }
    }

    /// Number of tracked particles.
    pub fn len(&self) -> usize {
        self.tracers.len()
    }

    /// Returns `true` if there are no particles.
    pub fn is_empty(&self) -> bool {
        self.tracers.is_empty()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 4  FINITE-TIME LYAPUNOV EXPONENT
// ══════════════════════════════════════════════════════════════════════════════

/// Computes the Finite-Time Lyapunov Exponent (FTLE) of tracer sets.
///
/// FTLE measures the average exponential rate of separation of nearby
/// trajectories over a finite advection window `T`.  Two estimators are
/// provided:
///
/// * [`ftle`](Self::ftle): the scalar rate of one tracer *pair*; it depends on
///   the orientation of the pair relative to the stretching direction.
/// * [`ftle_tetrad`](Self::ftle_tetrad): the orientation-independent exponent
///   `σ = ln √λ_max(FᵀF) / |T|` of the Cauchy–Green strain tensor of a
///   four-tracer tetrahedron.
///
/// Both compare the tracers' current positions with their recorded
/// [`initial_position`](TracerParticle::initial_position), so the tracers must
/// have been advected for the window `T` since they were created.
#[derive(Debug, Clone)]
pub struct FiniteTimeLyapunov {
    /// Integration window T \[s\] over which to measure divergence.
    pub advection_time: f64,
}

impl FiniteTimeLyapunov {
    /// Create with the given advection time T.
    pub fn new(advection_time: f64) -> Self {
        Self { advection_time }
    }

    /// Length `|T|` of the advection window, or `None` when it is zero or not
    /// finite (no rate can be defined).
    fn window(&self) -> Option<f64> {
        let window = self.advection_time.abs();
        (window.is_finite() && window >= f64::EPSILON).then_some(window)
    }

    /// Compute the FTLE of the first two tracers of `tracer_pair`.
    ///
    /// Uses the scalar formula:
    /// σ = (1/|T|) ln(|x₂(T) − x₁(T)| / |X₂ − X₁|)
    ///
    /// where `X` are the recorded initial positions.  A pair that is pulled
    /// apart gives σ > 0, a pair that is squeezed together σ < 0.
    ///
    /// Returns 0.0 when fewer than two tracers are supplied, when the advection
    /// time is zero or not finite, or when the initial separation is below
    /// machine epsilon.  A pair that has collapsed onto a single point has the
    /// exponent `-∞`.
    pub fn ftle(&self, tracer_pair: &[TracerParticle]) -> f64 {
        let (Some(p1), Some(p2)) = (tracer_pair.first(), tracer_pair.get(1)) else {
            return 0.0;
        };
        let Some(window) = self.window() else {
            return 0.0;
        };
        let sep_init = vec3_norm(vec3_sub(p2.initial_position, p1.initial_position));
        if sep_init < f64::EPSILON {
            return 0.0;
        }
        let sep_now = vec3_norm(vec3_sub(p2.position, p1.position));
        (sep_now / sep_init).ln() / window
    }

    /// FTLE `σ = ln √λ_max(FᵀF) / |T|` of a deformation gradient `f`
    /// (row-major 3×3), i.e. the logarithm of its largest principal stretch per
    /// unit time.
    ///
    /// Returns 0.0 when the advection time is zero or not finite.
    pub fn ftle_from_deformation_gradient(&self, f: &[f64; 9]) -> f64 {
        match self.window() {
            Some(window) => max_principal_stretch(f).ln() / window,
            None => 0.0,
        }
    }

    /// Orientation-independent FTLE of the tetrahedron spanned by the first
    /// four tracers of `tracers`.
    ///
    /// The deformation gradient is fitted with
    /// [`deformation_gradient_from_tetrad`] and reduced with
    /// [`ftle_from_deformation_gradient`](Self::ftle_from_deformation_gradient).
    ///
    /// Returns 0.0 when fewer than four tracers are supplied or when their
    /// initial configuration is degenerate (coplanar or coincident).
    pub fn ftle_tetrad(&self, tracers: &[TracerParticle]) -> f64 {
        match deformation_gradient_from_tetrad(tracers) {
            Some(f) => self.ftle_from_deformation_gradient(&f),
            None => 0.0,
        }
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 5  LCS DETECTOR
// ══════════════════════════════════════════════════════════════════════════════

/// Detects Lagrangian Coherent Structures via ridges of the FTLE field.
#[derive(Debug, Clone)]
pub struct LcsDetector {
    /// Flattened FTLE scalar field (row-major, x fastest).
    pub ftle_field: Vec<f64>,
    /// Grid dimensions `[nx, ny, nz]`.
    pub grid_size: [usize; 3],
}

impl LcsDetector {
    /// Create a detector with the given FTLE field and grid size.
    ///
    /// The field must have length `nx * ny * nz`.
    pub fn new(ftle_field: Vec<f64>, grid_size: [usize; 3]) -> Self {
        debug_assert_eq!(
            ftle_field.len(),
            grid_size[0] * grid_size[1] * grid_size[2],
            "FTLE field size mismatch"
        );
        Self {
            ftle_field,
            grid_size,
        }
    }

    /// Detect ridge points in the FTLE field.
    ///
    /// A grid point is classified as a ridge if its FTLE value exceeds the
    /// local neighbourhood average by more than `threshold` standard
    /// deviations.
    ///
    /// Returns a list of `[x, y, z]` grid-index triplets (cast to f64) of
    /// all ridge points found.
    pub fn detect_ridges(&self) -> Vec<[f64; 3]> {
        let [nx, ny, nz] = self.grid_size;
        if nx * ny * nz == 0 {
            return Vec::new();
        }

        // Compute global mean and std-dev
        let n = self.ftle_field.len() as f64;
        let mean = self.ftle_field.iter().sum::<f64>() / n;
        let var = self
            .ftle_field
            .iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>()
            / n;
        let std = var.sqrt();
        let threshold = mean + std; // one standard deviation above mean

        let mut ridges = Vec::new();
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let idx = iz * ny * nx + iy * nx + ix;
                    if self.ftle_field[idx] > threshold {
                        ridges.push([ix as f64, iy as f64, iz as f64]);
                    }
                }
            }
        }
        ridges
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 6  MIXING DIAGNOSTICS
// ══════════════════════════════════════════════════════════════════════════════

/// Diagnostics for Lagrangian mixing and stretching.
#[derive(Debug, Clone)]
pub struct MixingDiagnostic {
    /// Characteristic stretch rate λ \[1/s\] of the flow, a reference value
    /// supplied by the caller (see [`MixingDiagnostic::expected_stretch`]).
    ///
    /// It does not enter [`MixingDiagnostic::compute_deformation_gradient`],
    /// which is measured from the tracer positions alone.
    pub stretch_rate: f64,
}

impl MixingDiagnostic {
    /// Create with the given stretch rate.
    pub fn new(stretch_rate: f64) -> Self {
        Self { stretch_rate }
    }

    /// Deformation gradient `F = ∂x/∂X` of the material tetrahedron spanned by
    /// the first four tracers, as a 3×3 matrix stored row-major in a
    /// `[f64; 9]` array.
    ///
    /// With `Δx_k = x_k − x_0` (current positions) and `ΔX_k = X_k − X_0`
    /// (recorded [`initial_position`](TracerParticle::initial_position)s) for
    /// `k = 1, 2, 3`, the linear map that carries the reference edges onto the
    /// current ones is `F = [Δx₁ Δx₂ Δx₃] · [ΔX₁ ΔX₂ ΔX₃]⁻¹`.
    ///
    /// Returns the identity when fewer than four tracers are supplied or when
    /// the reference tetrahedron is degenerate (see
    /// [`deformation_gradient_from_tetrad`]).
    pub fn compute_deformation_gradient(&self, tracers: &[TracerParticle]) -> [f64; 9] {
        deformation_gradient_from_tetrad(tracers).unwrap_or(MAT3_IDENTITY)
    }

    /// Stretch factor `exp(λ t)` that a material line element reaches after
    /// `elapsed` seconds when it is stretched exponentially at the
    /// characteristic [`stretch_rate`](Self::stretch_rate) `λ`.
    ///
    /// Compare with [`max_principal_stretch`] of the measured deformation
    /// gradient to judge how close the flow is to the reference rate.
    pub fn expected_stretch(&self, elapsed: f64) -> f64 {
        (self.stretch_rate * elapsed).exp()
    }
}

// ══════════════════════════════════════════════════════════════════════════════
// § 7  FREE FUNCTIONS
// ══════════════════════════════════════════════════════════════════════════════

/// Compute the Mean Square Displacement (MSD) of a set of tracers relative to
/// their recorded initial positions:
///
/// `MSD = (1/N) Σ_i |x_i − X_i|²`
///
/// Returns the ensemble-averaged squared displacement \[m²\], or 0.0 for an
/// empty set.
pub fn particle_dispersion(tracers: &[TracerParticle]) -> f64 {
    if tracers.is_empty() {
        return 0.0;
    }
    let total: f64 = tracers
        .iter()
        .map(|p| {
            let displacement = vec3_sub(p.position, p.initial_position);
            vec3_dot(displacement, displacement)
        })
        .sum();
    total / tracers.len() as f64
}

/// Fit the deformation gradient `F = ∂x/∂X` to the first four tracers.
///
/// The tracers span a material tetrahedron with edges `ΔX_k = X_k − X_0`
/// (recorded initial positions) that has been carried to the edges
/// `Δx_k = x_k − x_0` (current positions), `k = 1, 2, 3`.  The unique linear
/// map between them is `F = [Δx₁ Δx₂ Δx₃] · [ΔX₁ ΔX₂ ΔX₃]⁻¹`, returned
/// row-major.  For a flow that is linear over the tetrahedron it equals the
/// flow-map gradient `exp(L T)` exactly.
///
/// Returns `None` when fewer than four tracers are supplied, when the reference
/// tetrahedron is degenerate (its volume is below `1e-12` of the volume of the
/// box spanned by its edge lengths, i.e. coplanar or coincident tracers), or
/// when the fitted matrix is not finite.
pub fn deformation_gradient_from_tetrad(tracers: &[TracerParticle]) -> Option<[f64; 9]> {
    let [t0, t1, t2, t3] = tracers.get(..4)? else {
        return None;
    };
    let reference_edges = [
        vec3_sub(t1.initial_position, t0.initial_position),
        vec3_sub(t2.initial_position, t0.initial_position),
        vec3_sub(t3.initial_position, t0.initial_position),
    ];
    let current_edges = [
        vec3_sub(t1.position, t0.position),
        vec3_sub(t2.position, t0.position),
        vec3_sub(t3.position, t0.position),
    ];
    let reference = mat3_from_columns(reference_edges[0], reference_edges[1], reference_edges[2]);
    let current = mat3_from_columns(current_edges[0], current_edges[1], current_edges[2]);

    // |det| can never exceed the product of the edge lengths (Hadamard), so the
    // ratio is a unit-free measure of how flat the reference tetrahedron is.
    let box_volume = vec3_norm(reference_edges[0])
        * vec3_norm(reference_edges[1])
        * vec3_norm(reference_edges[2]);
    let det = mat3_det(&reference);
    if !det.is_finite() || det.abs() <= 1e-12 * box_volume {
        return None;
    }
    let reference_inverse = mat3_inverse(&reference)?;
    let f = mat3_mul(&current, &reference_inverse);
    f.iter().all(|v| v.is_finite()).then_some(f)
}

/// Largest principal stretch `√λ_max(FᵀF)` (the largest singular value) of a
/// deformation gradient `f` stored row-major.
///
/// It equals 1 for the identity and for any pure rotation, and `e^{λT}` for a
/// hyperbolic flow with stretching rate `λ` observed over a time `T`.
pub fn max_principal_stretch(f: &[f64; 9]) -> f64 {
    let cauchy_green = mat3_mul(&mat3_transpose(f), f);
    sym3_max_eigenvalue(&cauchy_green).max(0.0).sqrt()
}

/// Compute the relative dispersion (separation distance) between two tracers.
///
/// Returns `|p2.position − p1.position|` \[m\].
pub fn relative_dispersion(p1: &TracerParticle, p2: &TracerParticle) -> f64 {
    vec3_norm(vec3_sub(p2.position, p1.position))
}

// ══════════════════════════════════════════════════════════════════════════════
// § 8  TESTS
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── RK4 accuracy ─────────────────────────────────────────────────────────

    /// Uniform flow: u = (1, 0, 0). Exact solution: x(t) = x0 + t.
    fn uniform_flow(pos: [f64; 3], _t: f64) -> [f64; 3] {
        let _ = pos;
        [1.0, 0.0, 0.0]
    }

    #[test]
    fn test_rk4_uniform_flow_exact() {
        let mut p = TracerParticle::new([0.0, 0.0, 0.0], 0);
        let dt = 0.1;
        let n = 10;
        for i in 0..n {
            p.integrate_rk4(&uniform_flow, dt, i as f64 * dt);
        }
        // After 10 steps of dt=0.1 the particle should be at x=1.0
        assert!((p.position[0] - 1.0).abs() < 1e-12);
        assert!(p.position[1].abs() < 1e-14);
        assert!(p.position[2].abs() < 1e-14);
    }

    #[test]
    fn test_rk4_age_accumulation() {
        let mut p = TracerParticle::new([0.0, 0.0, 0.0], 1);
        let dt = 0.25;
        for i in 0..4 {
            p.integrate_rk4(&uniform_flow, dt, i as f64 * dt);
        }
        assert!((p.age - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rk4_circular_flow() {
        // Solid-body rotation: u = (-y, x, 0), period 2π
        let flow = |pos: [f64; 3], _t: f64| -> [f64; 3] { [-pos[1], pos[0], 0.0] };
        let mut p = TracerParticle::new([1.0, 0.0, 0.0], 0);
        let n = 1000;
        let dt = 2.0 * PI / n as f64;
        for i in 0..n {
            p.integrate_rk4(&flow, dt, i as f64 * dt);
        }
        // After one full revolution, particle should be back near [1, 0, 0]
        assert!((p.position[0] - 1.0).abs() < 1e-6, "x={}", p.position[0]);
        assert!(p.position[1].abs() < 1e-6, "y={}", p.position[1]);
    }

    #[test]
    fn test_rk4_velocity_updated() {
        let flow = |_pos: [f64; 3], _t: f64| -> [f64; 3] { [2.0, 3.0, 0.0] };
        let mut p = TracerParticle::new([0.0, 0.0, 0.0], 0);
        p.integrate_rk4(&flow, 0.1, 0.0);
        assert!((p.velocity[0] - 2.0).abs() < 1e-12);
        assert!((p.velocity[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_rk4_stationary_field_no_move() {
        let zero_flow = |_pos: [f64; 3], _t: f64| -> [f64; 3] { [0.0, 0.0, 0.0] };
        let mut p = TracerParticle::new([1.0, 2.0, 3.0], 0);
        for i in 0..100 {
            p.integrate_rk4(&zero_flow, 0.01, i as f64 * 0.01);
        }
        assert!((p.position[0] - 1.0).abs() < 1e-14);
        assert!((p.position[1] - 2.0).abs() < 1e-14);
        assert!((p.position[2] - 3.0).abs() < 1e-14);
    }

    // ── Particle tracker ──────────────────────────────────────────────────────

    #[test]
    fn test_tracker_length() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tracker = ParticleTracker::new(positions);
        assert_eq!(tracker.len(), 3);
        assert!(!tracker.is_empty());
    }

    #[test]
    fn test_tracker_empty() {
        let tracker = ParticleTracker::new(Vec::new());
        assert!(tracker.is_empty());
        assert_eq!(tracker.len(), 0);
    }

    #[test]
    fn test_tracker_step_moves_all() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
        let mut tracker = ParticleTracker::new(positions);
        tracker.step(&uniform_flow, 0.1, 0.0);
        for p in &tracker.tracers {
            assert!((p.position[0] - 0.1).abs() < 1e-12 || (p.position[0] - 1.1).abs() < 1e-12);
        }
    }

    #[test]
    fn test_tracker_ids_preserved() {
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tracker = ParticleTracker::new(positions);
        for (i, p) in tracker.tracers.iter().enumerate() {
            assert_eq!(p.id, i);
        }
    }

    // ── FTLE ──────────────────────────────────────────────────────────────────

    #[test]
    fn test_ftle_zero_advection_time() {
        let ftl = FiniteTimeLyapunov::new(0.0);
        let p1 = TracerParticle::new([0.0, 0.0, 0.0], 0);
        let p2 = TracerParticle::new([0.1, 0.0, 0.0], 1);
        let result = ftl.ftle(&[p1, p2]);
        assert_eq!(result, 0.0);
    }

    /// Saddle flow `u = (λx, −λy, 0)` with stretching rate λ = 0.5.
    fn saddle_flow(pos: [f64; 3], _t: f64) -> [f64; 3] {
        [0.5 * pos[0], -0.5 * pos[1], 0.0]
    }

    /// Simple shear `u = (γy, 0, 0)` with γ = 1.
    fn shear_flow(pos: [f64; 3], _t: f64) -> [f64; 3] {
        [pos[1], 0.0, 0.0]
    }

    /// Solid-body rotation `u = (−ωy, ωx, 0)` with ω = 1.
    fn rotation_flow(pos: [f64; 3], _t: f64) -> [f64; 3] {
        [-pos[1], pos[0], 0.0]
    }

    /// Advect every tracer for `t_end` seconds with RK4 steps of size `dt`
    /// (both negative for a backward integration).
    fn advect(
        tracers: &mut [TracerParticle],
        flow: &dyn Fn([f64; 3], f64) -> [f64; 3],
        t_end: f64,
        dt: f64,
    ) {
        let steps = (t_end / dt).round() as usize;
        for step in 0..steps {
            let t = step as f64 * dt;
            for p in tracers.iter_mut() {
                p.integrate_rk4(flow, dt, t);
            }
        }
    }

    /// Tracer at `origin` plus three tracers displaced by `edges`.
    fn tetrad(origin: [f64; 3], edges: [[f64; 3]; 3]) -> Vec<TracerParticle> {
        let mut tracers = vec![TracerParticle::new(origin, 0)];
        for (k, edge) in edges.iter().enumerate() {
            tracers.push(TracerParticle::new(vec3_add(origin, *edge), k + 1));
        }
        tracers
    }

    /// Regular but skewed tetrad away from the origin.
    fn skewed_tetrad() -> Vec<TracerParticle> {
        tetrad(
            [1.0, -0.5, 0.25],
            [[0.1, 0.02, 0.0], [-0.03, 0.12, 0.01], [0.01, 0.0, 0.09]],
        )
    }

    fn assert_matrix_close(actual: &[f64; 9], expected: &[f64; 9], tol: f64) {
        for (k, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
            assert!(
                (a - e).abs() < tol,
                "entry {k}: got {a}, expected {e} (full matrix {actual:?})"
            );
        }
    }

    #[test]
    fn test_ftle_saddle_flow_stretching_axis() {
        // Pair separated along x (the stretching direction): the separation
        // grows by e^{λT} = e, so σ = λ = 0.5.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut pair = [
            TracerParticle::new([1.0, 0.5, 0.0], 0),
            TracerParticle::new([1.001, 0.5, 0.0], 1),
        ];
        advect(&mut pair, &saddle_flow, 2.0, 0.01);
        let sigma = ftl.ftle(&pair);
        assert!((sigma - 0.5).abs() < 1e-3, "FTLE = {sigma}, expected 0.5");
    }

    #[test]
    fn test_ftle_saddle_flow_contracting_axis() {
        // Pair separated along y (the contracting direction): σ = −λ.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut pair = [
            TracerParticle::new([0.5, 1.0, 0.0], 0),
            TracerParticle::new([0.5, 1.001, 0.0], 1),
        ];
        advect(&mut pair, &saddle_flow, 2.0, 0.01);
        let sigma = ftl.ftle(&pair);
        assert!((sigma + 0.5).abs() < 1e-3, "FTLE = {sigma}, expected -0.5");
    }

    #[test]
    fn test_ftle_uniform_translation_is_zero() {
        let ftl = FiniteTimeLyapunov::new(1.0);
        let mut pair = [
            TracerParticle::new([0.0, 0.0, 0.0], 0),
            TracerParticle::new([0.01, 0.0, 0.0], 1),
        ];
        advect(&mut pair, &uniform_flow, 1.0, 0.01);
        let sigma = ftl.ftle(&pair);
        assert!(sigma.abs() < 1e-12, "translation must not stretch: {sigma}");
    }

    #[test]
    fn test_ftle_matches_closed_form_for_prescribed_positions() {
        // Initial separation 0.1, final separation 0.5, T = 2.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let p1 = TracerParticle::new([0.0, 0.0, 0.0], 0);
        let mut p2 = TracerParticle::new([0.1, 0.0, 0.0], 1);
        p2.position = [0.0, 0.3, 0.4];
        let sigma = ftl.ftle(&[p1, p2]);
        let expected = 5.0_f64.ln() / 2.0;
        assert!(
            (sigma - expected).abs() < 1e-14,
            "FTLE = {sigma}, expected {expected}"
        );
    }

    #[test]
    fn test_ftle_does_not_depend_on_tracer_age() {
        // The previous implementation inferred the reference separation from
        // the window and returned 1.0 for any flow.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let p1 = TracerParticle::new([0.0, 0.0, 0.0], 0);
        let mut p2 = TracerParticle::new([0.1, 0.0, 0.0], 1);
        p2.position = [0.1, 0.0, 0.0];
        p2.age = 123.0;
        assert_eq!(ftl.ftle(&[p1, p2]), 0.0);
    }

    #[test]
    fn test_ftle_backward_window_uses_absolute_time() {
        // Going backwards in time along the stretching axis the separation
        // shrinks by e^{-1}: σ = ln(e^{-1}) / |T| = −0.5.
        let ftl = FiniteTimeLyapunov::new(-2.0);
        let mut pair = [
            TracerParticle::new([1.0, 0.5, 0.0], 0),
            TracerParticle::new([1.001, 0.5, 0.0], 1),
        ];
        advect(&mut pair, &saddle_flow, -2.0, -0.01);
        let sigma = ftl.ftle(&pair);
        assert!((sigma + 0.5).abs() < 1e-3, "FTLE = {sigma}, expected -0.5");
    }

    #[test]
    fn test_ftle_degenerate_inputs_return_zero_without_panic() {
        let ftl = FiniteTimeLyapunov::new(1.0);
        assert_eq!(ftl.ftle(&[]), 0.0);
        assert_eq!(ftl.ftle(&[TracerParticle::new([0.0; 3], 0)]), 0.0);

        // Coincident initial positions: no finite reference separation.
        let a = TracerParticle::new([1.0, 1.0, 1.0], 0);
        let mut b = TracerParticle::new([1.0, 1.0, 1.0], 1);
        b.position = [2.0, 1.0, 1.0];
        assert_eq!(ftl.ftle(&[a.clone(), b.clone()]), 0.0);

        // Non-finite advection windows define no rate.
        let mut c = TracerParticle::new([1.0, 1.0, 1.0], 2);
        c.position = [3.0, 1.0, 1.0];
        let d = TracerParticle::new([1.1, 1.0, 1.0], 3);
        assert_eq!(
            FiniteTimeLyapunov::new(f64::NAN).ftle(&[d.clone(), c.clone()]),
            0.0
        );
        assert_eq!(FiniteTimeLyapunov::new(f64::INFINITY).ftle(&[d, c]), 0.0);
    }

    #[test]
    fn test_ftle_collapsed_pair_is_negative_infinity() {
        let ftl = FiniteTimeLyapunov::new(1.0);
        let p1 = TracerParticle::new([0.0; 3], 0);
        let mut p2 = TracerParticle::new([0.1, 0.0, 0.0], 1);
        p2.position = [0.0; 3];
        assert_eq!(ftl.ftle(&[p1, p2]), f64::NEG_INFINITY);
    }

    #[test]
    fn test_integrate_rk4_preserves_initial_position() {
        let mut p = TracerParticle::new([0.3, -0.2, 0.7], 5);
        for i in 0..50 {
            p.integrate_rk4(&saddle_flow, 0.02, i as f64 * 0.02);
        }
        assert_eq!(p.initial_position, [0.3, -0.2, 0.7]);
        assert!(p.position[0] > 0.3, "tracer must have moved along +x");
    }

    #[test]
    fn test_tracker_records_initial_positions() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let mut tracker = ParticleTracker::new(positions.clone());
        tracker.step(&uniform_flow, 0.5, 0.0);
        for (p, x0) in tracker.tracers.iter().zip(positions.iter()) {
            assert_eq!(&p.initial_position, x0);
            assert!((p.position[0] - x0[0] - 0.5).abs() < 1e-14);
        }
    }

    // ── FTLE from the Cauchy–Green tensor of a tracer tetrad ───────────────────

    #[test]
    fn test_ftle_tetrad_saddle_flow() {
        // F = diag(e^{λT}, e^{-λT}, 1): the largest stretch is e^{λT}, σ = λ.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &saddle_flow, 2.0, 0.01);
        let sigma = ftl.ftle_tetrad(&tracers);
        assert!((sigma - 0.5).abs() < 1e-8, "FTLE = {sigma}, expected 0.5");
    }

    #[test]
    fn test_ftle_tetrad_shear_flow() {
        // F = [[1, γT, 0], [0, 1, 0], [0, 0, 1]] with γT = 2 has the largest
        // singular value 1 + √2, so σ = ln(1 + √2) / T.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &shear_flow, 2.0, 0.01);
        let sigma = ftl.ftle_tetrad(&tracers);
        let expected = (1.0 + 2.0_f64.sqrt()).ln() / 2.0;
        assert!(
            (sigma - expected).abs() < 1e-8,
            "FTLE = {sigma}, expected {expected}"
        );
    }

    #[test]
    fn test_ftle_tetrad_rigid_rotation_is_zero() {
        // Rotation changes no distances, so no exponential stretching occurs.
        let ftl = FiniteTimeLyapunov::new(PI);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &rotation_flow, PI, PI / 400.0);
        let sigma = ftl.ftle_tetrad(&tracers);
        assert!(sigma.abs() < 1e-8, "rotation FTLE = {sigma}");
    }

    #[test]
    fn test_ftle_tetrad_is_orientation_independent() {
        // The pair FTLE of a pair oriented at 45° to the stretching axis is
        // smaller than the true exponent, the tetrad value is not.
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut tracers = tetrad(
            [0.0, 0.0, 0.0],
            [[0.05, 0.05, 0.0], [-0.05, 0.05, 0.0], [0.0, 0.0, 0.1]],
        );
        advect(&mut tracers, &saddle_flow, 2.0, 0.01);
        let pair = ftl.ftle(&tracers[..2]);
        let tetrad_value = ftl.ftle_tetrad(&tracers);
        assert!((tetrad_value - 0.5).abs() < 1e-8, "tetrad {tetrad_value}");
        assert!(pair < 0.45, "oblique pair FTLE {pair} should underestimate");
        assert!(pair > 0.0);
    }

    #[test]
    fn test_ftle_tetrad_degenerate_inputs_return_zero() {
        let ftl = FiniteTimeLyapunov::new(2.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &saddle_flow, 2.0, 0.01);
        assert_eq!(ftl.ftle_tetrad(&tracers[..3]), 0.0);
        assert_eq!(ftl.ftle_tetrad(&[]), 0.0);
        // Zero advection window.
        assert_eq!(FiniteTimeLyapunov::new(0.0).ftle_tetrad(&tracers), 0.0);
        // Coplanar reference configuration.
        let flat = tetrad(
            [0.0; 3],
            [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.1, 0.1, 0.0]],
        );
        assert_eq!(ftl.ftle_tetrad(&flat), 0.0);
    }

    #[test]
    fn test_ftle_from_deformation_gradient_closed_form() {
        let ftl = FiniteTimeLyapunov::new(4.0);
        let f = [3.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 1.0];
        let sigma = ftl.ftle_from_deformation_gradient(&f);
        assert!((sigma - 3.0_f64.ln() / 4.0).abs() < 1e-14);
        assert_eq!(
            FiniteTimeLyapunov::new(0.0).ftle_from_deformation_gradient(&f),
            0.0
        );
    }

    // ── LCS detector ──────────────────────────────────────────────────────────

    #[test]
    fn test_lcs_empty_grid() {
        let det = LcsDetector::new(Vec::new(), [0, 0, 0]);
        assert!(det.detect_ridges().is_empty());
    }

    #[test]
    fn test_lcs_uniform_field_no_ridges() {
        // Uniform FTLE → no point exceeds one std-dev above mean (std=0)
        let n = 27usize;
        let field = vec![1.0; n];
        let det = LcsDetector::new(field, [3, 3, 3]);
        let ridges = det.detect_ridges();
        assert!(ridges.is_empty());
    }

    #[test]
    fn test_lcs_single_spike() {
        let mut field = vec![0.0_f64; 8];
        field[4] = 10.0; // one large value
        let det = LcsDetector::new(field, [2, 2, 2]);
        let ridges = det.detect_ridges();
        assert!(!ridges.is_empty());
    }

    #[test]
    fn test_lcs_grid_size_consistency() {
        let field: Vec<f64> = (0..64).map(|i| i as f64).collect();
        let det = LcsDetector::new(field, [4, 4, 4]);
        let ridges = det.detect_ridges();
        // All ridge indices should be within grid bounds
        for r in &ridges {
            assert!(r[0] < 4.0 && r[1] < 4.0 && r[2] < 4.0);
        }
    }

    #[test]
    fn test_lcs_2d_slice() {
        let field: Vec<f64> = (0..100).map(|i| (i as f64).sin().abs()).collect();
        let det = LcsDetector::new(field, [10, 10, 1]);
        let _ridges = det.detect_ridges();
        // Just check it runs without panic
    }

    // ── Mixing diagnostics ────────────────────────────────────────────────────

    #[test]
    fn test_deformation_gradient_identity_for_undeformed_tetrad() {
        let diag = MixingDiagnostic::new(3.0);
        let f = diag.compute_deformation_gradient(&skewed_tetrad());
        assert_matrix_close(&f, &MAT3_IDENTITY, 1e-14);
    }

    #[test]
    fn test_deformation_gradient_saddle_flow() {
        let diag = MixingDiagnostic::new(0.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &saddle_flow, 2.0, 0.01);
        let f = diag.compute_deformation_gradient(&tracers);
        let e = 1.0_f64.exp();
        assert_matrix_close(&f, &[e, 0.0, 0.0, 0.0, 1.0 / e, 0.0, 0.0, 0.0, 1.0], 1e-8);
    }

    #[test]
    fn test_deformation_gradient_shear_flow() {
        let diag = MixingDiagnostic::new(0.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &shear_flow, 2.0, 0.01);
        let f = diag.compute_deformation_gradient(&tracers);
        assert_matrix_close(&f, &[1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 1e-9);
    }

    #[test]
    fn test_deformation_gradient_rotation_flow() {
        // Quarter turn: F = R(π/2), a proper rotation (det = 1).
        let diag = MixingDiagnostic::new(0.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &rotation_flow, 0.5 * PI, PI / 800.0);
        let f = diag.compute_deformation_gradient(&tracers);
        assert_matrix_close(&f, &[0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], 1e-8);
        assert!((mat3_det(&f) - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_deformation_gradient_independent_of_stretch_rate() {
        // The deformation gradient is measured from positions; the reference
        // stretch rate must not leak into it.
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &shear_flow, 1.0, 0.01);
        let f_zero = MixingDiagnostic::new(0.0).compute_deformation_gradient(&tracers);
        let f_big = MixingDiagnostic::new(25.0).compute_deformation_gradient(&tracers);
        assert_eq!(f_zero, f_big);
    }

    #[test]
    fn test_deformation_gradient_needs_four_tracers() {
        let diag = MixingDiagnostic::new(1.0);
        let mut tracers = skewed_tetrad();
        advect(&mut tracers, &saddle_flow, 1.0, 0.01);
        for count in 0..4 {
            assert_eq!(
                diag.compute_deformation_gradient(&tracers[..count]),
                MAT3_IDENTITY,
                "{count} tracers cannot determine a 3-D deformation gradient"
            );
        }
        assert!(deformation_gradient_from_tetrad(&tracers[..3]).is_none());
        assert!(deformation_gradient_from_tetrad(&tracers).is_some());
    }

    #[test]
    fn test_deformation_gradient_degenerate_reference_returns_identity() {
        let diag = MixingDiagnostic::new(1.0);
        // Four coplanar reference positions: the edge matrix is singular.
        let mut flat = tetrad(
            [0.0; 3],
            [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.1, 0.1, 0.0]],
        );
        advect(&mut flat, &saddle_flow, 1.0, 0.01);
        assert_eq!(diag.compute_deformation_gradient(&flat), MAT3_IDENTITY);
        assert!(deformation_gradient_from_tetrad(&flat).is_none());

        // Coincident tracers.
        let same = vec![TracerParticle::new([1.0, 2.0, 3.0], 0); 4];
        assert!(deformation_gradient_from_tetrad(&same).is_none());
    }

    #[test]
    fn test_deformation_gradient_scale_invariant_degeneracy_test() {
        // A thin but genuine tetrahedron of micron size must not be mistaken
        // for a degenerate one just because its volume is tiny.
        let diag = MixingDiagnostic::new(0.0);
        let mut tracers = tetrad(
            [0.0; 3],
            [[1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 1e-6]],
        );
        advect(&mut tracers, &shear_flow, 1.0, 0.01);
        let f = diag.compute_deformation_gradient(&tracers);
        assert_matrix_close(&f, &[1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 1e-6);
    }

    #[test]
    fn test_expected_stretch_is_exponential_in_time() {
        let diag = MixingDiagnostic::new(0.5);
        assert!((diag.expected_stretch(2.0) - 1.0_f64.exp()).abs() < 1e-14);
        assert_eq!(diag.expected_stretch(0.0), 1.0);
        assert!(
            (MixingDiagnostic::new(-1.0).expected_stretch(1.0) - (-1.0_f64).exp()).abs() < 1e-14
        );
    }

    #[test]
    fn test_max_principal_stretch_known_matrices() {
        assert!((max_principal_stretch(&MAT3_IDENTITY) - 1.0).abs() < 1e-14);
        let diag = [3.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 2.0];
        assert!((max_principal_stretch(&diag) - 3.0).abs() < 1e-12);
        let shear = [1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert!((max_principal_stretch(&shear) - (1.0 + 2.0_f64.sqrt())).abs() < 1e-12);
        assert_eq!(max_principal_stretch(&[0.0; 9]), 0.0);
    }

    #[test]
    fn test_sym3_max_eigenvalue_known_matrices() {
        // Diagonal.
        assert_eq!(
            sym3_max_eigenvalue(&[1.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0, 3.0]),
            7.0
        );
        // [[2, 1, 0], [1, 2, 0], [0, 0, 1]] has eigenvalues 3, 1, 1.
        let a = [2.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0];
        assert!((sym3_max_eigenvalue(&a) - 3.0).abs() < 1e-12);
        // Discrete Laplacian [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]: 2 + √2.
        let lap = [2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
        assert!((sym3_max_eigenvalue(&lap) - (2.0 + 2.0_f64.sqrt())).abs() < 1e-12);
        // Multiple of the identity.
        assert!(
            (sym3_max_eigenvalue(&[4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 4.0]) - 4.0).abs()
                < 1e-14
        );
    }

    #[test]
    fn test_mat3_inverse_roundtrip_and_singular() {
        let m = [2.0, 1.0, 0.5, -1.0, 3.0, 0.25, 0.0, 1.5, 4.0];
        let inv = mat3_inverse(&m).expect("regular matrix");
        assert_matrix_close(&mat3_mul(&m, &inv), &MAT3_IDENTITY, 1e-14);
        assert_matrix_close(&mat3_mul(&inv, &m), &MAT3_IDENTITY, 1e-14);
        let singular = [1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.0, 1.0, 0.0];
        assert!(mat3_inverse(&singular).is_none());
        assert_eq!(
            mat3_transpose(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
            [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]
        );
    }

    // ── Dispersion metrics ─────────────────────────────────────────────────────

    #[test]
    fn test_particle_dispersion_empty() {
        assert_eq!(particle_dispersion(&[]), 0.0);
    }

    #[test]
    fn test_particle_dispersion_colocated() {
        // Particles at same position with no velocity: MSD should be ≈ 0
        let tracers = vec![
            TracerParticle::new([0.0; 3], 0),
            TracerParticle::new([0.0; 3], 1),
        ];
        let msd = particle_dispersion(&tracers);
        assert!(msd < 1e-28);
    }

    #[test]
    fn test_particle_dispersion_uniform_translation() {
        // Every tracer is displaced by U t = 5 regardless of where it started.
        let mut tracers = vec![
            TracerParticle::new([0.0; 3], 0),
            TracerParticle::new([1.0, 2.0, 3.0], 1),
            TracerParticle::new([-4.0, 0.5, 0.0], 2),
        ];
        advect(&mut tracers, &uniform_flow, 5.0, 0.1);
        let msd = particle_dispersion(&tracers);
        assert!((msd - 25.0).abs() < 1e-9, "MSD = {msd}, expected 25");
    }

    #[test]
    fn test_particle_dispersion_matches_saddle_flow_closed_form() {
        // x(t) = X e^{λt}, y(t) = Y e^{-λt}: |x − X|² = X²(e^{λt} − 1)² + Y²(e^{-λt} − 1)².
        let starts = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [1.0, 1.0, 0.5]];
        let mut tracers: Vec<TracerParticle> = starts
            .iter()
            .enumerate()
            .map(|(id, x)| TracerParticle::new(*x, id))
            .collect();
        advect(&mut tracers, &saddle_flow, 2.0, 0.01);
        let growth = 1.0_f64.exp() - 1.0;
        let decay = (-1.0_f64).exp() - 1.0;
        let expected = starts
            .iter()
            .map(|x| x[0] * x[0] * growth * growth + x[1] * x[1] * decay * decay)
            .sum::<f64>()
            / 3.0;
        let msd = particle_dispersion(&tracers);
        assert!(
            (msd - expected).abs() < 1e-8 * expected,
            "MSD = {msd}, expected {expected}"
        );
    }

    #[test]
    fn test_particle_dispersion_ignores_velocity_and_age() {
        // The old estimator extrapolated x0 ≈ x − v·age; the recorded initial
        // positions must be used instead.
        let mut p = TracerParticle::new([1.0, 1.0, 1.0], 0);
        p.position = [1.0, 4.0, 5.0];
        p.velocity = [100.0, -50.0, 7.0];
        p.age = 9.0;
        assert!((particle_dispersion(&[p]) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_dispersion_single_particle_at_rest_is_zero() {
        let tracers = vec![TracerParticle::new([3.0, 4.0, 5.0], 0)];
        assert_eq!(particle_dispersion(&tracers), 0.0);
    }

    #[test]
    fn test_relative_dispersion_zero_same_pos() {
        let p1 = TracerParticle::new([1.0, 2.0, 3.0], 0);
        let p2 = TracerParticle::new([1.0, 2.0, 3.0], 1);
        assert_eq!(relative_dispersion(&p1, &p2), 0.0);
    }

    #[test]
    fn test_relative_dispersion_unit_separation() {
        let p1 = TracerParticle::new([0.0; 3], 0);
        let p2 = TracerParticle::new([1.0, 0.0, 0.0], 1);
        assert!((relative_dispersion(&p1, &p2) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_relative_dispersion_3d() {
        let p1 = TracerParticle::new([0.0; 3], 0);
        let p2 = TracerParticle::new([1.0, 1.0, 1.0], 1);
        let d = relative_dispersion(&p1, &p2);
        assert!((d - 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_relative_dispersion_grows_in_hyperbolic_flow() {
        let hyp_flow = |pos: [f64; 3], _t: f64| -> [f64; 3] { [pos[0], -pos[1], 0.0] };
        let mut p1 = TracerParticle::new([1.0, 0.0, 0.0], 0);
        let mut p2 = TracerParticle::new([1.001, 0.0, 0.0], 1);
        let d0 = relative_dispersion(&p1, &p2);
        for i in 0..100 {
            let t = i as f64 * 0.01;
            p1.integrate_rk4(&hyp_flow, 0.01, t);
            p2.integrate_rk4(&hyp_flow, 0.01, t);
        }
        let d1 = relative_dispersion(&p1, &p2);
        assert!(d1 > d0, "expected growing separation in hyperbolic flow");
    }

    // ── Additional coverage ────────────────────────────────────────────────────

    #[test]
    fn test_tracer_new_zero_velocity() {
        let p = TracerParticle::new([5.0, 3.0, 1.0], 42);
        assert_eq!(p.velocity, [0.0; 3]);
        assert_eq!(p.age, 0.0);
        assert_eq!(p.id, 42);
    }

    #[test]
    fn test_tracker_step_many_particles() {
        let positions: Vec<[f64; 3]> = (0..100).map(|i| [i as f64, 0.0, 0.0]).collect();
        let mut tracker = ParticleTracker::new(positions);
        tracker.step(&uniform_flow, 0.01, 0.0);
        for (i, p) in tracker.tracers.iter().enumerate() {
            let expected_x = i as f64 + 0.01;
            assert!((p.position[0] - expected_x).abs() < 1e-12);
        }
    }

    #[test]
    fn test_mixing_diagnostic_new() {
        let diag = MixingDiagnostic::new(2.5);
        assert!((diag.stretch_rate - 2.5).abs() < 1e-14);
    }

    #[test]
    fn test_ftle_new() {
        let ftl = FiniteTimeLyapunov::new(3.125);
        assert!((ftl.advection_time - 3.125).abs() < 1e-14);
    }

    #[test]
    fn test_poincare_section_time_varying() {
        // Time-varying oscillating flow
        let flow = |_pos: [f64; 3], t: f64| -> [f64; 3] { [t.sin(), t.cos(), 0.0] };
        let mut p = TracerParticle::new([0.0; 3], 0);
        for i in 0..100 {
            let t = i as f64 * 0.1;
            p.integrate_rk4(&flow, 0.1, t);
        }
        // Just check the particle remains finite
        assert!(p.position[0].is_finite());
        assert!(p.position[1].is_finite());
    }
}
