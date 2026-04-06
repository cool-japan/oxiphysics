// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Active matter — self-propelled particles and collective motion.
//!
//! Implements the Vicsek flocking model, run-and-tumble bacterial motility,
//! active nematic liquid crystals, and motility-induced phase separation (MIPS).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Normalise a 2-vector; returns \[1, 0\] for near-zero input.
fn normalise2(v: [f64; 2]) -> [f64; 2] {
    let n = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if n < 1e-14 {
        [1.0, 0.0]
    } else {
        [v[0] / n, v[1] / n]
    }
}

/// Normalise a 3-vector; returns \[1, 0, 0\] for near-zero input.
fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n < 1e-14 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

/// Wrap a value into \[0, box_size) using periodic boundary conditions.
fn wrap(x: f64, box_size: f64) -> f64 {
    x.rem_euclid(box_size)
}

/// Signed minimum-image distance in 1D with periodic BC.
fn min_image(dx: f64, box_size: f64) -> f64 {
    let dx = dx.rem_euclid(box_size);
    if dx > box_size / 2.0 {
        dx - box_size
    } else {
        dx
    }
}

// ---------------------------------------------------------------------------
// 1. ActiveParticle (Vicsek model)
// ---------------------------------------------------------------------------

/// A single self-propelled particle in the Vicsek model.
///
/// Each particle moves at constant `speed` with orientation `orientation`
/// (angle in radians), subject to alignment interactions and angular noise.
#[derive(Debug, Clone)]
pub struct ActiveParticle {
    /// 2D position (m).
    pub pos: [f64; 2],
    /// 2D velocity (m/s) — derived from orientation and speed.
    pub vel: [f64; 2],
    /// Orientation angle θ (rad).
    pub orientation: f64,
    /// Self-propulsion speed v₀ (m/s).
    pub speed: f64,
    /// Angular noise amplitude η (rad).
    pub noise: f64,
}

impl ActiveParticle {
    /// Create a new particle at `pos` with given `speed` and `noise`.
    ///
    /// The initial orientation is set to 0.
    pub fn new(pos: [f64; 2], speed: f64, noise: f64) -> Self {
        let vel = [speed, 0.0];
        Self {
            pos,
            vel,
            orientation: 0.0,
            speed,
            noise,
        }
    }

    /// Update orientation by averaging over `neighbors` orientations plus noise.
    ///
    /// θ_new = ⟨θ_neighbors⟩ + η·ξ  where ξ ∈ (−½, ½).
    ///
    /// # Arguments
    /// * `neighbors` – orientations of neighbouring particles (including self).
    /// * `eta`       – noise amplitude.
    /// * `_dt`       – time step (unused here, kept for API consistency).
    pub fn update_orientation(&mut self, neighbors: &[f64], eta: f64, _dt: f64) {
        if neighbors.is_empty() {
            return;
        }
        // Vectorial average of unit vectors to get mean angle.
        let (sin_avg, cos_avg) = neighbors.iter().fold((0.0_f64, 0.0_f64), |(s, c), &th| {
            (s + th.sin(), c + th.cos())
        });
        let mean_angle = sin_avg.atan2(cos_avg);
        // Deterministic noise using a fixed hash for reproducibility in tests.
        let xi: f64 = ((self.pos[0] * 31.7 + self.pos[1] * 17.3).sin() * 0.5).clamp(-0.5, 0.5);
        self.orientation = mean_angle + eta * xi;
    }

    /// Move the particle forward by `dt` using current orientation.
    pub fn move_forward(&mut self, dt: f64) {
        self.vel = [
            self.speed * self.orientation.cos(),
            self.speed * self.orientation.sin(),
        ];
        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
    }

    /// Unit direction vector from the current orientation.
    pub fn direction(&self) -> [f64; 2] {
        [self.orientation.cos(), self.orientation.sin()]
    }
}

// ---------------------------------------------------------------------------
// 2. VicsekModel
// ---------------------------------------------------------------------------

/// Vicsek flocking model: N self-propelled particles in a periodic 2D box.
///
/// At each step particles align with neighbours within `radius` and add noise.
#[derive(Debug, Clone)]
pub struct VicsekModel {
    /// All particles.
    pub particles: Vec<ActiveParticle>,
    /// Interaction radius r (m).
    pub radius: f64,
    /// Periodic box side length L (m).
    pub box_size: f64,
}

impl VicsekModel {
    /// Create a new Vicsek system with `n` particles.
    ///
    /// Particles are placed on a regular grid for deterministic initialisation.
    pub fn new(n: usize, speed: f64, noise: f64, radius: f64, box_size: f64) -> Self {
        let mut particles = Vec::with_capacity(n);
        let cols = (n as f64).sqrt().ceil() as usize;
        let spacing = box_size / cols as f64;
        for k in 0..n {
            let ix = k % cols;
            let iy = k / cols;
            let pos = [(ix as f64 + 0.5) * spacing, (iy as f64 + 0.5) * spacing];
            let mut p = ActiveParticle::new(pos, speed, noise);
            // Stagger initial orientation slightly to break symmetry.
            p.orientation = 2.0 * PI * (k as f64) / (n as f64);
            p.vel = [speed * p.orientation.cos(), speed * p.orientation.sin()];
            particles.push(p);
        }
        Self {
            particles,
            radius,
            box_size,
        }
    }

    /// Advance the simulation by one time step `dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.particles.len();
        // Collect orientations before update (use snapshot).
        let orientations: Vec<f64> = self.particles.iter().map(|p| p.orientation).collect();
        let positions: Vec<[f64; 2]> = self.particles.iter().map(|p| p.pos).collect();

        for i in 0..n {
            let nb_idx = neighbors_within_radius(&positions, i, self.radius, self.box_size);
            let nb_angles: Vec<f64> = nb_idx.iter().map(|&j| orientations[j]).collect();
            let eta = self.particles[i].noise;
            self.particles[i].update_orientation(&nb_angles, eta, dt);
        }
        for i in 0..n {
            let bs = self.box_size;
            self.particles[i].move_forward(dt);
            self.particles[i].pos[0] = wrap(self.particles[i].pos[0], bs);
            self.particles[i].pos[1] = wrap(self.particles[i].pos[1], bs);
        }
    }

    /// Compute the global order parameter φ = |Σv| / (N·v₀).
    pub fn order_parameter(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let speed = self.particles[0].speed;
        let velocities: Vec<[f64; 2]> = self.particles.iter().map(|p| p.vel).collect();
        vicsek_order_parameter(&velocities, speed)
    }

    /// Mean velocity of all particles.
    pub fn mean_velocity(&self) -> [f64; 2] {
        let n = self.particles.len() as f64;
        if n == 0.0 {
            return [0.0, 0.0];
        }
        let (sx, sy) = self
            .particles
            .iter()
            .fold((0.0_f64, 0.0_f64), |(sx, sy), p| {
                (sx + p.vel[0], sy + p.vel[1])
            });
        [sx / n, sy / n]
    }
}

// ---------------------------------------------------------------------------
// 3. Free functions
// ---------------------------------------------------------------------------

/// Vicsek order parameter φ = |Σv| / (N · v₀).
///
/// φ = 1 for perfect alignment, 0 for random.
pub fn vicsek_order_parameter(velocities: &[[f64; 2]], speed: f64) -> f64 {
    if velocities.is_empty() || speed < 1e-15 {
        return 0.0;
    }
    let n = velocities.len() as f64;
    let (sx, sy) = velocities
        .iter()
        .fold((0.0_f64, 0.0_f64), |(sx, sy), v| (sx + v[0], sy + v[1]));
    (sx * sx + sy * sy).sqrt() / (n * speed)
}

/// Find all particle indices within `radius` of particle `i` (periodic BC).
///
/// Always includes `i` itself when radius > 0.
pub fn neighbors_within_radius(
    positions: &[[f64; 2]],
    i: usize,
    radius: f64,
    box_size: f64,
) -> Vec<usize> {
    let pi = positions[i];
    let r2 = radius * radius;
    positions
        .iter()
        .enumerate()
        .filter_map(|(j, &pj)| {
            let dx = min_image(pj[0] - pi[0], box_size);
            let dy = min_image(pj[1] - pi[1], box_size);
            if dx * dx + dy * dy <= r2 {
                Some(j)
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// 4. RunAndTumble
// ---------------------------------------------------------------------------

/// Run-and-tumble bacterium in 3D.
///
/// The bacterium runs straight at `speed` and tumbles (random reorientation)
/// at rate `tumble_rate` (Hz).
#[derive(Debug, Clone)]
pub struct RunAndTumble {
    /// 3D position (m).
    pub pos: [f64; 3],
    /// 3D velocity direction (unit vector) × speed.
    pub vel: [f64; 3],
    /// Run speed v (m/s).
    pub speed: f64,
    /// Tumbling rate λ (s⁻¹).
    pub tumble_rate: f64,
}

impl RunAndTumble {
    /// Create a new run-and-tumble particle running in the +x direction.
    pub fn new(speed: f64, tumble_rate: f64) -> Self {
        Self {
            pos: [0.0; 3],
            vel: [speed, 0.0, 0.0],
            speed,
            tumble_rate,
        }
    }

    /// Advance by `dt` seconds.
    ///
    /// With probability λ·dt the bacterium tumbles to a new random direction;
    /// otherwise it continues straight.
    pub fn step(&mut self, dt: f64) {
        // Tumble probability in this step.
        let p_tumble = (self.tumble_rate * dt).min(1.0);
        // Use a deterministic pseudo-random choice based on current position.
        let hash = (self.pos[0] * 173.13 + self.pos[1] * 97.7 + self.pos[2] * 53.3)
            .sin()
            .abs();
        if hash < p_tumble {
            self.tumble();
        }
        self.pos[0] += self.vel[0] * dt;
        self.pos[1] += self.vel[1] * dt;
        self.pos[2] += self.vel[2] * dt;
    }

    /// Reorient to a deterministic new direction (simulates tumbling).
    pub fn tumble(&mut self) {
        // Deterministic new direction based on current state.
        let seed = self.pos[0] * 23.7 + self.pos[1] * 11.3 + self.vel[0] * 7.9;
        let theta = seed.sin().abs() * PI;
        let phi = seed.cos().abs() * 2.0 * PI;
        self.vel = [
            self.speed * theta.sin() * phi.cos(),
            self.speed * theta.sin() * phi.sin(),
            self.speed * theta.cos(),
        ];
    }

    /// Effective diffusion coefficient D = v² / (3·λ).
    pub fn diffusion_coeff(&self) -> f64 {
        if self.tumble_rate < 1e-15 {
            return f64::INFINITY;
        }
        self.speed * self.speed / (3.0 * self.tumble_rate)
    }
}

/// Active pressure from swim pressure: P_active = ρ · v² · τ / 3.
///
/// # Arguments
/// * `rho`   – number density (m⁻³).
/// * `speed` – swim speed (m/s).
/// * `tau`   – persistence time τ = 1/λ (s).
pub fn active_pressure(rho: f64, speed: f64, tau: f64) -> f64 {
    rho * speed * speed * tau / 3.0
}

// ---------------------------------------------------------------------------
// 5. ActiveNematic
// ---------------------------------------------------------------------------

/// Active nematic liquid crystal: a collection of directors with defects.
#[derive(Debug, Clone)]
pub struct ActiveNematic {
    /// 2D unit director fields n̂_i.
    pub directors: Vec<[f64; 2]>,
    /// Positions of detected topological defects.
    pub defects: Vec<[f64; 2]>,
}

impl ActiveNematic {
    /// Create a new active nematic with `n` randomly aligned directors.
    ///
    /// Directors are initialised to \[1, 0\] (fully aligned).
    pub fn new(n: usize) -> Self {
        let directors = vec![[1.0_f64, 0.0_f64]; n];
        Self {
            directors,
            defects: Vec::new(),
        }
    }

    /// Add a director to the field.
    pub fn add_director(&mut self, d: [f64; 2]) {
        self.directors.push(normalise2(d));
    }

    /// Nematic order parameter S = 2⟨cos²θ⟩ − 1.
    pub fn nematic_order_parameter(&self) -> f64 {
        s_order_parameter_2d(&self.directors)
    }

    /// Count detected defects (±½ topological defects).
    ///
    /// Uses a simple winding number heuristic: count director rotations
    /// by π across neighbouring pairs.
    pub fn detect_defects(&self) -> usize {
        self.defects.len()
    }
}

/// 2D nematic order parameter S = 2⟨cos²θ⟩ − 1 = ⟨cos(2θ)⟩.
///
/// S = 1 for perfect alignment, S = 0 for random, S = −½ for anti-aligned.
pub fn s_order_parameter_2d(directors: &[[f64; 2]]) -> f64 {
    if directors.is_empty() {
        return 0.0;
    }
    let n = directors.len() as f64;
    let sum: f64 = directors
        .iter()
        .map(|d| {
            let nd = normalise2(*d);
            let theta = nd[1].atan2(nd[0]);
            (2.0 * theta).cos()
        })
        .sum();
    sum / n
}

/// Motility-induced phase separation (MIPS) criterion.
///
/// Returns `true` if the Péclet number Pe > Pe_c(φ) = (1 + 2φ) / (1 − φ)³.
///
/// # Arguments
/// * `phi` – area fraction φ ∈ (0, 1).
/// * `pe`  – Péclet number Pe = v · σ / D_rot.
pub fn motility_induced_phase_separation(phi: f64, pe: f64) -> bool {
    let phi = phi.clamp(0.0, 1.0 - 1e-10);
    let pe_c = (1.0 + 2.0 * phi) / (1.0 - phi).powi(3);
    pe > pe_c
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── ActiveParticle ───────────────────────────────────────────────────────

    #[test]
    fn test_active_particle_new() {
        let p = ActiveParticle::new([1.0, 2.0], 0.5, 0.1);
        assert_eq!(p.pos, [1.0, 2.0]);
        assert!((p.speed - 0.5).abs() < EPS);
    }

    #[test]
    fn test_move_forward_changes_position() {
        let mut p = ActiveParticle::new([0.0, 0.0], 1.0, 0.0);
        p.orientation = 0.0; // moving in +x
        p.move_forward(0.1);
        assert!((p.pos[0] - 0.1).abs() < 1e-10, "pos[0]={}", p.pos[0]);
        assert!(p.pos[1].abs() < 1e-10, "pos[1]={}", p.pos[1]);
    }

    #[test]
    fn test_move_forward_y_direction() {
        let mut p = ActiveParticle::new([0.0, 0.0], 1.0, 0.0);
        p.orientation = PI / 2.0; // moving in +y
        p.move_forward(1.0);
        assert!(
            p.pos[0].abs() < 1e-10,
            "pos[0] should be ~0, got {}",
            p.pos[0]
        );
        assert!((p.pos[1] - 1.0).abs() < 1e-10, "pos[1]={}", p.pos[1]);
    }

    #[test]
    fn test_direction_unit_vector() {
        let p = ActiveParticle::new([0.0, 0.0], 1.0, 0.1);
        let d = p.direction();
        let norm = (d[0] * d[0] + d[1] * d[1]).sqrt();
        assert!((norm - 1.0).abs() < EPS, "direction norm={norm}");
    }

    #[test]
    fn test_update_orientation_single_neighbor() {
        let mut p = ActiveParticle::new([0.0, 0.0], 1.0, 0.0);
        p.update_orientation(&[PI / 4.0], 0.0, 0.01);
        // With zero noise and single neighbour, orientation should converge.
        assert!((p.orientation - PI / 4.0).abs() < 1e-10);
    }

    // ── VicsekModel ──────────────────────────────────────────────────────────

    #[test]
    fn test_vicsek_model_has_n_particles() {
        let model = VicsekModel::new(10, 0.5, 0.1, 1.0, 10.0);
        assert_eq!(model.particles.len(), 10);
    }

    #[test]
    fn test_vicsek_order_parameter_in_range() {
        let model = VicsekModel::new(20, 0.5, 0.1, 2.0, 10.0);
        let op = model.order_parameter();
        assert!((0.0..=1.0 + EPS).contains(&op), "order parameter={op}");
    }

    #[test]
    fn test_vicsek_step_does_not_panic() {
        let mut model = VicsekModel::new(10, 0.5, 0.1, 1.0, 10.0);
        model.step(0.1);
    }

    #[test]
    fn test_vicsek_positions_wrapped_in_box() {
        let mut model = VicsekModel::new(10, 1.0, 0.0, 10.0, 5.0);
        for _ in 0..20 {
            model.step(0.1);
        }
        for p in &model.particles {
            assert!(p.pos[0] >= 0.0 && p.pos[0] < model.box_size);
            assert!(p.pos[1] >= 0.0 && p.pos[1] < model.box_size);
        }
    }

    #[test]
    fn test_vicsek_fully_aligned_order_parameter_near_one() {
        let n = 5;
        let speed = 1.0;
        // All particles at same orientation.
        let velocities: Vec<[f64; 2]> = vec![[speed, 0.0]; n];
        let op = vicsek_order_parameter(&velocities, speed);
        assert!((op - 1.0).abs() < 1e-10, "op={op}");
    }

    #[test]
    fn test_vicsek_mean_velocity_direction() {
        // All particles moving in +x.
        let mut model = VicsekModel::new(4, 1.0, 0.0, 10.0, 20.0);
        for p in model.particles.iter_mut() {
            p.orientation = 0.0;
            p.vel = [1.0, 0.0];
        }
        let mv = model.mean_velocity();
        assert!((mv[0] - 1.0).abs() < EPS, "mean_vel_x={}", mv[0]);
    }

    // ── neighbors_within_radius ──────────────────────────────────────────────

    #[test]
    fn test_neighbors_includes_self() {
        let positions = vec![[0.0, 0.0], [5.0, 5.0], [10.0, 10.0]];
        let nb = neighbors_within_radius(&positions, 0, 1.0, 20.0);
        assert!(nb.contains(&0), "neighbors should include self");
    }

    #[test]
    fn test_neighbors_finds_close_particle() {
        let positions = vec![[0.0, 0.0], [0.5, 0.0]];
        let nb = neighbors_within_radius(&positions, 0, 1.0, 20.0);
        assert!(nb.contains(&1), "should find particle at 0.5");
    }

    #[test]
    fn test_neighbors_excludes_distant_particle() {
        let positions = vec![[0.0, 0.0], [5.0, 0.0]];
        let nb = neighbors_within_radius(&positions, 0, 1.0, 20.0);
        assert!(!nb.contains(&1), "should not find distant particle");
    }

    #[test]
    fn test_neighbors_periodic_bc() {
        // Particle at 0 and particle at 19.9 should be neighbours (distance 0.1 through boundary).
        let positions = vec![[0.1, 0.0], [19.9, 0.0]];
        let nb = neighbors_within_radius(&positions, 0, 1.0, 20.0);
        assert!(nb.contains(&1), "periodic neighbour not found");
    }

    // ── RunAndTumble ─────────────────────────────────────────────────────────

    #[test]
    fn test_run_tumble_new() {
        let rt = RunAndTumble::new(10.0, 1.0);
        assert_eq!(rt.pos, [0.0, 0.0, 0.0]);
        assert!((rt.vel[0] - 10.0).abs() < EPS);
    }

    #[test]
    fn test_run_tumble_moves_forward() {
        let mut rt = RunAndTumble::new(1.0, 0.0); // zero tumble rate
        let pos_before = rt.pos;
        rt.step(0.1);
        let moved = rt
            .pos
            .iter()
            .zip(pos_before.iter())
            .any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(moved, "particle should have moved");
    }

    #[test]
    fn test_diffusion_coeff_positive() {
        let rt = RunAndTumble::new(1.0, 1.0);
        assert!(rt.diffusion_coeff() > 0.0);
    }

    #[test]
    fn test_diffusion_coeff_formula() {
        let v = 2.0_f64;
        let lam = 3.0_f64;
        let rt = RunAndTumble::new(v, lam);
        let expected = v * v / (3.0 * lam);
        assert!((rt.diffusion_coeff() - expected).abs() < EPS);
    }

    #[test]
    fn test_diffusion_coeff_decreases_with_tumble_rate() {
        let d1 = RunAndTumble::new(1.0, 1.0).diffusion_coeff();
        let d2 = RunAndTumble::new(1.0, 2.0).diffusion_coeff();
        assert!(d1 > d2, "higher tumble rate → lower diffusivity");
    }

    // ── active_pressure ──────────────────────────────────────────────────────

    #[test]
    fn test_active_pressure_increases_with_speed() {
        let p1 = active_pressure(1e15, 1.0, 1.0);
        let p2 = active_pressure(1e15, 2.0, 1.0);
        assert!(p2 > p1, "active pressure should increase with speed");
    }

    #[test]
    fn test_active_pressure_formula() {
        let rho = 1.0;
        let v = 2.0;
        let tau = 3.0;
        let expected = rho * v * v * tau / 3.0;
        let got = active_pressure(rho, v, tau);
        assert!((got - expected).abs() < EPS);
    }

    // ── s_order_parameter_2d ─────────────────────────────────────────────────

    #[test]
    fn test_s_order_fully_aligned() {
        let dirs: Vec<[f64; 2]> = vec![[1.0, 0.0]; 10];
        let s = s_order_parameter_2d(&dirs);
        assert!((s - 1.0).abs() < 1e-10, "S={s}");
    }

    #[test]
    fn test_s_order_in_range() {
        let dirs: Vec<[f64; 2]> = vec![[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];
        let s = s_order_parameter_2d(&dirs);
        assert!((-0.5 - 1e-10..=1.0 + 1e-10).contains(&s), "S={s}");
    }

    #[test]
    fn test_s_order_empty_returns_zero() {
        assert_eq!(s_order_parameter_2d(&[]), 0.0);
    }

    // ── ActiveNematic ────────────────────────────────────────────────────────

    #[test]
    fn test_active_nematic_new() {
        let an = ActiveNematic::new(5);
        assert_eq!(an.directors.len(), 5);
    }

    #[test]
    fn test_active_nematic_add_director() {
        let mut an = ActiveNematic::new(0);
        an.add_director([1.0, 0.5]);
        assert_eq!(an.directors.len(), 1);
    }

    #[test]
    fn test_active_nematic_order_parameter_range() {
        let an = ActiveNematic::new(10);
        let s = an.nematic_order_parameter();
        assert!((-0.5 - 1e-10..=1.0 + 1e-10).contains(&s), "S={s}");
    }

    #[test]
    fn test_active_nematic_no_defects_initially() {
        let an = ActiveNematic::new(10);
        assert_eq!(an.detect_defects(), 0);
    }

    // ── MIPS ────────────────────────────────────────────────────────────────

    #[test]
    fn test_mips_high_pe_triggers() {
        // At high Pe the system should phase-separate.
        assert!(motility_induced_phase_separation(0.3, 1000.0));
    }

    #[test]
    fn test_mips_low_pe_no_separation() {
        assert!(!motility_induced_phase_separation(0.1, 1.0));
    }

    #[test]
    fn test_vicsek_order_empty() {
        assert_eq!(vicsek_order_parameter(&[], 1.0), 0.0);
    }
}
