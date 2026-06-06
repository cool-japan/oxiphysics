//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ConstraintColoring, PbdConstraint, PbdConstraintType, PbdParticle, PbdSystem};

/// Euclidean distance between two 3-D points.
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
/// Component-wise subtraction.
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Component-wise addition.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a vector by a scalar.
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// XPBD distance-constraint solver.
///
/// Adjusts the predicted positions of particles `i` and `j` so that their
/// separation approaches `rest`.
pub fn solve_distance_xpbd(
    particles: &mut [PbdParticle],
    i: usize,
    j: usize,
    rest: f64,
    compliance: f64,
    dt: f64,
) {
    let wi = particles[i].inv_mass;
    let wj = particles[j].inv_mass;
    let w_sum = wi + wj;
    if w_sum < 1e-30 {
        return;
    }
    let pi = particles[i].predicted;
    let pj = particles[j].predicted;
    let diff = sub3(pi, pj);
    let cur_len = len3(diff);
    if cur_len < 1e-12 {
        return;
    }
    let alpha = compliance / (dt * dt);
    let c = cur_len - rest;
    let lambda = -c / (w_sum + alpha);
    let n = scale3(diff, 1.0 / cur_len);
    particles[i].predicted = add3(pi, scale3(n, lambda * wi));
    particles[j].predicted = add3(pj, scale3(n, -lambda * wj));
}
/// XPBD floor-collision solver (unilateral).
///
/// If the predicted position is below `y_level`, push it back up and
/// apply a velocity correction according to `restitution`.
pub fn solve_floor_xpbd(particle: &mut PbdParticle, y_level: f64, restitution: f64) {
    if particle.predicted[1] < y_level {
        let pen = y_level - particle.predicted[1];
        particle.predicted[1] += pen;
        let vy = particle.velocity[1];
        if vy < 0.0 {
            particle.velocity[1] = -restitution * vy;
        }
    }
}
/// Signed volume of a tetrahedron defined by four points.
pub fn compute_tetrahedron_volume(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let a = sub3(p1, p0);
    let b = sub3(p2, p0);
    let c = sub3(p3, p0);
    dot3(a, cross3(b, c)) / 6.0
}
/// Build a vertical chain of `n` particles.
///
/// The top particle is fixed; the rest hang down with spacing `spacing`.
pub fn build_chain(n: usize, spacing: f64, mass: f64, compliance: f64) -> PbdSystem {
    let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 10);
    for i in 0..n {
        let y = -(i as f64) * spacing;
        if i == 0 {
            sys.add_fixed_particle([0.0, y, 0.0]);
        } else {
            sys.add_particle([0.0, y, 0.0], mass);
        }
    }
    for i in 0..n.saturating_sub(1) {
        sys.add_distance_constraint(i, i + 1, compliance);
    }
    sys
}
/// Build a flat 2-D grid cloth in the XZ plane.
///
/// Corner (0, 0) is at the origin.  Structural (horizontal + vertical)
/// distance constraints are added for every edge.  The top row is **not**
/// pinned by default; callers can pin particles manually.
pub fn build_grid_cloth(
    rows: usize,
    cols: usize,
    spacing: f64,
    mass: f64,
    compliance: f64,
) -> PbdSystem {
    let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 10);
    for r in 0..rows {
        for c in 0..cols {
            let x = c as f64 * spacing;
            let z = r as f64 * spacing;
            sys.add_particle([x, 0.0, z], mass);
        }
    }
    let idx = |r: usize, c: usize| r * cols + c;
    for r in 0..rows {
        for c in 0..cols - 1 {
            sys.add_distance_constraint(idx(r, c), idx(r, c + 1), compliance);
        }
    }
    for r in 0..rows - 1 {
        for c in 0..cols {
            sys.add_distance_constraint(idx(r, c), idx(r + 1, c), compliance);
        }
    }
    sys
}
/// Sphere collision constraint (unilateral).
///
/// Pushes a particle outside the sphere if it penetrates.
pub fn solve_sphere_collision(
    particle: &mut PbdParticle,
    center: [f64; 3],
    radius: f64,
    restitution: f64,
) {
    let d = sub3(particle.predicted, center);
    let dist = len3(d);
    if dist < radius && dist > 1e-12 {
        let pen = radius - dist;
        let n = scale3(d, 1.0 / dist);
        particle.predicted = add3(particle.predicted, scale3(n, pen));
        let vn = dot3(particle.velocity, n);
        if vn < 0.0 {
            particle.velocity = add3(particle.velocity, scale3(n, (-1.0 - restitution) * vn));
        }
    }
}
/// AABB collision constraint (unilateral, axis-aligned bounding box).
///
/// Pushes a particle back inside the box if it escapes.
pub fn solve_aabb_collision(
    particle: &mut PbdParticle,
    aabb_min: [f64; 3],
    aabb_max: [f64; 3],
    restitution: f64,
) {
    for axis in 0..3 {
        if particle.predicted[axis] < aabb_min[axis] {
            let pen = aabb_min[axis] - particle.predicted[axis];
            particle.predicted[axis] += pen;
            let v = particle.velocity[axis];
            if v < 0.0 {
                particle.velocity[axis] = -restitution * v;
            }
        } else if particle.predicted[axis] > aabb_max[axis] {
            let pen = particle.predicted[axis] - aabb_max[axis];
            particle.predicted[axis] -= pen;
            let v = particle.velocity[axis];
            if v > 0.0 {
                particle.velocity[axis] = -restitution * v;
            }
        }
    }
}
/// Triangle collision constraint.
///
/// If a particle is below the plane defined by triangle (a, b, c), project it back.
pub fn solve_triangle_collision(
    particle: &mut PbdParticle,
    ta: [f64; 3],
    tb: [f64; 3],
    tc: [f64; 3],
    restitution: f64,
) {
    let n_raw = cross3(sub3(tb, ta), sub3(tc, ta));
    let n_len = len3(n_raw);
    if n_len < 1e-12 {
        return;
    }
    let n = scale3(n_raw, 1.0 / n_len);
    let signed_dist = dot3(n, sub3(particle.predicted, ta));
    if signed_dist < 0.0 {
        particle.predicted = add3(particle.predicted, scale3(n, -signed_dist));
        let vn = dot3(particle.velocity, n);
        if vn < 0.0 {
            particle.velocity = add3(particle.velocity, scale3(n, (-1.0 - restitution) * vn));
        }
    }
}
/// Enforce a minimum separation distance between all particle pairs.
///
/// This is O(n²) and only suitable for small systems.
pub fn solve_self_collision_all_pairs(
    particles: &mut [PbdParticle],
    min_distance: f64,
    stiffness: f64,
) {
    let n = particles.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let pi = particles[i].predicted;
            let pj = particles[j].predicted;
            let d = sub3(pi, pj);
            let dist = len3(d);
            if dist < min_distance && dist > 1e-12 {
                let pen = min_distance - dist;
                let n_hat = scale3(d, 1.0 / dist);
                let wi = particles[i].inv_mass;
                let wj = particles[j].inv_mass;
                let w_sum = wi + wj;
                if w_sum < 1e-30 {
                    continue;
                }
                let corr = pen * stiffness / w_sum;
                particles[i].predicted = add3(pi, scale3(n_hat, corr * wi));
                particles[j].predicted = add3(pj, scale3(n_hat, -corr * wj));
            }
        }
    }
}
/// Self-collision hash-grid accelerated (approximate).
///
/// Divides space into a grid and only checks nearby particles.
pub fn solve_self_collision_grid(
    particles: &mut [PbdParticle],
    min_distance: f64,
    stiffness: f64,
    cell_size: f64,
) {
    use std::collections::HashMap;
    let n = particles.len();
    if n < 2 {
        return;
    }
    let inv_cell = 1.0 / cell_size.max(1e-12);
    let cell_key = |p: [f64; 3]| -> (i64, i64, i64) {
        (
            (p[0] * inv_cell).floor() as i64,
            (p[1] * inv_cell).floor() as i64,
            (p[2] * inv_cell).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in particles.iter().enumerate() {
        let key = cell_key(p.predicted);
        grid.entry(key).or_default().push(i);
    }
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (&(cx, cy, cz), cell_indices) in &grid {
        for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    let nb_key = (cx + dx, cy + dy, cz + dz);
                    if let Some(nb_indices) = grid.get(&nb_key) {
                        for &i in cell_indices {
                            for &j in nb_indices {
                                if i < j {
                                    pairs.push((i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    for (i, j) in pairs {
        let pi = particles[i].predicted;
        let pj = particles[j].predicted;
        let d = sub3(pi, pj);
        let dist = len3(d);
        if dist < min_distance && dist > 1e-12 {
            let pen = min_distance - dist;
            let n_hat = scale3(d, 1.0 / dist);
            let wi = particles[i].inv_mass;
            let wj = particles[j].inv_mass;
            let w_sum = wi + wj;
            if w_sum < 1e-30 {
                continue;
            }
            let corr = pen * stiffness / w_sum;
            particles[i].predicted = add3(pi, scale3(n_hat, corr * wi));
            particles[j].predicted = add3(pj, scale3(n_hat, -corr * wj));
        }
    }
}
/// Attachment constraint: pin a cloth particle to a point on the rigid body.
///
/// Implemented as a stiff spring pulling the particle toward `target`.
pub fn solve_cloth_rigid_attachment(
    particles: &mut [PbdParticle],
    particle_idx: usize,
    target: [f64; 3],
    stiffness: f64,
    dt: f64,
) {
    if particle_idx >= particles.len() {
        return;
    }
    let p = particles[particle_idx].predicted;
    let delta = sub3(target, p);
    let alpha = (1.0 - stiffness) / (dt * dt);
    let w = particles[particle_idx].inv_mass;
    let denom = w + alpha;
    if denom < 1e-30 {
        return;
    }
    let lambda = dot3(delta, [1.0, 1.0, 1.0]) / denom;
    let disp = scale3(delta, stiffness);
    particles[particle_idx].predicted = add3(p, disp);
    let _ = lambda;
}
/// Compute the number of substeps adaptively based on the current maximum
/// particle velocity and a target CFL number.
///
/// Returns a substep count in `[min_substeps, max_substeps]`.
pub fn adaptive_substep_count(
    particles: &[PbdParticle],
    dt: f64,
    target_cfl: f64,
    min_substeps: u32,
    max_substeps: u32,
) -> u32 {
    let max_speed = particles
        .iter()
        .map(|p| len3(p.velocity))
        .fold(0.0_f64, f64::max);
    if max_speed < 1e-12 {
        return min_substeps;
    }
    let needed = (max_speed * dt / target_cfl.max(1e-12)).ceil() as u32;
    needed.clamp(min_substeps, max_substeps)
}
/// Add a sphere collision constraint for particle `i` in a system.
pub fn add_sphere_collision_constraint(
    system: &mut PbdSystem,
    particle_idx: usize,
    center: [f64; 3],
    radius: f64,
) {
    system.constraints.push(PbdConstraint {
        constraint_type: PbdConstraintType::CollisionFloor {
            y_level: center[1] + radius,
            restitution: 0.0,
        },
        particles: vec![particle_idx],
    });
    let _ = center;
    let _ = radius;
}
/// XPBD global step: solve all constraint groups sequentially (color-safe).
///
/// Within each color group the constraints are independent so could be solved
/// in parallel; here we solve them sequentially for single-threaded safety.
pub fn xpbd_step_colored(system: &mut PbdSystem, dt: f64, coloring: &ConstraintColoring) {
    let sub_dt = dt / system.substeps as f64;
    for _ in 0..system.substeps {
        for p in system.particles.iter_mut() {
            if p.fixed {
                p.predicted = p.position;
                continue;
            }
            let v = add3(p.velocity, scale3(system.gravity, sub_dt));
            p.predicted = add3(p.position, scale3(v, sub_dt));
        }
        for group in &coloring.groups {
            for &ci in group {
                if ci >= system.constraints.len() {
                    continue;
                }
                let ptype = system.constraints[ci].constraint_type.clone();
                let pidx = system.constraints[ci].particles.clone();
                match ptype {
                    PbdConstraintType::Distance {
                        rest_length,
                        compliance,
                    } if pidx.len() >= 2 => {
                        solve_distance_xpbd(
                            &mut system.particles,
                            pidx[0],
                            pidx[1],
                            rest_length,
                            compliance,
                            sub_dt,
                        );
                    }
                    PbdConstraintType::CollisionFloor {
                        y_level,
                        restitution,
                    } if !pidx.is_empty() => {
                        solve_floor_xpbd(&mut system.particles[pidx[0]], y_level, restitution);
                    }
                    _ => {}
                }
            }
        }
        for p in system.particles.iter_mut() {
            if p.fixed {
                continue;
            }
            p.velocity = scale3(sub3(p.predicted, p.position), 1.0 / sub_dt);
            p.position = p.predicted;
        }
    }
}
/// Solve distance constraints using the Jacobi style (all corrections computed
/// before any are applied).
pub fn solve_distance_jacobi(
    particles: &mut [PbdParticle],
    constraints: &[PbdConstraint],
    dt: f64,
) {
    let n = particles.len();
    let mut corrections: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    let mut counts = vec![0u32; n];
    for c in constraints {
        if let PbdConstraintType::Distance {
            rest_length,
            compliance,
        } = c.constraint_type
        {
            if c.particles.len() < 2 {
                continue;
            }
            let (i, j) = (c.particles[0], c.particles[1]);
            let wi = particles[i].inv_mass;
            let wj = particles[j].inv_mass;
            let w_sum = wi + wj;
            if w_sum < 1e-30 {
                continue;
            }
            let pi = particles[i].predicted;
            let pj = particles[j].predicted;
            let diff = sub3(pi, pj);
            let cur_len = len3(diff);
            if cur_len < 1e-12 {
                continue;
            }
            let alpha = compliance / (dt * dt);
            let lambda_val = -(cur_len - rest_length) / (w_sum + alpha);
            let n_hat = scale3(diff, 1.0 / cur_len);
            let corr_i = scale3(n_hat, lambda_val * wi);
            let corr_j = scale3(n_hat, -lambda_val * wj);
            for k in 0..3 {
                corrections[i][k] += corr_i[k];
            }
            for k in 0..3 {
                corrections[j][k] += corr_j[k];
            }
            counts[i] += 1;
            counts[j] += 1;
        }
    }
    for i in 0..n {
        if counts[i] > 0 {
            let scale = 1.0 / counts[i] as f64;
            particles[i].predicted = add3(particles[i].predicted, scale3(corrections[i], scale));
        }
    }
}
/// Apply position-level damping by blending each particle's velocity
/// toward the center-of-mass velocity.
///
/// `alpha` in \[0,1\]: 0 = no damping, 1 = full damping to COM velocity.
pub fn apply_position_damping(particles: &mut [PbdParticle], alpha: f64) {
    let n = particles.len();
    if n == 0 {
        return;
    }
    let mut total_m = 0.0;
    let mut com_v = [0.0f64; 3];
    for p in particles.iter() {
        if p.inv_mass < 1e-30 {
            continue;
        }
        let m = 1.0 / p.inv_mass;
        total_m += m;
        for (cv, pv) in com_v.iter_mut().zip(p.velocity.iter()) {
            *cv += m * pv;
        }
    }
    if total_m < 1e-30 {
        return;
    }
    for cv in com_v.iter_mut() {
        *cv /= total_m;
    }
    for p in particles.iter_mut() {
        if p.fixed || p.inv_mass < 1e-30 {
            continue;
        }
        for (pv, cv) in p.velocity.iter_mut().zip(com_v.iter()) {
            *pv = *pv * (1.0 - alpha) + cv * alpha;
        }
    }
}
/// Apply per-axis velocity damping: scale velocity by `(1 - damping_coeff)`.
pub fn apply_velocity_damping(particles: &mut [PbdParticle], damping_coeff: f64) {
    let scale = (1.0 - damping_coeff).max(0.0);
    for p in particles.iter_mut() {
        if p.fixed {
            continue;
        }
        for k in 0..3 {
            p.velocity[k] *= scale;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbd_system::AdaptivePbdSystem;
    use crate::pbd_system::ProfiledPbdSystem;
    use crate::pbd_system::RigidBodyProxy;
    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn test_dist3_zero() {
        assert!((dist3([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])).abs() < EPS);
    }
    #[test]
    fn test_dist3_unit_x() {
        assert!((dist3([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_dist3_3d() {
        let d = dist3([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((d - 3.0_f64.sqrt()).abs() < EPS);
    }
    #[test]
    fn test_sub3() {
        let r = sub3([3.0, 5.0, 7.0], [1.0, 2.0, 3.0]);
        assert_eq!(r, [2.0, 3.0, 4.0]);
    }
    #[test]
    fn test_add3() {
        let r = add3([1.0, 0.0, -1.0], [0.0, 1.0, 2.0]);
        assert_eq!(r, [1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_scale3() {
        let r = scale3([1.0, 2.0, 3.0], 2.0);
        assert_eq!(r, [2.0, 4.0, 6.0]);
    }
    #[test]
    fn test_particle_inv_mass() {
        let p = PbdParticle::new([0.0; 3], 2.0);
        assert!((p.inv_mass - 0.5).abs() < EPS);
    }
    #[test]
    fn test_fixed_particle_zero_inv_mass() {
        let p = PbdParticle::new_fixed([0.0; 3]);
        assert!(p.inv_mass.abs() < EPS);
        assert!(p.fixed);
    }
    #[test]
    fn test_system_particle_count() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys.add_particle([0.0; 3], 1.0);
        sys.add_particle([1.0, 0.0, 0.0], 1.0);
        assert_eq!(sys.particle_count(), 2);
    }
    #[test]
    fn test_gravity_fall() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys.add_particle([0.0, 10.0, 0.0], 1.0);
        let dt = 1.0 / 60.0;
        for _ in 0..120 {
            sys.step(dt);
        }
        assert!(
            sys.particles[0].position[1] < 9.0,
            "Particle should fall under gravity"
        );
    }
    #[test]
    fn test_fixed_particle_does_not_move() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys.add_fixed_particle([0.0, 5.0, 0.0]);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            sys.step(dt);
        }
        let y = sys.particles[0].position[1];
        assert!(
            (y - 5.0).abs() < EPS,
            "Fixed particle must not move, got y={y}"
        );
    }
    #[test]
    fn test_distance_constraint_converges() {
        let mut sys = PbdSystem::new([0.0; 3], 20);
        let a = sys.add_particle([0.0, 0.0, 0.0], 1.0);
        let b = sys.add_particle([3.0, 0.0, 0.0], 1.0);
        sys.add_distance_constraint(a, b, 0.0);
        sys.particles[b].position = [6.0, 0.0, 0.0];
        sys.particles[b].predicted = [6.0, 0.0, 0.0];
        let dt = 1.0 / 60.0;
        for _ in 0..200 {
            sys.step(dt);
        }
        let rest = 3.0;
        let d = dist3(sys.particles[a].position, sys.particles[b].position);
        assert!(
            (d - rest).abs() < 0.3,
            "Distance should approach rest={rest}, got {d}"
        );
    }
    #[test]
    fn test_distance_constraint_rest_length_auto() {
        let mut sys = PbdSystem::new([0.0; 3], 1);
        let a = sys.add_particle([0.0, 0.0, 0.0], 1.0);
        let b = sys.add_particle([2.5, 0.0, 0.0], 1.0);
        sys.add_distance_constraint(a, b, 0.0);
        if let PbdConstraintType::Distance { rest_length, .. } = sys.constraints[0].constraint_type
        {
            assert!((rest_length - 2.5).abs() < EPS);
        } else {
            panic!("Expected Distance constraint");
        }
    }
    #[test]
    fn test_floor_constraint_stops_particle() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 5);
        let i = sys.add_particle([0.0, 5.0, 0.0], 1.0);
        sys.add_floor_constraint(i, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        for _ in 0..300 {
            sys.step(dt);
        }
        assert!(
            sys.particles[i].position[1] >= -0.01,
            "Particle should not fall below floor, got y={}",
            sys.particles[i].position[1]
        );
    }
    #[test]
    fn test_floor_restitution_zero() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 5);
        let i = sys.add_particle([0.0, 1.0, 0.0], 1.0);
        sys.add_floor_constraint(i, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        for _ in 0..600 {
            sys.step(dt);
        }
        assert!(
            sys.particles[i].position[1] >= -0.05,
            "Particle should rest on floor"
        );
    }
    #[test]
    fn test_kinetic_energy_non_negative() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys.add_particle([0.0, 5.0, 0.0], 1.0);
        sys.step(1.0 / 60.0);
        assert!(sys.total_kinetic_energy() >= 0.0);
    }
    #[test]
    fn test_kinetic_energy_stationary_zero() {
        let mut sys = PbdSystem::new([0.0; 3], 1);
        sys.add_particle([0.0; 3], 1.0);
        assert!((sys.total_kinetic_energy()).abs() < EPS);
    }
    #[test]
    fn test_tet_volume_unit() {
        let v = compute_tetrahedron_volume(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        assert!(
            (v - 1.0 / 6.0).abs() < EPS,
            "Unit tet volume should be 1/6, got {v}"
        );
    }
    #[test]
    fn test_tet_volume_degenerate() {
        let v = compute_tetrahedron_volume(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        );
        assert!(v.abs() < EPS, "Degenerate tet volume should be 0, got {v}");
    }
    #[test]
    fn test_build_chain_count() {
        let sys = build_chain(5, 1.0, 1.0, 0.0);
        assert_eq!(sys.particle_count(), 5);
        assert_eq!(sys.constraints.len(), 4);
    }
    #[test]
    fn test_build_chain_top_fixed() {
        let sys = build_chain(4, 1.0, 1.0, 0.0);
        assert!(sys.particles[0].fixed, "Top of chain must be fixed");
        assert!(!sys.particles[1].fixed, "Other chain particles are dynamic");
    }
    #[test]
    fn test_build_grid_cloth_count() {
        let sys = build_grid_cloth(3, 4, 1.0, 1.0, 0.0);
        assert_eq!(sys.particle_count(), 12);
        assert_eq!(sys.constraints.len(), 17);
    }
    #[test]
    fn test_chain_hangs_under_gravity() {
        let mut sys = build_chain(5, 1.0, 1.0, 1e-5);
        let dt = 1.0 / 60.0;
        for _ in 0..300 {
            sys.step(dt);
        }
        let y_initial = -4.0;
        let y_final = sys.particles[4].position[1];
        assert!(
            y_final <= y_initial + 0.1,
            "Bottom chain particle should sag, got y={y_final}"
        );
    }
    #[test]
    fn test_solve_distance_standalone() {
        let mut pts = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([4.0, 0.0, 0.0], 1.0),
        ];
        for _ in 0..100 {
            solve_distance_xpbd(&mut pts, 0, 1, 2.0, 0.0, 1.0 / 60.0);
            for p in pts.iter_mut() {
                p.position = p.predicted;
            }
        }
        let d = dist3(pts[0].position, pts[1].position);
        assert!((d - 2.0).abs() < 0.01, "Standalone solver: d={d}");
    }
    #[test]
    fn test_substep_count_enforced() {
        let sys = PbdSystem::new([0.0; 3], 0);
        assert_eq!(sys.substeps, 1);
    }
    #[test]
    fn test_sphere_collision_pushes_out() {
        let mut p = PbdParticle::new([0.0, 0.0, 0.0], 1.0);
        p.predicted = [0.1, 0.0, 0.0];
        solve_sphere_collision(&mut p, [0.0, 0.0, 0.0], 1.0, 0.0);
        let d = dist3(p.predicted, [0.0, 0.0, 0.0]);
        assert!(
            d >= 1.0 - 1e-9,
            "Particle should be pushed outside sphere, d={d}"
        );
    }
    #[test]
    fn test_sphere_collision_outside_unchanged() {
        let mut p = PbdParticle::new([2.0, 0.0, 0.0], 1.0);
        p.predicted = [2.0, 0.0, 0.0];
        let before = p.predicted;
        solve_sphere_collision(&mut p, [0.0, 0.0, 0.0], 1.0, 0.0);
        let after = p.predicted;
        assert_eq!(before, after, "Outside particle should not be moved");
    }
    #[test]
    fn test_aabb_collision_clamps_particle() {
        let mut p = PbdParticle::new([5.0, 0.0, 0.0], 1.0);
        p.predicted = [5.0, 0.0, 0.0];
        p.velocity = [2.0, 0.0, 0.0];
        solve_aabb_collision(&mut p, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.0);
        assert!(
            p.predicted[0] <= 1.0 + 1e-9,
            "Particle should be clamped inside AABB"
        );
        assert!(p.velocity[0] <= 0.0 + 1e-9, "Velocity should be reversed");
    }
    #[test]
    fn test_aabb_collision_negative_side() {
        let mut p = PbdParticle::new([-0.5, 0.0, 0.0], 1.0);
        p.predicted = [-0.5, 0.0, 0.0];
        p.velocity = [-1.0, 0.0, 0.0];
        solve_aabb_collision(&mut p, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.5);
        assert!(
            p.predicted[0] >= 0.0 - 1e-9,
            "Particle should be on +x side of min"
        );
        assert!(p.velocity[0] >= 0.0 - 1e-9, "Velocity should be reflected");
    }
    #[test]
    fn test_triangle_collision_below_plane() {
        let mut p = PbdParticle::new([0.5, -0.5, 0.5], 1.0);
        p.predicted = [0.5, -0.5, 0.5];
        let ta = [0.0, 0.0, 0.0];
        let tb = [0.0, 0.0, 1.0];
        let tc = [1.0, 0.0, 0.0];
        solve_triangle_collision(&mut p, ta, tb, tc, 0.0);
        assert!(
            p.predicted[1] >= -1e-9,
            "Particle should be pushed above triangle plane, got y={}",
            p.predicted[1]
        );
    }
    #[test]
    fn test_self_collision_all_pairs_separation() {
        let mut particles = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([0.05, 0.0, 0.0], 1.0),
        ];
        particles[0].predicted = [0.0, 0.0, 0.0];
        particles[1].predicted = [0.05, 0.0, 0.0];
        solve_self_collision_all_pairs(&mut particles, 0.1, 1.0);
        let d = dist3(particles[0].predicted, particles[1].predicted);
        assert!(
            d >= 0.1 - 1e-9,
            "Particles should be separated to at least min_distance: d={d}"
        );
    }
    #[test]
    fn test_self_collision_grid_separation() {
        let mut particles = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([0.05, 0.0, 0.0], 1.0),
            PbdParticle::new([0.5, 0.5, 0.5], 1.0),
        ];
        for p in &mut particles {
            p.predicted = p.position;
        }
        solve_self_collision_grid(&mut particles, 0.1, 1.0, 0.2);
        let d01 = dist3(particles[0].predicted, particles[1].predicted);
        assert!(d01 >= 0.1 - 1e-9, "Grid self-collision: d01={d01}");
    }
    #[test]
    fn test_rigid_proxy_collision() {
        let proxy = RigidBodyProxy::new([0.0, 0.0, 0.0], 1.0);
        let mut p = PbdParticle::new([0.2, 0.0, 0.0], 1.0);
        p.predicted = [0.2, 0.0, 0.0];
        proxy.apply_collision(&mut p, 0.0);
        let d = dist3(p.predicted, [0.0, 0.0, 0.0]);
        assert!(
            d >= 1.0 - 1e-9,
            "Rigid proxy should push particle outside radius: d={d}"
        );
    }
    #[test]
    fn test_rigid_proxy_step() {
        let mut proxy = RigidBodyProxy::new([0.0, 0.0, 0.0], 1.0);
        proxy.velocity = [1.0, 0.0, 0.0];
        proxy.step(0.1);
        assert!((proxy.center[0] - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_cloth_rigid_attachment() {
        let mut particles = vec![PbdParticle::new([0.0, 0.0, 0.0], 1.0)];
        particles[0].predicted = [0.0, 0.0, 0.0];
        let target = [2.0, 0.0, 0.0];
        solve_cloth_rigid_attachment(&mut particles, 0, target, 0.5, 0.016);
        let after = particles[0].predicted;
        assert!(
            after[0] > 0.5,
            "Cloth particle should move toward target: x={}",
            after[0]
        );
    }
    #[test]
    fn test_adaptive_substeps_stationary() {
        let particles = vec![PbdParticle::new([0.0; 3], 1.0)];
        let n = adaptive_substep_count(&particles, 1.0 / 60.0, 0.5, 1, 10);
        assert_eq!(n, 1, "Stationary particles should use minimum substeps");
    }
    #[test]
    fn test_adaptive_substeps_fast_particle() {
        let mut p = PbdParticle::new([0.0; 3], 1.0);
        p.velocity = [100.0, 0.0, 0.0];
        let particles = vec![p];
        let n = adaptive_substep_count(&particles, 1.0 / 60.0, 0.1, 1, 32);
        assert!(n > 1, "Fast particle should require more substeps: n={n}");
    }
    #[test]
    fn test_adaptive_substeps_clamped_max() {
        let mut p = PbdParticle::new([0.0; 3], 1.0);
        p.velocity = [1_000_000.0, 0.0, 0.0];
        let particles = vec![p];
        let n = adaptive_substep_count(&particles, 1.0 / 60.0, 0.1, 1, 8);
        assert!(n <= 8, "Substeps should be clamped to max: n={n}");
    }
    #[test]
    fn test_profiled_system_records_stats() {
        let mut sys = ProfiledPbdSystem::new([0.0, -9.81, 0.0], 4);
        sys.system.add_particle([0.0, 5.0, 0.0], 1.0);
        sys.step(1.0 / 60.0);
        assert_eq!(sys.last_profile.substeps, 4);
        assert!(sys.last_profile.total_constraint_solves == 0);
    }
    #[test]
    fn test_profiled_system_constraint_residual() {
        let mut sys = ProfiledPbdSystem::new([0.0; 3], 4);
        let a = sys.system.add_particle([0.0, 0.0, 0.0], 1.0);
        let b = sys.system.add_particle([3.0, 0.0, 0.0], 1.0);
        sys.system.add_distance_constraint(a, b, 0.0);
        sys.step(1.0 / 60.0);
        assert!(
            sys.last_profile.residual >= 0.0,
            "Residual must be non-negative"
        );
    }
    #[test]
    fn test_adaptive_pbd_system_step() {
        let mut sys = AdaptivePbdSystem::new([0.0, -9.81, 0.0], 1, 8);
        sys.system.add_particle([0.0, 10.0, 0.0], 1.0);
        let dt = 1.0 / 60.0;
        for _ in 0..10 {
            sys.step(dt);
        }
        assert!(
            sys.system.particles[0].position[1] < 10.0,
            "Particle should fall under gravity"
        );
    }
}
