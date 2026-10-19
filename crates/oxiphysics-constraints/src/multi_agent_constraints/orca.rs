// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Optimal Reciprocal Collision Avoidance (ORCA) for disc-shaped agents.
//!
//! A faithful implementation of the agent-agent part of RVO2
//! (van den Berg, Guy, Lin, Manocha, "Reciprocal n-body collision avoidance",
//! 2011): every neighbour contributes one half-plane of permitted velocities,
//! and the new velocity is the permitted velocity closest to the preferred one
//! inside the maximum-speed disc, found with RVO2's incremental 2-D linear
//! programs (`linearProgram1/2/3`).
//!
//! # The half-plane of one neighbour
//!
//! With `p = p_B − p_A`, `v = v_A − v_B`, `R = r_A + r_B`, the *truncated
//! velocity obstacle* `VO^τ` is the set of relative velocities that cause a
//! collision within the time horizon `τ`.  Its boundary consists of the two
//! tangent legs of the cone from the origin to the disc `D(p, R)` and the
//! circular cut-off `D(p/τ, R/τ)`.  `u` is the shortest vector from `v` to the
//! boundary of `VO^τ`; each agent takes half of the responsibility for avoiding
//! the obstacle, so agent A is confined to the half-plane through `v_A + u/2`
//! whose normal is the outward direction of the boundary at the closest point.
//! Any pair of velocities taken from the two reciprocal half-planes is free of
//! collisions for at least `τ`.
//!
//! When the agents already overlap the same construction is applied to the
//! cut-off circle of the *time step*, `D(p/Δt, R/Δt)`, which makes the agents
//! separate within one step instead of waiting for the horizon.

use super::{vec2_add, vec2_dot, vec2_len, vec2_normalize, vec2_scale, vec2_sub};

/// Tolerance below which two unit directions count as parallel.
const PARALLEL_TOLERANCE: f64 = 1e-9;

/// A 2D half-plane constraint for ORCA (Optimal Reciprocal Collision Avoidance).
///
/// The permitted velocities `v` satisfy `n · (v − p) ≥ 0`, where `n` is the unit
/// normal pointing into the permitted side and `p` is a point on the boundary.
#[derive(Debug, Clone)]
pub struct HalfPlane2D {
    /// Point on the boundary of the half-plane.
    pub point: [f64; 2],
    /// Unit normal of the half-plane, pointing into the permitted side.
    pub normal: [f64; 2],
}

/// A disc-shaped agent as seen by its neighbours.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrcaDisc {
    /// Centre of the disc.
    pub position: [f64; 2],
    /// Current velocity.
    pub velocity: [f64; 2],
    /// Radius of the disc.
    pub radius: f64,
}

/// A directed line of RVO2; the permitted side is on its left.
#[derive(Debug, Clone, Copy)]
struct Line {
    point: [f64; 2],
    direction: [f64; 2],
}

fn det(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn line_from_half_plane(half_plane: &HalfPlane2D) -> Line {
    Line {
        point: half_plane.point,
        // The permitted side n·(v − p) ≥ 0 is the left side of the direction
        // d = (n_y, −n_x), since det(d, w) = n·w.
        direction: [half_plane.normal[1], -half_plane.normal[0]],
    }
}

/// Compute the ORCA half-plane for agent A (radius `r_a`, position `p_a`,
/// velocity `v_a`) with respect to agent B (radius `r_b`, position `p_b`,
/// velocity `v_b`) over a time horizon `tau`, for a simulation step `dt`.
///
/// This is the per-neighbour body of RVO2's `computeNewVelocity`.  The agent
/// should choose a new velocity inside the returned half-plane (and inside the
/// half-planes of all other neighbours).
///
/// * Agents that do not overlap get the half-plane that bounds the truncated
///   velocity obstacle `VO^τ` (cut-off circle or tangent leg, whichever is
///   nearest to the relative velocity).  The half-plane is returned even when
///   the relative velocity is already outside the obstacle: it is then the
///   limit up to which the velocity may still change without entering it.
/// * Agents that overlap are separated within `dt`: the normal points away from
///   the other agent, so the constraint never pushes the agents together.
///
/// `tau` and `dt` are floored at `f64::EPSILON`.  Coincident agents with a
/// zero relative velocity have no preferred separation direction; the
/// half-plane then pushes along +x.
pub fn orca_half_plane(
    p_a: [f64; 2],
    v_a: [f64; 2],
    r_a: f64,
    p_b: [f64; 2],
    v_b: [f64; 2],
    r_b: f64,
    tau: f64,
    dt: f64,
) -> HalfPlane2D {
    let relative_position = vec2_sub(p_b, p_a);
    let relative_velocity = vec2_sub(v_a, v_b);
    let dist_sq = vec2_dot(relative_position, relative_position);
    let combined_radius = r_a + r_b;
    let combined_radius_sq = combined_radius * combined_radius;

    // `direction` is the direction of the boundary line (the permitted side is
    // on its left) and `u` the smallest change of the relative velocity that
    // reaches the boundary of the obstacle.
    let (direction, u) = if dist_sq > combined_radius_sq {
        // No collision.
        let inv_tau = 1.0 / tau.max(f64::EPSILON);
        // Vector from the cut-off centre to the relative velocity.
        let w = vec2_sub(relative_velocity, vec2_scale(relative_position, inv_tau));
        let w_length_sq = vec2_dot(w, w);
        let dot_product_1 = vec2_dot(w, relative_position);
        if dot_product_1 < 0.0 && dot_product_1 * dot_product_1 > combined_radius_sq * w_length_sq {
            // Project on the cut-off circle.
            let w_length = w_length_sq.sqrt();
            let unit_w = vec2_scale(w, 1.0 / w_length);
            (
                [unit_w[1], -unit_w[0]],
                vec2_scale(unit_w, combined_radius * inv_tau - w_length),
            )
        } else {
            // Project on a leg.
            let leg = (dist_sq - combined_radius_sq).sqrt();
            let [px, py] = relative_position;
            let leg_direction = if det(relative_position, w) > 0.0 {
                // Left leg.
                [
                    (px * leg - py * combined_radius) / dist_sq,
                    (px * combined_radius + py * leg) / dist_sq,
                ]
            } else {
                // Right leg.
                [
                    -(px * leg + py * combined_radius) / dist_sq,
                    -(-px * combined_radius + py * leg) / dist_sq,
                ]
            };
            let dot_product_2 = vec2_dot(relative_velocity, leg_direction);
            (
                leg_direction,
                vec2_sub(vec2_scale(leg_direction, dot_product_2), relative_velocity),
            )
        }
    } else {
        // Collision: project on the cut-off circle of the time step.
        let inv_dt = 1.0 / dt.max(f64::EPSILON);
        let w = vec2_sub(relative_velocity, vec2_scale(relative_position, inv_dt));
        let w_length = vec2_len(w);
        let unit_w = if w_length > 1e-12 {
            vec2_scale(w, 1.0 / w_length)
        } else {
            // The relative velocity sits exactly on the cut-off centre: separate
            // along the line of centres, away from the other agent.
            let away = vec2_normalize(vec2_scale(relative_position, -1.0));
            if vec2_len(away) > 0.0 {
                away
            } else {
                [1.0, 0.0]
            }
        };
        (
            [unit_w[1], -unit_w[0]],
            vec2_scale(unit_w, combined_radius * inv_dt - w_length),
        )
    };

    HalfPlane2D {
        point: vec2_add(v_a, vec2_scale(u, 0.5)),
        normal: [-direction[1], direction[0]],
    }
}

/// Largest-speed circle feasibility test and 1-D optimisation along line
/// `line_no` (RVO2 `linearProgram1`).  Returns the optimum on the line or
/// `None` when the line is infeasible given the lines before it.
fn linear_program_1(
    lines: &[Line],
    line_no: usize,
    radius: f64,
    opt_velocity: [f64; 2],
    direction_opt: bool,
) -> Option<[f64; 2]> {
    let line = lines[line_no];
    let dot_product = vec2_dot(line.point, line.direction);
    let discriminant =
        dot_product * dot_product + radius * radius - vec2_dot(line.point, line.point);
    if discriminant < 0.0 {
        // The maximum-speed circle fully invalidates this line.
        return None;
    }
    let sqrt_discriminant = discriminant.sqrt();
    let mut t_left = -dot_product - sqrt_discriminant;
    let mut t_right = -dot_product + sqrt_discriminant;

    for earlier in &lines[..line_no] {
        let denominator = det(line.direction, earlier.direction);
        let numerator = det(earlier.direction, vec2_sub(line.point, earlier.point));
        if denominator.abs() <= PARALLEL_TOLERANCE {
            // The lines are (almost) parallel.
            if numerator < 0.0 {
                return None;
            }
            continue;
        }
        let t = numerator / denominator;
        if denominator >= 0.0 {
            // The earlier line bounds this line on the right.
            t_right = t_right.min(t);
        } else {
            // The earlier line bounds this line on the left.
            t_left = t_left.max(t);
        }
        if t_left > t_right {
            return None;
        }
    }

    let along = |t: f64| vec2_add(line.point, vec2_scale(line.direction, t));
    if direction_opt {
        // Optimise the direction.
        if vec2_dot(opt_velocity, line.direction) > 0.0 {
            Some(along(t_right))
        } else {
            Some(along(t_left))
        }
    } else {
        // Optimise the closest point.
        let t = vec2_dot(line.direction, vec2_sub(opt_velocity, line.point));
        // Explicit comparisons rather than `f64::clamp`, which panics on NaN bounds.
        if t < t_left {
            Some(along(t_left))
        } else if t > t_right {
            Some(along(t_right))
        } else {
            Some(along(t))
        }
    }
}

/// Incremental 2-D linear program (RVO2 `linearProgram2`): the point inside the
/// disc of radius `radius` that satisfies every line and is closest to
/// `opt_velocity` (or extreme in direction `opt_velocity`).  Returns the index
/// of the first line that cannot be satisfied (`lines.len()` when all can) and
/// the best point found before it.
fn linear_program_2(
    lines: &[Line],
    radius: f64,
    opt_velocity: [f64; 2],
    direction_opt: bool,
) -> (usize, [f64; 2]) {
    let mut result = if direction_opt {
        // The optimisation velocity has unit length in this case.
        vec2_scale(opt_velocity, radius)
    } else if vec2_dot(opt_velocity, opt_velocity) > radius * radius {
        vec2_scale(vec2_normalize(opt_velocity), radius)
    } else {
        opt_velocity
    };
    for (index, line) in lines.iter().enumerate() {
        if det(line.direction, vec2_sub(line.point, result)) > 0.0 {
            // The result violates line `index`: re-optimise on that line.
            match linear_program_1(lines, index, radius, opt_velocity, direction_opt) {
                Some(on_line) => result = on_line,
                None => return (index, result),
            }
        }
    }
    (lines.len(), result)
}

/// Fallback when the lines are infeasible (RVO2 `linearProgram3`): find the
/// point inside the disc that minimises the largest violation, continuing from
/// the partial solution `result` that satisfied the lines before `begin_line`.
fn linear_program_3(lines: &[Line], begin_line: usize, radius: f64, start: [f64; 2]) -> [f64; 2] {
    let mut result = start;
    let mut distance = 0.0;
    for (index, current) in lines.iter().enumerate().skip(begin_line) {
        if det(current.direction, vec2_sub(current.point, result)) > distance {
            // The result does not satisfy the current line: project the earlier
            // lines onto it and optimise the direction along the boundary.
            let mut projected: Vec<Line> = Vec::with_capacity(index);
            for earlier in &lines[..index] {
                let determinant = det(current.direction, earlier.direction);
                let point = if determinant.abs() <= PARALLEL_TOLERANCE {
                    // The current line and the earlier line are parallel.
                    if vec2_dot(current.direction, earlier.direction) > 0.0 {
                        // Same direction: the earlier line adds nothing.
                        continue;
                    }
                    // Opposite directions: the boundary halfway between them.
                    vec2_scale(vec2_add(current.point, earlier.point), 0.5)
                } else {
                    let offset = det(earlier.direction, vec2_sub(current.point, earlier.point))
                        / determinant;
                    vec2_add(current.point, vec2_scale(current.direction, offset))
                };
                projected.push(Line {
                    point,
                    direction: vec2_normalize(vec2_sub(earlier.direction, current.direction)),
                });
            }
            let optimisation_direction = [-current.direction[1], current.direction[0]];
            let (failed, candidate) =
                linear_program_2(&projected, radius, optimisation_direction, true);
            // A failure cannot happen in principle (the result is by definition
            // feasible for this program); it is floating point noise, so the
            // current result is kept.
            if failed >= projected.len() {
                result = candidate;
            }
            distance = det(current.direction, vec2_sub(current.point, result));
        }
    }
    result
}

/// Find the velocity closest to `v_pref` that satisfies all ORCA half-plane
/// `constraints` and has a speed of at most `max_speed`.
///
/// This is RVO2's `linearProgram2` followed, when the constraints are mutually
/// infeasible within the speed limit, by `linearProgram3`, which then returns
/// the velocity that minimises the largest violation of any constraint.  Always
/// returns a velocity inside the maximum-speed disc.
pub fn orca_solve_velocity(
    v_pref: [f64; 2],
    max_speed: f64,
    constraints: &[HalfPlane2D],
) -> [f64; 2] {
    let max_speed = max_speed.max(0.0);
    let lines: Vec<Line> = constraints.iter().map(line_from_half_plane).collect();
    let (failed, result) = linear_program_2(&lines, max_speed, v_pref, false);
    if failed < lines.len() {
        linear_program_3(&lines, failed, max_speed, result)
    } else {
        result
    }
}

/// Compute the new velocity of `agent` from its preferred velocity and its
/// `neighbors` (RVO2 `Agent::computeNewVelocity` without static obstacles).
///
/// * `preferred_velocity` – the velocity the agent would take if alone.
/// * `max_speed` – the speed limit of the agent.
/// * `time_horizon` – look-ahead time `τ` for the velocity obstacles.
/// * `time_step` – the simulation step `Δt`, used to separate overlapping agents.
///
/// All agents must be updated from the same snapshot of positions and
/// velocities (compute every new velocity first, then move the agents) for the
/// reciprocal guarantee to hold.
pub fn orca_compute_new_velocity(
    agent: &OrcaDisc,
    preferred_velocity: [f64; 2],
    max_speed: f64,
    time_horizon: f64,
    time_step: f64,
    neighbors: &[OrcaDisc],
) -> [f64; 2] {
    let half_planes: Vec<HalfPlane2D> = neighbors
        .iter()
        .map(|other| {
            orca_half_plane(
                agent.position,
                agent.velocity,
                agent.radius,
                other.position,
                other.velocity,
                other.radius,
                time_horizon,
                time_step,
            )
        })
        .collect();
    orca_solve_velocity(preferred_velocity, max_speed, &half_planes)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TAU: f64 = 5.0;
    const DT: f64 = 0.1;

    /// Deterministic xorshift generator in `[0, 1)` for the sampling tests.
    struct Xorshift(u64);

    impl Xorshift {
        fn next_unit(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        fn range(&mut self, low: f64, high: f64) -> f64 {
            low + (high - low) * self.next_unit()
        }
    }

    fn satisfies(half_plane: &HalfPlane2D, velocity: [f64; 2]) -> bool {
        vec2_dot(half_plane.normal, vec2_sub(velocity, half_plane.point)) >= -1e-9
    }

    /// Smallest distance between two agents moving with constant velocities
    /// during `[0, horizon]`.
    fn closest_approach(
        relative_position: [f64; 2],
        relative_velocity: [f64; 2],
        horizon: f64,
    ) -> f64 {
        let speed_sq = vec2_dot(relative_velocity, relative_velocity);
        let t = if speed_sq > 1e-18 {
            (vec2_dot(relative_position, relative_velocity) / speed_sq).clamp(0.0, horizon)
        } else {
            0.0
        };
        vec2_len(vec2_sub(
            relative_position,
            vec2_scale(relative_velocity, t),
        ))
    }

    #[test]
    fn head_on_cutoff_circle_matches_a_hand_computation() {
        // A at (-5, 0) moving +1, B at (5, 0) moving -1, radius 0.5 each, tau 2.
        // p = (10, 0), v = (2, 0), w = v - p/tau = (-3, 0): the relative
        // velocity lies before the cut-off circle (|w| = 3 > R/tau = 0.5), so
        // u = (R/tau - |w|)·w/|w| = (2.5, 0) and the line passes through
        // v_A + u/2 = (2.25, 0).
        let half_plane = orca_half_plane(
            [-5.0, 0.0],
            [1.0, 0.0],
            0.5,
            [5.0, 0.0],
            [-1.0, 0.0],
            0.5,
            2.0,
            DT,
        );
        assert!((half_plane.point[0] - 2.25).abs() < 1e-12);
        assert!(half_plane.point[1].abs() < 1e-12);
        assert!((half_plane.normal[0] + 1.0).abs() < 1e-12);
        assert!(half_plane.normal[1].abs() < 1e-12);
        // "vx <= 2.25": A may speed up to 2.25 m/s before the pair would hit
        // within the horizon; the preferred 1 m/s is allowed.
        assert!(satisfies(&half_plane, [1.0, 0.0]));
        assert!(!satisfies(&half_plane, [3.0, 0.0]));
    }

    #[test]
    fn leg_projection_matches_a_hand_computation() {
        // A at the origin, B straight ahead at (4, 0) with radius 1, B at rest,
        // A moving at (2, 0.1) is heading into the obstacle (tau = 5).
        // |p|^2 = 16, R = 1, leg = sqrt(15).  w = v - p/tau = (1.2, 0.1);
        // dot(w, p) > 0, so the relative velocity is projected on a leg.
        // det(p, w) = 4*0.1 - 0*1.2 = 0.4 > 0 selects the left leg with
        // direction (leg, R)/|p| = (sqrt(15)/4, 1/4).
        let half_plane = orca_half_plane(
            [0.0, 0.0],
            [2.0, 0.1],
            0.5,
            [4.0, 0.0],
            [0.0, 0.0],
            0.5,
            5.0,
            DT,
        );
        let leg = 15.0_f64.sqrt();
        let direction = [leg / 4.0, 0.25];
        let relative_velocity = [2.0, 0.1];
        let along = vec2_dot(relative_velocity, direction);
        let u = [
            direction[0] * along - relative_velocity[0],
            direction[1] * along - relative_velocity[1],
        ];
        assert!((half_plane.point[0] - (2.0 + 0.5 * u[0])).abs() < 1e-12);
        assert!((half_plane.point[1] - (0.1 + 0.5 * u[1])).abs() < 1e-12);
        // The normal is the left normal of the leg direction.
        assert!((half_plane.normal[0] + direction[1]).abs() < 1e-12);
        assert!((half_plane.normal[1] - direction[0]).abs() < 1e-12);
        assert!((vec2_len(half_plane.normal) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn far_away_neighbour_does_not_constrain_the_current_velocity() {
        let half_plane = orca_half_plane(
            [0.0, 0.0],
            [1.0, 0.0],
            0.3,
            [10.0, 5.0],
            [0.0, 0.0],
            0.3,
            2.0,
            DT,
        );
        assert!(satisfies(&half_plane, [1.0, 0.0]));
        assert!((vec2_len(half_plane.normal) - 1.0).abs() < 1e-12);
        // Within the maximum speed the agent keeps all of its options.
        for velocity in [[1.5, 0.0], [0.0, 1.5], [-1.5, 0.0], [1.0, 1.0]] {
            assert!(
                satisfies(&half_plane, velocity),
                "velocity {velocity:?} should be free"
            );
        }
    }

    #[test]
    fn velocities_from_reciprocal_half_planes_never_collide_within_the_horizon() {
        let mut rng = Xorshift(0x9E37_79B9_7F4A_7C15);
        let mut checked = 0;
        for _ in 0..400 {
            let p_a = [rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)];
            let angle = rng.range(0.0, std::f64::consts::TAU);
            let distance = rng.range(1.2, 6.0);
            let p_b = [
                p_a[0] + distance * angle.cos(),
                p_a[1] + distance * angle.sin(),
            ];
            let v_a = [rng.range(-1.5, 1.5), rng.range(-1.5, 1.5)];
            let v_b = [rng.range(-1.5, 1.5), rng.range(-1.5, 1.5)];
            let (r_a, r_b) = (0.5, 0.5);
            let hp_a = orca_half_plane(p_a, v_a, r_a, p_b, v_b, r_b, TAU, DT);
            let hp_b = orca_half_plane(p_b, v_b, r_b, p_a, v_a, r_a, TAU, DT);
            for _ in 0..40 {
                let new_a = [rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)];
                let new_b = [rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)];
                if !satisfies(&hp_a, new_a) || !satisfies(&hp_b, new_b) {
                    continue;
                }
                checked += 1;
                let gap = closest_approach(vec2_sub(p_b, p_a), vec2_sub(new_a, new_b), TAU);
                assert!(
                    gap >= r_a + r_b - 1e-9,
                    "collision within the horizon: gap {gap}, p_a {p_a:?} v_a {v_a:?} \
                     p_b {p_b:?} v_b {v_b:?} -> {new_a:?} / {new_b:?}"
                );
            }
        }
        assert!(
            checked > 1000,
            "the sampling must exercise many pairs: {checked}"
        );
    }

    #[test]
    fn overlapping_agents_are_separated_never_pushed_together() {
        // A at the origin, B 0.8 m ahead, combined radius 1.0: they overlap.
        let p_a = [0.0, 0.0];
        let p_b = [0.8, 0.0];
        let half_plane = orca_half_plane(p_a, [0.0, 0.0], 0.5, p_b, [0.0, 0.0], 0.5, TAU, DT);
        // The normal points away from B, so the constraint is "v_x <= -(overlap)/(2 dt)".
        assert!(
            half_plane.normal[0] < -0.99,
            "normal {:?}",
            half_plane.normal
        );
        // Heading for B is forbidden, retreating is allowed.
        assert!(!satisfies(&half_plane, [1.0, 0.0]));
        assert!(satisfies(&half_plane, [-2.0, 0.0]));
        // After a step with any velocity pair from the two half-planes the
        // agents no longer overlap.
        let mirrored = orca_half_plane(p_b, [0.0, 0.0], 0.5, p_a, [0.0, 0.0], 0.5, TAU, DT);
        let new_a = [-2.0, 0.0];
        let new_b = [2.0, 0.0];
        assert!(satisfies(&half_plane, new_a) && satisfies(&mirrored, new_b));
        let gap = vec2_len(vec2_sub(
            vec2_add(p_b, vec2_scale(new_b, DT)),
            vec2_add(p_a, vec2_scale(new_a, DT)),
        ));
        assert!(
            gap >= 1.0 - 1e-9,
            "agents still overlap after a step: {gap}"
        );
    }

    #[test]
    fn overlapping_agent_with_a_pushing_preferred_velocity_retreats() {
        let agent = OrcaDisc {
            position: [0.0, 0.0],
            velocity: [0.0, 0.0],
            radius: 0.5,
        };
        let other = OrcaDisc {
            position: [0.8, 0.0],
            velocity: [0.0, 0.0],
            radius: 0.5,
        };
        let new_velocity = orca_compute_new_velocity(&agent, [1.0, 0.0], 3.0, TAU, DT, &[other]);
        assert!(
            new_velocity[0] <= -0.1 / DT + 1e-9,
            "the agent must back off, got {new_velocity:?}"
        );
    }

    #[test]
    fn coincident_agents_still_get_a_finite_half_plane() {
        let half_plane = orca_half_plane(
            [1.0, 1.0],
            [0.0, 0.0],
            0.5,
            [1.0, 1.0],
            [0.0, 0.0],
            0.5,
            TAU,
            DT,
        );
        assert!(half_plane.point.iter().all(|c| c.is_finite()));
        assert!((vec2_len(half_plane.normal) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn solve_without_constraints_keeps_the_preferred_velocity_inside_the_disc() {
        assert_eq!(orca_solve_velocity([0.5, 0.25], 2.0, &[]), [0.5, 0.25]);
        let clipped = orca_solve_velocity([6.0, 8.0], 2.0, &[]);
        assert!((clipped[0] - 1.2).abs() < 1e-12 && (clipped[1] - 1.6).abs() < 1e-12);
    }

    #[test]
    fn solve_projects_onto_a_violated_half_plane() {
        // Permitted: v_x <= 1 (normal (-1, 0) through (1, 0)).
        let wall = HalfPlane2D {
            point: [1.0, 0.0],
            normal: [-1.0, 0.0],
        };
        let result = orca_solve_velocity([3.0, 0.5], 5.0, &[wall]);
        assert!((result[0] - 1.0).abs() < 1e-12 && (result[1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn solve_finds_the_vertex_of_two_half_planes() {
        // v_x <= 1 and v_y <= 2: the closest permitted point to (5, 5) is (1, 2).
        let planes = [
            HalfPlane2D {
                point: [1.0, 0.0],
                normal: [-1.0, 0.0],
            },
            HalfPlane2D {
                point: [0.0, 2.0],
                normal: [0.0, -1.0],
            },
        ];
        let result = orca_solve_velocity([5.0, 5.0], 10.0, &planes);
        assert!((result[0] - 1.0).abs() < 1e-12 && (result[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn solve_respects_the_speed_limit_when_constraints_cannot_be_met() {
        // v_x >= 5 is out of reach at max speed 2: the closest point of the
        // disc to the half-plane is (2, 0).
        let demand = HalfPlane2D {
            point: [5.0, 0.0],
            normal: [1.0, 0.0],
        };
        let result = orca_solve_velocity([0.0, 0.0], 2.0, &[demand]);
        assert!(vec2_len(result) <= 2.0 + 1e-9);
        assert!((result[0] - 2.0).abs() < 1e-9 && result[1].abs() < 1e-9);
    }

    #[test]
    fn solve_minimises_the_largest_violation_of_infeasible_constraints() {
        // v_x >= 1 and v_x <= -1 cannot both hold; the least-bad answer splits
        // the violation evenly at v_x = 0.
        let planes = [
            HalfPlane2D {
                point: [1.0, 0.0],
                normal: [1.0, 0.0],
            },
            HalfPlane2D {
                point: [-1.0, 0.0],
                normal: [-1.0, 0.0],
            },
        ];
        let result = orca_solve_velocity([0.3, 0.7], 3.0, &planes);
        assert!(result[0].abs() < 1e-9, "got {result:?}");
        assert!(vec2_len(result) <= 3.0 + 1e-9);
    }

    #[test]
    fn solved_velocity_satisfies_feasible_constraints_and_stays_closest() {
        let mut rng = Xorshift(0xD1B5_4A32_D192_ED03);
        for _ in 0..300 {
            let count = 1 + (rng.next_unit() * 5.0) as usize;
            let planes: Vec<HalfPlane2D> = (0..count)
                .map(|_| {
                    let angle = rng.range(0.0, std::f64::consts::TAU);
                    HalfPlane2D {
                        point: [rng.range(-0.5, 0.5), rng.range(-0.5, 0.5)],
                        normal: [angle.cos(), angle.sin()],
                    }
                })
                .collect();
            let preferred = [rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)];
            let result = orca_solve_velocity(preferred, 2.0, &planes);
            assert!(vec2_len(result) <= 2.0 + 1e-9);
            if planes.iter().all(|plane| satisfies(plane, result)) {
                // Brute force: no feasible sample is closer to the preferred velocity.
                let best = vec2_len(vec2_sub(result, preferred));
                for _ in 0..200 {
                    let sample = [rng.range(-2.0, 2.0), rng.range(-2.0, 2.0)];
                    if vec2_len(sample) <= 2.0
                        && planes.iter().all(|plane| satisfies(plane, sample))
                    {
                        assert!(
                            vec2_len(vec2_sub(sample, preferred)) >= best - 1e-9,
                            "sample {sample:?} beats the LP answer {result:?}"
                        );
                    }
                }
            }
        }
    }

    /// One simultaneous ORCA step for a set of discs with the given preferred
    /// velocities.
    fn step_agents(agents: &mut [OrcaDisc], preferred: &[[f64; 2]], max_speed: f64) {
        let snapshot: Vec<OrcaDisc> = agents.to_vec();
        for (index, agent) in agents.iter_mut().enumerate() {
            let neighbors: Vec<OrcaDisc> = snapshot
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != index)
                .map(|(_, disc)| *disc)
                .collect();
            agent.velocity = orca_compute_new_velocity(
                &snapshot[index],
                preferred[index],
                max_speed,
                TAU,
                DT,
                &neighbors,
            );
        }
        for agent in agents.iter_mut() {
            agent.position = vec2_add(agent.position, vec2_scale(agent.velocity, DT));
        }
    }

    fn minimum_pair_gap(agents: &[OrcaDisc]) -> f64 {
        let mut gap = f64::MAX;
        for (i, first) in agents.iter().enumerate() {
            for second in &agents[i + 1..] {
                let centre_distance = vec2_len(vec2_sub(first.position, second.position));
                gap = gap.min(centre_distance - first.radius - second.radius);
            }
        }
        gap
    }

    #[test]
    fn exactly_head_on_agents_never_overlap_and_slow_down_instead() {
        // With perfect symmetry the relative velocity lies on the axis of the
        // velocity obstacle, so ORCA keeps both agents on the axis and brakes
        // them: they creep towards each other but never touch.
        let mut agents = [
            OrcaDisc {
                position: [-5.0, 0.0],
                velocity: [0.0, 0.0],
                radius: 0.5,
            },
            OrcaDisc {
                position: [5.0, 0.0],
                velocity: [0.0, 0.0],
                radius: 0.5,
            },
        ];
        let preferred = [[1.0, 0.0], [-1.0, 0.0]];
        let mut smallest_gap = f64::MAX;
        for _ in 0..200 {
            step_agents(&mut agents, &preferred, 1.5);
            smallest_gap = smallest_gap.min(minimum_pair_gap(&agents));
        }
        assert!(
            smallest_gap >= -1e-6,
            "agents overlapped by {}",
            -smallest_gap
        );
        let separation = vec2_len(vec2_sub(agents[0].position, agents[1].position));
        assert!(separation < 3.0, "they must have closed in: {separation}");
        assert!(
            vec2_len(agents[0].velocity) < 0.5,
            "they must have braked: {:?}",
            agents[0].velocity
        );
    }

    #[test]
    fn head_on_agents_with_a_small_lateral_offset_pass_each_other() {
        for offset in [0.1, 0.3] {
            let mut agents = [
                OrcaDisc {
                    position: [-5.0, 0.0],
                    velocity: [0.0, 0.0],
                    radius: 0.5,
                },
                OrcaDisc {
                    position: [5.0, offset],
                    velocity: [0.0, 0.0],
                    radius: 0.5,
                },
            ];
            let preferred = [[1.0, 0.0], [-1.0, 0.0]];
            let mut smallest_gap = f64::MAX;
            for _ in 0..300 {
                step_agents(&mut agents, &preferred, 1.5);
                smallest_gap = smallest_gap.min(minimum_pair_gap(&agents));
            }
            assert!(
                smallest_gap >= -1e-6,
                "offset {offset}: agents overlapped by {}",
                -smallest_gap
            );
            assert!(
                agents[0].position[0] > agents[1].position[0] + 10.0,
                "offset {offset}: they must end up on each other's side: {:?} {:?}",
                agents[0].position,
                agents[1].position
            );
        }
    }

    /// Simulate agents walking to goals at unit speed and return the smallest
    /// gap between any two of them and the distances left to their goals.
    fn walk_to_goals(starts: &[[f64; 2]], goals: &[[f64; 2]], steps: usize) -> (f64, Vec<f64>) {
        let mut agents: Vec<OrcaDisc> = starts
            .iter()
            .map(|start| OrcaDisc {
                position: *start,
                velocity: [0.0, 0.0],
                radius: 0.5,
            })
            .collect();
        let mut smallest_gap = f64::MAX;
        for _ in 0..steps {
            let preferred: Vec<[f64; 2]> = agents
                .iter()
                .zip(goals)
                .map(|(agent, goal)| {
                    let to_goal = vec2_sub(*goal, agent.position);
                    let distance = vec2_len(to_goal);
                    if distance < 1e-9 {
                        [0.0, 0.0]
                    } else {
                        vec2_scale(to_goal, 1.0_f64.min(distance / DT) / distance)
                    }
                })
                .collect();
            step_agents(&mut agents, &preferred, 1.5);
            smallest_gap = smallest_gap.min(minimum_pair_gap(&agents));
        }
        let remaining = agents
            .iter()
            .zip(goals)
            .map(|(agent, goal)| vec2_len(vec2_sub(*goal, agent.position)))
            .collect();
        (smallest_gap, remaining)
    }

    #[test]
    fn asymmetric_ring_of_agents_reaches_antipodal_goals_without_collisions() {
        let mut starts = Vec::new();
        let mut goals = Vec::new();
        for index in 0..6 {
            let i = index as f64;
            let angle = std::f64::consts::TAU * i / 6.0 + 0.07 * i * i;
            let radius = 5.0 + 0.35 * i;
            let start = [radius * angle.cos(), radius * angle.sin()];
            starts.push(start);
            goals.push([-start[0], -start[1]]);
        }
        let (smallest_gap, remaining) = walk_to_goals(&starts, &goals, 900);
        assert!(
            smallest_gap >= -1e-6,
            "agents overlapped by {}",
            -smallest_gap
        );
        for left in remaining {
            assert!(left < 0.05, "an agent stopped {left} m from its goal");
        }
    }

    #[test]
    fn crossing_streams_of_agents_reach_their_goals_without_collisions() {
        let starts = [[-6.0, 0.0], [0.4, -6.0], [6.0, 0.9], [-0.5, 6.0]];
        let goals = [[6.0, 0.0], [0.4, 6.0], [-6.0, 0.9], [-0.5, -6.0]];
        let (smallest_gap, remaining) = walk_to_goals(&starts, &goals, 500);
        assert!(
            smallest_gap >= -1e-6,
            "agents overlapped by {}",
            -smallest_gap
        );
        for left in remaining {
            assert!(left < 0.05, "an agent stopped {left} m from its goal");
        }
    }
}
