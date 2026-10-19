//! Tests for the Bouzidi-Firdaouss-Lallemand interpolated bounce-back ([`InterpolatedBc`]).
//!
//! The first implementation used the coefficients `q` and `1 - q` (which do not reduce to the
//! half-way bounce-back at `q = 1/2`) and, for `q < 1/2`, sampled the wrong populations. These
//! tests pin the two BFL branches, their continuity at `q = 1/2` and the second-order accuracy of
//! a Poiseuille channel whose walls cut the lattice links at arbitrary distances.

use super::*;

fn c_x(i: usize) -> f64 {
    f64::from(D2Q9_VELOCITIES[i][0])
}

fn c_y(i: usize) -> f64 {
    f64::from(D2Q9_VELOCITIES[i][1])
}

#[test]
fn both_branches_reduce_to_halfway_bounce_back_at_q_one_half() {
    let bc = InterpolatedBc::new(0.5);
    for (f_i, upstream, opposite) in [(0.11, 0.9, 0.3), (0.0, 1.0, 1.0), (0.37, -0.4, 2.5)] {
        assert_eq!(bc.reflected_population(f_i, upstream, opposite), f_i);
    }
}

#[test]
fn reflected_population_is_continuous_across_q_one_half() {
    let (f_i, upstream, opposite) = (0.21, 0.13, 0.34);
    let below = InterpolatedBc::new(0.5 - 1e-9).reflected_population(f_i, upstream, opposite);
    let above = InterpolatedBc::new(0.5 + 1e-9).reflected_population(f_i, upstream, opposite);
    assert!((below - f_i).abs() < 1e-8, "below: {below}");
    assert!((above - f_i).abs() < 1e-8, "above: {above}");
    assert!((below - above).abs() < 1e-8);
}

#[test]
fn small_q_samples_the_upstream_population_and_large_q_the_opposite_one() {
    let (f_i, upstream, opposite) = (0.2, 0.6, 0.8);
    // q < 1/2: 2 q f_i(x_f) + (1 - 2 q) f_i(x_f - c_i); the opposite population is not used.
    let bc = InterpolatedBc::new(0.2);
    let expected = 0.4 * f_i + 0.6 * upstream;
    assert!((bc.reflected_population(f_i, upstream, opposite) - expected).abs() < 1e-15);
    assert_eq!(
        bc.reflected_population(f_i, upstream, opposite),
        bc.reflected_population(f_i, upstream, -5.0)
    );
    // q >= 1/2: f_i / (2 q) + (2 q - 1) / (2 q) f_ī(x_f); the upstream population is not used.
    let bc = InterpolatedBc::new(0.8);
    let expected = f_i / 1.6 + 0.6 / 1.6 * opposite;
    assert!((bc.reflected_population(f_i, upstream, opposite) - expected).abs() < 1e-15);
    assert_eq!(
        bc.reflected_population(f_i, upstream, opposite),
        bc.reflected_population(f_i, -5.0, opposite)
    );
    // q = 1: the wall sits on the neighbouring node, f_ī = (f_i + f_ī(x_f)) / 2.
    let bc = InterpolatedBc::new(1.0);
    assert!(
        (bc.reflected_population(f_i, upstream, opposite) - 0.5 * (f_i + opposite)).abs() < 1e-15
    );
}

#[test]
fn rest_state_is_preserved_for_every_distance_fraction() {
    // The coefficients of both branches sum to one, so uniform populations are reproduced.
    for k in 1..=20 {
        let q = f64::from(k) / 20.0;
        let bc = InterpolatedBc::new(q);
        let value = 0.0731;
        assert!(
            (bc.reflected_population(value, value, value) - value).abs() < 1e-15,
            "q = {q}"
        );
    }
}

#[test]
fn apply_overwrites_only_the_population_pointing_out_of_the_wall() {
    for (i, &i_opp) in D2Q9_OPPOSITES.iter().enumerate().skip(1) {
        let mut f = [0.0_f64; 9];
        for (k, f_k) in f.iter_mut().enumerate() {
            *f_k = 0.01 * (k as f64 + 1.0);
        }
        let before = f;
        let bc = InterpolatedBc::new(0.3);
        bc.apply(&mut f, 0.5, 0.25, 0.125, i);
        let expected = bc.reflected_population(0.5, 0.25, 0.125);
        for (k, (&after_k, &before_k)) in f.iter().zip(before.iter()).enumerate() {
            if k == i_opp {
                assert_eq!(after_k, expected, "direction {i}: slot {k}");
            } else {
                assert_eq!(after_k, before_k, "direction {i}: slot {k} must not change");
            }
        }
    }
}

#[test]
fn try_new_accepts_only_fractions_in_the_half_open_unit_interval() {
    for q in [1e-6, 0.25, 0.5, 0.999, 1.0] {
        let bc = InterpolatedBc::try_new(q);
        assert_eq!(bc, Ok(InterpolatedBc::new(q)), "q = {q}");
    }
    for q in [0.0, -0.1, 1.000_001, 2.0, f64::NAN, f64::INFINITY] {
        match InterpolatedBc::try_new(q) {
            Err(crate::Error::BoundaryConditionError { face, reason }) => {
                assert_eq!(face, "interpolated");
                assert!(reason.contains("(0, 1]"), "{reason}");
            }
            other => panic!("q = {q}: expected a boundary condition error, got {other:?}"),
        }
    }
}

/// Body-force driven channel between two flat walls that cut the links off the lattice:
/// the `n` fluid rows sit at `y = 0..n-1`, the lower wall at `y = -q_bottom` and the upper wall
/// at `y = n - 1 + q_top`. One node wide and periodic in x; BGK collision with Guo forcing and
/// [`InterpolatedBc`] on every link that crosses a wall.
struct Channel {
    n: usize,
    tau: f64,
    acceleration: f64,
    bottom: InterpolatedBc,
    top: InterpolatedBc,
    f: Vec<[f64; 9]>,
}

impl Channel {
    /// Start from the analytic parabola so that only the wall error has to relax.
    fn new(n: usize, tau: f64, acceleration: f64, q_bottom: f64, q_top: f64) -> Self {
        let nu = (tau - 0.5) / 3.0;
        let f = (0..n)
            .map(|j| {
                let y = j as f64;
                let u = acceleration / (2.0 * nu) * (y + q_bottom) * (n as f64 - 1.0 + q_top - y);
                let mut node = [0.0; 9];
                for (i, f_i) in node.iter_mut().enumerate() {
                    *f_i = feq_single(i, 1.0, u, 0.0);
                }
                node
            })
            .collect();
        Self {
            n,
            tau,
            acceleration,
            bottom: InterpolatedBc::new(q_bottom),
            top: InterpolatedBc::new(q_top),
            f,
        }
    }

    /// Macroscopic density and velocity of a node (with the Guo half-force correction).
    fn macroscopic(&self, node: &[f64; 9]) -> (f64, f64, f64) {
        let rho: f64 = node.iter().sum();
        let jx: f64 = node.iter().enumerate().map(|(i, f_i)| f_i * c_x(i)).sum();
        let jy: f64 = node.iter().enumerate().map(|(i, f_i)| f_i * c_y(i)).sum();
        (rho, (jx + 0.5 * rho * self.acceleration) / rho, jy / rho)
    }

    fn step(&mut self) {
        let omega = 1.0 / self.tau;
        let n = self.n;
        let mut post = vec![[0.0_f64; 9]; n];
        for (j, post_j) in post.iter_mut().enumerate() {
            let (rho, ux, uy) = self.macroscopic(&self.f[j]);
            for (i, post_ji) in post_j.iter_mut().enumerate() {
                let eu = c_x(i) * ux + c_y(i) * uy;
                let force = (1.0 - 0.5 * omega)
                    * D2Q9_WEIGHTS[i]
                    * (3.0 * (c_x(i) - ux) + 9.0 * eu * c_x(i))
                    * rho
                    * self.acceleration;
                *post_ji =
                    self.f[j][i] - omega * (self.f[j][i] - feq_single(i, rho, ux, uy)) + force;
            }
        }
        let mut next = vec![[0.0_f64; 9]; n];
        for (j, next_j) in next.iter_mut().enumerate() {
            for (i, next_ji) in next_j.iter_mut().enumerate() {
                let source = j as isize - D2Q9_VELOCITIES[i][1] as isize;
                if (0..n as isize).contains(&source) {
                    *next_ji = post[source as usize][i];
                }
            }
            // Populations entering from a wall: i is the unknown, `towards_wall` the link that
            // crosses the wall, `upstream` the row one node further away from the wall.
            let (wall_populations, bc): ([usize; 3], &InterpolatedBc) = if j == 0 {
                ([2, 5, 6], &self.bottom)
            } else if j == n - 1 {
                ([4, 7, 8], &self.top)
            } else {
                continue;
            };
            for i in wall_populations {
                let towards_wall = D2Q9_OPPOSITES[i];
                let upstream = (j as isize + D2Q9_VELOCITIES[i][1] as isize) as usize;
                bc.apply(
                    next_j,
                    post[j][towards_wall],
                    post[upstream][towards_wall],
                    post[j][i],
                    towards_wall,
                );
            }
        }
        self.f = next;
    }

    fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    fn velocity_profile(&self) -> Vec<f64> {
        self.f.iter().map(|node| self.macroscopic(node).1).collect()
    }
}

/// Relative L2 error of the converged channel against the Poiseuille parabola whose zeros sit at
/// the true wall positions.
fn poiseuille_error(n: usize, tau: f64, q_bottom: f64, q_top: f64) -> f64 {
    // Keep the peak velocity independent of the resolution.
    let acceleration = 1e-6 * (8.0 / n as f64).powi(2);
    let nu = (tau - 0.5) / 3.0;
    let height = n as f64 - 1.0 + q_bottom + q_top;
    let mut channel = Channel::new(n, tau, acceleration, q_bottom, q_top);
    // The diffusive relaxation time of the slowest mode is H^2 / (pi^2 nu); this is ~12 times it.
    channel.run((1.2 * height * height / nu) as usize + 100);
    let profile = channel.velocity_profile();
    let mut diff = 0.0;
    let mut norm = 0.0;
    for (j, &u) in profile.iter().enumerate() {
        let y = j as f64;
        let exact = acceleration / (2.0 * nu) * (y + q_bottom) * (n as f64 - 1.0 + q_top - y);
        diff += (u - exact).powi(2);
        norm += exact * exact;
    }
    (diff / norm).sqrt()
}

#[test]
fn poiseuille_flow_with_off_lattice_walls_converges_at_second_order() {
    let tau = 0.8;
    for (q_bottom, q_top) in [
        (0.25, 0.25),
        (0.5, 0.5),
        (0.75, 0.75),
        (1.0, 1.0),
        (0.3, 0.7),
    ] {
        let errors: Vec<f64> = [8usize, 16, 32]
            .into_iter()
            .map(|n| poiseuille_error(n, tau, q_bottom, q_top))
            .collect();
        assert!(
            errors[2] < 5e-3,
            "q = ({q_bottom}, {q_top}): error on the finest grid {}",
            errors[2]
        );
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                order > 1.68,
                "q = ({q_bottom}, {q_top}): observed order {order} ({} -> {})",
                pair[0],
                pair[1]
            );
        }
    }
}

#[test]
fn poiseuille_error_is_small_for_other_relaxation_times() {
    for tau in [1.0, 1.4] {
        let coarse = poiseuille_error(8, tau, 0.25, 0.25);
        let fine = poiseuille_error(16, tau, 0.25, 0.25);
        assert!(fine < 0.06, "tau = {tau}: {fine}");
        assert!(
            coarse / fine > 3.2,
            "tau = {tau}: {coarse} -> {fine} is not second order"
        );
    }
}
