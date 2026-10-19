// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Flow-level validation of the curved boundary schemes: plane Couette and
//! Poiseuille channels with off-lattice walls, and circular Couette flow
//! around a rotating cylinder.

use super::*;

/// Minimal D2Q9 BGK solver with Guo forcing driving a [`CurvedBoundarySet`].
struct Solver {
    lattice: CurvedLattice,
    domain: LatticeDomain,
    fluid: Vec<bool>,
    f: Vec<Vec<f64>>,
    omega: f64,
    force: [f64; 3],
}

impl Solver {
    fn new(domain: LatticeDomain, fluid: Vec<bool>, tau: f64, force: [f64; 3]) -> Self {
        let lattice = CurvedLattice::d2q9();
        let f = (0..domain.num_nodes())
            .map(|_| {
                (0..9)
                    .map(|i| lattice.equilibrium(i, 1.0, [0.0; 3]))
                    .collect()
            })
            .collect();
        Self {
            lattice,
            domain,
            fluid,
            f,
            omega: 1.0 / tau,
            force,
        }
    }

    fn velocity(&self, idx: usize) -> [f64; 3] {
        let (rho, j) = self.lattice.moments(&self.f[idx]).expect("moments");
        [
            (j[0] + 0.5 * self.force[0]) / rho,
            (j[1] + 0.5 * self.force[1]) / rho,
            0.0,
        ]
    }

    /// One collide–stream–boundary step; returns the max velocity change.
    fn step(&mut self, set: &CurvedBoundarySet) -> f64 {
        let n = self.domain.num_nodes();
        let u_old: Vec<[f64; 3]> = (0..n).map(|k| self.velocity(k)).collect();
        let mut f_post = self.f.clone();
        for k in 0..n {
            if !self.fluid[k] {
                continue;
            }
            let (rho, _) = self.lattice.moments(&self.f[k]).expect("moments");
            let u = u_old[k];
            for (i, &w) in self.lattice.weights().iter().enumerate() {
                let c = self.lattice.velocity_f64(i).expect("c");
                let cu = c[0] * u[0] + c[1] * u[1];
                let guo = w
                    * (1.0 - 0.5 * self.omega)
                    * ((3.0 * (c[0] - u[0]) + 9.0 * cu * c[0]) * self.force[0]
                        + (3.0 * (c[1] - u[1]) + 9.0 * cu * c[1]) * self.force[1]);
                let feq = self.lattice.equilibrium(i, rho, u);
                f_post[k][i] = self.f[k][i] - self.omega * (self.f[k][i] - feq) + guo;
            }
        }
        for k in 0..n {
            if !self.fluid[k] {
                continue;
            }
            let node = self.domain.coords(k);
            for (i, c) in self.lattice.velocities().iter().enumerate() {
                if let Some(src) = self.domain.neighbour(node, [-c[0], -c[1], 0]) {
                    let s = self.domain.linear_index(src[0], src[1], src[2]);
                    if self.fluid[s] {
                        self.f[k][i] = f_post[s][i];
                    }
                }
            }
        }
        let force = self.force;
        let bf = move |_: usize| force;
        set.apply(&f_post, &mut self.f, &self.lattice, self.omega, Some(&bf))
            .expect("apply curved BC");
        (0..n)
            .filter(|&k| self.fluid[k])
            .map(|k| {
                let u = self.velocity(k);
                (u[0] - u_old[k][0]).abs().max((u[1] - u_old[k][1]).abs())
            })
            .fold(0.0, f64::max)
    }

    fn run(&mut self, set: &CurvedBoundarySet, max_steps: usize, tol: f64) {
        for _ in 0..max_steps {
            if self.step(set) < tol {
                return;
            }
        }
    }
}

/// Channel with `n_fluid` fluid rows (y = 1..=n_fluid) between a bottom wall at
/// `y = 1 − q_bot` and a top wall at `y = n_fluid + q_top`; periodic in x.
struct Channel {
    solver: Solver,
    set: CurvedBoundarySet,
    y_bot: f64,
    y_top: f64,
}

fn channel(
    n_fluid: usize,
    q_bot: f64,
    q_top: f64,
    u_top: f64,
    tau: f64,
    force_x: f64,
    method: CurvedBcMethod,
) -> Channel {
    let nx = 2;
    let domain = LatticeDomain::new(nx, n_fluid + 2, 1, 1.0)
        .expect("domain")
        .with_periodic([true, false, false]);
    let y_bot = 1.0 - q_bot;
    let y_top = n_fluid as f64 + q_top;
    let sdf = move |p: [f64; 3]| (p[1] - y_bot).min(y_top - p[1]);
    let mid = 0.5 * (y_bot + y_top);
    let wall = move |p: [f64; 3]| {
        if p[1] > mid {
            [u_top, 0.0, 0.0]
        } else {
            [0.0; 3]
        }
    };
    let set =
        CurvedBoundarySet::from_geometry_with(&sdf, &wall, &domain, &CurvedLattice::d2q9(), method);
    let fluid = (0..domain.num_nodes())
        .map(|k| sdf(domain.position(domain.coords(k))) > 0.0)
        .collect();
    Channel {
        solver: Solver::new(domain, fluid, tau, [force_x, 0.0, 0.0]),
        set,
        y_bot,
        y_top,
    }
}

impl Channel {
    /// Max and RMS error of `u_x` against `exact(y)`, normalised by `scale`.
    fn errors(&self, exact: &dyn Fn(f64) -> f64, scale: f64) -> (f64, f64) {
        let mut max_err = 0.0_f64;
        let mut sum = 0.0;
        let mut count = 0usize;
        for k in 0..self.solver.domain.num_nodes() {
            if !self.solver.fluid[k] {
                continue;
            }
            let y = self.solver.domain.position(self.solver.domain.coords(k))[1];
            let e = (self.solver.velocity(k)[0] - exact(y)).abs() / scale;
            max_err = max_err.max(e);
            sum += e * e;
            count += 1;
        }
        (max_err, (sum / count as f64).sqrt())
    }
}

const ALL_METHODS: [CurvedBcMethod; 4] = [
    CurvedBcMethod::Bouzidi,
    CurvedBcMethod::QuadraticInterpolation,
    CurvedBcMethod::FilippovaHanel,
    CurvedBcMethod::MomentBased,
];

#[test]
fn test_couette_offset_walls_all_methods() {
    let u_top = 0.05;
    for &(q_bot, q_top) in &[(0.3, 0.7), (0.7, 0.3)] {
        for method in ALL_METHODS {
            let mut ch = channel(8, q_bot, q_top, u_top, 0.9, 0.0, method);
            ch.solver.run(&ch.set, 40_000, 1e-13);
            let (y_b, y_t) = (ch.y_bot, ch.y_top);
            let exact = move |y: f64| u_top * (y - y_b) / (y_t - y_b);
            let (max_err, _) = ch.errors(&exact, u_top);
            assert!(
                max_err < 1e-6,
                "{method:?} q=({q_bot},{q_top}): max rel error {max_err}"
            );
        }
    }
}

#[test]
fn test_couette_half_way_placement_is_wrong_for_offset_walls() {
    // Forcing q = 1/2 on walls that really sit at q = 0.3 / 0.7 must show the
    // wall-position error that the curved schemes remove.
    let u_top = 0.05;
    let mut ch = channel(8, 0.3, 0.7, u_top, 0.9, 0.0, CurvedBcMethod::Bouzidi);
    for node in &mut ch.set.nodes {
        node.q = 0.5;
    }
    ch.solver.run(&ch.set, 40_000, 1e-13);
    let (y_b, y_t) = (ch.y_bot, ch.y_top);
    let exact = move |y: f64| u_top * (y - y_b) / (y_t - y_b);
    let (max_err, _) = ch.errors(&exact, u_top);
    assert!(max_err > 1e-2, "half-way placement error only {max_err}");
}

#[test]
fn test_poiseuille_second_order_convergence() {
    let tau = 0.8;
    let nu = (tau - 0.5) / 3.0;
    let u_max = 0.01;
    for method in ALL_METHODS {
        let mut rms = Vec::new();
        for &n_fluid in &[8usize, 16] {
            let (q_bot, q_top) = (0.3, 0.6);
            let h = n_fluid as f64 - 1.0 + q_bot + q_top;
            let g = 8.0 * nu * u_max / (h * h);
            let mut ch = channel(n_fluid, q_bot, q_top, 0.0, tau, g, method);
            ch.solver.run(&ch.set, 200_000, 1e-11);
            let (y_b, y_t) = (ch.y_bot, ch.y_top);
            let exact = move |y: f64| g / (2.0 * nu) * (y - y_b) * (y_t - y);
            rms.push(ch.errors(&exact, u_max).1);
        }
        let order = (rms[0] / rms[1]).log2();
        assert!(rms[1] < 2e-2, "{method:?}: rms error {} too large", rms[1]);
        assert!(
            order > 1.7,
            "{method:?}: observed order {order} (errors {rms:?})"
        );
    }
}

#[test]
fn test_rotating_cylinder_circular_couette() {
    // Inner cylinder (radius r1) rotating at ω inside a stationary outer
    // cylinder (radius r2): u_θ(r) = A r + B / r.
    let (r1, r2) = (4.3, 11.2);
    let omega_rot = 0.004;
    let center = [13.1, 12.9];
    let n = 27;
    let domain = LatticeDomain::new(n, n, 1, 1.0).expect("domain");
    let lattice = CurvedLattice::d2q9();
    let inner = RotatingCylinder::new(center, r1, omega_rot);
    let mut set = inner.boundary_set_with(&domain, &lattice, CurvedBcMethod::Bouzidi);
    let outer_sdf = |p: [f64; 3]| r2 - (p[0] - center[0]).hypot(p[1] - center[1]);
    let outer = CurvedBoundarySet::from_geometry(&outer_sdf, &domain, &lattice);
    set.nodes.extend(outer.nodes);
    let fluid: Vec<bool> = (0..domain.num_nodes())
        .map(|k| {
            let p = domain.position(domain.coords(k));
            inner.signed_distance([p[0], p[1]]) > 0.0 && outer_sdf(p) > 0.0
        })
        .collect();
    let mut solver = Solver::new(domain, fluid, 0.8, [0.0; 3]);
    solver.run(&set, 30_000, 1e-12);
    let a = -omega_rot * r1 * r1 / (r2 * r2 - r1 * r1);
    let b = omega_rot * r1 * r1 * r2 * r2 / (r2 * r2 - r1 * r1);
    let u_wall = omega_rot * r1;
    let mut max_err = 0.0_f64;
    for k in 0..domain.num_nodes() {
        if !solver.fluid[k] {
            continue;
        }
        let p = domain.position(domain.coords(k));
        let (dx, dy) = (p[0] - center[0], p[1] - center[1]);
        let r = dx.hypot(dy);
        let u = solver.velocity(k);
        let u_theta = (-dy * u[0] + dx * u[1]) / r;
        let exact = a * r + b / r;
        max_err = max_err.max((u_theta - exact).abs() / u_wall);
    }
    assert!(max_err < 0.05, "circular Couette max rel error {max_err}");
}
