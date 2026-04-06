//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// A discrete Lagrangian L_d(q_k, q_{k+1}, h) that approximates the action
/// integral over one time step.
///
/// In the simplest case (midpoint rule):
///   L_d = h * L( (q_k + q_{k+1})/2, (q_{k+1} - q_k)/h )
///
/// where L(q, v) = T(v) - V(q) is the continuous Lagrangian.
#[derive(Debug, Clone)]
pub struct DiscreteLagrangian {
    /// Mass of the particle.
    pub mass: f64,
    /// Gravitational acceleration (applied in -z direction).
    pub gravity: f64,
    /// Spring stiffness for a harmonic potential V = 0.5*k*|q|^2.
    pub spring_k: f64,
}
impl DiscreteLagrangian {
    /// Create a new discrete Lagrangian.
    pub fn new(mass: f64, gravity: f64, spring_k: f64) -> Self {
        Self {
            mass,
            gravity,
            spring_k,
        }
    }
    /// Evaluate the discrete Lagrangian using the midpoint rule.
    ///
    /// L_d(q_k, q_{k+1}, h) = h * \[T(v_mid) - V(q_mid)\]
    /// where v_mid = (q_{k+1} - q_k)/h, q_mid = (q_k + q_{k+1})/2.
    #[allow(dead_code)]
    pub fn evaluate(&self, q_k: [f64; 3], q_kp1: [f64; 3], h: f64) -> f64 {
        let v_mid = scale3(sub3(q_kp1, q_k), 1.0 / h);
        let q_mid = scale3(add3(q_k, q_kp1), 0.5);
        let kinetic = 0.5 * self.mass * dot3(v_mid, v_mid);
        let potential_grav = self.mass * self.gravity * q_mid[2];
        let potential_spring = 0.5 * self.spring_k * dot3(q_mid, q_mid);
        h * (kinetic - potential_grav - potential_spring)
    }
    /// Compute D1 L_d: partial derivative with respect to q_k.
    ///
    /// Uses finite differences for generality.
    #[allow(dead_code)]
    pub fn d1(&self, q_k: [f64; 3], q_kp1: [f64; 3], h: f64) -> [f64; 3] {
        let eps = 1e-7;
        let mut result = [0.0; 3];
        for i in 0..3 {
            let mut q_plus = q_k;
            let mut q_minus = q_k;
            q_plus[i] += eps;
            q_minus[i] -= eps;
            result[i] =
                (self.evaluate(q_plus, q_kp1, h) - self.evaluate(q_minus, q_kp1, h)) / (2.0 * eps);
        }
        result
    }
    /// Compute D2 L_d: partial derivative with respect to q_{k+1}.
    #[allow(dead_code)]
    pub fn d2(&self, q_k: [f64; 3], q_kp1: [f64; 3], h: f64) -> [f64; 3] {
        let eps = 1e-7;
        let mut result = [0.0; 3];
        for i in 0..3 {
            let mut q_plus = q_kp1;
            let mut q_minus = q_kp1;
            q_plus[i] += eps;
            q_minus[i] -= eps;
            result[i] =
                (self.evaluate(q_k, q_plus, h) - self.evaluate(q_k, q_minus, h)) / (2.0 * eps);
        }
        result
    }
}
/// A holonomic constraint g(q) = 0.
///
/// Represents a scalar constraint on the configuration space. The constraint
/// manifold is defined by {q : g(q) = 0}. The null space of Dg(q) gives
/// the allowed directions of motion.
#[derive(Debug, Clone)]
pub struct HolonomicConstraint {
    /// Constraint type.
    pub kind: HolonomicKind,
    /// Parameters for the constraint.
    pub params: HolonomicParams,
}
impl HolonomicConstraint {
    /// Create a distance constraint.
    pub fn distance(body_a: usize, body_b: usize, dist: f64) -> Self {
        Self {
            kind: HolonomicKind::Distance,
            params: HolonomicParams {
                distance: dist,
                normal: [0.0; 3],
                center: [0.0; 3],
                body_a,
                body_b,
            },
        }
    }
    /// Create a plane constraint.
    pub fn plane(body: usize, normal: [f64; 3], offset: f64) -> Self {
        Self {
            kind: HolonomicKind::Plane,
            params: HolonomicParams {
                distance: offset,
                normal,
                center: [0.0; 3],
                body_a: body,
                body_b: 0,
            },
        }
    }
    /// Create a sphere constraint.
    pub fn sphere(body: usize, center: [f64; 3], radius: f64) -> Self {
        Self {
            kind: HolonomicKind::Sphere,
            params: HolonomicParams {
                distance: radius,
                normal: [0.0; 3],
                center,
                body_a: body,
                body_b: 0,
            },
        }
    }
    /// Evaluate the constraint value g(q) given positions.
    #[allow(dead_code)]
    pub fn evaluate(&self, positions: &[[f64; 3]]) -> f64 {
        match self.kind {
            HolonomicKind::Distance | HolonomicKind::RigidRod => {
                let qa = positions[self.params.body_a];
                let qb = positions[self.params.body_b];
                let diff = sub3(qa, qb);
                let dist = len3(diff);
                dist - self.params.distance
            }
            HolonomicKind::Plane => {
                let q = positions[self.params.body_a];
                dot3(self.params.normal, q) - self.params.distance
            }
            HolonomicKind::Sphere => {
                let q = positions[self.params.body_a];
                let diff = sub3(q, self.params.center);
                len3(diff) - self.params.distance
            }
        }
    }
    /// Compute the constraint Jacobian dg/dq for the involved bodies.
    ///
    /// Returns a vector of (body_index, gradient) pairs.
    #[allow(dead_code)]
    pub fn jacobian(&self, positions: &[[f64; 3]]) -> Vec<(usize, [f64; 3])> {
        match self.kind {
            HolonomicKind::Distance | HolonomicKind::RigidRod => {
                let qa = positions[self.params.body_a];
                let qb = positions[self.params.body_b];
                let diff = sub3(qa, qb);
                let dist = len3(diff);
                if dist < 1e-15 {
                    return vec![
                        (self.params.body_a, [1.0, 0.0, 0.0]),
                        (self.params.body_b, [-1.0, 0.0, 0.0]),
                    ];
                }
                let n = scale3(diff, 1.0 / dist);
                vec![
                    (self.params.body_a, n),
                    (self.params.body_b, scale3(n, -1.0)),
                ]
            }
            HolonomicKind::Plane => vec![(self.params.body_a, self.params.normal)],
            HolonomicKind::Sphere => {
                let q = positions[self.params.body_a];
                let diff = sub3(q, self.params.center);
                let dist = len3(diff);
                if dist < 1e-15 {
                    return vec![(self.params.body_a, [1.0, 0.0, 0.0])];
                }
                let n = scale3(diff, 1.0 / dist);
                vec![(self.params.body_a, n)]
            }
        }
    }
}
/// Discrete Euler-Lagrange (DEL) equation solver.
///
/// The DEL equation is:
///   D2 L_d(q_{k-1}, q_k, h) + D1 L_d(q_k, q_{k+1}, h) = 0
///
/// Given q_{k-1} and q_k, we solve for q_{k+1}.
#[derive(Debug, Clone)]
pub struct DiscreteEulerLagrange {
    /// The discrete Lagrangian.
    pub lagrangian: DiscreteLagrangian,
    /// Time step.
    pub dt: f64,
    /// Newton solver tolerance.
    pub tol: f64,
    /// Maximum Newton iterations.
    pub max_iter: usize,
}
impl DiscreteEulerLagrange {
    /// Create a new DEL solver.
    pub fn new(lagrangian: DiscreteLagrangian, dt: f64) -> Self {
        Self {
            lagrangian,
            dt,
            tol: 1e-10,
            max_iter: 50,
        }
    }
    /// Compute the DEL residual: D2 L_d(q_{k-1}, q_k, h) + D1 L_d(q_k, q_{k+1}, h).
    #[allow(dead_code)]
    pub fn residual(&self, q_km1: [f64; 3], q_k: [f64; 3], q_kp1: [f64; 3]) -> [f64; 3] {
        let d2_prev = self.lagrangian.d2(q_km1, q_k, self.dt);
        let d1_next = self.lagrangian.d1(q_k, q_kp1, self.dt);
        add3(d2_prev, d1_next)
    }
    /// Solve for q_{k+1} given q_{k-1} and q_k using Newton's method.
    ///
    /// Returns the converged q_{k+1} and the number of iterations used.
    #[allow(dead_code)]
    pub fn solve(&self, q_km1: [f64; 3], q_k: [f64; 3]) -> ([f64; 3], usize) {
        let mut q_kp1 = add3(q_k, sub3(q_k, q_km1));
        let eps = 1e-7;
        for iter in 0..self.max_iter {
            let res = self.residual(q_km1, q_k, q_kp1);
            let res_norm = len3(res);
            if res_norm < self.tol {
                return (q_kp1, iter);
            }
            let mut jac = [[0.0; 3]; 3];
            for j in 0..3 {
                let mut q_plus = q_kp1;
                let mut q_minus = q_kp1;
                q_plus[j] += eps;
                q_minus[j] -= eps;
                let r_plus = self.residual(q_km1, q_k, q_plus);
                let r_minus = self.residual(q_km1, q_k, q_minus);
                for i in 0..3 {
                    jac[i][j] = (r_plus[i] - r_minus[i]) / (2.0 * eps);
                }
            }
            let inv_jac = inv3(jac);
            let dq = mat_vec3(inv_jac, res);
            q_kp1 = sub3(q_kp1, dq);
        }
        (q_kp1, self.max_iter)
    }
    /// Step forward from two configurations: q_{k-1}, q_k -> q_{k+1}.
    #[allow(dead_code)]
    pub fn step(&self, q_km1: [f64; 3], q_k: [f64; 3]) -> [f64; 3] {
        self.solve(q_km1, q_k).0
    }
}
/// Symplectic partitioned Runge-Kutta (SPRK) integrator.
///
/// For separable Hamiltonians H(q,p) = T(p) + V(q), SPRK methods
/// update q and p in interleaved stages.
#[derive(Debug, Clone)]
pub struct SymplecticPRK {
    /// Mass.
    pub mass: f64,
    /// Time step.
    pub dt: f64,
    /// Coefficients b_i for momentum updates.
    pub b_coeffs: Vec<f64>,
    /// Coefficients B_i for position updates.
    pub b_hat_coeffs: Vec<f64>,
}
impl SymplecticPRK {
    /// Create a first-order symplectic Euler (1-stage SPRK).
    pub fn symplectic_euler(mass: f64, dt: f64) -> Self {
        Self {
            mass,
            dt,
            b_coeffs: vec![1.0],
            b_hat_coeffs: vec![1.0],
        }
    }
    /// Create a second-order Stormer-Verlet (2-stage SPRK).
    pub fn verlet(mass: f64, dt: f64) -> Self {
        Self {
            mass,
            dt,
            b_coeffs: vec![0.5, 0.5],
            b_hat_coeffs: vec![1.0, 0.0],
        }
    }
    /// Perform one step of the SPRK method.
    #[allow(dead_code)]
    pub fn step(
        &self,
        q: [f64; 3],
        p: [f64; 3],
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let stages = self.b_coeffs.len();
        let mut q_cur = q;
        let mut p_cur = p;
        for s in 0..stages {
            let f = force_fn(q_cur);
            p_cur = add3(p_cur, scale3(f, self.b_coeffs[s] * self.dt));
            q_cur = add3(
                q_cur,
                scale3(p_cur, self.b_hat_coeffs[s] * self.dt / self.mass),
            );
        }
        (q_cur, p_cur)
    }
}
/// Variational integrator for N-body gravitational systems.
///
/// Uses the discrete Euler-Lagrange equations for the N-body Lagrangian
/// with pairwise gravitational interactions.
#[derive(Debug, Clone)]
pub struct NBodyVariational {
    /// Masses of the bodies.
    pub masses: Vec<f64>,
    /// Gravitational constant.
    pub g_const: f64,
    /// Time step.
    pub dt: f64,
    /// Softening parameter to avoid singularities.
    pub softening: f64,
}
impl NBodyVariational {
    /// Create a new N-body variational integrator.
    pub fn new(masses: Vec<f64>, g_const: f64, dt: f64) -> Self {
        Self {
            masses,
            g_const,
            dt,
            softening: 1e-4,
        }
    }
    /// Compute gravitational force on particle i from all others.
    #[allow(dead_code)]
    pub fn gravitational_force(&self, i: usize, positions: &[[f64; 3]]) -> [f64; 3] {
        let mut force = [0.0; 3];
        let n = self.masses.len();
        for j in 0..n {
            if i == j {
                continue;
            }
            let diff = sub3(positions[j], positions[i]);
            let dist2 = dot3(diff, diff) + self.softening * self.softening;
            let dist = dist2.sqrt();
            let f_mag = self.g_const * self.masses[i] * self.masses[j] / dist2;
            let f_dir = scale3(diff, f_mag / dist);
            force = add3(force, f_dir);
        }
        force
    }
    /// Perform one Stormer-Verlet step for the N-body system.
    #[allow(dead_code)]
    pub fn step(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = self.masses.len();
        let mut v_half = Vec::with_capacity(n);
        for i in 0..n {
            let f = self.gravitational_force(i, positions);
            let a = scale3(f, 1.0 / self.masses[i]);
            v_half.push(add3(velocities[i], scale3(a, 0.5 * self.dt)));
        }
        let mut q_new = Vec::with_capacity(n);
        for i in 0..n {
            q_new.push(add3(positions[i], scale3(v_half[i], self.dt)));
        }
        let mut v_new = Vec::with_capacity(n);
        for i in 0..n {
            let f = self.gravitational_force(i, &q_new);
            let a = scale3(f, 1.0 / self.masses[i]);
            v_new.push(add3(v_half[i], scale3(a, 0.5 * self.dt)));
        }
        (q_new, v_new)
    }
    /// Compute total energy (kinetic + gravitational potential).
    #[allow(dead_code)]
    pub fn total_energy(&self, positions: &[[f64; 3]], velocities: &[[f64; 3]]) -> f64 {
        let n = self.masses.len();
        let mut ke = 0.0;
        for i in 0..n {
            ke += 0.5 * self.masses[i] * dot3(velocities[i], velocities[i]);
        }
        let mut pe = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let diff = sub3(positions[i], positions[j]);
                let dist = (dot3(diff, diff) + self.softening * self.softening).sqrt();
                pe -= self.g_const * self.masses[i] * self.masses[j] / dist;
            }
        }
        ke + pe
    }
    /// Compute total linear momentum.
    #[allow(dead_code)]
    pub fn total_momentum(&self, velocities: &[[f64; 3]]) -> [f64; 3] {
        let mut p = [0.0; 3];
        for i in 0..self.masses.len() {
            p = add3(p, scale3(velocities[i], self.masses[i]));
        }
        p
    }
    /// Compute total angular momentum.
    #[allow(dead_code)]
    pub fn total_angular_momentum(
        &self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> [f64; 3] {
        let mut l = [0.0; 3];
        for i in 0..self.masses.len() {
            let p = scale3(velocities[i], self.masses[i]);
            l = add3(l, cross3(positions[i], p));
        }
        l
    }
}
/// Stormer-Verlet (leapfrog) symplectic integrator for Hamiltonian systems.
///
/// State: position q, momentum p.
/// Equations: dq/dt = p/m, dp/dt = -dV/dq.
#[derive(Debug, Clone)]
pub struct SymplecticVerlet {
    /// Mass.
    pub mass: f64,
    /// Time step.
    pub dt: f64,
}
impl SymplecticVerlet {
    /// Create a new symplectic Verlet integrator.
    pub fn new(mass: f64, dt: f64) -> Self {
        Self { mass, dt }
    }
    /// Perform one step: (q, p) -> (q', p') given a force function.
    ///
    /// The force function takes a position and returns the force.
    #[allow(dead_code)]
    pub fn step(
        &self,
        q: [f64; 3],
        p: [f64; 3],
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let f0 = force_fn(q);
        let p_half = add3(p, scale3(f0, 0.5 * self.dt));
        let q_new = add3(q, scale3(p_half, self.dt / self.mass));
        let f1 = force_fn(q_new);
        let p_new = add3(p_half, scale3(f1, 0.5 * self.dt));
        (q_new, p_new)
    }
    /// Integrate over multiple steps and return the trajectory.
    #[allow(dead_code)]
    pub fn integrate(
        &self,
        q0: [f64; 3],
        p0: [f64; 3],
        steps: usize,
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> Vec<([f64; 3], [f64; 3])> {
        let mut trajectory = Vec::with_capacity(steps + 1);
        trajectory.push((q0, p0));
        let mut q = q0;
        let mut p = p0;
        for _ in 0..steps {
            let (qn, pn) = self.step(q, p, force_fn);
            q = qn;
            p = pn;
            trajectory.push((q, p));
        }
        trajectory
    }
    /// Compute the Hamiltonian H = |p|^2/(2m) + V(q) for energy monitoring.
    #[allow(dead_code)]
    pub fn hamiltonian(&self, p: [f64; 3], potential: f64) -> f64 {
        dot3(p, p) / (2.0 * self.mass) + potential
    }
}
/// Lie group integrator for rigid body rotation (RKMK-style).
///
/// Updates rotation matrix R and body angular velocity omega using
/// the Lie group structure of SO(3).
#[derive(Debug, Clone)]
pub struct LieGroupIntegrator {
    /// Time step.
    pub dt: f64,
}
impl LieGroupIntegrator {
    /// Create a new Lie group integrator.
    pub fn new(dt: f64) -> Self {
        Self { dt }
    }
    /// Perform one integration step for a free rigid body (no external torques).
    ///
    /// Uses the implicit midpoint rule on the Lie group:
    ///   R_{k+1} = R_k * exp(h * omega_mid)
    ///   I * omega_{k+1} = I * omega_k + h * (I*omega_mid) x omega_mid
    ///
    /// For torque-free motion, angular momentum is conserved.
    #[allow(dead_code)]
    pub fn step_free(&self, state: &LieGroupState) -> LieGroupState {
        let omega = state.omega_body;
        let inertia = state.inertia;
        let l = [
            inertia[0] * omega[0],
            inertia[1] * omega[1],
            inertia[2] * omega[2],
        ];
        let torque = cross3(l, omega);
        let omega_new = [
            omega[0] + self.dt * torque[0] / inertia[0],
            omega[1] + self.dt * torque[1] / inertia[1],
            omega[2] + self.dt * torque[2] / inertia[2],
        ];
        let omega_mid = scale3(add3(omega, omega_new), 0.5);
        let dtheta = scale3(omega_mid, self.dt);
        let delta_r = exp_so3(dtheta);
        let rotation_new = mat_mul3(state.rotation, delta_r);
        LieGroupState {
            rotation: rotation_new,
            omega_body: omega_new,
            inertia,
        }
    }
    /// Step with external torque (body frame).
    #[allow(dead_code)]
    pub fn step_with_torque(&self, state: &LieGroupState, torque_body: [f64; 3]) -> LieGroupState {
        let omega = state.omega_body;
        let inertia = state.inertia;
        let l = [
            inertia[0] * omega[0],
            inertia[1] * omega[1],
            inertia[2] * omega[2],
        ];
        let gyro_torque = cross3(l, omega);
        let total_torque = add3(gyro_torque, torque_body);
        let omega_new = [
            omega[0] + self.dt * total_torque[0] / inertia[0],
            omega[1] + self.dt * total_torque[1] / inertia[1],
            omega[2] + self.dt * total_torque[2] / inertia[2],
        ];
        let omega_mid = scale3(add3(omega, omega_new), 0.5);
        let dtheta = scale3(omega_mid, self.dt);
        let delta_r = exp_so3(dtheta);
        let rotation_new = mat_mul3(state.rotation, delta_r);
        LieGroupState {
            rotation: rotation_new,
            omega_body: omega_new,
            inertia,
        }
    }
    /// Integrate over multiple steps for a free rigid body.
    #[allow(dead_code)]
    pub fn integrate_free(&self, initial: &LieGroupState, steps: usize) -> Vec<LieGroupState> {
        let mut trajectory = Vec::with_capacity(steps + 1);
        trajectory.push(initial.clone());
        let mut state = initial.clone();
        for _ in 0..steps {
            state = self.step_free(&state);
            trajectory.push(state.clone());
        }
        trajectory
    }
}
/// Kinds of holonomic constraints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HolonomicKind {
    /// Distance constraint: |q_a - q_b| = d.
    Distance,
    /// Plane constraint: n . q = d.
    Plane,
    /// Sphere constraint: |q - c| = r.
    Sphere,
    /// Rigid rod constraint between two points.
    RigidRod,
}
/// Discrete null space method for projecting updates onto the constraint manifold.
///
/// Given a constraint g(q) = 0 with Jacobian G = dg/dq, the null space
/// projector is P = I - G^T (G G^T)^{-1} G. Updates are projected as
/// dq_projected = P * dq.
#[derive(Debug, Clone)]
pub struct DiscreteNullSpace {
    /// Constraint stabilization parameter (Baumgarte).
    pub beta: f64,
    /// Maximum projection iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}
impl DiscreteNullSpace {
    /// Create a new discrete null space projector.
    pub fn new(beta: f64) -> Self {
        Self {
            beta,
            max_iter: 20,
            tol: 1e-10,
        }
    }
    /// Project a position update onto the constraint manifold for a single
    /// distance constraint between two particles.
    ///
    /// Returns corrected positions for body_a and body_b.
    #[allow(dead_code)]
    pub fn project_distance(
        &self,
        q_a: [f64; 3],
        q_b: [f64; 3],
        target_dist: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
    ) -> ([f64; 3], [f64; 3]) {
        let mut qa = q_a;
        let mut qb = q_b;
        let w_total = inv_mass_a + inv_mass_b;
        if w_total < 1e-15 {
            return (qa, qb);
        }
        for _ in 0..self.max_iter {
            let diff = sub3(qa, qb);
            let dist = len3(diff);
            let error = dist - target_dist;
            if error.abs() < self.tol {
                break;
            }
            if dist < 1e-15 {
                break;
            }
            let n = scale3(diff, 1.0 / dist);
            let correction = error / w_total;
            qa = sub3(qa, scale3(n, inv_mass_a * correction));
            qb = add3(qb, scale3(n, inv_mass_b * correction));
        }
        (qa, qb)
    }
    /// Project positions onto a plane constraint n.q = d.
    #[allow(dead_code)]
    pub fn project_plane(&self, q: [f64; 3], normal: [f64; 3], offset: f64) -> [f64; 3] {
        let error = dot3(normal, q) - offset;
        sub3(q, scale3(normal, error))
    }
    /// Compute the null space projection matrix P = I - G^T (G G^T)^{-1} G
    /// for a single constraint gradient.
    #[allow(dead_code)]
    pub fn null_space_projector(&self, gradient: [f64; 3]) -> [[f64; 3]; 3] {
        let g_norm2 = dot3(gradient, gradient);
        if g_norm2 < 1e-30 {
            return identity3();
        }
        let outer = outer3(gradient, gradient);
        let inv_gnorm2 = 1.0 / g_norm2;
        let mut proj = identity3();
        for i in 0..3 {
            for j in 0..3 {
                proj[i][j] -= outer[i][j] * inv_gnorm2;
            }
        }
        proj
    }
}
/// Rotation represented as a 3x3 matrix for Lie group integration.
///
/// Uses the exponential map from so(3) to SO(3) for updating rotations.
#[derive(Debug, Clone)]
pub struct LieGroupState {
    /// Rotation matrix (row-major).
    pub rotation: [[f64; 3]; 3],
    /// Angular velocity in body frame.
    pub omega_body: [f64; 3],
    /// Inertia tensor in body frame (diagonal).
    pub inertia: [f64; 3],
}
impl LieGroupState {
    /// Create a new Lie group state with identity rotation.
    pub fn new(inertia: [f64; 3]) -> Self {
        Self {
            rotation: identity3(),
            omega_body: [0.0; 3],
            inertia,
        }
    }
    /// Create a state with given rotation and angular velocity.
    pub fn with_state(rotation: [[f64; 3]; 3], omega_body: [f64; 3], inertia: [f64; 3]) -> Self {
        Self {
            rotation,
            omega_body,
            inertia,
        }
    }
    /// Compute the body-frame angular momentum L = I * omega.
    #[allow(dead_code)]
    pub fn angular_momentum(&self) -> [f64; 3] {
        [
            self.inertia[0] * self.omega_body[0],
            self.inertia[1] * self.omega_body[1],
            self.inertia[2] * self.omega_body[2],
        ]
    }
    /// Compute the rotational kinetic energy T = 0.5 * omega^T I omega.
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * (self.inertia[0] * self.omega_body[0] * self.omega_body[0]
            + self.inertia[1] * self.omega_body[1] * self.omega_body[1]
            + self.inertia[2] * self.omega_body[2] * self.omega_body[2])
    }
}
/// Parameters for holonomic constraints.
#[derive(Debug, Clone)]
pub struct HolonomicParams {
    /// Target distance / offset.
    pub distance: f64,
    /// Normal vector (for plane constraint).
    pub normal: [f64; 3],
    /// Center point (for sphere constraint).
    pub center: [f64; 3],
    /// Index of first body.
    pub body_a: usize,
    /// Index of second body (if applicable).
    pub body_b: usize,
}
impl HolonomicParams {
    /// Create default parameters.
    pub fn default_params() -> Self {
        Self {
            distance: 1.0,
            normal: [0.0, 0.0, 1.0],
            center: [0.0; 3],
            body_a: 0,
            body_b: 1,
        }
    }
}
/// Momentum map computation for symmetry groups.
///
/// For a system with N particles, computes conserved quantities associated
/// with spatial symmetries (translation -> linear momentum, rotation -> angular momentum).
#[derive(Debug, Clone)]
pub struct MomentumMap {
    /// Number of particles.
    pub n_particles: usize,
}
impl MomentumMap {
    /// Create a new momentum map for n particles.
    pub fn new(n_particles: usize) -> Self {
        Self { n_particles }
    }
    /// Compute the total linear momentum.
    #[allow(dead_code)]
    pub fn linear_momentum(&self, momenta: &[[f64; 3]]) -> [f64; 3] {
        let mut total = [0.0; 3];
        for p in momenta.iter().take(self.n_particles) {
            total = add3(total, *p);
        }
        total
    }
    /// Compute the total angular momentum about the origin.
    ///
    /// L = sum_i (q_i x p_i)
    #[allow(dead_code)]
    pub fn angular_momentum(&self, positions: &[[f64; 3]], momenta: &[[f64; 3]]) -> [f64; 3] {
        let mut total = [0.0; 3];
        for i in 0..self.n_particles {
            total = add3(total, cross3(positions[i], momenta[i]));
        }
        total
    }
    /// Compute angular momentum about a given center.
    #[allow(dead_code)]
    pub fn angular_momentum_about(
        &self,
        positions: &[[f64; 3]],
        momenta: &[[f64; 3]],
        center: [f64; 3],
    ) -> [f64; 3] {
        let mut total = [0.0; 3];
        for i in 0..self.n_particles {
            let r = sub3(positions[i], center);
            total = add3(total, cross3(r, momenta[i]));
        }
        total
    }
    /// Compute the center of mass.
    #[allow(dead_code)]
    pub fn center_of_mass(&self, positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
        let mut com = [0.0; 3];
        let mut total_mass = 0.0;
        for i in 0..self.n_particles {
            com = add3(com, scale3(positions[i], masses[i]));
            total_mass += masses[i];
        }
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        scale3(com, 1.0 / total_mass)
    }
    /// Compute total kinetic energy.
    #[allow(dead_code)]
    pub fn kinetic_energy(&self, momenta: &[[f64; 3]], masses: &[f64]) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.n_particles {
            ke += dot3(momenta[i], momenta[i]) / (2.0 * masses[i]);
        }
        ke
    }
    /// Check if linear momentum is conserved between two snapshots.
    #[allow(dead_code)]
    pub fn is_linear_momentum_conserved(
        &self,
        momenta_before: &[[f64; 3]],
        momenta_after: &[[f64; 3]],
        tol: f64,
    ) -> bool {
        let p_before = self.linear_momentum(momenta_before);
        let p_after = self.linear_momentum(momenta_after);
        len3(sub3(p_before, p_after)) < tol
    }
    /// Check if angular momentum is conserved between two snapshots.
    #[allow(dead_code)]
    pub fn is_angular_momentum_conserved(
        &self,
        positions_before: &[[f64; 3]],
        momenta_before: &[[f64; 3]],
        positions_after: &[[f64; 3]],
        momenta_after: &[[f64; 3]],
        tol: f64,
    ) -> bool {
        let l_before = self.angular_momentum(positions_before, momenta_before);
        let l_after = self.angular_momentum(positions_after, momenta_after);
        len3(sub3(l_before, l_after)) < tol
    }
}
/// Fourth-order variational integrator using composition method.
///
/// Composes second-order steps with specific time-step fractions
/// (Yoshida's coefficients) to achieve fourth-order accuracy.
#[derive(Debug, Clone)]
pub struct FourthOrderVariational {
    /// Mass.
    pub mass: f64,
    /// Time step.
    pub dt: f64,
}
impl FourthOrderVariational {
    /// Create a new fourth-order variational integrator.
    pub fn new(mass: f64, dt: f64) -> Self {
        Self { mass, dt }
    }
    /// Yoshida coefficients for fourth-order composition.
    #[allow(dead_code)]
    fn yoshida_coefficients() -> [f64; 3] {
        let cbrt2 = 2.0_f64.cbrt();
        let w1 = 1.0 / (2.0 - cbrt2);
        let w0 = -cbrt2 / (2.0 - cbrt2);
        [w1, w0, w1]
    }
    /// Perform a leapfrog sub-step with given time step fraction.
    #[allow(dead_code)]
    fn leapfrog_substep(
        &self,
        q: [f64; 3],
        p: [f64; 3],
        sub_dt: f64,
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let f0 = force_fn(q);
        let p_half = add3(p, scale3(f0, 0.5 * sub_dt));
        let q_new = add3(q, scale3(p_half, sub_dt / self.mass));
        let f1 = force_fn(q_new);
        let p_new = add3(p_half, scale3(f1, 0.5 * sub_dt));
        (q_new, p_new)
    }
    /// Perform one fourth-order step.
    #[allow(dead_code)]
    pub fn step(
        &self,
        q: [f64; 3],
        p: [f64; 3],
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let coeffs = Self::yoshida_coefficients();
        let mut q_cur = q;
        let mut p_cur = p;
        for &w in &coeffs {
            let sub_dt = w * self.dt;
            let (qn, pn) = self.leapfrog_substep(q_cur, p_cur, sub_dt, force_fn);
            q_cur = qn;
            p_cur = pn;
        }
        (q_cur, p_cur)
    }
    /// Integrate over multiple steps.
    #[allow(dead_code)]
    pub fn integrate(
        &self,
        q0: [f64; 3],
        p0: [f64; 3],
        steps: usize,
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> Vec<([f64; 3], [f64; 3])> {
        let mut trajectory = Vec::with_capacity(steps + 1);
        trajectory.push((q0, p0));
        let mut q = q0;
        let mut p = p0;
        for _ in 0..steps {
            let (qn, pn) = self.step(q, p, force_fn);
            q = qn;
            p = pn;
            trajectory.push((q, p));
        }
        trajectory
    }
}
/// Variational collision response that preserves the symplectic structure.
///
/// Uses a discrete collision map that is derived from a variational principle
/// rather than ad-hoc velocity reflection.
#[derive(Debug, Clone)]
pub struct VariationalCollision {
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// Penetration tolerance.
    pub penetration_tol: f64,
}
impl VariationalCollision {
    /// Create a new variational collision handler.
    pub fn new(restitution: f64, friction: f64) -> Self {
        Self {
            restitution,
            friction,
            penetration_tol: 1e-6,
        }
    }
    /// Apply variational collision response for a particle hitting a plane.
    ///
    /// The plane is defined by normal `n` and offset `d` (n.q >= d).
    /// Returns (corrected_position, corrected_momentum).
    #[allow(dead_code)]
    pub fn collide_plane(
        &self,
        q: [f64; 3],
        p: [f64; 3],
        mass: f64,
        normal: [f64; 3],
        offset: f64,
    ) -> ([f64; 3], [f64; 3]) {
        let penetration = offset - dot3(normal, q);
        if penetration <= self.penetration_tol {
            return (q, p);
        }
        let q_corrected = add3(q, scale3(normal, penetration));
        let v = scale3(p, 1.0 / mass);
        let vn = dot3(v, normal);
        if vn >= 0.0 {
            return (q_corrected, p);
        }
        let impulse_n = -(1.0 + self.restitution) * vn;
        let vt = sub3(v, scale3(normal, vn));
        let vt_mag = len3(vt);
        let friction_impulse = if vt_mag > 1e-15 {
            let max_friction = self.friction * impulse_n.abs();
            let needed = vt_mag;
            let friction_mag = needed.min(max_friction);
            let tangent = normalize3(vt);
            scale3(tangent, -friction_mag)
        } else {
            [0.0; 3]
        };
        let total_impulse = add3(scale3(normal, impulse_n), friction_impulse);
        let p_corrected = add3(p, scale3(total_impulse, mass));
        (q_corrected, p_corrected)
    }
    /// Apply variational collision between two particles.
    ///
    /// Returns corrected (q_a, p_a, q_b, p_b).
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn collide_particles(
        &self,
        q_a: [f64; 3],
        p_a: [f64; 3],
        mass_a: f64,
        q_b: [f64; 3],
        p_b: [f64; 3],
        mass_b: f64,
        min_dist: f64,
    ) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
        let diff = sub3(q_a, q_b);
        let dist = len3(diff);
        if dist >= min_dist || dist < 1e-15 {
            return (q_a, p_a, q_b, p_b);
        }
        let normal = scale3(diff, 1.0 / dist);
        let penetration = min_dist - dist;
        let total_mass = mass_a + mass_b;
        let wa = mass_b / total_mass;
        let wb = mass_a / total_mass;
        let qa_new = add3(q_a, scale3(normal, wa * penetration));
        let qb_new = sub3(q_b, scale3(normal, wb * penetration));
        let va = scale3(p_a, 1.0 / mass_a);
        let vb = scale3(p_b, 1.0 / mass_b);
        let v_rel = dot3(sub3(va, vb), normal);
        if v_rel >= 0.0 {
            return (qa_new, p_a, qb_new, p_b);
        }
        let reduced_mass = (mass_a * mass_b) / total_mass;
        let impulse = (1.0 + self.restitution) * reduced_mass * v_rel;
        let pa_new = sub3(p_a, scale3(normal, impulse));
        let pb_new = add3(p_b, scale3(normal, impulse));
        (qa_new, pa_new, qb_new, pb_new)
    }
    /// Compute the energy change due to collision (should be <= 0 for valid restitution).
    #[allow(dead_code)]
    pub fn energy_change(&self, p_before: [f64; 3], p_after: [f64; 3], mass: f64) -> f64 {
        let ke_before = dot3(p_before, p_before) / (2.0 * mass);
        let ke_after = dot3(p_after, p_after) / (2.0 * mass);
        ke_after - ke_before
    }
}
/// Discrete Legendre transforms connecting the Lagrangian and Hamiltonian
/// pictures in discrete mechanics.
#[derive(Debug, Clone)]
pub struct DiscreteLegendreTransform {
    /// The underlying discrete Lagrangian.
    pub lagrangian: DiscreteLagrangian,
    /// Time step.
    pub dt: f64,
}
impl DiscreteLegendreTransform {
    /// Create a new discrete Legendre transform.
    pub fn new(lagrangian: DiscreteLagrangian, dt: f64) -> Self {
        Self { lagrangian, dt }
    }
    /// Left discrete Legendre transform: p_k^- = -D1 L_d(q_k, q_{k+1}).
    #[allow(dead_code)]
    pub fn left_transform(&self, q_k: [f64; 3], q_kp1: [f64; 3]) -> [f64; 3] {
        let d1 = self.lagrangian.d1(q_k, q_kp1, self.dt);
        scale3(d1, -1.0)
    }
    /// Right discrete Legendre transform: p_{k+1}^+ = D2 L_d(q_k, q_{k+1}).
    #[allow(dead_code)]
    pub fn right_transform(&self, q_k: [f64; 3], q_kp1: [f64; 3]) -> [f64; 3] {
        self.lagrangian.d2(q_k, q_kp1, self.dt)
    }
    /// Matching condition: at the DEL solution, p_k^+ = p_k^-.
    /// Returns the mismatch.
    #[allow(dead_code)]
    pub fn momentum_mismatch(&self, q_km1: [f64; 3], q_k: [f64; 3], q_kp1: [f64; 3]) -> [f64; 3] {
        let p_plus = self.right_transform(q_km1, q_k);
        let p_minus = self.left_transform(q_k, q_kp1);
        sub3(p_plus, p_minus)
    }
}
/// Variational midpoint integrator: a second-order symplectic method.
///
/// Uses the midpoint discrete Lagrangian and solves the resulting
/// discrete Euler-Lagrange equations.
#[derive(Debug, Clone)]
pub struct VariationalMidpoint {
    /// Mass.
    pub mass: f64,
    /// Time step.
    pub dt: f64,
}
impl VariationalMidpoint {
    /// Create a new variational midpoint integrator.
    pub fn new(mass: f64, dt: f64) -> Self {
        Self { mass, dt }
    }
    /// Step for a particle in a potential field.
    ///
    /// Given q_k and v_k, compute q_{k+1} and v_{k+1}.
    /// Uses the midpoint rule: force evaluated at midpoint.
    #[allow(dead_code)]
    pub fn step(
        &self,
        q: [f64; 3],
        v: [f64; 3],
        force_fn: &dyn Fn([f64; 3]) -> [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let q_mid = add3(q, scale3(v, 0.5 * self.dt));
        let f_mid = force_fn(q_mid);
        let a_mid = scale3(f_mid, 1.0 / self.mass);
        let v_new = add3(v, scale3(a_mid, self.dt));
        let q_new = add3(q, scale3(add3(v, v_new), 0.5 * self.dt));
        (q_new, v_new)
    }
    /// Compute the discrete action over a trajectory.
    #[allow(dead_code)]
    pub fn discrete_action(
        &self,
        trajectory: &[[f64; 3]],
        potential_fn: &dyn Fn([f64; 3]) -> f64,
    ) -> f64 {
        let mut action = 0.0;
        for i in 0..trajectory.len() - 1 {
            let q_k = trajectory[i];
            let q_kp1 = trajectory[i + 1];
            let v = scale3(sub3(q_kp1, q_k), 1.0 / self.dt);
            let q_mid = scale3(add3(q_k, q_kp1), 0.5);
            let kinetic = 0.5 * self.mass * dot3(v, v);
            let potential = potential_fn(q_mid);
            action += self.dt * (kinetic - potential);
        }
        action
    }
}
/// Constrained variational integrator combining DEL equations with
/// holonomic constraints.
///
/// Solves the augmented system:
///   D2 L_d(q_{k-1}, q_k) + D1 L_d(q_k, q_{k+1}) + G^T lambda = 0
///   g(q_{k+1}) = 0
#[derive(Debug, Clone)]
pub struct ConstrainedVariationalIntegrator {
    /// Mass of particles.
    pub masses: Vec<f64>,
    /// Time step.
    pub dt: f64,
    /// Constraints.
    pub constraints: Vec<HolonomicConstraint>,
    /// Solver tolerance.
    pub tol: f64,
    /// Maximum iterations.
    pub max_iter: usize,
}
impl ConstrainedVariationalIntegrator {
    /// Create a new constrained variational integrator.
    pub fn new(masses: Vec<f64>, dt: f64, constraints: Vec<HolonomicConstraint>) -> Self {
        Self {
            masses,
            dt,
            constraints,
            tol: 1e-8,
            max_iter: 100,
        }
    }
    /// Compute unconstrained update for particle i: q_{k+1} = 2*q_k - q_{k-1} + h^2*f/m.
    #[allow(dead_code)]
    fn unconstrained_step(
        &self,
        q_km1: [f64; 3],
        q_k: [f64; 3],
        force: [f64; 3],
        mass: f64,
    ) -> [f64; 3] {
        let dt2 = self.dt * self.dt;
        let accel = scale3(force, dt2 / mass);
        add3(sub3(scale3(q_k, 2.0), q_km1), accel)
    }
    /// Perform one step with constraint projection (SHAKE-like).
    ///
    /// 1. Compute unconstrained update.
    /// 2. Project onto constraint manifold iteratively.
    #[allow(dead_code)]
    pub fn step(
        &self,
        positions_km1: &[[f64; 3]],
        positions_k: &[[f64; 3]],
        forces: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let n = self.masses.len();
        let mut q_new: Vec<[f64; 3]> = Vec::with_capacity(n);
        for i in 0..n {
            q_new.push(self.unconstrained_step(
                positions_km1[i],
                positions_k[i],
                forces[i],
                self.masses[i],
            ));
        }
        for _ in 0..self.max_iter {
            let mut max_error = 0.0_f64;
            for constraint in &self.constraints {
                let error = constraint.evaluate(&q_new);
                max_error = max_error.max(error.abs());
                if error.abs() < self.tol {
                    continue;
                }
                let jac = constraint.jacobian(&q_new);
                let mut w_sum = 0.0;
                for &(idx, grad) in &jac {
                    w_sum += dot3(grad, grad) / self.masses[idx];
                }
                if w_sum < 1e-30 {
                    continue;
                }
                let lambda = -error / w_sum;
                for &(idx, grad) in &jac {
                    let correction = scale3(grad, lambda / self.masses[idx]);
                    q_new[idx] = add3(q_new[idx], correction);
                }
            }
            if max_error < self.tol {
                break;
            }
        }
        q_new
    }
    /// Compute velocity from position difference: v_k = (q_{k+1} - q_{k-1}) / (2*h).
    #[allow(dead_code)]
    pub fn velocity(&self, q_km1: [f64; 3], q_kp1: [f64; 3]) -> [f64; 3] {
        scale3(sub3(q_kp1, q_km1), 1.0 / (2.0 * self.dt))
    }
}
/// Energy monitor for tracking conservation properties of variational integrators.
#[derive(Debug, Clone)]
pub struct EnergyMonitor {
    /// Recorded energy values.
    pub energies: Vec<f64>,
    /// Recorded momentum magnitudes.
    pub momenta: Vec<f64>,
}
impl EnergyMonitor {
    /// Create a new empty energy monitor.
    pub fn new() -> Self {
        Self {
            energies: Vec::new(),
            momenta: Vec::new(),
        }
    }
    /// Record an energy value.
    #[allow(dead_code)]
    pub fn record_energy(&mut self, e: f64) {
        self.energies.push(e);
    }
    /// Record a momentum magnitude.
    #[allow(dead_code)]
    pub fn record_momentum(&mut self, p_mag: f64) {
        self.momenta.push(p_mag);
    }
    /// Compute the maximum energy drift.
    #[allow(dead_code)]
    pub fn max_energy_drift(&self) -> f64 {
        if self.energies.is_empty() {
            return 0.0;
        }
        let e0 = self.energies[0];
        self.energies
            .iter()
            .map(|&e| (e - e0).abs())
            .fold(0.0_f64, f64::max)
    }
    /// Compute the relative energy error.
    #[allow(dead_code)]
    pub fn relative_energy_error(&self) -> f64 {
        if self.energies.is_empty() {
            return 0.0;
        }
        let e0 = self.energies[0];
        if e0.abs() < 1e-30 {
            return self.max_energy_drift();
        }
        self.max_energy_drift() / e0.abs()
    }
    /// Compute the maximum momentum drift.
    #[allow(dead_code)]
    pub fn max_momentum_drift(&self) -> f64 {
        if self.momenta.is_empty() {
            return 0.0;
        }
        let p0 = self.momenta[0];
        self.momenta
            .iter()
            .map(|&p| (p - p0).abs())
            .fold(0.0_f64, f64::max)
    }
    /// Check if energy is bounded (not growing).
    #[allow(dead_code)]
    pub fn is_energy_bounded(&self, tol: f64) -> bool {
        self.max_energy_drift() < tol
    }
}
/// Backward error analysis for symplectic integrators.
///
/// A symplectic integrator exactly solves a modified Hamiltonian
/// H_mod = H + h^p * H_p + h^{p+1} * H_{p+1} + ...
///
/// This structure estimates the modified Hamiltonian error.
#[derive(Debug, Clone)]
pub struct BackwardErrorAnalysis {
    /// Time step.
    pub dt: f64,
    /// Order of the integrator.
    pub order: usize,
}
impl BackwardErrorAnalysis {
    /// Create a new backward error analyzer.
    pub fn new(dt: f64, order: usize) -> Self {
        Self { dt, order }
    }
    /// Estimate the modified Hamiltonian energy for a harmonic oscillator.
    ///
    /// For H = p^2/(2m) + k*q^2/2, the modified Hamiltonian of the
    /// symplectic Euler method is:
    ///   H_mod = H + (h/2) * k * p * q / m + O(h^2)
    #[allow(dead_code)]
    pub fn modified_hamiltonian_harmonic(&self, q: f64, p: f64, mass: f64, stiffness: f64) -> f64 {
        let h_original = p * p / (2.0 * mass) + 0.5 * stiffness * q * q;
        let h1_correction = (self.dt / 2.0) * stiffness * p * q / mass;
        h_original + h1_correction
    }
    /// Estimate the energy error bound for a p-th order symplectic integrator.
    ///
    /// |H(q_n, p_n) - H(q_0, p_0)| <= C * h^p * t * exp(gamma * t)
    ///
    /// For exponentially long times, the error grows very slowly.
    #[allow(dead_code)]
    pub fn energy_error_bound(&self, time: f64, _energy_scale: f64, gamma: f64) -> f64 {
        let h_p = self.dt.powi(self.order as i32);
        h_p * time * (gamma * time).exp()
    }
    /// Compute the shadow Hamiltonian error: H_mod - H.
    ///
    /// Approximated by tracking energy drift over a trajectory.
    #[allow(dead_code)]
    pub fn shadow_hamiltonian_drift(&self, energies: &[f64]) -> f64 {
        if energies.len() < 2 {
            return 0.0;
        }
        let e0 = energies[0];
        let mut max_drift = 0.0_f64;
        for &e in energies.iter().skip(1) {
            max_drift = max_drift.max((e - e0).abs());
        }
        max_drift
    }
    /// Estimate the effective order of the integrator from energy data.
    ///
    /// If we halve the time step and the energy error decreases by a factor
    /// of 2^p, then the order is p.
    #[allow(dead_code)]
    pub fn estimate_order(&self, error_h: f64, error_h_half: f64) -> f64 {
        if error_h_half.abs() < 1e-30 || error_h.abs() < 1e-30 {
            return 0.0;
        }
        (error_h / error_h_half).ln() / (2.0_f64).ln()
    }
    /// Compute the phase space volume preservation error.
    ///
    /// For a symplectic integrator, the phase space volume should be preserved.
    /// This measures the deviation from volume preservation using the Jacobian determinant.
    #[allow(dead_code)]
    pub fn volume_preservation_error(&self, jacobian_det: f64) -> f64 {
        (jacobian_det - 1.0).abs()
    }
}
/// Discrete Noether's theorem: if the discrete Lagrangian is invariant under
/// a symmetry group action, then the corresponding momentum map is preserved
/// by the discrete flow.
///
/// This structure verifies Noether conservation for discrete systems.
#[derive(Debug, Clone)]
pub struct DiscreteNoether {
    /// Tolerance for conservation check.
    pub tol: f64,
}
impl DiscreteNoether {
    /// Create a new discrete Noether checker.
    pub fn new(tol: f64) -> Self {
        Self { tol }
    }
    /// Check translational invariance of the discrete Lagrangian.
    ///
    /// Evaluates L_d(q_k + eps, q_{k+1} + eps) - L_d(q_k, q_{k+1}) for
    /// a small displacement eps in each direction.
    #[allow(dead_code)]
    pub fn check_translation_invariance(
        &self,
        lagrangian: &DiscreteLagrangian,
        q_k: [f64; 3],
        q_kp1: [f64; 3],
        h: f64,
    ) -> bool {
        let eps_val = 0.01;
        let l_original = lagrangian.evaluate(q_k, q_kp1, h);
        for axis in 0..3 {
            let mut shift = [0.0; 3];
            shift[axis] = eps_val;
            let q_k_shifted = add3(q_k, shift);
            let q_kp1_shifted = add3(q_kp1, shift);
            let l_shifted = lagrangian.evaluate(q_k_shifted, q_kp1_shifted, h);
            if (l_shifted - l_original).abs() > self.tol {
                return false;
            }
        }
        true
    }
    /// Check rotational invariance (about z-axis) of the discrete Lagrangian.
    #[allow(dead_code)]
    pub fn check_rotation_invariance_z(
        &self,
        lagrangian: &DiscreteLagrangian,
        q_k: [f64; 3],
        q_kp1: [f64; 3],
        h: f64,
    ) -> bool {
        let angle: f64 = 0.01;
        let c = angle.cos();
        let s = angle.sin();
        let rotate = |q: [f64; 3]| -> [f64; 3] { [c * q[0] - s * q[1], s * q[0] + c * q[1], q[2]] };
        let l_original = lagrangian.evaluate(q_k, q_kp1, h);
        let l_rotated = lagrangian.evaluate(rotate(q_k), rotate(q_kp1), h);
        (l_rotated - l_original).abs() < self.tol
    }
    /// Verify Noether conservation: given momentum map values at two time steps,
    /// check they are equal to within tolerance.
    #[allow(dead_code)]
    pub fn verify_conservation(&self, momentum_before: [f64; 3], momentum_after: [f64; 3]) -> bool {
        len3(sub3(momentum_before, momentum_after)) < self.tol
    }
    /// Compute the discrete momentum map from the discrete Lagrangian.
    ///
    /// p_k = -D1 L_d(q_k, q_{k+1}, h) (the left discrete Legendre transform).
    #[allow(dead_code)]
    pub fn discrete_momentum(
        &self,
        lagrangian: &DiscreteLagrangian,
        q_k: [f64; 3],
        q_kp1: [f64; 3],
        h: f64,
    ) -> [f64; 3] {
        let d1 = lagrangian.d1(q_k, q_kp1, h);
        scale3(d1, -1.0)
    }
}
