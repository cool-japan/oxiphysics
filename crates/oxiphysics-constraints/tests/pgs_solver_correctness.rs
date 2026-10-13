//! Regression tests for the projected Gauss-Seidel solvers
//! ([`PgsSolver`], [`SorPgsSolver`] and [`solve_block_pgs`]).
//!
//! They pin down four behaviours that used to be fabricated:
//!
//! * the residual and the `converged` flag describe the real state of the
//!   iteration (not a constraint count or a body-speed check) and the
//!   `Constraint` hooks `is_active`, `residual` and `post_solve` are honoured,
//! * warm starting applies `+Jᵀλ` to the first and `-Jᵀλ` to the second body
//!   (momentum conserving, angular part included) and round-trips through
//!   [`PgsState`],
//! * successive over-relaxation scales the impulse update before the
//!   accumulated impulse is projected, instead of scaling body velocities,
//! * the default relaxation of constraints that do not implement it never
//!   rescales absolute velocities.

use std::cell::Cell;
use std::rc::Rc;

use oxiphysics_constraints::pgs_solver::{
    ConstraintBlock, PgsSolver, PgsState, SolverStats, SorPgsSolver, solve_block_pgs,
};
use oxiphysics_constraints::{Constraint, ContactConstraint, DistanceJoint};
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::{RigidBody, RigidBodySet};

const DT: f64 = 1.0 / 60.0;
const TIGHT: f64 = 1e-12;

// ── Scene helpers ────────────────────────────────────────────────────────────

fn dynamic_body(mass: f64, position: Vec3, velocity: Vec3) -> RigidBody {
    let mut body = RigidBody::new(mass);
    body.transform.position = position;
    body.velocity = velocity;
    body
}

fn velocity_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies.get(handle).map(|b| b.velocity).unwrap_or_default()
}

fn angular_velocity_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies
        .get(handle)
        .map(|b| b.angular_velocity)
        .unwrap_or_default()
}

fn total_kinetic_energy(bodies: &RigidBodySet, handles: &[BodyHandle]) -> f64 {
    handles
        .iter()
        .filter_map(|&h| bodies.get(h))
        .map(RigidBody::kinetic_energy)
        .sum()
}

fn assert_vec_close(actual: Vec3, expected: Vec3, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).norm() <= tolerance,
        "{what}: expected {expected:?}, got {actual:?}"
    );
}

/// Contact between a body `above` and the body `below` with a shared contact
/// point and an upward normal.
fn upward_contact(
    above: BodyHandle,
    below: BodyHandle,
    point: Vec3,
    restitution: f64,
    friction: f64,
) -> ContactConstraint {
    ContactConstraint::new(
        above,
        below,
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        point,
        point,
        restitution,
        friction,
    )
}

/// A 2 kg box sliding over static ground.  The contact point lies straight
/// below its centre of mass and the approach speed (0.1 m/s) is below the
/// restitution threshold.
struct SlidingBox {
    bodies: RigidBodySet,
    box_handle: BodyHandle,
    constraints: Vec<Box<dyn Constraint>>,
}

fn sliding_box(slide_speed: f64, friction: f64) -> SlidingBox {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let box_handle = bodies.insert(dynamic_body(
        2.0,
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(slide_speed, -0.1, 0.0),
    ));
    let contact = upward_contact(box_handle, ground, Vec3::zeros(), 0.0, friction);
    SlidingBox {
        bodies,
        box_handle,
        constraints: vec![Box::new(contact)],
    }
}

/// Three unit boxes (unit mass) stacked on static ground, all falling at
/// `fall_speed` as they would after one gravity step.  All contact points lie
/// on the vertical line through the centres of mass, so no angular velocity
/// ever builds up.
struct Stack {
    bodies: RigidBodySet,
    boxes: [BodyHandle; 3],
    constraints: Vec<Box<dyn Constraint>>,
}

fn resting_stack(fall_speed: f64, friction: f64) -> Stack {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let mut boxes = [ground; 3];
    for (level, slot) in boxes.iter_mut().enumerate() {
        *slot = bodies.insert(dynamic_body(
            1.0,
            Vec3::new(0.0, 0.5 + level as f64, 0.0),
            Vec3::new(0.0, -fall_speed, 0.0),
        ));
    }
    let mut constraints: Vec<Box<dyn Constraint>> = Vec::new();
    let mut below = ground;
    for (level, &above) in boxes.iter().enumerate() {
        let point = Vec3::new(0.0, level as f64, 0.0);
        constraints.push(Box::new(upward_contact(above, below, point, 0.0, friction)));
        below = above;
    }
    Stack {
        bodies,
        boxes,
        constraints,
    }
}

// ── Probe constraint ─────────────────────────────────────────────────────────

/// Call counters shared between a [`Probe`] and the test (the solver owns the
/// boxed constraint, so the test observes it through this handle).
#[derive(Default)]
struct Counters {
    prepare: Cell<usize>,
    velocity: Cell<usize>,
    position: Cell<usize>,
    post_solve: Cell<usize>,
}

fn bump(cell: &Cell<usize>) {
    cell.set(cell.get() + 1);
}

/// What the velocity solve of a [`Probe`] does to its body.
#[derive(Clone, Copy)]
enum Behaviour {
    /// Leaves the body alone, like a satisfied constraint.
    Idle,
    /// Adds a fixed amount to the x velocity on every call (never converges).
    Drift(f64),
    /// Scales the x velocity by a factor on every call (converges
    /// geometrically).
    Decay(f64),
}

/// A custom constraint that counts the hook calls it receives.
struct Probe {
    handle: BodyHandle,
    active: bool,
    behaviour: Behaviour,
    /// Value returned by `residual()`; `None` keeps the trait default.
    reported: Option<f64>,
    counters: Rc<Counters>,
}

impl Constraint for Probe {
    fn prepare(&mut self, _bodies: &RigidBodySet, _dt: f64) {
        bump(&self.counters.prepare);
    }

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        bump(&self.counters.velocity);
        if let Some(body) = bodies.get_mut(self.handle) {
            match self.behaviour {
                Behaviour::Idle => {}
                Behaviour::Drift(delta) => body.velocity.x += delta,
                Behaviour::Decay(factor) => body.velocity.x *= factor,
            }
        }
    }

    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {
        bump(&self.counters.position);
    }

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.handle]
    }

    fn is_active(&self) -> bool {
        self.active
    }

    fn post_solve(&mut self) {
        bump(&self.counters.post_solve);
    }

    fn residual(&self) -> f64 {
        self.reported.unwrap_or(0.0)
    }
}

/// One body with x velocity `speed` and a [`Probe`] attached to it.
struct ProbeScene {
    bodies: RigidBodySet,
    body: BodyHandle,
    constraints: Vec<Box<dyn Constraint>>,
    counters: Rc<Counters>,
}

fn probe_scene(
    speed: f64,
    active: bool,
    behaviour: Behaviour,
    reported: Option<f64>,
) -> ProbeScene {
    let mut bodies = RigidBodySet::new();
    let body = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::new(speed, 0.0, 0.0)));
    let counters = Rc::new(Counters::default());
    let probe = Probe {
        handle: body,
        active,
        behaviour,
        reported,
        counters: Rc::clone(&counters),
    };
    ProbeScene {
        bodies,
        body,
        constraints: vec![Box::new(probe)],
        counters,
    }
}

// ── F0483: hooks and honest convergence reporting ────────────────────────────

fn run_pgs(constraints: &mut [Box<dyn Constraint>], bodies: &mut RigidBodySet) -> SolverStats {
    PgsSolver::new(10, 2).solve(constraints, bodies, DT)
}

fn run_sor(constraints: &mut [Box<dyn Constraint>], bodies: &mut RigidBodySet) -> SolverStats {
    SorPgsSolver::new(10, 2, 1.4).solve(constraints, bodies, DT)
}

fn run_warm(constraints: &mut [Box<dyn Constraint>], bodies: &mut RigidBodySet) -> SolverStats {
    let mut state = PgsState::new(constraints.len());
    PgsSolver::new(10, 2).solve_warm(constraints, bodies, DT, &mut state)
}

fn run_sor_warm(constraints: &mut [Box<dyn Constraint>], bodies: &mut RigidBodySet) -> SolverStats {
    let mut state = PgsState::new(constraints.len());
    SorPgsSolver::new(10, 2, 1.4).solve_warm(constraints, bodies, DT, &mut state)
}

fn run_block(constraints: &mut [Box<dyn Constraint>], bodies: &mut RigidBodySet) -> SolverStats {
    let blocks = [ConstraintBlock::new(vec![0], "probe")];
    solve_block_pgs(constraints, &blocks, bodies, DT, 10, 2)
}

type SolveFn = fn(&mut [Box<dyn Constraint>], &mut RigidBodySet) -> SolverStats;

const ALL_SOLVERS: [(&str, SolveFn); 5] = [
    ("PgsSolver::solve", run_pgs),
    ("SorPgsSolver::solve", run_sor),
    ("PgsSolver::solve_warm", run_warm),
    ("SorPgsSolver::solve_warm", run_sor_warm),
    ("solve_block_pgs", run_block),
];

#[test]
fn inactive_constraints_are_skipped_by_every_solver_phase() {
    for (name, solve) in ALL_SOLVERS {
        let mut scene = probe_scene(2.0, false, Behaviour::Drift(1.0), None);
        let stats = solve(&mut scene.constraints, &mut scene.bodies);
        let c = &scene.counters;
        assert_eq!(c.prepare.get(), 0, "{name}: prepare");
        assert_eq!(c.velocity.get(), 0, "{name}: solve_velocity");
        assert_eq!(c.position.get(), 0, "{name}: solve_position");
        assert_eq!(c.post_solve.get(), 0, "{name}: post_solve");
        assert_eq!(
            velocity_of(&scene.bodies, scene.body),
            Vec3::new(2.0, 0.0, 0.0),
            "{name}: an inactive constraint must apply no impulse"
        );
        assert!(stats.converged, "{name}: nothing active means converged");
        assert_eq!(stats.residual, 0.0, "{name}: residual");
    }
}

#[test]
fn active_constraints_receive_each_hook_the_documented_number_of_times() {
    let mut scene = probe_scene(2.0, true, Behaviour::Idle, None);
    let stats = PgsSolver::new(10, 3).solve(&mut scene.constraints, &mut scene.bodies, DT);
    let c = &scene.counters;
    assert_eq!(c.prepare.get(), 1);
    assert_eq!(c.velocity.get(), stats.iterations_used);
    assert_eq!(
        stats.iterations_used, 1,
        "an idle constraint is converged at once"
    );
    assert_eq!(c.position.get(), 3);
    assert_eq!(
        c.post_solve.get(),
        1,
        "post_solve runs once, after the loop"
    );

    let stats = PgsSolver::new(10, 3).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged);
    assert_eq!(c.post_solve.get(), 2, "once per solve");
}

#[test]
fn post_solve_runs_after_the_last_velocity_iteration() {
    // A constraint that never converges uses the whole iteration budget; the
    // hook still runs exactly once.
    let mut scene = probe_scene(0.0, true, Behaviour::Drift(1.0), None);
    let stats = SorPgsSolver::new(6, 0, 1.0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(stats.iterations_used, 6);
    assert_eq!(scene.counters.velocity.get(), 6);
    assert_eq!(scene.counters.post_solve.get(), 1);
}

#[test]
fn unreported_residual_is_measured_from_the_velocity_change() {
    // v <- v / 2 on every call: the update of call k has size 2^-k, so the
    // summed update first drops to the 1e-6 tolerance at k = 20.
    let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), None);
    let stats = PgsSolver::new(100, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged);
    assert_eq!(stats.iterations_used, 20);
    assert!(
        (stats.residual - 0.5_f64.powi(20)).abs() < TIGHT,
        "residual {}",
        stats.residual
    );
    assert!(stats.residual <= 1e-6);
}

#[test]
fn reported_residual_decides_convergence() {
    // The constraint changes nothing but reports an update of 0.25: the
    // solver must believe the constraint instead of the (unchanged) bodies.
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, Some(0.25));
    let stats = PgsSolver::new(7, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(!stats.converged);
    assert_eq!(stats.iterations_used, 7);
    assert_eq!(stats.residual, 0.25);

    // Reporting a residual below the tolerance converges at once.
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, Some(1e-9));
    let stats = PgsSolver::new(7, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged);
    assert_eq!(stats.iterations_used, 1);
    assert_eq!(stats.residual, 1e-9);
}

#[test]
fn residual_sums_the_updates_of_all_active_constraints() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, Some(0.25));
    let second = Rc::new(Counters::default());
    scene.constraints.push(Box::new(Probe {
        handle: scene.body,
        active: true,
        behaviour: Behaviour::Idle,
        reported: Some(0.5),
        counters: Rc::clone(&second),
    }));
    scene.constraints.push(Box::new(Probe {
        handle: scene.body,
        active: false,
        behaviour: Behaviour::Idle,
        reported: Some(100.0),
        counters: Rc::new(Counters::default()),
    }));
    let stats = PgsSolver::new(3, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(
        stats.residual, 0.75,
        "0.25 + 0.5, the inactive one is ignored"
    );
    assert_eq!(stats.iterations_used, 3);
    assert!(!stats.converged);
}

#[test]
fn constraint_that_keeps_changing_velocities_never_converges() {
    let mut scene = probe_scene(0.0, true, Behaviour::Drift(1.0), None);
    let stats = PgsSolver::new(9, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(!stats.converged);
    assert_eq!(stats.iterations_used, 9);
    assert_eq!(stats.residual, 1.0);
    assert_eq!(velocity_of(&scene.bodies, scene.body).x, 9.0);
}

#[test]
fn zero_iteration_budget_reports_not_converged() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, None);
    let stats = PgsSolver::new(0, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(stats.iterations_used, 0);
    assert!(!stats.converged);
    assert_eq!(stats.residual, f64::MAX);
    assert_eq!(scene.counters.post_solve.get(), 1);
}

#[test]
fn sliding_box_converges_while_still_sliding() {
    let mut scene = sliding_box(3.0, 0.3);
    let stats = PgsSolver::new(50, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);

    // The contact reaches its fixed point in two sweeps: the first applies
    // the impulses, the second finds nothing left to change.
    assert!(stats.converged, "stats {stats:?}");
    assert_eq!(stats.iterations_used, 2);
    assert_eq!(stats.residual, 0.0);

    // m = 2 kg, approach speed 0.1 m/s -> lambda_n = 0.2 N s; friction is
    // clamped at the Coulomb limit mu * lambda_n = 0.06 N s.
    let rows = scene.constraints[0].accumulated_impulses();
    assert!((rows[0] - 0.2).abs() < TIGHT, "normal impulse {}", rows[0]);
    let friction = rows[1].hypot(rows[2]);
    assert!(
        (friction - 0.06).abs() < TIGHT,
        "friction impulse {friction} must sit on the Coulomb limit"
    );

    let v = velocity_of(&scene.bodies, scene.box_handle);
    assert!(
        (v.x - (3.0 - 0.03)).abs() < TIGHT,
        "still sliding, vx = {}",
        v.x
    );
    assert!(v.y.abs() < TIGHT, "vy = {}", v.y);
    // Friction acts below the centre of mass: tau = r x J with r = (0, -0.5, 0)
    // and J = (-0.06, 0, 0) gives (0, 0, -0.03); the inertia is 2 kg m^2.
    let w = angular_velocity_of(&scene.bodies, scene.box_handle);
    assert_vec_close(w, Vec3::new(0.0, 0.0, -0.015), TIGHT, "angular velocity");
}

#[test]
fn sliding_box_converges_for_the_sor_and_block_solvers_too() {
    let mut scene = sliding_box(3.0, 0.3);
    let stats = SorPgsSolver::new(50, 0, 1.2).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged, "SOR stats {stats:?}");
    assert!(stats.iterations_used < 50);
    assert!(velocity_of(&scene.bodies, scene.box_handle).x > 2.9);

    let mut scene = sliding_box(3.0, 0.3);
    let blocks = [ConstraintBlock::new(vec![0], "contact")];
    let stats = solve_block_pgs(
        &mut scene.constraints,
        &blocks,
        &mut scene.bodies,
        DT,
        50,
        2,
    );
    assert!(stats.converged, "block stats {stats:?}");
    assert!(stats.iterations_used < 50);
    assert_eq!(stats.residual, 0.0);
}

#[test]
fn sticking_box_converges_inside_the_friction_cone() {
    // With enough friction the contact point sticks: the friction impulse
    // stays strictly inside the Coulomb cone and the tangential velocity of
    // the contact point (not of the centre of mass: the box pivots) vanishes.
    let mut scene = sliding_box(0.01, 5.0);
    let stats = PgsSolver::new(200, 0).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged, "stats {stats:?}");
    assert!(
        stats.iterations_used <= 3,
        "iterations {}",
        stats.iterations_used
    );
    assert!(stats.residual <= 1e-6);

    let rows = scene.constraints[0].accumulated_impulses();
    let friction = rows[1].hypot(rows[2]);
    let coulomb_limit = 5.0 * rows[0];
    assert!(
        friction > 1e-3 && friction < 0.5 * coulomb_limit,
        "friction impulse {friction} must be well inside the cone (limit {coulomb_limit})"
    );
    let lever = Vec3::new(0.0, -0.5, 0.0);
    let v = velocity_of(&scene.bodies, scene.box_handle);
    let w = angular_velocity_of(&scene.bodies, scene.box_handle);
    let contact_velocity = v + w.cross(&lever);
    assert!(
        contact_velocity.norm() < 1e-9,
        "the contact point must stick, its velocity is {contact_velocity:?}"
    );
}

#[test]
fn resting_stack_converges_to_rest() {
    let mut stack = resting_stack(9.81 * DT, 0.5);
    let stats = PgsSolver::new(400, 0).solve(&mut stack.constraints, &mut stack.bodies, DT);
    assert!(stats.converged, "stats {stats:?}");
    assert!(
        stats.iterations_used < 400,
        "iterations {}",
        stats.iterations_used
    );
    assert!(stats.residual <= 1e-6);
    for &handle in &stack.boxes {
        let v = velocity_of(&stack.bodies, handle);
        assert!(v.norm() < 1e-5, "box must come to rest, v = {v:?}");
    }
    // The impulses balance the weight of everything above each contact:
    // lambda = (3, 2, 1) * m * g * dt.
    let expected = [3.0, 2.0, 1.0].map(|n| n * 9.81 * DT);
    for (constraint, expected) in stack.constraints.iter().zip(expected) {
        let lambda = constraint.accumulated_impulses()[0];
        assert!(
            (lambda - expected).abs() < 1e-4,
            "normal impulse {lambda}, expected {expected}"
        );
    }
}

#[test]
fn small_iteration_budget_on_a_stack_is_reported_as_unconverged() {
    let mut stack = resting_stack(9.81 * DT, 0.5);
    let stats = PgsSolver::new(2, 0).solve(&mut stack.constraints, &mut stack.bodies, DT);
    assert!(
        !stats.converged,
        "two sweeps cannot settle a three box stack"
    );
    assert_eq!(stats.iterations_used, 2);
    assert!(stats.residual > 1e-6);
}

// ── F0484: warm starting ─────────────────────────────────────────────────────

/// Two dynamic bodies touching at an off-axis point with an oblique normal,
/// so that the linear and the angular parts of the warm start are exercised.
struct PairScene {
    bodies: RigidBodySet,
    a: BodyHandle,
    b: BodyHandle,
    contact_point: Vec3,
    constraint: ContactConstraint,
}

fn oblique_pair() -> PairScene {
    let mut bodies = RigidBodySet::new();
    let mut body_a = dynamic_body(2.0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(1.0, 0.5, -0.2));
    body_a.angular_velocity = Vec3::new(0.1, 0.2, 0.3);
    let mut body_b = dynamic_body(3.0, Vec3::new(0.5, 0.0, 0.5), Vec3::new(-0.3, 0.2, 0.4));
    body_b.angular_velocity = Vec3::new(0.0, -0.2, 0.1);
    let a = bodies.insert(body_a);
    let b = bodies.insert(body_b);
    let contact_point = Vec3::new(0.3, 0.5, -0.2);
    let constraint = ContactConstraint::new(
        a,
        b,
        Vec3::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
        0.0,
        contact_point,
        contact_point,
        0.0,
        0.5,
    );
    PairScene {
        bodies,
        a,
        b,
        contact_point,
        constraint,
    }
}

/// Velocities `(v_a, w_a, v_b, w_b)` of two bodies.
fn pair_velocities(
    bodies: &RigidBodySet,
    a: BodyHandle,
    b: BodyHandle,
) -> (Vec3, Vec3, Vec3, Vec3) {
    (
        velocity_of(bodies, a),
        angular_velocity_of(bodies, a),
        velocity_of(bodies, b),
        angular_velocity_of(bodies, b),
    )
}

/// Total linear momentum and angular momentum about the world origin.
fn momenta(bodies: &RigidBodySet, handles: &[BodyHandle]) -> (Vec3, Vec3) {
    let mut linear = Vec3::zeros();
    let mut angular = Vec3::zeros();
    for &handle in handles {
        if let Some(body) = bodies.get(handle) {
            let momentum = body.velocity * body.mass;
            linear += momentum;
            angular += body.transform.position.cross(&momentum)
                + body.world_inertia() * body.angular_velocity;
        }
    }
    (linear, angular)
}

#[test]
fn warm_start_applies_opposite_linear_and_angular_impulses() {
    let PairScene {
        mut bodies,
        a,
        b,
        contact_point,
        constraint,
    } = oblique_pair();
    let (normal, tangent1, tangent2) =
        (constraint.normal, constraint.tangent1, constraint.tangent2);
    let (lambda_n, lambda_t1, lambda_t2) = (0.7, 0.1, -0.05);

    let before = pair_velocities(&bodies, a, b);
    let momenta_before = momenta(&bodies, &[a, b]);
    let (inv_mass_a, inv_inertia_a, pos_a) = {
        let body = bodies.get(a).expect("body a");
        (
            body.inverse_mass,
            body.world_inverse_inertia,
            body.transform.position,
        )
    };
    let (inv_mass_b, inv_inertia_b, pos_b) = {
        let body = bodies.get(b).expect("body b");
        (
            body.inverse_mass,
            body.world_inverse_inertia,
            body.transform.position,
        )
    };

    let mut state = PgsState::new(1);
    state.record(0, vec![lambda_n, lambda_t1, lambda_t2]);
    state.prev_dt = DT;
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(constraint)];
    // No velocity or position iterations: the result is the warm start alone.
    let solver = PgsSolver::new(0, 0);
    solver.solve_warm(&mut constraints, &mut bodies, DT, &mut state);

    let impulse = normal * lambda_n + tangent1 * lambda_t1 + tangent2 * lambda_t2;
    let r_a = contact_point - pos_a;
    let r_b = contact_point - pos_b;
    let (v_a0, w_a0, v_b0, w_b0) = before;
    assert_vec_close(
        velocity_of(&bodies, a),
        v_a0 + impulse * inv_mass_a,
        TIGHT,
        "linear velocity of A gets +J",
    );
    assert_vec_close(
        velocity_of(&bodies, b),
        v_b0 - impulse * inv_mass_b,
        TIGHT,
        "linear velocity of B gets -J",
    );
    assert_vec_close(
        angular_velocity_of(&bodies, a),
        w_a0 + inv_inertia_a * r_a.cross(&impulse),
        TIGHT,
        "angular velocity of A gets +r x J",
    );
    assert_vec_close(
        angular_velocity_of(&bodies, b),
        w_b0 - inv_inertia_b * r_b.cross(&impulse),
        TIGHT,
        "angular velocity of B gets -r x J",
    );

    // The warm start is an internal impulse pair: both momenta are conserved.
    let momenta_after = momenta(&bodies, &[a, b]);
    assert_vec_close(momenta_after.0, momenta_before.0, TIGHT, "linear momentum");
    assert_vec_close(momenta_after.1, momenta_before.1, TIGHT, "angular momentum");

    // The state now holds exactly what the constraint accumulated.
    assert_eq!(state.rows(0), &[lambda_n, lambda_t1, lambda_t2]);
    assert_eq!(
        constraints[0].accumulated_impulses(),
        vec![lambda_n, lambda_t1, lambda_t2]
    );
}

#[test]
fn warm_start_scales_the_cached_impulses_by_the_time_step_ratio() {
    let initial_velocity = {
        let scene = oblique_pair();
        velocity_of(&scene.bodies, scene.a)
    };
    let run = |prev_dt: f64| {
        let PairScene {
            mut bodies,
            a,
            constraint,
            ..
        } = oblique_pair();
        let mut state = PgsState::new(1);
        state.record(0, vec![0.4, 0.05, 0.02]);
        state.prev_dt = prev_dt;
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(constraint)];
        PgsSolver::new(0, 0).solve_warm(&mut constraints, &mut bodies, DT, &mut state);
        (velocity_of(&bodies, a), state.rows(0).to_vec())
    };
    let (v_same, rows_same) = run(DT);
    let (v_half, rows_half) = run(DT / 2.0);
    assert_eq!(
        rows_same,
        vec![0.4, 0.05, 0.02],
        "equal time steps keep the impulses"
    );
    // prev_dt = dt / 2 doubles the impulses; the velocity change doubles with it.
    assert_vec_close(
        v_half - initial_velocity,
        (v_same - initial_velocity) * 2.0,
        TIGHT,
        "velocity change for dt ratio 2",
    );
    for (half, same) in rows_half.iter().zip(&rows_same) {
        assert!((half - 2.0 * same).abs() < TIGHT, "rows {rows_half:?}");
    }
}

#[test]
fn disabled_warm_starting_ignores_the_state_but_still_records_it() {
    let PairScene {
        mut bodies,
        a,
        b,
        constraint,
        ..
    } = oblique_pair();
    let before = pair_velocities(&bodies, a, b);
    let mut state = PgsState::new(1);
    state.record(0, vec![0.7, 0.1, -0.05]);
    state.prev_dt = DT;
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(constraint)];
    let solver = PgsSolver::with_options(0, 0, 1e-6, false);
    solver.solve_warm(&mut constraints, &mut bodies, DT, &mut state);
    assert_eq!(
        pair_velocities(&bodies, a, b),
        before,
        "nothing may be applied"
    );
    assert_eq!(state.rows(0), &[0.0, 0.0, 0.0]);
    assert_eq!(state.lambdas[0], 0.0);
    assert_eq!(state.cached_impulses[0], [0.0; 3]);
}

#[test]
fn warm_start_is_not_applied_twice() {
    // Frame 1 from a cold start.
    let mut first = sliding_box(0.0, 0.5);
    let mut state = PgsState::new(1);
    let stats =
        PgsSolver::new(20, 0).solve_warm(&mut first.constraints, &mut first.bodies, DT, &mut state);
    assert!(stats.converged);
    assert!(
        (state.rows(0)[0] - 0.2).abs() < TIGHT,
        "rows {:?}",
        state.rows(0)
    );
    assert!((state.lambdas[0] - 0.2).abs() < TIGHT);
    assert_eq!(
        state.rows(0),
        first.constraints[0].accumulated_impulses().as_slice()
    );

    // Frame 2: the same incoming velocity and a brand-new constraint (so its
    // own accumulators are zero) seeded from the state.  The warm start alone
    // already stops the box; the very first sweep then has nothing to change.
    let mut second = sliding_box(0.0, 0.5);
    let stats = PgsSolver::new(20, 0).solve_warm(
        &mut second.constraints,
        &mut second.bodies,
        DT,
        &mut state,
    );
    assert!(stats.converged);
    assert_eq!(
        stats.iterations_used, 1,
        "a double applied warm start would leave work for a second sweep"
    );
    assert_eq!(stats.residual, 0.0);
    assert!(velocity_of(&second.bodies, second.box_handle).norm() < TIGHT);

    // A solver without warm starting needs the full two sweeps again.
    let mut cold = sliding_box(0.0, 0.5);
    let mut stale = state.clone();
    let stats = PgsSolver::with_options(20, 0, 1e-6, false).solve_warm(
        &mut cold.constraints,
        &mut cold.bodies,
        DT,
        &mut stale,
    );
    assert_eq!(stats.iterations_used, 2);
}

#[test]
fn pgs_state_grows_to_the_number_of_constraints_and_keeps_rows_apart() {
    let mut stack = resting_stack(9.81 * DT, 0.5);
    let mut state = PgsState::new(0);
    let stats = PgsSolver::new(400, 0).solve_warm(
        &mut stack.constraints,
        &mut stack.bodies,
        DT,
        &mut state,
    );
    assert!(stats.converged);
    assert_eq!(state.lambdas.len(), 3);
    assert_eq!(state.row_impulses.len(), 3);
    for (index, constraint) in stack.constraints.iter().enumerate() {
        assert_eq!(
            state.rows(index),
            constraint.accumulated_impulses().as_slice()
        );
        assert!((state.lambdas[index] - state.rows(index)[0]).abs() < 1e-9);
        assert_eq!(state.cached_impulses[index][0], state.rows(index)[0]);
    }
    assert!(state.rows(0)[0] > state.rows(1)[0] && state.rows(1)[0] > state.rows(2)[0]);
}

#[test]
fn pgs_state_record_tolerates_inconsistent_vector_lengths() {
    // `PgsState` has public fields: recording must grow every vector instead
    // of indexing past the end of one that a caller left short.
    let mut state = PgsState::new(1);
    state.row_impulses.clear();
    state.cached_impulses.clear();
    state.record(2, vec![3.0, 4.0]);
    assert_eq!(state.lambdas.len(), 3);
    assert_eq!(state.cached_impulses.len(), 3);
    assert_eq!(state.row_impulses.len(), 3);
    assert_eq!(state.lambdas[2], 5.0, "Euclidean norm of the rows");
    assert_eq!(state.cached_impulses[2], [3.0, 4.0, 0.0]);
    assert_eq!(state.rows(2), &[3.0, 4.0]);
    assert!(state.rows(0).is_empty());
    assert!(state.rows(7).is_empty(), "unknown constraints have no rows");
}

#[test]
fn pgs_state_blend_mixes_the_impulse_rows_and_pads_missing_rows_with_zero() {
    let mut own = PgsState::new(2);
    own.record(0, vec![2.0, 0.0]);
    own.record(1, vec![1.0]);
    let mut other = PgsState::new(2);
    other.record(0, vec![0.0, 4.0, 8.0]);
    other.record(1, vec![]);
    own.blend(&other, 0.25);
    assert_eq!(own.rows(0), &[1.5, 1.0, 2.0]);
    assert_eq!(own.rows(1), &[0.75]);
    // alpha is clamped to [0, 1]: 0 keeps the state, anything above 1 takes the target.
    let mut kept = own.clone();
    kept.blend(&other, -3.0);
    assert_eq!(kept.rows(0), own.rows(0));
    own.blend(&other, 7.0);
    assert_eq!(own.rows(0), &[0.0, 4.0, 8.0]);
    assert_eq!(own.rows(1), &[0.0]);
}

#[test]
fn warm_started_stack_needs_fewer_sweeps_than_a_cold_stack() {
    let fall = 9.81 * DT;
    let mut state = PgsState::new(3);
    let mut frame_one = resting_stack(fall, 0.5);
    let cold = PgsSolver::new(400, 0).solve_warm(
        &mut frame_one.constraints,
        &mut frame_one.bodies,
        DT,
        &mut state,
    );
    assert!(cold.converged);

    let mut frame_two = resting_stack(fall, 0.5);
    let warm = PgsSolver::new(400, 0).solve_warm(
        &mut frame_two.constraints,
        &mut frame_two.bodies,
        DT,
        &mut state,
    );
    assert!(warm.converged);
    assert!(
        warm.iterations_used < cold.iterations_used,
        "warm {} vs cold {}",
        warm.iterations_used,
        cold.iterations_used
    );
    assert!(warm.iterations_used <= 2);
}

#[test]
fn constraints_without_impulse_rows_leave_the_warm_start_state_empty() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, None);
    let mut state = PgsState::new(1);
    state.record(0, vec![9.0, 9.0, 9.0]);
    state.prev_dt = DT;
    let stats =
        PgsSolver::new(5, 0).solve_warm(&mut scene.constraints, &mut scene.bodies, DT, &mut state);
    assert!(stats.converged);
    assert!(state.rows(0).is_empty(), "the default hook has no rows");
    assert_eq!(state.lambdas[0], 0.0);
    assert_eq!(
        velocity_of(&scene.bodies, scene.body),
        Vec3::new(1.0, 0.0, 0.0),
        "the default warm start applies nothing"
    );
}

#[test]
fn contact_warm_start_clamps_to_the_feasible_set() {
    let PairScene {
        mut bodies,
        a,
        b,
        constraint,
        ..
    } = oblique_pair();
    let before = pair_velocities(&bodies, a, b);
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(constraint)];
    // A pulling normal impulse is infeasible; friction beyond mu * lambda_n is
    // clamped (mu = 0.5, lambda_n = 0 -> no friction either).
    constraints[0].apply_warm_start(&mut bodies, &[-3.0, 2.0, -2.0]);
    assert_eq!(constraints[0].accumulated_impulses(), vec![0.0, 0.0, 0.0]);
    assert_eq!(pair_velocities(&bodies, a, b), before);

    constraints[0].apply_warm_start(&mut bodies, &[1.0, 2.0, -2.0]);
    assert_eq!(constraints[0].accumulated_impulses(), vec![1.0, 0.5, -0.5]);
}

// ── F0486 / F0487: successive over-relaxation ────────────────────────────────

/// Kinetic energy of the stack after `k` sweeps of SOR with factor `omega`.
fn stack_energy_after(sweeps: usize, omega: f64) -> f64 {
    let mut stack = resting_stack(9.81 * DT, 0.0);
    let mut solver = SorPgsSolver::new(sweeps, 0, omega);
    // A negative tolerance disables the early exit: exactly `sweeps` sweeps run.
    solver.tolerance = -1.0;
    solver.solve(&mut stack.constraints, &mut stack.bodies, DT);
    total_kinetic_energy(&stack.bodies, &stack.boxes)
}

#[test]
fn over_relaxation_never_increases_the_energy_of_a_resting_stack() {
    let initial = stack_energy_after(0, 1.5);
    assert!((initial - 1.5 * (9.81 * DT).powi(2)).abs() < TIGHT);
    let mut previous = initial;
    for sweeps in 1..=40 {
        let energy = stack_energy_after(sweeps, 1.5);
        assert!(
            energy <= previous + 1e-15,
            "sweep {sweeps}: energy rose from {previous} to {energy}"
        );
        assert!(
            energy <= initial,
            "sweep {sweeps}: energy {energy} above the initial one"
        );
        previous = energy;
    }
    assert!(
        previous < 1e-12 * initial,
        "the stack must settle, energy {previous} after 40 sweeps"
    );
}

#[test]
fn every_relaxation_factor_is_energy_stable_on_the_stack() {
    for omega in [0.3, 0.7, 1.0, 1.3, 1.7, 1.95] {
        let initial = stack_energy_after(0, omega);
        let mut previous = initial;
        for sweeps in 1..=25 {
            let energy = stack_energy_after(sweeps, omega);
            assert!(
                energy <= previous + 1e-15,
                "omega {omega}, sweep {sweeps}: energy rose from {previous} to {energy}"
            );
            previous = energy;
        }
    }
}

#[test]
fn over_relaxed_contact_follows_the_sor_recurrence() {
    // One unit mass falling onto static ground at 1 m/s.  The Gauss-Seidel
    // update removes the whole approach speed; SOR with omega scales that
    // update, so the speed error is multiplied by (1 - omega) every sweep and
    // the accumulated impulse approaches m * 1 m/s.
    let omega = 1.5;
    for sweeps in 1..=8 {
        let mut bodies = RigidBodySet::new();
        let ground = bodies.insert(RigidBody::new_static());
        let ball = bodies.insert(dynamic_body(
            1.0,
            Vec3::new(0.0, 0.5, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        ));
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(upward_contact(
            ball,
            ground,
            Vec3::zeros(),
            0.0,
            0.0,
        ))];
        let mut solver = SorPgsSolver::new(sweeps, 0, omega);
        solver.tolerance = -1.0;
        solver.solve(&mut constraints, &mut bodies, DT);

        let expected_vy = -(1.0 - omega).powi(sweeps as i32);
        let vy = velocity_of(&bodies, ball).y;
        assert!(
            (vy - expected_vy).abs() < TIGHT,
            "{sweeps} sweeps: vy = {vy}, expected {expected_vy}"
        );
        // Unit mass: the accumulated impulse is the velocity change.
        let lambda = constraints[0].accumulated_impulses()[0];
        assert!(
            (lambda - (vy + 1.0)).abs() < TIGHT,
            "{sweeps} sweeps: impulse {lambda} must equal the velocity change {}",
            vy + 1.0
        );
        assert!(lambda >= 0.0, "the accumulated impulse stays feasible");
    }
}

#[test]
fn over_relaxation_does_not_pull_on_a_separating_contact() {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let ball = bodies.insert(dynamic_body(
        1.0,
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(0.0, 2.0, 0.0),
    ));
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(upward_contact(
        ball,
        ground,
        Vec3::zeros(),
        0.0,
        0.0,
    ))];
    let stats = SorPgsSolver::new(10, 0, 1.5).solve(&mut constraints, &mut bodies, DT);
    assert!(stats.converged);
    assert_eq!(velocity_of(&bodies, ball), Vec3::new(0.0, 2.0, 0.0));
    assert_eq!(constraints[0].accumulated_impulses()[0], 0.0);
}

#[test]
fn over_relaxed_update_is_projected_onto_the_normal_impulse_bound() {
    // A seeded impulse of 1 pushes the ball away at 1 m/s.  The relaxed update
    // wants to remove 1.5 N s, but only 1 N s are available to take back.
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let ball = bodies.insert(dynamic_body(1.0, Vec3::new(0.0, 0.5, 0.0), Vec3::zeros()));
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(upward_contact(
        ball,
        ground,
        Vec3::zeros(),
        0.0,
        0.0,
    ))];
    let mut state = PgsState::new(1);
    state.record(0, vec![1.0, 0.0, 0.0]);
    state.prev_dt = DT;
    let mut solver = SorPgsSolver::new(1, 0, 1.5);
    solver.tolerance = -1.0;
    solver.solve_warm(&mut constraints, &mut bodies, DT, &mut state);
    assert_eq!(constraints[0].accumulated_impulses()[0], 0.0);
    assert!(velocity_of(&bodies, ball).y.abs() < TIGHT);
    assert_eq!(state.rows(0)[0], 0.0);
}

#[test]
fn contact_relaxed_solve_reports_the_size_of_its_update() {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let ball = bodies.insert(dynamic_body(
        1.0,
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
    ));
    let mut contact = upward_contact(ball, ground, Vec3::zeros(), 0.0, 0.0);
    contact.prepare(&bodies, DT);
    let update = contact.solve_velocity_relaxed(&mut bodies, DT, 1.5);
    assert!((update - 1.5).abs() < TIGHT, "update {update}");
    assert_eq!(contact.residual(), update);
    assert!((velocity_of(&bodies, ball).y - 0.5).abs() < TIGHT);
}

#[test]
fn contact_relaxation_factor_is_sanitised() {
    let solve_once = |omega: f64| {
        let mut bodies = RigidBodySet::new();
        let ground = bodies.insert(RigidBody::new_static());
        let ball = bodies.insert(dynamic_body(
            1.0,
            Vec3::new(0.0, 0.5, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        ));
        let mut contact = upward_contact(ball, ground, Vec3::zeros(), 0.0, 0.0);
        contact.prepare(&bodies, DT);
        contact.solve_velocity_relaxed(&mut bodies, DT, omega);
        velocity_of(&bodies, ball).y
    };
    // omega >= 2 would diverge: it is limited to 1.99.
    assert!((solve_once(5.0) - 0.99).abs() < TIGHT);
    // NaN falls back to the plain Gauss-Seidel update.
    assert!(solve_once(f64::NAN).abs() < TIGHT);
    // A tiny factor is raised to 0.01.
    assert!((solve_once(-3.0) - (-0.99)).abs() < TIGHT);
}

#[test]
fn sor_leaves_a_satisfied_joint_scene_untouched() {
    // Two bodies fly at 5 m/s.  The distance joint (anchors 2 m apart, allowed
    // range 1 m to 3 m) is satisfied.  The old "SOR" multiplied every
    // velocity by omega after every sweep: 1.5^10 = 57 times faster.
    let mut bodies = RigidBodySet::new();
    let velocity = Vec3::new(5.0, 0.0, 0.0);
    let a = bodies.insert(dynamic_body(1.0, Vec3::zeros(), velocity));
    let b = bodies.insert(dynamic_body(2.0, Vec3::new(2.0, 0.0, 0.0), velocity));
    let joint = DistanceJoint::new(a, b, Vec3::zeros(), Vec3::zeros(), 1.0, 3.0);
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(joint)];
    let stats = SorPgsSolver::new(10, 0, 1.5).solve(&mut constraints, &mut bodies, DT);
    assert!(stats.converged);
    assert_vec_close(velocity_of(&bodies, a), velocity, TIGHT, "body a");
    assert_vec_close(velocity_of(&bodies, b), velocity, TIGHT, "body b");
}

#[test]
fn sor_keeps_the_common_velocity_of_a_violated_joint() {
    // Same bodies, but the joint is violated (5 m apart, at most 3 m).  Both
    // bodies move together, so there is no relative velocity to remove and the
    // common 5 m/s must survive for any omega.
    for omega in [0.5, 1.0, 1.5, 1.9] {
        let mut bodies = RigidBodySet::new();
        let velocity = Vec3::new(5.0, 0.0, 0.0);
        let a = bodies.insert(dynamic_body(1.0, Vec3::zeros(), velocity));
        let b = bodies.insert(dynamic_body(2.0, Vec3::new(5.0, 0.0, 0.0), velocity));
        let joint = DistanceJoint::new(a, b, Vec3::zeros(), Vec3::zeros(), 1.0, 3.0);
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(joint)];
        SorPgsSolver::new(10, 0, omega).solve(&mut constraints, &mut bodies, DT);
        assert_vec_close(velocity_of(&bodies, a), velocity, TIGHT, "body a");
        assert_vec_close(velocity_of(&bodies, b), velocity, TIGHT, "body b");
    }
}

/// A scene mixing frictional contacts and a violated distance joint.
fn mixed_scene() -> (RigidBodySet, [BodyHandle; 3], Vec<Box<dyn Constraint>>) {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let first = bodies.insert(dynamic_body(
        1.0,
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(1.0, -0.2, 0.0),
    ));
    let second = bodies.insert(dynamic_body(
        2.0,
        Vec3::new(2.0, 0.5, 0.0),
        Vec3::new(-0.5, -0.2, 0.3),
    ));
    let constraints: Vec<Box<dyn Constraint>> = vec![
        Box::new(upward_contact(
            first,
            ground,
            Vec3::new(0.1, 0.0, 0.0),
            0.0,
            0.4,
        )),
        Box::new(upward_contact(
            second,
            ground,
            Vec3::new(2.0, 0.0, 0.1),
            0.0,
            0.4,
        )),
        Box::new(DistanceJoint::new(
            first,
            second,
            Vec3::zeros(),
            Vec3::zeros(),
            0.5,
            1.5,
        )),
    ];
    (bodies, [ground, first, second], constraints)
}

#[test]
fn sor_with_unit_omega_is_bit_for_bit_plain_pgs() {
    let (mut plain_bodies, handles, mut plain_constraints) = mixed_scene();
    let plain = PgsSolver::new(12, 3).solve(&mut plain_constraints, &mut plain_bodies, DT);

    let (mut sor_bodies, _, mut sor_constraints) = mixed_scene();
    let sor = SorPgsSolver::new(12, 3, 1.0).solve(&mut sor_constraints, &mut sor_bodies, DT);

    assert_eq!(plain.iterations_used, sor.iterations_used);
    assert_eq!(plain.converged, sor.converged);
    assert_eq!(plain.residual.to_bits(), sor.residual.to_bits());
    for handle in handles {
        let bits = |v: Vec3| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()];
        assert_eq!(
            bits(velocity_of(&plain_bodies, handle)),
            bits(velocity_of(&sor_bodies, handle))
        );
        assert_eq!(
            bits(angular_velocity_of(&plain_bodies, handle)),
            bits(angular_velocity_of(&sor_bodies, handle))
        );
    }
    for (p, s) in plain_constraints.iter().zip(&sor_constraints) {
        assert_eq!(p.accumulated_impulses(), s.accumulated_impulses());
    }
}

#[test]
fn constraints_without_relaxation_support_are_never_over_relaxed() {
    // The probe decays the x velocity by a factor 1/2 per call.  It does not
    // implement `solve_velocity_relaxed`, so omega > 1 must not amplify it
    // (blending after the fact could break bounds the constraint maintains).
    let after_one_sweep = |omega: f64| {
        let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), None);
        let stats =
            SorPgsSolver::new(1, 0, omega).solve(&mut scene.constraints, &mut scene.bodies, DT);
        (velocity_of(&scene.bodies, scene.body).x, stats.residual)
    };
    for omega in [1.0, 1.5, 1.9] {
        let (vx, residual) = after_one_sweep(omega);
        assert_eq!(vx, 0.5, "omega {omega}: plain update expected");
        assert_eq!(residual, 0.5);
    }
}

/// Bilateral constraint that equalises the x velocities of two bodies with an
/// impulse pair `+P` on body A and `-P` on body B.
struct MatchSpeed {
    a: BodyHandle,
    b: BodyHandle,
}

impl Constraint for MatchSpeed {
    fn prepare(&mut self, _bodies: &RigidBodySet, _dt: f64) {}

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let (Some(body_a), Some(body_b)) = (bodies.get(self.a), bodies.get(self.b)) else {
            return;
        };
        let stiffness = body_a.inverse_mass + body_b.inverse_mass;
        if stiffness <= 0.0 {
            return;
        }
        let impulse = (body_b.velocity.x - body_a.velocity.x) / stiffness;
        if let Some(body) = bodies.get_mut(self.a) {
            body.velocity.x += body.inverse_mass * impulse;
        }
        if let Some(body) = bodies.get_mut(self.b) {
            body.velocity.x -= body.inverse_mass * impulse;
        }
    }

    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.a, self.b]
    }
}

#[test]
fn two_body_constraint_without_a_residual_is_measured_on_both_bodies() {
    // Masses 1 kg and 3 kg, x velocities 0 and 2 m/s: the bilateral constraint
    // brings both to the common velocity 1.5 m/s, a change of +1.5 and -0.5.
    let mut bodies = RigidBodySet::new();
    let a = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::zeros()));
    let b = bodies.insert(dynamic_body(
        3.0,
        Vec3::new(4.0, 0.0, 0.0),
        Vec3::new(2.0, 0.0, 0.0),
    ));
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(MatchSpeed { a, b })];
    let first = PgsSolver::new(1, 0).solve(&mut constraints, &mut bodies, DT);
    assert!(!first.converged, "the first sweep still changes velocities");
    assert_eq!(first.iterations_used, 1);
    assert!(
        (first.residual - 2.5_f64.sqrt()).abs() < TIGHT,
        "update size sqrt(1.5^2 + 0.5^2), got {}",
        first.residual
    );
    assert!((velocity_of(&bodies, a).x - 1.5).abs() < TIGHT);
    assert!((velocity_of(&bodies, b).x - 1.5).abs() < TIGHT);

    // The constraint is now at its fixed point: nothing left to change.
    let second = PgsSolver::new(5, 0).solve(&mut constraints, &mut bodies, DT);
    assert!(second.converged);
    assert_eq!(second.iterations_used, 1);
    assert!(second.residual < TIGHT, "residual {}", second.residual);
}

#[test]
fn under_relaxation_blends_the_velocity_update_of_any_constraint() {
    // v: 1 -> 1/2 for the plain update; omega = 1/2 keeps 1/2 of the change:
    // v = 1 + (1/2 - 1) / 2 = 3/4.
    let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), None);
    let stats = SorPgsSolver::new(1, 0, 0.5).solve(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(velocity_of(&scene.bodies, scene.body).x, 0.75);
    assert_eq!(stats.residual, 0.25);

    // The same blend applied to a two body constraint removes exactly the
    // fraction omega of the relative velocity the plain solve removes, and it
    // conserves momentum.
    let after_one_sweep = |omega: f64| {
        let mut bodies = RigidBodySet::new();
        let a = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::zeros()));
        let b = bodies.insert(dynamic_body(
            3.0,
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ));
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(MatchSpeed { a, b })];
        let mut solver = SorPgsSolver::new(1, 0, omega);
        solver.tolerance = -1.0;
        solver.solve(&mut constraints, &mut bodies, DT);
        let (va, vb) = (velocity_of(&bodies, a).x, velocity_of(&bodies, b).x);
        (vb - va, va + 3.0 * vb)
    };
    for (omega, remaining) in [(1.0, 0.0), (0.5, 1.0), (0.25, 1.5)] {
        let (relative_speed, momentum) = after_one_sweep(omega);
        assert!(
            (relative_speed - remaining).abs() < TIGHT,
            "omega {omega}: relative speed {relative_speed}, expected {remaining}"
        );
        assert!((momentum - 6.0).abs() < TIGHT, "momentum {momentum}");
    }
    // omega > 1 is not applied to a constraint that cannot be relaxed exactly.
    let (relative_speed, _) = after_one_sweep(1.5);
    assert!(relative_speed.abs() < TIGHT);
}

#[test]
fn default_relaxed_solve_blends_only_for_factors_between_zero_and_one() {
    // Outside (0, 1) the default is the plain update, whatever the factor.
    for omega in [f64::NAN, 0.0, -2.0, 1.0, 3.0] {
        let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), Some(0.125));
        let reported = scene.constraints[0].solve_velocity_relaxed(&mut scene.bodies, DT, omega);
        assert_eq!(
            velocity_of(&scene.bodies, scene.body).x,
            0.5,
            "omega {omega}"
        );
        assert_eq!(
            reported, 0.125,
            "omega {omega}: the constraint's own residual"
        );
    }
    // Inside (0, 1): v = 1 + 0.25 (1/2 - 1) and the reported update shrinks
    // by the same factor.
    let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), Some(0.125));
    let reported = scene.constraints[0].solve_velocity_relaxed(&mut scene.bodies, DT, 0.25);
    assert_eq!(velocity_of(&scene.bodies, scene.body).x, 0.875);
    assert_eq!(reported, 0.03125);
}

#[test]
fn over_relaxation_accelerates_a_resting_stack_and_under_relaxation_slows_it() {
    let sweeps_needed = |omega: f64| {
        let mut stack = resting_stack(9.81 * DT, 0.0);
        let stats =
            SorPgsSolver::new(1000, 0, omega).solve(&mut stack.constraints, &mut stack.bodies, DT);
        assert!(stats.converged, "omega {omega}: {stats:?}");
        stats.iterations_used
    };
    let plain = sweeps_needed(1.0);
    let over = sweeps_needed(1.5);
    let under = sweeps_needed(0.5);
    assert!(
        over * 3 / 2 < plain,
        "SOR(1.5) needs {over} sweeps, plain PGS {plain}"
    );
    assert!(
        under > plain,
        "SOR(0.5) needs {under} sweeps, plain PGS {plain}"
    );
}

#[test]
fn block_pgs_reports_the_real_residual_and_runs_post_solve_once() {
    let mut scene = probe_scene(1.0, true, Behaviour::Decay(0.5), None);
    let blocks = [ConstraintBlock::new(vec![0], "probe")];
    let stats = solve_block_pgs(
        &mut scene.constraints,
        &blocks,
        &mut scene.bodies,
        DT,
        100,
        1,
    );
    assert!(stats.converged);
    assert_eq!(stats.iterations_used, 20);
    assert_eq!(scene.counters.post_solve.get(), 1);
    assert_eq!(scene.counters.prepare.get(), 1);

    let mut scene = probe_scene(0.0, true, Behaviour::Drift(1.0), None);
    let stats = solve_block_pgs(&mut scene.constraints, &blocks, &mut scene.bodies, DT, 4, 3);
    assert!(!stats.converged);
    assert_eq!(stats.iterations_used, 4);
    assert_eq!(stats.residual, 1.0, "residual of the last inner sweep");
    assert_eq!(velocity_of(&scene.bodies, scene.body).x, 12.0);
}
