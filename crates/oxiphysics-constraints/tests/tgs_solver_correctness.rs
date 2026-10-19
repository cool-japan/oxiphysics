//! Regression tests for [`TgsSolver`].
//!
//! They pin down two behaviours that used to be fabricated:
//!
//! * every substep hands the substep size `dt / substeps` (never `0.0`) to
//!   [`Constraint::solve_position`], so dt-dependent position solvers such as
//!   the XPBD joint really correct their positions under TGS,
//! * the returned [`SolverStats`] describe the real state of the iteration
//!   (sweeps run, residual of the worst substep, convergence) instead of a
//!   constant `converged: true`, and the `Constraint` hooks `is_active`,
//!   `residual` and `post_solve` are honoured exactly as the PGS solvers do.

use std::cell::{Cell, RefCell};
use std::rc::Rc;

use oxiphysics_constraints::joints::XpbdJoint;
use oxiphysics_constraints::pgs_solver::SolverStats;
use oxiphysics_constraints::{Constraint, ContactConstraint, TgsSolver};
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::{RigidBody, RigidBodySet};

const DT: f64 = 1.0 / 60.0;

// ── Scene helpers ────────────────────────────────────────────────────────────

fn dynamic_body(mass: f64, position: Vec3, velocity: Vec3) -> RigidBody {
    let mut body = RigidBody::new(mass);
    body.transform.position = position;
    body.velocity = velocity;
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body
}

fn velocity_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies.get(handle).map(|b| b.velocity).unwrap_or_default()
}

fn position_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies
        .get(handle)
        .map(|b| b.transform.position)
        .unwrap_or_default()
}

// ── Probe constraint ─────────────────────────────────────────────────────────

/// Call log shared between a [`Probe`] and the test (the solver owns the boxed
/// constraint, so the test observes it through this handle).
#[derive(Default)]
struct Log {
    prepare_dts: RefCell<Vec<f64>>,
    velocity_dts: RefCell<Vec<f64>>,
    position_dts: RefCell<Vec<f64>>,
    post_solve: Cell<usize>,
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

/// A custom constraint that records the hook calls it receives.
struct Probe {
    handle: BodyHandle,
    active: bool,
    behaviour: Behaviour,
    /// Value returned by `residual()`; `None` keeps the trait default.
    reported: Option<f64>,
    log: Rc<Log>,
}

impl Constraint for Probe {
    fn prepare(&mut self, _bodies: &RigidBodySet, dt: f64) {
        self.log.prepare_dts.borrow_mut().push(dt);
    }

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        self.log.velocity_dts.borrow_mut().push(dt);
        if let Some(body) = bodies.get_mut(self.handle) {
            match self.behaviour {
                Behaviour::Idle => {}
                Behaviour::Drift(delta) => body.velocity.x += delta,
                Behaviour::Decay(factor) => body.velocity.x *= factor,
            }
        }
    }

    fn solve_position(&mut self, _bodies: &mut RigidBodySet, dt: f64) {
        self.log.position_dts.borrow_mut().push(dt);
    }

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.handle]
    }

    fn is_active(&self) -> bool {
        self.active
    }

    fn post_solve(&mut self) {
        self.log.post_solve.set(self.log.post_solve.get() + 1);
    }

    fn residual(&self) -> f64 {
        self.reported.unwrap_or(0.0)
    }
}

/// One gravity-free body with x velocity `speed` and a [`Probe`] attached.
struct ProbeScene {
    bodies: RigidBodySet,
    body: BodyHandle,
    constraints: Vec<Box<dyn Constraint>>,
    log: Rc<Log>,
}

fn probe_scene(
    speed: f64,
    active: bool,
    behaviour: Behaviour,
    reported: Option<f64>,
) -> ProbeScene {
    let mut bodies = RigidBodySet::new();
    let body = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::new(speed, 0.0, 0.0)));
    let log = Rc::new(Log::default());
    let probe = Probe {
        handle: body,
        active,
        behaviour,
        reported,
        log: Rc::clone(&log),
    };
    ProbeScene {
        bodies,
        body,
        constraints: vec![Box::new(probe)],
        log,
    }
}

fn solver(substeps: usize, velocity_iterations: usize) -> TgsSolver {
    TgsSolver::new(substeps, velocity_iterations, Vec3::zeros())
}

// ── F0494: the substep size reaches solve_position ───────────────────────────

/// Two dynamic unit masses whose anchors are `gap` metres too far apart, joined
/// by a distance-1 XPBD joint with the given compliance.
fn stretched_xpbd_pair(
    gap: f64,
    compliance: f64,
) -> (RigidBodySet, [BodyHandle; 2], Vec<Box<dyn Constraint>>) {
    let mut bodies = RigidBodySet::new();
    let a = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::zeros()));
    let b = bodies.insert(dynamic_body(
        1.0,
        Vec3::new(1.0 + gap, 0.0, 0.0),
        Vec3::zeros(),
    ));
    let joint = XpbdJoint::new(a, b, Vec3::zeros(), Vec3::zeros(), 1.0, compliance, 0.0);
    (bodies, [a, b], vec![Box::new(joint)])
}

fn anchor_distance(bodies: &RigidBodySet, pair: [BodyHandle; 2]) -> f64 {
    (position_of(bodies, pair[1]) - position_of(bodies, pair[0])).norm()
}

#[test]
fn rigid_xpbd_joint_closes_its_gap_under_tgs() {
    let (mut bodies, pair, mut constraints) = stretched_xpbd_pair(0.1, 0.0);
    assert!((anchor_distance(&bodies, pair) - 1.1).abs() < 1e-12);

    solver(1, 4).step(&mut constraints, &mut bodies, DT);

    let error = (anchor_distance(&bodies, pair) - 1.0).abs();
    assert!(
        error < 1e-9,
        "a rigid XPBD joint must close a 0.1 m gap in one TGS step, error = {error}"
    );
}

#[test]
fn compliant_xpbd_joint_is_corrected_with_the_substep_size() {
    // XPBD with a constraint value C and unit masses (w = 2): the correction of
    // one pass leaves C * a / (w + a) with a = compliance / sub_dt^2, so the
    // remaining gap pins down the dt the joint received.
    let compliance = 1e-3;
    let gap = 0.1;
    let substeps = 4;
    let sub_dt = DT / substeps as f64;
    let alpha = compliance / (sub_dt * sub_dt);
    let expected = gap * (alpha / (2.0 + alpha)).powi(substeps as i32);

    let (mut bodies, pair, mut constraints) = stretched_xpbd_pair(gap, compliance);
    solver(substeps, 4).step(&mut constraints, &mut bodies, DT);

    let remaining = anchor_distance(&bodies, pair) - 1.0;
    assert!(
        (remaining - expected).abs() < 1e-9,
        "remaining gap {remaining}, expected {expected} for sub_dt = {sub_dt}"
    );
    assert!(remaining < gap, "the joint has to pull the bodies together");
}

#[test]
fn solve_position_receives_the_substep_size_on_every_pass() {
    let mut scene = probe_scene(0.0, true, Behaviour::Idle, None);
    let tgs = TgsSolver::with_options(4, 3, 2, Vec3::zeros(), 0.2, false);
    tgs.step(&mut scene.constraints, &mut scene.bodies, DT);

    let sub_dt = DT / 4.0;
    for (phase, dts) in [
        ("prepare", scene.log.prepare_dts.borrow()),
        ("solve_velocity", scene.log.velocity_dts.borrow()),
        ("solve_position", scene.log.position_dts.borrow()),
    ] {
        assert!(!dts.is_empty(), "{phase} was never called");
        assert!(
            dts.iter().all(|&dt| (dt - sub_dt).abs() < 1e-15),
            "{phase} must receive the substep size {sub_dt}, got {dts:?}"
        );
    }
    assert_eq!(scene.log.prepare_dts.borrow().len(), 4, "one per substep");
    assert_eq!(
        scene.log.position_dts.borrow().len(),
        8,
        "position_iterations (2) passes in each of the 4 substeps"
    );
}

#[test]
fn explicit_position_pass_uses_the_given_step_size() {
    let mut scene = probe_scene(0.0, true, Behaviour::Idle, None);
    let tgs = TgsSolver::with_options(2, 3, 3, Vec3::zeros(), 0.2, false);
    tgs.solve_position_constraints(&mut scene.constraints, &mut scene.bodies, 0.0125);
    assert_eq!(*scene.log.position_dts.borrow(), vec![0.0125; 3]);
}

// ── F0491: honest statistics ─────────────────────────────────────────────────

#[test]
fn satisfied_constraints_converge_after_one_sweep_per_substep() {
    let mut scene = probe_scene(2.0, true, Behaviour::Idle, None);
    let stats = solver(3, 8).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(stats.iterations_used, 3, "one sweep per substep suffices");
    assert!(stats.converged);
    assert_eq!(stats.residual, 0.0);
}

#[test]
fn constraint_that_keeps_changing_velocities_never_converges() {
    // The probe shifts the velocity by 1 m/s on every sweep and reports
    // nothing: the solver measures the change and cannot reach the tolerance.
    let mut scene = probe_scene(0.0, true, Behaviour::Drift(1.0), None);
    let stats = solver(2, 5).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(!stats.converged, "a drifting constraint is not converged");
    assert_eq!(stats.iterations_used, 10, "the whole budget is used");
    assert_eq!(stats.residual, 1.0);
    assert_eq!(velocity_of(&scene.bodies, scene.body).x, 10.0);
}

#[test]
fn reported_residual_decides_convergence() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, Some(0.25));
    let stats = solver(2, 3).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(!stats.converged);
    assert_eq!(stats.iterations_used, 6);
    assert_eq!(stats.residual, 0.25);

    let mut scene = probe_scene(1.0, true, Behaviour::Idle, Some(1e-9));
    let stats = solver(2, 3).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged, "a residual below the tolerance converges");
    assert_eq!(stats.iterations_used, 2);
    assert_eq!(stats.residual, 1e-9);
}

#[test]
fn tolerance_field_controls_when_a_substep_stops() {
    // v <- v / 2 on every sweep: the update of sweep k has size 2^-k.
    let make = || probe_scene(1.0, true, Behaviour::Decay(0.5), None);

    let mut scene = make();
    let stats = solver(1, 100).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged);
    assert_eq!(
        stats.iterations_used, 20,
        "2^-20 is the first update <= 1e-6"
    );
    assert!((stats.residual - 0.5_f64.powi(20)).abs() < 1e-15);

    let mut scene = make();
    let loose = solver(1, 100).with_tolerance(1e-3);
    let stats = loose.step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(stats.converged);
    assert_eq!(
        stats.iterations_used, 10,
        "2^-10 is the first update <= 1e-3"
    );

    let mut scene = make();
    let strict = solver(1, 20).with_tolerance(1e-9);
    let stats = strict.step(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(!stats.converged, "20 sweeps cannot reach 1e-9");
    assert_eq!(stats.iterations_used, 20);
}

#[test]
fn default_tolerance_is_one_micro_unit() {
    assert_eq!(TgsSolver::default().tolerance, 1e-6);
    assert_eq!(solver(2, 2).tolerance, 1e-6);
    assert_eq!(
        TgsSolver::with_options(2, 2, 1, Vec3::zeros(), 0.2, false).tolerance,
        1e-6
    );
}

#[test]
fn step_reports_the_worst_substep_not_the_last_one() {
    // The constraint reports a residual of 0.5 for its first two sweeps (the
    // whole first substep) and 0.0 afterwards; a stat computed from the last
    // substep alone would say "converged", the worst substep says otherwise.
    struct Fading {
        handle: BodyHandle,
        calls: usize,
    }
    impl Constraint for Fading {
        fn prepare(&mut self, _bodies: &RigidBodySet, _dt: f64) {}
        fn solve_velocity(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {
            self.calls += 1;
        }
        fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
        fn body_handles(&self) -> Vec<BodyHandle> {
            vec![self.handle]
        }
        fn residual(&self) -> f64 {
            if self.calls <= 2 { 0.5 } else { 0.0 }
        }
    }
    let mut bodies = RigidBodySet::new();
    let handle = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::zeros()));
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(Fading { handle, calls: 0 })];
    // Substep 1 uses both of its sweeps and ends at residual 0.5; substep 2
    // converges with its first sweep (residual 0.0).
    let stats = solver(2, 2).step(&mut constraints, &mut bodies, DT);
    assert_eq!(stats.iterations_used, 3);
    assert_eq!(stats.residual, 0.5);
    assert!(!stats.converged);
}

#[test]
fn nothing_to_iterate_reports_no_iteration_performed() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, None);
    let stats = solver(2, 0).step(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(stats.iterations_used, 0);
    assert!(!stats.converged);
    assert_eq!(stats.residual, f64::MAX);
}

#[test]
fn degenerate_steps_leave_the_bodies_untouched() {
    for (label, tgs, dt) in [
        (
            "zero substeps",
            TgsSolver::new(0, 4, Vec3::new(0.0, -9.81, 0.0)),
            DT,
        ),
        (
            "dt = 0",
            TgsSolver::new(2, 4, Vec3::new(0.0, -9.81, 0.0)),
            0.0,
        ),
        (
            "dt < 0",
            TgsSolver::new(2, 4, Vec3::new(0.0, -9.81, 0.0)),
            -DT,
        ),
        (
            "dt = NaN",
            TgsSolver::new(2, 4, Vec3::new(0.0, -9.81, 0.0)),
            f64::NAN,
        ),
        (
            "dt = inf",
            TgsSolver::new(2, 4, Vec3::new(0.0, -9.81, 0.0)),
            f64::INFINITY,
        ),
    ] {
        let mut scene = probe_scene(2.0, true, Behaviour::Drift(1.0), None);
        let stats = tgs.step(&mut scene.constraints, &mut scene.bodies, dt);
        assert_eq!(stats.iterations_used, 0, "{label}");
        assert!(!stats.converged, "{label}");
        assert_eq!(stats.residual, f64::MAX, "{label}");
        assert_eq!(
            velocity_of(&scene.bodies, scene.body),
            Vec3::new(2.0, 0.0, 0.0),
            "{label}: no gravity, no constraint impulse"
        );
        assert_eq!(
            position_of(&scene.bodies, scene.body),
            Vec3::zeros(),
            "{label}"
        );
        assert!(scene.log.prepare_dts.borrow().is_empty(), "{label}");
    }
}

// ── Hooks honoured by every TGS entry point ──────────────────────────────────

type RunFn = fn(&TgsSolver, &mut ProbeScene) -> SolverStats;

fn run_step(tgs: &TgsSolver, scene: &mut ProbeScene) -> SolverStats {
    tgs.step(&mut scene.constraints, &mut scene.bodies, DT)
}

fn run_solve(tgs: &TgsSolver, scene: &mut ProbeScene) -> SolverStats {
    tgs.solve(&mut scene.constraints, &mut scene.bodies, DT)
}

fn run_diagnostics(tgs: &TgsSolver, scene: &mut ProbeScene) -> SolverStats {
    let (stats, diagnostics) =
        tgs.step_with_diagnostics(&mut scene.constraints, &mut scene.bodies, DT);
    assert_eq!(diagnostics.len(), tgs.substeps);
    stats
}

fn run_island(tgs: &TgsSolver, scene: &mut ProbeScene) -> SolverStats {
    tgs.solve_island(&[scene.body], &mut scene.constraints, &mut scene.bodies, DT)
}

const ALL_ENTRY_POINTS: [(&str, RunFn); 4] = [
    ("step", run_step),
    ("solve", run_solve),
    ("step_with_diagnostics", run_diagnostics),
    ("solve_island", run_island),
];

#[test]
fn inactive_constraints_are_skipped_by_every_phase_of_every_entry_point() {
    for (name, run) in ALL_ENTRY_POINTS {
        let mut scene = probe_scene(2.0, false, Behaviour::Drift(1.0), None);
        let stats = run(
            &TgsSolver::with_options(3, 4, 2, Vec3::zeros(), 0.2, false),
            &mut scene,
        );
        let log = &scene.log;
        assert!(log.prepare_dts.borrow().is_empty(), "{name}: prepare");
        assert!(
            log.velocity_dts.borrow().is_empty(),
            "{name}: solve_velocity"
        );
        assert!(
            log.position_dts.borrow().is_empty(),
            "{name}: solve_position"
        );
        assert_eq!(log.post_solve.get(), 0, "{name}: post_solve");
        assert_eq!(
            velocity_of(&scene.bodies, scene.body),
            Vec3::new(2.0, 0.0, 0.0),
            "{name}: an inactive constraint must apply no impulse"
        );
        assert!(stats.converged, "{name}: nothing active means converged");
        assert_eq!(stats.residual, 0.0, "{name}");
        assert_eq!(
            stats.iterations_used, 3,
            "{name}: one empty sweep per substep"
        );
    }
}

#[test]
fn post_solve_runs_once_per_substep_after_the_sweeps() {
    for (name, run) in ALL_ENTRY_POINTS {
        let mut scene = probe_scene(0.0, true, Behaviour::Drift(1.0), None);
        let tgs = TgsSolver::with_options(3, 4, 1, Vec3::zeros(), 0.2, false);
        let stats = run(&tgs, &mut scene);
        assert_eq!(scene.log.post_solve.get(), 3, "{name}: once per substep");
        assert_eq!(scene.log.prepare_dts.borrow().len(), 3, "{name}: prepare");
        assert_eq!(scene.log.velocity_dts.borrow().len(), 12, "{name}: sweeps");
        assert_eq!(
            scene.log.position_dts.borrow().len(),
            3,
            "{name}: positions"
        );
        assert_eq!(stats.iterations_used, 12, "{name}");
        assert!(!stats.converged, "{name}");
    }
}

// ── Contacts under TGS ───────────────────────────────────────────────────────

/// A unit box resting on static ground under gravity (contact at its centre
/// line, no penetration).
fn resting_box(penetration: f64) -> (RigidBodySet, BodyHandle, Vec<Box<dyn Constraint>>) {
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let body = bodies.insert(dynamic_body(1.0, Vec3::new(0.0, 0.5, 0.0), Vec3::zeros()));
    let contact = ContactConstraint::new(
        body,
        ground,
        Vec3::new(0.0, 1.0, 0.0),
        penetration,
        Vec3::zeros(),
        Vec3::zeros(),
        0.0,
        0.5,
    );
    (bodies, body, vec![Box::new(contact)])
}

#[test]
fn resting_contact_under_gravity_converges() {
    let (mut bodies, body, mut constraints) = resting_box(0.0);
    let tgs = TgsSolver::new(4, 8, Vec3::new(0.0, -9.81, 0.0));
    let stats = tgs.step(&mut constraints, &mut bodies, DT);
    assert!(stats.converged, "residual = {}", stats.residual);
    assert!(stats.residual <= tgs.tolerance);
    assert!((4..=32).contains(&stats.iterations_used));
    let vy = velocity_of(&bodies, body).y;
    assert!(
        vy.abs() < 1e-9,
        "the contact must cancel gravity, vy = {vy}"
    );
}

#[test]
fn diagnostics_report_the_contact_penetration() {
    let (mut bodies, _body, mut constraints) = resting_box(0.03);
    let tgs = TgsSolver::new(3, 6, Vec3::zeros());
    let (stats, diagnostics) = tgs.step_with_diagnostics(&mut constraints, &mut bodies, DT);
    assert_eq!(diagnostics.len(), 3);
    assert!(stats.converged);
    for (substep, d) in diagnostics.iter().enumerate() {
        assert_eq!(
            d.max_penetration, 0.03,
            "substep {substep}: the largest contact penetration is reported"
        );
    }
}

#[test]
fn contact_reports_its_penetration_estimate() {
    let (_, _, overlapping) = resting_box(0.03);
    assert_eq!(overlapping[0].penetration_estimate(), 0.03);
    let (_, _, separated) = resting_box(-0.02);
    assert_eq!(
        separated[0].penetration_estimate(),
        0.0,
        "a separated contact does not penetrate"
    );
}

#[test]
fn diagnostics_without_penetrating_constraints_report_zero() {
    let mut scene = probe_scene(1.0, true, Behaviour::Idle, None);
    let (_, diagnostics) =
        solver(2, 2).step_with_diagnostics(&mut scene.constraints, &mut scene.bodies, DT);
    assert!(diagnostics.iter().all(|d| d.max_penetration == 0.0));
}

#[test]
fn inactive_contacts_do_not_count_towards_max_penetration() {
    struct Hidden(ContactConstraint);
    impl Constraint for Hidden {
        fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
            self.0.prepare(bodies, dt);
        }
        fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
            self.0.solve_velocity(bodies, dt);
        }
        fn solve_position(&mut self, bodies: &mut RigidBodySet, dt: f64) {
            self.0.solve_position(bodies, dt);
        }
        fn body_handles(&self) -> Vec<BodyHandle> {
            self.0.body_handles()
        }
        fn is_active(&self) -> bool {
            false
        }
        fn penetration_estimate(&self) -> f64 {
            self.0.penetration_estimate()
        }
    }
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let body = bodies.insert(dynamic_body(1.0, Vec3::new(0.0, 0.5, 0.0), Vec3::zeros()));
    let contact = ContactConstraint::new(
        body,
        ground,
        Vec3::new(0.0, 1.0, 0.0),
        0.2,
        Vec3::zeros(),
        Vec3::zeros(),
        0.0,
        0.0,
    );
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(Hidden(contact))];
    let (_, diagnostics) = solver(1, 2).step_with_diagnostics(&mut constraints, &mut bodies, DT);
    assert_eq!(diagnostics[0].max_penetration, 0.0);
}

// ── Island solve ─────────────────────────────────────────────────────────────

#[test]
fn island_solve_integrates_only_the_island_bodies() {
    let mut bodies = RigidBodySet::new();
    let inside = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0)));
    let outside = bodies.insert(dynamic_body(1.0, Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0)));
    let mut constraints: Vec<Box<dyn Constraint>> = Vec::new();
    let tgs = solver(4, 2);
    let stats = tgs.solve_island(&[inside], &mut constraints, &mut bodies, DT);
    assert!(stats.converged, "no constraints, nothing to resolve");
    assert!((position_of(&bodies, inside).x - DT).abs() < 1e-12);
    assert_eq!(position_of(&bodies, outside), Vec3::zeros());
}
