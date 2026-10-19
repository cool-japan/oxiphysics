use super::*;
use oxiphysics_rigid::RigidBody;

fn make_two_bodies() -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();

    let mut b_a = RigidBody::new(1.0);
    b_a.transform.position = Vec3::new(0.0, 0.0, 0.0);
    b_a.linear_damping = 0.0;
    b_a.angular_damping = 0.0;
    let ha = bodies.insert(b_a);

    let mut b_b = RigidBody::new(1.0);
    b_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
    b_b.linear_damping = 0.0;
    b_b.angular_damping = 0.0;
    let hb = bodies.insert(b_b);

    (bodies, ha, hb)
}

#[test]
fn test_six_dof_all_free_by_default() {
    let (_, ha, hb) = make_two_bodies();
    let c = SixDofConstraint::new(ha, hb);
    for i in 0..3 {
        assert!(!c.linear[i].locked);
        assert!(!c.angular[i].locked);
        assert!(c.linear[i].limited.is_none());
        assert!(c.angular[i].limited.is_none());
    }
}

#[test]
fn test_lock_linear() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(1);
    assert!(c.linear[1].locked);
    assert!(!c.linear[0].locked);
    assert!(!c.linear[2].locked);
}

#[test]
fn test_lock_angular() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_angular(2);
    assert!(c.angular[2].locked);
    assert!(!c.angular[0].locked);
    assert!(!c.angular[1].locked);
}

#[test]
fn test_set_linear_limit() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -1.0, 1.0).unwrap();
    assert_eq!(c.linear[0].limited, Some((-1.0, 1.0)));
}

#[test]
fn test_set_angular_limit() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_angular_limit(0, -0.5, 0.5).unwrap();
    assert_eq!(c.angular[0].limited, Some((-0.5, 0.5)));
}

#[test]
fn test_weld_locks_all_dof() {
    let (_, ha, hb) = make_two_bodies();
    let c = SixDofConstraint::weld(ha, hb);
    for i in 0..3 {
        assert!(c.linear[i].locked, "linear[{i}] should be locked");
        assert!(c.angular[i].locked, "angular[{i}] should be locked");
    }
}

#[test]
fn test_locked_linear_dof_produces_zero_velocity() {
    let mut bodies = RigidBodySet::new();

    // Place both bodies at the same position so there is no positional error;
    // the test checks only that the velocity constraint zeroes relative motion.
    let mut b_a = RigidBody::new(1.0);
    b_a.transform.position = Vec3::new(0.0, 0.0, 0.0);
    b_a.velocity = Vec3::new(5.0, 0.0, 0.0);
    b_a.linear_damping = 0.0;
    b_a.angular_damping = 0.0;
    let ha = bodies.insert(b_a);

    let mut b_b = RigidBody::new(1.0);
    b_b.transform.position = Vec3::new(0.0, 0.0, 0.0); // same position → no bias
    b_b.velocity = Vec3::new(-5.0, 0.0, 0.0);
    b_b.linear_damping = 0.0;
    b_b.angular_damping = 0.0;
    let hb = bodies.insert(b_b);

    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(0); // lock X translation

    let dt = 1.0 / 60.0;
    c.prepare(&bodies, dt);

    for _ in 0..20 {
        c.solve_velocity(&mut bodies, dt);
    }

    let vel_a = bodies.get(ha).unwrap().velocity;
    let vel_b = bodies.get(hb).unwrap().velocity;
    let rel_vel_x = (vel_a - vel_b).x;

    assert!(
        rel_vel_x.abs() < 0.5,
        "Locked X DOF should produce near-zero relative velocity along X, got {rel_vel_x}"
    );
}

#[test]
fn test_locked_angular_dof_produces_zero_ang_vel() {
    let (mut bodies, ha, hb) = make_two_bodies();

    // Give bodies opposing angular velocities about Z
    bodies.get_mut(ha).unwrap().angular_velocity = Vec3::new(0.0, 0.0, 5.0);
    bodies.get_mut(hb).unwrap().angular_velocity = Vec3::new(0.0, 0.0, -5.0);

    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_angular(2); // lock Z rotation

    let dt = 1.0 / 60.0;
    c.prepare(&bodies, dt);

    for _ in 0..20 {
        c.solve_velocity(&mut bodies, dt);
    }

    let ang_a = bodies.get(ha).unwrap().angular_velocity;
    let ang_b = bodies.get(hb).unwrap().angular_velocity;
    let rel_ang_z = (ang_a - ang_b).z;

    assert!(
        rel_ang_z.abs() < 1.0,
        "Locked Z angular DOF should produce near-zero relative angular velocity about Z, got {rel_ang_z}"
    );
}

#[test]
fn test_body_handles() {
    let (_, ha, hb) = make_two_bodies();
    let c = SixDofConstraint::new(ha, hb);
    let handles = c.body_handles();
    assert_eq!(handles.len(), 2);
    assert_eq!(handles[0], ha);
    assert_eq!(handles[1], hb);
}

#[test]
fn test_jacobian_linear_rows_orthogonality() {
    let r_a = Vec3::new(0.1, 0.2, 0.3);
    let r_b = Vec3::new(-0.1, 0.1, -0.2);
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];

    let j = compute_jacobian_6dof(r_a, r_b, &axes);

    // Rows 0, 1, 2 should correspond to X, Y, Z linear axes
    // The first 3 columns of each linear row should match the axis
    for (row, ax) in axes.iter().enumerate() {
        assert!((j[row][0] - ax.x).abs() < 1e-12);
        assert!((j[row][1] - ax.y).abs() < 1e-12);
        assert!((j[row][2] - ax.z).abs() < 1e-12);
    }

    // Angular rows (3,4,5) body-A cols should equal axes
    for (i, ax) in axes.iter().enumerate() {
        let row = 3 + i;
        assert!((j[row][3] - ax.x).abs() < 1e-12);
        assert!((j[row][4] - ax.y).abs() < 1e-12);
        assert!((j[row][5] - ax.z).abs() < 1e-12);
    }
}

#[test]
fn test_six_dof_with_motor() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    let motor = MotorConfig::new(2.0, 100.0, 50.0);
    c.set_linear_motor(0, motor).unwrap();
    assert!(c.linear[0].motor.is_some());
    assert_eq!(c.linear[0].motor.as_ref().unwrap().target_velocity, 2.0);
}

#[test]
#[should_panic]
fn test_lock_linear_out_of_bounds() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(3); // should panic
}

#[test]
#[should_panic]
fn test_lock_angular_out_of_bounds() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_angular(3); // should panic
}

// ── Spring/damper ─────────────────────────────────────────────────────────

#[test]
fn test_spring_damper_config_creation() {
    let sd = SpringDamperConfig::new(100.0, 10.0, 1.0, 0.0);
    assert!((sd.stiffness - 100.0).abs() < 1e-12);
    assert!((sd.damping - 10.0).abs() < 1e-12);
}

#[test]
fn test_spring_force_at_rest_position() {
    let sd = SpringDamperConfig::new(100.0, 10.0, 1.0, 0.0);
    let f = sd.force(1.0, 0.0); // at rest length, zero velocity → zero force
    assert!(f.abs() < 1e-12, "At rest: force = {f}");
}

#[test]
fn test_spring_force_compressed() {
    let sd = SpringDamperConfig::new(100.0, 0.0, 0.5, 0.0);
    let f = sd.force(0.0, 0.0); // displacement = rest_length from neutral
    // Extension = 0.0 - 0.5 = -0.5; F = -k * (-0.5) = +50
    assert!((f - 50.0).abs() < 1e-10, "Compressed spring: f = {f}");
}

#[test]
fn test_damper_force_velocity() {
    let sd = SpringDamperConfig::new(0.0, 50.0, 0.0, 0.0);
    let f = sd.force(0.0, 2.0); // velocity = 2 m/s
    // F = -c * vel = -100
    assert!((f + 100.0).abs() < 1e-10, "Damper: f = {f}");
}

// ── Velocity drive ────────────────────────────────────────────────────────

#[test]
fn test_velocity_drive_config_linear() {
    let vd = VelocityDriveConfig::new_linear(5.0, 1000.0);
    assert!((vd.target_velocity - 5.0).abs() < 1e-12);
    assert!((vd.max_force - 1000.0).abs() < 1e-12);
}

#[test]
fn test_velocity_drive_impulse_positive_error() {
    let vd = VelocityDriveConfig::new_linear(5.0, 1000.0);
    let imp = vd.compute_impulse(3.0, 2.0, 1.0 / 60.0); // actual_vel < target
    assert!(
        imp > 0.0,
        "Positive velocity error → positive impulse: {imp}"
    );
}

#[test]
fn test_velocity_drive_impulse_clamped() {
    let vd = VelocityDriveConfig::new_linear(100.0, 1.0); // small max force
    let imp = vd.compute_impulse(0.0, 10.0, 1.0 / 60.0);
    let limit = 1.0 * (1.0 / 60.0);
    assert!(
        imp.abs() <= limit + 1e-12,
        "Impulse should be clamped: {imp}"
    );
}

// ── Position drive ────────────────────────────────────────────────────────

#[test]
fn test_position_drive_config_creation() {
    let pd = PositionDriveConfig::new(1.0, 100.0, 10.0, 500.0);
    assert!((pd.target_position - 1.0).abs() < 1e-12);
}

#[test]
fn test_position_drive_force_toward_target() {
    let pd = PositionDriveConfig::new(2.0, 100.0, 10.0, 500.0);
    let f = pd.compute_force(0.0, 0.0); // position error = 2.0
    assert!(f > 0.0, "Force should drive toward target: {f}");
}

#[test]
fn test_position_drive_zero_error() {
    let pd = PositionDriveConfig::new(1.0, 100.0, 10.0, 500.0);
    let f = pd.compute_force(1.0, 0.0); // at target, zero velocity
    assert!(f.abs() < 1e-10, "Zero position error: f = {f}");
}

// ── Compliance / restitution ──────────────────────────────────────────────

#[test]
fn test_compliance_reduces_effective_mass() {
    let em_no_compliance = 1.0; // 1/k
    let em_with_compliance = apply_compliance(em_no_compliance, 0.01, 1.0 / 60.0);
    assert!(
        em_with_compliance < em_no_compliance,
        "Compliance should reduce effective mass: {em_with_compliance} vs {em_no_compliance}"
    );
}

#[test]
fn test_restitution_impulse_sign() {
    let imp = restitution_impulse(1.0, -2.0, 0.5); // approaching velocity -2
    assert!(
        imp > 0.0,
        "Restitution impulse should oppose approach: {imp}"
    );
}

#[test]
fn test_restitution_zero_coefficient() {
    let imp = restitution_impulse(1.0, -2.0, 0.0);
    // e=0: perfectly inelastic; impulse = eff_mass * |v_approach| = 2
    assert!((imp - 2.0).abs() < 1e-12, "Zero restitution: imp = {imp}");
}

// ── 6-DOF limits ─────────────────────────────────────────────────────────

#[test]
fn test_axis_config_free_no_limits() {
    let ac = AxisConfig::free();
    assert!(!ac.locked);
    assert!(ac.limited.is_none());
}

#[test]
fn test_axis_config_limited_clamps_correctly() {
    let ac = AxisConfig::limited(-1.0, 1.0);
    assert_eq!(ac.limited, Some((-1.0, 1.0)));
}

// ── Jacobian correctness ──────────────────────────────────────────────────

#[test]
fn test_jacobian_body_b_linear_is_negated() {
    let r_a = Vec3::new(0.0, 0.0, 0.0);
    let r_b = Vec3::new(0.0, 0.0, 0.0);
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let j = compute_jacobian_6dof(r_a, r_b, &axes);
    // For row 0 (X axis): columns 0..3 = [1,0,0], columns 6..9 = [-1,0,0]
    assert!(
        (j[0][6] + 1.0).abs() < 1e-12,
        "Body B X-linear should be -1: {}",
        j[0][6]
    );
    assert!((j[0][7]).abs() < 1e-12, "Body B Y should be 0: {}", j[0][7]);
}

#[test]
fn test_jacobian_angular_rows_body_b_negated() {
    let r_a = Vec3::zeros();
    let r_b = Vec3::zeros();
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let j = compute_jacobian_6dof(r_a, r_b, &axes);
    // Row 3 (X angular): body B angular should be -1
    assert!(
        (j[3][9] + 1.0).abs() < 1e-12,
        "Body B angular X row: {}",
        j[3][9]
    );
}

#[test]
fn test_jacobian_r_cross_axis_terms() {
    // Non-zero r_a: r_a = [0, 0, 1], axis = [1, 0, 0]
    // r_a × ax = [0,0,1] × [1,0,0] = [0*0-1*0, 1*1-0*0, 0*0-0*1] = [0,1,0]
    let r_a = Vec3::new(0.0, 0.0, 1.0);
    let r_b = Vec3::zeros();
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let j = compute_jacobian_6dof(r_a, r_b, &axes);
    // Row 0 (X axis), cols 3..6 = r_a × ax = [0,1,0]
    assert!((j[0][3]).abs() < 1e-12, "j[0][3] = {}", j[0][3]);
    assert!((j[0][4] - 1.0).abs() < 1e-12, "j[0][4] = {}", j[0][4]);
    assert!((j[0][5]).abs() < 1e-12, "j[0][5] = {}", j[0][5]);
}

// ── Motor angular ─────────────────────────────────────────────────────────

#[test]
fn test_angular_motor_config_attached() {
    let (_, ha, hb) = make_two_bodies();

    let mut c = SixDofConstraint::new(ha, hb);
    let motor = MotorConfig::new(5.0, 1000.0, 100.0);
    c.set_angular_motor(2, motor).unwrap();

    // Verify motor config was stored
    assert!(
        c.angular[2].motor.is_some(),
        "Motor should be attached to angular axis 2"
    );
    let m = c.angular[2].motor.as_ref().unwrap();
    assert!((m.target_velocity - 5.0).abs() < 1e-12);
    assert!((m.max_force - 1000.0).abs() < 1e-12);
}

// ── Test world helpers ───────────────────────────────────────────────────────

const DT: f64 = 1.0 / 60.0;

fn dynamic_body(mass: f64, position: Vec3) -> RigidBody {
    let mut body = RigidBody::new(mass);
    body.transform.position = position;
    body.linear_damping = 0.0;
    body.angular_damping = 0.0;
    body
}

/// A static body A at the origin and a dynamic unit-mass body B at the origin.
fn static_and_dynamic() -> (RigidBodySet, BodyHandle, BodyHandle) {
    let mut bodies = RigidBodySet::new();
    let ha = bodies.insert(RigidBody::new_static());
    let hb = bodies.insert(dynamic_body(1.0, Vec3::zeros()));
    (bodies, ha, hb)
}

/// Two dynamic unit-mass bodies at the origin.
fn coincident_bodies() -> (RigidBodySet, BodyHandle, BodyHandle) {
    let (mut bodies, ha, hb) = make_two_bodies();
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.position = Vec3::zeros();
    }
    (bodies, ha, hb)
}

/// Rotate a body about a world axis and refresh its world inertia.
fn rotate_body(bodies: &mut RigidBodySet, handle: BodyHandle, axis: Vec3, angle: f64) {
    if let Some(body) = bodies.get_mut(handle) {
        let turn = Quat::from_scaled_axis(axis.normalize() * angle);
        body.transform.rotation = turn * body.transform.rotation;
        body.update_world_inertia();
    }
}

fn set_velocity(bodies: &mut RigidBodySet, handle: BodyHandle, linear: Vec3, angular: Vec3) {
    if let Some(body) = bodies.get_mut(handle) {
        body.velocity = linear;
        body.angular_velocity = angular;
    }
}

fn position_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies
        .get(handle)
        .map(|b| b.transform.position)
        .unwrap_or_default()
}

fn velocity_of(bodies: &RigidBodySet, handle: BodyHandle) -> Vec3 {
    bodies.get(handle).map(|b| b.velocity).unwrap_or_default()
}

fn integrate_all(bodies: &mut RigidBodySet, dt: f64) {
    let handles: Vec<BodyHandle> = bodies.iter().map(|(handle, _)| handle).collect();
    for handle in handles {
        if let Some(body) = bodies.get_mut(handle) {
            body.integrate_velocity(dt);
        }
    }
}

/// One simulation step of a single constraint: the phases of a PGS solve
/// followed by the position integration.
fn step(
    bodies: &mut RigidBodySet,
    c: &mut SixDofConstraint,
    velocity_iterations: usize,
    position_iterations: usize,
) {
    c.prepare(bodies, DT);
    for _ in 0..velocity_iterations {
        c.solve_velocity(bodies, DT);
    }
    for _ in 0..position_iterations {
        c.solve_position(bodies, DT);
    }
    integrate_all(bodies, DT);
}

/// The joint coordinate `dof` of `c` for the current body state.
fn coordinate_of(c: &mut SixDofConstraint, bodies: &RigidBodySet, dof: usize) -> f64 {
    c.prepare(bodies, DT);
    c.coordinates()[dof]
}

/// Angle of the rotation of B relative to A.
fn relative_orientation_angle(bodies: &RigidBodySet, ha: BodyHandle, hb: BodyHandle) -> f64 {
    let qa = bodies.get(ha).map(|b| b.transform.rotation);
    let qb = bodies.get(hb).map(|b| b.transform.rotation);
    match (qa, qb) {
        (Some(qa), Some(qb)) => (qa.inverse() * qb).angle(),
        _ => f64::NAN,
    }
}

/// World-space distance between the two anchors of `c`.
fn anchor_gap(bodies: &RigidBodySet, c: &SixDofConstraint) -> f64 {
    match (bodies.get(c.body_a), bodies.get(c.body_b)) {
        (Some(a), Some(b)) => {
            let anchor_a = a.transform.position + a.transform.rotation * c.anchor_a;
            let anchor_b = b.transform.position + b.transform.rotation * c.anchor_b;
            (anchor_b - anchor_a).norm()
        }
        _ => f64::NAN,
    }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: expected {expected}, got {actual}"
    );
}

// ── F0481: position correction must close the anchor gap ─────────────────────

#[test]
fn test_solve_position_closes_gap_monotonically() {
    let (mut bodies, ha, hb) = coincident_bodies();
    // Anchor B sits 0.1 m further along +X than anchor A.
    let mut c =
        SixDofConstraint::weld(ha, hb).with_anchors(Vec3::zeros(), Vec3::new(0.1, 0.0, 0.0));

    let mut previous = anchor_gap(&bodies, &c);
    assert_close(previous, 0.1, 1e-12, "initial gap");
    for step in 0..60 {
        c.solve_position(&mut bodies, DT);
        let gap = anchor_gap(&bodies, &c);
        assert!(
            gap < previous,
            "step {step}: gap grew from {previous} to {gap}"
        );
        previous = gap;
    }
    // Each pass removes a fraction BAUMGARTE_6DOF of the remaining gap.
    assert_close(
        previous,
        0.1 * (1.0 - BAUMGARTE_6DOF).powi(60),
        1e-12,
        "gap after 60 passes",
    );
    assert!(previous < 1e-6);

    // Equal masses share the correction: A moved towards B's anchor (+X),
    // B moved towards A's anchor (-X), and the centre of mass stayed put.
    let pos_a = position_of(&bodies, ha);
    let pos_b = position_of(&bodies, hb);
    assert!(pos_a.x > 0.0, "body A must move along +X, got {}", pos_a.x);
    assert!(pos_b.x < 0.0, "body B must move along -X, got {}", pos_b.x);
    assert!((pos_a + pos_b).norm() < 1e-12, "centre of mass moved");
}

#[test]
fn test_solve_position_moves_only_the_dynamic_body() {
    let mut bodies = RigidBodySet::new();
    let ha = bodies.insert(RigidBody::new_static());
    let hb = bodies.insert(dynamic_body(2.0, Vec3::new(0.3, 0.0, 0.0)));
    let mut c = SixDofConstraint::weld(ha, hb);

    let before = anchor_gap(&bodies, &c);
    assert_close(before, 0.3, 1e-12, "initial gap");
    c.solve_position(&mut bodies, DT);
    assert_eq!(
        position_of(&bodies, ha),
        Vec3::zeros(),
        "the static body must not move"
    );
    assert_close(
        anchor_gap(&bodies, &c),
        before * (1.0 - BAUMGARTE_6DOF),
        1e-12,
        "one Baumgarte step",
    );
}

#[test]
fn test_solve_position_with_lever_arms_converges() {
    // The anchors are far from the centres of mass, so closing the gap also
    // turns the bodies; the angular DOFs are free.
    let (mut bodies, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb)
        .with_anchors(Vec3::new(0.5, 0.3, 0.1), Vec3::new(-0.2, 0.1, 0.4));
    for axis in 0..3 {
        c.lock_linear(axis);
    }
    let initial = anchor_gap(&bodies, &c);
    assert!(initial > 0.4);
    for _ in 0..300 {
        c.solve_position(&mut bodies, DT);
    }
    assert!(
        anchor_gap(&bodies, &c) < 1e-6,
        "gap = {}",
        anchor_gap(&bodies, &c)
    );
    assert!(
        relative_orientation_angle(&bodies, ha, hb) > 1e-3,
        "closing an off-centre gap must rotate the bodies"
    );
}

// ── F0478: angular locks measure the orientation error ───────────────────────

#[test]
fn test_weld_removes_an_injected_orientation_error() {
    // The reference rotation is captured by the first prepare; the error is
    // injected afterwards, so only a real orientation measurement can see it.
    for position_iterations in [0, 4] {
        let (mut bodies, ha, hb) = coincident_bodies();
        let mut c = SixDofConstraint::weld(ha, hb);
        step(&mut bodies, &mut c, 8, position_iterations);
        assert!(c.reference_rotation().is_some());

        rotate_body(&mut bodies, hb, Vec3::new(0.0, 0.0, 1.0), 0.2);
        let mut previous = relative_orientation_angle(&bodies, ha, hb);
        assert_close(previous, 0.2, 1e-9, "injected error");
        for k in 0..60 {
            step(&mut bodies, &mut c, 8, position_iterations);
            let angle = relative_orientation_angle(&bodies, ha, hb);
            // The error shrinks geometrically until it reaches the numerical
            // floor of the quaternion arithmetic.
            if previous > 1e-9 {
                assert!(
                    angle < previous,
                    "{position_iterations} position passes, step {k}: error grew {previous} -> {angle}"
                );
            }
            previous = angle;
        }
        assert!(
            previous < 1e-3,
            "{position_iterations} position passes: orientation error {previous} rad"
        );
    }
}

#[test]
fn test_weld_keeps_a_spinning_body_aligned() {
    // Body B starts spinning about Z; the weld has to make A follow.
    let (mut bodies, ha, hb) = coincident_bodies();
    set_velocity(&mut bodies, hb, Vec3::zeros(), Vec3::new(0.0, 0.0, 1.0));
    let mut c = SixDofConstraint::weld(ha, hb);
    for _ in 0..120 {
        step(&mut bodies, &mut c, 8, 2);
    }
    assert!(
        relative_orientation_angle(&bodies, ha, hb) < 1e-3,
        "relative angle = {}",
        relative_orientation_angle(&bodies, ha, hb)
    );
    // Equal inertias share the spin: both bodies turn at 0.5 rad/s.
    let spin_a = bodies
        .get(ha)
        .map(|b| b.angular_velocity)
        .unwrap_or_default();
    assert_close(spin_a.z, 0.5, 1e-9, "shared spin");
}

#[test]
fn test_locked_angular_axis_corrects_only_its_own_axis() {
    // Only rotation about Z is locked.
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(&mut bodies, hb, Vec3::zeros(), Vec3::new(0.7, 0.0, 0.4));
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_angular(2);
    c.prepare(&bodies, DT);
    c.solve_velocity(&mut bodies, DT);
    let spin = bodies
        .get(hb)
        .map(|b| b.angular_velocity)
        .unwrap_or_default();
    assert_close(
        spin.x,
        0.7,
        1e-12,
        "the free X rotation must not be touched",
    );
    assert_close(spin.z, 0.0, 1e-12, "the locked Z rotation must be removed");

    // A pure Z error is removed by the position pass without disturbing X.
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_angular(2);
    c.prepare(&bodies, DT);
    rotate_body(&mut bodies, hb, Vec3::new(0.0, 0.0, 1.0), 0.1);
    assert_close(
        coordinate_of(&mut c, &bodies, 5),
        0.1,
        1e-12,
        "injected error",
    );
    for _ in 0..80 {
        step(&mut bodies, &mut c, 8, 2);
    }
    assert!(
        coordinate_of(&mut c, &bodies, 5).abs() < 1e-3,
        "the locked Z rotation must be removed"
    );
    assert_close(
        coordinate_of(&mut c, &bodies, 3),
        0.0,
        1e-12,
        "X stays untouched",
    );
}

#[test]
fn test_angular_coordinates_are_measured_in_the_frame_of_body_a() {
    // A non-trivial reference rotation: A is turned 0.4 rad about X and B sits
    // 0.7 rad about Y relative to A.
    let (mut bodies, ha, hb) = coincident_bodies();
    rotate_body(&mut bodies, ha, Vec3::new(1.0, 0.0, 0.0), 0.4);
    let a_rotation = bodies.get(ha).map(|b| b.transform.rotation);
    let a_rotation = a_rotation.unwrap_or_else(Quat::identity);
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.rotation = a_rotation * Quat::from_scaled_axis(Vec3::new(0.0, 0.7, 0.0));
        b.update_world_inertia();
    }
    let mut c = SixDofConstraint::weld(ha, hb);
    c.prepare(&bodies, DT);
    for (dof, value) in c.coordinates().iter().enumerate() {
        assert!(
            value.abs() < 1e-12,
            "coordinate {dof} starts at zero: {value}"
        );
    }

    // Turn B by 0.2 rad about the world axis that is A's local Z.
    rotate_body(&mut bodies, hb, a_rotation * Vec3::new(0.0, 0.0, 1.0), 0.2);
    c.prepare(&bodies, DT);
    let coordinates = c.coordinates();
    assert_close(coordinates[3], 0.0, 1e-9, "angle about local X");
    assert_close(coordinates[4], 0.0, 1e-9, "angle about local Y");
    assert_close(coordinates[5], 0.2, 1e-9, "angle about local Z");

    // And the weld removes it again.
    for _ in 0..100 {
        step(&mut bodies, &mut c, 8, 2);
    }
    assert!(coordinate_of(&mut c, &bodies, 5).abs() < 1e-3);
}

#[test]
fn test_reference_rotation_is_captured_once_and_can_be_chosen() {
    let (mut bodies, ha, hb) = coincident_bodies();
    rotate_body(&mut bodies, hb, Vec3::new(0.0, 0.0, 1.0), 0.3);

    // Captured from the first prepare: the initial 0.3 rad is "zero".
    let mut captured = SixDofConstraint::weld(ha, hb);
    assert!(captured.reference_rotation().is_none());
    captured.prepare(&bodies, DT);
    assert!(captured.reference_rotation().is_some());
    assert_close(captured.coordinates()[5], 0.0, 1e-12, "captured reference");
    rotate_body(&mut bodies, hb, Vec3::new(0.0, 0.0, 1.0), 0.1);
    captured.prepare(&bodies, DT);
    assert_close(
        captured.coordinates()[5],
        0.1,
        1e-12,
        "reference stays fixed",
    );

    // Chosen explicitly: angles are measured from the aligned configuration.
    let mut aligned = SixDofConstraint::weld(ha, hb).with_reference_rotation(Quat::identity());
    aligned.prepare(&bodies, DT);
    assert_close(aligned.coordinates()[5], 0.4, 1e-12, "identity reference");
}

#[test]
fn test_rotation_vector_is_the_shortest_scaled_axis() {
    assert!(rotation_vector(&Quat::identity()).norm() < 1e-15);
    let scaled = Vec3::new(0.3, -0.2, 0.5);
    assert!((rotation_vector(&Quat::from_scaled_axis(scaled)) - scaled).norm() < 1e-12);
    let axis = Vec3::new(1.0, 2.0, 2.0) / 3.0;
    assert!((rotation_vector(&Quat::from_scaled_axis(axis * 3.0)) - axis * 3.0).norm() < 1e-12);
    // A turn of 4 rad about Z (w < 0) is the turn of 4 - 2π about Z.
    let long_way = rotation_vector(&Quat::from_scaled_axis(Vec3::new(0.0, 0.0, 4.0)));
    assert_close(
        long_way.z,
        4.0 - std::f64::consts::TAU,
        1e-12,
        "shortest way round",
    );
    // Tiny angles keep their direction.
    let tiny = Vec3::new(1e-14, 0.0, 0.0);
    assert!((rotation_vector(&Quat::from_scaled_axis(tiny)) - tiny).norm() < 1e-20);
}

// ── F0479: limits ────────────────────────────────────────────────────────────

#[test]
fn test_linear_limit_stops_a_pushed_body_at_the_upper_bound() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(&mut bodies, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::zeros());
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();

    let mut max_x = f64::MIN;
    for _ in 0..90 {
        step(&mut bodies, &mut c, 8, 2);
        max_x = max_x.max(position_of(&bodies, hb).x);
    }
    assert!(
        max_x <= 0.1 + 1e-9,
        "the body passed the limit: x = {max_x}"
    );
    assert_close(position_of(&bodies, hb).x, 0.1, 1e-6, "resting position");
    assert!(velocity_of(&bodies, hb).norm() < 1e-6, "the body must stop");
}

#[test]
fn test_linear_limit_stops_a_pushed_body_at_the_lower_bound() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(&mut bodies, hb, Vec3::new(-1.0, 0.0, 0.0), Vec3::zeros());
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();

    let mut min_x = f64::MAX;
    for _ in 0..90 {
        step(&mut bodies, &mut c, 8, 2);
        min_x = min_x.min(position_of(&bodies, hb).x);
    }
    assert!(
        min_x >= -0.1 - 1e-9,
        "the body passed the limit: x = {min_x}"
    );
    assert_close(position_of(&bodies, hb).x, -0.1, 1e-6, "resting position");
    assert!(velocity_of(&bodies, hb).norm() < 1e-6, "the body must stop");
}

#[test]
fn test_limit_of_body_a_frame_follows_the_axes_of_a() {
    // A is turned 90° about Z, so its local X axis is the world Y axis.
    let (mut bodies, ha, hb) = static_and_dynamic();
    rotate_body(
        &mut bodies,
        ha,
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::FRAC_PI_2,
    );
    set_velocity(&mut bodies, hb, Vec3::new(0.0, 1.0, 0.0), Vec3::zeros());
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();
    for _ in 0..60 {
        step(&mut bodies, &mut c, 8, 2);
    }
    assert_close(
        position_of(&bodies, hb).y,
        0.1,
        1e-6,
        "limit along A's local X",
    );
}

#[test]
fn test_limits_leave_motion_inside_the_range_untouched() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(
        &mut bodies,
        hb,
        Vec3::new(0.5, -0.25, 0.0),
        Vec3::new(0.0, 0.0, 0.5),
    );
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -1.0, 1.0).unwrap();
    c.set_linear_limit(1, -1.0, 1.0).unwrap();
    c.set_angular_limit(2, -1.0, 1.0).unwrap();

    c.prepare(&bodies, DT);
    c.solve_velocity(&mut bodies, DT);
    assert_eq!(velocity_of(&bodies, hb), Vec3::new(0.5, -0.25, 0.0));
    assert_eq!(
        bodies.get(hb).map(|b| b.angular_velocity),
        Some(Vec3::new(0.0, 0.0, 0.5))
    );
    assert_eq!(c.residual(), 0.0, "no row has anything to correct");

    // Position pass: inside the limits nothing moves either.
    c.solve_position(&mut bodies, DT);
    assert_eq!(position_of(&bodies, hb), Vec3::zeros());
}

#[test]
fn test_limited_body_held_against_a_constant_push() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();
    let mut max_x = f64::MIN;
    for _ in 0..180 {
        // A constant acceleration of 10 m/s² along +X.
        if let Some(b) = bodies.get_mut(hb) {
            b.velocity.x += 10.0 * DT;
        }
        step(&mut bodies, &mut c, 8, 2);
        max_x = max_x.max(position_of(&bodies, hb).x);
    }
    assert!(max_x <= 0.1 + 1e-9, "x = {max_x}");
    assert_close(
        position_of(&bodies, hb).x,
        0.1,
        1e-6,
        "pressed against the stop",
    );
}

#[test]
fn test_angular_limit_stops_rotation_at_both_bounds() {
    for sign in [1.0, -1.0] {
        let (mut bodies, ha, hb) = static_and_dynamic();
        set_velocity(
            &mut bodies,
            hb,
            Vec3::zeros(),
            Vec3::new(0.0, 0.0, 2.0 * sign),
        );
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_angular_limit(2, -0.1, 0.1).unwrap();
        let mut extreme = 0.0_f64;
        for _ in 0..60 {
            step(&mut bodies, &mut c, 8, 2);
            let angle = coordinate_of(&mut c, &bodies, 5);
            extreme = extreme.max(sign * angle);
        }
        assert!(
            extreme <= 0.1 + 1e-9,
            "sign {sign}: passed the limit, angle {extreme}"
        );
        assert_close(
            sign * coordinate_of(&mut c, &bodies, 5),
            0.1,
            1e-6,
            "resting angle",
        );
    }
}

#[test]
fn test_violated_limit_is_pushed_back_without_overshoot() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.position = Vec3::new(0.3, 0.0, 0.0);
    }
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();
    let mut lowest = f64::MAX;
    for _ in 0..40 {
        step(&mut bodies, &mut c, 8, 2);
        lowest = lowest.min(position_of(&bodies, hb).x);
    }
    assert!(
        lowest >= 0.1 - 1e-9,
        "the push-back overshot to x = {lowest}"
    );
    assert_close(position_of(&bodies, hb).x, 0.1, 1e-6, "back at the limit");
}

#[test]
fn test_violated_limit_only_stops_further_outward_motion_in_the_velocity_pass() {
    // (side of the violation, velocity along +X, expected velocity afterwards)
    let cases = [
        (0.3, 1.0, 0.0),   // beyond the upper bound, moving further out
        (0.3, -1.0, -1.0), // beyond the upper bound, moving back in
        (0.3, 0.0, 0.0),
        (-0.3, -1.0, 0.0), // beyond the lower bound, moving further out
        (-0.3, 1.0, 1.0),  // beyond the lower bound, moving back in
    ];
    for (x, velocity, expected) in cases {
        let (mut bodies, ha, hb) = static_and_dynamic();
        if let Some(b) = bodies.get_mut(hb) {
            b.transform.position = Vec3::new(x, 0.0, 0.0);
        }
        set_velocity(
            &mut bodies,
            hb,
            Vec3::new(velocity, 0.0, 0.0),
            Vec3::zeros(),
        );
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_linear_limit(0, -0.1, 0.1).unwrap();
        c.prepare(&bodies, DT);
        c.solve_velocity(&mut bodies, DT);
        assert_close(
            velocity_of(&bodies, hb).x,
            expected,
            1e-12,
            &format!("x = {x}, v = {velocity}"),
        );
    }
}

#[test]
fn test_position_pass_projects_onto_one_sided_limits() {
    // Beyond the upper bound of a linear and an angular limit: the pass moves
    // towards the bound and never past it.
    let (mut bodies, ha, hb) = static_and_dynamic();
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.position = Vec3::new(0.3, 0.0, 0.0);
    }
    rotate_body(&mut bodies, hb, Vec3::new(0.0, 0.0, 1.0), 0.3);
    // Angles are measured from the aligned configuration, so the 0.3 rad is a
    // real violation of the angular limit.
    let mut c = SixDofConstraint::new(ha, hb).with_reference_rotation(Quat::identity());
    c.set_linear_limit(0, -0.1, 0.1).unwrap();
    c.set_angular_limit(2, -0.1, 0.1).unwrap();
    for _ in 0..40 {
        c.solve_position(&mut bodies, DT);
    }
    let x = position_of(&bodies, hb).x;
    assert!((0.1 - 1e-9..0.1 + 1e-3).contains(&x), "x = {x}");
    let angle = coordinate_of(&mut c, &bodies, 5);
    assert!((0.1 - 1e-9..0.1 + 1e-3).contains(&angle), "angle = {angle}");

    // Below the lower bound the push goes the other way.
    let (mut bodies, ha, hb) = static_and_dynamic();
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.position = Vec3::new(-0.3, 0.0, 0.0);
    }
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.1, 0.1).unwrap();
    for _ in 0..40 {
        c.solve_position(&mut bodies, DT);
    }
    let x = position_of(&bodies, hb).x;
    assert!(x <= -0.1 + 1e-9 && x > -0.1 - 1e-3, "x = {x}");
}

#[test]
fn test_degenerate_limit_acts_like_a_lock() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(&mut bodies, hb, Vec3::new(2.0, 0.0, 0.0), Vec3::zeros());
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, 0.0, 0.0).unwrap();
    for _ in 0..30 {
        step(&mut bodies, &mut c, 8, 2);
    }
    assert!(position_of(&bodies, hb).x.abs() < 1e-3);
    assert!(velocity_of(&bodies, hb).x.abs() < 1e-3);
}

// ── F0479/F0480: configuration is validated, the hot path never panics ───────

#[test]
fn test_setters_validate_their_input() {
    let (_, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::new(ha, hb);

    assert!(matches!(
        c.set_linear_limit(0, 1.0, -1.0),
        Err(Error::InvalidJointLimits { .. })
    ));
    assert!(matches!(
        c.set_angular_limit(1, f64::NAN, 1.0),
        Err(Error::InvalidJointLimits { .. })
    ));
    assert!(matches!(
        c.set_linear_limit(3, -1.0, 1.0),
        Err(Error::ConfigurationError { .. })
    ));
    assert!(
        c.linear[0].limited.is_none(),
        "a rejected limit is not stored"
    );
    assert!(c.angular[1].limited.is_none());

    assert!(
        c.set_linear_limit(1, 0.5, 0.5).is_ok(),
        "min == max is allowed"
    );
    assert!(c.set_angular_limit(2, f64::NEG_INFINITY, 1.0).is_ok());

    assert!(matches!(
        c.set_linear_motor(0, MotorConfig::new(1.0, -5.0, 1.0)),
        Err(Error::ControllerError { .. })
    ));
    assert!(matches!(
        c.set_angular_motor(0, MotorConfig::new(f64::NAN, 5.0, 1.0)),
        Err(Error::ControllerError { .. })
    ));
    assert!(matches!(
        c.set_angular_motor(0, MotorConfig::new(1.0, f64::NAN, 1.0)),
        Err(Error::ControllerError { .. })
    ));
    assert!(matches!(
        c.set_angular_motor(5, MotorConfig::new(1.0, 5.0, 1.0)),
        Err(Error::ConfigurationError { .. })
    ));
    assert!(c.linear[0].motor.is_none() && c.angular[0].motor.is_none());
    assert!(
        c.set_linear_motor(0, MotorConfig::new(-1.0, 0.0, 1.0))
            .is_ok()
    );
}

#[test]
fn test_hot_path_survives_invalid_public_configuration() {
    // The configuration fields are public, so the setters cannot guard them.
    let (mut bodies, ha, hb) = static_and_dynamic();
    if let Some(b) = bodies.get_mut(hb) {
        b.transform.position = Vec3::new(1.5, 0.0, 0.0);
    }
    let mut c = SixDofConstraint::new(ha, hb);
    c.linear[0] = AxisConfig::limited(1.0, -1.0); // inverted: acts as [-1, 1]
    c.linear[1] = AxisConfig::limited(f64::NAN, 1.0); // unusable: ignored
    c.linear[2] = AxisConfig::limited(f64::NEG_INFINITY, f64::INFINITY);
    c.angular[0] = AxisConfig::free().with_motor(MotorConfig::new(f64::NAN, f64::NAN, 1.0));
    c.angular[2] = AxisConfig::free().with_motor(MotorConfig::new(1.0, -5.0, 1.0));

    for _ in 0..60 {
        step(&mut bodies, &mut c, 8, 2);
    }
    let body = bodies.get(hb);
    assert!(body.is_some_and(|b| b.velocity.iter().all(|v| v.is_finite())));
    assert!(body.is_some_and(|b| b.angular_velocity.iter().all(|v| v.is_finite())));
    let x = position_of(&bodies, hb).x;
    assert!(
        x > 0.99 && x < 1.01,
        "the inverted limit acts as [-1, 1], x = {x}"
    );
    let spin = body.map(|b| b.angular_velocity.z).unwrap_or_default();
    assert!(
        spin > 0.0,
        "a motor with a negative force limit uses |F|, got {spin}"
    );
}

#[test]
fn test_paused_or_invalid_time_step_keeps_everything_finite() {
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(dynamic_body(1.0, Vec3::zeros()));
        let hb = bodies.insert(dynamic_body(2.0, Vec3::new(1.0, 0.0, 0.0)));
        set_velocity(
            &mut bodies,
            ha,
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(0.1, 0.2, 0.3),
        );
        let mut c = SixDofConstraint::weld(ha, hb);
        c.set_linear_limit(0, -1.0, 1.0).unwrap();
        c.set_linear_motor(1, MotorConfig::new(1.0, 10.0, 1.0))
            .unwrap();
        c.prepare(&bodies, dt);
        for _ in 0..5 {
            c.solve_velocity(&mut bodies, dt);
        }
        c.solve_position(&mut bodies, dt);
        for handle in [ha, hb] {
            let body = bodies.get(handle);
            assert!(
                body.is_some_and(|b| b.velocity.iter().all(|v| v.is_finite())
                    && b.angular_velocity.iter().all(|v| v.is_finite())
                    && b.transform.position.iter().all(|v| v.is_finite())),
                "dt = {dt}: non-finite state"
            );
        }
    }
}

// ── Motors ───────────────────────────────────────────────────────────────────

#[test]
fn test_linear_motor_drives_the_coordinate_rate_to_its_target() {
    for target in [2.0, -2.0] {
        let (mut bodies, ha, hb) = static_and_dynamic();
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_linear_motor(0, MotorConfig::new(target, 1000.0, 1.0))
            .unwrap();
        c.prepare(&bodies, DT);
        c.solve_velocity(&mut bodies, DT);
        // The rate of the coordinate (velocity of B relative to A) reaches the target.
        assert_close(velocity_of(&bodies, hb).x, target, 1e-12, "motor target");
        assert_eq!(velocity_of(&bodies, hb).y, 0.0, "other axes are free");
    }
}

#[test]
fn test_motor_force_limit_bounds_the_impulse_of_a_step() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    // 30 N for 1/60 s is an impulse of 0.5 N·s: 0.5 m/s for a unit mass.
    c.set_linear_motor(0, MotorConfig::new(5.0, 30.0, 1.0))
        .unwrap();
    c.prepare(&bodies, DT);
    for _ in 0..10 {
        c.solve_velocity(&mut bodies, DT);
    }
    assert_close(
        velocity_of(&bodies, hb).x,
        0.5,
        1e-12,
        "force-limited speed",
    );
    // The next step (fresh accumulators) adds another 0.5 m/s.
    c.prepare(&bodies, DT);
    for _ in 0..10 {
        c.solve_velocity(&mut bodies, DT);
    }
    assert_close(velocity_of(&bodies, hb).x, 1.0, 1e-12, "second step");
}

#[test]
fn test_angular_motor_turns_body_b_in_the_positive_sense() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_angular_motor(2, MotorConfig::new(3.0, 1000.0, 1.0))
        .unwrap();
    for _ in 0..30 {
        step(&mut bodies, &mut c, 4, 0);
    }
    let spin = bodies
        .get(hb)
        .map(|b| b.angular_velocity)
        .unwrap_or_default();
    assert_close(spin.z, 3.0, 1e-9, "angular rate");
    // 0.5 s at 3 rad/s.
    assert_close(coordinate_of(&mut c, &bodies, 5), 1.5, 1e-6, "angle");
}

#[test]
fn test_motor_driving_into_a_limit_stops_at_the_limit() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_limit(0, -0.2, 0.2).unwrap();
    c.set_linear_motor(0, MotorConfig::new(5.0, 1000.0, 1.0))
        .unwrap();
    let mut max_x = f64::MIN;
    for _ in 0..150 {
        step(&mut bodies, &mut c, 8, 2);
        max_x = max_x.max(position_of(&bodies, hb).x);
    }
    assert!(
        max_x <= 0.2 + 1e-9,
        "the motor pushed through the limit: {max_x}"
    );
    assert_close(position_of(&bodies, hb).x, 0.2, 1e-6, "held at the limit");
    assert!(velocity_of(&bodies, hb).x.abs() < 1e-6);
}

#[test]
fn test_locked_axis_ignores_its_motor() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(0);
    c.linear[0].motor = Some(MotorConfig::new(4.0, 1000.0, 1.0));
    c.prepare(&bodies, DT);
    c.solve_velocity(&mut bodies, DT);
    assert_eq!(velocity_of(&bodies, hb).x, 0.0);
}

// ── Hooks: residual and relaxation ───────────────────────────────────────────

#[test]
fn test_residual_reports_the_velocity_update() {
    let (mut bodies, ha, hb) = static_and_dynamic();
    set_velocity(&mut bodies, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::zeros());
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(0);
    assert_eq!(c.residual(), 0.0, "nothing solved yet");
    c.prepare(&bodies, DT);
    c.solve_velocity(&mut bodies, DT);
    assert_close(c.residual(), 1.0, 1e-12, "the whole 1 m/s was removed");
    c.solve_velocity(&mut bodies, DT);
    assert!(
        c.residual() < 1e-12,
        "fixed point reached: {}",
        c.residual()
    );
}

#[test]
fn test_relaxed_solve_scales_the_impulse_update() {
    for (omega, expected_velocity) in [(0.5, 0.5), (1.0, 0.0), (1.5, -0.5)] {
        let (mut bodies, ha, hb) = static_and_dynamic();
        set_velocity(&mut bodies, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::zeros());
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_linear(0);
        c.prepare(&bodies, DT);
        let update = c.solve_velocity_relaxed(&mut bodies, DT, omega);
        assert_close(
            velocity_of(&bodies, hb).x,
            expected_velocity,
            1e-12,
            &format!("omega {omega}"),
        );
        assert_close(
            update,
            omega,
            1e-12,
            &format!("update size for omega {omega}"),
        );
    }
    // Out-of-range and NaN factors are sanitised, never propagated.
    for omega in [0.0, 100.0, f64::NAN] {
        let (mut bodies, ha, hb) = static_and_dynamic();
        set_velocity(&mut bodies, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::zeros());
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_linear(0);
        c.prepare(&bodies, DT);
        c.solve_velocity_relaxed(&mut bodies, DT, omega);
        assert!(velocity_of(&bodies, hb).x.is_finite());
    }
}

#[test]
fn test_prepare_resets_the_accumulated_impulses() {
    // A stale accumulator must not survive into the next step: with a motor
    // clamp it would shrink the usable impulse, with a one-sided limit it
    // would turn into a pulling impulse.
    let (mut bodies, ha, hb) = static_and_dynamic();
    let mut c = SixDofConstraint::new(ha, hb);
    c.set_linear_motor(0, MotorConfig::new(5.0, 30.0, 1.0))
        .unwrap();
    for _ in 0..4 {
        c.prepare(&bodies, DT);
        c.solve_velocity(&mut bodies, DT);
    }
    // Every step contributes the full 0.5 m/s.
    assert_close(
        velocity_of(&bodies, hb).x,
        2.0,
        1e-12,
        "four steps of 0.5 m/s",
    );
}

// ── Physical consistency ─────────────────────────────────────────────────────

#[test]
fn test_weld_impulses_conserve_momentum_and_satisfy_the_weld() {
    let mut bodies = RigidBodySet::new();
    let ha = bodies.insert(dynamic_body(1.0, Vec3::zeros()));
    let hb = bodies.insert(dynamic_body(3.0, Vec3::new(1.0, 0.0, 0.0)));
    set_velocity(
        &mut bodies,
        ha,
        Vec3::new(1.0, 2.0, -1.0),
        Vec3::new(0.3, -0.2, 0.5),
    );
    set_velocity(
        &mut bodies,
        hb,
        Vec3::new(-2.0, 0.5, 1.0),
        Vec3::new(-0.1, 0.4, 0.2),
    );
    // Both anchors are the world point (0.5, 0, 0).
    let mut c = SixDofConstraint::weld(ha, hb)
        .with_anchors(Vec3::new(0.5, 0.0, 0.0), Vec3::new(-0.5, 0.0, 0.0));

    let momenta = |bodies: &RigidBodySet| {
        let mut linear = Vec3::zeros();
        let mut angular = Vec3::zeros();
        for (_, body) in bodies.iter() {
            let mass = 1.0 / body.inverse_mass;
            linear += body.velocity * mass;
            angular += body.world_inertia() * body.angular_velocity
                + body.transform.position.cross(&(body.velocity * mass));
        }
        (linear, angular)
    };
    let (linear_before, angular_before) = momenta(&bodies);

    c.prepare(&bodies, DT);
    for _ in 0..300 {
        c.solve_velocity(&mut bodies, DT);
    }
    let (linear_after, angular_after) = momenta(&bodies);
    assert!(
        (linear_after - linear_before).norm() < 1e-10,
        "linear momentum"
    );
    assert!(
        (angular_after - angular_before).norm() < 1e-10,
        "angular momentum"
    );

    // The weld holds: no relative velocity at the anchor, no relative spin.
    let a = bodies.get(ha).map(|b| (b.velocity, b.angular_velocity));
    let b = bodies.get(hb).map(|b| (b.velocity, b.angular_velocity));
    let (va, wa) = a.unwrap_or_default();
    let (vb, wb) = b.unwrap_or_default();
    let anchor_a = va + wa.cross(&Vec3::new(0.5, 0.0, 0.0));
    let anchor_b = vb + wb.cross(&Vec3::new(-0.5, 0.0, 0.0));
    assert!(
        (anchor_a - anchor_b).norm() < 1e-9,
        "anchors still move apart"
    );
    assert!(
        (wa - wb).norm() < 1e-9,
        "bodies still spin against each other"
    );
}

#[test]
fn test_weld_carries_a_hanging_load() {
    let mut bodies = RigidBodySet::new();
    let ha = bodies.insert(RigidBody::new_static());
    let hb = bodies.insert(dynamic_body(2.0, Vec3::zeros()));
    let mut c = SixDofConstraint::weld(ha, hb);
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    for _ in 0..120 {
        if let Some(b) = bodies.get_mut(hb) {
            b.integrate_forces(DT, &gravity);
        }
        step(&mut bodies, &mut c, 8, 2);
    }
    assert!(position_of(&bodies, hb).norm() < 1e-6, "the load sagged");
    assert!(relative_orientation_angle(&bodies, ha, hb) < 1e-6);
}

#[test]
fn test_ball_joint_pendulum_keeps_its_anchors_together() {
    // Body B hangs from the origin by an anchor 0.5 m from its centre; only
    // the linear DOFs are locked, so it swings.
    let mut bodies = RigidBodySet::new();
    let ha = bodies.insert(RigidBody::new_static());
    let hb = bodies.insert(dynamic_body(1.0, Vec3::new(0.5, 0.0, 0.0)));
    let mut c =
        SixDofConstraint::new(ha, hb).with_anchors(Vec3::zeros(), Vec3::new(-0.5, 0.0, 0.0));
    for axis in 0..3 {
        c.lock_linear(axis);
    }
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let mut worst_gap = 0.0_f64;
    let mut lowest = 0.0_f64;
    let mut largest_swing = 0.0_f64;
    for _ in 0..240 {
        if let Some(b) = bodies.get_mut(hb) {
            b.integrate_forces(DT, &gravity);
        }
        step(&mut bodies, &mut c, 20, 4);
        worst_gap = worst_gap.max(anchor_gap(&bodies, &c));
        lowest = lowest.min(position_of(&bodies, hb).y);
        largest_swing = largest_swing.max(relative_orientation_angle(&bodies, ha, hb));
    }
    assert!(worst_gap < 2e-3, "the anchors drifted apart by {worst_gap}");
    assert!(
        lowest < -0.3,
        "the pendulum must swing down, lowest y = {lowest}"
    );
    assert!(
        largest_swing > 0.5,
        "the free angular DOFs let the body rotate, largest swing {largest_swing}"
    );
}

#[test]
fn test_pgs_solvers_converge_on_a_weld() {
    use crate::pgs_solver::{PgsSolver, SorPgsSolver};

    let make = || {
        let (mut bodies, ha, hb) = coincident_bodies();
        set_velocity(
            &mut bodies,
            ha,
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(0.1, 0.2, 0.3),
        );
        set_velocity(
            &mut bodies,
            hb,
            Vec3::new(-1.0, 0.0, 1.0),
            Vec3::new(-0.2, 0.1, 0.4),
        );
        let constraints: Vec<Box<dyn Constraint>> = vec![Box::new(SixDofConstraint::weld(ha, hb))];
        (bodies, ha, hb, constraints)
    };
    let welded = |bodies: &RigidBodySet, ha: BodyHandle, hb: BodyHandle| {
        let a = bodies.get(ha).map(|b| (b.velocity, b.angular_velocity));
        let b = bodies.get(hb).map(|b| (b.velocity, b.angular_velocity));
        let ((va, wa), (vb, wb)) = (a.unwrap_or_default(), b.unwrap_or_default());
        ((va - vb).norm(), (wa - wb).norm())
    };

    let (mut bodies, ha, hb, mut constraints) = make();
    let stats = PgsSolver::new(50, 2).solve(&mut constraints, &mut bodies, DT);
    assert!(stats.converged, "residual = {}", stats.residual);
    assert!(
        stats.iterations_used <= 3,
        "decoupled rows are solved by the first sweep, used {}",
        stats.iterations_used
    );
    let (linear, angular) = welded(&bodies, ha, hb);
    assert!(linear < 1e-9 && angular < 1e-9, "{linear}, {angular}");

    // The relaxed solver goes through `solve_velocity_relaxed`.
    let (mut bodies, ha, hb, mut constraints) = make();
    let stats = SorPgsSolver::new(80, 2, 1.5).solve(&mut constraints, &mut bodies, DT);
    assert!(stats.converged, "residual = {}", stats.residual);
    let (linear, angular) = welded(&bodies, ha, hb);
    assert!(linear < 1e-5 && angular < 1e-5, "{linear}, {angular}");
}

#[test]
fn test_dof_count_counts_locked_and_limited_axes() {
    let (_, ha, hb) = make_two_bodies();
    assert_eq!(SixDofConstraint::new(ha, hb).dof_count(), 0);
    assert_eq!(SixDofConstraint::weld(ha, hb).dof_count(), 6);
    let mut c = SixDofConstraint::new(ha, hb);
    c.lock_linear(1);
    c.set_angular_limit(0, -0.5, 0.5).unwrap();
    c.set_linear_motor(2, MotorConfig::new(1.0, 1.0, 1.0))
        .unwrap();
    assert_eq!(c.dof_count(), 2, "a motor alone constrains nothing");
}

#[test]
fn test_missing_bodies_are_ignored() {
    let (mut bodies, ha, hb) = make_two_bodies();
    let mut c = SixDofConstraint::weld(ha, hb);
    bodies.remove(hb);
    c.prepare(&bodies, DT);
    c.solve_velocity(&mut bodies, DT);
    c.solve_position(&mut bodies, DT);
    assert_eq!(c.residual(), 0.0);
    assert_eq!(c.coordinates(), [0.0; 6]);
}
