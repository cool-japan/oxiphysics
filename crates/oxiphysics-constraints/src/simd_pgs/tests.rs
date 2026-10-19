//! Tests of the SoA batch PGS solver.

use super::*;
use crate::contact::ContactConstraint;
use crate::pgs_solver::PgsSolver;
use crate::traits::Constraint;
use oxiphysics_core::math::{Mat3 as NaMat3, Vec3};
use oxiphysics_rigid::{RigidBody, RigidBodySet};

const DT: f64 = 1.0 / 60.0;

type BodyVelocity = ([f64; 3], [f64; 3]);

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn diag(x: f64, y: f64, z: f64) -> [[f64; 3]; 3] {
    [[x, 0.0, 0.0], [0.0, y, 0.0], [0.0, 0.0, z]]
}

fn unit_y() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: expected {expected}, got {actual}"
    );
}

fn assert_vec_close(actual: [f64; 3], expected: [f64; 3], tolerance: f64, what: &str) {
    for axis in 0..3 {
        assert_close(
            actual[axis],
            expected[axis],
            tolerance,
            &format!("{what}[{axis}]"),
        );
    }
}

fn build(rows: &[SoaConstraintRow], velocities: &[BodyVelocity]) -> SoaConstraints {
    let mut soa = SoaConstraints::new();
    for row in rows {
        soa.push(*row);
    }
    for (body, (linear, angular)) in velocities.iter().enumerate() {
        soa.set_body_velocity(body, *linear, *angular);
    }
    soa
}

fn velocities_of(soa: &SoaConstraints) -> Vec<BodyVelocity> {
    (0..soa.body_count())
        .filter_map(|body| soa.body_velocity(body))
        .collect()
}

/// Apply the impulse `delta_lambda` of `row` to `velocities` straight from
/// the physics: `P = Δλ n` acts on B at `r_b` and `-P` on A at `r_a`.
fn reference_apply(row: &SoaConstraintRow, delta_lambda: f64, velocities: &mut [BodyVelocity]) {
    let impulse = scale(row.normal, delta_lambda);
    let (v_b, w_b) = velocities[row.body_b];
    velocities[row.body_b] = (
        add(v_b, scale(impulse, row.inv_mass_b)),
        add(w_b, mat_vec(&row.inv_inertia_b, cross(row.r_b, impulse))),
    );
    let (v_a, w_a) = velocities[row.body_a];
    velocities[row.body_a] = (
        sub(v_a, scale(impulse, row.inv_mass_a)),
        sub(w_a, mat_vec(&row.inv_inertia_a, cross(row.r_a, impulse))),
    );
}

/// Independent Gauss-Seidel reference: velocities of the two contact
/// points from `v + ω × r`, the effective mass from its definition, rows
/// processed strictly one after the other.
fn reference_solve(
    rows: &[SoaConstraintRow],
    velocities: &mut [BodyVelocity],
    iterations: usize,
    omega: f64,
) -> Vec<f64> {
    let mut lambda: Vec<f64> = rows
        .iter()
        .map(|row| project(row.warm_lambda, row.lo, row.hi))
        .collect();
    for (row, &warm) in rows.iter().zip(&lambda) {
        reference_apply(row, warm, velocities);
    }
    for _ in 0..iterations {
        for (index, row) in rows.iter().enumerate() {
            let (v_a, w_a) = velocities[row.body_a];
            let (v_b, w_b) = velocities[row.body_b];
            let point_a = add(v_a, cross(w_a, row.r_a));
            let point_b = add(v_b, cross(w_b, row.r_b));
            let rel_vel = dot(row.normal, sub(point_b, point_a));
            let rxn_a = cross(row.r_a, row.normal);
            let rxn_b = cross(row.r_b, row.normal);
            let k = row.inv_mass_a
                + row.inv_mass_b
                + dot(rxn_a, mat_vec(&row.inv_inertia_a, rxn_a))
                + dot(rxn_b, mat_vec(&row.inv_inertia_b, rxn_b));
            let eff_mass = if k > 1e-15 { 1.0 / k } else { 0.0 };
            let old = lambda[index];
            let new = project(
                old - eff_mass * (rel_vel + row.bias) * omega,
                row.lo,
                row.hi,
            );
            lambda[index] = new;
            reference_apply(row, new - old, velocities);
        }
    }
    lambda
}

/// Solve `soa` with the batch solver and the reference and compare both
/// the velocities and the accumulated impulses.
fn assert_matches_reference(
    rows: &[SoaConstraintRow],
    velocities: &[BodyVelocity],
    solver: &BatchPgsSolver,
    what: &str,
) {
    let mut soa = build(rows, velocities);
    solver.solve(&mut soa, DT);
    let mut expected_velocities = velocities.to_vec();
    let expected_lambda = reference_solve(
        rows,
        &mut expected_velocities,
        solver.max_iterations,
        solver.omega,
    );
    for (body, (v, w)) in velocities_of(&soa).iter().enumerate() {
        assert_vec_close(
            *v,
            expected_velocities[body].0,
            1e-11,
            &format!("{what}: v[{body}]"),
        );
        assert_vec_close(
            *w,
            expected_velocities[body].1,
            1e-11,
            &format!("{what}: w[{body}]"),
        );
    }
    for (row, (&got, &want)) in soa.lambda.iter().zip(&expected_lambda).enumerate() {
        assert_close(got, want, 1e-11, &format!("{what}: lambda[{row}]"));
    }
}

/// A box (body 1) over a static ground (body 0) touching it in its four
/// bottom corners.
fn box_on_four_corners(box_velocity: BodyVelocity) -> (Vec<SoaConstraintRow>, Vec<BodyVelocity>) {
    let inv_inertia = diag(1.5, 2.0, 2.5);
    let rows = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
        .iter()
        .map(|&(x, z)| SoaConstraintRow {
            r_a: [x, 0.0, z],
            r_b: [x, -0.5, z],
            inv_inertia_b: inv_inertia,
            ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 0.5)
        })
        .collect();
    (rows, vec![([0.0; 3], [0.0; 3]), box_velocity])
}

// ── Basics (kept from the original suite) ─────────────────────────────────

#[test]
fn empty_constraints_converge_immediately() {
    let mut soa = SoaConstraints::new();
    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    assert!(stats.converged);
    assert_eq!(stats.iterations_used, 0);
    assert_eq!(stats.residual, 0.0);
}

#[test]
fn single_contact_impulse_nonnegative() {
    let mut builder = SoaContactBuilder::new();
    builder.add_contact(
        SoaContact::centered(0, 1, unit_y(), 0.05, 1.0, 1.0),
        DT,
        0.2,
    );
    let mut soa = builder.build();
    // Give body A a downward velocity so it violates the constraint
    soa.set_body_velocity(0, [0.0, -1.0, 0.0], [0.0; 3]);

    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    // Lambda must be non-negative for a contact constraint
    assert!(soa.lambda[0] >= 0.0);
    assert!(stats.iterations_used >= 1);
}

// ── F0474: lever arms and angular velocity ────────────────────────────────

#[test]
fn push_precomputes_the_angular_terms_and_effective_mass() {
    let mut soa = SoaConstraints::new();
    soa.push(SoaConstraintRow {
        r_a: [0.0, -1.0, 0.0],
        r_b: [1.0, 0.0, 0.0],
        inv_inertia_a: diag(1.0, 1.0, 3.0),
        inv_inertia_b: diag(2.0, 3.0, 4.0),
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.25, 0.5)
    });
    // r_a × n = (0,-1,0) × (0,1,0) = 0 ; r_b × n = (1,0,0) × (0,1,0) = (0,0,1)
    assert_vec_close(
        [soa.rxn_ax[0], soa.rxn_ay[0], soa.rxn_az[0]],
        [0.0; 3],
        1e-15,
        "r_a x n",
    );
    assert_vec_close(
        [soa.rxn_bx[0], soa.rxn_by[0], soa.rxn_bz[0]],
        [0.0, 0.0, 1.0],
        1e-15,
        "r_b x n",
    );
    assert_vec_close(
        [soa.irxn_bx[0], soa.irxn_by[0], soa.irxn_bz[0]],
        [0.0, 0.0, 4.0],
        1e-15,
        "I_b^-1 (r_b x n)",
    );
    // K = 0.25 + 0.5 + 0 + (0,0,1)·(0,0,4) = 4.75
    assert_close(soa.eff_mass[0], 1.0 / 4.75, 1e-15, "effective mass");
}

#[test]
fn off_centre_contact_spins_the_box_like_the_scalar_solver() {
    // A 2 kg box falls onto static ground and touches it at an offset of
    // (0.4, 0, 0.2) from its centre line.
    let inv_inertia = NaMat3::from_diagonal(&Vec3::new(1.5, 2.0, 2.5));
    let box_position = Vec3::new(0.0, 0.5, 0.0);
    let contact = Vec3::new(0.4, 0.0, 0.2);
    let velocity = Vec3::new(0.1, -1.0, 0.0);
    let spin = Vec3::new(0.3, 0.2, -0.1);

    // Scalar solver on a ContactConstraint (normal from B=ground to A=box).
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let mut box_body = RigidBody::new(2.0);
    box_body.transform.position = box_position;
    box_body.velocity = velocity;
    box_body.angular_velocity = spin;
    box_body.world_inverse_inertia = inv_inertia;
    let box_handle = bodies.insert(box_body);
    let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(ContactConstraint::new(
        box_handle,
        ground,
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        contact,
        contact,
        0.0,
        0.0,
    ))];
    PgsSolver::with_options(1, 0, 0.0, false).solve(&mut constraints, &mut bodies, DT);
    let reference = bodies
        .get(box_handle)
        .map(|b| (b.velocity, b.angular_velocity));
    let (ref_v, ref_w) = reference.unwrap_or_default();

    // The same contact in the SoA solver: A = ground, B = box, n = A→B.
    let r_box = contact - box_position;
    let row = SoaConstraintRow {
        r_a: [contact.x, contact.y, contact.z],
        r_b: [r_box.x, r_box.y, r_box.z],
        inv_inertia_b: [[1.5, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.5]],
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 0.5)
    };
    let mut soa = build(
        &[row],
        &[
            ([0.0; 3], [0.0; 3]),
            (
                [velocity.x, velocity.y, velocity.z],
                [spin.x, spin.y, spin.z],
            ),
        ],
    );
    BatchPgsSolver::new(1, 0.0, 1.0).solve(&mut soa, DT);

    let (v, w) = soa.body_velocity(1).unwrap_or_default();
    assert_vec_close(v, [ref_v.x, ref_v.y, ref_v.z], 1e-12, "linear velocity");
    assert_vec_close(w, [ref_w.x, ref_w.y, ref_w.z], 1e-12, "angular velocity");
    // The impulse at +x pushes the box up there, which spins it about +z.
    assert!(
        w[2] - spin.z > 0.01,
        "an off-centre contact must induce rotation, got ω_z = {}",
        w[2]
    );
    assert_vec_close(
        soa.body_velocity(0).unwrap_or_default().1,
        [0.0; 3],
        0.0,
        "the static ground must not rotate",
    );
}

#[test]
fn four_corner_contacts_match_the_scalar_solver_iteration_by_iteration() {
    let inv_inertia = NaMat3::from_diagonal(&Vec3::new(1.5, 2.0, 2.5));
    let velocity = [0.2, -1.5, 0.1];
    let spin = [0.3, 0.1, -0.2];
    let (rows, velocities) = box_on_four_corners((velocity, spin));

    for iterations in [1, 2, 5, 40] {
        // Scalar PGS with ContactConstraints in the same row order.
        let mut bodies = RigidBodySet::new();
        let ground = bodies.insert(RigidBody::new_static());
        let mut box_body = RigidBody::new(2.0);
        box_body.transform.position = Vec3::new(0.0, 0.5, 0.0);
        box_body.velocity = Vec3::new(velocity[0], velocity[1], velocity[2]);
        box_body.angular_velocity = Vec3::new(spin[0], spin[1], spin[2]);
        box_body.world_inverse_inertia = inv_inertia;
        let box_handle = bodies.insert(box_body);
        let mut constraints: Vec<Box<dyn Constraint>> =
            [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
                .iter()
                .map(|&(x, z)| {
                    let point = Vec3::new(x, 0.0, z);
                    Box::new(ContactConstraint::new(
                        box_handle,
                        ground,
                        Vec3::new(0.0, 1.0, 0.0),
                        0.0,
                        point,
                        point,
                        0.0,
                        0.0,
                    )) as Box<dyn Constraint>
                })
                .collect();
        PgsSolver::with_options(iterations, 0, 0.0, false).solve(&mut constraints, &mut bodies, DT);
        let (ref_v, ref_w) = bodies
            .get(box_handle)
            .map(|b| (b.velocity, b.angular_velocity))
            .unwrap_or_default();

        let mut soa = build(&rows, &velocities);
        assert!(
            !soa.lanes_are_independent(0),
            "four rows on one box share a moving body"
        );
        BatchPgsSolver::new(iterations, 0.0, 1.0).solve(&mut soa, DT);
        let (v, w) = soa.body_velocity(1).unwrap_or_default();
        assert_vec_close(
            v,
            [ref_v.x, ref_v.y, ref_v.z],
            1e-11,
            &format!("v after {iterations}"),
        );
        assert_vec_close(
            w,
            [ref_w.x, ref_w.y, ref_w.z],
            1e-11,
            &format!("w after {iterations}"),
        );
    }
}

#[test]
fn converged_corner_contacts_stop_the_box() {
    let (rows, velocities) = box_on_four_corners(([0.0, -1.5, 0.0], [0.4, 0.0, -0.3]));
    let mut soa = build(&rows, &velocities);
    let stats = BatchPgsSolver::new(500, 1e-10, 1.0).solve(&mut soa, DT);
    assert!(stats.converged, "residual = {}", stats.residual);
    let (v, w) = soa.body_velocity(1).unwrap_or_default();
    // Every corner is at rest: the whole box has stopped and stopped turning.
    assert!(v[1].abs() < 1e-8, "box still moving: vy = {}", v[1]);
    assert_vec_close(w, [0.0; 3], 1e-8, "angular velocity of the supported box");
    // Impulses are never pulling.
    assert!(soa.lambda.iter().all(|&l| l >= 0.0));
    // Momentum balance: the box (2 kg) lost 3 N·s along y.
    assert_close(soa.lambda.iter().sum::<f64>(), 3.0, 1e-7, "total impulse");
}

// ── Batches are plain Gauss-Seidel ────────────────────────────────────────

/// Eight independent pairs of dynamic bodies, rows with lever arms,
/// inertia and bias, all approaching.
fn independent_pairs() -> (Vec<SoaConstraintRow>, Vec<BodyVelocity>) {
    let mut rows = Vec::new();
    let mut velocities = Vec::new();
    let normals = [
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.6, 0.8, 0.0],
    ];
    for pair in 0..8 {
        let normal = normals[pair % 4];
        let s = pair as f64;
        rows.push(SoaConstraintRow {
            r_a: [0.1 * s, -0.2, 0.05],
            r_b: [-0.3, 0.1 * s, 0.2],
            inv_inertia_a: diag(1.0 + s, 2.0, 3.0),
            inv_inertia_b: diag(2.0, 1.0 + 0.5 * s, 1.5),
            bias: -0.05 * s,
            ..SoaConstraintRow::linear(normal, 2 * pair, 2 * pair + 1, 1.0, 0.5 + 0.25 * s)
        });
        velocities.push((scale(normal, 1.0 + s), [0.1, -0.2 * s, 0.3]));
        velocities.push((scale(normal, -0.5), [0.0, 0.4, -0.1 * s]));
    }
    (rows, velocities)
}

#[test]
fn independent_groups_are_batched_and_equal_to_plain_gauss_seidel() {
    let (rows, velocities) = independent_pairs();
    let soa = build(&rows, &velocities);
    assert!(soa.lanes_are_independent(0) && soa.lanes_are_independent(4));
    for omega in [1.0, 0.7, 1.4] {
        for iterations in [1, 3, 25] {
            assert_matches_reference(
                &rows,
                &velocities,
                &BatchPgsSolver::new(iterations, 0.0, omega),
                &format!("independent pairs, omega {omega}, {iterations} iterations"),
            );
        }
    }
}

#[test]
fn rows_sharing_a_moving_body_are_solved_sequentially() {
    // Eight redundant contacts between the same two dynamic bodies: a batch
    // evaluation (all rows read the same stale velocity) would apply the
    // impulse eight times.
    let rows: Vec<SoaConstraintRow> = (0..8)
        .map(|_| SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0))
        .collect();
    let velocities = [([0.0; 3], [0.0; 3]), ([0.0, -2.0, 0.0], [0.0; 3])];
    let mut soa = build(&rows, &velocities);
    assert!(!soa.lanes_are_independent(0));
    assert!(!soa.lanes_are_independent(4));

    let stats = BatchPgsSolver::new(1, 0.0, 1.0).solve(&mut soa, DT);
    // The first row removes the whole approach velocity (two equal masses
    // meet at the average velocity), the other seven see rv = 0 and add
    // nothing.
    assert_close(soa.lambda[0], 1.0, 1e-15, "first row impulse");
    assert!(soa.lambda[1..].iter().all(|&l| l == 0.0));
    let (va, _) = soa.body_velocity(0).unwrap_or_default();
    let (vb, _) = soa.body_velocity(1).unwrap_or_default();
    assert_vec_close(va, [0.0, -1.0, 0.0], 1e-15, "body A");
    assert_vec_close(vb, [0.0, -1.0, 0.0], 1e-15, "body B");
    assert_close(stats.residual, 1.0, 1e-15, "residual of the only iteration");
}

#[test]
fn rows_sharing_only_a_static_body_are_batched() {
    // Eight boxes resting on the same static ground (body 0).
    let rows: Vec<SoaConstraintRow> = (1..=8)
        .map(|body| SoaConstraintRow {
            r_b: [0.0, -0.5, 0.0],
            inv_inertia_b: diag(2.0, 2.0, 2.0),
            ..SoaConstraintRow::linear(unit_y(), 0, body, 0.0, 1.0)
        })
        .collect();
    let mut velocities = vec![([0.0; 3], [0.0; 3])];
    for body in 1..=8 {
        let s = body as f64;
        velocities.push(([0.1 * s, -s, 0.0], [0.0, 0.2 * s, 0.1]));
    }
    let soa = build(&rows, &velocities);
    assert!(soa.lanes_are_independent(0) && soa.lanes_are_independent(4));
    assert_matches_reference(
        &rows,
        &velocities,
        &BatchPgsSolver::new(4, 0.0, 1.0),
        "boxes on shared ground",
    );
}

#[test]
fn lane_independence_looks_at_moved_bodies_only() {
    // Rows 0..4: A = static 0, B = bodies 1..5 -> independent.
    let mut soa = SoaConstraints::new();
    for body in 1..=4 {
        soa.push(SoaConstraintRow::linear(unit_y(), 0, body, 0.0, 1.0));
    }
    assert!(soa.lanes_are_independent(0));

    // Row 3 now reads body 1, which row 0 moves.
    let mut soa = SoaConstraints::new();
    for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 4)] {
        soa.push(SoaConstraintRow::linear(unit_y(), a, b, 0.0, 1.0));
    }
    assert!(!soa.lanes_are_independent(0));

    // A row whose body A moves (positive inverse mass) conflicts with any
    // other row that touches that body, even if the other row treats it as
    // static.
    let mut soa = SoaConstraints::new();
    soa.push(SoaConstraintRow::linear(unit_y(), 5, 6, 1.0, 1.0));
    soa.push(SoaConstraintRow::linear(unit_y(), 7, 5, 1.0, 0.0));
    soa.push(SoaConstraintRow::linear(unit_y(), 8, 9, 1.0, 1.0));
    soa.push(SoaConstraintRow::linear(unit_y(), 10, 11, 1.0, 1.0));
    assert!(!soa.lanes_are_independent(0));

    // A body that only rotates counts as moving too.
    let mut soa = SoaConstraints::new();
    let spinner = |a, b| SoaConstraintRow {
        r_b: [1.0, 0.0, 0.0],
        inv_inertia_b: diag(1.0, 1.0, 1.0),
        ..SoaConstraintRow::linear(unit_y(), a, b, 0.0, 0.0)
    };
    soa.push(spinner(0, 1));
    soa.push(spinner(2, 1));
    soa.push(spinner(3, 4));
    soa.push(spinner(5, 6));
    assert!(!soa.lanes_are_independent(0));
}

#[test]
fn stack_of_ten_boxes_matches_the_scalar_solver() {
    let boxes = 10;
    let mut bodies = RigidBodySet::new();
    let ground = bodies.insert(RigidBody::new_static());
    let mut handles = vec![ground];
    for level in 0..boxes {
        let mut body = RigidBody::new(1.0);
        body.transform.position = Vec3::new(0.0, 0.5 + level as f64, 0.0);
        body.velocity = Vec3::new(0.0, -1.0, 0.0);
        handles.push(bodies.insert(body));
    }
    let mut constraints: Vec<Box<dyn Constraint>> = Vec::new();
    let mut rows = Vec::new();
    for level in 0..boxes {
        let point = Vec3::new(0.0, level as f64, 0.0);
        constraints.push(Box::new(ContactConstraint::new(
            handles[level + 1],
            handles[level],
            Vec3::new(0.0, 1.0, 0.0),
            0.0,
            point,
            point,
            0.0,
            0.0,
        )));
        let inv_mass_below = if level == 0 { 0.0 } else { 1.0 };
        rows.push(SoaConstraintRow {
            r_a: [0.0, 0.5, 0.0],
            r_b: [0.0, -0.5, 0.0],
            inv_inertia_a: if level == 0 {
                diag(0.0, 0.0, 0.0)
            } else {
                diag(1.0, 1.0, 1.0)
            },
            inv_inertia_b: diag(1.0, 1.0, 1.0),
            ..SoaConstraintRow::linear(unit_y(), level, level + 1, inv_mass_below, 1.0)
        });
    }
    let mut velocities = vec![([0.0; 3], [0.0; 3])];
    velocities.extend((0..boxes).map(|_| ([0.0, -1.0, 0.0], [0.0; 3])));

    let iterations = 4;
    PgsSolver::with_options(iterations, 0, 0.0, false).solve(&mut constraints, &mut bodies, DT);
    let mut soa = build(&rows, &velocities);
    BatchPgsSolver::new(iterations, 0.0, 1.0).solve(&mut soa, DT);
    for (level, &handle) in handles.iter().enumerate().skip(1) {
        let expected = bodies.get(handle).map(|b| b.velocity).unwrap_or_default();
        let (v, _) = soa.body_velocity(level).unwrap_or_default();
        assert_vec_close(
            v,
            [expected.x, expected.y, expected.z],
            1e-11,
            &format!("box {level}"),
        );
    }
    // After many more sweeps the stack is at rest and row k carries the
    // weight of the boxes above it.
    let mut soa = build(&rows, &velocities);
    let stats = BatchPgsSolver::new(3000, 1e-12, 1.0).solve(&mut soa, DT);
    assert!(stats.converged, "residual = {}", stats.residual);
    for level in 0..boxes {
        assert_close(
            soa.lambda[level],
            (boxes - level) as f64,
            1e-7,
            &format!("impulse of contact {level}"),
        );
    }
}

// ── F0476: warm start ─────────────────────────────────────────────────────

#[test]
fn warm_lambda_is_applied_to_the_bodies_before_iterating() {
    // The bodies already separate at 1 m/s; a warm start of 2 N·s is
    // applied first (A -2, B +2 m/s), the iteration removes it again and
    // the clamp stops at zero: the net impulse is zero, never negative.
    let row = SoaConstraintRow {
        warm_lambda: 2.0,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0)
    };
    let mut soa = build(&[row], &[([0.0; 3], [0.0; 3]), ([0.0, 1.0, 0.0], [0.0; 3])]);
    assert_eq!(soa.lambda[0], 2.0);
    BatchPgsSolver::new(5, 0.0, 1.0).solve(&mut soa, DT);

    assert_eq!(soa.lambda[0], 0.0, "separating bodies need no impulse");
    let (va, _) = soa.body_velocity(0).unwrap_or_default();
    let (vb, _) = soa.body_velocity(1).unwrap_or_default();
    assert_vec_close(va, [0.0; 3], 1e-15, "body A is left alone");
    assert_vec_close(vb, [0.0, 1.0, 0.0], 1e-15, "body B is left alone");
}

#[test]
fn net_applied_impulse_always_equals_the_accumulated_lambda() {
    for warm in [0.0, 0.5, 2.0, 7.5] {
        for approach in [-3.0, -0.5, 0.0, 0.8] {
            let row = SoaConstraintRow {
                warm_lambda: warm,
                ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.5, 2.0)
            };
            let mut soa = build(
                &[row],
                &[([0.0; 3], [0.0; 3]), ([0.0, approach, 0.0], [0.0; 3])],
            );
            BatchPgsSolver::new(4, 0.0, 1.0).solve(&mut soa, DT);
            let lambda = soa.lambda[0];
            assert!(
                lambda >= 0.0,
                "pulling impulse {lambda} (warm {warm}, approach {approach})"
            );
            let (va, _) = soa.body_velocity(0).unwrap_or_default();
            let (vb, _) = soa.body_velocity(1).unwrap_or_default();
            // The velocities changed by exactly the accumulated impulse.
            assert_close(va[1], -0.5 * lambda, 1e-12, "A");
            assert_close(vb[1] - approach, 2.0 * lambda, 1e-12, "B");
            // The relative velocity never ends up approaching.
            assert!(vb[1] - va[1] >= -1e-12);
        }
    }
}

#[test]
fn exact_warm_start_converges_without_changing_anything() {
    // B falls at 2 m/s onto static A; the exact impulse is 2 N·s (m = 1).
    let row = SoaConstraintRow {
        warm_lambda: 2.0,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 1.0)
    };
    let mut soa = build(
        &[row],
        &[([0.0; 3], [0.0; 3]), ([0.0, -2.0, 0.0], [0.0; 3])],
    );
    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    assert!(stats.converged);
    assert_eq!(stats.iterations_used, 1);
    assert_close(soa.lambda[0], 2.0, 1e-15, "lambda");
    assert_close(soa.body_vy[1], 0.0, 1e-15, "velocity");
}

#[test]
fn warm_lambda_is_projected_onto_the_bounds() {
    let negative = SoaConstraintRow {
        warm_lambda: -3.0,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0)
    };
    let mut soa = SoaConstraints::new();
    soa.push(negative);
    assert_eq!(
        soa.lambda[0], 0.0,
        "a contact cannot carry a pulling warm start"
    );
    let bounded = SoaConstraintRow {
        warm_lambda: 9.0,
        hi: 4.0,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0)
    };
    soa.push(bounded);
    assert_eq!(soa.lambda[1], 4.0);
}

#[test]
fn builder_cold_starts_and_folds_restitution_into_the_bias() {
    let dt = DT;
    let mut builder = SoaContactBuilder::new();
    let contact = SoaContact {
        restitution: 0.5,
        approach_vel: -2.0,
        ..SoaContact::centered(0, 1, unit_y(), 0.045 + 0.005, 0.0, 1.0)
    };
    builder.add_contact(contact, dt, 0.2);
    let soa = builder.build();
    assert_eq!(
        soa.lambda[0], 0.0,
        "the builder does not invent a warm start"
    );
    // slop 0.005: excess 0.045 -> -(0.2/dt)·0.045 ; restitution 0.5·(-2).
    let expected = -(0.2 / dt) * 0.045 - 1.0;
    assert_close(soa.bias[0], expected, 1e-12, "bias");
    assert_eq!(soa.lambda_lo[0], 0.0);
    assert!(soa.lambda_hi[0].is_infinite());
}

#[test]
fn builder_restitution_applies_only_to_approaching_contacts() {
    let mut builder = SoaContactBuilder::new();
    for approach_vel in [0.7, 0.0, -0.4] {
        let contact = SoaContact {
            restitution: 1.0,
            approach_vel,
            ..SoaContact::centered(0, 1, unit_y(), 0.0, 0.0, 1.0)
        };
        builder.add_contact(contact, DT, 0.2);
    }
    let soa = builder.build();
    assert_eq!(soa.bias[0], 0.0, "separating contacts do not bounce");
    assert_eq!(soa.bias[1], 0.0);
    assert_close(soa.bias[2], -0.4, 1e-15, "approaching contact");

    let mut gated = SoaContactBuilder::new().with_restitution_threshold(1.0);
    let contact = SoaContact {
        restitution: 1.0,
        approach_vel: -0.4,
        ..SoaContact::centered(0, 1, unit_y(), 0.0, 0.0, 1.0)
    };
    gated.add_contact(contact, DT, 0.2);
    assert_eq!(gated.build().bias[0], 0.0, "slow impacts are inelastic");
}

#[test]
fn restitution_bounces_the_body_back() {
    let mut builder = SoaContactBuilder::new();
    let contact = SoaContact {
        restitution: 0.5,
        approach_vel: -2.0,
        ..SoaContact::centered(0, 1, unit_y(), 0.0, 0.0, 1.0)
    };
    builder.add_contact(contact, DT, 0.2);
    let mut soa = builder.build();
    soa.set_body_velocity(1, [0.0, -2.0, 0.0], [0.0; 3]);
    BatchPgsSolver::default().solve(&mut soa, DT);
    // v' = -e v = 1 m/s, from an impulse of (1 + e)·m·|v| = 3.
    assert_close(soa.body_vy[1], 1.0, 1e-12, "rebound velocity");
    assert_close(soa.lambda[0], 3.0, 1e-12, "impulse");
}

#[test]
fn builder_position_correction_pushes_penetrating_bodies_apart() {
    let mut builder = SoaContactBuilder::new().with_slop(0.0);
    builder.add_contact(
        SoaContact::centered(0, 1, unit_y(), 0.05, 0.0, 1.0),
        DT,
        0.2,
    );
    let mut soa = builder.build();
    BatchPgsSolver::default().solve(&mut soa, DT);
    // The body leaves at (β/dt)·penetration = 0.2·60·0.05 = 0.6 m/s.
    assert_close(soa.body_vy[1], 0.6, 1e-12, "separation velocity");
}

#[test]
fn builder_ignores_a_non_positive_time_step() {
    let mut builder = SoaContactBuilder::new();
    builder.add_contact(
        SoaContact::centered(0, 1, unit_y(), 0.2, 0.0, 1.0),
        0.0,
        0.2,
    );
    let soa = builder.build();
    assert_eq!(soa.bias[0], 0.0);
    assert!(soa.bias[0].is_finite());
}

// ── Solver mechanics ──────────────────────────────────────────────────────

#[test]
fn solver_is_stable_for_inverted_and_nan_bounds() {
    let mut soa = SoaConstraints::new();
    soa.push(SoaConstraintRow {
        lo: 1.0,
        hi: -1.0,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0)
    });
    soa.push(SoaConstraintRow {
        lo: f64::NAN,
        hi: f64::NAN,
        ..SoaConstraintRow::linear(unit_y(), 0, 1, 1.0, 1.0)
    });
    soa.set_body_velocity(1, [0.0, -1.0, 0.0], [0.0; 3]);
    // No panic, whatever the bounds are.
    BatchPgsSolver::default().solve(&mut soa, DT);
}

#[test]
fn omega_scales_the_impulse_update() {
    let make = || {
        build(
            &[SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 1.0)],
            &[([0.0; 3], [0.0; 3]), ([0.0, -2.0, 0.0], [0.0; 3])],
        )
    };
    let mut plain = make();
    BatchPgsSolver::new(1, 0.0, 1.0).solve(&mut plain, DT);
    assert_close(plain.lambda[0], 2.0, 1e-15, "omega = 1");

    let mut under = make();
    BatchPgsSolver::new(1, 0.0, 0.5).solve(&mut under, DT);
    assert_close(
        under.lambda[0],
        1.0,
        1e-15,
        "omega = 0.5 takes half the step",
    );

    let mut over = make();
    BatchPgsSolver::new(1, 0.0, 1.5).solve(&mut over, DT);
    assert_close(over.lambda[0], 3.0, 1e-15, "omega = 1.5 overshoots");
    // ... and the unilateral row cannot pull the overshoot back in one
    // sweep, but it does so over the following iterations.
    let mut settled = make();
    BatchPgsSolver::new(60, 1e-12, 1.5).solve(&mut settled, DT);
    assert_close(settled.lambda[0], 2.0, 1e-9, "converged SOR impulse");
}

#[test]
fn relative_velocity_is_recorded_per_row() {
    let mut soa = build(
        &[SoaConstraintRow {
            r_b: [1.0, 0.0, 0.0],
            inv_inertia_b: diag(1.0, 1.0, 1.0),
            ..SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 1.0)
        }],
        &[([0.0; 3], [0.0; 3]), ([0.0, -1.0, 0.0], [0.0, 0.0, 2.0])],
    );
    BatchPgsSolver::new(1, 0.0, 1.0).solve(&mut soa, DT);
    // rv = n·v_b + (r_b × n)·ω_b = -1 + (0,0,1)·(0,0,2) = 1: separating at
    // the contact point, the row records it and applies nothing.
    assert_close(soa.rel_vel[0], 1.0, 1e-15, "recorded relative velocity");
    assert_eq!(soa.lambda[0], 0.0);
}

#[test]
fn inconsistent_or_undersized_sets_never_panic() {
    // Edited by hand: one per-row array is longer than the others.
    let mut soa = build(
        &[SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 1.0)],
        &[([0.0; 3], [0.0; 3]), ([0.0, -1.0, 0.0], [0.0; 3])],
    );
    soa.nx.push(1.0);
    assert!(!soa.is_consistent());
    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    assert_eq!(stats.iterations_used, 0);
    assert!(!stats.converged);
    assert_eq!(soa.body_vy[1], -1.0, "nothing is touched");

    // Body arrays that were shrunk below the used indices are regrown.
    let mut soa = build(
        &[SoaConstraintRow::linear(unit_y(), 0, 3, 0.0, 1.0)],
        &[([0.0; 3], [0.0; 3])],
    );
    soa.resize_bodies(1);
    assert!(soa.is_consistent());
    soa.body_a_idx[0] = 0;
    soa.body_b_idx[0] = 3;
    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    assert!(stats.converged);
    assert_eq!(soa.body_count(), 4);
}

#[test]
fn push_grows_the_body_arrays() {
    let mut soa = SoaConstraints::new();
    soa.push(SoaConstraintRow::linear(unit_y(), 3, 7, 1.0, 1.0));
    assert_eq!(soa.body_count(), 8);
    assert_eq!(soa.body_vx.len(), 8);
    assert_eq!(soa.body_wz.len(), 8);
    assert!(soa.body_velocity(8).is_none());
    soa.clear();
    assert!(soa.is_empty());
    assert_eq!(soa.body_count(), 8, "clear keeps the body velocities");
}

#[test]
fn rows_without_effective_mass_do_nothing() {
    // Two static bodies: K = 0, the row must neither divide by zero nor
    // move anything.
    let mut soa = build(
        &[SoaConstraintRow::linear(unit_y(), 0, 1, 0.0, 0.0)],
        &[([0.0, 1.0, 0.0], [0.0; 3]), ([0.0, -1.0, 0.0], [0.0; 3])],
    );
    assert_eq!(soa.eff_mass[0], 0.0);
    let stats = BatchPgsSolver::default().solve(&mut soa, DT);
    assert!(stats.converged);
    assert_eq!(soa.lambda[0], 0.0);
    assert_eq!(soa.body_vy[0], 1.0);
    assert_eq!(soa.body_vy[1], -1.0);
}
