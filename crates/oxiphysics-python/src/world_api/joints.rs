// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint solver of the rigid-body world: distance, point-to-point (ball-socket),
//! hinge and slider (prismatic) joints between two bodies.
//!
//! Every joint is a set of scalar rows `J` over the linear and angular
//! velocities `(v_A, ω_A, v_B, ω_B)` of its two bodies. Rows that belong
//! together form a block of up to three rows, which is solved jointly with the
//! block's effective-mass matrix `K = J M⁻¹ Jᵀ`:
//!
//! * the position pass removes the position error `C` of the rows by moving and
//!   rotating the bodies, `λ = -s K⁻¹ C`, without touching the velocities,
//! * the velocity pass removes the relative velocity `J v` along the rows by
//!   applying impulses, `λ = -s K⁻¹ J v`.
//!
//! Both passes are repeated for the configured number of solver iterations; the
//! relaxation factor `s` in `[0, 1]` is the stiffness of the joint (`1` removes
//! the whole violation of a row in one iteration). The bodies are rigid with a
//! scalar moment of inertia, anchors are body-local points, so a joint whose
//! anchor is off the centre of mass exerts a torque.
//!
//! The rows of the joints, with `r_X` the world offset from the centre of mass
//! of body `X` to its anchor, `d = (p_B + r_B) - (p_A + r_A)` the separation of
//! the anchors and `t` a direction:
//!
//! * a linear row along `t` is `(v_A, ω_A, v_B, ω_B) ↦ t·(v_B + ω_B×r_B -
//!   v_A - ω_A×r_A)`, with the angular part of `A` extended to
//!   `-((r_A + d) × t)` when `t` is fixed in the frame of `A` (slider),
//! * an angular row along `t` is `t·(ω_B - ω_A)`.
//!
//! Distance: one linear row along `d` with `C = |d| - L`. Point-to-point: three
//! linear rows along the world axes with `C = d`. Hinge: point-to-point plus two
//! angular rows perpendicular to the hinge axis of `A` with `C = t·(a_A × a_B)`
//! for the world axes `a_A`, `a_B` of the hinge on the two bodies. Slider: two
//! linear rows perpendicular to the slider axis with `C = t·d`, plus three
//! angular rows that keep the relative orientation of its creation.

use super::math_helpers::{quat_conjugate, quat_mul, quat_normalize, quat_rotate_vec};

/// Below this distance two anchors have no direction.
const MIN_DIRECTION_LENGTH: f64 = 1.0e-12;

/// A block whose effective-mass matrix has a pivot below this fraction of its
/// largest diagonal entry is skipped (both bodies immovable, or no freedom).
const SINGULAR_FRACTION: f64 = 1.0e-12;

/// State of a body as the joint solver sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct JointBody {
    /// Centre of mass position.
    pub position: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// Linear velocity.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
    /// Inverse mass (`0` for a body that cannot move).
    pub inv_mass: f64,
    /// Inverse of the scalar moment of inertia (`0` for a body that cannot turn).
    pub inv_inertia: f64,
}

impl JointBody {
    /// Add the impulse `(linear, angular)` to the velocities.
    fn apply_impulse(&mut self, linear: [f64; 3], angular: [f64; 3]) {
        for k in 0..3 {
            self.velocity[k] += self.inv_mass * linear[k];
            self.angular_velocity[k] += self.inv_inertia * angular[k];
        }
    }

    /// Move the body by `inv_mass * linear` and turn it by the rotation vector
    /// `inv_inertia * angular`.
    fn apply_displacement(&mut self, linear: [f64; 3], angular: [f64; 3]) {
        for (position, push) in self.position.iter_mut().zip(linear) {
            *position += self.inv_mass * push;
        }
        let rotation = scale(angular, self.inv_inertia);
        if dot(rotation, rotation) > 0.0 {
            self.orientation = rotate_orientation(self.orientation, rotation);
        }
    }
}

/// Kind of a joint with the data fixed when the joint is created.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum JointKind {
    /// Keep the anchors a fixed distance apart.
    Distance {
        /// Required distance of the anchors (m).
        length: f64,
    },
    /// Make the anchors coincide.
    PointToPoint,
    /// Make the anchors coincide and allow rotation only about the hinge axis.
    Hinge,
    /// Allow the anchors to separate only along the slider axis and keep the
    /// relative orientation.
    Slider {
        /// Orientation of body `B` relative to body `A` to be kept,
        /// `conj(q_A) * q_B` at creation, as `[x, y, z, w]`.
        reference_rotation: [f64; 4],
    },
}

/// A joint between two bodies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Joint {
    /// Kind of the joint.
    pub kind: JointKind,
    /// Anchor on body `A` in its local frame.
    pub anchor_a: [f64; 3],
    /// Anchor on body `B` in its local frame.
    pub anchor_b: [f64; 3],
    /// Unit axis (hinge, slider) in the local frame of body `A`.
    pub axis_a: [f64; 3],
    /// Unit axis (hinge) in the local frame of body `B`.
    pub axis_b: [f64; 3],
    /// Relaxation factor in `[0, 1]`: the fraction of a violation removed per
    /// solver iteration.
    pub stiffness: f64,
}

/// One scalar row of a joint: the Jacobian entries for `(v_A, ω_A, v_B, ω_B)`.
#[derive(Debug, Clone, Copy, Default)]
struct Row {
    lin_a: [f64; 3],
    ang_a: [f64; 3],
    lin_b: [f64; 3],
    ang_b: [f64; 3],
}

impl Row {
    /// `J v`: the velocity of the constrained quantity.
    fn velocity(&self, a: &JointBody, b: &JointBody) -> f64 {
        dot(self.lin_a, a.velocity)
            + dot(self.ang_a, a.angular_velocity)
            + dot(self.lin_b, b.velocity)
            + dot(self.ang_b, b.angular_velocity)
    }

    /// Entry of the effective-mass matrix, `J_i M⁻¹ J_jᵀ`.
    fn coupling(&self, other: &Row, a: &JointBody, b: &JointBody) -> f64 {
        a.inv_mass * dot(self.lin_a, other.lin_a)
            + a.inv_inertia * dot(self.ang_a, other.ang_a)
            + b.inv_mass * dot(self.lin_b, other.lin_b)
            + b.inv_inertia * dot(self.ang_b, other.ang_b)
    }
}

/// Up to three rows that are solved together, with their position errors.
#[derive(Debug, Clone, Copy)]
struct RowBlock {
    len: usize,
    rows: [Row; 3],
    error: [f64; 3],
}

impl RowBlock {
    fn new(len: usize) -> Self {
        Self {
            len,
            rows: [Row::default(); 3],
            error: [0.0; 3],
        }
    }

    /// Multipliers `λ` with `K λ = rhs`, or `None` if the block has no freedom.
    fn multipliers(&self, a: &JointBody, b: &JointBody, rhs: [f64; 3]) -> Option<[f64; 3]> {
        let mut matrix = [[0.0; 3]; 3];
        for (i, row_i) in self.rows.iter().enumerate().take(self.len) {
            for (j, row_j) in self.rows.iter().enumerate().take(self.len) {
                matrix[i][j] = row_i.coupling(row_j, a, b);
            }
        }
        // Rows beyond the block decouple: identity with a zero right-hand side.
        for (i, diagonal) in matrix.iter_mut().enumerate().skip(self.len) {
            diagonal[i] = 1.0;
        }
        solve3(matrix, rhs)
    }

    /// Sum of `λ_i * J_i` over the rows, for the four Jacobian blocks.
    fn combine(&self, lambda: [f64; 3]) -> Row {
        let mut sum = Row::default();
        for (row, &weight) in self.rows.iter().zip(&lambda).take(self.len) {
            for k in 0..3 {
                sum.lin_a[k] += weight * row.lin_a[k];
                sum.ang_a[k] += weight * row.ang_a[k];
                sum.lin_b[k] += weight * row.lin_b[k];
                sum.ang_b[k] += weight * row.ang_b[k];
            }
        }
        sum
    }

    /// Remove the fraction `relax` of the position error.
    fn solve_position(&self, a: &mut JointBody, b: &mut JointBody, relax: f64) {
        let rhs = scale(self.error, -relax);
        if let Some(lambda) = self.multipliers(a, b, rhs) {
            let push = self.combine(lambda);
            a.apply_displacement(push.lin_a, push.ang_a);
            b.apply_displacement(push.lin_b, push.ang_b);
        }
    }

    /// Remove the fraction `relax` of the relative velocity along the rows.
    fn solve_velocity(&self, a: &mut JointBody, b: &mut JointBody, relax: f64) {
        let mut rhs = [0.0; 3];
        for (slot, row) in rhs.iter_mut().zip(&self.rows).take(self.len) {
            *slot = -relax * row.velocity(a, b);
        }
        if let Some(lambda) = self.multipliers(a, b, rhs) {
            let impulse = self.combine(lambda);
            a.apply_impulse(impulse.lin_a, impulse.ang_a);
            b.apply_impulse(impulse.lin_b, impulse.ang_b);
        }
    }
}

/// World-space geometry of a joint for the current poses of its bodies.
struct Frame {
    /// Offset from the centre of mass of `A` to its anchor.
    r_a: [f64; 3],
    /// Offset from the centre of mass of `B` to its anchor.
    r_b: [f64; 3],
    /// Anchor of `B` minus anchor of `A`.
    separation: [f64; 3],
}

impl Joint {
    /// Number of blocks of the joint.
    fn block_count(&self) -> usize {
        match self.kind {
            JointKind::Distance { .. } | JointKind::PointToPoint => 1,
            JointKind::Hinge | JointKind::Slider { .. } => 2,
        }
    }

    fn frame(&self, a: &JointBody, b: &JointBody) -> Frame {
        let r_a = quat_rotate_vec(a.orientation, self.anchor_a);
        let r_b = quat_rotate_vec(b.orientation, self.anchor_b);
        let separation = sub(add(b.position, r_b), add(a.position, r_a));
        Frame {
            r_a,
            r_b,
            separation,
        }
    }

    /// The `index`-th block for the current poses, `None` if it is undefined
    /// (a distance joint whose anchors coincide).
    fn block(&self, index: usize, a: &JointBody, b: &JointBody) -> Option<RowBlock> {
        let frame = self.frame(a, b);
        match (self.kind, index) {
            (JointKind::Distance { length }, 0) => {
                let distance = norm(frame.separation);
                if distance < MIN_DIRECTION_LENGTH {
                    return None;
                }
                let direction = scale(frame.separation, 1.0 / distance);
                let mut block = RowBlock::new(1);
                block.rows[0] = linear_row(direction, &frame, [0.0; 3]);
                block.error[0] = distance - length;
                Some(block)
            }
            (JointKind::PointToPoint | JointKind::Hinge, 0) => Some(point_block(&frame)),
            (JointKind::Hinge, 1) => {
                let axis_a = quat_rotate_vec(a.orientation, self.axis_a);
                let axis_b = quat_rotate_vec(b.orientation, self.axis_b);
                let misalignment = cross(axis_a, axis_b);
                let (first, second) = perpendicular_pair(axis_a);
                let mut block = RowBlock::new(2);
                for (k, tangent) in [first, second].into_iter().enumerate() {
                    block.rows[k] = angular_row(tangent);
                    block.error[k] = dot(tangent, misalignment);
                }
                Some(block)
            }
            (JointKind::Slider { .. }, 0) => {
                let axis = quat_rotate_vec(a.orientation, self.axis_a);
                let (first, second) = perpendicular_pair(axis);
                let mut block = RowBlock::new(2);
                for (k, tangent) in [first, second].into_iter().enumerate() {
                    block.rows[k] = linear_row(tangent, &frame, frame.separation);
                    block.error[k] = dot(tangent, frame.separation);
                }
                Some(block)
            }
            (JointKind::Slider { reference_rotation }, 1) => {
                let wanted = quat_mul(a.orientation, reference_rotation);
                let mut offset = quat_mul(b.orientation, quat_conjugate(wanted));
                if offset[3] < 0.0 {
                    offset = [-offset[0], -offset[1], -offset[2], -offset[3]];
                }
                let mut block = RowBlock::new(3);
                for (k, axis) in WORLD_AXES.into_iter().enumerate() {
                    block.rows[k] = angular_row(axis);
                    block.error[k] = 2.0 * offset[k];
                }
                Some(block)
            }
            _ => None,
        }
    }

    /// Remove the fraction `stiffness` of the position error of every block.
    pub(super) fn solve_positions(&self, a: &mut JointBody, b: &mut JointBody) {
        for index in 0..self.block_count() {
            if let Some(block) = self.block(index, a, b) {
                block.solve_position(a, b, self.stiffness);
            }
        }
    }

    /// Remove the fraction `stiffness` of the relative velocity of every block.
    pub(super) fn solve_velocities(&self, a: &mut JointBody, b: &mut JointBody) {
        for index in 0..self.block_count() {
            if let Some(block) = self.block(index, a, b) {
                block.solve_velocity(a, b, self.stiffness);
            }
        }
    }

    /// Residual of the joint: the norm of the position error over all blocks.
    #[cfg(test)]
    pub(super) fn position_error(&self, a: &JointBody, b: &JointBody) -> f64 {
        (0..self.block_count())
            .filter_map(|index| self.block(index, a, b))
            .map(|block| {
                block
                    .error
                    .iter()
                    .take(block.len)
                    .map(|e| e * e)
                    .sum::<f64>()
            })
            .sum::<f64>()
            .sqrt()
    }
}

/// Solve the joints for the bodies they connect.
///
/// `joints` holds `(index of A, index of B, joint)` into `bodies`. The position
/// pass runs `iterations` sweeps over all joints, then the velocity pass runs
/// `iterations` sweeps. Joints that connect a body with itself are ignored.
pub(super) fn solve_joints(bodies: &mut [JointBody], joints: &[(usize, usize, Joint)], iterations: u32) {
    for pass in [Pass::Positions, Pass::Velocities] {
        for _ in 0..iterations {
            for &(index_a, index_b, joint) in joints {
                let Some((a, b)) = pair_mut(bodies, index_a, index_b) else {
                    continue;
                };
                match pass {
                    Pass::Positions => joint.solve_positions(a, b),
                    Pass::Velocities => joint.solve_velocities(a, b),
                }
            }
        }
    }
}

/// Which quantity a solver sweep corrects.
#[derive(Clone, Copy)]
enum Pass {
    Positions,
    Velocities,
}

/// Two distinct elements of a slice, mutably.
fn pair_mut(
    bodies: &mut [JointBody],
    first: usize,
    second: usize,
) -> Option<(&mut JointBody, &mut JointBody)> {
    if first == second || first >= bodies.len() || second >= bodies.len() {
        return None;
    }
    if first < second {
        let (head, tail) = bodies.split_at_mut(second);
        Some((&mut head[first], &mut tail[0]))
    } else {
        let (head, tail) = bodies.split_at_mut(first);
        Some((&mut tail[0], &mut head[second]))
    }
}

/// The world axes.
const WORLD_AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Linear row along `direction`; `offset_a` is added to `r_A` in the angular
/// part of body `A` (the anchor separation for a direction fixed in `A`).
fn linear_row(direction: [f64; 3], frame: &Frame, offset_a: [f64; 3]) -> Row {
    Row {
        lin_a: scale(direction, -1.0),
        ang_a: scale(cross(add(frame.r_a, offset_a), direction), -1.0),
        lin_b: direction,
        ang_b: cross(frame.r_b, direction),
    }
}

/// Angular row along `direction`.
fn angular_row(direction: [f64; 3]) -> Row {
    Row {
        lin_a: [0.0; 3],
        ang_a: scale(direction, -1.0),
        lin_b: [0.0; 3],
        ang_b: direction,
    }
}

/// The three linear rows of a point-to-point connection.
fn point_block(frame: &Frame) -> RowBlock {
    let mut block = RowBlock::new(3);
    for (k, axis) in WORLD_AXES.into_iter().enumerate() {
        block.rows[k] = linear_row(axis, frame, [0.0; 3]);
        block.error[k] = frame.separation[k];
    }
    block
}

/// Two unit vectors that complete the unit vector `axis` to an orthonormal basis.
fn perpendicular_pair(axis: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Cross with the world axis that is least aligned with `axis`.
    let helper = if axis[0].abs() <= axis[1].abs() && axis[0].abs() <= axis[2].abs() {
        WORLD_AXES[0]
    } else if axis[1].abs() <= axis[2].abs() {
        WORLD_AXES[1]
    } else {
        WORLD_AXES[2]
    };
    let first = normalized(cross(axis, helper));
    let second = normalized(cross(axis, first));
    (first, second)
}

/// Turn `orientation` by the world-frame rotation vector `rotation` (small angles).
fn rotate_orientation(orientation: [f64; 4], rotation: [f64; 3]) -> [f64; 4] {
    let half = scale(rotation, 0.5);
    let turned = quat_mul([half[0], half[1], half[2], 0.0], orientation);
    quat_normalize([
        orientation[0] + turned[0],
        orientation[1] + turned[1],
        orientation[2] + turned[2],
        orientation[3] + turned[3],
    ])
}

/// Solve the `3 x 3` system `m x = rhs` by Gaussian elimination with partial
/// pivoting; `None` if the matrix is singular.
fn solve3(mut m: [[f64; 3]; 3], mut rhs: [f64; 3]) -> Option<[f64; 3]> {
    let largest_diagonal = (0..3).map(|i| m[i][i].abs()).fold(0.0_f64, f64::max);
    if largest_diagonal.is_nan() || largest_diagonal <= 0.0 {
        return None;
    }
    let tolerance = SINGULAR_FRACTION * largest_diagonal;
    for col in 0..3 {
        let mut pivot_row = col;
        for row in (col + 1)..3 {
            if m[row][col].abs() > m[pivot_row][col].abs() {
                pivot_row = row;
            }
        }
        let pivot = m[pivot_row][col];
        if pivot.is_nan() || pivot.abs() <= tolerance {
            return None;
        }
        m.swap(col, pivot_row);
        rhs.swap(col, pivot_row);
        for row in (col + 1)..3 {
            let factor = m[row][col] / m[col][col];
            for k in col..3 {
                m[row][k] -= factor * m[col][k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        let tail: f64 = ((row + 1)..3).map(|k| m[row][k] * x[k]).sum();
        x[row] = (rhs[row] - tail) / m[row][row];
    }
    Some(x)
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], factor: f64) -> [f64; 3] {
    [a[0] * factor, a[1] * factor, a[2] * factor]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn normalized(a: [f64; 3]) -> [f64; 3] {
    let length = norm(a);
    if length > 0.0 {
        scale(a, 1.0 / length)
    } else {
        a
    }
}
