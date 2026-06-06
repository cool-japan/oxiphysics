// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Inverse kinematics — FABRIK and 2-bone analytic IK with joint cone limits.
//!
//! ## Types
//!
//! - `IkJoint`  — single joint with local offset and cone/twist limits
//! - `IkChain`  — kinematic chain: root + N joints + segment lengths
//! - `IkSolver` — solver variant (FABRIK or TwoBone analytic)
//! - `SolveReport` — result summary (iterations, convergence, residual)
//! - `IkError`  — error type (unused by `solve`, kept for the type inventory)
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::ik::{IkChain, IkSolver};
//!
//! // Three 1-metre segments hanging downward from origin
//! let mut chain = IkChain::new([0.0, 0.0, 0.0], &[
//!     (1.0, std::f64::consts::PI, 0.0),
//!     (1.0, std::f64::consts::PI, 0.0),
//!     (1.0, std::f64::consts::PI, 0.0),
//! ]);
//!
//! let solver = IkSolver::Fabrik { max_iterations: 30, tolerance: 1e-6 };
//! let report = solver.solve(&mut chain, [0.0, 3.0, 0.0]);
//! assert!(report.converged);
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Helper math (no external crates)
// ---------------------------------------------------------------------------

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn len_sq(a: [f64; 3]) -> f64 {
    dot(a, a)
}

#[inline]
fn len(a: [f64; 3]) -> f64 {
    len_sq(a).sqrt()
}

#[inline]
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    len(sub(b, a))
}

/// Returns `None` if `a` is degenerate (near-zero length).
#[inline]
fn normalize(a: [f64; 3]) -> Option<[f64; 3]> {
    let l = len(a);
    if l < 1e-12 {
        None
    } else {
        Some(scale(a, 1.0 / l))
    }
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Any unit vector not parallel to `v` (for building a perpendicular basis).
fn any_perpendicular(v: [f64; 3]) -> [f64; 3] {
    // Pick the world axis least aligned with v.
    let ax = v[0].abs();
    let ay = v[1].abs();
    let az = v[2].abs();
    let candidate = if ax <= ay && ax <= az {
        [1.0, 0.0, 0.0]
    } else if ay <= ax && ay <= az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize(cross(v, candidate)).unwrap_or([0.0, 1.0, 0.0])
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A single joint in the IK chain.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct IkJoint {
    /// Offset from this joint to the next joint (world-space after solving).
    pub local_offset: [f64; 3],
    /// Half-angle of the allowed cone (radians). 0 = fully locked, PI = unlimited.
    pub cone_limit_rad: f64,
    /// Allowed twist around the bone axis (radians). Ignored by FABRIK.
    pub twist_limit_rad: f64,
}

/// The kinematic chain (root anchor + N joints).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct IkChain {
    /// World-space position of the root (first joint base).
    pub root_world: [f64; 3],
    /// Joint definitions; `len()` equals the number of segments.
    pub joints: Vec<IkJoint>,
    /// Length of each segment.
    pub segment_lengths: Vec<f64>,
    /// Optional pole target (affects bend-plane for TwoBone; unused by FABRIK).
    pub pole_target: Option<[f64; 3]>,
}

impl IkChain {
    /// Build a chain from `root_world` and a list of
    /// `(segment_length, cone_limit_rad, twist_limit_rad)` tuples.
    ///
    /// `local_offset` is initialised as `seg_len * [0, -1, 0]` (default downward).
    pub fn new(root_world: [f64; 3], segments: &[(f64, f64, f64)]) -> Self {
        let joints: Vec<IkJoint> = segments
            .iter()
            .map(|&(seg_len, cone_limit_rad, twist_limit_rad)| IkJoint {
                local_offset: [0.0, -seg_len, 0.0],
                cone_limit_rad,
                twist_limit_rad,
            })
            .collect();
        let segment_lengths: Vec<f64> = segments.iter().map(|&(l, _, _)| l).collect();
        Self {
            root_world,
            joints,
            segment_lengths,
            pole_target: None,
        }
    }

    /// Total reachable length of the chain.
    pub fn total_reach(&self) -> f64 {
        self.segment_lengths.iter().sum()
    }
}

/// Result of a solve call.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SolveReport {
    pub iterations: usize,
    pub converged: bool,
    /// Distance from end-effector to target after solving.
    pub residual: f64,
}

/// Error type for IK operations (kept for the type inventory; not returned by `solve`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum IkError {
    TwoBoneRequiresTwoSegments,
    ChainTooShort,
}

/// The IK solver variant.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum IkSolver {
    Fabrik {
        max_iterations: usize,
        tolerance: f64,
    },
    /// Closed-form analytic IK for exactly 2 segments.
    TwoBone,
}

// ---------------------------------------------------------------------------
// Cone-limit helper
// ---------------------------------------------------------------------------

/// Clamp `child` so that the angle between the bone direction
/// (`joint → child`) and the reference direction (`parent → joint`)
/// does not exceed `cone_limit_rad`.
///
/// `parent` is the position *before* `joint`; skipped when there is no parent
/// (call with `parent = None` to skip).  Also skipped when the cone is
/// unlimited (`cone_limit_rad >= PI - 1e-6`).
fn clamp_cone(
    parent: Option<[f64; 3]>,
    joint: [f64; 3],
    child: &mut [f64; 3],
    cone_limit_rad: f64,
    seg_len: f64,
) {
    use std::f64::consts::PI;

    if cone_limit_rad >= PI - 1e-6 {
        return; // unlimited — nothing to clamp
    }

    let Some(parent_pos) = parent else { return };

    // Reference direction: parent → joint
    let ref_dir = match normalize(sub(joint, parent_pos)) {
        Some(d) => d,
        None => return, // degenerate, skip
    };

    // Bone direction: joint → child
    let bone_vec = sub(*child, joint);
    let bone_dir = match normalize(bone_vec) {
        Some(d) => d,
        None => return, // degenerate, skip
    };

    let cos_a = dot(ref_dir, bone_dir).clamp(-1.0, 1.0);
    let angle = cos_a.acos();

    if angle <= cone_limit_rad {
        return; // already within cone
    }

    // Rotate bone_dir back toward ref_dir so angle == cone_limit_rad.
    // axis = normalize(cross(ref_dir, bone_dir))
    let axis_raw = cross(ref_dir, bone_dir);
    let axis = match normalize(axis_raw) {
        Some(a) => a,
        None => return, // parallel or anti-parallel
    };

    // Rodrigues' rotation of ref_dir by cone_limit_rad around axis:
    let cos_c = cone_limit_rad.cos();
    let sin_c = cone_limit_rad.sin();
    let clamped_dir = add(
        add(scale(ref_dir, cos_c), scale(cross(axis, ref_dir), sin_c)),
        scale(axis, dot(axis, ref_dir) * (1.0 - cos_c)),
    );

    // Normalise and place child at the correct distance from joint.
    let clamped_dir = normalize(clamped_dir).unwrap_or(ref_dir);
    *child = add(joint, scale(clamped_dir, seg_len));
}

// ---------------------------------------------------------------------------
// Solver implementation
// ---------------------------------------------------------------------------

impl IkSolver {
    /// Solve the chain toward `target`, writing joint positions back into
    /// `chain.joints[i].local_offset`.
    pub fn solve(&self, chain: &mut IkChain, target: [f64; 3]) -> SolveReport {
        match self {
            IkSolver::Fabrik {
                max_iterations,
                tolerance,
            } => solve_fabrik(chain, target, *max_iterations, *tolerance),
            IkSolver::TwoBone => solve_two_bone(chain, target),
        }
    }
}

// ---------------------------------------------------------------------------
// FABRIK solver
// ---------------------------------------------------------------------------

fn solve_fabrik(
    chain: &mut IkChain,
    target: [f64; 3],
    max_iterations: usize,
    tolerance: f64,
) -> SolveReport {
    let n = chain.joints.len(); // number of segments
    if n == 0 {
        return SolveReport {
            iterations: 0,
            converged: false,
            residual: dist(chain.root_world, target),
        };
    }

    // Build initial positions from chain data.
    // positions[0] = root, positions[n] = tip.
    let mut positions: Vec<[f64; 3]> = Vec::with_capacity(n + 1);
    positions.push(chain.root_world);
    for i in 0..n {
        let prev = positions[i];
        positions.push(add(prev, chain.joints[i].local_offset));
    }

    let total_reach: f64 = chain.segment_lengths.iter().sum();
    let root_to_target = dist(chain.root_world, target);

    // Out of reach: fully extend toward target.
    if root_to_target > total_reach {
        let dir = normalize(sub(target, chain.root_world)).unwrap_or([1.0, 0.0, 0.0]);
        positions[0] = chain.root_world;
        for i in 0..n {
            positions[i + 1] = add(positions[i], scale(dir, chain.segment_lengths[i]));
        }
        // Write back
        for i in 0..n {
            chain.joints[i].local_offset = sub(positions[i + 1], positions[i]);
        }
        let residual = dist(positions[n], target);
        return SolveReport {
            iterations: 0,
            converged: false,
            residual,
        };
    }

    let mut converged = false;
    let mut iters = 0;

    for _ in 0..max_iterations {
        iters += 1;

        // --- Backward pass: effector → root ---
        positions[n] = target;
        for i in (0..n).rev() {
            let dir = normalize(sub(positions[i], positions[i + 1])).unwrap_or([0.0, 1.0, 0.0]);
            positions[i] = add(positions[i + 1], scale(dir, chain.segment_lengths[i]));

            // Cone clamp: parent of positions[i] is positions[i-1] (if it exists)
            if i > 0 {
                let parent = Some(positions[i - 1]);
                clamp_cone(
                    parent,
                    positions[i],
                    &mut positions[i + 1],
                    chain.joints[i].cone_limit_rad,
                    chain.segment_lengths[i],
                );
            }
        }

        // --- Forward pass: root → effector ---
        positions[0] = chain.root_world;
        for i in 0..n {
            let dir = normalize(sub(positions[i + 1], positions[i])).unwrap_or([0.0, -1.0, 0.0]);
            positions[i + 1] = add(positions[i], scale(dir, chain.segment_lengths[i]));

            // Cone clamp using next point as reference for child direction
            let child_ref = if i + 2 <= n { positions[i + 2] } else { target };
            let mut child_pos = child_ref;
            clamp_cone(
                Some(positions[i]),
                positions[i + 1],
                &mut child_pos,
                chain.joints[i].cone_limit_rad,
                chain.segment_lengths[i],
            );
            // Only update next position if it exists in our slice
            if i + 2 <= n {
                positions[i + 2] = child_pos;
            }
        }

        let residual = dist(positions[n], target);
        if residual < tolerance {
            converged = true;
            break;
        }
    }

    let residual = dist(positions[n], target);

    // Write back world offsets
    for i in 0..n {
        chain.joints[i].local_offset = sub(positions[i + 1], positions[i]);
    }

    SolveReport {
        iterations: iters,
        converged,
        residual,
    }
}

// ---------------------------------------------------------------------------
// TwoBone (analytic, closed-form) solver
// ---------------------------------------------------------------------------

fn solve_two_bone(chain: &mut IkChain, target: [f64; 3]) -> SolveReport {
    // Requires exactly 2 segments.
    if chain.joints.len() != 2 || chain.segment_lengths.len() != 2 {
        return SolveReport {
            iterations: 0,
            converged: false,
            residual: f64::INFINITY,
        };
    }

    let l1 = chain.segment_lengths[0];
    let l2 = chain.segment_lengths[1];

    let root = chain.root_world;
    let d_raw = dist(root, target);
    // Clamp d so the triangle inequality holds.
    let d = d_raw.clamp((l1 - l2).abs(), l1 + l2 - 1e-9);

    // Angle at the root (α), by law of cosines for the L1-d-L2 triangle:
    //   L2² = L1² + d² - 2·L1·d·cos α
    //   cos α = (L1² + d² - L2²) / (2·L1·d)
    let cos_alpha = ((l1 * l1 + d * d - l2 * l2) / (2.0 * l1 * d)).clamp(-1.0, 1.0);
    let sin_alpha = (1.0 - cos_alpha * cos_alpha).max(0.0).sqrt();

    // Direction from root to target.
    let along = normalize(sub(target, root)).unwrap_or([1.0, 0.0, 0.0]);

    // Build a perpendicular direction (bend-plane normal direction).
    let perp = if let Some(pole) = chain.pole_target {
        // n = cross(target-root, pole-root), then perp = cross(n, along)
        let pole_vec = sub(pole, root);
        let n = cross(sub(target, root), pole_vec);
        match normalize(n) {
            Some(n_norm) => {
                let p = cross(n_norm, along);
                normalize(p).unwrap_or_else(|| any_perpendicular(along))
            }
            None => any_perpendicular(along),
        }
    } else {
        any_perpendicular(along)
    };

    // Elbow position:
    //   project along root→target by  L1·cos α
    //   project perpendicular by      L1·sin α
    let elbow_pos = add(
        add(root, scale(along, l1 * cos_alpha)),
        scale(perp, l1 * sin_alpha),
    );

    // Tip position: place at distance L2 from elbow toward target.
    let to_target_from_elbow = normalize(sub(target, elbow_pos)).unwrap_or(along);
    let tip_pos = add(elbow_pos, scale(to_target_from_elbow, l2));

    // Write back.
    chain.joints[0].local_offset = sub(elbow_pos, root);
    chain.joints[1].local_offset = sub(tip_pos, elbow_pos);
    chain.segment_lengths[0] = l1;
    chain.segment_lengths[1] = l2;

    let residual = dist(tip_pos, target);
    SolveReport {
        iterations: 1,
        converged: residual < 1e-4,
        residual,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const PI: f64 = std::f64::consts::PI;

    // Helper: world-space tip position of the chain.
    fn tip(chain: &IkChain) -> [f64; 3] {
        let mut pos = chain.root_world;
        for j in &chain.joints {
            pos = add(pos, j.local_offset);
        }
        pos
    }

    // 1. FABRIK straight chain: 3 segments × 1 m, target straight up.
    #[test]
    fn test_fabrik_straight_chain() {
        let mut chain = IkChain::new(
            [0.0, 0.0, 0.0],
            &[(1.0, PI, 0.0), (1.0, PI, 0.0), (1.0, PI, 0.0)],
        );
        let solver = IkSolver::Fabrik {
            max_iterations: 50,
            tolerance: 1e-6,
        };
        let report = solver.solve(&mut chain, [0.0, 3.0, 0.0]);
        assert!(
            report.converged,
            "Should converge; residual={:.2e}",
            report.residual
        );
        assert!(
            report.residual < 1e-6,
            "Residual too large: {:.2e}",
            report.residual
        );
    }

    // 2. FABRIK unreachable: 2 × 1 m, target far away.
    #[test]
    fn test_fabrik_unreachable() {
        let mut chain = IkChain::new([0.0, 0.0, 0.0], &[(1.0, PI, 0.0), (1.0, PI, 0.0)]);
        let solver = IkSolver::Fabrik {
            max_iterations: 50,
            tolerance: 1e-6,
        };
        let report = solver.solve(&mut chain, [10.0, 0.0, 0.0]);
        assert!(!report.converged, "Should NOT converge");
        // Chain should be extended toward target (tip near x=2).
        let tip_pos = tip(&chain);
        assert!(
            tip_pos[0] > 1.8,
            "Chain should extend toward target; tip_x={:.3}",
            tip_pos[0]
        );
    }

    // 3. FABRIK partial reach: 3 × 1 m, target at [2.5, 0, 0].
    #[test]
    fn test_fabrik_partial_reach() {
        let mut chain = IkChain::new(
            [0.0, 0.0, 0.0],
            &[(1.0, PI, 0.0), (1.0, PI, 0.0), (1.0, PI, 0.0)],
        );
        let solver = IkSolver::Fabrik {
            max_iterations: 100,
            tolerance: 1e-4,
        };
        let report = solver.solve(&mut chain, [2.5, 0.0, 0.0]);
        assert!(
            report.residual < 1e-4,
            "Residual too large for partial-reach: {:.2e}",
            report.residual
        );
    }

    // 4. TwoBone 2-segment: L1=L2=1, target at [sqrt(2), 0, 0].
    #[test]
    fn test_two_bone_two_segment() {
        let mut chain = IkChain::new([0.0, 0.0, 0.0], &[(1.0, PI, 0.0), (1.0, PI, 0.0)]);
        let solver = IkSolver::TwoBone;
        let target = [2_f64.sqrt(), 0.0, 0.0];
        let report = solver.solve(&mut chain, target);
        assert!(
            report.residual < 1e-4,
            "TwoBone residual too large: {:.2e}",
            report.residual
        );
        // Verify segment lengths preserved.
        assert!(
            (len(chain.joints[0].local_offset) - 1.0).abs() < 1e-6,
            "Segment 0 length changed"
        );
        assert!(
            (len(chain.joints[1].local_offset) - 1.0).abs() < 1e-6,
            "Segment 1 length changed"
        );
    }

    // 5. TwoBone wrong length: 3-segment chain → residual=INFINITY.
    #[test]
    fn test_two_bone_wrong_length() {
        let mut chain = IkChain::new(
            [0.0, 0.0, 0.0],
            &[(1.0, PI, 0.0), (1.0, PI, 0.0), (1.0, PI, 0.0)],
        );
        let solver = IkSolver::TwoBone;
        let report = solver.solve(&mut chain, [1.0, 0.0, 0.0]);
        assert!(!report.converged);
        assert!(
            report.residual.is_infinite(),
            "Expected INFINITY residual, got {}",
            report.residual
        );
    }

    // 6. Cone limit: tight cone (0.1 rad) keeps bone angles small.
    #[test]
    fn test_cone_limit() {
        let cone = 0.1_f64; // very tight
        let mut chain = IkChain::new(
            [0.0, 0.0, 0.0],
            &[(1.0, cone, 0.0), (1.0, cone, 0.0), (1.0, cone, 0.0)],
        );
        let solver = IkSolver::Fabrik {
            max_iterations: 50,
            tolerance: 1e-6,
        };
        // Target off to the side — chain cannot fully reach it with tight cone.
        solver.solve(&mut chain, [3.0, 0.0, 0.0]);

        // Verify angles between successive bones stay within cone + small epsilon.
        let mut positions = vec![chain.root_world];
        for j in &chain.joints {
            let prev = *positions.last().expect("positions non-empty");
            positions.push(add(prev, j.local_offset));
        }

        for i in 1..chain.joints.len() {
            let d0 = normalize(sub(positions[i], positions[i - 1])).unwrap_or([0.0, 1.0, 0.0]);
            let d1 = normalize(sub(positions[i + 1], positions[i])).unwrap_or([0.0, 1.0, 0.0]);
            let angle = dot(d0, d1).clamp(-1.0, 1.0).acos();
            assert!(
                angle <= cone + 1e-3,
                "Bone {} angle {:.4} rad exceeds cone limit {:.4} rad",
                i,
                angle,
                cone
            );
        }
    }

    // 7. Serde round-trip for IkChain and IkSolver.
    #[test]
    fn test_serde_round_trip() {
        let chain = IkChain::new([1.0, 2.0, 3.0], &[(1.5, 0.5, 0.1), (2.0, PI, 0.0)]);

        let json = serde_json::to_string(&chain).expect("serialize IkChain");
        let chain2: IkChain = serde_json::from_str(&json).expect("deserialize IkChain");
        assert_eq!(chain, chain2);

        let solvers = [
            IkSolver::Fabrik {
                max_iterations: 20,
                tolerance: 1e-5,
            },
            IkSolver::TwoBone,
        ];
        for solver in &solvers {
            let json = serde_json::to_string(solver).expect("serialize IkSolver");
            let solver2: IkSolver = serde_json::from_str(&json).expect("deserialize IkSolver");
            assert_eq!(solver, &solver2);
        }

        // SolveReport round-trip
        let report = SolveReport {
            iterations: 5,
            converged: true,
            residual: 1e-7,
        };
        let json = serde_json::to_string(&report).expect("serialize SolveReport");
        let report2: SolveReport = serde_json::from_str(&json).expect("deserialize SolveReport");
        assert_eq!(report, report2);

        // IkError round-trip
        let err = IkError::TwoBoneRequiresTwoSegments;
        let json = serde_json::to_string(&err).expect("serialize IkError");
        let err2: IkError = serde_json::from_str(&json).expect("deserialize IkError");
        assert_eq!(err, err2);
    }
}
