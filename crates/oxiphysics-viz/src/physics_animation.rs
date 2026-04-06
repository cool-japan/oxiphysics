// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics-based animation and keyframe systems.
//!
//! Provides keyframe animation clips, blending, physics trajectory simulation,
//! morph targets, and FABRIK-based inverse kinematics.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// AnimTransform
// ---------------------------------------------------------------------------

/// A transform in 3-D space: position, rotation (quaternion), and scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnimTransform {
    /// World-space position \[x, y, z\].
    pub position: [f64; 3],
    /// Rotation as unit quaternion \[x, y, z, w\].
    pub rotation: [f64; 4],
    /// Per-axis scale \[sx, sy, sz\].
    pub scale: [f64; 3],
}

impl AnimTransform {
    /// Identity transform: no translation, no rotation, unit scale.
    pub fn identity() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
        }
    }

    /// Create a transform with only a translation.
    pub fn from_position(pos: [f64; 3]) -> Self {
        Self {
            position: pos,
            ..Self::identity()
        }
    }
}

impl Default for AnimTransform {
    fn default() -> Self {
        Self::identity()
    }
}

// ---------------------------------------------------------------------------
// EasingFn
// ---------------------------------------------------------------------------

/// Easing function for keyframe interpolation.
#[derive(Debug, Clone, Copy)]
pub enum EasingFn {
    /// Constant velocity — no acceleration or deceleration.
    Linear,
    /// Starts slow, ends fast (quadratic ease-in).
    EaseIn,
    /// Starts fast, ends slow (quadratic ease-out).
    EaseOut,
    /// Slow at both ends (cubic ease-in-out).
    EaseInOut,
    /// Cubic Bézier defined by two control-point x/y pairs `[x1,y1,x2,y2]`.
    Bezier([f64; 4]),
    /// Spring-based easing driven by physical spring parameters.
    Spring {
        /// Spring stiffness constant k.
        stiffness: f64,
        /// Damping coefficient c.
        damping: f64,
    },
}

impl EasingFn {
    /// Map a normalized time `t ∈ [0,1]` through the easing curve → `[0,1]`.
    pub fn apply(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            EasingFn::Linear => t,
            EasingFn::EaseIn => t * t,
            EasingFn::EaseOut => t * (2.0 - t),
            EasingFn::EaseInOut => {
                if t < 0.5 {
                    2.0 * t * t
                } else {
                    -1.0 + (4.0 - 2.0 * t) * t
                }
            }
            EasingFn::Bezier([x1, y1, x2, y2]) => {
                // Numerical solve for Bézier parameter using Newton's method, then evaluate y(u).
                let solve_t = |x_target: f64| -> f64 {
                    let mut u = x_target;
                    for _ in 0..8 {
                        let bx = 3.0 * (1.0 - u).powi(2) * u * x1
                            + 3.0 * (1.0 - u) * u * u * x2
                            + u * u * u;
                        let dx = 3.0 * (1.0 - u).powi(2) * x1
                            + 6.0 * (1.0 - u) * u * (x2 - x1)
                            + 3.0 * u * u * (1.0 - x2);
                        if dx.abs() < 1e-12 {
                            break;
                        }
                        u -= (bx - x_target) / dx;
                        u = u.clamp(0.0, 1.0);
                    }
                    u
                };
                let u = solve_t(t);
                3.0 * (1.0 - u).powi(2) * u * y1 + 3.0 * (1.0 - u) * u * u * y2 + u * u * u
            }
            EasingFn::Spring { stiffness, damping } => {
                // Spring settling toward 1 from 0 with zero initial velocity.
                let omega_n = stiffness.sqrt().max(1e-6);
                let zeta = damping / (2.0 * omega_n);
                let env = (-zeta * omega_n * t).exp();
                let disp = if zeta < 0.9999 {
                    // Under-damped
                    let omega_d = omega_n * (1.0 - zeta * zeta).sqrt();
                    env * ((omega_d * t).cos() + zeta / omega_d.max(1e-12) * (omega_d * t).sin())
                } else if zeta < 1.0001 {
                    // Critically damped: x(t) = (1 + ω_n * t) * e^(-ω_n * t)
                    (1.0 + omega_n * t) * (-omega_n * t).exp()
                } else {
                    // Over-damped
                    let disc = (zeta * zeta - 1.0).max(0.0).sqrt();
                    let r1 = -omega_n * (zeta - disc);
                    let r2 = -omega_n * (zeta + disc);
                    let denom = (r2 - r1).abs().max(1e-12);
                    let c1 = r2 / denom;
                    let c2 = -r1 / denom;
                    c1 * (r1 * t).exp() + c2 * (r2 * t).exp()
                };
                (1.0 - disp).clamp(0.0, 1.0)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Keyframe
// ---------------------------------------------------------------------------

/// A single keyframe in an animation clip.
#[derive(Debug, Clone)]
pub struct Keyframe {
    /// Time in seconds at which this keyframe occurs.
    pub time: f64,
    /// Transform value at this keyframe.
    pub transform: AnimTransform,
    /// Easing function applied between this keyframe and the next.
    pub easing: EasingFn,
}

impl Keyframe {
    /// Create a new keyframe at the given time with a linear easing.
    pub fn new(time: f64, transform: AnimTransform) -> Self {
        Self {
            time,
            transform,
            easing: EasingFn::Linear,
        }
    }

    /// Create a keyframe with a custom easing.
    pub fn with_easing(time: f64, transform: AnimTransform, easing: EasingFn) -> Self {
        Self {
            time,
            transform,
            easing,
        }
    }
}

// ---------------------------------------------------------------------------
// AnimationClip
// ---------------------------------------------------------------------------

/// A named animation clip composed of an ordered sequence of [`Keyframe`]s.
#[derive(Debug, Clone)]
pub struct AnimationClip {
    /// Human-readable clip name.
    pub name: String,
    /// Ordered keyframes (sorted by time).
    pub keyframes: Vec<Keyframe>,
    /// Total duration in seconds.
    pub duration: f64,
    /// Whether the clip loops back to the start at the end.
    pub looping: bool,
}

impl AnimationClip {
    /// Create a new, empty animation clip.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            keyframes: Vec::new(),
            duration: 0.0,
            looping: false,
        }
    }

    /// Add a keyframe, maintaining sorted order by time.
    pub fn add_keyframe(&mut self, kf: Keyframe) {
        let pos = self.keyframes.partition_point(|k| k.time < kf.time);
        if kf.time > self.duration {
            self.duration = kf.time;
        }
        self.keyframes.insert(pos, kf);
    }

    /// Sample the animation at time `t` (seconds).
    ///
    /// Returns the interpolated [`AnimTransform`] at that moment.
    pub fn sample(&self, t: f64) -> AnimTransform {
        if self.keyframes.is_empty() {
            return AnimTransform::identity();
        }

        let t = if self.looping && self.duration > 1e-12 {
            t % self.duration
        } else {
            t.clamp(0.0, self.duration)
        };

        // Binary search for the surrounding keyframes
        let pos = self.keyframes.partition_point(|k| k.time <= t);

        if pos == 0 {
            return self.keyframes[0].transform;
        }
        if pos >= self.keyframes.len() {
            return self
                .keyframes
                .last()
                .expect("collection should not be empty")
                .transform;
        }

        let kf_a = &self.keyframes[pos - 1];
        let kf_b = &self.keyframes[pos];
        let span = kf_b.time - kf_a.time;

        let local_t = if span > 1e-12 {
            (t - kf_a.time) / span
        } else {
            1.0
        };

        let eased = kf_a.easing.apply(local_t);
        lerp_transform(&kf_a.transform, &kf_b.transform, eased)
    }

    /// Return the total duration of the clip.
    pub fn duration(&self) -> f64 {
        self.duration
    }
}

// ---------------------------------------------------------------------------
// AnimationBlender
// ---------------------------------------------------------------------------

/// Blends multiple animation clips together using weighted mixing.
#[derive(Debug, Clone)]
pub struct AnimationBlender {
    /// Animation clips to blend.
    pub clips: Vec<AnimationClip>,
    /// Per-clip blend weights (should sum to 1 for best results).
    pub weights: Vec<f64>,
}

impl AnimationBlender {
    /// Create a new blender with the given clips (equal weights).
    pub fn new(clips: Vec<AnimationClip>) -> Self {
        let n = clips.len();
        let weights = if n > 0 {
            vec![1.0 / n as f64; n]
        } else {
            vec![]
        };
        Self { clips, weights }
    }

    /// Set explicit blend weights.
    pub fn set_weights(&mut self, weights: Vec<f64>) {
        self.weights = weights;
    }

    /// Sample all clips at time `t` and blend by weights.
    pub fn blend(&self, t: f64) -> AnimTransform {
        if self.clips.is_empty() {
            return AnimTransform::identity();
        }

        let w_sum: f64 = self.weights.iter().sum();
        let w_sum = if w_sum < 1e-12 { 1.0 } else { w_sum };

        let mut out = AnimTransform::identity();
        out.position = [0.0; 3];
        out.rotation = [0.0; 4];
        out.scale = [0.0; 3];

        let mut total_w = 0.0_f64;
        for (i, clip) in self.clips.iter().enumerate() {
            let w = self.weights.get(i).copied().unwrap_or(0.0) / w_sum;
            if w.abs() < 1e-12 {
                continue;
            }
            let s = clip.sample(t);

            for k in 0..3 {
                out.position[k] += w * s.position[k];
            }
            for k in 0..3 {
                out.scale[k] += w * s.scale[k];
            }
            // Accumulate rotation: first clip sets base, then slerp blend
            if total_w == 0.0 {
                out.rotation = s.rotation;
            } else {
                let slerp_t = w / (total_w + w);
                out.rotation = slerp_quat(out.rotation, s.rotation, slerp_t);
            }
            total_w += w;
        }
        out
    }

    /// Crossfade between clip `from` and clip `to` with the given progress `[0,1]`.
    pub fn crossfade(&self, from: usize, to: usize, progress: f64) -> AnimTransform {
        let t_from = self
            .clips
            .get(from)
            .map_or(AnimTransform::identity(), |c| c.sample(0.0));
        let t_to = self
            .clips
            .get(to)
            .map_or(AnimTransform::identity(), |c| c.sample(0.0));
        lerp_transform(&t_from, &t_to, progress.clamp(0.0, 1.0))
    }
}

// ---------------------------------------------------------------------------
// PhysicsAnimator
// ---------------------------------------------------------------------------

/// Simulates rigid-body trajectories with gravity and linear damping.
#[derive(Debug, Clone)]
pub struct PhysicsAnimator {
    /// Gravitational acceleration vector \[gx, gy, gz\] in m/s².
    pub gravity: [f64; 3],
    /// Linear damping coefficient (0 = no damping, 1 = critical).
    pub damping: f64,
}

impl PhysicsAnimator {
    /// Create a new `PhysicsAnimator` with standard Earth gravity (–y).
    pub fn new() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            damping: 0.0,
        }
    }

    /// Create a `PhysicsAnimator` with custom gravity.
    pub fn with_gravity(gravity: [f64; 3]) -> Self {
        Self {
            gravity,
            damping: 0.0,
        }
    }

    /// Simulate a ballistic trajectory using the Euler method.
    ///
    /// Returns the list of positions at each time step (length = `n_steps + 1`).
    pub fn simulate_trajectory(
        &self,
        initial_pos: [f64; 3],
        initial_vel: [f64; 3],
        dt: f64,
        n_steps: usize,
    ) -> Vec<[f64; 3]> {
        let mut positions = Vec::with_capacity(n_steps + 1);
        let mut pos = initial_pos;
        let mut vel = initial_vel;
        positions.push(pos);

        let d = 1.0 - self.damping.clamp(0.0, 0.9999) * dt;

        for _ in 0..n_steps {
            vel[0] = vel[0] * d + self.gravity[0] * dt;
            vel[1] = vel[1] * d + self.gravity[1] * dt;
            vel[2] = vel[2] * d + self.gravity[2] * dt;
            pos[0] += vel[0] * dt;
            pos[1] += vel[1] * dt;
            pos[2] += vel[2] * dt;
            positions.push(pos);
        }
        positions
    }
}

impl Default for PhysicsAnimator {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Interpolation helpers
// ---------------------------------------------------------------------------

/// Linearly interpolate between two [`AnimTransform`]s.
///
/// Position and scale are interpolated linearly; rotation uses [`slerp_quat`].
pub fn lerp_transform(a: &AnimTransform, b: &AnimTransform, t: f64) -> AnimTransform {
    let lerp1 = |x: f64, y: f64| x + (y - x) * t;
    AnimTransform {
        position: [
            lerp1(a.position[0], b.position[0]),
            lerp1(a.position[1], b.position[1]),
            lerp1(a.position[2], b.position[2]),
        ],
        rotation: slerp_quat(a.rotation, b.rotation, t),
        scale: [
            lerp1(a.scale[0], b.scale[0]),
            lerp1(a.scale[1], b.scale[1]),
            lerp1(a.scale[2], b.scale[2]),
        ],
    }
}

/// Spherical linear interpolation between two unit quaternions.
///
/// Returns the quaternion that rotates from `a` to `b` by fraction `t ∈ [0,1]`.
pub fn slerp_quat(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    let t = t.clamp(0.0, 1.0);

    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];

    // Ensure shortest path
    let (b, dot) = if dot < 0.0 {
        ([-b[0], -b[1], -b[2], -b[3]], -dot)
    } else {
        (b, dot)
    };

    let dot = dot.clamp(-1.0, 1.0);

    // If quaternions are very close, fall back to linear lerp + normalize
    if dot > 0.9995 {
        let r = [
            a[0] + (b[0] - a[0]) * t,
            a[1] + (b[1] - a[1]) * t,
            a[2] + (b[2] - a[2]) * t,
            a[3] + (b[3] - a[3]) * t,
        ];
        return normalize_quat(r);
    }

    let theta_0 = dot.acos(); // angle between a and b
    let theta = theta_0 * t; // angle for interpolation
    let sin_theta_0 = theta_0.sin().max(1e-12);

    let s0 = (theta_0 - theta).sin() / sin_theta_0;
    let s1 = theta.sin() / sin_theta_0;

    normalize_quat([
        s0 * a[0] + s1 * b[0],
        s0 * a[1] + s1 * b[1],
        s0 * a[2] + s1 * b[2],
        s0 * a[3] + s1 * b[3],
    ])
}

/// Normalize a quaternion to unit length.
fn normalize_quat(q: [f64; 4]) -> [f64; 4] {
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if len < 1e-12 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / len, q[1] / len, q[2] / len, q[3] / len]
    }
}

// ---------------------------------------------------------------------------
// MorphTarget
// ---------------------------------------------------------------------------

/// A morph target (blend shape) that stores per-vertex displacement deltas.
#[derive(Debug, Clone)]
pub struct MorphTarget {
    /// Name of this morph target.
    pub name: String,
    /// Per-vertex displacement vectors \[dx, dy, dz\].
    pub deltas: Vec<[f64; 3]>,
    /// Current blend weight in \[0, 1\].
    pub weight: f64,
}

impl MorphTarget {
    /// Create a new morph target with the given name and deltas.
    pub fn new(name: impl Into<String>, deltas: Vec<[f64; 3]>) -> Self {
        Self {
            name: name.into(),
            deltas,
            weight: 0.0,
        }
    }

    /// Apply this morph target to a base vertex list, returning the blended positions.
    pub fn apply_to(&self, base: &[[f64; 3]]) -> Vec<[f64; 3]> {
        base.iter()
            .enumerate()
            .map(|(i, &p)| {
                let d = self.deltas.get(i).copied().unwrap_or([0.0; 3]);
                [
                    p[0] + self.weight * d[0],
                    p[1] + self.weight * d[1],
                    p[2] + self.weight * d[2],
                ]
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// InverseKinematics (FABRIK)
// ---------------------------------------------------------------------------

/// A simple linear bone chain for FABRIK inverse kinematics.
#[derive(Debug, Clone)]
pub struct InverseKinematics {
    /// Joint positions in the chain (root first, end-effector last).
    pub chain: Vec<[f64; 3]>,
    /// IK target position for the end-effector.
    pub target: [f64; 3],
}

impl InverseKinematics {
    /// Create a new IK chain.
    pub fn new(chain: Vec<[f64; 3]>, target: [f64; 3]) -> Self {
        Self { chain, target }
    }

    /// Run FABRIK (Forward And Backward Reaching IK) for `iterations` passes.
    ///
    /// Returns the updated joint positions.
    #[allow(clippy::too_many_arguments)]
    pub fn fabrik_solve(&self, iterations: usize) -> Vec<[f64; 3]> {
        let n = self.chain.len();
        if n < 2 {
            return self.chain.clone();
        }

        // Precompute segment lengths
        let lengths: Vec<f64> = self
            .chain
            .windows(2)
            .map(|w| {
                let d = sub3(w[1], w[0]);
                len3(d)
            })
            .collect();

        // Check reachability
        let total_length: f64 = lengths.iter().sum();
        let root_to_target = len3(sub3(self.target, self.chain[0]));

        let mut joints = self.chain.clone();

        if root_to_target >= total_length {
            // Fully stretched toward target
            for i in 1..n {
                let dir = normalize3(sub3(self.target, joints[i - 1]));
                joints[i] = add3(joints[i - 1], scale3(dir, lengths[i - 1]));
            }
            return joints;
        }

        let root = joints[0];

        for _ in 0..iterations {
            // Forward pass: pull end-effector to target
            joints[n - 1] = self.target;
            for i in (0..n - 1).rev() {
                let dir = normalize3(sub3(joints[i], joints[i + 1]));
                joints[i] = add3(joints[i + 1], scale3(dir, lengths[i]));
            }

            // Backward pass: fix root
            joints[0] = root;
            for i in 0..n - 1 {
                let dir = normalize3(sub3(joints[i + 1], joints[i]));
                joints[i + 1] = add3(joints[i], scale3(dir, lengths[i]));
            }

            // Check convergence
            if len3(sub3(joints[n - 1], self.target)) < 1e-6 {
                break;
            }
        }

        joints
    }
}

// ---------------------------------------------------------------------------
// 3-D vector helpers (private)
// ---------------------------------------------------------------------------

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(v, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- EasingFn ---

    #[test]
    fn test_easing_linear_endpoints() {
        assert!((EasingFn::Linear.apply(0.0) - 0.0).abs() < 1e-12);
        assert!((EasingFn::Linear.apply(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_easing_linear_midpoint() {
        assert!((EasingFn::Linear.apply(0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_easing_ease_in_slower_than_linear() {
        // EaseIn at t=0.5 should be less than linear (0.5)
        let val = EasingFn::EaseIn.apply(0.5);
        assert!(val < 0.5, "ease_in(0.5) = {val}");
    }

    #[test]
    fn test_easing_ease_out_faster_than_linear() {
        let val = EasingFn::EaseOut.apply(0.5);
        assert!(val > 0.5, "ease_out(0.5) = {val}");
    }

    #[test]
    fn test_easing_ease_in_out_symmetric() {
        let a = EasingFn::EaseInOut.apply(0.25);
        let b = EasingFn::EaseInOut.apply(0.75);
        // Symmetry: f(0.25) + f(0.75) ≈ 1
        assert!((a + b - 1.0).abs() < 1e-10, "a+b = {}", a + b);
    }

    #[test]
    fn test_easing_spring_starts_at_zero() {
        let e = EasingFn::Spring {
            stiffness: 100.0,
            damping: 20.0,
        };
        let v = e.apply(0.0);
        assert!(v.abs() < 0.1, "spring(0) = {v}");
    }

    #[test]
    fn test_easing_spring_ends_near_one() {
        let e = EasingFn::Spring {
            stiffness: 100.0,
            damping: 20.0,
        };
        let v = e.apply(1.0);
        assert!(v > 0.5, "spring(1) = {v}");
    }

    #[test]
    fn test_easing_bezier_endpoints() {
        let e = EasingFn::Bezier([0.25, 0.1, 0.25, 1.0]);
        assert!((e.apply(0.0) - 0.0).abs() < 1e-6);
        // Note: Bézier at t=1.0 — the Bezier curve passes through (1,1)
        assert!(e.apply(1.0) > 0.9, "bezier(1.0) should be near 1");
    }

    // --- AnimTransform ---

    #[test]
    fn test_anim_transform_identity() {
        let id = AnimTransform::identity();
        assert_eq!(id.position, [0.0; 3]);
        assert_eq!(id.rotation, [0.0, 0.0, 0.0, 1.0]);
        assert_eq!(id.scale, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_anim_transform_from_position() {
        let t = AnimTransform::from_position([1.0, 2.0, 3.0]);
        assert_eq!(t.position, [1.0, 2.0, 3.0]);
        assert_eq!(t.rotation, [0.0, 0.0, 0.0, 1.0]);
    }

    // --- lerp_transform ---

    #[test]
    fn test_lerp_transform_at_zero() {
        let a = AnimTransform::from_position([0.0, 0.0, 0.0]);
        let b = AnimTransform::from_position([10.0, 0.0, 0.0]);
        let r = lerp_transform(&a, &b, 0.0);
        assert!((r.position[0]).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_transform_at_one() {
        let a = AnimTransform::from_position([0.0, 0.0, 0.0]);
        let b = AnimTransform::from_position([10.0, 0.0, 0.0]);
        let r = lerp_transform(&a, &b, 1.0);
        assert!((r.position[0] - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_lerp_transform_midpoint() {
        let a = AnimTransform::from_position([0.0, 0.0, 0.0]);
        let b = AnimTransform::from_position([4.0, 0.0, 0.0]);
        let r = lerp_transform(&a, &b, 0.5);
        assert!((r.position[0] - 2.0).abs() < 1e-12);
    }

    // --- slerp_quat ---

    #[test]
    fn test_slerp_quat_identity_to_identity() {
        let id = [0.0, 0.0, 0.0, 1.0_f64];
        let r = slerp_quat(id, id, 0.5);
        assert!((r[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_slerp_quat_at_zero_returns_a() {
        let a = [0.0_f64, 0.0, 0.0, 1.0];
        let b = [0.0_f64, 0.0, 1.0, 0.0];
        let r = slerp_quat(a, b, 0.0);
        assert!((r[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_slerp_quat_at_one_returns_b() {
        let a = [0.0_f64, 0.0, 0.0, 1.0];
        let b = [1.0_f64, 0.0, 0.0, 0.0];
        let r = slerp_quat(a, b, 1.0);
        assert!((r[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_slerp_quat_result_is_unit() {
        let a = [0.0_f64, 0.0, 0.0, 1.0];
        let b = [0.5_f64, 0.5, 0.5, 0.5]; // unit quaternion
        let r = slerp_quat(a, b, 0.3);
        let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "slerp result not unit: {len}");
    }

    // --- AnimationClip ---

    #[test]
    fn test_clip_empty_returns_identity() {
        let clip = AnimationClip::new("test");
        let s = clip.sample(0.5);
        assert_eq!(s, AnimTransform::identity());
    }

    #[test]
    fn test_clip_single_keyframe() {
        let mut clip = AnimationClip::new("test");
        let tf = AnimTransform::from_position([5.0, 0.0, 0.0]);
        clip.add_keyframe(Keyframe::new(1.0, tf));
        let s = clip.sample(0.5);
        assert_eq!(s.position, [5.0, 0.0, 0.0]);
    }

    #[test]
    fn test_clip_two_keyframes_midpoint() {
        let mut clip = AnimationClip::new("test");
        clip.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        clip.add_keyframe(Keyframe::new(
            2.0,
            AnimTransform::from_position([4.0, 0.0, 0.0]),
        ));
        let s = clip.sample(1.0);
        assert!(
            (s.position[0] - 2.0).abs() < 1e-10,
            "midpoint x = {}",
            s.position[0]
        );
    }

    #[test]
    fn test_clip_duration_grows_with_keyframes() {
        let mut clip = AnimationClip::new("test");
        assert!((clip.duration() - 0.0).abs() < 1e-12);
        clip.add_keyframe(Keyframe::new(3.0, AnimTransform::identity()));
        assert!((clip.duration() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_clip_sample_before_start_returns_first() {
        let mut clip = AnimationClip::new("test");
        clip.add_keyframe(Keyframe::new(
            1.0,
            AnimTransform::from_position([7.0, 0.0, 0.0]),
        ));
        let s = clip.sample(-1.0);
        assert_eq!(s.position, [7.0, 0.0, 0.0]);
    }

    #[test]
    fn test_clip_sample_after_end_returns_last() {
        let mut clip = AnimationClip::new("test");
        clip.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        clip.add_keyframe(Keyframe::new(
            1.0,
            AnimTransform::from_position([10.0, 0.0, 0.0]),
        ));
        let s = clip.sample(5.0);
        assert!((s.position[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_clip_looping() {
        let mut clip = AnimationClip::new("loop");
        clip.looping = true;
        clip.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        clip.add_keyframe(Keyframe::new(
            1.0,
            AnimTransform::from_position([1.0, 0.0, 0.0]),
        ));
        // Sample at t=1.5 should be same as t=0.5
        let s1 = clip.sample(0.5);
        let s2 = clip.sample(1.5);
        assert!(
            (s1.position[0] - s2.position[0]).abs() < 1e-10,
            "looped sample mismatch: {} vs {}",
            s1.position[0],
            s2.position[0]
        );
    }

    // --- AnimationBlender ---

    #[test]
    fn test_blender_empty_returns_identity() {
        let blender = AnimationBlender::new(vec![]);
        let s = blender.blend(0.0);
        assert_eq!(s, AnimTransform::identity());
    }

    #[test]
    fn test_blender_single_clip() {
        let mut clip = AnimationClip::new("c");
        clip.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([3.0, 0.0, 0.0]),
        ));
        let blender = AnimationBlender::new(vec![clip]);
        let s = blender.blend(0.0);
        assert!((s.position[0] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_blender_crossfade_at_zero() {
        let mut c1 = AnimationClip::new("a");
        c1.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        let mut c2 = AnimationClip::new("b");
        c2.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([10.0, 0.0, 0.0]),
        ));
        let blender = AnimationBlender::new(vec![c1, c2]);
        let r = blender.crossfade(0, 1, 0.0);
        assert!((r.position[0]).abs() < 1e-6);
    }

    #[test]
    fn test_blender_crossfade_at_one() {
        let mut c1 = AnimationClip::new("a");
        c1.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        let mut c2 = AnimationClip::new("b");
        c2.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([10.0, 0.0, 0.0]),
        ));
        let blender = AnimationBlender::new(vec![c1, c2]);
        let r = blender.crossfade(0, 1, 1.0);
        assert!((r.position[0] - 10.0).abs() < 1e-6);
    }

    // --- PhysicsAnimator ---

    #[test]
    fn test_trajectory_length() {
        let pa = PhysicsAnimator::new();
        let traj = pa.simulate_trajectory([0.0; 3], [1.0, 0.0, 0.0], 0.01, 100);
        assert_eq!(traj.len(), 101);
    }

    #[test]
    fn test_trajectory_first_point_is_initial() {
        let pa = PhysicsAnimator::new();
        let init = [1.0, 2.0, 3.0];
        let traj = pa.simulate_trajectory(init, [0.0; 3], 0.01, 10);
        assert_eq!(traj[0], init);
    }

    #[test]
    fn test_trajectory_falls_under_gravity() {
        let pa = PhysicsAnimator::new();
        let traj = pa.simulate_trajectory([0.0, 10.0, 0.0], [0.0; 3], 0.1, 50);
        // After several steps, y should decrease
        assert!(traj.last().unwrap()[1] < traj[0][1]);
    }

    #[test]
    fn test_trajectory_horizontal_motion() {
        let pa = PhysicsAnimator::with_gravity([0.0; 3]);
        let traj = pa.simulate_trajectory([0.0; 3], [1.0, 0.0, 0.0], 1.0, 5);
        // No gravity: x should increase linearly
        assert!((traj[1][0] - 1.0).abs() < 1e-10);
        assert!((traj[5][0] - 5.0).abs() < 1e-10);
    }

    // --- MorphTarget ---

    #[test]
    fn test_morph_target_zero_weight() {
        let base = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut mt = MorphTarget::new("smile", vec![[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]);
        mt.weight = 0.0;
        let result = mt.apply_to(&base);
        assert_eq!(result[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_morph_target_full_weight() {
        let base = vec![[0.0, 0.0, 0.0]];
        let mut mt = MorphTarget::new("brow", vec![[0.0, 0.5, 0.0]]);
        mt.weight = 1.0;
        let result = mt.apply_to(&base);
        assert!((result[0][1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_morph_target_half_weight() {
        let base = vec![[0.0, 0.0, 0.0]];
        let mut mt = MorphTarget::new("x", vec![[2.0, 0.0, 0.0]]);
        mt.weight = 0.5;
        let result = mt.apply_to(&base);
        assert!((result[0][0] - 1.0).abs() < 1e-12);
    }

    // --- InverseKinematics (FABRIK) ---

    #[test]
    fn test_fabrik_two_bone_reaches_target() {
        // Two-bone arm: lengths 1 + 1 = 2, target at (0.5, 1.0, 0).
        // The chain starts slightly off-axis so FABRIK can bend the elbow.
        let chain = vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [1.0, 0.0, 0.0]];
        let target = [0.5, 1.0, 0.0];
        let ik = InverseKinematics::new(chain, target);
        let result = ik.fabrik_solve(30);
        let end = result.last().unwrap();
        let dist = len3(sub3(*end, target));
        assert!(dist < 0.05, "end effector dist to target = {dist}");
    }

    #[test]
    fn test_fabrik_preserves_root() {
        let chain = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]];
        let ik = InverseKinematics::new(chain, [1.0, 1.0, 0.0]);
        let result = ik.fabrik_solve(10);
        let root = result[0];
        assert!(len3(root) < 1e-10, "root should stay at origin");
    }

    #[test]
    fn test_fabrik_fully_extended() {
        // Target beyond reach: chain should stretch toward target
        let chain = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let ik = InverseKinematics::new(chain, [100.0, 0.0, 0.0]);
        let result = ik.fabrik_solve(5);
        assert_eq!(result.len(), 2);
        // End effector should point toward target
        assert!(result[1][0] > 0.0);
    }

    #[test]
    fn test_fabrik_single_joint_returns_chain() {
        let chain = vec![[0.0, 0.0, 0.0]];
        let ik = InverseKinematics::new(chain.clone(), [1.0, 0.0, 0.0]);
        let result = ik.fabrik_solve(5);
        assert_eq!(result, chain);
    }

    // --- Vector helpers ---

    #[test]
    fn test_len3_zero() {
        assert!(len3([0.0; 3]) < 1e-12);
    }

    #[test]
    fn test_normalize3_unit_vector() {
        let v = normalize3([3.0, 0.0, 4.0]);
        let l = len3(v);
        assert!((l - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_keyframe_ordering() {
        let mut clip = AnimationClip::new("ord");
        clip.add_keyframe(Keyframe::new(
            2.0,
            AnimTransform::from_position([2.0, 0.0, 0.0]),
        ));
        clip.add_keyframe(Keyframe::new(
            0.0,
            AnimTransform::from_position([0.0, 0.0, 0.0]),
        ));
        clip.add_keyframe(Keyframe::new(
            1.0,
            AnimTransform::from_position([1.0, 0.0, 0.0]),
        ));
        // After sorting, times should be 0, 1, 2
        assert!((clip.keyframes[0].time - 0.0).abs() < 1e-12);
        assert!((clip.keyframes[1].time - 1.0).abs() < 1e-12);
        assert!((clip.keyframes[2].time - 2.0).abs() < 1e-12);
    }
}
