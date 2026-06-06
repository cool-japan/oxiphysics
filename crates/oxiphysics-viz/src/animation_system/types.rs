//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::functions::*;
use super::functions::{Mat4, Quat, Vec3};

/// LOD (level-of-detail) level for animation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AnimLod {
    /// Full quality: all bones, full keyframes.
    Full,
    /// Medium quality: reduced keyframes.
    Medium,
    /// Low quality: only root motion.
    Low,
    /// Disabled: no animation (static bind pose).
    Disabled,
}
/// Interpolation mode for a keyframe segment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InterpMode {
    /// Step (hold previous value until next key).
    Step,
    /// Linear interpolation.
    Linear,
    /// Cubic Hermite spline.
    Hermite,
    /// Cubic Bézier (two in/out tangent handles per segment).
    Bezier,
}
/// High-level animator combining a clip library, state machine, skeleton, and LOD.
#[derive(Debug, Clone)]
pub struct Animator {
    /// Registered clips by name.
    pub clips: HashMap<String, AnimationClip>,
    /// Animation state machine.
    pub state_machine: AnimStateMachine,
    /// Skeleton.
    pub skeleton: Skeleton,
    /// Current local bone matrices.
    pub local_matrices: Vec<Mat4>,
    /// Current world bone matrices.
    pub world_matrices: Vec<Mat4>,
    /// Current skinning matrices.
    pub skinning_matrices: Vec<Mat4>,
}
impl Animator {
    /// Create an animator from a skeleton.
    pub fn new(skeleton: Skeleton) -> Self {
        let n = skeleton.bones.len();
        let ident = identity_mat4();
        Self {
            clips: HashMap::new(),
            state_machine: AnimStateMachine::new(),
            local_matrices: vec![ident; n],
            world_matrices: vec![ident; n],
            skinning_matrices: vec![ident; n],
            skeleton,
        }
    }
    /// Register a clip.
    pub fn add_clip(&mut self, clip: AnimationClip) {
        self.clips.insert(clip.name.clone(), clip);
    }
    /// Advance the animation by `dt` seconds and recompute bone matrices.
    pub fn update(&mut self, dt: f64) {
        let clips_snapshot = self.clips.clone();
        self.state_machine.update(dt, &clips_snapshot);
        let pose = self.state_machine.evaluate(&self.clips);
        for (i, bone) in self.skeleton.bones.iter().enumerate() {
            if let Some(&mat) = pose.get(&bone.id) {
                self.local_matrices[i] = mat;
            }
        }
        self.world_matrices = self.skeleton.compute_world_matrices(&self.local_matrices);
        self.skinning_matrices = self
            .skeleton
            .compute_skinning_matrices(&self.local_matrices);
    }
    /// Fire a state-machine trigger.
    pub fn fire_trigger(&mut self, trigger: &str) {
        self.state_machine.fire_trigger(trigger);
    }
}
/// A single Vec3 keyframe.
#[derive(Debug, Clone, Copy)]
pub struct Vec3Key {
    /// Time of the key (seconds).
    pub time: f64,
    /// Value at this key.
    pub value: Vec3,
    /// Interpolation mode to use from this key to the next.
    pub interp: InterpMode,
    /// Optional out-tangent for Hermite.
    pub out_tangent: Vec3,
    /// Optional in-tangent for Hermite.
    pub in_tangent: Vec3,
}
impl Vec3Key {
    /// Construct a linear Vec3 keyframe.
    pub fn linear(time: f64, value: Vec3) -> Self {
        Self {
            time,
            value,
            interp: InterpMode::Linear,
            out_tangent: [0.0; 3],
            in_tangent: [0.0; 3],
        }
    }
}
/// A node in a blend tree.
#[derive(Debug, Clone)]
pub enum BlendNode {
    /// Leaf: play a single clip.
    Clip {
        /// Name of the clip to play.
        clip_name: String,
        /// Current playback time.
        time: f64,
    },
    /// Linear blend between two child nodes.
    Lerp {
        /// First child.
        child_a: Box<BlendNode>,
        /// Second child.
        child_b: Box<BlendNode>,
        /// Blend weight (0 → child_a, 1 → child_b).
        weight: f64,
    },
    /// Weighted additive blend of N children.
    Additive {
        /// Child nodes with associated weights.
        children: Vec<(Box<BlendNode>, f64)>,
    },
}
impl BlendNode {
    /// Evaluate this blend node against a clip library.
    /// Returns a map of node_id → Mat4.
    pub fn evaluate(&self, clips: &HashMap<String, AnimationClip>) -> HashMap<u32, Mat4> {
        match self {
            Self::Clip { clip_name, time } => {
                if let Some(clip) = clips.get(clip_name) {
                    clip.sample(*time)
                } else {
                    HashMap::new()
                }
            }
            Self::Lerp {
                child_a,
                child_b,
                weight,
            } => {
                let a = child_a.evaluate(clips);
                let b = child_b.evaluate(clips);
                let mut out = a.clone();
                for (&id, mb) in &b {
                    let ma = out.entry(id).or_insert(*mb);
                    *ma = lerp_mat4(*ma, *mb, *weight);
                }
                out
            }
            Self::Additive { children } => {
                let mut out: HashMap<u32, Mat4> = HashMap::new();
                for (child, w) in children {
                    let result = child.evaluate(clips);
                    for (id, m) in result {
                        let entry = out.entry(id).or_insert(identity_mat4());
                        *entry = lerp_mat4(*entry, m, *w);
                    }
                }
                out
            }
        }
    }
}
/// A single motion-capture frame: per-bone rotation quaternions.
#[derive(Debug, Clone)]
pub struct MocapFrame {
    /// Frame timestamp in seconds.
    pub time: f64,
    /// Bone rotations indexed by bone ID.
    pub bone_rotations: HashMap<u32, Quat>,
    /// Root position.
    pub root_position: Vec3,
}
impl MocapFrame {
    /// Create a new empty frame.
    pub fn new(time: f64, root_position: Vec3) -> Self {
        Self {
            time,
            bone_rotations: HashMap::new(),
            root_position,
        }
    }
    /// Set bone rotation.
    pub fn set_bone(&mut self, bone_id: u32, rotation: Quat) {
        self.bone_rotations.insert(bone_id, rotation);
    }
}
/// Animation LOD manager: tracks per-entity LOD.
#[derive(Debug, Clone, Default)]
pub struct AnimLodManager {
    /// Entity ID → current LOD level.
    pub lods: HashMap<u32, AnimLod>,
}
impl AnimLodManager {
    /// Create a new manager.
    pub fn new() -> Self {
        Self::default()
    }
    /// Update the LOD for an entity given its camera distance.
    pub fn update_entity(&mut self, entity_id: u32, distance: f64) {
        self.lods.insert(entity_id, lod_for_distance(distance));
    }
    /// Get the LOD for an entity (default: Full).
    pub fn get_lod(&self, entity_id: u32) -> AnimLod {
        self.lods.get(&entity_id).copied().unwrap_or(AnimLod::Full)
    }
    /// Return the update rate multiplier (how often to update) for a given LOD.
    pub fn update_rate(lod: AnimLod) -> f64 {
        match lod {
            AnimLod::Full => 1.0,
            AnimLod::Medium => 0.5,
            AnimLod::Low => 0.25,
            AnimLod::Disabled => 0.0,
        }
    }
}
/// A sorted sequence of Vec3 keyframes forming an animation curve.
#[derive(Debug, Clone, Default)]
pub struct Vec3Curve {
    /// Keyframes sorted by time.
    pub keys: Vec<Vec3Key>,
}
impl Vec3Curve {
    /// Create an empty curve.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a keyframe.
    pub fn push_key(&mut self, key: Vec3Key) {
        self.keys.push(key);
        self.keys.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Evaluate the curve at time `t`.
    pub fn sample(&self, t: f64) -> Vec3 {
        let n = self.keys.len();
        if n == 0 {
            return [0.0; 3];
        }
        if t <= self.keys[0].time {
            return self.keys[0].value;
        }
        if t >= self.keys[n - 1].time {
            return self.keys[n - 1].value;
        }
        let pos = self.keys.partition_point(|k| k.time <= t);
        let i = pos - 1;
        let j = i + 1;
        let ki = &self.keys[i];
        let kj = &self.keys[j];
        let dt = kj.time - ki.time;
        let raw_t = if dt.abs() < 1e-300 {
            0.0
        } else {
            (t - ki.time) / dt
        };
        match ki.interp {
            InterpMode::Step => ki.value,
            InterpMode::Linear => vec3_lerp(ki.value, kj.value, raw_t),
            InterpMode::Hermite => {
                let mut out = [0.0; 3];
                for (d, o) in out.iter_mut().enumerate() {
                    *o = cubic_hermite(
                        ki.value[d],
                        ki.out_tangent[d] * dt,
                        kj.value[d],
                        kj.in_tangent[d] * dt,
                        raw_t,
                    );
                }
                out
            }
            InterpMode::Bezier => vec3_lerp(ki.value, kj.value, raw_t),
        }
    }
}
/// A motion-capture sequence consisting of ordered frames.
#[derive(Debug, Clone, Default)]
pub struct MocapSequence {
    /// Ordered frames.
    pub frames: Vec<MocapFrame>,
    /// Frames per second (informational).
    pub fps: f64,
    /// Whether the sequence loops.
    pub looping: bool,
}
impl MocapSequence {
    /// Create a new empty sequence.
    pub fn new(fps: f64, looping: bool) -> Self {
        Self {
            frames: Vec::new(),
            fps,
            looping,
        }
    }
    /// Append a frame.
    pub fn push_frame(&mut self, frame: MocapFrame) {
        self.frames.push(frame);
    }
    /// Duration in seconds.
    pub fn duration(&self) -> f64 {
        self.frames.last().map(|f| f.time).unwrap_or(0.0)
    }
    /// Sample the sequence at time `t` with SLERP between adjacent frames.
    pub fn sample(&self, t: f64) -> Option<MocapFrame> {
        let n = self.frames.len();
        if n == 0 {
            return None;
        }
        let t = if self.looping {
            let dur = self.duration();
            if dur > 0.0 { t.rem_euclid(dur) } else { 0.0 }
        } else {
            t.clamp(0.0, self.duration())
        };
        if n == 1 || t <= self.frames[0].time {
            return Some(self.frames[0].clone());
        }
        if t >= self.frames[n - 1].time {
            return Some(self.frames[n - 1].clone());
        }
        let pos = self.frames.partition_point(|f| f.time <= t);
        let i = pos - 1;
        let j = i + 1;
        let fi = &self.frames[i];
        let fj = &self.frames[j];
        let dt = fj.time - fi.time;
        let alpha = if dt.abs() < 1e-300 {
            0.0
        } else {
            (t - fi.time) / dt
        };
        let root = vec3_lerp(fi.root_position, fj.root_position, alpha);
        let mut out = MocapFrame::new(t, root);
        let mut bone_ids: Vec<u32> = fi.bone_rotations.keys().cloned().collect();
        for id in fj.bone_rotations.keys() {
            if !bone_ids.contains(id) {
                bone_ids.push(*id);
            }
        }
        for bid in bone_ids {
            let qa = fi
                .bone_rotations
                .get(&bid)
                .copied()
                .unwrap_or([0.0, 0.0, 0.0, 1.0]);
            let qb = fj
                .bone_rotations
                .get(&bid)
                .copied()
                .unwrap_or([0.0, 0.0, 0.0, 1.0]);
            out.set_bone(bid, quat_slerp(qa, qb, alpha));
        }
        Some(out)
    }
}
/// Easing function type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Easing {
    /// Linear (no easing).
    Linear,
    /// Quadratic ease-in.
    EaseInQuad,
    /// Quadratic ease-out.
    EaseOutQuad,
    /// Quadratic ease-in-out.
    EaseInOutQuad,
    /// Cubic ease-in.
    EaseInCubic,
    /// Cubic ease-out.
    EaseOutCubic,
    /// Cubic ease-in-out.
    EaseInOutCubic,
    /// Sinusoidal ease-in.
    EaseInSine,
    /// Sinusoidal ease-out.
    EaseOutSine,
    /// Sinusoidal ease-in-out.
    EaseInOutSine,
    /// Exponential ease-in.
    EaseInExpo,
    /// Exponential ease-out.
    EaseOutExpo,
    /// Elastic ease-out.
    EaseOutElastic,
    /// Bounce ease-out.
    EaseOutBounce,
}
/// A single bone in a skeletal hierarchy.
#[derive(Debug, Clone)]
pub struct Bone {
    /// Unique bone identifier.
    pub id: u32,
    /// Bone name (e.g. `"spine"`).
    pub name: String,
    /// Index of the parent bone, or `None` for roots.
    pub parent: Option<u32>,
    /// Inverse bind-pose matrix (skin-space → bone-space).
    pub inverse_bind_pose: Mat4,
    /// Local rest transform (bind pose in parent space).
    pub rest_pose: Mat4,
}
impl Bone {
    /// Create a root bone with identity matrices.
    pub fn root(id: u32, name: impl Into<String>) -> Self {
        let identity = [
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
        Self {
            id,
            name: name.into(),
            parent: None,
            inverse_bind_pose: identity,
            rest_pose: identity,
        }
    }
}
/// A skeleton: an ordered list of bones forming a hierarchy.
#[derive(Debug, Clone, Default)]
pub struct Skeleton {
    /// All bones (index = position in array; id may differ).
    pub bones: Vec<Bone>,
}
impl Skeleton {
    /// Create an empty skeleton.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a bone.
    pub fn add_bone(&mut self, bone: Bone) -> usize {
        let idx = self.bones.len();
        self.bones.push(bone);
        idx
    }
    /// Find a bone by name.
    pub fn find_bone(&self, name: &str) -> Option<&Bone> {
        self.bones.iter().find(|b| b.name == name)
    }
    /// Compute world-space matrices from an array of local-space matrices.
    /// `local_mats` must be indexed by bone position (same ordering as `self.bones`).
    pub fn compute_world_matrices(&self, local_mats: &[Mat4]) -> Vec<Mat4> {
        let n = self.bones.len();
        let mut world = vec![[0.0_f64; 16]; n];
        let identity: Mat4 = [
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
        for (i, bone) in self.bones.iter().enumerate() {
            let local = if i < local_mats.len() {
                local_mats[i]
            } else {
                identity
            };
            world[i] = match bone.parent {
                None => local,
                Some(pid) => {
                    let pidx = self.bones.iter().position(|b| b.id == pid);
                    match pidx {
                        Some(pi) if pi < i => mat4_mul(world[pi], local),
                        _ => local,
                    }
                }
            };
        }
        world
    }
    /// Compute final skinning matrices: world * inverse_bind_pose.
    pub fn compute_skinning_matrices(&self, local_mats: &[Mat4]) -> Vec<Mat4> {
        let world = self.compute_world_matrices(local_mats);
        world
            .iter()
            .zip(self.bones.iter())
            .map(|(&w, b)| mat4_mul(w, b.inverse_bind_pose))
            .collect()
    }
}
/// A joint in an IK chain.
#[derive(Debug, Clone, Copy)]
pub struct IkJoint {
    /// Current position in world space.
    pub position: Vec3,
    /// Current orientation as quaternion.
    pub rotation: Quat,
    /// Maximum angular change per iteration (radians).
    pub angle_limit: f64,
}
impl IkJoint {
    /// Create a new IK joint.
    pub fn new(position: Vec3) -> Self {
        Self {
            position,
            rotation: [0.0, 0.0, 0.0, 1.0],
            angle_limit: std::f64::consts::PI,
        }
    }
}
/// A sorted sequence of scalar keyframes forming an animation curve.
#[derive(Debug, Clone, Default)]
pub struct ScalarCurve {
    /// Keyframes sorted by time.
    pub keys: Vec<ScalarKey>,
    /// Easing applied after interpolation.
    pub easing: Option<Easing>,
}
impl ScalarCurve {
    /// Create an empty curve.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a keyframe and keep the list sorted.
    pub fn push_key(&mut self, key: ScalarKey) {
        self.keys.push(key);
        self.keys.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Evaluate the curve at time `t`.
    pub fn sample(&self, t: f64) -> f64 {
        let n = self.keys.len();
        if n == 0 {
            return 0.0;
        }
        if t <= self.keys[0].time {
            return self.keys[0].value;
        }
        if t >= self.keys[n - 1].time {
            return self.keys[n - 1].value;
        }
        let pos = self.keys.partition_point(|k| k.time <= t);
        let i = pos - 1;
        let j = i + 1;
        let ki = &self.keys[i];
        let kj = &self.keys[j];
        let dt = kj.time - ki.time;
        let raw_t = if dt.abs() < 1e-300 {
            0.0
        } else {
            (t - ki.time) / dt
        };
        let et = match self.easing {
            Some(e) => apply_easing(e, raw_t),
            None => raw_t,
        };
        match ki.interp {
            InterpMode::Step => ki.value,
            InterpMode::Linear => ki.value + (kj.value - ki.value) * et,
            InterpMode::Hermite => cubic_hermite(
                ki.value,
                ki.out_tangent * dt,
                kj.value,
                kj.in_tangent * dt,
                et,
            ),
            InterpMode::Bezier => {
                let p1 = ki.value + ki.out_tangent * dt / 3.0;
                let p2 = kj.value - kj.in_tangent * dt / 3.0;
                cubic_bezier(ki.value, p1, p2, kj.value, et)
            }
        }
    }
}
/// A sorted sequence of quaternion keyframes.
#[derive(Debug, Clone, Default)]
pub struct QuatCurve {
    /// Keyframes sorted by time.
    pub keys: Vec<QuatKey>,
}
impl QuatCurve {
    /// Create an empty curve.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a keyframe.
    pub fn push_key(&mut self, key: QuatKey) {
        self.keys.push(key);
        self.keys.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Evaluate the curve at time `t` (SLERP or Step).
    pub fn sample(&self, t: f64) -> Quat {
        let n = self.keys.len();
        if n == 0 {
            return [0.0, 0.0, 0.0, 1.0];
        }
        if t <= self.keys[0].time {
            return self.keys[0].value;
        }
        if t >= self.keys[n - 1].time {
            return self.keys[n - 1].value;
        }
        let pos = self.keys.partition_point(|k| k.time <= t);
        let i = pos - 1;
        let j = i + 1;
        let ki = &self.keys[i];
        let kj = &self.keys[j];
        let dt = kj.time - ki.time;
        let alpha = if dt.abs() < 1e-300 {
            0.0
        } else {
            (t - ki.time) / dt
        };
        match ki.interp {
            InterpMode::Step => ki.value,
            _ => quat_slerp(ki.value, kj.value, alpha),
        }
    }
}
/// Per-node animation data (position + rotation + scale curves).
#[derive(Debug, Clone, Default)]
pub struct NodeAnimation {
    /// Node/bone identifier.
    pub node_id: u32,
    /// Position (translation) curve.
    pub position: Vec3Curve,
    /// Rotation curve.
    pub rotation: QuatCurve,
    /// Scale curve.
    pub scale: Vec3Curve,
}
impl NodeAnimation {
    /// Create node animation for the given node ID.
    pub fn new(node_id: u32) -> Self {
        Self {
            node_id,
            ..Self::default()
        }
    }
    /// Sample a TRS matrix at time `t`.
    pub fn sample_trs(&self, t: f64) -> Mat4 {
        let pos = self.position.sample(t);
        let rot = self.rotation.sample(t);
        let scl = self.scale.sample(t);
        let scl_safe = if vec3_len(scl) < 1e-300 {
            [1.0; 3]
        } else {
            scl
        };
        trs_matrix(pos, rot, scl_safe)
    }
}
/// A single scalar keyframe.
#[derive(Debug, Clone, Copy)]
pub struct ScalarKey {
    /// Time of the key (seconds).
    pub time: f64,
    /// Value at this key.
    pub value: f64,
    /// In-tangent for Hermite/Bézier interpolation.
    pub in_tangent: f64,
    /// Out-tangent for Hermite/Bézier interpolation.
    pub out_tangent: f64,
    /// Interpolation mode to use from this key to the next.
    pub interp: InterpMode,
}
impl ScalarKey {
    /// Construct a linear scalar keyframe.
    pub fn linear(time: f64, value: f64) -> Self {
        Self {
            time,
            value,
            in_tangent: 0.0,
            out_tangent: 0.0,
            interp: InterpMode::Linear,
        }
    }
    /// Construct a Hermite scalar keyframe.
    pub fn hermite(time: f64, value: f64, in_t: f64, out_t: f64) -> Self {
        Self {
            time,
            value,
            in_tangent: in_t,
            out_tangent: out_t,
            interp: InterpMode::Hermite,
        }
    }
}
/// Skinning weight entry: bone index + influence weight.
#[derive(Debug, Clone, Copy)]
pub struct SkinWeight {
    /// Bone index in the skeleton.
    pub bone_index: usize,
    /// Influence weight (should sum to 1 across all weights for a vertex).
    pub weight: f64,
}
impl SkinWeight {
    /// Create a new skin weight.
    pub fn new(bone_index: usize, weight: f64) -> Self {
        Self { bone_index, weight }
    }
}
/// An animation clip containing channel data for multiple nodes.
#[derive(Debug, Clone)]
pub struct AnimationClip {
    /// Clip name.
    pub name: String,
    /// Duration in seconds.
    pub duration: f64,
    /// Frames per second (informational).
    pub fps: f64,
    /// Whether this clip loops.
    pub looping: bool,
    /// Per-node animation channels.
    pub nodes: Vec<NodeAnimation>,
}
impl AnimationClip {
    /// Create a new clip.
    pub fn new(name: impl Into<String>, duration: f64, fps: f64, looping: bool) -> Self {
        Self {
            name: name.into(),
            duration,
            fps,
            looping,
            nodes: Vec::new(),
        }
    }
    /// Add a node animation channel.
    pub fn add_node(&mut self, node_anim: NodeAnimation) {
        self.nodes.push(node_anim);
    }
    /// Wrap or clamp `t` according to the loop flag.
    pub fn normalise_time(&self, t: f64) -> f64 {
        if self.duration <= 0.0 {
            return 0.0;
        }
        if self.looping {
            t.rem_euclid(self.duration)
        } else {
            t.clamp(0.0, self.duration)
        }
    }
    /// Sample all node TRS matrices at time `t`.
    /// Returns a map from node_id → Mat4.
    pub fn sample(&self, t: f64) -> HashMap<u32, Mat4> {
        let nt = self.normalise_time(t);
        self.nodes
            .iter()
            .map(|na| (na.node_id, na.sample_trs(nt)))
            .collect()
    }
}
/// A transition between two states.
#[derive(Debug, Clone)]
pub struct Transition {
    /// Source state name.
    pub from: String,
    /// Target state name.
    pub to: String,
    /// Blend duration (seconds) for cross-fading.
    pub blend_duration: f64,
    /// Optional trigger name that activates this transition.
    pub trigger: Option<String>,
}
impl Transition {
    /// Create a triggered transition.
    pub fn triggered(
        from: impl Into<String>,
        to: impl Into<String>,
        blend_duration: f64,
        trigger: impl Into<String>,
    ) -> Self {
        Self {
            from: from.into(),
            to: to.into(),
            blend_duration,
            trigger: Some(trigger.into()),
        }
    }
    /// Create an automatic (no trigger) transition.
    pub fn automatic(from: impl Into<String>, to: impl Into<String>, blend_duration: f64) -> Self {
        Self {
            from: from.into(),
            to: to.into(),
            blend_duration,
            trigger: None,
        }
    }
}
/// Animation state machine.
#[derive(Debug, Clone, Default)]
pub struct AnimStateMachine {
    /// All registered states.
    pub states: HashMap<String, AnimState>,
    /// All registered transitions.
    pub transitions: Vec<Transition>,
    /// Currently active state name.
    pub current: Option<String>,
    /// Current playback time within the active state.
    pub time: f64,
    /// Current blend weight (0 → old state, 1 → new state) during crossfade.
    pub blend_weight: f64,
    /// Previous state name (during crossfade).
    pub prev_state: Option<String>,
    /// Previous state time (during crossfade).
    pub prev_time: f64,
    /// Active triggers.
    pub triggers: Vec<String>,
}
impl AnimStateMachine {
    /// Create an empty state machine.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a state.
    pub fn add_state(&mut self, state: AnimState) {
        self.states.insert(state.name.clone(), state);
    }
    /// Register a transition.
    pub fn add_transition(&mut self, t: Transition) {
        self.transitions.push(t);
    }
    /// Set the initial (entry) state.
    pub fn set_entry(&mut self, name: &str) {
        self.current = Some(name.to_string());
        self.time = 0.0;
    }
    /// Fire a named trigger.
    pub fn fire_trigger(&mut self, trigger: &str) {
        self.triggers.push(trigger.to_string());
    }
    /// Advance time by `dt` seconds and process trigger-based transitions.
    pub fn update(&mut self, dt: f64, clips: &HashMap<String, AnimationClip>) {
        if let Some(ref cur) = self.current.clone() {
            let speed = self.states.get(cur).map(|s| s.speed).unwrap_or(1.0);
            self.time += dt * speed;
            if let Some(state) = self.states.get(cur)
                && let Some(clip) = clips.get(&state.clip_name)
                && clip.looping
            {
                self.time = self.time.rem_euclid(clip.duration.max(1e-9));
            }
        }
        if self.blend_weight < 1.0 {
            self.blend_weight = (self.blend_weight + dt / 0.3).min(1.0);
        }
        let triggers_snapshot = self.triggers.clone();
        for trigger in &triggers_snapshot {
            if let Some(ref cur) = self.current.clone() {
                for tr in self.transitions.clone().iter() {
                    if &tr.from == cur
                        && let Some(ref t) = tr.trigger
                        && t == trigger
                    {
                        self.prev_state = self.current.clone();
                        self.prev_time = self.time;
                        self.current = Some(tr.to.clone());
                        self.time = 0.0;
                        self.blend_weight = 0.0;
                        break;
                    }
                }
            }
        }
        self.triggers.clear();
    }
    /// Evaluate the current (blended) pose.
    pub fn evaluate(&self, clips: &HashMap<String, AnimationClip>) -> HashMap<u32, Mat4> {
        let cur_pose = self
            .current
            .as_ref()
            .and_then(|cur| {
                self.states
                    .get(cur)
                    .and_then(|s| clips.get(&s.clip_name).map(|c| c.sample(self.time)))
            })
            .unwrap_or_default();
        if self.blend_weight < 1.0
            && let Some(prev_pose) = self.prev_state.as_ref().and_then(|ps| {
                self.states
                    .get(ps)
                    .and_then(|s| clips.get(&s.clip_name).map(|c| c.sample(self.prev_time)))
            })
        {
            let w = self.blend_weight;
            let mut blended = prev_pose;
            for (id, m_cur) in &cur_pose {
                let entry = blended.entry(*id).or_insert(*m_cur);
                *entry = lerp_mat4(*entry, *m_cur, w);
            }
            return blended;
        }
        cur_pose
    }
}
/// A single state in the animation state machine.
#[derive(Debug, Clone)]
pub struct AnimState {
    /// State name.
    pub name: String,
    /// Clip to play in this state.
    pub clip_name: String,
    /// Playback speed multiplier.
    pub speed: f64,
}
impl AnimState {
    /// Create a new animation state.
    pub fn new(name: impl Into<String>, clip_name: impl Into<String>, speed: f64) -> Self {
        Self {
            name: name.into(),
            clip_name: clip_name.into(),
            speed,
        }
    }
}
/// A single quaternion keyframe.
#[derive(Debug, Clone, Copy)]
pub struct QuatKey {
    /// Time of the key (seconds).
    pub time: f64,
    /// Quaternion value at this key.
    pub value: Quat,
    /// Interpolation mode (Linear → SLERP; Step → step).
    pub interp: InterpMode,
}
impl QuatKey {
    /// Construct a linear (SLERP) quaternion keyframe.
    pub fn linear(time: f64, value: Quat) -> Self {
        Self {
            time,
            value,
            interp: InterpMode::Linear,
        }
    }
}
