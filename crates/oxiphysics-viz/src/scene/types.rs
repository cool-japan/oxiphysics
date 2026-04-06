//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{NodeId, mat4_mul};
use super::types_rendering::{DrawCall, LodEntry};

/// A single node in the scene graph.
pub struct SceneNode {
    /// Unique identifier (index in the owning `Scene::nodes` vector).
    pub id: NodeId,
    /// Human-readable name.
    pub name: String,
    /// Local transform relative to the parent node (or world if root).
    pub transform: Transform,
    /// Optional parent node.
    pub parent: Option<NodeId>,
    /// Direct children.
    pub children: Vec<NodeId>,
    /// Index into a mesh array, if this node carries geometry.
    pub mesh_index: Option<usize>,
    /// Whether this node (and its subtree) should be rendered.
    pub visible: bool,
}
/// Per-instance data for GPU instanced rendering.
#[derive(Debug, Clone)]
pub struct InstanceData {
    /// World transform matrix (column-major 4×4) for each instance.
    pub transforms: Vec<[[f32; 4]; 4]>,
    /// Optional per-instance colour modifiers.
    pub colors: Vec<[f32; 4]>,
}
impl InstanceData {
    /// Create empty instance data.
    pub fn new() -> Self {
        Self {
            transforms: Vec::new(),
            colors: Vec::new(),
        }
    }
    /// Add an instance with a 4×4 transform and a colour.
    pub fn push(&mut self, transform: [[f32; 4]; 4], color: [f32; 4]) {
        self.transforms.push(transform);
        self.colors.push(color);
    }
    /// Number of instances.
    pub fn len(&self) -> usize {
        self.transforms.len()
    }
    /// Returns `true` when there are no instances.
    pub fn is_empty(&self) -> bool {
        self.transforms.is_empty()
    }
    /// Identity transform as `[[f32; 4\]; 4]`.
    pub fn identity_transform() -> [[f32; 4]; 4] {
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}
/// Light manager: manages a list of lights and supports queries.
#[derive(Default)]
pub struct LightManager {
    pub(super) lights: Vec<(String, Light)>,
}
impl LightManager {
    /// Create an empty light manager.
    pub fn new() -> Self {
        Self { lights: Vec::new() }
    }
    /// Add a named light; returns its index.
    pub fn add(&mut self, name: &str, light: Light) -> usize {
        let idx = self.lights.len();
        self.lights.push((name.to_string(), light));
        idx
    }
    /// Remove a light by index.
    pub fn remove(&mut self, idx: usize) {
        if idx < self.lights.len() {
            self.lights.remove(idx);
        }
    }
    /// Find a light by name; returns its index.
    pub fn find(&self, name: &str) -> Option<usize> {
        self.lights.iter().position(|(n, _)| n == name)
    }
    /// Get a reference to a light by index.
    pub fn get(&self, idx: usize) -> Option<&Light> {
        self.lights.get(idx).map(|(_, l)| l)
    }
    /// Get a mutable reference to a light by index.
    pub fn get_mut(&mut self, idx: usize) -> Option<&mut Light> {
        self.lights.get_mut(idx).map(|(_, l)| l)
    }
    /// Number of lights.
    pub fn len(&self) -> usize {
        self.lights.len()
    }
    /// Returns `true` when there are no lights.
    pub fn is_empty(&self) -> bool {
        self.lights.is_empty()
    }
    /// Total intensity across all lights.
    pub fn total_intensity(&self) -> f32 {
        self.lights.iter().map(|(_, l)| l.intensity()).sum()
    }
    /// Iterate over all lights.
    pub fn iter(&self) -> impl Iterator<Item = &(String, Light)> {
        self.lights.iter()
    }
}
/// Atmospheric fog model selection.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FogModel {
    /// Uniform exponential decay: factor = exp(-density * distance).
    Exponential,
    /// Squared exponential: factor = exp(-(density * distance)²).
    ExponentialSquared,
    /// Linear fog between `near` and `far` clip distances.
    Linear {
        /// Distance at which fog begins (no fog).
        near: f64,
        /// Distance at which fog is fully opaque.
        far: f64,
    },
}
#[allow(dead_code)]
impl FogModel {
    /// Compute the fog factor (0 = fully fogged, 1 = no fog) at distance `d`.
    pub fn factor(&self, d: f64, density: f64) -> f64 {
        match self {
            FogModel::Exponential => (-density * d).exp().clamp(0.0, 1.0),
            FogModel::ExponentialSquared => {
                let x = density * d;
                (-(x * x)).exp().clamp(0.0, 1.0)
            }
            FogModel::Linear { near, far } => {
                let range = (far - near).abs().max(1e-15);
                (1.0 - (d - near) / range).clamp(0.0, 1.0)
            }
        }
    }
    /// Apply fog to a surface colour `surface` at distance `d`.
    ///
    /// Returns the blended colour `surface * f + fog_color * (1 - f)`.
    pub fn apply_to_color(
        &self,
        surface: [f32; 3],
        fog_color: [f32; 3],
        d: f64,
        density: f64,
    ) -> [f32; 3] {
        let f = self.factor(d, density) as f32;
        [
            surface[0] * f + fog_color[0] * (1.0 - f),
            surface[1] * f + fog_color[1] * (1.0 - f),
            surface[2] * f + fog_color[2] * (1.0 - f),
        ]
    }
}
/// Manages LOD selection for a set of objects with configurable distance tiers.
///
/// Each object has an ID and a centre position.  The manager computes the
/// appropriate `LodConfig` entry for each object based on its distance to the
/// camera position.
#[allow(dead_code)]
pub struct LodManager {
    /// LOD configuration (shared by all managed objects).
    pub config: LodConfig,
    /// Per-object centre positions `[x, y, z]`.
    pub positions: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl LodManager {
    /// Create a new LOD manager with the given config and object count.
    pub fn new(config: LodConfig) -> Self {
        Self {
            config,
            positions: Vec::new(),
        }
    }
    /// Register an object at world-space position `pos`.  Returns the object index.
    pub fn add_object(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.positions.len();
        self.positions.push(pos);
        idx
    }
    /// Compute the LOD mesh index for object `idx` from `camera_pos`.
    pub fn lod_for(&self, idx: usize, camera_pos: [f64; 3]) -> Option<usize> {
        let pos = self.positions.get(idx)?;
        let dx = pos[0] - camera_pos[0];
        let dy = pos[1] - camera_pos[1];
        let dz = pos[2] - camera_pos[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        self.config.select_lod(dist)
    }
    /// Batch-compute LOD indices for all objects.
    ///
    /// Returns a `Vec<Option`usize`>` where entry `i` holds the mesh index
    /// for object `i`.
    pub fn compute_all_lods(&self, camera_pos: [f64; 3]) -> Vec<Option<usize>> {
        (0..self.positions.len())
            .map(|i| self.lod_for(i, camera_pos))
            .collect()
    }
    /// Number of managed objects.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Returns `true` when there are no objects.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}
/// A serialized representation of a scene node.
#[derive(Debug, Clone)]
pub struct SerializedNode {
    /// Node id.
    pub id: NodeId,
    /// Name.
    pub name: String,
    /// Parent id (None for roots).
    pub parent: Option<NodeId>,
    /// Position \[x, y, z\].
    pub position: [f64; 3],
    /// Rotation quaternion \[x, y, z, w\].
    pub rotation: [f64; 4],
    /// Scale \[x, y, z\].
    pub scale: [f64; 3],
    /// Mesh index.
    pub mesh_index: Option<usize>,
    /// Visibility.
    pub visible: bool,
}
/// HDR tone-mapping parameters for scene post-processing.
///
/// Stores the complete set of parameters needed to configure a full HDR
/// rendering pipeline: exposure, knee, and colour-space conversion.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct HdrToneMapper {
    /// Exposure value in stops (EV).  0.0 = no change; positive = brighter.
    pub ev: f32,
    /// Tone-mapping knee (soft-clip inflection point).  Typical: 0.5–1.0.
    pub knee: f32,
    /// Maximum luminance of the display; values above this are clipped.
    pub max_luminance: f32,
    /// Whether to apply sRGB gamma encoding after tone-mapping.
    pub apply_srgb_gamma: bool,
}
#[allow(dead_code)]
impl HdrToneMapper {
    /// Create a default tone-mapper with neutral settings.
    pub fn default_settings() -> Self {
        Self {
            ev: 0.0,
            knee: 0.75,
            max_luminance: 1.0,
            apply_srgb_gamma: true,
        }
    }
    /// Apply exposure compensation to a single channel value.
    #[inline]
    pub fn apply_exposure(&self, v: f32) -> f32 {
        v * (2.0_f32).powf(self.ev)
    }
    /// Soft-clip around the knee point using a smooth-step function.
    ///
    /// Values below `knee` are passed through; values above are
    /// smoothly compressed toward `max_luminance`.
    #[inline]
    pub fn soft_clip(&self, v: f32) -> f32 {
        if v <= self.knee {
            return v;
        }
        let t = (v - self.knee) / (self.max_luminance - self.knee).max(1e-6);
        let t_clamped = t.clamp(0.0, 1.0);
        let smooth = t_clamped * t_clamped * (3.0 - 2.0 * t_clamped);
        self.knee + smooth * (self.max_luminance - self.knee)
    }
    /// Apply sRGB gamma encoding (IEC 61966-2-1) to a linear value.
    #[inline]
    pub fn srgb_encode(&self, v: f32) -> f32 {
        let c = v.clamp(0.0, 1.0);
        if c <= 0.003_130_8 {
            12.92 * c
        } else {
            1.055 * c.powf(1.0 / 2.4) - 0.055
        }
    }
    /// Process a single HDR channel: exposure → soft-clip → optional sRGB.
    pub fn process_channel(&self, v: f32) -> f32 {
        let exposed = self.apply_exposure(v);
        let clipped = self.soft_clip(exposed);
        if self.apply_srgb_gamma {
            self.srgb_encode(clipped)
        } else {
            clipped
        }
    }
    /// Compute the effective gain needed to normalise a scene with the given
    /// average luminance to the standard reference of 0.18 mid-grey.
    pub fn key_gain(&self, avg_scene_luminance: f32) -> f32 {
        if avg_scene_luminance < 1e-8 {
            return 1.0;
        }
        0.18 / avg_scene_luminance
    }
}
/// A flat-array scene graph.
pub struct Scene {
    /// All nodes, indexed by `NodeId`.
    pub nodes: Vec<SceneNode>,
    /// Top-level (parentless) nodes.
    pub root_nodes: Vec<NodeId>,
}
impl Scene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root_nodes: Vec::new(),
        }
    }
    /// Add a node and return its `NodeId`.
    pub fn add_node(&mut self, name: &str, transform: Transform, parent: Option<NodeId>) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(SceneNode {
            id,
            name: name.to_string(),
            transform,
            parent,
            children: Vec::new(),
            mesh_index: None,
            visible: true,
        });
        if let Some(pid) = parent {
            self.nodes[pid].children.push(id);
        } else {
            self.root_nodes.push(id);
        }
        id
    }
    /// Compute the world-space TRS matrix for `id` by composing the parent chain.
    pub fn world_transform(&self, id: NodeId) -> [[f64; 4]; 4] {
        let node = &self.nodes[id];
        let local = node.transform.to_matrix();
        match node.parent {
            None => local,
            Some(pid) => mat4_mul(self.world_transform(pid), local),
        }
    }
    /// Set the visibility of a single node.
    pub fn set_visible(&mut self, id: NodeId, visible: bool) {
        self.nodes[id].visible = visible;
    }
    /// Collect all node IDs whose `visible` flag is `true`.
    pub fn visible_nodes(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.visible)
            .map(|n| n.id)
            .collect()
    }
    /// Remove a node and unlink it from its parent and children.
    ///
    /// Children of the removed node are re-parented to the removed node's parent (or become roots).
    pub fn remove_node(&mut self, id: NodeId) {
        let parent = self.nodes[id].parent;
        let children: Vec<NodeId> = self.nodes[id].children.clone();
        for &child in &children {
            self.nodes[child].parent = parent;
            if let Some(pid) = parent {
                self.nodes[pid].children.push(child);
            } else {
                self.root_nodes.push(child);
            }
        }
        if let Some(pid) = parent {
            self.nodes[pid].children.retain(|&c| c != id);
        } else {
            self.root_nodes.retain(|&r| r != id);
        }
        self.nodes[id].children.clear();
        self.nodes[id].parent = None;
        self.nodes[id].visible = false;
    }
    /// Find the first node whose name equals `name`.
    pub fn find_by_name(&self, name: &str) -> Option<NodeId> {
        self.nodes.iter().find(|n| n.name == name).map(|n| n.id)
    }
    /// Return the depth of `id` in the scene tree (roots have depth 0).
    pub fn depth(&self, id: NodeId) -> usize {
        match self.nodes[id].parent {
            None => 0,
            Some(pid) => 1 + self.depth(pid),
        }
    }
}
impl Scene {
    /// Test which nodes pass a given frustum (sphere test using a provided
    /// per-node radius table).
    ///
    /// Returns a list of `(node_id, camera_distance)` for visible nodes.
    pub fn cull_with_frustum(
        &self,
        frustum: &Frustum,
        camera_pos: [f64; 3],
        node_radii: &[f64],
    ) -> Vec<(NodeId, f64)> {
        let mut visible = Vec::new();
        for node in &self.nodes {
            if !node.visible {
                continue;
            }
            let radius = if node.id < node_radii.len() {
                node_radii[node.id]
            } else {
                1.0
            };
            let wt = self.world_transform(node.id);
            let pos = [wt[3][0], wt[3][1], wt[3][2]];
            if frustum.contains_sphere(pos, radius) {
                let dx = pos[0] - camera_pos[0];
                let dy = pos[1] - camera_pos[1];
                let dz = pos[2] - camera_pos[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                visible.push((node.id, dist));
            }
        }
        visible
    }
    /// Collect all descendants of a node (depth-first).
    pub fn descendants(&self, id: NodeId) -> Vec<NodeId> {
        let mut result = Vec::new();
        let mut stack = vec![id];
        while let Some(cur) = stack.pop() {
            for &child in &self.nodes[cur].children {
                result.push(child);
                stack.push(child);
            }
        }
        result
    }
    /// Set the mesh_index for a node.
    pub fn set_mesh_index(&mut self, id: NodeId, mesh_index: Option<usize>) {
        self.nodes[id].mesh_index = mesh_index;
    }
    /// Return all nodes that carry geometry (have `mesh_index` set).
    pub fn geometry_nodes(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.visible && n.mesh_index.is_some())
            .map(|n| n.id)
            .collect()
    }
}
impl Scene {
    /// Serialize the scene to a list of `SerializedNode` records.
    pub fn serialize(&self) -> Vec<SerializedNode> {
        self.nodes
            .iter()
            .map(|n| SerializedNode {
                id: n.id,
                name: n.name.clone(),
                parent: n.parent,
                position: n.transform.position,
                rotation: n.transform.rotation,
                scale: n.transform.scale,
                mesh_index: n.mesh_index,
                visible: n.visible,
            })
            .collect()
    }
    /// Deserialize a scene from a list of `SerializedNode` records.
    ///
    /// Rebuilds parent/child links and root_nodes.
    pub fn deserialize(records: &[SerializedNode]) -> Self {
        let mut scene = Scene::new();
        for rec in records {
            let node = SceneNode {
                id: rec.id,
                name: rec.name.clone(),
                transform: Transform {
                    position: rec.position,
                    rotation: rec.rotation,
                    scale: rec.scale,
                },
                parent: rec.parent,
                children: Vec::new(),
                mesh_index: rec.mesh_index,
                visible: rec.visible,
            };
            scene.nodes.push(node);
        }
        for i in 0..scene.nodes.len() {
            if let Some(pid) = scene.nodes[i].parent {
                let child_id = scene.nodes[i].id;
                scene.nodes[pid].children.push(child_id);
            } else {
                scene.root_nodes.push(scene.nodes[i].id);
            }
        }
        scene
    }
}
impl Scene {
    /// Compute a tightly-fitting orthographic shadow-map frustum for a
    /// directional light specified by `light_dir` (unit vector toward the light).
    ///
    /// The frustum encapsulates all currently-visible nodes using their
    /// world-space positions derived from `world_transform`.  The returned
    /// `[[f64; 4\]; 4]` is a column-major orthographic projection matrix whose
    /// clip volume contains the scene bounding box when viewed along `light_dir`.
    ///
    /// Returns the identity matrix when the scene has no visible nodes.
    pub fn compute_shadow_map_frustum(&self, light_dir: [f64; 3]) -> [[f64; 4]; 4] {
        let visible: Vec<NodeId> = self.visible_nodes();
        if visible.is_empty() {
            return [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }
        let len = (light_dir[0].powi(2) + light_dir[1].powi(2) + light_dir[2].powi(2))
            .sqrt()
            .max(1e-15);
        let ld = [light_dir[0] / len, light_dir[1] / len, light_dir[2] / len];
        let up_hint = if ld[1].abs() < 0.9 {
            [0.0, 1.0, 0.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let right = [
            up_hint[1] * ld[2] - up_hint[2] * ld[1],
            up_hint[2] * ld[0] - up_hint[0] * ld[2],
            up_hint[0] * ld[1] - up_hint[1] * ld[0],
        ];
        let rlen = (right[0].powi(2) + right[1].powi(2) + right[2].powi(2))
            .sqrt()
            .max(1e-15);
        let right = [right[0] / rlen, right[1] / rlen, right[2] / rlen];
        let up = [
            ld[1] * right[2] - ld[2] * right[1],
            ld[2] * right[0] - ld[0] * right[2],
            ld[0] * right[1] - ld[1] * right[0],
        ];
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut min_z = f64::INFINITY;
        let mut max_z = f64::NEG_INFINITY;
        for id in &visible {
            let wt = self.world_transform(*id);
            let p = [wt[3][0], wt[3][1], wt[3][2]];
            let lx = p[0] * right[0] + p[1] * right[1] + p[2] * right[2];
            let ly = p[0] * up[0] + p[1] * up[1] + p[2] * up[2];
            let lz = p[0] * ld[0] + p[1] * ld[1] + p[2] * ld[2];
            if lx < min_x {
                min_x = lx;
            }
            if lx > max_x {
                max_x = lx;
            }
            if ly < min_y {
                min_y = ly;
            }
            if ly > max_y {
                max_y = ly;
            }
            if lz < min_z {
                min_z = lz;
            }
            if lz > max_z {
                max_z = lz;
            }
        }
        let w = (max_x - min_x).max(1e-6);
        let h = (max_y - min_y).max(1e-6);
        let d = (max_z - min_z).max(1e-6);
        let cx = (min_x + max_x) * 0.5;
        let cy = (min_y + max_y) * 0.5;
        let cz = min_z;
        [
            [2.0 / w, 0.0, 0.0, 0.0],
            [0.0, 2.0 / h, 0.0, 0.0],
            [0.0, 0.0, -2.0 / d, 0.0],
            [-2.0 * cx / w, -2.0 * cy / h, -(max_z + cz) / d, 1.0],
        ]
    }
    /// Set the `visible` flag to `false` for objects that are entirely hidden
    /// behind a coarse occluder set.
    ///
    /// This is a CPU-side approximation.  For each node we project its
    /// world-space position through `view_proj` and test whether its normalised
    /// device coordinate depth is behind the maximum depth of any occluder in
    /// `occluder_aabbs`.  Nodes that are definitively behind all occluders have
    /// their `visible` flag set to `false`; others are left unchanged.
    ///
    /// `occluder_aabbs` is a slice of `(min, max)` axis-aligned bounding
    /// boxes (world space) that act as solid blockers.  `view_proj` is a
    /// column-major combined view-projection matrix.
    pub fn cull_invisible_objects(
        &mut self,
        view_proj: [[f64; 4]; 4],
        occluder_aabbs: &[([f64; 3], [f64; 3])],
    ) {
        for node in self.nodes.iter_mut() {
            if !node.visible {
                continue;
            }
            let wt = { node.transform.position };
            let p = wt;
            let m = view_proj;
            let cx = m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0];
            let cy = m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1];
            let cz = m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2];
            let cw = m[0][3] * p[0] + m[1][3] * p[1] + m[2][3] * p[2] + m[3][3];
            if cw.abs() < 1e-12 {
                continue;
            }
            let ndcx = cx / cw;
            let ndcy = cy / cw;
            let ndcz = cz / cw;
            for &(omin, omax) in occluder_aabbs {
                let mut occ_max_z = f64::NEG_INFINITY;
                for bx in [omin[0], omax[0]] {
                    for by in [omin[1], omax[1]] {
                        for bz in [omin[2], omax[2]] {
                            let ocw = m[0][3] * bx + m[1][3] * by + m[2][3] * bz + m[3][3];
                            if ocw.abs() < 1e-12 {
                                continue;
                            }
                            let ocz = (m[0][2] * bx + m[1][2] * by + m[2][2] * bz + m[3][2]) / ocw;
                            if ocz > occ_max_z {
                                occ_max_z = ocz;
                            }
                        }
                    }
                }
                let occluder_cx = ((m[0][0] * omin[0]
                    + m[1][0] * omin[1]
                    + m[2][0] * omin[2]
                    + m[3][0])
                    + (m[0][0] * omax[0] + m[1][0] * omax[1] + m[2][0] * omax[2] + m[3][0]))
                    * 0.5
                    / (((m[0][3] * omin[0] + m[1][3] * omin[1] + m[2][3] * omin[2] + m[3][3])
                        + (m[0][3] * omax[0] + m[1][3] * omax[1] + m[2][3] * omax[2] + m[3][3]))
                        * 0.5)
                        .max(1e-12);
                let occluder_cy = ((m[0][1] * omin[0]
                    + m[1][1] * omin[1]
                    + m[2][1] * omin[2]
                    + m[3][1])
                    + (m[0][1] * omax[0] + m[1][1] * omax[1] + m[2][1] * omax[2] + m[3][1]))
                    * 0.5
                    / (((m[0][3] * omin[0] + m[1][3] * omin[1] + m[2][3] * omin[2] + m[3][3])
                        + (m[0][3] * omax[0] + m[1][3] * omax[1] + m[2][3] * omax[2] + m[3][3]))
                        * 0.5)
                        .max(1e-12);
                let half_ew = ((m[0][0] * (omax[0] - omin[0])
                    + m[1][0] * (omax[1] - omin[1])
                    + m[2][0] * (omax[2] - omin[2]))
                    / cw.abs()
                    * 0.5)
                    .abs();
                let half_eh = ((m[0][1] * (omax[0] - omin[0])
                    + m[1][1] * (omax[1] - omin[1])
                    + m[2][1] * (omax[2] - omin[2]))
                    / cw.abs()
                    * 0.5)
                    .abs();
                if ndcz > occ_max_z
                    && (ndcx - occluder_cx).abs() < half_ew + 0.1
                    && (ndcy - occluder_cy).abs() < half_eh + 0.1
                {
                    node.visible = false;
                    break;
                }
            }
            if !(-1.0..=1.0).contains(&ndcx)
                || !(-1.0..=1.0).contains(&ndcy)
                || !(-1.0..=1.0).contains(&ndcz)
            {
                node.visible = false;
            }
        }
    }
    /// Select the LOD level for a node based on its distance from the camera.
    ///
    /// Returns the mesh index from `lod_config` that corresponds to the
    /// appropriate detail level for `camera_distance`, or `None` when the
    /// config is empty.  This method does not modify the scene; callers
    /// should use the result to override `node.mesh_index` if desired.
    pub fn compute_lod_level(lod_config: &LodConfig, camera_distance: f64) -> Option<usize> {
        lod_config.select_lod(camera_distance)
    }
}
/// Sorted list of draw calls, ready for submission to a renderer.
///
/// Opaque objects are sorted front-to-back (ascending depth) to reduce
/// overdraw; transparent objects should be added via [`RenderQueue::push_transparent`]
/// and are sorted back-to-front (descending depth).
pub struct RenderQueue {
    pub(super) opaque: Vec<DrawCall>,
    pub(super) transparent: Vec<DrawCall>,
}
impl RenderQueue {
    /// Create an empty render queue.
    pub fn new() -> Self {
        Self {
            opaque: Vec::new(),
            transparent: Vec::new(),
        }
    }
    /// Clear all queued draw calls.
    pub fn clear(&mut self) {
        self.opaque.clear();
        self.transparent.clear();
    }
    /// Push an opaque draw call.
    pub fn push_opaque(&mut self, call: DrawCall) {
        self.opaque.push(call);
    }
    /// Push a transparent draw call.
    pub fn push_transparent(&mut self, call: DrawCall) {
        self.transparent.push(call);
    }
    /// Sort opaque calls front-to-back (smallest depth first).
    pub fn sort(&mut self) {
        self.opaque.sort_by(|a, b| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        self.transparent.sort_by(|a, b| {
            b.depth
                .partial_cmp(&a.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Returns the sorted opaque draw calls.
    pub fn opaque_calls(&self) -> &[DrawCall] {
        &self.opaque
    }
    /// Returns the sorted transparent draw calls.
    pub fn transparent_calls(&self) -> &[DrawCall] {
        &self.transparent
    }
    /// Total number of queued draw calls.
    pub fn len(&self) -> usize {
        self.opaque.len() + self.transparent.len()
    }
    /// Returns `true` when the queue is empty.
    pub fn is_empty(&self) -> bool {
        self.opaque.is_empty() && self.transparent.is_empty()
    }
}
/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl Aabb {
    /// Create an AABB from two corners.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Expand the AABB to include a point.
    pub fn expand(&mut self, p: [f64; 3]) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(p[i]);
            self.max[i] = self.max[i].max(p[i]);
        }
    }
    /// Return the centre of the AABB.
    pub fn centre(&self) -> [f64; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    /// Return the half-extents of the AABB.
    pub fn half_extents(&self) -> [f64; 3] {
        [
            (self.max[0] - self.min[0]) * 0.5,
            (self.max[1] - self.min[1]) * 0.5,
            (self.max[2] - self.min[2]) * 0.5,
        ]
    }
    /// Test whether the AABB intersects a sphere frustum (simplified
    /// sphere/AABB test).
    ///
    /// Returns `true` if the sphere overlaps the AABB.
    pub fn intersects_sphere(&self, centre: [f64; 3], radius: f64) -> bool {
        let mut d2 = 0.0f64;
        for i in 0..3 {
            let v = centre[i];
            if v < self.min[i] {
                d2 += (self.min[i] - v).powi(2);
            } else if v > self.max[i] {
                d2 += (v - self.max[i]).powi(2);
            }
        }
        d2 <= radius * radius
    }
    /// Merge two AABBs.
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }
}
/// A shadow map: a 2D depth texture generated by rendering the scene from a
/// light's perspective.  Used for PCF (Percentage Closer Filtering) shadow
/// testing.
#[allow(dead_code)]
pub struct ShadowMap {
    /// Width in texels.
    pub width: usize,
    /// Height in texels.
    pub height: usize,
    /// Depth values in `[0, 1]` stored row-major.
    pub depth: Vec<f32>,
    /// Light-space view-projection matrix (column-major).
    pub light_vp: [[f64; 4]; 4],
}
#[allow(dead_code)]
impl ShadowMap {
    /// Create an empty shadow map filled with 1.0 (maximum depth).
    pub fn new(width: usize, height: usize, light_vp: [[f64; 4]; 4]) -> Self {
        Self {
            width,
            height,
            depth: vec![1.0; width * height],
            light_vp,
        }
    }
    /// Read the depth at texel `(x, y)`.  Out-of-bounds returns 1.0.
    pub fn get_depth(&self, x: usize, y: usize) -> f32 {
        if x >= self.width || y >= self.height {
            return 1.0;
        }
        self.depth[y * self.width + x]
    }
    /// Write the depth at texel `(x, y)`.
    pub fn set_depth(&mut self, x: usize, y: usize, d: f32) {
        if x < self.width && y < self.height {
            self.depth[y * self.width + x] = d;
        }
    }
    /// Project a world-space point into shadow-map UV space.
    ///
    /// Returns `(u, v, depth)` where `u, v ∈ [0, 1]` and depth is the NDC z
    /// value.  Returns `None` when the point is behind the light.
    pub fn project(&self, world_pos: [f64; 3]) -> Option<(f32, f32, f32)> {
        let m = self.light_vp;
        let [px, py, pz] = world_pos;
        let cx = m[0][0] * px + m[1][0] * py + m[2][0] * pz + m[3][0];
        let cy = m[0][1] * px + m[1][1] * py + m[2][1] * pz + m[3][1];
        let cz = m[0][2] * px + m[1][2] * py + m[2][2] * pz + m[3][2];
        let cw = m[0][3] * px + m[1][3] * py + m[2][3] * pz + m[3][3];
        if cw.abs() < 1e-12 {
            return None;
        }
        let ndcx = (cx / cw * 0.5 + 0.5) as f32;
        let ndcy = (cy / cw * 0.5 + 0.5) as f32;
        let ndcz = (cz / cw) as f32;
        Some((ndcx.clamp(0.0, 1.0), ndcy.clamp(0.0, 1.0), ndcz))
    }
    /// Test if `world_pos` is in shadow using a 1-tap depth comparison.
    ///
    /// Returns `true` (in shadow) when the point's depth exceeds the
    /// stored shadow-map depth at the corresponding texel plus a bias.
    pub fn is_in_shadow(&self, world_pos: [f64; 3], bias: f32) -> bool {
        let Some((u, v, depth)) = self.project(world_pos) else {
            return true;
        };
        let tx = (u * (self.width - 1) as f32).round() as usize;
        let ty = (v * (self.height - 1) as f32).round() as usize;
        depth > self.get_depth(tx, ty) + bias
    }
    /// PCF 3×3 soft shadow: average 9 depth tests around the texel.
    pub fn pcf_shadow(&self, world_pos: [f64; 3], bias: f32) -> f32 {
        let Some((u, v, depth)) = self.project(world_pos) else {
            return 1.0;
        };
        let tx = (u * (self.width - 1) as f32).round() as isize;
        let ty = (v * (self.height - 1) as f32).round() as isize;
        let mut shadow_sum = 0.0_f32;
        let mut count = 0.0_f32;
        for dy in -1_isize..=1 {
            for dx in -1_isize..=1 {
                let sx = (tx + dx).clamp(0, self.width as isize - 1) as usize;
                let sy = (ty + dy).clamp(0, self.height as isize - 1) as usize;
                if depth > self.get_depth(sx, sy) + bias {
                    shadow_sum += 1.0;
                }
                count += 1.0;
            }
        }
        shadow_sum / count.max(1.0)
    }
}
/// An extended scene camera that supports orbit mode, fly mode, and camera shake.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SceneCamera {
    /// Eye position in world space.
    pub position: [f64; 3],
    /// Look-at target in world space.
    pub target: [f64; 3],
    /// World-space up vector.
    pub up: [f64; 3],
    /// Projection parameters.
    pub projection: SceneCameraProjection,
    /// Current shake offset applied to the eye position.
    pub shake_offset: [f64; 3],
    /// Current shake magnitude (decays over time).
    pub shake_magnitude: f64,
    /// Shake decay rate (fraction remaining per second, e.g. 0.9 = 90% after 1 s).
    pub shake_decay: f64,
}
impl SceneCamera {
    /// Create a new camera at `position` looking at `target`, with default projection.
    pub fn new(position: [f64; 3], target: [f64; 3]) -> Self {
        Self {
            position,
            target,
            up: [0.0, 1.0, 0.0],
            projection: SceneCameraProjection::default(),
            shake_offset: [0.0; 3],
            shake_magnitude: 0.0,
            shake_decay: 0.9,
        }
    }
    fn forward(&self) -> [f64; 3] {
        let [px, py, pz] = self.position;
        let [tx, ty, tz] = self.target;
        let fx = tx - px;
        let fy = ty - py;
        let fz = tz - pz;
        let len = (fx * fx + fy * fy + fz * fz).sqrt().max(1e-15);
        [fx / len, fy / len, fz / len]
    }
    fn right(&self) -> [f64; 3] {
        let [fx, fy, fz] = self.forward();
        let [ux, uy, uz] = self.up;
        let rx = fy * uz - fz * uy;
        let ry = fz * ux - fx * uz;
        let rz = fx * uy - fy * ux;
        let len = (rx * rx + ry * ry + rz * rz).sqrt().max(1e-15);
        [rx / len, ry / len, rz / len]
    }
    /// Orbit the camera around the target by `delta_yaw` (horizontal) and
    /// `delta_pitch` (vertical) radians.
    ///
    /// The camera maintains its distance to the target.
    pub fn orbit(&mut self, delta_yaw: f64, delta_pitch: f64) {
        let [tx, ty, tz] = self.target;
        let [px, py, pz] = self.position;
        let mut dx = px - tx;
        let mut dy = py - ty;
        let mut dz = pz - tz;
        let radius = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-15);
        let theta = dz.atan2(dx) + delta_yaw;
        let phi = (dy / radius).clamp(-1.0, 1.0).asin() + delta_pitch;
        let phi = phi.clamp(
            -std::f64::consts::FRAC_PI_2 + 0.01,
            std::f64::consts::FRAC_PI_2 - 0.01,
        );
        dx = radius * phi.cos() * theta.cos();
        dy = radius * phi.sin();
        dz = radius * phi.cos() * theta.sin();
        self.position = [tx + dx, ty + dy, tz + dz];
    }
    /// Zoom the camera toward/away from the target by `delta` world units.
    pub fn zoom(&mut self, delta: f64) {
        let [tx, ty, tz] = self.target;
        let [px, py, pz] = self.position;
        let dx = px - tx;
        let dy = py - ty;
        let dz = pz - tz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-15);
        let new_dist = (dist - delta).max(0.01);
        let scale = new_dist / dist;
        self.position = [tx + dx * scale, ty + dy * scale, tz + dz * scale];
    }
    /// Move the camera in fly mode.
    ///
    /// - `forward`: displacement along the view direction.
    /// - `right`: displacement along the camera-right axis.
    /// - `up`: displacement along the world-up axis.
    pub fn fly(&mut self, forward: f64, right: f64, up_amount: f64) {
        let [fx, fy, fz] = self.forward();
        let [rx, ry, rz] = self.right();
        let [ux, uy, uz] = self.up;
        let dx = fx * forward + rx * right + ux * up_amount;
        let dy = fy * forward + ry * right + uy * up_amount;
        let dz = fz * forward + rz * right + uz * up_amount;
        self.position[0] += dx;
        self.position[1] += dy;
        self.position[2] += dz;
        self.target[0] += dx;
        self.target[1] += dy;
        self.target[2] += dz;
    }
    /// Rotate the view direction (yaw/pitch) without moving the eye position (fly look).
    pub fn look(&mut self, delta_yaw: f64, delta_pitch: f64) {
        let [px, py, pz] = self.position;
        let [tx, ty, tz] = self.target;
        let mut fx = tx - px;
        let mut fy = ty - py;
        let mut fz = tz - pz;
        let dist = (fx * fx + fy * fy + fz * fz).sqrt().max(1e-15);
        let theta = fz.atan2(fx) + delta_yaw;
        let phi = (fy / dist).clamp(-1.0, 1.0).asin() + delta_pitch;
        let phi = phi.clamp(
            -std::f64::consts::FRAC_PI_2 + 0.01,
            std::f64::consts::FRAC_PI_2 - 0.01,
        );
        fx = dist * phi.cos() * theta.cos();
        fy = dist * phi.sin();
        fz = dist * phi.cos() * theta.sin();
        self.target = [px + fx, py + fy, pz + fz];
    }
    /// Apply a shake impulse with the given initial magnitude and decay rate.
    ///
    /// The shake offset is randomised each call to `update_shake`.
    pub fn apply_shake(&mut self, magnitude: f64, decay: f64) {
        self.shake_magnitude = magnitude;
        self.shake_decay = decay.clamp(0.0, 1.0);
    }
    /// Update the shake effect, advancing by `dt` seconds.
    ///
    /// Uses a simple deterministic pseudo-random offset (no external rand dep).
    pub fn update_shake(&mut self, dt: f64) {
        if self.shake_magnitude < 1e-6 {
            self.shake_offset = [0.0; 3];
            return;
        }
        let seed = (self.shake_magnitude * 1_000_007.0) as u64;
        let xr = ((seed.wrapping_mul(6364136223846793005).wrapping_add(1)) >> 33) as f64
            / u32::MAX as f64;
        let yr = ((seed.wrapping_mul(2862933555777941757).wrapping_add(1)) >> 33) as f64
            / u32::MAX as f64;
        let zr = ((seed.wrapping_mul(3202034522624059733).wrapping_add(1)) >> 33) as f64
            / u32::MAX as f64;
        let m = self.shake_magnitude;
        self.shake_offset = [
            (xr * 2.0 - 1.0) * m,
            (yr * 2.0 - 1.0) * m,
            (zr * 2.0 - 1.0) * m,
        ];
        self.shake_magnitude *= self.shake_decay.powf(dt);
        if self.shake_magnitude < 1e-6 {
            self.shake_magnitude = 0.0;
        }
    }
    /// Compute the view matrix incorporating the current shake offset.
    pub fn view_matrix(&self) -> [[f64; 4]; 4] {
        let [sx, sy, sz] = self.shake_offset;
        let ex = self.position[0] + sx;
        let ey = self.position[1] + sy;
        let ez = self.position[2] + sz;
        let [tx, ty, tz] = self.target;
        let [ux, uy, uz] = self.up;
        let fx = tx - ex;
        let fy = ty - ey;
        let fz = tz - ez;
        let fl = (fx * fx + fy * fy + fz * fz).sqrt().max(1e-15);
        let (fx, fy, fz) = (fx / fl, fy / fl, fz / fl);
        let rx = fy * uz - fz * uy;
        let ry = fz * ux - fx * uz;
        let rz = fx * uy - fy * ux;
        let rl = (rx * rx + ry * ry + rz * rz).sqrt().max(1e-15);
        let (rx, ry, rz) = (rx / rl, ry / rl, rz / rl);
        let upx = ry * fz - rz * fy;
        let upy = rz * fx - rx * fz;
        let upz = rx * fy - ry * fx;
        [
            [rx, upx, -fx, 0.0],
            [ry, upy, -fy, 0.0],
            [rz, upz, -fz, 0.0],
            [
                -(rx * ex + ry * ey + rz * ez),
                -(upx * ex + upy * ey + upz * ez),
                fx * ex + fy * ey + fz * ez,
                1.0,
            ],
        ]
    }
    /// Compute the projection matrix.
    pub fn projection_matrix(&self) -> [[f64; 4]; 4] {
        match self.projection {
            SceneCameraProjection::Perspective {
                fov_y,
                aspect,
                near,
                far,
            } => {
                let f = 1.0 / (fov_y * 0.5).tan();
                let ri = 1.0 / (near - far);
                [
                    [f / aspect, 0.0, 0.0, 0.0],
                    [0.0, f, 0.0, 0.0],
                    [0.0, 0.0, (near + far) * ri, -1.0],
                    [0.0, 0.0, 2.0 * near * far * ri, 0.0],
                ]
            }
            SceneCameraProjection::Orthographic {
                half_width,
                half_height,
                near,
                far,
            } => {
                let fmn_inv = 1.0 / (far - near);
                [
                    [1.0 / half_width, 0.0, 0.0, 0.0],
                    [0.0, 1.0 / half_height, 0.0, 0.0],
                    [0.0, 0.0, -2.0 * fmn_inv, 0.0],
                    [0.0, 0.0, -(far + near) * fmn_inv, 1.0],
                ]
            }
        }
    }
    /// Combined view-projection matrix.
    pub fn view_projection_matrix(&self) -> [[f64; 4]; 4] {
        mat4_mul(self.projection_matrix(), self.view_matrix())
    }
}
/// An extended scene light with attenuation coefficients and shadow support.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum SceneLight {
    /// Infinitely distant directional light.
    Directional {
        /// Normalised direction (toward the light source).
        direction: [f64; 3],
        /// RGB colour in linear space.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
    },
    /// Omnidirectional point light with physical attenuation.
    Point {
        /// World-space position.
        position: [f64; 3],
        /// RGB colour.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
        /// Constant attenuation term.
        attenuation_constant: f64,
        /// Linear attenuation term.
        attenuation_linear: f64,
        /// Quadratic attenuation term.
        attenuation_quadratic: f64,
    },
    /// Cone-shaped spotlight with inner and outer angles.
    Spotlight {
        /// World-space position.
        position: [f64; 3],
        /// Normalised direction.
        direction: [f64; 3],
        /// RGB colour.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
        /// Inner (full-bright) cone half-angle in radians.
        inner_angle: f64,
        /// Outer (fade-to-zero) cone half-angle in radians.
        outer_angle: f64,
        /// Constant attenuation.
        attenuation_constant: f64,
        /// Linear attenuation.
        attenuation_linear: f64,
        /// Quadratic attenuation.
        attenuation_quadratic: f64,
    },
}
impl SceneLight {
    /// Construct a directional light.
    pub fn directional(direction: [f64; 3], color: [f32; 3], intensity: f32) -> Self {
        let d = direction;
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-15);
        SceneLight::Directional {
            direction: [d[0] / len, d[1] / len, d[2] / len],
            color,
            intensity,
        }
    }
    /// Construct a point light with custom attenuation terms.
    pub fn point(
        position: [f64; 3],
        color: [f32; 3],
        intensity: f32,
        attenuation_constant: f64,
        attenuation_linear: f64,
        attenuation_quadratic: f64,
    ) -> Self {
        SceneLight::Point {
            position,
            color,
            intensity,
            attenuation_constant,
            attenuation_linear,
            attenuation_quadratic,
        }
    }
    /// Construct a spotlight.
    #[allow(clippy::too_many_arguments)]
    pub fn spotlight(
        position: [f64; 3],
        direction: [f64; 3],
        color: [f32; 3],
        intensity: f32,
        inner_angle: f64,
        outer_angle: f64,
        attenuation_constant: f64,
        attenuation_linear: f64,
        attenuation_quadratic: f64,
    ) -> Self {
        let d = direction;
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt().max(1e-15);
        SceneLight::Spotlight {
            position,
            direction: [d[0] / len, d[1] / len, d[2] / len],
            color,
            intensity,
            inner_angle,
            outer_angle,
            attenuation_constant,
            attenuation_linear,
            attenuation_quadratic,
        }
    }
    /// Compute physical attenuation at a given distance.
    ///
    /// Returns `1.0` for directional lights (no attenuation).
    /// For point/spot lights: `1 / (kc + kl*d + kq*d²)`.
    pub fn attenuation(&self, distance: f64) -> f64 {
        match self {
            SceneLight::Directional { .. } => 1.0,
            SceneLight::Point {
                attenuation_constant: kc,
                attenuation_linear: kl,
                attenuation_quadratic: kq,
                ..
            } => {
                let denom = kc + kl * distance + kq * distance * distance;
                1.0 / denom.max(1e-15)
            }
            SceneLight::Spotlight {
                attenuation_constant: kc,
                attenuation_linear: kl,
                attenuation_quadratic: kq,
                ..
            } => {
                let denom = kc + kl * distance + kq * distance * distance;
                1.0 / denom.max(1e-15)
            }
        }
    }
    /// Compute spotlight cone factor for a world-space point.
    ///
    /// Returns `1.0` for directional/point lights.
    /// For spotlights, interpolates between inner and outer cone angles.
    pub fn cone_factor(&self, point: [f64; 3]) -> f64 {
        match self {
            SceneLight::Directional { .. } | SceneLight::Point { .. } => 1.0,
            SceneLight::Spotlight {
                position,
                direction,
                inner_angle,
                outer_angle,
                ..
            } => {
                let lx = point[0] - position[0];
                let ly = point[1] - position[1];
                let lz = point[2] - position[2];
                let len = (lx * lx + ly * ly + lz * lz).sqrt().max(1e-15);
                let cos_angle = (lx * direction[0] + ly * direction[1] + lz * direction[2]) / len;
                let angle = cos_angle.clamp(-1.0, 1.0).acos();
                let inner = *inner_angle;
                let outer = *outer_angle;
                if angle <= inner {
                    1.0
                } else if angle >= outer {
                    0.0
                } else {
                    ((outer - angle) / (outer - inner)).clamp(0.0, 1.0)
                }
            }
        }
    }
    /// Compute an orthographic shadow frustum for this light centered at `scene_centre`
    /// with given `scene_radius`.
    ///
    /// Returns `None` for point lights (which require cube-map shadows).
    pub fn shadow_frustum(&self, scene_radius: f64) -> Option<ShadowFrustum> {
        match self {
            SceneLight::Directional { direction, .. } => {
                let [dx, dy, dz] = *direction;
                let pos = [
                    -dx * scene_radius * 2.0,
                    -dy * scene_radius * 2.0,
                    -dz * scene_radius * 2.0,
                ];
                Some(ShadowFrustum {
                    position: pos,
                    direction: *direction,
                    half_width: scene_radius,
                    near: 0.01,
                    far: scene_radius * 4.0,
                })
            }
            SceneLight::Spotlight {
                position,
                direction,
                outer_angle,
                ..
            } => Some(ShadowFrustum {
                position: *position,
                direction: *direction,
                half_width: outer_angle.tan() * scene_radius,
                near: 0.01,
                far: scene_radius * 2.0,
            }),
            SceneLight::Point { .. } => None,
        }
    }
    /// Return the RGB colour of this light.
    pub fn color(&self) -> [f32; 3] {
        match self {
            SceneLight::Directional { color, .. } => *color,
            SceneLight::Point { color, .. } => *color,
            SceneLight::Spotlight { color, .. } => *color,
        }
    }
    /// Return the intensity of this light.
    pub fn intensity(&self) -> f32 {
        match self {
            SceneLight::Directional { intensity, .. } => *intensity,
            SceneLight::Point { intensity, .. } => *intensity,
            SceneLight::Spotlight { intensity, .. } => *intensity,
        }
    }
}
/// Projection variant for [`SceneCamera`].
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SceneCameraProjection {
    /// Perspective projection with vertical FOV, aspect ratio, and clip planes.
    Perspective {
        /// Vertical field of view in radians.
        fov_y: f64,
        /// Aspect ratio (width / height).
        aspect: f64,
        /// Near clip plane.
        near: f64,
        /// Far clip plane.
        far: f64,
    },
    /// Orthographic projection.
    Orthographic {
        /// Half-width of the viewing volume.
        half_width: f64,
        /// Half-height of the viewing volume.
        half_height: f64,
        /// Near clip plane.
        near: f64,
        /// Far clip plane.
        far: f64,
    },
}
/// A view frustum represented by 6 planes (each `[nx, ny, nz, d]` with the
/// convention that a point is inside when `n·p + d >= 0`).
pub struct Frustum {
    pub(super) planes: [[f64; 4]; 6],
}
impl Frustum {
    /// Construct a frustum from a combined view-projection matrix (column-major).
    ///
    /// Uses the Gribb/Hartmann method to extract the 6 frustum planes.
    pub fn from_view_projection(m: [[f64; 4]; 4]) -> Self {
        let row = |r: usize| -> [f64; 4] { [m[0][r], m[1][r], m[2][r], m[3][r]] };
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);
        let add = |a: [f64; 4], b: [f64; 4]| -> [f64; 4] {
            [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
        };
        let sub = |a: [f64; 4], b: [f64; 4]| -> [f64; 4] {
            [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
        };
        let planes = [
            add(r3, r0),
            sub(r3, r0),
            add(r3, r1),
            sub(r3, r1),
            add(r3, r2),
            sub(r3, r2),
        ];
        let norm_plane = |p: [f64; 4]| -> [f64; 4] {
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if len < 1e-20 {
                return p;
            }
            [p[0] / len, p[1] / len, p[2] / len, p[3] / len]
        };
        Self {
            planes: planes.map(norm_plane),
        }
    }
    /// Test whether a point is inside (or on the boundary of) the frustum.
    pub fn contains_point(&self, p: [f64; 3]) -> bool {
        for plane in &self.planes {
            if plane[0] * p[0] + plane[1] * p[1] + plane[2] * p[2] + plane[3] < 0.0 {
                return false;
            }
        }
        true
    }
    /// Test whether a sphere is (possibly) inside the frustum.
    ///
    /// Returns `false` only if the sphere is entirely outside at least one plane.
    pub fn contains_sphere(&self, centre: [f64; 3], radius: f64) -> bool {
        for plane in &self.planes {
            let dist =
                plane[0] * centre[0] + plane[1] * centre[1] + plane[2] * centre[2] + plane[3];
            if dist < -radius {
                return false;
            }
        }
        true
    }
}
/// PBR metallic-roughness material.
#[derive(Debug, Clone)]
pub struct Material {
    /// Base colour (albedo) in linear RGB + alpha.
    pub albedo: [f32; 4],
    /// Perceptual roughness `[0, 1]`; 0 = mirror, 1 = fully diffuse.
    pub roughness: f32,
    /// Metalness `[0, 1]`; 0 = dielectric, 1 = metal.
    pub metalness: f32,
    /// Emissive colour multiplier (HDR allowed).
    pub emissive: [f32; 3],
    /// Human-readable name.
    pub name: String,
}
impl Material {
    /// A plain white PBR material with no emission.
    pub fn default_white() -> Self {
        Self {
            albedo: [1.0, 1.0, 1.0, 1.0],
            roughness: 0.5,
            metalness: 0.0,
            emissive: [0.0, 0.0, 0.0],
            name: "default_white".to_string(),
        }
    }
    /// A fully metallic material with the given albedo colour.
    pub fn metallic(albedo: [f32; 4]) -> Self {
        Self {
            albedo,
            roughness: 0.1,
            metalness: 1.0,
            emissive: [0.0, 0.0, 0.0],
            name: "metallic".to_string(),
        }
    }
    /// A purely emissive material.
    pub fn emissive(color: [f32; 3]) -> Self {
        Self {
            albedo: [0.0, 0.0, 0.0, 1.0],
            roughness: 1.0,
            metalness: 0.0,
            emissive: color,
            name: "emissive".to_string(),
        }
    }
    /// Returns `true` if this material has any non-zero emission.
    pub fn is_emissive(&self) -> bool {
        self.emissive[0] > 0.0 || self.emissive[1] > 0.0 || self.emissive[2] > 0.0
    }
}
/// Scene camera with view and projection matrices.
#[derive(Debug, Clone)]
pub struct Camera {
    /// Camera position in world space.
    pub position: [f64; 3],
    /// Point the camera looks at.
    pub target: [f64; 3],
    /// World-space up vector.
    pub up: [f64; 3],
    /// Projection parameters.
    pub projection: Projection,
}
impl Camera {
    /// Create a perspective camera at `position` looking at `target`.
    pub fn perspective(
        position: [f64; 3],
        target: [f64; 3],
        fov_y: f64,
        aspect: f64,
        near: f64,
        far: f64,
    ) -> Self {
        Self {
            position,
            target,
            up: [0.0, 1.0, 0.0],
            projection: Projection::Perspective {
                fov_y,
                aspect,
                near,
                far,
            },
        }
    }
    /// Create an orthographic camera at `position` looking at `target`.
    pub fn orthographic(
        position: [f64; 3],
        target: [f64; 3],
        half_width: f64,
        half_height: f64,
        near: f64,
        far: f64,
    ) -> Self {
        Self {
            position,
            target,
            up: [0.0, 1.0, 0.0],
            projection: Projection::Orthographic {
                half_width,
                half_height,
                near,
                far,
            },
        }
    }
    /// Compute the view matrix (column-major `[[f64; 4\]; 4]`).
    ///
    /// Uses a right-handed look-at convention.
    pub fn view_matrix(&self) -> [[f64; 4]; 4] {
        let [ex, ey, ez] = self.position;
        let [tx, ty, tz] = self.target;
        let [ux, uy, uz] = self.up;
        let fx = tx - ex;
        let fy = ty - ey;
        let fz = tz - ez;
        let fl = (fx * fx + fy * fy + fz * fz).sqrt().max(1e-15);
        let (fx, fy, fz) = (fx / fl, fy / fl, fz / fl);
        let rx = fy * uz - fz * uy;
        let ry = fz * ux - fx * uz;
        let rz = fx * uy - fy * ux;
        let rl = (rx * rx + ry * ry + rz * rz).sqrt().max(1e-15);
        let (rx, ry, rz) = (rx / rl, ry / rl, rz / rl);
        let upx = ry * fz - rz * fy;
        let upy = rz * fx - rx * fz;
        let upz = rx * fy - ry * fx;
        [
            [rx, upx, -fx, 0.0],
            [ry, upy, -fy, 0.0],
            [rz, upz, -fz, 0.0],
            [
                -(rx * ex + ry * ey + rz * ez),
                -(upx * ex + upy * ey + upz * ez),
                fx * ex + fy * ey + fz * ez,
                1.0,
            ],
        ]
    }
    /// Compute the projection matrix (column-major `[[f64; 4\]; 4]`).
    pub fn projection_matrix(&self) -> [[f64; 4]; 4] {
        match self.projection {
            Projection::Perspective {
                fov_y,
                aspect,
                near,
                far,
            } => {
                let f = 1.0 / (fov_y * 0.5).tan();
                let range_inv = 1.0 / (near - far);
                [
                    [f / aspect, 0.0, 0.0, 0.0],
                    [0.0, f, 0.0, 0.0],
                    [0.0, 0.0, (near + far) * range_inv, -1.0],
                    [0.0, 0.0, 2.0 * near * far * range_inv, 0.0],
                ]
            }
            Projection::Orthographic {
                half_width,
                half_height,
                near,
                far,
            } => {
                let rml = 2.0 * half_width;
                let tmb = 2.0 * half_height;
                let fmn = far - near;
                [
                    [2.0 / rml, 0.0, 0.0, 0.0],
                    [0.0, 2.0 / tmb, 0.0, 0.0],
                    [0.0, 0.0, -2.0 / fmn, 0.0],
                    [0.0, 0.0, -(far + near) / fmn, 1.0],
                ]
            }
        }
    }
}
/// Scene light source.
#[derive(Debug, Clone)]
pub enum Light {
    /// Infinitely distant directional light.
    Directional {
        /// Unit direction vector pointing *towards* the light.
        direction: [f64; 3],
        /// Linear-light colour/intensity (HDR allowed).
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
    },
    /// Omnidirectional point light.
    Point {
        /// World-space position.
        position: [f64; 3],
        /// Linear-light colour/intensity.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
        /// Range past which the light has no effect (0 = infinite).
        range: f32,
    },
    /// Cone-shaped spot light.
    Spot {
        /// World-space position.
        position: [f64; 3],
        /// Unit direction the spot points towards.
        direction: [f64; 3],
        /// Linear-light colour/intensity.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
        /// Inner cone half-angle in radians.
        inner_angle: f32,
        /// Outer cone half-angle in radians.
        outer_angle: f32,
        /// Effective range.
        range: f32,
    },
    /// Rectangular area light.
    Area {
        /// World-space centre of the rectangle.
        position: [f64; 3],
        /// Unit normal of the emitting face.
        normal: [f64; 3],
        /// Half-width and half-height of the rectangle.
        half_extents: [f32; 2],
        /// Linear-light colour/intensity.
        color: [f32; 3],
        /// Intensity multiplier.
        intensity: f32,
    },
}
impl Light {
    /// Return the colour (`[f32; 3]`) of this light.
    pub fn color(&self) -> [f32; 3] {
        match self {
            Light::Directional { color, .. } => *color,
            Light::Point { color, .. } => *color,
            Light::Spot { color, .. } => *color,
            Light::Area { color, .. } => *color,
        }
    }
    /// Return the intensity of this light.
    pub fn intensity(&self) -> f32 {
        match self {
            Light::Directional { intensity, .. } => *intensity,
            Light::Point { intensity, .. } => *intensity,
            Light::Spot { intensity, .. } => *intensity,
            Light::Area { intensity, .. } => *intensity,
        }
    }
}
/// Global scene environment settings: sky colour, ambient light, and fog.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SceneEnvironment {
    /// Sky (background) colour in linear RGB.
    pub sky_color: [f32; 3],
    /// Ambient light intensity multiplier.
    pub ambient_intensity: f32,
    /// Exponential fog density coefficient `λ`.
    ///
    /// Fog factor = `exp(-λ * d)`. Use 0 for no fog.
    pub fog_density: f32,
    /// Fog colour blended toward at high distances.
    pub fog_color: [f32; 3],
}
impl SceneEnvironment {
    /// Clear sky (light blue) with low ambient and no fog.
    pub fn clear_sky() -> Self {
        Self {
            sky_color: [0.5, 0.7, 1.0],
            ambient_intensity: 0.15,
            fog_density: 0.0,
            fog_color: [0.7, 0.8, 0.9],
        }
    }
    /// Overcast (grey) environment with moderate ambient.
    pub fn overcast() -> Self {
        Self {
            sky_color: [0.6, 0.6, 0.65],
            ambient_intensity: 0.4,
            fog_density: 0.01,
            fog_color: [0.6, 0.6, 0.65],
        }
    }
    /// Compute the exponential fog factor at distance `d`.
    ///
    /// Returns 1.0 (no fog) when `d == 0`, decays toward 0 with distance.
    pub fn fog_factor(&self, d: f64) -> f64 {
        let lambda = self.fog_density as f64;
        if lambda < 1e-20 {
            return 1.0;
        }
        (-lambda * d).exp()
    }
    /// Blend a surface colour with the fog colour based on distance.
    pub fn apply_fog(&self, surface: [f32; 3], distance: f64) -> [f32; 3] {
        let t = (1.0 - self.fog_factor(distance)) as f32;
        let [sr, sg, sb] = surface;
        let [fr, fg, fb] = self.fog_color;
        [
            sr * (1.0 - t) + fr * t,
            sg * (1.0 - t) + fg * t,
            sb * (1.0 - t) + fb * t,
        ]
    }
}
/// Cached frustum culler that stores per-object visibility as a `Vec`bool`.
///
/// After calling `cull()`, per-object visibility can be queried in O(1) via
/// `is_visible()`.  This avoids re-running the frustum test every time a
/// system queries visibility.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FrustumCuller {
    /// Visibility flags: `visibility\[i\]` is `true` when object `i` is visible.
    pub(super) visibility: Vec<bool>,
}
impl FrustumCuller {
    /// Create a culler pre-allocated for `capacity` objects (all initially invisible).
    pub fn new(capacity: usize) -> Self {
        Self {
            visibility: vec![false; capacity],
        }
    }
    /// Run the frustum test for all objects and cache the result.
    ///
    /// `centres\[i\]` and `radii\[i\]` describe the bounding sphere of object `i`.
    /// `centres` and `radii` must have the same length; the culler is resized
    /// to match.
    pub fn cull(&mut self, frustum: &Frustum, centres: &[[f64; 3]], radii: &[f64]) {
        let n = centres.len().min(radii.len());
        self.visibility.resize(n, false);
        for i in 0..n {
            self.visibility[i] = frustum.contains_sphere(centres[i], radii[i]);
        }
    }
    /// Returns `true` when object `index` passed the last `cull()` call.
    ///
    /// Out-of-bounds indices always return `false`.
    pub fn is_visible(&self, index: usize) -> bool {
        self.visibility.get(index).copied().unwrap_or(false)
    }
    /// Return an iterator over the indices of all visible objects.
    pub fn visible_indices(&self) -> Vec<usize> {
        self.visibility
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| if v { Some(i) } else { None })
            .collect()
    }
    /// Number of objects managed by this culler.
    pub fn len(&self) -> usize {
        self.visibility.len()
    }
    /// Returns `true` if there are no managed objects.
    pub fn is_empty(&self) -> bool {
        self.visibility.is_empty()
    }
    /// Count of currently visible objects.
    pub fn visible_count(&self) -> usize {
        self.visibility.iter().filter(|&&v| v).count()
    }
}
/// Affine transform composed of translation, rotation (unit quaternion), and scale.
pub struct Transform {
    /// World translation of this node relative to its parent.
    pub position: [f64; 3],
    /// Rotation represented as a unit quaternion `\[x, y, z, w\]`.
    pub rotation: [f64; 4],
    /// Non-uniform scale.
    pub scale: [f64; 3],
}
impl Transform {
    /// Identity transform: no translation, no rotation, unit scale.
    pub fn identity() -> Self {
        Self {
            position: [0.0, 0.0, 0.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
        }
    }
    /// Convert to a column-major 4×4 TRS matrix.
    pub fn to_matrix(&self) -> [[f64; 4]; 4] {
        let [x, y, z, w] = self.rotation;
        let [sx, sy, sz] = self.scale;
        let [tx, ty, tz] = self.position;
        let r00 = 1.0 - 2.0 * (y * y + z * z);
        let r10 = 2.0 * (x * y + z * w);
        let r20 = 2.0 * (x * z - y * w);
        let r01 = 2.0 * (x * y - z * w);
        let r11 = 1.0 - 2.0 * (x * x + z * z);
        let r21 = 2.0 * (y * z + x * w);
        let r02 = 2.0 * (x * z + y * w);
        let r12 = 2.0 * (y * z - x * w);
        let r22 = 1.0 - 2.0 * (x * x + y * y);
        [
            [r00 * sx, r10 * sx, r20 * sx, 0.0],
            [r01 * sy, r11 * sy, r21 * sy, 0.0],
            [r02 * sz, r12 * sz, r22 * sz, 0.0],
            [tx, ty, tz, 1.0],
        ]
    }
}
/// Node visibility state with LOD and culling.
#[derive(Debug, Clone)]
pub struct NodeRenderState {
    /// Whether the node passed the frustum culling test.
    pub in_frustum: bool,
    /// Selected LOD mesh index (None = use node's default mesh_index).
    pub active_lod: Option<usize>,
    /// Distance from the camera to this node.
    pub camera_distance: f64,
}
/// Shadow-map frustum parameters (for shadow rendering).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ShadowFrustum {
    /// World-space position of the shadow camera.
    pub position: [f64; 3],
    /// Direction the shadow camera looks.
    pub direction: [f64; 3],
    /// Half-width of the ortho shadow frustum.
    pub half_width: f64,
    /// Near clip of shadow frustum.
    pub near: f64,
    /// Far clip of shadow frustum.
    pub far: f64,
}
/// Projection type for a [`Camera`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Projection {
    /// Perspective projection.
    Perspective {
        /// Vertical field of view in radians.
        fov_y: f64,
        /// Aspect ratio (width / height).
        aspect: f64,
        /// Near clip plane distance.
        near: f64,
        /// Far clip plane distance.
        far: f64,
    },
    /// Orthographic projection.
    Orthographic {
        /// Half-width of the viewing volume.
        half_width: f64,
        /// Half-height of the viewing volume.
        half_height: f64,
        /// Near clip plane distance.
        near: f64,
        /// Far clip plane distance.
        far: f64,
    },
}
/// LOD configuration attached to a scene node.
#[derive(Debug, Clone)]
pub struct LodConfig {
    /// LOD entries ordered from nearest to farthest (ascending threshold).
    pub entries: Vec<LodEntry>,
}
impl LodConfig {
    /// Create a new empty LOD config.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }
    /// Add a LOD level. Entries should be added in order of increasing distance.
    pub fn add_level(&mut self, distance_threshold: f64, mesh_index: usize) {
        self.entries.push(LodEntry {
            distance_threshold,
            mesh_index,
        });
    }
    /// Select the appropriate mesh index for a given camera distance.
    ///
    /// Returns the mesh index of the highest-detail LOD that is within range.
    /// If the distance exceeds all thresholds, returns the last entry's mesh.
    pub fn select_lod(&self, distance: f64) -> Option<usize> {
        if self.entries.is_empty() {
            return None;
        }
        for entry in &self.entries {
            if distance <= entry.distance_threshold {
                return Some(entry.mesh_index);
            }
        }
        self.entries.last().map(|e| e.mesh_index)
    }
}
