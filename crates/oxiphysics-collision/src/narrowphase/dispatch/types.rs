//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::type_complexity)]
use crate::types::{CollisionPair, Contact, ContactManifold};
use oxiphysics_core::Transform;
use oxiphysics_geometry::Shape;
use std::collections::HashMap;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{NarrowPhaseFn, gjk_epa, shape_type_ordinal, try_specialized};

/// Identifies a shape category for dispatch table lookups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(dead_code)]
pub enum ShapeType {
    /// Sphere primitive.
    Sphere,
    /// Axis-aligned box primitive.
    Box,
    /// Capsule primitive.
    Capsule,
    /// Cylinder primitive.
    Cylinder,
    /// Cone primitive.
    Cone,
    /// Convex hull.
    ConvexHull,
    /// Triangle mesh.
    TriangleMesh,
    /// Compound shape.
    Compound,
    /// Height field.
    HeightField,
    /// Plane (infinite half-space).
    Plane,
    /// Custom user-defined shape.
    Custom(u32),
}
/// Configuration for the narrow-phase dispatcher.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DispatchConfig {
    /// Whether to use GJK+EPA as fallback when no specialized test is found.
    pub use_gjk_fallback: bool,
    /// Maximum GJK iterations for fallback.
    pub max_gjk_iterations: usize,
    /// Maximum EPA iterations for fallback.
    pub max_epa_iterations: usize,
    /// Tolerance for contact detection.
    pub contact_tolerance: f64,
}
/// A leaf triangle for mesh mid-phase.
#[derive(Debug, Clone, Copy)]
pub struct MeshTriangle {
    /// Vertex positions.
    pub vertices: [[f64; 3]; 3],
    /// Original triangle index.
    pub index: usize,
    /// Precomputed AABB.
    pub aabb: Aabb,
}
impl MeshTriangle {
    /// Create a mesh triangle and precompute its AABB.
    pub fn new(vertices: [[f64; 3]; 3], index: usize) -> Self {
        let min = [
            vertices[0][0].min(vertices[1][0]).min(vertices[2][0]),
            vertices[0][1].min(vertices[1][1]).min(vertices[2][1]),
            vertices[0][2].min(vertices[1][2]).min(vertices[2][2]),
        ];
        let max = [
            vertices[0][0].max(vertices[1][0]).max(vertices[2][0]),
            vertices[0][1].max(vertices[1][1]).max(vertices[2][1]),
            vertices[0][2].max(vertices[1][2]).max(vertices[2][2]),
        ];
        Self {
            vertices,
            index,
            aabb: Aabb::new(min, max),
        }
    }
}
/// A convex polygon contact patch in world space.
///
/// Produced by clipping an incident face against the reference face of an
/// opposing shape.  May contain up to 8 contact points.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ContactPatch {
    /// Contact points (world space).
    pub points: Vec<[f64; 3]>,
    /// Shared contact normal.
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
}
impl ContactPatch {
    /// Create an empty contact patch.
    #[allow(dead_code)]
    pub fn new(normal: [f64; 3], depth: f64) -> Self {
        Self {
            points: Vec::new(),
            normal,
            depth,
        }
    }
    /// Add a contact point.
    #[allow(dead_code)]
    pub fn add_point(&mut self, p: [f64; 3]) {
        self.points.push(p);
    }
    /// Number of contact points in this patch.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Whether the patch is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}
/// A convex piece from a convex decomposition.
#[derive(Debug)]
#[allow(dead_code)]
pub struct ConvexPiece {
    /// Vertices of the convex piece.
    pub vertices: Vec<[f64; 3]>,
    /// Shape type (always ConvexHull for pieces).
    pub shape_type: ShapeType,
    /// Local AABB for quick rejection.
    pub aabb_min: [f64; 3],
    /// Local AABB maximum.
    pub aabb_max: [f64; 3],
}
impl ConvexPiece {
    /// Create a convex piece from a vertex list.
    #[allow(dead_code)]
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        let aabb_min = vertices.iter().fold([f64::MAX; 3], |mut acc, v| {
            acc[0] = acc[0].min(v[0]);
            acc[1] = acc[1].min(v[1]);
            acc[2] = acc[2].min(v[2]);
            acc
        });
        let aabb_max = vertices.iter().fold([f64::MIN; 3], |mut acc, v| {
            acc[0] = acc[0].max(v[0]);
            acc[1] = acc[1].max(v[1]);
            acc[2] = acc[2].max(v[2]);
            acc
        });
        Self {
            vertices,
            shape_type: ShapeType::ConvexHull,
            aabb_min,
            aabb_max,
        }
    }
    /// Check AABB vs AABB overlap with another piece.
    #[allow(dead_code)]
    pub fn aabb_overlaps(&self, other: &ConvexPiece, expand: f64) -> bool {
        for i in 0..3 {
            if self.aabb_max[i] + expand < other.aabb_min[i] - expand {
                return false;
            }
            if other.aabb_max[i] + expand < self.aabb_min[i] - expand {
                return false;
            }
        }
        true
    }
}
/// Narrow-phase dispatcher with a pluggable algorithm registry.
pub struct NarrowPhaseDispatcher {
    /// Registered pair algorithms.
    pub(super) table: HashMap<DispatchKey, NarrowPhaseFn>,
    /// Configuration.
    #[allow(dead_code)]
    pub(super) config: DispatchConfig,
}
impl NarrowPhaseDispatcher {
    /// Create an empty dispatcher (no registered algorithms).
    #[allow(dead_code)]
    pub fn empty() -> Self {
        NarrowPhaseDispatcher {
            table: HashMap::new(),
            config: DispatchConfig::default(),
        }
    }
    /// Create an empty dispatcher with custom configuration.
    #[allow(dead_code)]
    pub fn with_config(config: DispatchConfig) -> Self {
        NarrowPhaseDispatcher {
            table: HashMap::new(),
            config,
        }
    }
    /// Register a collision algorithm for the given shape-type pair.
    #[allow(dead_code)]
    pub fn register_pair(&mut self, type_a: ShapeType, type_b: ShapeType, func: NarrowPhaseFn) {
        let key = DispatchKey::new(type_a, type_b);
        self.table.insert(key, func);
    }
    /// Unregister a collision algorithm for the given shape-type pair.
    /// Returns true if an entry was removed.
    #[allow(dead_code)]
    pub fn unregister_pair(&mut self, type_a: ShapeType, type_b: ShapeType) -> bool {
        let key = DispatchKey::new(type_a, type_b);
        self.table.remove(&key).is_some()
    }
    /// Check if a pair algorithm is registered for the given types.
    #[allow(dead_code)]
    pub fn has_pair(&self, type_a: ShapeType, type_b: ShapeType) -> bool {
        let key = DispatchKey::new(type_a, type_b);
        self.table.contains_key(&key)
    }
    /// Returns the number of registered pair algorithms.
    #[allow(dead_code)]
    pub fn registered_count(&self) -> usize {
        self.table.len()
    }
    /// Returns all registered dispatch keys.
    #[allow(dead_code)]
    pub fn registered_keys(&self) -> Vec<DispatchKey> {
        self.table.keys().copied().collect()
    }
    /// Dispatch a collision query between two shapes.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn dispatch(
        &self,
        shape_a: &dyn Shape,
        type_a: ShapeType,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        type_b: ShapeType,
        transform_b: &Transform,
        pair: CollisionPair,
    ) -> NarrowPhaseResult {
        let key = DispatchKey::new(type_a, type_b);
        if let Some(func) = self.table.get(&key) {
            return func(shape_a, transform_a, shape_b, transform_b, pair);
        }
        match gjk_epa(shape_a, transform_a, shape_b, transform_b, pair) {
            Some(manifold) => NarrowPhaseResult::contact(manifold),
            None => NarrowPhaseResult::separated(),
        }
    }
    /// Dispatch with swapped result: if the registered algorithm expects
    /// (A,B) but we have (B,A), swap the inputs and flip the contact normal.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn dispatch_symmetric(
        &self,
        shape_a: &dyn Shape,
        type_a: ShapeType,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        type_b: ShapeType,
        transform_b: &Transform,
        pair: CollisionPair,
    ) -> NarrowPhaseResult {
        let key = DispatchKey::new(type_a, type_b);
        if let Some(func) = self.table.get(&key) {
            let need_swap = shape_type_ordinal(type_a) > shape_type_ordinal(type_b);
            if need_swap {
                let swapped_pair = CollisionPair::new(pair.b, pair.a);
                let mut result = func(shape_b, transform_b, shape_a, transform_a, swapped_pair);
                if let Some(ref mut manifold) = result.manifold {
                    for contact in &mut manifold.contacts {
                        contact.normal = -contact.normal;
                        std::mem::swap(&mut contact.point_a, &mut contact.point_b);
                    }
                    manifold.pair = pair;
                }
                return result;
            }
            return func(shape_a, transform_a, shape_b, transform_b, pair);
        }
        match gjk_epa(shape_a, transform_a, shape_b, transform_b, pair) {
            Some(manifold) => NarrowPhaseResult::contact(manifold),
            None => NarrowPhaseResult::separated(),
        }
    }
    /// Generate contacts between two shapes (original API, no type hint needed).
    pub fn generate_contacts(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
        pair: CollisionPair,
    ) -> Option<ContactManifold> {
        if let Some(contact) = try_specialized(shape_a, transform_a, shape_b, transform_b) {
            let mut manifold = ContactManifold::new(pair);
            manifold.add_contact(contact);
            return Some(manifold);
        }
        gjk_epa(shape_a, transform_a, shape_b, transform_b, pair)
    }
    /// Batch dispatch: process multiple collision pairs at once.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn dispatch_batch(
        &self,
        pairs: &[(
            &dyn Shape,
            ShapeType,
            &Transform,
            &dyn Shape,
            ShapeType,
            &Transform,
            CollisionPair,
        )],
    ) -> Vec<NarrowPhaseResult> {
        pairs
            .iter()
            .map(|&(sa, ta_type, ta, sb, tb_type, tb, pair)| {
                self.dispatch(sa, ta_type, ta, sb, tb_type, tb, pair)
            })
            .collect()
    }
}
/// An axis-aligned bounding box for the mid-phase AABB tree.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl Aabb {
    /// Create a new AABB.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }
    /// Check AABB vs AABB overlap (closed intervals).
    pub fn overlaps(&self, other: &Aabb) -> bool {
        self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
            && self.min[2] <= other.max[2]
            && self.max[2] >= other.min[2]
    }
    /// Surface area (for SAH cost estimation).
    pub fn surface_area(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        2.0 * (dx * dy + dy * dz + dx * dz)
    }
}
/// Configuration for speculative contact generation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpeculativeConfig {
    /// Maximum distance for generating speculative contacts (collision margin).
    pub margin: f64,
    /// Velocity scale: if `v · n < -scale * margin` the contact is generated.
    pub velocity_scale: f64,
}
/// Statistics collected over a dispatch session.
#[derive(Debug, Clone, Default)]
pub struct DispatchStats {
    /// Total number of dispatch calls.
    pub total_dispatches: usize,
    /// Number of dispatches that resulted in contact.
    pub contacts_found: usize,
    /// Number of dispatches that used the GJK fallback.
    pub gjk_fallbacks: usize,
    /// Number of dispatches that used specialized algorithms.
    pub specialized_calls: usize,
    /// Number of dispatches pruned by broad phase.
    pub broad_phase_pruned: usize,
}
impl DispatchStats {
    /// Record a dispatch result.
    pub fn record(&mut self, had_contact: bool, used_gjk: bool) {
        self.total_dispatches += 1;
        if had_contact {
            self.contacts_found += 1;
        }
        if used_gjk {
            self.gjk_fallbacks += 1;
        } else {
            self.specialized_calls += 1;
        }
    }
    /// Record a broad-phase prune event.
    pub fn record_pruned(&mut self) {
        self.broad_phase_pruned += 1;
    }
    /// Contact rate: contacts_found / total_dispatches.
    pub fn contact_rate(&self) -> f64 {
        if self.total_dispatches == 0 {
            0.0
        } else {
            self.contacts_found as f64 / self.total_dispatches as f64
        }
    }
    /// Specialization rate: specialized_calls / total_dispatches.
    pub fn specialization_rate(&self) -> f64 {
        if self.total_dispatches == 0 {
            0.0
        } else {
            self.specialized_calls as f64 / self.total_dispatches as f64
        }
    }
}
/// Result of a heightfield sample at grid coordinates.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct HeightFieldSample {
    /// World-space position of the sample point.
    pub position: [f64; 3],
    /// Surface normal at this point.
    pub normal: [f64; 3],
    /// Height value.
    pub height: f64,
}
/// Simple heightfield representation for dispatch testing.
#[allow(dead_code)]
pub struct SimpleHeightField {
    /// Heights in row-major order (rows x cols).
    pub heights: Vec<f64>,
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Spacing between grid points.
    pub spacing: f64,
}
impl SimpleHeightField {
    /// Create a new heightfield.
    #[allow(dead_code)]
    pub fn new(heights: Vec<f64>, rows: usize, cols: usize, spacing: f64) -> Self {
        Self {
            heights,
            rows,
            cols,
            spacing,
        }
    }
    /// Sample height at grid coordinates (clamped).
    #[allow(dead_code)]
    pub fn height_at(&self, row: usize, col: usize) -> f64 {
        let r = row.min(self.rows - 1);
        let c = col.min(self.cols - 1);
        self.heights[r * self.cols + c]
    }
    /// Compute normal at grid coordinates using finite differences.
    #[allow(dead_code)]
    pub fn normal_at(&self, row: usize, col: usize) -> [f64; 3] {
        let h = self.height_at(row, col);
        let hx = if col + 1 < self.cols {
            self.height_at(row, col + 1)
        } else {
            h
        };
        let hz = if row + 1 < self.rows {
            self.height_at(row + 1, col)
        } else {
            h
        };
        let dx = hx - h;
        let dz = hz - h;
        let nx = -dx;
        let ny = self.spacing;
        let nz = -dz;
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len < 1e-12 {
            [0.0, 1.0, 0.0]
        } else {
            [nx / len, ny / len, nz / len]
        }
    }
    /// Test a sphere against the heightfield.
    #[allow(dead_code)]
    pub fn test_sphere(
        &self,
        sphere_center: [f64; 3],
        sphere_radius: f64,
    ) -> Option<([f64; 3], [f64; 3], f64)> {
        let col_f = sphere_center[0] / self.spacing;
        let row_f = sphere_center[2] / self.spacing;
        if col_f < 0.0 || row_f < 0.0 {
            return None;
        }
        let col = col_f as usize;
        let row = row_f as usize;
        if col >= self.cols || row >= self.rows {
            return None;
        }
        let h = self.height_at(row, col);
        let penetration = h + sphere_radius - sphere_center[1];
        if penetration > 0.0 {
            let normal = self.normal_at(row, col);
            let point = [sphere_center[0], h, sphere_center[2]];
            Some((point, normal, penetration))
        } else {
            None
        }
    }
}
/// A simple priority dispatch queue.
///
/// Entries with higher priority (deeper penetration estimate) are resolved first,
/// improving simulation stability when many collisions occur simultaneously.
#[allow(dead_code)]
pub struct DispatchQueue {
    pub(super) entries: Vec<DispatchQueueEntry>,
}
impl DispatchQueue {
    /// Create an empty dispatch queue.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }
    /// Push a new entry.
    #[allow(dead_code)]
    pub fn push(&mut self, pair: CollisionPair, priority: f64) {
        self.entries.push(DispatchQueueEntry { pair, priority });
    }
    /// Pop the highest-priority entry.
    #[allow(dead_code)]
    pub fn pop(&mut self) -> Option<DispatchQueueEntry> {
        if self.entries.is_empty() {
            return None;
        }
        let idx = self
            .entries
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.priority
                    .partial_cmp(&b.priority)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)?;
        Some(self.entries.swap_remove(idx))
    }
    /// Number of entries in the queue.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the queue is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the queue.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
/// A priority entry in a dispatch queue.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DispatchQueueEntry {
    /// Collision pair.
    pub pair: CollisionPair,
    /// Estimated penetration depth (higher = higher priority).
    pub priority: f64,
}
/// Dispatch result for a compound vs compound pair.
#[derive(Debug, Clone)]
pub struct CompoundDispatchResult {
    /// All contact manifolds found.
    pub manifolds: Vec<ContactManifold>,
    /// Number of narrowphase tests performed.
    pub tests_performed: usize,
    /// Number of pairs pruned by AABB overlap check.
    pub pairs_pruned: usize,
}
impl CompoundDispatchResult {
    /// Whether any contacts were found.
    pub fn has_contacts(&self) -> bool {
        !self.manifolds.is_empty()
    }
    /// Total contact count across all manifolds.
    pub fn total_contacts(&self) -> usize {
        self.manifolds.iter().map(|m| m.contacts.len()).sum()
    }
}
/// A contact patch built from manifold reduction.
///
/// After collision detection returns many raw contact points, the patch
/// reducer keeps only the most geometrically useful subset (typically ≤ 4).
#[allow(dead_code)]
pub struct ContactPatchReducer {
    /// Maximum number of contact points to keep.
    pub max_contacts: usize,
}
impl ContactPatchReducer {
    /// Create a reducer keeping at most `max_contacts` points.
    #[allow(dead_code)]
    pub fn new(max_contacts: usize) -> Self {
        Self { max_contacts }
    }
    /// Reduce a contact manifold to the most useful contacts.
    ///
    /// Strategy:
    /// 1. Keep the deepest contact.
    /// 2. Keep the contact farthest from the deepest.
    /// 3. Keep contacts that maximise triangle area of the current set.
    #[allow(dead_code)]
    pub fn reduce(&self, manifold: &mut ContactManifold) {
        if manifold.contacts.len() <= self.max_contacts {
            return;
        }
        let mut kept: Vec<Contact> = Vec::with_capacity(self.max_contacts);
        if let Some(idx) = manifold
            .contacts
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.depth
                    .partial_cmp(&b.depth)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
        {
            kept.push(manifold.contacts[idx].clone());
        }
        if kept.len() < self.max_contacts && !manifold.contacts.is_empty() {
            let ref_pt = kept[0].point_a;
            if let Some(idx) = manifold
                .contacts
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| {
                    let da = (a.point_a - ref_pt).norm_squared();
                    let db = (b.point_a - ref_pt).norm_squared();
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
            {
                kept.push(manifold.contacts[idx].clone());
            }
        }
        while kept.len() < self.max_contacts {
            let mut best_area = -1.0;
            let mut best_idx = 0;
            'outer: for (i, c) in manifold.contacts.iter().enumerate() {
                for k in &kept {
                    if (c.point_a - k.point_a).norm_squared() < 1e-10 {
                        continue 'outer;
                    }
                }
                let area = if kept.len() >= 2 {
                    let ab = kept[1].point_a - kept[0].point_a;
                    let ac = c.point_a - kept[0].point_a;
                    ab.cross(&ac).norm()
                } else {
                    (c.point_a - kept[0].point_a).norm()
                };
                if area > best_area {
                    best_area = area;
                    best_idx = i;
                }
            }
            if best_area <= 0.0 {
                break;
            }
            kept.push(manifold.contacts[best_idx].clone());
        }
        manifold.contacts = kept;
    }
}
/// A compound shape: a collection of sub-shapes with local transforms.
#[allow(dead_code)]
pub struct CompoundShape {
    /// Sub-shapes.
    pub children: Vec<Box<dyn Shape>>,
    /// Sub-shape types for dispatch.
    pub child_types: Vec<ShapeType>,
    /// Local transform of each sub-shape relative to the compound origin.
    pub local_transforms: Vec<Transform>,
}
impl CompoundShape {
    /// Create an empty compound shape.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            children: Vec::new(),
            child_types: Vec::new(),
            local_transforms: Vec::new(),
        }
    }
    /// Add a child shape.
    #[allow(dead_code)]
    pub fn add_child(
        &mut self,
        shape: Box<dyn Shape>,
        shape_type: ShapeType,
        local_transform: Transform,
    ) {
        self.children.push(shape);
        self.child_types.push(shape_type);
        self.local_transforms.push(local_transform);
    }
    /// Number of children.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.children.len()
    }
    /// Whether the compound has no children.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.children.is_empty()
    }
}
/// Result returned by a narrow-phase algorithm.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct NarrowPhaseResult {
    /// Optional contact manifold (None means separated).
    pub manifold: Option<ContactManifold>,
}
impl NarrowPhaseResult {
    /// Construct a result indicating no contact.
    #[allow(dead_code)]
    pub fn separated() -> Self {
        NarrowPhaseResult { manifold: None }
    }
    /// Construct a result with a contact manifold.
    #[allow(dead_code)]
    pub fn contact(manifold: ContactManifold) -> Self {
        NarrowPhaseResult {
            manifold: Some(manifold),
        }
    }
    /// Returns true if a contact was found.
    #[allow(dead_code)]
    pub fn has_contact(&self) -> bool {
        self.manifold.is_some()
    }
    /// Returns the penetration depth if contact exists.
    #[allow(dead_code)]
    pub fn penetration_depth(&self) -> Option<f64> {
        self.manifold
            .as_ref()
            .and_then(|m| m.contacts.first())
            .map(|c| c.depth)
    }
}
/// Concave mesh represented as a convex decomposition.
#[allow(dead_code)]
pub struct ConcaveMesh {
    /// Convex pieces of the decomposition.
    pub pieces: Vec<ConvexPiece>,
}
impl ConcaveMesh {
    /// Create a concave mesh from a list of convex pieces.
    #[allow(dead_code)]
    pub fn new(pieces: Vec<ConvexPiece>) -> Self {
        Self { pieces }
    }
    /// Number of convex pieces.
    #[allow(dead_code)]
    pub fn num_pieces(&self) -> usize {
        self.pieces.len()
    }
}
/// Cache for frequently reused shape geometric features.
#[derive(Debug, Clone)]
pub struct ShapeFeatureCache {
    /// Cached support point for a given direction hash.
    pub(super) support_cache: std::collections::HashMap<u64, [f64; 3]>,
    /// Cache hit count.
    pub hits: usize,
    /// Cache miss count.
    pub misses: usize,
}
impl ShapeFeatureCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self {
            support_cache: std::collections::HashMap::new(),
            hits: 0,
            misses: 0,
        }
    }
    /// Look up or compute the support point for a direction.
    ///
    /// `dir_key` is a hash of the direction vector (caller's responsibility).
    pub fn get_or_insert(&mut self, dir_key: u64, compute: impl FnOnce() -> [f64; 3]) -> [f64; 3] {
        if let Some(&cached) = self.support_cache.get(&dir_key) {
            self.hits += 1;
            cached
        } else {
            self.misses += 1;
            let val = compute();
            self.support_cache.insert(dir_key, val);
            val
        }
    }
    /// Cache hit ratio.
    pub fn hit_ratio(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            self.hits as f64 / total as f64
        }
    }
    /// Clear the cache.
    pub fn clear(&mut self) {
        self.support_cache.clear();
        self.hits = 0;
        self.misses = 0;
    }
}
/// Key for looking up a registered collision algorithm.
/// Order is canonicalised: the lower-ordinal type comes first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(dead_code)]
pub struct DispatchKey(pub ShapeType, pub ShapeType);
impl DispatchKey {
    /// Create a canonicalised dispatch key (order-independent).
    #[allow(dead_code)]
    pub fn new(a: ShapeType, b: ShapeType) -> Self {
        if shape_type_ordinal(a) <= shape_type_ordinal(b) {
            DispatchKey(a, b)
        } else {
            DispatchKey(b, a)
        }
    }
}
