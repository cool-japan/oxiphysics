//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::KD_LEAF_SIZE;
use super::functions::dist_sq;

/// A single point with optional attributes.
#[derive(Debug, Clone, Default)]
pub struct PointCloudPoint {
    /// Position `[x, y, z]`.
    pub position: [f64; 3],
    /// Normal vector `[nx, ny, nz]` (unit length, if present).
    pub normal: Option<[f64; 3]>,
    /// RGB color `[r, g, b]` in `[0, 255]`.
    pub color: Option<[u8; 3]>,
    /// Scalar intensity value.
    pub intensity: Option<f64>,
    /// Classification label.
    pub classification: Option<u8>,
}
/// A node in a KD-tree for 3D points.
#[derive(Debug)]
enum KdNode {
    /// Leaf node containing point indices.
    Leaf {
        /// Indices into the original point cloud.
        indices: Vec<usize>,
    },
    /// Internal node splitting on a given axis.
    Split {
        /// Axis (0=x, 1=y, 2=z).
        axis: usize,
        /// Split value.
        split_value: f64,
        /// Left subtree (points with coord <= split_value).
        left: Box<KdNode>,
        /// Right subtree.
        right: Box<KdNode>,
    },
}
/// A KD-tree for efficient nearest neighbor queries in 3D.
#[derive(Debug)]
pub struct KdTree {
    /// Root node.
    root: Option<KdNode>,
    /// Reference to positions (copied).
    positions: Vec<[f64; 3]>,
}
impl KdTree {
    /// Build a KD-tree from a point cloud.
    pub fn build(cloud: &PointCloud) -> Self {
        let positions = cloud.positions.clone();
        let indices: Vec<usize> = (0..positions.len()).collect();
        let root = if indices.is_empty() {
            None
        } else {
            Some(Self::build_recursive(&positions, indices, 0))
        };
        Self { root, positions }
    }
    /// Build from raw positions.
    pub fn from_positions(positions: &[[f64; 3]]) -> Self {
        let pos = positions.to_vec();
        let indices: Vec<usize> = (0..pos.len()).collect();
        let root = if indices.is_empty() {
            None
        } else {
            Some(Self::build_recursive(&pos, indices, 0))
        };
        Self {
            root,
            positions: pos,
        }
    }
    fn build_recursive(positions: &[[f64; 3]], mut indices: Vec<usize>, depth: usize) -> KdNode {
        if indices.len() <= KD_LEAF_SIZE {
            return KdNode::Leaf { indices };
        }
        let axis = depth % 3;
        indices.sort_by(|&a, &b| {
            positions[a][axis]
                .partial_cmp(&positions[b][axis])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = indices.len() / 2;
        let split_value = positions[indices[mid]][axis];
        let left_indices = indices[..mid].to_vec();
        let right_indices = indices[mid..].to_vec();
        KdNode::Split {
            axis,
            split_value,
            left: Box::new(Self::build_recursive(positions, left_indices, depth + 1)),
            right: Box::new(Self::build_recursive(positions, right_indices, depth + 1)),
        }
    }
    /// Find the nearest neighbor to a query point.
    ///
    /// Returns `(index, squared_distance)`.
    pub fn nearest(&self, query: [f64; 3]) -> Option<(usize, f64)> {
        let root = self.root.as_ref()?;
        let mut best_idx = 0;
        let mut best_dist_sq = f64::MAX;
        Self::nearest_recursive(
            root,
            &self.positions,
            query,
            &mut best_idx,
            &mut best_dist_sq,
        );
        Some((best_idx, best_dist_sq))
    }
    fn nearest_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        best_idx: &mut usize,
        best_dist_sq: &mut f64,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if d < *best_dist_sq {
                        *best_dist_sq = d;
                        *best_idx = i;
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::nearest_recursive(first, positions, query, best_idx, best_dist_sq);
                if diff * diff < *best_dist_sq {
                    Self::nearest_recursive(second, positions, query, best_idx, best_dist_sq);
                }
            }
        }
    }
    /// Find k nearest neighbors.
    ///
    /// Returns a vector of `(index, squared_distance)` sorted by distance.
    pub fn k_nearest(&self, query: [f64; 3], k: usize) -> Vec<(usize, f64)> {
        if self.root.is_none() || k == 0 {
            return Vec::new();
        }
        let mut results: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        Self::knn_recursive(
            self.root.as_ref().expect("value should be Some"),
            &self.positions,
            query,
            k,
            &mut results,
        );
        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        results.truncate(k);
        results
    }
    fn knn_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        k: usize,
        results: &mut Vec<(usize, f64)>,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if results.len() < k {
                        results.push((i, d));
                        results.sort_by(|a, b| {
                            a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                        });
                    } else if d < results[k - 1].1 {
                        results[k - 1] = (i, d);
                        results.sort_by(|a, b| {
                            a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
                        });
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::knn_recursive(first, positions, query, k, results);
                let worst_dist = if results.len() < k {
                    f64::MAX
                } else {
                    results[results.len() - 1].1
                };
                if diff * diff < worst_dist {
                    Self::knn_recursive(second, positions, query, k, results);
                }
            }
        }
    }
    /// Find all points within a given radius (squared).
    pub fn radius_search(&self, query: [f64; 3], radius_sq: f64) -> Vec<(usize, f64)> {
        let mut results = Vec::new();
        if let Some(ref root) = self.root {
            Self::radius_recursive(root, &self.positions, query, radius_sq, &mut results);
        }
        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        results
    }
    fn radius_recursive(
        node: &KdNode,
        positions: &[[f64; 3]],
        query: [f64; 3],
        radius_sq: f64,
        results: &mut Vec<(usize, f64)>,
    ) {
        match node {
            KdNode::Leaf { indices } => {
                for &i in indices {
                    let d = dist_sq(positions[i], query);
                    if d <= radius_sq {
                        results.push((i, d));
                    }
                }
            }
            KdNode::Split {
                axis,
                split_value,
                left,
                right,
            } => {
                let diff = query[*axis] - split_value;
                let (first, second) = if diff <= 0.0 {
                    (left.as_ref(), right.as_ref())
                } else {
                    (right.as_ref(), left.as_ref())
                };
                Self::radius_recursive(first, positions, query, radius_sq, results);
                if diff * diff <= radius_sq {
                    Self::radius_recursive(second, positions, query, radius_sq, results);
                }
            }
        }
    }
    /// Number of points in the tree.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether the tree is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}
/// E57 scan header (simplified).
#[derive(Debug, Clone)]
pub struct E57ScanHeader {
    /// Scan name/description.
    pub name: String,
    /// Number of points in the scan.
    pub num_points: u64,
    /// Whether the scan has color data.
    pub has_color: bool,
    /// Whether the scan has intensity data.
    pub has_intensity: bool,
    /// Scanner position `[x, y, z]`.
    pub scanner_position: [f64; 3],
    /// Scanner orientation (quaternion `[x, y, z, w]`).
    pub scanner_orientation: [f64; 4],
}
impl E57ScanHeader {
    /// Create a mock E57 scan header for testing.
    pub fn mock(name: &str, num_points: u64) -> Self {
        Self {
            name: name.to_string(),
            num_points,
            has_color: true,
            has_intensity: true,
            scanner_position: [0.0; 3],
            scanner_orientation: [0.0, 0.0, 0.0, 1.0],
        }
    }
}
/// PLY file encoding type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlyEncoding {
    /// ASCII text format.
    Ascii,
    /// Binary little-endian format.
    BinaryLittleEndian,
}
/// A point cloud: a collection of points with optional per-point attributes.
#[derive(Debug, Clone)]
pub struct PointCloud {
    /// Positions of all points.
    pub positions: Vec<[f64; 3]>,
    /// Per-point normals (if present, same length as `positions`).
    pub normals: Option<Vec<[f64; 3]>>,
    /// Per-point RGB colors (if present).
    pub colors: Option<Vec<[u8; 3]>>,
    /// Per-point intensity values (if present).
    pub intensities: Option<Vec<f64>>,
    /// Per-point classification labels (if present).
    pub classifications: Option<Vec<u8>>,
}
impl PointCloud {
    /// Create an empty point cloud.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            normals: None,
            colors: None,
            intensities: None,
            classifications: None,
        }
    }
    /// Create a point cloud from positions only.
    pub fn from_positions(positions: Vec<[f64; 3]>) -> Self {
        Self {
            positions,
            normals: None,
            colors: None,
            intensities: None,
            classifications: None,
        }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether the cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Whether normals are present.
    pub fn has_normals(&self) -> bool {
        self.normals.is_some()
    }
    /// Whether colors are present.
    pub fn has_colors(&self) -> bool {
        self.colors.is_some()
    }
    /// Whether intensities are present.
    pub fn has_intensities(&self) -> bool {
        self.intensities.is_some()
    }
    /// Add a point with position only.
    pub fn add_point(&mut self, pos: [f64; 3]) {
        self.positions.push(pos);
        if let Some(ref mut n) = self.normals {
            n.push([0.0; 3]);
        }
        if let Some(ref mut c) = self.colors {
            c.push([255, 255, 255]);
        }
        if let Some(ref mut i) = self.intensities {
            i.push(0.0);
        }
        if let Some(ref mut cl) = self.classifications {
            cl.push(0);
        }
    }
    /// Add a point with all attributes.
    pub fn add_full_point(&mut self, point: PointCloudPoint) {
        self.positions.push(point.position);
        if let Some(ref mut n) = self.normals {
            n.push(point.normal.unwrap_or([0.0; 3]));
        } else if let Some(n) = point.normal {
            let mut norms = vec![[0.0_f64; 3]; self.positions.len() - 1];
            norms.push(n);
            self.normals = Some(norms);
        }
        if let Some(ref mut c) = self.colors {
            c.push(point.color.unwrap_or([255, 255, 255]));
        } else if let Some(c) = point.color {
            let mut cols = vec![[255_u8; 3]; self.positions.len() - 1];
            cols.push(c);
            self.colors = Some(cols);
        }
        if let Some(ref mut i) = self.intensities {
            i.push(point.intensity.unwrap_or(0.0));
        } else if let Some(intensity) = point.intensity {
            let mut ints = vec![0.0_f64; self.positions.len() - 1];
            ints.push(intensity);
            self.intensities = Some(ints);
        }
        if let Some(ref mut cl) = self.classifications {
            cl.push(point.classification.unwrap_or(0));
        }
    }
    /// Get a single point as a `PointCloudPoint`.
    pub fn get_point(&self, idx: usize) -> Option<PointCloudPoint> {
        if idx >= self.positions.len() {
            return None;
        }
        Some(PointCloudPoint {
            position: self.positions[idx],
            normal: self.normals.as_ref().map(|n| n[idx]),
            color: self.colors.as_ref().map(|c| c[idx]),
            intensity: self.intensities.as_ref().map(|i| i[idx]),
            classification: self.classifications.as_ref().map(|c| c[idx]),
        })
    }
    /// Reserve capacity for `n` additional points.
    pub fn reserve(&mut self, n: usize) {
        self.positions.reserve(n);
        if let Some(ref mut normals) = self.normals {
            normals.reserve(n);
        }
        if let Some(ref mut colors) = self.colors {
            colors.reserve(n);
        }
        if let Some(ref mut ints) = self.intensities {
            ints.reserve(n);
        }
    }
}
// Default impl is in trait_impls.rs to avoid duplication
/// LAS file header structure (subset of fields).
#[derive(Debug, Clone)]
pub struct LasHeader {
    /// File signature (should be "LASF").
    pub file_signature: String,
    /// File source ID.
    pub file_source_id: u16,
    /// Version major.
    pub version_major: u8,
    /// Version minor.
    pub version_minor: u8,
    /// Point data record format.
    pub point_format: u8,
    /// Point data record length.
    pub point_record_length: u16,
    /// Number of point records.
    pub num_points: u64,
    /// Scale factors `[sx, sy, sz]`.
    pub scale: [f64; 3],
    /// Offset `[ox, oy, oz]`.
    pub offset: [f64; 3],
    /// Min bounds `[min_x, min_y, min_z]`.
    pub min_bounds: [f64; 3],
    /// Max bounds `[max_x, max_y, max_z]`.
    pub max_bounds: [f64; 3],
}
impl LasHeader {
    /// Create a mock LAS header for testing.
    pub fn mock(num_points: u64) -> Self {
        Self {
            file_signature: "LASF".to_string(),
            file_source_id: 0,
            version_major: 1,
            version_minor: 4,
            point_format: 0,
            point_record_length: 20,
            num_points,
            scale: [0.001, 0.001, 0.001],
            offset: [0.0, 0.0, 0.0],
            min_bounds: [0.0, 0.0, 0.0],
            max_bounds: [100.0, 100.0, 50.0],
        }
    }
    /// Parse a LAS header from raw bytes (minimal parsing).
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 227 {
            return None;
        }
        let sig = std::str::from_utf8(&data[0..4]).ok()?;
        if sig != "LASF" {
            return None;
        }
        let file_source_id = u16::from_le_bytes([data[4], data[5]]);
        let version_major = data[24];
        let version_minor = data[25];
        let point_format = data[104];
        let point_record_length = u16::from_le_bytes([data[105], data[106]]);
        let num_points_legacy =
            u32::from_le_bytes([data[107], data[108], data[109], data[110]]) as u64;
        let parse_f64 = |offset: usize| -> f64 {
            if offset + 8 > data.len() {
                return 0.0;
            }
            let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap_or([0; 8]);
            f64::from_le_bytes(bytes)
        };
        let scale = [parse_f64(131), parse_f64(139), parse_f64(147)];
        let offset = [parse_f64(155), parse_f64(163), parse_f64(171)];
        let max_x = parse_f64(179);
        let min_x = parse_f64(187);
        let max_y = parse_f64(195);
        let min_y = parse_f64(203);
        let max_z = parse_f64(211);
        let min_z = parse_f64(219);
        Some(Self {
            file_signature: sig.to_string(),
            file_source_id,
            version_major,
            version_minor,
            point_format,
            point_record_length,
            num_points: num_points_legacy,
            scale,
            offset,
            min_bounds: [min_x, min_y, min_z],
            max_bounds: [max_x, max_y, max_z],
        })
    }
    /// Bounding box volume.
    pub fn bounding_volume(&self) -> f64 {
        let dx = self.max_bounds[0] - self.min_bounds[0];
        let dy = self.max_bounds[1] - self.min_bounds[1];
        let dz = self.max_bounds[2] - self.min_bounds[2];
        dx.abs() * dy.abs() * dz.abs()
    }
}
/// E57 file structure (simplified).
#[derive(Debug, Clone)]
pub struct E57File {
    /// Major version.
    pub version_major: u32,
    /// Minor version.
    pub version_minor: u32,
    /// Scans in the file.
    pub scans: Vec<E57ScanHeader>,
}
impl E57File {
    /// Create a mock E57 file for testing.
    pub fn mock() -> Self {
        Self {
            version_major: 1,
            version_minor: 0,
            scans: vec![E57ScanHeader::mock("scan_001", 1000)],
        }
    }
    /// Total number of points across all scans.
    pub fn total_points(&self) -> u64 {
        self.scans.iter().map(|s| s.num_points).sum()
    }
    /// Number of scans.
    pub fn num_scans(&self) -> usize {
        self.scans.len()
    }
}
/// Axis-aligned bounding box in 3D.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}
impl BoundingBox {
    /// Compute the bounding box of a point cloud.
    pub fn from_point_cloud(cloud: &PointCloud) -> Option<Self> {
        if cloud.is_empty() {
            return None;
        }
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for p in &cloud.positions {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        Some(Self { min, max })
    }
    /// Center of the bounding box.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }
    /// Extents (half-widths) of the bounding box.
    pub fn extents(&self) -> [f64; 3] {
        [
            0.5 * (self.max[0] - self.min[0]),
            0.5 * (self.max[1] - self.min[1]),
            0.5 * (self.max[2] - self.min[2]),
        ]
    }
    /// Diagonal length.
    pub fn diagonal(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Volume of the bounding box.
    pub fn volume(&self) -> f64 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        dx * dy * dz
    }
    /// Check if a point is inside the bounding box.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        p[0] >= self.min[0]
            && p[0] <= self.max[0]
            && p[1] >= self.min[1]
            && p[1] <= self.max[1]
            && p[2] >= self.min[2]
            && p[2] <= self.max[2]
    }
    /// Merge two bounding boxes.
    pub fn merge(&self, other: &BoundingBox) -> BoundingBox {
        BoundingBox {
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
/// A triangle in a mesh (indices into a point cloud).
#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    /// Vertex indices.
    pub indices: [usize; 3],
}
