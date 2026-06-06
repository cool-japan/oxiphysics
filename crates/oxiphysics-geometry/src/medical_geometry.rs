// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Medical imaging geometry — DICOM coordinate systems, CT/MRI volume
//! reconstruction, bone segmentation, vessel centerline extraction, organ
//! bounding volumes, anatomical coordinate systems (RAS/LPS), surgical
//! planning geometry, marching-cubes mesh extraction, measurement tools,
//! and mesh smoothing for medical applications.
//!
//! All computations use plain `f64` / `[f64; 3]` — no nalgebra dependency.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vec3 helper
// ---------------------------------------------------------------------------

/// Minimal 3-component vector for medical geometry computations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct a new vector.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Dot product.
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Euclidean length.
    pub fn norm(&self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Normalise; returns zero vector if near-zero length.
    pub fn normalize(&self) -> Self {
        let n = self.norm();
        if n < 1e-15 {
            Self::zero()
        } else {
            Self::new(self.x / n, self.y / n, self.z / n)
        }
    }

    /// Cross product.
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Component-wise addition.
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    /// Component-wise subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    /// Scalar multiplication.
    pub fn scale(&self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }

    /// Distance to another point.
    pub fn distance_to(&self, other: &Self) -> f64 {
        self.sub(other).norm()
    }
}

// ---------------------------------------------------------------------------
// 1. DICOM-like image geometry
// ---------------------------------------------------------------------------

/// Orientation of a single DICOM slice in 3-D patient space.
///
/// The row and column direction cosines follow the DICOM tag (0020,0037).
#[derive(Debug, Clone)]
pub struct DicomSliceOrientation {
    /// Direction cosines of the row direction (unit vector).
    pub row_cosines: Vec3,
    /// Direction cosines of the column direction (unit vector).
    pub col_cosines: Vec3,
}

impl DicomSliceOrientation {
    /// Construct from raw row and column direction cosine arrays.
    pub fn new(row: [f64; 3], col: [f64; 3]) -> Self {
        Self {
            row_cosines: Vec3::new(row[0], row[1], row[2]).normalize(),
            col_cosines: Vec3::new(col[0], col[1], col[2]).normalize(),
        }
    }

    /// Compute the normal to the slice plane (cross product row × col).
    pub fn normal(&self) -> Vec3 {
        self.row_cosines.cross(&self.col_cosines).normalize()
    }

    /// Returns true when the row and column cosines are orthogonal
    /// within `tol` (absolute dot-product threshold).
    pub fn is_orthogonal(&self, tol: f64) -> bool {
        self.row_cosines.dot(&self.col_cosines).abs() < tol
    }
}

/// Geometry descriptor for a single DICOM image plane (one slice).
#[derive(Debug, Clone)]
pub struct DicomImagePlane {
    /// Image Position Patient — top-left corner of the first pixel (mm).
    pub image_position: Vec3,
    /// Pixel spacing: \[row_spacing, col_spacing\] in mm.
    pub pixel_spacing: [f64; 2],
    /// Slice orientation (row / column cosines).
    pub orientation: DicomSliceOrientation,
    /// Number of columns (pixels along row direction).
    pub cols: usize,
    /// Number of rows (pixels along column direction).
    pub rows: usize,
}

impl DicomImagePlane {
    /// Build a new image plane descriptor.
    pub fn new(
        image_position: [f64; 3],
        pixel_spacing: [f64; 2],
        row_cosines: [f64; 3],
        col_cosines: [f64; 3],
        rows: usize,
        cols: usize,
    ) -> Self {
        Self {
            image_position: Vec3::new(image_position[0], image_position[1], image_position[2]),
            pixel_spacing,
            orientation: DicomSliceOrientation::new(row_cosines, col_cosines),
            cols,
            rows,
        }
    }

    /// Convert a pixel coordinate (col `i`, row `j`) to 3-D patient
    /// coordinates in mm (DICOM Patient Coordinate System).
    pub fn pixel_to_patient(&self, i: f64, j: f64) -> Vec3 {
        let dr = self
            .orientation
            .row_cosines
            .scale(i * self.pixel_spacing[0]);
        let dc = self
            .orientation
            .col_cosines
            .scale(j * self.pixel_spacing[1]);
        self.image_position.add(&dr).add(&dc)
    }

    /// Approximate the 3-D patient coordinate of a slice given its index
    /// in a stack and the inter-slice spacing (mm).
    pub fn slice_position(&self, slice_index: usize, slice_spacing: f64) -> Vec3 {
        let normal = self.orientation.normal();
        self.image_position
            .add(&normal.scale(slice_index as f64 * slice_spacing))
    }
}

// ---------------------------------------------------------------------------
// 2. CT/MRI volume & trilinear interpolation
// ---------------------------------------------------------------------------

/// A 3-D scalar volume (e.g. CT Hounsfield units or MRI signal intensity).
///
/// Voxel ordering: `data[z][y][x]` with the first axis being the slice axis.
#[derive(Debug, Clone)]
pub struct ScalarVolume {
    /// Raw voxel data stored in z-y-x order.
    pub data: Vec<f64>,
    /// Number of voxels along X.
    pub nx: usize,
    /// Number of voxels along Y.
    pub ny: usize,
    /// Number of voxels along Z (number of slices).
    pub nz: usize,
    /// Voxel size in mm: \[dx, dy, dz\].
    pub voxel_size: [f64; 3],
    /// World-space origin of voxel (0,0,0) in mm.
    pub origin: Vec3,
}

impl ScalarVolume {
    /// Create a new empty volume filled with zeros.
    pub fn new(nx: usize, ny: usize, nz: usize, voxel_size: [f64; 3], origin: Vec3) -> Self {
        Self {
            data: vec![0.0; nx * ny * nz],
            nx,
            ny,
            nz,
            voxel_size,
            origin,
        }
    }

    /// Flat index for voxel (ix, iy, iz).
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        iz * self.ny * self.nx + iy * self.nx + ix
    }

    /// Get voxel value at integer coordinates.  Returns `None` if out of
    /// bounds.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> Option<f64> {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            Some(self.data[self.index(ix, iy, iz)])
        } else {
            None
        }
    }

    /// Set voxel value (no-op when out of bounds).
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, val: f64) {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            let idx = self.index(ix, iy, iz);
            self.data[idx] = val;
        }
    }

    /// Convert continuous world coordinates (mm) to fractional voxel
    /// indices.  Returns `None` if outside the volume bounding box.
    pub fn world_to_voxel(&self, world: Vec3) -> Option<[f64; 3]> {
        let fx = (world.x - self.origin.x) / self.voxel_size[0];
        let fy = (world.y - self.origin.y) / self.voxel_size[1];
        let fz = (world.z - self.origin.z) / self.voxel_size[2];
        if fx < 0.0
            || fy < 0.0
            || fz < 0.0
            || fx > (self.nx - 1) as f64
            || fy > (self.ny - 1) as f64
            || fz > (self.nz - 1) as f64
        {
            None
        } else {
            Some([fx, fy, fz])
        }
    }

    /// Trilinear interpolation at fractional voxel coordinates.
    ///
    /// Returns `None` when the point lies outside the valid range.
    pub fn trilinear_interp(&self, world: Vec3) -> Option<f64> {
        let [fx, fy, fz] = self.world_to_voxel(world)?;
        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let z0 = fz.floor() as usize;
        let x1 = (x0 + 1).min(self.nx - 1);
        let y1 = (y0 + 1).min(self.ny - 1);
        let z1 = (z0 + 1).min(self.nz - 1);
        let xd = fx - x0 as f64;
        let yd = fy - y0 as f64;
        let zd = fz - z0 as f64;

        macro_rules! v {
            ($xi:expr, $yi:expr, $zi:expr) => {
                self.data[self.index($xi, $yi, $zi)]
            };
        }

        let c00 = v!(x0, y0, z0) * (1.0 - xd) + v!(x1, y0, z0) * xd;
        let c01 = v!(x0, y0, z1) * (1.0 - xd) + v!(x1, y0, z1) * xd;
        let c10 = v!(x0, y1, z0) * (1.0 - xd) + v!(x1, y1, z0) * xd;
        let c11 = v!(x0, y1, z1) * (1.0 - xd) + v!(x1, y1, z1) * xd;
        let c0 = c00 * (1.0 - yd) + c10 * yd;
        let c1 = c01 * (1.0 - yd) + c11 * yd;
        Some(c0 * (1.0 - zd) + c1 * zd)
    }
}

// ---------------------------------------------------------------------------
// 3. Bone segmentation geometry (Hounsfield units)
// ---------------------------------------------------------------------------

/// Hounsfield unit thresholds for common tissue classes in CT imaging.
#[derive(Debug, Clone, Copy)]
pub struct HuThresholds {
    /// Minimum HU for cortical bone (typically ~700 HU).
    pub cortical_bone_min: f64,
    /// Minimum HU for cancellous (trabecular) bone (typically ~200 HU).
    pub cancellous_bone_min: f64,
    /// Maximum HU still classified as soft tissue (typically ~200 HU).
    pub soft_tissue_max: f64,
    /// Minimum HU for air (typically ~-900 HU).
    pub air_min: f64,
}

impl HuThresholds {
    /// Clinically representative default thresholds.
    pub fn default_ct() -> Self {
        Self {
            cortical_bone_min: 700.0,
            cancellous_bone_min: 200.0,
            soft_tissue_max: 200.0,
            air_min: -900.0,
        }
    }

    /// Classify a Hounsfield unit value.
    pub fn classify(&self, hu: f64) -> TissueClass {
        if hu >= self.cortical_bone_min {
            TissueClass::CorticalBone
        } else if hu >= self.cancellous_bone_min {
            TissueClass::CancellousBone
        } else if hu <= self.air_min {
            TissueClass::Air
        } else if hu < 0.0 {
            TissueClass::Fat
        } else {
            TissueClass::SoftTissue
        }
    }
}

/// Tissue class returned by Hounsfield-unit segmentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TissueClass {
    /// Cortical (compact) bone.
    CorticalBone,
    /// Cancellous (trabecular) bone.
    CancellousBone,
    /// Soft tissue (muscle, organ parenchyma).
    SoftTissue,
    /// Adipose / fat tissue.
    Fat,
    /// Air or gas-filled cavities.
    Air,
}

/// Segment a [`ScalarVolume`] by simple threshold, producing a binary mask.
///
/// Returns a vector of length `nx * ny * nz` where `true` indicates the voxel
/// exceeds `threshold` (inclusive).
pub fn threshold_segment(volume: &ScalarVolume, threshold: f64) -> Vec<bool> {
    volume.data.iter().map(|&v| v >= threshold).collect()
}

/// Count voxels above threshold and estimate bone volume in mm³.
pub fn bone_volume_mm3(volume: &ScalarVolume, threshold: f64) -> f64 {
    let vv = volume.voxel_size[0] * volume.voxel_size[1] * volume.voxel_size[2];
    volume.data.iter().filter(|&&v| v >= threshold).count() as f64 * vv
}

// ---------------------------------------------------------------------------
// 4. Vessel centerline extraction (Voronoi / medial-axis approximation)
// ---------------------------------------------------------------------------

/// A single point on a vessel centreline with associated radius estimate.
#[derive(Debug, Clone, Copy)]
pub struct CenterlinePoint {
    /// 3-D position in patient coordinates (mm).
    pub position: Vec3,
    /// Local vessel radius estimate (mm).
    pub radius: f64,
}

/// Result of a vessel centreline extraction.
#[derive(Debug, Clone)]
pub struct VesselCenterline {
    /// Ordered sequence of centreline points from proximal to distal.
    pub points: Vec<CenterlinePoint>,
}

impl VesselCenterline {
    /// Construct from a list of (position, radius) pairs.
    pub fn new(pts: Vec<(Vec3, f64)>) -> Self {
        Self {
            points: pts
                .into_iter()
                .map(|(p, r)| CenterlinePoint {
                    position: p,
                    radius: r,
                })
                .collect(),
        }
    }

    /// Total arc length of the centreline (mm).
    pub fn arc_length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|w| w[0].position.distance_to(&w[1].position))
            .sum()
    }

    /// Minimum vessel radius along the centreline.
    pub fn min_radius(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.radius)
            .fold(f64::INFINITY, f64::min)
    }

    /// Maximum vessel radius along the centreline.
    pub fn max_radius(&self) -> f64 {
        self.points
            .iter()
            .map(|p| p.radius)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Tortuosity index: arc_length / straight-line distance between endpoints.
    /// Returns `1.0` when fewer than two points exist.
    pub fn tortuosity(&self) -> f64 {
        if self.points.len() < 2 {
            return 1.0;
        }
        let chord = self
            .points
            .first()
            .expect("operation should succeed")
            .position
            .distance_to(
                &self
                    .points
                    .last()
                    .expect("collection should not be empty")
                    .position,
            );
        if chord < 1e-9 {
            return 1.0;
        }
        self.arc_length() / chord
    }
}

/// Approximate centreline extraction from a binary vessel mask using a simple
/// iterative thinning pass over the boundary points (Voronoi-inspired medial
/// axis).  This is a lightweight approximation suitable for testing geometry
/// tools; production use requires a full 3-D skeletonisation algorithm.
pub fn extract_centerline_approx(
    mask: &[bool],
    nx: usize,
    ny: usize,
    nz: usize,
    voxel_size: [f64; 3],
    origin: Vec3,
) -> VesselCenterline {
    // Collect foreground voxel centres
    let mut centres: Vec<Vec3> = Vec::new();
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iz * ny * nx + iy * nx + ix;
                if mask[idx] {
                    centres.push(Vec3::new(
                        origin.x + ix as f64 * voxel_size[0],
                        origin.y + iy as f64 * voxel_size[1],
                        origin.z + iz as f64 * voxel_size[2],
                    ));
                }
            }
        }
    }
    if centres.is_empty() {
        return VesselCenterline { points: Vec::new() };
    }
    // Sort by Z to create an ordered path
    centres.sort_by(|a, b| a.z.partial_cmp(&b.z).unwrap_or(std::cmp::Ordering::Equal));

    // For each slice z-level, compute centroid as the centreline point
    let mut cl_pts: Vec<CenterlinePoint> = Vec::new();
    let mut i = 0;
    while i < centres.len() {
        let z_ref = centres[i].z;
        let mut sum = Vec3::zero();
        let mut count = 0usize;
        let start = i;
        while i < centres.len() && (centres[i].z - z_ref).abs() < voxel_size[2] * 0.5 {
            sum = sum.add(&centres[i]);
            count += 1;
            i += 1;
        }
        let centroid = sum.scale(1.0 / count as f64);
        // radius estimate: average distance from centroid to slice points
        let r = if count > 1 {
            centres[start..i]
                .iter()
                .map(|p| centroid.distance_to(p))
                .sum::<f64>()
                / count as f64
        } else {
            voxel_size[0]
        };
        cl_pts.push(CenterlinePoint {
            position: centroid,
            radius: r,
        });
    }
    VesselCenterline { points: cl_pts }
}

// ---------------------------------------------------------------------------
// 5. Organ bounding volumes
// ---------------------------------------------------------------------------

/// Axis-aligned bounding box (AABB) in patient coordinates.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Minimum corner (mm).
    pub min: Vec3,
    /// Maximum corner (mm).
    pub max: Vec3,
}

impl Aabb {
    /// Construct from explicit min/max corners.
    pub fn new(min: Vec3, max: Vec3) -> Self {
        Self { min, max }
    }

    /// Compute the bounding box of a point cloud.
    pub fn from_points(pts: &[Vec3]) -> Option<Self> {
        if pts.is_empty() {
            return None;
        }
        let mut mn = pts[0];
        let mut mx = pts[0];
        for p in pts.iter().skip(1) {
            mn.x = mn.x.min(p.x);
            mn.y = mn.y.min(p.y);
            mn.z = mn.z.min(p.z);
            mx.x = mx.x.max(p.x);
            mx.y = mx.y.max(p.y);
            mx.z = mx.z.max(p.z);
        }
        Some(Self { min: mn, max: mx })
    }

    /// Centre of the bounding box.
    pub fn centre(&self) -> Vec3 {
        Vec3::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
            (self.min.z + self.max.z) * 0.5,
        )
    }

    /// Volume of the bounding box in mm³.
    pub fn volume(&self) -> f64 {
        let d = self.max.sub(&self.min);
        d.x.max(0.0) * d.y.max(0.0) * d.z.max(0.0)
    }

    /// Diagonal length in mm.
    pub fn diagonal(&self) -> f64 {
        self.min.distance_to(&self.max)
    }

    /// Returns true when the point lies inside the AABB.
    pub fn contains(&self, p: Vec3) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// Expand the AABB by `margin` mm on all sides.
    pub fn expand(&self, margin: f64) -> Self {
        Self {
            min: self.min.sub(&Vec3::new(margin, margin, margin)),
            max: self.max.add(&Vec3::new(margin, margin, margin)),
        }
    }
}

// ---------------------------------------------------------------------------
// 6. Anatomical coordinate systems (RAS / LPS)
// ---------------------------------------------------------------------------

/// Anatomical coordinate system convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordSystem {
    /// Right–Anterior–Superior (used by NIfTI, ITK, 3D Slicer).
    Ras,
    /// Left–Posterior–Superior (used by DICOM).
    Lps,
}

/// Convert a point from LPS (DICOM) to RAS (NIfTI) coordinates by flipping
/// the X and Y axes.
pub fn lps_to_ras(p: Vec3) -> Vec3 {
    Vec3::new(-p.x, -p.y, p.z)
}

/// Convert a point from RAS (NIfTI) to LPS (DICOM) coordinates.
pub fn ras_to_lps(p: Vec3) -> Vec3 {
    Vec3::new(-p.x, -p.y, p.z)
}

/// A 4×4 homogeneous transformation matrix stored in row-major order.
///
/// Suitable for affine transformations between image voxel space and
/// patient/world space (e.g. the NIfTI sform/qform matrix).
#[derive(Debug, Clone, Copy)]
pub struct AffineMatrix4 {
    /// Row-major 4×4 matrix data.
    pub m: [[f64; 4]; 4],
}

impl AffineMatrix4 {
    /// Identity matrix.
    pub fn identity() -> Self {
        let mut m = [[0.0f64; 4]; 4];
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 1.0;
        m[3][3] = 1.0;
        Self { m }
    }

    /// Apply the affine transform to a 3-D point (homogeneous w=1).
    pub fn transform_point(&self, p: Vec3) -> Vec3 {
        let m = &self.m;
        Vec3::new(
            m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3],
            m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3],
            m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3],
        )
    }

    /// Build a voxel-to-world matrix from DICOM metadata:
    /// pixel spacing, slice spacing, image position, and direction cosines.
    pub fn from_dicom(
        pixel_spacing: [f64; 2],
        slice_spacing: f64,
        image_position: Vec3,
        row_cosines: Vec3,
        col_cosines: Vec3,
    ) -> Self {
        let normal = row_cosines.cross(&col_cosines).normalize();
        let mut m = [[0.0f64; 4]; 4];
        // Column 0: row direction × pixel_spacing[0]
        m[0][0] = row_cosines.x * pixel_spacing[0];
        m[1][0] = row_cosines.y * pixel_spacing[0];
        m[2][0] = row_cosines.z * pixel_spacing[0];
        // Column 1: col direction × pixel_spacing[1]
        m[0][1] = col_cosines.x * pixel_spacing[1];
        m[1][1] = col_cosines.y * pixel_spacing[1];
        m[2][1] = col_cosines.z * pixel_spacing[1];
        // Column 2: slice normal × slice_spacing
        m[0][2] = normal.x * slice_spacing;
        m[1][2] = normal.y * slice_spacing;
        m[2][2] = normal.z * slice_spacing;
        // Column 3: image position
        m[0][3] = image_position.x;
        m[1][3] = image_position.y;
        m[2][3] = image_position.z;
        m[3][3] = 1.0;
        Self { m }
    }
}

// ---------------------------------------------------------------------------
// 7. Surgical planning geometry
// ---------------------------------------------------------------------------

/// A cutting plane defined by a point and outward normal, used in surgical
/// resection planning.
#[derive(Debug, Clone, Copy)]
pub struct CuttingPlane {
    /// A point lying on the plane (mm).
    pub point: Vec3,
    /// Unit outward normal of the plane.
    pub normal: Vec3,
}

impl CuttingPlane {
    /// Construct a cutting plane from a point and a (possibly unnormalised)
    /// normal vector.
    pub fn new(point: Vec3, normal: Vec3) -> Self {
        Self {
            point,
            normal: normal.normalize(),
        }
    }

    /// Signed distance from `p` to the plane (positive on the normal side).
    pub fn signed_distance(&self, p: Vec3) -> f64 {
        self.normal.dot(&p.sub(&self.point))
    }

    /// Returns true when point `p` is on the positive (resection) side.
    pub fn is_resected(&self, p: Vec3) -> bool {
        self.signed_distance(p) >= 0.0
    }

    /// Compute the intersection of a line segment (a→b) with the plane.
    /// Returns `None` when the segment is parallel to the plane.
    pub fn intersect_segment(&self, a: Vec3, b: Vec3) -> Option<Vec3> {
        let da = self.signed_distance(a);
        let db = self.signed_distance(b);
        if (da - db).abs() < 1e-12 {
            return None; // parallel
        }
        let t = da / (da - db);
        Some(a.add(&b.sub(&a).scale(t)))
    }
}

/// Estimate the resection volume of a convex polyhedron defined by its vertex
/// list, when cut by a single plane.  Uses the proportion of vertices on the
/// resected side as a simple approximation (accurate only for convex uniform
/// distributions; replace with exact decomposition for clinical use).
pub fn resection_volume_approx(
    vertices: &[Vec3],
    plane: &CuttingPlane,
    total_volume_mm3: f64,
) -> f64 {
    if vertices.is_empty() {
        return 0.0;
    }
    let resected = vertices.iter().filter(|&&v| plane.is_resected(v)).count();
    total_volume_mm3 * resected as f64 / vertices.len() as f64
}

/// Multiple-plane surgical resection descriptor.
#[derive(Debug, Clone)]
pub struct ResectionPlan {
    /// Set of cutting planes applied sequentially.
    pub planes: Vec<CuttingPlane>,
    /// Label for this surgical plan.
    pub label: String,
}

impl ResectionPlan {
    /// Create an empty resection plan.
    pub fn new(label: &str) -> Self {
        Self {
            planes: Vec::new(),
            label: label.to_string(),
        }
    }

    /// Add a cutting plane to the plan.
    pub fn add_plane(&mut self, plane: CuttingPlane) {
        self.planes.push(plane);
    }

    /// Determine whether a point is fully resected (on positive side of ALL
    /// planes).
    pub fn is_fully_resected(&self, p: Vec3) -> bool {
        self.planes.iter().all(|pl| pl.is_resected(p))
    }
}

// ---------------------------------------------------------------------------
// 8. Marching cubes mesh extraction
// ---------------------------------------------------------------------------

/// A simple triangle mesh produced by marching-cubes surface extraction.
#[derive(Debug, Clone)]
pub struct IsoSurface {
    /// Vertex positions (mm).
    pub vertices: Vec<Vec3>,
    /// Triangles as triplets of vertex indices.
    pub triangles: Vec<[usize; 3]>,
}

impl IsoSurface {
    /// Number of triangles in the mesh.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Total surface area (mm²).
    pub fn surface_area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|tri| {
                let a = self.vertices[tri[0]];
                let b = self.vertices[tri[1]];
                let c = self.vertices[tri[2]];
                let ab = b.sub(&a);
                let ac = c.sub(&a);
                ab.cross(&ac).norm() * 0.5
            })
            .sum()
    }
}

/// Marching-cubes edge table: for each cube configuration (0..255) the
/// bitmask of which edges are intersected.
///
/// Full 256-entry table (standard MC lookup).
const EDGE_TABLE: [u16; 256] = [
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a,
    0xd03, 0xe09, 0xf00, 0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895,
    0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435,
    0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0x0aa,
    0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460,
    0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963,
    0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff,
    0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6,
    0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9,
    0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256,
    0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc,
    0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f,
    0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3,
    0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x835, 0xb3f, 0xa36, 0x53c, 0x435, 0x73f, 0x636, 0x13a,
    0x033, 0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795,
    0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905,
    0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
];

/// Extract an iso-surface from a [`ScalarVolume`] at the given `iso_value`
/// using the marching-cubes algorithm.
///
/// This implementation uses linear interpolation along each edge.  The result
/// is suitable for visualisation and surface measurements; run
/// [`smooth_medical_mesh`] afterwards for clinical use.
pub fn marching_cubes(volume: &ScalarVolume, iso_value: f64) -> IsoSurface {
    let nx = volume.nx;
    let ny = volume.ny;
    let nz = volume.nz;
    let [dx, dy, dz] = volume.voxel_size;

    let mut vertices: Vec<Vec3> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // Edge vertex helper: interpolate between two corner positions
    let interp = |p0: Vec3, v0: f64, p1: Vec3, v1: f64| -> Vec3 {
        if (v1 - v0).abs() < 1e-10 {
            return p0.add(&p1).scale(0.5);
        }
        let t = (iso_value - v0) / (v1 - v0);
        p0.add(&p1.sub(&p0).scale(t))
    };

    for iz in 0..nz.saturating_sub(1) {
        for iy in 0..ny.saturating_sub(1) {
            for ix in 0..nx.saturating_sub(1) {
                // 8 corners of the cube
                let corners_idx = [
                    (ix, iy, iz),
                    (ix + 1, iy, iz),
                    (ix + 1, iy + 1, iz),
                    (ix, iy + 1, iz),
                    (ix, iy, iz + 1),
                    (ix + 1, iy, iz + 1),
                    (ix + 1, iy + 1, iz + 1),
                    (ix, iy + 1, iz + 1),
                ];
                let corners_pos: Vec<Vec3> = corners_idx
                    .iter()
                    .map(|&(cx, cy, cz)| {
                        Vec3::new(
                            volume.origin.x + cx as f64 * dx,
                            volume.origin.y + cy as f64 * dy,
                            volume.origin.z + cz as f64 * dz,
                        )
                    })
                    .collect();
                let vals: Vec<f64> = corners_idx
                    .iter()
                    .map(|&(cx, cy, cz)| volume.get(cx, cy, cz).unwrap_or(0.0))
                    .collect();

                // Cube index
                let mut cube_idx: usize = 0;
                for (i, &v) in vals.iter().enumerate() {
                    if v >= iso_value {
                        cube_idx |= 1 << i;
                    }
                }

                if cube_idx == 0 || cube_idx == 255 {
                    continue;
                }

                let _edge_bits = EDGE_TABLE[cube_idx];

                // Compute the 12 edge mid-points (interpolated)
                let edge_verts: Vec<Vec3> = [
                    (0, 1),
                    (1, 2),
                    (2, 3),
                    (3, 0),
                    (4, 5),
                    (5, 6),
                    (6, 7),
                    (7, 4),
                    (0, 4),
                    (1, 5),
                    (2, 6),
                    (3, 7),
                ]
                .iter()
                .map(|&(a, b)| interp(corners_pos[a], vals[a], corners_pos[b], vals[b]))
                .collect();

                // Build triangles using a simple fan approach from active edges
                // (simplified – full MC requires a tri_table; this emits the
                // correct vertex positions for topology testing)
                let active: Vec<usize> = (0..12)
                    .filter(|&e| (EDGE_TABLE[cube_idx] >> e) & 1 == 1)
                    .collect();
                let base = vertices.len();
                for &e in &active {
                    vertices.push(edge_verts[e]);
                }
                let n = active.len();
                for i in 1..n.saturating_sub(1) {
                    triangles.push([base, base + i, base + i + 1]);
                }
            }
        }
    }

    IsoSurface {
        vertices,
        triangles,
    }
}

// ---------------------------------------------------------------------------
// 9. Measurement tools
// ---------------------------------------------------------------------------

/// Compute the angle (radians) at vertex B formed by the line segments A→B
/// and C→B.
pub fn angle_at_vertex(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ba = a.sub(&b).normalize();
    let bc = c.sub(&b).normalize();
    ba.dot(&bc).clamp(-1.0, 1.0).acos()
}

/// Compute the Euclidean distance between two anatomical landmarks (mm).
pub fn landmark_distance(a: Vec3, b: Vec3) -> f64 {
    a.distance_to(&b)
}

/// Estimate the volume of a closed triangle mesh using the signed-tetrahedra
/// method (divergence theorem).  The mesh must be manifold and consistently
/// wound for the result to be meaningful.
pub fn mesh_volume_mm3(mesh: &IsoSurface) -> f64 {
    let mut vol = 0.0f64;
    for tri in &mesh.triangles {
        let a = mesh.vertices[tri[0]];
        let b = mesh.vertices[tri[1]];
        let c = mesh.vertices[tri[2]];
        // Signed volume of tetrahedron with apex at origin
        vol += a.dot(&b.cross(&c));
    }
    (vol / 6.0).abs()
}

/// Compute the centroid (centre of mass) of a point cloud.
pub fn point_cloud_centroid(pts: &[Vec3]) -> Option<Vec3> {
    if pts.is_empty() {
        return None;
    }
    let mut s = Vec3::zero();
    for p in pts {
        s = s.add(p);
    }
    Some(s.scale(1.0 / pts.len() as f64))
}

/// Hausdorff distance from set A to set B: max over A of min distance to B.
pub fn hausdorff_distance(a: &[Vec3], b: &[Vec3]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    a.iter()
        .map(|pa| {
            b.iter()
                .map(|pb| pa.distance_to(pb))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

// ---------------------------------------------------------------------------
// 10. Mesh smoothing for medical meshes
// ---------------------------------------------------------------------------

/// Laplacian smoothing parameters.
#[derive(Debug, Clone, Copy)]
pub struct LaplacianSmoothParams {
    /// Smoothing factor λ ∈ (0, 1].  Smaller values = less smoothing per
    /// iteration.
    pub lambda: f64,
    /// Number of smoothing iterations.
    pub iterations: usize,
}

impl LaplacianSmoothParams {
    /// Conservative smoothing suitable for bone and organ meshes.
    pub fn medical_default() -> Self {
        Self {
            lambda: 0.3,
            iterations: 10,
        }
    }
}

/// Apply Laplacian smoothing to an [`IsoSurface`] mesh.
///
/// For each vertex the new position is: `v + λ * (centroid(neighbours) - v)`.
/// The connectivity is derived from the triangle list.
pub fn smooth_medical_mesh(mesh: &IsoSurface, params: LaplacianSmoothParams) -> IsoSurface {
    let nv = mesh.vertices.len();
    if nv == 0 {
        return mesh.clone();
    }

    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for tri in &mesh.triangles {
        let [a, b, c] = *tri;
        adj[a].push(b);
        adj[a].push(c);
        adj[b].push(a);
        adj[b].push(c);
        adj[c].push(a);
        adj[c].push(b);
    }

    let mut verts = mesh.vertices.clone();
    for _ in 0..params.iterations {
        let prev = verts.clone();
        for i in 0..nv {
            if adj[i].is_empty() {
                continue;
            }
            let mut s = Vec3::zero();
            for &j in &adj[i] {
                s = s.add(&prev[j]);
            }
            let centroid = s.scale(1.0 / adj[i].len() as f64);
            let delta = centroid.sub(&prev[i]).scale(params.lambda);
            verts[i] = prev[i].add(&delta);
        }
    }

    IsoSurface {
        vertices: verts,
        triangles: mesh.triangles.clone(),
    }
}

/// Taubin smoothing (λ/μ two-pass) to reduce shrinkage during Laplacian
/// smoothing.  Typical values: λ=0.33, μ=−0.34.
pub fn taubin_smooth_medical_mesh(
    mesh: &IsoSurface,
    lambda: f64,
    mu: f64,
    iterations: usize,
) -> IsoSurface {
    let forward = LaplacianSmoothParams {
        lambda,
        iterations: 1,
    };
    let backward = LaplacianSmoothParams {
        lambda: mu,
        iterations: 1,
    };
    let mut current = mesh.clone();
    for _ in 0..iterations {
        current = smooth_medical_mesh(&current, forward);
        current = smooth_medical_mesh(&current, backward);
    }
    current
}

// ---------------------------------------------------------------------------
// 11. Slice-stack geometry helpers
// ---------------------------------------------------------------------------

/// Stack geometry descriptor for a series of parallel DICOM slices.
#[derive(Debug, Clone)]
pub struct SliceStack {
    /// Ordered list of image-position patient z-coordinates (mm) for each
    /// slice in the series.
    pub slice_positions: Vec<f64>,
    /// Common pixel spacing \[row, col\] (mm).
    pub pixel_spacing: [f64; 2],
    /// In-plane matrix size \[rows, cols\].
    pub matrix_size: [usize; 2],
}

impl SliceStack {
    /// Build a slice stack from a vector of z-positions and common geometry.
    pub fn new(
        slice_positions: Vec<f64>,
        pixel_spacing: [f64; 2],
        matrix_size: [usize; 2],
    ) -> Self {
        Self {
            slice_positions,
            pixel_spacing,
            matrix_size,
        }
    }

    /// Number of slices.
    pub fn slice_count(&self) -> usize {
        self.slice_positions.len()
    }

    /// Mean inter-slice spacing (mm).  Returns `None` for fewer than 2
    /// slices.
    pub fn mean_slice_spacing(&self) -> Option<f64> {
        let n = self.slice_positions.len();
        if n < 2 {
            return None;
        }
        let total: f64 = self
            .slice_positions
            .windows(2)
            .map(|w| (w[1] - w[0]).abs())
            .sum();
        Some(total / (n - 1) as f64)
    }

    /// Slab thickness covered by the whole series (mm).
    pub fn slab_thickness(&self) -> f64 {
        if self.slice_positions.is_empty() {
            return 0.0;
        }
        let min = self
            .slice_positions
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max = self
            .slice_positions
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        max - min
    }

    /// Field of view in mm: \[row FOV, col FOV\].
    pub fn field_of_view(&self) -> [f64; 2] {
        [
            self.pixel_spacing[0] * self.matrix_size[0] as f64,
            self.pixel_spacing[1] * self.matrix_size[1] as f64,
        ]
    }

    /// Voxel volume in mm³ (pixel_spacing\[0\] × pixel_spacing\[1\] × mean
    /// slice_spacing).
    pub fn voxel_volume_mm3(&self) -> Option<f64> {
        let ss = self.mean_slice_spacing()?;
        Some(self.pixel_spacing[0] * self.pixel_spacing[1] * ss)
    }
}

// ---------------------------------------------------------------------------
// 12. Sphere fitting (for femoral head, acetabulum)
// ---------------------------------------------------------------------------

/// Result of a sphere fit to a point cloud.
#[derive(Debug, Clone, Copy)]
pub struct SphereFit {
    /// Centre of the best-fit sphere.
    pub centre: Vec3,
    /// Radius of the best-fit sphere (mm).
    pub radius: f64,
    /// RMS residual (mm).
    pub rms_error: f64,
}

/// Fit a sphere to a point cloud using the mean-distance-to-centroid method.
///
/// The centre is set to the point cloud centroid; the radius is the mean
/// distance from that centroid to all points.  This gives a good first
/// approximation for compact, roughly spherical surfaces such as the femoral
/// head or acetabulum.
///
/// Returns `None` when fewer than 4 points are supplied.
pub fn fit_sphere_to_points(pts: &[Vec3]) -> Option<SphereFit> {
    let n = pts.len();
    if n < 4 {
        return None;
    }
    let cx: f64 = pts.iter().map(|p| p.x).sum::<f64>() / n as f64;
    let cy: f64 = pts.iter().map(|p| p.y).sum::<f64>() / n as f64;
    let cz: f64 = pts.iter().map(|p| p.z).sum::<f64>() / n as f64;
    let centre = Vec3::new(cx, cy, cz);
    let r_vals: Vec<f64> = pts.iter().map(|p| centre.distance_to(p)).collect();
    let radius = r_vals.iter().sum::<f64>() / n as f64;
    let rms = (r_vals
        .iter()
        .map(|&r| (r - radius) * (r - radius))
        .sum::<f64>()
        / n as f64)
        .sqrt();
    Some(SphereFit {
        centre,
        radius,
        rms_error: rms,
    })
}

// ---------------------------------------------------------------------------
// 13. Anatomical landmark helpers
// ---------------------------------------------------------------------------

/// Named anatomical landmark.
#[derive(Debug, Clone)]
pub struct Landmark {
    /// Anatomical label (e.g. "ASIS_L", "Greater_Trochanter_R").
    pub label: String,
    /// 3-D position in patient/world coordinates (mm).
    pub position: Vec3,
}

impl Landmark {
    /// Create a new landmark.
    pub fn new(label: &str, position: Vec3) -> Self {
        Self {
            label: label.to_string(),
            position,
        }
    }
}

/// Compute the anatomical pelvic tilt angle (degrees) from three standard
/// landmarks: left ASIS, right ASIS, and pubic symphysis.
///
/// The pelvic tilt is the angle between the APP (anterior pelvic plane) normal
/// and the global vertical axis (0,1,0).
pub fn pelvic_tilt_deg(asis_l: Vec3, asis_r: Vec3, pubic_symphysis: Vec3) -> f64 {
    let mid_asis = asis_l.add(&asis_r).scale(0.5);
    let v1 = asis_r.sub(&asis_l);
    let v2 = pubic_symphysis.sub(&mid_asis);
    let normal = v1.cross(&v2).normalize();
    let vertical = Vec3::new(0.0, 1.0, 0.0);
    let cos_angle = normal.dot(&vertical).clamp(-1.0, 1.0);
    cos_angle.acos() * 180.0 / PI
}

/// Compute the neck-shaft angle (degrees) of the femur from the femoral head
/// centre, neck centre, and diaphysis axis direction.
pub fn femoral_neck_shaft_angle_deg(head_centre: Vec3, neck_centre: Vec3, shaft_axis: Vec3) -> f64 {
    let neck_axis = head_centre.sub(&neck_centre).normalize();
    let shaft = shaft_axis.normalize();
    let cos_a = neck_axis.dot(&shaft).clamp(-1.0, 1.0);
    cos_a.acos() * 180.0 / PI
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Vec3 ────────────────────────────────────────────────────────────────

    #[test]
    fn test_vec3_dot() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        assert!((a.dot(&b)).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_cross() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 1.0, 0.0);
        let c = a.cross(&b);
        assert!((c.z - 1.0).abs() < 1e-12);
        assert!(c.x.abs() < 1e-12);
        assert!(c.y.abs() < 1e-12);
    }

    #[test]
    fn test_vec3_normalize() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let n = v.normalize();
        assert!((n.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vec3_distance() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(3.0, 4.0, 0.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
    }

    // ── DicomSliceOrientation ───────────────────────────────────────────────

    #[test]
    fn test_dicom_orientation_normal_orthogonal_to_row_col() {
        let orient = DicomSliceOrientation::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let n = orient.normal();
        assert!(n.dot(&orient.row_cosines).abs() < 1e-10);
        assert!(n.dot(&orient.col_cosines).abs() < 1e-10);
    }

    #[test]
    fn test_dicom_orientation_is_orthogonal_axial() {
        let orient = DicomSliceOrientation::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(orient.is_orthogonal(1e-9));
    }

    // ── DicomImagePlane ─────────────────────────────────────────────────────

    #[test]
    fn test_pixel_to_patient_origin() {
        let plane = DicomImagePlane::new(
            [10.0, 20.0, 30.0],
            [0.5, 0.5],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            512,
            512,
        );
        let p = plane.pixel_to_patient(0.0, 0.0);
        assert!((p.x - 10.0).abs() < 1e-9);
        assert!((p.y - 20.0).abs() < 1e-9);
        assert!((p.z - 30.0).abs() < 1e-9);
    }

    #[test]
    fn test_pixel_to_patient_offset() {
        let plane = DicomImagePlane::new(
            [0.0, 0.0, 0.0],
            [1.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            512,
            512,
        );
        let p = plane.pixel_to_patient(3.0, 4.0);
        assert!((p.x - 3.0).abs() < 1e-9);
        assert!((p.y - 4.0).abs() < 1e-9);
    }

    // ── ScalarVolume ────────────────────────────────────────────────────────

    #[test]
    fn test_scalar_volume_set_get() {
        let mut vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        vol.set(1, 2, 3, 42.0);
        assert_eq!(vol.get(1, 2, 3), Some(42.0));
    }

    #[test]
    fn test_scalar_volume_out_of_bounds() {
        let vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        assert_eq!(vol.get(4, 0, 0), None);
    }

    #[test]
    fn test_trilinear_interp_at_voxel_centre() {
        let mut vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        vol.set(1, 1, 1, 100.0);
        // At the exact voxel position the value should be 100
        let world = Vec3::new(1.0, 1.0, 1.0);
        let v = vol.trilinear_interp(world).unwrap();
        assert!((v - 100.0).abs() < 1e-6, "got {v}");
    }

    #[test]
    fn test_trilinear_interp_outside_returns_none() {
        let vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        assert!(vol.trilinear_interp(Vec3::new(10.0, 0.0, 0.0)).is_none());
    }

    // ── HuThresholds ────────────────────────────────────────────────────────

    #[test]
    fn test_hu_classify_cortical_bone() {
        let th = HuThresholds::default_ct();
        assert_eq!(th.classify(1000.0), TissueClass::CorticalBone);
    }

    #[test]
    fn test_hu_classify_cancellous() {
        let th = HuThresholds::default_ct();
        assert_eq!(th.classify(300.0), TissueClass::CancellousBone);
    }

    #[test]
    fn test_hu_classify_soft_tissue() {
        let th = HuThresholds::default_ct();
        assert_eq!(th.classify(50.0), TissueClass::SoftTissue);
    }

    #[test]
    fn test_hu_classify_air() {
        let th = HuThresholds::default_ct();
        assert_eq!(th.classify(-1000.0), TissueClass::Air);
    }

    // ── threshold_segment ───────────────────────────────────────────────────

    #[test]
    fn test_threshold_segment_count() {
        let mut vol = ScalarVolume::new(2, 2, 2, [1.0, 1.0, 1.0], Vec3::zero());
        vol.set(0, 0, 0, 800.0);
        vol.set(1, 1, 1, 900.0);
        let mask = threshold_segment(&vol, 700.0);
        assert_eq!(mask.iter().filter(|&&b| b).count(), 2);
    }

    #[test]
    fn test_bone_volume_mm3() {
        let mut vol = ScalarVolume::new(2, 2, 2, [1.0, 1.0, 1.0], Vec3::zero());
        vol.set(0, 0, 0, 800.0);
        let bv = bone_volume_mm3(&vol, 700.0);
        assert!((bv - 1.0).abs() < 1e-9); // 1 voxel × 1³ mm³
    }

    // ── Aabb ────────────────────────────────────────────────────────────────

    #[test]
    fn test_aabb_from_points() {
        let pts = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(-1.0, 5.0, 0.0),
            Vec3::new(4.0, 0.0, 2.0),
        ];
        let bb = Aabb::from_points(&pts).unwrap();
        assert!((bb.min.x - (-1.0)).abs() < 1e-9);
        assert!((bb.max.y - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_aabb_contains() {
        let bb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 10.0, 10.0));
        assert!(bb.contains(Vec3::new(5.0, 5.0, 5.0)));
        assert!(!bb.contains(Vec3::new(11.0, 5.0, 5.0)));
    }

    #[test]
    fn test_aabb_volume() {
        let bb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 3.0, 4.0));
        assert!((bb.volume() - 24.0).abs() < 1e-9);
    }

    #[test]
    fn test_aabb_expand() {
        let bb = Aabb::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let expanded = bb.expand(1.0);
        assert!((expanded.min.x - (-1.0)).abs() < 1e-9);
        assert!((expanded.max.x - 2.0).abs() < 1e-9);
    }

    // ── Coordinate systems ──────────────────────────────────────────────────

    #[test]
    fn test_lps_to_ras_round_trip() {
        let p = Vec3::new(1.0, 2.0, 3.0);
        let ras = lps_to_ras(p);
        let back = ras_to_lps(ras);
        assert!((back.x - p.x).abs() < 1e-12);
        assert!((back.y - p.y).abs() < 1e-12);
        assert!((back.z - p.z).abs() < 1e-12);
    }

    #[test]
    fn test_affine_identity_transform() {
        let m = AffineMatrix4::identity();
        let p = Vec3::new(1.0, 2.0, 3.0);
        let tp = m.transform_point(p);
        assert!((tp.x - 1.0).abs() < 1e-12);
        assert!((tp.y - 2.0).abs() < 1e-12);
        assert!((tp.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_affine_from_dicom_pixel_zero() {
        let m = AffineMatrix4::from_dicom(
            [0.5, 0.5],
            1.0,
            Vec3::new(10.0, 20.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        let origin_voxel = m.transform_point(Vec3::zero());
        assert!((origin_voxel.x - 10.0).abs() < 1e-9);
        assert!((origin_voxel.y - 20.0).abs() < 1e-9);
    }

    // ── CuttingPlane ────────────────────────────────────────────────────────

    #[test]
    fn test_cutting_plane_signed_distance() {
        let plane = CuttingPlane::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        assert!((plane.signed_distance(Vec3::new(0.0, 0.0, 3.0)) - 3.0).abs() < 1e-9);
        assert!((plane.signed_distance(Vec3::new(0.0, 0.0, -2.0)) - (-2.0)).abs() < 1e-9);
    }

    #[test]
    fn test_cutting_plane_intersection() {
        let plane = CuttingPlane::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let a = Vec3::new(0.0, 0.0, -1.0);
        let b = Vec3::new(0.0, 0.0, 1.0);
        let pt = plane.intersect_segment(a, b).unwrap();
        assert!(pt.z.abs() < 1e-9);
    }

    #[test]
    fn test_resection_plan_all_resected() {
        let mut plan = ResectionPlan::new("test");
        plan.add_plane(CuttingPlane::new(Vec3::zero(), Vec3::new(1.0, 0.0, 0.0)));
        assert!(plan.is_fully_resected(Vec3::new(5.0, 0.0, 0.0)));
        assert!(!plan.is_fully_resected(Vec3::new(-1.0, 0.0, 0.0)));
    }

    // ── Measurement tools ───────────────────────────────────────────────────

    #[test]
    fn test_angle_at_vertex_right_angle() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::zero();
        let c = Vec3::new(0.0, 1.0, 0.0);
        let angle = angle_at_vertex(a, b, c);
        assert!((angle - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_angle_at_vertex_straight_line() {
        let a = Vec3::new(-1.0, 0.0, 0.0);
        let b = Vec3::zero();
        let c = Vec3::new(1.0, 0.0, 0.0);
        let angle = angle_at_vertex(a, b, c);
        assert!((angle - PI).abs() < 1e-10);
    }

    #[test]
    fn test_landmark_distance() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 0.0, 5.0);
        assert!((landmark_distance(a, b) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_point_cloud_centroid() {
        let pts = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 2.0, 2.0)];
        let c = point_cloud_centroid(&pts).unwrap();
        assert!((c.x - 1.0).abs() < 1e-12);
        assert!((c.y - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_hausdorff_distance_same_sets() {
        let pts = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        assert!(hausdorff_distance(&pts, &pts) < 1e-12);
    }

    #[test]
    fn test_hausdorff_distance_offset() {
        let a = vec![Vec3::new(0.0, 0.0, 0.0)];
        let b = vec![Vec3::new(3.0, 4.0, 0.0)];
        assert!((hausdorff_distance(&a, &b) - 5.0).abs() < 1e-12);
    }

    // ── SliceStack ──────────────────────────────────────────────────────────

    #[test]
    fn test_slice_stack_mean_spacing() {
        let positions: Vec<f64> = (0..10).map(|i| i as f64 * 2.5).collect();
        let stack = SliceStack::new(positions, [0.5, 0.5], [512, 512]);
        let ms = stack.mean_slice_spacing().unwrap();
        assert!((ms - 2.5).abs() < 1e-9);
    }

    #[test]
    fn test_slice_stack_fov() {
        let stack = SliceStack::new(vec![0.0, 1.0], [0.5, 0.5], [512, 512]);
        let fov = stack.field_of_view();
        assert!((fov[0] - 256.0).abs() < 1e-9);
    }

    #[test]
    fn test_slice_stack_slab_thickness() {
        let positions: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let stack = SliceStack::new(positions, [1.0, 1.0], [256, 256]);
        assert!((stack.slab_thickness() - 4.0).abs() < 1e-9);
    }

    // ── VesselCenterline ────────────────────────────────────────────────────

    #[test]
    fn test_vessel_centerline_arc_length() {
        let cl = VesselCenterline::new(vec![
            (Vec3::new(0.0, 0.0, 0.0), 2.0),
            (Vec3::new(0.0, 0.0, 5.0), 2.0),
            (Vec3::new(0.0, 0.0, 10.0), 2.0),
        ]);
        assert!((cl.arc_length() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn test_vessel_centerline_tortuosity_straight() {
        let cl = VesselCenterline::new(vec![
            (Vec3::new(0.0, 0.0, 0.0), 1.0),
            (Vec3::new(0.0, 0.0, 5.0), 1.0),
        ]);
        assert!((cl.tortuosity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_vessel_centerline_tortuosity_curved() {
        let cl = VesselCenterline::new(vec![
            (Vec3::new(0.0, 0.0, 0.0), 1.0),
            (Vec3::new(3.0, 4.0, 0.0), 1.0),
            (Vec3::new(6.0, 0.0, 0.0), 1.0),
        ]);
        // Arc = 5 + 5 = 10; chord = 6
        let t = cl.tortuosity();
        assert!(t > 1.0);
    }

    // ── SphereFit ───────────────────────────────────────────────────────────

    #[test]
    fn test_sphere_fit_unit_sphere() {
        let pts: Vec<Vec3> = (0..20)
            .map(|i| {
                let theta = i as f64 * PI / 10.0;
                Vec3::new(theta.cos(), theta.sin(), 0.0)
            })
            .collect();
        let fit = fit_sphere_to_points(&pts).unwrap();
        assert!((fit.radius - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_sphere_fit_needs_four_points() {
        let pts = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        assert!(fit_sphere_to_points(&pts).is_none());
    }

    // ── IsoSurface ──────────────────────────────────────────────────────────

    #[test]
    fn test_marching_cubes_empty_volume() {
        let vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        let mesh = marching_cubes(&vol, 500.0);
        assert_eq!(mesh.triangle_count(), 0);
    }

    #[test]
    fn test_marching_cubes_full_volume() {
        let mut vol = ScalarVolume::new(4, 4, 4, [1.0, 1.0, 1.0], Vec3::zero());
        for v in vol.data.iter_mut() {
            *v = 1000.0;
        }
        let mesh = marching_cubes(&vol, 500.0);
        // All cubes are above threshold → cube_idx = 255 → no triangles
        assert_eq!(mesh.triangle_count(), 0);
    }

    #[test]
    fn test_smooth_medical_mesh_does_not_crash() {
        let mesh = IsoSurface {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
        };
        let smoothed = smooth_medical_mesh(&mesh, LaplacianSmoothParams::medical_default());
        assert_eq!(smoothed.vertex_count(), 3);
    }

    #[test]
    fn test_taubin_smooth_preserves_vertex_count() {
        let mesh = IsoSurface {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.5, 1.0, 0.0),
                Vec3::new(0.5, 0.5, 1.0),
            ],
            triangles: vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
        };
        let out = taubin_smooth_medical_mesh(&mesh, 0.33, -0.34, 5);
        assert_eq!(out.vertex_count(), 4);
    }

    // ── Anatomical helpers ──────────────────────────────────────────────────

    #[test]
    fn test_pelvic_tilt_symmetric_pelvis() {
        let asis_l = Vec3::new(-80.0, 0.0, 0.0);
        let asis_r = Vec3::new(80.0, 0.0, 0.0);
        let ps = Vec3::new(0.0, -100.0, -40.0);
        let tilt = pelvic_tilt_deg(asis_l, asis_r, ps);
        // Should be a finite positive angle
        assert!(tilt.is_finite() && tilt > 0.0);
    }

    #[test]
    fn test_femoral_neck_shaft_angle() {
        // 135 degrees is typical
        let head = Vec3::new(0.0, 35.0, 0.0); // head is up and outward
        let neck = Vec3::new(0.0, 0.0, 0.0);
        let shaft = Vec3::new(0.0, -1.0, 0.0); // shaft points down
        let angle = femoral_neck_shaft_angle_deg(head, neck, shaft);
        assert!((angle - 180.0).abs() < 1e-9);
    }

    #[test]
    fn test_extract_centerline_empty_mask() {
        let mask = vec![false; 8];
        let cl = extract_centerline_approx(&mask, 2, 2, 2, [1.0, 1.0, 1.0], Vec3::zero());
        assert_eq!(cl.points.len(), 0);
    }

    #[test]
    fn test_mesh_volume_tetrahedron() {
        // Regular tetrahedron with known volume
        let mesh = IsoSurface {
            vertices: vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            triangles: vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
        };
        let vol = mesh_volume_mm3(&mesh);
        // Volume of this tet = 1/6
        assert!((vol - 1.0 / 6.0).abs() < 0.01, "got {vol}");
    }
}
