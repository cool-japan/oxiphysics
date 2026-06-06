//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

use super::functions::{
    cross3_raw, dot3_raw, negate3_raw, obb_obb_test, project_obb_onto_axis, scale3_raw, sub3_raw,
};

/// Contact information produced by an OBB-OBB SAT query.
#[derive(Debug, Clone)]
pub struct ObbContact {
    /// Penetration depth (always positive when contact occurs).
    pub depth: f64,
    /// Contact normal pointing from OBB B toward OBB A.
    pub normal: Vec3,
    /// World-space contact point (midpoint of witness features).
    pub point: Vec3,
}
/// SAT-based OBB vs OBB collision detector.
pub struct ObbSat;
impl ObbSat {
    /// Project both OBBs onto `axis` and return the overlap amount.
    ///
    /// Returns `None` if the shapes are separated along this axis.
    /// A return value of `Some(overlap)` means they penetrate by `overlap` on this axis.
    pub fn test_separation_axis(obb_a: &ObbShape, obb_b: &ObbShape, axis: &Vec3) -> Option<f64> {
        let len_sq = axis.norm_squared();
        if len_sq < 1e-12 {
            return Some(0.0);
        }
        let axis_norm = axis / len_sq.sqrt();
        let ra = (obb_a.half_extents.x * obb_a.axes[0].dot(&axis_norm)).abs()
            + (obb_a.half_extents.y * obb_a.axes[1].dot(&axis_norm)).abs()
            + (obb_a.half_extents.z * obb_a.axes[2].dot(&axis_norm)).abs();
        let rb = (obb_b.half_extents.x * obb_b.axes[0].dot(&axis_norm)).abs()
            + (obb_b.half_extents.y * obb_b.axes[1].dot(&axis_norm)).abs()
            + (obb_b.half_extents.z * obb_b.axes[2].dot(&axis_norm)).abs();
        let center_diff = obb_b.center - obb_a.center;
        let d = center_diff.dot(&axis_norm).abs();
        let sep = d - (ra + rb);
        if sep > 0.0 { None } else { Some(-sep) }
    }
    /// Test OBB A against OBB B using SAT.
    ///
    /// Returns `None` if separated; returns `Some(ObbContact)` with the
    /// minimum-overlap contact if penetrating.
    pub fn query(obb_a: &ObbShape, obb_b: &ObbShape) -> Option<ObbContact> {
        let center_diff = obb_b.center - obb_a.center;
        let mut min_depth = f64::INFINITY;
        let mut best_axis = Vec3::new(0.0, 1.0, 0.0);
        for i in 0..3 {
            let axis = obb_a.axes[i];
            let overlap = Self::test_separation_axis(obb_a, obb_b, &axis)?;
            if overlap < min_depth {
                min_depth = overlap;
                best_axis = if center_diff.dot(&axis) >= 0.0 {
                    axis
                } else {
                    -axis
                };
            }
        }
        for j in 0..3 {
            let axis = obb_b.axes[j];
            let overlap = Self::test_separation_axis(obb_a, obb_b, &axis)?;
            if overlap < min_depth {
                min_depth = overlap;
                best_axis = if center_diff.dot(&axis) >= 0.0 {
                    axis
                } else {
                    -axis
                };
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let cross = obb_a.axes[i].cross(&obb_b.axes[j]);
                let cross_len_sq = cross.norm_squared();
                if cross_len_sq < 1e-12 {
                    continue;
                }
                let overlap = Self::test_separation_axis(obb_a, obb_b, &cross)?;
                if overlap < min_depth {
                    min_depth = overlap;
                    let axis = cross / cross_len_sq.sqrt();
                    best_axis = if center_diff.dot(&axis) >= 0.0 {
                        axis
                    } else {
                        -axis
                    };
                }
            }
        }
        let point = Self::find_contact_point(obb_a, obb_b, &best_axis, min_depth);
        Some(ObbContact {
            depth: min_depth,
            normal: best_axis,
            point,
        })
    }
    /// Compute an approximate contact point as the midpoint of the witness
    /// features on each OBB projected along the contact normal.
    pub fn find_contact_point(
        obb_a: &ObbShape,
        obb_b: &ObbShape,
        normal: &Vec3,
        depth: f64,
    ) -> Vec3 {
        let support_a = Self::support_point(obb_a, normal);
        let support_b = Self::support_point(obb_b, &-normal);
        (support_a + support_b) * 0.5 + normal * (depth * 0.5)
    }
    /// Compute the support point of `obb` along `direction` (world space).
    fn support_point(obb: &ObbShape, direction: &Vec3) -> Vec3 {
        let mut p = obb.center;
        for i in 0..3 {
            let dot = obb.axes[i].dot(direction);
            let he = if i == 0 {
                obb.half_extents.x
            } else if i == 1 {
                obb.half_extents.y
            } else {
                obb.half_extents.z
            };
            p += obb.axes[i] * (he * dot.signum());
        }
        p
    }
}
/// Result of a SAT collision test.
#[derive(Debug, Clone)]
pub struct SatResult {
    /// Penetration depth (positive means overlap).
    pub penetration_depth: f64,
    /// Contact normal (unit vector).
    pub contact_normal: [f64; 3],
    /// Approximate contact point in world space.
    pub contact_point: [f64; 3],
}
/// An oriented bounding box using raw `[f64;3]` arrays.
#[derive(Debug, Clone)]
pub struct Obb {
    /// World-space center.
    pub center: [f64; 3],
    /// Half-extents along local axes.
    pub half_extents: [f64; 3],
    /// Rotation matrix (columns are the local axes in world space).
    pub rotation: [[f64; 3]; 3],
}
impl Obb {
    /// Create a new OBB.
    pub fn new(center: [f64; 3], half_extents: [f64; 3], rotation: [[f64; 3]; 3]) -> Self {
        Self {
            center,
            half_extents,
            rotation,
        }
    }
    /// Create an axis-aligned OBB (identity rotation).
    pub fn axis_aligned(center: [f64; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            half_extents,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }
    /// Create an OBB from a center, local axes (columns of rotation matrix), and half-extents.
    pub fn from_center_axes(center: [f64; 3], axes: [[f64; 3]; 3], half_extents: [f64; 3]) -> Self {
        Self {
            center,
            half_extents,
            rotation: axes,
        }
    }
    /// Project this OBB onto `axis` and return the `(min, max)` interval.
    pub fn project_onto_axis(&self, axis: [f64; 3]) -> (f64, f64) {
        project_obb_onto_axis(self, axis)
    }
}
/// High-level SAT test facade with explicit edge-edge and contact-normal APIs.
///
/// All methods operate on the raw `Obb` type and return plain `[f64; 3]` data
/// so they are usable independently of the nalgebra-backed types in the rest
/// of the crate.
pub struct ObbSatTest;
impl ObbSatTest {
    /// Test the *edge-edge* separating axes between two OBBs.
    ///
    /// In the standard 15-axis SAT for OBB-OBB, 9 of the axes come from
    /// cross products of edge-direction pairs (`a.rotation[i] × b.rotation[j]`).
    /// This method tests only those 9 axes, ignoring the 6 face-normal axes.
    ///
    /// Returns the minimum overlap among all edge-edge axes, or `None` if any
    /// edge-edge axis separates the OBBs.
    ///
    /// A `Some(overlap)` with `overlap == 0.0` from a degenerate (near-parallel)
    /// cross product should be treated as a face-contact; the caller is responsible
    /// for choosing between face and edge contacts.
    pub fn test_edge_edge(a: &Obb, b: &Obb) -> Option<f64> {
        let diff = sub3_raw(b.center, a.center);
        let mut min_overlap = f64::INFINITY;
        for i in 0..3 {
            for j in 0..3 {
                let cross = cross3_raw(a.rotation[i], b.rotation[j]);
                let len_sq = dot3_raw(cross, cross);
                if len_sq < 1e-12 {
                    min_overlap = min_overlap.min(0.0);
                    continue;
                }
                let inv = 1.0 / len_sq.sqrt();
                let axis = scale3_raw(cross, inv);
                let (min_a, max_a) = project_obb_onto_axis(a, axis);
                let (min_b, max_b) = project_obb_onto_axis(b, axis);
                let overlap = max_a.min(max_b) - min_a.max(min_b);
                if overlap < 0.0 {
                    return None;
                }
                let oriented = if dot3_raw(diff, axis) >= 0.0 {
                    axis
                } else {
                    negate3_raw(axis)
                };
                let _ = oriented;
                min_overlap = min_overlap.min(overlap);
            }
        }
        Some(min_overlap)
    }
    /// Compute the contact normal from the minimum-penetration SAT axis.
    ///
    /// Tests all 15 axes (3 face normals of A, 3 of B, 9 edge cross products)
    /// and returns the axis with the smallest overlap as the contact normal,
    /// oriented from B toward A (consistent with impulse-resolution convention).
    ///
    /// Returns `None` if any axis separates the OBBs (no contact).
    pub fn compute_contact_normal(a: &Obb, b: &Obb) -> Option<[f64; 3]> {
        let diff = sub3_raw(b.center, a.center);
        let mut min_depth = f64::INFINITY;
        let mut best_axis = [0.0_f64, 1.0, 0.0];
        let mut test = |axis: [f64; 3]| -> bool {
            let len_sq = dot3_raw(axis, axis);
            if len_sq < 1e-12 {
                return true;
            }
            let inv = 1.0 / len_sq.sqrt();
            let n = scale3_raw(axis, inv);
            let (min_a, max_a) = project_obb_onto_axis(a, n);
            let (min_b, max_b) = project_obb_onto_axis(b, n);
            let overlap = max_a.min(max_b) - min_a.max(min_b);
            if overlap < 0.0 {
                return false;
            }
            if overlap < min_depth {
                min_depth = overlap;
                best_axis = if dot3_raw(diff, n) >= 0.0 {
                    n
                } else {
                    negate3_raw(n)
                };
            }
            true
        };
        for i in 0..3 {
            if !test(a.rotation[i]) {
                return None;
            }
        }
        for i in 0..3 {
            if !test(b.rotation[i]) {
                return None;
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let cross = cross3_raw(a.rotation[i], b.rotation[j]);
                if !test(cross) {
                    return None;
                }
            }
        }
        Some(best_axis)
    }
    /// Compute a full contact result (normal + depth + point) using the SAT.
    ///
    /// Convenience wrapper around `compute_contact_normal` that also computes
    /// penetration depth and an approximate contact point.
    pub fn contact(a: &Obb, b: &Obb) -> Option<SatResult> {
        obb_obb_test(a, b)
    }
}
/// An oriented bounding box defined in world space.
#[derive(Debug, Clone)]
pub struct ObbShape {
    /// Half-widths along each local axis.
    pub half_extents: Vec3,
    /// World-space center of the OBB.
    pub center: Vec3,
    /// World-space orientation axes (columns of the rotation matrix).
    pub axes: [Vec3; 3],
}
impl ObbShape {
    /// Create a new OBB.
    pub fn new(half_extents: Vec3, center: Vec3, axes: [Vec3; 3]) -> Self {
        Self {
            half_extents,
            center,
            axes,
        }
    }
    /// Create an axis-aligned OBB (identity orientation).
    pub fn axis_aligned(half_extents: Vec3, center: Vec3) -> Self {
        Self {
            half_extents,
            center,
            axes: [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
        }
    }
}
/// Whether the contact normal is primarily from a face or an edge-edge interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactFeatureType {
    /// Contact normal is aligned with a face normal (face-face or face-edge).
    FaceContact,
    /// Contact normal is a cross product of edges.
    EdgeEdgeContact,
}
