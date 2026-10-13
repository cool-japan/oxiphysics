//! Box-box contact manifold generation.
//!
//! [`BoxBoxManifold`] runs the full 15-axis separating-axis test (the three
//! face normals of each box and the nine edge-edge cross products) and then
//! builds the manifold for the axis of least penetration:
//!
//! - **face axis** - the incident face of the other box is clipped against the
//!   side planes of the reference face (Sutherland-Hodgman) and the clipped
//!   polygon is reduced to at most four points;
//! - **edge axis** - the two support edges are located and their closest points
//!   form a single contact.
//!
//! The contact normal points from box A toward box B. `pos_a` lies on box A and
//! `pos_b` on box B, and `pos_a - pos_b` is parallel to the normal. A contact
//! that penetrates has `pos_a - pos_b = depth * normal`. A contact that only
//! lies inside the margin (a speculative contact, for example at the raised
//! side of a slightly tilted face, or between boxes that do not overlap at all)
//! has `depth == 0` and `pos_a - pos_b = -gap * normal` with
//! `0 <= gap <= tolerance`.

use super::{
    ContactPoint, MAX_MANIFOLD_POINTS, ManifoldReduction, PolyhedronContact, add3, cross3, dot3,
    find_incident_face, neg3, norm3, normalize3, scale3, sub3,
};
use crate::narrowphase::dispatch::{BoxSatAxis, box_box_sat_with_margin};

/// A convex polygon clipped against four side planes has at most eight
/// vertices; the clipping stage keeps all of them before the reduction to
/// [`MAX_MANIFOLD_POINTS`].
const MAX_CLIPPED_POINTS: usize = 2 * MAX_MANIFOLD_POINTS;

/// Smallest axis length accepted as a valid (non-degenerate) box axis.
const MIN_AXIS_LENGTH: f64 = 1e-12;

/// Specialized box-box contact manifold generator.
///
/// Handles all 6 face pairs and the 9 distinct edge-edge axes of two oriented
/// boxes via the separating-axis test, then generates a manifold of up to
/// [`MAX_MANIFOLD_POINTS`] contacts. Face contacts come from clipping the
/// incident face against the reference face; edge contacts are a single point.
#[derive(Debug, Clone)]
pub struct BoxBoxManifold {
    /// Contact margin: boxes separated by at most this distance (on every
    /// separating-axis candidate) still produce contacts, with zero depth.
    pub tolerance: f64,
}

/// An oriented box in world space.
#[derive(Debug, Clone, Copy)]
struct OrientedBox {
    center: [f64; 3],
    half: [f64; 3],
    /// Orthonormal world-space axes.
    axes: [[f64; 3]; 3],
}

/// One face of an [`OrientedBox`].
#[derive(Debug, Clone, Copy)]
struct BoxFace {
    /// Outward unit normal.
    normal: [f64; 3],
    /// Corners, counter-clockwise when seen from outside (looking along
    /// `-normal`).
    corners: [[f64; 3]; 4],
}

impl OrientedBox {
    /// Build a box, normalising the axes. `None` for non-finite or degenerate
    /// input (zero-length axis, negative half-extent).
    fn new(center: [f64; 3], half: [f64; 3], axes: [[f64; 3]; 3]) -> Option<Self> {
        let finite = center.iter().chain(half.iter()).all(|v| v.is_finite())
            && axes.iter().flatten().all(|v| v.is_finite());
        if !finite || half.iter().any(|&h| h < 0.0) {
            return None;
        }
        let mut unit_axes = [[0.0_f64; 3]; 3];
        for (unit, axis) in unit_axes.iter_mut().zip(axes.iter()) {
            let length = norm3(*axis);
            if length < MIN_AXIS_LENGTH {
                return None;
            }
            *unit = scale3(*axis, 1.0 / length);
        }
        Some(Self {
            center,
            half,
            axes: unit_axes,
        })
    }

    /// Outward normals of the six faces: `+axis0, -axis0, +axis1, -axis1,
    /// +axis2, -axis2`.
    fn face_normals(&self) -> [[f64; 3]; 6] {
        let [a0, a1, a2] = self.axes;
        [a0, neg3(a0), a1, neg3(a1), a2, neg3(a2)]
    }

    /// Face `index` (numbered as in [`Self::face_normals`]).
    fn face(&self, index: usize) -> BoxFace {
        let k = (index / 2) % 3;
        let sign = if index.is_multiple_of(2) { 1.0 } else { -1.0 };
        let normal = scale3(self.axes[k], sign);
        let center = add3(self.center, scale3(normal, self.half[k]));
        let u = scale3(self.axes[(k + 1) % 3], self.half[(k + 1) % 3]);
        let v = scale3(self.axes[(k + 2) % 3], self.half[(k + 2) % 3]);
        let mut corners = [
            sub3(sub3(center, u), v),
            sub3(add3(center, u), v),
            add3(add3(center, u), v),
            add3(sub3(center, u), v),
        ];
        // The axes need not form a right-handed frame: enforce the winding.
        let winding = cross3(sub3(corners[1], corners[0]), sub3(corners[2], corners[1]));
        if dot3(winding, normal) < 0.0 {
            corners.swap(1, 3);
        }
        BoxFace { normal, corners }
    }

    /// The face whose outward normal is most aligned with `dir`.
    fn face_along(&self, dir: [f64; 3]) -> BoxFace {
        let (index, _) = find_incident_face(&self.face_normals(), neg3(dir));
        self.face(index)
    }

    /// End points of the edge parallel to local axis `axis` that is extreme in
    /// direction `dir`.
    fn support_edge(&self, axis: usize, dir: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let mut middle = self.center;
        for (k, (ax, &half)) in self.axes.iter().zip(self.half.iter()).enumerate() {
            if k == axis {
                continue;
            }
            let sign = if dot3(*ax, dir) >= 0.0 { 1.0 } else { -1.0 };
            middle = add3(middle, scale3(*ax, sign * half));
        }
        let along = scale3(self.axes[axis % 3], self.half[axis % 3]);
        (sub3(middle, along), add3(middle, along))
    }
}

/// How far `point`, lying in the plane of `face`, misses the face: zero when
/// it is inside the rectangle, otherwise the largest distance by which it lies
/// outside one of the four side planes.
fn outside_distance(face: &BoxFace, point: [f64; 3]) -> f64 {
    (0..face.corners.len())
        .map(|i| {
            let start = face.corners[i];
            let edge = sub3(face.corners[(i + 1) % face.corners.len()], start);
            let inward = normalize3(cross3(face.normal, edge));
            -dot3(sub3(point, start), inward)
        })
        .fold(0.0, f64::max)
}

/// The corner of `incident` whose projection onto the plane of `reference`
/// lies closest to the reference face (the deeper corner wins a tie).
fn nearest_corner(reference: &BoxFace, incident: &BoxFace) -> Option<[f64; 3]> {
    let key = |corner: [f64; 3]| {
        let depth = dot3(sub3(reference.corners[0], corner), reference.normal);
        let on_plane = add3(corner, scale3(reference.normal, depth));
        (outside_distance(reference, on_plane), depth)
    };
    incident.corners.iter().copied().min_by(|&a, &b| {
        let (miss_a, depth_a) = key(a);
        let (miss_b, depth_b) = key(b);
        miss_a.total_cmp(&miss_b).then(depth_b.total_cmp(&depth_a))
    })
}

impl BoxBoxManifold {
    /// Create a new box-box manifold generator.
    pub fn new(tolerance: f64) -> Self {
        Self { tolerance }
    }

    /// Compute the box-box contact manifold.
    ///
    /// Each box is given by its centre, half-extents and three orthogonal
    /// world-space axes (normalised internally). Returns up to
    /// [`MAX_MANIFOLD_POINTS`] contacts, or an empty list when the boxes are
    /// separated by more than [`BoxBoxManifold::tolerance`] along some
    /// separating axis or when the input is degenerate (non-finite values,
    /// negative half-extents, zero-length axes).
    ///
    /// The contact normal points from box A toward box B; `pos_a` is on A and
    /// `pos_b` on B, and `depth` is the (non-negative) penetration depth of
    /// that contact. A contact that only lies inside the margin has
    /// `depth == 0.0`; see the module documentation for the witness offsets.
    pub fn compute(
        &self,
        center_a: [f64; 3],
        half_a: [f64; 3],
        axes_a: [[f64; 3]; 3],
        center_b: [f64; 3],
        half_b: [f64; 3],
        axes_b: [[f64; 3]; 3],
    ) -> Vec<ContactPoint> {
        let (Some(box_a), Some(box_b)) = (
            OrientedBox::new(center_a, half_a, axes_a),
            OrientedBox::new(center_b, half_b, axes_b),
        ) else {
            return Vec::new();
        };
        let margin = if self.tolerance.is_finite() {
            self.tolerance.max(0.0)
        } else {
            0.0
        };
        let Some(sat) = box_box_sat_with_margin(
            box_a.center,
            box_a.axes,
            box_a.half,
            box_b.center,
            box_b.axes,
            box_b.half,
            margin,
        ) else {
            return Vec::new();
        };
        // The SAT helper reports the minimum translation direction from B
        // toward A; this generator reports normals from A toward B.
        let normal_ab = neg3(sat.normal);
        match sat.axis {
            BoxSatAxis::FaceA => Self::face_contacts(&box_a, &box_b, normal_ab, true, margin),
            BoxSatAxis::FaceB => Self::face_contacts(&box_b, &box_a, sat.normal, false, margin),
            BoxSatAxis::EdgeEdge(axis_a, axis_b) => {
                Self::edge_contact(&box_a, &box_b, axis_a, axis_b, normal_ab, sat.depth)
            }
        }
    }

    /// Face contact: clip the incident face of `incident` against the side
    /// planes of the reference face of `reference`.
    ///
    /// `outward` is the direction in which the reference face looks (away from
    /// the reference box); `reference_is_a` tells which of the two boxes is
    /// box A so the result can be expressed in the A-to-B convention.
    fn face_contacts(
        reference: &OrientedBox,
        incident: &OrientedBox,
        outward: [f64; 3],
        reference_is_a: bool,
        margin: f64,
    ) -> Vec<ContactPoint> {
        let ref_face = reference.face_along(outward);
        let (incident_index, _) = find_incident_face(&incident.face_normals(), ref_face.normal);
        let inc_face = incident.face(incident_index);

        let clipper = PolyhedronContact::new(MAX_CLIPPED_POINTS);
        let mut clipped = clipper.face_face_contact(
            &ref_face.corners,
            &inc_face.corners,
            ref_face.normal,
            margin,
        );
        if clipped.is_empty()
            && let Some(corner) = nearest_corner(&ref_face, &inc_face)
        {
            // Nothing of the incident face lies over the reference face: the
            // boxes only touch within the margin at an edge or corner of the
            // incident face. Use the incident corner closest to the reference
            // face, which is where the speculative contact belongs.
            let penetration = dot3(sub3(ref_face.corners[0], corner), ref_face.normal);
            if penetration > -margin {
                let on_reference = add3(corner, scale3(ref_face.normal, penetration));
                clipped.push(ContactPoint::new(
                    on_reference,
                    corner,
                    ref_face.normal,
                    penetration,
                ));
            }
        }

        let contacts: Vec<ContactPoint> = clipped
            .into_iter()
            .map(|c| {
                let depth = c.depth.max(0.0);
                if reference_is_a {
                    ContactPoint::new(c.pos_a, c.pos_b, c.normal, depth)
                } else {
                    // Reference is box B: its face point is `pos_b`, the
                    // incident point (on A) is `pos_a`, and the A-to-B normal
                    // is opposite to the reference face normal.
                    ContactPoint::new(c.pos_b, c.pos_a, neg3(c.normal), depth)
                }
            })
            .collect();
        ManifoldReduction::reduce(&contacts)
    }

    /// Edge-edge contact: the closest points of the two support edges.
    fn edge_contact(
        box_a: &OrientedBox,
        box_b: &OrientedBox,
        axis_a: usize,
        axis_b: usize,
        normal_ab: [f64; 3],
        depth: f64,
    ) -> Vec<ContactPoint> {
        let (a0, a1) = box_a.support_edge(axis_a, normal_ab);
        let (b0, b1) = box_b.support_edge(axis_b, neg3(normal_ab));
        let (on_a, on_b, _) =
            PolyhedronContact::new(1).edge_edge_closest(a0, sub3(a1, a0), b0, sub3(b1, b0));
        vec![ContactPoint::new(on_a, on_b, normal_ab, depth.max(0.0))]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::narrowphase::dispatch::box_box_sat_detailed;

    const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    const UNIT_HALF: [f64; 3] = [1.0; 3];

    /// Columns of the rotation by `angle` about `axis` (Rodrigues' formula).
    fn rotation_axes(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
        let k = normalize3(axis);
        let (sin, cos) = angle.sin_cos();
        let rotate = |v: [f64; 3]| {
            let kxv = cross3(k, v);
            let kdv = dot3(k, v);
            [
                v[0] * cos + kxv[0] * sin + k[0] * kdv * (1.0 - cos),
                v[1] * cos + kxv[1] * sin + k[1] * kdv * (1.0 - cos),
                v[2] * cos + kxv[2] * sin + k[2] * kdv * (1.0 - cos),
            ]
        };
        [
            rotate(IDENTITY[0]),
            rotate(IDENTITY[1]),
            rotate(IDENTITY[2]),
        ]
    }

    fn assert_close(actual: [f64; 3], expected: [f64; 3], tol: f64) {
        let diff = norm3(sub3(actual, expected));
        assert!(diff < tol, "{actual:?} differs from {expected:?} by {diff}");
    }

    /// The witness relation of a contact: `pos_a - pos_b` is parallel to the
    /// normal and equals `depth` along it for a penetrating contact, or minus
    /// the gap (at most `margin`) for a contact that only lies inside the
    /// margin.
    fn assert_witness_offset(contact: &ContactPoint, margin: f64) {
        let offset = sub3(contact.pos_a, contact.pos_b);
        let along = dot3(offset, contact.normal);
        assert_close(offset, scale3(contact.normal, along), 1e-7);
        if contact.depth > 0.0 {
            assert!(
                (along - contact.depth).abs() < 1e-7,
                "{contact:?}: offset {along}"
            );
        } else {
            assert!(
                along <= 1e-7 && along >= -margin - 1e-7,
                "{contact:?}: offset {along}"
            );
        }
    }

    fn generator() -> BoxBoxManifold {
        BoxBoxManifold::new(0.01)
    }

    #[test]
    fn separated_boxes_have_no_contacts() {
        let gen_ = generator();
        // Unit-half boxes: B at (-3, 0, 0) leaves a gap of 1.
        let contacts = gen_.compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [-3.0, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        assert!(contacts.is_empty(), "got {contacts:?}");
        for centre_b in [[3.0, 0.0, 0.0], [0.0, -2.5, 0.0], [0.0, 0.0, 2.2]] {
            let contacts =
                gen_.compute([0.0; 3], UNIT_HALF, IDENTITY, centre_b, UNIT_HALF, IDENTITY);
            assert!(contacts.is_empty(), "centre {centre_b:?}: {contacts:?}");
        }
    }

    #[test]
    fn face_contact_reports_axis_normal_and_true_depth() {
        let contacts = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [1.5, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        let mut corners_seen = Vec::new();
        for c in &contacts {
            assert_close(c.normal, [1.0, 0.0, 0.0], 1e-9);
            assert!((c.depth - 0.5).abs() < 1e-9, "depth {}", c.depth);
            assert!((c.pos_a[0] - 1.0).abs() < 1e-9, "pos_a {:?}", c.pos_a);
            assert!((c.pos_b[0] - 0.5).abs() < 1e-9, "pos_b {:?}", c.pos_b);
            assert_close(sub3(c.pos_a, c.pos_b), scale3(c.normal, c.depth), 1e-9);
            assert!((c.pos_a[1].abs() - 1.0).abs() < 1e-9);
            assert!((c.pos_a[2].abs() - 1.0).abs() < 1e-9);
            corners_seen.push((c.pos_a[1] > 0.0, c.pos_a[2] > 0.0));
        }
        corners_seen.sort_unstable();
        corners_seen.dedup();
        assert_eq!(corners_seen.len(), 4, "all four corners must be covered");
    }

    #[test]
    fn penetration_axis_is_the_least_penetrating_one() {
        // Deep overlap in x (1.9), shallow in y (0.2): the contact must be
        // along y with depth 0.2, not along the deepest axis.
        let contacts = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [0.1, 1.8, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        assert!(!contacts.is_empty());
        for c in &contacts {
            assert_close(c.normal, [0.0, 1.0, 0.0], 1e-9);
            assert!((c.depth - 0.2).abs() < 1e-9, "depth {}", c.depth);
        }
    }

    #[test]
    fn swapping_the_boxes_mirrors_the_manifold() {
        let gen_ = generator();
        let ab = gen_.compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [1.5, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        let ba = gen_.compute(
            [1.5, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
        );
        assert_eq!(ab.len(), ba.len());
        for c in &ab {
            let mirrored = ba
                .iter()
                .find(|d| {
                    norm3(sub3(d.pos_a, c.pos_b)) < 1e-9 && norm3(sub3(d.pos_b, c.pos_a)) < 1e-9
                })
                .expect("mirrored contact");
            assert_close(mirrored.normal, scale3(c.normal, -1.0), 1e-9);
            assert!((mirrored.depth - c.depth).abs() < 1e-9);
        }
    }

    #[test]
    fn small_box_on_large_box_contacts_its_own_footprint() {
        let large = [2.0, 1.0, 2.0];
        let small = [0.5; 3];
        let contacts = generator().compute(
            [0.0; 3],
            large,
            IDENTITY,
            [0.3, 1.45, -0.2],
            small,
            IDENTITY,
        );
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        for c in &contacts {
            assert_close(c.normal, [0.0, 1.0, 0.0], 1e-9);
            assert!((c.depth - 0.05).abs() < 1e-9, "depth {}", c.depth);
            assert!((c.pos_a[1] - 1.0).abs() < 1e-9);
            assert!((c.pos_b[1] - 0.95).abs() < 1e-9);
            // The footprint of the small box: x in {-0.2, 0.8}, z in {-0.7, 0.3}.
            assert!(
                ((c.pos_b[0] + 0.2).abs() < 1e-9) || ((c.pos_b[0] - 0.8).abs() < 1e-9),
                "x {}",
                c.pos_b[0]
            );
            assert!(
                ((c.pos_b[2] + 0.7).abs() < 1e-9) || ((c.pos_b[2] - 0.3).abs() < 1e-9),
                "z {}",
                c.pos_b[2]
            );
        }
    }

    #[test]
    fn large_incident_face_is_clipped_to_the_reference_face() {
        // Same scene with the roles reversed: the small box is the reference,
        // the large box's top face is clipped down to the small footprint.
        let contacts = generator().compute(
            [0.3, 1.45, -0.2],
            [0.5; 3],
            IDENTITY,
            [0.0; 3],
            [2.0, 1.0, 2.0],
            IDENTITY,
        );
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        for c in &contacts {
            assert_close(c.normal, [0.0, -1.0, 0.0], 1e-9);
            assert!((c.depth - 0.05).abs() < 1e-9);
            assert!((c.pos_a[1] - 0.95).abs() < 1e-9);
            assert!((c.pos_b[1] - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn tilted_box_on_ground_touches_with_its_lowest_edge() {
        // Ground = box B (top face at y = 0.5); box A is tilted 10 degrees about
        // x and sunk by 0.05 so only one bottom edge penetrates. The face axis
        // of B wins, so box B is the reference box.
        let tilt = 10.0_f64.to_radians();
        let axes_a = rotation_axes([1.0, 0.0, 0.0], tilt);
        let half_a = [0.5; 3];
        let centre_a = [0.0, 0.45 + 0.5 * (tilt.cos() + tilt.sin()), 0.0];
        let half_b = [5.0, 0.5, 5.0];
        let sat = box_box_sat_detailed([0.0; 3], IDENTITY, half_b, centre_a, axes_a, half_a)
            .expect("overlap");
        assert!(
            matches!(sat.axis, BoxSatAxis::FaceA),
            "ground is first here"
        );
        let sat = box_box_sat_detailed(centre_a, axes_a, half_a, [0.0; 3], IDENTITY, half_b)
            .expect("overlap");
        assert!(matches!(sat.axis, BoxSatAxis::FaceB), "reference must be B");

        let contacts = generator().compute(centre_a, half_a, axes_a, [0.0; 3], half_b, IDENTITY);
        assert_eq!(contacts.len(), 2, "got {contacts:?}");
        for c in &contacts {
            assert_close(c.normal, [0.0, -1.0, 0.0], 1e-9);
            assert!((c.depth - 0.05).abs() < 1e-9, "depth {}", c.depth);
            assert!((c.pos_a[1] - 0.45).abs() < 1e-9, "pos_a {:?}", c.pos_a);
            assert!((c.pos_b[1] - 0.5).abs() < 1e-9, "pos_b {:?}", c.pos_b);
            assert!((c.pos_a[0].abs() - 0.5).abs() < 1e-9);
            // The lowest edge of the tilted box: its two lowest corners have
            // z = 0.5 * (cos - sin) after the tilt about x.
            assert!((c.pos_a[2] - 0.5 * (tilt.cos() - tilt.sin())).abs() < 1e-9);
        }
    }

    #[test]
    fn slightly_tilted_box_gets_speculative_contacts_on_its_raised_edge() {
        // Ground = box B (top face at y = 0.5). Box A is tilted by 0.01 rad
        // about z: its low bottom edge sinks 5 mm into the ground while the
        // opposite bottom edge hovers about 5 mm above it, inside the 10 mm
        // margin. All four bottom corners are reported; only the sunk ones
        // penetrate, the raised ones carry zero depth and keep their gap in
        // the witness offset.
        let tilt = 0.01_f64;
        let axes_a = rotation_axes([0.0, 0.0, 1.0], tilt);
        let half_a = [0.5; 3];
        let half_b = [5.0, 0.5, 5.0];
        // Heights of the lowest and of the raised bottom corners relative to
        // the centre of A.
        let low = -0.5 * (tilt.sin() + tilt.cos());
        let high = 0.5 * (tilt.sin() - tilt.cos());
        let centre_a = [0.0, 0.495 - low, 0.0];
        let raised_gap = tilt.sin() - 0.005;
        assert!(raised_gap > 0.0 && raised_gap < 0.01);
        assert!((centre_a[1] + high - 0.5 - raised_gap).abs() < 1e-12);

        let contacts = generator().compute(centre_a, half_a, axes_a, [0.0; 3], half_b, IDENTITY);
        assert_eq!(contacts.len(), 4, "got {contacts:?}");
        let (mut sunk, mut raised) = (0, 0);
        for c in &contacts {
            assert_close(c.normal, [0.0, -1.0, 0.0], 1e-9);
            assert_witness_offset(c, 0.01);
            assert!((c.pos_b[1] - 0.5).abs() < 1e-9, "pos_b {:?}", c.pos_b);
            if c.depth > 0.0 {
                sunk += 1;
                assert!((c.depth - 0.005).abs() < 1e-9, "depth {}", c.depth);
                assert!((c.pos_a[1] - 0.495).abs() < 1e-9, "pos_a {:?}", c.pos_a);
            } else {
                raised += 1;
                assert_eq!(c.depth, 0.0);
                assert!(
                    (c.pos_a[1] - (0.5 + raised_gap)).abs() < 1e-9,
                    "pos_a {:?}",
                    c.pos_a
                );
            }
        }
        assert_eq!((sunk, raised), (2, 2));
    }

    #[test]
    fn edge_edge_contact_is_a_single_point_between_the_edges() {
        let half = [0.5; 3];
        let centre_b = [0.75, 0.75, 0.75];
        let axes_b = rotation_axes([1.0, 1.0, 1.0], std::f64::consts::FRAC_PI_4);
        let sat = box_box_sat_detailed([0.0; 3], IDENTITY, half, centre_b, axes_b, half)
            .expect("overlap");
        assert!(
            matches!(sat.axis, BoxSatAxis::EdgeEdge(_, _)),
            "scene must exercise the edge-edge axis"
        );

        let contacts = generator().compute([0.0; 3], half, IDENTITY, centre_b, half, axes_b);
        assert_eq!(contacts.len(), 1, "got {contacts:?}");
        let c = &contacts[0];
        assert!((norm3(c.normal) - 1.0).abs() < 1e-9);
        assert!(
            dot3(c.normal, sub3(centre_b, [0.0; 3])) > 0.0,
            "normal must point from A toward B"
        );
        assert!(c.depth > 0.0);
        assert!((c.depth - sat.depth).abs() < 1e-12);
        assert_close(sub3(c.pos_a, c.pos_b), scale3(c.normal, c.depth), 1e-6);
        // pos_a lies on an edge of A (two coordinates at the face planes).
        let on_faces_a = c
            .pos_a
            .iter()
            .filter(|v| (v.abs() - 0.5).abs() < 1e-9)
            .count();
        assert!(
            on_faces_a >= 2,
            "pos_a {:?} is not on an edge of A",
            c.pos_a
        );
        assert!(c.pos_a.iter().all(|v| v.abs() <= 0.5 + 1e-9));
        // pos_b lies on an edge of B (checked in B's frame).
        let local_b: Vec<f64> = axes_b
            .iter()
            .map(|axis| dot3(*axis, sub3(c.pos_b, centre_b)))
            .collect();
        let on_faces_b = local_b
            .iter()
            .filter(|v| (v.abs() - 0.5).abs() < 1e-9)
            .count();
        assert!(
            on_faces_b >= 2,
            "pos_b local {local_b:?} is not on an edge of B"
        );
        assert!(local_b.iter().all(|v| v.abs() <= 0.5 + 1e-9));
    }

    #[test]
    fn tolerance_acts_as_a_contact_margin() {
        let separated_by_5mm = [2.005, 0.0, 0.0];
        let near = BoxBoxManifold::new(0.01).compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            separated_by_5mm,
            UNIT_HALF,
            IDENTITY,
        );
        assert_eq!(near.len(), 4, "got {near:?}");
        for c in &near {
            assert_close(c.normal, [1.0, 0.0, 0.0], 1e-9);
            assert_eq!(c.depth, 0.0, "a gap is reported as zero depth");
        }
        let far = BoxBoxManifold::new(0.001).compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            separated_by_5mm,
            UNIT_HALF,
            IDENTITY,
        );
        assert!(far.is_empty(), "got {far:?}");
    }

    #[test]
    fn rotated_overlap_yields_valid_manifold() {
        let axes_b = rotation_axes([0.0, 1.0, 0.0], 0.6);
        let contacts = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [0.4, 1.4, 0.3],
            [0.8, 0.5, 0.6],
            axes_b,
        );
        assert!(!contacts.is_empty() && contacts.len() <= MAX_MANIFOLD_POINTS);
        for c in &contacts {
            assert!((norm3(c.normal) - 1.0).abs() < 1e-9);
            assert!(c.depth >= 0.0);
            assert!(
                c.normal[1] > 0.99,
                "stacked boxes touch along y: {:?}",
                c.normal
            );
            // A's top face is at y = 1, B's bottom face at y = 1.4 - 0.5.
            assert!((c.depth - 0.1).abs() < 1e-9, "depth {}", c.depth);
            assert_close(sub3(c.pos_a, c.pos_b), scale3(c.normal, c.depth), 1e-9);
        }
    }

    #[test]
    fn rotated_boxes_beyond_the_margin_do_not_touch() {
        let axes_b = rotation_axes([0.0, 1.0, 0.0], 0.6);
        // B's bottom face is at y = 1.2: a gap of 0.2 above A's top face.
        let contacts = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [0.4, 1.7, 0.3],
            [0.8, 0.5, 0.6],
            axes_b,
        );
        assert!(contacts.is_empty(), "got {contacts:?}");
    }

    #[test]
    fn degenerate_input_yields_no_contacts() {
        let gen_ = generator();
        let nan_centre = [f64::NAN, 0.0, 0.0];
        assert!(
            gen_.compute(
                nan_centre, UNIT_HALF, IDENTITY, [0.0; 3], UNIT_HALF, IDENTITY
            )
            .is_empty()
        );
        let flat_axes = [[0.0; 3], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(
            gen_.compute(
                [0.0; 3],
                UNIT_HALF,
                flat_axes,
                [0.5, 0.0, 0.0],
                UNIT_HALF,
                IDENTITY
            )
            .is_empty()
        );
        assert!(
            gen_.compute(
                [0.0; 3],
                [-1.0, 1.0, 1.0],
                IDENTITY,
                [0.5, 0.0, 0.0],
                UNIT_HALF,
                IDENTITY
            )
            .is_empty()
        );
    }

    #[test]
    fn boxes_touching_corner_to_corner_within_the_margin_get_one_contact() {
        // Two unit cubes offset diagonally by 4 mm in x and z: they lie inside
        // each other's margin, but their side faces do not overlap, so the clip
        // of the incident face is empty and the corner of the incident face
        // closest to the reference face is used.
        let contacts = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [2.004, 0.5, 2.004],
            UNIT_HALF,
            IDENTITY,
        );
        assert_eq!(contacts.len(), 1, "got {contacts:?}");
        let c = &contacts[0];
        assert_close(c.normal, [1.0, 0.0, 0.0], 1e-12);
        assert_eq!(c.depth, 0.0, "a gap is reported as zero depth");
        // The B corner that is nearest to A's +x face (y inside A's range).
        assert_close(c.pos_b, [1.004, -0.5, 1.004], 1e-12);
        assert_close(c.pos_a, [1.0, -0.5, 1.004], 1e-12);
    }

    /// Deterministic linear congruential generator for reproducible scenes.
    struct Lcg(u64);

    impl Lcg {
        fn new(seed: u64) -> Self {
            Self(seed)
        }

        /// Uniform sample in `[0, 1)`.
        fn unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) as f64 / (1u64 << 31) as f64
        }

        fn range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.unit()
        }

        fn direction(&mut self) -> [f64; 3] {
            normalize3([
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            ])
        }
    }

    fn box_corners(centre: [f64; 3], half: [f64; 3], axes: [[f64; 3]; 3]) -> Vec<[f64; 3]> {
        let mut corners = Vec::with_capacity(8);
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                for sz in [-1.0, 1.0] {
                    let offset = add3(
                        add3(scale3(axes[0], sx * half[0]), scale3(axes[1], sy * half[1])),
                        scale3(axes[2], sz * half[2]),
                    );
                    corners.push(add3(centre, offset));
                }
            }
        }
        corners
    }

    /// Signed gap between the projections of two vertex sets onto `axis`
    /// (negative when they overlap: minus the depth needed to separate them).
    fn projected_gap(a: &[[f64; 3]], b: &[[f64; 3]], axis: [f64; 3]) -> f64 {
        let extent = |points: &[[f64; 3]]| {
            points
                .iter()
                .map(|p| dot3(*p, axis))
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                    (lo.min(v), hi.max(v))
                })
        };
        let ((a_lo, a_hi), (b_lo, b_hi)) = (extent(a), extent(b));
        (b_lo - a_hi).max(a_lo - b_hi)
    }

    fn inside_box(point: [f64; 3], centre: [f64; 3], half: [f64; 3], axes: [[f64; 3]; 3]) -> bool {
        let offset = sub3(point, centre);
        axes.iter()
            .zip(half.iter())
            .all(|(axis, &h)| dot3(*axis, offset).abs() <= h + 1e-7)
    }

    #[test]
    fn random_box_pairs_satisfy_the_manifold_invariants() {
        let margin = 0.01;
        let generator = BoxBoxManifold::new(margin);
        let mut rng = Lcg::new(99);
        let (mut overlapping, mut apart) = (0, 0);
        for _ in 0..1500 {
            let half_a = [
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
            ];
            let half_b = [
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
                rng.range(0.3, 1.5),
            ];
            let axes_a = rotation_axes(rng.direction(), rng.range(0.0, std::f64::consts::PI));
            let axes_b = rotation_axes(rng.direction(), rng.range(0.0, std::f64::consts::PI));
            let centre_a = [0.0; 3];
            let reach = norm3(half_a) + norm3(half_b);
            let centre_b = scale3(rng.direction(), rng.range(0.1, 1.15 * reach));

            // Independent separating-axis test on explicit vertex projections.
            let corners_a = box_corners(centre_a, half_a, axes_a);
            let corners_b = box_corners(centre_b, half_b, axes_b);
            let mut axes: Vec<[f64; 3]> = axes_a.iter().chain(axes_b.iter()).copied().collect();
            for axis_a in &axes_a {
                for axis_b in &axes_b {
                    let cross = cross3(*axis_a, *axis_b);
                    if norm3(cross) > 1e-6 {
                        axes.push(normalize3(cross));
                    }
                }
            }
            let max_gap = axes
                .iter()
                .map(|axis| projected_gap(&corners_a, &corners_b, *axis))
                .fold(f64::NEG_INFINITY, f64::max);

            let contacts = generator.compute(centre_a, half_a, axes_a, centre_b, half_b, axes_b);
            if max_gap > margin + 1e-6 {
                apart += 1;
                assert!(contacts.is_empty(), "gap {max_gap}: {contacts:?}");
                continue;
            }
            if max_gap > -1e-6 {
                // Touching or inside the margin: only checked by dedicated tests.
                continue;
            }
            overlapping += 1;
            assert!(
                !contacts.is_empty() && contacts.len() <= MAX_MANIFOLD_POINTS,
                "overlap {}: {contacts:?}",
                -max_gap
            );
            let edge_bias = 1e-4
                * half_a
                    .iter()
                    .chain(half_b.iter())
                    .fold(1.0_f64, |m, h| m.max(*h));
            for c in &contacts {
                assert!((norm3(c.normal) - 1.0).abs() < 1e-9);
                assert!(dot3(c.normal, sub3(centre_b, centre_a)) > -1e-9);
                // No contact is deeper than the axis the SAT picked, which is
                // the axis of least penetration up to the edge bias.
                assert!(
                    c.depth >= 0.0 && c.depth <= -max_gap + edge_bias + 1e-9,
                    "depth {} vs least penetration {}",
                    c.depth,
                    -max_gap
                );
                assert!(inside_box(c.pos_a, centre_a, half_a, axes_a), "{c:?}");
                assert!(inside_box(c.pos_b, centre_b, half_b, axes_b), "{c:?}");
                assert_witness_offset(c, margin);
            }
        }
        assert!(overlapping > 500, "only {overlapping} overlapping pairs");
        assert!(apart > 300, "only {apart} separated pairs");
    }

    #[test]
    fn axes_are_normalised_before_use() {
        let stretched = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.5]];
        let a = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            stretched,
            [1.5, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        let b = generator().compute(
            [0.0; 3],
            UNIT_HALF,
            IDENTITY,
            [1.5, 0.0, 0.0],
            UNIT_HALF,
            IDENTITY,
        );
        assert_eq!(a.len(), b.len());
        assert!(a.iter().all(|c| (c.depth - 0.5).abs() < 1e-9));
    }
}
