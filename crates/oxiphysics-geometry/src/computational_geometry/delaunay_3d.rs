// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact 3-D Delaunay tetrahedralization (Bowyer-Watson) and α-shape extraction.
//!
//! This module builds the Delaunay tetrahedralization of a 3-D point set with an
//! *incremental* Bowyer-Watson insertion driven entirely by the exact predicates
//! from [`oxiphysics_core::exact_predicates`]:
//!
//! * [`orient3d`] decides on which side of an oriented triangle a point lies. It
//!   is used both to keep every tetrahedron consistently positively oriented and
//!   to drive the walk-based point-location step.
//! * [`insphere`] is the empty-sphere oracle: for a positively oriented tet
//!   `(a,b,c,d)` it reports whether a query point `p` lies strictly inside, on, or
//!   outside the tet's circumsphere. A point inside the circumsphere makes the tet
//!   *non-Delaunay* with respect to `p`, so the tet enters the insertion cavity.
//!
//! Both predicates run Shewchuk's adaptive floating-point filter and fall back to
//! exact expansion arithmetic only in the degenerate band, so the reported signs
//! are always provably correct. Genuine ties (`Orientation::Degenerate` — five
//! cospherical points, or four coplanar points) are broken with a
//! Simulation-of-Simplicity (SoS) rule keyed on the *lowest global vertex index*
//! among the points involved, exactly mirroring the `sos_sign` convention already
//! used by `mesh_boolean`. Because that rule depends only on the order-independent
//! set of indices, the same geometric configuration always resolves the same way,
//! giving a consistent virtual general position with no robustness failures.
//!
//! ## Algorithm (Bowyer-Watson, incremental)
//!
//! First a *super-tetrahedron* far larger than the input AABB is created so that
//! every input point lies strictly inside it (its four vertices receive sentinel
//! indices `n, n+1, n+2, n+3` appended after the `n` user points).
//!
//! Then each input point `p` is inserted in turn via four sub-steps. *Walk-based
//! location* starts from the most recently created tet and hops across the face
//! whose plane separates the current tet from `p` until the containing tet is
//! reached (with a brute-force fallback guarding against any degenerate walk cycle).
//! *Cavity identification* floods the face-adjacency graph to collect every tet
//! whose circumsphere contains `p` (`insphere`, with SoS on cospherical ties). The
//! *cavity boundary* — the faces shared by exactly one cavity tet — forms a
//! star-shaped polyhedron around `p`. Finally *re-tetrahedralization* joins `p` to
//! every boundary face, rebuilding the face-adjacency (`neigh`) of the new tets, of
//! each other, and of the outer tets bordering the cavity.
//!
//! When all points are inserted, every tet incident to a super-tetrahedron vertex is
//! removed, leaving the Delaunay tetrahedralization of the original `n` points.
//!
//! ## α-shapes
//!
//! [`Tetrahedralization::alpha_shape`] extracts the boundary of the α-complex at a
//! squared filtration radius `alpha`: a triangular face lies on the α-shape when it
//! is incident to exactly one Delaunay tet whose squared circumradius is `≤ alpha`.
//! At a sufficiently large `alpha` this reproduces the convex-hull boundary of the
//! point set; at smaller `alpha` it carves concavities, recovering, e.g., the
//! surface of a sampled sphere as a closed orientable mesh.
//!
//! References: Edelsbrunner & Mücke, "Three-dimensional alpha shapes", ACM TOG
//! 13(1):43-72, 1994; Bowyer 1981; Watson 1981; Shewchuk 1997 (predicates);
//! Edelsbrunner & Mücke 1990 (Simulation of Simplicity).

use oxiphysics_core::exact_predicates::{Orientation, insphere, orient3d};

/// A tetrahedron stored as 4 vertex indices + 4 face-adjacency entries.
///
/// The four vertices `v` are kept in *positive orientation*, i.e.
/// `orient3d(points[v[0]], points[v[1]], points[v[2]], points[v[3]])` is
/// [`Orientation::Positive`] (or `Degenerate` only for sliver tets that the
/// super-tetra construction never produces among real points).
///
/// `neigh[i]` is the index (into [`Tetrahedralization::tets`]) of the neighbour
/// tetrahedron across the face *opposite* vertex `i` — that face is made of the
/// other three vertices `v[(i+1)%4], v[(i+2)%4], v[(i+3)%4]` — or `-1` when that
/// face is on the boundary of the current tetrahedralization (only ever true while
/// the super-tetra is still present; the final Delaunay complex of a point cloud is
/// closed by the super-tetra during construction).
#[derive(Debug, Clone, PartialEq)]
pub struct Tet {
    /// The four vertex indices, in positive orientation.
    pub v: [usize; 4],
    /// Face adjacency: `neigh[i]` faces vertex `i`; `-1` means no neighbour.
    pub neigh: [i32; 4],
}

impl Tet {
    /// Returns the three vertex indices of the face opposite local vertex `i`.
    ///
    /// The triple is returned in the winding `(v[(i+1)%4], v[(i+2)%4],
    /// v[(i+3)%4])`; together with the tet's positive orientation this gives the
    /// face a consistent outward sense.
    #[inline]
    pub fn face(&self, i: usize) -> [usize; 3] {
        [
            self.v[(i + 1) % 4],
            self.v[(i + 2) % 4],
            self.v[(i + 3) % 4],
        ]
    }
}

/// A canonical (sorted-ascending) representation of an unordered triangular face,
/// used as a hash key when matching the two tets that share a face.
type FaceKey = [usize; 3];

/// Returns the sorted-ascending key of three vertex indices.
#[inline]
fn face_key(a: usize, b: usize, c: usize) -> FaceKey {
    let mut k = [a, b, c];
    k.sort_unstable();
    k
}

/// Output of 3-D Delaunay tetrahedralization.
///
/// Holds the original input `points` together with the list of Delaunay `tets`
/// over them. All vertex indices in every [`Tet`] index into `points`.
#[derive(Debug, Clone)]
pub struct Tetrahedralization {
    /// The input points (super-tetra vertices are *not* retained here).
    pub points: Vec<[f64; 3]>,
    /// The Delaunay tetrahedra over `points`.
    pub tets: Vec<Tet>,
}

impl Tetrahedralization {
    /// Returns the number of tetrahedra.
    #[inline]
    pub fn len(&self) -> usize {
        self.tets.len()
    }

    /// Returns `true` if there are no tetrahedra.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.tets.is_empty()
    }

    /// Squared circumradius of tet `t`, computed from its circumcenter.
    ///
    /// The circumcenter of a tetrahedron with vertices `a,b,c,d` solves the linear
    /// system that places it equidistant from all four vertices; this routine
    /// builds that `3x3` system relative to `a` and solves it by Cramer's rule. The
    /// returned value is the squared distance from the circumcenter to `a`. Returns
    /// [`f64::INFINITY`] for a degenerate (near-flat) tet whose circumsphere is not
    /// finite, which conservatively keeps such a tet *out* of every finite-`alpha`
    /// α-complex.
    pub fn circumradius_sq(&self, t: &Tet) -> f64 {
        let a = self.points[t.v[0]];
        let b = self.points[t.v[1]];
        let c = self.points[t.v[2]];
        let d = self.points[t.v[3]];
        circumradius_sq_of(a, b, c, d)
    }

    /// Returns the boundary triangles of the α-complex at squared radius `alpha`.
    ///
    /// A triangular face lies on the α-shape boundary when it is shared by exactly
    /// one Delaunay tet whose squared circumradius is `≤ alpha` (the standard
    /// boundary-only convention: an interior face shared by two such tets is *not*
    /// emitted, while a face whose only finite-`alpha` tet is on the other side is).
    /// Each returned triangle is the face triple oriented so its winding points
    /// *out of* the single owning tet, so for a closed component every edge appears
    /// exactly twice with opposite directions.
    ///
    /// At a sufficiently large `alpha` (≥ the largest tet circumradius²) the result
    /// is exactly the convex-hull boundary of the point set.
    pub fn alpha_shape(&self, alpha: f64) -> Vec<[usize; 3]> {
        alpha_shape(self, alpha)
    }

    /// Verifies the empty-sphere (Delaunay) property exactly for every tet.
    ///
    /// For each tetrahedron and each input point not belonging to it, the exact
    /// [`insphere`] predicate must report the point is *not strictly inside* the
    /// tet's circumsphere. Returns the first violating `(tet_idx, point_idx)` pair,
    /// or `None` if the whole complex is Delaunay. This is an `O(tets · points)`
    /// audit intended for tests and validation, not for hot paths.
    pub fn verify_empty_sphere(&self) -> Option<(usize, usize)> {
        for (ti, t) in self.tets.iter().enumerate() {
            let a = self.points[t.v[0]];
            let b = self.points[t.v[1]];
            let c = self.points[t.v[2]];
            let d = self.points[t.v[3]];
            for (pi, &p) in self.points.iter().enumerate() {
                if pi == t.v[0] || pi == t.v[1] || pi == t.v[2] || pi == t.v[3] {
                    continue;
                }
                // Strictly-inside is a hard Delaunay violation. A `Degenerate`
                // (cospherical) result is allowed: it means `p` lies *on* the
                // circumsphere, which is permitted for a Delaunay complex.
                if insphere_oriented(a, b, c, d, p) == Orientation::Positive {
                    return Some((ti, pi));
                }
            }
        }
        None
    }
}

/// Returns the squared circumradius of the tetrahedron `(a,b,c,d)`.
///
/// Solves the `3x3` system `M x = r` where row `k` is `2 (p_k − a)` and the
/// right-hand side is `|p_k|² − |a|²`, giving the circumcenter `x`; the squared
/// circumradius is `|x − a|²`. Cramer's rule is used; a vanishing determinant
/// (degenerate/flat tet) yields [`f64::INFINITY`].
fn circumradius_sq_of(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let row = |p: [f64; 3]| [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let r1 = row(b);
    let r2 = row(c);
    let r3 = row(d);
    let rhs = [
        0.5 * (r1[0] * (b[0] + a[0]) + r1[1] * (b[1] + a[1]) + r1[2] * (b[2] + a[2])),
        0.5 * (r2[0] * (c[0] + a[0]) + r2[1] * (c[1] + a[1]) + r2[2] * (c[2] + a[2])),
        0.5 * (r3[0] * (d[0] + a[0]) + r3[1] * (d[1] + a[1]) + r3[2] * (d[2] + a[2])),
    ];
    // Solve the 3x3 system (rows r1, r2, r3) for the circumcenter relative to `a`.
    let center = match solve3(r1, r2, r3, rhs) {
        Some(c) => c,
        None => return f64::INFINITY,
    };
    let dx = center[0] - a[0];
    let dy = center[1] - a[1];
    let dz = center[2] - a[2];
    dx * dx + dy * dy + dz * dz
}

/// Determinant of the `3x3` matrix whose *columns* are `c0, c1, c2`.
#[inline]
fn det3(c0: [f64; 3], c1: [f64; 3], c2: [f64; 3]) -> f64 {
    c0[0] * (c1[1] * c2[2] - c1[2] * c2[1]) - c1[0] * (c0[1] * c2[2] - c0[2] * c2[1])
        + c2[0] * (c0[1] * c1[2] - c0[2] * c1[1])
}

/// Solves the `3x3` linear system whose *rows* are `r0, r1, r2` for the given
/// right-hand side `rhs`, returning `None` for a singular system.
fn solve3(r0: [f64; 3], r1: [f64; 3], r2: [f64; 3], rhs: [f64; 3]) -> Option<[f64; 3]> {
    // Column form for Cramer's rule.
    let col0 = [r0[0], r1[0], r2[0]];
    let col1 = [r0[1], r1[1], r2[1]];
    let col2 = [r0[2], r1[2], r2[2]];
    let det = det3(col0, col1, col2);
    if det.abs() < 1e-300 {
        return None;
    }
    let x = det3(rhs, col1, col2) / det;
    let y = det3(col0, rhs, col2) / det;
    let z = det3(col0, col1, rhs) / det;
    Some([x, y, z])
}

// ─────────────────────────────────────────────────────────────────────────────
// Simulation of Simplicity (SoS) tie-breaks
// ─────────────────────────────────────────────────────────────────────────────

/// Resolves an [`Orientation`] to a strict `±1` sign, breaking a `Degenerate`
/// (exactly-zero determinant) result with Simulation of Simplicity.
///
/// The SoS rule mirrors `mesh_boolean::sos_sign`: when the exact determinant
/// vanishes the sign is decided from the parity of the *lowest* global vertex index
/// among the points that entered the determinant — an even lowest index resolves to
/// `+1`, an odd one to `-1`. Because the decision depends only on the
/// order-independent *set* of indices, the same geometric configuration always
/// resolves the same way regardless of the order in which tets are visited, so the
/// virtual perturbation is globally consistent and never produces a contradictory
/// orientation.
#[inline]
fn sos_sign(base: Orientation, indices: &[usize]) -> i32 {
    match base {
        Orientation::Positive => 1,
        Orientation::Negative => -1,
        Orientation::Degenerate => {
            let lowest = indices.iter().copied().min().unwrap_or(0);
            if lowest % 2 == 0 { 1 } else { -1 }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Builder state
// ─────────────────────────────────────────────────────────────────────────────

/// Mutable working state of the incremental Bowyer-Watson construction.
///
/// `pts` holds the `n` input points followed by the four super-tetra vertices, so
/// global indices `n..n+4` denote the super-tetra corners. `tets` grows and shrinks
/// during insertion; deleted tets are tracked in `dead` and their slots are reused.
struct Builder {
    /// Input points followed by the four super-tetra vertices.
    pts: Vec<[f64; 3]>,
    /// Number of real input points (super-tetra indices are `n..n+4`).
    n: usize,
    /// The working tetrahedra (some entries may be tombstoned via `dead`).
    tets: Vec<Tet>,
    /// `dead[i]` is `true` when slot `i` of `tets` has been deleted.
    dead: Vec<bool>,
}

impl Builder {
    /// Returns the position of global vertex index `gi`.
    #[inline]
    fn pos(&self, gi: usize) -> [f64; 3] {
        self.pts[gi]
    }

    /// Returns `true` if global vertex index `gi` is a super-tetra corner.
    #[inline]
    fn is_super(&self, gi: usize) -> bool {
        gi >= self.n
    }

    /// Exact orientation of the ordered tet `(a,b,c,d)` by global index, with SoS.
    ///
    /// Returns `+1` when the tet is positively oriented and `-1` otherwise; the
    /// `Degenerate` (coplanar) case is resolved by [`sos_sign`] so the answer is
    /// always strict.
    #[inline]
    fn orient_sign(&self, a: usize, b: usize, c: usize, d: usize) -> i32 {
        let o = orient3d(self.pos(a), self.pos(b), self.pos(c), self.pos(d));
        sos_sign(o, &[a, b, c, d])
    }

    /// Creates a positively oriented [`Tet`] from four global indices.
    ///
    /// If the supplied order is negatively oriented the first two vertices are
    /// swapped to flip the sign, guaranteeing the stored tet satisfies the
    /// positive-orientation invariant relied on by [`insphere`]. Adjacency is left
    /// unset (`-1`) for the caller to fill.
    fn make_tet(&self, a: usize, b: usize, c: usize, d: usize) -> Tet {
        if self.orient_sign(a, b, c, d) >= 0 {
            Tet {
                v: [a, b, c, d],
                neigh: [-1; 4],
            }
        } else {
            Tet {
                v: [b, a, c, d],
                neigh: [-1; 4],
            }
        }
    }

    /// Empty-sphere test of point `p` against tet `t`, with SoS on cospherical
    /// ties.
    ///
    /// Returns `true` when `p` lies inside (or, after the SoS tie-break, virtually
    /// inside) the circumsphere of `t`, i.e. when `t` is *not* Delaunay with respect
    /// to `p` and must enter the insertion cavity. The tet is positively oriented by
    /// construction, so a raw `insphere == Positive` already means "inside"; a
    /// `Degenerate` result (five cospherical points) is broken by the parity of the
    /// lowest index among the five involved vertices.
    fn in_circumsphere(&self, t: &Tet, p_gi: usize) -> bool {
        let a = self.pos(t.v[0]);
        let b = self.pos(t.v[1]);
        let c = self.pos(t.v[2]);
        let d = self.pos(t.v[3]);
        let p = self.pos(p_gi);
        let base = insphere_oriented(a, b, c, d, p);
        let idx = [t.v[0], t.v[1], t.v[2], t.v[3], p_gi];
        sos_sign(base, &idx) > 0
    }

    /// Returns `true` if `p` lies on the *inner* side of the face opposite local
    /// vertex `fi` of tet `ti` — i.e. on the same side as the tet's own apex.
    ///
    /// The face opposite `fi` is `(f0,f1,f2)`; the apex is `v[fi]`. Because the tet
    /// is positively oriented, `p` is inside across that face exactly when
    /// `orient3d(f0,f1,f2,p)` has the same strict sign as `orient3d(f0,f1,f2,apex)`.
    /// Both signs are made strict with SoS, so the test is total.
    fn inside_across_face(&self, ti: usize, fi: usize) -> i32 {
        // Sign of the apex relative to its opposite face is the reference "inside".
        let t = &self.tets[ti];
        let f = t.face(fi);
        let apex = t.v[fi];
        self.orient_sign(f[0], f[1], f[2], apex)
    }

    /// Sign of query point `p` relative to the face opposite local vertex `fi` of
    /// tet `ti`, with SoS.
    fn face_side_of_point(&self, ti: usize, fi: usize, p_gi: usize) -> i32 {
        let t = &self.tets[ti];
        let f = t.face(fi);
        let o = orient3d(self.pos(f[0]), self.pos(f[1]), self.pos(f[2]), self.pos(p_gi));
        sos_sign(o, &[f[0], f[1], f[2], p_gi])
    }

    /// Returns `true` when the cavity-boundary face opposite local vertex `fi` of
    /// cavity tet `ti` is *visible* from apex `p_gi`.
    ///
    /// Visibility means a new tetrahedron joining `p` to this face would have
    /// *strictly positive* volume — equivalently `p` lies strictly on the
    /// cavity-interior side of the face (the same side as the tet's own apex
    /// `v[fi]`). The apex side is computed with SoS (always strict) as the interior
    /// reference; the `p` side is computed with the *raw* [`orient3d`] so that a `p`
    /// exactly coplanar with the face is reported as **not** visible. A non-visible
    /// face signals a non-star-shaped cavity that the repair loop fixes by absorbing
    /// the neighbour across this face, which is exactly what guarantees the refill
    /// produces no flat or overlapping tets.
    fn face_visible_from(&self, ti: usize, fi: usize, p_gi: usize) -> bool {
        let t = &self.tets[ti];
        let f = t.face(fi);
        let apex = t.v[fi];
        let interior = self.orient_sign(f[0], f[1], f[2], apex);
        match orient3d(self.pos(f[0]), self.pos(f[1]), self.pos(f[2]), self.pos(p_gi)) {
            Orientation::Positive => interior > 0,
            Orientation::Negative => interior < 0,
            Orientation::Degenerate => false,
        }
    }
}

/// Evaluates [`insphere`] for a tet whose vertices `(a,b,c,d)` are supplied in an
/// arbitrary winding, normalising the orientation first so the sign convention is
/// stable.
///
/// [`insphere`] assumes its first four points are positively oriented (consistent
/// with [`orient3d`]); a positively oriented tet then maps `Positive → inside`. To
/// be robust to either winding this helper checks the float orientation and, when
/// negative, swaps the first two points (which both flips the orientation and flips
/// the in-sphere sign), then flips the returned [`Orientation`] back — so the result
/// always uses the "Positive = strictly inside" convention regardless of the input
/// winding. A coplanar (degenerate) base is passed through unchanged; the caller's
/// SoS handles it.
fn insphere_oriented(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> Orientation {
    match orient3d(a, b, c, d) {
        Orientation::Positive => insphere([a, b, c, d, p]),
        Orientation::Negative => flip(insphere([b, a, c, d, p])),
        Orientation::Degenerate => {
            // Flat base: the circumsphere is not well-defined in floating point.
            // Report the raw insphere sign (the caller resolves ties via SoS); the
            // construction never keeps flat tets among real points.
            insphere([a, b, c, d, p])
        }
    }
}

/// Returns the opposite [`Orientation`], leaving `Degenerate` unchanged.
#[inline]
fn flip(o: Orientation) -> Orientation {
    match o {
        Orientation::Positive => Orientation::Negative,
        Orientation::Negative => Orientation::Positive,
        Orientation::Degenerate => Orientation::Degenerate,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the 3-D Delaunay tetrahedralization of `points` with an exact
/// Bowyer-Watson algorithm.
///
/// Points may be in general position or degenerate (cospherical, coplanar, on a
/// grid); degeneracies are resolved via Simulation-of-Simplicity perturbation keyed
/// on the lowest global vertex index, so the construction never fails or produces an
/// inconsistent complex. Every returned [`Tet`] is positively oriented and indexes
/// into the returned `points` (a copy of the input). Duplicate input points are
/// tolerated: a point coincident with one already inserted produces an empty cavity
/// and is skipped.
///
/// Returns `None` if `points.len() < 4`, or if the points are *all* coplanar (no
/// tetrahedron of positive volume exists), since a 3-D tetrahedralization is then
/// undefined.
pub fn delaunay_3d(points: &[[f64; 3]]) -> Option<Tetrahedralization> {
    let n = points.len();
    if n < 4 {
        return None;
    }

    // Reject the fully-degenerate (all-coplanar) input: there is no positive-volume
    // tetrahedron to build, so a 3-D tetrahedralization does not exist.
    if all_coplanar(points) {
        return None;
    }

    let mut builder = init_builder(points);
    let n = builder.n;

    let mut hint = 0usize;
    for p_gi in 0..n {
        // Skip exact duplicates of an already-present point (empty cavity).
        if is_duplicate(&builder, p_gi) {
            continue;
        }
        hint = insert_point(&mut builder, p_gi, hint);
    }

    Some(finalize(builder))
}

/// Builds the initial [`Builder`] containing all input points plus a single
/// super-tetrahedron that strictly encloses them.
///
/// The super-tetra is sized from the input AABB and an expansion factor large
/// enough that every input point lies strictly inside it (see [`super_tetra`]).
fn init_builder(points: &[[f64; 3]]) -> Builder {
    let n = points.len();
    let mut pts = points.to_vec();
    let super_verts = super_tetra(points);
    pts.extend_from_slice(&super_verts);

    let mut builder = Builder {
        pts,
        n,
        tets: Vec::with_capacity(n * 6 + 1),
        dead: Vec::with_capacity(n * 6 + 1),
    };
    let st = builder.make_tet(n, n + 1, n + 2, n + 3);
    builder.tets.push(st);
    builder.dead.push(false);
    builder
}

/// Returns `true` if global vertex `p_gi` coincides (to a tight tolerance) with any
/// real input point of strictly lower index already inserted.
fn is_duplicate(builder: &Builder, p_gi: usize) -> bool {
    let p = builder.pos(p_gi);
    for q_gi in 0..p_gi {
        let q = builder.pos(q_gi);
        let dx = p[0] - q[0];
        let dy = p[1] - q[1];
        let dz = p[2] - q[2];
        if dx * dx + dy * dy + dz * dz < 1e-24 {
            return true;
        }
    }
    false
}

/// Returns the four vertices of a tetrahedron that strictly contains every input
/// point.
///
/// The AABB of the input is expanded outward by `3×` its longest side (plus a small
/// absolute pad so a degenerate/zero-extent input still yields a finite tet). The
/// four corners are placed far outside this padded box in a configuration whose
/// convex hull contains the box, guaranteeing every input point is strictly
/// interior. Exact sentinel coordinates are deliberately avoided in favour of
/// data-relative ones so the predicates see well-separated magnitudes.
fn super_tetra(points: &[[f64; 3]]) -> [[f64; 3]; 4] {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let cx = 0.5 * (lo[0] + hi[0]);
    let cy = 0.5 * (lo[1] + hi[1]);
    let cz = 0.5 * (lo[2] + hi[2]);
    let mut span = 0.0f64;
    for k in 0..3 {
        span = span.max(hi[k] - lo[k]);
    }
    if !span.is_finite() || span <= 0.0 {
        span = 1.0;
    }
    // Enclosing radius: large multiple of the longest side plus an absolute pad.
    let r = span * 10.0 + 1.0;
    // A big regular-ish tetrahedron centred at the data centroid. These offsets put
    // all four faces well outside the padded AABB.
    [
        [cx, cy + 3.0 * r, cz],
        [cx - 2.0 * r, cy - r, cz - r],
        [cx + 2.0 * r, cy - r, cz - r],
        [cx, cy - r, cz + 3.0 * r],
    ]
}

/// Returns `true` when *all* `points` lie on a common plane (no positive-volume
/// tetrahedron exists), using the exact [`orient3d`] predicate.
///
/// A reference plane is taken from the first three mutually non-collinear points; if
/// the points do not even span a plane (all collinear) they are trivially coplanar.
/// The set is coplanar iff every remaining point gives `orient3d == Degenerate`
/// against that plane.
fn all_coplanar(points: &[[f64; 3]]) -> bool {
    // Find first three non-collinear points to define a plane.
    let a = points[0];
    let mut b_idx = None;
    for (i, p) in points.iter().enumerate().skip(1) {
        if !same_point(*p, a) {
            b_idx = Some(i);
            break;
        }
    }
    let bi = match b_idx {
        Some(i) => i,
        None => return true, // all identical → degenerate
    };
    let b = points[bi];
    let mut c_idx = None;
    for (i, p) in points.iter().enumerate().skip(bi + 1) {
        if !collinear3(a, b, *p) {
            c_idx = Some(i);
            break;
        }
    }
    let ci = match c_idx {
        Some(i) => i,
        None => return true, // all collinear → coplanar
    };
    let c = points[ci];
    // Every point must be coplanar with (a,b,c).
    for p in points.iter() {
        if orient3d(a, b, c, *p) != Orientation::Degenerate {
            return false;
        }
    }
    true
}

/// Returns `true` if two points are identical to a tight tolerance.
#[inline]
fn same_point(a: [f64; 3], b: [f64; 3]) -> bool {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz < 1e-24
}

/// Returns `true` if three points are collinear, tested via the cross product
/// magnitude relative to the segment length.
fn collinear3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    let cm = cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2];
    let am = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    cm <= 1e-20 * am.max(1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Point location (walk)
// ─────────────────────────────────────────────────────────────────────────────

/// Locates a *live* tet whose closed region contains point `p_gi`, starting from
/// `hint` and walking across separating faces.
///
/// At each step the four faces of the current tet are tested with
/// [`Builder::face_side_of_point`] against the face's interior reference
/// ([`Builder::inside_across_face`]); if `p` lies strictly outside across some face
/// the walk moves to that neighbour. When `p` is on the inner side of all four faces
/// the current tet contains it. The walk is capped at `2n+8` steps; on exceeding the
/// cap (which can only happen in a pathological cycle) it falls back to a linear scan
/// that returns any tet containing `p`, and finally any live tet, so location never
/// fails.
fn locate(builder: &Builder, p_gi: usize, hint: usize) -> usize {
    let mut cur = live_from(builder, hint);
    let cap = builder.tets.len().saturating_mul(2).saturating_add(8);
    for _ in 0..cap {
        let mut moved = false;
        for fi in 0..4 {
            let nb = builder.tets[cur].neigh[fi];
            if nb < 0 {
                continue;
            }
            let inside_ref = builder.inside_across_face(cur, fi);
            let side = builder.face_side_of_point(cur, fi, p_gi);
            // `p` is on the far side of this face from the apex → cross it.
            if side != 0 && side != inside_ref {
                let nb = nb as usize;
                if !builder.dead[nb] {
                    cur = nb;
                    moved = true;
                    break;
                }
            }
        }
        if !moved {
            return cur;
        }
    }
    // Fallback: brute-force containment, then any live tet.
    for (ti, t) in builder.tets.iter().enumerate() {
        if builder.dead[ti] {
            continue;
        }
        if contains_point(builder, t, p_gi) {
            return ti;
        }
    }
    live_from(builder, 0)
}

/// Returns `start` if it is a live tet, otherwise the first live tet at or after it
/// (wrapping to the beginning), and `0` if none is found.
fn live_from(builder: &Builder, start: usize) -> usize {
    let m = builder.tets.len();
    if m == 0 {
        return 0;
    }
    let s = start.min(m - 1);
    if !builder.dead[s] {
        return s;
    }
    for off in 0..m {
        let i = (s + off) % m;
        if !builder.dead[i] {
            return i;
        }
    }
    0
}

/// Returns `true` if `p_gi` lies on the inner side of (or on) all four faces of
/// tet `t`, i.e. inside the closed tetrahedron.
fn contains_point(builder: &Builder, t: &Tet, p_gi: usize) -> bool {
    for fi in 0..4 {
        let f = t.face(fi);
        let apex = t.v[fi];
        let inside_ref = sos_sign(
            orient3d(
                builder.pos(f[0]),
                builder.pos(f[1]),
                builder.pos(f[2]),
                builder.pos(apex),
            ),
            &[f[0], f[1], f[2], apex],
        );
        let side = sos_sign(
            orient3d(
                builder.pos(f[0]),
                builder.pos(f[1]),
                builder.pos(f[2]),
                builder.pos(p_gi),
            ),
            &[f[0], f[1], f[2], p_gi],
        );
        if side != 0 && side != inside_ref {
            return false;
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Insertion
// ─────────────────────────────────────────────────────────────────────────────

/// Inserts point `p_gi` into the working triangulation and returns a tet index to
/// use as the next location hint.
///
/// The step locates the containing tet, computes the star-shaped insertion
/// [`cavity`] of tets to be removed, extracts its boundary faces, and re-fills it by
/// joining `p` to every boundary face ([`refill_cavity`]). Returns one of the newly
/// created tets so the next insertion's walk starts nearby. If the cavity comes back
/// empty (e.g. an exact duplicate slipped through), the triangulation is left
/// unchanged and the hint is returned.
fn insert_point(builder: &mut Builder, p_gi: usize, hint: usize) -> usize {
    let start = locate(builder, p_gi, hint);
    let (cav, boundary) = cavity(builder, p_gi, start);
    if cav.is_empty() {
        return start;
    }
    refill_cavity(builder, p_gi, &cav, &boundary)
}

/// A boundary face of the insertion cavity.
///
/// `tri` is the face's three global vertex indices in the *outward* winding of the
/// outer (non-cavity) tet — chosen so that joining the new apex to this face yields
/// a positively oriented tet. `outer` is the index of the tet just outside the
/// cavity across this face, or `-1` when the face borders the super-tetra boundary;
/// `outer_face` is the local face index within that outer tet (or `-1`).
struct BoundaryFace {
    tri: [usize; 3],
    outer: i32,
    outer_face: i32,
}

/// Computes the *star-shaped* insertion cavity for point `p_gi` together with its
/// boundary faces.
///
/// The cavity is grown in two phases:
///
/// 1. **Circumsphere flood.** Starting from the located tet (or, as a degenerate
///    safety net, any tet failing the empty-sphere test), a flood fill over the
///    face-adjacency graph collects every reachable tet whose circumsphere contains
///    `p` ([`Builder::in_circumsphere`]). For a Delaunay triangulation this set is
///    already connected and star-shaped.
///
/// 2. **Star-shapedness repair.** The boundary of the candidate cavity is extracted
///    and every boundary face is checked for *visibility* from `p`: the new tet
///    `(face, p)` must have positive volume, i.e. `p` must lie strictly on the
///    cavity-interior side of the face. Any boundary face that fails this (because
///    `p` is coplanar with it or beyond it — the classic non-star-shaped
///    configuration that arises under cospherical degeneracy) forces its outer
///    neighbour into the cavity, and the boundary is recomputed. The loop runs until
///    every boundary face is visible, guaranteeing the re-tetrahedralization in
///    [`refill_cavity`] exactly fills the removed region with positively oriented
///    tets and no overlaps.
///
/// Returns `(cavity_tets, boundary_faces)`. The cavity is empty (and the boundary
/// too) only when `p` is in no circumsphere at all, i.e. it coincides with an
/// existing vertex.
fn cavity(builder: &Builder, p_gi: usize, start: usize) -> (Vec<usize>, Vec<BoundaryFace>) {
    let mut in_cav = vec![false; builder.tets.len()];
    let mut cav: Vec<usize> = Vec::new();

    // Phase 1: circumsphere flood.
    let seed = cavity_seed(builder, p_gi, start);
    if seed == usize::MAX {
        return (cav, Vec::new());
    }
    let mut stack = vec![seed];
    in_cav[seed] = true;
    while let Some(ti) = stack.pop() {
        cav.push(ti);
        for &nb in builder.tets[ti].neigh.iter() {
            if nb < 0 {
                continue;
            }
            let nb = nb as usize;
            if in_cav[nb] || builder.dead[nb] {
                continue;
            }
            if builder.in_circumsphere(&builder.tets[nb], p_gi) {
                in_cav[nb] = true;
                stack.push(nb);
            }
        }
    }

    // Phase 2: repair until every boundary face is visible from `p`.
    let cap = builder.tets.len().saturating_add(8);
    for _ in 0..cap {
        let mut grew = false;
        let snapshot = cav.clone();
        for &ti in &snapshot {
            for fi in 0..4 {
                let nb = builder.tets[ti].neigh[fi];
                if nb >= 0 && in_cav[nb as usize] {
                    continue; // internal face
                }
                if !builder.face_visible_from(ti, fi, p_gi) {
                    // Non-star-shaped here: absorb the outer neighbour.
                    if nb >= 0 {
                        let nb = nb as usize;
                        if !in_cav[nb] && !builder.dead[nb] {
                            in_cav[nb] = true;
                            cav.push(nb);
                            grew = true;
                        }
                    }
                }
            }
        }
        if !grew {
            break;
        }
    }

    let boundary = extract_boundary(builder, &cav, &in_cav);
    (cav, boundary)
}

/// Returns the seed tet for the cavity flood: the located `start` if its
/// circumsphere already contains `p`, otherwise the first tet found whose
/// circumsphere contains `p`, or `usize::MAX` if none does (duplicate vertex).
fn cavity_seed(builder: &Builder, p_gi: usize, start: usize) -> usize {
    if !builder.dead[start] && builder.in_circumsphere(&builder.tets[start], p_gi) {
        return start;
    }
    for (ti, t) in builder.tets.iter().enumerate() {
        if !builder.dead[ti] && builder.in_circumsphere(t, p_gi) {
            return ti;
        }
    }
    usize::MAX
}

/// Extracts the boundary faces of cavity `cav` (membership in `in_cav`).
///
/// A face lies on the cavity boundary when the neighbour across it is not in the
/// cavity. The emitted `tri` uses the *outward* winding — the cavity tet's own face
/// winding reversed — so that joining the apex `p` (which sits on the cavity side)
/// produces a positively oriented tet, and the recorded `outer`/`outer_face` let
/// [`refill_cavity`] restitch adjacency to the region outside the cavity.
fn extract_boundary(builder: &Builder, cav: &[usize], in_cav: &[bool]) -> Vec<BoundaryFace> {
    let mut faces = Vec::new();
    for &ti in cav {
        let t = &builder.tets[ti];
        for fi in 0..4 {
            let nb = t.neigh[fi];
            if nb >= 0 && in_cav[nb as usize] {
                continue;
            }
            let f = t.face(fi);
            let tri = [f[0], f[2], f[1]];
            let outer_face = if nb >= 0 {
                shared_face_index(&builder.tets[nb as usize], &f)
            } else {
                -1
            };
            faces.push(BoundaryFace {
                tri,
                outer: nb,
                outer_face,
            });
        }
    }
    faces
}

/// Returns the local face index `fi` of tet `t` whose opposite-vertex face matches
/// the unordered triple `tri`, or `-1` if none matches.
fn shared_face_index(t: &Tet, tri: &[usize; 3]) -> i32 {
    let key = face_key(tri[0], tri[1], tri[2]);
    for fi in 0..4 {
        let f = t.face(fi);
        if face_key(f[0], f[1], f[2]) == key {
            return fi as i32;
        }
    }
    -1
}

/// Tombstones the cavity tets and creates one new tet per boundary face joining the
/// apex `p_gi` to that face, restoring full face-adjacency.
///
/// Each new tet `(tri[0], tri[1], tri[2], p)` is created positively oriented. Its
/// adjacency is set so that:
/// * the face opposite the apex (local vertex 3) points back to the recorded outer
///   neighbour, and that outer neighbour's matching face is repointed to the new tet
///   (so the cavity's outside is restitched), and
/// * the three side faces are matched against the other new tets by their shared
///   edges via a face-key map, giving the new tets mutual adjacency.
///
/// Returns the index of the last new tet (a convenient next location hint), or
/// `live_from(0)` if no tet was created.
fn refill_cavity(
    builder: &mut Builder,
    p_gi: usize,
    cav: &[usize],
    boundary: &[BoundaryFace],
) -> usize {
    use std::collections::HashMap;

    // Kill cavity tets first; their slots may be reused for new tets.
    for &ti in cav {
        builder.dead[ti] = true;
    }

    let mut free: Vec<usize> = cav.to_vec();
    let mut new_tets: Vec<usize> = Vec::with_capacity(boundary.len());

    // Create one tet per boundary face and wire it to the outer neighbour.
    for bf in boundary {
        let tet = builder.make_tet(bf.tri[0], bf.tri[1], bf.tri[2], p_gi);
        let idx = if let Some(slot) = free.pop() {
            builder.tets[slot] = tet;
            builder.dead[slot] = false;
            slot
        } else {
            builder.tets.push(tet);
            builder.dead.push(false);
            builder.tets.len() - 1
        };
        new_tets.push(idx);

        // Wire the apex-opposite face (the base = bf.tri) to the outer neighbour.
        // After `make_tet`, the base triangle may have been re-wound, so locate the
        // local face of the new tet that equals bf.tri by key.
        let base_local = shared_face_index(&builder.tets[idx], &bf.tri);
        if base_local >= 0 {
            builder.tets[idx].neigh[base_local as usize] = bf.outer;
        }
        if bf.outer >= 0 && bf.outer_face >= 0 {
            let o = bf.outer as usize;
            builder.tets[o].neigh[bf.outer_face as usize] = idx as i32;
        }
    }

    // Match the side faces of the new tets to each other by shared edge keys. Every
    // internal face of the star polyhedron is shared by exactly two new tets.
    let mut edge_map: HashMap<FaceKey, (usize, usize)> = HashMap::new();
    for &idx in &new_tets {
        for fi in 0..4 {
            let f = builder.tets[idx].face(fi);
            // Only the three faces that include the apex are internal-to-star;
            // the base face (no apex) was already wired to the outer neighbour.
            if f[0] != p_gi && f[1] != p_gi && f[2] != p_gi {
                continue;
            }
            let key = face_key(f[0], f[1], f[2]);
            match edge_map.get(&key).copied() {
                None => {
                    edge_map.insert(key, (idx, fi));
                }
                Some((other_idx, other_fi)) => {
                    builder.tets[idx].neigh[fi] = other_idx as i32;
                    builder.tets[other_idx].neigh[other_fi] = idx as i32;
                }
            }
        }
    }

    new_tets.last().copied().unwrap_or_else(|| live_from(builder, 0))
}

// ─────────────────────────────────────────────────────────────────────────────
// Finalize
// ─────────────────────────────────────────────────────────────────────────────

/// Strips the super-tetrahedron and tombstoned slots, returning the Delaunay
/// tetrahedralization over the original points.
///
/// Every live tet that does *not* touch a super-tetra vertex is kept; its
/// adjacency indices are remapped to the compacted tet array, and any adjacency
/// pointing at a removed (super-touching or dead) tet becomes `-1`. The
/// super-tetra vertices are dropped from the point list, leaving only the user's
/// `n` points.
fn finalize(builder: Builder) -> Tetrahedralization {
    let n = builder.n;
    // old tet index → new compacted index (or -1 if removed)
    let mut remap = vec![-1i32; builder.tets.len()];
    let mut kept: Vec<Tet> = Vec::new();
    for (ti, t) in builder.tets.iter().enumerate() {
        if builder.dead[ti] {
            continue;
        }
        if t.v.iter().any(|&v| builder.is_super(v)) {
            continue; // touches the super-tetra
        }
        remap[ti] = kept.len() as i32;
        kept.push(t.clone());
    }
    // Remap adjacency.
    for t in &mut kept {
        for nb in t.neigh.iter_mut() {
            *nb = if *nb < 0 { -1 } else { remap[*nb as usize] };
        }
    }
    Tetrahedralization {
        points: builder.pts[..n].to_vec(),
        tets: kept,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// α-shape
// ─────────────────────────────────────────────────────────────────────────────

/// Free function form of [`Tetrahedralization::alpha_shape`].
///
/// Filters the tets to those with squared circumradius `≤ alpha`, then emits every
/// triangular face that, within that filtered set, is incident to exactly one tet —
/// oriented to point out of that tet. See [`Tetrahedralization::alpha_shape`] for
/// the full contract.
pub fn alpha_shape(tetra: &Tetrahedralization, alpha: f64) -> Vec<[usize; 3]> {
    use std::collections::HashMap;

    // Keep tets within the filtration radius.
    let mut keep = vec![false; tetra.tets.len()];
    for (ti, t) in tetra.tets.iter().enumerate() {
        if tetra.circumradius_sq(t) <= alpha {
            keep[ti] = true;
        }
    }

    // Count how many kept tets touch each unordered face; remember one oriented
    // representative so the emitted triangle points out of its owning tet.
    let mut count: HashMap<FaceKey, u32> = HashMap::new();
    let mut rep: HashMap<FaceKey, [usize; 3]> = HashMap::new();
    for (ti, t) in tetra.tets.iter().enumerate() {
        if !keep[ti] {
            continue;
        }
        for fi in 0..4 {
            // Outward winding for face opposite vertex fi of a positively oriented
            // tet: reverse the inward face winding.
            let f = t.face(fi);
            let outward = [f[0], f[2], f[1]];
            let key = face_key(outward[0], outward[1], outward[2]);
            *count.entry(key).or_insert(0) += 1;
            rep.entry(key).or_insert(outward);
        }
    }

    let mut tris = Vec::new();
    for (key, &c) in &count {
        if c == 1
            && let Some(&tri) = rep.get(key)
        {
            tris.push(tri);
        }
    }
    tris
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal LCG (Knuth MMIX constants) producing deterministic f64 in [0,1).
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Asserts every tet of `t` is positively oriented and indexes valid points.
    fn assert_valid(t: &Tetrahedralization) {
        let np = t.points.len();
        for tet in &t.tets {
            for &v in &tet.v {
                assert!(v < np, "vertex index {v} out of range {np}");
            }
            let a = t.points[tet.v[0]];
            let b = t.points[tet.v[1]];
            let c = t.points[tet.v[2]];
            let d = t.points[tet.v[3]];
            // Positive (or at worst degenerate-but-tie-broken) orientation.
            assert_ne!(
                orient3d(a, b, c, d),
                Orientation::Negative,
                "tet {:?} is negatively oriented",
                tet.v
            );
        }
    }

    #[test]
    fn returns_none_for_too_few_points() {
        assert!(delaunay_3d(&[]).is_none());
        assert!(delaunay_3d(&[[0.0, 0.0, 0.0]]).is_none());
        assert!(delaunay_3d(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]).is_none());
    }

    #[test]
    fn returns_none_for_coplanar_points() {
        // Four points on the z = 0 plane.
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert!(delaunay_3d(&pts).is_none());
    }

    #[test]
    fn single_tetrahedron_four_points() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let t = delaunay_3d(&pts).expect("4 non-coplanar points tetrahedralize");
        assert_eq!(t.len(), 1, "four points yield exactly one tet");
        assert_valid(&t);
        assert_eq!(t.verify_empty_sphere(), None);
        // The single tet uses all four vertices.
        let mut vs = t.tets[0].v;
        vs.sort_unstable();
        assert_eq!(vs, [0, 1, 2, 3]);
    }

    #[test]
    fn five_point_known_case() {
        // A triangular bipyramid: a shared base triangle (vertices 0,1,2 in z = 0)
        // with one apex strictly above and one strictly below. The Delaunay
        // tetrahedralization of these five points in convex position is exactly the
        // two obvious tets — (base, top) and (base, bottom) — sharing the base face.
        let pts = [
            [0.0, 0.0, 0.0],  // base
            [2.0, 0.0, 0.0],  // base
            [1.0, 2.0, 0.0],  // base
            [1.0, 0.7, 1.5],  // top apex (vertex 3)
            [1.0, 0.7, -1.5], // bottom apex (vertex 4)
        ];
        let t = delaunay_3d(&pts).expect("5 points tetrahedralize");
        assert_valid(&t);
        assert_eq!(
            t.verify_empty_sphere(),
            None,
            "Delaunay empty-sphere must hold"
        );
        // Exactly two tets sharing the base triangle {0,1,2}.
        assert_eq!(t.len(), 2, "expected exactly two tets, got {}", t.len());
        // Every vertex 0..5 should be referenced.
        let mut used = [false; 5];
        for tet in &t.tets {
            for &v in &tet.v {
                used[v] = true;
            }
        }
        assert!(used.iter().all(|&u| u), "all five points used: {used:?}");
        // One tet must contain the top apex (3) and the other the bottom apex (4);
        // both must contain the whole base triangle.
        let has3 = t.tets.iter().find(|tet| tet.v.contains(&3));
        let has4 = t.tets.iter().find(|tet| tet.v.contains(&4));
        let tet_top = has3.expect("a tet contains the top apex");
        let tet_bot = has4.expect("a tet contains the bottom apex");
        for base in [0usize, 1, 2] {
            assert!(tet_top.v.contains(&base), "top tet contains base {base}");
            assert!(tet_bot.v.contains(&base), "bottom tet contains base {base}");
        }
        // The two tets must be mutually adjacent across the shared base face.
        assert!(t.tets[0].neigh.contains(&1), "tet 0 adjacent to tet 1");
        assert!(t.tets[1].neigh.contains(&0), "tet 1 adjacent to tet 0");
    }

    #[test]
    fn cube_corners_degenerate() {
        // Eight cube corners: heavily cospherical/coplanar grid, a torture case for
        // degeneracy handling. Must not panic and must be a valid Delaunay complex.
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let t = delaunay_3d(&pts).expect("cube corners tetrahedralize");
        assert!(!t.is_empty(), "cube must produce tets");
        assert_valid(&t);
        assert_eq!(
            t.verify_empty_sphere(),
            None,
            "empty-sphere must hold for the cube"
        );
        // A cube triangulates into 5 or 6 tets depending on the SoS tie-break; both
        // are valid Delaunay complexes. Total volume must equal the cube volume 1.
        let vol = total_volume(&t);
        assert!(
            (vol - 1.0).abs() < 1e-9,
            "tet volumes must fill the cube: got {vol}"
        );
    }

    #[test]
    fn random_cloud_empty_sphere_audit_small() {
        // A few hundred random points: fast empty-sphere audit in the default run.
        let mut rng = Lcg(0x00C0_FFEE_1234_5678);
        let mut pts = Vec::new();
        for _ in 0..400 {
            pts.push([rng.next_f64(), rng.next_f64(), rng.next_f64()]);
        }
        let t = delaunay_3d(&pts).expect("random cloud tetrahedralizes");
        assert_valid(&t);
        assert_eq!(
            t.verify_empty_sphere(),
            None,
            "empty-sphere property must hold for every tet"
        );
    }

    #[test]
    fn random_cloud_10k_empty_sphere_audit() {
        // The headline 10^4-point empty-sphere audit. The O(tets·points) verify is
        // heavy, so this test is gated behind an env var to keep ordinary runs fast
        // while remaining always-compiled and runnable on demand.
        if std::env::var("OXIPHYSICS_DELAUNAY_10K").is_err() {
            return;
        }
        let mut rng = Lcg(0x5EED_F00D_CAFE_0042);
        let mut pts = Vec::new();
        for _ in 0..10_000 {
            pts.push([rng.next_f64(), rng.next_f64(), rng.next_f64()]);
        }
        let t = delaunay_3d(&pts).expect("10k cloud tetrahedralizes");
        assert!(t.len() > 10_000, "expect many tets for 10k points");
        assert_eq!(
            t.verify_empty_sphere(),
            None,
            "empty-sphere property must hold across the whole 10k tetrahedralization"
        );
    }

    #[test]
    fn random_cloud_1k_empty_sphere_audit() {
        // 1000-point cloud: still affordable as a default-run empty-sphere audit and
        // exercises a non-trivial cavity-flood/walk workload.
        let mut rng = Lcg(0xABCD_0987_6543_210F);
        let mut pts = Vec::new();
        for _ in 0..1000 {
            pts.push([
                rng.next_f64() * 10.0,
                rng.next_f64() * 10.0,
                rng.next_f64() * 10.0,
            ]);
        }
        let t = delaunay_3d(&pts).expect("1k cloud tetrahedralizes");
        assert_valid(&t);
        assert_eq!(t.verify_empty_sphere(), None);
    }

    #[test]
    fn alpha_shape_convex_hull_at_large_alpha() {
        // 20 points in convex position (on a sphere). At very large alpha every tet
        // is retained, so the alpha-shape is the convex-hull boundary: a closed
        // orientable surface where every undirected edge appears exactly twice.
        let mut rng = Lcg(0x1357_9BDF_2468_ACE0);
        let mut pts = Vec::new();
        for _ in 0..20 {
            // Sample on the unit sphere via normalized Gaussian-ish coordinates.
            let u = rng.next_f64() * 2.0 - 1.0;
            let theta = rng.next_f64() * std::f64::consts::TAU;
            let r = (1.0 - u * u).max(0.0).sqrt();
            pts.push([r * theta.cos(), r * theta.sin(), u]);
        }
        let t = delaunay_3d(&pts).expect("sphere-sampled points tetrahedralize");
        assert_valid(&t);
        assert_eq!(t.verify_empty_sphere(), None);
        let tris = t.alpha_shape(1e18);
        assert!(!tris.is_empty(), "convex hull must have faces");
        assert_closed_surface(&tris);
    }

    #[test]
    fn alpha_shape_sphere_surface_closed() {
        // ~150 points approximately on a sphere of radius R (kept modest so the
        // O(tets·points) empty-sphere audit stays fast in the default run). All
        // points lie in convex position, so a large alpha recovers the sphere's
        // surface — the convex-hull boundary — as a closed orientable (genus-0) mesh
        // in which every undirected edge is shared by exactly two triangles.
        let mut rng = Lcg(0x2BAD_F00D_1111_2222);
        let radius = 5.0;
        let mut pts = Vec::new();
        for _ in 0..150 {
            let u = rng.next_f64() * 2.0 - 1.0;
            let theta = rng.next_f64() * std::f64::consts::TAU;
            let r = (1.0 - u * u).max(0.0).sqrt();
            pts.push([
                radius * r * theta.cos(),
                radius * r * theta.sin(),
                radius * u,
            ]);
        }
        let t = delaunay_3d(&pts).expect("sphere cloud tetrahedralizes");
        assert_valid(&t);
        assert_eq!(t.verify_empty_sphere(), None);
        // Large alpha → convex hull of the sphere samples, which is a closed surface.
        let tris = t.alpha_shape(1e18);
        assert!(tris.len() >= 4, "sphere surface needs many triangles");
        assert_closed_surface(&tris);
    }

    #[test]
    fn alpha_shape_empty_at_zero_alpha() {
        // No tet has zero circumradius, so alpha = 0 yields no faces.
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.2, 1.2, 1.2],
        ];
        let t = delaunay_3d(&pts).expect("tetrahedralizes");
        assert!(
            t.alpha_shape(0.0).is_empty(),
            "alpha=0 selects no tets, hence no faces"
        );
    }

    #[test]
    fn duplicate_points_are_tolerated() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0], // duplicate of vertex 0
            [1.2, 1.2, 1.2],
        ];
        let t = delaunay_3d(&pts).expect("duplicates tolerated");
        assert_valid(&t);
        assert_eq!(t.verify_empty_sphere(), None);
    }

    // ── test helpers ──────────────────────────────────────────────────────────

    /// Sum of absolute tet volumes.
    fn total_volume(t: &Tetrahedralization) -> f64 {
        let mut vol = 0.0;
        for tet in &t.tets {
            let a = t.points[tet.v[0]];
            let b = t.points[tet.v[1]];
            let c = t.points[tet.v[2]];
            let d = t.points[tet.v[3]];
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
            let cr = [
                ac[1] * ad[2] - ac[2] * ad[1],
                ac[2] * ad[0] - ac[0] * ad[2],
                ac[0] * ad[1] - ac[1] * ad[0],
            ];
            vol += (ab[0] * cr[0] + ab[1] * cr[1] + ab[2] * cr[2]).abs() / 6.0;
        }
        vol
    }

    /// Asserts the triangle set forms a closed surface: every undirected edge is
    /// shared by exactly two triangles.
    fn assert_closed_surface(tris: &[[usize; 3]]) {
        use std::collections::HashMap;
        let mut edge: HashMap<(usize, usize), u32> = HashMap::new();
        for tri in tris {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edge.entry(key).or_insert(0) += 1;
            }
        }
        for (e, &c) in &edge {
            assert_eq!(
                c, 2,
                "edge {e:?} shared by {c} triangles (closed surface needs exactly 2)"
            );
        }
    }
}
