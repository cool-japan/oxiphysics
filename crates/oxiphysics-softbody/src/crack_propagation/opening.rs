// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topological crack opening for [`SoftBodyMesh`].
//!
//! Splitting an edge at its midpoint only refines the mesh: `a - mid - b` is
//! still connected and the triangles still reference `a - b`.  A crack has to
//! *separate* the material, so [`open_crack`] cuts the mesh topology along a
//! crack segment:
//!
//! 1. Every edge whose end nodes lie on opposite sides of the crack and whose
//!    crossing point lies on the segment is a **cut edge**.
//! 2. Every triangle with a cut edge is **cut** and joins the side of the crack
//!    that holds its centroid.  Two cut triangles on opposite sides never stay
//!    joined along a shared edge, while two cut triangles on the *same* side (for
//!    example both triangles of a quad whose diagonal the crack crosses) and cut
//!    triangles next to untouched ones stay joined.
//! 3. Around every node on such a boundary the incident triangles fall into
//!    *sectors*: two triangles belong to the same sector when they share an edge
//!    through the node and are not on opposite sides of the crack.  A node whose
//!    star has several sectors is split: the sector with the largest area keeps
//!    the node, every other sector gets a duplicate with the same position and
//!    velocity, and the nodal mass is shared out by sector area, so no mass is
//!    created.  (Triangle fans that merely touch at the node, a bow-tie, are
//!    different sectors as well.)
//! 4. Triangles are rewired to the duplicates of their sector and the edge list is
//!    brought in line: every triangle edge follows its triangle (an edge shared by
//!    the two sides becomes one edge per side), and cut edges without a triangle
//!    are removed.
//!
//! After the cut the two crack faces share no node.  A crack that spans the whole
//! body therefore leaves two [`identify_fragments`](super::identify_fragments);
//! a crack that ends inside the body opens a slit whose tip stays closed.

use super::SoftBodyMesh;
use std::collections::{HashMap, HashSet};

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

fn is_finite(a: V3) -> bool {
    a.iter().all(|c| c.is_finite())
}

/// Undirected edge key with the smaller node index first.
fn edge_key(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

/// What [`open_crack`] changed in the mesh.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CrackOpening {
    /// Number of cut edges whose connectivity changed: edges without a triangle
    /// that were removed, and triangle edges that were duplicated for the two crack
    /// faces or re-attached to duplicated nodes.  An edge that an earlier call
    /// already opened is not counted again.
    pub severed_edges: usize,
    /// `(original, duplicate)` for every node that was split.  Per-node data that
    /// lives outside the mesh (velocities, constraints, ...) can be mirrored from
    /// these pairs.
    pub duplicated_nodes: Vec<(usize, usize)>,
}

/// The ruled surface swept by a crack segment moving along the surface normal.
struct CrackPlane {
    origin: V3,
    direction: V3,
    length: f64,
    normal: V3,
}

impl CrackPlane {
    fn new(mesh: &SoftBodyMesh, from: V3, to: V3, surface_normal: Option<V3>) -> Option<Self> {
        if !is_finite(from) || !is_finite(to) {
            return None;
        }
        let along = sub(to, from);
        let length = norm(along);
        if !length.is_finite() || length <= 1e-300 {
            return None;
        }
        let direction = scale(along, 1.0 / length);
        let surface = surface_normal
            .filter(|n| is_finite(*n) && norm(*n) > 1e-12)
            .or_else(|| estimate_surface_normal(mesh, from, to))
            .unwrap_or_else(|| default_surface_normal(direction));
        // Only the part of the surface normal perpendicular to the crack matters.
        let in_plane = sub(surface, scale(direction, dot(surface, direction)));
        let in_plane_len = norm(in_plane);
        if in_plane_len < 1e-9 * norm(surface).max(1e-300) {
            // The crack would run along the surface normal.
            return None;
        }
        let surface = scale(in_plane, 1.0 / in_plane_len);
        Some(Self {
            origin: from,
            direction,
            length,
            normal: cross(surface, direction),
        })
    }

    /// Signed distance of `x` from the crack surface; `>= 0` is the positive side.
    fn side(&self, x: V3) -> f64 {
        dot(sub(x, self.origin), self.normal)
    }

    /// `true` if the edge between `xa` and `xb` crosses the crack segment.
    fn crosses(&self, xa: V3, xb: V3) -> bool {
        let sa = self.side(xa);
        let sb = self.side(xb);
        if (sa >= 0.0) == (sb >= 0.0) {
            return false;
        }
        let edge = sub(xb, xa);
        let hit = add(xa, scale(edge, sa / (sa - sb)));
        let along = dot(sub(hit, self.origin), self.direction);
        let tolerance = 1e-12 * self.length.max(norm(edge));
        along >= -tolerance && along <= self.length + tolerance
    }
}

/// Distance from `point` to the segment `from`-`to`.
fn point_segment_distance(point: V3, from: V3, to: V3) -> f64 {
    let seg = sub(to, from);
    let len_sq = dot(seg, seg);
    let t = if len_sq > 0.0 {
        (dot(sub(point, from), seg) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm(sub(point, add(from, scale(seg, t))))
}

/// Normal of the triangle whose centroid is closest to the segment.
fn estimate_surface_normal(mesh: &SoftBodyMesh, from: V3, to: V3) -> Option<V3> {
    let n = mesh.nodes.len();
    let mut best: Option<(f64, V3)> = None;
    for tri in &mesh.triangles {
        if tri.iter().any(|&i| i >= n) {
            continue;
        }
        let p = tri.map(|i| mesh.nodes[i].position);
        let normal = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        let magnitude = norm(normal);
        if !magnitude.is_finite() || magnitude <= 1e-300 {
            continue;
        }
        let centroid = scale(add(add(p[0], p[1]), p[2]), 1.0 / 3.0);
        let distance = point_segment_distance(centroid, from, to);
        if best.is_none_or(|(d, _)| distance < d) {
            best = Some((distance, normal));
        }
    }
    best.map(|(_, normal)| normal)
}

/// Surface normal for meshes without triangles: the coordinate axis that is least
/// aligned with the crack, preferring z (then y) so that a crack in the xy plane
/// uses the z axis.
fn default_surface_normal(direction: V3) -> V3 {
    let mut axis = 2;
    let mut least = direction[2].abs();
    for candidate in [1usize, 0] {
        if direction[candidate].abs() < least {
            least = direction[candidate].abs();
            axis = candidate;
        }
    }
    let mut normal = [0.0; 3];
    normal[axis] = 1.0;
    normal
}

/// Every undirected edge of the mesh: the edge list followed by the triangle
/// edges that are not in it, each once.
fn candidate_edges(mesh: &SoftBodyMesh) -> Vec<(usize, usize)> {
    let n = mesh.nodes.len();
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    let mut edges = Vec::new();
    let mut push = |a: usize, b: usize| {
        if a != b && a < n && b < n && seen.insert(edge_key(a, b)) {
            edges.push(edge_key(a, b));
        }
    };
    for &[a, b] in &mesh.edges {
        push(a, b);
    }
    for &[a, b, c] in &mesh.triangles {
        push(a, b);
        push(b, c);
        push(c, a);
    }
    edges
}

fn triangle_area(mesh: &SoftBodyMesh, tri: usize) -> f64 {
    let p = mesh.triangles[tri].map(|i| mesh.nodes[i].position);
    0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0])))
}

fn find_root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Side of the crack (`true` = positive) of every **cut** triangle, i.e. every
/// triangle with at least one cut edge, judged by its centroid.
fn cut_triangle_sides(
    mesh: &SoftBodyMesh,
    plane: &CrackPlane,
    cut: &HashSet<(usize, usize)>,
) -> HashMap<usize, bool> {
    let n = mesh.nodes.len();
    let mut sides = HashMap::new();
    for (t, &[a, b, c]) in mesh.triangles.iter().enumerate() {
        if a >= n || b >= n || c >= n {
            continue;
        }
        let touched = [(a, b), (b, c), (c, a)]
            .iter()
            .any(|&(x, y)| cut.contains(&edge_key(x, y)));
        if !touched {
            continue;
        }
        let centroid = scale(
            add(
                add(mesh.nodes[a].position, mesh.nodes[b].position),
                mesh.nodes[c].position,
            ),
            1.0 / 3.0,
        );
        sides.insert(t, plane.side(centroid) >= 0.0);
    }
    sides
}

/// `true` if both triangles are cut and lie on opposite sides of the crack.
fn on_opposite_sides(sides: &HashMap<usize, bool>, t1: usize, t2: usize) -> bool {
    matches!((sides.get(&t1), sides.get(&t2)), (Some(a), Some(b)) if a != b)
}

/// Edge -> incident triangles.
fn incident_triangles(triangles: &[[usize; 3]]) -> HashMap<(usize, usize), Vec<usize>> {
    let mut edge_triangles: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (t, &[a, b, c]) in triangles.iter().enumerate() {
        for (x, y) in [(a, b), (b, c), (c, a)] {
            let list = edge_triangles.entry(edge_key(x, y)).or_default();
            if list.last() != Some(&t) {
                list.push(t);
            }
        }
    }
    edge_triangles
}

/// Nodes at the end of an edge that separates triangles on opposite sides of the
/// crack, sorted.
fn interface_nodes(
    edge_triangles: &HashMap<(usize, usize), Vec<usize>>,
    sides: &HashMap<usize, bool>,
) -> Vec<usize> {
    let mut nodes: Vec<usize> = Vec::new();
    for (&(a, b), tris) in edge_triangles {
        let separates = tris.iter().enumerate().any(|(i, &t1)| {
            tris[i + 1..]
                .iter()
                .any(|&t2| on_opposite_sides(sides, t1, t2))
        });
        if separates {
            nodes.push(a);
            nodes.push(b);
        }
    }
    nodes.sort_unstable();
    nodes.dedup();
    nodes
}

/// `true` if triangles `t1` and `t2` (both containing `v`) share an edge through
/// `v` and are not on opposite sides of the crack.
fn stay_joined_at(
    mesh: &SoftBodyMesh,
    sides: &HashMap<usize, bool>,
    t1: usize,
    t2: usize,
    v: usize,
) -> bool {
    !on_opposite_sides(sides, t1, t2)
        && mesh.triangles[t1]
            .iter()
            .any(|&u| u != v && mesh.triangles[t2].contains(&u))
}

/// Nodes duplicated by [`split_stars`] and the triangle corners that moved to them.
struct StarSplit {
    /// `(triangle, original node)` -> duplicate used by that triangle.
    remap: HashMap<(usize, usize), usize>,
    duplicated: Vec<(usize, usize)>,
}

/// Splits every interface node into one node per sector of its triangle star (see
/// the module documentation).
fn split_stars(
    mesh: &mut SoftBodyMesh,
    sides: &HashMap<usize, bool>,
    interface: &[usize],
) -> StarSplit {
    let n = mesh.nodes.len();
    let interface_set: HashSet<usize> = interface.iter().copied().collect();

    let mut stars: HashMap<usize, Vec<usize>> = HashMap::new();
    for (t, tri) in mesh.triangles.iter().enumerate() {
        if tri.iter().any(|&i| i >= n) {
            continue;
        }
        for &v in tri {
            if interface_set.contains(&v) {
                let star = stars.entry(v).or_default();
                if star.last() != Some(&t) {
                    star.push(t);
                }
            }
        }
    }

    let mut split = StarSplit {
        remap: HashMap::new(),
        duplicated: Vec::new(),
    };
    for &v in interface {
        let Some(star) = stars.get(&v) else { continue };
        if star.len() < 2 {
            continue;
        }
        // Sectors: connected components of the star under "stay joined".
        let mut parent: Vec<usize> = (0..star.len()).collect();
        for (i, &ti) in star.iter().enumerate() {
            for (j, &tj) in star.iter().enumerate().skip(i + 1) {
                if stay_joined_at(mesh, sides, ti, tj, v) {
                    let (ri, rj) = (find_root(&mut parent, i), find_root(&mut parent, j));
                    if ri != rj {
                        parent[rj.max(ri)] = rj.min(ri);
                    }
                }
            }
        }
        // Group the triangles by sector, in order of first appearance.
        let mut sector_of_root: HashMap<usize, usize> = HashMap::new();
        let mut sectors: Vec<Vec<usize>> = Vec::new();
        for (i, &tri) in star.iter().enumerate() {
            let root = find_root(&mut parent, i);
            let index = *sector_of_root.entry(root).or_insert_with(|| {
                sectors.push(Vec::new());
                sectors.len() - 1
            });
            sectors[index].push(tri);
        }
        if sectors.len() < 2 {
            continue;
        }
        let areas: Vec<f64> = sectors
            .iter()
            .map(|sector| sector.iter().map(|&t| triangle_area(mesh, t)).sum())
            .collect();
        let total_area: f64 = areas.iter().sum();
        let mut keeper = 0;
        for (i, &area) in areas.iter().enumerate() {
            if area > areas[keeper] {
                keeper = i;
            }
        }

        let position = mesh.nodes[v].position;
        let velocity = mesh.nodes[v].velocity;
        let mass = mesh.nodes[v].mass;
        let mut shared_out = 0.0;
        for (i, sector) in sectors.iter().enumerate() {
            if i == keeper {
                continue;
            }
            let share = if total_area > 1e-300 {
                areas[i] / total_area
            } else {
                1.0 / sectors.len() as f64
            };
            let duplicate = mesh.add_node(position, mass * share);
            mesh.nodes[duplicate].velocity = velocity;
            shared_out += mass * share;
            for &t in sector {
                split.remap.insert((t, v), duplicate);
            }
            split.duplicated.push((v, duplicate));
        }
        mesh.nodes[v].mass = mass - shared_out;
    }
    split
}

/// Opens the crack segment `from` -> `to` in `mesh` and reports what changed.
///
/// `surface_normal` is the normal of the surface the crack runs on; the crack
/// cuts along the plane through the segment that contains it.  With `None` the
/// normal of the triangle closest to the segment is used (triangles are assumed to
/// be consistently wound), and meshes without triangles use the coordinate axis
/// least aligned with the crack (z for a crack in the xy plane).  A crack that
/// would run along the surface normal, a zero-length or non-finite segment, and a
/// segment that crosses no edge leave the mesh untouched.
///
/// See the [module documentation](self) for the topological operation.  The
/// positions of all nodes stay unchanged: the new crack faces coincide until the
/// dynamics pull them apart.  Edge indices are only stable as long as no edge
/// without a triangle is cut (those edges are removed).
pub fn open_crack(
    mesh: &mut SoftBodyMesh,
    from: [f64; 3],
    to: [f64; 3],
    surface_normal: Option<[f64; 3]>,
) -> CrackOpening {
    let Some(plane) = CrackPlane::new(mesh, from, to, surface_normal) else {
        return CrackOpening::default();
    };
    let candidates = candidate_edges(mesh);
    let cut: HashSet<(usize, usize)> = candidates
        .iter()
        .copied()
        .filter(|&(a, b)| plane.crosses(mesh.nodes[a].position, mesh.nodes[b].position))
        .collect();
    if cut.is_empty() {
        return CrackOpening::default();
    }

    // Triangle edges and sides before the split.
    let original_triangles = mesh.triangles.clone();
    let edge_triangles = incident_triangles(&original_triangles);
    let sides = cut_triangle_sides(mesh, &plane, &cut);
    let interface = interface_nodes(&edge_triangles, &sides);

    let split = split_stars(mesh, &sides, &interface);
    for (&(t, v), &duplicate) in &split.remap {
        for corner in mesh.triangles[t].iter_mut() {
            if *corner == v {
                *corner = duplicate;
            }
        }
    }

    let image = |t: usize, x: usize| split.remap.get(&(t, x)).copied().unwrap_or(x);
    // Distinct images of the edge (a, b) in the triangles around it, in order.
    let images_of = |a: usize, b: usize| -> Vec<(usize, usize)> {
        let mut images: Vec<(usize, usize)> = Vec::new();
        if let Some(tris) = edge_triangles.get(&edge_key(a, b)) {
            for &t in tris {
                let pair = (image(t, a), image(t, b));
                if !images
                    .iter()
                    .any(|&other| edge_key(other.0, other.1) == edge_key(pair.0, pair.1))
                {
                    images.push(pair);
                }
            }
        }
        images
    };

    // Edge list: triangle edges follow their triangles, loose cut edges vanish.
    let mut kept: Vec<[usize; 2]> = Vec::with_capacity(mesh.edges.len());
    let mut appended: Vec<[usize; 2]> = Vec::new();
    for &[a, b] in &mesh.edges {
        let images = images_of(a, b);
        if images.is_empty() {
            if !cut.contains(&edge_key(a, b)) {
                kept.push([a, b]);
            }
            continue;
        }
        kept.push([images[0].0, images[0].1]);
        for pair in &images[1..] {
            appended.push([pair.0, pair.1]);
        }
    }
    kept.extend(appended);
    mesh.edges = kept;

    let severed_edges = cut
        .iter()
        .filter(|&&(a, b)| {
            let images = images_of(a, b);
            images.is_empty()
                || images.len() > 1
                || edge_key(images[0].0, images[0].1) != edge_key(a, b)
        })
        .count();
    CrackOpening {
        severed_edges,
        duplicated_nodes: split.duplicated,
    }
}

/// Opens the crack segment `from` -> `to` in `mesh` and returns the number of
/// edges whose connectivity it severed.
///
/// Shorthand for [`open_crack`] with the surface normal estimated from the mesh;
/// use [`open_crack`] to learn which nodes were duplicated.
pub fn open_crack_segment(mesh: &mut SoftBodyMesh, from: [f64; 3], to: [f64; 3]) -> usize {
    open_crack(mesh, from, to, None).severed_edges
}
