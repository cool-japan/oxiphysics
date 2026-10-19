// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conforming longest-edge (Rivara) bisection of triangle meshes.
//!
//! [`refine_conforming`] bisects the marked triangles of a mesh by their
//! longest edge and propagates the bisection to the neighbours across the split
//! edge, so the refined mesh never contains hanging nodes.
//!
//! # Algorithm
//!
//! Every triangle's *longest edge* is the maximum of a strict total order on
//! edges: longer squared length first, ties broken by the (smaller) pair of
//! global node indices.  Because the order is global, two triangles that share
//! an edge always agree on whether it is the longest edge of the other one.
//!
//! To bisect a triangle `t` the *longest-edge propagation path* is followed:
//! `t₀ = t`, `t₁` is a neighbour across the longest edge `E₀` of `t₀` whose own
//! longest edge differs from `E₀`, `t₂` a neighbour across the longest edge of
//! `t₁`, and so on.  The edges along the path strictly increase in the order,
//! so the path ends at a *terminal edge* `E`: every triangle carrying `E` has it
//! as its longest edge (or `E` lies on the boundary).  All triangles carrying
//! `E` are bisected at once about one shared midpoint node, which keeps the mesh
//! conforming after every single step.  The path is then followed again from `t`
//! until `t` itself has been bisected.
//!
//! Each bisection replaces the triangle `[n₀, n₁, n₂]` with longest edge
//! `(n_{i0}, n_{i1})` and opposite vertex `n_{i2}` by the children
//! `[n_{i0}, m, n_{i2}]` and `[n_{i1}, n_{i2}, m]`, where `m` is the midpoint
//! node.  `(i0, i1, i2)` is a cyclic permutation of `(0, 1, 2)`, so the
//! orientation of the parent is preserved.

use super::{TriElement, compute_h_element, sub3};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

/// An undirected edge as the ordered pair of its global node indices.
type EdgeKey = (usize, usize);

fn edge_key(a: usize, b: usize) -> EdgeKey {
    if a <= b { (a, b) } else { (b, a) }
}

/// Squared length of edge `k` (vertex `k` to vertex `(k + 1) % 3`).
fn edge_len2(elem: &TriElement, k: usize) -> f64 {
    let d = sub3(elem.coords[(k + 1) % 3], elem.coords[k]);
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

/// Global key of edge `k` of an element.
fn element_edge_key(elem: &TriElement, k: usize) -> EdgeKey {
    edge_key(elem.nodes[k], elem.nodes[(k + 1) % 3])
}

/// `true` if edge `a` of `elem` ranks above edge `b` in the longest-edge order.
fn outranks(elem: &TriElement, a: usize, b: usize) -> bool {
    match edge_len2(elem, a).total_cmp(&edge_len2(elem, b)) {
        Ordering::Greater => true,
        Ordering::Less => false,
        Ordering::Equal => element_edge_key(elem, a) < element_edge_key(elem, b),
    }
}

/// Local index `k` of the longest edge (vertex `k` to vertex `(k + 1) % 3`).
fn longest_edge_index(elem: &TriElement) -> usize {
    let mut best = 0;
    for k in 1..3 {
        if outranks(elem, k, best) {
            best = k;
        }
    }
    best
}

/// A triangle that repeats a node index cannot be bisected consistently.
fn is_degenerate(elem: &TriElement) -> bool {
    let [a, b, c] = elem.nodes;
    a == b || b == c || a == c
}

/// Mesh being refined together with the edge → triangles incidence map.
struct Refiner<'a> {
    nodes: &'a mut Vec<[f64; 3]>,
    elements: &'a mut Vec<TriElement>,
    on_edge: HashMap<EdgeKey, Vec<usize>>,
}

impl<'a> Refiner<'a> {
    fn new(nodes: &'a mut Vec<[f64; 3]>, elements: &'a mut Vec<TriElement>) -> Self {
        let mut refiner = Refiner {
            nodes,
            elements,
            on_edge: HashMap::new(),
        };
        for slot in 0..refiner.elements.len() {
            refiner.register(slot);
        }
        refiner
    }

    /// Record the edges of the triangle in `slot` (degenerate ones are skipped).
    fn register(&mut self, slot: usize) {
        let elem = &self.elements[slot];
        if is_degenerate(elem) {
            return;
        }
        let nodes = elem.nodes;
        for k in 0..3 {
            let list = self
                .on_edge
                .entry(edge_key(nodes[k], nodes[(k + 1) % 3]))
                .or_default();
            if !list.contains(&slot) {
                list.push(slot);
            }
        }
    }

    /// Forget the edges of the triangle currently stored in `slot`.
    fn unregister(&mut self, slot: usize) {
        let nodes = self.elements[slot].nodes;
        for k in 0..3 {
            let key = edge_key(nodes[k], nodes[(k + 1) % 3]);
            let now_empty = match self.on_edge.get_mut(&key) {
                Some(list) => {
                    list.retain(|&s| s != slot);
                    list.is_empty()
                }
                None => false,
            };
            if now_empty {
                self.on_edge.remove(&key);
            }
        }
    }

    fn longest_edge_key(&self, slot: usize) -> EdgeKey {
        let elem = &self.elements[slot];
        element_edge_key(elem, longest_edge_index(elem))
    }

    /// Follow the longest-edge propagation path from `start` to its terminal
    /// edge.  The path is acyclic; the bound only guards against corrupted
    /// input and still yields an edge that can be bisected conformingly.
    fn terminal_edge(&self, start: usize) -> EdgeKey {
        let mut current = start;
        let mut edge = self.longest_edge_key(current);
        for _ in 0..=self.elements.len() {
            let next = self.on_edge.get(&edge).and_then(|list| {
                list.iter()
                    .copied()
                    .find(|&s| s != current && self.longest_edge_key(s) != edge)
            });
            match next {
                Some(neighbour) => {
                    current = neighbour;
                    edge = self.longest_edge_key(current);
                }
                None => break,
            }
        }
        edge
    }

    /// Bisect every triangle that carries `edge` about one new midpoint node.
    /// Returns `false` if nothing could be split.
    fn split_edge(&mut self, edge: EdgeKey, bisected: &mut HashSet<usize>) -> bool {
        let slots: Vec<usize> = self.on_edge.get(&edge).cloned().unwrap_or_default();
        let Some(&first) = slots.first() else {
            return false;
        };
        let position = |elem: &TriElement, id: usize| {
            elem.nodes
                .iter()
                .position(|&n| n == id)
                .map(|k| elem.coords[k])
        };
        let (Some(pa), Some(pb)) = (
            position(&self.elements[first], edge.0),
            position(&self.elements[first], edge.1),
        ) else {
            return false;
        };
        let mid = [
            0.5 * (pa[0] + pb[0]),
            0.5 * (pa[1] + pb[1]),
            0.5 * (pa[2] + pb[2]),
        ];
        let mid_idx = self.nodes.len();
        self.nodes.push(mid);

        for slot in slots {
            self.unregister(slot);
            let parent = self.elements[slot].clone();
            let Some(k) = (0..3).find(|&k| element_edge_key(&parent, k) == edge) else {
                self.register(slot);
                continue;
            };
            let (i0, i1, i2) = (k, (k + 1) % 3, (k + 2) % 3);
            let coords_a = [parent.coords[i0], mid, parent.coords[i2]];
            let coords_b = [parent.coords[i1], parent.coords[i2], mid];
            self.elements[slot] = TriElement {
                nodes: [parent.nodes[i0], mid_idx, parent.nodes[i2]],
                coords: coords_a,
                degree: parent.degree,
                h_e: compute_h_element(&coords_a),
            };
            self.elements.push(TriElement {
                nodes: [parent.nodes[i1], parent.nodes[i2], mid_idx],
                coords: coords_b,
                degree: parent.degree,
                h_e: compute_h_element(&coords_b),
            });
            let new_slot = self.elements.len() - 1;
            self.register(slot);
            self.register(new_slot);
            bisected.insert(slot);
        }
        true
    }

    /// Make sure the triangle that sits in `start` is bisected (possibly as a
    /// side effect of bisecting another triangle's propagation path).
    fn bisect_element(&mut self, start: usize, bisected: &mut HashSet<usize>) {
        if is_degenerate(&self.elements[start]) {
            return;
        }
        while !bisected.contains(&start) {
            let edge = self.terminal_edge(start);
            if !self.split_edge(edge, bisected) {
                break;
            }
        }
    }
}

/// Bisect the triangles flagged in `marked` (one flag per element) by their
/// longest edge, appending the new midpoint nodes to `nodes` and keeping the
/// mesh conforming.
///
/// The parent triangle's slot receives its first child and the second child is
/// appended to `elements`; triangles that must be bisected only to close a
/// hanging node are refined as well, so more elements than the marked ones may
/// change.  Triangles that repeat a node index are left untouched.
pub(super) fn refine_conforming(
    nodes: &mut Vec<[f64; 3]>,
    elements: &mut Vec<TriElement>,
    marked: &[bool],
) {
    let mut refiner = Refiner::new(nodes, elements);
    let mut bisected: HashSet<usize> = HashSet::new();
    for (slot, &is_marked) in marked.iter().enumerate() {
        if is_marked {
            refiner.bisect_element(slot, &mut bisected);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::MeshRefinement;
    use super::*;

    /// Mesh of the unit square: `cells × cells` squares, each cut along the
    /// diagonal from its lower-left to its upper-right corner.  The vertex
    /// order of every triangle is rotated by `(i + j) % 3` so the longest edge
    /// is not always local edge 0.
    fn square_mesh(cells: usize) -> MeshRefinement {
        let mut nodes = Vec::new();
        for j in 0..=cells {
            for i in 0..=cells {
                // Division keeps the boundary coordinates exactly 0 and 1.
                nodes.push([i as f64 / cells as f64, j as f64 / cells as f64, 0.0]);
            }
        }
        let id = |i: usize, j: usize| j * (cells + 1) + i;
        let mut elements = Vec::new();
        for j in 0..cells {
            for i in 0..cells {
                let tris = [
                    [id(i, j), id(i + 1, j), id(i + 1, j + 1)],
                    [id(i, j), id(i + 1, j + 1), id(i, j + 1)],
                ];
                for tri in tris {
                    let r = (i + j) % 3;
                    let rotated = [tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]];
                    let coords = [nodes[rotated[0]], nodes[rotated[1]], nodes[rotated[2]]];
                    elements.push(TriElement::new(rotated, coords));
                }
            }
        }
        MeshRefinement::new(elements, nodes)
    }

    fn signed_area(e: &TriElement) -> f64 {
        let [p0, p1, p2] = e.coords;
        0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]))
    }

    /// Distance from point `p` to the segment `a-b`.
    fn dist_to_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
        let ab = sub3(b, a);
        let len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        let ap = sub3(p, a);
        let t = ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / len2).clamp(0.0, 1.0);
        let d = [
            p[0] - a[0] - t * ab[0],
            p[1] - a[1] - t * ab[1],
            p[2] - a[2],
        ];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    /// Check every invariant of a refined unit-square mesh.
    fn assert_valid_conforming_mesh(mesh: &MeshRefinement) {
        let mut area = 0.0;
        for (e, elem) in mesh.elements.iter().enumerate() {
            for k in 0..3 {
                assert_eq!(
                    mesh.nodes[elem.nodes[k]], elem.coords[k],
                    "element {e}: coords[{k}] disagree with node {}",
                    elem.nodes[k]
                );
            }
            assert!(signed_area(elem) > 0.0, "element {e} lost its orientation");
            assert_eq!(elem.h_e, compute_h_element(&elem.coords));
            area += elem.area();
        }
        assert!((area - 1.0).abs() < 1e-12, "total area {area}");

        let mut edge_count: HashMap<EdgeKey, usize> = HashMap::new();
        for elem in &mesh.elements {
            for k in 0..3 {
                *edge_count.entry(element_edge_key(elem, k)).or_insert(0) += 1;
            }
        }
        for (&(a, b), &count) in &edge_count {
            assert!(count <= 2, "edge ({a},{b}) shared by {count} triangles");
            if count == 1 {
                // A boundary edge must lie on one side of the unit square.
                let (pa, pb) = (mesh.nodes[a], mesh.nodes[b]);
                let on_side = (0..2).any(|c| [0.0, 1.0].iter().any(|&s| pa[c] == s && pb[c] == s));
                assert!(on_side, "interior edge ({a},{b}) has only one triangle");
            }
        }
        // No node may sit strictly inside an edge (hanging node).
        for (n, &p) in mesh.nodes.iter().enumerate() {
            for &(a, b) in edge_count.keys() {
                if n == a || n == b {
                    continue;
                }
                let d = dist_to_segment(p, mesh.nodes[a], mesh.nodes[b]);
                assert!(d > 1e-12, "node {n} hangs on edge ({a},{b})");
            }
        }
    }

    #[test]
    fn shared_edge_gets_one_midpoint_for_both_neighbours() {
        // A unit square cut into two triangles; marking one bisects both
        // (they share their longest edge), with a single new node.
        let mut mesh = square_mesh(1);
        assert_eq!((mesh.nodes.len(), mesh.elements.len()), (4, 2));
        mesh.refine_marked(&[true, false]);
        assert_eq!(mesh.nodes.len(), 5, "the diagonal midpoint is shared");
        assert_eq!(mesh.elements.len(), 4);
        assert_eq!(mesh.nodes[4], [0.5, 0.5, 0.0]);
        assert_valid_conforming_mesh(&mesh);
        // All four children meet at the new node.
        assert!(mesh.elements.iter().all(|e| e.nodes.contains(&4)));
    }

    #[test]
    fn children_use_the_bisected_longest_edge_for_every_rotation() {
        // For each cyclic rotation the longest edge sits at a different local
        // position; children must still be consistent with the node array.
        for rotation in 0..3 {
            let nodes = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
            let tri = [0usize, 1, 2];
            let rotated = [
                tri[rotation],
                tri[(rotation + 1) % 3],
                tri[(rotation + 2) % 3],
            ];
            let coords = [nodes[rotated[0]], nodes[rotated[1]], nodes[rotated[2]]];
            let elem = TriElement::new(rotated, coords);
            let mut mesh = MeshRefinement::new(vec![elem], nodes);
            mesh.refine_marked(&[true]);
            assert_eq!(mesh.elements.len(), 2);
            assert_eq!(mesh.nodes.len(), 4);
            // The longest edge is (0,0)-(3,0); its midpoint is (1.5, 0).
            assert_eq!(mesh.nodes[3], [1.5, 0.0, 0.0]);
            let area: f64 = mesh.elements.iter().map(|e| e.area()).sum();
            assert!((area - 1.5).abs() < 1e-14);
            for elem in &mesh.elements {
                assert!(elem.nodes.contains(&3));
                for k in 0..3 {
                    assert_eq!(mesh.nodes[elem.nodes[k]], elem.coords[k]);
                }
                assert!(signed_area(elem) > 0.0);
            }
        }
    }

    #[test]
    fn refinement_propagates_to_neighbours_whose_longest_edge_is_elsewhere() {
        // T (marked) has the long base (0,0)-(2,0) as its longest edge.  The
        // neighbour T' below the base owns a longer side, so it must be bisected
        // first; its neighbour T'' shares that longer side.
        let nodes = vec![
            [0.0, 0.0, 0.0],   // 0
            [2.0, 0.0, 0.0],   // 1
            [1.0, 1.0, 0.0],   // 2   apex of T
            [1.0, -3.0, 0.0],  // 3   apex of T'
            [-2.0, -3.0, 0.0], // 4  apex of T''
        ];
        let tri = |n: [usize; 3]| TriElement::new(n, [nodes[n[0]], nodes[n[1]], nodes[n[2]]]);
        let elements = vec![
            tri([0, 1, 2]), // T
            tri([1, 0, 3]), // T'  (longest edge 0-3, length sqrt(10))
            tri([0, 4, 3]), // T'' shares 0-3
        ];
        let total: f64 = elements.iter().map(|e| e.area()).sum();
        let mut mesh = MeshRefinement::new(elements, nodes);
        mesh.refine_marked(&[true, false, false]);
        assert!(mesh.elements.len() > 4, "closure must refine beyond T");
        let after: f64 = mesh.elements.iter().map(|e| e.area()).sum();
        assert!((after - total).abs() < 1e-12);
        for (e, elem) in mesh.elements.iter().enumerate() {
            for k in 0..3 {
                assert_eq!(mesh.nodes[elem.nodes[k]], elem.coords[k], "element {e}");
            }
            assert!(signed_area(elem) > 0.0, "element {e} flipped");
        }
        // Conformity: every edge is shared by at most two triangles and no node
        // hangs on an edge.
        let mut count: HashMap<EdgeKey, usize> = HashMap::new();
        for elem in &mesh.elements {
            for k in 0..3 {
                *count.entry(element_edge_key(elem, k)).or_insert(0) += 1;
            }
        }
        assert!(count.values().all(|&c| c <= 2));
        for (n, &p) in mesh.nodes.iter().enumerate() {
            for &(a, b) in count.keys() {
                if n != a && n != b {
                    assert!(
                        dist_to_segment(p, mesh.nodes[a], mesh.nodes[b]) > 1e-12,
                        "node {n} hangs on edge ({a},{b})"
                    );
                }
            }
        }
        // The shared base (0,0)-(2,0) was split, i.e. T no longer exists whole.
        assert!(
            !mesh.elements.iter().any(|e| {
                let mut ids = e.nodes;
                ids.sort_unstable();
                ids == [0, 1, 2]
            }),
            "the marked triangle must be bisected"
        );
    }

    #[test]
    fn repeated_random_marking_keeps_the_mesh_conforming() {
        let mut mesh = square_mesh(3);
        assert_valid_conforming_mesh(&mesh);
        let mut state: u64 = 0xDEAD_BEEF_1234_5678;
        for round in 0..6 {
            let marked: Vec<bool> = (0..mesh.elements.len())
                .map(|_| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    (state >> 60) < 4
                })
                .collect();
            let before = mesh.elements.len();
            let n_marked = marked.iter().filter(|&&m| m).count();
            mesh.refine_marked(&marked);
            assert!(
                mesh.elements.len() >= before + n_marked.min(1),
                "round {round}: nothing was refined"
            );
            assert_valid_conforming_mesh(&mesh);
        }
    }

    #[test]
    fn refining_every_triangle_halves_every_area() {
        let mut mesh = square_mesh(2);
        let before = mesh.elements.len();
        mesh.refine_marked(&vec![true; before]);
        // Each pair across an interior edge is bisected together, and every
        // triangle is bisected exactly once.
        assert_eq!(mesh.elements.len(), 2 * before);
        assert_valid_conforming_mesh(&mesh);
        let max_area = mesh
            .elements
            .iter()
            .map(|e| e.area())
            .fold(0.0_f64, f64::max);
        assert!((max_area - 1.0 / 16.0).abs() < 1e-14, "max area {max_area}");
    }

    #[test]
    fn unmarked_mesh_and_degenerate_triangles_are_left_alone() {
        let mut mesh = square_mesh(2);
        let snapshot: Vec<([usize; 3], [[f64; 3]; 3])> =
            mesh.elements.iter().map(|e| (e.nodes, e.coords)).collect();
        let n_nodes = mesh.nodes.len();
        mesh.refine_marked(&vec![false; mesh.elements.len()]);
        assert_eq!(mesh.nodes.len(), n_nodes);
        let after: Vec<([usize; 3], [[f64; 3]; 3])> =
            mesh.elements.iter().map(|e| (e.nodes, e.coords)).collect();
        assert_eq!(snapshot, after);

        // A triangle that repeats a node cannot be bisected.
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let degenerate = TriElement::new([0, 1, 1], [nodes[0], nodes[1], nodes[1]]);
        let mut mesh = MeshRefinement::new(vec![degenerate], nodes);
        mesh.refine_marked(&[true]);
        assert_eq!((mesh.elements.len(), mesh.nodes.len()), (1, 2));
    }

    #[test]
    fn coincident_triangles_are_bisected_together() {
        // Four triangles on the same three nodes (as the adaptive driver tests
        // build) share every edge; they must be split together about one node.
        let nodes = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let coords = [nodes[0], nodes[1], nodes[2]];
        let elements: Vec<TriElement> =
            (0..4).map(|_| TriElement::new([0, 1, 2], coords)).collect();
        let mut mesh = MeshRefinement::new(elements, nodes);
        mesh.refine_marked(&[true, false, false, false]);
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.elements.len(), 8);
        assert!(mesh.elements.iter().all(|e| e.nodes.contains(&3)));
    }
}
