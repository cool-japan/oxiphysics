// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Static 3-D k-d tree for nearest-neighbour and k-nearest-neighbour queries
//! over boundary samples. Queries may lie anywhere in space (the shrinking-ball
//! iteration probes points far outside the sample's bounding box).

use super::dot3;

/// Implicit k-d tree node.
#[derive(Debug, Clone, Copy)]
struct KdNode {
    /// Index of the point stored at this node.
    point: usize,
    /// Split axis (0, 1, 2).
    axis: usize,
    /// Left child (`usize::MAX` if none).
    left: usize,
    /// Right child (`usize::MAX` if none).
    right: usize,
}

/// A static k-d tree over a borrowed point set.
#[derive(Debug, Clone)]
pub(super) struct KdTree<'a> {
    points: &'a [[f64; 3]],
    nodes: Vec<KdNode>,
    root: usize,
}

const NONE: usize = usize::MAX;

fn sq_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    dot3(d, d)
}

impl<'a> KdTree<'a> {
    /// Build the tree over the subset `indices` of `points`.
    pub(super) fn build(points: &'a [[f64; 3]], indices: &[usize]) -> Self {
        let mut idx: Vec<usize> = indices.to_vec();
        let mut nodes = Vec::with_capacity(idx.len());
        let root = Self::build_rec(points, &mut idx, 0, &mut nodes);
        Self {
            points,
            nodes,
            root,
        }
    }

    fn build_rec(
        points: &[[f64; 3]],
        idx: &mut [usize],
        depth: usize,
        nodes: &mut Vec<KdNode>,
    ) -> usize {
        if idx.is_empty() {
            return NONE;
        }
        // Split on the axis of largest extent for balanced, well-shaped cells.
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &i in idx.iter() {
            for k in 0..3 {
                lo[k] = lo[k].min(points[i][k]);
                hi[k] = hi[k].max(points[i][k]);
            }
        }
        let mut axis = depth % 3;
        let mut best = -1.0;
        for k in 0..3 {
            let ext = hi[k] - lo[k];
            if ext > best {
                best = ext;
                axis = k;
            }
        }
        let mid = idx.len() / 2;
        idx.select_nth_unstable_by(mid, |&a, &b| points[a][axis].total_cmp(&points[b][axis]));
        let point = idx[mid];
        let slot = nodes.len();
        nodes.push(KdNode {
            point,
            axis,
            left: NONE,
            right: NONE,
        });
        let (left_part, rest) = idx.split_at_mut(mid);
        let right_part = &mut rest[1..];
        let left = Self::build_rec(points, left_part, depth + 1, nodes);
        let right = Self::build_rec(points, right_part, depth + 1, nodes);
        nodes[slot].left = left;
        nodes[slot].right = right;
        slot
    }

    /// Nearest stored point to `q` among those accepted by `accept`; returns
    /// `(index, squared distance)`.
    pub(super) fn nearest_filtered<F: Fn(usize) -> bool>(
        &self,
        q: [f64; 3],
        accept: F,
    ) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        self.nearest_rec(self.root, q, &accept, &mut best);
        best
    }

    fn nearest_rec<F: Fn(usize) -> bool>(
        &self,
        node: usize,
        q: [f64; 3],
        accept: &F,
        best: &mut Option<(usize, f64)>,
    ) {
        if node == NONE {
            return;
        }
        let nd = self.nodes[node];
        let p = self.points[nd.point];
        if accept(nd.point) {
            let d = sq_dist(p, q);
            if best.is_none_or(|(_, bd)| d < bd) {
                *best = Some((nd.point, d));
            }
        }
        let diff = q[nd.axis] - p[nd.axis];
        let (near, far) = if diff < 0.0 {
            (nd.left, nd.right)
        } else {
            (nd.right, nd.left)
        };
        self.nearest_rec(near, q, accept, best);
        if best.is_none_or(|(_, bd)| diff * diff < bd) {
            self.nearest_rec(far, q, accept, best);
        }
    }

    /// The `k` nearest stored points to `q` (excluding indices rejected by
    /// `accept`), sorted by increasing distance.
    pub(super) fn k_nearest<F: Fn(usize) -> bool>(
        &self,
        q: [f64; 3],
        k: usize,
        accept: F,
    ) -> Vec<(usize, f64)> {
        let mut found: Vec<(usize, f64)> = Vec::with_capacity(k + 1);
        if k > 0 {
            self.knn_rec(self.root, q, k, &accept, &mut found);
        }
        found
    }

    fn knn_rec<F: Fn(usize) -> bool>(
        &self,
        node: usize,
        q: [f64; 3],
        k: usize,
        accept: &F,
        found: &mut Vec<(usize, f64)>,
    ) {
        if node == NONE {
            return;
        }
        let nd = self.nodes[node];
        let p = self.points[nd.point];
        if accept(nd.point) {
            let d = sq_dist(p, q);
            if found.len() < k || found.last().is_some_and(|&(_, w)| d < w) {
                let pos = found.partition_point(|&(_, w)| w <= d);
                found.insert(pos, (nd.point, d));
                if found.len() > k {
                    found.pop();
                }
            }
        }
        let diff = q[nd.axis] - p[nd.axis];
        let (near, far) = if diff < 0.0 {
            (nd.left, nd.right)
        } else {
            (nd.right, nd.left)
        };
        self.knn_rec(near, q, k, accept, found);
        let worst = if found.len() < k {
            f64::INFINITY
        } else {
            found.last().map_or(f64::INFINITY, |&(_, w)| w)
        };
        if diff * diff < worst {
            self.knn_rec(far, q, k, accept, found);
        }
    }
}
