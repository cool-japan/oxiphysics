// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conversion of a weighted graph of medial balls into a curve skeleton
//! ([`MedialAxis`] with branches and junctions).
//!
//! The input is a set of medial balls (centre + radius) together with an
//! adjacency graph `G` that only connects balls lying in the interior of the
//! shape (Voronoi edges between inner Voronoi vertices, or surface-neighbour
//! relations between shrinking balls). Some balls are *anchors* (genuine
//! medial balls such as Voronoi poles); the others only provide connectivity.
//!
//! Pipeline:
//!
//! 1. **merge** balls joined by an edge shorter than `merge_tolerance`
//!    (co-spherical input yields many tets with the same circumcentre);
//! 2. **minimum spanning forest** (Kruskal) of the merged graph;
//! 3. **Steiner pruning**: repeatedly remove non-anchor leaves, leaving the
//!    smallest subtree that spans every anchor;
//! 4. **scale-axis twig pruning** (after Giesen et al., "The Scale Axis
//!    Transform", 2009): a leaf chain ending at a junction is removed when every
//!    ball on it is contained in some ball of the remaining skeleton scaled by
//!    `1 + prune_scale` — such a twig carries no shape information;
//! 5. **loop closure**: a graph edge whose endpoints hang off skeleton nodes
//!    `A`, `B` closes a genuine tunnel when the skeleton detour from `A` to `B`
//!    exceeds `|A − B| + 2π · loop_factor · max(r_A, r_B)`; such edges are
//!    added back, restoring loops (e.g. the circle of a torus) that the
//!    spanning tree cut, and the twig pruning is repeated on the closed graph;
//! 6. **branch extraction**: maximal chains between nodes of degree ≠ 2;
//!    nodes of degree ≥ 3 are junctions.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

use super::{MedialAxis, MedialAxisPoint, dist3, sub3};

/// One medial ball fed into the skeleton graph.
#[derive(Debug, Clone, Copy)]
pub(super) struct SkelNode {
    /// Ball centre.
    pub(super) position: [f64; 3],
    /// Ball radius.
    pub(super) radius: f64,
    /// Whether the ball must be kept by the Steiner pruning.
    pub(super) anchor: bool,
}

/// Tunable parameters of [`skeletonize`].
#[derive(Debug, Clone, Copy)]
pub(super) struct SkeletonParams {
    /// Edges shorter than this collapse their endpoints into one ball.
    pub(super) merge_tolerance: f64,
    /// Scale-axis factor `s − 1` used by the twig pruning.
    pub(super) prune_scale: f64,
    /// Multiplier of `2π r` for the loop-closure significance test.
    pub(super) loop_factor: f64,
}

/// Disjoint-set forest with path halving and union by size.
struct UnionFind {
    parent: Vec<usize>,
    size: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) -> bool {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return false;
        }
        let (big, small) = if self.size[ra] >= self.size[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[small] = big;
        self.size[big] += self.size[small];
        true
    }
}

/// Min-heap entry for the bounded Dijkstra.
#[derive(Debug, Clone, Copy, PartialEq)]
struct HeapItem {
    dist: f64,
    node: usize,
}

impl Eq for HeapItem {}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .dist
            .total_cmp(&self.dist)
            .then_with(|| other.node.cmp(&self.node))
    }
}

/// Merged graph of balls: nodes plus weighted, deduplicated edges.
struct MergedGraph {
    nodes: Vec<SkelNode>,
    edges: Vec<(usize, usize, f64)>,
}

/// Collapse balls joined by very short edges and deduplicate the edges.
fn merge_nodes(nodes: &[SkelNode], edges: &[(usize, usize)], tolerance: f64) -> MergedGraph {
    let n = nodes.len();
    let mut uf = UnionFind::new(n);
    for &(a, b) in edges {
        if a < n && b < n && dist3(nodes[a].position, nodes[b].position) <= tolerance {
            uf.union(a, b);
        }
    }
    let mut cluster_of = vec![usize::MAX; n];
    let mut sums: Vec<([f64; 3], f64, usize, bool)> = Vec::new();
    for i in 0..n {
        let root = uf.find(i);
        if cluster_of[root] == usize::MAX {
            cluster_of[root] = sums.len();
            sums.push(([0.0; 3], 0.0, 0, false));
        }
        let c = cluster_of[root];
        cluster_of[i] = c;
        let entry = &mut sums[c];
        for k in 0..3 {
            entry.0[k] += nodes[i].position[k];
        }
        entry.1 += nodes[i].radius;
        entry.2 += 1;
        entry.3 |= nodes[i].anchor;
    }
    let merged: Vec<SkelNode> = sums
        .iter()
        .map(|&(p, r, count, anchor)| {
            let inv = 1.0 / count as f64;
            SkelNode {
                position: [p[0] * inv, p[1] * inv, p[2] * inv],
                radius: r * inv,
                anchor,
            }
        })
        .collect();
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    let mut merged_edges = Vec::new();
    for &(a, b) in edges {
        if a >= n || b >= n {
            continue;
        }
        let (ca, cb) = (cluster_of[a], cluster_of[b]);
        if ca == cb {
            continue;
        }
        let key = (ca.min(cb), ca.max(cb));
        if seen.insert(key) {
            let w = dist3(merged[ca].position, merged[cb].position);
            merged_edges.push((key.0, key.1, w));
        }
    }
    MergedGraph {
        nodes: merged,
        edges: merged_edges,
    }
}

/// Kruskal minimum spanning forest; returns the forest adjacency lists.
fn spanning_forest(n: usize, edges: &[(usize, usize, f64)]) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..edges.len()).collect();
    order.sort_by(|&x, &y| {
        edges[x]
            .2
            .total_cmp(&edges[y].2)
            .then_with(|| (edges[x].0, edges[x].1).cmp(&(edges[y].0, edges[y].1)))
    });
    let mut uf = UnionFind::new(n);
    let mut adj = vec![Vec::new(); n];
    for e in order {
        let (a, b, _) = edges[e];
        if uf.union(a, b) {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    adj
}

/// Remove non-anchor leaves until only the anchor-spanning subtrees remain.
fn steiner_prune(adj: &[Vec<usize>], nodes: &[SkelNode]) -> Vec<bool> {
    let n = nodes.len();
    let mut alive = vec![true; n];
    let mut degree: Vec<usize> = adj.iter().map(Vec::len).collect();
    let mut queue: VecDeque<usize> = (0..n)
        .filter(|&i| !nodes[i].anchor && degree[i] <= 1)
        .collect();
    while let Some(i) = queue.pop_front() {
        if !alive[i] {
            continue;
        }
        alive[i] = false;
        for &j in &adj[i] {
            if alive[j] {
                degree[j] = degree[j].saturating_sub(1);
                if !nodes[j].anchor && degree[j] <= 1 {
                    queue.push_back(j);
                }
            }
        }
    }
    alive
}

/// Live degree of node `i` in `adj` restricted to `alive` nodes.
fn live_degree(adj: &[Vec<usize>], alive: &[bool], i: usize) -> usize {
    adj[i].iter().filter(|&&j| alive[j]).count()
}

/// A leaf chain of the skeleton ending next to a junction.
struct Twig {
    nodes: Vec<usize>,
    junction: usize,
    length: f64,
}

/// Collect every twig: walk from each leaf until a node of degree ≥ 3.
/// Chains that end at another leaf (path components) are not twigs.
fn collect_twigs(adj: &[Vec<usize>], alive: &[bool], nodes: &[SkelNode]) -> Vec<Twig> {
    let mut twigs = Vec::new();
    for leaf in 0..adj.len() {
        if !alive[leaf] || live_degree(adj, alive, leaf) != 1 {
            continue;
        }
        let mut chain = vec![leaf];
        let mut prev = leaf;
        let mut cur = leaf;
        let mut length = 0.0;
        let junction = loop {
            let next = adj[cur].iter().copied().find(|&j| alive[j] && j != prev);
            let Some(next) = next else { break None };
            length += dist3(nodes[cur].position, nodes[next].position);
            let deg = live_degree(adj, alive, next);
            if deg >= 3 {
                break Some(next);
            }
            if deg == 1 {
                break None;
            }
            chain.push(next);
            prev = cur;
            cur = next;
        };
        if let Some(junction) = junction {
            twigs.push(Twig {
                nodes: chain,
                junction,
                length,
            });
        }
    }
    twigs
}

/// Scale-axis twig pruning (see the module docs).
fn prune_twigs(adj: &[Vec<usize>], alive: &mut [bool], nodes: &[SkelNode], scale: f64) {
    let factor = 1.0 + scale.max(0.0);
    loop {
        let mut twigs = collect_twigs(adj, alive, nodes);
        twigs.sort_by(|a, b| a.length.total_cmp(&b.length));
        let mut pruned_any = false;
        let mut on_twig = vec![false; nodes.len()];
        for twig in &twigs {
            // Re-validate: earlier removals in this round may have turned the
            // junction into a regular chain node.
            if !alive[twig.junction] || live_degree(adj, alive, twig.junction) < 3 {
                continue;
            }
            if twig.nodes.iter().any(|&i| !alive[i]) {
                continue;
            }
            for &i in &twig.nodes {
                on_twig[i] = true;
            }
            let covered = twig.nodes.iter().all(|&x| {
                let bx = nodes[x];
                (0..nodes.len()).any(|y| {
                    alive[y]
                        && !on_twig[y]
                        && dist3(bx.position, nodes[y].position) + bx.radius
                            <= factor * nodes[y].radius
                })
            });
            for &i in &twig.nodes {
                on_twig[i] = false;
            }
            if covered {
                for &i in &twig.nodes {
                    alive[i] = false;
                }
                pruned_any = true;
            }
        }
        if !pruned_any {
            break;
        }
    }
}

/// Multi-source BFS over the forest: nearest surviving skeleton node of every
/// node (`None` for nodes whose tree contains no survivor).
fn attachments(adj: &[Vec<usize>], alive: &[bool]) -> Vec<Option<usize>> {
    let n = adj.len();
    let mut att: Vec<Option<usize>> = (0..n).map(|i| alive[i].then_some(i)).collect();
    let mut queue: VecDeque<usize> = (0..n).filter(|&i| alive[i]).collect();
    while let Some(i) = queue.pop_front() {
        for &j in &adj[i] {
            if att[j].is_none() {
                att[j] = att[i];
                queue.push_back(j);
            }
        }
    }
    att
}

/// Dijkstra from `src` in `skel` that gives up beyond `bound`; returns whether
/// `dst` is reachable within `bound`.
fn within_distance(
    skel: &[Vec<usize>],
    nodes: &[SkelNode],
    src: usize,
    dst: usize,
    bound: f64,
) -> bool {
    let mut best: HashMap<usize, f64> = HashMap::new();
    let mut heap = BinaryHeap::new();
    best.insert(src, 0.0);
    heap.push(HeapItem {
        dist: 0.0,
        node: src,
    });
    while let Some(HeapItem { dist, node }) = heap.pop() {
        if node == dst {
            return true;
        }
        if best.get(&node).is_some_and(|&d| dist > d) {
            continue;
        }
        for &next in &skel[node] {
            let nd = dist + dist3(nodes[node].position, nodes[next].position);
            if nd > bound {
                continue;
            }
            if best.get(&next).is_none_or(|&d| nd < d) {
                best.insert(next, nd);
                heap.push(HeapItem {
                    dist: nd,
                    node: next,
                });
            }
        }
    }
    false
}

/// Add back graph edges that close significant loops (see the module docs).
fn close_loops(
    graph: &MergedGraph,
    tree: &[Vec<usize>],
    alive: &[bool],
    loop_factor: f64,
) -> Vec<Vec<usize>> {
    let nodes = &graph.nodes;
    let n = nodes.len();
    let mut skel: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            if alive[i] {
                tree[i].iter().copied().filter(|&j| alive[j]).collect()
            } else {
                Vec::new()
            }
        })
        .collect();
    let att = attachments(tree, alive);
    let mut candidates: Vec<(usize, usize, f64)> = graph
        .edges
        .iter()
        .filter_map(|&(a, b, w)| match (att[a], att[b]) {
            (Some(x), Some(y)) if x != y => Some((x.min(y), x.max(y), w)),
            _ => None,
        })
        .collect();
    candidates.sort_by(|x, y| {
        x.2.total_cmp(&y.2)
            .then_with(|| (x.0, x.1).cmp(&(y.0, y.1)))
    });
    let mut tried: HashSet<(usize, usize)> = HashSet::new();
    for (a, b, _) in candidates {
        if !tried.insert((a, b)) || skel[a].contains(&b) {
            continue;
        }
        let gap = dist3(nodes[a].position, nodes[b].position);
        let bound =
            gap + 2.0 * std::f64::consts::PI * loop_factor * nodes[a].radius.max(nodes[b].radius);
        if !within_distance(&skel, nodes, a, b, bound) {
            skel[a].push(b);
            skel[b].push(a);
        }
    }
    skel
}

/// Split the skeleton graph into branches (maximal chains between nodes of
/// degree ≠ 2) and junctions (degree ≥ 3). Pure cycles become one closed
/// branch whose first and last index coincide; isolated nodes become
/// single-node branches.
fn extract_branches(skel: &[Vec<usize>], alive: &[bool]) -> (Vec<Vec<usize>>, Vec<usize>) {
    let n = skel.len();
    let degree: Vec<usize> = skel.iter().map(Vec::len).collect();
    let mut visited_edge: HashSet<(usize, usize)> = HashSet::new();
    let mut branches = Vec::new();
    let mut on_branch = vec![false; n];
    let walk = |start: usize, first: usize, visited_edge: &mut HashSet<(usize, usize)>| {
        let mut branch = vec![start];
        let mut prev = start;
        let mut cur = first;
        visited_edge.insert((start.min(first), start.max(first)));
        loop {
            branch.push(cur);
            if degree[cur] != 2 || cur == start {
                break;
            }
            let next = skel[cur]
                .iter()
                .copied()
                .find(|&j| j != prev && !visited_edge.contains(&(cur.min(j), cur.max(j))));
            let Some(next) = next else { break };
            visited_edge.insert((cur.min(next), cur.max(next)));
            prev = cur;
            cur = next;
        }
        branch
    };
    for i in 0..n {
        if !alive[i] || degree[i] == 2 {
            continue;
        }
        if degree[i] == 0 {
            branches.push(vec![i]);
            on_branch[i] = true;
            continue;
        }
        for &j in &skel[i] {
            if visited_edge.contains(&(i.min(j), i.max(j))) {
                continue;
            }
            let b = walk(i, j, &mut visited_edge);
            for &k in &b {
                on_branch[k] = true;
            }
            branches.push(b);
        }
    }
    // Remaining degree-2 nodes belong to pure cycles.
    for i in 0..n {
        if !alive[i] || on_branch[i] || degree[i] != 2 {
            continue;
        }
        let first = skel[i][0];
        let b = walk(i, first, &mut visited_edge);
        for &k in &b {
            on_branch[k] = true;
        }
        branches.push(b);
    }
    let junctions = (0..n).filter(|&i| alive[i] && degree[i] >= 3).collect();
    (branches, junctions)
}

/// Farthest-point sampling of `positions` with covering radius `radius`.
///
/// Seeds are picked greedily (first point, then always the point farthest from
/// all seeds) until every point lies within `radius` of a seed; each point is
/// assigned to its nearest seed. Returns the seed index of every point and the
/// number of seeds.
pub(super) fn farthest_point_clusters(positions: &[[f64; 3]], radius: f64) -> (Vec<usize>, usize) {
    let n = positions.len();
    let mut assign = vec![0usize; n];
    if n == 0 {
        return (assign, 0);
    }
    let mut dist = vec![f64::INFINITY; n];
    let mut seeds = 0usize;
    let mut next = 0usize;
    loop {
        let seed_pos = positions[next];
        for i in 0..n {
            let d = dist3(positions[i], seed_pos);
            if d < dist[i] {
                dist[i] = d;
                assign[i] = seeds;
            }
        }
        seeds += 1;
        let (far, far_d) =
            dist.iter().enumerate().fold(
                (0usize, -1.0_f64),
                |acc, (i, &d)| if d > acc.1 { (i, d) } else { acc },
            );
        if far_d <= radius {
            break;
        }
        next = far;
    }
    (assign, seeds)
}

/// Collapse balls into farthest-point clusters of covering radius `radius`:
/// each cluster becomes one ball at the centroid of its members' centres with
/// their mean radius, and clusters are adjacent when any of their members are.
pub(super) fn cluster_balls(
    nodes: &[SkelNode],
    edges: &[(usize, usize)],
    radius: f64,
) -> (Vec<SkelNode>, Vec<(usize, usize)>) {
    let positions: Vec<[f64; 3]> = nodes.iter().map(|b| b.position).collect();
    let (assign, m) = farthest_point_clusters(&positions, radius);
    let mut sum = vec![[0.0_f64; 3]; m];
    let mut rad = vec![0.0_f64; m];
    let mut count = vec![0usize; m];
    let mut anchor = vec![false; m];
    for (i, b) in nodes.iter().enumerate() {
        let c = assign[i];
        for (acc, x) in sum[c].iter_mut().zip(b.position) {
            *acc += x;
        }
        rad[c] += b.radius;
        count[c] += 1;
        anchor[c] |= b.anchor;
    }
    let clustered: Vec<SkelNode> = (0..m)
        .map(|c| {
            let inv = 1.0 / count[c].max(1) as f64;
            SkelNode {
                position: [sum[c][0] * inv, sum[c][1] * inv, sum[c][2] * inv],
                radius: rad[c] * inv,
                anchor: anchor[c],
            }
        })
        .collect();
    let mut cluster_edges: Vec<(usize, usize)> = edges
        .iter()
        .filter(|&&(a, b)| a < nodes.len() && b < nodes.len())
        .map(|&(a, b)| (assign[a], assign[b]))
        .filter(|&(a, b)| a != b)
        .map(|(a, b)| (a.min(b), a.max(b)))
        .collect();
    cluster_edges.sort_unstable();
    cluster_edges.dedup();
    (clustered, cluster_edges)
}

/// Run the whole pipeline and assemble a compact [`MedialAxis`].
pub(super) fn skeletonize(
    nodes: &[SkelNode],
    edges: &[(usize, usize)],
    params: SkeletonParams,
) -> MedialAxis {
    let graph = merge_nodes(nodes, edges, params.merge_tolerance);
    let n = graph.nodes.len();
    let tree = spanning_forest(n, &graph.edges);
    let mut alive = steiner_prune(&tree, &graph.nodes);
    prune_twigs(&tree, &mut alive, &graph.nodes, params.prune_scale);
    let skel = close_loops(&graph, &tree, &alive, params.loop_factor);
    // Closing a loop can turn the stub beyond the closing node into a twig.
    prune_twigs(&skel, &mut alive, &graph.nodes, params.prune_scale);
    let skel: Vec<Vec<usize>> = skel
        .iter()
        .enumerate()
        .map(|(i, nb)| {
            if alive[i] {
                nb.iter().copied().filter(|&j| alive[j]).collect()
            } else {
                Vec::new()
            }
        })
        .collect();
    let (branches, junctions) = extract_branches(&skel, &alive);

    // Compact re-indexing of the surviving nodes.
    let mut new_index = vec![usize::MAX; n];
    let mut points = Vec::new();
    for i in 0..n {
        if alive[i] {
            new_index[i] = points.len();
            points.push(MedialAxisPoint::new(
                graph.nodes[i].position,
                graph.nodes[i].radius,
                [0.0; 3],
                graph.nodes[i].radius,
            ));
        }
    }
    let branches: Vec<Vec<usize>> = branches
        .into_iter()
        .map(|b| b.into_iter().map(|i| new_index[i]).collect())
        .collect();
    let junctions: Vec<usize> = junctions.into_iter().map(|i| new_index[i]).collect();

    // Tangents from the first branch that traverses each point, as the chord
    // between the branch points half a medial radius of arc length before and
    // after it: the local direction at the scale of the shape rather than of
    // the sampling noise.
    let mut has_tangent = vec![false; points.len()];
    for branch in &branches {
        if branch.len() < 2 {
            continue;
        }
        let mut arc = vec![0.0_f64; branch.len()];
        for k in 1..branch.len() {
            arc[k] = arc[k - 1] + dist3(points[branch[k - 1]].position, points[branch[k]].position);
        }
        for k in 0..branch.len() {
            let idx = branch[k];
            if has_tangent[idx] {
                continue;
            }
            let window = 0.5 * points[idx].radius;
            let mut back = k;
            while back > 0 && arc[k] - arc[back] < window {
                back -= 1;
            }
            let mut fwd = k;
            while fwd + 1 < branch.len() && arc[fwd] - arc[k] < window {
                fwd += 1;
            }
            let d = sub3(points[branch[fwd]].position, points[branch[back]].position);
            if let Some(t) = super::poles::unit_or_none(d) {
                points[idx].tangent = t;
                has_tangent[idx] = true;
            }
        }
    }

    MedialAxis {
        points,
        branches,
        junctions,
    }
}
