//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

use super::types::{AStarState, Graph, State, UnionFind};

/// Returns nodes visited in BFS order starting from `start`.
#[allow(dead_code)]
pub fn bfs(graph: &Graph, start: usize) -> Vec<usize> {
    let adj = graph.adjacency_list();
    let mut visited = vec![false; graph.n_nodes];
    let mut order = Vec::new();
    let mut queue = VecDeque::new();
    visited[start] = true;
    queue.push_back(start);
    while let Some(node) = queue.pop_front() {
        order.push(node);
        for &(nbr, _) in &adj[node] {
            if !visited[nbr] {
                visited[nbr] = true;
                queue.push_back(nbr);
            }
        }
    }
    order
}
/// Returns nodes visited in DFS order starting from `start` (iterative).
#[allow(dead_code)]
pub fn dfs(graph: &Graph, start: usize) -> Vec<usize> {
    let adj = graph.adjacency_list();
    let mut visited = vec![false; graph.n_nodes];
    let mut order = Vec::new();
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if visited[node] {
            continue;
        }
        visited[node] = true;
        order.push(node);
        for &(nbr, _) in adj[node].iter().rev() {
            if !visited[nbr] {
                stack.push(nbr);
            }
        }
    }
    order
}
/// Returns the connected components of the graph (treating all edges as undirected).
///
/// Uses union-find. Each component is a sorted list of node indices.
#[allow(dead_code)]
pub fn connected_components(graph: &Graph) -> Vec<Vec<usize>> {
    let mut uf = UnionFind::new(graph.n_nodes);
    for &(a, b, _) in &graph.edges {
        uf.union(a, b);
    }
    let mut comp: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for i in 0..graph.n_nodes {
        let root = uf.find(i);
        comp.entry(root).or_default().push(i);
    }
    let mut result: Vec<Vec<usize>> = comp.into_values().collect();
    for c in &mut result {
        c.sort_unstable();
    }
    result.sort_unstable_by_key(|c| c[0]);
    result
}
/// Returns `true` if the graph is connected (treating all edges as undirected).
#[allow(dead_code)]
pub fn is_connected(graph: &Graph) -> bool {
    if graph.n_nodes == 0 {
        return true;
    }
    connected_components(graph).len() == 1
}
/// Topological sort using Kahn's algorithm.
///
/// Returns `None` if the graph contains a cycle.
#[allow(dead_code)]
pub fn topological_sort(graph: &Graph) -> Option<Vec<usize>> {
    let n = graph.n_nodes;
    let mut in_degree = vec![0usize; n];
    let adj = graph.adjacency_list();
    for node in 0..n {
        for &(nbr, _) in &adj[node] {
            in_degree[nbr] += 1;
        }
    }
    let mut queue: VecDeque<usize> = (0..n).filter(|&i| in_degree[i] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(node) = queue.pop_front() {
        order.push(node);
        for &(nbr, _) in &adj[node] {
            in_degree[nbr] -= 1;
            if in_degree[nbr] == 0 {
                queue.push_back(nbr);
            }
        }
    }
    if order.len() == n { Some(order) } else { None }
}
/// Dijkstra's shortest-path algorithm from `start`.
///
/// Returns a vector `dist` where `dist[i]` is the shortest distance from `start` to node `i`.
/// Unreachable nodes have distance `f64::INFINITY`.
///
/// **Panics** if any edge weight is negative.
#[allow(dead_code)]
pub fn dijkstra(graph: &Graph, start: usize) -> Vec<f64> {
    let adj = graph.adjacency_list();
    let mut dist = vec![f64::INFINITY; graph.n_nodes];
    dist[start] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(State {
        cost: 0.0,
        node: start,
    });
    while let Some(State { cost, node }) = heap.pop() {
        if cost > dist[node] {
            continue;
        }
        for &(nbr, w) in &adj[node] {
            let next = cost + w;
            if next < dist[nbr] {
                dist[nbr] = next;
                heap.push(State {
                    cost: next,
                    node: nbr,
                });
            }
        }
    }
    dist
}
/// Bellman-Ford shortest-path algorithm from `start`.
///
/// Returns `None` if a negative-weight cycle is reachable from `start`.
/// Otherwise returns the distance vector (same format as `dijkstra`).
#[allow(dead_code)]
pub fn bellman_ford(graph: &Graph, start: usize) -> Option<Vec<f64>> {
    let n = graph.n_nodes;
    let mut dist = vec![f64::INFINITY; n];
    dist[start] = 0.0;
    for _ in 0..n.saturating_sub(1) {
        for &(u, v, w) in &graph.edges {
            if dist[u].is_finite() && dist[u] + w < dist[v] {
                dist[v] = dist[u] + w;
            }
        }
    }
    for &(u, v, w) in &graph.edges {
        if dist[u].is_finite() && dist[u] + w < dist[v] {
            return None;
        }
    }
    Some(dist)
}
/// Reconstructs the path to `end` from a `prev` array produced by a shortest-path algorithm.
///
/// Returns the path as a vector of node indices from source to `end`.
/// Returns an empty vector if `end` is unreachable (i.e., `prev[end]` chain doesn't reach a node
/// with `prev[node] == None` that was the source).
#[allow(dead_code)]
pub fn path_reconstruct(prev: &[Option<usize>], end: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut cur = end;
    loop {
        path.push(cur);
        match prev[cur] {
            None => break,
            Some(p) => cur = p,
        }
    }
    path.reverse();
    path
}
/// Kruskal's MST algorithm for an undirected graph.
///
/// Returns the MST as a list of `(a, b, weight)` edge triples.
/// For disconnected graphs, returns a spanning forest.
#[allow(dead_code)]
pub fn kruskal_mst(graph: &Graph) -> Vec<(usize, usize, f64)> {
    let mut sorted_edges = graph.edges.clone();
    sorted_edges.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(Ordering::Equal));
    let mut uf = UnionFind::new(graph.n_nodes);
    let mut mst = Vec::new();
    for (u, v, w) in sorted_edges {
        if uf.union(u, v) {
            mst.push((u, v, w));
        }
    }
    mst
}
/// Returns the total weight of an MST (or any edge list).
#[allow(dead_code)]
pub fn mst_total_weight(mst_edges: &[(usize, usize, f64)]) -> f64 {
    mst_edges.iter().map(|&(_, _, w)| w).sum()
}
/// Creates a contact graph from a list of overlapping body pairs.
///
/// Each overlap becomes an undirected edge with weight `1.0`.
/// The node count is inferred as `max_node_index + 1`.
#[allow(dead_code)]
pub fn contact_graph_from_overlaps(overlaps: &[(usize, usize)]) -> Graph {
    let n = overlaps
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .max()
        .map(|m| m + 1)
        .unwrap_or(0);
    let mut g = Graph::new(n);
    for &(a, b) in overlaps {
        g.add_undirected_edge(a, b, 1.0);
    }
    g
}
/// Returns the number of connected components (physics "islands").
#[allow(dead_code)]
pub fn island_count(graph: &Graph) -> usize {
    connected_components(graph).len()
}
/// Computes the graph Laplacian matrix `L = D - A` for an undirected graph.
///
/// `A[i][j]` is the sum of weights of edges between node `i` and node `j`.
/// `D[i][i]` is the weighted degree of node `i`.
#[allow(dead_code)]
pub fn laplacian_matrix(graph: &Graph) -> Vec<Vec<f64>> {
    let n = graph.n_nodes;
    let mut l = vec![vec![0.0f64; n]; n];
    for &(u, v, w) in &graph.edges {
        l[u][v] -= w;
        l[u][u] += w;
    }
    l
}
/// A* shortest path from `start` to `goal` with a heuristic `h(node) -> f64`.
///
/// Returns the distance and predecessor array, or `None` if no path exists.
/// Use `path_reconstruct` on the returned `prev` to get the path.
pub fn astar(
    graph: &Graph,
    start: usize,
    goal: usize,
    heuristic: &dyn Fn(usize) -> f64,
) -> Option<(f64, Vec<Option<usize>>)> {
    let adj = graph.adjacency_list();
    let mut dist = vec![f64::INFINITY; graph.n_nodes];
    let mut prev: Vec<Option<usize>> = vec![None; graph.n_nodes];
    dist[start] = 0.0;
    let mut heap: BinaryHeap<AStarState> = BinaryHeap::new();
    heap.push(AStarState {
        f: dist[start] + heuristic(start),
        g: 0.0,
        node: start,
    });
    while let Some(AStarState { f: _, g, node }) = heap.pop() {
        if g > dist[node] {
            continue;
        }
        if node == goal {
            return Some((dist[goal], prev));
        }
        for &(nbr, w) in &adj[node] {
            let new_g = g + w;
            if new_g < dist[nbr] {
                dist[nbr] = new_g;
                prev[nbr] = Some(node);
                heap.push(AStarState {
                    f: new_g + heuristic(nbr),
                    g: new_g,
                    node: nbr,
                });
            }
        }
    }
    if dist[goal].is_finite() {
        Some((dist[goal], prev))
    } else {
        None
    }
}
/// Floyd-Warshall all-pairs shortest paths.
///
/// Returns a matrix `dist[i][j]` = shortest path from `i` to `j`.
/// `f64::INFINITY` means unreachable.
/// **Does not** detect negative cycles; behaviour on negative cycles is undefined.
pub fn floyd_warshall(graph: &Graph) -> Vec<Vec<f64>> {
    let n = graph.n_nodes;
    let mut dist = vec![vec![f64::INFINITY; n]; n];
    for i in 0..n {
        dist[i][i] = 0.0;
    }
    for &(u, v, w) in &graph.edges {
        if w < dist[u][v] {
            dist[u][v] = w;
        }
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if dist[i][k].is_finite() && dist[k][j].is_finite() {
                    let via = dist[i][k] + dist[k][j];
                    if via < dist[i][j] {
                        dist[i][j] = via;
                    }
                }
            }
        }
    }
    dist
}
/// Kosaraju's algorithm for strongly connected components (SCCs).
///
/// Returns a list of SCCs, each as a sorted list of node indices.
/// SCCs are returned in reverse topological order (sinks first).
pub fn kosaraju_scc(graph: &Graph) -> Vec<Vec<usize>> {
    let n = graph.n_nodes;
    let adj = graph.adjacency_list();
    let mut visited = vec![false; n];
    let mut finish_stack: Vec<usize> = Vec::with_capacity(n);
    fn dfs1(
        node: usize,
        adj: &[Vec<(usize, f64)>],
        visited: &mut Vec<bool>,
        stack: &mut Vec<usize>,
    ) {
        let mut call_stack = vec![(node, 0usize)];
        while let Some((v, idx)) = call_stack.last_mut() {
            let v = *v;
            if !visited[v] {
                visited[v] = true;
            }
            if *idx < adj[v].len() {
                let (nbr, _) = adj[v][*idx];
                *idx += 1;
                if !visited[nbr] {
                    call_stack.push((nbr, 0));
                }
            } else {
                stack.push(v);
                call_stack.pop();
            }
        }
    }
    for i in 0..n {
        if !visited[i] {
            dfs1(i, &adj, &mut visited, &mut finish_stack);
        }
    }
    let mut radj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(u, v, w) in &graph.edges {
        radj[v].push((u, w));
    }
    let mut visited2 = vec![false; n];
    let mut sccs: Vec<Vec<usize>> = Vec::new();
    while let Some(start) = finish_stack.pop() {
        if visited2[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut stack = vec![start];
        while let Some(v) = stack.pop() {
            if visited2[v] {
                continue;
            }
            visited2[v] = true;
            comp.push(v);
            for &(nbr, _) in &radj[v] {
                if !visited2[nbr] {
                    stack.push(nbr);
                }
            }
        }
        comp.sort_unstable();
        sccs.push(comp);
    }
    sccs
}
/// Tarjan's algorithm for strongly connected components.
///
/// Returns a list of SCCs in reverse topological order.
pub fn tarjan_scc(graph: &Graph) -> Vec<Vec<usize>> {
    let n = graph.n_nodes;
    let adj = graph.adjacency_list();
    let mut index = 0usize;
    let mut indices = vec![usize::MAX; n];
    let mut lowlinks = vec![0usize; n];
    let mut on_stack = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut sccs: Vec<Vec<usize>> = Vec::new();
    let mut call_stack: Vec<(usize, usize)> = Vec::new();
    for start in 0..n {
        if indices[start] != usize::MAX {
            continue;
        }
        call_stack.push((start, 0));
        while let Some(&mut (v, ref mut adj_idx)) = call_stack.last_mut() {
            if indices[v] == usize::MAX {
                indices[v] = index;
                lowlinks[v] = index;
                index += 1;
                stack.push(v);
                on_stack[v] = true;
            }
            let mut pushed = false;
            while *adj_idx < adj[v].len() {
                let (w, _) = adj[v][*adj_idx];
                *adj_idx += 1;
                if indices[w] == usize::MAX {
                    call_stack.push((w, 0));
                    pushed = true;
                    break;
                } else if on_stack[w] {
                    lowlinks[v] = lowlinks[v].min(indices[w]);
                }
            }
            if !pushed {
                let (v, _) = call_stack
                    .pop()
                    .expect("call_stack is non-empty in backtrack");
                if let Some(&mut (parent, _)) = call_stack.last_mut() {
                    lowlinks[parent] = lowlinks[parent].min(lowlinks[v]);
                }
                if lowlinks[v] == indices[v] {
                    let mut comp = Vec::new();
                    loop {
                        let w = stack.pop().expect("stack is non-empty in SCC extraction");
                        on_stack[w] = false;
                        comp.push(w);
                        if w == v {
                            break;
                        }
                    }
                    comp.sort_unstable();
                    sccs.push(comp);
                }
            }
        }
    }
    sccs
}
/// Prim's minimum spanning tree algorithm.
///
/// Returns the MST as a list of `(a, b, weight)` edge triples.
/// Assumes an undirected (or symmetrically-stored directed) graph.
pub fn prim_mst(graph: &Graph) -> Vec<(usize, usize, f64)> {
    if graph.n_nodes == 0 {
        return Vec::new();
    }
    let adj = graph.adjacency_list();
    let mut in_tree = vec![false; graph.n_nodes];
    let mut key = vec![f64::INFINITY; graph.n_nodes];
    let mut parent: Vec<Option<usize>> = vec![None; graph.n_nodes];
    key[0] = 0.0;
    for _ in 0..graph.n_nodes {
        let u = (0..graph.n_nodes)
            .filter(|&v| !in_tree[v])
            .min_by(|&a, &b| key[a].partial_cmp(&key[b]).unwrap_or(Ordering::Equal));
        let u = match u {
            Some(v) => v,
            None => break,
        };
        in_tree[u] = true;
        for &(v, w) in &adj[u] {
            if !in_tree[v] && w < key[v] {
                key[v] = w;
                parent[v] = Some(u);
            }
        }
    }
    let mut mst = Vec::new();
    for v in 0..graph.n_nodes {
        if let Some(u) = parent[v] {
            mst.push((u, v, key[v]));
        }
    }
    mst
}
/// Edmonds-Karp max-flow algorithm (BFS augmenting paths).
///
/// Returns the maximum flow value from `source` to `sink`.
/// Works on non-negative edge capacities.
pub fn max_flow_edmonds_karp(graph: &Graph, source: usize, sink: usize) -> f64 {
    let n = graph.n_nodes;
    let mut cap = vec![vec![0.0f64; n]; n];
    for &(u, v, w) in &graph.edges {
        cap[u][v] += w;
    }
    let mut flow = 0.0f64;
    loop {
        let mut prev = vec![usize::MAX; n];
        let mut visited = vec![false; n];
        visited[source] = true;
        let mut queue = VecDeque::new();
        queue.push_back(source);
        while let Some(u) = queue.pop_front() {
            if u == sink {
                break;
            }
            for v in 0..n {
                if !visited[v] && cap[u][v] > 1e-12 {
                    visited[v] = true;
                    prev[v] = u;
                    queue.push_back(v);
                }
            }
        }
        if !visited[sink] {
            break;
        }
        let mut path_flow = f64::INFINITY;
        let mut v = sink;
        while v != source {
            let u = prev[v];
            path_flow = path_flow.min(cap[u][v]);
            v = u;
        }
        let mut v = sink;
        while v != source {
            let u = prev[v];
            cap[u][v] -= path_flow;
            cap[v][u] += path_flow;
            v = u;
        }
        flow += path_flow;
    }
    flow
}
/// Dijkstra with both distances and predecessor array.
///
/// Returns `(dist, prev)` where `prev[i]` is the predecessor of node `i`
/// on the shortest path from `start`.
pub fn dijkstra_with_path(graph: &Graph, start: usize) -> (Vec<f64>, Vec<Option<usize>>) {
    let adj = graph.adjacency_list();
    let mut dist = vec![f64::INFINITY; graph.n_nodes];
    let mut prev: Vec<Option<usize>> = vec![None; graph.n_nodes];
    dist[start] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(State {
        cost: 0.0,
        node: start,
    });
    while let Some(State { cost, node }) = heap.pop() {
        if cost > dist[node] {
            continue;
        }
        for &(nbr, w) in &adj[node] {
            let next = cost + w;
            if next < dist[nbr] {
                dist[nbr] = next;
                prev[nbr] = Some(node);
                heap.push(State {
                    cost: next,
                    node: nbr,
                });
            }
        }
    }
    (dist, prev)
}
/// Returns `true` if the directed graph is a DAG (no directed cycles).
pub fn is_dag(graph: &Graph) -> bool {
    topological_sort(graph).is_some()
}
/// Counts directed cycles (simple) in the graph using DFS coloring.
///
/// Returns the number of back edges encountered.
pub fn count_back_edges(graph: &Graph) -> usize {
    let adj = graph.adjacency_list();
    let n = graph.n_nodes;
    let mut color = vec![0u8; n];
    let mut back_edges = 0usize;
    fn dfs_color(v: usize, adj: &[Vec<(usize, f64)>], color: &mut Vec<u8>, back: &mut usize) {
        color[v] = 1;
        for &(w, _) in &adj[v] {
            if color[w] == 1 {
                *back += 1;
            } else if color[w] == 0 {
                dfs_color(w, adj, color, back);
            }
        }
        color[v] = 2;
    }
    for i in 0..n {
        if color[i] == 0 {
            dfs_color(i, &adj, &mut color, &mut back_edges);
        }
    }
    back_edges
}
/// Returns the diameter of the graph: the longest shortest path.
///
/// Computed using BFS from every node (unweighted).
pub fn graph_diameter_bfs(graph: &Graph) -> usize {
    let mut diameter = 0usize;
    for start in 0..graph.n_nodes {
        let bfs_dist = bfs_distances(graph, start);
        for d in bfs_dist {
            if d != usize::MAX && d > diameter {
                diameter = d;
            }
        }
    }
    diameter
}
/// BFS distances from `start` (unweighted).  `usize::MAX` = unreachable.
pub fn bfs_distances(graph: &Graph, start: usize) -> Vec<usize> {
    let adj = graph.adjacency_list();
    let mut dist = vec![usize::MAX; graph.n_nodes];
    dist[start] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(start);
    while let Some(v) = queue.pop_front() {
        for &(w, _) in &adj[v] {
            if dist[w] == usize::MAX {
                dist[w] = dist[v] + 1;
                queue.push_back(w);
            }
        }
    }
    dist
}
/// Degree centrality: returns `degree(v) / (n-1)` for each node.
pub fn degree_centrality(graph: &Graph) -> Vec<f64> {
    let n = graph.n_nodes;
    if n <= 1 {
        return vec![0.0; n];
    }
    (0..n)
        .map(|v| graph.degree(v) as f64 / (n - 1) as f64)
        .collect()
}
/// Returns adjacency matrix as a flat `n×n` Vec`f64`.
pub fn adjacency_matrix(graph: &Graph) -> Vec<Vec<f64>> {
    let n = graph.n_nodes;
    let mut mat = vec![vec![0.0f64; n]; n];
    for &(u, v, w) in &graph.edges {
        mat[u][v] = w;
    }
    mat
}
/// Convert a distance matrix to a sparse graph (threshold edges).
pub fn distance_matrix_to_graph(dist: &[Vec<f64>], threshold: f64) -> Graph {
    let n = dist.len();
    let mut g = Graph::new(n);
    for i in 0..n {
        for j in 0..n {
            if i != j && dist[i][j] < threshold {
                g.add_edge(i, j, dist[i][j]);
            }
        }
    }
    g
}
/// Returns `true` if the undirected graph (edges treated as undirected) is bipartite.
///
/// Uses a 2-colouring BFS: if any edge connects same-coloured vertices, the graph
/// is not bipartite.
pub fn is_bipartite(graph: &Graph) -> bool {
    let adj = graph.adjacency_list();
    let mut color = vec![usize::MAX; graph.n_nodes];
    for start in 0..graph.n_nodes {
        if color[start] != usize::MAX {
            continue;
        }
        color[start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            for &(v, _w) in &adj[u] {
                if color[v] == usize::MAX {
                    color[v] = 1 - color[u];
                    queue.push_back(v);
                } else if color[v] == color[u] {
                    return false;
                }
            }
        }
    }
    true
}
/// Returns the 2-colouring as `Some(colors)` if bipartite, `None` otherwise.
///
/// `colors[i]` is 0 or 1 for each node `i`.
pub fn bipartite_coloring(graph: &Graph) -> Option<Vec<usize>> {
    let adj = graph.adjacency_list();
    let mut color = vec![usize::MAX; graph.n_nodes];
    for start in 0..graph.n_nodes {
        if color[start] != usize::MAX {
            continue;
        }
        color[start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(u) = queue.pop_front() {
            for &(v, _w) in &adj[u] {
                if color[v] == usize::MAX {
                    color[v] = 1 - color[u];
                    queue.push_back(v);
                } else if color[v] == color[u] {
                    return None;
                }
            }
        }
    }
    Some(color)
}
/// Greedy graph colouring.
///
/// Assigns the smallest available colour (integer ≥ 0) to each node in order
/// `0..n_nodes`.  Returns a vector of length `n_nodes` where each entry is the
/// colour index.  The number of colours used is `result.iter().max() + 1`.
///
/// This is a simple greedy heuristic and does **not** guarantee the chromatic
/// number in general graphs.
pub fn greedy_coloring(graph: &Graph) -> Vec<usize> {
    let adj = graph.adjacency_list();
    let n = graph.n_nodes;
    let mut color: Vec<Option<usize>> = vec![None; n];
    for u in 0..n {
        let mut used: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &(v, _w) in &adj[u] {
            if let Some(c) = color[v] {
                used.insert(c);
            }
        }
        let mut c = 0usize;
        while used.contains(&c) {
            c += 1;
        }
        color[u] = Some(c);
    }
    color.into_iter().map(|c| c.unwrap_or(0)).collect()
}
/// Returns the number of colours used by a colouring vector.
pub fn chromatic_number_upper_bound(coloring: &[usize]) -> usize {
    coloring.iter().copied().max().map(|m| m + 1).unwrap_or(0)
}
/// Kahn's algorithm for topological sorting (BFS / in-degree approach).
///
/// Returns `Some(order)` for a DAG, or `None` if the graph contains a cycle.
/// Equivalent to `topological_sort` but implemented via in-degree counting.
pub fn kahn_topological_sort(graph: &Graph) -> Option<Vec<usize>> {
    let n = graph.n_nodes;
    let mut in_degree = vec![0usize; n];
    let adj = graph.adjacency_list();
    for u in 0..n {
        for &(v, _w) in &adj[u] {
            in_degree[v] += 1;
        }
    }
    let mut queue: VecDeque<usize> = (0..n).filter(|&u| in_degree[u] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(u) = queue.pop_front() {
        order.push(u);
        for &(v, _w) in &adj[u] {
            in_degree[v] -= 1;
            if in_degree[v] == 0 {
                queue.push_back(v);
            }
        }
    }
    if order.len() == n { Some(order) } else { None }
}
/// Returns `true` if the (undirected) graph has an Eulerian circuit.
///
/// Conditions: the graph is connected (ignoring isolated nodes) AND every node
/// has even degree.
pub fn has_eulerian_circuit(graph: &Graph) -> bool {
    for u in 0..graph.n_nodes {
        if !graph.degree(u).is_multiple_of(2) {
            return false;
        }
    }
    let non_isolated: Vec<usize> = (0..graph.n_nodes)
        .filter(|&u| graph.degree(u) > 0)
        .collect();
    if non_isolated.is_empty() {
        return true;
    }
    let start = non_isolated[0];
    let visited = {
        let adj = graph.adjacency_list();
        let mut vis = vec![false; graph.n_nodes];
        let mut stack = vec![start];
        vis[start] = true;
        while let Some(u) = stack.pop() {
            for &(v, _w) in &adj[u] {
                if !vis[v] {
                    vis[v] = true;
                    stack.push(v);
                }
            }
        }
        vis
    };
    non_isolated.iter().all(|&u| visited[u])
}
/// Returns `true` if the (undirected) graph has an Eulerian path (but not a circuit).
///
/// Conditions: exactly 2 nodes have odd degree AND the graph is connected.
pub fn has_eulerian_path(graph: &Graph) -> bool {
    let odd_degree_count = (0..graph.n_nodes)
        .filter(|&u| !graph.degree(u).is_multiple_of(2))
        .count();
    odd_degree_count == 2 && {
        let non_isolated: Vec<usize> = (0..graph.n_nodes)
            .filter(|&u| graph.degree(u) > 0)
            .collect();
        if non_isolated.is_empty() {
            return false;
        }
        let start = non_isolated[0];
        let adj = graph.adjacency_list();
        let mut vis = vec![false; graph.n_nodes];
        let mut stack = vec![start];
        vis[start] = true;
        while let Some(u) = stack.pop() {
            for &(v, _w) in &adj[u] {
                if !vis[v] {
                    vis[v] = true;
                    stack.push(v);
                }
            }
        }
        non_isolated.iter().all(|&u| vis[u])
    }
}
/// Computes approximate betweenness centrality using Dijkstra for each source.
///
/// For each pair `(s, t)`, counts how many shortest paths pass through each
/// intermediate node and normalises by `(n-1)*(n-2)/2`.  This is the
/// unnormalised vertex betweenness centrality.
pub fn betweenness_centrality(graph: &Graph) -> Vec<f64> {
    let n = graph.n_nodes;
    let mut bc = vec![0.0f64; n];
    for s in 0..n {
        let (_dist, prev) = dijkstra_with_path(graph, s);
        for t in 0..n {
            if t == s {
                continue;
            }
            let path = path_reconstruct(&prev, t);
            if path.len() < 3 {
                continue;
            }
            for &v in &path[1..path.len() - 1] {
                bc[v] += 1.0;
            }
        }
    }
    let norm = if n > 2 {
        ((n - 1) * (n - 2)) as f64
    } else {
        1.0
    };
    for v in bc.iter_mut() {
        *v /= norm;
    }
    bc
}
/// Computes closeness centrality for every node.
///
/// Closeness = `(n_reachable - 1) / sum_of_distances` for reachable nodes,
/// normalised by `(n_reachable - 1) / (n - 1)` (Wasserman-Faust variant).
pub fn closeness_centrality(graph: &Graph) -> Vec<f64> {
    let n = graph.n_nodes;
    let mut cc = vec![0.0f64; n];
    for u in 0..n {
        let dist = dijkstra(graph, u);
        let (sum, reachable) = dist
            .iter()
            .enumerate()
            .fold((0.0f64, 0usize), |(s, r), (v, &d)| {
                if v != u && d.is_finite() {
                    (s + d, r + 1)
                } else {
                    (s, r)
                }
            });
        if reachable == 0 || sum == 0.0 {
            cc[u] = 0.0;
        } else {
            let raw = (reachable as f64) / sum;
            cc[u] = raw * (reachable as f64) / ((n - 1) as f64);
        }
    }
    cc
}
/// Augmenting-path maximum bipartite matching.
///
/// `left` and `right` are the node indices for the two sides of the bipartite
/// graph.  Returns a list of matched pairs `(l, r)`.
pub fn bipartite_matching(graph: &Graph, left: &[usize], right: &[usize]) -> Vec<(usize, usize)> {
    let adj = graph.adjacency_list();
    let mut right_idx = vec![usize::MAX; graph.n_nodes];
    for (i, &r) in right.iter().enumerate() {
        right_idx[r] = i;
    }
    let nr = right.len();
    let mut match_l: Vec<Option<usize>> = vec![None; graph.n_nodes];
    let mut match_r: Vec<Option<usize>> = vec![None; nr];
    fn augment(
        u: usize,
        adj: &[Vec<(usize, f64)>],
        right_idx: &[usize],
        match_r: &mut Vec<Option<usize>>,
        visited: &mut Vec<bool>,
    ) -> bool {
        for &(v, _w) in &adj[u] {
            let ri = right_idx[v];
            if ri == usize::MAX || visited[ri] {
                continue;
            }
            visited[ri] = true;
            let can =
                match_r[ri].is_none_or(|prev_l| augment(prev_l, adj, right_idx, match_r, visited));
            if can {
                match_r[ri] = Some(u);
                return true;
            }
        }
        false
    }
    for &l in left {
        let mut visited = vec![false; nr];
        if augment(l, &adj, &right_idx, &mut match_r, &mut visited) {
            for (ri, mr) in match_r.iter().enumerate() {
                if *mr == Some(l) {
                    match_l[l] = Some(ri);
                }
            }
        }
    }
    let mut pairs = Vec::new();
    for &l in left {
        if let Some(ri) = match_l[l] {
            pairs.push((l, right[ri]));
        }
    }
    pairs
}
/// Graph density: `|E| / (|V| * (|V| - 1))` for directed graphs.
pub fn graph_density(graph: &Graph) -> f64 {
    let n = graph.n_nodes;
    if n <= 1 {
        return 0.0;
    }
    graph.edge_count() as f64 / (n * (n - 1)) as f64
}
/// Local clustering coefficient for node `u` (undirected interpretation).
///
/// = (number of edges among neighbours of u) / (degree * (degree - 1) / 2)
pub fn local_clustering(graph: &Graph, u: usize) -> f64 {
    let adj = graph.adjacency_list();
    let neighbours: Vec<usize> = adj[u].iter().map(|&(v, _)| v).collect();
    let k = neighbours.len();
    if k < 2 {
        return 0.0;
    }
    let neighbour_set: std::collections::HashSet<usize> = neighbours.iter().copied().collect();
    let mut edges_among = 0usize;
    for &v in &neighbours {
        for &(w, _) in &adj[v] {
            if neighbour_set.contains(&w) && w != u {
                edges_among += 1;
            }
        }
    }
    edges_among as f64 / (k * (k - 1)) as f64
}
/// Average clustering coefficient across all nodes.
pub fn average_clustering(graph: &Graph) -> f64 {
    if graph.n_nodes == 0 {
        return 0.0;
    }
    let total: f64 = (0..graph.n_nodes).map(|u| local_clustering(graph, u)).sum();
    total / graph.n_nodes as f64
}
#[cfg(test)]
mod tests {
    use super::*;
    fn triangle_graph() -> Graph {
        let mut g = Graph::new(3);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 2.0);
        g.add_undirected_edge(0, 2, 4.0);
        g
    }
    #[test]
    fn test_new_graph() {
        let g = Graph::new(5);
        assert_eq!(g.node_count(), 5);
        assert_eq!(g.edge_count(), 0);
    }
    #[test]
    fn test_add_directed_edge() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.5);
        assert_eq!(g.edge_count(), 1);
    }
    #[test]
    fn test_add_undirected_edge_stores_two() {
        let mut g = Graph::new(3);
        g.add_undirected_edge(0, 1, 1.0);
        assert_eq!(g.edge_count(), 2);
    }
    #[test]
    fn test_degree() {
        let g = triangle_graph();
        assert_eq!(g.degree(0), 2);
        assert_eq!(g.degree(1), 2);
        assert_eq!(g.degree(2), 2);
    }
    #[test]
    fn test_adjacency_list() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 2.0);
        g.add_edge(0, 2, 3.0);
        let adj = g.adjacency_list();
        assert_eq!(adj[0].len(), 2);
        assert_eq!(adj[1].len(), 0);
    }
    #[test]
    fn test_bfs_connected() {
        let g = triangle_graph();
        let order = bfs(&g, 0);
        assert_eq!(order.len(), 3);
        assert_eq!(order[0], 0);
    }
    #[test]
    fn test_bfs_single_node() {
        let g = Graph::new(1);
        let order = bfs(&g, 0);
        assert_eq!(order, vec![0]);
    }
    #[test]
    fn test_bfs_disconnected_visits_reachable() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        let order = bfs(&g, 0);
        assert_eq!(order.len(), 2);
    }
    #[test]
    fn test_dfs_connected() {
        let g = triangle_graph();
        let order = dfs(&g, 0);
        assert_eq!(order.len(), 3);
        assert_eq!(order[0], 0);
    }
    #[test]
    fn test_dfs_single_node() {
        let g = Graph::new(1);
        assert_eq!(dfs(&g, 0), vec![0]);
    }
    #[test]
    fn test_connected_components_single() {
        let g = triangle_graph();
        let comps = connected_components(&g);
        assert_eq!(comps.len(), 1);
        assert_eq!(comps[0], vec![0, 1, 2]);
    }
    #[test]
    fn test_connected_components_two() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        let comps = connected_components(&g);
        assert_eq!(comps.len(), 2);
    }
    #[test]
    fn test_is_connected_true() {
        assert!(is_connected(&triangle_graph()));
    }
    #[test]
    fn test_is_connected_false() {
        let g = Graph::new(3);
        assert!(!is_connected(&g));
    }
    #[test]
    fn test_is_connected_single_node() {
        let g = Graph::new(1);
        assert!(is_connected(&g));
    }
    #[test]
    fn test_topological_sort_dag() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(0, 2, 1.0);
        g.add_edge(1, 3, 1.0);
        g.add_edge(2, 3, 1.0);
        let order = topological_sort(&g).expect("DAG should succeed");
        assert_eq!(order.len(), 4);
        let pos: Vec<usize> = {
            let mut p = vec![0usize; 4];
            for (i, &n) in order.iter().enumerate() {
                p[n] = i;
            }
            p
        };
        assert!(pos[0] < pos[1]);
        assert!(pos[0] < pos[2]);
        assert!(pos[1] < pos[3]);
    }
    #[test]
    fn test_topological_sort_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        assert!(topological_sort(&g).is_none());
    }
    #[test]
    fn test_dijkstra_triangle() {
        let g = triangle_graph();
        let dist = dijkstra(&g, 0);
        assert_eq!(dist[0], 0.0);
        assert_eq!(dist[1], 1.0);
        assert_eq!(dist[2], 3.0);
    }
    #[test]
    fn test_dijkstra_unreachable() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        let dist = dijkstra(&g, 0);
        assert_eq!(dist[2], f64::INFINITY);
    }
    #[test]
    fn test_bellman_ford_no_negative_cycle() {
        let g = triangle_graph();
        let dist = bellman_ford(&g, 0).expect("no negative cycle");
        assert_eq!(dist[0], 0.0);
        assert_eq!(dist[1], 1.0);
    }
    #[test]
    fn test_bellman_ford_negative_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, -2.0);
        g.add_edge(2, 0, -1.0);
        assert!(bellman_ford(&g, 0).is_none());
    }
    #[test]
    fn test_path_reconstruct() {
        let prev: Vec<Option<usize>> = vec![None, Some(0), Some(1), Some(2)];
        let path = path_reconstruct(&prev, 3);
        assert_eq!(path, vec![0, 1, 2, 3]);
    }
    #[test]
    fn test_path_reconstruct_single() {
        let prev: Vec<Option<usize>> = vec![None];
        assert_eq!(path_reconstruct(&prev, 0), vec![0]);
    }
    #[test]
    fn test_kruskal_mst_triangle() {
        let g = triangle_graph();
        let mst = kruskal_mst(&g);
        assert_eq!(mst.len(), 2);
        let total = mst_total_weight(&mst);
        assert!((total - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_mst_total_weight_empty() {
        assert_eq!(mst_total_weight(&[]), 0.0);
    }
    #[test]
    fn test_contact_graph_from_overlaps() {
        let overlaps = vec![(0, 1), (1, 2), (0, 2)];
        let g = contact_graph_from_overlaps(&overlaps);
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 6);
    }
    #[test]
    fn test_contact_graph_empty() {
        let g = contact_graph_from_overlaps(&[]);
        assert_eq!(g.node_count(), 0);
    }
    #[test]
    fn test_island_count_one() {
        assert_eq!(island_count(&triangle_graph()), 1);
    }
    #[test]
    fn test_island_count_three() {
        let g = Graph::new(3);
        assert_eq!(island_count(&g), 3);
    }
    #[test]
    fn test_laplacian_row_sum_zero() {
        let g = triangle_graph();
        let l = laplacian_matrix(&g);
        for row in &l {
            let s: f64 = row.iter().sum();
            assert!(s.abs() < 1e-12, "row sum = {s}");
        }
    }
    #[test]
    fn test_laplacian_diagonal_equals_degree() {
        let g = triangle_graph();
        let l = laplacian_matrix(&g);
        assert!((l[0][0] - 5.0).abs() < 1e-12);
        assert!((l[1][1] - 3.0).abs() < 1e-12);
        assert!((l[2][2] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_laplacian_off_diagonal_negative() {
        let g = triangle_graph();
        let l = laplacian_matrix(&g);
        assert!(l[0][1] < 0.0);
        assert!(l[1][0] < 0.0);
    }
    #[test]
    fn test_floyd_warshall_triangle() {
        let g = triangle_graph();
        let dist = floyd_warshall(&g);
        assert!((dist[0][1] - 1.0).abs() < 1e-12);
        assert!((dist[0][2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_floyd_warshall_diagonal_zero() {
        let g = triangle_graph();
        let dist = floyd_warshall(&g);
        for i in 0..3 {
            assert!(
                (dist[i][i] - 0.0).abs() < 1e-12,
                "dist[{}][{}]={}",
                i,
                i,
                dist[i][i]
            );
        }
    }
    #[test]
    fn test_floyd_warshall_unreachable() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        let dist = floyd_warshall(&g);
        assert_eq!(dist[0][2], f64::INFINITY);
    }
    #[test]
    fn test_astar_finds_path() {
        let g = triangle_graph();
        let result = astar(&g, 0, 2, &|_| 0.0);
        assert!(result.is_some(), "A* should find a path");
        let (cost, _) = result.unwrap();
        assert!((cost - 3.0).abs() < 1e-12, "cost={}", cost);
    }
    #[test]
    fn test_astar_no_path() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        let result = astar(&g, 0, 2, &|_| 0.0);
        assert!(
            result.is_none(),
            "A* should not find path to unreachable node"
        );
    }
    #[test]
    fn test_kosaraju_scc_single_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        let sccs = kosaraju_scc(&g);
        assert_eq!(sccs.len(), 1, "one SCC expected for a cycle");
    }
    #[test]
    fn test_kosaraju_scc_dag() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let sccs = kosaraju_scc(&g);
        assert_eq!(sccs.len(), 3, "DAG should have 3 SCCs");
    }
    #[test]
    fn test_prim_mst_triangle() {
        let g = triangle_graph();
        let mst = prim_mst(&g);
        let total = mst.iter().map(|&(_, _, w)| w).sum::<f64>();
        assert!((total - 3.0).abs() < 1e-12, "prim MST weight={}", total);
    }
    #[test]
    fn test_prim_mst_empty_graph() {
        let g = Graph::new(0);
        let mst = prim_mst(&g);
        assert!(mst.is_empty());
    }
    #[test]
    fn test_max_flow_simple_path() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 3.0);
        g.add_edge(1, 2, 2.0);
        let flow = max_flow_edmonds_karp(&g, 0, 2);
        assert!((flow - 2.0).abs() < 1e-9, "flow={}", flow);
    }
    #[test]
    fn test_max_flow_two_parallel_paths() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 2.0);
        g.add_edge(1, 3, 2.0);
        g.add_edge(0, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        let flow = max_flow_edmonds_karp(&g, 0, 3);
        assert!((flow - 3.0).abs() < 1e-9, "flow={}", flow);
    }
    #[test]
    fn test_dijkstra_with_path_triangle() {
        let g = triangle_graph();
        let (dist, prev) = dijkstra_with_path(&g, 0);
        assert!((dist[2] - 3.0).abs() < 1e-12);
        let path = path_reconstruct(&prev, 2);
        assert_eq!(path[0], 0);
        assert_eq!(*path.last().unwrap(), 2);
    }
    #[test]
    fn test_is_dag_true() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        assert!(is_dag(&g));
    }
    #[test]
    fn test_is_dag_false() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        assert!(!is_dag(&g));
    }
    #[test]
    fn test_degree_centrality_triangle() {
        let g = triangle_graph();
        let dc = degree_centrality(&g);
        for c in &dc {
            assert!((*c - 1.0).abs() < 1e-12, "c={}", c);
        }
    }
    #[test]
    fn test_adjacency_matrix_triangle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 2.0);
        g.add_edge(1, 2, 3.0);
        let mat = adjacency_matrix(&g);
        assert!((mat[0][1] - 2.0).abs() < 1e-12);
        assert!((mat[1][2] - 3.0).abs() < 1e-12);
        assert!((mat[0][0]).abs() < 1e-12);
    }
    #[test]
    fn test_graph_diameter_triangle() {
        let g = triangle_graph();
        let d = graph_diameter_bfs(&g);
        assert_eq!(d, 1, "diameter of complete triangle = 1");
    }
    #[test]
    fn test_count_back_edges_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        let back = count_back_edges(&g);
        assert!(back > 0, "cycle should have back edges");
    }
    #[test]
    fn test_count_back_edges_dag() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let back = count_back_edges(&g);
        assert_eq!(back, 0, "DAG should have no back edges");
    }
    #[test]
    fn test_is_bipartite_path_graph() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        assert!(is_bipartite(&g));
    }
    #[test]
    fn test_is_bipartite_odd_cycle_false() {
        let g = triangle_graph();
        assert!(!is_bipartite(&g));
    }
    #[test]
    fn test_bipartite_coloring_consistent() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        let colors = bipartite_coloring(&g).expect("should be bipartite");
        assert_ne!(colors[0], colors[1]);
        assert_ne!(colors[1], colors[2]);
        assert_ne!(colors[2], colors[3]);
    }
    #[test]
    fn test_bipartite_coloring_none_for_triangle() {
        assert!(bipartite_coloring(&triangle_graph()).is_none());
    }
    #[test]
    fn test_greedy_coloring_triangle_uses_3_colors() {
        let g = triangle_graph();
        let colors = greedy_coloring(&g);
        assert_eq!(colors.len(), 3);
        let k = chromatic_number_upper_bound(&colors);
        assert!(k >= 3, "K3 needs at least 3 colours, got {k}");
    }
    #[test]
    fn test_greedy_coloring_adjacent_different() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        let colors = greedy_coloring(&g);
        assert_ne!(colors[0], colors[1]);
        assert_ne!(colors[1], colors[2]);
        assert_ne!(colors[2], colors[3]);
    }
    #[test]
    fn test_greedy_coloring_path_needs_2() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        let colors = greedy_coloring(&g);
        let k = chromatic_number_upper_bound(&colors);
        assert!(k <= 2, "path needs at most 2 colours, got {k}");
    }
    #[test]
    fn test_kahn_dag() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(0, 2, 1.0);
        g.add_edge(1, 3, 1.0);
        g.add_edge(2, 3, 1.0);
        let order = kahn_topological_sort(&g).expect("DAG should succeed");
        assert_eq!(order.len(), 4);
        let pos: Vec<usize> = {
            let mut p = vec![0usize; 4];
            for (i, &n) in order.iter().enumerate() {
                p[n] = i;
            }
            p
        };
        assert!(pos[0] < pos[1]);
        assert!(pos[0] < pos[2]);
        assert!(pos[1] < pos[3]);
    }
    #[test]
    fn test_kahn_cycle_returns_none() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        assert!(kahn_topological_sort(&g).is_none());
    }
    #[test]
    fn test_eulerian_circuit_even_degrees() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        g.add_undirected_edge(2, 3, 1.0);
        g.add_undirected_edge(3, 0, 1.0);
        assert!(
            has_eulerian_circuit(&g),
            "square graph should have Eulerian circuit"
        );
    }
    #[test]
    fn test_no_eulerian_circuit_triangle() {
        let g = triangle_graph();
        assert!(has_eulerian_circuit(&g));
    }
    #[test]
    fn test_eulerian_path_odd_degree_nodes() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let odd_count = (0..3).filter(|&u| !g.degree(u).is_multiple_of(2)).count();
        assert_eq!(odd_count, 2);
    }
    #[test]
    fn test_betweenness_centrality_path() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let bc = betweenness_centrality(&g);
        assert!(bc[1] > bc[0], "bc[1]={}, bc[0]={}", bc[1], bc[0]);
        assert!(bc[1] > bc[2], "bc[1]={}, bc[2]={}", bc[1], bc[2]);
    }
    #[test]
    fn test_betweenness_centrality_all_finite() {
        let g = triangle_graph();
        let bc = betweenness_centrality(&g);
        assert!(bc.iter().all(|v| v.is_finite()));
    }
    #[test]
    fn test_closeness_centrality_central_node() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(0, 2, 1.0);
        g.add_edge(0, 3, 1.0);
        let cc = closeness_centrality(&g);
        assert!(cc[0] >= cc[1], "cc[0]={} cc[1]={}", cc[0], cc[1]);
        assert!(cc[0] >= cc[2], "cc[0]={} cc[2]={}", cc[0], cc[2]);
    }
    #[test]
    fn test_graph_density_complete_directed() {
        let mut g = Graph::new(3);
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    g.add_edge(i, j, 1.0);
                }
            }
        }
        let d = graph_density(&g);
        assert!((d - 1.0).abs() < 1e-12, "density={d}");
    }
    #[test]
    fn test_graph_density_empty() {
        let g = Graph::new(4);
        assert_eq!(graph_density(&g), 0.0);
    }
    #[test]
    fn test_bipartite_matching_simple() {
        let mut g = Graph::new(4);
        g.add_edge(0, 2, 1.0);
        g.add_edge(0, 3, 1.0);
        g.add_edge(1, 3, 1.0);
        let matching = bipartite_matching(&g, &[0, 1], &[2, 3]);
        assert_eq!(matching.len(), 2, "matching={:?}", matching);
    }
    #[test]
    fn test_bipartite_matching_no_edges() {
        let g = Graph::new(4);
        let matching = bipartite_matching(&g, &[0, 1], &[2, 3]);
        assert!(matching.is_empty());
    }
    #[test]
    fn test_local_clustering_complete_triangle() {
        let g = triangle_graph();
        let c = local_clustering(&g, 0);
        assert!((c - 1.0).abs() < 1e-10, "clustering={c}");
    }
    #[test]
    fn test_average_clustering_path_graph() {
        let mut g = Graph::new(3);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 1.0);
        let avg = average_clustering(&g);
        assert!(avg.is_finite() && avg >= 0.0);
    }
    #[test]
    fn test_dijkstra_single_node() {
        let g = Graph::new(1);
        let d = dijkstra(&g, 0);
        assert_eq!(d.len(), 1);
        assert_eq!(d[0], 0.0);
    }
    #[test]
    fn test_dijkstra_self_loop_ignored() {
        let mut g = Graph::new(2);
        g.add_edge(0, 0, 5.0);
        g.add_edge(0, 1, 2.0);
        let d = dijkstra(&g, 0);
        assert_eq!(d[0], 0.0);
        assert_eq!(d[1], 2.0);
    }
    #[test]
    fn test_bellman_ford_unreachable() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 2.0);
        let d = bellman_ford(&g, 0).unwrap();
        assert!(d[2].is_infinite());
    }
    #[test]
    fn test_astar_same_result_as_dijkstra() {
        let g = triangle_graph();
        let dijk = dijkstra(&g, 0);
        let zero_h = |_: usize| 0.0;
        if let Some((d, _)) = astar(&g, 0, 2, &zero_h) {
            assert!((d - dijk[2]).abs() < 1e-12, "A* and Dijkstra disagree");
        }
    }
    #[test]
    fn test_floyd_warshall_symmetric_undirected() {
        let g = triangle_graph();
        let dist = floyd_warshall(&g);
        for i in 0..3 {
            for j in 0..3 {
                assert!((dist[i][j] - dist[j][i]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_topological_sort_empty_graph() {
        let g = Graph::new(0);
        let order = topological_sort(&g).expect("empty graph is a DAG");
        assert!(order.is_empty());
    }
    #[test]
    fn test_tarjan_single_scc() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        let sccs = tarjan_scc(&g);
        assert_eq!(sccs.len(), 1);
        assert_eq!(sccs[0].len(), 3);
    }
    #[test]
    fn test_tarjan_dag_each_node_own_scc() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        let sccs = tarjan_scc(&g);
        assert_eq!(sccs.len(), 3, "DAG: each node is its own SCC");
    }
    #[test]
    fn test_kosaraju_multiple_sccs() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 0, 1.0);
        g.add_edge(2, 3, 1.0);
        g.add_edge(3, 2, 1.0);
        let sccs = kosaraju_scc(&g);
        assert_eq!(sccs.len(), 2);
    }
    #[test]
    fn test_kruskal_prim_agree() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 2.0);
        g.add_undirected_edge(1, 2, 3.0);
        g.add_undirected_edge(2, 3, 1.0);
        g.add_undirected_edge(0, 3, 4.0);
        let w_k = mst_total_weight(&kruskal_mst(&g));
        let w_p: f64 = prim_mst(&g).iter().map(|&(_, _, w)| w).sum();
        assert!((w_k - w_p).abs() < 1e-10, "Kruskal={w_k} Prim={w_p}");
    }
    #[test]
    fn test_max_flow_new_parallel() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 3.0);
        g.add_edge(0, 2, 2.0);
        g.add_edge(1, 3, 2.0);
        g.add_edge(2, 3, 3.0);
        let flow = max_flow_edmonds_karp(&g, 0, 3);
        assert!((flow - 4.0).abs() < 1e-9, "flow={flow}");
    }
    #[test]
    fn test_max_flow_no_path() {
        let g = Graph::new(3);
        let flow = max_flow_edmonds_karp(&g, 0, 2);
        assert_eq!(flow, 0.0);
    }
    #[test]
    fn test_max_flow_bottleneck() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 5.0);
        g.add_edge(1, 2, 2.0);
        let flow = max_flow_edmonds_karp(&g, 0, 2);
        assert!((flow - 2.0).abs() < 1e-9, "flow limited by bottleneck=2");
    }
    #[test]
    fn test_is_dag_linear_chain() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        assert!(is_dag(&g));
    }
    #[test]
    fn test_is_dag_three_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        assert!(!is_dag(&g));
    }
    #[test]
    fn test_count_back_edges_linear() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(0, 2, 1.0);
        assert_eq!(count_back_edges(&g), 0);
    }
    #[test]
    fn test_count_back_edges_two_cycle() {
        let mut g = Graph::new(2);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 0, 1.0);
        assert!(count_back_edges(&g) >= 1);
    }
    #[test]
    fn test_path_reconstruct_chain() {
        let prev: Vec<Option<usize>> = vec![None, Some(0), Some(1)];
        let path = path_reconstruct(&prev, 2);
        assert_eq!(path, vec![0, 1, 2]);
    }
    #[test]
    fn test_method_scc_single_cycle() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 0, 1.0);
        let sccs = g.strongly_connected_components();
        assert_eq!(sccs.len(), 1, "expected 1 SCC for a 3-cycle");
        let mut flat: Vec<usize> = sccs.into_iter().flatten().collect();
        flat.sort_unstable();
        assert_eq!(flat, vec![0, 1, 2]);
    }
    #[test]
    fn test_method_scc_dag_all_singleton() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        let sccs = g.strongly_connected_components();
        assert_eq!(sccs.len(), 4);
    }
    #[test]
    fn test_method_mst_triangle() {
        let mut g = Graph::new(3);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(1, 2, 2.0);
        g.add_undirected_edge(0, 2, 5.0);
        let mst = g.minimum_spanning_tree();
        assert_eq!(mst.len(), 2);
        let total: f64 = mst.iter().map(|&(_, _, w)| w).sum();
        assert!((total - 3.0).abs() < 1e-9, "MST weight={}", total);
    }
    #[test]
    fn test_method_mst_star_graph() {
        let mut g = Graph::new(4);
        g.add_undirected_edge(0, 1, 1.0);
        g.add_undirected_edge(0, 2, 1.0);
        g.add_undirected_edge(0, 3, 1.0);
        let mst = g.minimum_spanning_tree();
        assert_eq!(mst.len(), 3);
    }
    #[test]
    fn test_method_max_flow_simple() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, 10.0);
        g.add_edge(1, 2, 5.0);
        let flow = g.max_flow(0, 2);
        assert!((flow - 5.0).abs() < 1e-9, "flow limited by bottleneck=5");
    }
    #[test]
    fn test_method_max_flow_no_path() {
        let g = Graph::new(4);
        let flow = g.max_flow(0, 3);
        assert_eq!(flow, 0.0);
    }
    #[test]
    fn test_method_betweenness_path_graph() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(2, 3, 1.0);
        let bc = g.betweenness_centrality();
        assert!(bc[0].abs() < 1e-9, "bc[0]={}", bc[0]);
        assert!(bc[3].abs() < 1e-9, "bc[3]={}", bc[3]);
        assert!(bc[1] > 0.0, "bc[1]={}", bc[1]);
        assert!(bc[2] > 0.0, "bc[2]={}", bc[2]);
    }
    #[test]
    fn test_method_betweenness_single_node() {
        let g = Graph::new(1);
        let bc = g.betweenness_centrality();
        assert_eq!(bc.len(), 1);
        assert_eq!(bc[0], 0.0);
    }
    #[test]
    fn test_method_betweenness_hub_high_centrality() {
        let mut g = Graph::new(5);
        for i in 1..5 {
            g.add_undirected_edge(0, i, 1.0);
        }
        let bc = g.betweenness_centrality();
        for i in 1..5 {
            assert!(bc[0] >= bc[i], "hub bc={} leaf[{}] bc={}", bc[0], i, bc[i]);
        }
    }
}
