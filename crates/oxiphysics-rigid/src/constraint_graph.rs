#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constraint graph analysis for rigid body systems.
//!
//! Provides topology analysis of multi-body systems including connected
//! component (island) finding, Grübler mobility analysis, kinematic chain
//! ordering, closed-loop detection, and redundancy counting.

/// A node index representing a rigid body in the constraint graph.
pub type BodyIndex = usize;

/// An edge index representing a constraint between two bodies.
pub type ConstraintIndex = usize;

/// Represents a single constraint connecting two bodies.
///
/// Each constraint has a number of freedoms it removes (`dof_removed`),
/// which equals `6 - f_i` in Grübler's formula.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConstraintEdge {
    /// First body index.
    pub body_a: BodyIndex,
    /// Second body index.
    pub body_b: BodyIndex,
    /// Number of degrees of freedom removed by this constraint (1–6).
    pub dof_removed: u32,
}

/// Constraint graph for a rigid body system.
///
/// Bodies are nodes; constraints are edges. Supports topology queries such as
/// island detection, mobility analysis, and closed-loop detection.
#[derive(Debug, Default, Clone)]
pub struct ConstraintGraph {
    /// Number of bodies (nodes) in the graph.
    body_count: usize,
    /// List of constraints (edges).
    constraints: Vec<ConstraintEdge>,
}

impl ConstraintGraph {
    /// Create an empty constraint graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a body to the graph and return its index.
    pub fn add_body(&mut self) -> BodyIndex {
        let idx = self.body_count;
        self.body_count += 1;
        idx
    }

    /// Add a constraint between two bodies.
    ///
    /// # Arguments
    /// * `body_a` — index of the first body.
    /// * `body_b` — index of the second body.
    /// * `dof_removed` — degrees of freedom removed (1–6).
    ///
    /// # Panics
    /// Panics if either body index is out of range or `dof_removed` is 0.
    pub fn add_constraint(
        &mut self,
        body_a: BodyIndex,
        body_b: BodyIndex,
        dof_removed: u32,
    ) -> ConstraintIndex {
        assert!(body_a < self.body_count, "body_a index out of range");
        assert!(body_b < self.body_count, "body_b index out of range");
        assert!((1..=6).contains(&dof_removed), "dof_removed must be 1–6");
        let idx = self.constraints.len();
        self.constraints.push(ConstraintEdge {
            body_a,
            body_b,
            dof_removed,
        });
        idx
    }

    /// Find connected components (islands) using union-find.
    ///
    /// Returns a `Vec` of length `body_count` where each entry is the
    /// canonical root body index for that island.
    pub fn find_islands(&self) -> Vec<BodyIndex> {
        let n = self.body_count;
        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank: Vec<u32> = vec![0; n];

        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        fn union(parent: &mut Vec<usize>, rank: &mut Vec<u32>, a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra == rb {
                return;
            }
            match rank[ra].cmp(&rank[rb]) {
                std::cmp::Ordering::Less => parent[ra] = rb,
                std::cmp::Ordering::Greater => parent[rb] = ra,
                std::cmp::Ordering::Equal => {
                    parent[rb] = ra;
                    rank[ra] += 1;
                }
            }
        }

        for edge in &self.constraints {
            union(&mut parent, &mut rank, edge.body_a, edge.body_b);
        }

        // Path-compress all
        for i in 0..n {
            find(&mut parent, i);
        }

        parent
    }

    /// Count the number of distinct islands (connected components).
    pub fn island_count(&self) -> usize {
        if self.body_count == 0 {
            return 0;
        }
        let islands = self.find_islands();
        let mut roots: Vec<usize> = islands.clone();
        roots.sort_unstable();
        roots.dedup();
        roots.len()
    }

    /// Grübler mobility formula: `M = 6*(n-1) - Σ(dof_removed_i)`.
    ///
    /// Where `n` is the number of bodies. Returns the degrees of freedom of
    /// the mechanism. Negative values indicate over-constraint.
    pub fn mobility_analysis(&self) -> i64 {
        let n = self.body_count as i64;
        if n == 0 {
            return 0;
        }
        let constraint_sum: i64 = self.constraints.iter().map(|e| e.dof_removed as i64).sum();
        6 * (n - 1) - constraint_sum
    }

    /// Perform a topological sort (Kahn's algorithm) to order bodies in a
    /// kinematic chain.
    ///
    /// Returns `Some(order)` if the graph is a DAG (treating constraints as
    /// directed from `body_a` → `body_b`), or `None` if a cycle exists.
    pub fn kinematic_chain(&self) -> Option<Vec<BodyIndex>> {
        let n = self.body_count;
        let mut in_degree = vec![0usize; n];
        let mut adj: Vec<Vec<usize>> = vec![vec![]; n];

        for edge in &self.constraints {
            adj[edge.body_a].push(edge.body_b);
            in_degree[edge.body_b] += 1;
        }

        let mut queue: std::collections::VecDeque<usize> =
            (0..n).filter(|&i| in_degree[i] == 0).collect();

        let mut order = Vec::with_capacity(n);
        while let Some(node) = queue.pop_front() {
            order.push(node);
            for &next in &adj[node] {
                in_degree[next] -= 1;
                if in_degree[next] == 0 {
                    queue.push_back(next);
                }
            }
        }

        if order.len() == n { Some(order) } else { None }
    }

    /// Detect whether the constraint graph contains any cycles (closed loops).
    ///
    /// Returns `true` if at least one closed loop exists.
    pub fn closed_loop_detection(&self) -> bool {
        // A closed loop exists iff topological sort fails
        self.kinematic_chain().is_none()
    }

    /// Estimate the number of redundant constraints.
    ///
    /// A constraint is redundant when the mechanism mobility drops below zero,
    /// or equivalently: `redundancy = max(0, -mobility)`.
    ///
    /// This is the degree of static indeterminacy for a system of rigid bodies.
    pub fn redundancy_count(&self) -> u64 {
        let m = self.mobility_analysis();
        if m < 0 { (-m) as u64 } else { 0 }
    }

    /// Return the number of bodies in the graph.
    pub fn body_count(&self) -> usize {
        self.body_count
    }

    /// Return the number of constraints in the graph.
    pub fn constraint_count(&self) -> usize {
        self.constraints.len()
    }

    /// Iterate over all constraint edges.
    pub fn constraints(&self) -> &[ConstraintEdge] {
        &self.constraints
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_chain(n: usize, dof_removed: u32) -> ConstraintGraph {
        let mut g = ConstraintGraph::new();
        for _ in 0..n {
            g.add_body();
        }
        for i in 0..n - 1 {
            g.add_constraint(i, i + 1, dof_removed);
        }
        g
    }

    // ── add_body / add_constraint ─────────────────────────────────────────

    #[test]
    fn add_body_increments_count() {
        let mut g = ConstraintGraph::new();
        assert_eq!(g.add_body(), 0);
        assert_eq!(g.add_body(), 1);
        assert_eq!(g.add_body(), 2);
        assert_eq!(g.body_count(), 3);
    }

    #[test]
    fn add_constraint_increments_count() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 5);
        assert_eq!(g.constraint_count(), 1);
    }

    #[test]
    #[should_panic]
    fn add_constraint_invalid_body_panics() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_constraint(0, 99, 3);
    }

    #[test]
    #[should_panic]
    fn add_constraint_zero_dof_panics() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 0);
    }

    // ── find_islands ──────────────────────────────────────────────────────

    #[test]
    fn empty_graph_zero_islands() {
        let g = ConstraintGraph::new();
        assert_eq!(g.island_count(), 0);
    }

    #[test]
    fn single_body_one_island() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        assert_eq!(g.island_count(), 1);
    }

    #[test]
    fn two_disconnected_bodies_two_islands() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        assert_eq!(g.island_count(), 2);
    }

    #[test]
    fn connected_pair_one_island() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 3);
        assert_eq!(g.island_count(), 1);
    }

    #[test]
    fn chain_of_four_one_island() {
        let g = simple_chain(4, 5);
        assert_eq!(g.island_count(), 1);
    }

    #[test]
    fn two_separate_chains_two_islands() {
        let mut g = ConstraintGraph::new();
        for _ in 0..4 {
            g.add_body();
        }
        g.add_constraint(0, 1, 5);
        g.add_constraint(2, 3, 5);
        assert_eq!(g.island_count(), 2);
    }

    #[test]
    fn island_roots_correct_after_union() {
        let mut g = ConstraintGraph::new();
        for _ in 0..3 {
            g.add_body();
        }
        g.add_constraint(0, 1, 5);
        g.add_constraint(1, 2, 5);
        let islands = g.find_islands();
        // All three should share the same root
        assert_eq!(islands[0], islands[1]);
        assert_eq!(islands[1], islands[2]);
    }

    // ── mobility_analysis ─────────────────────────────────────────────────

    #[test]
    fn single_body_mobility_zero() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        // M = 6*(1-1) - 0 = 0
        assert_eq!(g.mobility_analysis(), 0);
    }

    #[test]
    fn two_bodies_no_constraints_mobility_six() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        // M = 6*(2-1) - 0 = 6
        assert_eq!(g.mobility_analysis(), 6);
    }

    #[test]
    fn revolute_joint_two_bodies_mobility_one() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        // Revolute removes 5 DOF (keeps 1)
        g.add_constraint(0, 1, 5);
        // M = 6*(2-1) - 5 = 1
        assert_eq!(g.mobility_analysis(), 1);
    }

    #[test]
    fn fully_fixed_two_bodies_mobility_zero() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 6);
        // M = 6*(2-1) - 6 = 0
        assert_eq!(g.mobility_analysis(), 0);
    }

    #[test]
    fn three_body_chain_four_revolutes_mobility_correct() {
        let mut g = ConstraintGraph::new();
        for _ in 0..3 {
            g.add_body();
        }
        // Two revolute joints (each removes 5)
        g.add_constraint(0, 1, 5);
        g.add_constraint(1, 2, 5);
        // M = 6*2 - 10 = 2
        assert_eq!(g.mobility_analysis(), 2);
    }

    #[test]
    fn over_constrained_negative_mobility() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        // Two fixed joints between same pair = over-constrained
        g.add_constraint(0, 1, 6);
        g.add_constraint(0, 1, 6);
        // M = 6*(2-1) - 12 = -6
        assert_eq!(g.mobility_analysis(), -6);
    }

    // ── kinematic_chain / closed_loop_detection ───────────────────────────

    #[test]
    fn acyclic_chain_topological_sort_succeeds() {
        let g = simple_chain(5, 5);
        let order = g.kinematic_chain();
        assert!(order.is_some());
        let order = order.unwrap();
        assert_eq!(order.len(), 5);
    }

    #[test]
    fn cyclic_graph_topological_sort_returns_none() {
        let mut g = ConstraintGraph::new();
        for _ in 0..3 {
            g.add_body();
        }
        g.add_constraint(0, 1, 5);
        g.add_constraint(1, 2, 5);
        g.add_constraint(2, 0, 5); // forms a cycle
        assert!(g.kinematic_chain().is_none());
    }

    #[test]
    fn no_closed_loop_in_tree() {
        let g = simple_chain(4, 5);
        assert!(!g.closed_loop_detection());
    }

    #[test]
    fn closed_loop_detected_in_cycle() {
        let mut g = ConstraintGraph::new();
        for _ in 0..4 {
            g.add_body();
        }
        g.add_constraint(0, 1, 5);
        g.add_constraint(1, 2, 5);
        g.add_constraint(2, 3, 5);
        g.add_constraint(3, 0, 5); // closed loop
        assert!(g.closed_loop_detection());
    }

    #[test]
    fn single_body_no_closed_loop() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        assert!(!g.closed_loop_detection());
    }

    // ── redundancy_count ──────────────────────────────────────────────────

    #[test]
    fn no_redundancy_for_positive_mobility() {
        let g = simple_chain(3, 5);
        // M = 6*2 - 10 = 2 (positive)
        assert_eq!(g.redundancy_count(), 0);
    }

    #[test]
    fn redundancy_for_negative_mobility() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 6);
        g.add_constraint(0, 1, 6);
        // M = -6, redundancy = 6
        assert_eq!(g.redundancy_count(), 6);
    }

    #[test]
    fn redundancy_zero_for_zero_mobility() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 6);
        // M = 0
        assert_eq!(g.redundancy_count(), 0);
    }

    // ── constraints accessor ──────────────────────────────────────────────

    #[test]
    fn constraints_slice_has_correct_length() {
        let g = simple_chain(5, 3);
        assert_eq!(g.constraints().len(), 4);
    }

    #[test]
    fn constraint_edge_fields_correct() {
        let mut g = ConstraintGraph::new();
        g.add_body();
        g.add_body();
        g.add_constraint(0, 1, 4);
        let edge = &g.constraints()[0];
        assert_eq!(edge.body_a, 0);
        assert_eq!(edge.body_b, 1);
        assert_eq!(edge.dof_removed, 4);
    }
}
