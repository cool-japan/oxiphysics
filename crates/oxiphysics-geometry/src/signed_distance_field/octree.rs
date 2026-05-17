// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive octree SDF storage.

// ─────────────────────────────────────────────────────────────────────────────
// Octree SDF Storage
// ─────────────────────────────────────────────────────────────────────────────

/// Node in an octree used to store SDF values adaptively.
#[derive(Debug, Clone)]
pub enum OctreeNode {
    /// Leaf node storing the SDF value at the cell centre.
    Leaf {
        /// SDF value at the centre.
        value: f64,
        /// Cell centre position.
        centre: [f64; 3],
        /// Half-size of the cell.
        half_size: f64,
    },
    /// Internal node with 8 children.
    Internal {
        /// Cell centre.
        centre: [f64; 3],
        /// Half-size.
        half_size: f64,
        /// Child nodes (octants).
        children: Box<[OctreeNode; 8]>,
    },
}

impl OctreeNode {
    /// Approximate SDF value at `p` by descending the octree.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        match self {
            OctreeNode::Leaf { value, .. } => *value,
            OctreeNode::Internal {
                centre,
                children,
                half_size: _,
            } => {
                let ix = if p[0] >= centre[0] { 1 } else { 0 };
                let iy = if p[1] >= centre[1] { 1 } else { 0 };
                let iz = if p[2] >= centre[2] { 1 } else { 0 };
                let child_idx = ix + 2 * iy + 4 * iz;
                children[child_idx].evaluate(p)
            }
        }
    }

    /// Half-size (radius) of this node's cell.
    pub fn half_size(&self) -> f64 {
        match self {
            OctreeNode::Leaf { half_size, .. } => *half_size,
            OctreeNode::Internal { half_size, .. } => *half_size,
        }
    }
}

/// Adaptive octree SDF storage.
#[derive(Debug, Clone)]
pub struct OctreeSdf {
    /// Root node.
    pub root: OctreeNode,
    /// Maximum subdivision depth.
    pub max_depth: usize,
    /// Subdivision threshold: refine if SDF value is below this.
    pub refine_threshold: f64,
}

impl OctreeSdf {
    /// Build an octree SDF from function `f` centred at `centre` with `half_size`.
    pub fn build<F>(
        f: &F,
        centre: [f64; 3],
        half_size: f64,
        max_depth: usize,
        refine_threshold: f64,
    ) -> Self
    where
        F: Fn([f64; 3]) -> f64,
    {
        let root = Self::build_node(f, centre, half_size, 0, max_depth, refine_threshold);
        Self {
            root,
            max_depth,
            refine_threshold,
        }
    }

    fn build_node<F>(
        f: &F,
        centre: [f64; 3],
        half_size: f64,
        depth: usize,
        max_depth: usize,
        threshold: f64,
    ) -> OctreeNode
    where
        F: Fn([f64; 3]) -> f64,
    {
        let value = f(centre);
        if depth >= max_depth || value.abs() > threshold {
            return OctreeNode::Leaf {
                value,
                centre,
                half_size,
            };
        }
        let hs = half_size * 0.5;
        let children: [OctreeNode; 8] = std::array::from_fn(|k| {
            let cx = centre[0] + if k & 1 != 0 { hs } else { -hs };
            let cy = centre[1] + if k & 2 != 0 { hs } else { -hs };
            let cz = centre[2] + if k & 4 != 0 { hs } else { -hs };
            Self::build_node(f, [cx, cy, cz], hs, depth + 1, max_depth, threshold)
        });
        OctreeNode::Internal {
            centre,
            half_size,
            children: Box::new(children),
        }
    }

    /// Evaluate the octree SDF at point `p`.
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        self.root.evaluate(p)
    }

    /// Count total number of leaf nodes.
    pub fn count_leaves(&self) -> usize {
        Self::count_leaves_node(&self.root)
    }

    fn count_leaves_node(node: &OctreeNode) -> usize {
        match node {
            OctreeNode::Leaf { .. } => 1,
            OctreeNode::Internal { children, .. } => {
                children.iter().map(Self::count_leaves_node).sum()
            }
        }
    }
}
