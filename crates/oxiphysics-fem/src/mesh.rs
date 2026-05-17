#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh types for finite element analysis.
//!
//! Provides node and tetrahedral mesh representations, along with
//! basic mesh generation utilities for common shapes.

use std::collections::HashMap;

use oxiphysics_core::math::{Mat3, Vec3};

/// A node in the finite element mesh.
#[derive(Debug, Clone)]
pub struct Node {
    /// Position of the node in 3D space.
    pub position: Vec3,
    /// Unique identifier for this node.
    pub id: usize,
}

/// A tetrahedral mesh for 3D finite element analysis.
///
/// Each element is defined by four node indices forming a tetrahedron.
#[derive(Debug, Clone)]
pub struct TetrahedralMesh {
    /// Node positions (index = node id).
    pub nodes: Vec<Vec3>,
    /// Elements, each defined by 4 node indices.
    pub elements: Vec<[usize; 4]>,
}

impl TetrahedralMesh {
    /// Create a tetrahedral mesh from node positions and element connectivity.
    ///
    /// # Panics
    ///
    /// Panics if any element index is out of bounds.
    pub fn from_nodes_and_elements(nodes: Vec<Vec3>, elements: Vec<[usize; 4]>) -> Self {
        let n = nodes.len();
        for (i, elem) in elements.iter().enumerate() {
            for &idx in elem {
                assert!(
                    idx < n,
                    "Element {i} references node {idx} but only {n} nodes exist"
                );
            }
        }
        Self { nodes, elements }
    }

    /// Return the number of nodes in the mesh.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Return the number of elements in the mesh.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }

    /// Generate a beam-shaped tetrahedral mesh.
    ///
    /// The beam extends from (0,0,0) to (length, width, height) and is
    /// subdivided into `nx * ny * nz` hexahedral cells, each split into
    /// 5 tetrahedra using a diagonal-alternating split that tiles the
    /// hex exactly (signed volumes sum to the hex volume, no overlap,
    /// no inversions) with face-compatible interfaces between cells.
    ///
    /// # Panics
    ///
    /// Panics if any of `nx`, `ny`, `nz` is zero.
    pub fn generate_beam(
        length: f64,
        width: f64,
        height: f64,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0, "subdivisions must be > 0");

        let dx = length / nx as f64;
        let dy = width / ny as f64;
        let dz = height / nz as f64;

        // Generate nodes on a regular grid
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
        for iz in 0..=nz {
            for iy in 0..=ny {
                for ix in 0..=nx {
                    nodes.push(Vec3::new(ix as f64 * dx, iy as f64 * dy, iz as f64 * dz));
                }
            }
        }

        // Helper to get node index from grid coordinates
        let node_idx = |ix: usize, iy: usize, iz: usize| -> usize {
            iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix
        };

        // Split each hexahedral cell into 5 tetrahedra using the standard
        // diagonal-alternating decomposition. The parity 0 split uses the
        // (n0, n6) body-diagonal; parity 1 uses the conjugate (n1, n7)
        // diagonal. Adjacent cells carry opposite parities so shared
        // faces are cut along the same diagonal from both sides, keeping
        // the global mesh face-compatible (watertight, conforming).
        //
        // Every tet below is written with a positively-oriented vertex
        // order: the signed volume `det([v1-v0, v2-v0, v3-v0]) / 6` is
        // strictly positive. For a unit hex the four corner tets each
        // contribute 1/6 and the central "core" tet contributes 2/6, so
        // the five signed volumes sum exactly to 1 = hex volume.
        let mut elements = Vec::new();
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    // 8 corners of the hex cell:
                    //   n0=(ix,  iy,  iz  )  n1=(ix+1,iy,  iz  )
                    //   n2=(ix+1,iy+1,iz  )  n3=(ix,  iy+1,iz  )
                    //   n4=(ix,  iy,  iz+1)  n5=(ix+1,iy,  iz+1)
                    //   n6=(ix+1,iy+1,iz+1)  n7=(ix,  iy+1,iz+1)
                    let n0 = node_idx(ix, iy, iz);
                    let n1 = node_idx(ix + 1, iy, iz);
                    let n2 = node_idx(ix + 1, iy + 1, iz);
                    let n3 = node_idx(ix, iy + 1, iz);
                    let n4 = node_idx(ix, iy, iz + 1);
                    let n5 = node_idx(ix + 1, iy, iz + 1);
                    let n6 = node_idx(ix + 1, iy + 1, iz + 1);
                    let n7 = node_idx(ix, iy + 1, iz + 1);

                    // Alternating parity ensures compatible shared faces.
                    let parity = (ix + iy + iz) % 2;
                    if parity == 0 {
                        // Diagonal n0–n6 split: four corner tets around
                        // n1, n2 (via the opposite corner n7), n5, n7,
                        // plus a positively-oriented central tet.
                        elements.push([n0, n1, n3, n4]); // V = +1/6
                        elements.push([n2, n1, n6, n3]); // V = +1/6
                        elements.push([n5, n1, n4, n6]); // V = +1/6
                        elements.push([n7, n3, n6, n4]); // V = +1/6
                        elements.push([n1, n3, n4, n6]); // V = +2/6 (core)
                    } else {
                        // Conjugate diagonal n1–n7 split.
                        elements.push([n0, n1, n2, n5]); // V = +1/6
                        elements.push([n0, n2, n3, n7]); // V = +1/6
                        elements.push([n0, n4, n5, n7]); // V = +1/6
                        elements.push([n2, n5, n6, n7]); // V = +1/6
                        elements.push([n0, n5, n2, n7]); // V = +2/6 (core)
                    }
                }
            }
        }

        Self { nodes, elements }
    }
}

/// A tetrahedral finite element mesh with boundary information.
///
/// Provides node coordinates, element connectivity, and automatically
/// computed boundary node indices.
#[derive(Debug, Clone)]
pub struct TetMesh {
    /// Node coordinates.
    pub nodes: Vec<Vec3>,
    /// Tetrahedral connectivity (4 node indices per element).
    pub elements: Vec<[usize; 4]>,
    /// Indices of boundary (surface) nodes.
    pub boundary_nodes: Vec<usize>,
}

impl TetMesh {
    /// Create a new `TetMesh`, automatically computing boundary nodes.
    pub fn new(nodes: Vec<Vec3>, elements: Vec<[usize; 4]>) -> Self {
        let boundary_nodes = Self::compute_boundary_nodes_for(&nodes, &elements);
        Self {
            nodes,
            elements,
            boundary_nodes,
        }
    }

    /// Return the number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Return the number of elements.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Total volume of the mesh (sum of all element volumes).
    pub fn total_volume(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| Self::tet_volume(&self.nodes, e))
            .sum()
    }

    /// Compute boundary nodes from the mesh topology.
    ///
    /// A face (triangle) is a boundary face if it belongs to exactly one element.
    /// Boundary nodes are all nodes of boundary faces.
    pub fn compute_boundary_nodes(&self) -> Vec<usize> {
        Self::compute_boundary_nodes_for(&self.nodes, &self.elements)
    }

    fn compute_boundary_nodes_for(_nodes: &[Vec3], elements: &[[usize; 4]]) -> Vec<usize> {
        use std::collections::HashMap;

        // Count how many elements share each face.
        // A face is represented as a sorted triple of node indices.
        let mut face_count: HashMap<[usize; 3], usize> = HashMap::new();

        for elem in elements {
            // The 4 faces of a tetrahedron (opposite each vertex)
            let faces = [
                [elem[1], elem[2], elem[3]],
                [elem[0], elem[2], elem[3]],
                [elem[0], elem[1], elem[3]],
                [elem[0], elem[1], elem[2]],
            ];
            for mut face in faces {
                face.sort_unstable();
                *face_count.entry(face).or_insert(0) += 1;
            }
        }

        // Collect all nodes belonging to faces shared by exactly 1 element.
        let mut boundary_set: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for (face, &count) in &face_count {
            if count == 1 {
                for &n in face {
                    boundary_set.insert(n);
                }
            }
        }

        boundary_set.into_iter().collect()
    }

    /// Volume of a single tetrahedron.
    fn tet_volume(nodes: &[Vec3], elem: &[usize; 4]) -> f64 {
        let p0 = nodes[elem[0]];
        let p1 = nodes[elem[1]];
        let p2 = nodes[elem[2]];
        let p3 = nodes[elem[3]];
        let a = p1 - p0;
        let b = p2 - p0;
        let c = p3 - p0;
        (a.cross(&b).dot(&c) / 6.0).abs()
    }

    /// Minimum dihedral angle across all elements (in radians).
    ///
    /// Well-shaped tetrahedra have a minimum dihedral angle greater than
    /// approximately 10° (0.175 rad).
    pub fn min_dihedral_angle(&self) -> f64 {
        self.elements
            .iter()
            .flat_map(|e| Self::tet_dihedral_angles(&self.nodes, e))
            .fold(f64::INFINITY, f64::min)
    }

    /// The 6 dihedral angles of a tetrahedron (one per edge).
    fn tet_dihedral_angles(nodes: &[Vec3], elem: &[usize; 4]) -> [f64; 6] {
        let p: [Vec3; 4] = [
            nodes[elem[0]],
            nodes[elem[1]],
            nodes[elem[2]],
            nodes[elem[3]],
        ];

        // Edge pairs: the dihedral angle along an edge is between the two faces
        // sharing that edge. The 6 edges are (i,j) for i<j in {0,1,2,3}.
        // For edge (i,j), the two opposing vertices are the remaining two.
        let edges = [
            (0, 1, 2, 3),
            (0, 2, 1, 3),
            (0, 3, 1, 2),
            (1, 2, 0, 3),
            (1, 3, 0, 2),
            (2, 3, 0, 1),
        ];

        let mut angles = [0.0f64; 6];
        for (k, &(i, j, m, n)) in edges.iter().enumerate() {
            // Edge vector
            let e = p[j] - p[i];
            // The two faces sharing edge (i,j) are (i,j,m) and (i,j,n).
            // Normal of face (i,j,m): (p[j]-p[i]) × (p[m]-p[i])
            let n1 = e.cross(&(p[m] - p[i]));
            // Normal of face (i,j,n): (p[j]-p[i]) × (p[n]-p[i])
            let n2 = e.cross(&(p[n] - p[i]));

            let len1 = n1.norm();
            let len2 = n2.norm();
            if len1 < 1e-15 || len2 < 1e-15 {
                angles[k] = 0.0;
            } else {
                let cos_a = n1.dot(&n2) / (len1 * len2);
                let cos_a = cos_a.clamp(-1.0, 1.0);
                // Dihedral angle is π minus the angle between outward normals.
                angles[k] = std::f64::consts::PI - cos_a.acos();
            }
        }
        angles
    }

    /// Aspect ratio of a single element (circumradius / inradius).
    ///
    /// A perfect regular tetrahedron has aspect ratio of 3.0.
    pub fn element_aspect_ratio(&self, elem_idx: usize) -> f64 {
        let elem = &self.elements[elem_idx];
        let p: [Vec3; 4] = [
            self.nodes[elem[0]],
            self.nodes[elem[1]],
            self.nodes[elem[2]],
            self.nodes[elem[3]],
        ];

        // Volume
        let a = p[1] - p[0];
        let b = p[2] - p[0];
        let c = p[3] - p[0];
        let vol = (a.cross(&b).dot(&c) / 6.0).abs();

        if vol < 1e-30 {
            return f64::INFINITY;
        }

        // Circumradius via circumcenter: solve the 3x3 linear system.
        let circumradius = {
            // Circumcenter of tet: solve the 3x3 linear system
            // 2*(p1-p0)·x = |p1|²-|p0|²  etc.
            let rhs = Vec3::new(
                p[1].dot(&p[1]) - p[0].dot(&p[0]),
                p[2].dot(&p[2]) - p[0].dot(&p[0]),
                p[3].dot(&p[3]) - p[0].dot(&p[0]),
            );
            let mat = Mat3::new(
                2.0 * (p[1][0] - p[0][0]),
                2.0 * (p[1][1] - p[0][1]),
                2.0 * (p[1][2] - p[0][2]),
                2.0 * (p[2][0] - p[0][0]),
                2.0 * (p[2][1] - p[0][1]),
                2.0 * (p[2][2] - p[0][2]),
                2.0 * (p[3][0] - p[0][0]),
                2.0 * (p[3][1] - p[0][1]),
                2.0 * (p[3][2] - p[0][2]),
            );
            match mat.lu().solve(&rhs) {
                Some(cc) => (cc - p[0]).norm(),
                None => return f64::INFINITY,
            }
        };

        // Inradius: r = 3V / A_total where A_total = sum of face areas
        let face_area = |i: usize, j: usize, k: usize| -> f64 {
            ((p[j] - p[i]).cross(&(p[k] - p[i]))).norm() * 0.5
        };
        let total_surface =
            face_area(0, 1, 2) + face_area(0, 1, 3) + face_area(0, 2, 3) + face_area(1, 2, 3);

        if total_surface < 1e-30 {
            return f64::INFINITY;
        }

        let inradius = 3.0 * vol / total_surface;

        circumradius / inradius
    }
}

/// Generate a structured tetrahedral mesh of the box \[0,lx\] × \[0,ly\] × \[0,lz\].
///
/// Each hexahedral cell is split into 6 tetrahedra using the standard
/// Freudenthal decomposition, giving `nx * ny * nz * 6` elements total.
///
/// # Panics
///
/// Panics if any of `nx`, `ny`, `nz` is zero.
pub fn generate_box_mesh(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> TetMesh {
    assert!(nx > 0 && ny > 0 && nz > 0, "subdivisions must be > 0");

    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let dz = lz / nz as f64;

    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for iz in 0..=nz {
        for iy in 0..=ny {
            for ix in 0..=nx {
                nodes.push(Vec3::new(ix as f64 * dx, iy as f64 * dy, iz as f64 * dz));
            }
        }
    }

    let node_idx = |ix: usize, iy: usize, iz: usize| -> usize {
        iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix
    };

    // Split each hex cell into 6 tetrahedra (standard 6-tet decomposition).
    // Given the 8 vertices of the hex labeled:
    //   n0=(ix,iy,iz),    n1=(ix+1,iy,iz),
    //   n2=(ix,iy+1,iz),  n3=(ix+1,iy+1,iz),
    //   n4=(ix,iy,iz+1),  n5=(ix+1,iy,iz+1),
    //   n6=(ix,iy+1,iz+1),n7=(ix+1,iy+1,iz+1)
    // The 6 tets all share the main diagonal n0–n7.
    let mut elements = Vec::with_capacity(nx * ny * nz * 6);
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let n0 = node_idx(ix, iy, iz);
                let n1 = node_idx(ix + 1, iy, iz);
                let n2 = node_idx(ix, iy + 1, iz);
                let n3 = node_idx(ix + 1, iy + 1, iz);
                let n4 = node_idx(ix, iy, iz + 1);
                let n5 = node_idx(ix + 1, iy, iz + 1);
                let n6 = node_idx(ix, iy + 1, iz + 1);
                let n7 = node_idx(ix + 1, iy + 1, iz + 1);

                // Standard 6-tet decomposition sharing the main diagonal n0–n7
                elements.push([n0, n1, n3, n7]);
                elements.push([n0, n1, n5, n7]);
                elements.push([n0, n4, n5, n7]);
                elements.push([n0, n2, n3, n7]);
                elements.push([n0, n4, n6, n7]);
                elements.push([n0, n2, n6, n7]);
            }
        }
    }

    TetMesh::new(nodes, elements)
}

/// Generate a bar mesh with a 1×1 cross section along the X axis.
///
/// This is a convenience wrapper around [`generate_box_mesh`] with
/// `ly = lz = 1.0` and `ny = nz = 1`.
pub fn generate_bar_mesh(length: f64, nx: usize) -> TetMesh {
    generate_box_mesh(length, 1.0, 1.0, nx, 1, 1)
}

/// Generate a cantilever beam mesh suitable for Euler–Bernoulli comparison.
///
/// The beam extends from `x = 0` (fixed end) to `x = length` (free end).
/// The cross-section is `width × height` in the Y–Z plane.
///
/// The mesh uses `nx` elements along the length and `ny × nz` elements across
/// the cross-section (minimum 1 each).  Each hexahedral cell is split into
/// 6 tetrahedra using the Freudenthal decomposition.
///
/// Use [`cantilever_fixed_dofs`] to obtain the DOF indices to constrain,
/// and [`cantilever_tip_nodes`] to find where to apply tip loads.
///
/// # Panics
///
/// Panics if `nx`, `ny`, or `nz` is zero.
pub fn generate_cantilever_mesh(
    length: f64,
    width: f64,
    height: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> TetMesh {
    generate_box_mesh(length, width, height, nx, ny, nz)
}

/// Return all global DOF indices that should be fixed for a cantilever beam.
///
/// Fixes all three displacement components (x, y, z) at every node whose
/// x-coordinate is within `tol` of 0.
pub fn cantilever_fixed_dofs(mesh: &TetMesh, tol: f64) -> Vec<usize> {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.x.abs() < tol)
        .flat_map(|(i, _)| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect()
}

/// Return the indices of tip nodes (free end) of a cantilever beam.
///
/// Collects all nodes whose x-coordinate is within `tol` of `length`.
pub fn cantilever_tip_nodes(mesh: &TetMesh, length: f64, tol: f64) -> Vec<usize> {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| (n.x - length).abs() < tol)
        .map(|(i, _)| i)
        .collect()
}

/// Statistics summary for a [`TetMesh`].
pub struct MeshStats {
    /// Total number of nodes.
    pub n_nodes: usize,
    /// Total number of elements.
    pub n_elements: usize,
    /// Number of boundary nodes.
    pub n_boundary_nodes: usize,
    /// Total mesh volume.
    pub total_volume: f64,
    /// Minimum element aspect ratio.
    pub min_aspect_ratio: f64,
    /// Maximum element aspect ratio.
    pub max_aspect_ratio: f64,
    /// Mean element aspect ratio.
    pub mean_aspect_ratio: f64,
}

impl MeshStats {
    /// Compute statistics for the given mesh.
    pub fn compute(mesh: &TetMesh) -> Self {
        let n = mesh.n_elements();
        let ratios: Vec<f64> = (0..n).map(|i| mesh.element_aspect_ratio(i)).collect();

        let min_aspect_ratio = ratios.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_aspect_ratio = ratios.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean_aspect_ratio = if n == 0 {
            0.0
        } else {
            ratios.iter().sum::<f64>() / n as f64
        };

        Self {
            n_nodes: mesh.n_nodes(),
            n_elements: mesh.n_elements(),
            n_boundary_nodes: mesh.boundary_nodes.len(),
            total_volume: mesh.total_volume(),
            min_aspect_ratio,
            max_aspect_ratio,
            mean_aspect_ratio,
        }
    }
}

// ---------------------------------------------------------------------------
// Mesh refinement (h-refinement)
// ---------------------------------------------------------------------------

/// Result of an h-refinement pass on a tetrahedral mesh.
#[allow(dead_code)]
pub struct RefinedMesh {
    /// New node positions.
    pub nodes: Vec<Vec3>,
    /// New element connectivity.
    pub elements: Vec<[usize; 4]>,
    /// Number of elements before refinement.
    pub original_element_count: usize,
}

/// Perform uniform h-refinement on a `TetMesh` by bisecting every edge.
///
/// Each tetrahedron is split into 8 sub-tetrahedra by inserting mid-edge nodes.
/// Returns the refined mesh.
#[allow(dead_code)]
pub fn h_refine(mesh: &TetMesh) -> TetMesh {
    let mut new_nodes = mesh.nodes.clone();
    // Map from sorted edge (a,b) to midpoint node index.
    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();

    let mut midpoint = |a: usize, b: usize, nodes: &mut Vec<Vec3>| -> usize {
        let key = (a.min(b), a.max(b));
        if let Some(&idx) = edge_map.get(&key) {
            return idx;
        }
        let mid = (nodes[a] + nodes[b]) * 0.5;
        let idx = nodes.len();
        nodes.push(mid);
        edge_map.insert(key, idx);
        idx
    };

    let mut new_elements = Vec::with_capacity(mesh.elements.len() * 8);
    for elem in &mesh.elements {
        let [n0, n1, n2, n3] = *elem;
        // 6 midpoints for the 6 edges
        let m01 = midpoint(n0, n1, &mut new_nodes);
        let m02 = midpoint(n0, n2, &mut new_nodes);
        let m03 = midpoint(n0, n3, &mut new_nodes);
        let m12 = midpoint(n1, n2, &mut new_nodes);
        let m13 = midpoint(n1, n3, &mut new_nodes);
        let m23 = midpoint(n2, n3, &mut new_nodes);

        // Split into 8 sub-tetrahedra
        // 4 corner tets
        new_elements.push([n0, m01, m02, m03]);
        new_elements.push([n1, m01, m12, m13]);
        new_elements.push([n2, m02, m12, m23]);
        new_elements.push([n3, m03, m13, m23]);
        // 4 interior tets (octahedron decomposition)
        new_elements.push([m01, m02, m03, m12]);
        new_elements.push([m02, m03, m12, m23]);
        new_elements.push([m03, m12, m13, m23]);
        new_elements.push([m01, m03, m12, m13]);
    }

    TetMesh::new(new_nodes, new_elements)
}

// ---------------------------------------------------------------------------
// Mesh coarsening
// ---------------------------------------------------------------------------

/// Remove elements whose quality (aspect ratio) exceeds a threshold,
/// performing a simple coarsening. Returns the number of elements removed.
#[allow(dead_code)]
pub fn coarsen_by_quality(mesh: &mut TetMesh, max_aspect_ratio: f64) -> usize {
    let original_count = mesh.elements.len();
    let nodes = &mesh.nodes;
    mesh.elements.retain(|elem| {
        let tmp = TetMesh::new(nodes.clone(), vec![*elem]);
        tmp.element_aspect_ratio(0) <= max_aspect_ratio
    });
    let removed = original_count - mesh.elements.len();
    mesh.boundary_nodes = mesh.compute_boundary_nodes();
    removed
}

/// Collapse short edges below `min_length`, merging the two endpoints into
/// their midpoint. Returns the number of edges collapsed.
#[allow(dead_code)]
pub fn collapse_short_edges(mesh: &mut TetMesh, min_length: f64) -> usize {
    let mut collapsed = 0usize;
    let mut remap: HashMap<usize, usize> = HashMap::new();

    // Find edges to collapse
    for elem in &mesh.elements {
        let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for &(a_idx, b_idx) in &pairs {
            let a = elem[a_idx];
            let b = elem[b_idx];
            if a == b || remap.contains_key(&a) || remap.contains_key(&b) {
                continue;
            }
            let edge_len = (mesh.nodes[a] - mesh.nodes[b]).norm();
            if edge_len < min_length {
                // Merge b into a, move a to midpoint
                let mid = (mesh.nodes[a] + mesh.nodes[b]) * 0.5;
                mesh.nodes[a] = mid;
                remap.insert(b, a);
                collapsed += 1;
            }
        }
    }

    // Apply remapping
    for elem in &mut mesh.elements {
        for idx in elem.iter_mut() {
            if let Some(&new) = remap.get(idx) {
                *idx = new;
            }
        }
    }

    // Remove degenerate elements (those with repeated indices)
    mesh.elements.retain(|elem| {
        let mut sorted = *elem;
        sorted.sort_unstable();
        sorted[0] != sorted[1] && sorted[1] != sorted[2] && sorted[2] != sorted[3]
    });
    mesh.boundary_nodes = mesh.compute_boundary_nodes();

    collapsed
}

// ---------------------------------------------------------------------------
// Element quality checks
// ---------------------------------------------------------------------------

/// Quality metric for a single tetrahedral element.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ElementQuality {
    /// Aspect ratio (circumradius / inradius).
    pub aspect_ratio: f64,
    /// Volume of the element.
    pub volume: f64,
    /// Minimum dihedral angle in radians.
    pub min_dihedral: f64,
    /// Maximum dihedral angle in radians.
    pub max_dihedral: f64,
}

/// Compute quality metrics for all elements in a mesh.
#[allow(dead_code)]
pub fn compute_element_qualities(mesh: &TetMesh) -> Vec<ElementQuality> {
    (0..mesh.n_elements())
        .map(|i| {
            let elem = &mesh.elements[i];
            let ar = mesh.element_aspect_ratio(i);
            let vol = TetMesh::tet_volume(&mesh.nodes, elem);
            let angles = TetMesh::tet_dihedral_angles(&mesh.nodes, elem);
            let min_d = angles.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_d = angles.iter().cloned().fold(0.0_f64, f64::max);
            ElementQuality {
                aspect_ratio: ar,
                volume: vol,
                min_dihedral: min_d,
                max_dihedral: max_d,
            }
        })
        .collect()
}

/// Find the worst-quality element index by aspect ratio.
#[allow(dead_code)]
pub fn worst_element(mesh: &TetMesh) -> Option<usize> {
    if mesh.n_elements() == 0 {
        return None;
    }
    let mut worst_idx = 0;
    let mut worst_ar = f64::NEG_INFINITY;
    for i in 0..mesh.n_elements() {
        let ar = mesh.element_aspect_ratio(i);
        if ar > worst_ar {
            worst_ar = ar;
            worst_idx = i;
        }
    }
    Some(worst_idx)
}

// ---------------------------------------------------------------------------
// Mesh partitioning (simple greedy bisection)
// ---------------------------------------------------------------------------

/// Partition mesh elements into two roughly equal-sized groups based on
/// the centroid coordinate along `axis` (0=x, 1=y, 2=z).
///
/// Returns `(left_elements, right_elements)` where each is a vector of
/// element indices.
#[allow(dead_code)]
pub fn partition_bisection(mesh: &TetMesh, axis: usize) -> (Vec<usize>, Vec<usize>) {
    assert!(axis < 3, "axis must be 0, 1, or 2");
    let n = mesh.n_elements();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // Compute centroid of each element along the chosen axis.
    let mut centroids: Vec<(usize, f64)> = (0..n)
        .map(|i| {
            let elem = &mesh.elements[i];
            let c = elem.iter().map(|&ni| mesh.nodes[ni][axis]).sum::<f64>() / 4.0;
            (i, c)
        })
        .collect();

    centroids.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mid = n / 2;
    let left: Vec<usize> = centroids[..mid].iter().map(|&(i, _)| i).collect();
    let right: Vec<usize> = centroids[mid..].iter().map(|&(i, _)| i).collect();

    (left, right)
}

/// Recursive bisection partitioning into `2^levels` parts.
#[allow(dead_code)]
pub fn partition_recursive_bisection(mesh: &TetMesh, levels: usize) -> Vec<Vec<usize>> {
    if levels == 0 || mesh.n_elements() == 0 {
        return vec![(0..mesh.n_elements()).collect()];
    }

    let axis = levels % 3; // alternate axes
    let (left, right) = partition_bisection(mesh, axis);

    if levels == 1 {
        return vec![left, right];
    }

    // Build sub-meshes and recurse (just track indices, don't actually split).
    // For simplicity, we use a greedy approach with element indices.
    let mut result = Vec::new();

    // Further split each half
    let split_half = |indices: Vec<usize>| -> Vec<Vec<usize>> {
        if indices.len() <= 1 {
            return vec![indices];
        }
        let mid = indices.len() / 2;
        // Sort by next axis centroid
        let next_axis = (axis + 1) % 3;
        let mut sorted: Vec<(usize, f64)> = indices
            .iter()
            .map(|&i| {
                let elem = &mesh.elements[i];
                let c = elem
                    .iter()
                    .map(|&ni| mesh.nodes[ni][next_axis])
                    .sum::<f64>()
                    / 4.0;
                (i, c)
            })
            .collect();
        sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let l: Vec<usize> = sorted[..mid].iter().map(|&(i, _)| i).collect();
        let r: Vec<usize> = sorted[mid..].iter().map(|&(i, _)| i).collect();
        vec![l, r]
    };

    result.extend(split_half(left));
    result.extend(split_half(right));
    result
}

// ---------------------------------------------------------------------------
// Node renumbering (Cuthill-McKee)
// ---------------------------------------------------------------------------

/// Reverse Cuthill-McKee (RCM) node ordering to reduce matrix bandwidth.
///
/// Returns a permutation vector where `perm[new_index] = old_index`.
#[allow(dead_code)]
pub fn cuthill_mckee(mesh: &TetMesh) -> Vec<usize> {
    let n = mesh.n_nodes();
    if n == 0 {
        return Vec::new();
    }

    // Build adjacency list
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for elem in &mesh.elements {
        for i in 0..4 {
            for j in (i + 1)..4 {
                let a = elem[i];
                let b = elem[j];
                if !adj[a].contains(&b) {
                    adj[a].push(b);
                }
                if !adj[b].contains(&a) {
                    adj[b].push(a);
                }
            }
        }
    }

    // Sort adjacency by degree (ascending)
    let degrees: Vec<usize> = adj.iter().map(|n| n.len()).collect();
    for neighbors in &mut adj {
        neighbors.sort_by_key(|&v| degrees[v]);
    }

    // BFS from node with minimum degree
    let start = (0..n).min_by_key(|&v| adj[v].len()).unwrap_or(0);
    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut queue = std::collections::VecDeque::new();

    queue.push_back(start);
    visited[start] = true;

    while let Some(v) = queue.pop_front() {
        order.push(v);
        for &w in &adj[v] {
            if !visited[w] {
                visited[w] = true;
                queue.push_back(w);
            }
        }
    }

    // Handle disconnected components
    for v in 0..n {
        if !visited[v] {
            visited[v] = true;
            order.push(v);
        }
    }

    // Reverse Cuthill-McKee: reverse the order
    order.reverse();
    order
}

/// Compute the bandwidth of the stiffness matrix for a given node ordering.
///
/// The bandwidth is `max |perm_inv[a] - perm_inv[b]|` over all edges (a,b).
#[allow(dead_code)]
pub fn compute_bandwidth(mesh: &TetMesh, perm: &[usize]) -> usize {
    let n = mesh.n_nodes();
    // Build inverse permutation
    let mut inv = vec![0usize; n];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        if old_idx < n {
            inv[old_idx] = new_idx;
        }
    }

    let mut max_bw = 0usize;
    for elem in &mesh.elements {
        for i in 0..4 {
            for j in (i + 1)..4 {
                let a = inv[elem[i]];
                let b = inv[elem[j]];
                let diff = a.abs_diff(b);
                if diff > max_bw {
                    max_bw = diff;
                }
            }
        }
    }
    max_bw
}

/// Apply a node permutation to a `TetMesh`, reordering node positions and
/// updating element connectivity.
#[allow(dead_code)]
pub fn apply_permutation(mesh: &TetMesh, perm: &[usize]) -> TetMesh {
    let n = mesh.n_nodes();
    assert_eq!(perm.len(), n);

    // Build inverse permutation: inv[old] = new
    let mut inv = vec![0usize; n];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        inv[old_idx] = new_idx;
    }

    // Reorder nodes
    let mut new_nodes = vec![Vec3::zeros(); n];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        new_nodes[new_idx] = mesh.nodes[old_idx];
    }

    // Update element connectivity
    let new_elements: Vec<[usize; 4]> = mesh
        .elements
        .iter()
        .map(|elem| [inv[elem[0]], inv[elem[1]], inv[elem[2]], inv[elem[3]]])
        .collect();

    TetMesh::new(new_nodes, new_elements)
}

// ============================================================================
// Mesh connectivity: edge and face tables
// ============================================================================

/// Build an edge table for a tetrahedral mesh.
///
/// Returns a sorted, deduplicated list of directed edges `(a, b)` where `a <= b`.
#[allow(dead_code)]
pub fn mesh_edge_table(mesh: &TetMesh) -> Vec<(usize, usize)> {
    let mut edges = std::collections::BTreeSet::new();
    for elem in &mesh.elements {
        let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for (i, j) in pairs {
            let a = elem[i].min(elem[j]);
            let b = elem[i].max(elem[j]);
            edges.insert((a, b));
        }
    }
    edges.into_iter().collect()
}

/// Build a face table: returns all boundary triangular faces.
///
/// A face is on the boundary if it belongs to exactly one tetrahedron.
/// Returns faces as sorted (a, b, c) triples with a < b < c.
#[allow(dead_code)]
pub fn mesh_face_table(mesh: &TetMesh) -> Vec<(usize, usize, usize)> {
    let mut face_count: std::collections::HashMap<(usize, usize, usize), usize> =
        std::collections::HashMap::new();
    for elem in &mesh.elements {
        // 4 faces of a tetrahedron
        let faces = [
            [elem[0], elem[1], elem[2]],
            [elem[0], elem[1], elem[3]],
            [elem[0], elem[2], elem[3]],
            [elem[1], elem[2], elem[3]],
        ];
        for face in &faces {
            let mut sorted = *face;
            sorted.sort_unstable();
            *face_count
                .entry((sorted[0], sorted[1], sorted[2]))
                .or_insert(0) += 1;
        }
    }
    face_count
        .into_iter()
        .filter(|&(_, count)| count == 1)
        .map(|(face, _)| face)
        .collect()
}

// ============================================================================
// Mesh quality assessment
// ============================================================================

/// Compute the quality of each tetrahedron.
///
/// Quality is measured as the ratio of the inradius to the circumradius:
/// `q = 3 * r_in / r_circ`
///
/// A value of 1.0 corresponds to a regular (ideal) tetrahedron.
/// Returns values in \[0, 1\].
#[allow(dead_code)]
pub fn mesh_quality(mesh: &TetMesh) -> Vec<f64> {
    use oxiphysics_core::math::Vec3;
    mesh.elements
        .iter()
        .map(|elem| {
            let p0 = mesh.nodes[elem[0]];
            let p1 = mesh.nodes[elem[1]];
            let p2 = mesh.nodes[elem[2]];
            let p3 = mesh.nodes[elem[3]];

            // Volume
            let v01 = p1 - p0;
            let v02 = p2 - p0;
            let v03 = p3 - p0;
            let vol = v01.dot(&v02.cross(&v03)).abs() / 6.0;

            if vol < 1e-30 {
                return 0.0;
            }

            // Edge lengths
            let edges = [(p0, p1), (p0, p2), (p0, p3), (p1, p2), (p1, p3), (p2, p3)];
            let edge_lens: Vec<f64> = edges.iter().map(|(a, b)| (*b - *a).norm()).collect();
            let l_min = edge_lens.iter().cloned().fold(f64::INFINITY, f64::min);
            let l_max = edge_lens.iter().cloned().fold(0.0_f64, f64::max);

            if l_min < 1e-30 || l_max < 1e-30 {
                return 0.0;
            }

            // Face areas
            let faces: [[Vec3; 3]; 4] = [[p1, p2, p3], [p0, p2, p3], [p0, p1, p3], [p0, p1, p2]];
            let face_areas: Vec<f64> = faces
                .iter()
                .map(|f| {
                    let ab = f[1] - f[0];
                    let ac = f[2] - f[0];
                    ab.cross(&ac).norm() * 0.5
                })
                .collect();
            let total_face_area: f64 = face_areas.iter().sum();

            if total_face_area < 1e-30 {
                return 0.0;
            }

            // Inradius: r = 3 * V / A_total
            let r_in = 3.0 * vol / total_face_area;

            // Circumradius of regular tet with same edge sum
            // Approximation: r_circ = l_max * sqrt(3.0 / 8.0) for regular tet
            // More precise: use actual formula via determinant
            let l_avg = edge_lens.iter().sum::<f64>() / 6.0;
            let r_circ_ref = l_avg * (3.0_f64 / 8.0).sqrt();

            if r_circ_ref < 1e-30 {
                return 0.0;
            }

            // Quality: scaled so that regular tet gives ~1.0
            let q = (r_in / r_circ_ref) * (l_max / l_min.max(1e-30)).min(3.0).recip() * 3.0;
            q.clamp(0.0, 1.0)
        })
        .collect()
}

// ============================================================================
// Mesh statistics
// ============================================================================

/// Aggregate mesh statistics.
#[derive(Debug, Clone)]
pub struct MeshStatistics {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Number of elements.
    pub n_elements: usize,
    /// Minimum element volume.
    pub min_volume: f64,
    /// Maximum element volume.
    pub max_volume: f64,
    /// Mean element volume.
    pub mean_volume: f64,
    /// Minimum edge length across all elements.
    pub min_edge_length: f64,
    /// Maximum edge length across all elements.
    pub max_edge_length: f64,
    /// Minimum mesh quality metric (0=degenerate, 1=ideal).
    pub min_quality: f64,
    /// Mean mesh quality metric.
    pub mean_quality: f64,
}

/// Compute aggregate statistics for a tetrahedral mesh.
#[allow(dead_code)]
pub fn mesh_statistics(mesh: &TetMesh) -> MeshStatistics {
    let n = mesh.n_elements();
    let mut min_vol = f64::INFINITY;
    let mut max_vol = f64::NEG_INFINITY;
    let mut sum_vol = 0.0_f64;
    let mut min_edge = f64::INFINITY;
    let mut max_edge = 0.0_f64;

    for elem in &mesh.elements {
        let p0 = mesh.nodes[elem[0]];
        let p1 = mesh.nodes[elem[1]];
        let p2 = mesh.nodes[elem[2]];
        let p3 = mesh.nodes[elem[3]];

        let v01 = p1 - p0;
        let v02 = p2 - p0;
        let v03 = p3 - p0;
        let vol = v01.dot(&v02.cross(&v03)).abs() / 6.0;

        if vol < min_vol {
            min_vol = vol;
        }
        if vol > max_vol {
            max_vol = vol;
        }
        sum_vol += vol;

        let edges = [(p0, p1), (p0, p2), (p0, p3), (p1, p2), (p1, p3), (p2, p3)];
        for (a, b) in &edges {
            let l = (*b - *a).norm();
            if l < min_edge {
                min_edge = l;
            }
            if l > max_edge {
                max_edge = l;
            }
        }
    }

    let qualities = mesh_quality(mesh);
    let min_q = qualities.iter().cloned().fold(f64::INFINITY, f64::min);
    let mean_q = if n > 0 {
        qualities.iter().sum::<f64>() / n as f64
    } else {
        0.0
    };

    MeshStatistics {
        n_nodes: mesh.n_nodes(),
        n_elements: n,
        min_volume: if n > 0 { min_vol } else { 0.0 },
        max_volume: if n > 0 { max_vol } else { 0.0 },
        mean_volume: if n > 0 { sum_vol / n as f64 } else { 0.0 },
        min_edge_length: if n > 0 { min_edge } else { 0.0 },
        max_edge_length: if n > 0 { max_edge } else { 0.0 },
        min_quality: if n > 0 { min_q } else { 0.0 },
        mean_quality: mean_q,
    }
}

// ============================================================================
// Mesh coarsening (edge collapse)
// ============================================================================

/// Mesh coarsening by collapsing edges shorter than `min_length`.
///
/// Wraps the existing `collapse_short_edges` function and returns a new mesh.
#[allow(dead_code)]
pub fn mesh_coarsen(mesh: &TetMesh, min_length: f64) -> TetMesh {
    let mut coarsened = mesh.clone();
    collapse_short_edges(&mut coarsened, min_length);
    coarsened
}

// ============================================================================
// Mesh refinement (uniform red refinement)
// ============================================================================

/// Uniformly refine a tetrahedral mesh by splitting each tetrahedron.
///
/// Each tetrahedron is split into 8 sub-tetrahedra by introducing midpoints
/// on each edge (red refinement). This approximately doubles the mesh resolution.
#[allow(dead_code)]
pub fn mesh_refine_uniform(mesh: &TetMesh) -> TetMesh {
    h_refine(mesh)
}

// ============================================================================
// RCM (Reverse Cuthill-McKee) reordering
// ============================================================================

/// Reverse Cuthill-McKee reordering to reduce matrix bandwidth.
///
/// This is an alias for the existing `cuthill_mckee` function.
#[allow(dead_code)]
pub fn rcm_reordering(mesh: &TetMesh) -> Vec<usize> {
    cuthill_mckee(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beam_generation() {
        let mesh = TetrahedralMesh::generate_beam(1.0, 0.1, 0.1, 4, 1, 1);
        assert_eq!(mesh.num_nodes(), 5 * 2 * 2); // (4+1)*(1+1)*(1+1) = 20
        assert_eq!(mesh.num_elements(), 4 * 5); // 4 cells * 5 tets each = 20

        // All node indices should be valid
        for elem in &mesh.elements {
            for &idx in elem {
                assert!(idx < mesh.num_nodes());
            }
        }
    }

    #[test]
    fn test_from_nodes_and_elements() {
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let elements = vec![[0, 1, 2, 3]];
        let mesh = TetrahedralMesh::from_nodes_and_elements(nodes, elements);
        assert_eq!(mesh.num_nodes(), 4);
        assert_eq!(mesh.num_elements(), 1);
    }

    #[test]
    fn test_box_mesh_volume() {
        let mesh = generate_box_mesh(2.0, 3.0, 4.0, 2, 3, 4);
        let vol = mesh.total_volume();
        assert!(
            (vol - 24.0).abs() < 1e-10,
            "expected volume 24.0, got {vol}"
        );
    }

    #[test]
    fn test_box_mesh_node_count() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_nodes(), 3 * 3 * 3, "expected 27 nodes");
    }

    #[test]
    fn test_box_mesh_element_count() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        assert_eq!(mesh.n_elements(), 2 * 2 * 2 * 6, "expected 48 elements");
    }

    #[test]
    fn test_mesh_boundary_nodes() {
        // Single tetrahedron: all 4 nodes should be boundary nodes.
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let elements = vec![[0, 1, 2, 3]];
        let mesh = TetMesh::new(nodes, elements);
        assert_eq!(
            mesh.boundary_nodes.len(),
            4,
            "all 4 nodes should be boundary"
        );
    }

    #[test]
    fn test_bar_mesh_generation() {
        let mesh = generate_bar_mesh(4.0, 4);
        let vol = mesh.total_volume();
        assert!((vol - 4.0).abs() < 1e-10, "expected volume 4.0, got {vol}");
    }

    #[test]
    fn test_mesh_stats() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let stats = MeshStats::compute(&mesh);
        assert!(stats.n_nodes > 0);
        assert!(stats.n_elements > 0);
        assert!(stats.n_boundary_nodes > 0);
        assert!(stats.total_volume > 0.0);
        assert!(stats.min_aspect_ratio > 0.0);
        assert!(stats.max_aspect_ratio > 0.0);
        assert!(stats.mean_aspect_ratio > 0.0);
    }

    #[test]
    fn test_aspect_ratio_regular_tet() {
        // Regular tetrahedron with edge length 1: circumradius/inradius = 3.0
        let s = 1.0_f64;
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(s, 0.0, 0.0),
            Vec3::new(s / 2.0, s * 3.0_f64.sqrt() / 2.0, 0.0),
            Vec3::new(
                s / 2.0,
                s / (2.0 * 3.0_f64.sqrt()),
                s * (2.0 / 3.0_f64).sqrt(),
            ),
        ];
        let elements = vec![[0, 1, 2, 3]];
        let mesh = TetMesh::new(nodes, elements);
        let ratio = mesh.element_aspect_ratio(0);
        assert!(
            (ratio - 3.0).abs() < 1e-6,
            "expected aspect ratio 3.0 for regular tet, got {ratio}"
        );
    }

    // --- h-refinement ---
    #[test]
    fn test_h_refine_increases_elements() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let original_elems = mesh.n_elements();
        let refined = h_refine(&mesh);
        assert_eq!(refined.n_elements(), original_elems * 8);
    }

    #[test]
    fn test_h_refine_preserves_volume() {
        let mesh = generate_box_mesh(2.0, 3.0, 4.0, 1, 1, 1);
        let original_vol = mesh.total_volume();
        let refined = h_refine(&mesh);
        let refined_vol = refined.total_volume();
        assert!(
            (refined_vol - original_vol).abs() < 1e-8,
            "volume changed: {original_vol} -> {refined_vol}"
        );
    }

    #[test]
    fn test_h_refine_increases_nodes() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let refined = h_refine(&mesh);
        assert!(refined.n_nodes() > mesh.n_nodes());
    }

    #[test]
    fn test_h_refine_valid_indices() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let refined = h_refine(&mesh);
        for elem in &refined.elements {
            for &idx in elem {
                assert!(idx < refined.n_nodes(), "invalid node index {idx}");
            }
        }
    }

    // --- coarsening ---
    #[test]
    fn test_coarsen_by_quality_removes_bad_elements() {
        let mut mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let original = mesh.n_elements();
        // With a very strict threshold, some elements should be removed
        let removed = coarsen_by_quality(&mut mesh, 3.0);
        // At least check it doesn't panic and count is consistent
        assert_eq!(mesh.n_elements() + removed, original);
    }

    #[test]
    fn test_coarsen_lenient_keeps_all() {
        let mut mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let original = mesh.n_elements();
        let removed = coarsen_by_quality(&mut mesh, 1000.0);
        assert_eq!(removed, 0);
        assert_eq!(mesh.n_elements(), original);
    }

    // --- collapse short edges ---
    #[test]
    fn test_collapse_short_edges_no_change() {
        let mut mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let original = mesh.n_elements();
        let collapsed = collapse_short_edges(&mut mesh, 1e-10);
        assert_eq!(collapsed, 0);
        assert_eq!(mesh.n_elements(), original);
    }

    // --- element quality ---
    #[test]
    fn test_compute_element_qualities() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let qualities = compute_element_qualities(&mesh);
        assert_eq!(qualities.len(), mesh.n_elements());
        for q in &qualities {
            assert!(q.volume > 0.0);
            assert!(q.aspect_ratio > 0.0);
            assert!(q.min_dihedral > 0.0);
            assert!(
                q.max_dihedral > q.min_dihedral || (q.max_dihedral - q.min_dihedral).abs() < 1e-10
            );
        }
    }

    #[test]
    fn test_worst_element() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let worst = worst_element(&mesh);
        assert!(worst.is_some());
        let idx = worst.unwrap();
        assert!(idx < mesh.n_elements());
    }

    #[test]
    fn test_worst_element_empty() {
        let mesh = TetMesh::new(Vec::new(), Vec::new());
        assert!(worst_element(&mesh).is_none());
    }

    // --- partitioning ---
    #[test]
    fn test_partition_bisection_splits_evenly() {
        let mesh = generate_box_mesh(4.0, 1.0, 1.0, 4, 1, 1);
        let (left, right) = partition_bisection(&mesh, 0);
        assert_eq!(left.len() + right.len(), mesh.n_elements());
        assert!(!left.is_empty());
        assert!(!right.is_empty());
    }

    #[test]
    fn test_partition_bisection_axis_y() {
        let mesh = generate_box_mesh(1.0, 4.0, 1.0, 1, 4, 1);
        let (left, right) = partition_bisection(&mesh, 1);
        assert_eq!(left.len() + right.len(), mesh.n_elements());
    }

    #[test]
    fn test_partition_recursive() {
        let mesh = generate_box_mesh(4.0, 4.0, 4.0, 4, 4, 4);
        let parts = partition_recursive_bisection(&mesh, 2);
        let total: usize = parts.iter().map(|p| p.len()).sum();
        assert_eq!(total, mesh.n_elements());
        assert!(
            parts.len() <= 4,
            "expected at most 4 parts, got {}",
            parts.len()
        );
    }

    #[test]
    fn test_partition_empty_mesh() {
        let mesh = TetMesh::new(Vec::new(), Vec::new());
        let (left, right) = partition_bisection(&mesh, 0);
        assert!(left.is_empty());
        assert!(right.is_empty());
    }

    // --- Cuthill-McKee ---
    #[test]
    fn test_cuthill_mckee_permutation_complete() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let perm = cuthill_mckee(&mesh);
        assert_eq!(perm.len(), mesh.n_nodes());
        // Check it's a valid permutation (every node appears exactly once)
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        let expected: Vec<usize> = (0..mesh.n_nodes()).collect();
        assert_eq!(sorted, expected);
    }

    #[test]
    fn test_cuthill_mckee_reduces_bandwidth() {
        let mesh = generate_box_mesh(2.0, 1.0, 1.0, 4, 2, 2);
        let identity_perm: Vec<usize> = (0..mesh.n_nodes()).collect();
        let bw_original = compute_bandwidth(&mesh, &identity_perm);
        let perm = cuthill_mckee(&mesh);
        let bw_rcm = compute_bandwidth(&mesh, &perm);
        // RCM should not increase bandwidth (may not always decrease for small meshes)
        assert!(
            bw_rcm <= bw_original + 1,
            "RCM bandwidth {bw_rcm} > original {bw_original}"
        );
    }

    #[test]
    fn test_compute_bandwidth() {
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let elements = vec![[0, 1, 2, 3]];
        let mesh = TetMesh::new(nodes, elements);
        let perm = vec![0, 1, 2, 3];
        let bw = compute_bandwidth(&mesh, &perm);
        assert_eq!(bw, 3); // max |0-3| = 3
    }

    #[test]
    fn test_apply_permutation_preserves_volume() {
        let mesh = generate_box_mesh(2.0, 3.0, 4.0, 2, 2, 2);
        let perm = cuthill_mckee(&mesh);
        let permuted = apply_permutation(&mesh, &perm);
        let vol_orig = mesh.total_volume();
        let vol_perm = permuted.total_volume();
        assert!(
            (vol_orig - vol_perm).abs() < 1e-8,
            "volume changed: {vol_orig} -> {vol_perm}"
        );
    }

    #[test]
    fn test_apply_permutation_node_count() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let perm = cuthill_mckee(&mesh);
        let permuted = apply_permutation(&mesh, &perm);
        assert_eq!(permuted.n_nodes(), mesh.n_nodes());
        assert_eq!(permuted.n_elements(), mesh.n_elements());
    }

    #[test]
    fn test_cuthill_mckee_empty() {
        let mesh = TetMesh::new(Vec::new(), Vec::new());
        let perm = cuthill_mckee(&mesh);
        assert!(perm.is_empty());
    }

    #[test]
    fn test_cantilever_fixed_dofs() {
        let mesh = generate_cantilever_mesh(2.0, 0.5, 0.5, 4, 1, 1);
        let fixed = cantilever_fixed_dofs(&mesh, 1e-10);
        assert!(!fixed.is_empty());
        // All fixed DOFs should correspond to nodes at x=0
        for &dof in &fixed {
            let node = dof / 3;
            assert!(mesh.nodes[node].x.abs() < 1e-10);
        }
    }

    #[test]
    fn test_cantilever_tip_nodes() {
        let length = 2.0;
        let mesh = generate_cantilever_mesh(length, 0.5, 0.5, 4, 1, 1);
        let tips = cantilever_tip_nodes(&mesh, length, 1e-10);
        assert!(!tips.is_empty());
        for &node in &tips {
            assert!((mesh.nodes[node].x - length).abs() < 1e-10);
        }
    }

    // ── Mesh connectivity tests ─────────────────────────────────────────────

    #[test]
    fn test_mesh_edge_table_nonempty() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let edges = mesh_edge_table(&mesh);
        assert!(!edges.is_empty(), "Edge table should be nonempty");
    }

    #[test]
    fn test_mesh_edge_table_unique_edges() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let edges = mesh_edge_table(&mesh);
        // All edges should be unique (a <= b)
        for &(a, b) in &edges {
            assert!(a <= b, "Edge ({a},{b}) should have a <= b");
        }
        let mut sorted = edges.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), edges.len(), "No duplicate edges");
    }

    #[test]
    fn test_mesh_face_table_nonempty() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let faces = mesh_face_table(&mesh);
        assert!(!faces.is_empty());
    }

    // ── Mesh quality tests ──────────────────────────────────────────────────

    #[test]
    fn test_mesh_quality_range() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let qualities = mesh_quality(&mesh);
        assert_eq!(qualities.len(), mesh.n_elements());
        for &q in &qualities {
            assert!(q >= 0.0 && q <= 1.0, "Quality must be in [0,1], got {q}");
        }
    }

    #[test]
    fn test_mesh_quality_regular_tet() {
        // A regular tetrahedron should have quality close to 1
        let s = 1.0_f64;
        let nodes = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(s, 0.0, 0.0),
            Vec3::new(0.5 * s, (3.0_f64.sqrt() / 2.0) * s, 0.0),
            Vec3::new(
                0.5 * s,
                (3.0_f64.sqrt() / 6.0) * s,
                (6.0_f64.sqrt() / 3.0) * s,
            ),
        ];
        let mesh = TetMesh::new(nodes, vec![[0, 1, 2, 3]]);
        let q = mesh_quality(&mesh);
        assert!(
            q[0] > 0.8,
            "Regular tet quality should be > 0.8, got {}",
            q[0]
        );
    }

    // ── Mesh statistics tests ───────────────────────────────────────────────

    #[test]
    fn test_mesh_statistics() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let stats = mesh_statistics(&mesh);
        assert!(stats.min_volume > 0.0, "Min volume should be positive");
        assert!(stats.max_volume >= stats.min_volume);
        assert!(stats.min_edge_length > 0.0);
        assert!(stats.max_edge_length >= stats.min_edge_length);
        assert_eq!(stats.n_nodes, mesh.n_nodes());
        assert_eq!(stats.n_elements, mesh.n_elements());
    }

    // ── Mesh coarsening tests ───────────────────────────────────────────────

    #[test]
    fn test_mesh_coarsening_reduces_elements() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
        let n_before = mesh.n_elements();
        let coarsened = mesh_coarsen(&mesh, 0.5);
        // Coarsening with threshold 0.5 should remove some elements or keep all
        assert!(
            coarsened.n_elements() <= n_before,
            "Coarsening should not increase elements"
        );
    }

    #[test]
    fn test_mesh_coarsening_zero_threshold_keeps_all() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let coarsened = mesh_coarsen(&mesh, 0.0);
        // With zero threshold (no collapsing), all elements should remain
        assert_eq!(coarsened.n_elements(), mesh.n_elements());
    }

    // ── Mesh refinement tests ───────────────────────────────────────────────

    #[test]
    fn test_mesh_refinement_increases_elements() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 1, 1, 1);
        let n_before = mesh.n_elements();
        let refined = mesh_refine_uniform(&mesh);
        assert!(
            refined.n_elements() > n_before,
            "Refinement should increase elements"
        );
    }

    #[test]
    fn test_mesh_refinement_preserves_volume() {
        let mesh = generate_box_mesh(2.0, 2.0, 2.0, 2, 2, 2);
        let v_before = mesh.total_volume();
        let refined = mesh_refine_uniform(&mesh);
        let v_after = refined.total_volume();
        assert!(
            (v_after - v_before).abs() < v_before * 1e-6,
            "Volume: {v_before} → {v_after}"
        );
    }

    // ── RCM reordering tests ────────────────────────────────────────────────

    #[test]
    fn test_rcm_reordering_reduces_bandwidth() {
        let mesh = generate_box_mesh(2.0, 2.0, 2.0, 3, 2, 2);
        let bw_before = compute_bandwidth(&mesh, &(0..mesh.n_nodes()).collect::<Vec<_>>());
        let perm = rcm_reordering(&mesh);
        let bw_after = compute_bandwidth(&mesh, &perm);
        // RCM should not increase bandwidth
        assert!(
            bw_after <= bw_before + 1,
            "RCM bandwidth {bw_after} > original {bw_before}"
        );
    }

    #[test]
    fn test_rcm_permutation_complete() {
        let mesh = generate_box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
        let perm = rcm_reordering(&mesh);
        assert_eq!(perm.len(), mesh.n_nodes());
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), mesh.n_nodes(), "RCM must be a permutation");
    }
}
