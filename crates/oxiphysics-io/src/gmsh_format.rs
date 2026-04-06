// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gmsh mesh format (.msh) version 2.2 ASCII support.
//!
//! Implements reading and writing of the Gmsh 2.2 ASCII `.msh` format,
//! including node coordinates, element connectivity, and physical groups.

use std::collections::HashMap;

// ── Data types ────────────────────────────────────────────────────────────────

/// A node in a Gmsh mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshNode {
    /// 1-based node number as stored in the file.
    pub id: usize,
    /// Node coordinates `[x, y, z]`.
    pub coords: [f64; 3],
}

/// Gmsh element types (subset used in practice).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum GmshElementType {
    /// 2-node line.
    Line2 = 1,
    /// 3-node triangle.
    Triangle3 = 2,
    /// 4-node quadrilateral.
    Quad4 = 3,
    /// 4-node tetrahedron.
    Tet4 = 4,
    /// 8-node hexahedron.
    Hex8 = 5,
    /// 1-node point.
    Point1 = 15,
}

impl GmshElementType {
    /// Try to construct from the integer type code.
    pub fn from_code(code: usize) -> Option<Self> {
        match code {
            1 => Some(Self::Line2),
            2 => Some(Self::Triangle3),
            3 => Some(Self::Quad4),
            4 => Some(Self::Tet4),
            5 => Some(Self::Hex8),
            15 => Some(Self::Point1),
            _ => None,
        }
    }

    /// Number of nodes for this element type.
    pub fn node_count(self) -> usize {
        match self {
            Self::Point1 => 1,
            Self::Line2 => 2,
            Self::Triangle3 => 3,
            Self::Quad4 => 4,
            Self::Tet4 => 4,
            Self::Hex8 => 8,
        }
    }
}

/// A single Gmsh mesh element.
#[derive(Debug, Clone, PartialEq)]
pub struct GmshElement {
    /// 1-based element number.
    pub id: usize,
    /// Element type code.
    pub element_type: GmshElementType,
    /// Physical group tag (first tag in the tag list, 0 if absent).
    pub physical_tag: usize,
    /// Geometric entity tag (second tag, 0 if absent).
    pub geometric_tag: usize,
    /// 1-based node connectivity.
    pub nodes: Vec<usize>,
}

/// A physical group: a named set of elements sharing a tag.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicalGroup {
    /// Spatial dimension (0–3).
    pub dimension: usize,
    /// Integer tag matching `GmshElement::physical_tag`.
    pub tag: usize,
    /// Human-readable name.
    pub name: String,
}

/// A complete Gmsh 2.2 mesh.
#[derive(Debug, Clone, Default)]
pub struct GmshMesh {
    /// Node list (1-based IDs preserved).
    pub nodes: Vec<GmshNode>,
    /// Element list.
    pub elements: Vec<GmshElement>,
    /// Physical groups declared in `$PhysicalNames`.
    pub physical_groups: Vec<PhysicalGroup>,
}

impl GmshMesh {
    /// Create an empty Gmsh mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Return the node with the given 1-based ID, or `None`.
    pub fn node_by_id(&self, id: usize) -> Option<&GmshNode> {
        self.nodes.iter().find(|n| n.id == id)
    }

    /// Count of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Count of elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
}

/// Generic FEM mesh connectivity.
#[derive(Debug, Clone, Default)]
pub struct FemMesh {
    /// Node coordinates `[x, y, z]` in node-index order (0-based).
    pub coords: Vec<[f64; 3]>,
    /// Element connectivity as lists of 0-based node indices.
    pub connectivity: Vec<Vec<usize>>,
}

// ── Parsing ───────────────────────────────────────────────────────────────────

/// Parse a Gmsh 2.2 ASCII `.msh` file from a string.
///
/// Sections recognised: `$MeshFormat`, `$PhysicalNames`, `$Nodes`, `$Elements`.
pub fn read_gmsh_v2(src: &str) -> GmshMesh {
    let mut mesh = GmshMesh::new();
    let lines: Vec<&str> = src.lines().collect();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();
        match line {
            "$MeshFormat" => {
                // Skip until $EndMeshFormat
                while i < lines.len() && lines[i].trim() != "$EndMeshFormat" {
                    i += 1;
                }
            }
            "$PhysicalNames" => {
                i += 1;
                if i < lines.len() {
                    let _count: usize = lines[i].trim().parse().unwrap_or(0);
                    i += 1;
                }
                while i < lines.len() && lines[i].trim() != "$EndPhysicalNames" {
                    let parts: Vec<&str> = lines[i].split_whitespace().collect();
                    if parts.len() >= 3 {
                        let dim = parts[0].parse::<usize>().unwrap_or(0);
                        let tag = parts[1].parse::<usize>().unwrap_or(0);
                        // Name may be quoted
                        let name = parts[2..].join(" ").trim_matches('"').to_string();
                        mesh.physical_groups.push(PhysicalGroup {
                            dimension: dim,
                            tag,
                            name,
                        });
                    }
                    i += 1;
                }
            }
            "$Nodes" => {
                i += 1;
                if i < lines.len() {
                    let _n: usize = lines[i].trim().parse().unwrap_or(0);
                    i += 1;
                }
                while i < lines.len() && lines[i].trim() != "$EndNodes" {
                    let parts: Vec<&str> = lines[i].split_whitespace().collect();
                    if parts.len() >= 4
                        && let (Ok(id), Ok(x), Ok(y), Ok(z)) = (
                            parts[0].parse::<usize>(),
                            parts[1].parse::<f64>(),
                            parts[2].parse::<f64>(),
                            parts[3].parse::<f64>(),
                        )
                    {
                        mesh.nodes.push(GmshNode {
                            id,
                            coords: [x, y, z],
                        });
                    }
                    i += 1;
                }
            }
            "$Elements" => {
                i += 1;
                if i < lines.len() {
                    let _n: usize = lines[i].trim().parse().unwrap_or(0);
                    i += 1;
                }
                while i < lines.len() && lines[i].trim() != "$EndElements" {
                    let parts: Vec<&str> = lines[i].split_whitespace().collect();
                    if parts.len() >= 4
                        && let (Ok(id), Ok(type_code), Ok(n_tags)) = (
                            parts[0].parse::<usize>(),
                            parts[1].parse::<usize>(),
                            parts[2].parse::<usize>(),
                        )
                        && let Some(etype) = GmshElementType::from_code(type_code)
                    {
                        let tag_start = 3;
                        let node_start = tag_start + n_tags;
                        let physical_tag = if n_tags >= 1 {
                            parts
                                .get(tag_start)
                                .and_then(|s| s.parse().ok())
                                .unwrap_or(0)
                        } else {
                            0
                        };
                        let geometric_tag = if n_tags >= 2 {
                            parts
                                .get(tag_start + 1)
                                .and_then(|s| s.parse().ok())
                                .unwrap_or(0)
                        } else {
                            0
                        };
                        let nodes: Vec<usize> = parts[node_start..]
                            .iter()
                            .filter_map(|s| s.parse::<usize>().ok())
                            .collect();
                        mesh.elements.push(GmshElement {
                            id,
                            element_type: etype,
                            physical_tag,
                            geometric_tag,
                            nodes,
                        });
                    }
                    i += 1;
                }
            }
            _ => {}
        }
        i += 1;
    }

    mesh
}

/// Serialise a `GmshMesh` to Gmsh 2.2 ASCII format.
pub fn write_gmsh_v2(mesh: &GmshMesh) -> String {
    let mut out = String::new();

    // Header
    out.push_str("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    // Physical names
    if !mesh.physical_groups.is_empty() {
        out.push_str("$PhysicalNames\n");
        out.push_str(&format!("{}\n", mesh.physical_groups.len()));
        for pg in &mesh.physical_groups {
            out.push_str(&format!("{} {} \"{}\"\n", pg.dimension, pg.tag, pg.name));
        }
        out.push_str("$EndPhysicalNames\n");
    }

    // Nodes
    out.push_str("$Nodes\n");
    out.push_str(&format!("{}\n", mesh.nodes.len()));
    for node in &mesh.nodes {
        out.push_str(&format!(
            "{} {} {} {}\n",
            node.id, node.coords[0], node.coords[1], node.coords[2]
        ));
    }
    out.push_str("$EndNodes\n");

    // Elements
    out.push_str("$Elements\n");
    out.push_str(&format!("{}\n", mesh.elements.len()));
    for el in &mesh.elements {
        let type_code = el.element_type as usize;
        // Always write 2 tags: physical + geometric
        let node_str: Vec<String> = el.nodes.iter().map(|n| n.to_string()).collect();
        out.push_str(&format!(
            "{} {} 2 {} {} {}\n",
            el.id,
            type_code,
            el.physical_tag,
            el.geometric_tag,
            node_str.join(" ")
        ));
    }
    out.push_str("$EndElements\n");

    out
}

/// Collect the 1-based node IDs belonging to a specific physical group tag.
///
/// An element belongs to the group if `element.physical_tag == tag`.
/// All node IDs referenced by matching elements are returned (deduplicated,
/// sorted).
pub fn physical_group_nodes(mesh: &GmshMesh, tag: usize) -> Vec<usize> {
    let mut node_set: std::collections::BTreeSet<usize> = Default::default();
    for el in &mesh.elements {
        if el.physical_tag == tag {
            for &nid in &el.nodes {
                node_set.insert(nid);
            }
        }
    }
    node_set.into_iter().collect()
}

/// Convert a `GmshMesh` to a generic `FemMesh`.
///
/// Nodes are re-indexed to 0-based.  Only elements with a known element type
/// are included.
pub fn convert_to_fem_mesh(mesh: &GmshMesh) -> FemMesh {
    // Build a map from 1-based Gmsh ID → 0-based FEM index
    let mut id_to_idx: HashMap<usize, usize> = HashMap::new();
    let mut coords: Vec<[f64; 3]> = Vec::with_capacity(mesh.nodes.len());

    for (idx, node) in mesh.nodes.iter().enumerate() {
        id_to_idx.insert(node.id, idx);
        coords.push(node.coords);
    }

    let mut connectivity: Vec<Vec<usize>> = Vec::with_capacity(mesh.elements.len());
    for el in &mesh.elements {
        let zero_based: Vec<usize> = el
            .nodes
            .iter()
            .filter_map(|nid| id_to_idx.get(nid).copied())
            .collect();
        connectivity.push(zero_based);
    }

    FemMesh {
        coords,
        connectivity,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle_msh_src() -> &'static str {
        "$MeshFormat\n\
         2.2 0 8\n\
         $EndMeshFormat\n\
         $PhysicalNames\n\
         1\n\
         2 1 \"Surface\"\n\
         $EndPhysicalNames\n\
         $Nodes\n\
         3\n\
         1 0.0 0.0 0.0\n\
         2 1.0 0.0 0.0\n\
         3 0.0 1.0 0.0\n\
         $EndNodes\n\
         $Elements\n\
         1\n\
         1 2 2 1 0 1 2 3\n\
         $EndElements\n"
    }

    // --- read_gmsh_v2 ---

    #[test]
    fn test_read_node_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.node_count(), 3);
    }

    #[test]
    fn test_read_node_coords() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let n1 = mesh.node_by_id(1).unwrap();
        assert_eq!(n1.coords, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_read_node_id_preserved() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert!(mesh.node_by_id(2).is_some());
        assert!(mesh.node_by_id(99).is_none());
    }

    #[test]
    fn test_read_element_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.element_count(), 1);
    }

    #[test]
    fn test_read_element_type() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.elements[0].element_type, GmshElementType::Triangle3);
    }

    #[test]
    fn test_read_element_physical_tag() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.elements[0].physical_tag, 1);
    }

    #[test]
    fn test_read_element_nodes() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.elements[0].nodes, vec![1, 2, 3]);
    }

    #[test]
    fn test_read_physical_group_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.physical_groups.len(), 1);
    }

    #[test]
    fn test_read_physical_group_name() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.physical_groups[0].name, "Surface");
    }

    #[test]
    fn test_read_physical_group_dimension() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        assert_eq!(mesh.physical_groups[0].dimension, 2);
    }

    #[test]
    fn test_read_empty_string() {
        let mesh = read_gmsh_v2("");
        assert_eq!(mesh.node_count(), 0);
        assert_eq!(mesh.element_count(), 0);
    }

    // --- GmshElementType ---

    #[test]
    fn test_element_type_from_code_triangle() {
        assert_eq!(
            GmshElementType::from_code(2),
            Some(GmshElementType::Triangle3)
        );
    }

    #[test]
    fn test_element_type_from_code_unknown() {
        assert!(GmshElementType::from_code(999).is_none());
    }

    #[test]
    fn test_element_type_node_count_tet() {
        assert_eq!(GmshElementType::Tet4.node_count(), 4);
    }

    #[test]
    fn test_element_type_node_count_triangle() {
        assert_eq!(GmshElementType::Triangle3.node_count(), 3);
    }

    // --- write_gmsh_v2 ---

    #[test]
    fn test_write_contains_mesh_format() {
        let mesh = GmshMesh::new();
        let out = write_gmsh_v2(&mesh);
        assert!(out.contains("$MeshFormat"));
        assert!(out.contains("2.2 0 8"));
    }

    #[test]
    fn test_write_contains_nodes_section() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        assert!(out.contains("$Nodes"));
        assert!(out.contains("$EndNodes"));
    }

    #[test]
    fn test_write_contains_elements_section() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        assert!(out.contains("$Elements"));
        assert!(out.contains("$EndElements"));
    }

    // --- Roundtrip ---

    #[test]
    fn test_roundtrip_node_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        let mesh2 = read_gmsh_v2(&out);
        assert_eq!(mesh2.node_count(), mesh.node_count());
    }

    #[test]
    fn test_roundtrip_element_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        let mesh2 = read_gmsh_v2(&out);
        assert_eq!(mesh2.element_count(), mesh.element_count());
    }

    #[test]
    fn test_roundtrip_node_coords() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        let mesh2 = read_gmsh_v2(&out);
        let n2 = mesh2.node_by_id(2).unwrap();
        assert!((n2.coords[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_roundtrip_physical_group_name() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        let mesh2 = read_gmsh_v2(&out);
        assert_eq!(mesh2.physical_groups[0].name, "Surface");
    }

    #[test]
    fn test_roundtrip_element_type() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let out = write_gmsh_v2(&mesh);
        let mesh2 = read_gmsh_v2(&out);
        assert_eq!(mesh2.elements[0].element_type, GmshElementType::Triangle3);
    }

    // --- physical_group_nodes ---

    #[test]
    fn test_physical_group_nodes_correct() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let nodes = physical_group_nodes(&mesh, 1);
        assert_eq!(nodes, vec![1, 2, 3]);
    }

    #[test]
    fn test_physical_group_nodes_missing_tag() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let nodes = physical_group_nodes(&mesh, 999);
        assert!(nodes.is_empty());
    }

    #[test]
    fn test_physical_group_nodes_deduplicated() {
        // Two elements sharing node 1
        let src = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n\
                   $Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 1 1 0\n$EndNodes\n\
                   $Elements\n2\n\
                   1 2 2 1 0 1 2 3\n\
                   2 2 2 1 0 2 3 4\n\
                   $EndElements\n";
        let mesh = read_gmsh_v2(src);
        let nodes = physical_group_nodes(&mesh, 1);
        // Nodes 1,2,3,4 — no duplicates
        assert_eq!(nodes.len(), 4);
    }

    // --- convert_to_fem_mesh ---

    #[test]
    fn test_fem_mesh_node_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let fem = convert_to_fem_mesh(&mesh);
        assert_eq!(fem.coords.len(), 3);
    }

    #[test]
    fn test_fem_mesh_connectivity_count() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let fem = convert_to_fem_mesh(&mesh);
        assert_eq!(fem.connectivity.len(), 1);
    }

    #[test]
    fn test_fem_mesh_connectivity_zero_based() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let fem = convert_to_fem_mesh(&mesh);
        // 1-based [1,2,3] → 0-based [0,1,2]
        assert_eq!(fem.connectivity[0], vec![0, 1, 2]);
    }

    #[test]
    fn test_fem_mesh_coords_correct() {
        let mesh = read_gmsh_v2(triangle_msh_src());
        let fem = convert_to_fem_mesh(&mesh);
        // Node 2 has coords [1,0,0]
        assert!((fem.coords[1][0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_fem_mesh_empty() {
        let mesh = GmshMesh::new();
        let fem = convert_to_fem_mesh(&mesh);
        assert!(fem.coords.is_empty());
        assert!(fem.connectivity.is_empty());
    }
}
