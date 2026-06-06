//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{DirichletBc, ElementTypeMapping, FeElementType, MultiPointConstraint};
use std::fmt::Write as _;

/// Degree of freedom index (1=Tx, 2=Ty, 3=Tz, 4=Rx, 5=Ry, 6=Rz).
pub type DofIndex = u8;
/// Full element type mapping table.
pub fn element_type_table() -> Vec<ElementTypeMapping> {
    vec![
        ElementTypeMapping {
            fe_type: FeElementType::Line2,
            abaqus: "T3D2",
            nastran: "CBAR",
            gmsh_tag: 1,
            vtk_code: 3,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Tri3,
            abaqus: "S3",
            nastran: "CTRIA3",
            gmsh_tag: 2,
            vtk_code: 5,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Quad4,
            abaqus: "S4",
            nastran: "CQUAD4",
            gmsh_tag: 3,
            vtk_code: 9,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Tet4,
            abaqus: "C3D4",
            nastran: "CTETRA",
            gmsh_tag: 4,
            vtk_code: 10,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Hex8,
            abaqus: "C3D8",
            nastran: "CHEXA",
            gmsh_tag: 5,
            vtk_code: 12,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Tri6,
            abaqus: "S6",
            nastran: "CTRIA6",
            gmsh_tag: 9,
            vtk_code: 22,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Tet10,
            abaqus: "C3D10",
            nastran: "CTETRA10",
            gmsh_tag: 11,
            vtk_code: 24,
        },
        ElementTypeMapping {
            fe_type: FeElementType::Hex20,
            abaqus: "C3D20",
            nastran: "CHEXA20",
            gmsh_tag: 17,
            vtk_code: 25,
        },
    ]
}
/// Export Dirichlet BCs to a simple text format.
pub fn export_dirichlet_bcs(bcs: &[DirichletBc]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# Dirichlet boundary conditions");
    let _ = writeln!(out, "# node_id dof value");
    for bc in bcs {
        let _ = writeln!(out, "{} {} {:.10e}", bc.node_id, bc.dof, bc.value);
    }
    out
}
/// Export MPCs to a simple text format.
pub fn export_mpcs(mpcs: &[MultiPointConstraint]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# Multi-point constraints");
    for mpc in mpcs {
        let _ = write!(
            out,
            "DEPENDENT {} {} =",
            mpc.dependent_node, mpc.dependent_dof
        );
        for (n, d, c) in &mpc.terms {
            let _ = write!(out, " {:.6e}*N{}DOF{}", c, n, d);
        }
        let _ = writeln!(out, " + {:.6e}", mpc.rhs);
    }
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_element_io::AbaqusInpIo;
    use crate::finite_element_io::ExodusLikeMesh;
    use crate::finite_element_io::FeElement;
    use crate::finite_element_io::FeMesh;
    use crate::finite_element_io::FeNode;
    use crate::finite_element_io::GmshMesh;
    use crate::finite_element_io::GmshVersion;
    use crate::finite_element_io::LinearElasticMaterial;
    use crate::finite_element_io::MapdlDeck;
    use crate::finite_element_io::MeshPartition;
    use crate::finite_element_io::NastranBulkData;
    use crate::finite_element_io::NastranGrid;
    use crate::finite_element_io::NodalDisplacements;
    use crate::finite_element_io::NodalStresses;
    use crate::finite_element_io::RestartCheckpoint;
    use crate::finite_element_io::RoundRobinPartitioner;
    use crate::finite_element_io::VtkResultExporter;
    use crate::finite_element_io::types::*;
    #[test]
    fn test_fe_node_new() {
        let n = FeNode::new(1, [1.0, 2.0, 3.0]);
        assert_eq!(n.id, 1);
        assert_eq!(n.coords, [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_fe_node_distance() {
        let a = FeNode::new(1, [0.0, 0.0, 0.0]);
        let b = FeNode::new(2, [3.0, 4.0, 0.0]);
        assert!((a.distance(&b) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_element_type_node_count() {
        assert_eq!(FeElementType::Tet4.node_count(), Some(4));
        assert_eq!(FeElementType::Hex8.node_count(), Some(8));
        assert_eq!(FeElementType::Tri3.node_count(), Some(3));
        assert_eq!(FeElementType::Line2.node_count(), Some(2));
        assert_eq!(FeElementType::Unknown("X".into()).node_count(), None);
    }
    #[test]
    fn test_element_type_volumetric() {
        assert!(FeElementType::Tet4.is_volumetric());
        assert!(FeElementType::Hex8.is_volumetric());
        assert!(!FeElementType::Tri3.is_volumetric());
        assert!(!FeElementType::Line2.is_volumetric());
    }
    #[test]
    fn test_material_shear_modulus() {
        let mat = LinearElasticMaterial::new(1, "Steel", 200e9, 0.3, 7800.0);
        let g = mat.shear_modulus();
        let expected = 200e9 / (2.0 * 1.3);
        assert!((g - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_material_bulk_modulus() {
        let mat = LinearElasticMaterial::new(1, "Steel", 200e9, 0.3, 7800.0);
        let k = mat.bulk_modulus();
        let expected = 200e9 / (3.0 * 0.4);
        assert!((k - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_dirichlet_bc_fixed() {
        let bc = DirichletBc::fixed(10, 3);
        assert_eq!(bc.node_id, 10);
        assert_eq!(bc.dof, 3);
        assert_eq!(bc.value, 0.0);
    }
    #[test]
    fn test_dirichlet_bc_prescribed() {
        let bc = DirichletBc::prescribed(5, 1, 0.001);
        assert!((bc.value - 0.001).abs() < 1e-15);
    }
    #[test]
    fn test_fe_mesh_bounding_box_empty() {
        let mesh = FeMesh::new();
        assert!(mesh.bounding_box().is_none());
    }
    #[test]
    fn test_fe_mesh_bounding_box() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        mesh.nodes.push(FeNode::new(2, [1.0, 2.0, 3.0]));
        let (lo, hi) = mesh.bounding_box().unwrap();
        assert_eq!(lo, [0.0, 0.0, 0.0]);
        assert_eq!(hi, [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_fe_mesh_find_node() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(42, [1.0, 2.0, 3.0]));
        assert!(mesh.find_node(42).is_some());
        assert!(mesh.find_node(99).is_none());
    }
    #[test]
    fn test_abaqus_parse_nodes() {
        let inp = "*NODE\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n3, 0.0, 1.0, 0.0\n";
        let mesh = AbaqusInpIo::parse(inp);
        assert_eq!(mesh.nodes.len(), 3);
        assert_eq!(mesh.nodes[0].id, 1);
        assert!((mesh.nodes[1].coords[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_abaqus_parse_elements() {
        let inp = "*ELEMENT, TYPE=C3D4\n1, 1, 2, 3, 4\n2, 5, 6, 7, 8\n";
        let mesh = AbaqusInpIo::parse(inp);
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Tet4);
        assert_eq!(mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
    }
    #[test]
    fn test_abaqus_parse_material() {
        let inp = "*MATERIAL, NAME=Steel\n*ELASTIC\n200e9, 0.3\n";
        let mesh = AbaqusInpIo::parse(inp);
        assert_eq!(mesh.materials.len(), 1);
        assert!((mesh.materials[0].young_modulus - 200e9).abs() < 1.0);
    }
    #[test]
    fn test_abaqus_write_contains_node_keyword() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        let out = AbaqusInpIo::write(&mesh);
        assert!(out.contains("*NODE"));
        assert!(out.contains("1,"));
    }
    #[test]
    fn test_abaqus_parse_nset() {
        let inp = "*NSET, NSET=Fixed\n1, 2, 3\n";
        let mesh = AbaqusInpIo::parse(inp);
        let ids = mesh.node_sets.get("Fixed").cloned().unwrap_or_default();
        assert!(ids.contains(&1));
        assert!(ids.contains(&3));
    }
    #[test]
    fn test_nastran_parse_grid() {
        let bdf = "BEGIN BULK\nGRID,1,0,1.0,2.0,3.0,0\nENDDATA\n";
        let data = NastranBulkData::parse(bdf);
        assert_eq!(data.grids.len(), 1);
        assert_eq!(data.grids[0].id, 1);
        assert!((data.grids[0].coords[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_nastran_parse_ctetra() {
        let bdf = "BEGIN BULK\nCTETRA,1,1,10,20,30,40\nENDDATA\n";
        let data = NastranBulkData::parse(bdf);
        assert_eq!(data.elements.len(), 1);
        assert_eq!(data.elements[0].element_type, FeElementType::Tet4);
    }
    #[test]
    fn test_nastran_parse_mat1() {
        let bdf = "BEGIN BULK\nMAT1,1,200E9,77E9,0.3,7800\nENDDATA\n";
        let data = NastranBulkData::parse(bdf);
        assert_eq!(data.materials.len(), 1);
        assert_eq!(data.materials[0].mid, 1);
        assert!((data.materials[0].e - 200e9).abs() < 1.0);
    }
    #[test]
    fn test_nastran_parse_spc() {
        let bdf = "BEGIN BULK\nSPC,1,100,123456,0.0\nENDDATA\n";
        let data = NastranBulkData::parse(bdf);
        assert_eq!(data.spcs.len(), 1);
        assert_eq!(data.spcs[0].components, "123456");
    }
    #[test]
    fn test_nastran_write_roundtrip() {
        let mut data = NastranBulkData::new();
        data.grids.push(NastranGrid::new(1, [1.0, 2.0, 3.0]));
        let text = data.write();
        assert!(text.contains("BEGIN BULK"));
        assert!(text.contains("GRID"));
        assert!(text.contains("ENDDATA"));
    }
    #[test]
    fn test_gmsh_parse_v2_nodes() {
        let msh = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n\
                   $Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n\
                   $Elements\n0\n$EndElements\n";
        let mesh = GmshMesh::parse(msh);
        assert_eq!(mesh.version, Some(GmshVersion::V2));
        assert_eq!(mesh.nodes.len(), 3);
    }
    #[test]
    fn test_gmsh_parse_elements_tet() {
        let msh = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n\
                   $Nodes\n0\n$EndNodes\n\
                   $Elements\n1\n1 4 0 1 2 3 4\n$EndElements\n";
        let mesh = GmshMesh::parse(msh);
        assert_eq!(mesh.elements.len(), 1);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Tet4);
    }
    #[test]
    fn test_gmsh_write_v2_roundtrip() {
        let mut mesh = GmshMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        mesh.elements
            .push(FeElement::new(1, FeElementType::Tet4, vec![1, 2, 3, 4]));
        let text = mesh.write_v2();
        assert!(text.contains("$MeshFormat"));
        assert!(text.contains("$Nodes"));
        assert!(text.contains("$Elements"));
    }
    #[test]
    fn test_exodus_coords() {
        let mut ex = ExodusLikeMesh::new("Test");
        ex.num_nodes = 2;
        ex.coordinates = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert!((ex.x(0) - 1.0).abs() < 1e-10);
        assert!((ex.y(0) - 2.0).abs() < 1e-10);
        assert!((ex.z(0) - 3.0).abs() < 1e-10);
        assert!((ex.x(1) - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_exodus_total_elements() {
        let mut ex = ExodusLikeMesh::new("Test");
        ex.element_blocks
            .insert(1, (FeElementType::Tet4, vec![1, 2, 3, 4, 5, 6, 7, 8]));
        assert_eq!(ex.total_elements(), 2);
    }
    #[test]
    fn test_exodus_serialize_text() {
        let mut ex = ExodusLikeMesh::new("MyMesh");
        ex.num_nodes = 1;
        ex.coordinates = vec![0.0, 0.0, 0.0];
        let text = ex.serialize_text();
        assert!(text.contains("MyMesh"));
        assert!(text.contains("NUM_NODES: 1"));
    }
    #[test]
    fn test_nodal_displacements_magnitude() {
        let mut nd = NodalDisplacements::zeros(vec![1, 2]);
        nd.displacements[0] = [3.0, 4.0, 0.0];
        assert!((nd.magnitude(0) - 5.0).abs() < 1e-10);
        assert_eq!(nd.magnitude(1), 0.0);
    }
    #[test]
    fn test_nodal_displacements_max_magnitude() {
        let mut nd = NodalDisplacements::zeros(vec![1, 2, 3]);
        nd.displacements[0] = [1.0, 0.0, 0.0];
        nd.displacements[1] = [0.0, 2.0, 0.0];
        nd.displacements[2] = [0.0, 0.0, 0.5];
        assert!((nd.max_magnitude() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_nodal_stresses_von_mises_uniaxial() {
        let mut ns = NodalStresses::zeros(vec![1]);
        ns.stresses[0] = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = ns.von_mises(0);
        assert!((vm - 100.0).abs() < 1e-8);
    }
    #[test]
    fn test_nodal_stresses_von_mises_hydrostatic() {
        let mut ns = NodalStresses::zeros(vec![1]);
        ns.stresses[0] = [50.0, 50.0, 50.0, 0.0, 0.0, 0.0];
        let vm = ns.von_mises(0);
        assert!(
            vm < 1e-8,
            "hydrostatic stress should have von Mises ≈ 0, got {}",
            vm
        );
    }
    #[test]
    fn test_vtk_export_contains_header() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        mesh.elements
            .push(FeElement::new(1, FeElementType::Tet4, vec![1, 2, 3, 4]));
        let out = VtkResultExporter::export_legacy_ascii(&mesh, None, None);
        assert!(out.contains("# vtk DataFile Version 3.0"));
        assert!(out.contains("DATASET UNSTRUCTURED_GRID"));
    }
    #[test]
    fn test_vtk_export_with_displacements() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        mesh.elements
            .push(FeElement::new(1, FeElementType::Tet4, vec![1]));
        let disp = NodalDisplacements {
            node_ids: vec![1],
            displacements: vec![[0.001, 0.002, 0.003]],
        };
        let out = VtkResultExporter::export_legacy_ascii(&mesh, Some(&disp), None);
        assert!(out.contains("VECTORS displacement double"));
    }
    #[test]
    fn test_round_robin_partition_node_count() {
        let mut mesh = FeMesh::new();
        for i in 1..=6 {
            mesh.nodes.push(FeNode::new(i, [0.0, 0.0, 0.0]));
        }
        let parts = RoundRobinPartitioner::partition(&mesh, 3);
        assert_eq!(parts.len(), 3);
        let total_owned: usize = parts.iter().map(|p| p.owned_nodes.len()).sum();
        assert_eq!(total_owned, 6);
    }
    #[test]
    fn test_round_robin_partition_single() {
        let mut mesh = FeMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        let parts = RoundRobinPartitioner::partition(&mesh, 1);
        assert_eq!(parts[0].owned_nodes.len(), 1);
    }
    #[test]
    fn test_element_type_table_count() {
        let table = element_type_table();
        assert!(table.len() >= 8);
    }
    #[test]
    fn test_element_type_table_tet4_vtk() {
        let table = element_type_table();
        let entry = table
            .iter()
            .find(|e| e.fe_type == FeElementType::Tet4)
            .expect("Tet4 must be in the table");
        assert_eq!(entry.vtk_code, 10);
        assert_eq!(entry.abaqus, "C3D4");
        assert_eq!(entry.gmsh_tag, 4);
    }
    #[test]
    fn test_export_dirichlet_bcs() {
        let bcs = vec![
            DirichletBc::fixed(1, 1),
            DirichletBc::prescribed(2, 2, 0.005),
        ];
        let text = export_dirichlet_bcs(&bcs);
        assert!(text.contains("1 1 0"));
        assert!(text.contains("2 2"));
    }
    #[test]
    fn test_export_mpcs() {
        let mut mpc = MultiPointConstraint::new(10, 3);
        mpc.add_term(20, 3, -1.0);
        let text = export_mpcs(&[mpc]);
        assert!(text.contains("DEPENDENT 10 3"));
        assert!(text.contains("N20DOF3"));
    }
    #[test]
    fn test_restart_serialize_deserialize() {
        let mut ckpt = RestartCheckpoint::new("Beam simulation", 5, 0.05);
        ckpt.displacements = vec![0.001, 0.002, 0.003];
        ckpt.velocities = vec![0.1, 0.2, 0.3];
        ckpt.residual_norm = 1e-8;
        ckpt.metadata.insert("solver".to_string(), "CG".to_string());
        let text = ckpt.serialize();
        let loaded = RestartCheckpoint::deserialize(&text).unwrap();
        assert_eq!(loaded.title, "Beam simulation");
        assert_eq!(loaded.step, 5);
        assert!((loaded.time - 0.05).abs() < 1e-12);
        assert_eq!(loaded.displacements.len(), 3);
        assert!((loaded.displacements[0] - 0.001).abs() < 1e-12);
        assert_eq!(
            loaded.metadata.get("solver").map(|s| s.as_str()),
            Some("CG")
        );
    }
    #[test]
    fn test_restart_empty_displacements() {
        let ckpt = RestartCheckpoint::new("Empty", 0, 0.0);
        let text = ckpt.serialize();
        let loaded = RestartCheckpoint::deserialize(&text).unwrap();
        assert!(loaded.displacements.is_empty());
    }
    #[test]
    fn test_mapdl_parse_material_ex() {
        let deck_str = "MP,EX,1,200E9\nMP,NUXY,1,0.3\n";
        let deck = MapdlDeck::parse(deck_str);
        let mat = deck.materials.get(&1).copied().unwrap_or_default();
        assert!((mat.0 - 200e9).abs() < 1.0);
        assert!((mat.1 - 0.3).abs() < 1e-10);
    }
    #[test]
    fn test_mapdl_parse_displacement_bc() {
        let deck_str = "D,100,UX,0.0\nD,100,UY,0.0\n";
        let deck = MapdlDeck::parse(deck_str);
        let constraints = deck.disp_constraints.get(&100).cloned().unwrap_or_default();
        assert_eq!(constraints.len(), 2);
    }
    #[test]
    fn test_mapdl_parse_force() {
        let deck_str = "F,50,FZ,1000.0\n";
        let deck = MapdlDeck::parse(deck_str);
        let forces = deck.nodal_forces.get(&50).cloned().unwrap_or_default();
        assert_eq!(forces.len(), 1);
        assert!((forces[0].1 - 1000.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_partition_total_nodes() {
        let mut p = MeshPartition::new(0, 4);
        p.owned_nodes = vec![1, 2, 3];
        p.ghost_nodes = vec![4, 5];
        assert_eq!(p.total_nodes(), 5);
    }
    #[test]
    fn test_mesh_partition_neighbour_fraction() {
        let mut p = MeshPartition::new(0, 4);
        p.neighbours = vec![1, 2];
        let frac = p.neighbour_fraction();
        assert!((frac - 2.0 / 3.0).abs() < 1e-10);
    }
}
