//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BrepSolid, TriangleMesh};

/// Cross product of two 3-vectors.
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Length of a 3-vector.
pub(super) fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Normalize a 3-vector.
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / l, v[1] / l, v[2] / l]
    }
}
/// Subtract two 3-vectors.
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3-vectors.
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector.
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Linearly interpolate between two points.
pub(super) fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}
/// Tessellate a BREP solid into a triangle mesh.
///
/// Uses simple fan triangulation for planar faces with straight edges.
pub fn tessellate_brep(solid: &BrepSolid, _segments_per_edge: usize) -> TriangleMesh {
    let mut mesh = TriangleMesh::new();
    let mut vert_map: Vec<usize> = Vec::with_capacity(solid.vertices.len());
    for v in &solid.vertices {
        vert_map.push(mesh.add_vertex(v.position));
    }
    for face in &solid.faces {
        if face.edge_loops.is_empty() {
            continue;
        }
        let outer = &face.edge_loops[0];
        let mut face_verts: Vec<usize> = Vec::new();
        for &eidx in outer {
            let edge = &solid.edges[eidx];
            let sv = vert_map[edge.start_vertex];
            if face_verts.is_empty()
                || *face_verts.last().expect("collection should not be empty") != sv
            {
                face_verts.push(sv);
            }
        }
        if face_verts.len() >= 3 {
            for i in 1..face_verts.len() - 1 {
                mesh.add_triangle(face_verts[0], face_verts[i], face_verts[i + 1]);
            }
        }
    }
    mesh.compute_normals();
    mesh
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cad_io::Assembly;
    use crate::cad_io::AssemblyComponent;
    use crate::cad_io::AssemblyTransform;
    use crate::cad_io::BoundingBox;
    use crate::cad_io::IgesCircularArc;
    use crate::cad_io::IgesEntityType;
    use crate::cad_io::LengthUnit;
    use crate::cad_io::StepParser;
    use crate::cad_io::StlExporter;
    use crate::cad_io::UnitConverter;
    #[test]
    fn test_unit_mm_to_m() {
        let val = LengthUnit::Millimeter.convert(1000.0, LengthUnit::Meter);
        assert!((val - 1.0).abs() < 1e-10, "val={val}");
    }
    #[test]
    fn test_unit_inch_to_mm() {
        let val = LengthUnit::Inch.convert(1.0, LengthUnit::Millimeter);
        assert!((val - 25.4).abs() < 1e-10, "val={val}");
    }
    #[test]
    fn test_unit_converter_point() {
        let conv = UnitConverter::new(LengthUnit::Meter, LengthUnit::Millimeter);
        let p = conv.convert_point([1.0, 2.0, 3.0]);
        assert!((p[0] - 1000.0).abs() < 1e-6);
        assert!((p[1] - 2000.0).abs() < 1e-6);
    }
    #[test]
    fn test_unit_converter_points() {
        let conv = UnitConverter::new(LengthUnit::Centimeter, LengthUnit::Meter);
        let pts = conv.convert_points(&[[100.0, 200.0, 300.0]]);
        assert!((pts[0][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_from_points() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]];
        let bb = BoundingBox::from_points(&points);
        assert!((bb.min[0] - (-1.0)).abs() < 1e-10);
        assert!((bb.max[0] - 1.0).abs() < 1e-10);
        assert!((bb.max[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_center() {
        let bb = BoundingBox::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = bb.center();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_volume() {
        let bb = BoundingBox::new([0.0; 3], [2.0, 3.0, 4.0]);
        assert!((bb.volume() - 24.0).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_surface_area() {
        let bb = BoundingBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        assert!((bb.surface_area() - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_bbox_contains_point() {
        let bb = BoundingBox::new([0.0; 3], [1.0, 1.0, 1.0]);
        assert!(bb.contains_point([0.5, 0.5, 0.5]));
        assert!(!bb.contains_point([2.0, 0.5, 0.5]));
    }
    #[test]
    fn test_bbox_intersects() {
        let a = BoundingBox::new([0.0; 3], [2.0, 2.0, 2.0]);
        let b = BoundingBox::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        assert!(a.intersects(&b));
        let c = BoundingBox::new([5.0; 3], [6.0; 3]);
        assert!(!a.intersects(&c));
    }
    #[test]
    fn test_bbox_diagonal() {
        let bb = BoundingBox::new([0.0; 3], [3.0, 4.0, 0.0]);
        assert!((bb.diagonal() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_brep_box_topology() {
        let solid = BrepSolid::create_box("test_box", 1.0, 1.0, 1.0);
        assert_eq!(solid.vertices.len(), 8);
        assert_eq!(solid.edges.len(), 12);
        assert_eq!(solid.faces.len(), 6);
        assert_eq!(solid.euler_characteristic(), 2);
    }
    #[test]
    fn test_brep_validate() {
        let solid = BrepSolid::create_box("test", 1.0, 1.0, 1.0);
        assert!(solid.validate().is_ok());
    }
    #[test]
    fn test_brep_bounding_box() {
        let solid = BrepSolid::create_box("test", 2.0, 3.0, 4.0);
        let bb = solid.bounding_box();
        assert!((bb.min[0]).abs() < 1e-10);
        assert!((bb.max[0] - 2.0).abs() < 1e-10);
        assert!((bb.max[1] - 3.0).abs() < 1e-10);
        assert!((bb.max[2] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_brep_edge_length() {
        let mut solid = BrepSolid::new("test");
        let v0 = solid.add_vertex([0.0, 0.0, 0.0]);
        let v1 = solid.add_vertex([3.0, 4.0, 0.0]);
        solid.add_line_edge(v0, v1);
        let length = solid.edges[0].approximate_length(&solid.vertices);
        assert!((length - 5.0).abs() < 1e-6, "length={length}");
    }
    #[test]
    fn test_tessellate_box() {
        let solid = BrepSolid::create_box("test", 1.0, 1.0, 1.0);
        let mesh = tessellate_brep(&solid, 10);
        assert!(mesh.triangle_count() > 0);
        assert!(mesh.vertex_count() > 0);
    }
    #[test]
    fn test_mesh_surface_area() {
        let mut mesh = TriangleMesh::new();
        mesh.add_vertex([0.0, 0.0, 0.0]);
        mesh.add_vertex([1.0, 0.0, 0.0]);
        mesh.add_vertex([0.0, 1.0, 0.0]);
        mesh.add_triangle(0, 1, 2);
        assert!((mesh.surface_area() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_merge() {
        let mut m1 = TriangleMesh::new();
        m1.add_vertex([0.0; 3]);
        m1.add_vertex([1.0, 0.0, 0.0]);
        m1.add_vertex([0.0, 1.0, 0.0]);
        m1.add_triangle(0, 1, 2);
        let mut m2 = TriangleMesh::new();
        m2.add_vertex([2.0, 0.0, 0.0]);
        m2.add_vertex([3.0, 0.0, 0.0]);
        m2.add_vertex([2.0, 1.0, 0.0]);
        m2.add_triangle(0, 1, 2);
        m1.merge(&m2);
        assert_eq!(m1.vertex_count(), 6);
        assert_eq!(m1.triangle_count(), 2);
        assert_eq!(m1.triangles[1], [3, 4, 5]);
    }
    #[test]
    fn test_mesh_normals() {
        let mut mesh = TriangleMesh::new();
        mesh.add_vertex([0.0, 0.0, 0.0]);
        mesh.add_vertex([1.0, 0.0, 0.0]);
        mesh.add_vertex([0.0, 1.0, 0.0]);
        mesh.add_triangle(0, 1, 2);
        mesh.compute_normals();
        assert_eq!(mesh.normals.len(), 1);
        assert!((mesh.normals[0][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_stl_ascii_export() {
        let mut mesh = TriangleMesh::new();
        mesh.add_vertex([0.0, 0.0, 0.0]);
        mesh.add_vertex([1.0, 0.0, 0.0]);
        mesh.add_vertex([0.0, 1.0, 0.0]);
        mesh.add_triangle(0, 1, 2);
        mesh.compute_normals();
        let exporter = StlExporter::new("test");
        let stl = exporter.to_ascii_stl(&mesh);
        assert!(stl.contains("solid test"));
        assert!(stl.contains("endsolid test"));
        assert!(stl.contains("facet normal"));
        assert!(stl.contains("vertex"));
    }
    #[test]
    fn test_stl_binary_export() {
        let mut mesh = TriangleMesh::new();
        mesh.add_vertex([0.0, 0.0, 0.0]);
        mesh.add_vertex([1.0, 0.0, 0.0]);
        mesh.add_vertex([0.0, 1.0, 0.0]);
        mesh.add_triangle(0, 1, 2);
        mesh.compute_normals();
        let exporter = StlExporter::new("test");
        let data = exporter.to_binary_stl(&mesh);
        assert_eq!(data.len(), 80 + 4 + 50);
    }
    #[test]
    fn test_step_parse_entities() {
        let content = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('test'),'2;1');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('origin',(0.0,0.0,0.0));
#2=CARTESIAN_POINT('p1',(1.0,2.0,3.0));
#3=DIRECTION('dir',(0.0,0.0,1.0));
ENDSEC;
END-ISO-10303-21;
"#;
        let mut parser = StepParser::new();
        assert!(parser.parse(content).is_ok());
        assert_eq!(parser.entity_count(), 3);
    }
    #[test]
    fn test_step_find_entities() {
        let content = r#"DATA;
#1=CARTESIAN_POINT('a',(1.0,2.0,3.0));
#2=DIRECTION('d',(0.0,0.0,1.0));
#3=CARTESIAN_POINT('b',(4.0,5.0,6.0));
ENDSEC;
"#;
        let mut parser = StepParser::new();
        parser.parse(content).unwrap();
        let points = parser.find_entities("CARTESIAN_POINT");
        assert_eq!(points.len(), 2);
    }
    #[test]
    fn test_step_find_by_id() {
        let content = "DATA;\n#42=CARTESIAN_POINT('p',(1.0,2.0,3.0));\nENDSEC;\n";
        let mut parser = StepParser::new();
        parser.parse(content).unwrap();
        let e = parser.find_by_id(42);
        assert!(e.is_some());
        assert_eq!(e.unwrap().entity_type, "CARTESIAN_POINT");
    }
    #[test]
    fn test_step_extract_points() {
        let content = "DATA;\n#1=CARTESIAN_POINT('p',(10.0,20.0,30.0));\nENDSEC;\n";
        let mut parser = StepParser::new();
        parser.parse(content).unwrap();
        let points = parser.extract_cartesian_points();
        assert_eq!(points.len(), 1);
        assert!((points[0].1[0] - 10.0).abs() < 1e-10);
        assert!((points[0].1[1] - 20.0).abs() < 1e-10);
        assert!((points[0].1[2] - 30.0).abs() < 1e-10);
    }
    #[test]
    fn test_iges_entity_type_roundtrip() {
        let etype = IgesEntityType::from_type_number(110);
        assert_eq!(etype, IgesEntityType::Line);
        assert_eq!(etype.to_type_number(), 110);
    }
    #[test]
    fn test_iges_entity_type_curve() {
        assert!(IgesEntityType::Line.is_curve());
        assert!(IgesEntityType::CircularArc.is_curve());
        assert!(!IgesEntityType::Plane.is_curve());
    }
    #[test]
    fn test_iges_entity_type_surface() {
        assert!(IgesEntityType::Plane.is_surface());
        assert!(IgesEntityType::RuledSurface.is_surface());
        assert!(!IgesEntityType::Line.is_surface());
    }
    #[test]
    fn test_iges_circular_arc_radius() {
        let arc = IgesCircularArc {
            z_displacement: 0.0,
            center: [0.0, 0.0],
            start: [3.0, 4.0],
            end: [5.0, 0.0],
        };
        assert!((arc.radius() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_iges_circular_arc_center_3d() {
        let arc = IgesCircularArc {
            z_displacement: 5.0,
            center: [1.0, 2.0],
            start: [0.0, 0.0],
            end: [0.0, 0.0],
        };
        let c = arc.center_3d();
        assert!((c[2] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_assembly_transform_identity() {
        let t = AssemblyTransform::identity();
        let p = t.apply([1.0, 2.0, 3.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_assembly_transform_translation() {
        let t = AssemblyTransform::translation([10.0, 20.0, 30.0]);
        let p = t.apply([1.0, 2.0, 3.0]);
        assert!((p[0] - 11.0).abs() < 1e-10);
        assert!((p[1] - 22.0).abs() < 1e-10);
    }
    #[test]
    fn test_assembly_transform_inverse() {
        let t = AssemblyTransform::translation([5.0, 0.0, 0.0]);
        let inv = t.inverse();
        let p = t.apply([0.0; 3]);
        let p_back = inv.apply(p);
        for (i, &v) in p_back.iter().enumerate() {
            assert!(v.abs() < 1e-10, "p_back[{i}]={v}");
        }
    }
    #[test]
    fn test_assembly_transform_compose() {
        let t1 = AssemblyTransform::translation([1.0, 0.0, 0.0]);
        let t2 = AssemblyTransform::translation([0.0, 2.0, 0.0]);
        let composed = t1.compose(&t2);
        let p = composed.apply([0.0; 3]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_assembly_component_count() {
        let solid = BrepSolid::new("root");
        let mut root = AssemblyComponent::new("root", solid, AssemblyTransform::identity());
        let child1 = AssemblyComponent::new(
            "child1",
            BrepSolid::new("c1"),
            AssemblyTransform::identity(),
        );
        let child2 = AssemblyComponent::new(
            "child2",
            BrepSolid::new("c2"),
            AssemblyTransform::identity(),
        );
        root.add_child(child1);
        root.add_child(child2);
        assert_eq!(root.total_components(), 3);
    }
    #[test]
    fn test_assembly_leaf_names() {
        let solid = BrepSolid::new("root");
        let mut root = AssemblyComponent::new("root", solid, AssemblyTransform::identity());
        root.add_child(AssemblyComponent::new(
            "leaf_a",
            BrepSolid::new("a"),
            AssemblyTransform::identity(),
        ));
        root.add_child(AssemblyComponent::new(
            "leaf_b",
            BrepSolid::new("b"),
            AssemblyTransform::identity(),
        ));
        let names = root.leaf_names();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&"leaf_a".to_string()));
        assert!(names.contains(&"leaf_b".to_string()));
    }
    #[test]
    fn test_assembly_bbox() {
        let solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        let root = AssemblyComponent::new(
            "root",
            solid,
            AssemblyTransform::translation([5.0, 0.0, 0.0]),
        );
        let assembly = Assembly::new("test_asm", root, LengthUnit::Millimeter);
        let bb = assembly.bounding_box();
        assert!((bb.min[0] - 5.0).abs() < 1e-10);
        assert!((bb.max[0] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_assembly_stl_export() {
        let solid = BrepSolid::create_box("box", 1.0, 1.0, 1.0);
        let root = AssemblyComponent::new("root", solid, AssemblyTransform::identity());
        let assembly = Assembly::new("test", root, LengthUnit::Millimeter);
        let stl = assembly.to_stl();
        assert!(stl.contains("solid test"));
        assert!(stl.contains("endsolid test"));
    }
}
