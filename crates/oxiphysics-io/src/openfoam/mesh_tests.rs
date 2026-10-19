// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the structured box mesh (face winding and ordering), the polyMesh
//! writers and the Green-Gauss gradient.

use super::functions::{
    face_area, face_area_vector, face_centroid, face_normal, foam_header, foam_header_with_note,
    green_gauss_gradient, parse_foam_faces, parse_foam_label_list, parse_foam_points,
};
use super::types::FoamMesh;
use crate::openfoam_format;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Cell centres of `FoamMesh::box_mesh(l.., n..)`, in cell-index order.
fn box_centres(l: [f64; 3], n: [usize; 3]) -> Vec<[f64; 3]> {
    let mut centres = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                centres.push([
                    (i as f64 + 0.5) * l[0] / n[0] as f64,
                    (j as f64 + 0.5) * l[1] / n[1] as f64,
                    (k as f64 + 0.5) * l[2] / n[2] as f64,
                ]);
            }
        }
    }
    centres
}

fn face_vertices(mesh: &FoamMesh, face: usize) -> Vec<[f64; 3]> {
    mesh.faces[face].iter().map(|&v| mesh.points[v]).collect()
}

fn gradient(mesh: &FoamMesh, values: &[f64], centres: &[[f64; 3]]) -> Vec<[f64; 3]> {
    green_gauss_gradient(
        values,
        centres,
        &mesh.faces,
        &mesh.points,
        &mesh.owner,
        &mesh.neighbour,
        mesh.n_cells,
    )
}

const MESHES: [([f64; 3], [usize; 3]); 3] = [
    ([1.0, 2.0, 3.0], [3, 3, 3]),
    ([2.0, 3.0, 4.0], [4, 3, 2]),
    ([1.0, 1.0, 1.0], [1, 2, 3]),
];

#[test]
fn internal_face_normals_point_from_owner_to_neighbour() {
    for (l, n) in MESHES {
        let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
        let centres = box_centres(l, n);
        let n_internal = mesh.n_internal_faces();
        assert!(n_internal > 0);
        for face in 0..n_internal {
            let normal = face_normal(&face_vertices(&mesh, face));
            let owner = mesh.owner[face];
            let neighbour = mesh.neighbour[face] as usize;
            let towards_neighbour = sub(centres[neighbour], centres[owner]);
            assert!(
                dot(normal, towards_neighbour) > 0.0,
                "mesh {n:?}, internal face {face}: normal {normal:?} points back to the owner"
            );
        }
    }
}

#[test]
fn y_faces_have_a_positive_y_normal() {
    let mesh = FoamMesh::box_mesh(1.0, 2.0, 3.0, 3, 3, 3);
    let mut y_faces = 0;
    for face in 0..mesh.n_internal_faces() {
        let normal = face_normal(&face_vertices(&mesh, face));
        if normal[1].abs() > 0.5 {
            y_faces += 1;
            assert!(normal[1] > 0.99, "face {face}: {normal:?}");
        }
    }
    assert_eq!(y_faces, 3 * 2 * 3, "internal y faces of a 3x3x3 box");
}

#[test]
fn boundary_face_normals_point_out_of_the_box() {
    for (l, n) in MESHES {
        let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
        let centres = box_centres(l, n);
        for face in mesh.n_internal_faces()..mesh.n_faces() {
            let verts = face_vertices(&mesh, face);
            let normal = face_normal(&verts);
            let outward = sub(face_centroid(&verts), centres[mesh.owner[face]]);
            assert!(
                dot(normal, outward) > 0.0,
                "mesh {n:?}, boundary face {face}: normal {normal:?} points into the cell"
            );
        }
    }
}

#[test]
fn internal_faces_are_in_upper_triangular_order() {
    for (l, n) in MESHES {
        let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
        let n_internal = mesh.n_internal_faces();
        let mut previous: Option<(usize, usize)> = None;
        for face in 0..n_internal {
            let owner = mesh.owner[face];
            let neighbour = mesh.neighbour[face] as usize;
            assert!(owner < neighbour, "face {face}: {owner} !< {neighbour}");
            if let Some(prev) = previous {
                assert!((owner, neighbour) > prev, "face {face} out of order");
            }
            previous = Some((owner, neighbour));
        }
        assert!(
            mesh.neighbour[n_internal..].iter().all(|&nb| nb == -1),
            "boundary faces must follow the internal ones"
        );
        assert!(mesh.check_topology(), "mesh {n:?}");
    }
}

#[test]
fn internal_face_counts_match_the_structured_grid() {
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
    assert_eq!(mesh.n_internal_faces(), 12);
    assert_eq!(mesh.n_boundary_faces(), 24);
    assert_eq!(mesh.n_faces(), 36);
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 4, 3, 2);
    assert_eq!(mesh.n_internal_faces(), 3 * 3 * 2 + 4 * 2 * 2 + 4 * 3);
    assert_eq!(mesh.n_boundary_faces(), 2 * (3 * 2 + 4 * 2 + 4 * 3));
}

#[test]
fn neighbour_file_lists_only_the_internal_faces() {
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 2, 2);
    let text = mesh.write_neighbour();
    let written = parse_foam_label_list(&text);
    assert_eq!(written.len(), mesh.n_internal_faces());
    assert!(
        written.iter().all(|&n| n >= 0),
        "no -1 entries: {written:?}"
    );
    let internal: Vec<i64> = mesh.neighbour.iter().copied().filter(|&n| n >= 0).collect();
    assert_eq!(written, internal);
    // The leading size is nInternalFaces, not nFaces.
    let count_line = text
        .lines()
        .skip_while(|l| !l.starts_with("// * *"))
        .skip(1)
        .find(|l| !l.trim().is_empty())
        .unwrap();
    assert_eq!(count_line.trim(), mesh.n_internal_faces().to_string());
    assert!(!text.contains("-1"));
}

#[test]
fn owner_file_lists_every_face_and_records_the_mesh_sizes() {
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
    let text = mesh.write_owner();
    assert_eq!(parse_foam_label_list(&text).len(), mesh.n_faces());
    assert!(
        text.contains("note        \"nPoints:27 nCells:8 nFaces:36 nInternalFaces:12\";"),
        "{text}"
    );
    let neighbour = mesh.write_neighbour();
    assert!(neighbour.contains("nInternalFaces:12"));
    // Other files do not get a note.
    assert!(!mesh.write_faces().contains("note"));
    assert!(!mesh.write_points().contains("note"));
}

#[test]
fn header_note_is_optional_and_quoted() {
    let plain = foam_header("labelList", "owner");
    assert!(!plain.contains("note"));
    let noted = foam_header_with_note("labelList", "owner", Some("nPoints:1 nCells:1"));
    assert!(noted.contains("note        \"nPoints:1 nCells:1\";\n"));
    assert!(noted.contains("class       labelList;\n"));
    assert!(noted.contains("object      owner;\n"));
    assert_eq!(foam_header_with_note("labelList", "owner", None), plain);
}

#[test]
fn polymesh_files_read_back_through_the_simplified_reader() {
    for (l, n) in MESHES {
        let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
        let read = openfoam_format::read_foam_mesh(
            &mesh.write_points(),
            &mesh.write_faces(),
            &mesh.write_owner(),
            &mesh.write_neighbour(),
        )
        .unwrap_or_else(|e| panic!("mesh {n:?} is not readable: {e}"));
        assert_eq!(read.n_points(), mesh.n_points());
        assert_eq!(read.n_faces(), mesh.n_faces());
        assert_eq!(read.n_internal_faces(), mesh.n_internal_faces());
        assert_eq!(read.faces, mesh.faces);
        assert_eq!(read.owner, mesh.owner);
        let internal: Vec<usize> = mesh
            .neighbour
            .iter()
            .filter_map(|&nb| usize::try_from(nb).ok())
            .collect();
        assert_eq!(read.neighbour, internal);
        for (read_p, p) in read.points.chunks(3).zip(&mesh.points) {
            assert_eq!(read_p, p.as_slice());
        }
    }
}

#[test]
fn points_and_faces_parse_back_with_the_openfoam_module_readers() {
    let mesh = FoamMesh::box_mesh(1.0, 2.0, 3.0, 2, 2, 2);
    assert_eq!(parse_foam_points(&mesh.write_points()), mesh.points);
    assert_eq!(parse_foam_faces(&mesh.write_faces()), mesh.faces);
}

#[test]
fn check_topology_accepts_both_neighbour_layouts() {
    let mut mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    assert!(mesh.check_topology());
    let n_internal = mesh.n_internal_faces();
    let n_boundary = mesh.n_boundary_faces();
    mesh.neighbour.truncate(n_internal);
    assert!(mesh.check_topology(), "standard OpenFOAM neighbour list");
    assert_eq!(mesh.n_internal_faces(), n_internal);
    assert_eq!(mesh.n_boundary_faces(), n_boundary);
}

#[test]
fn check_topology_rejects_inconsistent_meshes() {
    let base = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    assert!(base.check_topology());

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.owner.pop();
    assert!(!m.check_topology(), "owner shorter than faces");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.neighbour.push(-1);
    assert!(!m.check_topology(), "neighbour longer than faces");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.faces.swap(0, 1);
    m.owner.swap(0, 1);
    m.neighbour.swap(0, 1);
    assert!(!m.check_topology(), "internal faces out of order");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    let (own, nbr) = (m.owner[2], m.neighbour[2]);
    m.owner[2] = nbr as usize;
    m.neighbour[2] = own as i64;
    assert!(!m.check_topology(), "owner above neighbour");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.neighbour[0] = -1;
    assert!(!m.check_topology(), "boundary entry before internal faces");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    let last = m.neighbour.len() - 1;
    m.neighbour[last] = 4;
    assert!(!m.check_topology(), "internal face after boundary faces");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.neighbour[3] = m.n_cells as i64;
    assert!(!m.check_topology(), "neighbour beyond the last cell");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.faces[5][1] = m.points.len();
    assert!(!m.check_topology(), "face refers to a missing point");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.boundary_patches[0].start_face = 0;
    assert!(!m.check_topology(), "patch overlaps the internal faces");

    let mut m = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    m.boundary_patches[5].n_faces += 1;
    assert!(!m.check_topology(), "patch runs past the last face");
}

#[test]
fn degenerate_box_sizes_give_an_empty_mesh() {
    for (nx, ny, nz) in [(0, 2, 2), (2, 0, 2), (2, 2, 0), (0, 0, 0)] {
        let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, nx, ny, nz);
        assert_eq!(mesh.n_cells, 0);
        assert!(mesh.points.is_empty() && mesh.faces.is_empty());
        assert!(mesh.check_topology());
    }
}

#[test]
fn single_cell_box_has_only_boundary_faces() {
    let mesh = FoamMesh::box_mesh(2.0, 3.0, 4.0, 1, 1, 1);
    assert_eq!(mesh.n_internal_faces(), 0);
    assert_eq!(mesh.n_boundary_faces(), 6);
    assert!(mesh.neighbour.iter().all(|&n| n == -1));
    assert!(mesh.check_topology());
    let text = mesh.write_neighbour();
    assert!(parse_foam_label_list(&text).is_empty());
}

fn linear_field(centres: &[[f64; 3]], g: [f64; 3]) -> Vec<f64> {
    centres.iter().map(|c| dot(g, *c)).collect()
}

/// Cells of the box whose six neighbours all exist.
fn interior_cells(n: [usize; 3]) -> Vec<usize> {
    let mut cells = Vec::new();
    for k in 1..n[2].saturating_sub(1) {
        for j in 1..n[1].saturating_sub(1) {
            for i in 1..n[0].saturating_sub(1) {
                cells.push(k * n[1] * n[0] + j * n[0] + i);
            }
        }
    }
    cells
}

#[test]
fn green_gauss_recovers_a_linear_gradient_in_interior_cells() {
    let target = [2.0, 3.0, -1.0];
    for (l, n) in [([1.0, 1.0, 1.0], [4, 4, 4]), ([2.0, 3.0, 4.0], [5, 4, 3])] {
        let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
        let centres = box_centres(l, n);
        let values = linear_field(&centres, target);
        let grad = gradient(&mesh, &values, &centres);
        assert_eq!(grad.len(), mesh.n_cells);
        let interior = interior_cells(n);
        assert!(!interior.is_empty());
        for cell in interior {
            for (axis, (got, want)) in grad[cell].iter().zip(&target).enumerate() {
                assert!(
                    (got - want).abs() < 1e-10,
                    "mesh {n:?} cell {cell} axis {axis}: {got} vs {want}"
                );
            }
        }
    }
}

#[test]
fn green_gauss_gradient_is_not_a_third_of_the_true_value() {
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 3, 3, 3);
    let centres = box_centres([1.0; 3], [3, 3, 3]);
    let values: Vec<f64> = centres.iter().map(|c| c[0]).collect();
    let grad = gradient(&mesh, &values, &centres);
    // The only fully interior cell of a 3x3x3 grid is the middle one.
    let middle = 13;
    assert!((grad[middle][0] - 1.0).abs() < 1e-12, "{:?}", grad[middle]);
    assert!(grad[middle][1].abs() < 1e-12 && grad[middle][2].abs() < 1e-12);
}

#[test]
fn green_gauss_of_a_constant_field_vanishes_everywhere() {
    let mesh = FoamMesh::box_mesh(2.0, 3.0, 4.0, 4, 3, 2);
    let centres = box_centres([2.0, 3.0, 4.0], [4, 3, 2]);
    let values = vec![7.25; mesh.n_cells];
    for g in gradient(&mesh, &values, &centres) {
        assert!(g.iter().all(|c| c.abs() < 1e-12), "{g:?}");
    }
}

#[test]
fn green_gauss_accepts_the_standard_neighbour_list() {
    let (l, n) = ([1.0, 2.0, 1.5], [4, 3, 3]);
    let mesh = FoamMesh::box_mesh(l[0], l[1], l[2], n[0], n[1], n[2]);
    let centres = box_centres(l, n);
    let values = linear_field(&centres, [1.0, -2.0, 0.5]);
    let padded = gradient(&mesh, &values, &centres);
    let standard = green_gauss_gradient(
        &values,
        &centres,
        &mesh.faces,
        &mesh.points,
        &mesh.owner,
        &mesh.neighbour[..mesh.n_internal_faces()],
        mesh.n_cells,
    );
    assert_eq!(padded, standard);
    // And with no neighbour list at all every face is a boundary face.
    let none = green_gauss_gradient(
        &values,
        &centres,
        &mesh.faces,
        &mesh.points,
        &mesh.owner,
        &[],
        mesh.n_cells,
    );
    assert_eq!(none.len(), mesh.n_cells);
}

#[test]
fn green_gauss_tolerates_short_or_invalid_topology() {
    let mesh = FoamMesh::box_mesh(1.0, 1.0, 1.0, 2, 2, 2);
    let centres = box_centres([1.0; 3], [2, 2, 2]);
    let values = vec![1.0; mesh.n_cells];
    // An owner list that stops early skips the remaining faces.
    let short_owner = green_gauss_gradient(
        &values,
        &centres,
        &mesh.faces,
        &mesh.points,
        &mesh.owner[..5],
        &mesh.neighbour,
        mesh.n_cells,
    );
    assert_eq!(short_owner.len(), 8);
    // Out-of-range owner and neighbour cells and unknown points are ignored.
    let mut faces = mesh.faces.clone();
    faces[0] = vec![0, 1, 10_000];
    let mut owner = mesh.owner.clone();
    owner[1] = 999;
    let mut neighbour = mesh.neighbour.clone();
    neighbour[2] = 999;
    let grad = green_gauss_gradient(
        &values,
        &centres,
        &faces,
        &mesh.points,
        &owner,
        &neighbour,
        mesh.n_cells,
    );
    assert!(grad.iter().flatten().all(|c| c.is_finite()));
    // No cells at all.
    assert!(green_gauss_gradient(&[], &[], &[], &[], &[], &[], 0).is_empty());
}

#[test]
fn degenerate_cells_get_a_zero_gradient() {
    let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
    let grad = green_gauss_gradient(
        &[5.0],
        &[[0.5, 0.3, 0.0]],
        &[vec![0, 1, 2]],
        &points,
        &[0],
        &[],
        1,
    );
    assert_eq!(grad, vec![[0.0; 3]]);
}

#[test]
fn face_area_vector_follows_the_vertex_order() {
    let square = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 3.0, 0.0],
        [0.0, 3.0, 0.0],
    ];
    assert_eq!(face_area_vector(&square), [0.0, 0.0, 6.0]);
    let mut reversed = square;
    reversed.reverse();
    assert_eq!(face_area_vector(&reversed), [0.0, 0.0, -6.0]);
    assert!((face_area(&square) - 6.0).abs() < 1e-14);
    assert_eq!(face_area_vector(&square[..2]), [0.0; 3]);
    assert_eq!(face_area_vector(&[]), [0.0; 3]);
}

#[test]
fn face_area_vector_of_a_warped_quad_is_the_fan_sum() {
    // Bilinear patch with one corner lifted: S = 1/2 (d1 x d2) over the diagonals.
    let quad = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.5],
        [0.0, 1.0, 0.0],
    ];
    let sf = face_area_vector(&quad);
    let d1 = sub(quad[2], quad[0]);
    let d2 = sub(quad[3], quad[1]);
    let cross = [
        d1[1] * d2[2] - d1[2] * d2[1],
        d1[2] * d2[0] - d1[0] * d2[2],
        d1[0] * d2[1] - d1[1] * d2[0],
    ];
    for (got, want) in sf.iter().zip(cross.map(|c| 0.5 * c)) {
        assert!((got - want).abs() < 1e-14, "{sf:?}");
    }
}
