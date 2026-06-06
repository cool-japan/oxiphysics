//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{VtuAnnotation, VtuCellQuality, VtuWriter};

/// VTK cell type: vertex (1 point).
pub const VTK_VERTEX: u8 = 1;
/// VTK cell type: poly vertex (variable points).
pub const VTK_POLY_VERTEX: u8 = 2;
/// VTK cell type: line (2 points).
pub const VTK_LINE: u8 = 3;
/// VTK cell type: triangle (3 points).
pub const VTK_TRIANGLE: u8 = 5;
/// VTK cell type: triangle strip (variable points).
pub const VTK_TRIANGLE_STRIP: u8 = 6;
/// VTK cell type: polygon (variable points).
pub const VTK_POLYGON: u8 = 7;
/// VTK cell type: quad (4 points).
pub const VTK_QUAD: u8 = 9;
/// VTK cell type: tetrahedron (4 points).
pub const VTK_TET: u8 = 10;
/// VTK cell type: voxel (8 points, axis-aligned hex).
pub const VTK_VOXEL: u8 = 11;
/// VTK cell type: hexahedron (8 points).
pub const VTK_HEX: u8 = 12;
/// VTK cell type: wedge / pentahedron (6 points).
pub const VTK_WEDGE: u8 = 13;
/// VTK cell type: pyramid (5 points).
pub const VTK_PYRAMID: u8 = 14;
/// Return the expected number of nodes for a given VTK cell type.
/// Returns 0 for variable-size cells (polygon, poly-vertex, etc.).
pub fn cell_type_num_nodes(cell_type: u8) -> usize {
    match cell_type {
        VTK_VERTEX => 1,
        VTK_POLY_VERTEX => 0,
        VTK_LINE => 2,
        VTK_TRIANGLE => 3,
        VTK_TRIANGLE_STRIP => 0,
        VTK_POLYGON => 0,
        VTK_QUAD => 4,
        VTK_TET => 4,
        VTK_VOXEL => 8,
        VTK_HEX => 8,
        VTK_WEDGE => 6,
        VTK_PYRAMID => 5,
        _ => 0,
    }
}
/// VTK cell type: quadratic triangle (6 points).
pub const VTK_QUADRATIC_TRIANGLE: u8 = 22;
/// VTK cell type: quadratic quad (8 points).
pub const VTK_QUADRATIC_QUAD: u8 = 23;
/// VTK cell type: quadratic tetrahedron (10 points).
pub const VTK_QUADRATIC_TET: u8 = 24;
/// VTK cell type: quadratic hexahedron (20 points).
pub const VTK_QUADRATIC_HEX: u8 = 25;
/// Return the number of nodes for common higher-order VTK cell types.
pub fn higher_order_cell_num_nodes(cell_type: u8) -> usize {
    match cell_type {
        VTK_QUADRATIC_TRIANGLE => 6,
        VTK_QUADRATIC_QUAD => 8,
        VTK_QUADRATIC_TET => 10,
        VTK_QUADRATIC_HEX => 20,
        other => cell_type_num_nodes(other),
    }
}
/// Compute per-cell quality metrics for all cells in a `VtuWriter`.
pub fn compute_cell_quality(writer: &VtuWriter) -> Vec<VtuCellQuality> {
    let mut result = Vec::new();
    let pts = writer.points();
    let types = writer.cell_types();
    let cells = writer.cells_per_cell();
    for (ci, (&ct, node_ids)) in types.iter().zip(cells.iter()).enumerate() {
        let positions: Vec<[f64; 3]> = node_ids
            .iter()
            .filter_map(|&ni| if ni < pts.len() { Some(pts[ni]) } else { None })
            .collect();
        if positions.len() < 2 {
            result.push(VtuCellQuality {
                cell_index: ci,
                cell_type: ct,
                aspect_ratio: 0.0,
                min_edge: 0.0,
                max_edge: 0.0,
                scaled_jacobian: 0.0,
            });
            continue;
        }
        let mut edges = Vec::new();
        for i in 0..positions.len() {
            let j = (i + 1) % positions.len();
            let p = positions[i];
            let q = positions[j];
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            edges.push(d);
        }
        let min_e = edges.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_e = edges.iter().cloned().fold(0.0_f64, f64::max);
        let ar = if max_e > 1e-30 { min_e / max_e } else { 0.0 };
        let sj = if positions.len() >= 4 {
            let e1 = sub3(positions[1], positions[0]);
            let e2 = sub3(positions[2], positions[0]);
            let e3 = sub3(positions[3], positions[0]);
            let jdet = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
                - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
                + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
            let l1 = norm3(e1);
            let l2 = norm3(e2);
            let l3 = norm3(e3);
            let den = l1 * l2 * l3;
            if den > 1e-30 { jdet / den } else { 0.0 }
        } else {
            ar
        };
        result.push(VtuCellQuality {
            cell_index: ci,
            cell_type: ct,
            aspect_ratio: ar,
            min_edge: min_e,
            max_edge: max_e,
            scaled_jacobian: sj,
        });
    }
    result
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Convert a structured grid (Nx × Ny × Nz) to VTU unstructured HEX8 mesh.
///
/// `coords` is a flat array of `(Nx+1) × (Ny+1) × (Nz+1)` node positions
/// in (x, y, z) row-major order (z-fastest).
pub fn structured_to_vtu(nx: usize, ny: usize, nz: usize, coords: &[[f64; 3]]) -> VtuWriter {
    let mut writer = VtuWriter::new();
    writer.add_points(coords);
    let idx =
        |i: usize, j: usize, k: usize| -> usize { i * (ny + 1) * (nz + 1) + j * (nz + 1) + k };
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let n = [
                    idx(i, j, k),
                    idx(i + 1, j, k),
                    idx(i + 1, j + 1, k),
                    idx(i, j + 1, k),
                    idx(i, j, k + 1),
                    idx(i + 1, j, k + 1),
                    idx(i + 1, j + 1, k + 1),
                    idx(i, j + 1, k + 1),
                ];
                writer.add_cell(n.to_vec(), VTK_HEX);
            }
        }
    }
    writer
}
/// Generate a uniform Cartesian structured grid as VTU.
///
/// Domain: \[x0, x1\] × \[y0, y1\] × \[z0, z1\] with `nx`×`ny`×`nz` cells.
pub fn uniform_grid_to_vtu(
    nx: usize,
    ny: usize,
    nz: usize,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
) -> VtuWriter {
    let dx = (x1 - x0) / nx as f64;
    let dy = (y1 - y0) / ny as f64;
    let dz = (z1 - z0) / nz as f64;
    let mut coords = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
    for i in 0..=(nx) {
        for j in 0..=(ny) {
            for k in 0..=(nz) {
                coords.push([x0 + i as f64 * dx, y0 + j as f64 * dy, z0 + k as f64 * dz]);
            }
        }
    }
    structured_to_vtu(nx, ny, nz, &coords)
}
/// Interpolate point (node) scalar data to cell-centre values by averaging.
///
/// For each cell the value is the arithmetic mean of its node values.
/// Returns a vec with one value per cell.
pub fn point_to_cell_data(writer: &VtuWriter, point_data: &[f64]) -> Vec<f64> {
    let cells = writer.cells_per_cell();
    let mut result = Vec::with_capacity(cells.len());
    for node_ids in cells {
        if node_ids.is_empty() {
            result.push(0.0);
            continue;
        }
        let sum: f64 = node_ids
            .iter()
            .map(|&ni| *point_data.get(ni).unwrap_or(&0.0))
            .sum();
        result.push(sum / node_ids.len() as f64);
    }
    result
}
/// Interpolate cell-centre scalar data to node values by averaging over
/// all cells that share each node (unweighted).
///
/// Returns a vec with one value per point.
pub fn cell_to_point_data(writer: &VtuWriter, cell_data: &[f64]) -> Vec<f64> {
    let n_pts = writer.num_points();
    let mut sums = vec![0.0_f64; n_pts];
    let mut counts = vec![0usize; n_pts];
    let cells = writer.cells_per_cell();
    for (ci, node_ids) in cells.iter().enumerate() {
        let cv = *cell_data.get(ci).unwrap_or(&0.0);
        for &ni in node_ids {
            if ni < n_pts {
                sums[ni] += cv;
                counts[ni] += 1;
            }
        }
    }
    sums.iter()
        .zip(counts.iter())
        .map(|(s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
        .collect()
}
/// Estimate the gradient of a scalar field at each node using a weighted
/// sum over node neighbours (based on connectivity).
///
/// This is a simple O(n*k) implementation for unstructured meshes.
/// Returns interleaved `[gx, gy, gz, gx, gy, gz, ...]`.
pub fn estimate_gradient(writer: &VtuWriter, scalars: &[f64]) -> Vec<f64> {
    let n_pts = writer.num_points();
    let pts = writer.points();
    let cells = writer.cells_per_cell();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_pts];
    for node_ids in cells {
        for &ni in node_ids {
            for &nj in node_ids {
                if ni != nj && ni < n_pts && nj < n_pts && !adj[ni].contains(&nj) {
                    adj[ni].push(nj);
                }
            }
        }
    }
    let mut grad = vec![0.0_f64; n_pts * 3];
    for i in 0..n_pts {
        if i >= scalars.len() {
            continue;
        }
        let pi = pts[i];
        let si = scalars[i];
        let mut gx = 0.0_f64;
        let mut gy = 0.0_f64;
        let mut gz = 0.0_f64;
        let mut w_sum = 0.0_f64;
        for &j in &adj[i] {
            if j >= scalars.len() {
                continue;
            }
            let pj = pts[j];
            let sj = scalars[j];
            let dx = pj[0] - pi[0];
            let dy = pj[1] - pi[1];
            let dz = pj[2] - pi[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < 1e-30 {
                continue;
            }
            let ds = sj - si;
            let w = 1.0 / d2;
            gx += w * ds * dx / d2;
            gy += w * ds * dy / d2;
            gz += w * ds * dz / d2;
            w_sum += w;
        }
        if w_sum > 1e-30 {
            grad[i * 3] = gx / w_sum;
            grad[i * 3 + 1] = gy / w_sum;
            grad[i * 3 + 2] = gz / w_sum;
        }
    }
    grad
}
/// Write a VTU file with metadata annotations in the XML header comment.
pub fn write_vtu_with_annotations(writer: &VtuWriter, annotations: &[VtuAnnotation]) -> String {
    let xml = writer.to_xml();
    let mut header = String::from("<!-- VTU Annotations\n");
    for ann in annotations {
        header.push_str(&format!("  {}: {}\n", ann.key, ann.value));
    }
    header.push_str("-->\n");
    format!("{}{}", header, xml)
}
/// Compute the axis-aligned bounding box of all points.
///
/// Returns `(min, max)` or `None` if there are no points.
pub fn compute_bounding_box(points: &[[f64; 3]]) -> Option<([f64; 3], [f64; 3])> {
    if points.is_empty() {
        return None;
    }
    let mut min = points[0];
    let mut max = points[0];
    for &p in &points[1..] {
        for k in 0..3 {
            if p[k] < min[k] {
                min[k] = p[k];
            }
            if p[k] > max[k] {
                max[k] = p[k];
            }
        }
    }
    Some((min, max))
}
/// Compute the diagonal length of a bounding box.
pub fn bbox_diagonal(min: [f64; 3], max: [f64; 3]) -> f64 {
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Compute the volume (Lx × Ly × Lz) of a bounding box.
pub fn bbox_volume(min: [f64; 3], max: [f64; 3]) -> f64 {
    (max[0] - min[0]).abs() * (max[1] - min[1]).abs() * (max[2] - min[2]).abs()
}
/// Filter cells whose average scalar value (computed from point data) is
/// above a threshold.  Returns the point indices of qualifying cells.
pub fn threshold_cells(writer: &VtuWriter, point_data: &[f64], threshold: f64) -> Vec<usize> {
    let cell_vals = point_to_cell_data(writer, point_data);
    cell_vals
        .iter()
        .enumerate()
        .filter(|&(_, &v)| v > threshold)
        .map(|(i, _)| i)
        .collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::vtu::PvdWriter;
    use crate::vtu::PvtuWriter;
    use crate::vtu::VtuReader;
    use crate::vtu::VtuTimeSeries;
    use crate::vtu::VtuWriterLegacy;
    use oxiphysics_core::math::Vec3;
    #[test]
    fn test_vtu_write() {
        let path = std::env::temp_dir().join("oxiphy_test.vtu");
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let cells = vec![[0, 1, 2, 3]];
        let pvals = vec![100.0, 200.0, 300.0, 400.0];
        let cvals = vec![0.166];
        let pdata = vec![("temperature", pvals.as_slice())];
        let cdata = vec![("volume", cvals.as_slice())];
        VtuWriterLegacy::write(path.to_str().unwrap_or(""), &pts, &cells, &pdata, &cdata).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<VTKFile"));
        assert!(content.contains("UnstructuredGrid"));
        assert!(content.contains("temperature"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_vtu_builder_empty_grid_xml() {
        let writer = VtuWriter::new();
        let xml = writer.to_xml();
        assert!(xml.contains("<VTKFile"), "missing VTKFile tag: {xml}");
        assert!(xml.contains("NumberOfPoints=\"0\""));
        assert!(xml.contains("NumberOfCells=\"0\""));
    }
    #[test]
    fn test_vtu_builder_single_point() {
        let mut writer = VtuWriter::new();
        writer.add_point([1.0, 2.0, 3.0]);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfPoints=\"1\""));
        assert!(
            xml.contains("1 2 3") || xml.contains("1.0 2.0 3.0"),
            "point coords missing: {xml}"
        );
    }
    #[test]
    fn test_vtu_builder_triangle_cell() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point([0.0, 1.0, 0.0]);
        writer.add_cell(vec![0, 1, 2], VTK_TRIANGLE);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfCells=\"1\""));
        assert!(xml.contains('5'), "triangle cell type 5 missing: {xml}");
    }
    #[test]
    fn test_vtu_builder_point_data_included() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point_data_scalar("pressure", &[101.0, 202.0]);
        let xml = writer.to_xml();
        assert!(xml.contains("pressure"), "scalar field name missing: {xml}");
        assert!(xml.contains("101"), "scalar value missing: {xml}");
    }
    #[test]
    fn test_vtu_builder_vector_data_included() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point_data_vector("velocity", &[[1.0, 0.0, 0.0]]);
        let xml = writer.to_xml();
        assert!(xml.contains("velocity"), "vector field name missing: {xml}");
        assert!(xml.contains("NumberOfComponents=\"3\""));
    }
    #[test]
    fn test_vtu_reader_parses_points() {
        let mut writer = VtuWriter::new();
        writer.add_point([1.5, 2.5, 3.5]);
        writer.add_point([4.0, 5.0, 6.0]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.points().len(), 2, "expected 2 points");
        assert!((reader.points()[0][0] - 1.5).abs() < 1e-6);
        assert!((reader.points()[1][2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_pvd_format() {
        let pvd = PvdWriter::write(&[(0.0, "frame_0.vtu"), (0.1, "frame_1.vtu")]);
        assert!(pvd.contains("Collection"), "PVD missing Collection: {pvd}");
        assert!(pvd.contains("frame_0.vtu"));
        assert!(pvd.contains("timestep=\"0\"") || pvd.contains("timestep=\"0.0\""));
    }
    #[test]
    fn test_vtu_write_to_file() {
        let path = std::env::temp_dir().join("oxiphy_builder_test.vtu");
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_point([1.0, 0.0, 0.0]);
        writer.add_point([0.0, 1.0, 0.0]);
        writer.add_cell(vec![0, 1, 2], VTK_TRIANGLE);
        writer.write_to_file(path.to_str().unwrap_or("")).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("<VTKFile"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_cell_type_num_nodes() {
        assert_eq!(cell_type_num_nodes(VTK_VERTEX), 1);
        assert_eq!(cell_type_num_nodes(VTK_LINE), 2);
        assert_eq!(cell_type_num_nodes(VTK_TRIANGLE), 3);
        assert_eq!(cell_type_num_nodes(VTK_QUAD), 4);
        assert_eq!(cell_type_num_nodes(VTK_TET), 4);
        assert_eq!(cell_type_num_nodes(VTK_HEX), 8);
        assert_eq!(cell_type_num_nodes(VTK_WEDGE), 6);
        assert_eq!(cell_type_num_nodes(VTK_PYRAMID), 5);
        assert_eq!(cell_type_num_nodes(VTK_POLYGON), 0);
        assert_eq!(cell_type_num_nodes(255), 0);
    }
    #[test]
    fn test_add_points_bulk() {
        let mut writer = VtuWriter::new();
        let first = writer.add_points(&[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(first, 0);
        assert_eq!(writer.num_points(), 3);
    }
    #[test]
    fn test_add_triangle_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains('5'));
    }
    #[test]
    fn test_add_quad_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_quad(0, 1, 2, 3);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains('9'));
    }
    #[test]
    fn test_add_tet_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_tet(0, 1, 2, 3);
        assert_eq!(writer.num_cells(), 1);
    }
    #[test]
    fn test_add_hex_shorthand() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 8]);
        writer.add_hex(0, 1, 2, 3, 4, 5, 6, 7);
        assert_eq!(writer.num_cells(), 1);
        let xml = writer.to_xml();
        assert!(xml.contains("12"));
    }
    #[test]
    fn test_cell_data_scalar_in_xml() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("volume", &[0.5]);
        let xml = writer.to_xml();
        assert!(
            xml.contains("<CellData>"),
            "CellData section missing: {xml}"
        );
        assert!(xml.contains("volume"), "cell scalar name missing: {xml}");
        assert!(xml.contains("0.5"), "cell scalar value missing: {xml}");
    }
    #[test]
    fn test_cell_data_vector_in_xml() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_vector("normal", &[[0.0, 0.0, 1.0]]);
        let xml = writer.to_xml();
        assert!(xml.contains("normal"));
        assert!(xml.contains("<CellData>"));
    }
    #[test]
    fn test_bounding_box_empty() {
        let writer = VtuWriter::new();
        assert!(writer.bounding_box().is_none());
    }
    #[test]
    fn test_bounding_box_multiple_points() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[1.0, -2.0, 3.0], [-1.0, 5.0, 0.0], [0.0, 0.0, 10.0]]);
        let (min, max) = writer.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-12);
        assert!((min[1] - (-2.0)).abs() < 1e-12);
        assert!((min[2] - 0.0).abs() < 1e-12);
        assert!((max[0] - 1.0).abs() < 1e-12);
        assert!((max[1] - 5.0).abs() < 1e-12);
        assert!((max[2] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_validate_valid_writer() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        writer.add_triangle(0, 1, 2);
        writer.add_point_data_scalar("temp", &[1.0, 2.0, 3.0]);
        writer.add_cell_data_scalar("vol", &[0.5]);
        assert!(writer.validate().is_ok());
    }
    #[test]
    fn test_validate_bad_connectivity() {
        let mut writer = VtuWriter::new();
        writer.add_point([0.0, 0.0, 0.0]);
        writer.add_cell(vec![0, 99], VTK_LINE);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_validate_mismatched_point_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_point_data_scalar("temp", &[1.0, 2.0]);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_validate_mismatched_cell_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("vol", &[0.5, 1.0]);
        assert!(writer.validate().is_err());
    }
    #[test]
    fn test_reader_extracts_num_cells() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell(vec![2, 3], VTK_LINE);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.num_cells(), 2);
    }
    #[test]
    fn test_reader_cell_types() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 4]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell(vec![2, 3], VTK_LINE);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.cell_types().len(), 2);
        assert_eq!(reader.cell_types()[0], VTK_TRIANGLE);
        assert_eq!(reader.cell_types()[1], VTK_LINE);
    }
    #[test]
    fn test_reader_point_data_names() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 2]);
        writer.add_point_data_scalar("pressure", &[1.0, 2.0]);
        writer.add_point_data_scalar("temperature", &[300.0, 310.0]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let names = reader.point_data_names();
        assert!(names.iter().any(|n| n == "pressure"));
        assert!(names.iter().any(|n| n == "temperature"));
    }
    #[test]
    fn test_reader_cell_data_names() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 3]);
        writer.add_triangle(0, 1, 2);
        writer.add_cell_data_scalar("volume", &[0.5]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let names = reader.cell_data_names();
        assert!(names.iter().any(|n| n == "volume"));
    }
    #[test]
    fn test_reader_bounding_box() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[-1.0, -2.0, -3.0], [4.0, 5.0, 6.0]]);
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        let (min, max) = reader.bounding_box().unwrap();
        assert!((min[0] - (-1.0)).abs() < 1e-6);
        assert!((max[2] - 6.0).abs() < 1e-6);
    }
    #[test]
    fn test_pvd_parallel() {
        let pvd = PvdWriter::write_parallel(&[
            (0.0, vec!["frame_0_0.vtu", "frame_0_1.vtu"]),
            (0.1, vec!["frame_1_0.vtu", "frame_1_1.vtu"]),
        ]);
        assert!(pvd.contains("frame_0_0.vtu"));
        assert!(pvd.contains("frame_0_1.vtu"));
        assert!(pvd.contains("part=\"0\""));
        assert!(pvd.contains("part=\"1\""));
    }
    #[test]
    fn test_pvtu_writer() {
        let pvtu = PvtuWriter::write(
            &["piece_0.vtu", "piece_1.vtu"],
            &["temperature"],
            &["volume"],
        );
        assert!(pvtu.contains("PUnstructuredGrid"));
        assert!(pvtu.contains("piece_0.vtu"));
        assert!(pvtu.contains("piece_1.vtu"));
        assert!(pvtu.contains("temperature"));
        assert!(pvtu.contains("volume"));
        assert!(pvtu.contains("PPoints"));
    }
    #[test]
    fn test_roundtrip_write_read_with_cell_data() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ]);
        writer.add_tet(0, 1, 2, 3);
        writer.add_point_data_scalar("density", &[1.0, 1.1, 1.2, 1.3]);
        writer.add_cell_data_scalar("quality", &[0.95]);
        assert!(writer.validate().is_ok());
        let xml = writer.to_xml();
        let reader = VtuReader::from_xml(&xml).unwrap();
        assert_eq!(reader.points().len(), 4);
        assert_eq!(reader.num_cells(), 1);
        assert!(reader.point_data_names().iter().any(|n| n == "density"));
        assert!(reader.cell_data_names().iter().any(|n| n == "quality"));
    }
    #[test]
    fn test_default_trait() {
        let w = VtuWriter::default();
        assert_eq!(w.num_points(), 0);
        assert_eq!(w.num_cells(), 0);
    }
    #[test]
    fn test_multiple_cells_mixed_types() {
        let mut writer = VtuWriter::new();
        writer.add_points(&[[0.0; 3]; 10]);
        writer.add_triangle(0, 1, 2);
        writer.add_quad(3, 4, 5, 6);
        writer.add_cell(vec![7, 8], VTK_LINE);
        writer.add_cell(vec![9], VTK_VERTEX);
        assert_eq!(writer.num_cells(), 4);
        let xml = writer.to_xml();
        assert!(xml.contains("NumberOfCells=\"4\""));
    }
    #[test]
    fn test_time_series_num_frames() {
        let ts = VtuTimeSeries::new("/tmp", "test_ts");
        assert_eq!(ts.num_frames(), 0);
        assert!(ts.times().is_empty());
    }
}
