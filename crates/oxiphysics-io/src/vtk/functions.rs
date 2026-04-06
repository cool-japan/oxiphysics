//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::io::Write;

use super::types::{PvdEntry, VtkCellType, VtkLegacyData, VtuGrid};

#[cfg(test)]
mod tests {

    use crate::VtkWriter;

    use oxiphysics_core::math::Vec3;
    use std::path::Path;
    #[test]
    fn test_vtk_point_cloud_write() {
        let path = "/tmp/oxiphy_test_points.vtk";
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        VtkWriter::write_points(path, &pts).unwrap();
        assert!(Path::new(path).exists());
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("POINTS 3 float"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_vtk_unstructured_grid_write() {
        let path = "/tmp/oxiphy_test_ugrid.vtk";
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let cells = vec![[0, 1, 2, 3]];
        let scalars = vec![1.0, 2.0, 3.0, 4.0];
        VtkWriter::write_unstructured_grid(path, &pts, &cells, Some(("pressure", &scalars)), None)
            .unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("UNSTRUCTURED_GRID"));
        assert!(content.contains("CELL_TYPES 1"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_vtk_header_format() {
        let path = "/tmp/oxiphy_test_header.vtk";
        let pts = vec![Vec3::new(0.0, 0.0, 0.0)];
        VtkWriter::write_points(path, &pts).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], "# vtk DataFile Version 3.0");
        assert_eq!(lines[2], "ASCII");
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_vtk_polydata_write() {
        let path = "/tmp/oxiphy_test_polydata.vtk";
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        VtkWriter::write_polydata(path, &pts, &tris).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("POLYGONS 1 4"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_vtk_write_read_roundtrip() {
        let path = "/tmp/oxiphy_test_roundtrip.vtk";
        let positions = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
            Vec3::new(7.0, 8.0, 9.0),
            Vec3::new(0.5, 1.5, 2.5),
        ];
        let scalars = vec![10.0, 20.0, 30.0, 40.0];
        let cells: Vec<[usize; 4]> = vec![];
        VtkWriter::write_unstructured_grid(
            path,
            &positions,
            &cells,
            Some(("density", &scalars)),
            None,
        )
        .unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(
            content.contains("POINTS 4 float"),
            "expected POINTS 4 float header"
        );
        assert!(
            content.contains("SCALARS density float 1"),
            "expected scalars section"
        );
        let mut parsed: Vec<[f64; 3]> = Vec::new();
        let mut in_points = false;
        for line in content.lines() {
            if line.starts_with("POINTS") {
                in_points = true;
                continue;
            }
            if in_points {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() == 3 {
                    if let (Ok(x), Ok(y), Ok(z)) = (
                        parts[0].parse::<f64>(),
                        parts[1].parse::<f64>(),
                        parts[2].parse::<f64>(),
                    ) {
                        parsed.push([x, y, z]);
                    }
                } else {
                    in_points = false;
                }
            }
        }
        assert_eq!(parsed.len(), 4, "expected 4 parsed points");
        let tol = 1e-4;
        assert!((parsed[0][0] - 1.0).abs() < tol);
        assert!((parsed[0][1] - 2.0).abs() < tol);
        assert!((parsed[0][2] - 3.0).abs() < tol);
        assert!((parsed[1][0] - 4.0).abs() < tol);
        assert!((parsed[2][1] - 8.0).abs() < tol);
        assert!((parsed[3][2] - 2.5).abs() < tol);
        std::fs::remove_file(path).ok();
    }
}
/// Validate a `VtuGrid` for consistency.
///
/// Checks that all cell connectivity indices are valid and all point data
/// arrays have the correct length.
#[allow(dead_code)]
pub fn validate_vtu_grid(grid: &VtuGrid) -> Vec<String> {
    let mut errors = Vec::new();
    let n_pts = grid.n_points();
    for (ci, conn) in grid.cells.iter().enumerate() {
        for &idx in conn {
            if idx >= n_pts {
                errors.push(format!("Cell {ci}: point index {idx} >= n_points {n_pts}"));
            }
        }
    }
    for arr in &grid.point_data {
        if arr.len() != n_pts {
            errors.push(format!(
                "Point data '{}': len {} != n_points {}",
                arr.name(),
                arr.len(),
                n_pts
            ));
        }
    }
    let n_cells = grid.n_cells();
    for arr in &grid.cell_data {
        if arr.len() != n_cells {
            errors.push(format!(
                "Cell data '{}': len {} != n_cells {}",
                arr.name(),
                arr.len(),
                n_cells
            ));
        }
    }
    errors
}
#[cfg(test)]
mod vtu_grid_tests {
    use super::*;

    use crate::vtk::VtkPolyDataGrid;
    use crate::vtk::VtkRectilinearGrid;
    use crate::vtk::VtkStructuredGrid;
    use crate::vtk::VtkTimeSeries;
    use crate::vtk::types::*;

    #[test]
    fn test_vtu_empty_grid() {
        let grid = VtuGrid::new();
        let xml = grid.to_vtu_string();
        assert!(
            xml.starts_with("<?xml version=\"1.0\"?>"),
            "missing XML header"
        );
        assert!(xml.contains("<VTKFile"), "missing VTKFile element");
        assert!(xml.contains("UnstructuredGrid"), "missing UnstructuredGrid");
    }
    #[test]
    fn test_vtu_single_point() {
        let mut grid = VtuGrid::new();
        grid.add_point([1.0, 2.0, 3.0]);
        let xml = grid.to_vtu_string();
        assert!(
            xml.contains("NumberOfPoints=\"1\""),
            "expected NumberOfPoints=\"1\""
        );
    }
    #[test]
    fn test_vtu_tetrahedron() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point([0.0, 1.0, 0.0]);
        grid.add_point([0.0, 0.0, 1.0]);
        grid.add_cell(vec![0, 1, 2, 3], VtkCellType::Tetra);
        let xml = grid.to_vtu_string();
        assert!(
            xml.contains(">10<") || xml.contains("          10"),
            "expected cell type 10 in output; got:\n{}",
            xml
        );
    }
    #[test]
    fn test_vtu_point_data() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point([0.0, 1.0, 0.0]);
        grid.add_cell(vec![0, 1, 2], VtkCellType::Triangle);
        grid.add_point_scalar("pressure", vec![1.0, 2.0, 3.0]);
        let xml = grid.to_vtu_string();
        assert!(
            xml.contains("Name=\"pressure\""),
            "pressure field missing from PointData"
        );
        assert!(xml.contains("<PointData>"), "PointData section missing");
    }
    #[test]
    fn test_vtu_from_points() {
        let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let grid = VtuGrid::from_points(&pts);
        assert_eq!(grid.n_points(), 3);
        assert_eq!(grid.n_cells(), 3);
        let xml = grid.to_vtu_string();
        assert!(xml.contains("NumberOfCells=\"3\""));
    }
    #[test]
    fn test_vtu_from_tet_mesh() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let elems = [[0, 1, 2, 3]];
        let grid = VtuGrid::from_tet_mesh(&nodes, &elems);
        assert_eq!(grid.n_points(), 4);
        assert_eq!(grid.n_cells(), 1);
        let xml = grid.to_vtu_string();
        assert!(xml.contains("NumberOfCells=\"1\""));
        assert!(xml.contains("10"), "expected cell type 10 for tet");
    }
    #[test]
    fn test_pvd_collection() {
        let steps = vec![
            (0.0_f64, "sim_0000.vtu".to_owned()),
            (1.0_f64, "sim_0001.vtu".to_owned()),
        ];
        let pvd = VtuGrid::write_pvd_collection("simulation", &steps);
        assert!(
            pvd.contains("timestep=\"0\"")
                || pvd.contains("timestep=\"0.0\"")
                || pvd.contains("timestep=\"0\""),
            "missing first timestep entry"
        );
        let count = pvd.matches("<DataSet").count();
        assert_eq!(count, 2, "expected 2 DataSet entries, got {}", count);
    }
    #[test]
    fn test_vtu_offsets() {
        let mut grid = VtuGrid::new();
        for i in 0..8 {
            grid.add_point([i as f64, 0.0, 0.0]);
        }
        grid.add_cell(vec![0, 1, 2, 3], VtkCellType::Tetra);
        grid.add_cell(vec![4, 5, 6, 7], VtkCellType::Tetra);
        let xml = grid.to_vtu_string();
        assert!(
            xml.contains("4 8"),
            "expected offsets '4 8' in output;\n{}",
            xml
        );
    }
    #[test]
    fn test_structured_grid_uniform_point_count() {
        let g = VtkStructuredGrid::uniform(0.0, 1.0, 3, 0.0, 1.0, 3, 0.0, 1.0, 3);
        assert_eq!(g.n_points(), 27);
        assert_eq!(g.points.len(), 27);
    }
    #[test]
    fn test_structured_grid_single_point() {
        let g = VtkStructuredGrid::uniform(0.0, 0.0, 1, 0.0, 0.0, 1, 0.0, 0.0, 1);
        assert_eq!(g.n_points(), 1);
    }
    #[test]
    fn test_structured_grid_to_vts_contains_header() {
        let g = VtkStructuredGrid::uniform(0.0, 1.0, 2, 0.0, 1.0, 2, 0.0, 1.0, 2);
        let s = g.to_vts_string();
        assert!(s.contains("<VTKFile type=\"StructuredGrid\""));
        assert!(s.contains("<Piece"));
    }
    #[test]
    fn test_structured_grid_with_scalar() {
        let mut g = VtkStructuredGrid::uniform(0.0, 1.0, 2, 0.0, 1.0, 2, 0.0, 1.0, 2);
        let vals: Vec<f64> = (0..g.n_points()).map(|i| i as f64).collect();
        g.add_point_scalar("index", vals);
        let s = g.to_vts_string();
        assert!(s.contains("Name=\"index\""));
    }
    #[test]
    fn test_rectilinear_grid_n_points() {
        let g = VtkRectilinearGrid::new(vec![0.0, 0.5, 1.0], vec![0.0, 1.0], vec![0.0]);
        assert_eq!(g.n_points(), 6);
    }
    #[test]
    fn test_rectilinear_grid_dims() {
        let g = VtkRectilinearGrid::new(vec![0.0, 1.0], vec![0.0, 1.0, 2.0], vec![0.0, 0.5, 1.0]);
        assert_eq!(g.dims(), [2, 3, 3]);
    }
    #[test]
    fn test_rectilinear_grid_vtr_string() {
        let g = VtkRectilinearGrid::new(vec![0.0, 1.0], vec![0.0, 1.0], vec![0.0, 1.0]);
        let s = g.to_vtr_string();
        assert!(s.contains("<VTKFile type=\"RectilinearGrid\""));
        assert!(s.contains("<Coordinates>"));
        assert!(s.contains("Name=\"x\""));
        assert!(s.contains("Name=\"y\""));
        assert!(s.contains("Name=\"z\""));
    }
    #[test]
    fn test_polydata_add_triangle() {
        let mut g = VtkPolyDataGrid::new();
        g.add_point([0.0, 0.0, 0.0]);
        g.add_point([1.0, 0.0, 0.0]);
        g.add_point([0.0, 1.0, 0.0]);
        g.add_triangle(0, 1, 2);
        assert_eq!(g.n_triangles(), 1);
        assert_eq!(g.n_points(), 3);
    }
    #[test]
    fn test_polydata_triangle_normals() {
        let mut g = VtkPolyDataGrid::new();
        g.add_point([0.0, 0.0, 0.0]);
        g.add_point([1.0, 0.0, 0.0]);
        g.add_point([0.0, 1.0, 0.0]);
        g.add_triangle(0, 1, 2);
        let normals = g.compute_triangle_normals();
        assert_eq!(normals.len(), 1);
        assert!(
            (normals[0][2]).abs() > 0.9,
            "Normal z-component should be ~1"
        );
    }
    #[test]
    fn test_polydata_to_vtu_grid() {
        let mut g = VtkPolyDataGrid::new();
        g.add_point([0.0, 0.0, 0.0]);
        g.add_point([1.0, 0.0, 0.0]);
        g.add_point([0.0, 1.0, 0.0]);
        g.add_triangle(0, 1, 2);
        let vtu = g.to_vtu_grid();
        assert_eq!(vtu.n_points(), 3);
        assert_eq!(vtu.n_cells(), 1);
    }
    #[test]
    fn test_time_series_push_and_count() {
        let mut ts = VtkTimeSeries::new("sim");
        let g1 = VtuGrid::from_points(&[[0.0, 0.0, 0.0]]);
        let g2 = VtuGrid::from_points(&[[1.0, 0.0, 0.0]]);
        ts.push(0.0, g1);
        ts.push(0.1, g2);
        assert_eq!(ts.n_steps(), 2);
    }
    #[test]
    fn test_time_series_pvd_string() {
        let mut ts = VtkTimeSeries::new("run");
        ts.push(0.0, VtuGrid::new());
        ts.push(1.0, VtuGrid::new());
        let pvd = ts.to_pvd_string();
        assert!(
            pvd.contains("<DataSet"),
            "PVD should contain DataSet entries"
        );
        let count = pvd.matches("<DataSet").count();
        assert_eq!(count, 2);
    }
    #[test]
    fn test_time_series_vtu_string() {
        let mut ts = VtkTimeSeries::new("s");
        let mut g = VtuGrid::new();
        g.add_point([1.0, 2.0, 3.0]);
        ts.push(0.5, g);
        let vtu = ts.vtu_string(0).expect("should have step 0");
        assert!(
            vtu.contains("1 2 3") || vtu.contains("1.0 2.0 3.0") || vtu.contains("1 2 3"),
            "VTU should contain coordinates"
        );
    }
    #[test]
    fn test_time_series_none_out_of_range() {
        let ts = VtkTimeSeries::new("s");
        assert!(ts.vtu_string(0).is_none());
    }
    #[test]
    fn test_validate_vtu_grid_valid() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point([0.0, 1.0, 0.0]);
        grid.add_cell(vec![0, 1, 2], VtkCellType::Triangle);
        grid.add_point_scalar("p", vec![1.0, 2.0, 3.0]);
        let errors = validate_vtu_grid(&grid);
        assert!(
            errors.is_empty(),
            "Valid grid should have no errors: {:?}",
            errors
        );
    }
    #[test]
    fn test_validate_vtu_grid_bad_connectivity() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.cells.push(vec![0, 5]);
        grid.cell_types.push(VtkCellType::Line);
        let errors = validate_vtu_grid(&grid);
        assert!(
            !errors.is_empty(),
            "Out-of-range index should produce an error"
        );
    }
    #[test]
    fn test_validate_vtu_grid_mismatched_point_data() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point_scalar("bad", vec![1.0]);
        let errors = validate_vtu_grid(&grid);
        assert!(
            !errors.is_empty(),
            "Mismatched point data should produce error"
        );
    }
}
/// Write a legacy VTK unstructured grid in binary (big-endian) format.
///
/// VTK legacy binary files use big-endian byte order.
/// Supports arbitrary cell types and float64 point coordinates.
#[allow(dead_code)]
pub fn write_vtk_binary_unstructured(
    path: &str,
    points: &[[f64; 3]],
    cells: &[Vec<usize>],
    cell_types: &[u8],
) -> crate::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "OxiPhysics binary unstructured grid")?;
    writeln!(w, "BINARY")?;
    writeln!(w, "DATASET UNSTRUCTURED_GRID")?;
    writeln!(w, "POINTS {} double", points.len())?;
    for p in points {
        for &coord in p {
            w.write_all(&coord.to_be_bytes())?;
        }
    }
    writeln!(w)?;
    let total_entries: usize = cells.iter().map(|c| 1 + c.len()).sum();
    writeln!(w, "CELLS {} {}", cells.len(), total_entries)?;
    for cell in cells {
        let n = cell.len() as i32;
        w.write_all(&n.to_be_bytes())?;
        for &idx in cell {
            let i = idx as i32;
            w.write_all(&i.to_be_bytes())?;
        }
    }
    writeln!(w)?;
    writeln!(w, "CELL_TYPES {}", cell_types.len())?;
    for &ct in cell_types {
        let ct_i32 = ct as i32;
        w.write_all(&ct_i32.to_be_bytes())?;
    }
    writeln!(w)?;
    w.flush()?;
    Ok(())
}
/// Parse a legacy VTK ASCII file, extracting points and point scalars.
///
/// This is a simplified reader that handles POINTS and POINT_DATA/SCALARS sections.
#[allow(dead_code)]
pub fn read_vtk_legacy_ascii(content: &str) -> std::result::Result<VtkLegacyData, String> {
    let mut lines = content.lines();
    let magic = lines.next().ok_or("Missing magic line")?;
    if !magic.starts_with("# vtk DataFile") {
        return Err(format!("Not a VTK file: {magic}"));
    }
    let title = lines.next().ok_or("Missing title")?.to_string();
    let format = lines.next().ok_or("Missing format line")?;
    if format.trim() != "ASCII" {
        return Err(format!("Only ASCII format supported, got: {format}"));
    }
    let ds_line = lines.next().ok_or("Missing DATASET line")?;
    let dataset_type = ds_line
        .strip_prefix("DATASET ")
        .ok_or("Expected DATASET line")?
        .trim()
        .to_string();
    let mut result = VtkLegacyData {
        title,
        dataset_type,
        points: Vec::new(),
        point_scalars: Vec::new(),
    };
    let remaining: Vec<&str> = lines.collect();
    let mut i = 0;
    while i < remaining.len() {
        let line = remaining[i].trim();
        if line.starts_with("POINTS ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let n: usize = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            i += 1;
            let mut pts = Vec::with_capacity(n);
            let mut coords_read = 0;
            let mut current = [0.0f64; 3];
            while coords_read < n * 3 && i < remaining.len() {
                for tok in remaining[i].split_whitespace() {
                    if let Ok(v) = tok.parse::<f64>() {
                        let idx = coords_read % 3;
                        current[idx] = v;
                        if idx == 2 {
                            pts.push(current);
                            current = [0.0; 3];
                        }
                        coords_read += 1;
                    }
                }
                i += 1;
            }
            result.points = pts;
            continue;
        }
        if line.starts_with("SCALARS ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let name = parts.get(1).unwrap_or(&"scalar").to_string();
            i += 1;
            if i < remaining.len() && remaining[i].trim().starts_with("LOOKUP_TABLE") {
                i += 1;
            }
            let n = result.points.len();
            let mut vals = Vec::with_capacity(n);
            while vals.len() < n && i < remaining.len() {
                for tok in remaining[i].split_whitespace() {
                    if let Ok(v) = tok.parse::<f64>() {
                        vals.push(v);
                    }
                }
                i += 1;
            }
            result.point_scalars.push((name, vals));
            continue;
        }
        i += 1;
    }
    Ok(result)
}
/// Write a ParaView Data (PVD) collection file pointing to a sequence of VTU files.
#[allow(dead_code)]
pub fn write_pvd_file(path: &str, entries: &[PvdEntry]) -> crate::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        w,
        r#"<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">"#
    )?;
    writeln!(w, r#"  <Collection>"#)?;
    for entry in entries {
        writeln!(
            w,
            r#"    <DataSet timestep="{:.6e}" group="{}" part="{}" file="{}"/>"#,
            entry.time, entry.group, entry.part, entry.filename
        )?;
    }
    writeln!(w, r#"  </Collection>"#)?;
    writeln!(w, r#"</VTKFile>"#)?;
    w.flush()?;
    Ok(())
}
/// Generate a sequence of PVD entries for evenly spaced time steps.
#[allow(dead_code)]
pub fn pvd_entries_uniform(
    base_name: &str,
    extension: &str,
    t_start: f64,
    dt: f64,
    n_steps: usize,
) -> Vec<PvdEntry> {
    (0..n_steps)
        .map(|i| {
            let t = t_start + i as f64 * dt;
            let fname = format!("{base_name}_{i:06}.{extension}");
            PvdEntry::new(t, fname)
        })
        .collect()
}
/// Write a triangle surface mesh as a VTK POLYDATA file (ASCII legacy).
#[allow(dead_code)]
pub fn write_vtk_polydata(
    path: &str,
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    normals: Option<&[[f64; 3]]>,
) -> crate::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "OxiPhysics surface mesh")?;
    writeln!(w, "ASCII")?;
    writeln!(w, "DATASET POLYDATA")?;
    writeln!(w, "POINTS {} double", vertices.len())?;
    for v in vertices {
        writeln!(w, "{:.10e} {:.10e} {:.10e}", v[0], v[1], v[2])?;
    }
    let n_tris = triangles.len();
    writeln!(w, "POLYGONS {} {}", n_tris, n_tris * 4)?;
    for tri in triangles {
        writeln!(w, "3 {} {} {}", tri[0], tri[1], tri[2])?;
    }
    if let Some(nrms) = normals {
        writeln!(w, "POINT_DATA {}", vertices.len())?;
        writeln!(w, "NORMALS normals double")?;
        for n in nrms {
            writeln!(w, "{:.10e} {:.10e} {:.10e}", n[0], n[1], n[2])?;
        }
    }
    w.flush()?;
    Ok(())
}
/// Write VTK CELL_DATA with a symmetric 3×3 stress tensor per cell.
///
/// VTK represents tensors as 9-component arrays (row-major).
#[allow(dead_code)]
pub fn write_vtk_cell_stress(
    path: &str,
    points: &[[f64; 3]],
    cells: &[Vec<usize>],
    cell_types: &[u8],
    stress_tensors: &[[[f64; 3]; 3]],
) -> crate::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "OxiPhysics cell stress data")?;
    writeln!(w, "ASCII")?;
    writeln!(w, "DATASET UNSTRUCTURED_GRID")?;
    writeln!(w, "POINTS {} double", points.len())?;
    for p in points {
        writeln!(w, "{:.10e} {:.10e} {:.10e}", p[0], p[1], p[2])?;
    }
    let total_entries: usize = cells.iter().map(|c| 1 + c.len()).sum();
    writeln!(w, "CELLS {} {}", cells.len(), total_entries)?;
    for cell in cells {
        let indices: Vec<String> = cell.iter().map(|i| i.to_string()).collect();
        writeln!(w, "{} {}", cell.len(), indices.join(" "))?;
    }
    writeln!(w, "CELL_TYPES {}", cell_types.len())?;
    for &ct in cell_types {
        writeln!(w, "{ct}")?;
    }
    if !stress_tensors.is_empty() {
        writeln!(w, "CELL_DATA {}", stress_tensors.len())?;
        writeln!(w, "TENSORS stress double")?;
        for sigma in stress_tensors {
            for row in sigma {
                writeln!(w, "{:.10e} {:.10e} {:.10e}", row[0], row[1], row[2])?;
            }
            writeln!(w)?;
        }
    }
    w.flush()?;
    Ok(())
}
/// Write VTK POINT_DATA with velocity vectors for each point.
#[allow(dead_code)]
pub fn write_vtk_point_velocity(
    path: &str,
    points: &[[f64; 3]],
    velocities: &[[f64; 3]],
) -> crate::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "OxiPhysics velocity field")?;
    writeln!(w, "ASCII")?;
    writeln!(w, "DATASET POLYDATA")?;
    writeln!(w, "POINTS {} double", points.len())?;
    for p in points {
        writeln!(w, "{:.10e} {:.10e} {:.10e}", p[0], p[1], p[2])?;
    }
    writeln!(w, "POINT_DATA {}", points.len())?;
    writeln!(w, "VECTORS velocity double")?;
    for v in velocities {
        writeln!(w, "{:.10e} {:.10e} {:.10e}", v[0], v[1], v[2])?;
    }
    w.flush()?;
    Ok(())
}
#[cfg(test)]
mod tests_vtk_extended {
    use super::*;

    use crate::vtk::types::*;

    #[test]
    fn test_vtk_legacy_data_empty() {
        let d = VtkLegacyData::empty();
        assert!(d.points.is_empty());
        assert!(d.point_scalars.is_empty());
    }
    #[test]
    fn test_read_vtk_legacy_ascii_points() {
        let content = "# vtk DataFile Version 3.0\nTest\nASCII\nDATASET POLYDATA\n\
            POINTS 3 float\n0.0 0.0 0.0\n1.0 0.0 0.0\n0.0 1.0 0.0\n";
        let data = read_vtk_legacy_ascii(content).unwrap();
        assert_eq!(data.points.len(), 3);
        assert_eq!(data.title, "Test");
        assert_eq!(data.dataset_type, "POLYDATA");
    }
    #[test]
    fn test_read_vtk_legacy_ascii_point_values() {
        let content = "# vtk DataFile Version 3.0\nTest\nASCII\nDATASET POLYDATA\n\
            POINTS 2 float\n0.0 0.0 0.0\n1.0 0.0 0.0\n";
        let data = read_vtk_legacy_ascii(content).unwrap();
        assert!((data.points[1][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_read_vtk_legacy_ascii_scalars() {
        let content = "# vtk DataFile Version 3.0\nTest\nASCII\nDATASET POLYDATA\n\
            POINTS 2 float\n0.0 0.0 0.0\n1.0 0.0 0.0\n\
            POINT_DATA 2\nSCALARS pressure float 1\nLOOKUP_TABLE default\n0.5 1.5\n";
        let data = read_vtk_legacy_ascii(content).unwrap();
        assert_eq!(data.point_scalars.len(), 1);
        assert_eq!(data.point_scalars[0].0, "pressure");
        assert!((data.point_scalars[0].1[0] - 0.5).abs() < 1e-10);
        assert!((data.point_scalars[0].1[1] - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_read_vtk_legacy_ascii_bad_magic() {
        let content = "NOTAVTK\nTitle\nASCII\nDATASET POLYDATA\n";
        assert!(read_vtk_legacy_ascii(content).is_err());
    }
    #[test]
    fn test_read_vtk_legacy_ascii_binary_unsupported() {
        let content = "# vtk DataFile Version 3.0\nTitle\nBINARY\nDATASET POLYDATA\n";
        assert!(read_vtk_legacy_ascii(content).is_err());
    }
    #[test]
    fn test_pvd_entry_new() {
        let e = PvdEntry::new(0.5, "step_000000.vtu");
        assert!((e.time - 0.5).abs() < 1e-12);
        assert_eq!(e.filename, "step_000000.vtu");
        assert_eq!(e.part, 0);
    }
    #[test]
    fn test_pvd_entries_uniform_count() {
        let entries = pvd_entries_uniform("sim", "vtu", 0.0, 0.01, 5);
        assert_eq!(entries.len(), 5);
    }
    #[test]
    fn test_pvd_entries_uniform_times() {
        let entries = pvd_entries_uniform("sim", "vtu", 1.0, 0.1, 3);
        assert!((entries[0].time - 1.0).abs() < 1e-12);
        assert!((entries[1].time - 1.1).abs() < 1e-12);
        assert!((entries[2].time - 1.2).abs() < 1e-12);
    }
    #[test]
    fn test_pvd_entries_uniform_filenames() {
        let entries = pvd_entries_uniform("flow", "vtu", 0.0, 1.0, 2);
        assert_eq!(entries[0].filename, "flow_000000.vtu");
        assert_eq!(entries[1].filename, "flow_000001.vtu");
    }
    #[test]
    fn test_write_pvd_file() {
        let entries = pvd_entries_uniform("sim", "vtu", 0.0, 0.01, 3);
        let path = "/tmp/test_oxiphysics_pvd.pvd";
        write_pvd_file(path, &entries).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("Collection"));
        assert!(content.contains("timestep"));
        assert!(content.contains("sim_000000.vtu"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_polydata_no_normals() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let path = "/tmp/test_oxiphysics_polydata.vtk";
        write_vtk_polydata(path, &verts, &tris, None).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("POLYDATA"));
        assert!(content.contains("POINTS 3"));
        assert!(content.contains("POLYGONS 1"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_polydata_with_normals() {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 1, 2]];
        let normals = vec![[0.0, 0.0, 1.0]; 3];
        let path = "/tmp/test_oxiphysics_polydata_normals.vtk";
        write_vtk_polydata(path, &verts, &tris, Some(&normals)).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("NORMALS"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_cell_stress() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![5_u8];
        let stress = vec![[[1.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 0.5]]];
        let path = "/tmp/test_oxiphysics_stress.vtk";
        write_vtk_cell_stress(path, &points, &cells, &cell_types, &stress).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("TENSORS stress"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_point_velocity() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let velocities = vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        let path = "/tmp/test_oxiphysics_velocity.vtk";
        write_vtk_point_velocity(path, &points, &velocities).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("VECTORS velocity"));
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_binary_unstructured() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![5_u8];
        let path = "/tmp/test_oxiphysics_binary.vtk";
        write_vtk_binary_unstructured(path, &points, &cells, &cell_types).unwrap();
        let meta = std::fs::metadata(path).unwrap();
        assert!(meta.len() > 0);
        std::fs::remove_file(path).ok();
    }
    #[test]
    fn test_write_vtk_binary_header_ascii() {
        let points = vec![[0.5, 0.5, 0.5]];
        let cells = vec![vec![0_usize]];
        let cell_types = vec![1_u8];
        let path = "/tmp/test_oxiphysics_binary2.vtk";
        write_vtk_binary_unstructured(path, &points, &cells, &cell_types).unwrap();
        let bytes = std::fs::read(path).unwrap();
        let header_str = std::str::from_utf8(&bytes[..5]).unwrap_or("");
        assert_eq!(header_str, "# vtk");
        std::fs::remove_file(path).ok();
    }
}
/// Returns the expected number of nodes for a given [`VtkCellType`].
///
/// Returns `0` for unknown types.
#[allow(dead_code)]
pub fn cell_type_node_count(ct: VtkCellType) -> usize {
    match ct {
        VtkCellType::Vertex => 1,
        VtkCellType::Line => 2,
        VtkCellType::Triangle => 3,
        VtkCellType::Quad => 4,
        VtkCellType::Tetra => 4,
        VtkCellType::Hexahedron => 8,
        VtkCellType::Wedge => 6,
        VtkCellType::Pyramid => 5,
    }
}
/// Return a human-readable name for a [`VtkCellType`].
#[allow(dead_code)]
pub fn cell_type_name(ct: VtkCellType) -> &'static str {
    match ct {
        VtkCellType::Vertex => "Vertex",
        VtkCellType::Line => "Line",
        VtkCellType::Triangle => "Triangle",
        VtkCellType::Quad => "Quad",
        VtkCellType::Tetra => "Tetra",
        VtkCellType::Hexahedron => "Hexahedron",
        VtkCellType::Wedge => "Wedge",
        VtkCellType::Pyramid => "Pyramid",
    }
}
/// Compute axis-aligned bounding box of a [`VtuGrid`]'s points.
///
/// Returns `([min_x, min_y, min_z], [max_x, max_y, max_z])`, or `None` if
/// the grid has no points.
#[allow(dead_code)]
pub fn grid_bounding_box(grid: &VtuGrid) -> Option<([f64; 3], [f64; 3])> {
    if grid.points.is_empty() {
        return None;
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in &grid.points {
        for d in 0..3 {
            if p[d] < lo[d] {
                lo[d] = p[d];
            }
            if p[d] > hi[d] {
                hi[d] = p[d];
            }
        }
    }
    Some((lo, hi))
}
/// Compute centroid of all points in a [`VtuGrid`].
#[allow(dead_code)]
pub fn grid_centroid(grid: &VtuGrid) -> Option<[f64; 3]> {
    let n = grid.points.len();
    if n == 0 {
        return None;
    }
    let sum = grid.points.iter().fold([0.0_f64; 3], |acc, p| {
        [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
    });
    Some([sum[0] / n as f64, sum[1] / n as f64, sum[2] / n as f64])
}
/// Merge two [`VtuGrid`]s into one (points and cells are concatenated).
///
/// Point indices in `b`'s cells are offset by `a.n_points()`.
#[allow(dead_code)]
pub fn merge_vtu_grids(a: &VtuGrid, b: &VtuGrid) -> VtuGrid {
    let mut out = VtuGrid::new();
    for &p in &a.points {
        out.add_point(p);
    }
    for &p in &b.points {
        out.add_point(p);
    }
    for (conn, &ct) in a.cells.iter().zip(a.cell_types.iter()) {
        out.add_cell(conn.clone(), ct);
    }
    let offset = a.n_points();
    for (conn, &ct) in b.cells.iter().zip(b.cell_types.iter()) {
        let shifted: Vec<usize> = conn.iter().map(|&i| i + offset).collect();
        out.add_cell(shifted, ct);
    }
    out
}
/// Compute per-cell volume for tetrahedral cells in a [`VtuGrid`].
///
/// Non-tetrahedral cells get volume `0.0`.
#[allow(dead_code)]
pub fn compute_cell_volumes(grid: &VtuGrid) -> Vec<f64> {
    grid.cells
        .iter()
        .zip(grid.cell_types.iter())
        .map(|(conn, &ct)| {
            if matches!(ct, VtkCellType::Tetra) && conn.len() == 4 {
                let p = |i: usize| grid.points[conn[i]];
                let [a, b, c, d] = [p(0), p(1), p(2), p(3)];
                let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let e3 = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
                let cross = [
                    e2[1] * e3[2] - e2[2] * e3[1],
                    e2[2] * e3[0] - e2[0] * e3[2],
                    e2[0] * e3[1] - e2[1] * e3[0],
                ];
                (e1[0] * cross[0] + e1[1] * cross[1] + e1[2] * cross[2]).abs() / 6.0
            } else {
                0.0
            }
        })
        .collect()
}
/// Compute the Euclidean distance between two 3-D points.
#[allow(dead_code)]
pub fn point_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
#[cfg(test)]
mod vtk_extended_tests {
    use super::*;
    use crate::VtkDataArray;

    use crate::vtk::VtkFieldData;
    use crate::vtk::VtkMultiBlock;

    use crate::vtk::VtkTimeSeries;
    use crate::vtk::types::*;

    #[test]
    fn test_field_data_add_and_get() {
        let mut fd = VtkFieldData::new();
        fd.add("time", vec![1.5, 2.5]);
        let rec = fd.get("time").expect("should find 'time'");
        assert_eq!(rec.values, vec![1.5, 2.5]);
    }
    #[test]
    fn test_field_data_add_scalar() {
        let mut fd = VtkFieldData::new();
        fd.add_scalar("step", 42.0);
        let rec = fd.get("step").expect("should find 'step'");
        assert_eq!(rec.values.len(), 1);
        assert!((rec.values[0] - 42.0).abs() < 1e-12);
    }
    #[test]
    fn test_field_data_get_missing_returns_none() {
        let fd = VtkFieldData::new();
        assert!(fd.get("missing").is_none());
    }
    #[test]
    fn test_field_data_to_vtk_string_contains_name() {
        let mut fd = VtkFieldData::new();
        fd.add("pressure", vec![101325.0]);
        let s = fd.to_vtk_field_string();
        assert!(
            s.contains("pressure"),
            "field string should contain field name"
        );
        assert!(
            s.contains("FIELD"),
            "field string should contain FIELD keyword"
        );
    }
    #[test]
    fn test_field_data_empty_string_empty() {
        let fd = VtkFieldData::new();
        assert!(fd.to_vtk_field_string().is_empty());
    }
    #[test]
    fn test_field_data_multiple_records() {
        let mut fd = VtkFieldData::new();
        fd.add_scalar("dt", 0.001);
        fd.add_scalar("iter", 100.0);
        let s = fd.to_vtk_field_string();
        assert!(s.contains("dt"));
        assert!(s.contains("iter"));
        assert!(s.contains("FieldData 2"));
    }
    #[test]
    fn test_multi_block_empty() {
        let mb = VtkMultiBlock::new("empty");
        assert_eq!(mb.n_blocks(), 0);
        assert_eq!(mb.total_points(), 0);
    }
    #[test]
    fn test_multi_block_add_block() {
        let mut mb = VtkMultiBlock::new("test");
        let g = VtuGrid::from_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        mb.add_block("domain_a", g);
        assert_eq!(mb.n_blocks(), 1);
        assert_eq!(mb.total_points(), 2);
    }
    #[test]
    fn test_multi_block_total_cells() {
        let mut mb = VtkMultiBlock::new("m");
        let mut g1 = VtuGrid::new();
        g1.add_point([0.0, 0.0, 0.0]);
        g1.add_point([1.0, 0.0, 0.0]);
        g1.add_cell(vec![0, 1], VtkCellType::Line);
        let g2 = VtuGrid::from_points(&[[5.0, 0.0, 0.0]]);
        mb.add_block("a", g1);
        mb.add_block("b", g2);
        assert_eq!(mb.total_cells(), 2);
    }
    #[test]
    fn test_multi_block_vtm_string_contains_blocks() {
        let mut mb = VtkMultiBlock::new("sim");
        mb.add_block("fluid", VtuGrid::new());
        mb.add_block("solid", VtuGrid::new());
        let vtm = mb.to_vtm_string();
        assert!(vtm.contains("fluid"));
        assert!(vtm.contains("solid"));
        assert!(vtm.contains("vtkMultiBlockDataSet"));
    }
    #[test]
    fn test_multi_block_write_to_dir() {
        let mut mb = VtkMultiBlock::new("test_mb");
        mb.add_block("part0", VtuGrid::from_points(&[[0.0, 0.0, 0.0]]));
        let written = mb.write_to_dir("/tmp").unwrap();
        assert!(!written.is_empty());
        for f in &written {
            assert!(std::path::Path::new(f).exists(), "file should exist: {f}");
            std::fs::remove_file(f).ok();
        }
    }
    #[test]
    fn test_cell_type_node_counts() {
        assert_eq!(cell_type_node_count(VtkCellType::Vertex), 1);
        assert_eq!(cell_type_node_count(VtkCellType::Line), 2);
        assert_eq!(cell_type_node_count(VtkCellType::Triangle), 3);
        assert_eq!(cell_type_node_count(VtkCellType::Quad), 4);
        assert_eq!(cell_type_node_count(VtkCellType::Tetra), 4);
        assert_eq!(cell_type_node_count(VtkCellType::Hexahedron), 8);
        assert_eq!(cell_type_node_count(VtkCellType::Wedge), 6);
        assert_eq!(cell_type_node_count(VtkCellType::Pyramid), 5);
    }
    #[test]
    fn test_cell_type_names() {
        assert_eq!(cell_type_name(VtkCellType::Triangle), "Triangle");
        assert_eq!(cell_type_name(VtkCellType::Tetra), "Tetra");
        assert_eq!(cell_type_name(VtkCellType::Hexahedron), "Hexahedron");
    }
    #[test]
    fn test_grid_bounding_box_empty_is_none() {
        let grid = VtuGrid::new();
        assert!(grid_bounding_box(&grid).is_none());
    }
    #[test]
    fn test_grid_bounding_box_values() {
        let mut g = VtuGrid::new();
        g.add_point([1.0, 2.0, 3.0]);
        g.add_point([4.0, 5.0, 6.0]);
        let (lo, hi) = grid_bounding_box(&g).unwrap();
        assert!((lo[0] - 1.0).abs() < 1e-12);
        assert!((hi[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_grid_centroid_single_point() {
        let mut g = VtuGrid::new();
        g.add_point([3.0, 4.0, 5.0]);
        let c = grid_centroid(&g).unwrap();
        assert!((c[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_grid_centroid_two_points() {
        let mut g = VtuGrid::new();
        g.add_point([0.0, 0.0, 0.0]);
        g.add_point([2.0, 4.0, 6.0]);
        let c = grid_centroid(&g).unwrap();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_grid_centroid_empty_is_none() {
        let g = VtuGrid::new();
        assert!(grid_centroid(&g).is_none());
    }
    #[test]
    fn test_merge_vtu_grids_point_count() {
        let a = VtuGrid::from_points(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let b = VtuGrid::from_points(&[[2.0, 0.0, 0.0]]);
        let merged = merge_vtu_grids(&a, &b);
        assert_eq!(merged.n_points(), 3);
    }
    #[test]
    fn test_merge_vtu_grids_cell_count() {
        let a = VtuGrid::from_points(&[[0.0, 0.0, 0.0]]);
        let b = VtuGrid::from_points(&[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let merged = merge_vtu_grids(&a, &b);
        assert_eq!(merged.n_cells(), 3);
    }
    #[test]
    fn test_merge_vtu_grids_indices_offset() {
        let mut a = VtuGrid::new();
        a.add_point([0.0, 0.0, 0.0]);
        a.add_point([1.0, 0.0, 0.0]);
        a.add_cell(vec![0, 1], VtkCellType::Line);
        let mut b = VtuGrid::new();
        b.add_point([2.0, 0.0, 0.0]);
        b.add_point([3.0, 0.0, 0.0]);
        b.add_cell(vec![0, 1], VtkCellType::Line);
        let merged = merge_vtu_grids(&a, &b);
        assert_eq!(merged.cells[1], vec![2, 3]);
    }
    #[test]
    fn test_compute_cell_volumes_tet() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let grid = VtuGrid::from_tet_mesh(&nodes, &[[0, 1, 2, 3]]);
        let vols = compute_cell_volumes(&grid);
        assert_eq!(vols.len(), 1);
        assert!(
            (vols[0] - 1.0 / 6.0).abs() < 1e-10,
            "tet volume should be 1/6, got {}",
            vols[0]
        );
    }
    #[test]
    fn test_compute_cell_volumes_vertex_zero() {
        let grid = VtuGrid::from_points(&[[0.0, 0.0, 0.0]]);
        let vols = compute_cell_volumes(&grid);
        assert_eq!(vols[0], 0.0);
    }
    #[test]
    fn test_compute_cell_volumes_multiple_tets() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ];
        let elements = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let grid = VtuGrid::from_tet_mesh(&nodes, &elements);
        let vols = compute_cell_volumes(&grid);
        assert_eq!(vols.len(), 2);
        for v in &vols {
            assert!(
                (v - 1.0 / 6.0).abs() < 1e-10,
                "each tet should have volume 1/6"
            );
        }
    }
    #[test]
    fn test_point_distance_known() {
        let d = point_distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_point_distance_same_point_is_zero() {
        let d = point_distance([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(d.abs() < 1e-12);
    }
    #[test]
    fn test_vtk_data_array_scalar_len() {
        let arr = VtkDataArray::Scalar {
            name: "p".into(),
            values: vec![1.0, 2.0, 3.0],
        };
        assert_eq!(arr.len(), 3);
        assert_eq!(arr.n_components(), 1);
        assert!(!arr.is_empty());
    }
    #[test]
    fn test_vtk_data_array_vector3_n_components() {
        let arr = VtkDataArray::Vector3 {
            name: "v".into(),
            values: vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        };
        assert_eq!(arr.n_components(), 3);
        assert_eq!(arr.len(), 2);
    }
    #[test]
    fn test_vtk_data_array_integer() {
        let arr = VtkDataArray::Integer {
            name: "id".into(),
            values: vec![1, 2, 3, 4],
        };
        assert_eq!(arr.len(), 4);
        assert_eq!(arr.name(), "id");
    }
    #[test]
    fn test_vtk_data_array_empty() {
        let arr = VtkDataArray::Scalar {
            name: "empty".into(),
            values: vec![],
        };
        assert!(arr.is_empty());
    }
    #[test]
    fn test_vtu_cell_scalar_data_in_xml() {
        let mut grid = VtuGrid::new();
        for i in 0..4 {
            grid.add_point([i as f64, 0.0, 0.0]);
        }
        grid.add_cell(vec![0, 1, 2, 3], VtkCellType::Quad);
        grid.add_cell_scalar("material_id", vec![7.0]);
        let xml = grid.to_vtu_string();
        assert!(xml.contains("<CellData>"), "CellData section missing");
        assert!(xml.contains("material_id"), "material_id field missing");
    }
    #[test]
    fn test_vtu_vector_point_data_in_xml() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point([0.0, 1.0, 0.0]);
        grid.add_cell(vec![0, 1, 2], VtkCellType::Triangle);
        grid.add_point_vector("velocity", vec![[1.0, 0.0, 0.0]; 3]);
        let xml = grid.to_vtu_string();
        assert!(xml.contains("velocity"), "velocity field should be in XML");
        assert!(xml.contains("NumberOfComponents=\"3\""));
    }
    #[test]
    fn test_time_series_estimated_size() {
        let mut ts = VtkTimeSeries::new("s");
        let mut g = VtuGrid::new();
        for i in 0..10 {
            g.add_point([i as f64, 0.0, 0.0]);
        }
        ts.push(0.0, g);
        let sz = ts.estimated_size_bytes();
        assert!(sz > 0);
    }
    #[test]
    fn test_time_series_increasing_times() {
        let mut ts = VtkTimeSeries::new("flow");
        ts.push(0.0, VtuGrid::new());
        ts.push(1.0, VtuGrid::new());
        ts.push(2.5, VtuGrid::new());
        assert_eq!(ts.n_steps(), 3);
        assert!((ts.times[2] - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_validate_vtu_grid_cell_data_mismatch() {
        let mut grid = VtuGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_point([0.0, 1.0, 0.0]);
        grid.add_cell(vec![0, 1, 2], VtkCellType::Triangle);
        grid.cell_data.push(VtkDataArray::Scalar {
            name: "bad_cell".into(),
            values: vec![1.0, 2.0],
        });
        let errors = validate_vtu_grid(&grid);
        assert!(
            !errors.is_empty(),
            "mismatched cell data should produce error"
        );
    }
    #[test]
    fn test_validate_vtu_grid_empty_valid() {
        let grid = VtuGrid::new();
        let errors = validate_vtu_grid(&grid);
        assert!(errors.is_empty(), "empty grid should be valid");
    }
}
