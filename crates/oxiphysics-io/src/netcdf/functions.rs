//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::type_complexity)]
use crate::netcdf::types::*;

/// Create an AMBER-convention NetCDF trajectory file.
///
/// AMBER trajectories use:
/// - dimensions: frame (unlimited), atom, spatial (3)
/// - variables: coordinates(frame, atom, spatial), cell_lengths(frame, spatial),
///   cell_angles(frame, spatial), time(frame)
#[allow(dead_code)]
pub fn create_amber_trajectory(
    n_atoms: usize,
    n_frames: usize,
    coordinates: Vec<f64>,
    time: Vec<f64>,
    cell_lengths: Option<Vec<f64>>,
    cell_angles: Option<Vec<f64>>,
) -> NetCdfFile {
    let mut file = NetCdfFile::new();
    file.add_unlimited_dimension("frame", n_frames);
    file.add_dimension("atom", n_atoms);
    file.add_dimension("spatial", 3);
    file.add_attribute("Conventions", "AMBER");
    file.add_attribute("ConventionVersion", "1.0");
    file.add_attribute("program", "OxiPhysics");
    let mut coords_var = NetCdfVariable::new(
        "coordinates",
        vec!["frame".into(), "atom".into(), "spatial".into()],
        coordinates,
    );
    coords_var.add_attribute("units", "angstrom");
    file.add_variable(coords_var);
    let mut time_var = NetCdfVariable::new("time", vec!["frame".into()], time);
    time_var.add_attribute("units", "picosecond");
    file.add_variable(time_var);
    if let Some(lengths) = cell_lengths {
        let mut var = NetCdfVariable::new(
            "cell_lengths",
            vec!["frame".into(), "spatial".into()],
            lengths,
        );
        var.add_attribute("units", "angstrom");
        file.add_variable(var);
    }
    if let Some(angles) = cell_angles {
        let mut var = NetCdfVariable::new(
            "cell_angles",
            vec!["frame".into(), "spatial".into()],
            angles,
        );
        var.add_attribute("units", "degree");
        file.add_variable(var);
    }
    file
}
/// Extract coordinates for a specific frame from an AMBER trajectory.
#[allow(dead_code)]
pub fn extract_frame_coordinates(file: &NetCdfFile, frame_idx: usize) -> Option<Vec<[f64; 3]>> {
    let coords_var = file.get_variable("coordinates")?;
    let n_atoms = file.get_dimension_size("atom")?;
    let start = frame_idx * n_atoms * 3;
    let end = start + n_atoms * 3;
    if end > coords_var.data.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(n_atoms);
    for i in 0..n_atoms {
        let base = start + i * 3;
        positions.push([
            coords_var.data[base],
            coords_var.data[base + 1],
            coords_var.data[base + 2],
        ]);
    }
    Some(positions)
}
/// Create a coordinate variable (a 1D variable with the same name as its dimension).
#[allow(dead_code)]
pub fn create_coordinate_variable(name: &str, values: Vec<f64>, units: &str) -> NetCdfVariable {
    let mut var = NetCdfVariable::new(name, vec![name.to_string()], values);
    var.add_attribute("units", units);
    var
}
/// Write a [`NetCdfFile`] as human-readable CDL text.
#[allow(dead_code)]
pub fn write_ncf_text(file: &NetCdfFile) -> String {
    let mut out = String::new();
    out.push_str("netcdf data {\n");
    out.push_str("dimensions:\n");
    for (name, size) in &file.dimensions {
        if file.unlimited_dim.as_deref() == Some(name.as_str()) {
            out.push_str(&format!("\t{name} = UNLIMITED ; // ({size} currently)\n"));
        } else {
            out.push_str(&format!("\t{name} = {size} ;\n"));
        }
    }
    out.push_str("variables:\n");
    for var in &file.variables {
        let dims = var.dims.join(", ");
        out.push_str(&format!("\tdouble {}({}) ;\n", var.name, dims));
        for attr in &var.attributes {
            out.push_str(&format!(
                "\t\t{}:{} = \"{}\" ;\n",
                var.name, attr.key, attr.value
            ));
        }
    }
    if !file.attributes.is_empty() {
        out.push_str("// global attributes:\n");
        for (key, value) in &file.attributes {
            out.push_str(&format!("\t\t:data:{key} = \"{value}\" ;\n"));
        }
    }
    out.push_str("data:\n");
    for var in &file.variables {
        let values: Vec<String> = var.data.iter().map(|v| format!("{v}")).collect();
        out.push_str(&format!("\t{} = {} ;\n", var.name, values.join(", ")));
    }
    out.push_str("}\n");
    out
}
/// Parse a CDL text string produced by [`write_ncf_text`] into a [`NetCdfFile`].
///
/// Returns `Err(String)` on the first parse failure.
#[allow(dead_code)]
pub fn parse_ncf_text(cdl: &str) -> Result<NetCdfFile, String> {
    let mut dimensions: Vec<(String, usize)> = Vec::new();
    let mut variable_headers: Vec<(String, Vec<String>)> = Vec::new();
    let mut variable_attrs: std::collections::HashMap<String, Vec<VariableAttribute>> =
        std::collections::HashMap::new();
    let mut attributes: Vec<(String, String)> = Vec::new();
    let mut data_map: std::collections::HashMap<String, Vec<f64>> =
        std::collections::HashMap::new();
    let mut unlimited_dim: Option<String> = None;
    #[derive(PartialEq)]
    pub(super) enum Section {
        None,
        Dimensions,
        Variables,
        Attributes,
        Data,
    }
    let mut section = Section::None;
    for raw in cdl.lines() {
        let line = raw.trim();
        if line.is_empty() || line == "{" || line == "}" || line.starts_with("netcdf ") {
            continue;
        }
        if line == "dimensions:" {
            section = Section::Dimensions;
            continue;
        }
        if line == "variables:" {
            section = Section::Variables;
            continue;
        }
        if line.starts_with("// global attributes") {
            section = Section::Attributes;
            continue;
        }
        if line == "data:" {
            section = Section::Data;
            continue;
        }
        let line_no_comment = if let Some(pos) = line.find("//") {
            line[..pos].trim()
        } else {
            line
        };
        let line_clean = line_no_comment.trim_end_matches(';').trim();
        match section {
            Section::Dimensions => {
                let parts: Vec<&str> = line_clean.splitn(2, '=').collect();
                if parts.len() != 2 {
                    return Err(format!("Bad dimension line: {raw}"));
                }
                let dim_name = parts[0].trim().to_string();
                let size_str = parts[1].trim();
                if size_str == "UNLIMITED" {
                    dimensions.push((dim_name.clone(), 0));
                    unlimited_dim = Some(dim_name);
                } else {
                    let dim_size = size_str
                        .parse::<usize>()
                        .map_err(|_| format!("Bad dimension size in: {raw}"))?;
                    dimensions.push((dim_name, dim_size));
                }
            }
            Section::Variables => {
                if line_clean.contains(':') && !line_clean.starts_with("double ") {
                    let parts: Vec<&str> = line_clean.splitn(2, ':').collect();
                    if parts.len() == 2 {
                        let var_name = parts[0].trim().to_string();
                        let rest: Vec<&str> = parts[1].splitn(2, '=').collect();
                        if rest.len() == 2 {
                            let key = rest[0].trim().to_string();
                            let val = rest[1].trim().trim_matches('"').to_string();
                            variable_attrs
                                .entry(var_name)
                                .or_default()
                                .push(VariableAttribute { key, value: val });
                        }
                    }
                    continue;
                }
                if !line_clean.starts_with("double ") {
                    continue;
                }
                let rest = line_clean["double ".len()..].trim();
                if let Some(paren) = rest.find('(') {
                    let var_name = rest[..paren].trim().to_string();
                    let dims_str = rest[paren + 1..].trim_end_matches(')').trim();
                    let dims: Vec<String> = if dims_str.is_empty() {
                        Vec::new()
                    } else {
                        dims_str.split(',').map(|s| s.trim().to_string()).collect()
                    };
                    variable_headers.push((var_name, dims));
                } else {
                    return Err(format!("Bad variable line: {raw}"));
                }
            }
            Section::Attributes => {
                if !line_clean.starts_with(":data:") {
                    continue;
                }
                let rest = &line_clean[":data:".len()..];
                let parts: Vec<&str> = rest.splitn(2, '=').collect();
                if parts.len() != 2 {
                    return Err(format!("Bad attribute line: {raw}"));
                }
                let key = parts[0].trim().to_string();
                let val = parts[1].trim().trim_matches('"').to_string();
                attributes.push((key, val));
            }
            Section::Data => {
                let parts: Vec<&str> = line_clean.splitn(2, '=').collect();
                if parts.len() != 2 {
                    return Err(format!("Bad data line: {raw}"));
                }
                let var_name = parts[0].trim().to_string();
                let values: Result<Vec<f64>, _> = parts[1]
                    .split(',')
                    .map(|s| s.trim().parse::<f64>())
                    .collect();
                let values =
                    values.map_err(|_| format!("Bad numeric data in variable '{var_name}'"))?;
                data_map.insert(var_name, values);
            }
            Section::None => {}
        }
    }
    let variables: Vec<NetCdfVariable> = variable_headers
        .into_iter()
        .map(|(name, dims)| {
            let data = data_map.remove(&name).unwrap_or_default();
            let attrs = variable_attrs.remove(&name).unwrap_or_default();
            NetCdfVariable {
                name,
                dims,
                data,
                attributes: attrs,
            }
        })
        .collect();
    Ok(NetCdfFile {
        dimensions,
        variables,
        attributes,
        unlimited_dim,
    })
}
/// Reshape flat data into a 2D array (row-major) given row and column sizes.
#[allow(dead_code)]
pub fn reshape_2d(data: &[f64], n_rows: usize, n_cols: usize) -> Vec<Vec<f64>> {
    let mut result = Vec::with_capacity(n_rows);
    for r in 0..n_rows {
        let start = r * n_cols;
        let end = (start + n_cols).min(data.len());
        if start < data.len() {
            result.push(data[start..end].to_vec());
        }
    }
    result
}
/// Compute statistics for a variable's data.
#[allow(dead_code)]
pub fn variable_stats(var: &NetCdfVariable) -> Option<(f64, f64, f64)> {
    if var.data.is_empty() {
        return None;
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for &v in &var.data {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    let mean = sum / var.data.len() as f64;
    Some((min, max, mean))
}
/// Create a NetCDF string variable for atom names, residue names, etc.
#[allow(dead_code)]
pub fn create_string_variable(name: &str, dims: Vec<String>, values: Vec<String>) -> Nc4Variable {
    let mut var = Nc4Variable::string_var(name, dims, values);
    var.add_attribute("dtype", "string");
    var
}
/// Copy global attributes from a source file into each variable of a target file,
/// unless the variable already has an attribute with the same key.
#[allow(dead_code)]
pub fn inherit_global_attributes(source: &NetCdfFile, target: &mut NetCdfFile) {
    for (key, value) in &source.attributes {
        for var in target.variables.iter_mut() {
            if var.get_attribute(key).is_none() {
                var.add_attribute(key, value);
            }
        }
    }
}
/// Merge two `NetCdfFile` datasets by combining variables (source variables
/// overwrite target variables with the same name).
#[allow(dead_code)]
pub fn merge_netcdf_files(base: &NetCdfFile, overlay: &NetCdfFile) -> NetCdfFile {
    let mut result = base.clone();
    for var in &overlay.variables {
        if let Some(existing) = result.variables.iter_mut().find(|v| v.name == var.name) {
            *existing = var.clone();
        } else {
            result.variables.push(var.clone());
        }
    }
    for (key, value) in &overlay.attributes {
        if let Some(existing) = result.attributes.iter_mut().find(|(k, _)| k == key) {
            existing.1 = value.clone();
        } else {
            result.attributes.push((key.clone(), value.clone()));
        }
    }
    result
}
/// Parse dimension declarations from a CDL string.
///
/// Returns a list of `(name, size)` pairs for lines matching `name` = `size` ;`.
#[allow(dead_code)]
pub fn parse_cdl_dimensions(cdl: &str) -> Vec<(String, usize)> {
    let mut result = Vec::new();
    let mut in_dimensions = false;
    for line in cdl.lines() {
        let trimmed = line.trim();
        if trimmed == "dimensions:" {
            in_dimensions = true;
            continue;
        }
        if in_dimensions {
            if trimmed == "variables:"
                || trimmed == "data:"
                || trimmed.starts_with("//")
                || trimmed == "}"
            {
                break;
            }
            let clean = trimmed.trim_end_matches(';').trim();
            let parts: Vec<&str> = clean.splitn(2, '=').collect();
            if parts.len() == 2 {
                let name = parts[0].trim().to_string();
                let size_str = parts[1].trim();
                if let Ok(size) = size_str.parse::<usize>() {
                    result.push((name, size));
                }
            }
        }
    }
    result
}
/// Compute trajectory statistics.
///
/// Returns `(total_time, avg_displacement, n_frames)`.
///
/// `total_time` = last time – first time (0 if fewer than 2 frames).
/// `avg_displacement` = mean Euclidean distance from frame 0 to the last frame
/// (averaged over atoms; 0 if no atoms or only 1 frame).
#[allow(dead_code)]
pub fn trajectory_stats(times: &[f64], positions: &[Vec<[f64; 3]>]) -> (f64, f64, usize) {
    let n_frames = times.len();
    if n_frames == 0 {
        return (0.0, 0.0, 0);
    }
    let total_time = if n_frames >= 2 {
        times[n_frames - 1] - times[0]
    } else {
        0.0
    };
    let avg_displacement = if n_frames >= 2 && !positions.is_empty() {
        let first = &positions[0];
        let last = &positions[positions.len() - 1];
        let n = first.len().min(last.len());
        if n == 0 {
            0.0
        } else {
            let sum: f64 = (0..n)
                .map(|i| {
                    let dx = last[i][0] - first[i][0];
                    let dy = last[i][1] - first[i][1];
                    let dz = last[i][2] - first[i][2];
                    (dx * dx + dy * dy + dz * dz).sqrt()
                })
                .sum();
            sum / n as f64
        }
    } else {
        0.0
    };
    (total_time, avg_displacement, n_frames)
}
#[cfg(test)]
mod tests {
    use super::*;

    fn make_file() -> NetCdfFile {
        NetCdfFile {
            dimensions: vec![("time".to_string(), 3), ("space".to_string(), 2)],
            variables: vec![
                NetCdfVariable {
                    name: "temperature".to_string(),
                    dims: vec!["time".to_string(), "space".to_string()],
                    data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                    attributes: Vec::new(),
                },
                NetCdfVariable {
                    name: "pressure".to_string(),
                    dims: vec!["time".to_string()],
                    data: vec![101.3, 100.8, 99.5],
                    attributes: Vec::new(),
                },
            ],
            attributes: vec![
                ("title".to_string(), "Test dataset".to_string()),
                ("units".to_string(), "SI".to_string()),
            ],
            unlimited_dim: None,
        }
    }
    #[test]
    fn test_write_ncf_text_structure() {
        let file = make_file();
        let cdl = write_ncf_text(&file);
        assert!(cdl.contains("dimensions:"));
        assert!(cdl.contains("time = 3"));
        assert!(cdl.contains("space = 2"));
        assert!(cdl.contains("variables:"));
        assert!(cdl.contains("double temperature(time, space)"));
        assert!(cdl.contains("double pressure(time)"));
        assert!(cdl.contains("data:"));
        assert!(cdl.contains("temperature ="));
        assert!(cdl.contains("pressure ="));
        assert!(cdl.contains(r#":data:title = "Test dataset""#));
    }
    #[test]
    fn test_ncf_roundtrip() {
        let original = make_file();
        let cdl = write_ncf_text(&original);
        let parsed = parse_ncf_text(&cdl).expect("parse failed");
        assert_eq!(parsed.dimensions.len(), original.dimensions.len());
        for (orig, got) in original.dimensions.iter().zip(parsed.dimensions.iter()) {
            assert_eq!(orig.0, got.0, "dimension name mismatch");
            assert_eq!(orig.1, got.1, "dimension size mismatch");
        }
        assert_eq!(parsed.variables.len(), original.variables.len());
        for (orig, got) in original.variables.iter().zip(parsed.variables.iter()) {
            assert_eq!(orig.name, got.name, "variable name mismatch");
            assert_eq!(orig.dims, got.dims, "variable dims mismatch");
            assert_eq!(orig.data.len(), got.data.len(), "data length mismatch");
            for (a, b) in orig.data.iter().zip(got.data.iter()) {
                assert!((a - b).abs() < 1e-9, "data value mismatch: {a} vs {b}");
            }
        }
        assert_eq!(parsed.attributes.len(), original.attributes.len());
        for (orig, got) in original.attributes.iter().zip(parsed.attributes.iter()) {
            assert_eq!(orig.0, got.0, "attribute key mismatch");
            assert_eq!(orig.1, got.1, "attribute value mismatch");
        }
    }
    #[test]
    fn test_ncf_empty_dataset() {
        let file = NetCdfFile {
            dimensions: vec![],
            variables: vec![],
            attributes: vec![],
            unlimited_dim: None,
        };
        let cdl = write_ncf_text(&file);
        let parsed = parse_ncf_text(&cdl).expect("parse failed");
        assert!(parsed.dimensions.is_empty());
        assert!(parsed.variables.is_empty());
        assert!(parsed.attributes.is_empty());
    }
    #[test]
    fn test_netcdf_file_new() {
        let file = NetCdfFile::new();
        assert!(file.dimensions.is_empty());
        assert!(file.variables.is_empty());
        assert!(file.attributes.is_empty());
        assert!(file.unlimited_dim.is_none());
    }
    #[test]
    fn test_add_dimension() {
        let mut file = NetCdfFile::new();
        file.add_dimension("time", 10);
        assert_eq!(file.get_dimension_size("time"), Some(10));
        assert_eq!(file.get_dimension_size("missing"), None);
    }
    #[test]
    fn test_add_unlimited_dimension() {
        let mut file = NetCdfFile::new();
        file.add_unlimited_dimension("frame", 5);
        assert!(file.is_unlimited_dimension("frame"));
        assert!(!file.is_unlimited_dimension("atom"));
        assert_eq!(file.get_dimension_size("frame"), Some(5));
    }
    #[test]
    fn test_variable_new_and_attributes() {
        let mut var = NetCdfVariable::new("temp", vec!["time".into()], vec![1.0, 2.0, 3.0]);
        assert_eq!(var.len(), 3);
        assert!(!var.is_empty());
        var.add_attribute("units", "K");
        assert_eq!(var.get_attribute("units"), Some("K"));
        assert_eq!(var.get_attribute("missing"), None);
    }
    #[test]
    fn test_variable_empty() {
        let var = NetCdfVariable::new("empty", vec![], vec![]);
        assert!(var.is_empty());
        assert_eq!(var.len(), 0);
    }
    #[test]
    fn test_get_variable() {
        let mut file = NetCdfFile::new();
        file.add_variable(NetCdfVariable::new("x", vec![], vec![1.0]));
        file.add_variable(NetCdfVariable::new("y", vec![], vec![2.0]));
        assert!(file.get_variable("x").is_some());
        assert!(file.get_variable("z").is_none());
        assert_eq!(file.get_variable("y").unwrap().data[0], 2.0);
    }
    #[test]
    fn test_get_variable_mut() {
        let mut file = NetCdfFile::new();
        file.add_variable(NetCdfVariable::new("x", vec![], vec![1.0]));
        file.get_variable_mut("x").unwrap().data.push(2.0);
        assert_eq!(file.get_variable("x").unwrap().data.len(), 2);
    }
    #[test]
    fn test_global_attribute() {
        let mut file = NetCdfFile::new();
        file.add_attribute("title", "test");
        assert_eq!(file.get_attribute("title"), Some("test"));
        assert_eq!(file.get_attribute("missing"), None);
    }
    #[test]
    fn test_variable_names() {
        let mut file = NetCdfFile::new();
        file.add_variable(NetCdfVariable::new("a", vec![], vec![]));
        file.add_variable(NetCdfVariable::new("b", vec![], vec![]));
        let names = file.variable_names();
        assert_eq!(names, vec!["a", "b"]);
    }
    #[test]
    fn test_dimension_names() {
        let mut file = NetCdfFile::new();
        file.add_dimension("time", 10);
        file.add_dimension("space", 3);
        let names = file.dimension_names();
        assert_eq!(names, vec!["time", "space"]);
    }
    #[test]
    fn test_dimension_descriptor() {
        let unlimited = NetCdfDimension {
            name: "frame".into(),
            size: None,
        };
        assert!(unlimited.is_unlimited());
        assert_eq!(unlimited.effective_size(), 0);
        let fixed = NetCdfDimension {
            name: "atom".into(),
            size: Some(100),
        };
        assert!(!fixed.is_unlimited());
        assert_eq!(fixed.effective_size(), 100);
    }
    #[test]
    fn test_amber_trajectory_creation() {
        let n_atoms = 3;
        let n_frames = 2;
        let coords: Vec<f64> = (0..18).map(|i| i as f64).collect();
        let time = vec![0.0, 1.0];
        let cell_lengths = Some(vec![10.0, 10.0, 10.0, 10.0, 10.0, 10.0]);
        let cell_angles = Some(vec![90.0, 90.0, 90.0, 90.0, 90.0, 90.0]);
        let file =
            create_amber_trajectory(n_atoms, n_frames, coords, time, cell_lengths, cell_angles);
        assert_eq!(file.get_attribute("Conventions"), Some("AMBER"));
        assert_eq!(file.get_dimension_size("atom"), Some(3));
        assert_eq!(file.get_dimension_size("spatial"), Some(3));
        assert!(file.is_unlimited_dimension("frame"));
        assert!(file.get_variable("coordinates").is_some());
        assert!(file.get_variable("time").is_some());
        assert!(file.get_variable("cell_lengths").is_some());
        assert!(file.get_variable("cell_angles").is_some());
    }
    #[test]
    fn test_amber_trajectory_no_cell() {
        let file = create_amber_trajectory(2, 1, vec![0.0; 6], vec![0.0], None, None);
        assert!(file.get_variable("cell_lengths").is_none());
        assert!(file.get_variable("cell_angles").is_none());
    }
    #[test]
    fn test_extract_frame_coordinates() {
        let coords = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let file = create_amber_trajectory(2, 2, coords, vec![0.0, 1.0], None, None);
        let frame0 = extract_frame_coordinates(&file, 0).unwrap();
        assert_eq!(frame0.len(), 2);
        assert!((frame0[0][0] - 1.0).abs() < 1e-10);
        assert!((frame0[1][2] - 6.0).abs() < 1e-10);
        let frame1 = extract_frame_coordinates(&file, 1).unwrap();
        assert!((frame1[0][0] - 7.0).abs() < 1e-10);
        assert!(extract_frame_coordinates(&file, 5).is_none());
    }
    #[test]
    fn test_coordinate_variable() {
        let var = create_coordinate_variable("time", vec![0.0, 0.5, 1.0], "seconds");
        assert_eq!(var.name, "time");
        assert_eq!(var.dims, vec!["time"]);
        assert_eq!(var.get_attribute("units"), Some("seconds"));
        assert_eq!(var.data.len(), 3);
    }
    #[test]
    fn test_reshape_2d() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let matrix = reshape_2d(&data, 2, 3);
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix[0], vec![1.0, 2.0, 3.0]);
        assert_eq!(matrix[1], vec![4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_reshape_2d_empty() {
        let matrix = reshape_2d(&[], 0, 3);
        assert!(matrix.is_empty());
    }
    #[test]
    fn test_variable_stats() {
        let var = NetCdfVariable::new("t", vec![], vec![1.0, 5.0, 3.0]);
        let (min, max, mean) = variable_stats(&var).unwrap();
        assert!((min - 1.0).abs() < 1e-10);
        assert!((max - 5.0).abs() < 1e-10);
        assert!((mean - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_variable_stats_empty() {
        let var = NetCdfVariable::new("empty", vec![], vec![]);
        assert!(variable_stats(&var).is_none());
    }
    #[test]
    fn test_unlimited_dim_roundtrip() {
        let mut file = NetCdfFile::new();
        file.add_unlimited_dimension("frame", 5);
        file.add_dimension("atom", 3);
        file.add_variable(NetCdfVariable::new(
            "x",
            vec!["frame".into(), "atom".into()],
            vec![1.0; 15],
        ));
        let cdl = write_ncf_text(&file);
        assert!(cdl.contains("UNLIMITED"));
        let parsed = parse_ncf_text(&cdl).unwrap();
        assert!(parsed.unlimited_dim.is_some());
        assert_eq!(parsed.unlimited_dim.as_deref(), Some("frame"));
    }
    #[test]
    fn test_variable_attributes_roundtrip() {
        let mut file = NetCdfFile::new();
        file.add_dimension("x", 3);
        let mut var = NetCdfVariable::new("temp", vec!["x".into()], vec![1.0, 2.0, 3.0]);
        var.add_attribute("units", "kelvin");
        var.add_attribute("long_name", "temperature");
        file.add_variable(var);
        let cdl = write_ncf_text(&file);
        let parsed = parse_ncf_text(&cdl).unwrap();
        let parsed_var = parsed.get_variable("temp").unwrap();
        assert_eq!(parsed_var.get_attribute("units"), Some("kelvin"));
        assert_eq!(parsed_var.get_attribute("long_name"), Some("temperature"));
    }
    #[test]
    fn test_amber_trajectory_roundtrip_cdl() {
        let file = create_amber_trajectory(
            2,
            1,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![0.0],
            None,
            None,
        );
        let cdl = write_ncf_text(&file);
        let parsed = parse_ncf_text(&cdl).unwrap();
        assert_eq!(parsed.get_attribute("Conventions"), Some("AMBER"));
        assert!(parsed.get_variable("coordinates").is_some());
    }
    #[test]
    fn test_full_amber_with_cell() {
        let file = create_amber_trajectory(
            1,
            2,
            vec![0.0; 6],
            vec![0.0, 0.5],
            Some(vec![10.0, 10.0, 10.0, 10.0, 10.0, 10.0]),
            Some(vec![90.0, 90.0, 90.0, 90.0, 90.0, 90.0]),
        );
        let coords_var = file.get_variable("coordinates").unwrap();
        assert_eq!(coords_var.get_attribute("units"), Some("angstrom"));
        let time_var = file.get_variable("time").unwrap();
        assert_eq!(time_var.get_attribute("units"), Some("picosecond"));
    }
    #[test]
    fn test_nc4_data_type_byte_width() {
        assert_eq!(Nc4DataType::Float64.byte_width(), 8);
        assert_eq!(Nc4DataType::Float32.byte_width(), 4);
        assert_eq!(Nc4DataType::Int32.byte_width(), 4);
        assert_eq!(Nc4DataType::UInt8.byte_width(), 1);
        assert_eq!(Nc4DataType::String.byte_width(), 0);
    }
    #[test]
    fn test_nc4_data_type_name() {
        assert_eq!(Nc4DataType::Float64.type_name(), "double");
        assert_eq!(Nc4DataType::String.type_name(), "string");
        assert_eq!(Nc4DataType::Int32.type_name(), "int");
    }
    #[test]
    fn test_nc4_variable_float64() {
        let var = Nc4Variable::float64("temp", vec!["time".into()], vec![1.0, 2.0, 3.0]);
        assert_eq!(var.data_type, Nc4DataType::Float64);
        assert_eq!(var.len(), 3);
        assert!(!var.is_empty());
    }
    #[test]
    fn test_nc4_variable_string() {
        let var = Nc4Variable::string_var(
            "names",
            vec!["n".into()],
            vec!["ALA".to_string(), "GLY".to_string()],
        );
        assert_eq!(var.data_type, Nc4DataType::String);
        assert_eq!(var.string_data.len(), 2);
        assert_eq!(var.string_data[0], "ALA");
        assert_eq!(var.len(), 2);
    }
    #[test]
    fn test_nc4_variable_int32() {
        let var = Nc4Variable::int32("step", vec!["frame".into()], vec![0, 100, 200]);
        assert_eq!(var.data_type, Nc4DataType::Int32);
        assert_eq!(var.int_data.len(), 3);
        assert_eq!(var.int_data[1], 100);
    }
    #[test]
    fn test_nc4_variable_attributes() {
        let mut var = Nc4Variable::float64("temp", vec![], vec![]);
        var.add_attribute("units", "K");
        var.add_attribute("long_name", "temperature");
        assert_eq!(var.get_attribute("units"), Some("K"));
        assert_eq!(var.get_attribute("missing"), None);
    }
    #[test]
    fn test_nc4_group_add_variable() {
        let mut group = Nc4Group::new("trajectory");
        group.add_variable(Nc4Variable::float64("x", vec![], vec![1.0, 2.0]));
        group.add_variable(Nc4Variable::float64("y", vec![], vec![3.0, 4.0]));
        assert_eq!(group.variables.len(), 2);
        let xv = group.get_variable("x").unwrap();
        assert!((xv.float_data[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_nc4_group_add_subgroup() {
        let mut root = Nc4Group::new("/");
        let mut child = Nc4Group::new("data");
        child.add_variable(Nc4Variable::float64("t", vec![], vec![0.0]));
        root.add_subgroup(child);
        assert_eq!(root.subgroups.len(), 1);
        assert_eq!(root.subgroups[0].name, "data");
    }
    #[test]
    fn test_nc4_group_total_variable_count() {
        let mut root = Nc4Group::new("/");
        root.add_variable(Nc4Variable::float64("a", vec![], vec![]));
        let mut child = Nc4Group::new("child");
        child.add_variable(Nc4Variable::float64("b", vec![], vec![]));
        child.add_variable(Nc4Variable::float64("c", vec![], vec![]));
        root.add_subgroup(child);
        assert_eq!(root.total_variable_count(), 3);
    }
    #[test]
    fn test_nc4_group_all_variables() {
        let mut root = Nc4Group::new("/");
        root.add_variable(Nc4Variable::float64("a", vec![], vec![]));
        let mut child = Nc4Group::new("child");
        child.add_variable(Nc4Variable::float64("b", vec![], vec![]));
        root.add_subgroup(child);
        let all = root.all_variables();
        assert_eq!(all.len(), 2);
    }
    #[test]
    fn test_nc4_group_attributes() {
        let mut group = Nc4Group::new("meta");
        group.add_attribute("author", "test");
        assert_eq!(group.get_attribute("author"), Some("test"));
        assert_eq!(group.get_attribute("missing"), None);
    }
    #[test]
    fn test_nc4_file_add_dimension() {
        let mut f = Nc4File::new();
        f.add_dimension("atom", 100);
        f.add_unlimited_dimension("frame");
        assert_eq!(f.get_dimension_size("atom"), Some(100));
        assert_eq!(f.get_dimension_size("frame"), None);
    }
    #[test]
    fn test_nc4_file_extend_unlimited() {
        let mut f = Nc4File::new();
        f.add_unlimited_dimension("frame");
        f.extend_unlimited();
        f.extend_unlimited();
        assert_eq!(f.unlimited_size, 2);
    }
    #[test]
    fn test_nc4_file_global_attributes() {
        let mut f = Nc4File::new();
        f.add_global_attribute("Conventions", "CF-1.8");
        f.add_global_attribute("title", "Test");
        assert_eq!(f.get_global_attribute("Conventions"), Some("CF-1.8"));
        assert_eq!(f.get_global_attribute("missing"), None);
    }
    #[test]
    fn test_nc4_file_add_get_variable() {
        let mut f = Nc4File::new();
        f.add_variable(Nc4Variable::float64(
            "pressure",
            vec!["time".into()],
            vec![101.3, 100.5],
        ));
        let v = f.get_variable("pressure").unwrap();
        assert_eq!(v.float_data.len(), 2);
        assert!(f.get_variable("missing").is_none());
    }
    #[test]
    fn test_nc4_file_add_group() {
        let mut f = Nc4File::new();
        let mut g = Nc4Group::new("simdata");
        g.add_variable(Nc4Variable::int32("step", vec![], vec![0, 1, 2]));
        f.add_group(g);
        let retrieved = f.get_group("simdata").unwrap();
        assert_eq!(retrieved.variables.len(), 1);
        assert!(f.get_group("missing").is_none());
    }
    #[test]
    fn test_nc4_file_total_variable_count() {
        let mut f = Nc4File::new();
        f.add_variable(Nc4Variable::float64("a", vec![], vec![]));
        let mut g = Nc4Group::new("sub");
        g.add_variable(Nc4Variable::float64("b", vec![], vec![]));
        f.add_group(g);
        assert_eq!(f.total_variable_count(), 2);
    }
    #[test]
    fn test_nc4_file_to_bytes_non_empty() {
        let mut f = Nc4File::new();
        f.add_global_attribute("test", "value");
        f.add_variable(Nc4Variable::float64("x", vec![], vec![1.0, 2.0, 3.0]));
        let bytes = f.to_bytes();
        assert!(!bytes.is_empty());
        assert_eq!(&bytes[0..5], b"OXNC4");
    }
    #[test]
    fn test_create_string_variable() {
        let var = create_string_variable(
            "residue_name",
            vec!["atom".into()],
            vec!["ALA".to_string(), "GLY".to_string(), "PRO".to_string()],
        );
        assert_eq!(var.name, "residue_name");
        assert_eq!(var.string_data.len(), 3);
        assert_eq!(var.string_data[1], "GLY");
        assert_eq!(var.get_attribute("dtype"), Some("string"));
    }
    #[test]
    fn test_inherit_global_attributes() {
        let mut source = NetCdfFile::new();
        source.add_attribute("Conventions", "AMBER");
        source.add_attribute("program", "OxiPhysics");
        let mut target = NetCdfFile::new();
        let mut var = NetCdfVariable::new("x", vec![], vec![1.0]);
        var.add_attribute("Conventions", "CF");
        target.add_variable(var);
        inherit_global_attributes(&source, &mut target);
        let xv = target.get_variable("x").unwrap();
        assert_eq!(xv.get_attribute("Conventions"), Some("CF"));
        assert_eq!(xv.get_attribute("program"), Some("OxiPhysics"));
    }
    #[test]
    fn test_merge_netcdf_files() {
        let mut base = NetCdfFile::new();
        base.add_dimension("time", 3);
        base.add_variable(NetCdfVariable::new(
            "a",
            vec!["time".into()],
            vec![1.0, 2.0, 3.0],
        ));
        base.add_variable(NetCdfVariable::new(
            "b",
            vec!["time".into()],
            vec![4.0, 5.0, 6.0],
        ));
        base.add_attribute("source", "base");
        let mut overlay = NetCdfFile::new();
        overlay.add_variable(NetCdfVariable::new(
            "a",
            vec!["time".into()],
            vec![10.0, 20.0, 30.0],
        ));
        overlay.add_variable(NetCdfVariable::new("c", vec![], vec![99.0]));
        overlay.add_attribute("source", "overlay");
        let merged = merge_netcdf_files(&base, &overlay);
        assert_eq!(merged.variable_names().len(), 3);
        let av = merged.get_variable("a").unwrap();
        assert!((av.data[0] - 10.0).abs() < 1e-10);
        let bv = merged.get_variable("b").unwrap();
        assert!((bv.data[0] - 4.0).abs() < 1e-10);
        assert_eq!(merged.get_attribute("source"), Some("overlay"));
    }
    #[test]
    fn test_netcdf_file_new_empty() {
        let f = NetcdfFile::new();
        assert!(f.dimensions.is_empty());
        assert!(f.variables.is_empty());
        assert!(f.global_attrs.is_empty());
    }
    #[test]
    fn test_netcdf_file_add_dimension() {
        let mut f = NetcdfFile::new();
        f.add_dimension("time", 10);
        f.add_dimension("space", 3);
        assert_eq!(f.dimensions.get("time"), Some(&10));
        assert_eq!(f.dimensions.get("space"), Some(&3));
    }
    #[test]
    fn test_netcdf_file_add_variable_increases_count() {
        let mut f = NetcdfFile::new();
        assert_eq!(f.variables.len(), 0);
        f.add_variable(NetcdfVariable {
            name: "temp".to_string(),
            dimensions: vec!["time".to_string()],
            data: vec![1.0, 2.0, 3.0],
            units: "K".to_string(),
            long_name: "temperature".to_string(),
        });
        assert_eq!(f.variables.len(), 1);
        f.add_variable(NetcdfVariable {
            name: "pressure".to_string(),
            dimensions: vec!["time".to_string()],
            data: vec![100.0],
            units: "Pa".to_string(),
            long_name: "pressure".to_string(),
        });
        assert_eq!(f.variables.len(), 2);
    }
    #[test]
    fn test_netcdf_file_add_global_attr() {
        let mut f = NetcdfFile::new();
        f.add_global_attr("Conventions", "CF-1.8");
        f.add_global_attr("title", "Test");
        assert_eq!(f.global_attrs.len(), 2);
        assert!(
            f.global_attrs
                .iter()
                .any(|(k, v)| k == "Conventions" && v == "CF-1.8")
        );
    }
    #[test]
    fn test_netcdf_file_write_cdl_contains_dimension_name() {
        let mut f = NetcdfFile::new();
        f.add_dimension("latitude", 180);
        f.add_dimension("longitude", 360);
        let cdl = f.write_cdl();
        assert!(cdl.contains("latitude"), "CDL missing dimension 'latitude'");
        assert!(cdl.contains("180"), "CDL missing dimension size 180");
        assert!(cdl.contains("longitude"));
        assert!(cdl.contains("360"));
    }
    #[test]
    fn test_netcdf_file_write_cdl_contains_variable() {
        let mut f = NetcdfFile::new();
        f.add_dimension("t", 5);
        f.add_variable(NetcdfVariable {
            name: "velocity".to_string(),
            dimensions: vec!["t".to_string()],
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            units: "m/s".to_string(),
            long_name: "wind velocity".to_string(),
        });
        let cdl = f.write_cdl();
        assert!(cdl.contains("velocity"));
        assert!(cdl.contains("m/s"));
        assert!(cdl.contains("wind velocity"));
    }
    #[test]
    fn test_parse_cdl_dimensions() {
        let cdl = "netcdf test {\ndimensions:\n\ttime = 10 ;\n\tatom = 100 ;\nvariables:\n}\n";
        let dims = parse_cdl_dimensions(cdl);
        assert_eq!(dims.len(), 2, "expected 2 dimensions, got {:?}", dims);
        assert!(dims.iter().any(|(n, s)| n == "time" && *s == 10));
        assert!(dims.iter().any(|(n, s)| n == "atom" && *s == 100));
    }
    #[test]
    fn test_trajectory_stats_n_frames() {
        let times = vec![0.0, 0.5, 1.0, 1.5, 2.0];
        let positions: Vec<Vec<[f64; 3]>> = times.iter().map(|_| vec![[0.0, 0.0, 0.0]]).collect();
        let (_total_time, _avg_disp, n_frames) = trajectory_stats(&times, &positions);
        assert_eq!(n_frames, 5);
    }
    #[test]
    fn test_trajectory_stats_total_time() {
        let times = vec![0.0, 1.0, 2.0];
        let positions: Vec<Vec<[f64; 3]>> = times.iter().map(|_| vec![[0.0, 0.0, 0.0]]).collect();
        let (total_time, _, _) = trajectory_stats(&times, &positions);
        assert!((total_time - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_trajectory_stats_avg_displacement() {
        let times = vec![0.0, 1.0];
        let positions = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[3.0, 4.0, 0.0], [1.0, 0.0, 0.0]],
        ];
        let (_total_time, avg_disp, _) = trajectory_stats(&times, &positions);
        assert!((avg_disp - 2.5).abs() < 1e-10, "avg_disp={}", avg_disp);
    }
    #[test]
    fn test_trajectory_write_creates_file() {
        let times = vec![0.0, 0.5, 1.0];
        let positions = vec![
            vec![[1.0, 2.0, 3.0]],
            vec![[1.1, 2.1, 3.1]],
            vec![[1.2, 2.2, 3.2]],
        ];
        let path = "/tmp/test_trajectory.cdl";
        NetcdfFile::trajectory_write(path, &times, &positions).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("frame"), "missing 'frame' dimension");
        assert!(content.contains("atom"));
        assert!(content.contains("coordinates"));
        assert!(content.contains("time"));
    }
}
/// Write a simple particle trajectory in CDL format using [`NetcdfWriter`].
///
/// Dimensions: `time` (unlimited), `atom`, `spatial` (3).
/// Variables: `time(time)`, `x(time, atom)`, `y(time, atom)`, `z(time, atom)`.
#[allow(dead_code)]
pub fn write_particle_trajectory_cdl(
    times: &[f64],
    positions: &[Vec<[f64; 3]>],
    units_time: &str,
    units_pos: &str,
) -> String {
    let n_frames = times.len();
    let n_atoms = positions.first().map(|f| f.len()).unwrap_or(0);
    let mut w = NetcdfWriter::new("trajectory");
    w.add_unlimited_dimension("time", n_frames);
    w.add_dimension("atom", n_atoms);
    w.add_dimension("spatial", 3);
    w.add_global_attribute("Conventions", "OxiPhysics");
    w.add_global_attribute("n_frames", &format!("{}", n_frames));
    w.add_global_attribute("n_atoms", &format!("{}", n_atoms));
    w.add_coordinate("time", times.to_vec(), units_time);
    let x_data: Vec<f64> = positions
        .iter()
        .flat_map(|f| f.iter().map(|p| p[0]))
        .collect();
    let y_data: Vec<f64> = positions
        .iter()
        .flat_map(|f| f.iter().map(|p| p[1]))
        .collect();
    let z_data: Vec<f64> = positions
        .iter()
        .flat_map(|f| f.iter().map(|p| p[2]))
        .collect();
    w.add_variable("x", &["time", "atom"], x_data);
    w.add_variable_attribute("units", units_pos);
    w.add_variable_attribute("long_name", "x coordinate");
    w.add_variable("y", &["time", "atom"], y_data);
    w.add_variable_attribute("units", units_pos);
    w.add_variable_attribute("long_name", "y coordinate");
    w.add_variable("z", &["time", "atom"], z_data);
    w.add_variable_attribute("units", units_pos);
    w.add_variable_attribute("long_name", "z coordinate");
    w.finish()
}
/// Parse a particle trajectory CDL back into time and position arrays.
///
/// Returns `(times, positions)` where `positions\[frame\]\[atom\] = \[x, y, z\]`.
#[allow(dead_code)]
pub fn read_particle_trajectory_cdl(cdl: &str) -> Result<(Vec<f64>, Vec<Vec<[f64; 3]>>), String> {
    let reader = NetcdfReader::from_cdl(cdl)?;
    let times = reader
        .get_data("time")
        .ok_or("missing 'time' variable")?
        .to_vec();
    let n_frames = times.len();
    let n_atoms = reader.get_dimension("atom").unwrap_or(0);
    let x_data = reader.get_data("x").ok_or("missing 'x' variable")?.to_vec();
    let y_data = reader.get_data("y").ok_or("missing 'y' variable")?.to_vec();
    let z_data = reader.get_data("z").ok_or("missing 'z' variable")?.to_vec();
    let mut positions = Vec::with_capacity(n_frames);
    for frame in 0..n_frames {
        let mut frame_pos = Vec::with_capacity(n_atoms);
        for atom in 0..n_atoms {
            let idx = frame * n_atoms + atom;
            let x = x_data.get(idx).copied().unwrap_or(0.0);
            let y = y_data.get(idx).copied().unwrap_or(0.0);
            let z = z_data.get(idx).copied().unwrap_or(0.0);
            frame_pos.push([x, y, z]);
        }
        positions.push(frame_pos);
    }
    Ok((times, positions))
}
/// Build a regularly-spaced time coordinate (start, step, n).
#[allow(dead_code)]
pub fn linspace_time(start: f64, step: f64, n: usize) -> Vec<f64> {
    (0..n).map(|i| start + i as f64 * step).collect()
}
/// Build X/Y/Z Cartesian coordinate arrays for a uniform grid.
#[allow(dead_code)]
pub fn uniform_grid_coords(
    nx: usize,
    ny: usize,
    nz: usize,
    origin: [f64; 3],
    spacing: [f64; 3],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = (0..nx).map(|i| origin[0] + i as f64 * spacing[0]).collect();
    let y: Vec<f64> = (0..ny).map(|j| origin[1] + j as f64 * spacing[1]).collect();
    let z: Vec<f64> = (0..nz).map(|k| origin[2] + k as f64 * spacing[2]).collect();
    (x, y, z)
}
/// Common CF (Climate and Forecast) convention global attributes.
#[allow(dead_code)]
pub fn cf_global_attributes(title: &str, institution: &str, source: &str) -> Vec<(String, String)> {
    vec![
        ("Conventions".to_string(), "CF-1.8".to_string()),
        ("title".to_string(), title.to_string()),
        ("institution".to_string(), institution.to_string()),
        ("source".to_string(), source.to_string()),
        ("history".to_string(), "Created by OxiPhysics".to_string()),
    ]
}
/// Add CF convention attributes to a [`NetCdfFile`].
#[allow(dead_code)]
pub fn apply_cf_conventions(file: &mut NetCdfFile, title: &str) {
    file.add_attribute("Conventions", "CF-1.8");
    file.add_attribute("title", title);
    file.add_attribute("institution", "OxiPhysics");
}
/// Compute extended statistics for a variable.
#[allow(dead_code)]
pub fn extended_variable_stats(var: &NetCdfVariable) -> Option<NetcdfVariableStats> {
    if var.data.is_empty() {
        return None;
    }
    let n = var.data.len();
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0;
    for &v in &var.data {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
        sum += v;
    }
    let mean = sum / n as f64;
    let var_sq = var.data.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n as f64;
    Some(NetcdfVariableStats {
        name: var.name.clone(),
        min,
        max,
        mean,
        std_dev: var_sq.sqrt(),
        count: n,
    })
}
#[cfg(test)]
mod tests_netcdf_new {
    use super::*;

    #[test]
    fn writer_basic_cdl() {
        let mut w = NetcdfWriter::new("test");
        w.add_dimension("time", 3);
        w.add_variable("temperature", &["time"], vec![300.0, 301.0, 302.0]);
        w.add_variable_attribute("units", "K");
        let cdl = w.finish();
        assert!(cdl.starts_with("netcdf test {"));
        assert!(cdl.contains("time = 3"));
        assert!(cdl.contains("double temperature(time)"));
        assert!(cdl.contains("temperature:units = \"K\""));
        assert!(cdl.contains("300"));
    }
    #[test]
    fn writer_unlimited_dimension() {
        let mut w = NetcdfWriter::new("sim");
        w.add_unlimited_dimension("frame", 5);
        w.add_dimension("atom", 10);
        let cdl = w.to_cdl();
        assert!(cdl.contains("UNLIMITED"));
        assert!(cdl.contains("5 currently"));
    }
    #[test]
    fn writer_global_attributes() {
        let mut w = NetcdfWriter::new("data");
        w.add_global_attribute("Conventions", "CF-1.8");
        w.add_global_attribute("title", "Test");
        let cdl = w.to_cdl();
        assert!(cdl.contains(":Conventions = \"CF-1.8\""));
        assert!(cdl.contains(":title = \"Test\""));
    }
    #[test]
    fn writer_coordinate_variable() {
        let mut w = NetcdfWriter::new("grid");
        w.add_dimension("x", 4);
        w.add_coordinate("x", vec![0.0, 1.0, 2.0, 3.0], "m");
        let cdl = w.to_cdl();
        assert!(cdl.contains("x:units = \"m\""));
        assert!(cdl.contains("0, 1, 2, 3"));
    }
    #[test]
    fn writer_n_dimensions_and_variables() {
        let mut w = NetcdfWriter::new("d");
        w.add_dimension("a", 1);
        w.add_dimension("b", 2);
        w.add_variable("v", &["a"], vec![1.0]);
        assert_eq!(w.n_dimensions(), 2);
        assert_eq!(w.n_variables(), 1);
    }
    #[test]
    fn writer_write_to_file() {
        let mut w = NetcdfWriter::new("out");
        w.add_dimension("t", 2);
        w.add_variable("time", &["t"], vec![0.0, 1.0]);
        let path = "/tmp/test_netcdf_writer.cdl";
        w.write_to_file(path).unwrap();
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("netcdf out"));
        assert!(content.contains("time = 0, 1"));
    }
    #[test]
    fn reader_from_cdl_roundtrip() {
        let mut w = NetcdfWriter::new("rt");
        w.add_dimension("n", 3);
        w.add_variable("pressure", &["n"], vec![100.0, 101.0, 102.0]);
        w.add_variable_attribute("units", "Pa");
        w.add_global_attribute("source", "test");
        let cdl = w.finish();
        let reader = NetcdfReader::from_cdl(&cdl).unwrap();
        assert_eq!(reader.get_dimension("n"), Some(3));
        let p = reader.get_data("pressure").unwrap();
        assert_eq!(p.len(), 3);
        assert!((p[1] - 101.0).abs() < 1e-9);
    }
    #[test]
    fn reader_variable_names() {
        let mut w = NetcdfWriter::new("x");
        w.add_dimension("k", 2);
        w.add_variable("alpha", &["k"], vec![1.0, 2.0]);
        w.add_variable("beta", &["k"], vec![3.0, 4.0]);
        let cdl = w.finish();
        let reader = NetcdfReader::from_cdl(&cdl).unwrap();
        let names = reader.variable_names();
        assert!(names.contains(&"alpha"));
        assert!(names.contains(&"beta"));
    }
    #[test]
    fn reader_missing_variable_returns_none() {
        let w = NetcdfWriter::new("empty");
        let cdl = w.finish();
        let reader = NetcdfReader::from_cdl(&cdl).unwrap();
        assert!(reader.get_variable("nonexistent").is_none());
        assert!(reader.get_data("nonexistent").is_none());
    }
    #[test]
    fn particle_trajectory_write_read_roundtrip() {
        let times = vec![0.0, 1.0, 2.0];
        let positions = vec![
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            vec![[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]],
            vec![[1.2, 2.2, 3.2], [4.2, 5.2, 6.2]],
        ];
        let cdl = write_particle_trajectory_cdl(&times, &positions, "ps", "angstrom");
        assert!(cdl.contains("n_frames"));
        assert!(cdl.contains("n_atoms"));
        assert!(cdl.contains("OxiPhysics"));
        let (rt_times, rt_pos) = read_particle_trajectory_cdl(&cdl).unwrap();
        assert_eq!(rt_times.len(), 3);
        assert!((rt_times[2] - 2.0).abs() < 1e-9);
        assert_eq!(rt_pos.len(), 3);
        assert_eq!(rt_pos[0].len(), 2);
        assert!((rt_pos[1][0][0] - 1.1).abs() < 1e-9);
        assert!((rt_pos[2][1][2] - 6.2).abs() < 1e-9);
    }
    #[test]
    fn particle_trajectory_empty() {
        let cdl = write_particle_trajectory_cdl(&[], &[], "ps", "nm");
        assert!(cdl.contains("netcdf trajectory"));
        assert!(cdl.contains("dimensions:"));
        assert!(cdl.contains("variables:"));
    }
    #[test]
    fn linspace_time_basic() {
        let t = linspace_time(0.0, 0.5, 5);
        assert_eq!(t.len(), 5);
        assert!((t[0] - 0.0).abs() < 1e-12);
        assert!((t[4] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn linspace_time_single() {
        let t = linspace_time(3.0, 1.0, 1);
        assert_eq!(t, vec![3.0]);
    }
    #[test]
    fn uniform_grid_coords_sizes() {
        let (x, y, z) = uniform_grid_coords(4, 5, 6, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(x.len(), 4);
        assert_eq!(y.len(), 5);
        assert_eq!(z.len(), 6);
    }
    #[test]
    fn uniform_grid_coords_values() {
        let (x, _y, _z) = uniform_grid_coords(3, 1, 1, [1.0, 0.0, 0.0], [0.5, 1.0, 1.0]);
        assert!((x[0] - 1.0).abs() < 1e-12);
        assert!((x[1] - 1.5).abs() < 1e-12);
        assert!((x[2] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn cf_attrs_content() {
        let attrs = cf_global_attributes("My Title", "ACME Corp", "OxiPhysics 1.0");
        assert!(
            attrs
                .iter()
                .any(|(k, v)| k == "Conventions" && v.contains("CF"))
        );
        assert!(attrs.iter().any(|(k, _)| k == "title"));
        assert!(attrs.iter().any(|(k, _)| k == "institution"));
    }
    #[test]
    fn apply_cf_adds_attributes() {
        let mut file = NetCdfFile::new();
        apply_cf_conventions(&mut file, "My Dataset");
        assert_eq!(file.get_attribute("Conventions"), Some("CF-1.8"));
        assert_eq!(file.get_attribute("title"), Some("My Dataset"));
    }
    #[test]
    fn extended_stats_basic() {
        let var = NetCdfVariable::new("v", vec![], vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let s = extended_variable_stats(&var).unwrap();
        assert_eq!(s.count, 5);
        assert!((s.min - 1.0).abs() < 1e-12);
        assert!((s.max - 5.0).abs() < 1e-12);
        assert!((s.mean - 3.0).abs() < 1e-12);
        assert!((s.std_dev - 2.0f64.sqrt()).abs() < 1e-9);
    }
    #[test]
    fn extended_stats_range() {
        let var = NetCdfVariable::new("v", vec![], vec![-5.0, 5.0]);
        let s = extended_variable_stats(&var).unwrap();
        assert!((s.range() - 10.0).abs() < 1e-12);
    }
    #[test]
    fn extended_stats_empty() {
        let var = NetCdfVariable::new("empty", vec![], vec![]);
        assert!(extended_variable_stats(&var).is_none());
    }
    #[test]
    fn extended_stats_cv_zero_mean() {
        let var = NetCdfVariable::new("z", vec![], vec![0.0, 0.0, 0.0]);
        let s = extended_variable_stats(&var).unwrap();
        assert!((s.cv() - 0.0).abs() < f64::EPSILON * 10.0);
    }
    #[test]
    fn dimension_unlimited_check() {
        let d_fixed = NetCdfDimension {
            name: "x".to_string(),
            size: Some(10),
        };
        let d_unlim = NetCdfDimension {
            name: "t".to_string(),
            size: None,
        };
        assert!(!d_fixed.is_unlimited());
        assert!(d_unlim.is_unlimited());
        assert_eq!(d_fixed.effective_size(), 10);
        assert_eq!(d_unlim.effective_size(), 0);
    }
    #[test]
    fn merge_netcdf_files_basic() {
        let mut base = NetCdfFile::new();
        base.add_variable(NetCdfVariable::new("a", vec![], vec![1.0]));
        base.add_attribute("source", "base");
        let mut overlay = NetCdfFile::new();
        overlay.add_variable(NetCdfVariable::new("b", vec![], vec![2.0]));
        overlay.add_attribute("source", "overlay");
        let merged = merge_netcdf_files(&base, &overlay);
        assert!(merged.get_variable("a").is_some());
        assert!(merged.get_variable("b").is_some());
        assert_eq!(merged.get_attribute("source"), Some("overlay"));
    }
    #[test]
    fn inherit_global_attributes_basic() {
        let mut source = NetCdfFile::new();
        source.add_attribute("units", "SI");
        let mut target = NetCdfFile::new();
        target.add_variable(NetCdfVariable::new("T", vec![], vec![300.0]));
        inherit_global_attributes(&source, &mut target);
        let attr = target.get_variable("T").unwrap().get_attribute("units");
        assert_eq!(attr, Some("SI"));
    }
    #[test]
    fn reshape_2d_basic() {
        let data: Vec<f64> = (1..=6).map(|i| i as f64).collect();
        let r = reshape_2d(&data, 2, 3);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], vec![1.0, 2.0, 3.0]);
        assert_eq!(r[1], vec![4.0, 5.0, 6.0]);
    }
    #[test]
    fn variable_stats_basic() {
        let var = NetCdfVariable::new("v", vec![], vec![2.0, 4.0, 6.0]);
        let (min, max, mean) = variable_stats(&var).unwrap();
        assert!((min - 2.0).abs() < 1e-12);
        assert!((max - 6.0).abs() < 1e-12);
        assert!((mean - 4.0).abs() < 1e-12);
    }
    #[test]
    fn variable_stats_empty() {
        let var = NetCdfVariable::new("empty", vec![], vec![]);
        assert!(variable_stats(&var).is_none());
    }
    #[test]
    fn amber_trajectory_structure() {
        let n_atoms = 5;
        let n_frames = 3;
        let coords: Vec<f64> = (0..n_atoms * n_frames * 3).map(|i| i as f64).collect();
        let times = vec![0.0, 1.0, 2.0];
        let file = create_amber_trajectory(n_atoms, n_frames, coords, times, None, None);
        assert_eq!(file.get_dimension_size("atom"), Some(n_atoms));
        assert_eq!(file.get_dimension_size("spatial"), Some(3));
        assert_eq!(file.get_attribute("Conventions"), Some("AMBER"));
        let coords_var = file.get_variable("coordinates").unwrap();
        assert_eq!(coords_var.data.len(), n_atoms * n_frames * 3);
    }
    #[test]
    fn amber_extract_frame_coordinates() {
        let n_atoms = 2;
        let n_frames = 2;
        let coords = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let file = create_amber_trajectory(n_atoms, n_frames, coords, vec![0.0, 1.0], None, None);
        let frame0 = extract_frame_coordinates(&file, 0).unwrap();
        assert_eq!(frame0.len(), 2);
        assert!((frame0[0][0] - 1.0).abs() < 1e-9);
        assert!((frame0[1][2] - 6.0).abs() < 1e-9);
        let frame1 = extract_frame_coordinates(&file, 1).unwrap();
        assert!((frame1[0][0] - 7.0).abs() < 1e-9);
    }
    #[test]
    fn amber_extract_out_of_range() {
        let file = create_amber_trajectory(2, 1, vec![0.0; 6], vec![0.0], None, None);
        assert!(extract_frame_coordinates(&file, 99).is_none());
    }
    #[test]
    fn coordinate_variable_units_attribute() {
        let var = create_coordinate_variable("time", vec![0.0, 1.0, 2.0], "ps");
        assert_eq!(var.get_attribute("units"), Some("ps"));
        assert_eq!(var.dims, vec!["time"]);
    }
    #[test]
    fn nc4_file_to_bytes_contains_magic() {
        let mut f = Nc4File::new();
        f.add_global_attribute("title", "Test");
        f.add_dimension("x", 10);
        let bytes = f.to_bytes();
        assert_eq!(&bytes[..5], b"OXNC4");
    }
    #[test]
    fn nc4_file_unlimited_extend() {
        let mut f = Nc4File::new();
        f.add_unlimited_dimension("time");
        f.extend_unlimited();
        f.extend_unlimited();
        assert_eq!(f.unlimited_size, 2);
        assert_eq!(f.get_dimension_size("time"), Some(2));
    }
    #[test]
    fn nc4_group_hierarchy() {
        let mut root = Nc4Group::new("root");
        let mut child = Nc4Group::new("child");
        child.add_variable(Nc4Variable::float64("v", vec![], vec![1.0, 2.0]));
        root.add_subgroup(child);
        assert_eq!(root.total_variable_count(), 1);
        let all = root.all_variables();
        assert_eq!(all.len(), 1);
        assert_eq!(all[0].name, "v");
    }
}
#[cfg(test)]
#[allow(dead_code)]
mod tests_netcdf_extra {
    use super::*;

    #[test]
    fn ncdf_variable_new_empty_data() {
        let var = NetCdfVariable::new("pressure", vec!["time".to_string()], vec![]);
        assert!(var.is_empty());
        assert_eq!(var.len(), 0);
    }
    #[test]
    fn ncdf_variable_add_and_get_attribute() {
        let mut var = NetCdfVariable::new("temp", vec![], vec![1.0]);
        var.add_attribute("units", "K");
        assert_eq!(var.get_attribute("units"), Some("K"));
        assert_eq!(var.get_attribute("missing"), None);
    }
    #[test]
    fn ncdf_variable_multiple_attributes() {
        let mut var = NetCdfVariable::new("vel", vec!["atom".to_string()], vec![1.0, 2.0]);
        var.add_attribute("units", "angstrom/ps");
        var.add_attribute("long_name", "velocity");
        assert_eq!(var.get_attribute("long_name"), Some("velocity"));
        assert_eq!(var.attributes.len(), 2);
    }
    #[test]
    fn ncdf_dimension_unlimited() {
        let dim = NetCdfDimension {
            name: "time".to_string(),
            size: None,
        };
        assert!(dim.is_unlimited());
        assert_eq!(dim.effective_size(), 0);
    }
    #[test]
    fn ncdf_dimension_fixed() {
        let dim = NetCdfDimension {
            name: "atom".to_string(),
            size: Some(100),
        };
        assert!(!dim.is_unlimited());
        assert_eq!(dim.effective_size(), 100);
    }
    #[test]
    fn netcdf_file_add_and_get_variable() {
        let mut f = NetCdfFile::new();
        let var = NetCdfVariable::new("coords", vec!["atom".to_string()], vec![1.0, 2.0, 3.0]);
        f.add_variable(var);
        let got = f.get_variable("coords").unwrap();
        assert_eq!(got.data, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn netcdf_file_get_variable_missing() {
        let f = NetCdfFile::new();
        assert!(f.get_variable("missing").is_none());
    }
    #[test]
    fn netcdf_file_add_dimension() {
        let mut f = NetCdfFile::new();
        f.add_dimension("atom", 50);
        assert_eq!(f.get_dimension_size("atom"), Some(50));
    }
    #[test]
    fn netcdf_file_unlimited_dimension() {
        let mut f = NetCdfFile::new();
        f.add_unlimited_dimension("time", 0);
        assert!(f.is_unlimited_dimension("time"));
        assert!(!f.is_unlimited_dimension("atom"));
    }
    #[test]
    fn netcdf_file_add_global_attribute() {
        let mut f = NetCdfFile::new();
        f.add_attribute("Conventions", "AMBER");
        assert_eq!(f.get_attribute("Conventions"), Some("AMBER"));
    }
    #[test]
    fn netcdf_file_variable_names() {
        let mut f = NetCdfFile::new();
        f.add_variable(NetCdfVariable::new("time", vec![], vec![]));
        f.add_variable(NetCdfVariable::new("coords", vec![], vec![]));
        let names = f.variable_names();
        assert!(names.contains(&"time"));
        assert!(names.contains(&"coords"));
    }
    #[test]
    fn netcdf_file_dimension_names() {
        let mut f = NetCdfFile::new();
        f.add_dimension("frame", 10);
        f.add_dimension("atom", 5);
        let names = f.dimension_names();
        assert_eq!(names.len(), 2);
    }
    #[test]
    fn ncf_text_roundtrip_basic() {
        let mut f = NetCdfFile::new();
        f.add_dimension("time", 3);
        let mut var = NetCdfVariable::new(
            "temperature",
            vec!["time".to_string()],
            vec![300.0, 310.0, 320.0],
        );
        var.add_attribute("units", "K");
        f.add_variable(var);
        let cdl = write_ncf_text(&f);
        let parsed = parse_ncf_text(&cdl).unwrap();
        let v = parsed.get_variable("temperature").unwrap();
        assert_eq!(v.data, vec![300.0, 310.0, 320.0]);
    }
    #[test]
    fn ncf_text_contains_dimension() {
        let mut f = NetCdfFile::new();
        f.add_dimension("atom", 42);
        let cdl = write_ncf_text(&f);
        assert!(cdl.contains("atom"));
    }
    #[test]
    fn ncf_text_contains_variable_attribute() {
        let mut f = NetCdfFile::new();
        let mut var = NetCdfVariable::new("v", vec![], vec![1.0]);
        var.add_attribute("units", "m/s");
        f.add_variable(var);
        let cdl = write_ncf_text(&f);
        assert!(cdl.contains("m/s"));
    }
    #[test]
    fn reshape_2d_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mat = reshape_2d(&data, 2, 3);
        assert_eq!(mat.len(), 2);
        assert_eq!(mat[0], vec![1.0, 2.0, 3.0]);
        assert_eq!(mat[1], vec![4.0, 5.0, 6.0]);
    }
    #[test]
    fn reshape_2d_single_row() {
        let data = vec![1.0, 2.0, 3.0];
        let mat = reshape_2d(&data, 1, 3);
        assert_eq!(mat, vec![vec![1.0, 2.0, 3.0]]);
    }
    #[test]
    fn reshape_2d_zero_rows() {
        let mat = reshape_2d(&[], 0, 3);
        assert!(mat.is_empty());
    }
    #[test]
    fn variable_stats_basic() {
        let var = NetCdfVariable::new("t", vec![], vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let (min, max, mean) = variable_stats(&var).unwrap();
        assert!((min - 1.0).abs() < 1e-12);
        assert!((max - 5.0).abs() < 1e-12);
        assert!((mean - 3.0).abs() < 1e-12);
    }
    #[test]
    fn variable_stats_empty_returns_none() {
        let var = NetCdfVariable::new("t", vec![], vec![]);
        assert!(variable_stats(&var).is_none());
    }
    #[test]
    fn amber_trajectory_n_frames_dimension() {
        let file = create_amber_trajectory(3, 2, vec![0.0; 18], vec![0.0, 1.0], None, None);
        assert_eq!(file.get_dimension_size("frame"), Some(2));
        assert_eq!(file.get_dimension_size("atom"), Some(3));
    }
    #[test]
    fn amber_trajectory_has_time_variable() {
        let file = create_amber_trajectory(2, 3, vec![0.0; 18], vec![0.0, 1.0, 2.0], None, None);
        assert!(file.get_variable("time").is_some());
    }
    #[test]
    fn amber_trajectory_has_coordinates_variable() {
        let file = create_amber_trajectory(
            2,
            1,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![0.0],
            None,
            None,
        );
        assert!(file.get_variable("coordinates").is_some());
    }
    #[test]
    fn amber_trajectory_amber_conventions_attribute() {
        let file = create_amber_trajectory(1, 1, vec![0.0; 3], vec![0.0], None, None);
        let conv = file.get_attribute("Conventions");
        assert!(conv.is_some());
    }
    #[test]
    fn coordinate_variable_data() {
        let var = create_coordinate_variable("x", vec![0.0, 1.0, 2.0], "nm");
        assert_eq!(var.data, vec![0.0, 1.0, 2.0]);
        assert_eq!(var.name, "x");
    }
    #[test]
    fn coordinate_variable_units_attr() {
        let var = create_coordinate_variable("time", vec![0.0, 0.5], "ps");
        assert_eq!(var.get_attribute("units"), Some("ps"));
    }
    #[test]
    fn nc4_data_type_float64_width() {
        assert_eq!(Nc4DataType::Float64.byte_width(), 8);
    }
    #[test]
    fn nc4_data_type_int32_width() {
        assert_eq!(Nc4DataType::Int32.byte_width(), 4);
    }
    #[test]
    fn nc4_data_type_type_names() {
        assert_eq!(Nc4DataType::Float64.type_name(), "double");
        assert_eq!(Nc4DataType::Int32.type_name(), "int");
    }
    #[test]
    fn nc4_variable_float64_constructor() {
        let v = Nc4Variable::float64("coords", vec!["atom".to_string()], vec![1.0, 2.0]);
        assert_eq!(v.name, "coords");
        assert_eq!(v.data_type, Nc4DataType::Float64);
        assert_eq!(v.float_data, vec![1.0, 2.0]);
    }
    #[test]
    fn nc4_variable_int32_constructor() {
        let v = Nc4Variable::int32("n", vec![], vec![1, 2, 3]);
        assert_eq!(v.data_type, Nc4DataType::Int32);
        assert_eq!(v.int_data, vec![1_i64, 2, 3]);
    }
    #[test]
    fn nc4_variable_add_get_attribute() {
        let mut v = Nc4Variable::float64("t", vec![], vec![]);
        v.add_attribute("units", "K");
        assert_eq!(v.get_attribute("units"), Some("K"));
    }
    #[test]
    fn nc4_variable_len_is_empty() {
        let v = Nc4Variable::float64("x", vec![], vec![1.0, 2.0, 3.0]);
        assert_eq!(v.len(), 3);
        assert!(!v.is_empty());
        let empty = Nc4Variable::float64("y", vec![], vec![]);
        assert!(empty.is_empty());
    }
    #[test]
    fn nc4_file_add_and_get_variable() {
        let mut f = Nc4File::new();
        f.add_variable(Nc4Variable::float64("temperature", vec![], vec![300.0]));
        let v = f.get_variable("temperature").unwrap();
        assert_eq!(v.float_data, vec![300.0]);
    }
    #[test]
    fn nc4_file_total_variable_count() {
        let mut f = Nc4File::new();
        f.add_variable(Nc4Variable::float64("a", vec![], vec![]));
        f.add_variable(Nc4Variable::float64("b", vec![], vec![]));
        assert_eq!(f.total_variable_count(), 2);
    }
    #[test]
    fn nc4_file_add_group_and_retrieve() {
        let mut f = Nc4File::new();
        let mut g = Nc4Group::new("traj");
        g.add_variable(Nc4Variable::float64("v", vec![], vec![1.0]));
        f.add_group(g);
        assert!(f.get_group("traj").is_some());
        assert!(f.get_group("missing").is_none());
    }
    #[test]
    fn nc4_file_global_attribute_roundtrip() {
        let mut f = Nc4File::new();
        f.add_global_attribute("title", "Test Run");
        assert_eq!(f.get_global_attribute("title"), Some("Test Run"));
    }
    #[test]
    fn nc4_file_is_unlimited_fixed_dim() {
        let mut f = Nc4File::new();
        f.add_dimension("atom", 50);
        assert!(!f.is_unlimited("atom"));
    }
    #[test]
    fn merge_netcdf_files_overlay_wins() {
        let mut base = NetCdfFile::new();
        base.add_variable(NetCdfVariable::new("temp", vec![], vec![100.0]));
        let mut overlay = NetCdfFile::new();
        overlay.add_variable(NetCdfVariable::new("temp", vec![], vec![200.0]));
        overlay.add_variable(NetCdfVariable::new("pressure", vec![], vec![1.0]));
        let merged = merge_netcdf_files(&base, &overlay);
        let temp = merged.get_variable("temp").unwrap();
        assert_eq!(temp.data, vec![200.0]);
        assert!(merged.get_variable("pressure").is_some());
    }
    #[test]
    fn parse_cdl_dimensions_basic() {
        let cdl = "netcdf x {\ndimensions:\n\tatom = 5 ;\n\tframe = 10 ;\nvariables:\n}";
        let dims = parse_cdl_dimensions(cdl);
        assert_eq!(dims.len(), 2);
        let names: Vec<&str> = dims.iter().map(|(n, _)| n.as_str()).collect();
        assert!(names.contains(&"atom"));
        assert!(names.contains(&"frame"));
    }
    #[test]
    fn parse_cdl_dimensions_empty() {
        let cdl = "netcdf x {\ndimensions:\nvariables:\n}";
        let dims = parse_cdl_dimensions(cdl);
        assert!(dims.is_empty());
    }
    #[test]
    fn netcdffile_hashmap_add_dimension() {
        let mut f = NetcdfFile::new();
        f.add_dimension("atom", 10);
        assert_eq!(f.dimensions.get("atom"), Some(&10));
    }
    #[test]
    fn netcdffile_write_cdl_contains_dims_and_vars() {
        let mut f = NetcdfFile::new();
        f.add_dimension("time", 5);
        f.variables.push(NetcdfVariable {
            name: "temperature".to_string(),
            dimensions: vec!["time".to_string()],
            data: vec![300.0; 5],
            units: "K".to_string(),
            long_name: "Temperature".to_string(),
        });
        let cdl = f.write_cdl();
        assert!(cdl.contains("time = 5"));
        assert!(cdl.contains("temperature"));
        assert!(cdl.contains("K"));
    }
    #[test]
    fn linspace_time_basic() {
        let t = linspace_time(0.0, 0.5, 5);
        assert_eq!(t.len(), 5);
        assert!((t[0] - 0.0).abs() < 1e-12);
        assert!((t[1] - 0.5).abs() < 1e-12);
        assert!((t[4] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn linspace_time_empty() {
        let t = linspace_time(0.0, 1.0, 0);
        assert!(t.is_empty());
    }
    #[test]
    fn cf_global_attributes_has_required_keys() {
        let attrs = cf_global_attributes("My Run", "UniX", "MD");
        let keys: Vec<&str> = attrs.iter().map(|(k, _)| k.as_str()).collect();
        assert!(keys.contains(&"title") || attrs.iter().any(|(_, v)| v.contains("My Run")));
    }
}
/// Write a CDL dimensions block from a list of `NetcdfDimSpec`.
#[allow(dead_code)]
pub fn write_cdl_dimensions<W: std::io::Write>(
    writer: &mut W,
    dims: &[NetcdfDimSpec],
) -> std::io::Result<()> {
    writeln!(writer, "dimensions:")?;
    for d in dims {
        writeln!(writer, "{}", d.to_cdl())?;
    }
    Ok(())
}
/// Wrap positions from Å → nm (divide by 10).
#[allow(dead_code)]
pub fn angstroms_to_nm(positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|p| [p[0] / 10.0, p[1] / 10.0, p[2] / 10.0])
        .collect()
}
/// Wrap positions from nm → Å (multiply by 10).
#[allow(dead_code)]
pub fn nm_to_angstroms(positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    positions
        .iter()
        .map(|p| [p[0] * 10.0, p[1] * 10.0, p[2] * 10.0])
        .collect()
}
/// Generate a linearly spaced time vector from `t0` to `t0 + (n-1)*dt`.
/// (Same as linspace_time but with explicit start time.)
#[allow(dead_code)]
pub fn time_vector(t0: f64, dt: f64, n: usize) -> Vec<f64> {
    (0..n).map(|i| t0 + i as f64 * dt).collect()
}
/// Compute mean-square displacement (MSD) for a 1-D position series.
/// `positions\[i\]` is the 1-D coordinate at frame `i`.
#[allow(dead_code)]
pub fn msd_1d(positions: &[f64]) -> Vec<f64> {
    if positions.is_empty() {
        return vec![];
    }
    let n = positions.len();
    let p0 = positions[0];
    positions
        .iter()
        .map(|&p| (p - p0).powi(2))
        .collect::<Vec<_>>()[..n]
        .to_vec()
}
/// Compute the mean (average) of a slice.
#[allow(dead_code)]
pub fn slice_mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    data.iter().sum::<f64>() / data.len() as f64
}
/// Compute the variance of a slice.
#[allow(dead_code)]
pub fn slice_variance(data: &[f64]) -> f64 {
    if data.len() < 2 {
        return 0.0;
    }
    let m = slice_mean(data);
    data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / data.len() as f64
}
/// Compute the autocorrelation at lag `lag` for a centred signal.
#[allow(dead_code)]
pub fn autocorrelation(data: &[f64], lag: usize) -> f64 {
    let n = data.len();
    if lag >= n {
        return 0.0;
    }
    let m = slice_mean(data);
    let var = slice_variance(data);
    if var == 0.0 {
        return 0.0;
    }
    let cov: f64 = (0..n - lag)
        .map(|i| (data[i] - m) * (data[i + lag] - m))
        .sum::<f64>()
        / (n - lag) as f64;
    cov / var
}
