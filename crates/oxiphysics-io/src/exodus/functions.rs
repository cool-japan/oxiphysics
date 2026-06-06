//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ElementQuality, ExodusMesh, ExodusResult, ExodusSideSet, ExodusTimeStep, ExodusVariable,
    FieldStats, MeshPartition, MeshValidationReport,
};

/// Write an `ExodusMesh` to a simplified text representation.
///
/// Format:
/// ```text
/// EXODUS_MESH
/// TITLE `title`
/// NODES `count`
/// `x` `y` `z`
/// ...
/// BLOCKS `count`
/// BLOCK `id` `name` `etype` <n_nodes_per_elem> `n_elements`
/// `conn0` `conn1` ...
/// ...
/// NODE_SETS `count`
/// NODE_SET `id` `name` `count`
/// `nid0` `nid1` ...
/// SIDE_SETS `count`
/// SIDE_SET `id` `name` `count`
/// `eid0` `sid0` `eid1` `sid1` ...
/// END
/// ```
pub fn write_text(mesh: &ExodusMesh) -> String {
    let mut s = String::new();
    s.push_str("EXODUS_MESH\n");
    s.push_str(&format!("TITLE {}\n", mesh.title));
    s.push_str(&format!("NODES {}\n", mesh.coordinates.len()));
    for coord in &mesh.coordinates {
        s.push_str(&format!("{} {} {}\n", coord[0], coord[1], coord[2]));
    }
    s.push_str(&format!("BLOCKS {}\n", mesh.blocks.len()));
    for block in &mesh.blocks {
        let n_elems = block
            .connectivity
            .len()
            .checked_div(block.n_nodes_per_elem)
            .unwrap_or(0);
        s.push_str(&format!(
            "BLOCK {} {} {} {} {}\n",
            block.id, block.name, block.element_type, block.n_nodes_per_elem, n_elems
        ));
        if !block.connectivity.is_empty() {
            let strs: Vec<String> = block.connectivity.iter().map(|v| v.to_string()).collect();
            s.push_str(&strs.join(" "));
            s.push('\n');
        }
    }
    s.push_str(&format!("NODE_SETS {}\n", mesh.node_sets.len()));
    for ns in &mesh.node_sets {
        s.push_str(&format!(
            "NODE_SET {} {} {}\n",
            ns.id,
            ns.name,
            ns.node_ids.len()
        ));
        if !ns.node_ids.is_empty() {
            let strs: Vec<String> = ns.node_ids.iter().map(|v| v.to_string()).collect();
            s.push_str(&strs.join(" "));
            s.push('\n');
        }
    }
    s.push_str(&format!("SIDE_SETS {}\n", mesh.side_sets.len()));
    for ss in &mesh.side_sets {
        s.push_str(&format!(
            "SIDE_SET {} {} {}\n",
            ss.id,
            ss.name,
            ss.element_ids.len()
        ));
        if !ss.element_ids.is_empty() {
            let mut pairs = Vec::new();
            for (eid, sid) in ss.element_ids.iter().zip(ss.side_ids.iter()) {
                pairs.push(format!("{} {}", eid, sid));
            }
            s.push_str(&pairs.join(" "));
            s.push('\n');
        }
    }
    s.push_str("END\n");
    s
}
/// Parse a text-format Exodus mesh string back into an `ExodusMesh`.
pub fn parse_text(s: &str) -> Result<ExodusMesh, String> {
    let mut lines = s.lines().peekable();
    match lines.next() {
        Some("EXODUS_MESH") => {}
        other => return Err(format!("Expected EXODUS_MESH header, got: {:?}", other)),
    }
    let title = match lines.next() {
        Some(line) => match line.strip_prefix("TITLE ") {
            Some(rest) => rest.trim().to_string(),
            None => return Err(format!("Expected TITLE, got: {:?}", line)),
        },
        None => return Err("Expected TITLE, got: None".to_string()),
    };
    let mut mesh = ExodusMesh::new(&title);
    let n_nodes: usize = match lines.next() {
        Some(line) => match line.strip_prefix("NODES ") {
            Some(rest) => rest
                .trim()
                .parse()
                .map_err(|e| format!("Bad node count: {}", e))?,
            None => return Err(format!("Expected NODES, got: {:?}", line)),
        },
        None => return Err("Expected NODES, got: None".to_string()),
    };
    for _ in 0..n_nodes {
        let line = lines.next().ok_or("Unexpected EOF reading nodes")?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(format!("Expected 3 coordinates, got: {}", line));
        }
        let x: f64 = parts[0].parse().map_err(|e| format!("Bad x: {}", e))?;
        let y: f64 = parts[1].parse().map_err(|e| format!("Bad y: {}", e))?;
        let z: f64 = parts[2].parse().map_err(|e| format!("Bad z: {}", e))?;
        mesh.add_node(x, y, z);
    }
    let n_blocks: usize = match lines.next() {
        Some(line) => match line.strip_prefix("BLOCKS ") {
            Some(rest) => rest
                .trim()
                .parse()
                .map_err(|e| format!("Bad block count: {}", e))?,
            None => return Err(format!("Expected BLOCKS, got: {:?}", line)),
        },
        None => return Err("Expected BLOCKS, got: None".to_string()),
    };
    for _ in 0..n_blocks {
        let header = lines.next().ok_or("Unexpected EOF reading block header")?;
        let parts: Vec<&str> = header.split_whitespace().collect();
        if parts.len() < 6 || parts[0] != "BLOCK" {
            return Err(format!("Expected BLOCK header, got: {}", header));
        }
        let id: u32 = parts[1]
            .parse()
            .map_err(|e| format!("Bad block id: {}", e))?;
        let name = parts[2].to_string();
        let etype = parts[3].to_string();
        let n_nodes_per_elem: usize = parts[4]
            .parse()
            .map_err(|e| format!("Bad n_nodes_per_elem: {}", e))?;
        let n_elems: usize = parts[5]
            .parse()
            .map_err(|e| format!("Bad n_elems: {}", e))?;
        let block_idx = mesh.add_block(id, &name, &etype, n_nodes_per_elem);
        if n_elems > 0 {
            let conn_line = lines.next().ok_or("Unexpected EOF reading connectivity")?;
            let conn: Result<Vec<usize>, _> = conn_line
                .split_whitespace()
                .map(|v| v.parse::<usize>())
                .collect();
            let conn = conn.map_err(|e| format!("Bad connectivity: {}", e))?;
            mesh.add_element_to_block(block_idx, &conn);
        }
    }
    let n_node_sets: usize = match lines.next() {
        Some(line) => match line.strip_prefix("NODE_SETS ") {
            Some(rest) => rest
                .trim()
                .parse()
                .map_err(|e| format!("Bad node set count: {}", e))?,
            None => return Err(format!("Expected NODE_SETS, got: {:?}", line)),
        },
        None => return Err("Expected NODE_SETS, got: None".to_string()),
    };
    for _ in 0..n_node_sets {
        let header = lines
            .next()
            .ok_or("Unexpected EOF reading node set header")?;
        let parts: Vec<&str> = header.split_whitespace().collect();
        if parts.len() < 4 || parts[0] != "NODE_SET" {
            return Err(format!("Expected NODE_SET header, got: {}", header));
        }
        let id: u32 = parts[1]
            .parse()
            .map_err(|e| format!("Bad node set id: {}", e))?;
        let name = parts[2].to_string();
        let count: usize = parts[3]
            .parse()
            .map_err(|e| format!("Bad node set count: {}", e))?;
        let nodes = if count > 0 {
            let data_line = lines.next().ok_or("Unexpected EOF reading node set data")?;
            let parsed: Result<Vec<usize>, _> = data_line
                .split_whitespace()
                .map(|v| v.parse::<usize>())
                .collect();
            parsed.map_err(|e| format!("Bad node set data: {}", e))?
        } else {
            Vec::new()
        };
        mesh.add_node_set(id, &name, nodes);
    }
    let n_side_sets: usize = match lines.next() {
        Some(line) => match line.strip_prefix("SIDE_SETS ") {
            Some(rest) => rest
                .trim()
                .parse()
                .map_err(|e| format!("Bad side set count: {}", e))?,
            None => return Err(format!("Expected SIDE_SETS, got: {:?}", line)),
        },
        None => return Err("Expected SIDE_SETS, got: None".to_string()),
    };
    for _ in 0..n_side_sets {
        let header = lines
            .next()
            .ok_or("Unexpected EOF reading side set header")?;
        let parts: Vec<&str> = header.split_whitespace().collect();
        if parts.len() < 4 || parts[0] != "SIDE_SET" {
            return Err(format!("Expected SIDE_SET header, got: {}", header));
        }
        let id: u32 = parts[1]
            .parse()
            .map_err(|e| format!("Bad side set id: {}", e))?;
        let name = parts[2].to_string();
        let count: usize = parts[3]
            .parse()
            .map_err(|e| format!("Bad side set count: {}", e))?;
        let (element_ids, side_ids) = if count > 0 {
            let data_line = lines.next().ok_or("Unexpected EOF reading side set data")?;
            let nums: Result<Vec<usize>, _> = data_line
                .split_whitespace()
                .map(|v| v.parse::<usize>())
                .collect();
            let nums = nums.map_err(|e| format!("Bad side set data: {}", e))?;
            let mut eids = Vec::with_capacity(count);
            let mut sids = Vec::with_capacity(count);
            for chunk in nums.chunks(2) {
                if chunk.len() == 2 {
                    eids.push(chunk[0]);
                    sids.push(chunk[1]);
                }
            }
            (eids, sids)
        } else {
            (Vec::new(), Vec::new())
        };
        mesh.side_sets.push(ExodusSideSet {
            id,
            name,
            element_ids,
            side_ids,
        });
    }
    match lines.next() {
        Some("END") => {}
        other => return Err(format!("Expected END, got: {:?}", other)),
    }
    Ok(mesh)
}
/// Compute the axis-aligned bounding box of a mesh.
///
/// Returns `(min, max)` corner points.
pub fn bounding_box(mesh: &ExodusMesh) -> ([f64; 3], [f64; 3]) {
    if mesh.coordinates.is_empty() {
        return ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
    }
    let first = mesh.coordinates[0];
    let mut min = first;
    let mut max = first;
    for &coord in &mesh.coordinates {
        for i in 0..3 {
            if coord[i] < min[i] {
                min[i] = coord[i];
            }
            if coord[i] > max[i] {
                max[i] = coord[i];
            }
        }
    }
    (min, max)
}
/// Translate all nodes in the mesh by a constant offset.
pub fn translate_mesh(mesh: &mut ExodusMesh, offset: [f64; 3]) {
    for coord in &mut mesh.coordinates {
        coord[0] += offset[0];
        coord[1] += offset[1];
        coord[2] += offset[2];
    }
}
/// Compute the time values for `n_steps` linearly spaced from `t_start` to `t_end`.
pub fn linear_time_steps(t_start: f64, t_end: f64, n_steps: usize) -> Vec<f64> {
    if n_steps == 0 {
        return vec![];
    }
    if n_steps == 1 {
        return vec![t_start];
    }
    let dt = (t_end - t_start) / (n_steps - 1) as f64;
    (0..n_steps).map(|i| t_start + i as f64 * dt).collect()
}
/// Compute the time values for `n_steps` logarithmically spaced from `t_start` to `t_end`.
///
/// Both `t_start` and `t_end` must be positive.
pub fn log_time_steps(t_start: f64, t_end: f64, n_steps: usize) -> Vec<f64> {
    if n_steps == 0 || t_start <= 0.0 || t_end <= 0.0 {
        return vec![];
    }
    if n_steps == 1 {
        return vec![t_start];
    }
    let log_start = t_start.ln();
    let log_end = t_end.ln();
    let d = (log_end - log_start) / (n_steps - 1) as f64;
    (0..n_steps)
        .map(|i| (log_start + i as f64 * d).exp())
        .collect()
}
/// Add a time step with a scalar node variable to an `ExodusResult`.
pub fn add_node_variable_step(
    result: &mut ExodusResult,
    time: f64,
    var_name: &str,
    values: Vec<f64>,
) {
    result.time_steps.push(ExodusTimeStep {
        time,
        variables: vec![ExodusVariable {
            name: var_name.to_string(),
            entity_type: "node".to_string(),
            values,
        }],
    });
}
/// Add a time step with an element variable to an `ExodusResult`.
pub fn add_element_variable_step(
    result: &mut ExodusResult,
    time: f64,
    var_name: &str,
    values: Vec<f64>,
) {
    result.time_steps.push(ExodusTimeStep {
        time,
        variables: vec![ExodusVariable {
            name: var_name.to_string(),
            entity_type: "element".to_string(),
            values,
        }],
    });
}
/// Validate an `ExodusMesh` for common consistency issues.
///
/// Checks:
/// - Node connectivity: all node IDs in blocks are valid (1-based, ≤ n_nodes).
/// - Block completeness: connectivity length is divisible by `n_nodes_per_elem`.
/// - Node set node IDs are in valid range.
/// - Side set element IDs are in valid range.
/// - No duplicate block IDs.
/// - No duplicate node-set IDs.
pub fn validate_mesh(mesh: &ExodusMesh) -> MeshValidationReport {
    let mut report = MeshValidationReport {
        is_valid: true,
        errors: Vec::new(),
        warnings: Vec::new(),
    };
    let n_nodes = mesh.node_count();
    let total_elements = mesh.element_count();
    if n_nodes == 0 {
        report.add_warning("Mesh has no nodes".to_string());
    }
    if total_elements == 0 {
        report.add_warning("Mesh has no elements".to_string());
    }
    let mut seen_block_ids: Vec<u32> = Vec::new();
    for (bi, block) in mesh.blocks.iter().enumerate() {
        if seen_block_ids.contains(&block.id) {
            report.add_error(format!(
                "Duplicate block id {} at block index {}",
                block.id, bi
            ));
        } else {
            seen_block_ids.push(block.id);
        }
        if block.n_nodes_per_elem == 0 {
            report.add_error(format!("Block {} has zero nodes_per_elem", block.id));
            continue;
        }
        if block.connectivity.len() % block.n_nodes_per_elem != 0 {
            report.add_error(format!(
                "Block {}: connectivity len {} not divisible by n_nodes_per_elem {}",
                block.id,
                block.connectivity.len(),
                block.n_nodes_per_elem
            ));
        }
        for &nid in &block.connectivity {
            if nid == 0 || nid > n_nodes {
                report.add_error(format!(
                    "Block {}: node id {} out of range [1, {}]",
                    block.id, nid, n_nodes
                ));
            }
        }
    }
    let mut seen_ns_ids: Vec<u32> = Vec::new();
    for ns in &mesh.node_sets {
        if seen_ns_ids.contains(&ns.id) {
            report.add_error(format!("Duplicate node set id {}", ns.id));
        } else {
            seen_ns_ids.push(ns.id);
        }
        for &nid in &ns.node_ids {
            if nid == 0 || nid > n_nodes {
                report.add_error(format!(
                    "Node set {}: node id {} out of range [1, {}]",
                    ns.id, nid, n_nodes
                ));
            }
        }
    }
    for ss in &mesh.side_sets {
        if ss.element_ids.len() != ss.side_ids.len() {
            report.add_error(format!(
                "Side set {}: element_ids and side_ids have different lengths ({} vs {})",
                ss.id,
                ss.element_ids.len(),
                ss.side_ids.len()
            ));
        }
        for &eid in &ss.element_ids {
            if eid == 0 || eid > total_elements {
                report.add_error(format!(
                    "Side set {}: element id {} out of range [1, {}]",
                    ss.id, eid, total_elements
                ));
            }
        }
    }
    report
}
/// Scale all node coordinates by a uniform factor.
pub fn scale_mesh(mesh: &mut ExodusMesh, factor: f64) {
    for coord in &mut mesh.coordinates {
        coord[0] *= factor;
        coord[1] *= factor;
        coord[2] *= factor;
    }
}
/// Compute the centroid of a mesh (arithmetic mean of node coordinates).
pub fn mesh_centroid(mesh: &ExodusMesh) -> [f64; 3] {
    let n = mesh.coordinates.len();
    if n == 0 {
        return [0.0, 0.0, 0.0];
    }
    let mut sum = [0.0_f64; 3];
    for &c in &mesh.coordinates {
        sum[0] += c[0];
        sum[1] += c[1];
        sum[2] += c[2];
    }
    [sum[0] / n as f64, sum[1] / n as f64, sum[2] / n as f64]
}
/// Build a unit-cube HEX8 mesh with 8 nodes and 1 element.
pub fn build_unit_hex_mesh(title: &str) -> ExodusMesh {
    let mut mesh = ExodusMesh::new(title);
    let coords = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ];
    for (x, y, z) in coords {
        mesh.add_node(x, y, z);
    }
    let bi = mesh.add_block(1, "hex_block", "HEX8", 8);
    mesh.add_element_to_block(bi, &[1, 2, 3, 4, 5, 6, 7, 8]);
    mesh
}
/// Build a unit tetrahedron mesh with 4 nodes and 1 TET4 element.
pub fn build_unit_tet_mesh(title: &str) -> ExodusMesh {
    let mut mesh = ExodusMesh::new(title);
    mesh.add_node(0.0, 0.0, 0.0);
    mesh.add_node(1.0, 0.0, 0.0);
    mesh.add_node(0.0, 1.0, 0.0);
    mesh.add_node(0.0, 0.0, 1.0);
    let bi = mesh.add_block(1, "tet_block", "TET4", 4);
    mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
    mesh
}
/// Merge two meshes together (renumber the second mesh's node and element IDs).
///
/// Returns a new mesh combining both. Node IDs in the second mesh are offset
/// by the number of nodes in the first mesh.
pub fn merge_meshes(a: &ExodusMesh, b: &ExodusMesh) -> ExodusMesh {
    let mut out = ExodusMesh::new(&format!("merged_{}_{}", a.title, b.title));
    for &c in &a.coordinates {
        out.coordinates.push(c);
    }
    let offset = a.node_count();
    for &c in &b.coordinates {
        out.coordinates.push(c);
    }
    for block in &a.blocks {
        let bi = out.add_block(
            block.id,
            &block.name,
            &block.element_type,
            block.n_nodes_per_elem,
        );
        out.add_element_to_block(bi, &block.connectivity);
    }
    let max_bid = a.blocks.iter().map(|b| b.id).max().unwrap_or(0);
    for block in &b.blocks {
        let new_id = block.id + max_bid;
        let bi = out.add_block(
            new_id,
            &block.name,
            &block.element_type,
            block.n_nodes_per_elem,
        );
        let shifted: Vec<usize> = block.connectivity.iter().map(|&n| n + offset).collect();
        out.add_element_to_block(bi, &shifted);
    }
    for ns in &a.node_sets {
        out.add_node_set(ns.id, &ns.name, ns.node_ids.clone());
    }
    for ns in &b.node_sets {
        let shifted: Vec<usize> = ns.node_ids.iter().map(|&n| n + offset).collect();
        out.add_node_set(ns.id + a.node_sets.len() as u32, &ns.name, shifted);
    }
    out
}
/// Compute the Euclidean distance between two nodes.
pub fn node_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Compute quality metrics for all TET4 elements in a mesh.
///
/// Elements that do not belong to TET4 blocks are skipped.
pub fn tet4_quality(mesh: &ExodusMesh) -> Vec<ElementQuality> {
    let mut result = Vec::new();
    for block in &mesh.blocks {
        if block.element_type != "TET4" || block.n_nodes_per_elem != 4 {
            continue;
        }
        let n_elem = block.connectivity.len() / 4;
        for ei in 0..n_elem {
            let base = ei * 4;
            let n0 = mesh.coordinates[block.connectivity[base] - 1];
            let n1 = mesh.coordinates[block.connectivity[base + 1] - 1];
            let n2 = mesh.coordinates[block.connectivity[base + 2] - 1];
            let n3 = mesh.coordinates[block.connectivity[base + 3] - 1];
            let edges = [
                node_distance(n0, n1),
                node_distance(n0, n2),
                node_distance(n0, n3),
                node_distance(n1, n2),
                node_distance(n1, n3),
                node_distance(n2, n3),
            ];
            let min_e = edges.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_e = edges.iter().cloned().fold(0.0_f64, f64::max);
            let aspect = if max_e > 1e-30 { min_e / max_e } else { 0.0 };
            let e1 = [n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]];
            let e2 = [n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]];
            let e3 = [n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]];
            let jdet = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
                - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
                + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
            let len1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
            let len2 = (e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]).sqrt();
            let len3 = (e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]).sqrt();
            let denom = len1 * len2 * len3;
            let scaled_jac = if denom > 1e-30 { jdet / denom } else { 0.0 };
            result.push(ElementQuality {
                block_id: block.id,
                element_index: ei,
                aspect_ratio: aspect,
                scaled_jacobian: scaled_jac,
                min_edge: min_e,
                max_edge: max_e,
            });
        }
    }
    result
}
/// Compute quality metrics for all HEX8 elements in a mesh.
pub fn hex8_quality(mesh: &ExodusMesh) -> Vec<ElementQuality> {
    let mut result = Vec::new();
    for block in &mesh.blocks {
        if block.element_type != "HEX8" || block.n_nodes_per_elem != 8 {
            continue;
        }
        let n_elem = block.connectivity.len() / 8;
        for ei in 0..n_elem {
            let base = ei * 8;
            let nodes: Vec<[f64; 3]> = (0..8)
                .map(|k| mesh.coordinates[block.connectivity[base + k] - 1])
                .collect();
            let edge_pairs = [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 0),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 4),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ];
            let edges: Vec<f64> = edge_pairs
                .iter()
                .map(|&(a, b)| node_distance(nodes[a], nodes[b]))
                .collect();
            let min_e = edges.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_e = edges.iter().cloned().fold(0.0_f64, f64::max);
            let aspect = if max_e > 1e-30 { min_e / max_e } else { 0.0 };
            let e1 = [
                nodes[1][0] - nodes[0][0],
                nodes[1][1] - nodes[0][1],
                nodes[1][2] - nodes[0][2],
            ];
            let e2 = [
                nodes[3][0] - nodes[0][0],
                nodes[3][1] - nodes[0][1],
                nodes[3][2] - nodes[0][2],
            ];
            let e3 = [
                nodes[4][0] - nodes[0][0],
                nodes[4][1] - nodes[0][1],
                nodes[4][2] - nodes[0][2],
            ];
            let jdet = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
                - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
                + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
            let len1 = (e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]).sqrt();
            let len2 = (e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]).sqrt();
            let len3 = (e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]).sqrt();
            let denom = len1 * len2 * len3;
            let scaled_jac = if denom > 1e-30 { jdet / denom } else { 0.0 };
            result.push(ElementQuality {
                block_id: block.id,
                element_index: ei,
                aspect_ratio: aspect,
                scaled_jacobian: scaled_jac,
                min_edge: min_e,
                max_edge: max_e,
            });
        }
    }
    result
}
/// Interpolate a scalar field at an arbitrary point using nearest-node
/// inverse-distance weighting.
///
/// `field` must have one value per node.
pub fn idw_interpolate(mesh: &ExodusMesh, field: &[f64], query: [f64; 3], power: f64) -> f64 {
    if field.is_empty() || mesh.coordinates.is_empty() {
        return 0.0;
    }
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (i, &c) in mesh.coordinates.iter().enumerate() {
        if i >= field.len() {
            break;
        }
        let d = node_distance(c, query);
        if d < 1e-15 {
            return field[i];
        }
        let w = 1.0 / d.powf(power);
        num += w * field[i];
        den += w;
    }
    if den.abs() < 1e-30 { 0.0 } else { num / den }
}
/// Find the index of the node closest to a query point.
pub fn find_nearest_node(mesh: &ExodusMesh, query: [f64; 3]) -> Option<usize> {
    if mesh.coordinates.is_empty() {
        return None;
    }
    let mut best_idx = 0;
    let mut best_dist = f64::INFINITY;
    for (i, &c) in mesh.coordinates.iter().enumerate() {
        let d = node_distance(c, query);
        if d < best_dist {
            best_dist = d;
            best_idx = i;
        }
    }
    Some(best_idx)
}
/// Extract the subset of nodes within a spherical region.
///
/// Returns 0-based indices of nodes whose distance from `center` is ≤ `radius`.
pub fn nodes_in_sphere(mesh: &ExodusMesh, center: [f64; 3], radius: f64) -> Vec<usize> {
    mesh.coordinates
        .iter()
        .enumerate()
        .filter(|&(_, c)| node_distance(*c, center) <= radius)
        .map(|(i, _)| i)
        .collect()
}
/// Extract a sub-mesh containing only elements in the specified block.
pub fn extract_block(mesh: &ExodusMesh, block_id: u32) -> Option<ExodusMesh> {
    let block = mesh.blocks.iter().find(|b| b.id == block_id)?;
    let mut out = ExodusMesh::new(&format!("{}_block_{}", mesh.title, block_id));
    let mut used_nodes: Vec<usize> = block.connectivity.clone();
    used_nodes.sort_unstable();
    used_nodes.dedup();
    let mut remap = std::collections::HashMap::new();
    for (new_idx, &old_id) in used_nodes.iter().enumerate() {
        remap.insert(old_id, new_idx + 1);
    }
    for &old_id in &used_nodes {
        if old_id <= mesh.coordinates.len() {
            let c = mesh.coordinates[old_id - 1];
            out.coordinates.push(c);
        }
    }
    let new_conn: Vec<usize> = block
        .connectivity
        .iter()
        .map(|&n| *remap.get(&n).unwrap_or(&n))
        .collect();
    let bi = out.add_block(
        block.id,
        &block.name,
        &block.element_type,
        block.n_nodes_per_elem,
    );
    out.add_element_to_block(bi, &new_conn);
    Some(out)
}
/// Extract a sub-mesh containing only the nodes in a node set.
pub fn extract_node_set_coords(mesh: &ExodusMesh, ns_id: u32) -> Option<Vec<[f64; 3]>> {
    let ns = mesh.node_sets.iter().find(|ns| ns.id == ns_id)?;
    let coords: Vec<[f64; 3]> = ns
        .node_ids
        .iter()
        .filter_map(|&nid| {
            if nid >= 1 && nid <= mesh.coordinates.len() {
                Some(mesh.coordinates[nid - 1])
            } else {
                None
            }
        })
        .collect();
    Some(coords)
}
/// Refine a mesh by midpoint subdivision of TET4 elements.
///
/// Each TET4 is split into 8 smaller tetrahedra by inserting 6 midpoint nodes
/// on every edge (regular octasection). This produces a conforming mesh.
pub fn refine_tet4_mesh(mesh: &ExodusMesh) -> ExodusMesh {
    let mut out = ExodusMesh::new(&format!("{}_refined", mesh.title));
    for &c in &mesh.coordinates {
        out.coordinates.push(c);
    }
    let mut edge_mid: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    let add_midpoint = |coords: &mut Vec<[f64; 3]>, a: usize, b: usize| -> usize {
        let ca = coords[a - 1];
        let cb = coords[b - 1];
        let mid = [
            (ca[0] + cb[0]) * 0.5,
            (ca[1] + cb[1]) * 0.5,
            (ca[2] + cb[2]) * 0.5,
        ];
        coords.push(mid);
        coords.len()
    };
    for block in &mesh.blocks {
        if block.element_type != "TET4" || block.n_nodes_per_elem != 4 {
            let bi = out.add_block(
                block.id,
                &block.name,
                &block.element_type,
                block.n_nodes_per_elem,
            );
            out.add_element_to_block(bi, &block.connectivity);
            continue;
        }
        let bi = out.add_block(block.id, &block.name, "TET4", 4);
        let n_elem = block.connectivity.len() / 4;
        for ei in 0..n_elem {
            let base = ei * 4;
            let (v0, v1, v2, v3) = (
                block.connectivity[base],
                block.connectivity[base + 1],
                block.connectivity[base + 2],
                block.connectivity[base + 3],
            );
            let get_mid = |a: usize,
                           b: usize,
                           coords: &mut Vec<[f64; 3]>,
                           cache: &mut std::collections::HashMap<(usize, usize), usize>|
             -> usize {
                let key = if a < b { (a, b) } else { (b, a) };
                if let Some(&mid) = cache.get(&key) {
                    mid
                } else {
                    let mid = add_midpoint(coords, a, b);
                    cache.insert(key, mid);
                    mid
                }
            };
            let m01 = get_mid(v0, v1, &mut out.coordinates, &mut edge_mid);
            let m02 = get_mid(v0, v2, &mut out.coordinates, &mut edge_mid);
            let m03 = get_mid(v0, v3, &mut out.coordinates, &mut edge_mid);
            let m12 = get_mid(v1, v2, &mut out.coordinates, &mut edge_mid);
            let m13 = get_mid(v1, v3, &mut out.coordinates, &mut edge_mid);
            let m23 = get_mid(v2, v3, &mut out.coordinates, &mut edge_mid);
            let sub_tets = [
                [v0, m01, m02, m03],
                [v1, m01, m12, m13],
                [v2, m02, m12, m23],
                [v3, m03, m13, m23],
                [m01, m02, m03, m13],
                [m01, m02, m12, m13],
                [m02, m03, m13, m23],
                [m02, m12, m13, m23],
            ];
            for tet in &sub_tets {
                out.add_element_to_block(bi, tet);
            }
        }
    }
    out
}
/// Extract all unique node IDs referenced in node sets as a sorted vec.
pub fn all_boundary_nodes(mesh: &ExodusMesh) -> Vec<usize> {
    let mut ids: Vec<usize> = mesh
        .node_sets
        .iter()
        .flat_map(|ns| ns.node_ids.iter().cloned())
        .collect();
    ids.sort_unstable();
    ids.dedup();
    ids
}
/// Count how many times each node appears in element connectivity.
///
/// Returns a vector of length `n_nodes` with valence counts (0-based index).
pub fn node_valence(mesh: &ExodusMesh) -> Vec<usize> {
    let n = mesh.node_count();
    let mut valence = vec![0usize; n];
    for block in &mesh.blocks {
        for &nid in &block.connectivity {
            if nid >= 1 && nid <= n {
                valence[nid - 1] += 1;
            }
        }
    }
    valence
}
/// Write an `ExodusResult` (mesh + time steps) to a text representation.
///
/// Variables are written after the mesh section with their time stamps.
pub fn write_result_text(result: &ExodusResult) -> String {
    let mut s = write_text(&result.mesh);
    if s.ends_with("END\n") {
        s.truncate(s.len() - 4);
    }
    s.push_str(&format!("TIME_STEPS {}\n", result.time_steps.len()));
    for ts in &result.time_steps {
        s.push_str(&format!("TIME_STEP {:.15e}\n", ts.time));
        s.push_str(&format!("VARIABLES {}\n", ts.variables.len()));
        for var in &ts.variables {
            s.push_str(&format!(
                "VAR {} {} {}\n",
                var.entity_type,
                var.name,
                var.values.len()
            ));
            let vals: Vec<String> = var.values.iter().map(|v| format!("{:.15e}", v)).collect();
            s.push_str(&vals.join(" "));
            s.push('\n');
        }
    }
    s.push_str("END\n");
    s
}
/// Parse an `ExodusResult` from the text format written by `write_result_text`.
pub fn parse_result_text(data: &str) -> std::result::Result<ExodusResult, String> {
    let ts_marker = "TIME_STEPS ";
    if let Some(pos) = data.find(ts_marker) {
        let mesh_part = format!("{}END\n", &data[..pos]);
        let mesh = parse_text(&mesh_part)?;
        let mut result = ExodusResult {
            mesh,
            time_steps: Vec::new(),
        };
        let rest = &data[pos..];
        let mut lines = rest.lines().peekable();
        let n_steps: usize = match lines.next() {
            Some(line) => match line.strip_prefix("TIME_STEPS ") {
                Some(rest) => rest
                    .trim()
                    .parse()
                    .map_err(|e| format!("Bad TIME_STEPS: {}", e))?,
                None => return Err(format!("Expected TIME_STEPS, got: {:?}", line)),
            },
            None => return Err("Expected TIME_STEPS, got: None".to_string()),
        };
        for _ in 0..n_steps {
            let time: f64 = match lines.next() {
                Some(line) => match line.strip_prefix("TIME_STEP ") {
                    Some(rest) => rest
                        .trim()
                        .parse()
                        .map_err(|e| format!("Bad TIME_STEP: {}", e))?,
                    None => return Err(format!("Expected TIME_STEP, got: {:?}", line)),
                },
                None => return Err("Expected TIME_STEP, got: None".to_string()),
            };
            let n_vars: usize = match lines.next() {
                Some(line) => match line.strip_prefix("VARIABLES ") {
                    Some(rest) => rest
                        .trim()
                        .parse()
                        .map_err(|e| format!("Bad VARIABLES: {}", e))?,
                    None => return Err(format!("Expected VARIABLES, got: {:?}", line)),
                },
                None => return Err("Expected VARIABLES, got: None".to_string()),
            };
            let mut variables = Vec::new();
            for _ in 0..n_vars {
                let hdr = lines.next().ok_or("Unexpected EOF reading VAR")?;
                let parts: Vec<&str> = hdr.split_whitespace().collect();
                if parts.len() < 4 || parts[0] != "VAR" {
                    return Err(format!("Expected VAR header, got: {}", hdr));
                }
                let entity_type = parts[1].to_string();
                let name = parts[2].to_string();
                let n_vals: usize = parts[3]
                    .parse()
                    .map_err(|e| format!("Bad VAR count: {}", e))?;
                let vals_line = lines.next().ok_or("Unexpected EOF reading VAR data")?;
                let values: std::result::Result<Vec<f64>, _> = vals_line
                    .split_whitespace()
                    .take(n_vals)
                    .map(|v| v.parse::<f64>())
                    .collect();
                let values = values.map_err(|e| format!("Bad VAR values: {}", e))?;
                variables.push(ExodusVariable {
                    name,
                    entity_type,
                    values,
                });
            }
            result.time_steps.push(ExodusTimeStep { time, variables });
        }
        Ok(result)
    } else {
        let mesh = parse_text(data)?;
        Ok(ExodusResult {
            mesh,
            time_steps: Vec::new(),
        })
    }
}
/// Partition a mesh into `n_parts` using a round-robin element assignment.
///
/// Only the first block is partitioned; subsequent blocks are copied to part 0.
pub fn partition_mesh_round_robin(mesh: &ExodusMesh, n_parts: usize) -> Vec<MeshPartition> {
    if n_parts == 0 || mesh.blocks.is_empty() {
        return Vec::new();
    }
    let mut parts: Vec<ExodusMesh> = (0..n_parts)
        .map(|i| ExodusMesh::new(&format!("{}_part_{}", mesh.title, i)))
        .collect();
    for p in &mut parts {
        for &c in &mesh.coordinates {
            p.coordinates.push(c);
        }
    }
    let block = &mesh.blocks[0];
    let n_elem = block
        .connectivity
        .len()
        .checked_div(block.n_nodes_per_elem)
        .unwrap_or(0);
    let block_idxs: Vec<usize> = (0..n_parts)
        .map(|pi| {
            parts[pi].add_block(
                block.id,
                &block.name,
                &block.element_type,
                block.n_nodes_per_elem,
            )
        })
        .collect();
    for ei in 0..n_elem {
        let part = ei % n_parts;
        let base = ei * block.n_nodes_per_elem;
        let conn = &block.connectivity[base..base + block.n_nodes_per_elem];
        parts[part].add_element_to_block(block_idxs[part], conn);
    }
    for ns in &mesh.node_sets {
        for p in &mut parts {
            p.add_node_set(ns.id, &ns.name, ns.node_ids.clone());
        }
    }
    parts
        .into_iter()
        .enumerate()
        .map(|(i, m)| MeshPartition {
            partition_id: i,
            mesh: m,
        })
        .collect()
}
/// Compute statistics for a scalar field.
pub fn field_stats(values: &[f64]) -> FieldStats {
    if values.is_empty() {
        return FieldStats {
            min: 0.0,
            max: 0.0,
            mean: 0.0,
            std_dev: 0.0,
            count: 0,
        };
    }
    let count = values.len();
    let mut min_v = values[0];
    let mut max_v = values[0];
    let mut sum = 0.0_f64;
    for &v in values {
        if v < min_v {
            min_v = v;
        }
        if v > max_v {
            max_v = v;
        }
        sum += v;
    }
    let mean = sum / count as f64;
    let var: f64 = values
        .iter()
        .map(|&v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / count as f64;
    FieldStats {
        min: min_v,
        max: max_v,
        mean,
        std_dev: var.sqrt(),
        count,
    }
}
/// Normalize a scalar field to the range \[0, 1\].
pub fn normalize_field(values: &[f64]) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let stats = field_stats(values);
    let range = stats.max - stats.min;
    if range < 1e-30 {
        return vec![0.5; values.len()];
    }
    values.iter().map(|&v| (v - stats.min) / range).collect()
}
/// Compute the L2 norm of a field.
pub fn field_l2_norm(values: &[f64]) -> f64 {
    let sum_sq: f64 = values.iter().map(|&v| v * v).sum();
    sum_sq.sqrt()
}
/// Export node coordinates to a CSV string.
pub fn export_nodes_csv(mesh: &ExodusMesh) -> String {
    let mut s = String::from("node_id,x,y,z\n");
    for (i, &c) in mesh.coordinates.iter().enumerate() {
        s.push_str(&format!("{},{},{},{}\n", i + 1, c[0], c[1], c[2]));
    }
    s
}
/// Export node field values to a CSV string alongside node coordinates.
pub fn export_field_csv(mesh: &ExodusMesh, field_name: &str, values: &[f64]) -> String {
    let mut s = format!("node_id,x,y,z,{}\n", field_name);
    for (i, &c) in mesh.coordinates.iter().enumerate() {
        let v = values.get(i).copied().unwrap_or(0.0);
        s.push_str(&format!("{},{},{},{},{}\n", i + 1, c[0], c[1], c[2], v));
    }
    s
}
/// Export element connectivity to a CSV string.
pub fn export_connectivity_csv(mesh: &ExodusMesh) -> String {
    let mut s = String::from("block_id,element_id");
    let max_npe = mesh
        .blocks
        .iter()
        .map(|b| b.n_nodes_per_elem)
        .max()
        .unwrap_or(0);
    for k in 0..max_npe {
        s.push_str(&format!(",n{}", k));
    }
    s.push('\n');
    let mut global_elem = 1;
    for block in &mesh.blocks {
        if block.n_nodes_per_elem == 0 {
            continue;
        }
        let n_elem = block.connectivity.len() / block.n_nodes_per_elem;
        for ei in 0..n_elem {
            let base = ei * block.n_nodes_per_elem;
            s.push_str(&format!("{},{}", block.id, global_elem));
            for k in 0..block.n_nodes_per_elem {
                s.push_str(&format!(",{}", block.connectivity[base + k]));
            }
            for _ in block.n_nodes_per_elem..max_npe {
                s.push_str(",-1");
            }
            s.push('\n');
            global_elem += 1;
        }
    }
    s
}
#[cfg(test)]
mod tests {
    use super::*;
    fn make_simple_mesh() -> ExodusMesh {
        let mut mesh = ExodusMesh::new("test_mesh");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node(1.0, 1.0, 0.0);
        mesh.add_node(0.0, 1.0, 0.0);
        let bi = mesh.add_block(1, "solid", "QUAD4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        mesh
    }
    #[test]
    fn test_new_mesh_title() {
        let mesh = ExodusMesh::new("MyMesh");
        assert_eq!(mesh.title, "MyMesh");
    }
    #[test]
    fn test_new_mesh_empty() {
        let mesh = ExodusMesh::new("empty");
        assert_eq!(mesh.node_count(), 0);
        assert_eq!(mesh.element_count(), 0);
        assert!(mesh.blocks.is_empty());
        assert!(mesh.node_sets.is_empty());
        assert!(mesh.side_sets.is_empty());
    }
    #[test]
    fn test_add_node_returns_one_based_id() {
        let mut mesh = ExodusMesh::new("t");
        let id1 = mesh.add_node(1.0, 2.0, 3.0);
        let id2 = mesh.add_node(4.0, 5.0, 6.0);
        assert_eq!(id1, 1);
        assert_eq!(id2, 2);
    }
    #[test]
    fn test_node_count() {
        let mesh = make_simple_mesh();
        assert_eq!(mesh.node_count(), 4);
    }
    #[test]
    fn test_element_count_quad() {
        let mesh = make_simple_mesh();
        assert_eq!(mesh.element_count(), 1);
    }
    #[test]
    fn test_add_block_returns_index() {
        let mut mesh = ExodusMesh::new("t");
        let idx0 = mesh.add_block(1, "a", "TET4", 4);
        let idx1 = mesh.add_block(2, "b", "HEX8", 8);
        assert_eq!(idx0, 0);
        assert_eq!(idx1, 1);
    }
    #[test]
    fn test_block_properties() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_block(42, "myblock", "HEX8", 8);
        assert_eq!(mesh.blocks[0].id, 42);
        assert_eq!(mesh.blocks[0].name, "myblock");
        assert_eq!(mesh.blocks[0].element_type, "HEX8");
        assert_eq!(mesh.blocks[0].n_nodes_per_elem, 8);
    }
    #[test]
    fn test_add_node_set() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node_set(10, "bottom", vec![1, 2]);
        assert_eq!(mesh.node_sets.len(), 1);
        assert_eq!(mesh.node_sets[0].id, 10);
        assert_eq!(mesh.node_sets[0].name, "bottom");
        assert_eq!(mesh.node_sets[0].node_ids, vec![1, 2]);
    }
    #[test]
    fn test_element_count_multiple_blocks() {
        let mut mesh = ExodusMesh::new("t");
        for i in 0..8 {
            mesh.add_node(i as f64, 0.0, 0.0);
        }
        let b0 = mesh.add_block(1, "b0", "QUAD4", 4);
        mesh.add_element_to_block(b0, &[1, 2, 3, 4]);
        mesh.add_element_to_block(b0, &[5, 6, 7, 8]);
        let b1 = mesh.add_block(2, "b1", "TET4", 4);
        mesh.add_element_to_block(b1, &[1, 2, 3, 4]);
        assert_eq!(mesh.element_count(), 3);
    }
    #[test]
    fn test_bounding_box_single_node() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(3.0, 4.0, 5.0);
        let (min, max) = bounding_box(&mesh);
        assert_eq!(min, [3.0, 4.0, 5.0]);
        assert_eq!(max, [3.0, 4.0, 5.0]);
    }
    #[test]
    fn test_bounding_box_multiple_nodes() {
        let mesh = make_simple_mesh();
        let (min, max) = bounding_box(&mesh);
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0, 0.0]);
    }
    #[test]
    fn test_bounding_box_empty_mesh() {
        let mesh = ExodusMesh::new("empty");
        let (min, max) = bounding_box(&mesh);
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_bounding_box_negative_coords() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(-5.0, -3.0, -1.0);
        mesh.add_node(2.0, 4.0, 6.0);
        let (min, max) = bounding_box(&mesh);
        assert_eq!(min, [-5.0, -3.0, -1.0]);
        assert_eq!(max, [2.0, 4.0, 6.0]);
    }
    #[test]
    fn test_translate_mesh() {
        let mut mesh = make_simple_mesh();
        translate_mesh(&mut mesh, [1.0, 2.0, 3.0]);
        assert_eq!(mesh.coordinates[0], [1.0, 2.0, 3.0]);
        assert_eq!(mesh.coordinates[1], [2.0, 2.0, 3.0]);
        assert_eq!(mesh.coordinates[2], [2.0, 3.0, 3.0]);
        assert_eq!(mesh.coordinates[3], [1.0, 3.0, 3.0]);
    }
    #[test]
    fn test_translate_mesh_zero() {
        let mut mesh = make_simple_mesh();
        let orig: Vec<[f64; 3]> = mesh.coordinates.clone();
        translate_mesh(&mut mesh, [0.0, 0.0, 0.0]);
        assert_eq!(mesh.coordinates, orig);
    }
    #[test]
    fn test_write_text_roundtrip() {
        let mesh = make_simple_mesh();
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse failed");
        assert_eq!(parsed.title, mesh.title);
        assert_eq!(parsed.node_count(), mesh.node_count());
        assert_eq!(parsed.element_count(), mesh.element_count());
    }
    #[test]
    fn test_write_text_contains_header() {
        let mesh = make_simple_mesh();
        let text = write_text(&mesh);
        assert!(text.starts_with("EXODUS_MESH\n"));
        assert!(text.contains("TITLE test_mesh"));
        assert!(text.contains("END"));
    }
    #[test]
    fn test_parse_text_coordinates() {
        let mesh = make_simple_mesh();
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse failed");
        assert!((parsed.coordinates[0][0] - 0.0).abs() < 1e-12);
        assert!((parsed.coordinates[1][0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_parse_text_block_type() {
        let mesh = make_simple_mesh();
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse failed");
        assert_eq!(parsed.blocks[0].element_type, "QUAD4");
        assert_eq!(parsed.blocks[0].n_nodes_per_elem, 4);
    }
    #[test]
    fn test_parse_text_bad_header() {
        let result = parse_text("BAD_HEADER\n");
        assert!(result.is_err());
    }
    #[test]
    fn test_parse_text_node_sets_roundtrip() {
        let mut mesh = ExodusMesh::new("ns_test");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node_set(5, "left", vec![1]);
        mesh.add_node_set(6, "right", vec![2]);
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse failed");
        assert_eq!(parsed.node_sets.len(), 2);
        assert_eq!(parsed.node_sets[0].id, 5);
        assert_eq!(parsed.node_sets[0].name, "left");
        assert_eq!(parsed.node_sets[0].node_ids, vec![1]);
        assert_eq!(parsed.node_sets[1].id, 6);
    }
    #[test]
    fn test_side_sets_roundtrip() {
        let mut mesh = ExodusMesh::new("ss_test");
        for i in 0..4 {
            mesh.add_node(i as f64, 0.0, 0.0);
        }
        let bi = mesh.add_block(1, "b", "QUAD4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        mesh.side_sets.push(ExodusSideSet {
            id: 10,
            name: "top".to_string(),
            element_ids: vec![1],
            side_ids: vec![3],
        });
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse failed");
        assert_eq!(parsed.side_sets.len(), 1);
        assert_eq!(parsed.side_sets[0].id, 10);
        assert_eq!(parsed.side_sets[0].element_ids, vec![1]);
        assert_eq!(parsed.side_sets[0].side_ids, vec![3]);
    }
    #[test]
    fn test_exodus_result_write_summary() {
        let mut mesh = ExodusMesh::new("result_test");
        mesh.add_node(0.0, 0.0, 0.0);
        let bi = mesh.add_block(1, "b", "TET4", 4);
        mesh.add_element_to_block(bi, &[1, 1, 1, 1]);
        let result = ExodusResult {
            mesh,
            time_steps: vec![
                ExodusTimeStep {
                    time: 0.0,
                    variables: vec![ExodusVariable {
                        name: "displacement".to_string(),
                        entity_type: "node".to_string(),
                        values: vec![0.0],
                    }],
                },
                ExodusTimeStep {
                    time: 1.0,
                    variables: vec![],
                },
            ],
        };
        let summary = result.write_summary();
        assert!(summary.contains("Title: result_test"));
        assert!(summary.contains("Nodes: 1"));
        assert!(summary.contains("Time Steps: 2"));
        assert!(summary.contains("t=0.000000"));
        assert!(summary.contains("t=1.000000"));
    }
    #[test]
    fn test_hex8_mesh() {
        let mut mesh = ExodusMesh::new("hex_mesh");
        let coords = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 1.0),
        ];
        for (x, y, z) in coords {
            mesh.add_node(x, y, z);
        }
        let bi = mesh.add_block(1, "vol", "HEX8", 8);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(mesh.node_count(), 8);
        assert_eq!(mesh.element_count(), 1);
        let (min, max) = bounding_box(&mesh);
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0, 1.0]);
    }
    #[test]
    fn test_write_text_empty_mesh() {
        let mesh = ExodusMesh::new("empty");
        let text = write_text(&mesh);
        let parsed = parse_text(&text).expect("parse of empty mesh failed");
        assert_eq!(parsed.node_count(), 0);
        assert_eq!(parsed.element_count(), 0);
    }
    #[test]
    fn test_linear_time_steps_count() {
        let ts = linear_time_steps(0.0, 1.0, 5);
        assert_eq!(ts.len(), 5);
        assert!((ts[0] - 0.0).abs() < 1e-12);
        assert!((ts[4] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_linear_time_steps_spacing() {
        let ts = linear_time_steps(0.0, 4.0, 5);
        for w in ts.windows(2) {
            assert!((w[1] - w[0] - 1.0).abs() < 1e-12, "Spacing should be 1.0");
        }
    }
    #[test]
    fn test_linear_time_steps_empty() {
        let ts = linear_time_steps(0.0, 1.0, 0);
        assert!(ts.is_empty());
    }
    #[test]
    fn test_linear_time_steps_single() {
        let ts = linear_time_steps(3.0, 7.0, 1);
        assert_eq!(ts.len(), 1);
        assert!((ts[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_log_time_steps_count() {
        let ts = log_time_steps(1.0, 100.0, 3);
        assert_eq!(ts.len(), 3);
        assert!((ts[0] - 1.0).abs() < 1e-9);
        assert!((ts[2] - 100.0).abs() < 1e-6);
    }
    #[test]
    fn test_log_time_steps_monotone() {
        let ts = log_time_steps(0.1, 10.0, 10);
        for w in ts.windows(2) {
            assert!(w[1] > w[0], "Log steps should be monotonically increasing");
        }
    }
    #[test]
    fn test_add_node_variable_step() {
        let mesh = ExodusMesh::new("t");
        let mut result = ExodusResult {
            mesh,
            time_steps: vec![],
        };
        add_node_variable_step(&mut result, 1.5, "temperature", vec![300.0, 310.0]);
        assert_eq!(result.time_steps.len(), 1);
        assert!((result.time_steps[0].time - 1.5).abs() < 1e-12);
        assert_eq!(result.time_steps[0].variables[0].name, "temperature");
        assert_eq!(result.time_steps[0].variables[0].values.len(), 2);
    }
    #[test]
    fn test_add_element_variable_step() {
        let mesh = ExodusMesh::new("t");
        let mut result = ExodusResult {
            mesh,
            time_steps: vec![],
        };
        add_element_variable_step(&mut result, 2.0, "stress", vec![1e6]);
        assert_eq!(result.time_steps[0].variables[0].entity_type, "element");
    }
    #[test]
    fn test_validate_valid_mesh() {
        let mesh = make_simple_mesh();
        let report = validate_mesh(&mesh);
        assert!(
            report.is_valid,
            "Simple mesh should be valid: {:?}",
            report.errors
        );
    }
    #[test]
    fn test_validate_empty_mesh_has_warnings() {
        let mesh = ExodusMesh::new("empty");
        let report = validate_mesh(&mesh);
        assert!(report.is_valid);
        assert!(
            !report.warnings.is_empty(),
            "Empty mesh should have warnings"
        );
    }
    #[test]
    fn test_validate_bad_node_id() {
        let mut mesh = ExodusMesh::new("bad");
        mesh.add_node(0.0, 0.0, 0.0);
        let bi = mesh.add_block(1, "b", "QUAD4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        let report = validate_mesh(&mesh);
        assert!(!report.is_valid, "Invalid node IDs should fail validation");
        assert!(!report.errors.is_empty());
    }
    #[test]
    fn test_validate_duplicate_block_id() {
        let mut mesh = ExodusMesh::new("dup");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_block(1, "a", "TET4", 4);
        mesh.add_block(1, "b", "TET4", 4);
        let report = validate_mesh(&mesh);
        assert!(!report.is_valid);
        assert!(
            report
                .errors
                .iter()
                .any(|e| e.contains("Duplicate block id"))
        );
    }
    #[test]
    fn test_validate_mismatched_connectivity_length() {
        let mut mesh = ExodusMesh::new("mis");
        for _ in 0..4 {
            mesh.add_node(0.0, 0.0, 0.0);
        }
        let bi = mesh.add_block(1, "b", "QUAD4", 4);
        mesh.blocks[bi].connectivity = vec![1, 2, 3];
        let report = validate_mesh(&mesh);
        assert!(!report.is_valid);
        assert!(report.errors.iter().any(|e| e.contains("not divisible")));
    }
    #[test]
    fn test_scale_mesh() {
        let mut mesh = make_simple_mesh();
        scale_mesh(&mut mesh, 2.0);
        assert!((mesh.coordinates[1][0] - 2.0).abs() < 1e-12);
        assert!((mesh.coordinates[2][0] - 2.0).abs() < 1e-12);
        assert!((mesh.coordinates[2][1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_mesh_centroid() {
        let mesh = make_simple_mesh();
        let c = mesh_centroid(&mesh);
        assert!((c[0] - 0.5).abs() < 1e-12);
        assert!((c[1] - 0.5).abs() < 1e-12);
        assert!(c[2].abs() < 1e-12);
    }
    #[test]
    fn test_mesh_centroid_empty() {
        let mesh = ExodusMesh::new("empty");
        let c = mesh_centroid(&mesh);
        assert_eq!(c, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_build_unit_hex_mesh() {
        let mesh = build_unit_hex_mesh("test_hex");
        assert_eq!(mesh.node_count(), 8);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.blocks[0].element_type, "HEX8");
        let report = validate_mesh(&mesh);
        assert!(
            report.is_valid,
            "Unit HEX8 mesh should be valid: {:?}",
            report.errors
        );
    }
    #[test]
    fn test_build_unit_tet_mesh() {
        let mesh = build_unit_tet_mesh("test_tet");
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.blocks[0].element_type, "TET4");
        let report = validate_mesh(&mesh);
        assert!(report.is_valid);
    }
    #[test]
    fn test_merge_meshes_node_count() {
        let a = make_simple_mesh();
        let b = build_unit_tet_mesh("b");
        let merged = merge_meshes(&a, &b);
        assert_eq!(merged.node_count(), 8);
    }
    #[test]
    fn test_merge_meshes_element_count() {
        let a = make_simple_mesh();
        let b = build_unit_tet_mesh("b");
        let merged = merge_meshes(&a, &b);
        assert_eq!(merged.element_count(), 2);
    }
    #[test]
    fn test_unit_hex_mesh_validation() {
        let mesh = build_unit_hex_mesh("validation_test");
        let report = validate_mesh(&mesh);
        assert!(report.is_valid);
        assert!(report.errors.is_empty());
    }
}
