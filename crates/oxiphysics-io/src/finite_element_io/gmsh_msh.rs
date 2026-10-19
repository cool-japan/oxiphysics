//! Reader for Gmsh `.msh` files in the ASCII formats 2.x and 4.1.
//!
//! Both layouts are section based (`$Nodes` ... `$EndNodes`) but differ in how
//! a section is organised:
//!
//! * **MSH 2.x** lists `node-number x y z` and
//!   `elm-number elm-type number-of-tags <tag>... node-number-list` lines
//!   behind a count line.
//! * **MSH 4.1** groups nodes and elements into entity blocks.  A `$Nodes`
//!   block is `entityDim entityTag parametric numNodes`, followed by the
//!   `numNodes` node tags (one per line) and then the `numNodes` coordinate
//!   lines (`x y z`, plus parametric coordinates when `parametric` is 1).  An
//!   `$Elements` block is `entityDim entityTag elementType numElements`,
//!   followed by `elementTag nodeTag...` lines.
//!
//! MSH 4.0 uses yet another node layout and binary files cannot be read from a
//! string, so both are rejected with an error rather than misread.

use super::types::{FeElement, FeElementType, FeNode, GmshMesh, GmshVersion};
use crate::Error;

/// Which section layout the file uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// MSH 2.x.
    V2,
    /// MSH 4.1.
    V41,
}

/// Element type for a Gmsh element-type number.
///
/// Types without a [`FeElementType`] counterpart are kept as
/// `Unknown("gmsh:<number>")` so that they survive a round trip.
pub(super) fn element_type_from_tag(tag: usize) -> FeElementType {
    match tag {
        1 => FeElementType::Line2,
        2 => FeElementType::Tri3,
        3 => FeElementType::Quad4,
        4 => FeElementType::Tet4,
        5 => FeElementType::Hex8,
        9 => FeElementType::Tri6,
        11 => FeElementType::Tet10,
        16 => FeElementType::Quad8,
        17 => FeElementType::Hex20,
        _ => FeElementType::Unknown(format!("gmsh:{tag}")),
    }
}

/// Gmsh element-type number for an element type; `0` when there is none.
pub(super) fn tag_of_element_type(element_type: &FeElementType) -> usize {
    match element_type {
        FeElementType::Line2 => 1,
        FeElementType::Tri3 => 2,
        FeElementType::Quad4 => 3,
        FeElementType::Tet4 => 4,
        FeElementType::Hex8 => 5,
        FeElementType::Tri6 => 9,
        FeElementType::Tet10 => 11,
        FeElementType::Quad8 => 16,
        FeElementType::Hex20 => 17,
        FeElementType::Unknown(name) => name
            .strip_prefix("gmsh:")
            .and_then(|tag| tag.parse().ok())
            .unwrap_or(0),
    }
}

/// Number of nodes of the Gmsh element types that map to a
/// [`FeElementType`]; used to catch truncated element lines.
fn node_count_of_tag(tag: usize) -> Option<usize> {
    element_type_from_tag(tag).node_count()
}

fn error(line: usize, message: impl AsRef<str>) -> Error {
    Error::Parse(format!("Gmsh MSH line {line}: {}", message.as_ref()))
}

/// Non-blank lines of the file with their 1-based numbers.
struct Lines<'a> {
    inner: std::str::Lines<'a>,
    number: usize,
}

impl<'a> Lines<'a> {
    fn new(source: &'a str) -> Self {
        Self {
            inner: source.lines(),
            number: 0,
        }
    }

    /// The next line that is not blank, trimmed, with its number.
    fn next_line(&mut self) -> Option<(usize, &'a str)> {
        for line in self.inner.by_ref() {
            self.number += 1;
            let text = line.trim();
            if !text.is_empty() {
                return Some((self.number, text));
            }
        }
        None
    }

    /// Like [`Lines::next_line`], but running out of lines is an error.
    fn require(&mut self, what: &str) -> Result<(usize, &'a str), Error> {
        self.next_line().ok_or_else(|| {
            error(
                self.number,
                format!("unexpected end of file while reading {what}"),
            )
        })
    }

    /// Consume the line that has to be exactly `marker`.
    fn expect_line(&mut self, marker: &str) -> Result<(), Error> {
        let (number, line) = self.require(marker)?;
        if line == marker {
            Ok(())
        } else {
            Err(error(number, format!("expected {marker}, found `{line}`")))
        }
    }
}

fn parse_usize(token: &str, what: &str, line: usize) -> Result<usize, Error> {
    token
        .parse::<usize>()
        .map_err(|_| error(line, format!("`{token}` is not a valid {what}")))
}

fn parse_f64(token: &str, what: &str, line: usize) -> Result<f64, Error> {
    token
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| error(line, format!("`{token}` is not a valid {what}")))
}

/// The first `N` whitespace separated fields of `line` as counts.
fn parse_counts<const N: usize>(
    line: &str,
    number: usize,
    what: &str,
) -> Result<[usize; N], Error> {
    let mut values = [0usize; N];
    let mut tokens = line.split_whitespace();
    for slot in &mut values {
        let token = tokens
            .next()
            .ok_or_else(|| error(number, format!("{what} needs {N} fields, found `{line}`")))?;
        *slot = parse_usize(token, what, number)?;
    }
    Ok(values)
}

/// What `$MeshFormat` declares.
#[derive(Debug, Clone, Copy)]
struct MeshFormat {
    layout: Layout,
}

impl MeshFormat {
    fn version(self) -> GmshVersion {
        match self.layout {
            Layout::V2 => GmshVersion::V2,
            Layout::V41 => GmshVersion::V4,
        }
    }
}

/// Read a Gmsh ASCII mesh into `mesh`.
///
/// # Errors
///
/// Fails with [`Error::Parse`] for binary files and for versions other than
/// 2.x and 4.1, for sections that are cut off or whose counts disagree with
/// their content, and for malformed numbers and element lines.
pub(super) fn parse_into(mesh: &mut GmshMesh, source: &str) -> Result<(), Error> {
    let mut lines = Lines::new(source);
    let mut format: Option<MeshFormat> = None;
    while let Some((number, line)) = lines.next_line() {
        match line {
            "$MeshFormat" => {
                let parsed = read_mesh_format(&mut lines)?;
                mesh.version = Some(parsed.version());
                format = Some(parsed);
            }
            "$PhysicalNames" => read_physical_names(&mut lines, mesh)?,
            "$Nodes" => {
                let layout = declared_layout(format, number, "$Nodes")?;
                match layout {
                    Layout::V2 => read_nodes_v2(&mut lines, mesh)?,
                    Layout::V41 => read_nodes_v41(&mut lines, mesh)?,
                }
            }
            "$Elements" => {
                let layout = declared_layout(format, number, "$Elements")?;
                match layout {
                    Layout::V2 => read_elements_v2(&mut lines, mesh)?,
                    Layout::V41 => read_elements_v41(&mut lines, mesh)?,
                }
            }
            section if section.starts_with('$') && !section.starts_with("$End") => {
                skip_section(&mut lines, section)?;
            }
            _ => {}
        }
    }
    if format.is_none() {
        return Err(Error::Parse(
            "Gmsh MSH: the file has no $MeshFormat section".to_string(),
        ));
    }
    Ok(())
}

fn declared_layout(
    format: Option<MeshFormat>,
    number: usize,
    section: &str,
) -> Result<Layout, Error> {
    format
        .map(|f| f.layout)
        .ok_or_else(|| error(number, format!("{section} appears before $MeshFormat")))
}

/// Skip a section this reader does not use, such as `$Entities`.
fn skip_section(lines: &mut Lines<'_>, header: &str) -> Result<(), Error> {
    let end = format!("$End{}", header.get(1..).unwrap_or(""));
    loop {
        let (_, line) = lines.require(&format!("section {header}"))?;
        if line == end {
            return Ok(());
        }
    }
}

fn read_mesh_format(lines: &mut Lines<'_>) -> Result<MeshFormat, Error> {
    let (number, line) = lines.require("the $MeshFormat line")?;
    let mut tokens = line.split_whitespace();
    let (Some(version), Some(file_type)) = (tokens.next(), tokens.next()) else {
        return Err(error(
            number,
            format!("$MeshFormat needs `version file-type data-size`, found `{line}`"),
        ));
    };
    let layout = match split_version(version) {
        Some((2, _)) => Layout::V2,
        Some((4, 1)) => Layout::V41,
        Some((4, minor)) => {
            return Err(error(
                number,
                format!(
                    "MSH version 4.{minor} is not supported; only 4.1 is (save with Mesh.MshFileVersion = 4.1)"
                ),
            ));
        }
        _ => {
            return Err(error(
                number,
                format!("MSH version {version} is not supported; only 2.x and 4.1 are"),
            ));
        }
    };
    match file_type {
        "0" => {}
        "1" => {
            return Err(error(
                number,
                "binary MSH files are not supported; save the mesh in ASCII format (Mesh.Binary = 0)",
            ));
        }
        other => {
            return Err(error(
                number,
                format!("`{other}` is not a valid MSH file type (0 = ASCII, 1 = binary)"),
            ));
        }
    }
    lines.expect_line("$EndMeshFormat")?;
    Ok(MeshFormat { layout })
}

/// `(major, minor)` of a version string such as `2.2` or `4.1`.
fn split_version(text: &str) -> Option<(u32, u32)> {
    let mut parts = text.split('.');
    let major = parts.next()?.parse().ok()?;
    let minor = match parts.next() {
        Some(minor) => minor.parse().ok()?,
        None => 0,
    };
    parts.next().is_none().then_some((major, minor))
}

/// Split `text` after its first white-space delimited token.
fn split_first_token(text: &str) -> Option<(&str, &str)> {
    let text = text.trim_start();
    if text.is_empty() {
        return None;
    }
    Some(match text.find(char::is_whitespace) {
        Some(end) => text.split_at(end),
        None => (text, ""),
    })
}

fn read_physical_names(lines: &mut Lines<'_>, mesh: &mut GmshMesh) -> Result<(), Error> {
    let (count_line, line) = lines.require("the physical name count")?;
    let declared = parse_usize(
        line.split_whitespace().next().unwrap_or(""),
        "physical name count",
        count_line,
    )?;
    let mut read = 0usize;
    loop {
        let (number, line) = lines.require("$PhysicalNames")?;
        if line == "$EndPhysicalNames" {
            break;
        }
        let malformed = || error(number, format!("malformed physical name `{line}`"));
        let (dim, rest) = split_first_token(line).ok_or_else(malformed)?;
        let (tag, rest) = split_first_token(rest).ok_or_else(malformed)?;
        let dim = dim.parse::<i32>().map_err(|_| malformed())?;
        let tag = tag.parse::<i32>().map_err(|_| malformed())?;
        mesh.physical_names
            .insert((dim, tag), rest.trim().trim_matches('"').to_string());
        read += 1;
    }
    if read != declared {
        return Err(error(
            count_line,
            format!("$PhysicalNames declares {declared} names but {read} are present"),
        ));
    }
    Ok(())
}

fn read_nodes_v2(lines: &mut Lines<'_>, mesh: &mut GmshMesh) -> Result<(), Error> {
    let (count_line, line) = lines.require("the node count")?;
    let declared = parse_usize(
        line.split_whitespace().next().unwrap_or(""),
        "node count",
        count_line,
    )?;
    let mut read = 0usize;
    loop {
        let (number, line) = lines.require("$Nodes")?;
        if line == "$EndNodes" {
            break;
        }
        let mut tokens = line.split_whitespace();
        let (Some(id), Some(x), Some(y), Some(z)) =
            (tokens.next(), tokens.next(), tokens.next(), tokens.next())
        else {
            return Err(error(
                number,
                format!("a node line needs `node-number x y z`, found `{line}`"),
            ));
        };
        mesh.nodes.push(FeNode::new(
            parse_usize(id, "node number", number)?,
            [
                parse_f64(x, "coordinate", number)?,
                parse_f64(y, "coordinate", number)?,
                parse_f64(z, "coordinate", number)?,
            ],
        ));
        read += 1;
    }
    if read != declared {
        return Err(error(
            count_line,
            format!("$Nodes declares {declared} nodes but {read} are present"),
        ));
    }
    Ok(())
}

fn read_nodes_v41(lines: &mut Lines<'_>, mesh: &mut GmshMesh) -> Result<(), Error> {
    let (header_line, line) = lines.require("the $Nodes header")?;
    let [blocks, declared, _min_tag, _max_tag] = parse_counts::<4>(
        line,
        header_line,
        "$Nodes header `numEntityBlocks numNodes minNodeTag maxNodeTag`",
    )?;
    let mut read = 0usize;
    for _ in 0..blocks {
        let (number, line) = lines.require("a node block header")?;
        let [_dim, _tag, parametric, in_block] = parse_counts::<4>(
            line,
            number,
            "node block header `entityDim entityTag parametric numNodesInBlock`",
        )?;
        if parametric > 1 {
            return Err(error(
                number,
                format!("`{parametric}` is not a valid parametric flag (0 or 1)"),
            ));
        }
        // Tags and coordinates are stored one after the other, so the tags are
        // collected first.  Nothing is reserved from the declared count: a
        // corrupt count just runs into the end of the file.
        let mut tags = Vec::new();
        for _ in 0..in_block {
            let (number, line) = lines.require("a node tag")?;
            let token = line.split_whitespace().next().unwrap_or("");
            tags.push(parse_usize(token, "node tag", number)?);
        }
        for tag in tags {
            let (number, line) = lines.require("node coordinates")?;
            let mut tokens = line.split_whitespace();
            let (Some(x), Some(y), Some(z)) = (tokens.next(), tokens.next(), tokens.next()) else {
                return Err(error(
                    number,
                    format!("a node needs the coordinates `x y z`, found `{line}`"),
                ));
            };
            // Parametric coordinates, if any, follow x, y and z and are not
            // part of a node's position.
            mesh.nodes.push(FeNode::new(
                tag,
                [
                    parse_f64(x, "coordinate", number)?,
                    parse_f64(y, "coordinate", number)?,
                    parse_f64(z, "coordinate", number)?,
                ],
            ));
            read += 1;
        }
    }
    lines.expect_line("$EndNodes")?;
    if read != declared {
        return Err(error(
            header_line,
            format!("$Nodes declares {declared} nodes but {read} are present"),
        ));
    }
    Ok(())
}

/// Check an element's connectivity against what its type needs.
fn check_element(id: usize, type_tag: usize, nodes: &[usize], number: usize) -> Result<(), Error> {
    if nodes.is_empty() {
        return Err(error(number, format!("element {id} has no nodes")));
    }
    match node_count_of_tag(type_tag) {
        Some(expected) if expected != nodes.len() => Err(error(
            number,
            format!(
                "element {id} of Gmsh type {type_tag} needs {expected} nodes, found {}",
                nodes.len()
            ),
        )),
        _ => Ok(()),
    }
}

fn parse_node_list(tokens: &[&str], number: usize) -> Result<Vec<usize>, Error> {
    tokens
        .iter()
        .map(|token| parse_usize(token, "node number", number))
        .collect()
}

fn read_elements_v2(lines: &mut Lines<'_>, mesh: &mut GmshMesh) -> Result<(), Error> {
    let (count_line, line) = lines.require("the element count")?;
    let declared = parse_usize(
        line.split_whitespace().next().unwrap_or(""),
        "element count",
        count_line,
    )?;
    let mut read = 0usize;
    loop {
        let (number, line) = lines.require("$Elements")?;
        if line == "$EndElements" {
            break;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        let [id, type_tag, tag_count] = match parts.as_slice() {
            [id, type_tag, tag_count, ..] => [*id, *type_tag, *tag_count],
            _ => {
                return Err(error(
                    number,
                    format!(
                        "an element line needs `elm-number elm-type number-of-tags ...`, found `{line}`"
                    ),
                ));
            }
        };
        let id = parse_usize(id, "element number", number)?;
        let type_tag = parse_usize(type_tag, "element type", number)?;
        let tag_count = parse_usize(tag_count, "number of tags", number)?;
        let node_tokens = 3usize
            .checked_add(tag_count)
            .and_then(|start| parts.get(start..))
            .ok_or_else(|| {
                error(
                    number,
                    format!(
                        "element {id} declares {tag_count} tags but its line has only {} fields",
                        parts.len()
                    ),
                )
            })?;
        let nodes = parse_node_list(node_tokens, number)?;
        check_element(id, type_tag, &nodes, number)?;
        mesh.elements
            .push(FeElement::new(id, element_type_from_tag(type_tag), nodes));
        read += 1;
    }
    if read != declared {
        return Err(error(
            count_line,
            format!("$Elements declares {declared} elements but {read} are present"),
        ));
    }
    Ok(())
}

fn read_elements_v41(lines: &mut Lines<'_>, mesh: &mut GmshMesh) -> Result<(), Error> {
    let (header_line, line) = lines.require("the $Elements header")?;
    let [blocks, declared, _min_tag, _max_tag] = parse_counts::<4>(
        line,
        header_line,
        "$Elements header `numEntityBlocks numElements minElementTag maxElementTag`",
    )?;
    let mut read = 0usize;
    for _ in 0..blocks {
        let (number, line) = lines.require("an element block header")?;
        let [_dim, _tag, type_tag, in_block] = parse_counts::<4>(
            line,
            number,
            "element block header `entityDim entityTag elementType numElementsInBlock`",
        )?;
        let element_type = element_type_from_tag(type_tag);
        for _ in 0..in_block {
            let (number, line) = lines.require("an element")?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            let Some((id, node_tokens)) = parts.split_first() else {
                return Err(error(
                    number,
                    "an element line needs `elementTag nodeTag...`",
                ));
            };
            let id = parse_usize(id, "element tag", number)?;
            let nodes = parse_node_list(node_tokens, number)?;
            check_element(id, type_tag, &nodes, number)?;
            mesh.elements
                .push(FeElement::new(id, element_type.clone(), nodes));
            read += 1;
        }
    }
    lines.expect_line("$EndElements")?;
    if read != declared {
        return Err(error(
            header_line,
            format!("$Elements declares {declared} elements but {read} are present"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(source: &str) -> Result<GmshMesh, Error> {
        GmshMesh::parse(source)
    }

    fn message(source: &str) -> String {
        parse(source).unwrap_err().to_string()
    }

    /// A unit cube meshed with one hexahedron, as `gmsh -3 -format msh41`
    /// writes it: eight point blocks, one quadrilateral surface block and
    /// the volume block.
    const UNIT_CUBE_41: &str = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 1 \"bottom face\"
3 2 \"cube\"
$EndPhysicalNames
$Entities
8 0 0 1
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
5 0 0 1 0
6 1 0 1 0
7 1 1 1 0
8 0 1 1 0
1 0 0 0 1 1 1 1 2 0 0
$EndEntities
$Nodes
8 8 1 8
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
0 5 0 1
5
0 0 1
0 6 0 1
6
1 0 1
0 7 0 1
7
1 1 1
0 8 0 1
8
0 1 1
$EndNodes
$Elements
2 2 1 2
2 1 3 1
1 1 2 3 4
3 1 5 1
2 1 2 3 4 5 6 7 8
$EndElements
";

    #[test]
    fn msh41_unit_cube_is_read_block_by_block() {
        let mesh = parse(UNIT_CUBE_41).unwrap();
        assert_eq!(mesh.version, Some(GmshVersion::V4));
        assert_eq!(mesh.nodes.len(), 8);
        for (index, node) in mesh.nodes.iter().enumerate() {
            assert_eq!(node.id, index + 1);
        }
        assert_eq!(mesh.nodes[0].coords, [0.0, 0.0, 0.0]);
        assert_eq!(mesh.nodes[2].coords, [1.0, 1.0, 0.0]);
        assert_eq!(mesh.nodes[6].coords, [1.0, 1.0, 1.0]);
        assert_eq!(mesh.nodes[7].coords, [0.0, 1.0, 1.0]);
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[0].id, 1);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Quad4);
        assert_eq!(mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
        assert_eq!(mesh.elements[1].id, 2);
        assert_eq!(mesh.elements[1].element_type, FeElementType::Hex8);
        assert_eq!(mesh.elements[1].connectivity, vec![1, 2, 3, 4, 5, 6, 7, 8]);
        assert_eq!(
            mesh.physical_names.get(&(3, 2)).map(String::as_str),
            Some("cube")
        );
        assert_eq!(
            mesh.physical_names.get(&(2, 1)).map(String::as_str),
            Some("bottom face")
        );
    }

    #[test]
    fn msh41_blocks_with_many_nodes_and_sparse_tags() {
        let source = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
2 5 10 50
1 7 0 2
10
20
0.5 0 0
1.5 0 0
3 9 0 3
30
40
50
-1.0e-1 2.5E+1 3
4 5 6
7 8 9
$EndNodes
$Elements
1 2 100 101
3 9 4 2
100 10 20 30 40
101 20 30 40 50
$EndElements
";
        let mesh = parse(source).unwrap();
        let ids: Vec<usize> = mesh.nodes.iter().map(|n| n.id).collect();
        assert_eq!(ids, vec![10, 20, 30, 40, 50]);
        assert_eq!(mesh.nodes[1].coords, [1.5, 0.0, 0.0]);
        assert_eq!(mesh.nodes[2].coords, [-0.1, 25.0, 3.0]);
        assert_eq!(mesh.nodes[4].coords, [7.0, 8.0, 9.0]);
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[1].id, 101);
        assert_eq!(mesh.elements[1].element_type, FeElementType::Tet4);
        assert_eq!(mesh.elements[1].connectivity, vec![20, 30, 40, 50]);
    }

    #[test]
    fn msh41_parametric_coordinates_are_not_part_of_the_position() {
        let source = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 2 1 2
1 3 1 2
1
2
0 0 0 0.25
1 0 0 0.75
$EndNodes
$Elements
1 1 1 1
1 3 1 1
1 1 2
$EndElements
";
        let mesh = parse(source).unwrap();
        assert_eq!(mesh.nodes[0].coords, [0.0, 0.0, 0.0]);
        assert_eq!(mesh.nodes[1].coords, [1.0, 0.0, 0.0]);
        assert_eq!(mesh.elements[0].element_type, FeElementType::Line2);
    }

    #[test]
    fn the_old_v4_element_panic_input_is_now_read() {
        // A block header `3 1 4 2` has four fields; the v2 reader treated the
        // fourth as a tag count and sliced past the end of the line.
        let source = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 2 1 2
3 1 4 2
1 1 2 3 4
2 4 3 2 1
$EndElements
";
        let mesh = parse(source).unwrap();
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.elements.len(), 2);
        assert_eq!(mesh.elements[1].connectivity, vec![4, 3, 2, 1]);
    }

    #[test]
    fn msh2_reads_tags_and_second_order_types() {
        let source = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
3
1 3 2 11 7 1 2 3 4
2 16 0 1 2 3 4 1 2 3 4
3 10 1 5 1 2 3 4 1 2 3 4 1
$EndElements
";
        let mesh = parse(source).unwrap();
        assert_eq!(mesh.version, Some(GmshVersion::V2));
        assert_eq!(mesh.elements[0].element_type, FeElementType::Quad4);
        assert_eq!(mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
        // Gmsh type 16 is the 8-node quadrangle, type 10 the 9-node one.
        assert_eq!(mesh.elements[1].element_type, FeElementType::Quad8);
        assert_eq!(
            mesh.elements[2].element_type,
            FeElementType::Unknown("gmsh:10".to_string())
        );
        assert_eq!(mesh.elements[2].connectivity.len(), 9);
    }

    #[test]
    fn msh2_tag_count_beyond_the_line_is_an_error_not_a_panic() {
        let source = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
0
$EndNodes
$Elements
1
1 4 9 1 2 3 4
$EndElements
";
        let text = message(source);
        assert!(
            text.contains("declares 9 tags") && text.contains("line 9"),
            "{text}"
        );
    }

    #[test]
    fn element_with_too_few_nodes_for_its_type_is_rejected() {
        let source = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
0
$EndNodes
$Elements
1
1 4 0 1 2 3
$EndElements
";
        let text = message(source);
        assert!(text.contains("needs 4 nodes, found 3"), "{text}");
    }

    #[test]
    fn binary_files_and_other_versions_are_rejected() {
        let binary = "$MeshFormat\n4.1 1 8\n$EndMeshFormat\n";
        assert!(message(binary).contains("binary"));
        let binary2 = "$MeshFormat\n2.2 1 8\n$EndMeshFormat\n";
        assert!(message(binary2).contains("binary"));
        let v40 = "$MeshFormat\n4.0 0 8\n$EndMeshFormat\n";
        assert!(message(v40).contains("4.0"));
        let v1 = "$MeshFormat\n1.0 0 8\n$EndMeshFormat\n";
        assert!(message(v1).contains("version 1.0"));
        let nonsense = "$MeshFormat\n4.1 7 8\n$EndMeshFormat\n";
        assert!(message(nonsense).contains("file type"));
    }

    #[test]
    fn missing_or_late_mesh_format_is_an_error() {
        assert!(message("").contains("$MeshFormat"));
        assert!(message("$Nodes\n0\n$EndNodes\n").contains("before $MeshFormat"));
        assert!(message("$Elements\n0\n$EndElements\n").contains("before $MeshFormat"));
    }

    #[test]
    fn declared_counts_have_to_match_the_content() {
        let v2 = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n1 0 0 0\n2 1 0 0\n$EndNodes\n";
        assert!(message(v2).contains("declares 3 nodes but 2"));
        let v2_elements =
            "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Elements\n2\n1 15 0 1\n$EndElements\n";
        assert!(message(v2_elements).contains("declares 2 elements but 1"));
        let v41 = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 3 1 3\n0 1 0 2\n1\n2\n0 0 0\n1 0 0\n$EndNodes\n";
        assert!(message(v41).contains("declares 3 nodes but 2"));
        let names = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n2\n1 1 \"a\"\n$EndPhysicalNames\n";
        assert!(message(names).contains("declares 2 names but 1"));
    }

    #[test]
    fn truncated_files_are_errors_at_every_depth() {
        let cut_nodes =
            "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 2 1 2\n0 1 0 2\n1\n2\n0 0 0\n";
        assert!(message(cut_nodes).contains("end of file"));
        let cut_section = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Entities\n1 2 3\n";
        assert!(message(cut_section).contains("$Entities"));
        let cut_format = "$MeshFormat\n2.2 0 8\n";
        assert!(message(cut_format).contains("$EndMeshFormat"));
        let huge_block =
            "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 1 1 1\n0 1 0 99999999999\n1\n";
        assert!(message(huge_block).contains("end of file"));
    }

    #[test]
    fn malformed_numbers_name_their_line() {
        let source = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n1\n1 0 zero 0\n$EndNodes\n";
        let text = message(source);
        assert!(text.contains("line 6") && text.contains("`zero`"), "{text}");
        let nan = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n1\n1 0 NaN 0\n$EndNodes\n";
        assert!(message(nan).contains("`NaN`"));
    }

    #[test]
    fn unknown_sections_and_stray_lines_are_skipped() {
        let source = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Comments
some text with $Nodes inside
$EndComments
this line belongs to no section
$NodeData
1
\"temperature\"
$EndNodeData
$Nodes
1
1 1 2 3
$EndNodes
";
        let mesh = parse(source).unwrap();
        assert_eq!(mesh.nodes.len(), 1);
        assert_eq!(mesh.nodes[0].coords, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn physical_names_keep_spaces_and_irregular_blanks() {
        let source = "\
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
2   7   \"inlet  plane\"
3 1 plain
$EndPhysicalNames
";
        let mesh = parse(source).unwrap();
        assert_eq!(
            mesh.physical_names.get(&(2, 7)).map(String::as_str),
            Some("inlet  plane")
        );
        assert_eq!(
            mesh.physical_names.get(&(3, 1)).map(String::as_str),
            Some("plain")
        );
        let bad = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n2 x \"a\"\n$EndPhysicalNames\n";
        assert!(message(bad).contains("malformed physical name"));
    }

    #[test]
    fn write_v2_round_trips_through_the_reader() {
        let mut mesh = GmshMesh::new();
        mesh.nodes.push(FeNode::new(1, [0.0, 0.0, 0.0]));
        mesh.nodes.push(FeNode::new(2, [1.0, 0.5, -0.25]));
        mesh.nodes.push(FeNode::new(3, [0.0, 1.0, 2.0e-3]));
        mesh.elements
            .push(FeElement::new(1, FeElementType::Tri3, vec![1, 2, 3]));
        mesh.elements.push(FeElement::new(
            2,
            FeElementType::Quad8,
            vec![1, 2, 3, 1, 2, 3, 1, 2],
        ));
        mesh.elements.push(FeElement::new(
            3,
            FeElementType::Unknown("gmsh:10".to_string()),
            vec![1, 2, 3, 1, 2, 3, 1, 2, 3],
        ));
        let text = mesh.write_v2();
        let back = parse(&text).unwrap();
        assert_eq!(back.nodes, mesh.nodes);
        assert_eq!(back.elements, mesh.elements);
        assert_eq!(back.version, Some(GmshVersion::V2));
    }

    #[test]
    fn element_type_numbers_map_both_ways() {
        for tag in [1, 2, 3, 4, 5, 9, 11, 16, 17] {
            assert_eq!(tag_of_element_type(&element_type_from_tag(tag)), tag);
        }
        assert_eq!(
            tag_of_element_type(&FeElementType::Unknown("gmsh:21".to_string())),
            21
        );
        assert_eq!(
            tag_of_element_type(&FeElementType::Unknown("x".to_string())),
            0
        );
    }
}
