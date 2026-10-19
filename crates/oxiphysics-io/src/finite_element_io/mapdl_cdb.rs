//! The fixed-format blocks of ANSYS MAPDL CDB decks.
//!
//! `CDWRITE` and the workbench solver input write nodes and elements as
//! `NBLOCK` / `EBLOCK` blocks.  The line after the command carries a Fortran
//! format (`(3i9,6e21.13e3)`, `(19i9)`, ...) and every row is cut into fields
//! of exactly those widths; a block ends at a `-1` row or, for older decks, at
//! the `N,R5.3,LOC,-1` line.
//!
//! * `NBLOCK` rows hold `node, solid-model entity, line location, x, y, z`
//!   followed by optional rotation angles.  Trailing fields that a short row
//!   leaves out are zero.
//! * `EBLOCK` rows of the `SOLID` layout hold eleven attribute fields
//!   (material, element type, real constant, section, coordinate system,
//!   birth/death, solid model reference, shape flag, **node count**, unused,
//!   **element number**) followed by the node numbers.  Elements with more
//!   nodes than fit in the row continue on further rows that hold node numbers
//!   only.
//! * `ET,<id>,<name>` commands tell which ANSYS element type each type number
//!   stands for, which decides whether a 4-node element is a tetrahedron or a
//!   shell.

use std::collections::HashMap;

use super::fortran_format::{FieldKind, FortranFormat, read_integer, read_real};
use super::types::{FeElement, FeElementType, FeNode};

/// Number of integer attribute fields in front of the nodes of a `SOLID`
/// `EBLOCK` row.
const EBLOCK_HEADER_FIELDS: usize = 11;
/// Position (0-based) of the material number in a `SOLID` `EBLOCK` row.
const EBLOCK_MATERIAL: usize = 0;
/// Position (0-based) of the element type number in a `SOLID` `EBLOCK` row.
const EBLOCK_TYPE: usize = 1;
/// Position (0-based) of the node count in a `SOLID` `EBLOCK` row.
const EBLOCK_NODE_COUNT: usize = 8;
/// Position (0-based) of the element number in a `SOLID` `EBLOCK` row.
const EBLOCK_ELEMENT_ID: usize = 10;

/// `ET` table: element type number of the deck to ANSYS element number
/// (`ET,1,185` and `ET,1,SOLID185` both map `1` to `185`).
pub(super) type ElementTypeTable = HashMap<usize, u32>;

/// Is `upper` (an upper-cased, trimmed line) the APDL command `name`?
///
/// APDL commands may be abbreviated to their first four letters.
pub(super) fn is_command(upper: &str, name: &str) -> bool {
    let word = upper
        .split(|c: char| c == ',' || c.is_whitespace())
        .next()
        .unwrap_or("");
    word.len() >= 4 && name.starts_with(word)
}

/// Record an `ET,<id>,<name>` command in `table`.
pub(super) fn record_element_type(line: &str, table: &mut ElementTypeTable) {
    let mut parts = line.split(',');
    parts.next();
    let type_id = parts
        .next()
        .and_then(|text| text.trim().parse::<usize>().ok());
    let number = parts.next().and_then(ansys_element_number);
    if let (Some(type_id), Some(number)) = (type_id, number) {
        table.insert(type_id, number);
    }
}

/// The number in an element name: `185` for `185` and for `SOLID185`.
fn ansys_element_number(name: &str) -> Option<u32> {
    let digits: String = name
        .trim()
        .chars()
        .skip_while(|c| !c.is_ascii_digit())
        .take_while(char::is_ascii_digit)
        .collect();
    digits.parse().ok()
}

/// Read the rows of an `NBLOCK`.
///
/// `first` is the index in `lines` of the line after the command.  Nodes are
/// appended to `nodes`, problems to `issues`; the return value is the index
/// of the first line after the block.
pub(super) fn read_nblock(
    lines: &[&str],
    first: usize,
    nodes: &mut Vec<FeNode>,
    issues: &mut Vec<String>,
) -> usize {
    let Some((format, start)) = read_format_line(lines, first, "NBLOCK", issues) else {
        return skip_block(lines, first);
    };
    let layout = match NodeLayout::new(format) {
        Ok(layout) => layout,
        Err(message) => {
            issues.push(format!("line {}: NBLOCK {message}", start));
            return skip_block(lines, start);
        }
    };
    let mut at = start;
    while let Some(&line) = lines.get(at) {
        match classify_row(line) {
            RowKind::Blank => {}
            RowKind::Terminator { consume } => {
                return if consume { at + 1 } else { at };
            }
            RowKind::Data => match layout.read_node(line) {
                Ok(node) => nodes.push(node),
                Err(message) => issues.push(format!("line {}: NBLOCK row: {message}", at + 1)),
            },
        }
        at += 1;
    }
    at
}

/// Read the rows of an `EBLOCK`.
///
/// `first` is the index in `lines` of the line after the command and
/// `command` the command line itself, whose third field is the layout key.
/// Only the `SOLID` layout is understood.  Elements are appended to
/// `elements`, problems to `issues`; the return value is the index of the first
/// line after the block.
pub(super) fn read_eblock(
    lines: &[&str],
    first: usize,
    command: &str,
    types: &ElementTypeTable,
    elements: &mut Vec<FeElement>,
    issues: &mut Vec<String>,
) -> usize {
    let key = command
        .split(',')
        .nth(2)
        .map(|text| text.trim().to_ascii_uppercase())
        .unwrap_or_default();
    let Some((format, start)) = read_format_line(lines, first, "EBLOCK", issues) else {
        return skip_block(lines, first);
    };
    if key != "SOLID" {
        issues.push(format!(
            "line {first}: EBLOCK layout `{key}` is not supported; only the SOLID layout is read"
        ));
        return skip_block(lines, start);
    }
    let fields = format.fields();
    if fields.len() <= EBLOCK_HEADER_FIELDS || fields.iter().any(|f| f.kind != FieldKind::Integer) {
        issues.push(format!(
            "line {}: an EBLOCK format must consist of more than {EBLOCK_HEADER_FIELDS} integer fields",
            start
        ));
        return skip_block(lines, start);
    }
    let mut at = start;
    while let Some(&line) = lines.get(at) {
        match classify_row(line) {
            RowKind::Blank => at += 1,
            RowKind::Terminator { consume } => {
                return if consume { at + 1 } else { at };
            }
            RowKind::Data => match read_solid_element(&format, lines, at, types) {
                Ok((element, next)) => {
                    elements.push(element);
                    at = next;
                }
                Err(message) => {
                    // The rows can no longer be told apart from continuation
                    // rows, so the rest of the block is abandoned.
                    issues.push(message);
                    return skip_block(lines, at + 1);
                }
            },
        }
    }
    at
}

/// What one line inside a block is.
enum RowKind {
    /// Nothing to read.
    Blank,
    /// The block ends here; `consume` tells whether the line belongs to the
    /// block (`-1`, `N,R5.3,LOC,-1`) or is the next command.
    Terminator { consume: bool },
    /// A row of the block.
    Data,
}

fn classify_row(line: &str) -> RowKind {
    let text = line.trim();
    if text.is_empty() || text.starts_with('!') {
        return RowKind::Blank;
    }
    if text.split_whitespace().next() == Some("-1")
        || text.to_ascii_uppercase().starts_with("N,R5.3,LOC")
    {
        return RowKind::Terminator { consume: true };
    }
    // Rows are numbers; anything that starts like an APDL command (`ET`, `/PREP7`,
    // `*SET`) is the next command of a block that lacks its terminator.
    if text.starts_with(|c: char| c.is_ascii_alphabetic() || c == '/' || c == '*') {
        return RowKind::Terminator { consume: false };
    }
    RowKind::Data
}

/// Index of the line after the block that starts at `from`.
fn skip_block(lines: &[&str], from: usize) -> usize {
    let mut at = from;
    while let Some(&line) = lines.get(at) {
        match classify_row(line) {
            RowKind::Terminator { consume } => return if consume { at + 1 } else { at },
            RowKind::Blank | RowKind::Data => at += 1,
        }
    }
    at
}

/// Parse the format line that has to follow a block command.
///
/// Returns the format and the index of the first row, or records an issue and
/// returns `None`.
fn read_format_line(
    lines: &[&str],
    first: usize,
    block: &str,
    issues: &mut Vec<String>,
) -> Option<(FortranFormat, usize)> {
    let mut at = first;
    while lines.get(at).is_some_and(|line| line.trim().is_empty()) {
        at += 1;
    }
    let Some(line) = lines.get(at) else {
        issues.push(format!(
            "line {first}: {block} is not followed by a format line"
        ));
        return None;
    };
    let text = line.trim();
    if !text.starts_with('(') {
        issues.push(format!(
            "line {}: {block} has to be followed by a Fortran format line such as `(3i9,6e21.13e3)`, found `{text}`",
            at + 1
        ));
        return None;
    }
    match FortranFormat::parse(text) {
        Ok(format) => Some((format, at + 1)),
        Err(message) => {
            issues.push(format!(
                "line {}: bad {block} format `{text}`: {message}",
                at + 1
            ));
            None
        }
    }
}

/// Where the node number and the coordinates sit in an `NBLOCK` row.
struct NodeLayout {
    format: FortranFormat,
    /// Field index of x, y and z; `None` when the format ends before it.
    coordinates: [Option<usize>; 3],
}

impl NodeLayout {
    fn new(format: FortranFormat) -> Result<Self, String> {
        if format.fields().first().map(|f| f.kind) != Some(FieldKind::Integer) {
            return Err("rows have to start with the integer node number".to_string());
        }
        let mut reals = format
            .fields()
            .iter()
            .enumerate()
            .filter(|(_, spec)| spec.kind == FieldKind::Real)
            .map(|(index, _)| index);
        let coordinates = [reals.next(), reals.next(), reals.next()];
        if coordinates[0].is_none() {
            return Err("format has no coordinate fields".to_string());
        }
        Ok(Self {
            format,
            coordinates,
        })
    }

    fn read_node(&self, line: &str) -> Result<FeNode, String> {
        let fields = self.format.split(line)?;
        let number = read_integer(fields.first().copied().unwrap_or(""))?;
        let id = usize::try_from(number)
            .ok()
            .filter(|id| *id > 0)
            .ok_or_else(|| format!("`{number}` is not a valid node number"))?;
        let mut coords = [0.0; 3];
        for (slot, index) in coords.iter_mut().zip(self.coordinates) {
            if let Some(index) = index {
                *slot = read_real(fields.get(index).copied().unwrap_or(""))?;
            }
        }
        Ok(FeNode::new(id, coords))
    }
}

fn row_integers(format: &FortranFormat, line: &str) -> Result<Vec<i64>, String> {
    format.split(line)?.into_iter().map(read_integer).collect()
}

/// Append the node numbers in `values` to `nodes`, up to `wanted` of them.
fn take_nodes(values: &[i64], wanted: usize, nodes: &mut Vec<usize>) -> Result<(), String> {
    for &value in values {
        if nodes.len() == wanted {
            break;
        }
        let id = usize::try_from(value)
            .ok()
            .filter(|id| *id > 0)
            .ok_or_else(|| format!("`{value}` is not a valid node number"))?;
        nodes.push(id);
    }
    Ok(())
}

/// Read the element that starts at row `at` of a `SOLID` `EBLOCK`, together
/// with its continuation rows; returns it and the index of the next row.
fn read_solid_element(
    format: &FortranFormat,
    lines: &[&str],
    at: usize,
    types: &ElementTypeTable,
) -> Result<(FeElement, usize), String> {
    let row_text = lines.get(at).copied().unwrap_or("");
    let header = row_integers(format, row_text)
        .map_err(|message| format!("line {}: EBLOCK row: {message}", at + 1))?;
    let field = |index: usize| header.get(index).copied().unwrap_or(0);

    let element_id = usize::try_from(field(EBLOCK_ELEMENT_ID))
        .ok()
        .filter(|id| *id > 0)
        .ok_or_else(|| {
            format!(
                "line {}: `{}` is not a valid element number",
                at + 1,
                field(EBLOCK_ELEMENT_ID)
            )
        })?;
    let wanted = usize::try_from(field(EBLOCK_NODE_COUNT))
        .ok()
        .filter(|count| *count > 0)
        .ok_or_else(|| {
            format!(
                "line {}: element {element_id} declares {} nodes",
                at + 1,
                field(EBLOCK_NODE_COUNT)
            )
        })?;

    let mut nodes = Vec::new();
    take_nodes(
        header.get(EBLOCK_HEADER_FIELDS..).unwrap_or(&[]),
        wanted,
        &mut nodes,
    )
    .map_err(|message| format!("line {}: EBLOCK row: {message}", at + 1))?;
    let mut next = at + 1;
    while nodes.len() < wanted {
        let row = lines
            .get(next)
            .copied()
            .filter(|row| matches!(classify_row(row), RowKind::Data));
        let Some(row) = row else {
            return Err(format!(
                "line {}: element {element_id} declares {wanted} nodes but only {} are given",
                at + 1,
                nodes.len()
            ));
        };
        let values = row_integers(format, row)
            .map_err(|message| format!("line {}: EBLOCK row: {message}", next + 1))?;
        take_nodes(&values, wanted, &mut nodes)
            .map_err(|message| format!("line {}: EBLOCK row: {message}", next + 1))?;
        next += 1;
    }

    let ansys_number = usize::try_from(field(EBLOCK_TYPE))
        .ok()
        .and_then(|type_id| types.get(&type_id).copied());
    let mut element = FeElement::new(
        element_id,
        classify_element(ansys_number, nodes.len()),
        nodes,
    );
    element.material_id = usize::try_from(field(EBLOCK_MATERIAL))
        .ok()
        .filter(|material| *material > 0);
    Ok((element, next))
}

/// Topology family of an ANSYS element type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Family {
    /// 3-D solid (tetrahedron, brick).
    Solid,
    /// Shell or 2-D solid (triangle, quadrilateral).
    Surface,
    /// Truss or beam.
    Line,
}

/// Family of the ANSYS element with the given number, when it is one of the
/// well-known structural and thermal element types listed here.
fn family_of(number: u32) -> Option<Family> {
    match number {
        // SOLID45/65/70/87/90/92/95/185/186/187/226/227/278/279/285.
        45 | 65 | 70 | 87 | 90 | 92 | 95 | 185 | 186 | 187 | 226 | 227 | 278 | 279 | 285 => {
            Some(Family::Solid)
        }
        // PLANE2/35/42/55/77/82/182/183 and SHELL41/43/63/91/93/99/143/181/281.
        2 | 35 | 41 | 42 | 43 | 55 | 63 | 77 | 82 | 91 | 93 | 99 | 143 | 181 | 182 | 183 | 281 => {
            Some(Family::Surface)
        }
        // LINK1/8/10/180, BEAM3/4/44/54/188 and PIPE16/288.
        1 | 3 | 4 | 8 | 10 | 16 | 44 | 54 | 180 | 188 | 288 => Some(Family::Line),
        _ => None,
    }
}

/// Element type for an element with `node_count` nodes.
///
/// `ansys_number` is the ANSYS element number the deck's `ET` table gives for
/// the element's type, when there is one.  Without an `ET` entry the element
/// is taken to be a solid, which is what the node count of a tetrahedron or
/// brick says; an element whose type is registered but not one of the known
/// families, such as a contact or mass element, stays `Unknown` instead of
/// being mistaken for a solid.
fn classify_element(ansys_number: Option<u32>, node_count: usize) -> FeElementType {
    let family = match ansys_number {
        None => Some(Family::Solid),
        Some(number) => family_of(number),
    };
    let known = match family {
        Some(Family::Solid) => match node_count {
            4 => Some(FeElementType::Tet4),
            8 => Some(FeElementType::Hex8),
            10 => Some(FeElementType::Tet10),
            20 => Some(FeElementType::Hex20),
            _ => None,
        },
        Some(Family::Surface) => match node_count {
            3 => Some(FeElementType::Tri3),
            4 => Some(FeElementType::Quad4),
            6 => Some(FeElementType::Tri6),
            8 => Some(FeElementType::Quad8),
            _ => None,
        },
        Some(Family::Line) => (node_count == 2).then_some(FeElementType::Line2),
        None => None,
    };
    known.unwrap_or_else(|| match ansys_number {
        Some(number) => FeElementType::Unknown(format!("mapdl:{number}:{node_count}")),
        None => FeElementType::Unknown(format!("mapdl:{node_count}-node")),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_element_io::types::MapdlDeck;

    /// Fortran `E21.13E3`: a 21-column field whose sign column is blank for
    /// positive values, so negative values touch the previous field.
    fn e21(value: f64) -> String {
        let magnitude = format!("{:.13E}", value.abs());
        let (mantissa, exponent) = magnitude.split_once('E').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let sign = if value.is_sign_negative() { '-' } else { ' ' };
        let exponent_sign = if exponent < 0 { '-' } else { '+' };
        format!("{sign}{mantissa}E{exponent_sign}{:03}", exponent.abs())
    }

    fn nblock_row(id: usize, xyz: [f64; 3]) -> String {
        format!(
            "{id:>9}{:>9}{:>9}{}{}{}",
            0,
            0,
            e21(xyz[0]),
            e21(xyz[1]),
            e21(xyz[2])
        )
    }

    fn eblock_row(values: &[i64]) -> String {
        values.iter().map(|value| format!("{value:>9}")).collect()
    }

    fn solid_row(material: i64, et: i64, nodes: &[i64], element: i64) -> String {
        let mut values = vec![
            material,
            et,
            1,
            1,
            0,
            0,
            0,
            0,
            nodes.len() as i64,
            0,
            element,
        ];
        values.extend(nodes.iter().take(8));
        eblock_row(&values)
    }

    fn deck(lines: &[String]) -> String {
        let mut text = lines.join("\n");
        text.push('\n');
        text
    }

    fn strs(lines: &[&str]) -> Vec<String> {
        lines.iter().map(|line| line.to_string()).collect()
    }

    #[test]
    fn e21_helper_matches_the_cdb_layout() {
        assert_eq!(e21(0.0), " 0.0000000000000E+000");
        assert_eq!(e21(-0.25), "-2.5000000000000E-001");
        assert_eq!(e21(1234.5), " 1.2345000000000E+003");
        assert_eq!(e21(0.0).len(), 21);
        assert_eq!(e21(-1.0e-100).len(), 21);
    }

    #[test]
    fn nblock_rows_are_read_by_column() {
        let mut lines = strs(&["NBLOCK,6,SOLID,         3,         3", "(3i9,6e21.13e3)"]);
        lines.push(nblock_row(1, [0.0, 0.0, 0.0]));
        lines.push(nblock_row(7, [1.5, -2.25, 3.0e-3]));
        lines.push(nblock_row(12, [-1.0, -2.0, -3.0]));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.nodes.len(), 3);
        assert_eq!(deck.nodes[0], FeNode::new(1, [0.0, 0.0, 0.0]));
        assert_eq!(deck.nodes[1], FeNode::new(7, [1.5, -2.25, 3.0e-3]));
        assert_eq!(deck.nodes[2], FeNode::new(12, [-1.0, -2.0, -3.0]));
    }

    #[test]
    fn literal_cdb_rows_with_touching_negative_values() {
        // Rows as `CDWRITE` prints them: the negative y and z coordinates fill
        // their 21 columns and abut the preceding field.
        let text = [
            "/PREP7",
            "NBLOCK,6,SOLID,         2,         2",
            "(3i9,6e21.13e3)",
            "        1        0        0 1.0000000000000E+000-2.0000000000000E+000-3.5000000000000E-002",
            "        2        0        0-4.0000000000000E+001 5.0000000000000E+000 6.0000000000000E+002",
            "-1",
        ]
        .join("\n");
        let deck = MapdlDeck::try_parse(&text).unwrap();
        assert_eq!(deck.nodes.len(), 2);
        assert_eq!(deck.nodes[0].id, 1);
        assert_eq!(deck.nodes[0].coords, [1.0, -2.0, -0.035]);
        assert_eq!(deck.nodes[1].id, 2);
        assert_eq!(deck.nodes[1].coords, [-40.0, 5.0, 600.0]);
    }

    #[test]
    fn old_style_narrow_nblock_format() {
        let text = [
            "NBLOCK,6,SOLID,         1,         1",
            "(3i8,6e16.9)",
            "       1       0       0 0.100000000E+01-0.200000000E+01 0.300000000E+01",
            "N,R5.3,LOC,       -1,",
        ]
        .join("\n");
        let deck = MapdlDeck::try_parse(&text).unwrap();
        assert_eq!(deck.nodes, vec![FeNode::new(1, [1.0, -2.0, 3.0])]);
    }

    #[test]
    fn short_nblock_rows_leave_trailing_coordinates_zero() {
        let mut lines = strs(&["NBLOCK,6,SOLID,         2,         2", "(3i9,6e21.13e3)"]);
        lines.push(format!("{:>9}{:>9}{:>9}{}{}", 1, 0, 0, e21(2.0), e21(-4.0)));
        lines.push(format!("{:>9}", 2));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.nodes[0].coords, [2.0, -4.0, 0.0]);
        assert_eq!(deck.nodes[1].coords, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn nblock_stops_at_either_terminator_and_commands_after_it_still_apply() {
        for terminator in ["-1", "N,R5.3,LOC,       -1,"] {
            let mut lines = strs(&["NBLOCK,6,SOLID,         1,         1", "(3i9,6e21.13e3)"]);
            lines.push(nblock_row(1, [1.0, 2.0, 3.0]));
            lines.push(terminator.to_string());
            lines.extend(strs(&["MP,EX,1,200E9", "D,1,UX,0.0", "F,1,FZ,50.0"]));
            let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
            assert_eq!(deck.nodes.len(), 1, "terminator {terminator}");
            assert_eq!(deck.materials.get(&1).map(|m| m.0), Some(200e9));
            assert_eq!(deck.disp_constraints.get(&1).map(Vec::len), Some(1));
            assert_eq!(deck.nodal_forces.get(&1).map(Vec::len), Some(1));
        }
    }

    #[test]
    fn nblock_without_a_terminator_ends_at_the_next_command() {
        let mut lines = strs(&["NBLOCK,6,SOLID,         1,         1", "(3i9,6e21.13e3)"]);
        lines.push(nblock_row(4, [1.0, 2.0, 3.0]));
        lines.push("MP,EX,2,70E9".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.nodes.len(), 1);
        assert_eq!(deck.materials.get(&2).map(|m| m.0), Some(70e9));
    }

    #[test]
    fn apdl_star_commands_end_a_block_that_lacks_its_terminator() {
        let mut lines = strs(&["NBLOCK,6,SOLID,         1,         1", "(3i9,6e21.13e3)"]);
        lines.push(nblock_row(1, [1.0, 2.0, 3.0]));
        lines.push("*SET,LOAD,5.0".to_string());
        lines.push("MP,EX,1,1E9".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.nodes.len(), 1);
        assert_eq!(deck.materials.get(&1).map(|m| m.0), Some(1e9));
    }

    #[test]
    fn eblock_reads_attribute_columns_and_the_element_number() {
        let mut lines = strs(&[
            "ET,       1,185",
            "EBLOCK,19,SOLID,         2,         2",
            "(19i9)",
        ]);
        lines.push(solid_row(3, 1, &[11, 12, 13, 14, 15, 16, 17, 18], 41));
        lines.push(solid_row(0, 1, &[12, 13, 19, 20, 16, 17, 21, 22], 42));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.elements.len(), 2);
        let first = &deck.elements[0];
        assert_eq!(first.id, 41);
        assert_eq!(first.element_type, FeElementType::Hex8);
        assert_eq!(first.connectivity, vec![11, 12, 13, 14, 15, 16, 17, 18]);
        assert_eq!(first.material_id, Some(3));
        let second = &deck.elements[1];
        assert_eq!(second.id, 42);
        assert_eq!(second.connectivity, vec![12, 13, 19, 20, 16, 17, 21, 22]);
        assert_eq!(second.material_id, None);
    }

    #[test]
    fn eblock_continuation_row_completes_a_ten_node_tetrahedron() {
        let mut lines = strs(&[
            "ET,       2,187",
            "EBLOCK,19,SOLID,         2,         2",
            "(19i9)",
        ]);
        lines.push(solid_row(1, 2, &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 5));
        lines.push(eblock_row(&[9, 10]));
        lines.push(solid_row(1, 2, &[2, 3, 4, 11, 12, 13, 14, 15, 16, 17], 6));
        lines.push(eblock_row(&[16, 17]));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.elements.len(), 2);
        assert_eq!(deck.elements[0].id, 5);
        assert_eq!(deck.elements[0].element_type, FeElementType::Tet10);
        assert_eq!(
            deck.elements[0].connectivity,
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        );
        assert_eq!(deck.elements[1].id, 6);
        assert_eq!(deck.elements[1].connectivity.len(), 10);
        assert_eq!(deck.elements[1].connectivity[8..], [16, 17]);
    }

    #[test]
    fn twenty_node_brick_spans_two_rows() {
        let mut lines = strs(&["EBLOCK,19,SOLID,         1,         1", "(19i9)"]);
        let nodes: Vec<i64> = (101..121).collect();
        lines.push(solid_row(1, 1, &nodes, 9));
        lines.push(eblock_row(&nodes[8..]));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.elements.len(), 1);
        assert_eq!(deck.elements[0].element_type, FeElementType::Hex20);
        assert_eq!(deck.elements[0].connectivity.len(), 20);
        assert_eq!(deck.elements[0].connectivity[19], 120);
    }

    #[test]
    fn solid_node_counts_without_an_et_table() {
        let mut lines = strs(&["EBLOCK,19,SOLID,         3,         3", "(19i9)"]);
        lines.push(solid_row(1, 1, &[1, 2, 3, 4], 1));
        lines.push(solid_row(1, 1, &[1, 2, 3, 4, 5, 6, 7, 8], 2));
        lines.push(solid_row(1, 1, &[1, 2, 3, 4, 5, 6], 3));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.elements[0].element_type, FeElementType::Tet4);
        assert_eq!(deck.elements[1].element_type, FeElementType::Hex8);
        assert_eq!(
            deck.elements[2].element_type,
            FeElementType::Unknown("mapdl:6-node".to_string())
        );
    }

    #[test]
    fn element_type_table_separates_shells_beams_and_contact_from_solids() {
        let mut lines = strs(&[
            "ET,       1,185",
            "ET,       2,SHELL181",
            "ET,       3,188",
            "ET,       4,174",
            "EBLOCK,19,SOLID,         4,         4",
            "(19i9)",
        ]);
        lines.push(solid_row(1, 1, &[1, 2, 3, 4, 5, 6, 7, 8], 1));
        lines.push(solid_row(1, 2, &[1, 2, 3, 4], 2));
        lines.push(solid_row(1, 3, &[5, 6], 3));
        lines.push(solid_row(1, 4, &[1, 2, 3, 4], 4));
        lines.push("-1".to_string());
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.elements[0].element_type, FeElementType::Hex8);
        assert_eq!(deck.elements[1].element_type, FeElementType::Quad4);
        assert_eq!(deck.elements[2].element_type, FeElementType::Line2);
        assert_eq!(
            deck.elements[3].element_type,
            FeElementType::Unknown("mapdl:174:4".to_string())
        );
    }

    #[test]
    fn eblock_with_wide_integer_fields() {
        let text = [
            "EBLOCK,19,SOLID,         1,         1".to_string(),
            "(19i10)".to_string(),
            [1i64, 1, 1, 1, 0, 0, 0, 0, 4, 0, 77, 100, 200, 300, 400]
                .iter()
                .map(|v| format!("{v:>10}"))
                .collect::<String>(),
            "-1".to_string(),
        ]
        .join("\n");
        let deck = MapdlDeck::try_parse(&text).unwrap();
        assert_eq!(deck.elements[0].id, 77);
        assert_eq!(deck.elements[0].connectivity, vec![100, 200, 300, 400]);
    }

    #[test]
    fn unsupported_eblock_layout_is_reported_and_skipped() {
        let mut lines = strs(&["EBLOCK,9,,         1,         1", "(9i7)"]);
        lines.push(eblock_row(&[1, 2, 3, 4, 5, 6, 7, 8, 9]));
        lines.push("-1".to_string());
        lines.push("MP,EX,1,1E9".to_string());
        let text = deck(&lines);
        let error = MapdlDeck::try_parse(&text).unwrap_err().to_string();
        assert!(
            error.contains("EBLOCK") && error.contains("SOLID"),
            "{error}"
        );
        let lenient = MapdlDeck::parse(&text);
        assert!(lenient.elements.is_empty());
        assert_eq!(lenient.materials.get(&1).map(|m| m.0), Some(1e9));
    }

    #[test]
    fn blocks_without_a_usable_format_line_are_reported() {
        for text in [
            "NBLOCK,6,SOLID,         1,         1\n        1        0        0\n-1\n",
            "NBLOCK,6,SOLID,         1,         1\n(3i9,6x\n-1\n",
            "EBLOCK,19,SOLID,         1,         1\n",
            "NBLOCK,6,SOLID,         1,         1\n(3i9)\n        1        0        0\n-1\n",
        ] {
            assert!(
                MapdlDeck::try_parse(text).is_err(),
                "must be reported: {text}"
            );
            assert!(MapdlDeck::parse(text).nodes.is_empty());
        }
    }

    #[test]
    fn malformed_rows_are_reported_and_skipped_by_the_lenient_parser() {
        let mut lines = strs(&["NBLOCK,6,SOLID,         3,         3", "(3i9,6e21.13e3)"]);
        lines.push(nblock_row(1, [1.0, 2.0, 3.0]));
        lines.push(format!(
            "{:>9}{:>9}{:>9}{:>21}{}{}",
            2,
            0,
            0,
            "*********************",
            e21(0.0),
            e21(0.0)
        ));
        lines.push(nblock_row(3, [4.0, 5.0, 6.0]));
        lines.push("-1".to_string());
        let text = deck(&lines);
        let error = MapdlDeck::try_parse(&text).unwrap_err().to_string();
        assert!(error.contains("NBLOCK row"), "{error}");
        let ids: Vec<usize> = MapdlDeck::parse(&text).nodes.iter().map(|n| n.id).collect();
        assert_eq!(ids, vec![1, 3]);
    }

    #[test]
    fn truncated_element_node_list_is_reported() {
        let mut lines = strs(&["EBLOCK,19,SOLID,         1,         1", "(19i9)"]);
        lines.push(solid_row(1, 1, &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 5));
        lines.push("-1".to_string());
        let error = MapdlDeck::try_parse(&deck(&lines)).unwrap_err().to_string();
        assert!(
            error.contains("element 5") && error.contains("10 nodes"),
            "{error}"
        );
    }

    #[test]
    fn element_type_commands_and_abbreviations() {
        let mut table = ElementTypeTable::new();
        record_element_type("ET,       1,185", &mut table);
        record_element_type("ET,2,SOLID187", &mut table);
        record_element_type("ET,x,185", &mut table);
        record_element_type("ET,3", &mut table);
        assert_eq!(table.get(&1), Some(&185));
        assert_eq!(table.get(&2), Some(&187));
        assert_eq!(table.len(), 2);
        assert!(is_command("NBLOCK,6,SOLID", "NBLOCK"));
        assert!(is_command("NBLO,6", "NBLOCK"));
        assert!(is_command("EBLOCK", "EBLOCK"));
        assert!(!is_command("NBLOCKS,6", "NBLOCK"));
        assert!(!is_command("N,1,0,0,0", "NBLOCK"));
        assert!(!is_command("EN,1,2,3,4", "EBLOCK"));
    }

    #[test]
    fn a_complete_deck_reads_blocks_and_commands_together() {
        let mut lines = strs(&[
            "/PREP7",
            "ET,       1,285",
            "MP,EX,1,2.0E11",
            "MP,NUXY,1,0.3",
            "NBLOCK,6,SOLID,         4,         4",
            "(3i9,6e21.13e3)",
        ]);
        lines.push(nblock_row(1, [0.0, 0.0, 0.0]));
        lines.push(nblock_row(2, [1.0, 0.0, 0.0]));
        lines.push(nblock_row(3, [0.0, -1.0, 0.0]));
        lines.push(nblock_row(4, [0.0, 0.0, 1.0]));
        lines.push("-1".to_string());
        lines.extend(strs(&["EBLOCK,19,SOLID,         1,         1", "(19i9)"]));
        lines.push(solid_row(1, 1, &[1, 2, 3, 4], 1));
        lines.push("-1".to_string());
        lines.extend(strs(&["D,2,UY,0.5", "F,4,FZ,-100.0", "FINISH"]));
        let deck = MapdlDeck::try_parse(&deck(&lines)).unwrap();
        assert_eq!(deck.nodes.len(), 4);
        assert_eq!(deck.elements.len(), 1);
        assert_eq!(deck.elements[0].element_type, FeElementType::Tet4);
        assert_eq!(deck.materials.get(&1).copied(), Some((2.0e11, 0.3, 0.0)));
        assert_eq!(deck.disp_constraints.get(&2), Some(&vec![(2, 0.5)]));
        assert_eq!(deck.nodal_forces.get(&4), Some(&vec![(3, -100.0)]));
    }
}
