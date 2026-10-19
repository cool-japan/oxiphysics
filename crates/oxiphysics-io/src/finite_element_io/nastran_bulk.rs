//! Card assembly for NASTRAN Bulk Data.
//!
//! A Bulk Data card is not a line: a card that needs more than one line of
//! fields continues on the following lines, and only the joined card can be
//! interpreted.  `CHEXA` is the typical case: its first line holds the
//! element and property ids and six grid points, `G7` and `G8` sit on a
//! continuation line, and the 20-node form needs a second one.
//!
//! Three layouts are understood.
//!
//! * **Small field, fixed format**: the card name in columns 1-8, then eight
//!   data fields of eight columns, then the continuation field in columns
//!   73-80.  A continuation line starts with a continuation marker (`+...`)
//!   or with a blank name field.
//! * **Large field, fixed format**: the card name ends in `*`, the data
//!   fields are sixteen columns wide and there are four per line.  Its
//!   continuation lines usually start with `*`; they have the field width of
//!   the card whatever marker they start with.
//! * **Free field**: comma separated; a continuation line starts with a
//!   marker, or with a comma when its name field is blank.
//!
//! Every layout is reduced to the same logical field sequence: field 0 is the
//! card name followed by the data fields of all the card's lines, with
//! missing trailing fields of a line filled up with empty ones so that the
//! fields of a continuation line keep their position.
//!
//! The accessors that follow the card assembly read the fields of such a
//! sequence: ids, reals in Fortran notation (`2.1+11` is `2.1E+11`) and
//! fields that may be left blank.

use std::fmt::Write as _;

use super::fortran_format::parse_fortran_real;
use super::types::{FeElementType, NastranElement, NastranSpc};

/// Data fields per line of a small-field card.
const SMALL_FIELD_STRIDE: usize = 8;
/// Data fields per line of a large-field card.
const LARGE_FIELD_STRIDE: usize = 4;
/// Width of a small data field in columns.
const SMALL_FIELD_WIDTH: usize = 8;
/// Width of a large data field in columns.
const LARGE_FIELD_WIDTH: usize = 16;
/// Column after which a fixed-format line carries only the continuation field.
const DATA_END_COLUMN: usize = 72;
/// Last column of a fixed-format line.
const LINE_END_COLUMN: usize = 80;

/// Remove the `$` comment and the trailing white space of a line.
///
/// Leading white space is kept: in the fixed format it distinguishes a
/// continuation line with a blank name field from a new card.
pub(super) fn strip_comment(line: &str) -> &str {
    line.split('$').next().unwrap_or("").trim_end()
}

/// A line that delimits the Bulk Data section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SectionMarker {
    /// `BEGIN BULK`: the Bulk Data section starts after this line.
    BeginBulk,
    /// `ENDDATA`: the Bulk Data section ends here.
    EndData,
}

/// Recognise `BEGIN BULK` and `ENDDATA` lines (case-insensitively).
pub(super) fn section_marker(line: &str) -> Option<SectionMarker> {
    let mut words = line.split_whitespace();
    let first = words.next()?;
    if first.eq_ignore_ascii_case("ENDDATA") {
        return Some(SectionMarker::EndData);
    }
    let begins_bulk = first.eq_ignore_ascii_case("BEGIN")
        && words
            .next()
            .is_some_and(|word| word.eq_ignore_ascii_case("BULK"));
    begins_bulk.then_some(SectionMarker::BeginBulk)
}

/// Does `line` continue the previous card instead of starting a new one?
fn is_continuation(line: &str) -> bool {
    match line.chars().next() {
        Some('+' | '*' | ',') => true,
        Some(first) if first.is_whitespace() => {
            let blanks = line.chars().take_while(|c| c.is_whitespace()).count();
            let starts_like_a_name = line
                .trim_start()
                .chars()
                .next()
                .is_some_and(|c| c.is_ascii_alphabetic());
            // A blank name field spans eight columns.  An indented card name
            // does not, and is a card of its own; so is a comma separated line
            // that starts with a name however deeply it is indented, because a
            // free-field continuation leaves its name out by starting with the
            // comma.
            !starts_like_a_name || (blanks >= SMALL_FIELD_WIDTH && !line.contains(','))
        }
        _ => false,
    }
}

/// Is the card name of the fixed-format `line` marked as large field?
fn is_large_field(line: &str) -> bool {
    line.chars()
        .take(SMALL_FIELD_WIDTH)
        .collect::<String>()
        .trim_end()
        .ends_with('*')
}

/// Cut a fixed-format line into its name field, its data fields and the
/// continuation field.
fn split_fixed(line: &str, large: bool) -> Vec<String> {
    let chars: Vec<char> = line.chars().collect();
    let field = |from: usize, to: usize| -> String {
        let to = to.min(chars.len());
        match chars.get(from..to) {
            Some(columns) => columns.iter().collect::<String>().trim().to_string(),
            None => String::new(),
        }
    };
    let width = if large {
        LARGE_FIELD_WIDTH
    } else {
        SMALL_FIELD_WIDTH
    };
    let mut fields = vec![field(0, SMALL_FIELD_WIDTH)];
    let mut column = SMALL_FIELD_WIDTH;
    while column < DATA_END_COLUMN {
        fields.push(field(column, column + width));
        column += width;
    }
    fields.push(field(DATA_END_COLUMN, LINE_END_COLUMN));
    fields
}

fn split_free(line: &str) -> Vec<String> {
    line.split(',').map(|f| f.trim().to_string()).collect()
}

/// The data fields of a continuation line, including its own continuation
/// field if it has one.
///
/// A fixed-format continuation line has the field width of the card it
/// continues (`large`), whatever character its marker starts with.
fn continuation_data(line: &str, large: bool) -> Vec<String> {
    if line.contains(',') {
        split_free(line).into_iter().skip(1).collect()
    } else {
        split_fixed(line, large).into_iter().skip(1).collect()
    }
}

/// A card whose lines are still being collected.
struct OpenCard {
    /// Card name followed by the data fields collected so far.
    fields: Vec<String>,
    /// Data fields per line: eight, or four for a large-field card.
    stride: usize,
    /// Number of lines of the card so far.
    lines: usize,
}

impl OpenCard {
    fn append(&mut self, data: Vec<String>) {
        // The previous line's data fields end at this position.  Whatever sits
        // behind it is that line's continuation marker, and a short line is
        // filled up so that this line's first field lands where it belongs.
        let target = 1 + self.stride * self.lines;
        self.fields.truncate(target);
        self.fields.resize(target, String::new());
        self.fields.extend(data);
        self.lines += 1;
    }
}

/// Collects the lines of a Bulk Data section into whole cards.
#[derive(Default)]
pub(super) struct CardJoiner {
    cards: Vec<Vec<String>>,
    open: Option<OpenCard>,
}

impl CardJoiner {
    pub(super) fn new() -> Self {
        Self::default()
    }

    /// Add the next line of the section.
    ///
    /// `line` must have its comment stripped (see [`strip_comment`]) and must
    /// not be blank, but its leading white space has to be intact.  A
    /// continuation line without a card in front of it is dropped.
    pub(super) fn push_line(&mut self, line: &str) {
        if is_continuation(line) {
            if let Some(card) = self.open.as_mut() {
                let large = card.stride == LARGE_FIELD_STRIDE;
                card.append(continuation_data(line, large));
            }
            return;
        }
        self.close();
        let fields = if line.contains(',') {
            split_free(line)
        } else {
            split_fixed(line, is_large_field(line))
        };
        let large = fields.first().is_some_and(|name| name.ends_with('*'));
        self.open = Some(OpenCard {
            fields,
            stride: if large {
                LARGE_FIELD_STRIDE
            } else {
                SMALL_FIELD_STRIDE
            },
            lines: 1,
        });
    }

    fn close(&mut self) {
        if let Some(mut card) = self.open.take() {
            while card.fields.last().is_some_and(String::is_empty) {
                card.fields.pop();
            }
            self.cards.push(card.fields);
        }
    }

    /// The finished cards: card name first, then the data fields in order,
    /// without trailing empty fields.
    pub(super) fn finish(mut self) -> Vec<Vec<String>> {
        self.close();
        self.cards
    }
}

/// The grid point ids of a `CTETRA` / `CHEXA` card: every field behind the
/// element and property ids.
///
/// `None` when a field is not a positive integer.  An empty field in the list
/// is an omitted mid-side node, which no element type here can represent.
pub(super) fn element_node_ids(fields: &[String]) -> Option<Vec<usize>> {
    fields
        .get(3..)?
        .iter()
        .map(|field| field.parse::<usize>().ok().filter(|id| *id > 0))
        .collect()
}

/// The property id of a `CTETRA` / `CHEXA` card, which defaults to the
/// element id when the field is blank.
pub(super) fn property_id(fields: &[String], element_id: usize) -> Option<usize> {
    match fields.get(2) {
        Some(field) if !field.is_empty() => field.parse().ok(),
        _ => Some(element_id),
    }
}

/// The field at `index` as an id: `None` when it is absent, blank or not a
/// non-negative integer.
pub(super) fn required_id(fields: &[String], index: usize) -> Option<usize> {
    fields.get(index)?.parse().ok()
}

/// The field at `index` as an id, `default` when it is absent or blank, and
/// `None` when it holds something that is not a non-negative integer.
pub(super) fn optional_id(fields: &[String], index: usize, default: usize) -> Option<usize> {
    match fields.get(index).map(String::as_str) {
        None | Some("") => Some(default),
        Some(text) => text.parse().ok(),
    }
}

/// The field at `index` as a real number, in Fortran notation (`2.1+11` is
/// `2.1E+11`): `None` when it is absent, blank or not a number.
pub(super) fn required_real(fields: &[String], index: usize) -> Option<f64> {
    parse_fortran_real(fields.get(index)?)
}

/// A real field that may be left blank: `Some(None)` when it is absent or
/// blank, `Some(Some(value))` for a number and `None` for anything else.
pub(super) fn blank_or_real(fields: &[String], index: usize) -> Option<Option<f64>> {
    match fields.get(index).map(String::as_str) {
        None | Some("") => Some(None),
        Some(text) => parse_fortran_real(text).map(Some),
    }
}

/// The field at `index` as a real number, `default` when it is absent or
/// blank, and `None` when it holds something that is not a number.
pub(super) fn optional_real(fields: &[String], index: usize, default: f64) -> Option<f64> {
    blank_or_real(fields, index).map(|value| value.unwrap_or(default))
}

/// Poisson's ratio of a `MAT1` card that gives too little to derive it.
const DEFAULT_POISSON_RATIO: f64 = 0.3;

/// Complete the isotropic constants `(E, G, nu)` of a `MAT1` card from the
/// ones it gives, using `G = E / (2 (1 + nu))`.
///
/// Whatever is blank is derived from the other two.  With a single constant
/// Poisson's ratio is [`DEFAULT_POISSON_RATIO`]; a card without `E` and `G`
/// has no stiffness and keeps zero for both.  `None` when the shear modulus
/// would have to be derived from a Poisson's ratio of `-1` or less.
pub(super) fn complete_moduli(
    young: Option<f64>,
    shear: Option<f64>,
    poisson: Option<f64>,
) -> Option<(f64, f64, f64)> {
    match (young, shear, poisson) {
        (Some(e), Some(g), Some(nu)) => Some((e, g, nu)),
        (Some(e), Some(g), None) => {
            let nu = if g == 0.0 {
                DEFAULT_POISSON_RATIO
            } else {
                e / (2.0 * g) - 1.0
            };
            Some((e, g, nu))
        }
        (Some(e), None, nu) => {
            let nu = nu.unwrap_or(DEFAULT_POISSON_RATIO);
            let denominator = 2.0 * (1.0 + nu);
            (denominator > 0.0).then(|| (e, e / denominator, nu))
        }
        (None, Some(g), nu) => {
            let nu = nu.unwrap_or(DEFAULT_POISSON_RATIO);
            Some((2.0 * g * (1.0 + nu), g, nu))
        }
        (None, None, nu) => Some((0.0, 0.0, nu.unwrap_or(DEFAULT_POISSON_RATIO))),
    }
}

/// The constraints of an `SPC` card: `SPC, SID, G1, C1, D1 [, G2, C2, D2]`.
///
/// A card holds one or two constraints; the second one is optional but has to
/// be complete when any of its fields is given.  `D` defaults to `0.0`.  An
/// invalid card yields no constraints at all rather than half of them.
pub(super) fn spc_constraints(fields: &[String]) -> Vec<NastranSpc> {
    let Some(sid) = required_id(fields, 1) else {
        return Vec::new();
    };
    let mut constraints = Vec::with_capacity(2);
    for first in [2usize, 5] {
        let given = (first..first + 3)
            .any(|index| fields.get(index).is_some_and(|field| !field.is_empty()));
        if !given && first == 5 {
            break;
        }
        let gid = required_id(fields, first);
        let components = fields
            .get(first + 1)
            .filter(|field| !field.is_empty())
            .cloned();
        let value = optional_real(fields, first + 2, 0.0);
        let (Some(gid), Some(components), Some(value)) = (gid, components, value) else {
            return Vec::new();
        };
        constraints.push(NastranSpc {
            sid,
            gid,
            components,
            value,
        });
    }
    constraints
}

/// Write `element` as a free-field `CTETRA` or `CHEXA` card.
///
/// The first line carries the element id, the property id and the first six
/// grid points; the remaining ones follow on continuation lines of up to
/// eight.  `markers` numbers the continuation markers so that every pair is
/// unique within the file.  An element that has no card here is written as a
/// comment.
pub(super) fn write_element(out: &mut String, element: &NastranElement, markers: &mut usize) {
    let (card, expected) = match element.element_type {
        FeElementType::Tet4 => ("CTETRA", 4),
        FeElementType::Tet10 => ("CTETRA", 10),
        FeElementType::Hex8 => ("CHEXA", 8),
        FeElementType::Hex20 => ("CHEXA", 20),
        _ => {
            let _ = writeln!(out, "$ unsupported element {}", element.eid);
            return;
        }
    };
    if element.nodes.len() != expected {
        let _ = writeln!(
            out,
            "$ unsupported element {}: {card} needs {expected} nodes, has {}",
            element.eid,
            element.nodes.len()
        );
        return;
    }
    let (first, rest) = element.nodes.split_at(6.min(element.nodes.len()));
    let mut groups: Vec<&[usize]> = vec![first];
    groups.extend(rest.chunks(SMALL_FIELD_STRIDE));
    let first_marker = *markers + 1;
    let line_count = groups.len();
    for (index, nodes) in groups.into_iter().enumerate() {
        let mut line = if index == 0 {
            format!("{card},{},{}", element.eid, element.pid)
        } else {
            format!("+{}", first_marker + index - 1)
        };
        for node in nodes {
            let _ = write!(line, ",{node}");
        }
        if index + 1 < line_count {
            let _ = write!(line, ",+{}", first_marker + index);
        }
        let _ = writeln!(out, "{line}");
    }
    *markers += line_count - 1;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_element_io::types::{NastranBulkData, NastranGrid};

    fn parse(cards: &[&str]) -> NastranBulkData {
        let mut text = String::from("BEGIN BULK\n");
        for card in cards {
            text.push_str(card);
            text.push('\n');
        }
        text.push_str("ENDDATA\n");
        NastranBulkData::parse(&text)
    }

    fn element(eid: usize, pid: usize, kind: FeElementType, nodes: &[usize]) -> NastranElement {
        NastranElement {
            eid,
            pid,
            element_type: kind,
            nodes: nodes.to_vec(),
        }
    }

    #[test]
    fn fixed_format_chexa_joins_its_continuation_line() {
        let data = parse(&[
            "CHEXA   1       1       1       2       3       4       5       6       +C1",
            "+C1     7       8",
        ]);
        assert_eq!(
            data.elements,
            vec![element(
                1,
                1,
                FeElementType::Hex8,
                &[1, 2, 3, 4, 5, 6, 7, 8]
            )]
        );
    }

    #[test]
    fn free_field_chexa_joins_its_continuation_line() {
        let data = parse(&["CHEXA,10,3,1,2,3,4,5,6,+C1", "+C1,7,8"]);
        assert_eq!(data.elements.len(), 1);
        assert_eq!(data.elements[0].eid, 10);
        assert_eq!(data.elements[0].pid, 3);
        assert_eq!(data.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn continuation_with_a_blank_name_field_keeps_its_columns() {
        // The continuation has no marker: its first eight columns are blank,
        // so the grid points start in column 9.  Trimming the line would
        // shift them into the name field.
        let data = parse(&[
            "CHEXA   2       1       11      12      13      14      15      16      +",
            "        17      18",
        ]);
        assert_eq!(data.elements.len(), 1);
        assert_eq!(data.elements[0].nodes, vec![11, 12, 13, 14, 15, 16, 17, 18]);
        let free = parse(&["CHEXA,3,1,1,2,3,4,5,6,", ",7,8"]);
        assert_eq!(free.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    fn joined(lines: &[&str]) -> Vec<Vec<String>> {
        let mut joiner = CardJoiner::new();
        for line in lines {
            joiner.push_line(line);
        }
        joiner.finish()
    }

    fn strings(fields: &[&str]) -> Vec<String> {
        fields.iter().map(|field| field.to_string()).collect()
    }

    #[test]
    fn a_short_parent_line_is_padded_to_eight_data_fields() {
        // Only three data fields on the first line: the continuation line's
        // fields are fields 9, 10, ... of the card, not fields 4, 5, ...
        let expected = vec![strings(&[
            "CARD", "1", "2", "3", "", "", "", "", "", "A", "B",
        ])];
        assert_eq!(joined(&["CARD,1,2,3", "+,A,B"]), expected);
        assert_eq!(
            joined(&["CARD    1       2       3", "+       A       B"]),
            expected
        );
    }

    #[test]
    fn the_continuation_marker_field_is_not_data() {
        let expected = vec![strings(&[
            "CARD", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        ])];
        assert_eq!(joined(&["CARD,1,2,3,4,5,6,7,8,+M1", "+M1,9,10"]), expected);
        assert_eq!(
            joined(&[
                "CARD    1       2       3       4       5       6       7       8       +M1",
                "+M1     9       10",
            ]),
            expected
        );
    }

    #[test]
    fn a_hole_left_by_a_short_parent_line_is_not_a_valid_element() {
        // G4..G6 are missing from the first line, so the node list has holes.
        assert!(
            parse(&["CTETRA,4,1,1,2,3", "+,4,5,6,7,8,9,10"])
                .elements
                .is_empty()
        );
    }

    #[test]
    fn a_chexa_line_with_ten_fields_returns_no_element_and_does_not_panic() {
        // name, EID, PID, G1..G6 and the continuation marker: ten fields.
        assert!(parse(&["CHEXA,1,1,1,2,3,4,5,6,+C1"]).elements.is_empty());
        assert!(
            parse(&["CHEXA   1       1       1       2       3       4       5       6       +C1"])
                .elements
                .is_empty()
        );
        // Shorter and longer cards fail the same way.
        assert!(parse(&["CHEXA,1,1,1,2,3"]).elements.is_empty());
        assert!(parse(&["CHEXA,1"]).elements.is_empty());
        assert!(parse(&["CHEXA"]).elements.is_empty());
        assert!(parse(&["CHEXA,1,1,1,2,3,4,5,6,7,8,9"]).elements.is_empty());
    }

    #[test]
    fn single_line_chexa_with_all_nodes_is_still_accepted() {
        let data = parse(&["CHEXA,5,2,1,2,3,4,5,6,7,8"]);
        assert_eq!(data.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn quadratic_elements_use_all_their_continuation_lines() {
        let tet10 = parse(&["CTETRA,1,1,1,2,3,4,5,6,+A", "+A,7,8,9,10"]);
        assert_eq!(tet10.elements[0].element_type, FeElementType::Tet10);
        assert_eq!(tet10.elements[0].nodes, (1..=10).collect::<Vec<_>>());

        let hex20 = parse(&[
            "CHEXA,2,1,1,2,3,4,5,6,+B",
            "+B,7,8,9,10,11,12,13,14,+C",
            "+C,15,16,17,18,19,20",
        ]);
        assert_eq!(hex20.elements[0].element_type, FeElementType::Hex20);
        assert_eq!(hex20.elements[0].nodes, (1..=20).collect::<Vec<_>>());
    }

    #[test]
    fn element_cards_with_missing_or_surplus_nodes_are_dropped() {
        // Six nodes would be neither a Tet4 nor a Tet10; a hole in the
        // mid-side nodes (G8 omitted) cannot be represented either.
        assert!(parse(&["CTETRA,1,1,1,2,3,4,5,6"]).elements.is_empty());
        assert!(
            parse(&["CTETRA,1,1,1,2,3,4,5,6,+", "+,7,,9,10"])
                .elements
                .is_empty()
        );
        assert!(parse(&["CTETRA,1,1,1,2,3"]).elements.is_empty());
        assert!(parse(&["CTETRA,1,1,1,2,x,4"]).elements.is_empty());
        assert!(parse(&["CTETRA,1,1,1,2,0,4"]).elements.is_empty());
        // 8 + 1 nodes is no hexahedron.
        assert!(
            parse(&["CHEXA,1,1,1,2,3,4,5,6,+", "+,7,8,9"])
                .elements
                .is_empty()
        );
    }

    #[test]
    fn blank_property_id_defaults_to_the_element_id() {
        let data = parse(&["CTETRA,7,,1,2,3,4"]);
        assert_eq!(data.elements[0].eid, 7);
        assert_eq!(data.elements[0].pid, 7);
    }

    #[test]
    fn shorthand_exponents_of_fixed_format_decks_are_read() {
        // `2.1+11` is 2.1E+11 and `-1.5-3` is -1.5E-3: the exponent letter is
        // left out to fit eight columns.
        let data = parse(&[
            "GRID    1               1.-3    -1.5-3  2.5+1",
            "MAT1    1       2.1+11          .3      7.8+3",
        ]);
        assert_eq!(data.grids.len(), 1);
        assert_eq!(data.grids[0].coords, [1.0e-3, -1.5e-3, 25.0]);
        assert_eq!(data.materials.len(), 1);
        let material = &data.materials[0];
        assert_eq!(material.e, 2.1e11);
        assert_eq!(material.nu, 0.3);
        assert_eq!(material.rho, 7800.0);
        assert!((material.g - 2.1e11 / 2.6).abs() < 1.0);
    }

    #[test]
    fn blank_grid_fields_take_their_defaults() {
        let data = parse(&["GRID,1,,1.5,-2.0", "GRID,2", "GRID,3,4,,,7.0,5"]);
        assert_eq!(data.grids.len(), 3);
        assert_eq!(data.grids[0].coords, [1.5, -2.0, 0.0]);
        assert_eq!(data.grids[0].cp, 0);
        assert_eq!(data.grids[1].coords, [0.0, 0.0, 0.0]);
        assert_eq!(data.grids[2].coords, [0.0, 0.0, 7.0]);
        assert_eq!((data.grids[2].cp, data.grids[2].cd), (4, 5));
    }

    #[test]
    fn cards_with_fields_that_are_not_numbers_are_skipped_not_zeroed() {
        let data = parse(&[
            "GRID,1,0,abc,0.0,0.0",
            "GRID,x,0,1.0,0.0,0.0",
            "GRID,2,0,1.0,0.0,nan",
            "GRID,3,0.5,1.0,0.0,0.0",
            "MAT1,1,200E9,77E9,0.3,heavy",
            "MAT1,2,stiff",
            "FORCE,1,5,0,x,0.0,0.0,1.0",
            "FORCE,1,5,0,100.0,0.0,up,1.0",
            "FORCE,1,five,0,100.0,0.0,0.0,1.0",
            "SPC,1,10,123,zzz",
            "SPC,1,10,123,0.0,5",
            "SPC,1,10,,0.0",
            "SPC,1,,123,0.0",
        ]);
        assert!(data.grids.is_empty(), "{:?}", data.grids);
        assert!(data.materials.is_empty(), "{:?}", data.materials);
        assert!(data.forces.is_empty(), "{:?}", data.forces);
        assert!(data.spcs.is_empty(), "{:?}", data.spcs);
    }

    #[test]
    fn mat1_constants_are_completed_from_the_ones_that_are_given() {
        let data = parse(&[
            "MAT1,1,200.0E9,,0.25",
            "MAT1,2,,80.0E9,0.25",
            "MAT1,3,200.0E9,80.0E9,",
            "MAT1,4,100.0E9",
            "MAT1,5,200.0E9,80.0E9,0.25,7850.0",
            "MAT1,6,,,0.2",
            // A shear modulus cannot be derived from nu = -1.
            "MAT1,7,100.0E9,,-1.0",
        ]);
        assert_eq!(data.materials.len(), 6);
        let constants = |index: usize| {
            let m = &data.materials[index];
            (m.e, m.g, m.nu)
        };
        assert_eq!(constants(0), (200.0e9, 80.0e9, 0.25));
        assert_eq!(constants(1), (200.0e9, 80.0e9, 0.25));
        assert_eq!(constants(2), (200.0e9, 80.0e9, 0.25));
        let (e, g, nu) = constants(3);
        assert_eq!((e, nu), (100.0e9, 0.3));
        assert!((g - 100.0e9 / 2.6).abs() < 1.0);
        assert_eq!(constants(4), (200.0e9, 80.0e9, 0.25));
        assert_eq!(data.materials[4].rho, 7850.0);
        assert_eq!(constants(5), (0.0, 0.0, 0.2));
    }

    #[test]
    fn spc_cards_carry_one_or_two_constraints() {
        let data = parse(&[
            "SPC     2       32      3       -2.6    1       3       2.9",
            "SPC,7,10,123456",
            "SPC,8,11,12,0.5",
        ]);
        let spcs = &data.spcs;
        assert_eq!(spcs.len(), 4);
        assert_eq!((spcs[0].sid, spcs[0].gid), (2, 32));
        assert_eq!(spcs[0].components, "3");
        assert_eq!(spcs[0].value, -2.6);
        assert_eq!((spcs[1].sid, spcs[1].gid), (2, 1));
        assert_eq!(spcs[1].components, "3");
        assert_eq!(spcs[1].value, 2.9);
        assert_eq!(spcs[2].components, "123456");
        assert_eq!(spcs[2].value, 0.0);
        assert_eq!((spcs[3].gid, spcs[3].value), (11, 0.5));
    }

    #[test]
    fn force_direction_components_may_be_blank() {
        let data = parse(&[
            "FORCE,2,5,,100.0,1.0,0.0,",
            "FORCE,3,6,,-2.5E+1,,,1.0",
            "FORCE   4       7               1.+2    0.      -1.     0.",
        ]);
        assert_eq!(data.forces.len(), 3);
        assert_eq!(data.forces[0].direction, [1.0, 0.0, 0.0]);
        assert_eq!(data.forces[1].scale, -25.0);
        assert_eq!(data.forces[1].direction, [0.0, 0.0, 1.0]);
        assert_eq!(data.forces[2].scale, 100.0);
        assert_eq!(data.forces[2].direction, [0.0, -1.0, 0.0]);
    }

    #[test]
    fn large_field_cards_are_split_at_sixteen_columns() {
        let first = format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}{:<8}",
            "CHEXA*", 21, 4, 1, 2, "*K1"
        );
        let second = format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}{:<8}",
            "*K1", 3, 4, 5, 6, "*K2"
        );
        let third = format!("{:<8}{:<16}{:<16}", "*K2", 7, 8);
        let data = parse(&[&first, &second, &third]);
        assert_eq!(data.elements.len(), 1);
        assert_eq!(data.elements[0].eid, 21);
        assert_eq!(data.elements[0].pid, 4);
        assert_eq!(data.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn continuation_lines_have_the_field_width_of_their_card() {
        // A small-field card may use `*` in its continuation marker, and a
        // large-field card `+`; the widths follow the card, not the marker.
        let small = parse(&[
            "CHEXA   1       1       1       2       3       4       5       6       *C1",
            "*C1     7       8",
        ]);
        assert_eq!(small.elements.len(), 1);
        assert_eq!(small.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
        let first = format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}{:<8}",
            "CHEXA*", 7, 1, 1, 2, "+L1"
        );
        let second = format!(
            "{:<8}{:<16}{:<16}{:<16}{:<16}{:<8}",
            "+L1", 3, 4, 5, 6, "+L2"
        );
        let third = format!("{:<8}{:<16}{:<16}", "+L2", 7, 8);
        let large = parse(&[&first, &second, &third]);
        assert_eq!(large.elements.len(), 1);
        assert_eq!(large.elements[0].eid, 7);
        assert_eq!(large.elements[0].nodes, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn continuations_attach_to_their_own_card_only() {
        let data = parse(&[
            "GRID,1,0,0.0,0.0,0.0",
            "CHEXA,1,1,1,2,3,4,5,6,+A",
            "+A,7,8",
            "GRID,2,0,1.0,0.0,0.0",
            "CTETRA,2,1,1,2,3,4",
            "MAT1,1,200E9,77E9,0.3,7800",
        ]);
        assert_eq!(data.grids.len(), 2);
        assert_eq!(data.grids[1].id, 2);
        assert_eq!(data.elements.len(), 2);
        assert_eq!(data.elements[0].element_type, FeElementType::Hex8);
        assert_eq!(data.elements[1].element_type, FeElementType::Tet4);
        assert_eq!(data.materials.len(), 1);
    }

    #[test]
    fn card_names_must_match_exactly() {
        let data = parse(&[
            "FORCE1,1,5,100.0,1.0,2.0",
            "GRIDX,1,0,0.0,0.0,0.0",
            "MAT10,1,200E9,77E9,0.3,7800",
            "SPC1,1,123,5",
            "FORCE,1,5,0,100.0,0.0,0.0,1.0",
        ]);
        assert!(data.grids.is_empty());
        assert!(data.materials.is_empty());
        assert!(data.spcs.is_empty());
        assert_eq!(data.forces.len(), 1);
        assert_eq!(data.forces[0].gid, 5);
    }

    #[test]
    fn comments_and_orphan_continuations_are_ignored() {
        let data = parse(&[
            "+ORPHAN,1,2,3",
            "$ a full line comment",
            "GRID,1,0,1.0,2.0,3.0 $ trailing comment",
            "",
            "CTETRA,1,1,1,2,3,4 $ another one",
        ]);
        assert_eq!(data.grids.len(), 1);
        assert_eq!(data.grids[0].coords, [1.0, 2.0, 3.0]);
        assert_eq!(data.elements.len(), 1);
    }

    #[test]
    fn only_the_bulk_data_section_is_read() {
        let text = "SOL 101\nCEND\nSPC = 1\nGRID,9,0,0.0,0.0,0.0\nBEGIN BULK\nGRID,1,0,1.0,0.0,0.0\nENDDATA\nGRID,2,0,2.0,0.0,0.0\n";
        let data = NastranBulkData::parse(text);
        assert_eq!(data.grids.len(), 1);
        assert_eq!(data.grids[0].id, 1);
    }

    #[test]
    fn indented_free_field_cards_are_not_continuations() {
        let data = parse(&["  GRID,1,0,1.0,2.0,3.0", "  GRID,2,0,4.0,5.0,6.0"]);
        assert_eq!(data.grids.len(), 2);
        // Indenting by a whole field width or more does not turn a card into a
        // continuation of the card before it.
        let deep = parse(&[
            "        GRID,1,0,1.0,2.0,3.0",
            "        GRID,2,0,4.0,5.0,6.0",
            "        CTETRA,1,1,1,2,3,4",
        ]);
        assert_eq!(deep.grids.len(), 2);
        assert_eq!(deep.grids[1].coords, [4.0, 5.0, 6.0]);
        assert_eq!(deep.elements.len(), 1);
    }

    #[test]
    fn write_emits_every_supported_element_and_reads_back_identically() {
        let mut data = NastranBulkData::new();
        data.grids.push(NastranGrid::new(1, [0.0, 0.0, 0.0]));
        data.elements = vec![
            element(1, 1, FeElementType::Tet4, &[1, 2, 3, 4]),
            element(2, 7, FeElementType::Hex8, &[1, 2, 3, 4, 5, 6, 7, 8]),
            element(3, 2, FeElementType::Tet10, &(11..=20).collect::<Vec<_>>()),
            element(4, 2, FeElementType::Hex20, &(21..=40).collect::<Vec<_>>()),
        ];
        let text = data.write();
        assert!(text.contains("CTETRA,1,1,1,2,3,4\n"), "{text}");
        assert!(
            text.contains("CHEXA,2,7,1,2,3,4,5,6,+1\n+1,7,8\n"),
            "{text}"
        );
        assert_eq!(text.matches("CTETRA").count(), 2);
        assert_eq!(text.matches("CHEXA").count(), 2);
        let back = NastranBulkData::parse(&text);
        assert_eq!(back.elements, data.elements);
        assert_eq!(back.grids.len(), 1);
    }

    #[test]
    fn write_marks_unsupported_and_inconsistent_elements() {
        let mut data = NastranBulkData::new();
        data.elements = vec![
            element(5, 1, FeElementType::Tri3, &[1, 2, 3]),
            element(6, 1, FeElementType::Tet4, &[1, 2, 3]),
            element(8, 1, FeElementType::Tet4, &[1, 2, 3, 4]),
        ];
        let text = data.write();
        assert!(text.contains("$ unsupported element 5\n"), "{text}");
        assert!(text.contains("$ unsupported element 6"), "{text}");
        let back = NastranBulkData::parse(&text);
        assert_eq!(back.elements.len(), 1);
        assert_eq!(back.elements[0].eid, 8);
    }

    #[test]
    fn written_continuation_markers_are_unique() {
        let mut data = NastranBulkData::new();
        data.elements = vec![
            element(1, 1, FeElementType::Hex8, &[1, 2, 3, 4, 5, 6, 7, 8]),
            element(2, 1, FeElementType::Hex20, &(1..=20).collect::<Vec<_>>()),
            element(3, 1, FeElementType::Hex8, &[9, 10, 11, 12, 13, 14, 15, 16]),
        ];
        let text = data.write();
        let mut leading: Vec<&str> = text
            .lines()
            .filter(|line| line.starts_with('+'))
            .map(|line| line.split(',').next().unwrap())
            .collect();
        let total = leading.len();
        leading.sort_unstable();
        leading.dedup();
        assert_eq!(total, 4);
        assert_eq!(leading.len(), total, "{text}");
    }

    #[test]
    fn section_markers_are_recognised_case_insensitively() {
        assert_eq!(section_marker("BEGIN BULK"), Some(SectionMarker::BeginBulk));
        assert_eq!(
            section_marker("begin   bulk"),
            Some(SectionMarker::BeginBulk)
        );
        assert_eq!(section_marker("  ENDDATA"), Some(SectionMarker::EndData));
        assert_eq!(section_marker("enddata"), Some(SectionMarker::EndData));
        assert_eq!(section_marker("BEGIN SUPER=1"), None);
        assert_eq!(section_marker("GRID,1"), None);
        assert_eq!(section_marker(""), None);
    }

    #[test]
    fn line_classification() {
        assert!(is_continuation("+C1     7       8"));
        assert!(is_continuation("*K1"));
        assert!(is_continuation(",7,8"));
        assert!(is_continuation("        7       8"));
        assert!(!is_continuation("CHEXA,1"));
        assert!(!is_continuation("  GRID,1"));
        assert!(is_continuation("  7,8"));
        assert!(is_continuation("        ABC     1.0"));
        assert!(!is_continuation("        GRID,1,0,1.0,2.0,3.0"));
        assert!(is_continuation("        ,7,8"));
        assert!(is_large_field("CHEXA*  rest"));
        assert!(!is_large_field("CHEXA   rest"));
        assert_eq!(strip_comment("  GRID,1 $ x"), "  GRID,1");
        assert_eq!(strip_comment("$"), "");
    }
}
