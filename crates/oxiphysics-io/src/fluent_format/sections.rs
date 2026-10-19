// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The sections of an ASCII Fluent `.msh` file.
//!
//! | Section | Content | Handling |
//! |---------|---------|----------|
//! | `(0 …)`, `(1 …)` | comment, title | skipped |
//! | `(2 ND)` | dimension | kept as the default coordinate count of node zones |
//! | `(10 (zone first last type [nd])(x y z …))` | nodes | read |
//! | `(12 (zone first last type element-type)[(types …)])` | cells | read |
//! | `(13 (zone first last bc-type face-type)(faces …))` | faces | read |
//! | `(39 …)`, `(45 (zone-id type name)(…))` | zone names and types | read |
//! | `(2000 …)` to `(3999 …)` | binary forms of the above | rejected |
//! | any other | periodicity, trees, … | skipped |
//!
//! All numbers in the headers and in the face and cell data are hexadecimal,
//! the coordinates are decimal.  A section whose zone number is `0` only
//! declares the total count of its kind and carries no data.
//!
//! The data of a section is a list of numbers, not of lines: it is read token
//! by token, so a writer may break the lines wherever it likes.  A face row is
//! `n0 n1 … c0 c1`: the nodes of the face followed by the cell on either side
//! (`0` for none).  In a zone of face type `0` (mixed) or `5` (polygonal) each
//! row starts with its node count; in a zone of face type `2`, `3` or `4` every
//! face has that many nodes and the row has no prefix.

use super::scan::{Scanner, Tokens};
use super::{
    FluentCell, FluentCellType, FluentFace, FluentFaceType, FluentMesh, FluentNode, FluentZoneType,
};
use crate::Error;

const SECTION_DIMENSION: u32 = 2;
const SECTION_NODES: u32 = 10;
const SECTION_CELLS: u32 = 12;
const SECTION_FACES: u32 = 13;
const SECTION_ZONE_ATTRIBUTES: u32 = 39;
const SECTION_ZONE: u32 = 45;

/// Sections 2000 to 3999 are the single and double precision binary forms of
/// the ASCII sections; their data is raw bytes that no parenthesis counting
/// can skip.
const BINARY_SECTIONS: std::ops::RangeInclusive<u32> = 2000..=3999;

/// Most cells a mesh may declare in uniform cell zones.
///
/// A zone of uniform cells is described by its index range alone, so the
/// header, not the data, would decide how much memory the reader asks for.
const MAX_UNIFORM_CELLS: usize = 1 << 27;

/// Parse the bytes of an ASCII Fluent mesh file.
///
/// # Errors
///
/// Returns [`Error::Parse`], naming the line, for binary sections, for
/// sections that are not closed, for malformed numbers, for data that does not
/// match the index range of its header, and for face and cell types that
/// [`FluentMesh`] cannot represent.
pub(super) fn parse_mesh(bytes: &[u8]) -> Result<FluentMesh, Error> {
    let mut reader = SectionReader {
        scanner: Scanner::new(bytes),
        mesh: FluentMesh::new(),
        dimension: None,
    };
    reader.read_all()?;
    Ok(reader.mesh)
}

/// The `zone first last` triple that heads node, cell and face sections.
#[derive(Debug, Clone, Copy)]
struct IndexRange {
    zone: usize,
    first: usize,
    last: usize,
}

impl IndexRange {
    /// Number of indices; `last >= first` is checked when the range is made.
    fn count(&self) -> usize {
        self.last - self.first + 1
    }
}

struct SectionReader<'a> {
    scanner: Scanner<'a>,
    mesh: FluentMesh,
    /// Coordinates per node declared by a `(2 ND)` section.
    dimension: Option<usize>,
}

impl<'a> SectionReader<'a> {
    fn read_all(&mut self) -> Result<(), Error> {
        loop {
            self.scanner.skip_whitespace();
            let start = self.scanner.position();
            match self.scanner.peek() {
                None => return Ok(()),
                Some(b'(') => self.scanner.bump(),
                Some(other) => {
                    return Err(self.scanner.error_at(
                        start,
                        format!(
                            "expected '(' at the start of a section, found {:?}",
                            char::from(other)
                        ),
                    ));
                }
            }
            self.scanner.skip_whitespace();
            let Some(number) = self.scanner.read_decimal() else {
                return Err(self
                    .scanner
                    .error_at(start, "a section has to start with its number"));
            };
            self.read_section(number, start)?;
        }
    }

    fn read_section(&mut self, number: u32, start: usize) -> Result<(), Error> {
        match number {
            SECTION_DIMENSION => self.dimension_section(start),
            SECTION_NODES => self.node_section(start),
            SECTION_CELLS => self.cell_section(start),
            SECTION_FACES => self.face_section(start),
            SECTION_ZONE_ATTRIBUTES | SECTION_ZONE => self.zone_section(),
            n if BINARY_SECTIONS.contains(&n) => Err(self.scanner.error_at(
                start,
                format!(
                    "section {n} holds binary data; only ASCII Fluent meshes can be read \
                     (write the mesh with the ASCII option)"
                ),
            )),
            _ => self.scanner.skip_to_section_end(),
        }
    }

    fn dimension_section(&mut self, start: usize) -> Result<(), Error> {
        self.scanner.skip_whitespace();
        let at = self.scanner.position();
        match self.scanner.read_decimal() {
            Some(dimension @ (2 | 3)) => {
                self.dimension = usize::try_from(dimension).ok();
            }
            _ => {
                return Err(self
                    .scanner
                    .error_at(at, "the dimension section has to say 2 or 3"));
            }
        }
        self.scanner.expect_section_end(start)
    }

    /// Read the `(zone first last type …)` header group.
    fn read_header(&mut self, what: &str, min: usize, max: usize) -> Result<Vec<usize>, Error> {
        self.scanner.skip_whitespace();
        let at = self.scanner.position();
        if self.scanner.peek() != Some(b'(') {
            return Err(self.scanner.error_at(
                at,
                format!("the {what} section has no (zone first last …) header"),
            ));
        }
        let text = self.scanner.read_group()?;
        let mut values = Vec::new();
        for token in text.split_ascii_whitespace() {
            let value = parse_hex(token).map_err(|message| {
                self.scanner
                    .error_at(at, format!("{what} section header: {message}"))
            })?;
            values.push(value);
        }
        if !(min..=max).contains(&values.len()) {
            return Err(self.scanner.error_at(
                at,
                format!(
                    "the {what} section header has {} fields; {min} to {max} are expected",
                    values.len()
                ),
            ));
        }
        Ok(values)
    }

    /// The index range of a header; `None` for the declaration section of
    /// zone `0`, which has no data.
    fn index_range(
        &self,
        what: &str,
        header: &[usize],
        start: usize,
    ) -> Result<Option<IndexRange>, Error> {
        let &[zone, first, last, ..] = header else {
            return Err(self
                .scanner
                .error_at(start, format!("the {what} section header is too short")));
        };
        if zone == 0 {
            return Ok(None);
        }
        if first == 0 || last < first {
            return Err(self.scanner.error_at(
                start,
                format!(
                    "{what} zone {zone:x} has the index range {first:x} to {last:x}; indices \
                     start at 1 and the last one cannot precede the first"
                ),
            ));
        }
        Ok(Some(IndexRange { zone, first, last }))
    }

    /// Read the `( … )` data group that follows a header.
    fn data_group(
        &mut self,
        what: &str,
        zone: usize,
        start: usize,
    ) -> Result<(usize, &'a str), Error> {
        self.scanner.skip_whitespace();
        let at = self.scanner.position();
        if self.scanner.peek() != Some(b'(') {
            return Err(self.scanner.error_at(
                at,
                format!(
                    "{what} zone {zone:x} has no data block; the section that starts on line {} \
                     ends without one",
                    self.scanner.line_at(start)
                ),
            ));
        }
        let text = self.scanner.read_group()?;
        Ok((at + 1, text))
    }

    /// The next token of `tokens` as a hexadecimal number, with its offset.
    fn next_hex(
        &self,
        tokens: &mut Tokens<'_>,
        missing: impl FnOnce() -> String,
    ) -> Result<(usize, usize), Error> {
        match tokens.next_token() {
            Some((offset, token)) => parse_hex(token)
                .map(|value| (offset, value))
                .map_err(|message| self.scanner.error_at(offset, message)),
            None => Err(self.scanner.error_at(tokens.end_offset(), missing())),
        }
    }

    fn node_section(&mut self, start: usize) -> Result<(), Error> {
        let header = self.read_header("node", 4, 5)?;
        let Some(range) = self.index_range("node", &header, start)? else {
            return self.scanner.expect_section_end(start);
        };
        let (base, data) = self.data_group("node", range.zone, start)?;
        let count = range.count();
        let nd = self.node_dimension(header.get(4).copied(), count, data, base, start)?;
        let mut tokens = Tokens::new(data, base);
        for k in 0..count {
            let mut coords = [0.0_f64; 3];
            for slot in coords.iter_mut().take(nd) {
                let Some((offset, token)) = tokens.next_token() else {
                    return Err(self.scanner.error_at(
                        tokens.end_offset(),
                        format!(
                            "node zone {:x} declares {count} nodes but its data ends after {k}",
                            range.zone
                        ),
                    ));
                };
                *slot = parse_coordinate(token)
                    .map_err(|message| self.scanner.error_at(offset, message))?;
            }
            self.mesh
                .nodes
                .push(FluentNode::new(range.first + k, coords));
        }
        if let Some((offset, _)) = tokens.next_token() {
            return Err(self.scanner.error_at(
                offset,
                format!(
                    "node zone {:x} holds more coordinates than the {count} nodes it declares",
                    range.zone
                ),
            ));
        }
        self.scanner.expect_section_end(start)
    }

    /// Coordinates per node: the header's `nd`, else the `(2 ND)` section,
    /// else whatever count divides the data evenly.
    fn node_dimension(
        &self,
        header_nd: Option<usize>,
        count: usize,
        data: &str,
        base: usize,
        start: usize,
    ) -> Result<usize, Error> {
        let declared = header_nd.or(self.dimension);
        if let Some(nd) = declared {
            return if nd == 2 || nd == 3 {
                Ok(nd)
            } else {
                Err(self.scanner.error_at(
                    start,
                    format!("a node zone cannot have {nd} coordinates per node; 2 or 3 are valid"),
                ))
            };
        }
        let tokens = Tokens::new(data, base).count_remaining();
        [3usize, 2]
            .into_iter()
            .find(|nd| count.checked_mul(*nd) == Some(tokens))
            .ok_or_else(|| {
                self.scanner.error_at(
                    start,
                    format!(
                        "the node data holds {tokens} numbers, which is neither 2 nor 3 per node \
                         for {count} nodes, and neither the header nor a dimension section \
                         says which"
                    ),
                )
            })
    }

    fn cell_section(&mut self, start: usize) -> Result<(), Error> {
        let header = self.read_header("cell", 4, 5)?;
        let Some(range) = self.index_range("cell", &header, start)? else {
            return self.scanner.expect_section_end(start);
        };
        let Some(&element_type) = header.get(4) else {
            return Err(self.scanner.error_at(
                start,
                format!(
                    "cell zone {:x} has no element type in its header",
                    range.zone
                ),
            ));
        };
        if element_type == 0 {
            // Mixed zone: the data lists the element type of every cell.
            let (base, data) = self.data_group("cell", range.zone, start)?;
            let mut tokens = Tokens::new(data, base);
            for k in 0..range.count() {
                let (offset, code) = self.next_hex(&mut tokens, || {
                    format!(
                        "cell zone {:x} is mixed but lists the types of only {k} cells",
                        range.zone
                    )
                })?;
                let cell_type = self.cell_type(code, offset)?;
                self.mesh
                    .cells
                    .push(FluentCell::new(range.first + k, cell_type, range.zone));
            }
            if let Some((offset, _)) = tokens.next_token() {
                return Err(self.scanner.error_at(
                    offset,
                    format!(
                        "cell zone {:x} lists more cell types than the {} cells it declares",
                        range.zone,
                        range.count()
                    ),
                ));
            }
        } else {
            let cell_type = self.cell_type(element_type, start)?;
            let count = range.count();
            if self.mesh.cells.len().saturating_add(count) > MAX_UNIFORM_CELLS {
                return Err(self.scanner.error_at(
                    start,
                    format!(
                        "cell zone {:x} declares {count} cells, which takes the mesh past the \
                         limit of {MAX_UNIFORM_CELLS}",
                        range.zone
                    ),
                ));
            }
            self.mesh.cells.try_reserve(count).map_err(|_| {
                self.scanner.error_at(
                    start,
                    format!("cell zone {:x}: cannot allocate {count} cells", range.zone),
                )
            })?;
            for id in range.first..=range.last {
                self.mesh
                    .cells
                    .push(FluentCell::new(id, cell_type.clone(), range.zone));
            }
        }
        self.scanner.expect_section_end(start)
    }

    fn cell_type(&self, code: usize, offset: usize) -> Result<FluentCellType, Error> {
        u32::try_from(code)
            .ok()
            .and_then(FluentCellType::try_from_code)
            .ok_or_else(|| {
                let reason = if code == 7 {
                    "polyhedral cells (element type 7) are not supported".to_string()
                } else {
                    format!("{code} is not a known cell element type")
                };
                self.scanner.error_at(offset, reason)
            })
    }

    fn face_section(&mut self, start: usize) -> Result<(), Error> {
        let header = self.read_header("face", 4, 5)?;
        let Some(range) = self.index_range("face", &header, start)? else {
            return self.scanner.expect_section_end(start);
        };
        let Some(&face_type) = header.get(4) else {
            return Err(self.scanner.error_at(
                start,
                format!("face zone {:x} has no face type in its header", range.zone),
            ));
        };
        // Face type 0 (mixed) and 5 (polygonal) put the node count in front of
        // every face; 2, 3 and 4 fix it for the whole zone.
        let fixed_nodes = match face_type {
            0 | 5 => None,
            2..=4 => Some(face_type),
            other => {
                return Err(self.scanner.error_at(
                    start,
                    format!(
                        "face type {other} is not supported; expected 0 (mixed), 2 (line), \
                         3 (triangle), 4 (quadrilateral) or 5 (polygon)"
                    ),
                ));
            }
        };
        let (base, data) = self.data_group("face", range.zone, start)?;
        let mut tokens = Tokens::new(data, base);
        let zone = range.zone;
        for id in range.first..=range.last {
            let missing = || format!("face zone {zone:x} ends before face {id}");
            let (count_offset, node_count) = match fixed_nodes {
                Some(count) => (base, count),
                None => self.next_hex(&mut tokens, missing)?,
            };
            let face_type = FluentFaceType::from_node_count(node_count).ok_or_else(|| {
                self.scanner.error_at(
                    count_offset,
                    format!(
                        "face {id} has {node_count} nodes; only faces with 2 to 4 nodes (lines, \
                         triangles and quadrilaterals) can be represented"
                    ),
                )
            })?;
            let mut node_ids = Vec::with_capacity(node_count);
            for _ in 0..node_count {
                node_ids.push(self.next_hex(&mut tokens, missing)?.1);
            }
            let left_cell = self.next_hex(&mut tokens, missing)?.1;
            let right_cell = self.next_hex(&mut tokens, missing)?.1;
            self.mesh.faces.push(FluentFace::new(
                id, face_type, node_ids, left_cell, right_cell,
            ));
        }
        if let Some((offset, _)) = tokens.next_token() {
            return Err(self.scanner.error_at(
                offset,
                format!(
                    "face zone {zone:x} holds more data than the {} faces it declares",
                    range.count()
                ),
            ));
        }
        self.scanner.expect_section_end(start)
    }

    /// `(45 (zone-id zone-type zone-name)(…))`, or the same under number 39.
    fn zone_section(&mut self) -> Result<(), Error> {
        self.scanner.skip_whitespace();
        let at = self.scanner.position();
        if self.scanner.peek() != Some(b'(') {
            return Err(self.scanner.error_at(
                at,
                "a zone section has no (zone-id zone-type zone-name) header",
            ));
        }
        let text = self.scanner.read_group()?;
        let mut words = text.split_ascii_whitespace();
        let (Some(id_text), Some(type_text)) = (words.next(), words.next()) else {
            return Err(self
                .scanner
                .error_at(at, "a zone header needs a zone id and a zone type"));
        };
        // The id is decimal in these sections; hexadecimal is accepted too for
        // writers that use the notation of the other sections.
        let zone_id = id_text
            .parse::<usize>()
            .or_else(|_| usize::from_str_radix(id_text, 16))
            .map_err(|_| {
                self.scanner
                    .error_at(at, format!("{id_text:?} is not a zone id"))
            })?;
        if self.mesh.zone_type(zone_id).is_none() {
            self.mesh
                .add_zone(zone_id, FluentZoneType::from_keyword(type_text));
        }
        // The data group holds zone attributes, which are not interpreted.
        self.scanner.skip_to_section_end()
    }
}

fn parse_hex(token: &str) -> Result<usize, String> {
    usize::from_str_radix(token, 16).map_err(|_| format!("{token:?} is not a hexadecimal number"))
}

fn parse_coordinate(token: &str) -> Result<f64, String> {
    token
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| format!("{token:?} is not a finite number"))
}
