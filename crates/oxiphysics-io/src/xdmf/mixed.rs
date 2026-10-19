// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! XDMF `Mixed` topology streams.
//!
//! A `Mixed` topology has no fixed element size, so its connectivity is one
//! flat stream in which every element is introduced by its XDMF cell-type
//! code:
//!
//! ```text
//! 4 0 1 2        Triangle      code, 3 node ids
//! 5 0 1 2 3      Quadrilateral code, 4 node ids
//! 2 3 7 8 9      Polyline      code, node count, node ids
//! ```
//!
//! The fixed-size cell types are `Triangle = 4`, `Quadrilateral = 5`,
//! `Tetrahedron = 6`, `Pyramid = 7`, `Wedge = 8` and `Hexahedron = 9`.
//! `Polyvertex = 1`, `Polyline = 2` and `Polygon = 3` have a variable size, so
//! the code is followed by an explicit node count. Any other code is rejected.
//!
//! The helpers here walk such a stream, so the writers can count the elements
//! (the `NumberOfElements` attribute) and lay the stream out one element per
//! line instead of chunking it by a node count that does not exist.

use std::io;

const POLYVERTEX: usize = 1;
const POLYLINE: usize = 2;
const POLYGON: usize = 3;
const TRIANGLE: usize = 4;
const QUADRILATERAL: usize = 5;
const TETRAHEDRON: usize = 6;
const PYRAMID: usize = 7;
const WEDGE: usize = 8;
const HEXAHEDRON: usize = 9;

fn invalid_input(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message)
}

/// Length in stream entries (cell-type code, optional node count and node
/// ids) of the element that starts at `offset`.
fn element_len(stream: &[usize], offset: usize) -> io::Result<usize> {
    let overrun = |code: usize| {
        invalid_input(format!(
            "XDMF mixed topology: element of cell type {code} at connectivity offset {offset} runs past the end of the {}-entry stream",
            stream.len()
        ))
    };
    let code = *stream.get(offset).ok_or_else(|| {
        invalid_input(format!(
            "XDMF mixed topology: offset {offset} is past the end of the stream"
        ))
    })?;
    let (header, nodes): (usize, usize) = match code {
        POLYVERTEX | POLYLINE | POLYGON => {
            let count = *stream.get(offset + 1).ok_or_else(|| overrun(code))?;
            let minimum = match code {
                POLYVERTEX => 1,
                POLYLINE => 2,
                _ => 3,
            };
            if count < minimum {
                return Err(invalid_input(format!(
                    "XDMF mixed topology: cell type {code} at connectivity offset {offset} needs at least {minimum} nodes, found {count}"
                )));
            }
            (2, count)
        }
        TRIANGLE => (1, 3),
        QUADRILATERAL | TETRAHEDRON => (1, 4),
        PYRAMID => (1, 5),
        WEDGE => (1, 6),
        HEXAHEDRON => (1, 8),
        other => {
            return Err(invalid_input(format!(
                "XDMF mixed topology: unknown cell type code {other} at connectivity offset {offset}"
            )));
        }
    };
    let len = header.checked_add(nodes).ok_or_else(|| overrun(code))?;
    match offset.checked_add(len) {
        Some(end) if end <= stream.len() => Ok(len),
        _ => Err(overrun(code)),
    }
}

/// Iterator over the elements of a mixed stream; each item is the element's
/// slice (cell-type code, optional node count, node ids). After the first
/// malformed element it yields one `Err` and stops.
struct MixedElements<'a> {
    stream: &'a [usize],
    offset: usize,
}

impl<'a> Iterator for MixedElements<'a> {
    type Item = io::Result<&'a [usize]>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.offset >= self.stream.len() {
            return None;
        }
        match element_len(self.stream, self.offset) {
            Ok(len) => {
                let entry = self.stream.get(self.offset..self.offset + len)?;
                self.offset += len;
                Some(Ok(entry))
            }
            Err(error) => {
                self.offset = self.stream.len();
                Some(Err(error))
            }
        }
    }
}

fn elements(stream: &[usize]) -> MixedElements<'_> {
    MixedElements { stream, offset: 0 }
}

/// Number of elements in an XDMF `Mixed` connectivity stream.
///
/// Walks the cell-type codes of the stream (see the [module
/// documentation](self) for the layout). This is the value of the
/// `NumberOfElements` attribute of a `Mixed` topology.
///
/// Returns an [`io::ErrorKind::InvalidInput`] error for an unknown cell-type
/// code, a polytype with too few nodes, or an element that runs past the end
/// of the stream.
pub fn xdmf_mixed_element_count(connectivity: &[usize]) -> io::Result<usize> {
    let mut count = 0usize;
    for element in elements(connectivity) {
        element?;
        count += 1;
    }
    Ok(count)
}

/// One text row per element of a mixed stream: the cell-type code (and node
/// count) followed by the node ids, separated by spaces. Validates the stream
/// like [`xdmf_mixed_element_count`].
pub(super) fn mixed_rows(connectivity: &[usize]) -> io::Result<Vec<String>> {
    let mut rows = Vec::new();
    for element in elements(connectivity) {
        rows.push(join_ids(element?));
    }
    Ok(rows)
}

/// One text row per element of a fixed-size topology with `npe` nodes per
/// element (`npe` must be non-zero; the last row of a ragged connectivity is
/// shorter).
pub(super) fn fixed_rows(connectivity: &[usize], npe: usize) -> Vec<String> {
    connectivity.chunks(npe.max(1)).map(join_ids).collect()
}

/// Space-separated node ids.
pub(super) fn join_ids(ids: &[usize]) -> String {
    ids.iter()
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}
