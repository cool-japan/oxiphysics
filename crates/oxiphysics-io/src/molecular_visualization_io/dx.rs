// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reader and writer for scalar volumetric grids in the OpenDX (`.dx`) text
//! format as produced by APBS and read by VMD.
//!
//! The layout written by APBS is
//!
//! ```text
//! object 1 class gridpositions counts nx ny nz
//! origin ox oy oz
//! delta dxx dxy dxz
//! delta dyx dyy dyz
//! delta dzx dzy dzz
//! object 2 class gridconnections counts nx ny nz
//! object 3 class array type double rank 0 items N data follows
//! <N values, the z index varying fastest>
//! attribute "dep" string "positions"
//! object "regular positions regular connections" class field
//! component "positions" value 1
//! component "connections" value 2
//! component "data" value 3
//! ```
//!
//! Objects are recognised by their `class`, never by their number, so files
//! that place the data array elsewhere are read as well.

use super::functions::{Result, Vec3};
use super::types::{DxGrid, MolVizError};
use std::iter::Peekable;
use std::str::Lines;

fn parse_error(message: String) -> MolVizError {
    MolVizError::ParseError { message }
}

fn dimension_error(message: String) -> MolVizError {
    MolVizError::DimensionMismatch { message }
}

fn missing_section(section: &str) -> MolVizError {
    MolVizError::MissingSection {
        section: section.to_string(),
    }
}

/// The word that follows `class` on an `object` line.
fn object_class(line: &str) -> Option<&str> {
    let mut tokens = line.split_whitespace();
    tokens.find(|token| *token == "class")?;
    tokens.next()
}

/// The token that follows `keyword` on `line`.
fn token_after<'a>(line: &'a str, keyword: &str) -> Option<&'a str> {
    let mut tokens = line.split_whitespace();
    tokens.find(|token| *token == keyword)?;
    tokens.next()
}

/// The three integers after `counts` on a `gridpositions` line.
fn parse_counts(line: &str) -> Result<[usize; 3]> {
    let mut tokens = line.split_whitespace();
    tokens
        .find(|token| *token == "counts")
        .ok_or_else(|| parse_error(format!("DX gridpositions line without counts: '{line}'")))?;
    let counts = tokens
        .map(|token| {
            token
                .parse::<usize>()
                .map_err(|_| parse_error(format!("bad DX grid count '{token}'")))
        })
        .collect::<Result<Vec<usize>>>()?;
    match counts.as_slice() {
        [nx, ny, nz] => Ok([*nx, *ny, *nz]),
        other => Err(dimension_error(format!(
            "a DX grid needs three counts, found {}",
            other.len()
        ))),
    }
}

/// The three numbers after the keyword of an `origin` or `delta` line.
fn parse_vec3(line: &str, what: &str) -> Result<Vec3> {
    let values = line
        .split_whitespace()
        .skip(1)
        .map(|token| {
            token
                .parse::<f64>()
                .map_err(|_| parse_error(format!("bad DX {what} component '{token}'")))
        })
        .collect::<Result<Vec<f64>>>()?;
    match values.as_slice() {
        [x, y, z] => Ok([*x, *y, *z]),
        _ => Err(parse_error(format!(
            "a DX {what} line needs three numbers: '{line}'"
        ))),
    }
}

/// Validate the header of a data array and return its declared item count.
///
/// Only text (`data follows`) scalar arrays are supported.
fn parse_array_header(line: &str) -> Result<Option<usize>> {
    if line
        .split_whitespace()
        .any(|token| matches!(token, "binary" | "msb" | "lsb" | "ieee"))
    {
        return Err(parse_error(
            "binary DX data arrays are not supported, only text arrays".to_string(),
        ));
    }
    if let Some(rank) = token_after(line, "rank")
        && rank != "0"
    {
        return Err(parse_error(format!(
            "only scalar DX data (rank 0) can be stored in a grid, found rank {rank}"
        )));
    }
    token_after(line, "items")
        .map(|token| {
            token
                .parse::<usize>()
                .map_err(|_| parse_error(format!("bad DX item count '{token}'")))
        })
        .transpose()
}

/// Read the numbers of a data block.  The block ends at the first line that
/// opens another DX statement (`attribute`, `object`, `component`, `end`).
fn read_data_block(lines: &mut Peekable<Lines<'_>>) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    while let Some(next) = lines.peek() {
        let line = next.trim();
        if line.is_empty() || line.starts_with('#') {
            lines.next();
            continue;
        }
        if line
            .split_whitespace()
            .next()
            .is_some_and(|word| matches!(word, "attribute" | "object" | "component" | "end"))
        {
            break;
        }
        lines.next();
        for token in line.split_whitespace() {
            let value = token
                .parse::<f64>()
                .map_err(|_| parse_error(format!("bad DX data value '{token}'")))?;
            values.push(value);
        }
    }
    Ok(values)
}

impl DxGrid {
    /// Parse a DX file from its full text content.
    ///
    /// The grid geometry is taken from the `gridpositions` object (`counts`,
    /// `origin` and three `delta` lines).  The values are read from the first
    /// text array of rank 0 (`class array ... data follows`), whichever object
    /// number it carries, up to the next `attribute`, `object`, `component` or
    /// `end` line.  The last comment line mentioning `units` is stored in
    /// [`DxGrid::units`].
    ///
    /// # Errors
    ///
    /// * [`MolVizError::MissingSection`] when the file has no `gridpositions`
    ///   object, no `origin` or `delta` lines, or no data array.
    /// * [`MolVizError::ParseError`] for a malformed number, a binary or
    ///   non-scalar array.
    /// * [`MolVizError::DimensionMismatch`] when a count is zero or overflows,
    ///   or when the number of values (or the declared `items`) differs from
    ///   `nx * ny * nz`.  A short data block is never padded with zeros.
    pub fn parse(content: &str) -> Result<Self> {
        let mut counts: Option<[usize; 3]> = None;
        let mut origin: Option<Vec3> = None;
        let mut deltas: Vec<Vec3> = Vec::new();
        let mut units = String::new();
        let mut data: Option<(Option<usize>, Vec<f64>)> = None;
        // True while the lines belong to the first gridpositions object.
        let mut in_positions = false;

        let mut lines = content.lines().peekable();
        while let Some(raw_line) = lines.next() {
            let line = raw_line.trim();
            if line.starts_with('#') {
                if line.contains("units") {
                    units = line.to_string();
                }
                continue;
            }
            match line.split_whitespace().next() {
                Some("object") => {
                    in_positions = false;
                    match object_class(line) {
                        Some("gridpositions") if counts.is_none() => {
                            counts = Some(parse_counts(line)?);
                            in_positions = true;
                        }
                        Some("array") if data.is_none() && line.contains("data follows") => {
                            let declared = parse_array_header(line)?;
                            data = Some((declared, read_data_block(&mut lines)?));
                        }
                        _ => {}
                    }
                }
                Some("origin") if in_positions => origin = Some(parse_vec3(line, "origin")?),
                Some("delta") if in_positions => deltas.push(parse_vec3(line, "delta")?),
                _ => {}
            }
        }

        let [nx, ny, nz] = counts.ok_or_else(|| missing_section("DX gridpositions object"))?;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(dimension_error(format!(
                "DX grid counts {nx} x {ny} x {nz} contain a zero"
            )));
        }
        let expected = nx
            .checked_mul(ny)
            .and_then(|n| n.checked_mul(nz))
            .ok_or_else(|| {
                dimension_error(format!("DX grid counts {nx} x {ny} x {nz} overflow"))
            })?;
        let origin = origin.ok_or_else(|| missing_section("DX origin"))?;
        let [delta_x, delta_y, delta_z] = match deltas.as_slice() {
            [] => return Err(missing_section("DX delta")),
            [x, y, z] => [*x, *y, *z],
            other => {
                return Err(parse_error(format!(
                    "a DX grid needs exactly three delta lines, found {}",
                    other.len()
                )));
            }
        };
        let (declared, values) =
            data.ok_or_else(|| missing_section("DX data array (class array, data follows)"))?;
        if let Some(items) = declared
            && items != expected
        {
            return Err(dimension_error(format!(
                "DX array declares {items} items but a {nx} x {ny} x {nz} grid has {expected} points"
            )));
        }
        if values.len() != expected {
            return Err(dimension_error(format!(
                "DX grid {nx} x {ny} x {nz} needs {expected} values but the data block holds {}",
                values.len()
            )));
        }
        Ok(DxGrid {
            nx,
            ny,
            nz,
            origin,
            delta_x,
            delta_y,
            delta_z,
            data: values,
            quantity: String::from("potential"),
            units,
        })
    }

    /// Serialize this grid to DX text in the layout APBS writes (see the module
    /// documentation).
    ///
    /// Numbers are written with the shortest representation that parses back
    /// to the identical `f64`, so `DxGrid::parse(&grid.to_dx_string())`
    /// reproduces the geometry and every value exactly.  The data array is
    /// written in the order of [`DxGrid::data`] (z index fastest), three values
    /// per line; the item count is the actual length of `data`.
    pub fn to_dx_string(&self) -> String {
        let mut out = String::new();
        out.push_str("# Data from OxiPhysics\n");
        if !self.units.is_empty() {
            let units = self.units.replace(['\n', '\r'], " ");
            if !units.starts_with('#') {
                out.push_str("# ");
            }
            out.push_str(&units);
            out.push('\n');
        }
        out.push_str(&format!(
            "object 1 class gridpositions counts {} {} {}\n",
            self.nx, self.ny, self.nz
        ));
        out.push_str(&format!(
            "origin {:e} {:e} {:e}\n",
            self.origin[0], self.origin[1], self.origin[2]
        ));
        for delta in [self.delta_x, self.delta_y, self.delta_z] {
            out.push_str(&format!(
                "delta {:e} {:e} {:e}\n",
                delta[0], delta[1], delta[2]
            ));
        }
        out.push_str(&format!(
            "object 2 class gridconnections counts {} {} {}\n",
            self.nx, self.ny, self.nz
        ));
        out.push_str(&format!(
            "object 3 class array type double rank 0 items {} data follows\n",
            self.data.len()
        ));
        for row in self.data.chunks(3) {
            let formatted: Vec<String> = row.iter().map(|v| format!("{v:e}")).collect();
            out.push_str(&formatted.join(" "));
            out.push('\n');
        }
        out.push_str("attribute \"dep\" string \"positions\"\n");
        out.push_str("object \"regular positions regular connections\" class field\n");
        out.push_str("component \"positions\" value 1\n");
        out.push_str("component \"connections\" value 2\n");
        out.push_str("component \"data\" value 3\n");
        out
    }
}
