// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crystallographic Information File (CIF 1.1) reader.
//!
//! The file is split into tokens exactly as the CIF syntax defines them, and the
//! tokens are then folded into data blocks:
//!
//! * `#` starts a comment when it begins a token and runs to the end of the line.
//! * `'...'` and `"..."` are quoted strings.  A quote ends the string only when it
//!   is followed by white space, so `'O'Brien'` is one value and `'x, y, z'` keeps
//!   its blanks.  Quoted strings never span lines.
//! * A line that starts with `;` opens a text field that runs to the next line
//!   starting with `;`.
//! * `_name` is a data name, `loop_` starts a table, `data_<name>` starts a block.
//!   A table lasts until the next data name, `loop_`, `data_`, `save_`, `global_`,
//!   `stop_` or the end of the input.  Its values are one stream that is cut into
//!   rows by the number of data names, so a row may span lines and a line may hold
//!   several rows.
//!
//! Numbers carry their standard uncertainty in parentheses (`0.2500(3)`); see
//! [`parse_cif_number`].

use std::collections::HashMap;

use super::functions::Result;
use super::types::{CifAtomSite, CifBlock, MolVizError};

/// Data name of the fractional x coordinate of an atom site.
const FRACT_X: &str = "_atom_site_fract_x";
/// Data name of the fractional y coordinate of an atom site.
const FRACT_Y: &str = "_atom_site_fract_y";
/// Data name of the fractional z coordinate of an atom site.
const FRACT_Z: &str = "_atom_site_fract_z";

/// One lexical item of a CIF file.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// `data_<name>`: opens a data block called `<name>`.
    DataBlock(String),
    /// `loop_`.
    Loop,
    /// `_name`: a data name.
    Tag(String),
    /// `save_...`, `global_` or `stop_`: reserved words that end a table.
    Reserved,
    /// A data value with its quotes or text-field delimiters removed.
    Value(String),
}

/// A token and the 1-based line it starts on.
#[derive(Debug, Clone, PartialEq)]
struct Lexed {
    token: Token,
    line: usize,
}

fn parse_error(message: String) -> MolVizError {
    MolVizError::ParseError { message }
}

/// The CIF white-space characters.
fn is_blank(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n')
}

/// `word` without the ASCII case-insensitive `prefix`, when it starts with it.
fn strip_prefix_ignore_case<'a>(word: &'a str, prefix: &str) -> Option<&'a str> {
    let head = word.get(..prefix.len())?;
    if head.eq_ignore_ascii_case(prefix) {
        word.get(prefix.len()..)
    } else {
        None
    }
}

/// Classify an unquoted word.
fn classify_word(word: &str) -> Token {
    if word.starts_with('_') {
        Token::Tag(word.to_string())
    } else if let Some(name) = strip_prefix_ignore_case(word, "data_") {
        Token::DataBlock(name.to_string())
    } else if word.eq_ignore_ascii_case("loop_") {
        Token::Loop
    } else if strip_prefix_ignore_case(word, "save_").is_some()
        || word.eq_ignore_ascii_case("global_")
        || word.eq_ignore_ascii_case("stop_")
    {
        Token::Reserved
    } else {
        Token::Value(word.to_string())
    }
}

/// Split `content` into CIF tokens.
fn lex(content: &str) -> Result<Vec<Lexed>> {
    let bytes = content.as_bytes();
    let len = bytes.len();
    let mut tokens = Vec::new();
    let mut pos = 0usize;
    let mut line = 1usize;
    let mut line_start = 0usize;
    while pos < len {
        let byte = bytes[pos];
        match byte {
            b'\n' => {
                pos += 1;
                line += 1;
                line_start = pos;
            }
            b' ' | b'\t' | b'\r' => pos += 1,
            b'#' => {
                while pos < len && bytes[pos] != b'\n' {
                    pos += 1;
                }
            }
            b';' if pos == line_start => {
                let body_start = pos + 1;
                let Some(offset) = content[body_start..].find("\n;") else {
                    return Err(parse_error(format!(
                        "line {line}: text field opened by ';' is never closed"
                    )));
                };
                let body_end = body_start + offset;
                let body = &content[body_start..body_end];
                let value = body.strip_suffix('\r').unwrap_or(body);
                let value = value
                    .strip_prefix("\r\n")
                    .or_else(|| value.strip_prefix('\n'))
                    .unwrap_or(value);
                tokens.push(Lexed {
                    token: Token::Value(value.to_string()),
                    line,
                });
                line += body.matches('\n').count() + 1;
                line_start = body_end + 1;
                pos = body_end + 2;
            }
            b'\'' | b'"' => {
                let start = pos + 1;
                let mut cursor = start;
                let mut end = None;
                while cursor < len && bytes[cursor] != b'\n' {
                    if bytes[cursor] == byte && bytes.get(cursor + 1).is_none_or(|&b| is_blank(b)) {
                        end = Some(cursor);
                        break;
                    }
                    cursor += 1;
                }
                let Some(end) = end else {
                    return Err(parse_error(format!(
                        "line {line}: quoted string is not closed before the end of the line"
                    )));
                };
                tokens.push(Lexed {
                    token: Token::Value(content[start..end].to_string()),
                    line,
                });
                pos = end + 1;
            }
            _ => {
                let start = pos;
                while pos < len && !is_blank(bytes[pos]) {
                    pos += 1;
                }
                tokens.push(Lexed {
                    token: classify_word(&content[start..pos]),
                    line,
                });
            }
        }
    }
    Ok(tokens)
}

/// Data name in canonical form: lower case, `.` category separators as `_`.
fn normalize_tag(tag: &str) -> String {
    tag.chars()
        .map(|c| {
            if c == '.' {
                '_'
            } else {
                c.to_ascii_lowercase()
            }
        })
        .collect()
}

/// Whether a table with these data names lists atom sites with fractional coordinates.
fn is_atom_site_loop(tags: &[String]) -> bool {
    let names: Vec<String> = tags.iter().map(|t| normalize_tag(t)).collect();
    [FRACT_X, FRACT_Y, FRACT_Z]
        .iter()
        .all(|wanted| names.iter().any(|n| n == wanted))
}

/// What the parser is waiting for.
#[derive(Debug)]
enum State {
    /// Between items.
    Idle,
    /// A data name was read; its value comes next.
    AwaitValue { name: String, line: usize },
    /// `loop_` was read; the data names are being collected.
    LoopHeader { tags: Vec<String>, line: usize },
    /// The values of a table are being collected (`keep` is false for tables that
    /// are not atom-site tables, whose values are skipped).
    LoopBody {
        tags: Vec<String>,
        values: Vec<String>,
        keep: bool,
    },
}

/// Folds the token stream into data blocks.
#[derive(Debug)]
struct Parser {
    blocks: Vec<CifBlock>,
    current: Option<CifBlock>,
    state: State,
}

impl Parser {
    fn new() -> Self {
        Self {
            blocks: Vec::new(),
            current: None,
            state: State::Idle,
        }
    }

    /// Close whatever is open: a dangling data name is an error, a finished table
    /// is converted.
    fn finish_state(&mut self) -> Result<()> {
        let state = std::mem::replace(&mut self.state, State::Idle);
        match state {
            State::Idle | State::LoopHeader { .. } => Ok(()),
            State::AwaitValue { name, line } => Err(parse_error(format!(
                "line {line}: data name '{name}' has no value"
            ))),
            State::LoopBody { tags, values, keep } => {
                if keep && let Some(blk) = self.current.as_mut() {
                    process_cif_loop(&tags, &values, blk)?;
                }
                Ok(())
            }
        }
    }

    fn feed(&mut self, item: Lexed) -> Result<()> {
        let Lexed { token, line } = item;
        if let Token::DataBlock(name) = token {
            self.finish_state()?;
            if let Some(blk) = self.current.take() {
                self.blocks.push(blk);
            }
            self.current = Some(CifBlock {
                name,
                tags: HashMap::new(),
                atom_sites: Vec::new(),
            });
            return Ok(());
        }
        if self.current.is_none() {
            // Anything before the first data block is not part of a block.
            return Ok(());
        }
        match token {
            // A data block header was handled above.
            Token::DataBlock(_) => Ok(()),
            Token::Tag(name) => self.feed_tag(name, line),
            Token::Loop => {
                if matches!(&self.state, State::LoopHeader { tags, .. } if tags.is_empty()) {
                    return Err(parse_error(format!(
                        "line {line}: loop_ is followed by another loop_ instead of data names"
                    )));
                }
                self.finish_state()?;
                self.state = State::LoopHeader {
                    tags: Vec::new(),
                    line,
                };
                Ok(())
            }
            Token::Reserved => {
                self.check_loop_header_has_tags(line)?;
                self.finish_state()
            }
            Token::Value(value) => self.feed_value(value, line),
        }
    }

    fn check_loop_header_has_tags(&self, line: usize) -> Result<()> {
        if matches!(&self.state, State::LoopHeader { tags, .. } if tags.is_empty()) {
            return Err(parse_error(format!("line {line}: loop_ has no data names")));
        }
        Ok(())
    }

    fn feed_tag(&mut self, name: String, line: usize) -> Result<()> {
        match &mut self.state {
            State::LoopHeader { tags, .. } => {
                tags.push(name);
                Ok(())
            }
            _ => {
                self.finish_state()?;
                self.state = State::AwaitValue { name, line };
                Ok(())
            }
        }
    }

    fn feed_value(&mut self, value: String, line: usize) -> Result<()> {
        let state = std::mem::replace(&mut self.state, State::Idle);
        match state {
            State::Idle => Err(parse_error(format!(
                "line {line}: value '{value}' does not belong to a data name"
            ))),
            State::AwaitValue { name, .. } => {
                if let Some(blk) = self.current.as_mut() {
                    blk.tags.insert(name, value);
                }
                Ok(())
            }
            State::LoopHeader {
                tags,
                line: loop_line,
            } => {
                if tags.is_empty() {
                    return Err(parse_error(format!(
                        "line {loop_line}: loop_ has no data names"
                    )));
                }
                let keep = is_atom_site_loop(&tags);
                let values = if keep { vec![value] } else { Vec::new() };
                self.state = State::LoopBody { tags, values, keep };
                Ok(())
            }
            State::LoopBody {
                tags,
                mut values,
                keep,
            } => {
                if keep {
                    values.push(value);
                }
                self.state = State::LoopBody { tags, values, keep };
                Ok(())
            }
        }
    }

    fn finish(mut self) -> Result<Vec<CifBlock>> {
        self.finish_state()?;
        if let Some(blk) = self.current.take() {
            self.blocks.push(blk);
        }
        Ok(self.blocks)
    }
}

/// Parse a CIF file into its data blocks.
///
/// Every `_name value` pair of a block is stored in [`CifBlock::tags`] (quotes and
/// text-field delimiters removed, standard uncertainties left in place; use
/// [`parse_cif_number`] to read them).  Tables whose data names include
/// `_atom_site_fract_x`, `_atom_site_fract_y` and `_atom_site_fract_z` are read into
/// [`CifBlock::atom_sites`]; other tables (symmetry operators,
/// `_atom_site_aniso_*` displacement parameters, bonds, ...) are skipped.
/// Content before the first `data_` line is ignored.
///
/// # Errors
///
/// [`MolVizError::ParseError`] for a quoted string or text field that is not closed,
/// a data name without a value, a value that belongs to no data name, a `loop_`
/// without data names, an atom-site table whose number of values is not a multiple
/// of its number of data names, and an atom site whose coordinates or occupancy are
/// not numbers (or whose coordinates are `?` / `.`).
pub fn parse_cif(content: &str) -> Result<Vec<CifBlock>> {
    let content = content.strip_prefix('\u{feff}').unwrap_or(content);
    let mut parser = Parser::new();
    for item in lex(content)? {
        parser.feed(item)?;
    }
    parser.finish()
}

/// Parse a CIF number, dropping its standard uncertainty.
///
/// `4.2500(3)` is `4.25`, `-0.0012(4)` is `-0.0012`, `1.5e-3(2)` is `0.0015`.
/// The CIF markers for *inapplicable* (`.`) and *unknown* (`?`) values give
/// `Ok(None)`.
///
/// # Errors
///
/// [`MolVizError::ParseError`] when the text is neither a marker nor a finite
/// decimal number optionally followed by `(digits)`.
pub fn parse_cif_number(text: &str) -> Result<Option<f64>> {
    let trimmed = text.trim();
    if trimmed == "." || trimmed == "?" {
        return Ok(None);
    }
    match strip_standard_uncertainty(trimmed).parse::<f64>() {
        Ok(value) if value.is_finite() => Ok(Some(value)),
        _ => Err(parse_error(format!("invalid CIF number '{trimmed}'"))),
    }
}

/// `4.25(3)` without its `(3)`; text without a trailing `(digits)` is unchanged.
fn strip_standard_uncertainty(text: &str) -> &str {
    let Some(body) = text.strip_suffix(')') else {
        return text;
    };
    let Some(open) = body.rfind('(') else {
        return text;
    };
    let digits = &body[open + 1..];
    if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        &body[..open]
    } else {
        text
    }
}

/// The element symbol a site label starts with: an initial letter and, when it is
/// followed by a lower-case letter, that letter (`Fe1` is `Fe`, `O2` is `O`).
fn element_from_label(label: &str) -> String {
    let mut chars = label.chars();
    match chars.next() {
        Some(first) if first.is_ascii_alphabetic() => {
            let mut symbol = String::from(first);
            if let Some(second) = chars.next()
                && second.is_ascii_lowercase()
            {
                symbol.push(second);
            }
            symbol
        }
        _ => "X".to_string(),
    }
}

/// Read the atom sites of one table into `blk`.
///
/// Tables without `_atom_site_fract_x/y/z` columns are not atom-site tables and
/// leave `blk` untouched.  A missing occupancy (no column, `?` or `.`) is `1.0`.
fn process_cif_loop(headers: &[String], values: &[String], blk: &mut CifBlock) -> Result<()> {
    let names: Vec<String> = headers.iter().map(|h| normalize_tag(h)).collect();
    let column = |wanted: &str| names.iter().position(|n| n == wanted);
    let (Some(idx_x), Some(idx_y), Some(idx_z)) =
        (column(FRACT_X), column(FRACT_Y), column(FRACT_Z))
    else {
        return Ok(());
    };
    let idx_label = column("_atom_site_label");
    let idx_type = column("_atom_site_type_symbol").or_else(|| column("_atom_site_element"));
    let idx_occ = column("_atom_site_occupancy");
    let width = headers.len();
    if !values.len().is_multiple_of(width) {
        return Err(parse_error(format!(
            "atom-site table has {} values for {width} data names; the last row is incomplete",
            values.len()
        )));
    }
    for row in values.chunks_exact(width) {
        let label = idx_label.map_or_else(|| "X".to_string(), |i| row[i].clone());
        let element = idx_type.map_or_else(|| element_from_label(&label), |i| row[i].clone());
        let coordinate = |index: usize, tag: &str| -> Result<f64> {
            match parse_cif_number(&row[index]) {
                Ok(Some(value)) => Ok(value),
                Ok(None) => Err(parse_error(format!(
                    "atom site '{label}': {tag} is '{}', the position is unknown",
                    row[index].trim()
                ))),
                Err(e) => Err(parse_error(format!("atom site '{label}': {tag}: {e}"))),
            }
        };
        let fx = coordinate(idx_x, FRACT_X)?;
        let fy = coordinate(idx_y, FRACT_Y)?;
        let fz = coordinate(idx_z, FRACT_Z)?;
        let occupancy = match idx_occ {
            Some(i) => parse_cif_number(&row[i])
                .map_err(|e| parse_error(format!("atom site '{label}': occupancy: {e}")))?
                .unwrap_or(1.0),
            None => 1.0,
        };
        blk.atom_sites.push(CifAtomSite {
            label,
            element,
            fx,
            fy,
            fz,
            occupancy,
        });
    }
    Ok(())
}
