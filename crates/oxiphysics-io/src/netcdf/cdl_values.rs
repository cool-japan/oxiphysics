// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Typed declarations and value lists of the CDL text written by
//! [`write_ncf_text`](super::write_ncf_text) and read by
//! [`parse_ncf_text`](super::parse_ncf_text).

use super::classic::NcType;
use super::types::NetCdfVariable;

/// Render `var.data` for the CDL `data:` section.
///
/// Numeric variables are written as `", "`-separated numbers.  `NC_CHAR`
/// variables are written as quoted strings, one per row of length `row_len`
/// (the last dimension), so that `ncgen` and [`parse_cdl_values`] restore the
/// exact bytes; non-printable bytes use `\xHH` escapes.
pub(super) fn format_cdl_values(var: &NetCdfVariable, row_len: Option<usize>) -> String {
    if var.nc_type == NcType::Char {
        let bytes: Vec<u8> = var
            .data
            .iter()
            .map(|&v| {
                if (0.0..=255.0).contains(&v) {
                    v as u8
                } else {
                    0
                }
            })
            .collect();
        let row = row_len.filter(|&r| r > 0).unwrap_or(bytes.len().max(1));
        return bytes
            .chunks(row)
            .map(quote_cdl_bytes)
            .collect::<Vec<_>>()
            .join(", ");
    }
    var.data
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

fn quote_cdl_bytes(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() + 2);
    out.push('"');
    for &b in bytes {
        match b {
            b'"' => out.push_str("\\\""),
            b'\\' => out.push_str("\\\\"),
            // `//` would start a CDL comment.
            b'/' => out.push_str("\\x2f"),
            b'\n' => out.push_str("\\n"),
            b'\t' => out.push_str("\\t"),
            0x20..=0x7E => out.push(char::from(b)),
            other => out.push_str(&format!("\\x{other:02x}")),
        }
    }
    out.push('"');
    out
}

/// Parse the value list of a CDL data statement (the text after `=`).
///
/// Quoted strings (as written for `NC_CHAR` variables) yield their byte
/// codes, concatenated; anything else must be a comma-separated list of
/// numbers.
pub(super) fn parse_cdl_values(text: &str) -> Result<Vec<f64>, String> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(Vec::new());
    }
    if !text.starts_with('"') {
        return text
            .split(',')
            .map(|s| s.trim().parse::<f64>().map_err(|e| format!("{s:?}: {e}")))
            .collect();
    }
    let mut out = Vec::new();
    let mut chars = text.chars().peekable();
    loop {
        while chars.peek().is_some_and(|c| c.is_whitespace()) {
            chars.next();
        }
        match chars.next() {
            Some('"') => {}
            None => break,
            Some(c) => return Err(format!("expected a quoted string, found {c:?}")),
        }
        loop {
            match chars.next() {
                None => return Err("unterminated string".to_string()),
                Some('"') => break,
                Some('\\') => {
                    let byte = match chars.next() {
                        Some('n') => b'\n',
                        Some('t') => b'\t',
                        Some('r') => b'\r',
                        Some('0') => 0,
                        Some('"') => b'"',
                        Some('\\') => b'\\',
                        Some('\'') => b'\'',
                        Some('x') => {
                            let hex: String = (0..2).filter_map(|_| chars.next()).collect();
                            u8::from_str_radix(&hex, 16)
                                .map_err(|_| format!("bad \\x escape {hex:?}"))?
                        }
                        other => return Err(format!("unknown escape {other:?}")),
                    };
                    out.push(f64::from(byte));
                }
                Some(c) => {
                    let mut buf = [0u8; 4];
                    out.extend(c.encode_utf8(&mut buf).bytes().map(f64::from));
                }
            }
        }
        while chars.peek().is_some_and(|c| c.is_whitespace()) {
            chars.next();
        }
        match chars.next() {
            Some(',') | None => {}
            Some(c) => return Err(format!("expected ',' between strings, found {c:?}")),
        }
    }
    Ok(out)
}

/// A parsed CDL variable declaration.
pub(super) struct CdlDeclaration {
    /// Declared external type.
    pub(super) nc_type: NcType,
    /// Variable name.
    pub(super) name: String,
    /// Dimension names (empty for a scalar).
    pub(super) dims: Vec<String>,
}

/// Parse a CDL variable declaration such as `float coordinates(frame, atom)`
/// or the scalar form `double t`.  Returns `None` when the line is not a
/// declaration (its first word is not a type keyword).
pub(super) fn parse_cdl_declaration(line: &str) -> Option<Result<CdlDeclaration, String>> {
    let line = line.trim();
    let (keyword, rest) = line.split_once(char::is_whitespace)?;
    let nc_type = NcType::from_cdl_name(keyword)?;
    let rest = rest.trim();
    let (name, dims) = match rest.find('(') {
        Some(paren) => {
            let Some(inner) = rest[paren + 1..].trim_end().strip_suffix(')') else {
                return Some(Err(format!("unterminated dimension list in {line:?}")));
            };
            let dims: Vec<String> = inner
                .split(',')
                .map(|d| d.trim().to_string())
                .filter(|d| !d.is_empty())
                .collect();
            (rest[..paren].trim().to_string(), dims)
        }
        None => (rest.to_string(), Vec::new()),
    };
    if name.is_empty() || name.contains(char::is_whitespace) {
        return Some(Err(format!("bad variable declaration {line:?}")));
    }
    Some(Ok(CdlDeclaration {
        nc_type,
        name,
        dims,
    }))
}
