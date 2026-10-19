// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Token-level reader for the part of the CIF 1.1 syntax that crystal
//! structure files use.
//!
//! CIF is a token stream, not a line format: a data name may be followed by
//! its value on the next line, a `loop_` table runs until the next data name,
//! `loop_` or `data_` token, and multi-line values are written as
//! semicolon-delimited text fields.  This module turns a file into the items
//! and loops of its first data block; interpreting them is left to the caller.

use crate::{Error, Result};
use std::collections::HashMap;

/// One lexical item of a CIF file.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// A data name such as `_cell_length_a`, stored in lower case.
    Tag(String),
    /// A value: bare, quoted, or a semicolon-delimited text field.
    Value(String),
    /// The `loop_` keyword.
    Loop,
    /// A `data_<name>` block header.
    Data,
    /// `save_...`, `global_` or `stop_`: structural keywords that end a loop
    /// and carry nothing this reader uses.
    Keyword,
}

/// A `loop_` table: its column data names and the flat, row-major values.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct CifLoop {
    /// Column data names in file order, lower case.
    pub(super) headers: Vec<String>,
    /// Values in file order; row `r`, column `c` is `values[r * headers.len() + c]`.
    pub(super) values: Vec<String>,
}

impl CifLoop {
    /// Index of the first column named by any of `tags`, trying `tags` in
    /// order of preference.
    pub(super) fn column(&self, tags: &[&str]) -> Option<usize> {
        tags.iter()
            .find_map(|tag| self.headers.iter().position(|h| h == tag))
    }

    /// The table rows.
    ///
    /// Fails when the number of values is not a whole number of rows, which
    /// means the table is truncated or corrupt.
    pub(super) fn rows(&self) -> Result<Vec<&[String]>> {
        let width = self.headers.len();
        if width == 0 || !self.values.len().is_multiple_of(width) {
            return Err(Error::Parse(format!(
                "CIF: loop over {} has {} values, which is not a multiple of its {} columns",
                self.headers.first().map(String::as_str).unwrap_or("?"),
                self.values.len(),
                width
            )));
        }
        Ok(self.values.chunks_exact(width).collect())
    }
}

/// The items and loops of the first data block of a CIF file.
#[derive(Debug, Clone, Default)]
pub(super) struct CifBlock {
    /// Data name (lower case) to value, for `name value` items outside loops.
    pub(super) items: HashMap<String, String>,
    /// Every `loop_` table in file order.
    pub(super) loops: Vec<CifLoop>,
}

impl CifBlock {
    /// The first loop that has a column named by any of `tags`.
    pub(super) fn loop_with_any(&self, tags: &[&str]) -> Option<&CifLoop> {
        self.loops.iter().find(|l| l.column(tags).is_some())
    }
}

/// Split one physical line into tokens, appending them to `out`.
///
/// A `#` that starts a token begins a comment that runs to the end of the
/// line.  A quoted value ends at the first matching quote that is followed by
/// whitespace or the end of the line, so quotes inside a value are allowed.
fn lex_line(line: &str, out: &mut Vec<Token>) {
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0usize;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        if c == '#' {
            return;
        }
        if c == '\'' || c == '"' {
            let mut end = None;
            let mut j = i + 1;
            while j < chars.len() {
                if chars[j] == c && chars.get(j + 1).is_none_or(|next| next.is_whitespace()) {
                    end = Some(j);
                    break;
                }
                j += 1;
            }
            match end {
                Some(e) => {
                    out.push(Token::Value(chars[i + 1..e].iter().collect()));
                    i = e + 1;
                }
                None => {
                    // Unterminated quote: take the rest of the line.
                    out.push(Token::Value(chars[i + 1..].iter().collect()));
                    i = chars.len();
                }
            }
            continue;
        }
        let start = i;
        while i < chars.len() && !chars[i].is_whitespace() {
            i += 1;
        }
        out.push(classify_bare(chars[start..i].iter().collect()));
    }
}

/// Classify an unquoted token as a data name, a keyword or a value.
fn classify_bare(text: String) -> Token {
    if text.starts_with('_') {
        return Token::Tag(text.to_ascii_lowercase());
    }
    let lower = text.to_ascii_lowercase();
    if lower == "loop_" {
        Token::Loop
    } else if lower.starts_with("data_") {
        Token::Data
    } else if lower.starts_with("save_") || lower == "global_" || lower == "stop_" {
        Token::Keyword
    } else {
        Token::Value(text)
    }
}

/// Turn the whole file into a token stream.
fn lex(content: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let mut lines = content.lines();
    while let Some(line) = lines.next() {
        let Some(first_text_line) = line.strip_prefix(';') else {
            lex_line(line, &mut tokens);
            continue;
        };
        // Semicolon-delimited text field: runs to the next line that starts
        // with ';' and is a single value, whatever it contains.
        let mut text = first_text_line.to_string();
        let mut remainder = None;
        for next in lines.by_ref() {
            if let Some(after) = next.strip_prefix(';') {
                remainder = Some(after);
                break;
            }
            text.push('\n');
            text.push_str(next);
        }
        let Some(after_field) = remainder else {
            return Err(Error::Parse(
                "CIF: unterminated semicolon text field".to_string(),
            ));
        };
        tokens.push(Token::Value(text));
        lex_line(after_field, &mut tokens);
    }
    Ok(tokens)
}

/// Parse the first data block of `content` into items and loops.
///
/// A `loop_` ends at the next data name, `loop_`, `data_` or structural
/// keyword, so items that follow a table are read as ordinary items again.
/// Later data blocks are ignored.
pub(super) fn parse_block(content: &str) -> Result<CifBlock> {
    let tokens = lex(content)?;
    let mut block = CifBlock::default();
    let mut seen_data_header = false;
    let mut i = 0usize;
    while i < tokens.len() {
        match &tokens[i] {
            Token::Data => {
                if seen_data_header {
                    break;
                }
                seen_data_header = true;
                i += 1;
            }
            Token::Loop => {
                i += 1;
                let mut headers = Vec::new();
                while let Some(Token::Tag(name)) = tokens.get(i) {
                    headers.push(name.clone());
                    i += 1;
                }
                if headers.is_empty() {
                    return Err(Error::Parse(
                        "CIF: loop_ is not followed by any data name".to_string(),
                    ));
                }
                let mut values = Vec::new();
                while let Some(Token::Value(value)) = tokens.get(i) {
                    values.push(value.clone());
                    i += 1;
                }
                block.loops.push(CifLoop { headers, values });
            }
            Token::Tag(name) => {
                if let Some(Token::Value(value)) = tokens.get(i + 1) {
                    block.items.insert(name.clone(), value.clone());
                    i += 2;
                } else {
                    // A data name without a value (truncated item).
                    block.items.insert(name.clone(), String::new());
                    i += 1;
                }
            }
            Token::Value(_) | Token::Keyword => i += 1,
        }
    }
    Ok(block)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn values_of(line: &str) -> Vec<String> {
        let mut tokens = Vec::new();
        lex_line(line, &mut tokens);
        tokens
            .into_iter()
            .filter_map(|t| match t {
                Token::Value(v) => Some(v),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn quoted_values_keep_their_spaces() {
        assert_eq!(values_of("'hello world' 3.14"), ["hello world", "3.14"]);
        assert_eq!(values_of("\"a b\" c"), ["a b", "c"]);
    }

    #[test]
    fn a_quote_only_closes_before_whitespace() {
        assert_eq!(values_of("'it's here' next"), ["it's here", "next"]);
        assert_eq!(values_of("O1' bare"), ["O1'", "bare"]);
    }

    #[test]
    fn unterminated_quote_takes_the_rest_of_the_line() {
        assert_eq!(values_of("'never closed here"), ["never closed here"]);
    }

    #[test]
    fn trailing_comments_are_ignored() {
        assert_eq!(values_of("1.0 2.0 # not 3.0"), ["1.0", "2.0"]);
        assert_eq!(values_of("# only a comment"), Vec::<String>::new());
        // A '#' inside a bare token is data, not a comment.
        assert_eq!(values_of("a#b c"), ["a#b", "c"]);
    }

    #[test]
    fn keywords_are_case_insensitive_and_quoted_text_is_never_a_keyword() {
        let mut tokens = Vec::new();
        lex_line("LOOP_ Data_Foo 'loop_' _Tag", &mut tokens);
        assert_eq!(
            tokens,
            vec![
                Token::Loop,
                Token::Data,
                Token::Value("loop_".to_string()),
                Token::Tag("_tag".to_string()),
            ]
        );
    }

    #[test]
    fn value_on_the_line_after_its_tag_is_read() {
        let block = parse_block("data_x\n_cell_length_a\n  5.43\n_cell_length_b 6.0\n").unwrap();
        assert_eq!(block.items["_cell_length_a"], "5.43");
        assert_eq!(block.items["_cell_length_b"], "6.0");
    }

    #[test]
    fn every_loop_is_kept_and_the_next_tag_closes_a_loop() {
        let cif = "\
data_x
loop_
_a_col1
_a_col2
1 2
3 4
loop_
_b_col
x
y
z
_after_loop value
";
        let block = parse_block(cif).unwrap();
        assert_eq!(block.loops.len(), 2);
        assert_eq!(block.loops[0].headers, ["_a_col1", "_a_col2"]);
        assert_eq!(block.loops[0].values, ["1", "2", "3", "4"]);
        assert_eq!(block.loops[1].headers, ["_b_col"]);
        assert_eq!(block.loops[1].values, ["x", "y", "z"]);
        assert_eq!(block.items["_after_loop"], "value");
    }

    #[test]
    fn data_header_ends_a_loop_and_only_the_first_block_is_read() {
        let cif = "data_one\n_k 1\nloop_\n_c\n1\n2\ndata_two\n_k 2\nloop_\n_d\n9\n";
        let block = parse_block(cif).unwrap();
        assert_eq!(block.items["_k"], "1");
        assert_eq!(block.loops.len(), 1);
        assert_eq!(block.loops[0].values, ["1", "2"]);
    }

    #[test]
    fn semicolon_text_fields_are_single_values() {
        let cif = "\
data_x
_title
;
loop_
_not_a_tag a very long line
;
_after 1
loop_
_t
_u
;multi
line
; plain
";
        let block = parse_block(cif).unwrap();
        assert_eq!(
            block.items["_title"],
            "\nloop_\n_not_a_tag a very long line"
        );
        assert_eq!(block.items["_after"], "1");
        assert_eq!(block.loops.len(), 1);
        assert_eq!(block.loops[0].values, ["multi\nline", "plain"]);
    }

    #[test]
    fn unterminated_text_field_is_an_error() {
        let err = parse_block("data_x\n_a\n;never closed\n").unwrap_err();
        assert!(err.to_string().contains("semicolon"), "{err}");
    }

    #[test]
    fn loop_without_data_names_is_an_error() {
        assert!(parse_block("data_x\nloop_\n1 2 3\n").is_err());
    }

    #[test]
    fn rows_require_whole_rows() {
        let ok = CifLoop {
            headers: vec!["_a".into(), "_b".into()],
            values: vec!["1".into(), "2".into(), "3".into(), "4".into()],
        };
        let rows = ok.rows().unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[1], ["3", "4"]);
        let ragged = CifLoop {
            headers: vec!["_a".into(), "_b".into()],
            values: vec!["1".into(), "2".into(), "3".into()],
        };
        assert!(ragged.rows().is_err());
    }

    #[test]
    fn column_lookup_honours_the_order_of_preference() {
        let lp = CifLoop {
            headers: vec!["_x".into(), "_old".into(), "_new".into()],
            values: Vec::new(),
        };
        assert_eq!(lp.column(&["_new", "_old"]), Some(2));
        assert_eq!(lp.column(&["_missing", "_old"]), Some(1));
        assert_eq!(lp.column(&["_missing"]), None);
    }
}
