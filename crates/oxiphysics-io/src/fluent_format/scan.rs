// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parenthesis scanner for ANSYS Fluent `.msh` files.
//!
//! A Fluent mesh is a sequence of parenthesised sections:
//!
//! ```text
//! (10 (1 1 2d4 1 3)(
//! 0.0 0.0 0.0
//! ...
//! ))
//! ```
//!
//! Where the parentheses fall on the lines differs between writers: Fluent
//! puts the `(` that opens the data on the line of the header, others put it
//! on a line of its own, and either may hold several records per line.  The
//! scanner therefore works on the parenthesis structure only and never on
//! lines.  Double quotes protect the parentheses inside comment strings such
//! as `(0 "Zone Sections (cells)")`.

use crate::Error;

/// UTF-8 byte order mark that some tools put in front of text files.
const BYTE_ORDER_MARK: [u8; 3] = [0xEF, 0xBB, 0xBF];

/// Cursor over the bytes of a Fluent mesh file.
pub(super) struct Scanner<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Scanner<'a> {
    /// Start at the beginning of `bytes`, behind a byte order mark if there is
    /// one.
    pub(super) fn new(bytes: &'a [u8]) -> Self {
        let pos = if bytes.starts_with(&BYTE_ORDER_MARK) {
            BYTE_ORDER_MARK.len()
        } else {
            0
        };
        Self { bytes, pos }
    }

    /// Byte offset of the cursor.
    pub(super) fn position(&self) -> usize {
        self.pos
    }

    /// The byte under the cursor.
    pub(super) fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    /// Move the cursor one byte forward.
    pub(super) fn bump(&mut self) {
        if self.pos < self.bytes.len() {
            self.pos += 1;
        }
    }

    pub(super) fn skip_whitespace(&mut self) {
        while self
            .bytes
            .get(self.pos)
            .is_some_and(u8::is_ascii_whitespace)
        {
            self.pos += 1;
        }
    }

    /// Read a run of decimal digits, such as a section number.
    pub(super) fn read_decimal(&mut self) -> Option<u32> {
        let start = self.pos;
        let mut value: u32 = 0;
        while let Some(digit) = self.bytes.get(self.pos).filter(|b| b.is_ascii_digit()) {
            value = value
                .checked_mul(10)?
                .checked_add(u32::from(digit - b'0'))?;
            self.pos += 1;
        }
        (self.pos > start).then_some(value)
    }

    /// 1-based number of the line that holds byte offset `pos`.
    pub(super) fn line_at(&self, pos: usize) -> usize {
        let before = self.bytes.get(..pos.min(self.bytes.len())).unwrap_or(&[]);
        1 + before.iter().filter(|&&byte| byte == b'\n').count()
    }

    /// A parse error that names the line of byte offset `pos`.
    pub(super) fn error_at(&self, pos: usize, message: impl std::fmt::Display) -> Error {
        Error::Parse(format!("Fluent mesh line {}: {message}", self.line_at(pos)))
    }

    /// Read the parenthesised group the cursor is on and return the text
    /// between its parentheses.
    ///
    /// # Errors
    ///
    /// Fails when the cursor is not on a `(`, when the group is not closed, and
    /// when its text is not valid UTF-8.
    pub(super) fn read_group(&mut self) -> Result<&'a str, Error> {
        let start = self.pos;
        let end = self.group_end()?;
        let inner = self.bytes.get(start + 1..end).unwrap_or(&[]);
        self.pos = end + 1;
        std::str::from_utf8(inner)
            .map_err(|_| self.error_at(start, "the text of a section is not valid UTF-8"))
    }

    /// Skip the rest of the section the cursor is in, up to and including the
    /// `)` that closes it.
    ///
    /// # Errors
    ///
    /// Fails when the section is not closed.
    pub(super) fn skip_to_section_end(&mut self) -> Result<(), Error> {
        let end = self.find_close(self.pos, 1)?;
        self.pos = end + 1;
        Ok(())
    }

    /// Consume the `)` that closes the current section, which has to be next
    /// apart from white space.
    ///
    /// # Errors
    ///
    /// Fails when anything else follows; `section_start` names the section in
    /// the message.
    pub(super) fn expect_section_end(&mut self, section_start: usize) -> Result<(), Error> {
        self.skip_whitespace();
        match self.peek() {
            Some(b')') => {
                self.bump();
                Ok(())
            }
            Some(other) => Err(self.error_at(
                self.pos,
                format!(
                    "unexpected {:?} where the section that starts on line {} should end",
                    char::from(other),
                    self.line_at(section_start)
                ),
            )),
            None => Err(self.error_at(
                self.pos,
                format!(
                    "the section that starts on line {} is not closed",
                    self.line_at(section_start)
                ),
            )),
        }
    }

    /// Offset of the `)` that matches the `(` under the cursor.
    fn group_end(&self) -> Result<usize, Error> {
        if self.peek() != Some(b'(') {
            return Err(self.error_at(self.pos, "expected '('"));
        }
        self.find_close(self.pos + 1, 1)
    }

    /// Offset of the `)` that brings the nesting depth, which is `depth` at
    /// `from`, down to zero.
    fn find_close(&self, from: usize, mut depth: usize) -> Result<usize, Error> {
        let mut quoted = false;
        let mut escaped = false;
        let tail = self.bytes.get(from..).unwrap_or(&[]);
        for (offset, &byte) in tail.iter().enumerate() {
            if quoted {
                match byte {
                    _ if escaped => escaped = false,
                    b'\\' => escaped = true,
                    b'"' => quoted = false,
                    _ => {}
                }
                continue;
            }
            match byte {
                b'"' => quoted = true,
                b'(' => depth += 1,
                b')' => {
                    depth -= 1;
                    if depth == 0 {
                        return Ok(from + offset);
                    }
                }
                _ => {}
            }
        }
        Err(self.error_at(
            from.saturating_sub(1),
            "unterminated section: the parenthesis opened here is never closed",
        ))
    }
}

/// White-space separated tokens of a text together with their byte offsets in
/// the file, for error messages that name the line of a bad token.
pub(super) struct Tokens<'a> {
    text: &'a str,
    pos: usize,
    base: usize,
}

impl<'a> Tokens<'a> {
    /// Tokens of `text`, which starts at byte offset `base` of the file.
    pub(super) fn new(text: &'a str, base: usize) -> Self {
        Self { text, pos: 0, base }
    }

    /// File offset of the end of the text.
    pub(super) fn end_offset(&self) -> usize {
        self.base + self.text.len()
    }

    /// The next token and its byte offset in the file.
    pub(super) fn next_token(&mut self) -> Option<(usize, &'a str)> {
        let bytes = self.text.as_bytes();
        while bytes.get(self.pos).is_some_and(u8::is_ascii_whitespace) {
            self.pos += 1;
        }
        let start = self.pos;
        while bytes
            .get(self.pos)
            .is_some_and(|b| !b.is_ascii_whitespace())
        {
            self.pos += 1;
        }
        if self.pos == start {
            None
        } else {
            self.text
                .get(start..self.pos)
                .map(|token| (self.base + start, token))
        }
    }

    /// Number of tokens left.
    pub(super) fn count_remaining(&self) -> usize {
        self.text
            .get(self.pos..)
            .map_or(0, |rest| rest.split_ascii_whitespace().count())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_are_balanced_and_nested_groups_stay_inside() {
        let mut scanner = Scanner::new(b"(a (b c) d)(e)");
        assert_eq!(scanner.read_group().unwrap(), "a (b c) d");
        assert_eq!(scanner.read_group().unwrap(), "e");
        assert_eq!(scanner.peek(), None);
    }

    #[test]
    fn quoted_parentheses_do_not_count() {
        let mut scanner = Scanner::new(br#"(0 "Zone Sections (cells")(next)"#);
        assert_eq!(scanner.read_group().unwrap(), r#"0 "Zone Sections (cells""#);
        assert_eq!(scanner.read_group().unwrap(), "next");
        let mut escaped = Scanner::new(br#"(0 "a \" ) b")"#);
        assert_eq!(escaped.read_group().unwrap(), r#"0 "a \" ) b""#);
    }

    #[test]
    fn unterminated_groups_are_errors_that_name_the_line() {
        let mut scanner = Scanner::new(b"\n\n(10 (1 1 2 1 3)(\n0 0 0\n");
        scanner.skip_whitespace();
        scanner.bump();
        let message = scanner.skip_to_section_end().unwrap_err().to_string();
        assert!(message.contains("unterminated"), "{message}");
        let mut quoted = Scanner::new(b"(0 \"never closed)");
        let message = quoted.read_group().unwrap_err().to_string();
        assert!(message.contains("line 1"), "{message}");
    }

    #[test]
    fn read_group_requires_an_opening_parenthesis() {
        let mut scanner = Scanner::new(b"x(1)");
        assert!(scanner.read_group().is_err());
    }

    #[test]
    fn invalid_utf8_inside_a_group_is_an_error() {
        let mut scanner = Scanner::new(b"(1 \xff 2)");
        let message = scanner.read_group().unwrap_err().to_string();
        assert!(message.contains("UTF-8"), "{message}");
        // The group is still skipped cleanly.
        let mut skipped = Scanner::new(b"(\xff)(2)");
        skipped.bump();
        skipped.skip_to_section_end().unwrap();
        assert_eq!(skipped.read_group().unwrap(), "2");
    }

    #[test]
    fn decimal_numbers_stop_at_the_first_other_byte() {
        let mut scanner = Scanner::new(b"2013 (");
        assert_eq!(scanner.read_decimal(), Some(2013));
        assert_eq!(scanner.peek(), Some(b' '));
        assert_eq!(scanner.read_decimal(), None);
        let mut huge = Scanner::new(b"99999999999");
        assert_eq!(huge.read_decimal(), None);
    }

    #[test]
    fn byte_order_mark_is_skipped() {
        let mut scanner = Scanner::new(b"\xEF\xBB\xBF(2 3)");
        assert_eq!(scanner.read_group().unwrap(), "2 3");
    }

    #[test]
    fn line_numbers_count_newlines_before_the_offset() {
        let scanner = Scanner::new(b"a\nb\r\nc");
        assert_eq!(scanner.line_at(0), 1);
        assert_eq!(scanner.line_at(2), 2);
        assert_eq!(scanner.line_at(5), 3);
        assert_eq!(scanner.line_at(1000), 3);
        let message = scanner.error_at(5, "boom").to_string();
        assert!(
            message.contains("line 3") && message.contains("boom"),
            "{message}"
        );
    }

    #[test]
    fn section_end_must_follow_directly() {
        let mut scanner = Scanner::new(b"  )rest");
        scanner.expect_section_end(0).unwrap();
        assert_eq!(scanner.peek(), Some(b'r'));
        let mut wrong = Scanner::new(b" x)");
        assert!(wrong.expect_section_end(0).is_err());
        let mut eof = Scanner::new(b"   ");
        let message = eof.expect_section_end(0).unwrap_err().to_string();
        assert!(message.contains("not closed"), "{message}");
    }

    #[test]
    fn tokens_report_file_offsets() {
        let mut tokens = Tokens::new(" 1  ab\n cd ", 10);
        assert_eq!(tokens.next_token(), Some((11, "1")));
        assert_eq!(tokens.count_remaining(), 2);
        assert_eq!(tokens.next_token(), Some((14, "ab")));
        assert_eq!(tokens.next_token(), Some((18, "cd")));
        assert_eq!(tokens.next_token(), None);
        assert_eq!(tokens.end_offset(), 21);
    }
}
