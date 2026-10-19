//! Fixed-width Fortran records, as used by the `NBLOCK` and `EBLOCK` rows of
//! ANSYS MAPDL CDB decks.
//!
//! The row layout of a CDB block is not implied by the block type: the line
//! after the `NBLOCK` / `EBLOCK` command carries a Fortran `FORMAT`
//! specification such as `(3i9,6e21.13e3)` or `(19i9)`, and every row has to
//! be cut into fields of exactly those widths.  White-space splitting cannot
//! work, because a negative `E21.13E3` value fills its whole field and abuts
//! the previous one (`1.0000000000000E+000-2.0000000000000E+000`).

/// Largest number of fields one format may expand to.
const MAX_FORMAT_FIELDS: usize = 1024;
/// Largest width of a single field.
const MAX_FIELD_WIDTH: usize = 1024;
/// Deepest nesting of parenthesised groups, e.g. `(2(i9,e21.13))`.
const MAX_GROUP_DEPTH: usize = 8;

/// What one edit descriptor reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum FieldKind {
    /// `Iw` / `Iw.m`.
    Integer,
    /// `Fw.d`, `Ew.d[Ee]`, `Dw.d`, `Gw.d`.
    Real,
    /// `Aw`.
    Text,
    /// `nX`: columns that are skipped.
    Skip,
}

/// One column group of a record: what it holds and how wide it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct FieldSpec {
    pub(super) kind: FieldKind,
    pub(super) width: usize,
}

/// A parsed Fortran `FORMAT` specification: the ordered column groups of one
/// record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FortranFormat {
    fields: Vec<FieldSpec>,
}

impl FortranFormat {
    /// Parse a specification such as `(3i9,6e21.13e3)`.
    ///
    /// Supported are the `I`, `F`, `E`, `D`, `G`, `A` and `X` edit
    /// descriptors, repeat counts (`3i9`) and parenthesised groups with
    /// repeat counts (`2(i9,f10.3)`).  Blanks and letter case are ignored.
    pub(super) fn parse(spec: &str) -> Result<Self, String> {
        let text: Vec<u8> = spec
            .bytes()
            .filter(|byte| !byte.is_ascii_whitespace())
            .map(|byte| byte.to_ascii_lowercase())
            .collect();
        let mut parser = FormatParser {
            text: &text,
            pos: 0,
        };
        let mut fields = Vec::new();
        parser.parse_group(&mut fields, 0)?;
        if parser.pos != text.len() {
            return Err(format!(
                "unexpected text after the closing parenthesis of format `{}`",
                spec.trim()
            ));
        }
        Ok(Self { fields })
    }

    /// The column groups of one record, in order.
    pub(super) fn fields(&self) -> &[FieldSpec] {
        &self.fields
    }

    /// Cut `line` into one trimmed slice per column group.
    ///
    /// Groups that lie (partly or wholly) beyond the end of a short line come
    /// back truncated or empty, which the readers treat as blank fields.
    ///
    /// # Errors
    ///
    /// Fails when a group boundary falls inside a multi-byte character.
    pub(super) fn split<'a>(&self, line: &'a str) -> Result<Vec<&'a str>, String> {
        let mut out = Vec::with_capacity(self.fields.len());
        let mut start = 0usize;
        for spec in &self.fields {
            let end = start + spec.width;
            let slice = if start >= line.len() {
                ""
            } else {
                line.get(start..end.min(line.len())).ok_or_else(|| {
                    format!(
                        "column {} of the record lies inside a multi-byte character",
                        start + 1
                    )
                })?
            };
            out.push(slice.trim());
            start = end;
        }
        Ok(out)
    }
}

/// Read an integer field; a blank field reads as zero, as in Fortran.
pub(super) fn read_integer(field: &str) -> Result<i64, String> {
    let text = field.trim();
    if text.is_empty() {
        return Ok(0);
    }
    text.parse::<i64>()
        .map_err(|_| format!("`{text}` is not an integer"))
}

/// Read a real field; a blank field reads as zero, as in Fortran.
pub(super) fn read_real(field: &str) -> Result<f64, String> {
    let text = field.trim();
    if text.is_empty() {
        return Ok(0.0);
    }
    parse_fortran_real(text).ok_or_else(|| format!("`{text}` is not a real number"))
}

/// Parse a real number written the Fortran way.
///
/// Besides the forms Rust understands this accepts the exponent letters `D`
/// and `d` (`1.0D+03`) and the implied exponent that Fortran uses when the
/// exponent letter is dropped (`1.5-100` is `1.5E-100`).  Non-finite values
/// are rejected: a coordinate that reads `NaN` is a corrupt file, not data.
pub(super) fn parse_fortran_real(text: &str) -> Option<f64> {
    let bytes = text.trim().as_bytes();
    if bytes.is_empty() {
        return None;
    }
    let mut normalised = String::with_capacity(bytes.len() + 1);
    let mut exponent_seen = false;
    for (index, &byte) in bytes.iter().enumerate() {
        match byte {
            b'd' | b'D' | b'e' | b'E' => {
                normalised.push('e');
                exponent_seen = true;
            }
            b'+' | b'-'
                if index > 0
                    && !exponent_seen
                    && (bytes[index - 1].is_ascii_digit() || bytes[index - 1] == b'.') =>
            {
                normalised.push('e');
                normalised.push(char::from(byte));
                exponent_seen = true;
            }
            _ => normalised.push(char::from(byte)),
        }
    }
    normalised
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
}

/// Recursive-descent parser over the blank-free, lower-cased specification.
struct FormatParser<'a> {
    text: &'a [u8],
    pos: usize,
}

impl FormatParser<'_> {
    fn peek(&self) -> Option<u8> {
        self.text.get(self.pos).copied()
    }

    fn bump(&mut self) -> Option<u8> {
        let byte = self.peek()?;
        self.pos += 1;
        Some(byte)
    }

    /// Read a run of decimal digits, if there is one.
    fn number(&mut self) -> Result<Option<usize>, String> {
        let start = self.pos;
        let mut value = 0usize;
        while let Some(byte) = self.peek().filter(u8::is_ascii_digit) {
            value = value
                .checked_mul(10)
                .and_then(|v| v.checked_add(usize::from(byte - b'0')))
                .filter(|v| *v <= MAX_FORMAT_FIELDS.max(MAX_FIELD_WIDTH))
                .ok_or_else(|| "a number in the format is too large".to_string())?;
            self.pos += 1;
        }
        Ok((self.pos > start).then_some(value))
    }

    /// A width that must be present and non-zero.
    fn width(&mut self, descriptor: char) -> Result<usize, String> {
        match self.number()? {
            Some(width) if (1..=MAX_FIELD_WIDTH).contains(&width) => Ok(width),
            Some(_) => Err(format!(
                "the width of the `{descriptor}` descriptor must be between 1 and {MAX_FIELD_WIDTH}"
            )),
            None => Err(format!("the `{descriptor}` descriptor needs a width")),
        }
    }

    fn parse_group(&mut self, out: &mut Vec<FieldSpec>, depth: usize) -> Result<(), String> {
        if depth > MAX_GROUP_DEPTH {
            return Err("the format nests parentheses too deeply".to_string());
        }
        if self.bump() != Some(b'(') {
            return Err("a format has to start with `(`".to_string());
        }
        loop {
            self.parse_item(out, depth)?;
            match self.bump() {
                Some(b',') => {}
                Some(b')') => return Ok(()),
                Some(other) => {
                    return Err(format!(
                        "unexpected `{}` in the format; expected `,` or `)`",
                        char::from(other)
                    ));
                }
                None => return Err("the format is missing its closing `)`".to_string()),
            }
        }
    }

    fn parse_item(&mut self, out: &mut Vec<FieldSpec>, depth: usize) -> Result<(), String> {
        let repeat = self.number()?;
        match self.peek() {
            Some(b'(') => {
                let mut inner = Vec::new();
                self.parse_group(&mut inner, depth + 1)?;
                for _ in 0..repeat.unwrap_or(1) {
                    push_all(out, &inner)?;
                }
                Ok(())
            }
            Some(letter) if letter.is_ascii_alphabetic() => {
                self.pos += 1;
                self.parse_descriptor(letter, repeat, out)
            }
            Some(other) => Err(format!(
                "unexpected `{}` in the format; expected an edit descriptor",
                char::from(other)
            )),
            None => Err("the format ends where an edit descriptor was expected".to_string()),
        }
    }

    fn parse_descriptor(
        &mut self,
        letter: u8,
        repeat: Option<usize>,
        out: &mut Vec<FieldSpec>,
    ) -> Result<(), String> {
        let (kind, width) = match letter {
            b'x' => {
                let width = repeat.unwrap_or(1);
                if !(1..=MAX_FIELD_WIDTH).contains(&width) {
                    return Err("the `X` descriptor needs a column count of at least 1".to_string());
                }
                push_all(
                    out,
                    &[FieldSpec {
                        kind: FieldKind::Skip,
                        width,
                    }],
                )?;
                return Ok(());
            }
            b'i' => {
                let width = self.width('I')?;
                if self.peek() == Some(b'.') {
                    self.pos += 1;
                    self.number()?;
                }
                (FieldKind::Integer, width)
            }
            b'f' | b'e' | b'd' | b'g' => {
                let width = self.width(char::from(letter).to_ascii_uppercase())?;
                if self.peek() == Some(b'.') {
                    self.pos += 1;
                    self.number()?;
                    if self.peek() == Some(b'e') {
                        self.pos += 1;
                        self.number()?;
                    }
                }
                (FieldKind::Real, width)
            }
            b'a' => (FieldKind::Text, self.width('A')?),
            other => {
                return Err(format!(
                    "the `{}` edit descriptor is not supported",
                    char::from(other).to_ascii_uppercase()
                ));
            }
        };
        let spec = FieldSpec { kind, width };
        for _ in 0..repeat.unwrap_or(1) {
            push_all(out, &[spec])?;
        }
        Ok(())
    }
}

/// Append `items` to `out`, refusing formats that expand to absurd sizes.
fn push_all(out: &mut Vec<FieldSpec>, items: &[FieldSpec]) -> Result<(), String> {
    if out.len() + items.len() > MAX_FORMAT_FIELDS {
        return Err(format!(
            "the format expands to more than {MAX_FORMAT_FIELDS} fields"
        ));
    }
    out.extend_from_slice(items);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kinds(format: &FortranFormat) -> Vec<(FieldKind, usize)> {
        format
            .fields()
            .iter()
            .map(|spec| (spec.kind, spec.width))
            .collect()
    }

    #[test]
    fn nblock_format_expands_repeat_counts() {
        let format = FortranFormat::parse("(3i9,6e21.13e3)").unwrap();
        let mut expected = vec![(FieldKind::Integer, 9); 3];
        expected.extend(vec![(FieldKind::Real, 21); 6]);
        assert_eq!(kinds(&format), expected);
    }

    #[test]
    fn eblock_format_is_nineteen_integer_fields() {
        let format = FortranFormat::parse("(19i9)").unwrap();
        assert_eq!(kinds(&format), vec![(FieldKind::Integer, 9); 19]);
    }

    #[test]
    fn older_narrow_formats_are_understood() {
        let format = FortranFormat::parse("(3i8,6e16.9)").unwrap();
        assert_eq!(format.fields().len(), 9);
        assert_eq!(format.fields()[0].width, 8);
        assert_eq!(format.fields()[3].width, 16);
        let wide = FortranFormat::parse("(19I10)").unwrap();
        assert_eq!(wide.fields()[18].width, 10);
    }

    #[test]
    fn letter_case_and_blanks_are_ignored() {
        let compact = FortranFormat::parse("(3i9,6e21.13e3)").unwrap();
        let spaced = FortranFormat::parse(" ( 3 I9 , 6 E21.13E3 ) ").unwrap();
        assert_eq!(compact, spaced);
    }

    #[test]
    fn groups_skips_and_other_descriptors() {
        let format = FortranFormat::parse("(2(i4,f6.2),3x,a5,d12.4,g10.3,i5.3)").unwrap();
        assert_eq!(
            kinds(&format),
            vec![
                (FieldKind::Integer, 4),
                (FieldKind::Real, 6),
                (FieldKind::Integer, 4),
                (FieldKind::Real, 6),
                (FieldKind::Skip, 3),
                (FieldKind::Text, 5),
                (FieldKind::Real, 12),
                (FieldKind::Real, 10),
                (FieldKind::Integer, 5),
            ]
        );
    }

    #[test]
    fn malformed_formats_are_rejected() {
        for bad in [
            "",
            "3i9",
            "(3i9",
            "(3i9))",
            "(i)",
            "(i0)",
            "(3q9)",
            "(3i9,)",
            "(/)",
            "(99999999999999999999i9)",
            "(2000i9)",
            "((((((((((i9))))))))))",
        ] {
            assert!(
                FortranFormat::parse(bad).is_err(),
                "format `{bad}` must be rejected"
            );
        }
    }

    #[test]
    fn split_cuts_fixed_widths_where_numbers_abut() {
        let format = FortranFormat::parse("(3i9,6e21.13e3)").unwrap();
        // The second coordinate is negative, so it fills its field and touches
        // the first one: there is no blank to split on.
        let row = "        7        0        0 1.0000000000000E+000-2.5000000000000E-001 3.0000000000000E+001";
        let fields = format.split(row).unwrap();
        assert_eq!(fields.len(), 9);
        assert_eq!(fields[0], "7");
        assert_eq!(fields[3], "1.0000000000000E+000");
        assert_eq!(fields[4], "-2.5000000000000E-001");
        assert_eq!(fields[5], "3.0000000000000E+001");
        assert_eq!(fields[6], "");
        assert_eq!(read_real(fields[4]).unwrap(), -0.25);
        assert_eq!(read_real(fields[5]).unwrap(), 30.0);
        assert_eq!(read_real(fields[8]).unwrap(), 0.0);
    }

    #[test]
    fn split_gives_truncated_and_empty_fields_for_short_lines() {
        let format = FortranFormat::parse("(3i4,2f6.1)").unwrap();
        let fields = format.split("  12   3").unwrap();
        assert_eq!(fields, vec!["12", "3", "", "", ""]);
        let partial = format.split("   1   2   3 1.5").unwrap();
        assert_eq!(partial, vec!["1", "2", "3", "1.5", ""]);
    }

    #[test]
    fn split_reports_a_group_boundary_inside_a_multibyte_character() {
        let format = FortranFormat::parse("(2i2)").unwrap();
        assert!(format.split("1\u{e9}34").is_err());
    }

    #[test]
    fn integer_fields_read_blank_as_zero_and_reject_junk() {
        assert_eq!(read_integer("   ").unwrap(), 0);
        assert_eq!(read_integer(" -1 ").unwrap(), -1);
        assert_eq!(read_integer("+42").unwrap(), 42);
        assert!(read_integer("1.5").is_err());
        assert!(read_integer("*********").is_err());
    }

    #[test]
    fn fortran_reals_in_every_notation() {
        let expected = [
            ("1.0E+000", 1.0),
            ("-2.5000000000000E-001", -0.25),
            ("1.0D+03", 1000.0),
            ("1.0d-3", 0.001),
            ("1.5-100", 1.5e-100),
            ("2.5+3", 2500.0),
            ("  .5 ", 0.5),
            ("7.", 7.0),
            ("-3", -3.0),
            ("+4.25e2", 425.0),
        ];
        for (text, value) in expected {
            let parsed = parse_fortran_real(text).unwrap_or_else(|| panic!("`{text}` must parse"));
            assert!(
                (parsed - value).abs() <= value.abs() * 1e-15,
                "`{text}` parsed as {parsed}, expected {value}"
            );
        }
    }

    #[test]
    fn non_numeric_and_non_finite_reals_are_rejected() {
        for text in [
            "", "abc", "1.0E", "--1", "nan", "inf", "-inf", "1e999", "1,5",
        ] {
            assert!(
                parse_fortran_real(text).is_none(),
                "`{text}` must be rejected"
            );
        }
    }
}
