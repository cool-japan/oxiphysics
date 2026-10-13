// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Readers for the fixed-width Fortran records used by the AMBER text formats.
//!
//! AMBER writes every record with a fixed column width: `10F8.3` for `mdcrd`
//! trajectories, `6F12.7` for `inpcrd` / `rst7` restart files and `10I8`,
//! `5E16.8` and `20a4` for the `prmtop` topology.  A value that fills its whole
//! column has no blank in front of it, so two neighbouring fields can run
//! together (`1.000-100.000` in an `F8.3` record), which defeats plain
//! whitespace splitting.
//!
//! The helpers in this module accept both layouts.  Blank-separated records are
//! split on whitespace, and a record whose fields touch is cut at multiples of
//! the column width instead.

/// Cut `line` into consecutive chunks of `width` characters.
///
/// Every chunk is returned untrimmed and blank chunks are kept; only the final
/// chunk may be shorter than `width`.  Columns are counted in characters, never
/// bytes, so a record containing multi-byte characters cannot be split inside a
/// character.  A `width` of zero yields no chunks.
pub(crate) fn fixed_chunks(line: &str, width: usize) -> Vec<&str> {
    let mut chunks = Vec::new();
    if width == 0 {
        return chunks;
    }
    let mut start = 0usize;
    let mut columns = 0usize;
    for (offset, _) in line.char_indices() {
        if columns == width {
            chunks.push(&line[start..offset]);
            start = offset;
            columns = 0;
        }
        columns += 1;
    }
    if start < line.len() {
        chunks.push(&line[start..]);
    }
    chunks
}

/// Cut `line` into `width`-character columns, trimming each one and dropping
/// the blank columns.
pub(crate) fn fixed_columns(line: &str, width: usize) -> Vec<&str> {
    fixed_chunks(line, width)
        .into_iter()
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .collect()
}

/// Parse one record of real numbers.
///
/// Blank-separated numbers are returned as they are.  When some blank-separated
/// token is not a number on its own, neighbouring fields touch (a glued pair
/// such as `1.000-100.000` always contains a second decimal point or sign), and
/// the record is re-read as consecutive `width`-character columns.
///
/// # Errors
///
/// Returns a description of the first field that is not a number.
pub(crate) fn parse_real_record(line: &str, width: usize) -> Result<Vec<f64>, String> {
    let blank_separated: Result<Vec<f64>, _> =
        line.split_whitespace().map(str::parse::<f64>).collect();
    if let Ok(values) = blank_separated {
        return Ok(values);
    }
    fixed_columns(line, width)
        .into_iter()
        .map(|field| {
            field
                .parse::<f64>()
                .map_err(|_| format!("`{field}` is not a real number"))
        })
        .collect()
}

/// Parse one record of integers.
///
/// Integers that fill their column run into their left neighbour without any
/// blank, and such a glued pair is itself a valid (larger) integer.  The only
/// reliable signal is the token length: a token wider than a column must hold
/// at least two fields, so the record is cut at multiples of `width`
/// characters in that case.  Otherwise it is split on whitespace.
///
/// # Errors
///
/// Returns a description of the first field that is not an integer.
pub(crate) fn parse_int_record(line: &str, width: usize) -> Result<Vec<i64>, String> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let fields = if tokens.iter().all(|token| token.chars().count() <= width) {
        tokens
    } else {
        fixed_columns(line, width)
    };
    fields
        .into_iter()
        .map(|field| {
            field
                .parse::<i64>()
                .map_err(|_| format!("`{field}` is not an integer"))
        })
        .collect()
}

/// Split one record of `width`-character text fields (`20a4` and similar).
///
/// Blank full-width columns are kept because they are legitimate empty strings;
/// a blank partial column at the end of the line is trailing padding and is
/// dropped.  Every field is returned trimmed.
pub(crate) fn parse_text_record(line: &str, width: usize) -> Vec<String> {
    fixed_chunks(line, width)
        .into_iter()
        .filter(|chunk| !chunk.trim().is_empty() || chunk.chars().count() == width)
        .map(|chunk| chunk.trim().to_string())
        .collect()
}

/// Deterministic helpers for the mutation sweeps in the test modules of the
/// text readers.
///
/// A sweep takes a valid file, damages it in a few random places (replace,
/// delete or insert a character, including multi-byte ones) and feeds it to a
/// reader.  The generator is seeded explicitly, so every run replays the same
/// damaged inputs.
#[cfg(test)]
pub(crate) mod test_support {
    /// SplitMix64 pseudo-random generator.
    pub(crate) struct SplitMix64(u64);

    impl SplitMix64 {
        pub(crate) fn new(seed: u64) -> Self {
            Self(seed)
        }

        pub(crate) fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// A value in `0..bound`; `bound` must be positive.
        pub(crate) fn below(&mut self, bound: usize) -> usize {
            let bound = u64::try_from(bound).unwrap_or(u64::MAX).max(1);
            usize::try_from(self.next_u64() % bound).unwrap_or(0)
        }
    }

    /// Characters inserted by [`mutate_text`]: digits, signs and separators of
    /// the numeric fields, section markers, and multi-byte characters.
    const ALPHABET: [char; 15] = [
        '0',
        '1',
        '9',
        '-',
        '+',
        ' ',
        '.',
        'E',
        '%',
        '\n',
        'x',
        '\u{e9}',
        '\u{3042}',
        '\u{1F600}',
        '\t',
    ];

    /// Apply `edits` random single-character edits to `text`.
    pub(crate) fn mutate_text(text: &str, rng: &mut SplitMix64, edits: usize) -> String {
        let mut chars: Vec<char> = text.chars().collect();
        for _ in 0..edits {
            let replacement = ALPHABET[rng.below(ALPHABET.len())];
            match rng.below(3) {
                0 if !chars.is_empty() => {
                    let at = rng.below(chars.len());
                    chars[at] = replacement;
                }
                1 if !chars.is_empty() => {
                    let at = rng.below(chars.len());
                    chars.remove(at);
                }
                _ => {
                    let at = rng.below(chars.len() + 1);
                    chars.insert(at, replacement);
                }
            }
        }
        chars.into_iter().collect()
    }

    #[test]
    fn generator_is_deterministic_and_bounded() {
        let mut first = SplitMix64::new(7);
        let mut second = SplitMix64::new(7);
        for _ in 0..100 {
            assert_eq!(first.next_u64(), second.next_u64());
        }
        let mut rng = SplitMix64::new(11);
        let mut seen = [false; 5];
        for _ in 0..500 {
            let value = rng.below(5);
            assert!(value < 5);
            seen[value] = true;
        }
        assert!(seen.iter().all(|&hit| hit), "every residue must occur");
    }

    #[test]
    fn mutation_changes_text_without_breaking_characters() {
        let original = "%FLAG POINTERS\n%FORMAT(10I8)\n       3\n";
        let mut rng = SplitMix64::new(3);
        let mut changed = 0usize;
        for _ in 0..200 {
            let mutated = mutate_text(original, &mut rng, 3);
            // `String` guarantees valid UTF-8, so every edit landed on a
            // character boundary.
            if mutated != original {
                changed += 1;
            }
        }
        assert!(
            changed > 150,
            "only {changed} of 200 mutations changed the text"
        );
        assert_eq!(mutate_text("", &mut rng, 0), "");
        // On an empty text the only possible edit is an insertion.
        assert_eq!(mutate_text("", &mut rng, 1).chars().count(), 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_chunks_cuts_at_multiples_of_width() {
        assert_eq!(fixed_chunks("abcdefghij", 4), vec!["abcd", "efgh", "ij"]);
        assert_eq!(fixed_chunks("abcd", 4), vec!["abcd"]);
        assert_eq!(fixed_chunks("abcde", 4), vec!["abcd", "e"]);
        assert!(fixed_chunks("", 4).is_empty());
        assert!(fixed_chunks("abc", 0).is_empty());
    }

    #[test]
    fn fixed_chunks_counts_characters_not_bytes() {
        // Every character below is two bytes long; a byte-based cut would land
        // inside a character and panic.
        let chunks = fixed_chunks("\u{e9}\u{e9}\u{e9}\u{e9}\u{e9}\u{e9}", 4);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].chars().count(), 4);
        assert_eq!(chunks[1].chars().count(), 2);
    }

    #[test]
    fn fixed_columns_trims_and_drops_blank_columns() {
        assert_eq!(fixed_columns("  1.5        -2.5", 6), vec!["1.5", "-2.5"]);
        assert!(fixed_columns("      ", 3).is_empty());
    }

    #[test]
    fn real_record_blank_separated() {
        let values = parse_real_record("  1.5  -2.25 3e2\n", 8).unwrap();
        assert_eq!(values, vec![1.5, -2.25, 300.0]);
        assert!(parse_real_record("   \n", 8).unwrap().is_empty());
    }

    #[test]
    fn real_record_adjacent_negative_fields() {
        // F8.3 fields: "   1.000", "-100.000" and "-200.500" touch each other.
        let line = "   1.000-100.000-200.500   3.000";
        let values = parse_real_record(line, 8).unwrap();
        assert_eq!(values, vec![1.0, -100.0, -200.5, 3.0]);
    }

    #[test]
    fn real_record_adjacent_positive_fields() {
        // Two full-width positive F8.3 fields: "1000.000" "1000.000".
        let values = parse_real_record("1000.0001000.000   2.500", 8).unwrap();
        assert_eq!(values, vec![1000.0, 1000.0, 2.5]);
    }

    #[test]
    fn real_record_adjacent_exponent_fields() {
        // E16.8 fields that fill their column.
        let line = "-1.23456789E-100-9.87654321E+002";
        let values = parse_real_record(line, 16).unwrap();
        assert_eq!(values.len(), 2);
        assert!((values[0] + 1.23456789e-100).abs() < 1e-108);
        assert!((values[1] + 9.87654321e2).abs() < 1e-6);
    }

    #[test]
    fn real_record_rejects_garbage() {
        let err = parse_real_record("   1.000  abc", 8).unwrap_err();
        assert!(err.contains("abc"), "unexpected message: {err}");
        assert!(parse_real_record("********", 8).is_err());
    }

    #[test]
    fn int_record_blank_separated() {
        assert_eq!(
            parse_int_record("       1       2      -3", 8).unwrap(),
            vec![1, 2, -3]
        );
        assert!(parse_int_record("", 8).unwrap().is_empty());
    }

    #[test]
    fn int_record_full_width_fields_are_cut_by_column() {
        // "12345678" fills an I8 column, so it touches the field before it and
        // the whole line is one blank-free token.
        let line = "       312345678       9";
        assert_eq!(parse_int_record(line, 8).unwrap(), vec![3, 12_345_678, 9]);
    }

    #[test]
    fn int_record_rejects_garbage() {
        assert!(parse_int_record("   1   x", 8).is_err());
        assert!(parse_int_record("   1.5", 8).is_err());
    }

    #[test]
    fn text_record_keeps_blank_full_columns() {
        assert_eq!(
            parse_text_record("N   H       CA  ", 4),
            vec!["N", "H", "", "CA"]
        );
    }

    #[test]
    fn text_record_handles_trimmed_trailing_padding() {
        // Trailing blanks were stripped, so the last name is a partial column.
        assert_eq!(parse_text_record("N   H", 4), vec!["N", "H"]);
        // A partial blank column carries no name.
        assert_eq!(parse_text_record("N   H     ", 4), vec!["N", "H"]);
        assert!(parse_text_record("", 4).is_empty());
    }

    #[test]
    fn mutated_records_never_panic_and_chunks_cover_the_line() {
        use super::test_support::{SplitMix64, mutate_text};
        let lines = [
            "   1.000-100.000-200.500   3.000",
            "       312345678       9",
            "N   H       CA  ",
            "-1.23456789E-100-9.87654321E+002",
            "\u{e9}\u{e9}\u{e9}\u{e9}\u{e9}\u{e9}",
        ];
        let mut rng = SplitMix64::new(0xF027_0001);
        let mut parsed = 0usize;
        let mut rejected = 0usize;
        for line in lines {
            for round in 0..400 {
                let damaged = mutate_text(line, &mut rng, 1 + round % 4);
                for width in [1usize, 3, 4, 8, 12, 16] {
                    let chunks = fixed_chunks(&damaged, width);
                    assert_eq!(chunks.concat(), damaged, "width {width}");
                    assert!(chunks.iter().all(|c| c.chars().count() <= width));
                    for field in parse_text_record(&damaged, width) {
                        assert!(field.chars().count() <= width);
                    }
                    match parse_real_record(&damaged, width) {
                        Ok(values) => {
                            parsed += 1;
                            assert!(values.len() <= damaged.chars().count());
                        }
                        Err(message) => {
                            rejected += 1;
                            assert!(message.contains("not a real number"), "{message}");
                        }
                    }
                    match parse_int_record(&damaged, width) {
                        Ok(values) => assert!(values.len() <= damaged.chars().count()),
                        Err(message) => assert!(message.contains("not an integer"), "{message}"),
                    }
                }
            }
        }
        assert!(
            parsed > 0 && rejected > 0,
            "parsed {parsed}, rejected {rejected}"
        );
    }
}
