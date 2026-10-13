// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Text serialisation of [`Rst7File`] (`rst7` / `inpcrd` restart files).
//!
//! Layout of the ASCII restart file:
//!
//! ```text
//! title
//! natom [time]
//! x1 y1 z1 x2 y2 z2 ...     coordinates, 6F12.7 (six values per line)
//! vx1 vy1 vz1 ...           velocities (optional, same format)
//! a b c alpha beta gamma    periodic box (optional, 6F12.7)
//! ```
//!
//! The body has no markers saying which optional blocks are present; they are
//! recognised by the number of values.  With `N` atoms a file holds `3N`
//! (coordinates), `3N + 6` (coordinates and box), `6N` (coordinates and
//! velocities) or `6N + 6` (all three) values.  The only count shared by two
//! layouts is `N = 2` with twelve values, which [`Rst7File::read_with_hint`]
//! resolves.

use super::{AmberIoError, Rst7File};
use crate::fortran_fields::parse_real_record;
use std::io::{BufRead, BufReader, Read, Write};

/// Column width of an `F12.7` field.
const FIELD_WIDTH: usize = 12;
/// Fields per line of the coordinate, velocity and box records.
const FIELDS_PER_LINE: usize = 6;
/// Number of values in the box record (`a b c alpha beta gamma`).
const BOX_VALUES: usize = 6;

fn parse_error(message: impl Into<String>) -> AmberIoError {
    AmberIoError::ParseError(message.into())
}

fn invalid_data(message: impl Into<String>) -> AmberIoError {
    AmberIoError::InvalidData(message.into())
}

/// Format `value` as an `F12.7` field.
fn format_field(value: f64) -> Result<String, AmberIoError> {
    let field = format!("{value:FIELD_WIDTH$.7}");
    if field.len() > FIELD_WIDTH {
        return Err(invalid_data(format!(
            "value {value} does not fit a 12.7 field, so the columns of the file would shift"
        )));
    }
    Ok(field)
}

/// Append `values` to `text`, [`FIELDS_PER_LINE`] per line.
fn push_block(text: &mut String, values: &[f64]) -> Result<(), AmberIoError> {
    for line in values.chunks(FIELDS_PER_LINE) {
        for &value in line {
            text.push_str(&format_field(value)?);
        }
        text.push('\n');
    }
    Ok(())
}

/// Optional blocks that follow the coordinates.
#[derive(Clone, Copy)]
struct Layout {
    velocities: bool,
    box_record: bool,
}

impl Layout {
    /// Number of values a body with this layout holds for `n_coord` coordinates.
    fn value_count(self, n_coord: usize) -> usize {
        let blocks = if self.velocities { 2 } else { 1 };
        let box_values = if self.box_record { BOX_VALUES } else { 0 };
        n_coord.saturating_mul(blocks).saturating_add(box_values)
    }

    /// Pick the layout of a body of `count` values.
    ///
    /// The candidates are tried in the order coordinates, coordinates + box,
    /// coordinates + velocities, coordinates + velocities + box, so the
    /// ambiguous two-atom case resolves to coordinates + box.  A
    /// `has_velocities` hint restricts the candidates.
    fn classify(count: usize, n_coord: usize, has_velocities: Option<bool>) -> Option<Self> {
        [(false, false), (false, true), (true, false), (true, true)]
            .into_iter()
            .map(|(velocities, box_record)| Layout {
                velocities,
                box_record,
            })
            .filter(|layout| has_velocities.is_none_or(|wanted| wanted == layout.velocities))
            .find(|layout| layout.value_count(n_coord) == count)
    }
}

impl Rst7File {
    /// Write the restart file to any [`Write`] sink.
    ///
    /// Coordinates, velocities and the box are written in `F12.7` fixed-point
    /// format, six values per line (two atoms per line, matching the AMBER
    /// specification).  Nothing is written when the file cannot be
    /// represented.
    ///
    /// # Errors
    ///
    /// Returns [`AmberIoError::InvalidData`] when the title spans several
    /// lines, the coordinate count is not a multiple of three, the velocities
    /// are neither empty nor as long as the coordinates, or a value does not
    /// fit its 12-character column (which would shift every later field).
    /// [`AmberIoError::IoError`] reports a failing sink.
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), AmberIoError> {
        if self.title.contains(['\n', '\r']) {
            return Err(invalid_data("the title must be a single line"));
        }
        if !self.coordinates.len().is_multiple_of(3) {
            return Err(invalid_data(format!(
                "{} coordinate values do not form whole atoms",
                self.coordinates.len()
            )));
        }
        if !self.velocities.is_empty() && self.velocities.len() != self.coordinates.len() {
            return Err(invalid_data(format!(
                "{} velocity values for {} coordinate values",
                self.velocities.len(),
                self.coordinates.len()
            )));
        }

        let mut text = String::new();
        text.push_str(&self.title);
        text.push('\n');
        match self.time {
            Some(time) => text.push_str(&format!("{:6}{:15.7e}\n", self.natom(), time)),
            None => text.push_str(&format!("{:6}\n", self.natom())),
        }
        push_block(&mut text, &self.coordinates)?;
        push_block(&mut text, &self.velocities)?;
        if let Some(box_dimensions) = self.box_dimensions {
            push_block(&mut text, &box_dimensions)?;
        }
        writer.write_all(text.as_bytes())?;
        writer.flush()?;
        Ok(())
    }

    /// Parse an AMBER restart file from any [`Read`] source.
    ///
    /// Equivalent to [`Rst7File::read_with_hint`] with no hint: a file whose
    /// value count fits both coordinates + box and coordinates + velocities
    /// (two atoms, twelve values) is read as coordinates + box.
    ///
    /// # Errors
    ///
    /// See [`Rst7File::read_with_hint`].
    pub fn read<R: Read>(reader: R) -> Result<Self, AmberIoError> {
        Self::read_with_hint(reader, None)
    }

    /// Parse an AMBER restart file, telling the reader whether the file holds
    /// velocities.
    ///
    /// The body is read as `F12.7` columns, which also separates fields that
    /// touch each other (`-100.0000000` fills its whole column); blank
    /// separated values are accepted as well.  The number of values selects
    /// the layout (see the module documentation).  `has_velocities` is only
    /// needed for the one ambiguous case, two atoms with twelve values, where
    /// `Some(true)` selects coordinates + velocities and `None` or
    /// `Some(false)` selects coordinates + box.  For every other atom count
    /// the hint is checked against the layout the value count implies.
    ///
    /// # Errors
    ///
    /// Returns [`AmberIoError::UnexpectedEof`] when the title or the atom
    /// count line is missing and [`AmberIoError::ParseError`] when the atom
    /// count, the time or a value is not a number, when the value count fits
    /// none of the layouts above, or when it contradicts `has_velocities`.
    pub fn read_with_hint<R: Read>(
        reader: R,
        has_velocities: Option<bool>,
    ) -> Result<Self, AmberIoError> {
        let mut lines = BufReader::new(reader).lines();
        let mut rst = Rst7File::new();

        rst.title = lines.next().ok_or(AmberIoError::UnexpectedEof)??;

        let header = lines.next().ok_or(AmberIoError::UnexpectedEof)??;
        let mut header_fields = header.split_whitespace();
        let natom: usize = header_fields
            .next()
            .ok_or_else(|| parse_error("missing natom"))?
            .parse()
            .map_err(|e: std::num::ParseIntError| parse_error(format!("natom: {e}")))?;
        if let Some(time) = header_fields.next() {
            rst.time = Some(
                time.parse::<f64>()
                    .map_err(|e| parse_error(format!("time: {e}")))?,
            );
        }

        let mut values: Vec<f64> = Vec::new();
        for line in lines {
            let line = line?;
            values.extend(
                parse_real_record(&line, FIELD_WIDTH)
                    .map_err(|e| parse_error(format!("restart body: {e}")))?,
            );
        }

        let n_coord = natom.saturating_mul(3);
        let layout = Layout::classify(values.len(), n_coord, has_velocities).ok_or_else(|| {
            let hint = match has_velocities {
                Some(true) => ", with velocities",
                Some(false) => ", without velocities",
                None => "",
            };
            parse_error(format!(
                "restart body holds {} values for {natom} atoms{hint}; expected {} (coordinates), {} (coordinates and box), {} (coordinates and velocities) or {} (coordinates, velocities and box)",
                values.len(),
                n_coord,
                n_coord.saturating_add(BOX_VALUES),
                n_coord.saturating_mul(2),
                n_coord.saturating_mul(2).saturating_add(BOX_VALUES),
            ))
        })?;

        let short = || parse_error("the restart body is shorter than its layout requires");
        let (coordinates, rest) = values.split_at_checked(n_coord).ok_or_else(short)?;
        rst.coordinates = coordinates.to_vec();
        let rest = if layout.velocities {
            let (velocities, rest) = rest.split_at_checked(n_coord).ok_or_else(short)?;
            rst.velocities = velocities.to_vec();
            rest
        } else {
            rest
        };
        if layout.box_record {
            rst.box_dimensions = Some(
                <[f64; BOX_VALUES]>::try_from(rest)
                    .map_err(|_| parse_error("the box record does not hold six values"))?,
            );
        }
        Ok(rst)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const BOX: [f64; 6] = [20.0, 21.5, 22.25, 90.0, 90.0, 109.4712206];

    fn fields(values: &[f64]) -> String {
        values
            .chunks(FIELDS_PER_LINE)
            .map(|row| row.iter().map(|v| format!("{v:12.7}")).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Restart text with the given atom count and body blocks.
    fn restart(natom: usize, blocks: &[&[f64]]) -> String {
        let mut text = format!("test restart\n{natom:6}\n");
        for block in blocks {
            text.push_str(&fields(block));
            text.push('\n');
        }
        text
    }

    fn read_str(text: &str) -> Result<Rst7File, AmberIoError> {
        Rst7File::read(Cursor::new(text.as_bytes()))
    }

    fn read_hint(text: &str, hint: Option<bool>) -> Result<Rst7File, AmberIoError> {
        Rst7File::read_with_hint(Cursor::new(text.as_bytes()), hint)
    }

    fn series(start: f64, count: usize) -> Vec<f64> {
        (0..count).map(|i| start + i as f64 * 0.25).collect()
    }

    // ── layouts ───────────────────────────────────────────────────────────────

    fn assert_all_layouts(natom: usize) {
        let coords = series(1.0, 3 * natom);
        let vels = series(-0.5, 3 * natom);

        let only_coords = read_str(&restart(natom, &[&coords])).unwrap();
        assert_eq!(only_coords.coordinates, coords);
        assert!(only_coords.velocities.is_empty());
        assert_eq!(only_coords.box_dimensions, None);
        assert_eq!(only_coords.natom(), natom);

        let with_box = read_str(&restart(natom, &[&coords, &BOX])).unwrap();
        assert_eq!(with_box.coordinates, coords);
        assert!(with_box.velocities.is_empty(), "the box is not velocities");
        assert_eq!(with_box.box_dimensions, Some(BOX));

        let with_vel = read_str(&restart(natom, &[&coords, &vels])).unwrap();
        assert_eq!(with_vel.coordinates, coords);
        assert_eq!(with_vel.velocities, vels);
        assert_eq!(with_vel.box_dimensions, None);

        let everything = read_str(&restart(natom, &[&coords, &vels, &BOX])).unwrap();
        assert_eq!(everything.coordinates, coords);
        assert_eq!(everything.velocities, vels);
        assert_eq!(everything.box_dimensions, Some(BOX));
    }

    #[test]
    fn all_four_layouts_for_one_atom() {
        assert_all_layouts(1);
    }

    #[test]
    fn all_four_layouts_for_three_atoms() {
        assert_all_layouts(3);
    }

    #[test]
    fn all_four_layouts_for_other_atom_counts() {
        for natom in [4, 5, 10, 17] {
            assert_all_layouts(natom);
        }
    }

    #[test]
    fn solvated_inpcrd_box_is_not_mistaken_for_velocities() {
        // tleap writes coordinates followed by the box: the old reader turned
        // the box into velocities for a single atom.
        let text = restart(
            1,
            &[&[1.0, 2.0, 3.0], &[50.0, 50.0, 50.0, 90.0, 90.0, 90.0]],
        );
        let rst = read_str(&text).unwrap();
        assert_eq!(
            rst.box_dimensions,
            Some([50.0, 50.0, 50.0, 90.0, 90.0, 90.0])
        );
        assert!(rst.velocities.is_empty());
        assert_eq!(rst.coordinates, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn two_atoms_with_twelve_values_follow_the_hint() {
        let first = series(1.0, 6);
        let second = series(-3.0, 6);
        let text = restart(2, &[&first, &second]);

        for hint in [None, Some(false)] {
            let rst = read_hint(&text, hint).unwrap();
            assert_eq!(rst.coordinates, first);
            assert!(rst.velocities.is_empty());
            assert_eq!(
                rst.box_dimensions,
                Some(<[f64; 6]>::try_from(&second[..]).unwrap())
            );
        }
        let rst = read_hint(&text, Some(true)).unwrap();
        assert_eq!(rst.coordinates, first);
        assert_eq!(rst.velocities, second);
        assert_eq!(rst.box_dimensions, None);

        // `read` is the unhinted entry point.
        assert_eq!(
            read_str(&text).unwrap().box_dimensions,
            read_hint(&text, None).unwrap().box_dimensions
        );
    }

    #[test]
    fn two_atoms_unambiguous_counts_accept_matching_hints() {
        let coords = series(1.0, 6);
        let vels = series(0.5, 6);
        let only = restart(2, &[&coords]);
        assert!(read_hint(&only, Some(false)).is_ok());
        assert!(read_hint(&only, None).is_ok());
        let both = restart(2, &[&coords, &vels, &BOX]);
        let rst = read_hint(&both, Some(true)).unwrap();
        assert_eq!(rst.velocities, vels);
        assert_eq!(rst.box_dimensions, Some(BOX));
        let rst = read_hint(&restart(2, &[&coords, &vels]), Some(true)).unwrap();
        assert_eq!(rst.velocities, vels);
    }

    #[test]
    fn hint_contradicting_the_value_count_is_an_error() {
        let coords = series(1.0, 6);
        let vels = series(0.5, 6);
        // coordinates only, but the caller claims velocities
        let err = read_hint(&restart(2, &[&coords]), Some(true)).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
        assert!(err.to_string().contains("with velocities"), "{err}");
        // coordinates + velocities + box, but the caller claims none
        let err = read_hint(&restart(2, &[&coords, &vels, &BOX]), Some(false)).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
        // one atom: coordinates + box cannot hold velocities
        let one = restart(1, &[&[1.0, 2.0, 3.0], &BOX]);
        assert!(read_hint(&one, Some(true)).is_err());
        assert!(read_hint(&one, Some(false)).is_ok());
    }

    #[test]
    fn value_count_fitting_no_layout_is_an_error() {
        for count in [1usize, 2, 4, 5, 7, 8, 11, 13, 17, 19] {
            let text = restart(2, &[&series(0.0, count)]);
            let err = read_str(&text).unwrap_err();
            assert!(
                matches!(err, AmberIoError::ParseError(_)),
                "count {count}: {err}"
            );
            assert!(err.to_string().contains("expected 6"), "{err}");
        }
        // three atoms announced, four coordinates present
        let err = read_str(&restart(3, &[&series(0.0, 4)])).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
    }

    #[test]
    fn zero_atoms() {
        let rst = read_str("empty\n     0\n").unwrap();
        assert!(rst.coordinates.is_empty());
        assert_eq!(rst.natom(), 0);
        assert_eq!(rst.box_dimensions, None);
        let rst = read_str(&restart(0, &[&BOX])).unwrap();
        assert!(rst.coordinates.is_empty());
        assert_eq!(rst.box_dimensions, Some(BOX));
    }

    // ── touching fields ───────────────────────────────────────────────────────

    #[test]
    fn negative_hundreds_round_trip() {
        let mut original = Rst7File::new();
        original.title = "touching fields".to_string();
        original.coordinates = vec![-150.0, 1.5, -999.9999999, 120.25, -100.0, 1000.0];
        original.velocities = vec![-100.5, 0.25, 3.0, -0.0000001, 1.0, -2.0];
        original.box_dimensions = Some([-120.0, 130.0, 140.0, 90.0, 90.0, 90.0]);
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(
            text.lines().any(|line| line
                .split_whitespace()
                .any(|token| token.len() > FIELD_WIDTH)),
            "the test needs touching fields:\n{text}"
        );
        let restored = Rst7File::read_with_hint(Cursor::new(text.as_bytes()), Some(true)).unwrap();
        assert_eq!(restored.title, original.title);
        assert_eq!(restored.coordinates, original.coordinates);
        assert_eq!(restored.velocities, original.velocities);
        assert_eq!(restored.box_dimensions, original.box_dimensions);
    }

    #[test]
    fn minus_150_round_trips_through_the_default_reader() {
        let mut original = Rst7File::new();
        original.title = "t".to_string();
        original.coordinates = vec![-150.0, -150.0, -150.0];
        original.box_dimensions = Some([-150.0; 6]);
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert_eq!(restored.coordinates, vec![-150.0; 3]);
        assert_eq!(restored.box_dimensions, Some([-150.0; 6]));
    }

    #[test]
    fn blank_separated_files_are_still_accepted() {
        let text = "hand written\n    2\n 1.0 2.0 3.0\n 4.0 5.0 6.0\n";
        let rst = read_str(text).unwrap();
        assert_eq!(rst.coordinates, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert!(rst.velocities.is_empty());
    }

    #[test]
    fn blank_lines_in_the_body_are_ignored() {
        let text = "t\n    1\n\n   1.0000000   2.0000000   3.0000000\n\n";
        assert_eq!(read_str(text).unwrap().coordinates, vec![1.0, 2.0, 3.0]);
    }

    // ── header ────────────────────────────────────────────────────────────────

    #[test]
    fn time_is_read_from_the_atom_count_line() {
        let text = "t\n     1  5.0000000e-1\n   1.0000000   2.0000000   3.0000000\n";
        let rst = read_str(text).unwrap();
        assert_eq!(rst.time, Some(0.5));
        assert_eq!(rst.natom(), 1);
    }

    #[test]
    fn missing_header_lines_are_unexpected_eof() {
        assert!(matches!(read_str(""), Err(AmberIoError::UnexpectedEof)));
        assert!(matches!(
            read_str("only a title\n"),
            Err(AmberIoError::UnexpectedEof)
        ));
    }

    #[test]
    fn malformed_header_and_body_are_parse_errors() {
        assert!(matches!(
            read_str("t\n\n"),
            Err(AmberIoError::ParseError(_))
        ));
        assert!(matches!(
            read_str("t\nabc\n"),
            Err(AmberIoError::ParseError(_))
        ));
        assert!(matches!(
            read_str("t\n    1 later\n"),
            Err(AmberIoError::ParseError(_))
        ));
        let err = read_str("t\n    1\n   1.0000000   x.0000000   3.0000000\n").unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
    }

    // ── writing ───────────────────────────────────────────────────────────────

    fn sample(natom: usize) -> Rst7File {
        let mut rst = Rst7File::new();
        rst.title = "sample".to_string();
        rst.time = Some(2.5);
        rst.coordinates = series(-10.0, 3 * natom);
        rst
    }

    #[test]
    fn written_file_uses_six_values_per_line() {
        let mut rst = sample(3);
        rst.box_dimensions = Some(BOX);
        let mut buf = Vec::new();
        rst.write(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "sample");
        assert_eq!(lines[1].split_whitespace().next(), Some("3"));
        // 9 coordinates -> two lines, then the box line
        assert_eq!(lines.len(), 2 + 2 + 1);
        assert_eq!(lines[2].len(), 6 * FIELD_WIDTH);
        assert_eq!(lines[3].len(), 3 * FIELD_WIDTH);
        assert_eq!(lines[4].len(), 6 * FIELD_WIDTH);
        assert_eq!(read_str(&text).unwrap().box_dimensions, Some(BOX));
    }

    /// A sink that rejects every write.
    struct ClosedSink;

    impl Write for ClosedSink {
        fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("sink closed"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("sink closed"))
        }
    }

    #[test]
    fn writer_reports_a_failing_sink() {
        let err = sample(2).write(&mut ClosedSink).unwrap_err();
        assert!(matches!(err, AmberIoError::IoError(_)), "{err}");
        assert!(err.to_string().contains("sink closed"), "{err}");
    }

    #[test]
    fn writer_rejects_inconsistent_arrays_without_writing() {
        let mut rst = sample(2);
        rst.coordinates.pop();
        let mut buf = Vec::new();
        let err = rst.write(&mut buf).unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
        assert!(buf.is_empty());

        let mut rst = sample(2);
        rst.velocities = vec![0.0; 3];
        let err = rst.write(&mut buf).unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
        assert!(buf.is_empty());

        let mut rst = sample(1);
        rst.title = "two\nlines".to_string();
        assert!(matches!(
            rst.write(&mut buf),
            Err(AmberIoError::InvalidData(_))
        ));
        assert!(buf.is_empty());
    }

    #[test]
    fn writer_rejects_values_that_overflow_the_column() {
        for bad in [12_345.678, -1000.0, 1.0e9] {
            let mut rst = sample(1);
            rst.coordinates[1] = bad;
            let mut buf = Vec::new();
            let err = rst.write(&mut buf).unwrap_err();
            assert!(matches!(err, AmberIoError::InvalidData(_)), "{bad}: {err}");
            assert!(buf.is_empty());
        }
        // The widest values that still fit.
        let mut rst = sample(1);
        rst.coordinates = vec![9999.9999999, -999.9999999, 0.0];
        let mut buf = Vec::new();
        rst.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert_eq!(restored.coordinates, rst.coordinates);
    }

    #[test]
    fn large_restart_round_trips() {
        let mut rst = sample(1000);
        rst.velocities = series(0.0, 3000);
        rst.box_dimensions = Some(BOX);
        let mut buf = Vec::new();
        rst.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert_eq!(restored.natom(), 1000);
        assert_eq!(restored.coordinates, rst.coordinates);
        assert_eq!(restored.velocities, rst.velocities);
        assert_eq!(restored.box_dimensions, Some(BOX));
        assert_eq!(restored.time, Some(2.5));
    }

    /// `(value count, coordinate count, hint, expected (velocities, box))`.
    type ClassifyCase = (usize, usize, Option<bool>, Option<(bool, bool)>);

    #[test]
    fn layout_classification_table() {
        let cases: [ClassifyCase; 8] = [
            (3, 3, None, Some((false, false))),
            (9, 3, None, Some((false, true))),
            (6, 3, None, Some((true, false))),
            (12, 3, None, Some((true, true))),
            (12, 6, None, Some((false, true))),
            (12, 6, Some(true), Some((true, false))),
            (12, 6, Some(false), Some((false, true))),
            (10, 6, None, None),
        ];
        for (count, n_coord, hint, expected) in cases {
            let got = Layout::classify(count, n_coord, hint).map(|l| (l.velocities, l.box_record));
            assert_eq!(
                got, expected,
                "count {count}, coordinates {n_coord}, hint {hint:?}"
            );
        }
    }

    #[test]
    fn mutated_restarts_never_panic_and_stay_consistent() {
        use crate::fortran_fields::test_support::{SplitMix64, mutate_text};
        let mut full = sample(3);
        full.velocities = series(-1.0, 9);
        full.box_dimensions = Some(BOX);
        let mut buf = Vec::new();
        full.write(&mut buf).unwrap();
        let seeds = [
            String::from_utf8(buf).unwrap(),
            restart(1, &[&[1.0, -150.0, 3.0], &BOX]),
            "hand written\n    2\n 1.0 2.0 3.0\n 4.0 5.0 6.0\n".to_string(),
        ];
        let hints = [None, Some(true), Some(false)];
        let mut rng = SplitMix64::new(0x7257_0001);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        for seed in &seeds {
            for round in 0..1500 {
                let text = mutate_text(seed, &mut rng, 1 + round % 4);
                match read_hint(&text, hints[round % hints.len()]) {
                    Ok(rst) => {
                        accepted += 1;
                        assert_eq!(rst.coordinates.len() % 3, 0, "{text}");
                        assert!(
                            rst.velocities.is_empty()
                                || rst.velocities.len() == rst.coordinates.len(),
                            "{text}"
                        );
                    }
                    Err(_) => rejected += 1,
                }
            }
        }
        assert!(accepted > 0, "no damaged restart was accepted");
        assert!(rejected > 0, "no damaged restart was rejected");
    }
}
