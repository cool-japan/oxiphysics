// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parsing of AMBER MDCRD trajectory text (`10F8.3` coordinate records).

use super::super::types::MdcrdFrame;
use crate::fortran_fields::parse_real_record;

/// Number of coordinate values per line of an MDCRD record (`10F8.3`).
const MDCRD_VALUES_PER_LINE: usize = 10;
/// Column width of an MDCRD field (`F8.3`).
const MDCRD_COLUMN_WIDTH: usize = 8;
/// Parse the periodic-box line of an MDCRD frame (`a b c`, extra values are
/// ignored).
fn parse_mdcrd_box(line: &str) -> Option<[f64; 3]> {
    match parse_real_record(line, MDCRD_COLUMN_WIDTH).ok()?.as_slice() {
        [a, b, c, ..] => Some([*a, *b, *c]),
        _ => None,
    }
}
/// Parse AMBER MDCRD trajectory data from a string.
///
/// MDCRD format: a title line, then for every frame the coordinates in
/// `10F8.3` records (ten values per line, fewer on the last line of a frame),
/// followed by one `3F8.3` box line when `has_box` is set.
///
/// Fields are cut into 8-character columns when they touch each other
/// (`1.000-100.000`) and split on blanks otherwise, so no value is ever
/// dropped from a frame.  Only complete frames are returned: parsing stops
/// silently at the first frame whose coordinate lines are truncated, do not
/// hold exactly `3 * n_atoms` numbers, or contain something that is not a
/// number.  A final frame without its box line, or with a box line of fewer
/// than three numbers, is returned with `box_dims == None`.  `n_atoms == 0`
/// yields no frames.
///
/// Use [`parse_mdcrd_checked`] to be told about such defects instead.
pub fn parse_mdcrd(s: &str, n_atoms: usize, has_box: bool) -> Vec<MdcrdFrame> {
    scan_mdcrd(s, n_atoms, has_box, false).unwrap_or_default()
}
/// Parse AMBER MDCRD trajectory data from a string, reporting every defect.
///
/// The format and the handling of touching fields are those of
/// [`parse_mdcrd`], but where that function stops silently, this one fails
/// with a message naming the frame (counted from 0) and, where there is one,
/// the line (counted from 1, the title being line 1).
///
/// # Errors
///
/// Returns a message when `n_atoms` is zero or too large, a coordinate line
/// holds something that is not a number, a frame does not hold exactly
/// `3 * n_atoms` values, the box record of a frame (`has_box`) is missing or
/// does not hold three numbers, or lines remain that do not make up a whole
/// frame.  Blank lines after the last frame are allowed.
pub fn parse_mdcrd_checked(
    s: &str,
    n_atoms: usize,
    has_box: bool,
) -> Result<Vec<MdcrdFrame>, String> {
    if n_atoms == 0 || n_atoms.checked_mul(3).is_none() {
        return Err(format!(
            "{n_atoms} atoms cannot be read from an MDCRD file: a frame needs 3 * n_atoms coordinates"
        ));
    }
    scan_mdcrd(s, n_atoms, has_box, true)
}
/// Core of [`parse_mdcrd`] and [`parse_mdcrd_checked`].
///
/// A defect ends the scan: with `strict` it is an error, otherwise the frames
/// read so far are the result.
fn scan_mdcrd(
    s: &str,
    n_atoms: usize,
    has_box: bool,
    strict: bool,
) -> Result<Vec<MdcrdFrame>, String> {
    let Some(n_coords) = n_atoms.checked_mul(3).filter(|&n| n > 0) else {
        return Ok(Vec::new());
    };
    let coords_lines = n_coords.div_ceil(MDCRD_VALUES_PER_LINE);
    let frame_lines = coords_lines + usize::from(has_box);
    let all_lines: Vec<&str> = s.lines().skip(1).collect();
    // Blank lines after the last record belong to no frame.
    let used = all_lines
        .iter()
        .rposition(|line| !line.trim().is_empty())
        .map_or(0, |last| last + 1);
    let lines = &all_lines[..used];
    let mut frames = Vec::new();
    let mut first = 0usize;
    while first + coords_lines <= lines.len() {
        let frame_index = frames.len();
        let box_line = first + coords_lines;
        let mut coords: Vec<f64> = Vec::with_capacity(n_coords);
        for (offset, line) in lines[first..box_line].iter().enumerate() {
            match parse_real_record(line, MDCRD_COLUMN_WIDTH) {
                Ok(values) => coords.extend(values),
                Err(message) if strict => {
                    return Err(format!(
                        "frame {frame_index}, line {}: {message}",
                        first + offset + 2
                    ));
                }
                Err(_) => return Ok(frames),
            }
        }
        if coords.len() != n_coords {
            if strict {
                return Err(format!(
                    "frame {frame_index}: the coordinate lines hold {} values, expected {n_coords}",
                    coords.len()
                ));
            }
            return Ok(frames);
        }
        let box_dims = if has_box {
            match lines.get(box_line).and_then(|line| parse_mdcrd_box(line)) {
                Some(dims) => Some(dims),
                None if strict => {
                    return Err(format!(
                        "frame {frame_index}, line {}: the box record is missing or does not hold three numbers",
                        box_line + 2
                    ));
                }
                None => None,
            }
        } else {
            None
        };
        frames.push(MdcrdFrame { coords, box_dims });
        first += frame_lines;
    }
    if strict && first < lines.len() {
        return Err(format!(
            "frame {}: {} trailing lines do not make up a whole frame of {coords_lines} coordinate lines",
            frames.len(),
            lines.len() - first
        ));
    }
    Ok(frames)
}
