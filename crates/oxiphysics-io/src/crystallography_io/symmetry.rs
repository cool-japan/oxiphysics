// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Space-group symmetry operations in the CIF `x,y,z` notation, and expansion
//! of an asymmetric unit into the full unit cell.

use super::CrystalAtom;
use crate::{Error, Result};
use std::collections::HashMap;

/// Fractional-coordinate tolerance under which two atoms of the same element
/// occupy the same site.
pub(super) const SITE_TOLERANCE: f64 = 1e-4;

/// Number of hash cells per axis used to find coincident sites.  The cell
/// edge (`1 / CELLS_PER_AXIS`) is larger than [`SITE_TOLERANCE`], so a site
/// can only coincide with atoms in its own or an adjacent cell.
const CELLS_PER_AXIS: i32 = 1000;

/// One space-group operation `x' = R . x + t` acting on fractional coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct SymmetryOp {
    /// Rotation part `R`; row `i` holds the coefficients of `x, y, z` in `x'_i`.
    pub(super) rotation: [[f64; 3]; 3],
    /// Translation part `t`.
    pub(super) translation: [f64; 3],
}

impl SymmetryOp {
    /// Apply the operation to a fractional position.
    pub(super) fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let mut out = self.translation;
        for (row, value) in self.rotation.iter().zip(out.iter_mut()) {
            *value += row[0] * p[0] + row[1] * p[1] + row[2] * p[2];
        }
        out
    }
}

/// Map a fractional coordinate into `[0, 1)`.
///
/// `x - floor(x)` rounds to exactly `1.0` for tiny negative values such as
/// `-1e-17`, which would otherwise leak out of the half-open interval, so that
/// case is mapped to `0.0`.
pub(super) fn wrap_unit(x: f64) -> f64 {
    let wrapped = x - x.floor();
    if wrapped >= 1.0 { 0.0 } else { wrapped }
}

fn axis_index(c: char) -> Option<usize> {
    match c {
        'x' => Some(0),
        'y' => Some(1),
        'z' => Some(2),
        _ => None,
    }
}

/// Read an unsigned decimal number (`12`, `0.5`, `.25`, `3.`) starting at
/// `start`; returns the value and the index just past it.
fn parse_decimal(chars: &[char], start: usize) -> Option<(f64, usize)> {
    let mut i = start;
    let mut seen_dot = false;
    let mut text = String::new();
    while let Some(&c) = chars.get(i) {
        if c.is_ascii_digit() {
            text.push(c);
        } else if c == '.' && !seen_dot {
            seen_dot = true;
            text.push(c);
        } else {
            break;
        }
        i += 1;
    }
    if text.is_empty() || text == "." {
        return None;
    }
    text.parse::<f64>().ok().map(|value| (value, i))
}

/// Read an unsigned decimal or a fraction `n/d` starting at `start`.
fn parse_number(chars: &[char], start: usize) -> Option<(f64, usize)> {
    let (numerator, i) = parse_decimal(chars, start)?;
    if chars.get(i) == Some(&'/') {
        let (denominator, next) = parse_decimal(chars, i + 1)?;
        if denominator == 0.0 {
            return None;
        }
        return Some((numerator / denominator, next));
    }
    Some((numerator, i))
}

/// Parse one coordinate expression such as `-x+y+1/2` into the coefficients
/// of `x, y, z` and the constant term.
fn parse_component(component: &str, whole: &str) -> Result<([f64; 3], f64)> {
    let bad = |why: &str| Error::Parse(format!("CIF: invalid symmetry operation '{whole}': {why}"));
    let chars: Vec<char> = component
        .chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| c.to_ascii_lowercase())
        .collect();
    if chars.is_empty() {
        return Err(bad("empty coordinate expression"));
    }
    let mut coefficients = [0.0_f64; 3];
    let mut shift = 0.0_f64;
    let mut i = 0usize;
    let mut first_term = true;
    while i < chars.len() {
        let mut sign = 1.0_f64;
        match chars[i] {
            '+' => i += 1,
            '-' => {
                sign = -1.0;
                i += 1;
            }
            _ if first_term => {}
            _ => return Err(bad("terms must be separated by '+' or '-'")),
        }
        first_term = false;
        let Some(&c) = chars.get(i) else {
            return Err(bad("sign without a term"));
        };
        if let Some(axis) = axis_index(c) {
            coefficients[axis] += sign;
            i += 1;
        } else if c.is_ascii_digit() || c == '.' {
            let (value, next) = parse_number(&chars, i).ok_or_else(|| bad("malformed number"))?;
            i = next;
            if let Some(axis) = chars.get(i).copied().and_then(axis_index) {
                coefficients[axis] += sign * value;
                i += 1;
            } else {
                shift += sign * value;
            }
        } else {
            return Err(bad(&format!("unexpected character '{c}'")));
        }
    }
    Ok((coefficients, shift))
}

/// Parse an operation in CIF notation, e.g. `x, y, z`, `-x+1/2, y+1/2, -z`,
/// `x-y, x, z+1/6` or `1/2+x,1/2+y,z`.  Letters may be upper or lower case,
/// whitespace is ignored, and translations may be fractions or decimals.
pub(super) fn parse_symmetry_op(text: &str) -> Result<SymmetryOp> {
    let parts: Vec<&str> = text.split(',').collect();
    if parts.len() != 3 {
        return Err(Error::Parse(format!(
            "CIF: invalid symmetry operation '{text}': expected 3 comma-separated coordinates, found {}",
            parts.len()
        )));
    }
    let mut rotation = [[0.0_f64; 3]; 3];
    let mut translation = [0.0_f64; 3];
    for (row, part) in parts.iter().enumerate() {
        let (coefficients, shift) = parse_component(part, text)?;
        rotation[row] = coefficients;
        translation[row] = shift;
    }
    Ok(SymmetryOp {
        rotation,
        translation,
    })
}

/// Hash cell containing the fractional position `p` (components in `[0, 1)`).
fn cell_of(p: [f64; 3]) -> [i32; 3] {
    p.map(|x| ((x * f64::from(CELLS_PER_AXIS)) as i32).clamp(0, CELLS_PER_AXIS - 1))
}

/// Whether two fractional positions coincide within [`SITE_TOLERANCE`] on
/// every axis, taking periodic images into account (minimum-image distance).
fn same_site(a: [f64; 3], b: [f64; 3]) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| {
        let d = (x - y).abs();
        d.min(1.0 - d) < SITE_TOLERANCE
    })
}

/// Expand an asymmetric unit into the full unit cell.
///
/// Every operation is applied to every atom, the images are wrapped into
/// `[0, 1)`, and an image is dropped when an atom of the same element already
/// sits at that site (within [`SITE_TOLERANCE`], periodic images included).
/// Atoms on special positions therefore appear once, and a file that already
/// lists the whole cell is not doubled.  Occupancy and displacement parameters
/// of the first atom found at a site are kept.
pub(super) fn expand_atoms(atoms: &[CrystalAtom], ops: &[SymmetryOp]) -> Vec<CrystalAtom> {
    let mut out: Vec<CrystalAtom> = Vec::new();
    let mut cells: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
    for atom in atoms {
        for op in ops {
            let position = op.apply(atom.fractional_pos).map(wrap_unit);
            let home = cell_of(position);
            let mut duplicate = false;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbour = [
                            (home[0] + dx).rem_euclid(CELLS_PER_AXIS),
                            (home[1] + dy).rem_euclid(CELLS_PER_AXIS),
                            (home[2] + dz).rem_euclid(CELLS_PER_AXIS),
                        ];
                        let Some(members) = cells.get(&neighbour) else {
                            continue;
                        };
                        if members.iter().any(|&idx| {
                            out[idx].element == atom.element
                                && same_site(out[idx].fractional_pos, position)
                        }) {
                            duplicate = true;
                            break 'search;
                        }
                    }
                }
            }
            if duplicate {
                continue;
            }
            cells.entry(home).or_default().push(out.len());
            out.push(CrystalAtom {
                element: atom.element.clone(),
                fractional_pos: position,
                occupancy: atom.occupancy,
                b_factor: atom.b_factor,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(element: &str, pos: [f64; 3]) -> CrystalAtom {
        CrystalAtom {
            element: element.to_string(),
            fractional_pos: pos,
            occupancy: 1.0,
            b_factor: 0.0,
        }
    }

    fn op(text: &str) -> SymmetryOp {
        parse_symmetry_op(text).expect("valid operation")
    }

    fn assert_close(a: [f64; 3], b: [f64; 3]) {
        for (x, y) in a.iter().zip(b.iter()) {
            assert!((x - y).abs() < 1e-12, "{a:?} vs {b:?}");
        }
    }

    #[test]
    fn identity_and_inversion_parse() {
        let id = op("x, y, z");
        assert_eq!(
            id.rotation,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(id.translation, [0.0; 3]);
        let inv = op("-x,-y,-z");
        assert_eq!(
            inv.rotation,
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]
        );
    }

    #[test]
    fn upper_case_signs_and_whitespace_are_accepted() {
        let a = op("+X, +Y , Z");
        assert_eq!(a, op("x,y,z"));
        let b = op("  -X + 1/2 ,Y+1/2, -Z  ");
        assert_eq!(b, op("-x+1/2,y+1/2,-z"));
    }

    #[test]
    fn translation_may_come_first_and_may_be_decimal() {
        let a = op("1/2+x, 1/2+y, z");
        assert_eq!(a.translation, [0.5, 0.5, 0.0]);
        assert_eq!(a.rotation[0], [1.0, 0.0, 0.0]);
        let b = op("x+0.25, y-.5, z+3/4");
        assert_eq!(b.translation, [0.25, -0.5, 0.75]);
    }

    #[test]
    fn hexagonal_operations_mix_axes() {
        let a = op("x-y, x, z+1/6");
        assert_eq!(
            a.rotation,
            [[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert!((a.translation[2] - 1.0 / 6.0).abs() < 1e-15);
        let b = op("-x+y, -x, z+1/3");
        assert_eq!(
            b.rotation,
            [[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        );
    }

    #[test]
    fn number_directly_before_an_axis_is_a_coefficient() {
        let a = op("2x, y, z");
        assert_eq!(a.rotation[0], [2.0, 0.0, 0.0]);
        let b = op("x+2y, y, z");
        assert_eq!(b.rotation[0], [1.0, 2.0, 0.0]);
    }

    #[test]
    fn malformed_operations_are_errors() {
        for bad in [
            "", "x,y", "x,y,z,x", "x,y,w", "x y z", "x,,z", "x,y,z+", "x,y,-", "xy,y,z", "x,y,1/0",
            "x,y,1/", "x,y,.", "x,y,z+*", "x+-y,y,z",
        ] {
            let err = parse_symmetry_op(bad).unwrap_err();
            assert!(
                err.to_string().contains("symmetry operation"),
                "{bad}: {err}"
            );
        }
    }

    #[test]
    fn apply_combines_rotation_and_translation() {
        let p = [0.1, 0.2, 0.3];
        assert_close(op("x,y,z").apply(p), p);
        assert_close(op("-x,-y,-z").apply(p), [-0.1, -0.2, -0.3]);
        assert_close(op("x+1/2,y+1/2,z").apply(p), [0.6, 0.7, 0.3]);
        assert_close(op("x-y,x,z+1/6").apply(p), [-0.1, 0.1, 0.3 + 1.0 / 6.0]);
        assert_close(op("z,x,y").apply(p), [0.3, 0.1, 0.2]);
    }

    #[test]
    fn wrap_unit_stays_in_the_half_open_interval() {
        assert_eq!(wrap_unit(0.25), 0.25);
        assert_eq!(wrap_unit(1.25), 0.25);
        assert!((wrap_unit(-0.25) - 0.75).abs() < 1e-15);
        assert_eq!(wrap_unit(1.0), 0.0);
        assert_eq!(wrap_unit(0.0), 0.0);
        // x - floor(x) evaluates to exactly 1.0 here; it must map to 0.
        assert_eq!(wrap_unit(-1e-17), 0.0);
        for x in [
            -3.7,
            -1e-300,
            1.0 - f64::EPSILON / 2.0,
            2.999999,
            1e15 + 0.5,
        ] {
            let w = wrap_unit(x);
            assert!((0.0..1.0).contains(&w), "wrap_unit({x}) = {w}");
        }
    }

    #[test]
    fn centering_translations_generate_the_conventional_cell() {
        let ops = [
            op("x,y,z"),
            op("x,y+1/2,z+1/2"),
            op("x+1/2,y,z+1/2"),
            op("x+1/2,y+1/2,z"),
        ];
        let cell = expand_atoms(&[atom("Cu", [0.0, 0.0, 0.0])], &ops);
        assert_eq!(cell.len(), 4);
        let expected = [
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ];
        for (got, want) in cell.iter().zip(expected) {
            assert_eq!(got.element, "Cu");
            assert_close(got.fractional_pos, want);
        }
    }

    #[test]
    fn images_landing_on_the_same_site_are_merged() {
        // Inversion maps the origin to itself, and x+1 to the same site as x.
        let ops = [op("x,y,z"), op("-x,-y,-z"), op("x+1,y,z")];
        let cell = expand_atoms(&[atom("O", [0.0, 0.0, 0.0])], &ops);
        assert_eq!(cell.len(), 1);
    }

    #[test]
    fn periodic_images_across_the_cell_boundary_are_merged() {
        // The inversion image of x = 1e-5 wraps to 0.99999: 2e-5 away through
        // the cell boundary, which is inside the tolerance.
        let ops = [op("x,y,z"), op("-x,-y,-z")];
        let cell = expand_atoms(&[atom("Si", [1e-5, 0.0, 0.0])], &ops);
        assert_eq!(cell.len(), 1);
        let apart = expand_atoms(&[atom("Si", [1e-3, 0.0, 0.0])], &ops);
        assert_eq!(apart.len(), 2);
    }

    #[test]
    fn different_elements_at_one_site_are_kept() {
        let ops = [op("x,y,z"), op("x,y,z")];
        let atoms = [
            atom("Fe", [0.25, 0.25, 0.25]),
            atom("Ni", [0.25, 0.25, 0.25]),
        ];
        let cell = expand_atoms(&atoms, &ops);
        assert_eq!(cell.len(), 2);
        assert_eq!(cell[0].element, "Fe");
        assert_eq!(cell[1].element, "Ni");
    }

    #[test]
    fn an_already_expanded_cell_is_not_doubled() {
        let ops = [op("x,y,z"), op("x+1/2,y+1/2,z+1/2")];
        let listed = [atom("Fe", [0.0, 0.0, 0.0]), atom("Fe", [0.5, 0.5, 0.5])];
        assert_eq!(expand_atoms(&listed, &ops).len(), 2);
    }

    #[test]
    fn expansion_keeps_occupancy_and_displacement_of_the_source_atom() {
        let mut source = atom("Na", [0.1, 0.2, 0.3]);
        source.occupancy = 0.5;
        source.b_factor = 1.25;
        let cell = expand_atoms(&[source], &[op("x,y,z"), op("-x,-y,-z")]);
        assert_eq!(cell.len(), 2);
        for a in &cell {
            assert_eq!(a.occupancy, 0.5);
            assert_eq!(a.b_factor, 1.25);
        }
        assert_close(cell[1].fractional_pos, [0.9, 0.8, 0.7]);
    }

    #[test]
    fn large_expansions_stay_fast_and_correct() {
        // 8 x 8 x 8 grid of distinct sites, expanded with 4 translations that
        // each map the grid onto itself: no atom may be duplicated.
        let mut atoms = Vec::new();
        for i in 0..8 {
            for j in 0..8 {
                for k in 0..8 {
                    atoms.push(atom("Ar", [i as f64 / 8.0, j as f64 / 8.0, k as f64 / 8.0]));
                }
            }
        }
        let ops = [
            op("x,y,z"),
            op("x+1/8,y,z"),
            op("x,y+1/8,z"),
            op("x,y,z+1/8"),
        ];
        assert_eq!(expand_atoms(&atoms, &ops).len(), 512);
    }
}
