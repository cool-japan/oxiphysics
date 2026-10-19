// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! 3x3 linear-algebra helpers shared by the crystallography readers.
//!
//! Lattices are stored as row matrices: `lattice[i]` is the i-th basis vector,
//! so a Cartesian position is `r = f . L` for the row vector of fractional
//! coordinates `f`, and `f = r . L^-1` in the other direction.

/// Lattices whose determinant (cell volume in cubic Angstrom) has a smaller
/// magnitude than this are treated as singular.
pub(super) const MIN_LATTICE_DET: f64 = 1e-12;

/// Determinant of a 3x3 matrix.
pub(super) fn determinant(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Inverse of a 3x3 matrix via its adjugate.
///
/// Returns `None` when the determinant is not finite or its magnitude is below
/// [`MIN_LATTICE_DET`].
pub(super) fn inverse(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = determinant(m);
    if !det.is_finite() || det.abs() < MIN_LATTICE_DET {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}

/// Row vector times matrix: `r . m`.
pub(super) fn row_times_matrix(r: [f64; 3], m: &[[f64; 3]; 3]) -> [f64; 3] {
    [
        r[0] * m[0][0] + r[1] * m[1][0] + r[2] * m[2][0],
        r[0] * m[0][1] + r[1] * m[1][1] + r[2] * m[2][1],
        r[0] * m[0][2] + r[1] * m[1][2] + r[2] * m[2][2],
    ]
}

/// Convert a Cartesian position (Angstrom) to fractional coordinates of the
/// lattice whose inverse is `lattice_inverse`.
pub(super) fn cartesian_to_fractional(r: [f64; 3], lattice_inverse: &[[f64; 3]; 3]) -> [f64; 3] {
    row_times_matrix(r, lattice_inverse)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRICLINIC: [[f64; 3]; 3] = [[4.1, 0.3, -0.2], [0.7, 5.2, 0.4], [-0.5, 0.9, 6.3]];

    #[test]
    fn determinant_of_diagonal_and_triclinic_matrices() {
        let diag = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        assert!((determinant(&diag) - 24.0).abs() < 1e-14);
        // Cofactor expansion along the first row, written out by hand.
        let expected = 4.1 * (5.2 * 6.3 - 0.4 * 0.9) - 0.3 * (0.7 * 6.3 - 0.4 * -0.5)
            + -0.2 * (0.7 * 0.9 - 5.2 * -0.5);
        assert!((determinant(&TRICLINIC) - expected).abs() < 1e-12);
    }

    #[test]
    fn inverse_times_matrix_is_identity() {
        let inv = inverse(&TRICLINIC).expect("triclinic lattice is regular");
        for (i, lattice_row) in TRICLINIC.iter().enumerate() {
            let row = row_times_matrix(*lattice_row, &inv);
            for (j, value) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 1e-13, "({i},{j}) = {value}");
            }
        }
    }

    #[test]
    fn singular_and_non_finite_matrices_have_no_inverse() {
        let singular = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 1.0]];
        assert!(inverse(&singular).is_none());
        let nan = [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(inverse(&nan).is_none());
        let tiny = [[1e-5, 0.0, 0.0], [0.0, 1e-5, 0.0], [0.0, 0.0, 1e-5]];
        assert!(inverse(&tiny).is_none());
    }

    #[test]
    fn cartesian_conversion_inverts_the_fractional_to_cartesian_map() {
        let frac = [0.25, 0.5, 0.75];
        let cart = row_times_matrix(frac, &TRICLINIC);
        let inv = inverse(&TRICLINIC).expect("regular");
        let back = cartesian_to_fractional(cart, &inv);
        for (got, want) in back.iter().zip(frac.iter()) {
            assert!((got - want).abs() < 1e-13);
        }
    }
}
