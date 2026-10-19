// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rest-state geometry of linear tetrahedra for the real-time neural soft body.
//!
//! A linear tetrahedron with rest corners `X0..X3` interpolates the
//! displacement field with the barycentric shape functions `N0..N3`.  With the
//! rest edge matrix `Dm = [X1 - X0 | X2 - X0 | X3 - X0]` (edge vectors as
//! columns):
//!
//! * the gradients of the shape functions are the rows of `Dm^-1`:
//!   `grad N_a = (Dm^-1)^T e_(a-1)` for `a = 1..3` and
//!   `grad N_0 = -(grad N_1 + grad N_2 + grad N_3)`,
//! * the displacement gradient is `grad u = Du Dm^-1` with
//!   `Du = [u1 - u0 | u2 - u0 | u3 - u0]`, and the small strain is
//!   `eps = sym(grad u)`,
//! * a (symmetric) stress `sigma` acting on the element exerts the nodal
//!   forces `f_a = -V0 sigma grad N_a` with `V0 = |det Dm| / 6`.  These sum
//!   to zero, so the element conserves linear momentum, and they vanish for a
//!   rigid translation.

/// Elements whose rest volume is below this fraction of the cube of their
/// longest rest edge are treated as degenerate.
const DEGENERATE_RELATIVE_VOLUME: f64 = 1e-12;

/// Rest-state data of one linear tetrahedron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetRestShape {
    /// Inverse `Dm^-1` of the rest edge matrix `Dm = [X1 - X0 | X2 - X0 | X3 - X0]`
    /// (row-major).  Row `a` is the gradient of the shape function `N_(a+1)`.
    pub dm_inv: [[f64; 3]; 3],
    /// Rest volume `V0 = |det Dm| / 6`.
    pub volume: f64,
}

impl TetRestShape {
    /// Build the rest data of the tetrahedron `elem` from the rest positions.
    ///
    /// Returns `None` if a corner index is out of range or the tetrahedron is
    /// degenerate (its rest volume vanishes relative to its edge lengths).
    pub fn from_rest(rest: &[[f64; 3]], elem: &[usize; 4]) -> Option<Self> {
        let x0 = *rest.get(elem[0])?;
        let mut edges = [[0.0f64; 3]; 3];
        for (edge, &corner) in edges.iter_mut().zip(elem[1..].iter()) {
            let x = *rest.get(corner)?;
            *edge = [x[0] - x0[0], x[1] - x0[1], x[2] - x0[2]];
        }
        let dm = [
            [edges[0][0], edges[1][0], edges[2][0]],
            [edges[0][1], edges[1][1], edges[2][1]],
            [edges[0][2], edges[1][2], edges[2][2]],
        ];
        let det = determinant(&dm);
        let longest = edges
            .iter()
            .map(|e| (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt())
            .fold(0.0f64, f64::max);
        // Written positively so that a NaN determinant is rejected as well.
        let non_degenerate = det.abs() > DEGENERATE_RELATIVE_VOLUME * longest.powi(3);
        if !non_degenerate {
            return None;
        }
        Some(Self {
            dm_inv: inverse(&dm, det),
            volume: det.abs() / 6.0,
        })
    }

    /// Displacement gradient `grad u = Du Dm^-1` (`[i][j] = d u_i / d X_j`)
    /// from the four nodal displacements.
    pub fn displacement_gradient(&self, u: &[[f64; 3]; 4]) -> [[f64; 3]; 3] {
        let mut gradient = [[0.0f64; 3]; 3];
        for (i, gradient_row) in gradient.iter_mut().enumerate() {
            for (j, entry) in gradient_row.iter_mut().enumerate() {
                *entry = (0..3)
                    .map(|c| (u[c + 1][i] - u[0][i]) * self.dm_inv[c][j])
                    .sum();
            }
        }
        gradient
    }

    /// Small strain `eps = sym(grad u)` in Voigt order
    /// `[xx, yy, zz, xy, xz, yz]` (tensor shear components, not engineering
    /// shear).
    pub fn voigt_strain(&self, u: &[[f64; 3]; 4]) -> [f64; 6] {
        let g = self.displacement_gradient(u);
        [
            g[0][0],
            g[1][1],
            g[2][2],
            0.5 * (g[0][1] + g[1][0]),
            0.5 * (g[0][2] + g[2][0]),
            0.5 * (g[1][2] + g[2][1]),
        ]
    }

    /// Nodal forces `f_a = -V0 sigma grad N_a` exerted on the four corners by
    /// the Voigt stress `[xx, yy, zz, xy, xz, yz]`.
    ///
    /// `f_1..f_3` are the columns of `-V0 sigma Dm^-T` and
    /// `f_0 = -(f_1 + f_2 + f_3)`, so the four forces sum to zero.
    pub fn nodal_forces(&self, stress: &[f64; 6]) -> [[f64; 3]; 4] {
        let sigma = [
            [stress[0], stress[3], stress[4]],
            [stress[3], stress[1], stress[5]],
            [stress[4], stress[5], stress[2]],
        ];
        let mut forces = [[0.0f64; 3]; 4];
        for a in 0..3 {
            let grad_n = self.dm_inv[a];
            for (i, sigma_row) in sigma.iter().enumerate() {
                forces[a + 1][i] = -self.volume
                    * (sigma_row[0] * grad_n[0]
                        + sigma_row[1] * grad_n[1]
                        + sigma_row[2] * grad_n[2]);
            }
        }
        forces[0] = [
            -(forces[1][0] + forces[2][0] + forces[3][0]),
            -(forces[1][1] + forces[2][1] + forces[3][1]),
            -(forces[1][2] + forces[2][2] + forces[3][2]),
        ];
        forces
    }
}

/// Determinant of a 3x3 matrix.
fn determinant(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Inverse of a 3x3 matrix with the given (non-zero) determinant.
fn inverse(m: &[[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = 1.0 / det;
    [
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
    ]
}
