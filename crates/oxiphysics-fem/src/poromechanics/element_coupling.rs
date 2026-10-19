// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Element-level Biot coupling integral of the u-p formulation.
//!
//! The coupling matrix of an element `Ω_e` with equal-order interpolation of
//! displacement and pore pressure (`N_p = N_u = N`) is
//!
//! ```text
//! Q_e = α ∫_{Ω_e} Bᵀ m N dV,        m = (1, 1, 1, 0, 0, 0)ᵀ,
//! ```
//!
//! and since `Bᵀ m` collects the gradient of the shape functions, the entry for
//! the displacement component `i` of node `a` and the pressure of node `b` is
//!
//! ```text
//! Q_e[3a + i, b] = α ∫ (∂N_a/∂x_i) N_b dV.
//! ```
//!
//! The integrals are evaluated exactly for the 4-node tetrahedron (constant
//! gradients, `∫ N_b dV = V/4`) and with the 2×2×2 Gauss rule for the 8-node
//! trilinear hexahedron.  The resulting `3n × n` row-major matrix can be added
//! to a [`ConsolidationSolver`](super::ConsolidationSolver) through
//! [`ConsolidationSolver::add_element_block`](super::ConsolidationSolver::add_element_block).

use crate::error::{Error, Result};

/// Abscissa `1/√3` of the two-point Gauss-Legendre rule on `[-1, 1]` (weights 1).
const GAUSS_2: f64 = 0.577_350_269_189_625_8;

/// Natural coordinates of the hexahedron nodes: the corners of the bottom face
/// counter-clockwise, then the corners of the top face in the same order.
const HEX8_NATURAL: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];

/// Smallest accepted ratio `|det J| / (|J₁|·|J₂|·|J₃|)` of a Jacobian with
/// row norms `|J_k|` (1 for an orthogonal map, 0 for a flat element).
const MIN_JACOBIAN_RATIO: f64 = 1.0e-12;

/// Cofactor matrix of a 3×3 matrix.
fn cofactors(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
        ],
        [
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
        ],
        [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ]
}

/// Inverse and (signed) determinant of the Jacobian `J[i][k] = ∂x_i/∂ξ_k`.
///
/// # Errors
/// [`Error::InvalidInput`] if the map is numerically singular.
fn invert_jacobian(j: &[[f64; 3]; 3], what: &str) -> Result<([[f64; 3]; 3], f64)> {
    let cof = cofactors(j);
    let det = j[0][0] * cof[0][0] + j[0][1] * cof[0][1] + j[0][2] * cof[0][2];
    let row_norm = |row: &[f64; 3]| row.iter().map(|v| v * v).sum::<f64>().sqrt();
    let hadamard = row_norm(&j[0]) * row_norm(&j[1]) * row_norm(&j[2]);
    if !det.is_finite() || det.abs() <= MIN_JACOBIAN_RATIO * hadamard {
        return Err(Error::InvalidInput(format!(
            "{what} is degenerate (zero or negligible volume)"
        )));
    }
    let mut inverse = [[0.0; 3]; 3];
    for (i, row) in inverse.iter_mut().enumerate() {
        for (k, entry) in row.iter_mut().enumerate() {
            *entry = cof[k][i] / det;
        }
    }
    Ok((inverse, det))
}

/// Gradients `∇N_a` of the linear tetrahedron shape functions (constant over
/// the element) and the element volume.
///
/// The node order is arbitrary; a mirrored (negatively oriented) element gives
/// the same gradients and a positive volume.
///
/// # Errors
/// [`Error::InvalidInput`] for a degenerate (flat) tetrahedron.
pub(super) fn tet4_gradients(coords: &[[f64; 3]; 4]) -> Result<([[f64; 3]; 4], f64)> {
    // J[i][k] = x_{k+1, i} - x_{0, i}; N_{k+1} = ξ_k, N_0 = 1 - ξ_0 - ξ_1 - ξ_2.
    let mut jacobian = [[0.0; 3]; 3];
    for (i, row) in jacobian.iter_mut().enumerate() {
        for (k, entry) in row.iter_mut().enumerate() {
            *entry = coords[k + 1][i] - coords[0][i];
        }
    }
    let (inverse, det) = invert_jacobian(&jacobian, "the tetrahedron")?;
    let mut gradients = [[0.0; 3]; 4];
    gradients[1..4].copy_from_slice(&inverse);
    // Partition of unity: the four gradients add up to zero.
    gradients[0] = std::array::from_fn(|i| -(gradients[1][i] + gradients[2][i] + gradients[3][i]));
    Ok((gradients, det.abs() / 6.0))
}

/// Coupling matrix of a 4-node tetrahedron (`12 × 4`, row-major).
fn tet4_coupling(coords: &[[f64; 3]; 4], alpha: f64) -> Result<Vec<f64>> {
    let (gradients, volume) = tet4_gradients(coords)?;
    // ∫ N_b dV = V/4 for every node b.
    let weight = alpha * volume / 4.0;
    let mut q = vec![0.0; 12 * 4];
    for (a, gradient) in gradients.iter().enumerate() {
        for (i, &g) in gradient.iter().enumerate() {
            for b in 0..4 {
                q[(3 * a + i) * 4 + b] = weight * g;
            }
        }
    }
    Ok(q)
}

/// Trilinear shape functions and their natural derivatives at `(ξ, η, ζ)`.
fn hex8_shape(natural: [f64; 3]) -> ([f64; 8], [[f64; 3]; 8]) {
    let mut n = [0.0; 8];
    let mut dn = [[0.0; 3]; 8];
    for (a, node) in HEX8_NATURAL.iter().enumerate() {
        let f = [
            1.0 + node[0] * natural[0],
            1.0 + node[1] * natural[1],
            1.0 + node[2] * natural[2],
        ];
        n[a] = 0.125 * f[0] * f[1] * f[2];
        dn[a] = [
            0.125 * node[0] * f[1] * f[2],
            0.125 * f[0] * node[1] * f[2],
            0.125 * f[0] * f[1] * node[2],
        ];
    }
    (n, dn)
}

/// Coupling matrix of an 8-node hexahedron (`24 × 8`, row-major) by 2×2×2 Gauss
/// quadrature.
fn hex8_coupling(coords: &[[f64; 3]; 8], alpha: f64) -> Result<Vec<f64>> {
    let mut q = vec![0.0; 24 * 8];
    let mut orientation = 0.0_f64;
    for &xi in &[-GAUSS_2, GAUSS_2] {
        for &eta in &[-GAUSS_2, GAUSS_2] {
            for &zeta in &[-GAUSS_2, GAUSS_2] {
                let (n, dn) = hex8_shape([xi, eta, zeta]);
                // J[i][k] = Σ_a x_{a,i} ∂N_a/∂ξ_k
                let mut jacobian = [[0.0; 3]; 3];
                for (x, d) in coords.iter().zip(dn.iter()) {
                    for i in 0..3 {
                        for k in 0..3 {
                            jacobian[i][k] += x[i] * d[k];
                        }
                    }
                }
                let (inverse, det) = invert_jacobian(&jacobian, "the hexahedron")?;
                if orientation == 0.0 {
                    orientation = det.signum();
                } else if det.signum() != orientation {
                    return Err(Error::InvalidInput(
                        "the hexahedron is twisted (the Jacobian changes sign)".to_string(),
                    ));
                }
                let weight = alpha * det.abs();
                for (a, d) in dn.iter().enumerate() {
                    for i in 0..3 {
                        // (∇N_a)_i = Σ_k ∂N_a/∂ξ_k ∂ξ_k/∂x_i
                        let gradient: f64 = (0..3).map(|k| d[k] * inverse[k][i]).sum();
                        for (b, &n_b) in n.iter().enumerate() {
                            q[(3 * a + i) * 8 + b] += weight * gradient * n_b;
                        }
                    }
                }
            }
        }
    }
    Ok(q)
}

/// Row-major `3n × n` coupling matrix `Q_e = α ∫ Bᵀ m N dV` of an element given
/// by the coordinates of its nodes (4 nodes: linear tetrahedron, 8 nodes:
/// trilinear hexahedron in the bottom-face-then-top-face order of
/// [`HexahedralElement`](crate::element::HexahedralElement)).
///
/// # Errors
/// [`Error::InvalidInput`] for a node count other than 4 or 8, a non-finite
/// coordinate or `alpha`, and degenerate or twisted elements.
pub(super) fn coupling_matrix(node_coords: &[[f64; 3]], alpha: f64) -> Result<Vec<f64>> {
    if !alpha.is_finite() {
        return Err(Error::InvalidInput(format!(
            "the Biot coefficient must be finite, got {alpha}"
        )));
    }
    if node_coords.iter().flatten().any(|c| !c.is_finite()) {
        return Err(Error::InvalidInput(
            "node coordinates must be finite".to_string(),
        ));
    }
    if let Ok(tet) = <[[f64; 3]; 4]>::try_from(node_coords) {
        return tet4_coupling(&tet, alpha);
    }
    if let Ok(hex) = <[[f64; 3]; 8]>::try_from(node_coords) {
        return hex8_coupling(&hex, alpha);
    }
    Err(Error::InvalidInput(format!(
        "the coupling matrix is implemented for 4-node tetrahedra and 8-node hexahedra, \
         got {} nodes",
        node_coords.len()
    )))
}

#[cfg(test)]
mod tests {
    use super::super::{
        BiotBlock, ConsolidationSolver, PorousElement, TerzaghiConsolidation, pore_pressure_ratio,
    };
    use super::*;
    use crate::constitutive::LinearElasticMaterial;
    use crate::element::{HexahedralElement, LinearTetrahedron};
    use crate::mesh::TetrahedralMesh;
    use oxiphysics_core::math::Vec3;

    const REFERENCE_TET: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    /// A sheared, rotated tetrahedron without any symmetry.
    fn general_tet() -> [[f64; 3]; 4] {
        [
            [0.3, -0.2, 0.1],
            [1.4, 0.1, 0.3],
            [0.2, 1.1, -0.2],
            [0.5, 0.4, 1.3],
        ]
    }

    const UNIT_CUBE: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];

    /// A hexahedron with non-parallel faces.
    fn distorted_hex() -> [[f64; 3]; 8] {
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.1, 0.0],
            [1.1, 1.3, 0.1],
            [-0.1, 1.0, 0.0],
            [0.1, -0.1, 1.0],
            [1.3, 0.0, 1.2],
            [1.2, 1.2, 1.1],
            [0.0, 1.1, 0.9],
        ]
    }

    /// `Σ_{a,i} Q[3a + i, b] u_{a,i}` for the displacement field `u(x) = f(x)`.
    fn coupled_pressure_load(
        q: &[f64],
        coords: &[[f64; 3]],
        field: impl Fn(&[f64; 3]) -> [f64; 3],
    ) -> Vec<f64> {
        let n = coords.len();
        (0..n)
            .map(|b| {
                coords
                    .iter()
                    .enumerate()
                    .map(|(a, x)| {
                        let u = field(x);
                        (0..3).map(|i| q[(3 * a + i) * n + b] * u[i]).sum::<f64>()
                    })
                    .sum()
            })
            .collect()
    }

    #[test]
    fn reference_tetrahedron_matches_hand_computed_entries() {
        // V = 1/6, ∇N_0 = (-1,-1,-1), ∇N_k = e_k:  Q[3a + i, b] = α ∇N_a,i / 24.
        let alpha = 0.8;
        let q = coupling_matrix(&REFERENCE_TET, alpha).expect("coupling");
        assert_eq!(q.len(), 48);
        let gradients = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        for (a, gradient) in gradients.iter().enumerate() {
            for (i, &g) in gradient.iter().enumerate() {
                for b in 0..4 {
                    let expected = alpha * g / 24.0;
                    assert!(
                        (q[(3 * a + i) * 4 + b] - expected).abs() < 1e-15,
                        "Q[{a},{i};{b}] = {} vs {expected}",
                        q[(3 * a + i) * 4 + b]
                    );
                }
            }
        }
    }

    #[test]
    fn tetrahedron_gradients_match_the_solid_element_b_matrix() {
        let coords = general_tet();
        let (gradients, volume) = tet4_gradients(&coords).expect("gradients");
        let nodes = coords.map(|c| Vec3::new(c[0], c[1], c[2]));
        let (b, signed_volume) = LinearTetrahedron::b_matrix(&nodes);
        assert!((volume - signed_volume.abs()).abs() < 1e-14);
        for (a, gradient) in gradients.iter().enumerate() {
            // Rows 0..3 of B hold ∂N_a/∂x, ∂N_a/∂y, ∂N_a/∂z on the diagonal blocks.
            let from_b = [b[3 * a], b[12 + 3 * a + 1], b[24 + 3 * a + 2]];
            for i in 0..3 {
                assert!(
                    (gradient[i] - from_b[i]).abs() < 1e-13,
                    "node {a}, component {i}: {} vs {}",
                    gradient[i],
                    from_b[i]
                );
            }
        }
        // Linear fields are reproduced: Σ_a x_{a,j} ∂N_a/∂x_i = δ_ij.
        for i in 0..3 {
            for j in 0..3 {
                let value: f64 = coords
                    .iter()
                    .zip(gradients.iter())
                    .map(|(x, g)| x[j] * g[i])
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 1e-13);
            }
        }
    }

    #[test]
    fn hexahedron_shape_functions_follow_the_node_order_of_the_solid_element() {
        for natural in [[0.3, -0.7, 0.1], [-1.0, 1.0, 0.5], [0.0, 0.0, 0.0]] {
            let (n, dn) = hex8_shape(natural);
            let reference_n =
                HexahedralElement::shape_functions(natural[0], natural[1], natural[2]);
            let (dxi, deta, dzeta) =
                HexahedralElement::shape_function_derivatives(natural[0], natural[1], natural[2]);
            for a in 0..8 {
                assert!((n[a] - reference_n[a]).abs() < 1e-15, "N_{a}");
                assert!((dn[a][0] - dxi[a]).abs() < 1e-15, "dN_{a}/dxi");
                assert!((dn[a][1] - deta[a]).abs() < 1e-15, "dN_{a}/deta");
                assert!((dn[a][2] - dzeta[a]).abs() < 1e-15, "dN_{a}/dzeta");
            }
            // Partition of unity.
            assert!((n.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn rigid_translation_causes_no_coupling() {
        let translation = |_: &[f64; 3]| [0.7, -1.3, 2.1];
        let tet = general_tet();
        let q = coupling_matrix(&tet, 1.0).expect("tet");
        for load in coupled_pressure_load(&q, &tet, translation) {
            assert!(load.abs() < 1e-13, "tet: {load}");
        }
        let hex = distorted_hex();
        let q = coupling_matrix(&hex, 1.0).expect("hex");
        for load in coupled_pressure_load(&q, &hex, translation) {
            assert!(load.abs() < 1e-13, "hex: {load}");
        }
    }

    #[test]
    fn uniform_dilation_loads_every_node_with_its_share_of_the_volume() {
        // u = ε x has divergence 3ε, so (Qᵀu)_b = 3 α ε ∫ N_b dV.
        let (alpha, eps) = (0.9, 1e-3);
        let dilation = |x: &[f64; 3]| [eps * x[0], eps * x[1], eps * x[2]];

        let tet = general_tet();
        let (_, volume) = tet4_gradients(&tet).expect("gradients");
        let q = coupling_matrix(&tet, alpha).expect("tet");
        for load in coupled_pressure_load(&q, &tet, dilation) {
            let expected = 3.0 * alpha * eps * volume / 4.0;
            assert!((load - expected).abs() < 1e-15, "tet: {load} vs {expected}");
        }

        // Unit cube: ∫ N_b dV = 1/8.
        let q = coupling_matrix(&UNIT_CUBE, alpha).expect("cube");
        for load in coupled_pressure_load(&q, &UNIT_CUBE, dilation) {
            let expected = 3.0 * alpha * eps / 8.0;
            assert!(
                (load - expected).abs() < 1e-15,
                "cube: {load} vs {expected}"
            );
        }

        // Distorted cell: the loads add up to 3 α ε V with V from the same quadrature.
        let hex = distorted_hex();
        let q = coupling_matrix(&hex, alpha).expect("hex");
        let total: f64 = coupled_pressure_load(&q, &hex, dilation).iter().sum();
        let expected = 3.0 * alpha * eps * HexahedralElement::volume(&hex);
        assert!(
            (total - expected).abs() < 1e-13,
            "distorted hex: {total} vs {expected}"
        );
    }

    #[test]
    fn coupling_depends_on_the_geometry_not_on_the_node_count_alone() {
        // Q scales with (length)^2: gradients ~ 1/h and volume ~ h^3.
        let scaled = |c: &[[f64; 3]], s: f64| -> Vec<[f64; 3]> {
            c.iter().map(|x| [s * x[0], s * x[1], s * x[2]]).collect()
        };
        for coords in [general_tet().to_vec(), distorted_hex().to_vec()] {
            let q1 = coupling_matrix(&coords, 1.0).expect("unit size");
            let q2 = coupling_matrix(&scaled(&coords, 2.0), 1.0).expect("double size");
            assert!(q1.iter().any(|v| v.abs() > 1e-3));
            for (a, b) in q1.iter().zip(q2.iter()) {
                assert!((4.0 * a - b).abs() < 1e-12, "{a} vs {b}");
            }
        }
        // A stretched tetrahedron couples differently from the reference one.
        let stretched = [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let reference = coupling_matrix(&REFERENCE_TET, 1.0).expect("reference");
        let other = coupling_matrix(&stretched, 1.0).expect("stretched");
        assert!(
            reference
                .iter()
                .zip(other.iter())
                .any(|(a, b)| (a - b).abs() > 1e-3)
        );
    }

    #[test]
    fn coupling_is_translation_invariant_and_rotates_like_a_vector() {
        let shift = [4.0, -2.5, 0.75];
        // Rotation by 90 degrees about z:  (x, y, z) -> (-y, x, z).
        let rotate = |x: &[f64; 3]| [-x[1], x[0], x[2]];
        for coords in [general_tet().to_vec(), distorted_hex().to_vec()] {
            let n = coords.len();
            let q = coupling_matrix(&coords, 1.0).expect("coupling");
            let moved: Vec<[f64; 3]> = coords
                .iter()
                .map(|x| [x[0] + shift[0], x[1] + shift[1], x[2] + shift[2]])
                .collect();
            let q_moved = coupling_matrix(&moved, 1.0).expect("shifted");
            for (a, b) in q.iter().zip(q_moved.iter()) {
                assert!((a - b).abs() < 1e-12);
            }
            let rotated: Vec<[f64; 3]> = coords.iter().map(rotate).collect();
            let q_rot = coupling_matrix(&rotated, 1.0).expect("rotated");
            // Row (a, i) of the rotated matrix is R applied to the rows (a, ·).
            for a in 0..n {
                for b in 0..n {
                    let old = [
                        q[(3 * a) * n + b],
                        q[(3 * a + 1) * n + b],
                        q[(3 * a + 2) * n + b],
                    ];
                    let expected = rotate(&old);
                    for i in 0..3 {
                        let got = q_rot[(3 * a + i) * n + b];
                        assert!(
                            (got - expected[i]).abs() < 1e-12,
                            "{got} vs {}",
                            expected[i]
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn mirrored_element_gives_the_same_coupling() {
        // Swapping two nodes reverses the orientation but not the integral of the
        // (relabelled) basis functions: Q(permuted nodes) is Q with permuted rows/columns.
        let tet = general_tet();
        let swapped = [tet[0], tet[2], tet[1], tet[3]];
        let q = coupling_matrix(&tet, 1.0).expect("tet");
        let q_swapped = coupling_matrix(&swapped, 1.0).expect("swapped");
        let permutation = [0, 2, 1, 3];
        for (new_a, &old_a) in permutation.iter().enumerate() {
            for i in 0..3 {
                for (new_b, &old_b) in permutation.iter().enumerate() {
                    let a = q_swapped[(3 * new_a + i) * 4 + new_b];
                    let b = q[(3 * old_a + i) * 4 + old_b];
                    assert!((a - b).abs() < 1e-13);
                }
            }
        }
    }

    #[test]
    fn porous_element_wraps_the_integral_and_checks_its_input() {
        let k = [1e-12, 0.0, 0.0, 0.0, 1e-12, 0.0, 0.0, 0.0, 1e-12];
        let element = PorousElement::new(vec![10, 11, 12, 13], k, 0.3, 1e7);
        let q = element
            .compute_coupling_matrix(&REFERENCE_TET, 0.5)
            .expect("coupling");
        assert_eq!(q, coupling_matrix(&REFERENCE_TET, 0.5).expect("direct"));
        // Coordinates must match the node list.
        assert!(matches!(
            element.compute_coupling_matrix(&UNIT_CUBE, 1.0),
            Err(Error::InvalidInput(_))
        ));
        // Unsupported topology.
        let five = PorousElement::new(vec![0, 1, 2, 3, 4], k, 0.3, 1e7);
        let coords = [
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0; 3],
        ];
        assert!(matches!(
            five.compute_coupling_matrix(&coords, 1.0),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn degenerate_twisted_and_non_finite_input_is_rejected() {
        let flat = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert!(matches!(
            coupling_matrix(&flat, 1.0),
            Err(Error::InvalidInput(_))
        ));
        // Two top nodes exchanged: the cell folds over itself.
        let mut twisted = UNIT_CUBE;
        twisted.swap(6, 7);
        assert!(matches!(
            coupling_matrix(&twisted, 1.0),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            coupling_matrix(&REFERENCE_TET, f64::NAN),
            Err(Error::InvalidInput(_))
        ));
        let mut bad = REFERENCE_TET;
        bad[2][1] = f64::INFINITY;
        assert!(matches!(
            coupling_matrix(&bad, 1.0),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            coupling_matrix(&REFERENCE_TET[..3], 1.0),
            Err(Error::InvalidInput(_))
        ));
    }

    /// One-dimensional consolidation of a column discretised with tetrahedra:
    /// element stiffness, permeability and the geometric coupling matrix are all
    /// assembled from the node coordinates and compared with Terzaghi's solution.
    #[test]
    fn terzaghi_column_of_tetrahedra_matches_the_analytic_solution() {
        let (height, n_cells) = (2.0, 20);
        let (young, poisson, mobility, load) = (6.0e6, 0.3, 2.0e-9, 1.0e4);
        // The column axis is x; all lateral displacements are suppressed, so the
        // skeleton is loaded in uniaxial strain with the oedometric modulus.
        let mesh = TetrahedralMesh::generate_beam(height, 1.0, 1.0, n_cells, 1, 1);
        let d_matrix = LinearElasticMaterial::new(young, poisson).constitutive_matrix();
        let oedometric = young * (1.0 - poisson) / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        let n_nodes = mesh.num_nodes();
        let mut solver = ConsolidationSolver::new(3 * n_nodes, n_nodes);
        let permeability = [mobility, 0.0, 0.0, 0.0, mobility, 0.0, 0.0, 0.0, mobility];

        for connectivity in &mesh.elements {
            let nodes = connectivity.map(|n| mesh.nodes[n]);
            let coords = connectivity.map(|n| {
                let p = mesh.nodes[n];
                [p.x, p.y, p.z]
            });
            let u_dofs: Vec<usize> = connectivity
                .iter()
                .flat_map(|&n| (0..3).map(move |c| 3 * n + c))
                .collect();
            let stiffness: Vec<f64> = LinearTetrahedron::element_stiffness(&nodes, &d_matrix)
                .iter()
                .flat_map(|row| row.iter().copied())
                .collect();
            solver
                .add_element_block(BiotBlock::Stiffness, &u_dofs, &u_dofs, &stiffness)
                .expect("stiffness");

            let element = PorousElement::new(connectivity.to_vec(), permeability, 0.3, young);
            let coupling = element
                .compute_coupling_matrix(&coords, 1.0)
                .expect("coupling");
            solver
                .add_element_block(BiotBlock::Coupling, &u_dofs, connectivity, &coupling)
                .expect("coupling block");

            let (gradients, volume) = tet4_gradients(&coords).expect("gradients");
            let mut flow = Vec::with_capacity(16);
            for ga in &gradients {
                for gb in &gradients {
                    let dot: f64 = ga.iter().zip(gb.iter()).map(|(x, y)| x * y).sum();
                    flow.push(mobility * volume * dot);
                }
            }
            solver
                .add_element_block(BiotBlock::Permeability, connectivity, connectivity, &flow)
                .expect("flow block");
        }

        // Boundary conditions: base fixed, lateral motion suppressed, top drained.
        let on_top = |n: usize| (mesh.nodes[n].x - height).abs() < 1e-9;
        for n in 0..n_nodes {
            solver.fix_displacement(3 * n + 1, 0.0).expect("u_y");
            solver.fix_displacement(3 * n + 2, 0.0).expect("u_z");
            if mesh.nodes[n].x.abs() < 1e-9 {
                solver.fix_displacement(3 * n, 0.0).expect("base");
            }
            if on_top(n) {
                solver.fix_pressure(n, 0.0).expect("drain");
            }
        }
        // Uniform pressure on the top face: consistent nodal forces (area / 3 per
        // face triangle and node).
        let mut f_u = vec![0.0; 3 * n_nodes];
        for connectivity in &mesh.elements {
            for face in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
                let ids = face.map(|k| connectivity[k]);
                if ids.iter().all(|&n| on_top(n)) {
                    let p = ids.map(|n| mesh.nodes[n]);
                    let area = 0.5 * (p[1] - p[0]).cross(&(p[2] - p[0])).norm();
                    for &n in &ids {
                        f_u[3 * n] -= load * area / 3.0;
                    }
                }
            }
        }
        let total_force: f64 = f_u.iter().sum();
        assert!((total_force + load).abs() < 1e-9 * load, "{total_force}");

        let consolidation = mobility * oedometric;
        let tv_end = 0.2;
        let t_end = tv_end * height * height / consolidation;
        let n_steps = 400;
        let dt = t_end / n_steps as f64;
        let f_p = vec![0.0; n_nodes];
        for _ in 0..n_steps {
            solver.step(dt, &f_u, &f_p).expect("Biot step");
        }

        let final_settlement = load * height / oedometric;
        let top: Vec<usize> = (0..n_nodes).filter(|&n| on_top(n)).collect();
        assert_eq!(top.len(), 4);
        let mean_settlement =
            -top.iter().map(|&n| solver.u_dofs[3 * n]).sum::<f64>() / top.len() as f64;
        let degree = mean_settlement / final_settlement;
        let exact =
            TerzaghiConsolidation::new(height, consolidation).degree_of_consolidation(t_end);
        assert!(
            ((degree - exact) / exact).abs() < 0.01,
            "U_num = {degree}, U_exact = {exact}"
        );
        // The deformation stays (nearly) one-dimensional: the four nodes of a cross
        // section move together up to the effect of the alternating tetrahedral
        // split (about 0.15 %).
        for &n in &top {
            let spread = (solver.u_dofs[3 * n] + mean_settlement).abs();
            assert!(spread < 5e-3 * mean_settlement, "spread {spread}");
        }
        // Excess pore pressure at the impervious base, a distance H from the drain.
        let base: Vec<usize> = (0..n_nodes)
            .filter(|&n| mesh.nodes[n].x.abs() < 1e-9)
            .collect();
        let p_exact = pore_pressure_ratio(height, height, tv_end);
        for &n in &base {
            let p = solver.p_dofs[n] / load;
            assert!(
                (p - p_exact).abs() < 0.02,
                "p_base/q = {p}, exact = {p_exact}"
            );
        }
    }
}
