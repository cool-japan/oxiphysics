// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense LU factorisation and the Dirichlet-to-Neumann condensation behind
//! [`super::BemFemCoupling::condense`].
//!
//! The collocation BEM for the Laplace equation produces the dense system
//! `H·u = G·q` (potential `u`, flux `q = ∂u/∂n`, one DOF per boundary element).
//! Eliminating the flux gives the *Dirichlet-to-Neumann map*
//!
//! ```text
//! q = G⁻¹ · H · u        ⇒        K = G⁻¹ · H
//! ```
//!
//! which is what an FEM domain sees of the BEM region.  `G` is factorised
//! **once** (LU with partial pivoting, `O(n³)`) and every column of `H` is then
//! back-substituted (`O(n²)` per column), so the whole condensation costs
//! `O(n³)` instead of the `O(n⁴)` of inverting `G` column by column.

use super::{BemFemCoupling, DenseMatrix, assemble_laplace_bem};
use crate::error::{Error, Result};

impl BemFemCoupling {
    /// Condense the BEM influence matrices into the Dirichlet-to-Neumann map
    /// used to couple the BEM region with an FEM stiffness.
    ///
    /// The collocation system of the Laplace BEM is `H·u = G·q` (potential `u`
    /// and flux `q = ∂u/∂n`, one DOF per BEM element).  Eliminating the flux
    /// gives `q = G⁻¹·H·u`, so the condensed matrix is
    ///
    /// `K = G⁻¹ · H`
    ///
    /// (note the order: `G⁻¹` multiplies `H` from the left).  It is indexed by
    /// BEM element, i.e. by collocation point, and maps element potentials to
    /// element fluxes.  Because the rigid-body condition makes `H·1 = 0`,
    /// constant potentials lie in its null space: `K·1 = 0`.
    ///
    /// `G` is LU-factorised once with partial pivoting and all columns of `H`
    /// are back-substituted, `O(n³)` in total (the column-by-column inverse of
    /// `G` cost `O(n⁴)`).
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidInput`] for an empty BEM mesh or a non-finite
    ///   influence coefficient (for example a zero-area element);
    /// - [`Error::SingularSystem`] if `G` is singular to working precision.
    ///
    /// On error `condensed_stiffness` keeps its previous value.
    pub fn condense(&mut self) -> Result<()> {
        let (h, g) = assemble_laplace_bem(&self.bem_mesh);
        self.condensed_stiffness = Some(dirichlet_to_neumann_map(&h, &g)?);
        Ok(())
    }

    /// The condensed map scaled row-wise by the BEM element areas,
    /// `diag(A)·G⁻¹·H`.
    ///
    /// For constant elements `f_i = A_i·q_i` is the flux resultant integrated
    /// over element `i`, so this is the force-like stiffness (flux resultant per
    /// unit potential) to assemble next to an FEM stiffness.  Returns `None`
    /// until [`condense`](Self::condense) has succeeded.
    pub fn flux_resultant_stiffness(&self) -> Option<DenseMatrix> {
        let condensed = self.condensed_stiffness.as_ref()?;
        let mut scaled = condensed.clone();
        for (i, element) in self.bem_mesh.elements.iter().enumerate().take(scaled.rows) {
            for j in 0..scaled.cols {
                scaled.set(i, j, element.area * condensed.get(i, j));
            }
        }
        Some(scaled)
    }
}

/// LU factorisation `P·A = L·U` of a square matrix with partial (row) pivoting.
struct DenseLu {
    n: usize,
    /// Packed factors, row-major: the strict lower triangle holds `L` (its unit
    /// diagonal is implicit), the upper triangle including the diagonal holds `U`.
    lu: Vec<f64>,
    /// `perm[i]` is the row of `A` that ended up in row `i` of `P·A`.
    perm: Vec<usize>,
}

impl DenseLu {
    /// Factorise `a`.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidInput`] if `a` is not square or has a non-finite entry;
    /// - [`Error::SingularSystem`] if a pivot is not larger than
    ///   `n · ε · max|aᵢⱼ|`, i.e. `a` is singular to working precision.  The
    ///   threshold is relative to the matrix scale, so it behaves identically
    ///   for small and large elements.
    fn factor(a: &DenseMatrix) -> Result<Self> {
        if a.rows != a.cols {
            return Err(Error::InvalidInput(format!(
                "cannot factorise a non-square {} x {} matrix",
                a.rows, a.cols
            )));
        }
        let n = a.rows;
        let mut lu = a.data.clone();
        if lu.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput(
                "influence matrix contains a non-finite entry (degenerate boundary element?)"
                    .to_string(),
            ));
        }
        let scale = lu
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        let pivot_tolerance = n as f64 * f64::EPSILON * scale;
        let mut perm: Vec<usize> = (0..n).collect();

        for k in 0..n {
            let mut pivot_row = k;
            let mut pivot_abs = lu[k * n + k].abs();
            for i in (k + 1)..n {
                let candidate = lu[i * n + k].abs();
                if candidate > pivot_abs {
                    pivot_abs = candidate;
                    pivot_row = i;
                }
            }
            if pivot_abs <= pivot_tolerance {
                return Err(Error::SingularSystem);
            }
            if pivot_row != k {
                for j in 0..n {
                    lu.swap(k * n + j, pivot_row * n + j);
                }
                perm.swap(k, pivot_row);
            }
            let pivot = lu[k * n + k];
            for i in (k + 1)..n {
                let multiplier = lu[i * n + k] / pivot;
                lu[i * n + k] = multiplier;
                if multiplier != 0.0 {
                    for j in (k + 1)..n {
                        let update = multiplier * lu[k * n + j];
                        lu[i * n + j] -= update;
                    }
                }
            }
        }
        Ok(Self { n, lu, perm })
    }

    /// Solve `A·X = B` for all columns of `B` at once.
    ///
    /// The substitutions run over whole rows of `X`, so the inner loops are
    /// contiguous in memory.  `B` must have `n` rows.
    fn solve_matrix(&self, b: &DenseMatrix) -> Result<DenseMatrix> {
        let n = self.n;
        if b.rows != n {
            return Err(Error::InvalidInput(format!(
                "right-hand side has {} rows but the factorised matrix has {n}",
                b.rows
            )));
        }
        let m = b.cols;
        let mut x = vec![0.0; n * m];

        // Forward substitution: L·Y = P·B.
        for i in 0..n {
            let source = self.perm[i];
            x[i * m..(i + 1) * m].copy_from_slice(&b.data[source * m..(source + 1) * m]);
            for j in 0..i {
                let multiplier = self.lu[i * n + j];
                if multiplier != 0.0 {
                    let (solved, current) = x.split_at_mut(i * m);
                    let row_j = &solved[j * m..(j + 1) * m];
                    for (xi, xj) in current[..m].iter_mut().zip(row_j) {
                        *xi -= multiplier * xj;
                    }
                }
            }
        }

        // Back substitution: U·X = Y.
        for i in (0..n).rev() {
            for j in (i + 1)..n {
                let entry = self.lu[i * n + j];
                if entry != 0.0 {
                    let (upper, lower) = x.split_at_mut(j * m);
                    let row_i = &mut upper[i * m..(i + 1) * m];
                    for (xi, xj) in row_i.iter_mut().zip(&lower[..m]) {
                        *xi -= entry * xj;
                    }
                }
            }
            let diagonal = self.lu[i * n + i];
            for xi in &mut x[i * m..(i + 1) * m] {
                *xi /= diagonal;
            }
        }

        Ok(DenseMatrix {
            rows: n,
            cols: m,
            data: x,
        })
    }
}

/// Dirichlet-to-Neumann map `K = G⁻¹·H` of the collocation system `H·u = G·q`.
///
/// `G` is factorised once and all `n` columns of `H` are back-substituted:
/// `O(n³)` in total.  Note the order: `q = G⁻¹·(H·u)`, so `G⁻¹` multiplies `H`
/// from the **left** (`H·G⁻¹` would map a different space).
///
/// # Errors
///
/// - [`Error::InvalidInput`] if the matrices are empty, not square, differently
///   sized, or contain a non-finite entry;
/// - [`Error::SingularSystem`] if `G` is singular to working precision.
pub(super) fn dirichlet_to_neumann_map(h: &DenseMatrix, g: &DenseMatrix) -> Result<DenseMatrix> {
    if g.rows == 0 || g.rows != g.cols || h.rows != g.rows || h.cols != g.rows {
        return Err(Error::InvalidInput(format!(
            "H ({} x {}) and G ({} x {}) must be non-empty square matrices of equal size",
            h.rows, h.cols, g.rows, g.cols
        )));
    }
    if h.data.iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput(
            "influence matrix H contains a non-finite entry (degenerate boundary element?)"
                .to_string(),
        ));
    }
    let k = DenseLu::factor(g)?.solve_matrix(h)?;
    if k.data.iter().any(|value| !value.is_finite()) {
        return Err(Error::SingularSystem);
    }
    Ok(k)
}

#[cfg(test)]
mod tests {
    use super::super::{
        BoundaryElement, BoundaryMesh, ElementOrder, element_normal, triangle_area,
    };
    use super::*;

    fn matrix(rows: usize, cols: usize, values: &[f64]) -> DenseMatrix {
        assert_eq!(values.len(), rows * cols);
        DenseMatrix {
            rows,
            cols,
            data: values.to_vec(),
        }
    }

    fn multiply(a: &DenseMatrix, b: &DenseMatrix) -> DenseMatrix {
        let mut c = DenseMatrix::zeros(a.rows, b.cols);
        for i in 0..a.rows {
            for j in 0..b.cols {
                let mut sum = 0.0;
                for m in 0..a.cols {
                    sum += a.get(i, m) * b.get(m, j);
                }
                c.set(i, j, sum);
            }
        }
        c
    }

    fn max_abs_difference(a: &DenseMatrix, b: &DenseMatrix) -> f64 {
        assert_eq!((a.rows, a.cols), (b.rows, b.cols));
        a.data
            .iter()
            .zip(&b.data)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }

    /// Explicit inverse of a 3 x 3 matrix by cofactors.
    fn inverse_3x3(a: &DenseMatrix) -> DenseMatrix {
        let m = |i: usize, j: usize| a.get(i, j);
        let cof = |r0: usize, r1: usize, c0: usize, c1: usize| {
            m(r0, c0) * m(r1, c1) - m(r0, c1) * m(r1, c0)
        };
        let det = m(0, 0) * cof(1, 2, 1, 2) - m(0, 1) * cof(1, 2, 0, 2) + m(0, 2) * cof(1, 2, 0, 1);
        let adjugate = [
            cof(1, 2, 1, 2),
            -cof(0, 2, 1, 2),
            cof(0, 1, 1, 2),
            -cof(1, 2, 0, 2),
            cof(0, 2, 0, 2),
            -cof(0, 1, 0, 2),
            cof(1, 2, 0, 1),
            -cof(0, 2, 0, 1),
            cof(0, 1, 0, 1),
        ];
        matrix(3, 3, &adjugate.map(|v| v / det))
    }

    /// Three non-coplanar triangles forming an open "tent".
    fn tent_mesh() -> BoundaryMesh {
        let nodes = vec![
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [-0.5, 0.9, 0.0],
            [-0.5, -0.9, 0.0],
        ];
        let corners = [[0, 1, 2], [0, 2, 3], [0, 3, 1]];
        let elements = corners
            .iter()
            .map(|c| BoundaryElement {
                node_ids: c.to_vec(),
                normal: element_normal(&nodes[c[0]], &nodes[c[1]], &nodes[c[2]]),
                area: triangle_area(&nodes[c[0]], &nodes[c[1]], &nodes[c[2]]),
                order: ElementOrder::Constant,
            })
            .collect();
        BoundaryMesh { nodes, elements }
    }

    #[test]
    fn dirichlet_to_neumann_map_is_g_inverse_times_h_not_h_times_g_inverse() {
        let g = matrix(3, 3, &[4.0, 1.0, 0.5, 1.0, 5.0, 1.0, 0.5, 2.0, 6.0]);
        let h = matrix(3, 3, &[1.0, -1.0, 0.0, 0.0, 2.0, -2.0, 3.0, 0.0, -3.0]);
        let k = dirichlet_to_neumann_map(&h, &g).expect("G is regular");

        let g_inv = inverse_3x3(&g);
        let expected = multiply(&g_inv, &h);
        assert!(
            max_abs_difference(&k, &expected) < 1e-13,
            "K = {:?}, expected {:?}",
            k.data,
            expected.data
        );
        // G·K must reproduce H ...
        assert!(max_abs_difference(&multiply(&g, &k), &h) < 1e-12);
        // ... and the pre-fix product H·G⁻¹ is a genuinely different matrix.
        let wrong_order = multiply(&h, &g_inv);
        assert!(max_abs_difference(&k, &wrong_order) > 1e-2);
    }

    #[test]
    fn lu_pivots_around_a_zero_leading_entry_and_solves_many_right_hand_sides() {
        // A zero in the (0, 0) position forces a row exchange.
        let a = matrix(3, 3, &[0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0]);
        let b = matrix(3, 2, &[1.0, 0.0, 2.0, 1.0, 3.0, -1.0]);
        let x = DenseLu::factor(&a)
            .and_then(|lu| lu.solve_matrix(&b))
            .expect("A is regular");
        assert!(max_abs_difference(&multiply(&a, &x), &b) < 1e-13);
    }

    #[test]
    fn singular_g_is_reported_instead_of_silently_skipped() {
        // Second row = 2 * first row.
        let g = matrix(2, 2, &[1.0, 2.0, 2.0, 4.0]);
        let h = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let result = dirichlet_to_neumann_map(&h, &g);
        assert!(matches!(result, Err(Error::SingularSystem)), "{result:?}");

        // The threshold is relative: a regular but tiny-valued G is accepted ...
        let tiny = matrix(2, 2, &[1.0e-30, 0.0, 0.0, 2.0e-30]);
        assert!(dirichlet_to_neumann_map(&h, &tiny).is_ok());
        // ... and a huge singular one is still rejected.
        let huge = matrix(2, 2, &[1.0e30, 2.0e30, 2.0e30, 4.0e30]);
        assert!(matches!(
            dirichlet_to_neumann_map(&h, &huge),
            Err(Error::SingularSystem)
        ));
        // The zero matrix is singular too.
        let zero = DenseMatrix::zeros(2, 2);
        assert!(matches!(
            dirichlet_to_neumann_map(&h, &zero),
            Err(Error::SingularSystem)
        ));
    }

    #[test]
    fn malformed_input_is_rejected() {
        let g = matrix(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let wrong_size = matrix(3, 3, &[1.0; 9]);
        assert!(matches!(
            dirichlet_to_neumann_map(&wrong_size, &g),
            Err(Error::InvalidInput(_))
        ));
        let empty = DenseMatrix::zeros(0, 0);
        assert!(matches!(
            dirichlet_to_neumann_map(&empty, &empty),
            Err(Error::InvalidInput(_))
        ));
        let nan_h = matrix(2, 2, &[f64::NAN, 0.0, 0.0, 1.0]);
        assert!(matches!(
            dirichlet_to_neumann_map(&nan_h, &g),
            Err(Error::InvalidInput(_))
        ));
        let nan_g = matrix(2, 2, &[1.0, 0.0, 0.0, f64::INFINITY]);
        assert!(matches!(
            dirichlet_to_neumann_map(&g, &nan_g),
            Err(Error::InvalidInput(_))
        ));
        let rectangular = matrix(2, 3, &[1.0; 6]);
        assert!(matches!(
            DenseLu::factor(&rectangular),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn condense_three_element_mesh_matches_explicit_g_inverse_h() {
        let mesh = tent_mesh();
        let (h, g) = assemble_laplace_bem(&mesh);
        assert!(
            h.data.iter().any(|v| v.abs() > 1e-6),
            "the tent has non-zero flux kernels"
        );

        let mut coupling = BemFemCoupling::new(vec![0, 1, 2], mesh);
        coupling.condense().expect("tent G is regular");
        let k = coupling.condensed_stiffness.as_ref().expect("stored");

        let expected = multiply(&inverse_3x3(&g), &h);
        let scale = expected.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(
            max_abs_difference(k, &expected) < 1e-10 * scale,
            "K = {:?}\nexpected {:?}",
            k.data,
            expected.data
        );
        // G·K = H, the defining relation.
        let h_scale = h.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(max_abs_difference(&multiply(&g, k), &h) < 1e-10 * h_scale);
    }

    #[test]
    fn condensed_stiffness_of_a_closed_surface_annihilates_constant_potentials() {
        // H·1 = 0 (rigid-body condition), hence K·1 = G⁻¹·H·1 = 0.
        let mesh = BoundaryMesh::sphere(1.0, 4, 8);
        let n = mesh.num_elements();
        assert_eq!(n, 48);
        let (h, _) = assemble_laplace_bem(&mesh);
        assert!(
            h.data.iter().any(|v| v.abs() > 1e-3),
            "H must be non-trivial"
        );

        let mut coupling = BemFemCoupling::new((0..mesh.num_nodes()).collect(), mesh);
        coupling.condense().expect("sphere G is regular");
        let k = coupling.condensed_stiffness.as_ref().expect("stored");
        assert_eq!((k.rows, k.cols), (n, n));

        let k_ones = k.matvec(&vec![1.0; n]);
        let k_norm = k.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let residual = k_ones.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(k_norm > 1e-6, "K must not vanish: {k_norm}");
        assert!(
            residual <= 1e-10 * k_norm,
            "‖K·1‖ = {residual:e} vs ‖K‖ = {k_norm:e}"
        );
    }

    #[test]
    fn condense_scales_cubically_not_quartically() {
        // 624 elements: the column-by-column inverse needed ≈ 5·10¹⁰ flops
        // (about a minute); factor-once needs ≈ 8·10⁷ (milliseconds).
        let mesh = BoundaryMesh::sphere(1.0, 14, 24);
        let n = mesh.num_elements();
        assert_eq!(n, 624);
        let mut coupling = BemFemCoupling::new(vec![0], mesh);
        let start = std::time::Instant::now();
        coupling.condense().expect("sphere G is regular");
        let elapsed = start.elapsed();
        assert!(
            elapsed < std::time::Duration::from_secs(30),
            "condense took {elapsed:?}"
        );
        let k = coupling.condensed_stiffness.as_ref().expect("stored");
        let k_norm = k.data.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let residual = k
            .matvec(&vec![1.0; n])
            .iter()
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(
            residual <= 1e-9 * k_norm,
            "‖K·1‖ = {residual:e}, ‖K‖ = {k_norm:e}"
        );
    }

    #[test]
    fn flux_resultant_stiffness_scales_rows_by_the_element_areas() {
        let mesh = tent_mesh();
        let areas: Vec<f64> = mesh.elements.iter().map(|e| e.area).collect();
        let mut coupling = BemFemCoupling::new(vec![0, 1, 2], mesh);
        assert!(coupling.flux_resultant_stiffness().is_none());
        coupling.condense().expect("regular");
        let k = coupling.condensed_stiffness.clone().expect("stored");
        let k_f = coupling.flux_resultant_stiffness().expect("condensed");
        for (i, &area) in areas.iter().enumerate() {
            for j in 0..3 {
                let expected = area * k.get(i, j);
                assert!((k_f.get(i, j) - expected).abs() <= 1e-12 * (1.0 + expected.abs()));
            }
        }
    }

    #[test]
    fn condense_reports_empty_and_degenerate_meshes_and_keeps_the_previous_result() {
        // Empty mesh.
        let mut coupling = BemFemCoupling::new(Vec::new(), BoundaryMesh::new());
        assert!(matches!(coupling.condense(), Err(Error::InvalidInput(_))));
        assert!(coupling.condensed_stiffness.is_none());

        // A zero-area element makes G's self term 0/0.
        let mut mesh = tent_mesh();
        mesh.elements[1].area = 0.0;
        let mut coupling = BemFemCoupling::new(vec![0, 1, 2], mesh);
        let previous = DenseMatrix::zeros(3, 3);
        coupling.condensed_stiffness = Some(previous.clone());
        assert!(matches!(coupling.condense(), Err(Error::InvalidInput(_))));
        let kept = coupling.condensed_stiffness.as_ref().expect("kept");
        assert_eq!(kept.data, previous.data);
    }
}
