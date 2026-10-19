// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sparse neighbourhood filters for element-based topology optimisation.
//!
//! The classical density and sensitivity filters average a field over all
//! elements whose centroid lies within a radius `r_min` of the element under
//! consideration, with a radial weight such as the hat function
//! `w = max(0, r_min − d)`.  Evaluating that with an all-pairs loop costs `O(n²)`
//! distance evaluations per application and a dense filter matrix needs `n²`
//! entries (a 200 × 100 element mesh already needs 3.2 GB).
//!
//! [`SparseFilter`] removes both problems.  It is built **once per mesh**: the
//! centroids are binned into a uniform cell list whose cell size equals the
//! filter radius, so only the elements in the 3 × 3 surrounding cells can be
//! neighbours; squared distances are compared against `r_min²` before any square
//! root is taken.  The weights are stored as a [`CsrMatrix`] `H` together with the
//! row sums `Σⱼ Hᵢⱼ`, and every filter is then a single sparse matrix–vector
//! product (row-parallel with rayon for large meshes).  Time and memory are
//! `O(n · k)` with `k` the number of neighbours inside the radius.
//!
//! All results are bit-compatible with the former dense loops: the stored
//! weights are exactly the non-zero weights of the dense kernel and every row is
//! summed in ascending column order.

use crate::sparse::CsrMatrix;
use rayon::prelude::*;
use std::collections::HashMap;

/// Meshes with at least this many elements are built / applied row-parallel.
const PARALLEL_ROW_THRESHOLD: usize = 2048;

/// Rows whose weight sum does not exceed this value are treated as empty (the
/// filtered value falls back to the unfiltered one).
const MIN_WEIGHT_SUM: f64 = 1.0e-30;

/// Relative slack of the squared-distance prune.  Any pair whose kernel weight
/// can be non-zero (`d < r`) is always kept; the exact weight decides the rest.
const PRUNE_SLACK: f64 = 1.0 + 8.0 * f64::EPSILON;

/// Radial weight function of a [`SparseFilter`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterKernel {
    /// Linear cone (hat) kernel `w = max(0, r − d)` (Bruns & Tortorelli).
    Hat,
    /// Gaussian kernel `w = exp(−d² / (2σ²))` with `σ = r / 3`, truncated at
    /// `d = r` (elements with `d ≥ r` get zero weight).
    Gaussian,
}

/// Sparse, row-normalisable neighbourhood filter of a 2-D element mesh.
///
/// `h` holds the raw weights `Hᵢⱼ` (the diagonal included) and `row_sum[i]` is
/// `Σⱼ Hᵢⱼ`.  The filtered value of a field `v` is `(H v)ᵢ / row_sum[i]`.
#[derive(Debug, Clone)]
pub struct SparseFilter {
    /// Raw weights `H` (`n × n`, CSR, columns ascending within a row).
    pub h: CsrMatrix,
    /// Row sums `Σⱼ Hᵢⱼ`.
    pub row_sum: Vec<f64>,
}

/// Uniform cell list over 2-D points (cell edge = filter radius).
struct CellList {
    cells: HashMap<(i64, i64), Vec<usize>>,
}

impl CellList {
    fn cell_of(x: f64, y: f64, cell: f64) -> (i64, i64) {
        // `as` saturates for huge ratios; points that saturate share a cell and are
        // still separated by the exact distance test.
        ((x / cell).floor() as i64, (y / cell).floor() as i64)
    }

    fn new(points: &[(f64, f64)], cell: f64) -> Self {
        let mut cells: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, &(x, y)) in points.iter().enumerate() {
            if x.is_finite() && y.is_finite() {
                cells.entry(Self::cell_of(x, y, cell)).or_default().push(i);
            }
        }
        Self { cells }
    }

    /// Non-zero kernel weights `(j, w_ij)` of row `i`, columns ascending.
    fn row(
        &self,
        points: &[(f64, f64)],
        i: usize,
        radius: f64,
        kernel: FilterKernel,
    ) -> Vec<(usize, f64)> {
        let (xi, yi) = points[i];
        if !(xi.is_finite() && yi.is_finite()) {
            return Vec::new();
        }
        let (cx, cy) = Self::cell_of(xi, yi, radius);
        let mut candidates: Vec<usize> = Vec::new();
        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                let (Some(kx), Some(ky)) = (cx.checked_add(dx), cy.checked_add(dy)) else {
                    continue;
                };
                if let Some(list) = self.cells.get(&(kx, ky)) {
                    candidates.extend_from_slice(list);
                }
            }
        }
        candidates.sort_unstable();

        let r2 = radius * radius;
        let sigma = radius / 3.0;
        let mut row = Vec::with_capacity(candidates.len());
        for j in candidates {
            let (xj, yj) = points[j];
            let d2 = (xi - xj).powi(2) + (yi - yj).powi(2);
            let weight = match kernel {
                FilterKernel::Hat => {
                    if d2 > r2 * PRUNE_SLACK {
                        continue;
                    }
                    (radius - d2.sqrt()).max(0.0)
                }
                FilterKernel::Gaussian => {
                    if d2 < r2 {
                        (-d2 / (2.0 * sigma * sigma)).exp()
                    } else {
                        0.0
                    }
                }
            };
            if weight > 0.0 {
                row.push((j, weight));
            }
        }
        row
    }
}

impl SparseFilter {
    /// Build the filter of `centroids` for the given `radius` and `kernel`.
    ///
    /// A radius that is not positive and finite yields an empty filter (no
    /// element has any weight, every filtered value falls back to the unfiltered
    /// one).  Elements with a non-finite centroid have no neighbours.
    pub fn new(centroids: &[(f64, f64)], radius: f64, kernel: FilterKernel) -> Self {
        let n = centroids.len();
        let rows: Vec<Vec<(usize, f64)>> = if radius.is_finite() && radius > 0.0 {
            let cells = CellList::new(centroids, radius);
            let build_row = |i: usize| cells.row(centroids, i, radius, kernel);
            if n >= PARALLEL_ROW_THRESHOLD {
                (0..n).into_par_iter().map(build_row).collect()
            } else {
                (0..n).map(build_row).collect()
            }
        } else {
            vec![Vec::new(); n]
        };

        let nnz: usize = rows.iter().map(Vec::len).sum();
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::with_capacity(nnz);
        let mut values = Vec::with_capacity(nnz);
        let mut row_sum = Vec::with_capacity(n);
        row_ptr.push(0);
        for row in &rows {
            let mut sum = 0.0_f64;
            for &(j, w) in row {
                col_indices.push(j);
                values.push(w);
                sum += w;
            }
            row_sum.push(sum);
            row_ptr.push(col_indices.len());
        }
        Self {
            h: CsrMatrix {
                row_ptr,
                col_indices,
                values,
                nrows: n,
                ncols: n,
            },
            row_sum,
        }
    }

    /// Hat-kernel filter (`w = max(0, r_min − d)`), the standard choice.
    pub fn hat(centroids: &[(f64, f64)], r_min: f64) -> Self {
        Self::new(centroids, r_min, FilterKernel::Hat)
    }

    /// Number of elements the filter was built for.
    pub fn len(&self) -> usize {
        self.row_sum.len()
    }

    /// Whether the filter covers no element.
    pub fn is_empty(&self) -> bool {
        self.row_sum.is_empty()
    }

    /// Number of stored weights (`nnz` of `H`).
    pub fn nnz(&self) -> usize {
        self.h.nnz()
    }

    /// Evaluate `row_fn` for every row, in parallel for large meshes.
    fn map_rows<F>(&self, row_fn: F) -> Vec<f64>
    where
        F: Fn(usize) -> f64 + Sync + Send,
    {
        let n = self.len();
        if n >= PARALLEL_ROW_THRESHOLD {
            (0..n).into_par_iter().map(row_fn).collect()
        } else {
            (0..n).map(row_fn).collect()
        }
    }

    /// Sparse matrix–vector product `y = H x` (row-parallel for large meshes).
    ///
    /// # Panics
    /// Panics if `x.len()` differs from the number of elements.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(
            x.len(),
            self.len(),
            "field length must equal the element count"
        );
        let h = &self.h;
        self.map_rows(|i| {
            let mut sum = 0.0_f64;
            for idx in h.row_ptr[i]..h.row_ptr[i + 1] {
                sum += h.values[idx] * x[h.col_indices[idx]];
            }
            sum
        })
    }

    /// Weighted average `(H v)ᵢ / Σⱼ Hᵢⱼ` of `values`.
    ///
    /// Rows whose weight sum does not exceed `1e-30` keep `values[i]`.
    ///
    /// # Panics
    /// Panics if `values.len()` differs from the number of elements.
    pub fn average(&self, values: &[f64]) -> Vec<f64> {
        let numerator = self.apply(values);
        numerator
            .iter()
            .zip(self.row_sum.iter())
            .zip(values.iter())
            .map(|((&num, &den), &own)| if den > MIN_WEIGHT_SUM { num / den } else { own })
            .collect()
    }

    /// Density-weighted average `Σⱼ Hᵢⱼ ρⱼ vⱼ / Σⱼ Hᵢⱼ ρⱼ` of `values`.
    ///
    /// Rows whose weighted density sum does not exceed `1e-30` in magnitude keep
    /// `values[i]`.
    ///
    /// # Panics
    /// Panics if the slice lengths differ from the number of elements.
    pub fn density_weighted_average(&self, values: &[f64], densities: &[f64]) -> Vec<f64> {
        let n = self.len();
        assert_eq!(
            values.len(),
            n,
            "values length must equal the element count"
        );
        assert_eq!(
            densities.len(),
            n,
            "densities length must equal the element count"
        );
        let h = &self.h;
        let numerators_and_denominators: Vec<(f64, f64)> = {
            let row_pair = |i: usize| {
                let mut num = 0.0_f64;
                let mut den = 0.0_f64;
                for idx in h.row_ptr[i]..h.row_ptr[i + 1] {
                    let j = h.col_indices[idx];
                    let w = h.values[idx];
                    num += w * densities[j] * values[j];
                    den += w * densities[j];
                }
                (num, den)
            };
            if n >= PARALLEL_ROW_THRESHOLD {
                (0..n).into_par_iter().map(row_pair).collect()
            } else {
                (0..n).map(row_pair).collect()
            }
        };
        numerators_and_denominators
            .iter()
            .zip(values.iter())
            .map(|(&(num, den), &own)| {
                if den.abs() > MIN_WEIGHT_SUM {
                    num / den
                } else {
                    own
                }
            })
            .collect()
    }

    /// Row-normalised copy of `H` (every row sums to one).
    ///
    /// Rows whose weight sum does not exceed `1e-30` are left unnormalised.
    pub fn normalised(&self) -> CsrMatrix {
        let mut matrix = self.h.clone();
        for (i, &sum) in self.row_sum.iter().enumerate() {
            if sum > MIN_WEIGHT_SUM {
                for idx in matrix.row_ptr[i]..matrix.row_ptr[i + 1] {
                    matrix.values[idx] /= sum;
                }
            }
        }
        matrix
    }
}

#[cfg(test)]
mod tests;
