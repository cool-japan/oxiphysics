// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Incremental double-description enumeration of the extreme rays of a pointed
//! polyhedral cone.
//!
//! The cone is given by inequalities,
//!
//! ```text
//! K = { x in R^d : a_i . x <= 0 for every row a_i }
//! ```
//!
//! and [`ConeEnumerator`] maintains its extreme rays while half-spaces are
//! intersected one at a time (Motzkin et al.; Fukuda & Prodon, "Double
//! description method revisited", 1996).  It starts from the simplicial cone
//! of `d` linearly independent rows.  Inserting a row discards the rays
//! strictly outside the new half-space and creates one new ray on the bounding
//! hyperplane for every *adjacent* pair of rays (one outside, one inside).
//!
//! Adjacency is decided combinatorially by the zero sets
//! `Z(r) = { i : a_i . r = 0 }`: two rays are adjacent iff no third ray
//! satisfies at least the constraints their zero sets share.  The test stays
//! valid on highly degenerate input (many rows through one ray, repeated
//! rows), which is what contact wrench sets look like: all the cone edges of
//! one contact lie in a common plane.
//!
//! Because the rows can be inserted in any order and the enumerator can be
//! inspected between insertions, callers may add only the rows they need
//! (cutting-plane methods).  Rows are normalised to unit length and every ray
//! is kept at unit length, so one absolute tolerance serves all comparisons.
//! A [`Budget`] bounds the work; an insertion that would exceed it is refused
//! and leaves the enumerator unchanged.

/// Dot products below this magnitude count as zero (rows and rays have unit
/// length).
pub(super) const ZERO_TOL: f64 = 1e-9;

/// A row whose residual after orthogonalisation against the chosen rows is
/// below this is linearly dependent on them.
const RANK_TOL: f64 = 1e-9;

/// Smallest pivot accepted when inverting the initial basis.
const PIVOT_TOL: f64 = 1e-13;

/// Work limits of one enumeration.
#[derive(Debug, Clone, Copy)]
pub(super) struct Budget {
    /// Largest number of rays held at any time.
    pub(super) max_rays: usize,
    /// Largest amount of work in abstract units of roughly three nanoseconds
    /// (one comparison of a pair of zero sets costs a few units).
    pub(super) max_work: u64,
}

impl Budget {
    /// Limit used by the public grasp metrics: generous enough for every
    /// realistic grasp, small enough to bound the run time to a few seconds on
    /// pathological input.  Measured with an optimised build on random sphere
    /// grasps with linearised friction cones: up to about 150 wrenches take
    /// milliseconds to a fraction of a second, a few hundred wrenches can take
    /// one to two seconds, and the budget caps the run time beyond that.
    pub(super) const STANDARD: Self = Self {
        max_rays: 200_000,
        max_work: 800_000_000,
    };
}

/// An insertion was refused because the [`Budget`] would be exceeded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct BudgetExhausted;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Scale `v` to unit length; `false` if it is zero or not finite.
fn normalise(v: &mut [f64]) -> bool {
    let norm = dot(v, v).sqrt();
    if norm.is_finite() && norm > 0.0 {
        for component in v.iter_mut() {
            *component /= norm;
        }
        true
    } else {
        false
    }
}

fn set_bit(bits: &mut [u64], index: usize) {
    if let Some(word) = bits.get_mut(index / 64) {
        *word |= 1_u64 << (index % 64);
    }
}

/// `true` iff every bit of `inner` is also set in `outer`.
fn is_subset(inner: &[u64], outer: &[u64]) -> bool {
    inner.iter().zip(outer.iter()).all(|(i, o)| i & !o == 0)
}

/// Indices of the set bits of a bit set, in increasing order.
fn set_bits(bits: &[u64]) -> impl Iterator<Item = usize> + '_ {
    bits.iter().enumerate().flat_map(|(word_index, &word)| {
        let mut rest = word;
        std::iter::from_fn(move || {
            if rest == 0 {
                None
            } else {
                let bit = rest.trailing_zeros() as usize;
                rest &= rest - 1;
                Some(word_index * 64 + bit)
            }
        })
    })
}

/// Choose `dim` well-conditioned linearly independent rows (Gram-Schmidt with
/// the largest residual as pivot), starting with the `seed` rows.  `None` if
/// the rows have rank below `dim` (the cone is not pointed) or the seed rows
/// are dependent or out of range.
fn select_basis(rows: &[Vec<f64>], dim: usize, seed: &[usize]) -> Option<Vec<usize>> {
    let mut residual: Vec<Vec<f64>> = rows.to_vec();
    let mut taken = vec![false; rows.len()];
    let mut chosen = Vec::with_capacity(dim);
    for step in 0..dim {
        let (index, norm) = match seed.get(step) {
            Some(&index) => {
                let row = residual.get(index)?;
                if taken[index] {
                    return None;
                }
                (index, dot(row, row).sqrt())
            }
            None => {
                let mut best: Option<(usize, f64)> = None;
                for (index, row) in residual.iter().enumerate() {
                    if taken[index] {
                        continue;
                    }
                    let norm = dot(row, row).sqrt();
                    if best.is_none_or(|(_, b)| norm > b) {
                        best = Some((index, norm));
                    }
                }
                best?
            }
        };
        if norm <= RANK_TOL {
            return None;
        }
        let direction: Vec<f64> = residual[index].iter().map(|v| v / norm).collect();
        for row in residual.iter_mut() {
            let projection = dot(row, &direction);
            for (value, d) in row.iter_mut().zip(direction.iter()) {
                *value -= projection * d;
            }
        }
        taken[index] = true;
        chosen.push(index);
    }
    Some(chosen)
}

/// Inverse of a square matrix by Gauss-Jordan elimination with partial
/// pivoting; `None` if it is numerically singular.
fn invert(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut work: Vec<Vec<f64>> = matrix
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut extended = row.clone();
            extended.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            extended
        })
        .collect();
    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_abs = 0.0_f64;
        for (r, row) in work.iter().enumerate().skip(col) {
            if row[col].abs() > pivot_abs {
                pivot_row = r;
                pivot_abs = row[col].abs();
            }
        }
        if pivot_abs <= PIVOT_TOL {
            return None;
        }
        work.swap(col, pivot_row);
        let pivot = work[col][col];
        for value in work[col].iter_mut() {
            *value /= pivot;
        }
        let pivot_values = work[col].clone();
        for (r, row) in work.iter_mut().enumerate() {
            if r == col {
                continue;
            }
            let factor = row[col];
            if factor != 0.0 {
                for (value, p) in row.iter_mut().zip(pivot_values.iter()) {
                    *value -= factor * p;
                }
            }
        }
    }
    Some(work.into_iter().map(|row| row[n..].to_vec()).collect())
}

/// Rays together with their zero sets (flat storage).
struct RaySet {
    dim: usize,
    words: usize,
    xs: Vec<f64>,
    zeros: Vec<u64>,
}

impl RaySet {
    fn new(dim: usize, words: usize) -> Self {
        Self {
            dim,
            words,
            xs: Vec::new(),
            zeros: Vec::new(),
        }
    }

    fn len(&self) -> usize {
        self.xs.len() / self.dim
    }

    fn push(&mut self, ray: &[f64], zero: &[u64]) {
        self.xs.extend_from_slice(ray);
        self.zeros.extend_from_slice(zero);
    }

    fn ray(&self, index: usize) -> &[f64] {
        &self.xs[index * self.dim..(index + 1) * self.dim]
    }

    fn zero(&self, index: usize) -> &[u64] {
        &self.zeros[index * self.words..(index + 1) * self.words]
    }

    /// For every row, the rays that satisfy it with equality.
    fn rays_by_row(&self, rows: usize) -> Vec<Vec<u32>> {
        let mut lists: Vec<Vec<u32>> = vec![Vec::new(); rows];
        for (k, zero) in self.zeros.chunks_exact(self.words).enumerate() {
            for row in set_bits(zero) {
                if let Some(list) = lists.get_mut(row) {
                    list.push(k as u32);
                }
            }
        }
        lists
    }
}

/// The shortest of the per-row ray lists over the rows in `common`, if any.
fn shortest_list<'a>(lists: &'a [Vec<u32>], common: &[u64]) -> Option<&'a [u32]> {
    set_bits(common)
        .filter_map(|row| lists.get(row))
        .min_by_key(|list| list.len())
        .map(|list| list.as_slice())
}

/// Extreme rays of the pointed cone `{ x : a . x <= 0 for the processed rows }`,
/// maintained incrementally.
pub(super) struct ConeEnumerator {
    dim: usize,
    words: usize,
    rows: Vec<Vec<f64>>,
    processed: Vec<bool>,
    set: RaySet,
    work: u64,
    budget: Budget,
    common: Vec<u64>,
    values: Vec<f64>,
}

impl ConeEnumerator {
    /// Start from `dim` linearly independent rows, which always include the
    /// rows listed in `seed` (in that order).
    ///
    /// `rows` are inequality normals in `R^dim` (any non-zero length; they are
    /// normalised, zero rows are never inserted).  Returns `None` for malformed
    /// input (wrong row length, non-finite entry, `dim = 0`), when the rows
    /// have rank below `dim` (the cone contains a line and has no extreme
    /// rays), or when the seed is invalid.
    pub(super) fn new(
        rows: &[Vec<f64>],
        dim: usize,
        seed: &[usize],
        budget: Budget,
    ) -> Option<Self> {
        if dim == 0 || seed.len() > dim {
            return None;
        }
        let mut unit: Vec<Vec<f64>> = Vec::with_capacity(rows.len());
        let mut processed = Vec::with_capacity(rows.len());
        for row in rows {
            if row.len() != dim || row.iter().any(|v| !v.is_finite()) {
                return None;
            }
            let mut normalised = row.clone();
            // A zero row (`0 <= 0`) is vacuous: keep its index, never use it.
            processed.push(!normalise(&mut normalised));
            unit.push(normalised);
        }
        let basis = select_basis(&unit, dim, seed)?;
        let basis_matrix: Vec<Vec<f64>> = basis.iter().map(|&i| unit[i].clone()).collect();
        let inverse = invert(&basis_matrix)?;
        let words = unit.len().div_ceil(64);

        // The cone of the basis rows alone is simplicial: B x = -e_j gives its
        // j-th extreme ray, which satisfies every basis row except the j-th
        // with equality.
        let mut set = RaySet::new(dim, words);
        for j in 0..dim {
            let mut ray: Vec<f64> = inverse.iter().map(|inverse_row| -inverse_row[j]).collect();
            if !normalise(&mut ray) {
                return None;
            }
            let mut zero = vec![0_u64; words];
            for (k, &row_index) in basis.iter().enumerate() {
                if k != j {
                    set_bit(&mut zero, row_index);
                }
            }
            set.push(&ray, &zero);
        }
        for &row_index in &basis {
            processed[row_index] = true;
        }
        Some(Self {
            dim,
            words,
            rows: unit,
            processed,
            set,
            work: 0,
            budget,
            common: vec![0_u64; words],
            values: Vec::new(),
        })
    }

    /// Number of rows (including vacuous ones).
    pub(super) fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// The unit-length row with the given index.
    pub(super) fn row(&self, index: usize) -> &[f64] {
        &self.rows[index]
    }

    /// Whether the row has been intersected already (or is vacuous).
    pub(super) fn is_processed(&self, index: usize) -> bool {
        self.processed[index]
    }

    /// The unit-length extreme rays of the current cone.
    pub(super) fn rays(&self) -> impl Iterator<Item = &[f64]> {
        self.set.xs.chunks_exact(self.dim)
    }

    /// Intersect the cone with the half-space of row `index`.
    ///
    /// Inserting a processed row is a no-op.  If the work or ray budget would
    /// be exceeded the call is refused and the enumerator stays unchanged.
    pub(super) fn insert(&mut self, index: usize) -> Result<(), BudgetExhausted> {
        if self.processed.get(index).is_none_or(|done| *done) {
            return Ok(());
        }
        let dim = self.dim;
        let words = self.words;
        let row = &self.rows[index];
        self.values.clear();
        self.values
            .extend(self.set.xs.chunks_exact(dim).map(|ray| dot(row, ray)));
        let outside: Vec<usize> = (0..self.values.len())
            .filter(|&k| self.values[k] > ZERO_TOL)
            .collect();
        let inside: Vec<usize> = (0..self.values.len())
            .filter(|&k| self.values[k] < -ZERO_TOL)
            .collect();

        if outside.is_empty() {
            // Redundant half-space: only the zero sets change.
            for (k, zero) in self.set.zeros.chunks_exact_mut(words).enumerate() {
                if self.values[k].abs() <= ZERO_TOL {
                    set_bit(zero, index);
                }
            }
            self.processed[index] = true;
            return Ok(());
        }

        let total = self.set.len();
        let kept = total - outside.len();
        // Every insertion touches each ray a few times (classification, index,
        // rebuild).
        let mut work = self.work + (total * (dim + words)) as u64;
        let max_work = self.budget.max_work;
        let candidate_cost = 2 + words as u64;
        let scan_cost = 1 + words as u64 / 2;
        let rows = self.rows.len();
        // Adjacent rays share at least `needed` active constraints.
        let needed = dim.saturating_sub(2);
        let lists = self.set.rays_by_row(rows);
        let mut is_inside = vec![false; total];
        for &q in &inside {
            is_inside[q] = true;
        }
        let mut seen = vec![usize::MAX; total];
        let mut candidates: Vec<usize> = Vec::new();
        let mut fresh = RaySet::new(dim, words);
        for &p in &outside {
            // Every partner of `p` has at least `needed` of the active rows of
            // `p` in common, so it contains one of any `active - needed + 1` of
            // them: take the rows with the shortest lists.
            candidates.clear();
            if needed == 0 {
                candidates.extend(inside.iter().copied());
            } else {
                let mut active: Vec<usize> = set_bits(self.set.zero(p)).collect();
                if active.len() < needed {
                    continue;
                }
                active.sort_by_key(|&row| lists.get(row).map_or(0, |list| list.len()));
                let take = active.len() - needed + 1;
                for &row in active.iter().take(take) {
                    let Some(list) = lists.get(row) else {
                        continue;
                    };
                    work += list.len() as u64;
                    for &q in list {
                        let q = q as usize;
                        if is_inside[q] && seen[q] != p {
                            seen[q] = p;
                            candidates.push(q);
                        }
                    }
                }
                if work > max_work {
                    return Err(BudgetExhausted);
                }
            }
            for &q in &candidates {
                work += candidate_cost;
                let mut shared = 0_usize;
                for ((c, a), b) in self
                    .common
                    .iter_mut()
                    .zip(self.set.zero(p).iter())
                    .zip(self.set.zero(q).iter())
                {
                    *c = a & b;
                    shared += c.count_ones() as usize;
                }
                if shared < needed {
                    continue;
                }
                // Adjacent iff no third ray satisfies all constraints that `p`
                // and `q` have in common.
                let blocked = match shortest_list(&lists, &self.common) {
                    Some(list) => {
                        work += list.len() as u64 * scan_cost;
                        list.iter().any(|&k| {
                            let k = k as usize;
                            k != p && k != q && is_subset(&self.common, self.set.zero(k))
                        })
                    }
                    None => {
                        work += total as u64 * scan_cost;
                        (0..total).any(|k| k != p && k != q)
                    }
                };
                if work > max_work {
                    return Err(BudgetExhausted);
                }
                if blocked {
                    continue;
                }
                // values[p] > 0 > values[q]: a positive combination of the two
                // rays lying on the new hyperplane.
                let mut ray: Vec<f64> = self
                    .set
                    .ray(p)
                    .iter()
                    .zip(self.set.ray(q).iter())
                    .map(|(xp, xq)| self.values[p] * xq - self.values[q] * xp)
                    .collect();
                if !normalise(&mut ray) {
                    continue;
                }
                let mut zero = self.common.clone();
                set_bit(&mut zero, index);
                fresh.push(&ray, &zero);
                if kept + fresh.len() > self.budget.max_rays {
                    return Err(BudgetExhausted);
                }
            }
            if work > max_work {
                return Err(BudgetExhausted);
            }
        }

        let mut next = RaySet::new(dim, words);
        for (k, (ray, zero)) in self
            .set
            .xs
            .chunks_exact(dim)
            .zip(self.set.zeros.chunks_exact(words))
            .enumerate()
        {
            if self.values[k] > ZERO_TOL {
                continue;
            }
            if self.values[k] < -ZERO_TOL {
                next.push(ray, zero);
            } else {
                let mut on_plane = zero.to_vec();
                set_bit(&mut on_plane, index);
                next.push(ray, &on_plane);
            }
        }
        next.xs.extend_from_slice(&fresh.xs);
        next.zeros.extend_from_slice(&fresh.zeros);
        self.set = next;
        self.work = work;
        self.processed[index] = true;
        Ok(())
    }
}
