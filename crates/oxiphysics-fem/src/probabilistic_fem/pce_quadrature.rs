// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multi-dimensional Gauss-Hermite quadrature and pseudo-spectral projection
//! onto the multivariate probabilists' Hermite polynomial chaos.
//!
//! For `n_dim` independent standard normal variables `ξ = (ξ₁, …, ξ_d)` the
//! chaos basis is `Ψ_α(ξ) = Π_j He_{α_j}(ξ_j)` with `E[Ψ_α²] = Π_j α_j!`.  The
//! non-intrusive (pseudo-spectral) coefficient of a response `f(ξ)` is
//!
//! ```text
//! c_α = E[f(ξ) Ψ_α(ξ)] / α!  ≈  (1 / α!) Σ_q w_q f(ξ_q) Ψ_α(ξ_q),
//! ```
//!
//! evaluated on a quadrature rule `{(ξ_q, w_q)}` for the standard normal
//! measure.  [`NormalQuadrature`] provides the two classical constructions:
//!
//! * the **full tensor grid** of the `n`-point Gauss-Hermite rule (`nᵈ` nodes,
//!   exact for polynomials of degree ≤ `2n − 1` in every variable), and
//! * the **Smolyak sparse grid** of level `k` built with the combination
//!   technique from the same one-dimensional rules (exact for polynomials of
//!   total degree ≤ `2k + 1`, with a node count that grows only polynomially
//!   in `d` for fixed `k`).
//!
//! [`pce_project`] evaluates the response once per node and projects it onto
//! any list of multi-indices.

use std::collections::HashMap;
use std::f64::consts::{PI, SQRT_2};

use crate::error::{Error, Result};
use crate::stochastic_fem::gauss_hermite_quadrature;

/// Largest number of nodes a [`NormalQuadrature`] may contain (counted before
/// coincident sparse-grid nodes are merged); every node is one evaluation of
/// the response.
pub const MAX_QUADRATURE_NODES: usize = 1_000_000;

/// Largest number of stored node coordinates (`nodes × dimensions`) of a
/// [`NormalQuadrature`], which bounds its memory to about 0.5 GB.
pub const MAX_QUADRATURE_COORDINATES: usize = 64_000_000;

/// Largest number of points of the one-dimensional Gauss-Hermite rules used
/// here.  The outermost weights of an `n`-point rule decay like `e^(-2n)`, so
/// larger rules add nothing for the smooth responses of a polynomial chaos
/// expansion while their cost grows quickly.
pub const MAX_GAUSS_HERMITE_POINTS: usize = 100;

/// Size of the tensor grid above which the default rule of
/// [`SpectralStochasticFem`](super::SpectralStochasticFem) switches from the
/// full tensor grid to a Smolyak sparse grid (see
/// [`NormalQuadrature::for_pce`]).
pub const DEFAULT_MAX_TENSOR_NODES: usize = 100_000;

/// Quadrature rule for the standard normal measure on `ℝᵈ`.
///
/// `Σ_q w_q g(ξ_q)` approximates `E[g(ξ)]` for `ξ ~ N(0, I_d)`.  The weights
/// sum to one; sparse-grid rules contain negative weights.
#[derive(Debug, Clone)]
pub struct NormalQuadrature {
    n_dim: usize,
    /// Row-major node coordinates, `len() × n_dim`.
    nodes: Vec<f64>,
    weights: Vec<f64>,
}

impl NormalQuadrature {
    /// Number of random dimensions `d`.
    pub fn n_dim(&self) -> usize {
        self.n_dim
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.weights.len()
    }

    /// Whether the rule has no node.
    pub fn is_empty(&self) -> bool {
        self.weights.is_empty()
    }

    /// Coordinates of node `index` (a slice of length [`Self::n_dim`]).
    pub fn node(&self, index: usize) -> Option<&[f64]> {
        if index >= self.len() {
            return None;
        }
        let start = index.checked_mul(self.n_dim)?;
        let end = start.checked_add(self.n_dim)?;
        self.nodes.get(start..end)
    }

    /// Weight of node `index`.
    pub fn weight(&self, index: usize) -> Option<f64> {
        self.weights.get(index).copied()
    }

    /// Iterate over `(node coordinates, weight)` pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&[f64], f64)> + '_ {
        let d = self.n_dim;
        self.weights
            .iter()
            .enumerate()
            .map(move |(q, &w)| (&self.nodes[q * d..(q + 1) * d], w))
    }

    /// Sum of all weights (one up to rounding).
    pub fn weight_sum(&self) -> f64 {
        self.weights.iter().sum()
    }

    /// Approximate `E[f(ξ)]` for `ξ ~ N(0, I_d)`.
    pub fn integrate<F: Fn(&[f64]) -> f64>(&self, f: F) -> f64 {
        self.iter().map(|(xi, w)| w * f(xi)).sum()
    }

    /// Rule with a single node at the origin, the only rule in zero dimensions.
    fn origin(n_dim: usize) -> Self {
        Self {
            n_dim,
            nodes: vec![0.0; n_dim],
            weights: vec![1.0],
        }
    }

    /// Full tensor-product grid of the `points_per_dim`-point Gauss-Hermite
    /// rule: `points_per_dim^n_dim` nodes, exact for polynomials of degree
    /// ≤ `2·points_per_dim − 1` in every variable.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `points_per_dim` is zero, the grid would
    /// contain more than [`MAX_QUADRATURE_NODES`] nodes, or the one-dimensional
    /// rule cannot be computed.
    pub fn tensor(n_dim: usize, points_per_dim: usize) -> Result<Self> {
        if points_per_dim == 0 {
            return Err(Error::InvalidInput(
                "a quadrature rule needs at least one point per dimension".to_string(),
            ));
        }
        let total = tensor_node_count(n_dim, points_per_dim)
            .filter(|&count| {
                count <= MAX_QUADRATURE_NODES
                    && count
                        .checked_mul(n_dim)
                        .is_some_and(|coordinates| coordinates <= MAX_QUADRATURE_COORDINATES)
            })
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "the tensor grid of {points_per_dim} points in {n_dim} dimensions exceeds \
                     the limits of {MAX_QUADRATURE_NODES} nodes and \
                     {MAX_QUADRATURE_COORDINATES} coordinates; use a Smolyak sparse grid"
                ))
            })?;
        let rule = probabilists_gauss_hermite(points_per_dim)?;
        let mut nodes = Vec::with_capacity(total * n_dim);
        let mut weights = Vec::with_capacity(total);
        let mut digits = vec![0_usize; n_dim];
        for _ in 0..total {
            let mut weight = 1.0;
            for &digit in &digits {
                let (x, w) = rule[digit];
                nodes.push(x);
                weight *= w;
            }
            weights.push(weight);
            // Odometer increment of the multi-index.
            for digit in digits.iter_mut() {
                *digit += 1;
                if *digit < points_per_dim {
                    break;
                }
                *digit = 0;
            }
        }
        Ok(Self {
            n_dim,
            nodes,
            weights,
        })
    }

    /// Smolyak sparse grid of level `level` for the standard normal measure,
    /// assembled by the combination technique
    ///
    /// ```text
    /// A(k, d) = Σ_{k-d+1 ≤ |m| ≤ k} (-1)^(k-|m|) C(d-1, k-|m|) · ⊗_j U_{m_j + 1},
    /// ```
    ///
    /// where `m ∈ ℕᵈ`, `|m| = Σ m_j` and `U_n` is the `n`-point Gauss-Hermite
    /// rule.  The grid is exact for every polynomial of total degree
    /// ≤ `2k + 1`; level 0 is the single point at the origin and in one
    /// dimension the level-`k` grid is the `(k + 1)`-point rule.  Because the
    /// one-dimensional rules are not nested only nodes with identical
    /// coordinates (such as the origin) are merged.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if the grid would be larger than
    /// [`MAX_QUADRATURE_NODES`] or a one-dimensional rule cannot be computed.
    pub fn smolyak(n_dim: usize, level: usize) -> Result<Self> {
        if n_dim == 0 {
            return Ok(Self::origin(0));
        }
        let too_large = || {
            Error::InvalidInput(format!(
                "the Smolyak grid of level {level} in {n_dim} dimensions exceeds the limits of \
                 {MAX_QUADRATURE_NODES} nodes and {MAX_QUADRATURE_COORDINATES} coordinates"
            ))
        };
        if level >= MAX_GAUSS_HERMITE_POINTS {
            return Err(Error::InvalidInput(format!(
                "a Smolyak grid of level {level} needs a Gauss-Hermite rule of {} points, more \
                 than the supported {MAX_GAUSS_HERMITE_POINTS}",
                level + 1
            )));
        }
        // Number of product-rule nodes before merging, in closed form: the
        // multi-indices with |m| = q contribute Σ Π (m_j + 1) = C(q + 2d - 1, q)
        // nodes, and only q >= level - (d - 1) has a non-zero coefficient.  The
        // limits are checked before anything is allocated.
        let generated: f64 = (level.saturating_sub(n_dim - 1)..=level)
            .map(|q| binomial(q + 2 * n_dim - 1, q))
            .sum();
        if generated > MAX_QUADRATURE_NODES as f64
            || generated * n_dim as f64 > MAX_QUADRATURE_COORDINATES as f64
        {
            return Err(too_large());
        }
        let rules: Vec<Vec<(f64, f64)>> = (0..=level)
            .map(|m| probabilists_gauss_hermite(m + 1))
            .collect::<Result<_>>()?;
        let mut builder = SmolyakBuilder {
            n_dim,
            rules: &rules,
            multi_index: vec![0; n_dim],
            lookup: HashMap::new(),
            nodes: Vec::new(),
            weights: Vec::new(),
            generated: 0,
        };
        builder.visit(0, level).map_err(|_| too_large())?;
        Ok(builder.finish())
    }

    /// Rule used by default for a pseudo-spectral projection with
    /// `points_per_dim` Gauss-Hermite points per dimension: the full tensor
    /// grid while it has at most `max_tensor_nodes` nodes, otherwise the
    /// Smolyak grid of level `points_per_dim - 1`, which is exact up to the
    /// same total degree `2·points_per_dim − 1` as the tensor rule is in each
    /// variable.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `points_per_dim` is zero or the selected
    /// grid is too large (see [`NormalQuadrature::tensor`] and
    /// [`NormalQuadrature::smolyak`]).
    pub fn for_pce(n_dim: usize, points_per_dim: usize, max_tensor_nodes: usize) -> Result<Self> {
        if points_per_dim == 0 {
            return Err(Error::InvalidInput(
                "a quadrature rule needs at least one point per dimension".to_string(),
            ));
        }
        match tensor_node_count(n_dim, points_per_dim) {
            Some(count) if count <= max_tensor_nodes => Self::tensor(n_dim, points_per_dim),
            _ => Self::smolyak(n_dim, points_per_dim - 1),
        }
    }
}

/// Number of nodes `n^d` of a tensor grid, `None` on overflow.
fn tensor_node_count(n_dim: usize, points_per_dim: usize) -> Option<usize> {
    let exponent = u32::try_from(n_dim).ok()?;
    points_per_dim.checked_pow(exponent)
}

/// 64-bit hash of the IEEE bit patterns of a coordinate tuple.
fn coordinate_hash(coordinates: &[f64]) -> u64 {
    let mut hash = 0xCBF2_9CE4_8422_2325_u64;
    for &x in coordinates {
        hash = (hash ^ x.to_bits()).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        hash ^= hash >> 32;
    }
    hash
}

/// Binomial coefficient `C(n, k)` in floating point (exact while it is below
/// 2⁵³, a sound magnitude estimate beyond).
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    (1..=k).fold(1.0, |acc, i| acc * (n - k + i) as f64 / i as f64)
}

/// Gauss-Hermite rule for the standard normal density: the physicists' rule
/// with `x = ξ/√2` and weights divided by `√π`.  Weights sum to one.
fn probabilists_gauss_hermite(n: usize) -> Result<Vec<(f64, f64)>> {
    if n == 0 || n > MAX_GAUSS_HERMITE_POINTS {
        return Err(Error::InvalidInput(format!(
            "Gauss-Hermite rules with 1 to {MAX_GAUSS_HERMITE_POINTS} points are supported, \
             got {n}"
        )));
    }
    let (nodes, weights) = gauss_hermite_quadrature(n);
    if nodes.len() != n || weights.len() != n {
        return Err(Error::InvalidInput(format!(
            "the {n}-point Gauss-Hermite rule could not be computed"
        )));
    }
    let normalisation = 1.0 / PI.sqrt();
    Ok(nodes
        .iter()
        .zip(weights.iter())
        .map(|(&x, &w)| (SQRT_2 * x, w * normalisation))
        .collect())
}

/// Accumulates the terms of the Smolyak combination formula.
struct SmolyakBuilder<'a> {
    n_dim: usize,
    /// `rules[m]` is the `(m + 1)`-point probabilists' Gauss-Hermite rule.
    rules: &'a [Vec<(f64, f64)>],
    multi_index: Vec<usize>,
    /// Hash of the coordinates of a node to the positions of the nodes with that
    /// hash (candidates, confirmed by comparing the stored coordinates).
    lookup: HashMap<u64, Vec<usize>>,
    nodes: Vec<f64>,
    weights: Vec<f64>,
    /// Product-rule nodes generated so far, before merging.
    generated: usize,
}

impl SmolyakBuilder<'_> {
    /// Enumerate every multi-index `m` with `|m| <= level` exactly once, as a
    /// set of positive entries at increasing positions (so the recursion depth
    /// is at most `min(level, d)`), and add the terms whose coefficient is
    /// non-zero.  `multi_index` holds the entries placed before `first_dim`
    /// and `remaining = level - |m|`.
    fn visit(&mut self, first_dim: usize, remaining: usize) -> Result<()> {
        // The combination coefficient vanishes unless level - |m| <= d - 1.
        if remaining < self.n_dim {
            self.add_term(remaining)?;
        }
        if remaining == 0 {
            return Ok(());
        }
        for dim in first_dim..self.n_dim {
            for m in 1..=remaining {
                self.multi_index[dim] = m;
                self.visit(dim + 1, remaining - m)?;
            }
            self.multi_index[dim] = 0;
        }
        Ok(())
    }

    /// Add `(-1)^deficit C(d-1, deficit) ⊗_j U_{m_j + 1}`.
    fn add_term(&mut self, deficit: usize) -> Result<()> {
        let sign = if deficit.is_multiple_of(2) { 1.0 } else { -1.0 };
        let coefficient = sign * binomial(self.n_dim - 1, deficit);
        let mut count = 1_usize;
        for &m in &self.multi_index {
            count = count
                .checked_mul(m + 1)
                .ok_or_else(|| Error::InvalidInput("quadrature size overflow".to_string()))?;
        }
        self.generated = self.generated.saturating_add(count);
        if self.generated > MAX_QUADRATURE_NODES
            || self
                .generated
                .checked_mul(self.n_dim)
                .is_none_or(|coordinates| coordinates > MAX_QUADRATURE_COORDINATES)
        {
            return Err(Error::InvalidInput("too many quadrature nodes".to_string()));
        }
        let d = self.n_dim;
        let mut digits = vec![0_usize; d];
        let mut coordinates = vec![0.0_f64; d];
        for _ in 0..count {
            let mut weight = coefficient;
            for (j, &digit) in digits.iter().enumerate() {
                let (x, w) = self.rules[self.multi_index[j]][digit];
                // `+ 0.0` turns a negative zero into a positive one.
                coordinates[j] = x + 0.0;
                weight *= w;
            }
            let hash = coordinate_hash(&coordinates);
            let existing = self.lookup.get(&hash).and_then(|candidates| {
                candidates
                    .iter()
                    .copied()
                    .find(|&q| self.nodes[q * d..(q + 1) * d] == coordinates[..])
            });
            match existing {
                Some(position) => self.weights[position] += weight,
                None => {
                    self.lookup
                        .entry(hash)
                        .or_default()
                        .push(self.weights.len());
                    self.nodes.extend_from_slice(&coordinates);
                    self.weights.push(weight);
                }
            }
            for (j, digit) in digits.iter_mut().enumerate() {
                *digit += 1;
                if *digit < self.multi_index[j] + 1 {
                    break;
                }
                *digit = 0;
            }
        }
        Ok(())
    }

    /// Drop nodes whose weights cancelled exactly and build the rule.
    fn finish(self) -> NormalQuadrature {
        let d = self.n_dim;
        let mut nodes = Vec::with_capacity(self.nodes.len());
        let mut weights = Vec::with_capacity(self.weights.len());
        for (q, &weight) in self.weights.iter().enumerate() {
            if weight != 0.0 {
                nodes.extend_from_slice(&self.nodes[q * d..(q + 1) * d]);
                weights.push(weight);
            }
        }
        NormalQuadrature {
            n_dim: d,
            nodes,
            weights,
        }
    }
}

/// Pseudo-spectral projection of `f` onto the Hermite chaos.
///
/// Returns, for every multi-index of `alphas`, the coefficient
/// `c_α = E[f(ξ) Ψ_α(ξ)] / α!` with the expectation evaluated on `rule`.
/// `f` is called exactly once per node of the rule.
///
/// The result is exact whenever `f·Ψ_α` is integrated exactly by the rule; in
/// particular a polynomial response of total degree `p` is reproduced
/// exactly by a tensor rule with `n ≥ p + 1` points per dimension and by a
/// Smolyak grid of level `k` with `2k + 1 ≥ 2p`.
///
/// # Errors
/// [`Error::InvalidInput`] if a multi-index does not have `rule.n_dim()`
/// entries or `f` returns a non-finite value.
pub fn pce_project<F>(f: F, alphas: &[Vec<usize>], rule: &NormalQuadrature) -> Result<Vec<f64>>
where
    F: Fn(&[f64]) -> f64,
{
    let n_dim = rule.n_dim();
    if let Some(bad) = alphas.iter().find(|alpha| alpha.len() != n_dim) {
        return Err(Error::InvalidInput(format!(
            "multi-index of length {} does not match the {n_dim} random dimensions",
            bad.len()
        )));
    }
    let max_degree = alphas
        .iter()
        .flat_map(|alpha| alpha.iter().copied())
        .max()
        .unwrap_or(0);
    let stride = max_degree + 1;
    // table[j * stride + k] = He_k(ξ_j) at the current node.
    let mut table = vec![0.0_f64; n_dim * stride];
    let mut coefficients = vec![0.0_f64; alphas.len()];
    for (q, (xi, weight)) in rule.iter().enumerate() {
        let value = f(xi);
        if !value.is_finite() {
            return Err(Error::InvalidInput(format!(
                "the response is not finite ({value}) at quadrature node {q}"
            )));
        }
        for (row, &x) in table.chunks_exact_mut(stride).zip(xi.iter()) {
            row[0] = 1.0;
            if stride > 1 {
                row[1] = x;
            }
            for k in 1..max_degree {
                row[k + 1] = x * row[k] - k as f64 * row[k - 1];
            }
        }
        let weighted = weight * value;
        for (coefficient, alpha) in coefficients.iter_mut().zip(alphas.iter()) {
            let psi: f64 = alpha
                .iter()
                .enumerate()
                .map(|(j, &degree)| table[j * stride + degree])
                .product();
            *coefficient += weighted * psi;
        }
    }
    for (coefficient, alpha) in coefficients.iter_mut().zip(alphas.iter()) {
        *coefficient /= alpha_factorial(alpha);
    }
    Ok(coefficients)
}

/// `α! = Π_j α_j!`, the squared norm of `Ψ_α`.
fn alpha_factorial(alpha: &[usize]) -> f64 {
    alpha
        .iter()
        .map(|&a| (1..=a).fold(1.0_f64, |acc, k| acc * k as f64))
        .product()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `E[Π ξ_j^{a_j}]` for independent standard normals.
    fn gaussian_moment(exponents: &[usize]) -> f64 {
        exponents
            .iter()
            .map(|&a| {
                if a % 2 == 1 {
                    0.0
                } else {
                    // (a - 1)!! = 1·3·5···(a - 1)
                    (1..a).step_by(2).fold(1.0_f64, |acc, k| acc * k as f64)
                }
            })
            .product()
    }

    /// All exponent tuples of length `d` with total degree at most `max_degree`.
    fn exponent_tuples(d: usize, max_degree: usize) -> Vec<Vec<usize>> {
        fn extend(d: usize, remaining: usize, current: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
            if current.len() == d {
                out.push(current.clone());
                return;
            }
            for a in 0..=remaining {
                current.push(a);
                extend(d, remaining - a, current, out);
                current.pop();
            }
        }
        let mut out = Vec::new();
        extend(d, max_degree, &mut Vec::new(), &mut out);
        out
    }

    fn monomial(exponents: &[usize], xi: &[f64]) -> f64 {
        exponents
            .iter()
            .zip(xi.iter())
            .map(|(&a, &x)| x.powi(a as i32))
            .product()
    }

    #[test]
    fn one_dimensional_rule_reproduces_normal_moments() {
        for n in 1..=8 {
            let rule = NormalQuadrature::tensor(1, n).expect("rule");
            assert_eq!(rule.len(), n);
            assert!((rule.weight_sum() - 1.0).abs() < 1e-13, "n = {n}");
            // Exact up to degree 2n - 1.
            for degree in 0..2 * n {
                let numeric = rule.integrate(|xi| xi[0].powi(degree as i32));
                let exact = gaussian_moment(&[degree]);
                assert!(
                    (numeric - exact).abs() < 1e-11 * exact.abs().max(1.0),
                    "n = {n}, degree {degree}: {numeric} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn three_point_rule_is_the_classical_one() {
        let rule = NormalQuadrature::tensor(1, 3).expect("rule");
        let sqrt3 = 3.0_f64.sqrt();
        let mut pairs: Vec<(f64, f64)> = rule.iter().map(|(xi, w)| (xi[0], w)).collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert!((pairs[0].0 + sqrt3).abs() < 1e-13 && (pairs[0].1 - 1.0 / 6.0).abs() < 1e-13);
        assert!(pairs[1].0.abs() < 1e-13 && (pairs[1].1 - 2.0 / 3.0).abs() < 1e-13);
        assert!((pairs[2].0 - sqrt3).abs() < 1e-13 && (pairs[2].1 - 1.0 / 6.0).abs() < 1e-13);
    }

    #[test]
    fn tensor_grid_has_full_node_count_and_is_exact_per_variable() {
        let rule = NormalQuadrature::tensor(3, 4).expect("rule");
        assert_eq!(rule.n_dim(), 3);
        assert_eq!(rule.len(), 64);
        assert!((rule.weight_sum() - 1.0).abs() < 1e-13);
        // Degree <= 7 in every variable.
        for exponents in [[7, 0, 0], [6, 4, 2], [2, 6, 6], [4, 4, 4], [1, 3, 5]] {
            let numeric = rule.integrate(|xi| monomial(&exponents, xi));
            let exact = gaussian_moment(&exponents);
            assert!(
                (numeric - exact).abs() < 1e-10 * exact.abs().max(1.0),
                "{exponents:?}: {numeric} vs {exact}"
            );
        }
    }

    #[test]
    fn tensor_grid_covers_the_off_diagonal_points() {
        // The former implementation evaluated only the points xi_1 = ... = xi_d.
        let rule = NormalQuadrature::tensor(2, 3).expect("rule");
        assert_eq!(rule.len(), 9);
        let off_diagonal = rule
            .iter()
            .filter(|(xi, _)| (xi[0] - xi[1]).abs() > 1e-9)
            .count();
        assert_eq!(off_diagonal, 6);
        let mixed = rule.integrate(|xi| xi[0] * xi[0] * xi[1] * xi[1]);
        assert!((mixed - 1.0).abs() < 1e-13);
    }

    #[test]
    fn zero_dimensional_rules_are_a_single_point() {
        for rule in [
            NormalQuadrature::tensor(0, 3).expect("tensor"),
            NormalQuadrature::smolyak(0, 4).expect("smolyak"),
        ] {
            assert_eq!(rule.len(), 1);
            assert_eq!(rule.n_dim(), 0);
            assert_eq!(rule.node(0), Some(&[][..]));
            assert!((rule.weight(0).expect("weight") - 1.0).abs() < 1e-15);
            assert!((rule.integrate(|_| 2.5) - 2.5).abs() < 1e-15);
        }
    }

    #[test]
    fn smolyak_grid_is_exact_for_total_degree_two_k_plus_one() {
        for (d, level) in [(2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (5, 1)] {
            let rule = NormalQuadrature::smolyak(d, level).expect("sparse grid");
            assert!(
                (rule.weight_sum() - 1.0).abs() < 1e-12,
                "d = {d}, k = {level}"
            );
            for exponents in exponent_tuples(d, 2 * level + 1) {
                let numeric = rule.integrate(|xi| monomial(&exponents, xi));
                let exact = gaussian_moment(&exponents);
                assert!(
                    (numeric - exact).abs() < 1e-9 * exact.abs().max(1.0),
                    "d = {d}, k = {level}, {exponents:?}: {numeric} vs {exact}"
                );
            }
        }
    }

    #[test]
    fn smolyak_grid_is_not_exact_beyond_its_degree() {
        // Total degree 2k + 2 is the first one the level-k grid cannot reach in general.
        let rule = NormalQuadrature::smolyak(2, 1).expect("sparse grid");
        let numeric = rule.integrate(|xi| xi[0] * xi[0] * xi[1] * xi[1]);
        assert!((numeric - 1.0).abs() > 1e-3, "got {numeric}");
    }

    #[test]
    fn smolyak_node_counts() {
        assert_eq!(NormalQuadrature::smolyak(3, 0).expect("k = 0").len(), 1);
        // 2d + 1 nodes: the origin and two points on every axis.
        assert_eq!(NormalQuadrature::smolyak(3, 1).expect("k = 1").len(), 7);
        assert_eq!(NormalQuadrature::smolyak(10, 1).expect("k = 1").len(), 21);
        // Origin + 6 two-point axis nodes + 6 three-point axis nodes + 12 pair nodes.
        assert_eq!(NormalQuadrature::smolyak(3, 2).expect("k = 2").len(), 25);
        // One dimension: the (k + 1)-point rule.
        for k in 0..6 {
            assert_eq!(NormalQuadrature::smolyak(1, k).expect("1-D").len(), k + 1);
        }
    }

    #[test]
    fn smolyak_equals_the_gauss_hermite_rule_in_one_dimension() {
        let sparse = NormalQuadrature::smolyak(1, 4).expect("sparse");
        let full = NormalQuadrature::tensor(1, 5).expect("tensor");
        assert_eq!(sparse.len(), full.len());
        let mut a: Vec<(f64, f64)> = sparse.iter().map(|(xi, w)| (xi[0], w)).collect();
        let mut b: Vec<(f64, f64)> = full.iter().map(|(xi, w)| (xi[0], w)).collect();
        a.sort_by(|x, y| x.0.total_cmp(&y.0));
        b.sort_by(|x, y| x.0.total_cmp(&y.0));
        for (p, q) in a.iter().zip(b.iter()) {
            assert!((p.0 - q.0).abs() < 1e-13 && (p.1 - q.1).abs() < 1e-13);
        }
    }

    #[test]
    fn sparse_grid_is_far_smaller_than_the_tensor_grid() {
        let sparse = NormalQuadrature::smolyak(10, 3).expect("sparse");
        let tensor_count = tensor_node_count(10, 4).expect("count");
        assert!(sparse.len() < 5_000, "{} nodes", sparse.len());
        assert!(tensor_count > 1_000_000);
        assert!((sparse.weight_sum() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn oversized_rules_are_rejected() {
        assert!(matches!(
            NormalQuadrature::tensor(30, 5),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            NormalQuadrature::tensor(3, 0),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            NormalQuadrature::smolyak(500, 5),
            Err(Error::InvalidInput(_))
        ));
        // Few nodes but far too many stored coordinates (2d + 1 nodes of d entries).
        assert!(matches!(
            NormalQuadrature::smolyak(100_000, 1),
            Err(Error::InvalidInput(_))
        ));
        // One-dimensional rules are limited to MAX_GAUSS_HERMITE_POINTS points.
        assert!(matches!(
            NormalQuadrature::tensor(1, MAX_GAUSS_HERMITE_POINTS + 1),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            NormalQuadrature::smolyak(1, MAX_GAUSS_HERMITE_POINTS),
            Err(Error::InvalidInput(_))
        ));
        assert!(NormalQuadrature::tensor(1, MAX_GAUSS_HERMITE_POINTS).is_ok());
    }

    #[test]
    fn level_one_grid_in_many_dimensions_has_two_d_plus_one_nodes() {
        // The origin and two points on each axis: 2001 nodes of 1000 coordinates.
        let rule = NormalQuadrature::smolyak(1_000, 1).expect("sparse grid");
        assert_eq!(rule.len(), 2 * 1_000 + 1);
        assert_eq!(rule.n_dim(), 1_000);
        assert!((rule.weight_sum() - 1.0).abs() < 1e-9);
        // Exact for the quadratic E[xi_17^2] = 1 (degree 2 <= 2k + 1 = 3).
        assert!((rule.integrate(|xi| xi[17] * xi[17]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn node_count_estimate_matches_the_generated_nodes() {
        // Before merging the grid has Σ_q C(q + 2d - 1, q) nodes; the origin is the
        // only coordinate tuple that several terms share.
        for (d, level, merged) in [(3, 2, 25), (2, 1, 5), (4, 2, 41), (1, 3, 4)] {
            let rule = NormalQuadrature::smolyak(d, level).expect("sparse grid");
            let generated: f64 = (level.saturating_sub(d - 1)..=level)
                .map(|q| binomial(q + 2 * d - 1, q))
                .sum();
            assert!(rule.len() as f64 <= generated, "d = {d}, k = {level}");
            assert_eq!(rule.len(), merged, "d = {d}, k = {level}");
        }
        assert_eq!(binomial(5, 0) + binomial(6, 1) + binomial(7, 2), 28.0);
    }

    #[test]
    fn for_pce_switches_to_the_sparse_grid_above_the_threshold() {
        // 4^3 = 64 nodes fit a threshold of 64 ...
        let tensor = NormalQuadrature::for_pce(3, 4, 64).expect("tensor");
        assert_eq!(tensor.len(), 64);
        // ... but not a threshold of 63: level 3 sparse grid instead.
        let sparse = NormalQuadrature::for_pce(3, 4, 63).expect("sparse");
        let reference = NormalQuadrature::smolyak(3, 3).expect("reference");
        assert_eq!(sparse.len(), reference.len());
        assert!(sparse.len() != 64);
        // Overflowing tensor sizes also fall back to the sparse grid.
        let huge = NormalQuadrature::for_pce(40, 3, DEFAULT_MAX_TENSOR_NODES).expect("sparse");
        assert_eq!(huge.n_dim(), 40);
        assert!((huge.weight_sum() - 1.0).abs() < 1e-9);
    }

    fn total_degree_indices(d: usize, order: usize) -> Vec<Vec<usize>> {
        exponent_tuples(d, order)
    }

    #[test]
    fn projection_recovers_the_chaos_coefficients_of_a_polynomial() {
        // f = 1 + 2 xi_1 + 3 xi_1 xi_2 = He_0 + 2 He_1(xi_1) + 3 He_1(xi_1) He_1(xi_2).
        let alphas = total_degree_indices(2, 2);
        let f = |xi: &[f64]| 1.0 + 2.0 * xi[0] + 3.0 * xi[0] * xi[1];
        let tensor = NormalQuadrature::tensor(2, 3).expect("tensor");
        let sparse = NormalQuadrature::smolyak(2, 2).expect("sparse");
        for (name, rule) in [("tensor", &tensor), ("sparse", &sparse)] {
            let coefficients = pce_project(f, &alphas, rule).expect("projection");
            for (alpha, c) in alphas.iter().zip(coefficients.iter()) {
                let expected = match alpha.as_slice() {
                    [0, 0] => 1.0,
                    [1, 0] => 2.0,
                    [1, 1] => 3.0,
                    _ => 0.0,
                };
                assert!(
                    (c - expected).abs() < 1e-13,
                    "{name}: coefficient of {alpha:?} is {c}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn projection_uses_the_hermite_not_the_monomial_basis() {
        // xi^3 = He_3 + 3 He_1 and xi_1^2 xi_2 = (He_2 + 1) He_1(xi_2) in three variables.
        let alphas = total_degree_indices(3, 3);
        let f = |xi: &[f64]| xi[2].powi(3) + xi[0] * xi[0] * xi[1];
        let rule = NormalQuadrature::tensor(3, 4).expect("tensor");
        let coefficients = pce_project(f, &alphas, &rule).expect("projection");
        for (alpha, c) in alphas.iter().zip(coefficients.iter()) {
            let expected = match alpha.as_slice() {
                [0, 0, 3] => 1.0,
                [0, 0, 1] => 3.0,
                [2, 1, 0] => 1.0,
                [0, 1, 0] => 1.0,
                _ => 0.0,
            };
            assert!(
                (c - expected).abs() < 1e-12,
                "coefficient of {alpha:?} is {c}, expected {expected}"
            );
        }
        // Parseval: E[f^2] equals the sum of c_alpha^2 alpha! for a polynomial in the span.
        let second_moment = rule.integrate(|xi| f(xi) * f(xi));
        let spectral: f64 = coefficients
            .iter()
            .zip(alphas.iter())
            .map(|(c, alpha)| c * c * alpha_factorial(alpha))
            .sum();
        assert!((second_moment - spectral).abs() < 1e-9 * second_moment);
    }

    #[test]
    fn sparse_and_tensor_projections_agree_for_a_smooth_response() {
        let alphas = total_degree_indices(3, 2);
        let f = |xi: &[f64]| (0.3 * xi[0] - 0.2 * xi[1] + 0.1 * xi[2]).exp();
        let tensor = NormalQuadrature::tensor(3, 5).expect("tensor");
        let sparse = NormalQuadrature::smolyak(3, 4).expect("sparse");
        let a = pce_project(f, &alphas, &tensor).expect("tensor projection");
        let b = pce_project(f, &alphas, &sparse).expect("sparse projection");
        // exp(c . xi) = exp(|c|^2 / 2) Prod_j sum_n c_j^n He_n(xi_j) / n!, hence
        // the coefficient of Psi_alpha is exp(|c|^2 / 2) Prod_j c_j^alpha_j / alpha_j!.
        let c = [0.3, -0.2, 0.1];
        let scale = (0.5_f64 * c.iter().map(|v| v * v).sum::<f64>()).exp();
        for ((alpha, x), y) in alphas.iter().zip(a.iter()).zip(b.iter()) {
            let analytic = scale
                * alpha
                    .iter()
                    .zip(c.iter())
                    .map(|(&n, &cj)| cj.powi(n as i32) / alpha_factorial(&[n]))
                    .product::<f64>();
            assert!(
                (x - analytic).abs() < 3e-6,
                "tensor {alpha:?}: {x} vs {analytic}"
            );
            assert!(
                (y - analytic).abs() < 3e-6,
                "sparse {alpha:?}: {y} vs {analytic}"
            );
            assert!((x - y).abs() < 6e-6, "{alpha:?}: {x} vs {y}");
        }
    }

    #[test]
    fn projection_validates_its_input() {
        let rule = NormalQuadrature::tensor(2, 3).expect("tensor");
        assert!(matches!(
            pce_project(|_| 1.0, &[vec![0, 0, 1]], &rule),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            pce_project(|_| f64::NAN, &[vec![0, 0]], &rule),
            Err(Error::InvalidInput(_))
        ));
        // An empty list of multi-indices has an empty projection.
        let none = pce_project(|_| 1.0, &[], &rule).expect("empty projection");
        assert!(none.is_empty());
    }

    #[test]
    fn node_accessors_are_bounds_checked() {
        let rule = NormalQuadrature::tensor(2, 2).expect("tensor");
        assert_eq!(rule.len(), 4);
        assert!(!rule.is_empty());
        assert!(rule.node(3).is_some());
        assert!(rule.node(4).is_none());
        assert!(rule.weight(4).is_none());
        assert_eq!(rule.iter().count(), 4);
    }
}
