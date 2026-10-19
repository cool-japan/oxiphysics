// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The Fiedler vector of a [`Subgraph`] by restarted Lanczos iteration.
//!
//! The Fiedler vector is the eigenvector of the Laplacian `L = D - A` that
//! belongs to its second-smallest eigenvalue `lambda_2`.  The smallest
//! eigenvalue is `0`, with the constant vector as eigenvector.  The constant
//! vector is *deflated*: the start vector and every product `L q` are projected
//! onto the space of mean-free vectors, in which `lambda_2` is the smallest
//! eigenvalue.
//!
//! # Why Lanczos and not a power iteration
//!
//! A power iteration on `L` itself converges to the *largest* eigenvector,
//! which is the most oscillating one and carries no community structure.  The
//! smallest eigenpairs are the largest ones of `M = c I - L` for any
//! `c >= lambda_max`, and `c = 2 * max_degree` is such a bound (Gershgorin).
//! The power iteration on `M` is correct but converges at the rate
//! `(c - lambda_3) / (c - lambda_2)`.  A graph with several communities has a
//! whole cluster of very small eigenvalues, so the rate is within `1e-5` of `1`
//! and the iteration needs on the order of `10^5` products: with the iteration
//! capped at `10^4` products, bisecting a graph of 2000 nodes and 20
//! communities took about a second.
//!
//! Lanczos iteration builds an orthonormal basis of the Krylov space
//! `span{q, L q, L^2 q, ...}` in which `L` is tridiagonal, and takes the
//! smallest eigenpair of that small matrix (a *Ritz pair*) as the
//! approximation.  The number of products it needs scales with the square root
//! of the one above: tens to a few hundred for graphs with community
//! structure.  Every basis vector is orthogonalised against all previous ones
//! (full reorthogonalisation), so round-off cannot spoil the basis, and the
//! basis stays short (a cycle of [`KRYLOV_DIMENSION`] steps); the iteration is
//! then restarted from the current Ritz vector.  Graphs with a huge diameter
//! (long paths, thin lattices) have a tiny `lambda_2` and converge slowly under
//! any Krylov method; [`MAX_PRODUCTS`] bounds the work for them.

use super::Subgraph;
use super::jacobi::symmetric_eigen;

/// Number of Lanczos steps (and so basis vectors) per cycle.
const KRYLOV_DIMENSION: usize = 40;

/// Bound on the Laplacian products of one Fiedler vector computation.  A
/// cluster of nearly equal eigenvalues cannot be resolved cheaply, but any
/// vector in it separates the communities equally well, so the best Ritz
/// vector found within the budget is returned.
const MAX_PRODUCTS: usize = 2_000;

/// A Ritz pair is accepted once its residual `||L y - theta y||` is below this
/// fraction of the spectral radius bound `2 * max_degree`.
///
/// The same threshold decides when the Krylov space is exhausted: the residual
/// of every Ritz pair is at most `beta`, the length of the part of `L q`
/// outside the space, so once `beta` is below the tolerance all Ritz pairs are
/// eigenpairs to that accuracy (graphs with repeated eigenvalues, such as
/// those with many structurally identical nodes, reach this after fewer than
/// `size - 1` steps).
const RESIDUAL_TOLERANCE: f64 = 1e-8;

/// A vector is orthogonalised a second time when the first pass shrank it by
/// more than this factor (Daniel-Gragg-Kaufman-Stewart criterion).
const REORTHOGONALISE_BELOW: f64 = 0.7;

/// Fractional part of the golden ratio.  Its multiples give a deterministic
/// start vector that is not symmetric under any permutation of the nodes, so
/// it is not orthogonal to the Fiedler vector of a symmetric graph.
const GOLDEN_RATIO_CONJUGATE: f64 = 0.618_033_988_749_894_9;

/// Result of one Lanczos cycle.
struct Cycle {
    /// Unit-length, mean-free Ritz vector of the smallest Ritz value.
    vector: Vec<f64>,
    /// Norm of the residual `L y - theta y` of that Ritz pair.
    residual: f64,
    /// Laplacian products spent.
    products: usize,
    /// Whether the Krylov space was exhausted, i.e. the pair is exact.
    exact: bool,
}

/// Approximate Fiedler vector of `sub`: a unit vector orthogonal to the
/// constant vector that (nearly) minimises the Rayleigh quotient
/// `x^T L x / x^T x`.
///
/// A graph without edges has no spectral structure; its (deterministic,
/// normalised, mean-free) start vector is returned.  The result is
/// deterministic.
pub(super) fn fiedler_vector(sub: &Subgraph) -> Vec<f64> {
    let mut products = 0;
    fiedler_vector_counting(sub, &mut products)
}

/// [`fiedler_vector`] that also reports, in `products`, how many Laplacian
/// products it spent.
pub(super) fn fiedler_vector_counting(sub: &Subgraph, products: &mut usize) -> Vec<f64> {
    *products = 0;
    let size = sub.len();
    if size < 2 {
        return vec![0.0; size];
    }
    let mut vector: Vec<f64> = (0..size)
        .map(|i| ((i + 1) as f64 * GOLDEN_RATIO_CONJUGATE).fract() - 0.5)
        .collect();
    if !deflate_and_normalise(&mut vector) {
        return vector;
    }
    let max_degree = sub.degree.iter().copied().fold(0.0_f64, f64::max);
    if max_degree <= 0.0 {
        return vector;
    }
    let scale = 2.0 * max_degree;
    // The mean-free vectors form a space of dimension `size - 1`.
    let dimension = KRYLOV_DIMENSION.min(size - 1);
    loop {
        let cycle = lanczos_cycle(sub, &vector, dimension, scale);
        *products += cycle.products;
        vector = cycle.vector;
        let converged = cycle.exact || cycle.residual <= RESIDUAL_TOLERANCE * scale;
        if converged || *products + dimension > MAX_PRODUCTS {
            return vector;
        }
    }
}

/// One Lanczos cycle of at most `dimension` steps started from the unit,
/// mean-free vector `start`.  `scale` bounds the spectral radius of `L`.
fn lanczos_cycle(sub: &Subgraph, start: &[f64], dimension: usize, scale: f64) -> Cycle {
    let size = sub.len();
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(dimension);
    let mut alphas: Vec<f64> = Vec::with_capacity(dimension);
    let mut betas: Vec<f64> = Vec::with_capacity(dimension);
    // Norm of the part of `L q_last` outside the Krylov space: the residual of
    // every Ritz pair is this times the last component of its eigenvector.
    let mut boundary = 0.0_f64;
    let mut exact = false;
    let mut vector = start.to_vec();
    let mut work = vec![0.0_f64; size];
    for step in 0..dimension {
        sub.laplacian_product(&vector, &mut work);
        alphas.push(dot(&vector, &work));
        basis.push(vector);
        // Subtracting the projections on all basis vectors performs the
        // three-term Lanczos recurrence and reorthogonalises in one go.  The
        // constant vector, which is deflated, takes part in every pass: a
        // component along it would otherwise grow with the 1 / beta scaling
        // of a short vector and show up as a bogus Ritz value of 0.
        let before = norm(&work);
        orthogonalise(&mut work, &basis);
        if norm(&work) < REORTHOGONALISE_BELOW * before {
            orthogonalise(&mut work, &basis);
        }
        let beta = norm(&work);
        if beta <= RESIDUAL_TOLERANCE * scale {
            exact = true;
            boundary = beta;
            break;
        }
        boundary = beta;
        if step + 1 == dimension {
            break;
        }
        betas.push(beta);
        vector = work.iter().map(|entry| entry / beta).collect();
    }

    // Smallest eigenpair of the tridiagonal matrix T = Q^T L Q.
    let k = alphas.len();
    let mut tridiagonal = vec![vec![0.0_f64; k]; k];
    for (j, &alpha) in alphas.iter().enumerate() {
        tridiagonal[j][j] = alpha;
    }
    for (j, &beta) in betas.iter().enumerate() {
        tridiagonal[j][j + 1] = beta;
        tridiagonal[j + 1][j] = beta;
    }
    let (_, eigenvectors) = symmetric_eigen(tridiagonal);
    let coefficients: &[f64] = match eigenvectors.first() {
        Some(eigenvector) => eigenvector.as_slice(),
        None => &[],
    };
    let mut ritz = vec![0.0_f64; size];
    for (&coefficient, basis_vector) in coefficients.iter().zip(&basis) {
        for (slot, &entry) in ritz.iter_mut().zip(basis_vector) {
            *slot += coefficient * entry;
        }
    }
    let usable = deflate_and_normalise(&mut ritz);
    let residual = boundary * coefficients.last().map_or(0.0, |c| c.abs());
    Cycle {
        vector: if usable { ritz } else { start.to_vec() },
        residual: if usable { residual } else { f64::INFINITY },
        products: k,
        exact: exact && usable,
    }
}

/// Subtract the mean of `v`.
fn remove_mean(v: &mut [f64]) {
    let mean = v.iter().sum::<f64>() / v.len().max(1) as f64;
    for entry in v.iter_mut() {
        *entry -= mean;
    }
}

/// Vectors that keep less than this fraction of their length after the mean
/// is removed are numerically constant: what is left is rounding noise.
const MIN_MEAN_FREE_FRACTION: f64 = 1e-8;

/// Subtract the mean of `v` and scale it to unit length.  Returns `false`
/// (leaving `v` mean-free but unscaled) when nothing is left, i.e. `v` was
/// (numerically) constant or not finite.
fn deflate_and_normalise(v: &mut [f64]) -> bool {
    let original = norm(v);
    remove_mean(v);
    let length = norm(v);
    if !length.is_finite()
        || !original.is_finite()
        || length <= f64::MIN_POSITIVE
        || length < MIN_MEAN_FREE_FRACTION * original
    {
        return false;
    }
    for entry in v.iter_mut() {
        *entry /= length;
    }
    true
}

/// Remove from `w` its components along the constant vector and along the
/// (orthonormal, mean-free) `basis` vectors.
fn orthogonalise(w: &mut [f64], basis: &[Vec<f64>]) {
    remove_mean(w);
    for b in basis {
        let coefficient = dot(b, w);
        for (entry, &component) in w.iter_mut().zip(b) {
            *entry -= coefficient * component;
        }
    }
    remove_mean(w);
}

/// Euclidean inner product.
///
/// Four independent partial sums let the processor overlap the additions of a
/// long reduction (a single running sum is a dependency chain); the summation
/// order is fixed, so the result is deterministic.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut lanes = [0.0_f64; 4];
    let mut a_chunks = a.chunks_exact(4);
    let mut b_chunks = b.chunks_exact(4);
    for (x, y) in (&mut a_chunks).zip(&mut b_chunks) {
        for lane in 0..4 {
            lanes[lane] += x[lane] * y[lane];
        }
    }
    let tail: f64 = a_chunks
        .remainder()
        .iter()
        .zip(b_chunks.remainder())
        .map(|(x, y)| x * y)
        .sum();
    (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]) + tail
}

/// Euclidean norm.
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
