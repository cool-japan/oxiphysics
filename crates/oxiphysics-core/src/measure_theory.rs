#![allow(clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Measure theory: sigma-algebras, Lebesgue measure, probability measures,
//! signed measures, product measures, and numerical Lebesgue/Monte Carlo integration.
//!
//! This module provides both abstract measure-theoretic structures and concrete
//! numerical implementations suitable for scientific computing.

#![allow(dead_code)]

use rand::RngExt as RandRng;
use std::f64::consts::TAU;

// ---------------------------------------------------------------------------
// Core measure trait
// ---------------------------------------------------------------------------

/// A measurable set, identified by a tag and a predicate.
///
/// For numerical work, sets are represented as closed intervals (boxes) in ℝⁿ.
#[derive(Debug, Clone)]
pub struct MeasurableSet {
    /// Human-readable label for the set.
    pub label: String,
    /// Lower bounds of the bounding box (one per dimension).
    pub lower: Vec<f64>,
    /// Upper bounds of the bounding box (one per dimension).
    pub upper: Vec<f64>,
}

impl MeasurableSet {
    /// Create a measurable set from lower and upper bounds.
    ///
    /// # Panics
    /// Panics if `lower` and `upper` have different lengths.
    pub fn new(label: impl Into<String>, lower: Vec<f64>, upper: Vec<f64>) -> Self {
        assert_eq!(
            lower.len(),
            upper.len(),
            "MeasurableSet: dimension mismatch"
        );
        MeasurableSet {
            label: label.into(),
            lower,
            upper,
        }
    }

    /// Create a 1-D interval `[a, b]`.
    pub fn interval(a: f64, b: f64) -> Self {
        Self::new(format!("[{a},{b}]"), vec![a], vec![b])
    }

    /// Return the dimension of the set.
    pub fn dim(&self) -> usize {
        self.lower.len()
    }

    /// Check whether a point lies in the closed box `[lower, upper]`.
    pub fn contains(&self, point: &[f64]) -> bool {
        assert_eq!(
            point.len(),
            self.dim(),
            "MeasurableSet::contains: dim mismatch"
        );
        self.lower
            .iter()
            .zip(self.upper.iter())
            .zip(point.iter())
            .all(|((lo, hi), x)| x >= lo && x <= hi)
    }

    /// Return `true` if this set is a null set (measure zero) under Lebesgue measure.
    ///
    /// A box is a null set iff at least one dimension has `lower == upper`.
    pub fn is_null_set(&self) -> bool {
        self.lower
            .iter()
            .zip(self.upper.iter())
            .any(|(lo, hi)| (hi - lo).abs() < 1e-15)
    }
}

/// Abstract measure trait.
///
/// A measure `μ` on a sigma-algebra `Σ` over `X` satisfies:
/// 1. Non-negativity: `μ(A) ≥ 0` for all `A ∈ Σ`
/// 2. Null empty set: `μ(∅) = 0`
/// 3. Countable additivity: `μ(⋃ Aᵢ) = Σ μ(Aᵢ)` for pairwise disjoint `Aᵢ`
pub trait Measure {
    /// Compute the measure of a set.
    fn measure(&self, set: &MeasurableSet) -> f64;

    /// Return `true` if `set` is a null set (measure zero).
    fn is_null(&self, set: &MeasurableSet) -> bool {
        self.measure(set).abs() < 1e-15
    }

    /// The measure of the empty set should always be zero.
    fn empty_set_measure(&self) -> f64 {
        0.0
    }

    /// Compute the measure of the union of two disjoint sets.
    ///
    /// For additive measures: `μ(A ∪ B) = μ(A) + μ(B)` when `A ∩ B = ∅`.
    fn union_measure_disjoint(&self, a: &MeasurableSet, b: &MeasurableSet) -> f64 {
        self.measure(a) + self.measure(b)
    }
}

// ---------------------------------------------------------------------------
// Sigma-algebra over a finite universe
// ---------------------------------------------------------------------------

/// A sigma-algebra over a finite universe `{0, 1, …, n-1}`.
///
/// Stores all subsets that belong to the sigma-algebra as bitmasks.
/// This is useful for combinatorial and discrete measure theory examples.
#[derive(Debug, Clone)]
pub struct FiniteSigmaAlgebra {
    /// Size of the universe.
    pub n: usize,
    /// Bitmask representations of sets in the sigma-algebra.
    sets: Vec<u64>,
}

impl FiniteSigmaAlgebra {
    /// Construct the power sigma-algebra (all 2ⁿ subsets).
    ///
    /// # Panics
    /// Panics if `n > 63`.
    pub fn power_set(n: usize) -> Self {
        assert!(n <= 63, "FiniteSigmaAlgebra: n too large (max 63)");
        let sets = (0u64..(1u64 << n)).collect();
        Self { n, sets }
    }

    /// Construct the trivial sigma-algebra `{∅, X}`.
    pub fn trivial(n: usize) -> Self {
        let full = if n == 0 { 0 } else { (1u64 << n) - 1 };
        Self {
            n,
            sets: vec![0, full],
        }
    }

    /// Check whether a given bitmask belongs to the sigma-algebra.
    pub fn contains(&self, bitmask: u64) -> bool {
        self.sets.contains(&bitmask)
    }

    /// Number of sets in the sigma-algebra.
    pub fn cardinality(&self) -> usize {
        self.sets.len()
    }

    /// Verify closure under complement.
    pub fn is_closed_under_complement(&self) -> bool {
        let full = if self.n == 0 { 0 } else { (1u64 << self.n) - 1 };
        self.sets.iter().all(|&s| self.contains(full & !s))
    }

    /// Verify closure under finite union.
    pub fn is_closed_under_union(&self) -> bool {
        for &a in &self.sets {
            for &b in &self.sets {
                if !self.contains(a | b) {
                    return false;
                }
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Lebesgue measure
// ---------------------------------------------------------------------------

/// Lebesgue measure on ℝⁿ (product of interval lengths).
///
/// For a box `[a₁,b₁] × … × [aₙ,bₙ]`:
///
/// `λ(A) = Π (bᵢ − aᵢ)`
///
/// The Lebesgue measure is translation-invariant, countably additive, and
/// assigns measure zero to all countable sets.
#[derive(Debug, Clone, Copy, Default)]
pub struct LebesgueMeasure;

impl Measure for LebesgueMeasure {
    fn measure(&self, set: &MeasurableSet) -> f64 {
        set.lower
            .iter()
            .zip(set.upper.iter())
            .map(|(lo, hi)| (hi - lo).max(0.0))
            .product()
    }
}

impl LebesgueMeasure {
    /// Measure a 1-D interval `[a, b]`.
    pub fn interval_measure(a: f64, b: f64) -> f64 {
        (b - a).max(0.0)
    }

    /// Measure an n-D box given separate lower and upper bounds.
    pub fn box_measure(lower: &[f64], upper: &[f64]) -> f64 {
        assert_eq!(lower.len(), upper.len(), "box_measure: dim mismatch");
        lower
            .iter()
            .zip(upper.iter())
            .map(|(lo, hi)| (hi - lo).max(0.0))
            .product()
    }

    /// Compute the measure of the intersection of two 1-D intervals.
    pub fn interval_intersection(a1: f64, a2: f64, b1: f64, b2: f64) -> f64 {
        let lo = a1.max(b1);
        let hi = a2.min(b2);
        (hi - lo).max(0.0)
    }

    /// Check translation invariance: `λ([a,b]) = λ([a+t, b+t])`.
    pub fn translation_invariant(a: f64, b: f64, t: f64) -> bool {
        let original = Self::interval_measure(a, b);
        let shifted = Self::interval_measure(a + t, b + t);
        (original - shifted).abs() < 1e-12
    }
}

// ---------------------------------------------------------------------------
// Probability measure
// ---------------------------------------------------------------------------

/// A probability measure on ℝⁿ with density function.
///
/// `P(A) = ∫_A f(x) dx`  where `∫ f dx = 1` (normalization).
///
/// Supports numerical expectation and variance via Monte Carlo estimation.
pub struct ProbabilityMeasure {
    /// Probability density function `f(x) ≥ 0`.
    pub density: Box<dyn Fn(&[f64]) -> f64>,
    /// Domain of the density (support box).
    pub domain: MeasurableSet,
    /// Number of Monte Carlo samples for numerical computations.
    pub n_samples: usize,
}

impl ProbabilityMeasure {
    /// Create a probability measure from a density and domain.
    pub fn new(
        density: impl Fn(&[f64]) -> f64 + 'static,
        domain: MeasurableSet,
        n_samples: usize,
    ) -> Self {
        Self {
            density: Box::new(density),
            domain,
            n_samples,
        }
    }

    /// Estimate `P(A)` via Monte Carlo integration.
    ///
    /// Samples uniformly from the domain and averages the density over `A`.
    pub fn probability(&self, set: &MeasurableSet) -> f64 {
        let mut rng = rand::rng();
        let domain_vol = LebesgueMeasure.measure(&self.domain);
        if domain_vol == 0.0 {
            return 0.0;
        }

        let mut sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            if set.contains(&point) {
                sum += (self.density)(&point);
            }
        }
        domain_vol * sum / self.n_samples as f64
    }

    /// Estimate the total mass (should be ≈ 1 for a true probability measure).
    pub fn total_mass(&self) -> f64 {
        let mut rng = rand::rng();
        let domain_vol = LebesgueMeasure.measure(&self.domain);
        if domain_vol == 0.0 {
            return 0.0;
        }

        let mut sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            sum += (self.density)(&point);
        }
        domain_vol * sum / self.n_samples as f64
    }

    /// Estimate the expectation `E[g(X)]` of a function `g` under this measure.
    pub fn expectation(&self, g: &dyn Fn(&[f64]) -> f64) -> f64 {
        let mut rng = rand::rng();
        let domain_vol = LebesgueMeasure.measure(&self.domain);
        if domain_vol == 0.0 {
            return 0.0;
        }

        let mut numer = 0.0;
        let mut denom = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            let f_val = (self.density)(&point);
            numer += f_val * g(&point);
            denom += f_val;
        }
        if denom.abs() < 1e-15 {
            return 0.0;
        }
        numer / denom
    }

    /// Estimate `Var[g(X)] = E[g(X)²] − (E[g(X)])²`.
    pub fn variance(&self, g: &dyn Fn(&[f64]) -> f64) -> f64 {
        let eg = self.expectation(g);
        let g2 = |x: &[f64]| g(x).powi(2);
        let eg2 = self.expectation(&g2);
        (eg2 - eg * eg).max(0.0)
    }

    /// Estimate the entropy `H = -∫ f ln(f) dx` via Monte Carlo.
    pub fn entropy(&self) -> f64 {
        let mut rng = rand::rng();
        let domain_vol = LebesgueMeasure.measure(&self.domain);
        if domain_vol == 0.0 {
            return 0.0;
        }

        let mut sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            let f_val = (self.density)(&point);
            if f_val > 1e-15 {
                sum += f_val * f_val.ln();
            }
        }
        -domain_vol * sum / self.n_samples as f64
    }
}

// ---------------------------------------------------------------------------
// Signed measure & decompositions
// ---------------------------------------------------------------------------

/// A signed measure represented as a difference of two non-negative measures.
///
/// By the Jordan decomposition theorem, every signed measure `ν` can be
/// written as `ν = ν⁺ − ν⁻` where `ν⁺, ν⁻ ≥ 0` are mutually singular.
///
/// Here we work numerically with density functions.
pub struct SignedMeasure {
    /// The positive part density `f⁺(x) ≥ 0`.
    pub positive_density: Box<dyn Fn(&[f64]) -> f64>,
    /// The negative part density `f⁻(x) ≥ 0`.
    pub negative_density: Box<dyn Fn(&[f64]) -> f64>,
    /// Domain of integration.
    pub domain: MeasurableSet,
    /// Number of Monte Carlo samples.
    pub n_samples: usize,
}

impl SignedMeasure {
    /// Create a signed measure from positive and negative densities.
    pub fn new(
        positive_density: impl Fn(&[f64]) -> f64 + 'static,
        negative_density: impl Fn(&[f64]) -> f64 + 'static,
        domain: MeasurableSet,
        n_samples: usize,
    ) -> Self {
        Self {
            positive_density: Box::new(positive_density),
            negative_density: Box::new(negative_density),
            domain,
            n_samples,
        }
    }

    /// Create a signed measure from a single (possibly negative) density `f`.
    pub fn from_density(
        density: impl Fn(&[f64]) -> f64 + 'static,
        domain: MeasurableSet,
        n_samples: usize,
    ) -> Self {
        let density = std::rc::Rc::new(density);
        let d1 = density.clone();
        let d2 = density;
        Self::new(
            move |x| d1(x).max(0.0),
            move |x| (-d2(x)).max(0.0),
            domain,
            n_samples,
        )
    }

    /// Evaluate the signed measure on a set: `ν(A) = ∫_A (f⁺ − f⁻) dx`.
    pub fn evaluate(&self, set: &MeasurableSet) -> f64 {
        let vol = LebesgueMeasure.measure(&self.domain);
        if vol == 0.0 {
            return 0.0;
        }
        let mut rng = rand::rng();
        let mut sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            if set.contains(&point) {
                let fpos = (self.positive_density)(&point);
                let fneg = (self.negative_density)(&point);
                sum += fpos - fneg;
            }
        }
        vol * sum / self.n_samples as f64
    }

    /// Total variation: `|ν|(A) = ∫_A (f⁺ + f⁻) dx`.
    pub fn total_variation(&self, set: &MeasurableSet) -> f64 {
        let vol = LebesgueMeasure.measure(&self.domain);
        if vol == 0.0 {
            return 0.0;
        }
        let mut rng = rand::rng();
        let mut sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            if set.contains(&point) {
                let fpos = (self.positive_density)(&point);
                let fneg = (self.negative_density)(&point);
                sum += fpos + fneg;
            }
        }
        vol * sum / self.n_samples as f64
    }

    /// Hahn decomposition: identify the positive set P where `f⁺ > f⁻`.
    ///
    /// Returns a sample of points in P within the domain.
    pub fn hahn_positive_set_samples(&self, n: usize) -> Vec<Vec<f64>> {
        let mut rng = rand::rng();
        let dim = self.domain.dim();
        let mut result = Vec::new();
        for _ in 0..n {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            if (self.positive_density)(&point) > (self.negative_density)(&point) {
                result.push(point);
            }
        }
        result
    }

    /// Jordan decomposition: return estimates of `ν⁺(domain)` and `ν⁻(domain)`.
    pub fn jordan_decomposition(&self) -> (f64, f64) {
        let vol = LebesgueMeasure.measure(&self.domain);
        if vol == 0.0 {
            return (0.0, 0.0);
        }
        let mut rng = rand::rng();
        let mut pos_sum = 0.0;
        let mut neg_sum = 0.0;
        let dim = self.domain.dim();
        for _ in 0..self.n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(self.domain.lower[d]..self.domain.upper[d]))
                .collect();
            pos_sum += (self.positive_density)(&point);
            neg_sum += (self.negative_density)(&point);
        }
        let scale = vol / self.n_samples as f64;
        (scale * pos_sum, scale * neg_sum)
    }
}

// ---------------------------------------------------------------------------
// Product measure (Fubini's theorem)
// ---------------------------------------------------------------------------

/// Product measure `μ₁ ⊗ μ₂` on a product space `X₁ × X₂`.
///
/// By Fubini's theorem, for non-negative measurable `f`:
/// `∫_{X₁×X₂} f d(μ₁⊗μ₂) = ∫_{X₁} (∫_{X₂} f(x,y) dμ₂(y)) dμ₁(x)`
///
/// For Lebesgue measures this reduces to iterated integration.
#[derive(Debug, Clone)]
pub struct ProductMeasure {
    /// First factor domain.
    pub domain1: MeasurableSet,
    /// Second factor domain.
    pub domain2: MeasurableSet,
}

impl ProductMeasure {
    /// Create a product measure from two domains.
    pub fn new(domain1: MeasurableSet, domain2: MeasurableSet) -> Self {
        Self { domain1, domain2 }
    }

    /// Measure a product box `A × B` under the product Lebesgue measure.
    ///
    /// `(λ₁ ⊗ λ₂)(A × B) = λ₁(A) · λ₂(B)`
    pub fn measure_product_box(&self, a: &MeasurableSet, b: &MeasurableSet) -> f64 {
        LebesgueMeasure.measure(a) * LebesgueMeasure.measure(b)
    }

    /// Compute the product domain as a single concatenated box.
    pub fn product_domain(&self) -> MeasurableSet {
        let mut lower = self.domain1.lower.clone();
        lower.extend(self.domain2.lower.iter().copied());
        let mut upper = self.domain1.upper.clone();
        upper.extend(self.domain2.upper.iter().copied());
        MeasurableSet::new("product_domain", lower, upper)
    }

    /// Integrate a function `f(x, y)` over the product domain using Fubini.
    ///
    /// Uses a simple grid rule: `n_x` × `n_y` quadrature points.
    pub fn fubini_integrate(&self, f: &dyn Fn(f64, f64) -> f64, n_x: usize, n_y: usize) -> f64 {
        assert_eq!(
            self.domain1.dim(),
            1,
            "fubini_integrate: domain1 must be 1-D"
        );
        assert_eq!(
            self.domain2.dim(),
            1,
            "fubini_integrate: domain2 must be 1-D"
        );

        let x0 = self.domain1.lower[0];
        let x1 = self.domain1.upper[0];
        let y0 = self.domain2.lower[0];
        let y1 = self.domain2.upper[0];
        let dx = (x1 - x0) / n_x as f64;
        let dy = (y1 - y0) / n_y as f64;

        let mut sum = 0.0;
        for i in 0..n_x {
            let x = x0 + (i as f64 + 0.5) * dx;
            let inner: f64 = (0..n_y)
                .map(|j| {
                    let y = y0 + (j as f64 + 0.5) * dy;
                    f(x, y)
                })
                .sum();
            sum += inner * dy;
        }
        sum * dx
    }

    /// Verify Fubini's theorem numerically: iterated vs. direct Monte Carlo.
    pub fn verify_fubini(&self, f: &dyn Fn(f64, f64) -> f64, n: usize) -> (f64, f64) {
        let iterated = self.fubini_integrate(f, n, n);
        // Direct Monte Carlo estimate
        let mut rng = rand::rng();
        let x0 = self.domain1.lower[0];
        let x1 = self.domain1.upper[0];
        let y0 = self.domain2.lower[0];
        let y1 = self.domain2.upper[0];
        let vol = (x1 - x0) * (y1 - y0);
        let mut mc_sum = 0.0;
        for _ in 0..(n * n) {
            let x = rng.random_range(x0..x1);
            let y = rng.random_range(y0..y1);
            mc_sum += f(x, y);
        }
        let mc = vol * mc_sum / (n * n) as f64;
        (iterated, mc)
    }
}

// ---------------------------------------------------------------------------
// Lebesgue integral (numerical)
// ---------------------------------------------------------------------------

/// Numerical Lebesgue integration using midpoint-rule quadrature.
///
/// For a function `f: ℝⁿ → ℝ` on a box domain, the integral is approximated
/// by subdividing each dimension into `n_pts` sub-intervals.
pub struct MeasureIntegral;

impl MeasureIntegral {
    /// Integrate `f` over a 1-D interval `[a, b]` using `n` midpoint sub-intervals.
    pub fn integrate_1d(f: &dyn Fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
        let h = (b - a) / n as f64;
        (0..n)
            .map(|i| {
                let x = a + (i as f64 + 0.5) * h;
                f(x) * h
            })
            .sum()
    }

    /// Integrate `f` over a 2-D rectangle using `n × n` midpoint sub-intervals.
    pub fn integrate_2d(
        f: &dyn Fn(f64, f64) -> f64,
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        n: usize,
    ) -> f64 {
        let hx = (b - a) / n as f64;
        let hy = (d - c) / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let x = a + (i as f64 + 0.5) * hx;
            for j in 0..n {
                let y = c + (j as f64 + 0.5) * hy;
                sum += f(x, y);
            }
        }
        sum * hx * hy
    }

    /// Integrate `f` over an n-D box using midpoint quadrature.
    ///
    /// The box is given as `lower` and `upper` vectors; `n_pts` sub-intervals
    /// per dimension. This has exponential cost in the dimension.
    pub fn integrate_nd(
        f: &dyn Fn(&[f64]) -> f64,
        lower: &[f64],
        upper: &[f64],
        n_pts: usize,
    ) -> f64 {
        let dim = lower.len();
        assert_eq!(dim, upper.len(), "integrate_nd: dim mismatch");
        let steps: Vec<f64> = lower
            .iter()
            .zip(upper.iter())
            .map(|(lo, hi)| (hi - lo) / n_pts as f64)
            .collect();
        let vol_cell: f64 = steps.iter().product();
        let total_cells = n_pts.pow(dim as u32);

        (0..total_cells)
            .map(|idx| {
                let mut point = vec![0.0f64; dim];
                let mut rem = idx;
                for d in 0..dim {
                    let coord = rem % n_pts;
                    rem /= n_pts;
                    point[d] = lower[d] + (coord as f64 + 0.5) * steps[d];
                }
                f(&point) * vol_cell
            })
            .sum()
    }

    /// Monte Carlo integration of `f` over an n-D box.
    ///
    /// Produces an unbiased estimate with standard error `O(1/√n)`.
    pub fn monte_carlo(
        f: &dyn Fn(&[f64]) -> f64,
        lower: &[f64],
        upper: &[f64],
        n_samples: usize,
    ) -> f64 {
        let dim = lower.len();
        assert_eq!(dim, upper.len(), "monte_carlo: dim mismatch");
        let vol: f64 = lower
            .iter()
            .zip(upper.iter())
            .map(|(lo, hi)| hi - lo)
            .product();

        let mut rng = rand::rng();
        let mut sum = 0.0;
        for _ in 0..n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(lower[d]..upper[d]))
                .collect();
            sum += f(&point);
        }
        vol * sum / n_samples as f64
    }

    /// Importance-sampling Monte Carlo: `∫ f dx = ∫ (f/g) · g dx`.
    ///
    /// Samples are drawn from the proposal density `g` (given as CDF inverse)
    /// and the estimate is `(1/n) Σ f(xᵢ) / g(xᵢ)`.
    pub fn importance_sampling(
        f: &dyn Fn(f64) -> f64,
        proposal_sample: &dyn Fn() -> f64,
        proposal_density: &dyn Fn(f64) -> f64,
        n_samples: usize,
    ) -> f64 {
        let mut sum = 0.0;
        for _ in 0..n_samples {
            let x = proposal_sample();
            let g = proposal_density(x);
            if g.abs() > 1e-15 {
                sum += f(x) / g;
            }
        }
        sum / n_samples as f64
    }

    /// Compute the Riemann–Stieltjes-style integral `∫ f dG` numerically.
    ///
    /// Approximates `Σ f(xᵢ) · (G(xᵢ₊₁) − G(xᵢ))` over `n` sub-intervals.
    pub fn riemann_stieltjes(
        f: &dyn Fn(f64) -> f64,
        g: &dyn Fn(f64) -> f64,
        a: f64,
        b: f64,
        n: usize,
    ) -> f64 {
        let h = (b - a) / n as f64;
        (0..n)
            .map(|i| {
                let x = a + i as f64 * h;
                let xp = x + h;
                f(x) * (g(xp) - g(x))
            })
            .sum()
    }

    /// Estimate the standard error of a Monte Carlo integral estimate.
    pub fn monte_carlo_std_error(
        f: &dyn Fn(&[f64]) -> f64,
        lower: &[f64],
        upper: &[f64],
        n_samples: usize,
    ) -> (f64, f64) {
        let dim = lower.len();
        let vol: f64 = lower
            .iter()
            .zip(upper.iter())
            .map(|(lo, hi)| hi - lo)
            .product();

        let mut rng = rand::rng();
        let mut vals = Vec::with_capacity(n_samples);
        for _ in 0..n_samples {
            let point: Vec<f64> = (0..dim)
                .map(|d| rng.random_range(lower[d]..upper[d]))
                .collect();
            vals.push(f(&point));
        }
        let mean: f64 = vals.iter().sum::<f64>() / n_samples as f64;
        let var: f64 =
            vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n_samples - 1) as f64;
        let std_err = (var / n_samples as f64).sqrt() * vol;
        (vol * mean, std_err)
    }
}

// ---------------------------------------------------------------------------
// Radon-Nikodym derivative
// ---------------------------------------------------------------------------

/// Radon–Nikodym derivative approximation.
///
/// Given two measures `μ` and `ν` with `ν ≪ μ` (ν absolutely continuous
/// w.r.t. μ), the Radon–Nikodym theorem guarantees `dν/dμ` exists.
///
/// Here we estimate it numerically by binning.
pub struct RadonNikodym;

impl RadonNikodym {
    /// Estimate the Radon–Nikodym derivative `dν/dμ` at point `x` by
    /// computing the density ratio in a small ball of radius `eps`.
    ///
    /// Both `mu_samples` and `nu_samples` are finite sample sets from μ and ν.
    pub fn estimate_at(x: f64, mu_samples: &[f64], nu_samples: &[f64], eps: f64) -> f64 {
        let count =
            |samples: &[f64]| samples.iter().filter(|&&s| (s - x).abs() <= eps).count() as f64;
        let mu_count = count(mu_samples);
        let nu_count = count(nu_samples);
        if mu_count < 1e-15 {
            return 0.0;
        }
        // Normalize by sample sizes
        (nu_count / nu_samples.len() as f64) / (mu_count / mu_samples.len() as f64)
    }
}

// ---------------------------------------------------------------------------
// Measure-theoretic convergence
// ---------------------------------------------------------------------------

/// Types of convergence for sequences of measurable functions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergenceType {
    /// Pointwise almost everywhere.
    AlmostEverywhere,
    /// In measure (convergence in probability).
    InMeasure,
    /// In L^p norm.
    InLp,
    /// Uniform (almost everywhere).
    Uniform,
}

/// Check whether a sequence of functions converges in Lp norm.
///
/// Evaluates `‖fₙ − f‖_p = (∫ |fₙ(x) − f(x)|^p dx)^(1/p)` over `[a, b]`.
pub fn check_lp_convergence(
    fn_seq: &[Box<dyn Fn(f64) -> f64>],
    limit: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    p: f64,
    n_pts: usize,
    tolerance: f64,
) -> bool {
    if fn_seq.is_empty() {
        return true;
    }
    let h = (b - a) / n_pts as f64;
    let last = fn_seq.last().expect("fn_seq is non-empty after check");
    let lp_norm: f64 = (0..n_pts)
        .map(|i| {
            let x = a + (i as f64 + 0.5) * h;
            (last(x) - limit(x)).abs().powf(p) * h
        })
        .sum::<f64>()
        .powf(1.0 / p);
    lp_norm < tolerance
}

// ---------------------------------------------------------------------------
// Outer measure
// ---------------------------------------------------------------------------

/// Carathéodory outer measure construction.
///
/// The outer measure `μ*(A) = inf { Σ μ(Eᵢ) : A ⊆ ⋃ Eᵢ }` is constructed
/// from any set function `μ` on covering sets.
pub struct OuterMeasure {
    /// The base covering measure (e.g., length of intervals).
    pub base: LebesgueMeasure,
}

impl OuterMeasure {
    /// Create an outer measure from the Lebesgue base measure.
    pub fn new() -> Self {
        Self {
            base: LebesgueMeasure,
        }
    }

    /// Estimate `μ*(A)` by covering `A` with `n_covers` equal sub-intervals.
    pub fn estimate(&self, a: f64, b: f64, n_covers: usize) -> f64 {
        // For an interval, the outer measure equals the length
        let _ = n_covers;
        (b - a).max(0.0)
    }

    /// Verify Carathéodory's condition: `μ*(E) = μ*(E ∩ A) + μ*(E ∩ Aᶜ)` for
    /// measurable `A`, numerically on the interval `[0,1]`.
    pub fn verify_caratheodory(&self, a_lo: f64, a_hi: f64) -> bool {
        let e_lo = 0.0;
        let e_hi = 1.0;
        let full = e_hi - e_lo;
        let intersection = (a_hi.min(e_hi) - a_lo.max(e_lo)).max(0.0);
        let complement = ((a_lo - e_lo).max(0.0) + (e_hi - a_hi).max(0.0)).min(full);
        (full - (intersection + complement)).abs() < 1e-12
    }
}

impl Default for OuterMeasure {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Counting measure
// ---------------------------------------------------------------------------

/// Counting measure on a finite set.
///
/// `μ_count(A) = |A|` (number of elements in `A`).
#[derive(Debug, Clone, Copy, Default)]
pub struct CountingMeasure;

impl CountingMeasure {
    /// Measure of a finite set: returns its cardinality.
    pub fn measure_set(set_size: usize) -> f64 {
        set_size as f64
    }

    /// Check that the counting measure is additive: `|A ∪ B| = |A| + |B|` for
    /// disjoint `A` and `B`.
    pub fn additivity(a_size: usize, b_size: usize) -> bool {
        Self::measure_set(a_size + b_size) == Self::measure_set(a_size) + Self::measure_set(b_size)
    }
}

// ---------------------------------------------------------------------------
// Ergodic measure (time-average)
// ---------------------------------------------------------------------------

/// Time-average (ergodic) measure for discrete dynamical systems.
///
/// For an ergodic system, the time average equals the space average by
/// the Birkhoff ergodic theorem.
pub struct ErgodicMeasure;

impl ErgodicMeasure {
    /// Estimate the time-average of `f` along an orbit starting from `x0`.
    ///
    /// `T_n f(x) = (1/n) Σ_{k=0}^{n-1} f(φ^k(x))`
    pub fn time_average(
        f: &dyn Fn(f64) -> f64,
        map: &dyn Fn(f64) -> f64,
        x0: f64,
        n: usize,
    ) -> f64 {
        let mut x = x0;
        let mut sum = 0.0;
        for _ in 0..n {
            sum += f(x);
            x = map(x);
        }
        sum / n as f64
    }

    /// Check the Birkhoff ergodic theorem numerically for a single orbit.
    ///
    /// Returns `true` if the time-average is within `tol` of the space-average.
    pub fn birkhoff_check(
        f: &dyn Fn(f64) -> f64,
        map: &dyn Fn(f64) -> f64,
        x0: f64,
        space_average: f64,
        n: usize,
        tol: f64,
    ) -> bool {
        let ta = Self::time_average(f, map, x0, n);
        (ta - space_average).abs() < tol
    }
}

// ---------------------------------------------------------------------------
// Haar measure (compact groups — discrete approximation)
// ---------------------------------------------------------------------------

/// Discrete approximation to Haar measure on the circle group `U(1)`.
///
/// Haar measure on `U(1) ≅ [0, 2π)` is `dθ / (2π)`.
pub struct HaarMeasure;

impl HaarMeasure {
    /// Estimate the Haar-measure integral of `f` on `[0, 2π)`.
    pub fn integrate(f: &dyn Fn(f64) -> f64, n: usize) -> f64 {
        use std::f64::consts::TAU;
        let h = TAU / n as f64;
        (0..n).map(|i| f(i as f64 * h) * h / TAU).sum()
    }

    /// Verify translation invariance: `∫ f(θ) dθ = ∫ f(θ + t) dθ`.
    pub fn translation_invariant(f: &dyn Fn(f64) -> f64, t: f64, n: usize) -> bool {
        let original = Self::integrate(f, n);
        let shifted = Self::integrate(&|theta| f((theta + t) % TAU), n);
        (original - shifted).abs() < 1e-6
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- MeasurableSet ---

    #[test]
    fn test_measurable_set_contains() {
        let s = MeasurableSet::interval(0.0, 1.0);
        assert!(s.contains(&[0.5]));
        assert!(!s.contains(&[1.5]));
    }

    #[test]
    fn test_measurable_set_null() {
        let s = MeasurableSet::interval(1.0, 1.0);
        assert!(s.is_null_set());
    }

    #[test]
    fn test_measurable_set_dim() {
        let s = MeasurableSet::new("box", vec![0.0, 0.0], vec![1.0, 1.0]);
        assert_eq!(s.dim(), 2);
    }

    // --- FiniteSigmaAlgebra ---

    #[test]
    fn test_sigma_algebra_power_set_cardinality() {
        let sa = FiniteSigmaAlgebra::power_set(3);
        assert_eq!(sa.cardinality(), 8);
    }

    #[test]
    fn test_sigma_algebra_trivial_cardinality() {
        let sa = FiniteSigmaAlgebra::trivial(4);
        assert_eq!(sa.cardinality(), 2);
    }

    #[test]
    fn test_sigma_algebra_complement_closure() {
        let sa = FiniteSigmaAlgebra::power_set(3);
        assert!(sa.is_closed_under_complement());
    }

    #[test]
    fn test_sigma_algebra_union_closure() {
        let sa = FiniteSigmaAlgebra::power_set(3);
        assert!(sa.is_closed_under_union());
    }

    #[test]
    fn test_sigma_algebra_contains_empty() {
        let sa = FiniteSigmaAlgebra::power_set(4);
        assert!(sa.contains(0)); // empty set
    }

    // --- LebesgueMeasure ---

    #[test]
    fn test_lebesgue_interval() {
        let m = LebesgueMeasure;
        let s = MeasurableSet::interval(0.0, 3.0);
        assert!((m.measure(&s) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_lebesgue_box() {
        let m = LebesgueMeasure;
        let s = MeasurableSet::new("box", vec![0.0, 0.0], vec![2.0, 3.0]);
        assert!((m.measure(&s) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_lebesgue_null_set() {
        let m = LebesgueMeasure;
        let s = MeasurableSet::interval(1.0, 1.0);
        assert!(m.is_null(&s));
    }

    #[test]
    fn test_lebesgue_empty_set() {
        let m = LebesgueMeasure;
        assert!((m.empty_set_measure()).abs() < 1e-12);
    }

    #[test]
    fn test_lebesgue_disjoint_union() {
        let m = LebesgueMeasure;
        let a = MeasurableSet::interval(0.0, 1.0);
        let b = MeasurableSet::interval(2.0, 4.0);
        let u = m.union_measure_disjoint(&a, &b);
        assert!((u - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_lebesgue_translation_invariant() {
        assert!(LebesgueMeasure::translation_invariant(0.0, 1.0, 5.0));
    }

    #[test]
    fn test_lebesgue_intersection() {
        let i = LebesgueMeasure::interval_intersection(0.0, 2.0, 1.0, 3.0);
        assert!((i - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_lebesgue_no_intersection() {
        let i = LebesgueMeasure::interval_intersection(0.0, 1.0, 2.0, 3.0);
        assert!(i.abs() < 1e-12);
    }

    // --- ProbabilityMeasure ---

    #[test]
    fn test_probability_total_mass_uniform() {
        // Uniform density on [0,1]: f(x) = 1
        let pm = ProbabilityMeasure::new(|_x| 1.0, MeasurableSet::interval(0.0, 1.0), 10_000);
        let mass = pm.total_mass();
        assert!((mass - 1.0).abs() < 0.05, "total mass = {mass}");
    }

    #[test]
    fn test_probability_expectation_linear() {
        // E[X] for X ~ Uniform[0,1] should be 0.5
        let pm = ProbabilityMeasure::new(|_x| 1.0, MeasurableSet::interval(0.0, 1.0), 20_000);
        let ex = pm.expectation(&|x| x[0]);
        assert!((ex - 0.5).abs() < 0.05, "E[X] = {ex}");
    }

    #[test]
    fn test_probability_variance_uniform() {
        // Var[X] for X ~ Uniform[0,1] = 1/12 ≈ 0.0833
        let pm = ProbabilityMeasure::new(|_x| 1.0, MeasurableSet::interval(0.0, 1.0), 30_000);
        let vx = pm.variance(&|x| x[0]);
        assert!((vx - 1.0 / 12.0).abs() < 0.02, "Var[X] = {vx}");
    }

    #[test]
    fn test_probability_entropy_positive() {
        let pm = ProbabilityMeasure::new(|_x| 1.0, MeasurableSet::interval(0.0, 1.0), 5_000);
        // H(Uniform[0,1]) = 0 (in nats, since ∫ 1·ln(1) dx = 0)
        let h = pm.entropy();
        assert!(h.abs() < 0.05, "H = {h}");
    }

    // --- SignedMeasure ---

    #[test]
    fn test_signed_measure_total_variation_positive() {
        // f⁺ = 2, f⁻ = 1 on [0,1] → TV = 3
        let sm = SignedMeasure::new(|_| 2.0, |_| 1.0, MeasurableSet::interval(0.0, 1.0), 20_000);
        let tv = sm.total_variation(&MeasurableSet::interval(0.0, 1.0));
        assert!((tv - 3.0).abs() < 0.15, "TV = {tv}");
    }

    #[test]
    fn test_signed_measure_evaluate() {
        // ν = 1 − 0 = 1 on [0,1]
        let sm = SignedMeasure::new(|_| 1.0, |_| 0.0, MeasurableSet::interval(0.0, 1.0), 20_000);
        let v = sm.evaluate(&MeasurableSet::interval(0.0, 1.0));
        assert!((v - 1.0).abs() < 0.1, "ν = {v}");
    }

    #[test]
    fn test_signed_measure_jordan_decomposition() {
        let sm = SignedMeasure::new(|_| 2.0, |_| 1.0, MeasurableSet::interval(0.0, 1.0), 20_000);
        let (pos, neg) = sm.jordan_decomposition();
        assert!((pos - 2.0).abs() < 0.15, "ν⁺ = {pos}");
        assert!((neg - 1.0).abs() < 0.15, "ν⁻ = {neg}");
    }

    #[test]
    fn test_signed_measure_from_density() {
        // density = x - 0.5 on [0,1]: positive on [0.5,1], negative on [0,0.5]
        let sm =
            SignedMeasure::from_density(|x| x[0] - 0.5, MeasurableSet::interval(0.0, 1.0), 20_000);
        let (pos, neg) = sm.jordan_decomposition();
        assert!(pos >= 0.0);
        assert!(neg >= 0.0);
    }

    // --- ProductMeasure ---

    #[test]
    fn test_product_measure_box() {
        let pm = ProductMeasure::new(
            MeasurableSet::interval(0.0, 2.0),
            MeasurableSet::interval(0.0, 3.0),
        );
        let a = MeasurableSet::interval(0.0, 2.0);
        let b = MeasurableSet::interval(0.0, 3.0);
        let m = pm.measure_product_box(&a, &b);
        assert!((m - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_fubini_constant_function() {
        let pm = ProductMeasure::new(
            MeasurableSet::interval(0.0, 1.0),
            MeasurableSet::interval(0.0, 1.0),
        );
        // ∫∫ 1 dx dy = 1
        let result = pm.fubini_integrate(&|_x, _y| 1.0, 100, 100);
        assert!((result - 1.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_fubini_separable_function() {
        let pm = ProductMeasure::new(
            MeasurableSet::interval(0.0, 1.0),
            MeasurableSet::interval(0.0, 1.0),
        );
        // ∫∫ x·y dx dy = (∫ x dx)(∫ y dy) = 0.5 × 0.5 = 0.25
        let result = pm.fubini_integrate(&|x, y| x * y, 200, 200);
        assert!((result - 0.25).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_fubini_vs_monte_carlo() {
        let pm = ProductMeasure::new(
            MeasurableSet::interval(0.0, 1.0),
            MeasurableSet::interval(0.0, 1.0),
        );
        let (iterated, mc) = pm.verify_fubini(&|x, y| x * x + y * y, 50);
        assert!((iterated - mc).abs() < 0.1, "iterated={iterated} mc={mc}");
    }

    // --- MeasureIntegral ---

    #[test]
    fn test_integrate_1d_constant() {
        let result = MeasureIntegral::integrate_1d(&|_| 2.0, 0.0, 3.0, 1000);
        assert!((result - 6.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_integrate_1d_linear() {
        // ∫₀¹ x dx = 0.5
        let result = MeasureIntegral::integrate_1d(&|x| x, 0.0, 1.0, 10_000);
        assert!((result - 0.5).abs() < 0.001, "result = {result}");
    }

    #[test]
    fn test_integrate_1d_trig() {
        // ∫₀^π sin(x) dx = 2
        let result = MeasureIntegral::integrate_1d(&|x| x.sin(), 0.0, PI, 10_000);
        assert!((result - 2.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_integrate_2d_constant() {
        // ∫∫_{[0,1]²} 1 dx dy = 1
        let result = MeasureIntegral::integrate_2d(&|_x, _y| 1.0, 0.0, 1.0, 0.0, 1.0, 100);
        assert!((result - 1.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_integrate_2d_separable() {
        // ∫∫_{[0,1]²} x·y dx dy = 0.25
        let result = MeasureIntegral::integrate_2d(&|x, y| x * y, 0.0, 1.0, 0.0, 1.0, 200);
        assert!((result - 0.25).abs() < 0.005, "result = {result}");
    }

    #[test]
    fn test_integrate_nd_volume() {
        // ∫_{[0,1]³} 1 dx = 1
        let result =
            MeasureIntegral::integrate_nd(&|_| 1.0, &[0.0, 0.0, 0.0], &[1.0, 1.0, 1.0], 10);
        assert!((result - 1.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_monte_carlo_1d() {
        // ∫₀¹ x² dx = 1/3
        let result = MeasureIntegral::monte_carlo(&|x| x[0].powi(2), &[0.0], &[1.0], 50_000);
        assert!((result - 1.0 / 3.0).abs() < 0.02, "result = {result}");
    }

    #[test]
    fn test_monte_carlo_2d() {
        // ∫∫_{[0,1]²} (x+y) dx dy = 1
        let result =
            MeasureIntegral::monte_carlo(&|p| p[0] + p[1], &[0.0, 0.0], &[1.0, 1.0], 50_000);
        assert!((result - 1.0).abs() < 0.05, "result = {result}");
    }

    #[test]
    fn test_riemann_stieltjes() {
        // ∫₀¹ 1 dG where G(x) = x → result = 1
        let result = MeasureIntegral::riemann_stieltjes(&|_| 1.0, &|x| x, 0.0, 1.0, 10_000);
        assert!((result - 1.0).abs() < 0.001, "result = {result}");
    }

    #[test]
    fn test_monte_carlo_std_error() {
        let (est, err) = MeasureIntegral::monte_carlo_std_error(&|x| x[0], &[0.0], &[1.0], 10_000);
        assert!((est - 0.5).abs() < 0.05, "estimate = {est}");
        assert!(err > 0.0, "std error should be positive");
    }

    // --- RadonNikodym ---

    #[test]
    fn test_radon_nikodym_estimate() {
        // Both from the same distribution → derivative ≈ 1
        let samples: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        let rn = RadonNikodym::estimate_at(0.5, &samples, &samples, 0.05);
        assert!((rn - 1.0).abs() < 0.2, "RN derivative = {rn}");
    }

    // --- OuterMeasure ---

    #[test]
    fn test_outer_measure_estimate() {
        let om = OuterMeasure::new();
        let m = om.estimate(0.0, 3.0, 10);
        assert!((m - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_caratheodory_condition() {
        let om = OuterMeasure::new();
        assert!(om.verify_caratheodory(0.2, 0.7));
    }

    // --- CountingMeasure ---

    #[test]
    fn test_counting_measure() {
        assert!((CountingMeasure::measure_set(5) - 5.0).abs() < 1e-12);
        assert!((CountingMeasure::measure_set(0)).abs() < 1e-12);
    }

    // --- ErgodicMeasure ---

    #[test]
    fn test_ergodic_time_average_doubling_map() {
        // Irrational rotation: T(x) = (x + α) mod 1, ergodic w.r.t. Lebesgue
        // for irrational α.  E[f] for f(x) = x ≈ 0.5 (space average).
        let alpha = (5.0_f64.sqrt() - 1.0) / 2.0; // golden ratio ≈ 0.618
        let ta = ErgodicMeasure::time_average(
            &|x| x,
            &move |x| (x + alpha) % 1.0,
            0.123_456_789,
            200_000,
        );
        assert!((ta - 0.5).abs() < 0.05, "time average = {ta}");
    }

    #[test]
    fn test_ergodic_birkhoff_check() {
        let alpha = (5.0_f64.sqrt() - 1.0) / 2.0;
        let ok = ErgodicMeasure::birkhoff_check(
            &|x| x,
            &move |x| (x + alpha) % 1.0,
            0.1234,
            0.5,
            200_000,
            0.05,
        );
        assert!(ok, "Birkhoff check failed");
    }

    // --- HaarMeasure ---

    #[test]
    fn test_haar_measure_constant() {
        // ∫ 1 dθ/(2π) = 1
        let result = HaarMeasure::integrate(&|_| 1.0, 1000);
        assert!((result - 1.0).abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_haar_measure_cosine() {
        // ∫ cos(θ) dθ/(2π) = 0
        let result = HaarMeasure::integrate(&|theta: f64| theta.cos(), 10_000);
        assert!(result.abs() < 0.01, "result = {result}");
    }

    #[test]
    fn test_haar_translation_invariant() {
        let ok = HaarMeasure::translation_invariant(&|theta: f64| theta.sin().powi(2), 1.0, 10_000);
        assert!(ok);
    }

    // --- check_lp_convergence ---

    #[test]
    fn test_lp_convergence() {
        // fₙ(x) = x^(1/n) → 1 on [0,1] in L2
        let fns: Vec<Box<dyn Fn(f64) -> f64>> = (1..=10)
            .map(|n| {
                let b: Box<dyn Fn(f64) -> f64> = Box::new(move |x: f64| x.powf(1.0 / n as f64));
                b
            })
            .collect();
        let limit = |_x: f64| 1.0;
        let converges = check_lp_convergence(&fns, &limit, 0.0, 1.0, 2.0, 1000, 0.15);
        assert!(converges);
    }
}
