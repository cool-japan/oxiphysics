// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-collision.
//!
//! This module defines the error hierarchy used throughout the collision
//! detection subsystem, including classification helpers, recovery hints,
//! and a severity scale for runtime triage.

use thiserror::Error;

// ── Error ─────────────────────────────────────────────────────────────────────

/// Main error type for the collision module.
#[derive(Debug, Error)]
pub enum Error {
    /// Generic unclassified error message.
    #[error("{0}")]
    General(String),

    /// A GJK query did not converge within the allowed iteration budget.
    #[error("GJK did not converge after {iterations} iterations (dist={dist:.6})")]
    GjkNoConverge {
        /// Number of iterations attempted before giving up.
        iterations: usize,
        /// Remaining distance at termination.
        dist: f64,
    },

    /// The EPA polytope expansion failed to find a valid separating face.
    #[error("EPA failed: {reason}")]
    EpaFailure {
        /// Human-readable reason.
        reason: String,
    },

    /// Conservative-advancement CCD did not converge within the allowed budget.
    #[error("CCD did not converge: toi={toi:.6} after {iters} iters")]
    CcdNoConverge {
        /// Best TOI estimate when iterations were exhausted.
        toi: f64,
        /// Number of iterations attempted.
        iters: usize,
    },

    /// A shape is degenerate (zero volume, zero radius, etc.).
    #[error("degenerate shape: {detail}")]
    DegenerateShape {
        /// Description of the degeneracy.
        detail: String,
    },

    /// An index (body, shape, node …) was out of range.
    #[error("index {index} out of range (len={len})")]
    IndexOutOfRange {
        /// The offending index.
        index: usize,
        /// The valid range length.
        len: usize,
    },

    /// The broadphase structure has exceeded its capacity.
    #[error("broadphase capacity exceeded: {capacity} slots are full")]
    BroadphaseCapacity {
        /// Maximum capacity that was exceeded.
        capacity: usize,
    },

    /// A numerical computation produced a non-finite value (NaN / Inf).
    #[error("non-finite value encountered in {context}")]
    NonFinite {
        /// Where the non-finite was detected.
        context: String,
    },

    /// The requested feature is not yet implemented for this shape combination.
    #[error("unsupported shape pair: {shape_a} vs {shape_b}")]
    UnsupportedShapePair {
        /// Name of shape A.
        shape_a: String,
        /// Name of shape B.
        shape_b: String,
    },
}

/// Result type alias.
pub type Result<T> = std::result::Result<T, Error>;

// ── ErrorKind ────────────────────────────────────────────────────────────────

/// Coarse classification of [`enum@Error`] variants for programmatic branching.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ErrorKind {
    /// Algorithm failed to converge.
    Convergence,
    /// Numerical / degenerate geometry.
    Numerical,
    /// Out-of-range index.
    Bounds,
    /// Capacity / resource limit exceeded.
    Capacity,
    /// Feature not available for this shape combination.
    Unsupported,
    /// Any other error.
    Other,
}

impl Error {
    /// Return the coarse [`ErrorKind`] for this error.
    pub fn kind(&self) -> ErrorKind {
        match self {
            Error::GjkNoConverge { .. } | Error::CcdNoConverge { .. } => ErrorKind::Convergence,
            Error::EpaFailure { .. } | Error::DegenerateShape { .. } | Error::NonFinite { .. } => {
                ErrorKind::Numerical
            }
            Error::IndexOutOfRange { .. } => ErrorKind::Bounds,
            Error::BroadphaseCapacity { .. } => ErrorKind::Capacity,
            Error::UnsupportedShapePair { .. } => ErrorKind::Unsupported,
            Error::General(_) => ErrorKind::Other,
        }
    }

    /// Returns `true` when this error indicates a convergence failure.
    pub fn is_convergence(&self) -> bool {
        self.kind() == ErrorKind::Convergence
    }

    /// Returns `true` when this error indicates a numerical problem.
    pub fn is_numerical(&self) -> bool {
        self.kind() == ErrorKind::Numerical
    }

    /// Returns `true` when this error indicates an out-of-bounds access.
    pub fn is_bounds(&self) -> bool {
        self.kind() == ErrorKind::Bounds
    }
}

// ── Severity ─────────────────────────────────────────────────────────────────

/// Severity level for runtime triage and logging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Severity {
    /// Informational — simulation continues normally.
    Info = 0,
    /// Warning — result may be approximate.
    Warning = 1,
    /// Error — a contact or query was skipped.
    Error = 2,
    /// Fatal — the simulation state is likely corrupt.
    Fatal = 3,
}

impl Error {
    /// Suggest a [`Severity`] level for this error.
    pub fn severity(&self) -> Severity {
        match self {
            Error::GjkNoConverge { .. } | Error::CcdNoConverge { .. } => Severity::Warning,
            Error::EpaFailure { .. } => Severity::Warning,
            Error::DegenerateShape { .. } => Severity::Warning,
            Error::NonFinite { .. } => Severity::Fatal,
            Error::IndexOutOfRange { .. } => Severity::Error,
            Error::BroadphaseCapacity { .. } => Severity::Error,
            Error::UnsupportedShapePair { .. } => Severity::Warning,
            Error::General(_) => Severity::Info,
        }
    }

    /// Returns `true` if the severity is at least [`Severity::Error`].
    pub fn is_severe(&self) -> bool {
        self.severity() >= Severity::Error
    }
}

// ── RecoveryHint ─────────────────────────────────────────────────────────────

/// Suggested recovery action for a failed collision query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecoveryHint {
    /// Skip this pair this frame; it will be retried next frame.
    SkipPair,
    /// Treat the bodies as separated (no contact generated).
    TreatAsSeparated,
    /// Fall back to a simpler, less accurate algorithm.
    UseFallback,
    /// Halt the simulation — the state is likely corrupt.
    Abort,
}

impl Error {
    /// Suggest what the caller should do when encountering this error.
    pub fn recovery_hint(&self) -> RecoveryHint {
        match self {
            Error::GjkNoConverge { .. }
            | Error::EpaFailure { .. }
            | Error::CcdNoConverge { .. } => RecoveryHint::TreatAsSeparated,
            Error::DegenerateShape { .. } | Error::UnsupportedShapePair { .. } => {
                RecoveryHint::SkipPair
            }
            Error::NonFinite { .. } => RecoveryHint::Abort,
            Error::IndexOutOfRange { .. } | Error::BroadphaseCapacity { .. } => {
                RecoveryHint::UseFallback
            }
            Error::General(_) => RecoveryHint::SkipPair,
        }
    }
}

// ── ErrorCollector ────────────────────────────────────────────────────────────

/// Collects multiple non-fatal errors so the caller can inspect them in bulk
/// after a simulation step, rather than short-circuiting on the first failure.
#[derive(Debug, Default)]
pub struct ErrorCollector {
    errors: Vec<Error>,
}

impl ErrorCollector {
    /// Create an empty collector.
    pub fn new() -> Self {
        Self { errors: Vec::new() }
    }

    /// Push an error into the collector.
    ///
    /// Returns `Err(error)` if the error is [`Severity::Fatal`] (caller must abort).
    /// Otherwise stores it and returns `Ok(())`.
    pub fn push(&mut self, err: Error) -> Result<()> {
        if err.severity() == Severity::Fatal {
            return Err(err);
        }
        self.errors.push(err);
        Ok(())
    }

    /// Returns `true` if any errors were collected.
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Number of collected errors.
    pub fn len(&self) -> usize {
        self.errors.len()
    }

    /// Returns `true` when no errors have been collected.
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }

    /// Drain and return all collected errors.
    pub fn take_errors(&mut self) -> Vec<Error> {
        std::mem::take(&mut self.errors)
    }

    /// Count how many errors belong to a given [`ErrorKind`].
    pub fn count_by_kind(&self, kind: ErrorKind) -> usize {
        self.errors.iter().filter(|e| e.kind() == kind).count()
    }

    /// Returns all errors of at least the given severity.
    pub fn errors_above(&self, min: Severity) -> Vec<&Error> {
        self.errors.iter().filter(|e| e.severity() >= min).collect()
    }
}

// ── ErrorBudget ───────────────────────────────────────────────────────────────

/// Tracks per-category error budgets for rate-limited triage.
///
/// A budget tracks how many errors of each [`ErrorKind`] are tolerated before
/// an action (e.g., logging) should be triggered.  Budgets are decremented on
/// each call to [`ErrorBudget::check`]; when they reach zero the call returns
/// `true` indicating the caller should act.
#[derive(Debug, Clone)]
pub struct ErrorBudget {
    convergence: u32,
    numerical: u32,
    bounds: u32,
    capacity: u32,
    unsupported: u32,
    other: u32,
}

impl Default for ErrorBudget {
    fn default() -> Self {
        Self {
            convergence: 16,
            numerical: 8,
            bounds: 4,
            capacity: 2,
            unsupported: 32,
            other: 64,
        }
    }
}

impl ErrorBudget {
    /// Create a new budget with default thresholds.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a budget where every category has the same threshold.
    pub fn uniform(threshold: u32) -> Self {
        Self {
            convergence: threshold,
            numerical: threshold,
            bounds: threshold,
            capacity: threshold,
            unsupported: threshold,
            other: threshold,
        }
    }

    /// Decrement the budget for the given error kind.
    ///
    /// Returns `true` if the budget was already zero (caller should act),
    /// or if it just reached zero this call.
    pub fn check(&mut self, kind: ErrorKind) -> bool {
        let slot = match kind {
            ErrorKind::Convergence => &mut self.convergence,
            ErrorKind::Numerical => &mut self.numerical,
            ErrorKind::Bounds => &mut self.bounds,
            ErrorKind::Capacity => &mut self.capacity,
            ErrorKind::Unsupported => &mut self.unsupported,
            ErrorKind::Other => &mut self.other,
        };
        if *slot == 0 {
            return true;
        }
        *slot -= 1;
        *slot == 0
    }

    /// Reset all budgets to their default values.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Returns `true` if the budget for this kind is exhausted (== 0).
    pub fn is_exhausted(&self, kind: ErrorKind) -> bool {
        match kind {
            ErrorKind::Convergence => self.convergence == 0,
            ErrorKind::Numerical => self.numerical == 0,
            ErrorKind::Bounds => self.bounds == 0,
            ErrorKind::Capacity => self.capacity == 0,
            ErrorKind::Unsupported => self.unsupported == 0,
            ErrorKind::Other => self.other == 0,
        }
    }
}

// ── ErrorStats ────────────────────────────────────────────────────────────────

/// Accumulates error counts per kind and severity for a simulation step.
///
/// Useful for telemetry: after each physics step the engine can drain the
/// collector, tally the stats here, then emit a single structured log line
/// instead of one log per error.
#[derive(Debug, Default, Clone)]
pub struct ErrorStats {
    /// Counts per `ErrorKind`.
    pub by_kind: [u64; 6],
    /// Counts per `Severity`.
    pub by_severity: [u64; 4],
    /// Total error count.
    pub total: u64,
}

impl ErrorStats {
    /// Create a zeroed stats block.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a single error observation.
    pub fn record(&mut self, err: &Error) {
        self.total += 1;
        let ki = match err.kind() {
            ErrorKind::Convergence => 0,
            ErrorKind::Numerical => 1,
            ErrorKind::Bounds => 2,
            ErrorKind::Capacity => 3,
            ErrorKind::Unsupported => 4,
            ErrorKind::Other => 5,
        };
        self.by_kind[ki] += 1;
        let si = match err.severity() {
            Severity::Info => 0,
            Severity::Warning => 1,
            Severity::Error => 2,
            Severity::Fatal => 3,
        };
        self.by_severity[si] += 1;
    }

    /// Absorb all errors from a [`ErrorCollector`], recording each one.
    pub fn absorb(&mut self, collector: &mut ErrorCollector) {
        for err in collector.take_errors() {
            self.record(&err);
        }
    }

    /// Returns the count for a specific [`ErrorKind`].
    pub fn count_kind(&self, kind: ErrorKind) -> u64 {
        let ki = match kind {
            ErrorKind::Convergence => 0,
            ErrorKind::Numerical => 1,
            ErrorKind::Bounds => 2,
            ErrorKind::Capacity => 3,
            ErrorKind::Unsupported => 4,
            ErrorKind::Other => 5,
        };
        self.by_kind[ki]
    }

    /// Returns the count for a specific [`Severity`].
    pub fn count_severity(&self, sev: Severity) -> u64 {
        let si = match sev {
            Severity::Info => 0,
            Severity::Warning => 1,
            Severity::Error => 2,
            Severity::Fatal => 3,
        };
        self.by_severity[si]
    }

    /// Reset all counters to zero.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Returns `true` if any fatal errors were recorded.
    pub fn has_fatals(&self) -> bool {
        self.by_severity[3] > 0
    }

    /// Merge another stats block into this one.
    pub fn merge(&mut self, other: &ErrorStats) {
        self.total += other.total;
        for i in 0..6 {
            self.by_kind[i] += other.by_kind[i];
        }
        for i in 0..4 {
            self.by_severity[i] += other.by_severity[i];
        }
    }
}

// ── ContextualError ───────────────────────────────────────────────────────────

/// An [`enum@Error`] decorated with the body-pair indices that triggered it.
///
/// This lets the engine log "pair (3, 7) produced a GJK convergence failure"
/// without baking body indices into every error variant.
#[derive(Debug)]
pub struct ContextualError {
    /// The underlying error.
    pub error: Error,
    /// Index of the first body involved.
    pub body_a: usize,
    /// Index of the second body involved.
    pub body_b: usize,
    /// Optional free-form context tag (algorithm name, phase, etc.).
    pub tag: Option<&'static str>,
}

impl ContextualError {
    /// Wrap an [`enum@Error`] with body indices.
    pub fn new(error: Error, body_a: usize, body_b: usize) -> Self {
        Self {
            error,
            body_a,
            body_b,
            tag: None,
        }
    }

    /// Attach an optional tag and return `self`.
    pub fn with_tag(mut self, tag: &'static str) -> Self {
        self.tag = Some(tag);
        self
    }

    /// The kind of the wrapped error.
    pub fn kind(&self) -> ErrorKind {
        self.error.kind()
    }

    /// The severity of the wrapped error.
    pub fn severity(&self) -> Severity {
        self.error.severity()
    }

    /// Recovery hint for the wrapped error.
    pub fn recovery_hint(&self) -> RecoveryHint {
        self.error.recovery_hint()
    }
}

impl std::fmt::Display for ContextualError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.tag {
            Some(tag) => write!(
                f,
                "[{}] bodies ({},{}) — {}",
                tag, self.body_a, self.body_b, self.error
            ),
            None => write!(
                f,
                "bodies ({},{}) — {}",
                self.body_a, self.body_b, self.error
            ),
        }
    }
}

// ── ContextualCollector ────────────────────────────────────────────────────────

/// Like [`ErrorCollector`] but stores [`ContextualError`] values.
#[derive(Debug, Default)]
pub struct ContextualCollector {
    errors: Vec<ContextualError>,
}

impl ContextualCollector {
    /// Create an empty collector.
    pub fn new() -> Self {
        Self::default()
    }

    /// Push a contextual error.
    ///
    /// Fatal errors are NOT stored — they are returned directly so the caller
    /// can abort.
    pub fn push(&mut self, err: ContextualError) -> std::result::Result<(), ContextualError> {
        if err.severity() == Severity::Fatal {
            return Err(err);
        }
        self.errors.push(err);
        Ok(())
    }

    /// Number of collected errors.
    pub fn len(&self) -> usize {
        self.errors.len()
    }

    /// Returns `true` when no errors have been collected.
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }

    /// Drain and return all collected errors.
    pub fn take_errors(&mut self) -> Vec<ContextualError> {
        std::mem::take(&mut self.errors)
    }

    /// Count errors involving a specific body index.
    pub fn count_for_body(&self, body: usize) -> usize {
        self.errors
            .iter()
            .filter(|e| e.body_a == body || e.body_b == body)
            .count()
    }

    /// Collect all body pairs that produced errors.
    pub fn error_pairs(&self) -> Vec<(usize, usize)> {
        self.errors.iter().map(|e| (e.body_a, e.body_b)).collect()
    }

    /// Drain into an [`ErrorStats`] block.
    pub fn drain_into_stats(&mut self, stats: &mut ErrorStats) {
        for ce in self.take_errors() {
            stats.record(&ce.error);
        }
    }
}

// ── Retry policy ──────────────────────────────────────────────────────────────

/// Specifies how many times a failed collision query should be retried.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RetryPolicy {
    /// Maximum number of retries for convergence failures.
    pub max_convergence_retries: u32,
    /// Maximum number of retries for numerical failures.
    pub max_numerical_retries: u32,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        Self {
            max_convergence_retries: 2,
            max_numerical_retries: 1,
        }
    }
}

impl RetryPolicy {
    /// Number of retries for an error of the given kind.
    pub fn retries_for(&self, kind: ErrorKind) -> u32 {
        match kind {
            ErrorKind::Convergence => self.max_convergence_retries,
            ErrorKind::Numerical => self.max_numerical_retries,
            _ => 0,
        }
    }

    /// Returns `true` if the policy allows any retries for the kind.
    pub fn should_retry(&self, kind: ErrorKind) -> bool {
        self.retries_for(kind) > 0
    }
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general_error_kind_is_other() {
        let e = Error::General("oops".into());
        assert_eq!(e.kind(), ErrorKind::Other);
    }

    #[test]
    fn gjk_no_converge_kind() {
        let e = Error::GjkNoConverge {
            iterations: 64,
            dist: 0.001,
        };
        assert_eq!(e.kind(), ErrorKind::Convergence);
        assert!(e.is_convergence());
    }

    #[test]
    fn ccd_no_converge_kind() {
        let e = Error::CcdNoConverge {
            toi: 0.5,
            iters: 32,
        };
        assert_eq!(e.kind(), ErrorKind::Convergence);
        assert!(e.is_convergence());
    }

    #[test]
    fn non_finite_is_fatal_and_severe() {
        let e = Error::NonFinite {
            context: "EPA".into(),
        };
        assert_eq!(e.severity(), Severity::Fatal);
        assert!(e.is_severe());
        assert_eq!(e.recovery_hint(), RecoveryHint::Abort);
    }

    #[test]
    fn degenerate_shape_warning() {
        let e = Error::DegenerateShape {
            detail: "radius=0".into(),
        };
        assert_eq!(e.severity(), Severity::Warning);
        assert!(!e.is_severe());
        assert!(e.is_numerical());
        assert_eq!(e.recovery_hint(), RecoveryHint::SkipPair);
    }

    #[test]
    fn index_out_of_range_is_bounds() {
        let e = Error::IndexOutOfRange { index: 99, len: 10 };
        assert_eq!(e.kind(), ErrorKind::Bounds);
        assert!(e.is_bounds());
        assert!(e.is_severe());
    }

    #[test]
    fn collector_push_non_fatal_ok() {
        let mut col = ErrorCollector::new();
        col.push(Error::General("warn".into())).unwrap();
        assert_eq!(col.len(), 1);
        assert!(col.has_errors());
    }

    #[test]
    fn collector_push_fatal_returns_err() {
        let mut col = ErrorCollector::new();
        let result = col.push(Error::NonFinite {
            context: "test".into(),
        });
        assert!(result.is_err());
        // Fatal errors are NOT stored in the collector — they are propagated
        assert!(col.is_empty());
    }

    #[test]
    fn collector_count_by_kind() {
        let mut col = ErrorCollector::new();
        col.push(Error::GjkNoConverge {
            iterations: 10,
            dist: 0.01,
        })
        .unwrap();
        col.push(Error::CcdNoConverge { toi: 0.3, iters: 5 })
            .unwrap();
        col.push(Error::General("misc".into())).unwrap();
        assert_eq!(col.count_by_kind(ErrorKind::Convergence), 2);
        assert_eq!(col.count_by_kind(ErrorKind::Other), 1);
    }

    #[test]
    fn collector_errors_above_error_severity() {
        let mut col = ErrorCollector::new();
        col.push(Error::GjkNoConverge {
            iterations: 8,
            dist: 0.0,
        })
        .unwrap();
        col.push(Error::BroadphaseCapacity { capacity: 1024 })
            .unwrap();
        let severe = col.errors_above(Severity::Error);
        assert_eq!(severe.len(), 1); // only BroadphaseCapacity is Error-level
    }

    #[test]
    fn collector_take_errors_clears() {
        let mut col = ErrorCollector::new();
        col.push(Error::General("a".into())).unwrap();
        col.push(Error::General("b".into())).unwrap();
        let taken = col.take_errors();
        assert_eq!(taken.len(), 2);
        assert!(col.is_empty());
    }

    #[test]
    fn severity_ordering() {
        assert!(Severity::Fatal > Severity::Error);
        assert!(Severity::Error > Severity::Warning);
        assert!(Severity::Warning > Severity::Info);
    }

    #[test]
    fn unsupported_shape_pair() {
        let e = Error::UnsupportedShapePair {
            shape_a: "Torus".into(),
            shape_b: "Cone".into(),
        };
        assert_eq!(e.kind(), ErrorKind::Unsupported);
        assert_eq!(e.recovery_hint(), RecoveryHint::SkipPair);
        assert!(!e.is_severe());
    }

    #[test]
    fn broadphase_capacity_error_level() {
        let e = Error::BroadphaseCapacity { capacity: 4096 };
        assert_eq!(e.kind(), ErrorKind::Capacity);
        assert_eq!(e.severity(), Severity::Error);
        assert_eq!(e.recovery_hint(), RecoveryHint::UseFallback);
    }

    #[test]
    fn epa_failure_convergence_and_warn() {
        let e = Error::EpaFailure {
            reason: "degenerate polytope".into(),
        };
        assert_eq!(e.kind(), ErrorKind::Numerical);
        assert_eq!(e.severity(), Severity::Warning);
        assert_eq!(e.recovery_hint(), RecoveryHint::TreatAsSeparated);
    }

    #[test]
    fn general_error_display() {
        let e = Error::General("something went wrong".into());
        assert_eq!(format!("{e}"), "something went wrong");
    }

    #[test]
    fn gjk_no_converge_display_contains_iterations() {
        let e = Error::GjkNoConverge {
            iterations: 64,
            dist: 0.001_234,
        };
        let s = format!("{e}");
        assert!(
            s.contains("64"),
            "display should contain iteration count: {s}"
        );
    }

    // ── ErrorBudget tests ────────────────────────────────────────────────────

    #[test]
    fn budget_uniform_starts_full() {
        let b = ErrorBudget::uniform(5);
        assert!(!b.is_exhausted(ErrorKind::Convergence));
        assert!(!b.is_exhausted(ErrorKind::Numerical));
        assert!(!b.is_exhausted(ErrorKind::Bounds));
    }

    #[test]
    fn budget_check_decrements_to_zero() {
        let mut b = ErrorBudget::uniform(2);
        let first = b.check(ErrorKind::Convergence);
        assert!(!first, "first check should not signal exhaustion yet");
        let second = b.check(ErrorKind::Convergence);
        assert!(second, "second check should signal exhaustion (reached 0)");
        assert!(b.is_exhausted(ErrorKind::Convergence));
    }

    #[test]
    fn budget_check_already_zero_returns_true() {
        let mut b = ErrorBudget::uniform(0);
        // Budget already zero — every check should return true immediately
        assert!(b.check(ErrorKind::Numerical));
        assert!(b.check(ErrorKind::Numerical));
    }

    #[test]
    fn budget_reset_restores_defaults() {
        let mut b = ErrorBudget::uniform(1);
        b.check(ErrorKind::Convergence);
        assert!(b.is_exhausted(ErrorKind::Convergence));
        b.reset();
        // After reset the budget should no longer be exhausted (default > 0)
        assert!(!b.is_exhausted(ErrorKind::Convergence));
    }

    #[test]
    fn budget_each_kind_independent() {
        let mut b = ErrorBudget::uniform(1);
        b.check(ErrorKind::Convergence);
        assert!(b.is_exhausted(ErrorKind::Convergence));
        assert!(!b.is_exhausted(ErrorKind::Numerical));
        assert!(!b.is_exhausted(ErrorKind::Bounds));
    }

    // ── ErrorStats tests ─────────────────────────────────────────────────────

    #[test]
    fn stats_starts_zeroed() {
        let s = ErrorStats::new();
        assert_eq!(s.total, 0);
        assert_eq!(s.count_kind(ErrorKind::Convergence), 0);
        assert_eq!(s.count_severity(Severity::Fatal), 0);
    }

    #[test]
    fn stats_record_convergence() {
        let mut s = ErrorStats::new();
        let e = Error::GjkNoConverge {
            iterations: 10,
            dist: 0.0,
        };
        s.record(&e);
        assert_eq!(s.total, 1);
        assert_eq!(s.count_kind(ErrorKind::Convergence), 1);
        assert_eq!(s.count_severity(Severity::Warning), 1);
        assert!(!s.has_fatals());
    }

    #[test]
    fn stats_record_fatal() {
        let mut s = ErrorStats::new();
        let e = Error::NonFinite {
            context: "test".into(),
        };
        s.record(&e);
        assert_eq!(s.total, 1);
        assert!(s.has_fatals());
        assert_eq!(s.count_severity(Severity::Fatal), 1);
    }

    #[test]
    fn stats_merge() {
        let mut a = ErrorStats::new();
        a.record(&Error::General("x".into()));
        let mut b = ErrorStats::new();
        b.record(&Error::General("y".into()));
        b.record(&Error::GjkNoConverge {
            iterations: 1,
            dist: 0.0,
        });
        a.merge(&b);
        assert_eq!(a.total, 3);
        assert_eq!(a.count_kind(ErrorKind::Other), 2);
        assert_eq!(a.count_kind(ErrorKind::Convergence), 1);
    }

    #[test]
    fn stats_absorb_from_collector() {
        let mut col = ErrorCollector::new();
        col.push(Error::GjkNoConverge {
            iterations: 5,
            dist: 0.1,
        })
        .unwrap();
        col.push(Error::CcdNoConverge { toi: 0.5, iters: 3 })
            .unwrap();
        let mut stats = ErrorStats::new();
        stats.absorb(&mut col);
        assert_eq!(stats.total, 2);
        assert_eq!(stats.count_kind(ErrorKind::Convergence), 2);
        // Collector should be empty now
        assert!(col.is_empty());
    }

    #[test]
    fn stats_reset_clears_all() {
        let mut s = ErrorStats::new();
        s.record(&Error::General("a".into()));
        s.record(&Error::GjkNoConverge {
            iterations: 1,
            dist: 0.0,
        });
        s.reset();
        assert_eq!(s.total, 0);
        assert_eq!(s.count_kind(ErrorKind::Other), 0);
        assert_eq!(s.count_kind(ErrorKind::Convergence), 0);
    }

    #[test]
    fn stats_all_kinds_counted() {
        let mut s = ErrorStats::new();
        s.record(&Error::GjkNoConverge {
            iterations: 1,
            dist: 0.0,
        }); // Convergence
        s.record(&Error::EpaFailure { reason: "x".into() }); // Numerical
        s.record(&Error::IndexOutOfRange { index: 1, len: 0 }); // Bounds
        s.record(&Error::BroadphaseCapacity { capacity: 10 }); // Capacity
        s.record(&Error::UnsupportedShapePair {
            shape_a: "A".into(),
            shape_b: "B".into(),
        }); // Unsupported
        s.record(&Error::General("misc".into())); // Other
        assert_eq!(s.total, 6);
        assert_eq!(s.count_kind(ErrorKind::Convergence), 1);
        assert_eq!(s.count_kind(ErrorKind::Numerical), 1);
        assert_eq!(s.count_kind(ErrorKind::Bounds), 1);
        assert_eq!(s.count_kind(ErrorKind::Capacity), 1);
        assert_eq!(s.count_kind(ErrorKind::Unsupported), 1);
        assert_eq!(s.count_kind(ErrorKind::Other), 1);
    }

    // ── ContextualError tests ────────────────────────────────────────────────

    #[test]
    fn contextual_error_display_includes_bodies() {
        let ce = ContextualError::new(
            Error::GjkNoConverge {
                iterations: 5,
                dist: 0.01,
            },
            3,
            7,
        );
        let s = format!("{ce}");
        assert!(s.contains('3'), "should contain body_a: {s}");
        assert!(s.contains('7'), "should contain body_b: {s}");
    }

    #[test]
    fn contextual_error_with_tag_display() {
        let ce = ContextualError::new(Error::General("fail".into()), 0, 1).with_tag("GJK");
        let s = format!("{ce}");
        assert!(s.contains("GJK"), "display should include tag: {s}");
    }

    #[test]
    fn contextual_error_delegates_kind_and_severity() {
        let ce = ContextualError::new(
            Error::NonFinite {
                context: "test".into(),
            },
            0,
            1,
        );
        assert_eq!(ce.kind(), ErrorKind::Numerical);
        assert_eq!(ce.severity(), Severity::Fatal);
        assert_eq!(ce.recovery_hint(), RecoveryHint::Abort);
    }

    // ── ContextualCollector tests ─────────────────────────────────────────────

    #[test]
    fn contextual_collector_push_and_len() {
        let mut col = ContextualCollector::new();
        col.push(ContextualError::new(Error::General("x".into()), 0, 1))
            .unwrap();
        col.push(ContextualError::new(
            Error::GjkNoConverge {
                iterations: 1,
                dist: 0.0,
            },
            2,
            3,
        ))
        .unwrap();
        assert_eq!(col.len(), 2);
        assert!(!col.is_empty());
    }

    #[test]
    fn contextual_collector_fatal_not_stored() {
        let mut col = ContextualCollector::new();
        let result = col.push(ContextualError::new(
            Error::NonFinite {
                context: "epa".into(),
            },
            0,
            1,
        ));
        assert!(result.is_err(), "fatal should not be stored");
        assert!(col.is_empty());
    }

    #[test]
    fn contextual_collector_count_for_body() {
        let mut col = ContextualCollector::new();
        col.push(ContextualError::new(Error::General("a".into()), 5, 3))
            .unwrap();
        col.push(ContextualError::new(Error::General("b".into()), 5, 7))
            .unwrap();
        col.push(ContextualError::new(Error::General("c".into()), 1, 2))
            .unwrap();
        assert_eq!(col.count_for_body(5), 2);
        assert_eq!(col.count_for_body(3), 1);
        assert_eq!(col.count_for_body(99), 0);
    }

    #[test]
    fn contextual_collector_error_pairs() {
        let mut col = ContextualCollector::new();
        col.push(ContextualError::new(Error::General("a".into()), 1, 2))
            .unwrap();
        col.push(ContextualError::new(Error::General("b".into()), 3, 4))
            .unwrap();
        let pairs = col.error_pairs();
        assert_eq!(pairs.len(), 2);
        assert!(pairs.contains(&(1, 2)));
        assert!(pairs.contains(&(3, 4)));
    }

    #[test]
    fn contextual_collector_drain_into_stats() {
        let mut col = ContextualCollector::new();
        col.push(ContextualError::new(
            Error::GjkNoConverge {
                iterations: 3,
                dist: 0.01,
            },
            0,
            1,
        ))
        .unwrap();
        col.push(ContextualError::new(Error::General("misc".into()), 2, 3))
            .unwrap();
        let mut stats = ErrorStats::new();
        col.drain_into_stats(&mut stats);
        assert_eq!(stats.total, 2);
        assert!(col.is_empty());
    }

    #[test]
    fn contextual_collector_take_clears() {
        let mut col = ContextualCollector::new();
        col.push(ContextualError::new(Error::General("x".into()), 0, 0))
            .unwrap();
        let taken = col.take_errors();
        assert_eq!(taken.len(), 1);
        assert!(col.is_empty());
    }

    // ── RetryPolicy tests ────────────────────────────────────────────────────

    #[test]
    fn retry_policy_convergence_defaults() {
        let p = RetryPolicy::default();
        assert!(p.should_retry(ErrorKind::Convergence));
        assert!(p.should_retry(ErrorKind::Numerical));
        assert!(!p.should_retry(ErrorKind::Bounds));
        assert!(!p.should_retry(ErrorKind::Capacity));
        assert!(!p.should_retry(ErrorKind::Unsupported));
    }

    #[test]
    fn retry_policy_retries_for_kind() {
        let p = RetryPolicy::default();
        assert_eq!(p.retries_for(ErrorKind::Convergence), 2);
        assert_eq!(p.retries_for(ErrorKind::Numerical), 1);
        assert_eq!(p.retries_for(ErrorKind::Other), 0);
    }

    #[test]
    fn retry_policy_custom() {
        let p = RetryPolicy {
            max_convergence_retries: 5,
            max_numerical_retries: 3,
        };
        assert_eq!(p.retries_for(ErrorKind::Convergence), 5);
        assert_eq!(p.retries_for(ErrorKind::Numerical), 3);
        assert!(p.should_retry(ErrorKind::Convergence));
        assert!(p.should_retry(ErrorKind::Numerical));
    }

    // ── Edge-case / integration tests ────────────────────────────────────────

    #[test]
    fn collector_is_empty_initially() {
        let col = ErrorCollector::new();
        assert!(col.is_empty());
        assert_eq!(col.len(), 0);
    }

    #[test]
    fn error_kind_all_variants_have_recovery_hints() {
        // Smoke-test: all error variants produce a RecoveryHint without panicking
        let errors: Vec<Error> = vec![
            Error::General("g".into()),
            Error::GjkNoConverge {
                iterations: 1,
                dist: 0.0,
            },
            Error::EpaFailure { reason: "r".into() },
            Error::CcdNoConverge { toi: 0.0, iters: 0 },
            Error::DegenerateShape { detail: "d".into() },
            Error::IndexOutOfRange { index: 0, len: 0 },
            Error::BroadphaseCapacity { capacity: 0 },
            Error::NonFinite {
                context: "c".into(),
            },
            Error::UnsupportedShapePair {
                shape_a: "A".into(),
                shape_b: "B".into(),
            },
        ];
        for e in &errors {
            let _ = e.recovery_hint(); // must not panic
            let _ = e.kind();
            let _ = e.severity();
        }
    }

    #[test]
    fn collector_errors_above_info_includes_all() {
        let mut col = ErrorCollector::new();
        col.push(Error::General("a".into())).unwrap();
        col.push(Error::GjkNoConverge {
            iterations: 1,
            dist: 0.0,
        })
        .unwrap();
        col.push(Error::BroadphaseCapacity { capacity: 10 })
            .unwrap();
        // Info is the lowest severity — all non-fatal errors should appear
        let above_info = col.errors_above(Severity::Info);
        assert_eq!(above_info.len(), 3);
    }

    #[test]
    fn stats_count_severity_info() {
        let mut s = ErrorStats::new();
        s.record(&Error::General("x".into())); // severity = Info
        assert_eq!(s.count_severity(Severity::Info), 1);
        assert_eq!(s.count_severity(Severity::Warning), 0);
    }
}
