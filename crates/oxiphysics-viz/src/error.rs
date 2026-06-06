// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-viz.
//!
//! This module defines all error variants produced by the visualization subsystem,
//! along with conversion helpers and a typed `Result` alias.
//!
//! # Design
//!
//! Errors are grouped into several categories:
//! - **Geometry** – mesh/primitive construction failures.
//! - **Shader** – shader source or parameter validation failures.
//! - **Volume** – volume rendering / transfer function errors.
//! - **Colormap** – colormap parameter errors.
//! - **Scene** – scene-graph management errors.
//! - **IO** – file-system-related errors (stubbed, no actual IO dependency).

use std::fmt;
use thiserror::Error;

// ---------------------------------------------------------------------------
// Primary error type
// ---------------------------------------------------------------------------

/// All error variants produced by `oxiphysics-viz`.
#[derive(Debug, Clone, Error)]
pub enum Error {
    /// Generic string error — used when no specific variant applies.
    #[error("{0}")]
    General(String),

    /// A geometry operation (mesh building, primitive generation) failed.
    #[error("geometry error: {message}")]
    Geometry {
        /// Human-readable description of what went wrong.
        message: String,
    },

    /// A shader validation or parameter-binding error.
    #[error("shader error in '{shader}': {message}")]
    Shader {
        /// Name of the shader (e.g. `"phong_vs"` or `"pbr_fs"`).
        shader: String,
        /// Human-readable description.
        message: String,
    },

    /// A volume-rendering error (e.g. empty volume data, bad grid dimensions).
    #[error("volume rendering error: {message}")]
    Volume {
        /// Human-readable description.
        message: String,
    },

    /// A colormap parameter error (e.g. `min >= max`).
    #[error("colormap error: {message}")]
    Colormap {
        /// Human-readable description.
        message: String,
    },

    /// A scene-graph management error (e.g. node not found, invalid transform).
    #[error("scene error: {message}")]
    Scene {
        /// Human-readable description.
        message: String,
    },

    /// A post-processing pipeline error.
    #[error("post-processing error in stage '{stage}': {message}")]
    PostProcess {
        /// Name of the failing stage (e.g. `"bloom"`, `"ssao"`).
        stage: String,
        /// Human-readable description.
        message: String,
    },

    /// An attempted operation on an out-of-range index.
    #[error("index {index} out of range [0, {length})")]
    OutOfRange {
        /// The offending index.
        index: usize,
        /// The valid length.
        length: usize,
    },

    /// A resource (mesh, texture, shader program) was not found by name.
    #[error("resource not found: '{name}'")]
    NotFound {
        /// Name of the missing resource.
        name: String,
    },

    /// Duplicate resource registration.
    #[error("resource already exists: '{name}'")]
    AlreadyExists {
        /// Name of the duplicate resource.
        name: String,
    },

    /// Invalid parameter value.
    #[error("invalid parameter '{param}': {reason}")]
    InvalidParameter {
        /// Parameter name.
        param: String,
        /// Why it is invalid.
        reason: String,
    },
}

/// Result type alias for `oxiphysics-viz` operations.
pub type Result<T> = std::result::Result<T, Error>;

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

impl Error {
    /// Construct a [`Error::General`] from any `Display`-able value.
    pub fn general(msg: impl fmt::Display) -> Self {
        Error::General(msg.to_string())
    }

    /// Construct a [`Error::Geometry`] error.
    pub fn geometry(message: impl Into<String>) -> Self {
        Error::Geometry {
            message: message.into(),
        }
    }

    /// Construct a [`Error::Shader`] error.
    pub fn shader(shader: impl Into<String>, message: impl Into<String>) -> Self {
        Error::Shader {
            shader: shader.into(),
            message: message.into(),
        }
    }

    /// Construct a [`Error::Volume`] error.
    pub fn volume(message: impl Into<String>) -> Self {
        Error::Volume {
            message: message.into(),
        }
    }

    /// Construct a [`Error::Colormap`] error.
    pub fn colormap(message: impl Into<String>) -> Self {
        Error::Colormap {
            message: message.into(),
        }
    }

    /// Construct a [`Error::Scene`] error.
    pub fn scene(message: impl Into<String>) -> Self {
        Error::Scene {
            message: message.into(),
        }
    }

    /// Construct a [`Error::PostProcess`] error.
    pub fn post_process(stage: impl Into<String>, message: impl Into<String>) -> Self {
        Error::PostProcess {
            stage: stage.into(),
            message: message.into(),
        }
    }

    /// Construct an [`Error::OutOfRange`] error.
    pub fn out_of_range(index: usize, length: usize) -> Self {
        Error::OutOfRange { index, length }
    }

    /// Construct a [`Error::NotFound`] error.
    pub fn not_found(name: impl Into<String>) -> Self {
        Error::NotFound { name: name.into() }
    }

    /// Construct an [`Error::AlreadyExists`] error.
    pub fn already_exists(name: impl Into<String>) -> Self {
        Error::AlreadyExists { name: name.into() }
    }

    /// Construct an [`Error::InvalidParameter`] error.
    pub fn invalid_parameter(param: impl Into<String>, reason: impl Into<String>) -> Self {
        Error::InvalidParameter {
            param: param.into(),
            reason: reason.into(),
        }
    }

    /// Return `true` if this is a geometry error.
    pub fn is_geometry(&self) -> bool {
        matches!(self, Error::Geometry { .. })
    }

    /// Return `true` if this is a shader error.
    pub fn is_shader(&self) -> bool {
        matches!(self, Error::Shader { .. })
    }

    /// Return `true` if this is a not-found error.
    pub fn is_not_found(&self) -> bool {
        matches!(self, Error::NotFound { .. })
    }

    /// Return `true` if this is an out-of-range error.
    pub fn is_out_of_range(&self) -> bool {
        matches!(self, Error::OutOfRange { .. })
    }
}

// ---------------------------------------------------------------------------
// ErrorContext — wraps a Result with contextual information
// ---------------------------------------------------------------------------

/// A helper to attach additional context to any error.
///
/// Inspired by the `anyhow` / `context` pattern.
#[derive(Debug)]
pub struct ErrorContext {
    /// The underlying error.
    pub source: Error,
    /// Additional human-readable context.
    pub context: String,
}

impl fmt::Display for ErrorContext {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.context, self.source)
    }
}

impl ErrorContext {
    /// Wrap an error with a context string.
    pub fn new(source: Error, context: impl Into<String>) -> Self {
        Self {
            source,
            context: context.into(),
        }
    }
}

/// Extension trait that adds `.context(msg)` to `Result<T, Error>`.
pub trait ResultExt<T> {
    /// Attach a context string to an `Err` variant.
    fn context(self, ctx: &str) -> std::result::Result<T, ErrorContext>;
}

impl<T> ResultExt<T> for Result<T> {
    fn context(self, ctx: &str) -> std::result::Result<T, ErrorContext> {
        self.map_err(|e| ErrorContext::new(e, ctx))
    }
}

// ---------------------------------------------------------------------------
// ErrorCollection — accumulate multiple errors
// ---------------------------------------------------------------------------

/// A container that collects multiple errors during batch operations.
#[derive(Debug, Default)]
pub struct ErrorCollection {
    errors: Vec<Error>,
}

impl ErrorCollection {
    /// Create an empty collection.
    pub fn new() -> Self {
        Self { errors: Vec::new() }
    }

    /// Push an error into the collection.
    pub fn push(&mut self, e: Error) {
        self.errors.push(e);
    }

    /// Return `true` if no errors have been collected.
    pub fn is_empty(&self) -> bool {
        self.errors.is_empty()
    }

    /// Number of accumulated errors.
    pub fn len(&self) -> usize {
        self.errors.len()
    }

    /// Drain the collection into a `Vec`.
    pub fn into_vec(self) -> Vec<Error> {
        self.errors
    }

    /// Convert to a single `Result`: `Ok(())` if empty, otherwise `Err(General(...))`.
    pub fn into_result(self) -> Result<()> {
        if self.errors.is_empty() {
            Ok(())
        } else {
            let msg = self
                .errors
                .iter()
                .enumerate()
                .map(|(i, e)| format!("[{}] {}", i, e))
                .collect::<Vec<_>>()
                .join("; ");
            Err(Error::general(msg))
        }
    }
}

impl fmt::Display for ErrorCollection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} error(s)", self.errors.len())?;
        for (i, e) in self.errors.iter().enumerate() {
            write!(f, "\n  [{}] {}", i, e)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// VizErrorCode — machine-readable error codes
// ---------------------------------------------------------------------------

/// Machine-readable code associated with a viz error, for programmatic handling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VizErrorCode {
    /// General unclassified error.
    General = 0,
    /// Geometry failure.
    Geometry = 1,
    /// Shader failure.
    Shader = 2,
    /// Volume rendering failure.
    Volume = 3,
    /// Colormap failure.
    Colormap = 4,
    /// Scene management failure.
    Scene = 5,
    /// Post-processing failure.
    PostProcess = 6,
    /// Out-of-range index.
    OutOfRange = 7,
    /// Resource not found.
    NotFound = 8,
    /// Resource already exists.
    AlreadyExists = 9,
    /// Invalid parameter.
    InvalidParameter = 10,
}

impl VizErrorCode {
    /// Map a [`enum@Error`] variant to its numeric code.
    pub fn from_error(e: &Error) -> Self {
        match e {
            Error::General(_) => VizErrorCode::General,
            Error::Geometry { .. } => VizErrorCode::Geometry,
            Error::Shader { .. } => VizErrorCode::Shader,
            Error::Volume { .. } => VizErrorCode::Volume,
            Error::Colormap { .. } => VizErrorCode::Colormap,
            Error::Scene { .. } => VizErrorCode::Scene,
            Error::PostProcess { .. } => VizErrorCode::PostProcess,
            Error::OutOfRange { .. } => VizErrorCode::OutOfRange,
            Error::NotFound { .. } => VizErrorCode::NotFound,
            Error::AlreadyExists { .. } => VizErrorCode::AlreadyExists,
            Error::InvalidParameter { .. } => VizErrorCode::InvalidParameter,
        }
    }

    /// Numeric code value.
    pub fn code(self) -> u32 {
        self as u32
    }
}

impl fmt::Display for VizErrorCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "VIZ-{:04}", self.code())
    }
}

// ---------------------------------------------------------------------------
// ErrorSeverity — classify how severe an error is
// ---------------------------------------------------------------------------

/// Severity level of a visualization error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ErrorSeverity {
    /// Informational — not a real error, but worth recording.
    Info = 0,
    /// Warning — operation succeeded but something unexpected happened.
    Warning = 1,
    /// Error — operation failed; the caller must handle this.
    Error = 2,
    /// Fatal — unrecoverable; the subsystem must be torn down.
    Fatal = 3,
}

impl ErrorSeverity {
    /// Returns `true` if this severity indicates a genuine failure.
    pub fn is_failure(&self) -> bool {
        *self >= ErrorSeverity::Error
    }

    /// Human-readable label.
    pub fn label(&self) -> &'static str {
        match self {
            ErrorSeverity::Info => "INFO",
            ErrorSeverity::Warning => "WARN",
            ErrorSeverity::Error => "ERROR",
            ErrorSeverity::Fatal => "FATAL",
        }
    }

    /// Map an [`enum@Error`] variant to its default severity.
    pub fn from_error(e: &Error) -> Self {
        match e {
            Error::General(_) => ErrorSeverity::Error,
            Error::Geometry { .. } => ErrorSeverity::Error,
            Error::Shader { .. } => ErrorSeverity::Error,
            Error::Volume { .. } => ErrorSeverity::Error,
            Error::Colormap { .. } => ErrorSeverity::Warning,
            Error::Scene { .. } => ErrorSeverity::Error,
            Error::PostProcess { .. } => ErrorSeverity::Warning,
            Error::OutOfRange { .. } => ErrorSeverity::Error,
            Error::NotFound { .. } => ErrorSeverity::Error,
            Error::AlreadyExists { .. } => ErrorSeverity::Warning,
            Error::InvalidParameter { .. } => ErrorSeverity::Error,
        }
    }
}

impl fmt::Display for ErrorSeverity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label())
    }
}

// ---------------------------------------------------------------------------
// RetryHint — can this operation be retried?
// ---------------------------------------------------------------------------

/// Hint indicating whether an operation that returned an error can be retried.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RetryHint {
    /// The error is permanent; retrying will always fail.
    Never,
    /// The error might be transient; caller may retry immediately.
    Immediately,
    /// The error is transient but the caller should wait before retrying.
    AfterDelay,
}

impl Error {
    /// Suggest a retry strategy for this error.
    pub fn retry_hint(&self) -> RetryHint {
        match self {
            Error::General(_) => RetryHint::Never,
            Error::Geometry { .. } => RetryHint::Never,
            Error::Shader { .. } => RetryHint::Never,
            Error::Volume { .. } => RetryHint::Never,
            Error::Colormap { .. } => RetryHint::Never,
            Error::Scene { .. } => RetryHint::Immediately,
            Error::PostProcess { .. } => RetryHint::Immediately,
            Error::OutOfRange { .. } => RetryHint::Never,
            Error::NotFound { .. } => RetryHint::AfterDelay,
            Error::AlreadyExists { .. } => RetryHint::Never,
            Error::InvalidParameter { .. } => RetryHint::Never,
        }
    }

    /// Return the severity of this error.
    pub fn severity(&self) -> ErrorSeverity {
        ErrorSeverity::from_error(self)
    }

    /// Return a concise fix suggestion for this error.
    pub fn fix_hint(&self) -> &'static str {
        match self {
            Error::General(_) => "Check the operation parameters.",
            Error::Geometry { .. } => "Verify mesh input (no degenerate triangles, valid indices).",
            Error::Shader { .. } => "Review shader source for syntax/uniform errors.",
            Error::Volume { .. } => "Ensure volume data is non-empty and dimensions are valid.",
            Error::Colormap { .. } => "Ensure scalar_min < scalar_max.",
            Error::Scene { .. } => "Check scene node existence before operating on it.",
            Error::PostProcess { .. } => "Verify post-process stage parameters.",
            Error::OutOfRange { .. } => "Check the index against the container length.",
            Error::NotFound { .. } => "Register the resource before using it.",
            Error::AlreadyExists { .. } => {
                "Use a unique resource name or remove the existing one first."
            }
            Error::InvalidParameter { .. } => "Check parameter constraints in the documentation.",
        }
    }
}

// ---------------------------------------------------------------------------
// ErrorRecord — a single timestamped error entry
// ---------------------------------------------------------------------------

/// A single error entry with optional source location metadata.
#[derive(Debug, Clone)]
pub struct ErrorRecord {
    /// The underlying error.
    pub error: Error,
    /// Severity of the error.
    pub severity: ErrorSeverity,
    /// Optional source file name.
    pub file: Option<&'static str>,
    /// Optional source line number.
    pub line: Option<u32>,
    /// Optional module path.
    pub module: Option<&'static str>,
    /// Sequence number (insertion order).
    pub seq: usize,
}

impl ErrorRecord {
    /// Construct an `ErrorRecord` with full source location.
    pub fn with_location(
        error: Error,
        file: Option<&'static str>,
        line: Option<u32>,
        module: Option<&'static str>,
        seq: usize,
    ) -> Self {
        let severity = error.severity();
        Self {
            error,
            severity,
            file,
            line,
            module,
            seq,
        }
    }

    /// Construct a basic `ErrorRecord` without source location.
    pub fn new(error: Error, seq: usize) -> Self {
        let severity = error.severity();
        Self {
            error,
            severity,
            file: None,
            line: None,
            module: None,
            seq,
        }
    }
}

impl fmt::Display for ErrorRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}][{}] {}", self.seq, self.severity, self.error)?;
        if let (Some(file), Some(line)) = (self.file, self.line) {
            write!(f, " ({}:{})", file, line)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// ErrorReport — a rich diagnostic report
// ---------------------------------------------------------------------------

/// A structured diagnostic report that aggregates `ErrorRecord`s and
/// provides filtering, severity summaries, and formatted output.
#[derive(Debug, Default)]
pub struct ErrorReport {
    records: Vec<ErrorRecord>,
    next_seq: usize,
}

impl ErrorReport {
    /// Create an empty report.
    pub fn new() -> Self {
        Self {
            records: Vec::new(),
            next_seq: 0,
        }
    }

    /// Push an error into the report.
    pub fn push(&mut self, error: Error) {
        let seq = self.next_seq;
        self.next_seq += 1;
        self.records.push(ErrorRecord::new(error, seq));
    }

    /// Push an error with source location.
    pub fn push_with_location(
        &mut self,
        error: Error,
        file: &'static str,
        line: u32,
        module: &'static str,
    ) {
        let seq = self.next_seq;
        self.next_seq += 1;
        self.records.push(ErrorRecord::with_location(
            error,
            Some(file),
            Some(line),
            Some(module),
            seq,
        ));
    }

    /// Number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Return `true` if no records exist.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Count records at or above the given severity.
    pub fn count_at_least(&self, severity: ErrorSeverity) -> usize {
        self.records
            .iter()
            .filter(|r| r.severity >= severity)
            .count()
    }

    /// Return `true` if any record is an [`ErrorSeverity::Error`] or higher.
    pub fn has_errors(&self) -> bool {
        self.count_at_least(ErrorSeverity::Error) > 0
    }

    /// Return `true` if any record is [`ErrorSeverity::Fatal`].
    pub fn has_fatal(&self) -> bool {
        self.count_at_least(ErrorSeverity::Fatal) > 0
    }

    /// Iterate over all records.
    pub fn records(&self) -> &[ErrorRecord] {
        &self.records
    }

    /// Filter records by severity level (at least `min_severity`).
    pub fn filter_by_severity(&self, min_severity: ErrorSeverity) -> Vec<&ErrorRecord> {
        self.records
            .iter()
            .filter(|r| r.severity >= min_severity)
            .collect()
    }

    /// Convert to a single `Result<()>`:
    /// - `Ok(())` if no errors at `Error` level or above.
    /// - `Err(General(...))` listing all error/fatal messages.
    pub fn into_result(self) -> Result<()> {
        let failures: Vec<_> = self
            .records
            .iter()
            .filter(|r| r.severity >= ErrorSeverity::Error)
            .collect();
        if failures.is_empty() {
            return Ok(());
        }
        let msg = failures
            .iter()
            .map(|r| format!("{}", r))
            .collect::<Vec<_>>()
            .join("; ");
        Err(Error::general(msg))
    }

    /// Format the full report as a multi-line string.
    pub fn format_full(&self) -> String {
        let mut out = format!("ErrorReport ({} record(s)):\n", self.records.len());
        for r in &self.records {
            out.push_str(&format!("  {}\n", r));
            out.push_str(&format!("    hint: {}\n", r.error.fix_hint()));
        }
        out
    }

    /// Produce a compact one-line summary.
    pub fn summary(&self) -> String {
        let fatals = self.count_at_least(ErrorSeverity::Fatal);
        let errors = self.count_at_least(ErrorSeverity::Error) - fatals;
        let warnings = self.count_at_least(ErrorSeverity::Warning) - errors - fatals;
        format!(
            "{} fatal, {} error(s), {} warning(s)",
            fatals, errors, warnings
        )
    }
}

impl fmt::Display for ErrorReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.format_full())
    }
}

// ---------------------------------------------------------------------------
// ErrorChain — walk a chain of ErrorContext wrappers
// ---------------------------------------------------------------------------

/// A chain of context-wrapped errors, from innermost to outermost.
#[derive(Debug)]
pub struct ErrorChain {
    /// Ordered list: `[innermost_context, ..., outermost_context]`.
    pub contexts: Vec<String>,
    /// The root cause error.
    pub root: Error,
}

impl ErrorChain {
    /// Build an `ErrorChain` by unwrapping nested [`ErrorContext`]s.
    pub fn from_context(ctx: ErrorContext) -> Self {
        let contexts = vec![ctx.context.clone()];
        // In the current design ErrorContext holds a single Error (not nested).
        // If the source were another ErrorContext we would recurse; here we just
        // record the single context layer.
        Self {
            contexts,
            root: ctx.source,
        }
    }

    /// Return `true` if the root cause matches the given error code.
    pub fn root_code(&self) -> VizErrorCode {
        VizErrorCode::from_error(&self.root)
    }

    /// Format the chain in a breadcrumb style.
    pub fn format_breadcrumb(&self) -> String {
        let breadcrumb = self.contexts.join(" → ");
        format!("{}: {}", breadcrumb, self.root)
    }
}

impl fmt::Display for ErrorChain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.format_breadcrumb())
    }
}

// ---------------------------------------------------------------------------
// LogAdapter — format errors for different logging targets
// ---------------------------------------------------------------------------

/// Target for structured log output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogTarget {
    /// Plain human-readable text.
    Text,
    /// JSON object (single line).
    Json,
    /// CSV row (seq, severity, code, message).
    Csv,
}

/// Formats an `ErrorRecord` for a specified logging target.
pub fn format_log_entry(record: &ErrorRecord, target: LogTarget) -> String {
    let code = VizErrorCode::from_error(&record.error);
    match target {
        LogTarget::Text => {
            format!(
                "[{}] {} {} | {} | hint: {}",
                record.seq,
                record.severity,
                code,
                record.error,
                record.error.fix_hint(),
            )
        }
        LogTarget::Json => {
            let file = record.file.unwrap_or("unknown");
            let line = record.line.unwrap_or(0);
            let module = record.module.unwrap_or("unknown");
            format!(
                "{{\"seq\":{},\"severity\":\"{}\",\"code\":\"{}\",\"message\":\"{}\",\"file\":\"{}\",\"line\":{},\"module\":\"{}\"}}",
                record.seq,
                record.severity,
                code,
                record.error.to_string().replace('"', "\\\""),
                file,
                line,
                module,
            )
        }
        LogTarget::Csv => {
            format!(
                "{},{},{},\"{}\"",
                record.seq,
                record.severity,
                code,
                record.error.to_string().replace('"', "\"\""),
            )
        }
    }
}

// ---------------------------------------------------------------------------
// ErrorFilter — configurable severity filter
// ---------------------------------------------------------------------------

/// A filter that only passes errors at or above a minimum severity.
#[derive(Debug, Clone)]
pub struct ErrorFilter {
    /// Minimum severity to pass through.
    pub min_severity: ErrorSeverity,
    /// If `true`, also suppress errors matching any of the listed codes.
    pub suppressed_codes: Vec<VizErrorCode>,
}

impl Default for ErrorFilter {
    fn default() -> Self {
        Self {
            min_severity: ErrorSeverity::Warning,
            suppressed_codes: Vec::new(),
        }
    }
}

impl ErrorFilter {
    /// Construct a filter that passes all severities.
    pub fn pass_all() -> Self {
        Self {
            min_severity: ErrorSeverity::Info,
            suppressed_codes: Vec::new(),
        }
    }

    /// Suppress a specific error code.
    pub fn suppress(mut self, code: VizErrorCode) -> Self {
        self.suppressed_codes.push(code);
        self
    }

    /// Return `true` if the given error passes through this filter.
    pub fn passes(&self, error: &Error) -> bool {
        let sev = ErrorSeverity::from_error(error);
        if sev < self.min_severity {
            return false;
        }
        let code = VizErrorCode::from_error(error);
        !self.suppressed_codes.contains(&code)
    }

    /// Filter a slice of errors, returning only those that pass.
    pub fn apply<'a>(&self, errors: &'a [Error]) -> Vec<&'a Error> {
        errors.iter().filter(|e| self.passes(e)).collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_general_error_display() {
        let e = Error::general("something broke");
        assert!(e.to_string().contains("something broke"));
    }

    #[test]
    fn test_geometry_error_display() {
        let e = Error::geometry("degenerate triangle");
        assert!(e.to_string().contains("geometry error"));
        assert!(e.to_string().contains("degenerate triangle"));
        assert!(e.is_geometry());
        assert!(!e.is_shader());
    }

    #[test]
    fn test_shader_error_display() {
        let e = Error::shader("phong_vs", "missing uniform u_mvp");
        let s = e.to_string();
        assert!(s.contains("phong_vs"));
        assert!(s.contains("missing uniform u_mvp"));
        assert!(e.is_shader());
    }

    #[test]
    fn test_volume_error() {
        let e = Error::volume("grid is empty");
        assert!(e.to_string().contains("grid is empty"));
        assert!(!e.is_geometry());
    }

    #[test]
    fn test_out_of_range_error() {
        let e = Error::out_of_range(10, 5);
        let s = e.to_string();
        assert!(s.contains("10"));
        assert!(s.contains("5"));
        assert!(e.is_out_of_range());
    }

    #[test]
    fn test_not_found_error() {
        let e = Error::not_found("my_mesh");
        assert!(e.to_string().contains("my_mesh"));
        assert!(e.is_not_found());
        assert!(!e.is_geometry());
    }

    #[test]
    fn test_already_exists_error() {
        let e = Error::already_exists("my_texture");
        assert!(e.to_string().contains("my_texture"));
    }

    #[test]
    fn test_invalid_parameter_error() {
        let e = Error::invalid_parameter("roughness", "must be in [0, 1]");
        let s = e.to_string();
        assert!(s.contains("roughness"));
        assert!(s.contains("must be in [0, 1]"));
    }

    #[test]
    fn test_post_process_error() {
        let e = Error::post_process("bloom", "radius too large");
        let s = e.to_string();
        assert!(s.contains("bloom"));
        assert!(s.contains("radius too large"));
    }

    #[test]
    fn test_scene_error() {
        let e = Error::scene("node not found");
        assert!(e.to_string().contains("node not found"));
    }

    #[test]
    fn test_error_collection_empty() {
        let col = ErrorCollection::new();
        assert!(col.is_empty());
        assert_eq!(col.len(), 0);
        assert!(col.into_result().is_ok());
    }

    #[test]
    fn test_error_collection_with_errors() {
        let mut col = ErrorCollection::new();
        col.push(Error::geometry("bad vertex"));
        col.push(Error::shader("vs", "compile error"));
        assert_eq!(col.len(), 2);
        assert!(!col.is_empty());
        let result = col.into_result();
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("[0]"));
        assert!(msg.contains("[1]"));
    }

    #[test]
    fn test_error_collection_display() {
        let mut col = ErrorCollection::new();
        col.push(Error::general("err1"));
        let s = col.to_string();
        assert!(s.contains("1 error"));
    }

    #[test]
    fn test_error_context_wrapper() {
        let result: Result<i32> = Err(Error::geometry("bad normal"));
        let ctx_result = result.context("while building sphere");
        let ctx_err = ctx_result.unwrap_err();
        let s = ctx_err.to_string();
        assert!(s.contains("while building sphere"));
        assert!(s.contains("bad normal"));
    }

    #[test]
    fn test_error_context_ok_passthrough() {
        let result: Result<i32> = Ok(42);
        let ctx_result = result.context("should not appear");
        assert_eq!(ctx_result.unwrap(), 42);
    }

    #[test]
    fn test_viz_error_code_mapping() {
        let e = Error::geometry("test");
        assert_eq!(VizErrorCode::from_error(&e), VizErrorCode::Geometry);

        let e2 = Error::not_found("x");
        assert_eq!(VizErrorCode::from_error(&e2), VizErrorCode::NotFound);

        let e3 = Error::general("x");
        assert_eq!(VizErrorCode::from_error(&e3), VizErrorCode::General);
    }

    #[test]
    fn test_viz_error_code_display() {
        let code = VizErrorCode::Geometry;
        assert_eq!(code.to_string(), "VIZ-0001");

        let code2 = VizErrorCode::General;
        assert_eq!(code2.to_string(), "VIZ-0000");
    }

    #[test]
    fn test_viz_error_code_numeric() {
        assert_eq!(VizErrorCode::Shader.code(), 2);
        assert_eq!(VizErrorCode::NotFound.code(), 8);
        assert_eq!(VizErrorCode::InvalidParameter.code(), 10);
    }

    #[test]
    fn test_colormap_error() {
        let e = Error::colormap("min must be less than max");
        assert!(e.to_string().contains("min must be less than max"));
    }

    #[test]
    fn test_error_collection_into_vec() {
        let mut col = ErrorCollection::new();
        col.push(Error::general("a"));
        col.push(Error::general("b"));
        let v = col.into_vec();
        assert_eq!(v.len(), 2);
    }

    // ---------------------------------------------------------------------------
    // ErrorSeverity tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_severity_ordering() {
        assert!(ErrorSeverity::Info < ErrorSeverity::Warning);
        assert!(ErrorSeverity::Warning < ErrorSeverity::Error);
        assert!(ErrorSeverity::Error < ErrorSeverity::Fatal);
    }

    #[test]
    fn test_severity_is_failure() {
        assert!(!ErrorSeverity::Info.is_failure());
        assert!(!ErrorSeverity::Warning.is_failure());
        assert!(ErrorSeverity::Error.is_failure());
        assert!(ErrorSeverity::Fatal.is_failure());
    }

    #[test]
    fn test_severity_labels() {
        assert_eq!(ErrorSeverity::Info.label(), "INFO");
        assert_eq!(ErrorSeverity::Warning.label(), "WARN");
        assert_eq!(ErrorSeverity::Error.label(), "ERROR");
        assert_eq!(ErrorSeverity::Fatal.label(), "FATAL");
    }

    #[test]
    fn test_severity_from_error_geometry() {
        let e = Error::geometry("bad mesh");
        assert_eq!(ErrorSeverity::from_error(&e), ErrorSeverity::Error);
    }

    #[test]
    fn test_severity_from_error_colormap_is_warning() {
        let e = Error::colormap("bad range");
        assert_eq!(ErrorSeverity::from_error(&e), ErrorSeverity::Warning);
    }

    #[test]
    fn test_severity_display() {
        let s = ErrorSeverity::Fatal.to_string();
        assert_eq!(s, "FATAL");
    }

    // ---------------------------------------------------------------------------
    // RetryHint tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_retry_hint_geometry_never() {
        let e = Error::geometry("degenerate");
        assert_eq!(e.retry_hint(), RetryHint::Never);
    }

    #[test]
    fn test_retry_hint_not_found_after_delay() {
        let e = Error::not_found("resource");
        assert_eq!(e.retry_hint(), RetryHint::AfterDelay);
    }

    #[test]
    fn test_retry_hint_scene_immediately() {
        let e = Error::scene("temp glitch");
        assert_eq!(e.retry_hint(), RetryHint::Immediately);
    }

    // ---------------------------------------------------------------------------
    // ErrorRecord tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_error_record_new() {
        let rec = ErrorRecord::new(Error::geometry("oops"), 0);
        assert_eq!(rec.seq, 0);
        assert_eq!(rec.severity, ErrorSeverity::Error);
        assert!(rec.file.is_none());
    }

    #[test]
    fn test_error_record_with_location() {
        let rec = ErrorRecord::with_location(
            Error::shader("vs", "bad"),
            Some("foo.rs"),
            Some(42),
            Some("my_module"),
            1,
        );
        assert_eq!(rec.seq, 1);
        assert_eq!(rec.file, Some("foo.rs"));
        assert_eq!(rec.line, Some(42));
        let s = rec.to_string();
        assert!(s.contains("foo.rs"));
        assert!(s.contains("42"));
    }

    // ---------------------------------------------------------------------------
    // ErrorReport tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_error_report_empty() {
        let rep = ErrorReport::new();
        assert!(rep.is_empty());
        assert_eq!(rep.len(), 0);
        assert!(!rep.has_errors());
        assert!(rep.into_result().is_ok());
    }

    #[test]
    fn test_error_report_push_and_has_errors() {
        let mut rep = ErrorReport::new();
        rep.push(Error::geometry("fail"));
        assert!(!rep.is_empty());
        assert!(rep.has_errors());
    }

    #[test]
    fn test_error_report_only_warnings_not_errors() {
        let mut rep = ErrorReport::new();
        rep.push(Error::colormap("range")); // Warning severity
        rep.push(Error::already_exists("x")); // Warning severity
        assert!(!rep.has_errors());
        assert!(rep.into_result().is_ok());
    }

    #[test]
    fn test_error_report_filter_by_severity() {
        let mut rep = ErrorReport::new();
        rep.push(Error::colormap("warn level"));
        rep.push(Error::geometry("error level"));
        let errors_only = rep.filter_by_severity(ErrorSeverity::Error);
        assert_eq!(errors_only.len(), 1);
    }

    #[test]
    fn test_error_report_summary() {
        let mut rep = ErrorReport::new();
        rep.push(Error::geometry("e1"));
        rep.push(Error::geometry("e2"));
        rep.push(Error::colormap("w1"));
        let s = rep.summary();
        assert!(s.contains("error") || s.contains("warn"), "summary: {s}");
    }

    #[test]
    fn test_error_report_into_result_includes_messages() {
        let mut rep = ErrorReport::new();
        rep.push(Error::geometry("bad input"));
        let res = rep.into_result();
        let msg = res.unwrap_err().to_string();
        assert!(msg.contains("bad input"));
    }

    #[test]
    fn test_error_report_push_with_location() {
        let mut rep = ErrorReport::new();
        rep.push_with_location(
            Error::shader("vs", "compile fail"),
            "shader.rs",
            100,
            "renderer",
        );
        assert_eq!(rep.len(), 1);
        let full = rep.format_full();
        assert!(full.contains("shader.rs"));
    }

    // ---------------------------------------------------------------------------
    // ErrorChain tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_error_chain_breadcrumb() {
        let inner = Error::geometry("bad normal");
        let ctx = ErrorContext::new(inner, "while building capsule");
        let chain = ErrorChain::from_context(ctx);
        let s = chain.format_breadcrumb();
        assert!(s.contains("while building capsule"));
        assert!(s.contains("bad normal"));
        assert_eq!(chain.root_code(), VizErrorCode::Geometry);
    }

    #[test]
    fn test_error_chain_display() {
        let ctx = ErrorContext::new(Error::volume("empty"), "ray march setup");
        let chain = ErrorChain::from_context(ctx);
        let s = chain.to_string();
        assert!(s.contains("ray march setup"));
        assert!(s.contains("empty"));
    }

    // ---------------------------------------------------------------------------
    // LogAdapter tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_log_format_text() {
        let rec = ErrorRecord::new(Error::geometry("bad"), 3);
        let s = format_log_entry(&rec, LogTarget::Text);
        assert!(s.contains("[3]"));
        assert!(s.contains("bad"));
    }

    #[test]
    fn test_log_format_json() {
        let rec = ErrorRecord::new(Error::shader("ps", "error"), 0);
        let s = format_log_entry(&rec, LogTarget::Json);
        assert!(s.starts_with('{'));
        assert!(s.contains("\"seq\":0"));
        assert!(s.contains("\"severity\""));
    }

    #[test]
    fn test_log_format_csv() {
        let rec = ErrorRecord::new(Error::volume("data empty"), 7);
        let s = format_log_entry(&rec, LogTarget::Csv);
        // CSV: seq,severity,code,"message"
        assert!(s.starts_with("7,"));
        assert!(s.contains("data empty"));
    }

    // ---------------------------------------------------------------------------
    // ErrorFilter tests
    // ---------------------------------------------------------------------------

    #[test]
    fn test_error_filter_default_passes_warnings() {
        let f = ErrorFilter::default();
        assert!(f.passes(&Error::geometry("e")));
        assert!(f.passes(&Error::colormap("w")));
    }

    #[test]
    fn test_error_filter_suppress_code() {
        let f = ErrorFilter::default().suppress(VizErrorCode::Colormap);
        // Colormap errors are warning level AND suppressed → should not pass
        assert!(!f.passes(&Error::colormap("suppressed")));
        // Geometry errors pass (not suppressed)
        assert!(f.passes(&Error::geometry("ok")));
    }

    #[test]
    fn test_error_filter_pass_all() {
        let f = ErrorFilter::pass_all();
        // Even Info-level events pass
        assert!(f.passes(&Error::general("any")));
    }

    #[test]
    fn test_error_filter_apply_slice() {
        let errors = vec![
            Error::geometry("e1"),
            Error::colormap("w1"), // Warning
            Error::geometry("e2"),
        ];
        let f = ErrorFilter::default(); // passes Warning+
        let passing = f.apply(&errors);
        assert_eq!(passing.len(), 3); // all pass at Warning+ threshold
    }

    #[test]
    fn test_error_filter_high_threshold() {
        let f = ErrorFilter {
            min_severity: ErrorSeverity::Fatal,
            ..Default::default()
        };
        let errors = vec![
            Error::geometry("error level"),
            Error::colormap("warning level"),
        ];
        let passing = f.apply(&errors);
        // Neither geometry (Error) nor colormap (Warning) reaches Fatal
        assert_eq!(passing.len(), 0);
    }

    // ---------------------------------------------------------------------------
    // fix_hint / severity coherence
    // ---------------------------------------------------------------------------

    #[test]
    fn test_fix_hint_non_empty() {
        let errors = vec![
            Error::general("g"),
            Error::geometry("g"),
            Error::shader("s", "m"),
            Error::volume("v"),
            Error::colormap("c"),
            Error::scene("s"),
            Error::post_process("p", "m"),
            Error::out_of_range(0, 5),
            Error::not_found("n"),
            Error::already_exists("a"),
            Error::invalid_parameter("p", "r"),
        ];
        for e in &errors {
            assert!(
                !e.fix_hint().is_empty(),
                "fix_hint should not be empty for {:?}",
                e
            );
        }
    }

    #[test]
    fn test_error_severity_method() {
        let e = Error::geometry("test");
        assert_eq!(e.severity(), ErrorSeverity::Error);
        let w = Error::colormap("test");
        assert_eq!(w.severity(), ErrorSeverity::Warning);
    }
}
