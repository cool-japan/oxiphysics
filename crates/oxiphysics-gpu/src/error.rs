// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-gpu

use thiserror::Error;

/// Main error type for the gpu module.
#[derive(Debug, Error)]
pub enum Error {
    /// Generic error with a free-form message.
    #[error("{0}")]
    General(String),

    /// A GPU buffer allocation failed.
    #[error(
        "buffer allocation failed: requested {requested_bytes} bytes (available {available_bytes})"
    )]
    BufferAllocationFailed {
        /// Bytes requested.
        requested_bytes: usize,
        /// Bytes actually available.
        available_bytes: usize,
    },

    /// An invalid buffer handle was used.
    #[error("invalid buffer handle: {0}")]
    InvalidBufferHandle(usize),

    /// A shader compilation error (mock).
    #[error("shader compilation error in '{shader}': {message}")]
    ShaderCompilationError {
        /// Name of the offending shader.
        shader: String,
        /// Compiler message.
        message: String,
    },

    /// A dispatch exceeded the hardware work-group limit.
    #[error("dispatch size {dispatch_size} exceeds hardware limit {limit}")]
    DispatchLimitExceeded {
        /// Requested dispatch size (number of work-groups).
        dispatch_size: usize,
        /// Hardware maximum.
        limit: usize,
    },

    /// Out-of-bounds grid access.
    #[error("grid index ({i}, {j}, {k}) out of bounds for grid ({nx}, {ny}, {nz})")]
    GridIndexOutOfBounds {
        /// Requested x index.
        i: usize,
        /// Requested y index.
        j: usize,
        /// Requested z index.
        k: usize,
        /// Grid x dimension.
        nx: usize,
        /// Grid y dimension.
        ny: usize,
        /// Grid z dimension.
        nz: usize,
    },

    /// A kernel argument count mismatch.
    #[error("kernel '{kernel}' expects {expected} arguments but got {got}")]
    KernelArgCountMismatch {
        /// Kernel name.
        kernel: String,
        /// Expected number of arguments.
        expected: usize,
        /// Provided number of arguments.
        got: usize,
    },

    /// An unsupported backend feature was requested.
    #[error("unsupported feature: {feature}")]
    UnsupportedFeature {
        /// Description of the unsupported feature.
        feature: String,
    },
}

/// A pipeline-stage error: carries the stage name plus the underlying cause.
#[derive(Debug, Error)]
#[error("pipeline stage '{stage}' failed: {source}")]
pub struct PipelineStageError {
    /// Name of the pipeline stage (e.g. `"vertex_fetch"`, `"sph_density"`).
    pub stage: String,
    /// Root cause.
    pub source: Box<Error>,
}

/// Severity level for a GPU error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    /// Informational — execution can continue.
    Info,
    /// Warning — partial results may be degraded.
    Warning,
    /// Fatal — must abort current dispatch.
    Fatal,
}

impl std::fmt::Display for ErrorSeverity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ErrorSeverity::Info => write!(f, "INFO"),
            ErrorSeverity::Warning => write!(f, "WARNING"),
            ErrorSeverity::Fatal => write!(f, "FATAL"),
        }
    }
}

/// An error annotated with severity and an optional kernel name.
#[derive(Debug)]
pub struct AnnotatedError {
    /// The underlying error.
    pub error: Error,
    /// Severity classification.
    pub severity: ErrorSeverity,
    /// Optional kernel that triggered the error.
    pub kernel: Option<String>,
}

impl std::fmt::Display for AnnotatedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(ref k) = self.kernel {
            write!(f, "[{}] kernel '{}': {}", self.severity, k, self.error)
        } else {
            write!(f, "[{}] {}", self.severity, self.error)
        }
    }
}

impl AnnotatedError {
    /// Wrap an error as fatal with an optional kernel label.
    pub fn fatal(error: Error, kernel: Option<&str>) -> Self {
        Self {
            error,
            severity: ErrorSeverity::Fatal,
            kernel: kernel.map(str::to_string),
        }
    }

    /// Wrap an error as a warning.
    pub fn warning(error: Error, kernel: Option<&str>) -> Self {
        Self {
            error,
            severity: ErrorSeverity::Warning,
            kernel: kernel.map(str::to_string),
        }
    }
}

/// Result type alias
pub type Result<T> = std::result::Result<T, Error>;

/// GPU operation error (used by `LbmGpuSolver`, `BvhGpuTraverser`, etc.).
#[derive(Debug, Error)]
pub enum GpuError {
    /// GPU backend initialisation failed.
    #[error("GPU backend init failed: {0}")]
    BackendInit(String),

    /// Shader dispatch / pipeline error.
    #[error("shader dispatch error: {0}")]
    ShaderDispatch(String),

    /// Read-back from the GPU buffer failed.
    #[error("GPU buffer read-back failed: {0}")]
    ReadBack(String),

    /// An invalid buffer handle was used.
    #[error("invalid GPU buffer handle: {0}")]
    InvalidHandle(usize),
}

impl Error {
    /// Construct a [`Error::General`] from any `Display`-able value.
    pub fn general(msg: impl std::fmt::Display) -> Self {
        Error::General(msg.to_string())
    }

    /// True when this is a recoverable allocation error.
    pub fn is_allocation_error(&self) -> bool {
        matches!(self, Error::BufferAllocationFailed { .. })
    }

    /// True when this is a shader compilation error.
    pub fn is_shader_error(&self) -> bool {
        matches!(self, Error::ShaderCompilationError { .. })
    }

    /// True when this is an out-of-bounds grid error.
    pub fn is_grid_error(&self) -> bool {
        matches!(self, Error::GridIndexOutOfBounds { .. })
    }

    /// True when this is a kernel argument mismatch error.
    pub fn is_arg_mismatch(&self) -> bool {
        matches!(self, Error::KernelArgCountMismatch { .. })
    }

    /// True when this is an unsupported feature error.
    pub fn is_unsupported(&self) -> bool {
        matches!(self, Error::UnsupportedFeature { .. })
    }

    /// Wrap this error in a [`PipelineStageError`].
    pub fn in_stage(self, stage: impl Into<String>) -> PipelineStageError {
        PipelineStageError {
            stage: stage.into(),
            source: Box::new(self),
        }
    }

    /// Annotate with fatal severity.
    pub fn fatal(self, kernel: Option<&str>) -> AnnotatedError {
        AnnotatedError::fatal(self, kernel)
    }

    /// Annotate with warning severity.
    pub fn warning(self, kernel: Option<&str>) -> AnnotatedError {
        AnnotatedError::warning(self, kernel)
    }

    /// Convert into a `Result::Err`.
    pub fn into_err<T>(self) -> Result<T> {
        Err(self)
    }
}

// ── Convenience constructors ─────────────────────────────────────────────────

/// Build a [`Error::BufferAllocationFailed`] error.
pub fn alloc_err(requested_bytes: usize, available_bytes: usize) -> Error {
    Error::BufferAllocationFailed {
        requested_bytes,
        available_bytes,
    }
}

/// Build a [`Error::KernelArgCountMismatch`] error.
pub fn arg_mismatch_err(kernel: impl Into<String>, expected: usize, got: usize) -> Error {
    Error::KernelArgCountMismatch {
        kernel: kernel.into(),
        expected,
        got,
    }
}

/// Build a [`Error::GridIndexOutOfBounds`] error.
pub fn grid_oob_err(i: usize, j: usize, k: usize, nx: usize, ny: usize, nz: usize) -> Error {
    Error::GridIndexOutOfBounds {
        i,
        j,
        k,
        nx,
        ny,
        nz,
    }
}

/// Build a [`Error::DispatchLimitExceeded`] error.
pub fn dispatch_limit_err(dispatch_size: usize, limit: usize) -> Error {
    Error::DispatchLimitExceeded {
        dispatch_size,
        limit,
    }
}

/// Build a [`Error::ShaderCompilationError`].
pub fn shader_err(shader: impl Into<String>, message: impl Into<String>) -> Error {
    Error::ShaderCompilationError {
        shader: shader.into(),
        message: message.into(),
    }
}

/// Build an [`Error::UnsupportedFeature`].
pub fn unsupported_err(feature: impl Into<String>) -> Error {
    Error::UnsupportedFeature {
        feature: feature.into(),
    }
}

// ── Error collection ─────────────────────────────────────────────────────────

/// Collect multiple errors from a batch dispatch.  Returns `Ok(())` if the
/// vec is empty, or `Err` containing the first error otherwise.
pub fn collect_errors(errors: Vec<Error>) -> Result<()> {
    errors.into_iter().next().map_or(Ok(()), Err)
}

/// Check a boolean condition; return `Err(Error::General(msg))` if false.
pub fn check(condition: bool, msg: impl std::fmt::Display) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(Error::general(msg))
    }
}

#[cfg(test)]
mod error_tests {
    use super::*;

    #[test]
    fn test_general_error_message() {
        let e = Error::general("something went wrong");
        assert_eq!(e.to_string(), "something went wrong");
    }

    #[test]
    fn test_buffer_allocation_failed_message() {
        let e = Error::BufferAllocationFailed {
            requested_bytes: 1024,
            available_bytes: 512,
        };
        let msg = e.to_string();
        assert!(msg.contains("1024"), "should mention requested bytes");
        assert!(msg.contains("512"), "should mention available bytes");
        assert!(e.is_allocation_error());
    }

    #[test]
    fn test_invalid_buffer_handle() {
        let e = Error::InvalidBufferHandle(42);
        assert!(e.to_string().contains("42"));
    }

    #[test]
    fn test_shader_compilation_error() {
        let e = Error::ShaderCompilationError {
            shader: "sph_density".to_string(),
            message: "undefined symbol".to_string(),
        };
        let msg = e.to_string();
        assert!(msg.contains("sph_density"));
        assert!(msg.contains("undefined symbol"));
        assert!(e.is_shader_error());
    }

    #[test]
    fn test_dispatch_limit_exceeded() {
        let e = Error::DispatchLimitExceeded {
            dispatch_size: 100_000,
            limit: 65535,
        };
        let msg = e.to_string();
        assert!(msg.contains("100000"));
        assert!(msg.contains("65535"));
    }

    #[test]
    fn test_grid_index_out_of_bounds() {
        let e = Error::GridIndexOutOfBounds {
            i: 10,
            j: 5,
            k: 3,
            nx: 8,
            ny: 8,
            nz: 8,
        };
        let msg = e.to_string();
        assert!(msg.contains("10"));
        assert!(msg.contains('8'.to_string().as_str()));
    }

    #[test]
    fn test_is_not_shader_error() {
        let e = Error::general("not a shader error");
        assert!(!e.is_shader_error());
    }

    #[test]
    fn test_unsupported_feature() {
        let e = Error::UnsupportedFeature {
            feature: "ray_tracing".to_string(),
        };
        assert!(e.to_string().contains("ray_tracing"));
    }

    // ── New error variant / helper tests ─────────────────────────────────

    #[test]
    fn test_is_grid_error() {
        let e = grid_oob_err(1, 2, 3, 4, 5, 6);
        assert!(e.is_grid_error());
        assert!(!e.is_allocation_error());
    }

    #[test]
    fn test_is_arg_mismatch() {
        let e = arg_mismatch_err("test_kernel", 3, 2);
        assert!(e.is_arg_mismatch());
        assert!(!e.is_shader_error());
    }

    #[test]
    fn test_is_unsupported() {
        let e = unsupported_err("ray_tracing");
        assert!(e.is_unsupported());
    }

    #[test]
    fn test_in_stage_wraps_error() {
        let e = Error::general("boom");
        let wrapped = e.in_stage("sph_density");
        assert!(wrapped.to_string().contains("sph_density"));
        assert!(wrapped.to_string().contains("boom"));
    }

    #[test]
    fn test_alloc_err_convenience() {
        let e = alloc_err(512, 256);
        assert!(e.is_allocation_error());
        assert!(e.to_string().contains("512"));
    }

    #[test]
    fn test_dispatch_limit_err_convenience() {
        let e = dispatch_limit_err(99999, 65535);
        assert!(e.to_string().contains("99999"));
    }

    #[test]
    fn test_shader_err_convenience() {
        let e = shader_err("my_shader", "syntax error");
        assert!(e.is_shader_error());
        assert!(e.to_string().contains("syntax error"));
    }

    #[test]
    fn test_into_err() {
        let result: Result<i32> = Error::general("nope").into_err();
        assert!(result.is_err());
    }

    #[test]
    fn test_collect_errors_empty() {
        assert!(collect_errors(vec![]).is_ok());
    }

    #[test]
    fn test_collect_errors_nonempty() {
        let errs = vec![Error::general("first"), Error::general("second")];
        let result = collect_errors(errs);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("first"));
    }

    #[test]
    fn test_check_passes() {
        assert!(check(true, "should not fail").is_ok());
    }

    #[test]
    fn test_check_fails() {
        let r = check(false, "condition violated");
        assert!(r.is_err());
        assert!(r.unwrap_err().to_string().contains("condition violated"));
    }

    #[test]
    fn test_annotated_error_fatal_display() {
        let e = Error::general("crash");
        let ann = e.fatal(Some("sph_kernel"));
        let s = ann.to_string();
        assert!(s.contains("FATAL"));
        assert!(s.contains("sph_kernel"));
        assert!(s.contains("crash"));
    }

    #[test]
    fn test_annotated_error_warning_no_kernel() {
        let e = Error::general("degraded");
        let ann = e.warning(None);
        let s = ann.to_string();
        assert!(s.contains("WARNING"));
        assert!(s.contains("degraded"));
    }

    #[test]
    fn test_error_severity_ordering() {
        assert!(ErrorSeverity::Info < ErrorSeverity::Warning);
        assert!(ErrorSeverity::Warning < ErrorSeverity::Fatal);
    }

    #[test]
    fn test_error_severity_display() {
        assert_eq!(ErrorSeverity::Info.to_string(), "INFO");
        assert_eq!(ErrorSeverity::Warning.to_string(), "WARNING");
        assert_eq!(ErrorSeverity::Fatal.to_string(), "FATAL");
    }
}
