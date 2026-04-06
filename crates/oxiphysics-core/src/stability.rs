// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! API stability markers for OxiPhysics.
//!
//! This module provides types and traits to document the stability level
//! of public APIs. Consumers can query stability at runtime via the
//! [`HasStability`] trait, and library authors annotate items with the
//! marker structs ([`Stable`], [`Unstable`], [`Experimental`],
//! [`Deprecated`]) in documentation or associated types.
//!
//! # Stability Policy
//!
//! | Level          | Guarantee |
//! |----------------|-----------|
//! | **Stable**     | Follows semver strictly. Breaking changes only in major versions. |
//! | **Unstable**   | May change in minor versions with a deprecation notice. |
//! | **Experimental** | May change or be removed at any time without notice. |
//! | **Deprecated** | Scheduled for removal in a future version. |
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_core::stability::{StabilityLevel, HasStability};
//!
//! struct MyApi;
//!
//! impl HasStability for MyApi {
//!     fn stability() -> StabilityLevel {
//!         StabilityLevel::Stable
//!     }
//! }
//!
//! assert_eq!(MyApi::stability(), StabilityLevel::Stable);
//! assert_eq!(MyApi::stability().to_string(), "stable");
//! ```

/// Marker type for stable APIs.
///
/// Stable APIs follow semver strictly: breaking changes happen only
/// across major version boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Stable;

/// Marker type for unstable APIs.
///
/// Unstable APIs may change in minor releases. A deprecation notice
/// will be issued before removal whenever possible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Unstable;

/// Marker type for experimental APIs.
///
/// Experimental APIs may change or be removed at any time. They are
/// provided for early feedback and should not be relied upon in
/// production code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Experimental;

/// Marker type for deprecated APIs.
///
/// Deprecated APIs are scheduled for removal in a future version.
/// Migration guidance is provided in the item-level documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Deprecated;

/// Runtime-queryable API stability level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StabilityLevel {
    /// The API is stable and follows semver.
    Stable,
    /// The API is unstable and may change in minor versions.
    Unstable,
    /// The API is experimental and may change or be removed at any time.
    Experimental,
    /// The API is deprecated and will be removed in a future version.
    Deprecated,
}

impl StabilityLevel {
    /// Returns `true` if the API is considered safe for production use.
    ///
    /// Only [`StabilityLevel::Stable`] qualifies.
    pub fn is_production_ready(self) -> bool {
        matches!(self, Self::Stable)
    }

    /// Returns `true` if the API may change without a major version bump.
    pub fn may_change(self) -> bool {
        matches!(self, Self::Unstable | Self::Experimental)
    }

    /// Returns `true` if the API is deprecated.
    pub fn is_deprecated(self) -> bool {
        matches!(self, Self::Deprecated)
    }

    /// Returns a human-readable label suitable for documentation badges.
    pub fn badge_label(self) -> &'static str {
        match self {
            Self::Stable => "stability: stable",
            Self::Unstable => "stability: unstable",
            Self::Experimental => "stability: experimental",
            Self::Deprecated => "stability: deprecated",
        }
    }
}

impl std::fmt::Display for StabilityLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let label = match self {
            Self::Stable => "stable",
            Self::Unstable => "unstable",
            Self::Experimental => "experimental",
            Self::Deprecated => "deprecated",
        };
        f.write_str(label)
    }
}

impl std::str::FromStr for StabilityLevel {
    type Err = StabilityParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "stable" => Ok(Self::Stable),
            "unstable" => Ok(Self::Unstable),
            "experimental" => Ok(Self::Experimental),
            "deprecated" => Ok(Self::Deprecated),
            _ => Err(StabilityParseError {
                input: s.to_owned(),
            }),
        }
    }
}

/// Error returned when parsing an invalid stability level string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StabilityParseError {
    /// The invalid input string.
    pub input: String,
}

impl std::fmt::Display for StabilityParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "unknown stability level '{}': expected one of stable, unstable, experimental, deprecated",
            self.input
        )
    }
}

impl std::error::Error for StabilityParseError {}

/// Trait for types that declare their API stability level.
///
/// Implement this on public API entry-point types so consumers can
/// programmatically query how stable an API is.
///
/// ```no_run
/// use oxiphysics_core::stability::{StabilityLevel, HasStability};
///
/// struct RigidBodySolver;
///
/// impl HasStability for RigidBodySolver {
///     fn stability() -> StabilityLevel {
///         StabilityLevel::Stable
///     }
/// }
/// ```
pub trait HasStability {
    /// Returns the stability level of this API.
    fn stability() -> StabilityLevel;

    /// Returns `true` if this API is production-ready.
    fn is_production_ready() -> bool {
        Self::stability().is_production_ready()
    }

    /// Returns `true` if this API may change without a major version bump.
    fn may_change() -> bool {
        Self::stability().may_change()
    }
}

/// Returns the version of the `oxiphysics-core` crate at build time.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Returns the minimum supported Rust version (MSRV) for OxiPhysics.
pub fn msrv() -> &'static str {
    "nightly (edition 2024)"
}

/// Returns a summary of the stability policy as a static string slice.
pub fn stability_policy_summary() -> &'static str {
    "OxiPhysics follows semver for Stable APIs. \
     Unstable APIs may change in minor versions. \
     Experimental APIs may change at any time. \
     Deprecated APIs will be removed in a future major version."
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_stable() {
        assert_eq!(StabilityLevel::Stable.to_string(), "stable");
    }

    #[test]
    fn display_unstable() {
        assert_eq!(StabilityLevel::Unstable.to_string(), "unstable");
    }

    #[test]
    fn display_experimental() {
        assert_eq!(StabilityLevel::Experimental.to_string(), "experimental");
    }

    #[test]
    fn display_deprecated() {
        assert_eq!(StabilityLevel::Deprecated.to_string(), "deprecated");
    }

    #[test]
    fn parse_round_trip() {
        for level in [
            StabilityLevel::Stable,
            StabilityLevel::Unstable,
            StabilityLevel::Experimental,
            StabilityLevel::Deprecated,
        ] {
            let parsed: StabilityLevel = level.to_string().parse().expect("round-trip parse");
            assert_eq!(parsed, level);
        }
    }

    #[test]
    fn parse_case_insensitive() {
        assert_eq!(
            "STABLE".parse::<StabilityLevel>().expect("upper case"),
            StabilityLevel::Stable
        );
        assert_eq!(
            "Experimental"
                .parse::<StabilityLevel>()
                .expect("mixed case"),
            StabilityLevel::Experimental
        );
    }

    #[test]
    fn parse_invalid() {
        let err = "bogus".parse::<StabilityLevel>().expect_err("should fail");
        assert_eq!(err.input, "bogus");
        assert!(err.to_string().contains("bogus"));
    }

    #[test]
    fn production_ready() {
        assert!(StabilityLevel::Stable.is_production_ready());
        assert!(!StabilityLevel::Unstable.is_production_ready());
        assert!(!StabilityLevel::Experimental.is_production_ready());
        assert!(!StabilityLevel::Deprecated.is_production_ready());
    }

    #[test]
    fn may_change() {
        assert!(!StabilityLevel::Stable.may_change());
        assert!(StabilityLevel::Unstable.may_change());
        assert!(StabilityLevel::Experimental.may_change());
        assert!(!StabilityLevel::Deprecated.may_change());
    }

    #[test]
    fn is_deprecated() {
        assert!(StabilityLevel::Deprecated.is_deprecated());
        assert!(!StabilityLevel::Stable.is_deprecated());
    }

    #[test]
    fn badge_labels() {
        assert_eq!(StabilityLevel::Stable.badge_label(), "stability: stable");
        assert_eq!(
            StabilityLevel::Experimental.badge_label(),
            "stability: experimental"
        );
    }

    #[test]
    fn version_not_empty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn msrv_contains_nightly() {
        assert!(msrv().contains("nightly"));
    }

    #[test]
    fn policy_summary_not_empty() {
        assert!(!stability_policy_summary().is_empty());
    }

    #[test]
    fn has_stability_trait() {
        struct TestApi;
        impl HasStability for TestApi {
            fn stability() -> StabilityLevel {
                StabilityLevel::Unstable
            }
        }
        assert_eq!(TestApi::stability(), StabilityLevel::Unstable);
        assert!(!TestApi::is_production_ready());
        assert!(TestApi::may_change());
    }

    #[test]
    fn marker_structs_are_eq() {
        assert_eq!(Stable, Stable);
        assert_eq!(Unstable, Unstable);
        assert_eq!(Experimental, Experimental);
        assert_eq!(Deprecated, Deprecated);
    }
}
