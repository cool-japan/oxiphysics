//! # StabilityEnvelopeMonitor - Trait Implementations
//!
//! This module contains trait implementations for `StabilityEnvelopeMonitor`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::StabilityEnvelopeMonitor;

impl Default for StabilityEnvelopeMonitor {
    fn default() -> Self {
        Self {
            ssf_warning: 1.2,
            ssf_critical: 0.9,
            slip_warning: 0.2,
            slip_critical: 0.4,
            sideslip_warning_deg: 8.0,
            sideslip_critical_deg: 15.0,
        }
    }
}
