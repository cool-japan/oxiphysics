//! # IpcBarrierParams - Trait Implementations
//!
//! This module contains trait implementations for `IpcBarrierParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::IpcBarrierParams;

impl Default for IpcBarrierParams {
    fn default() -> Self {
        IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        }
    }
}
