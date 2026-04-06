// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-softbody

use thiserror::Error;

/// Main error type for the softbody module
#[derive(Debug, Error)]
pub enum Error {
    /// Generic error
    #[error("{0}")]
    General(String),
}

/// Result type alias
pub type Result<T> = std::result::Result<T, Error>;
