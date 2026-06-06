//! # SpectrumSource - Trait Implementations
//!
//! This module contains trait implementations for `SpectrumSource`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SpectrumSource;

impl Default for SpectrumSource {
    fn default() -> Self {
        SpectrumSource::Unknown(String::new())
    }
}

impl std::fmt::Display for SpectrumSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpectrumSource::Ftir => write!(f, "FTIR"),
            SpectrumSource::Nir => write!(f, "NIR"),
            SpectrumSource::Raman => write!(f, "Raman"),
            SpectrumSource::UvVis => write!(f, "UV-Vis"),
            SpectrumSource::MassSpec => write!(f, "MS"),
            SpectrumSource::Nmr => write!(f, "NMR"),
            SpectrumSource::Unknown(s) => write!(f, "{}", s),
        }
    }
}
