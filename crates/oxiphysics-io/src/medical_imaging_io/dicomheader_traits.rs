//! # DicomHeader - Trait Implementations
//!
//! This module contains trait implementations for `DicomHeader`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DicomHeader;
use std::fmt;

impl fmt::Display for DicomHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "DicomHeader(patient={}, modality={}, {}x{}, slice={:.2}mm)",
            self.patient_name, self.modality, self.columns, self.rows, self.slice_thickness
        )
    }
}
