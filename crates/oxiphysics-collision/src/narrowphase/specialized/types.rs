//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use crate::types::Contact;

/// Precise cylinder-cylinder contact test with mode detection.
///
/// Detects three contact modes:
/// 1. End-cap vs end-cap
/// 2. Side vs side (parallel axes)
/// 3. Side vs end-cap (crossing axes)
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum CylinderContactMode {
    /// Two side surfaces are closest.
    SideToSide,
    /// An end-cap of one cylinder is closest to the side of the other.
    CapToSide,
    /// Two end-caps are closest.
    CapToCap,
    /// Axes are parallel; side-to-side with full overlap.
    Parallel,
}
/// Result of a precise cylinder-cylinder contact query.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CylinderContactResult {
    /// Contact information (if any).
    pub contact: Option<Contact>,
    /// Contact mode detected.
    pub mode: CylinderContactMode,
}
