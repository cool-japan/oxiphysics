//! # NoBarostat - Trait Implementations
//!
//! This module contains trait implementations for `NoBarostat`.
//!
//! ## Implemented Traits
//!
//! - `Barostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::Barostat;
use super::types::NoBarostat;

impl Barostat for NoBarostat {
    fn apply(
        &mut self,
        _atoms: &mut AtomSet,
        _pbox: &mut PeriodicBox,
        _target_pressure: f64,
        _dt: f64,
    ) {
    }
}
