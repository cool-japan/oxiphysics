//! # ConstraintKind - Trait Implementations
//!
//! This module contains trait implementations for `ConstraintKind`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::ConstraintKind;

impl std::fmt::Display for ConstraintKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConstraintKind::Fixed => write!(f, "Fixed"),
            ConstraintKind::Ball => write!(f, "Ball"),
            ConstraintKind::Revolute => write!(f, "Revolute"),
            ConstraintKind::Prismatic => write!(f, "Prismatic"),
            ConstraintKind::Spring => write!(f, "Spring"),
            ConstraintKind::Contact => write!(f, "Contact"),
            ConstraintKind::SixDof => write!(f, "SixDof"),
            ConstraintKind::Gear => write!(f, "Gear"),
            ConstraintKind::Pulley => write!(f, "Pulley"),
            ConstraintKind::RackPinion => write!(f, "RackPinion"),
            ConstraintKind::Motor => write!(f, "Motor"),
            ConstraintKind::Pbd => write!(f, "Pbd"),
            ConstraintKind::Custom(id) => write!(f, "Custom({id})"),
        }
    }
}
