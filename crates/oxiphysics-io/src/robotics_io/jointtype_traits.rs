//! # JointType - Trait Implementations
//!
//! This module contains trait implementations for `JointType`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::JointType;
use std::fmt;

impl fmt::Display for JointType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            JointType::Fixed => "fixed",
            JointType::Revolute => "revolute",
            JointType::Prismatic => "prismatic",
            JointType::Continuous => "continuous",
            JointType::Planar => "planar",
            JointType::Floating => "floating",
        };
        write!(f, "{s}")
    }
}
