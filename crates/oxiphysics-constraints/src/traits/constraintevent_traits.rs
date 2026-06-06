//! # ConstraintEvent - Trait Implementations
//!
//! This module contains trait implementations for `ConstraintEvent`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::ConstraintEvent;

impl std::fmt::Display for ConstraintEvent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConstraintEvent::Broken {
                constraint_index,
                force,
            } => {
                write!(f, "Broken(idx={constraint_index}, force={force:.3e})")
            }
            ConstraintEvent::ContactBegin { constraint_index } => {
                write!(f, "ContactBegin(idx={constraint_index})")
            }
            ConstraintEvent::ContactEnd { constraint_index } => {
                write!(f, "ContactEnd(idx={constraint_index})")
            }
            ConstraintEvent::LimitHit {
                constraint_index,
                is_lower,
            } => {
                write!(f, "LimitHit(idx={constraint_index}, lower={is_lower})")
            }
            ConstraintEvent::ConvergenceWarning {
                constraint_index,
                residual,
            } => {
                write!(
                    f,
                    "ConvergenceWarning(idx={constraint_index}, res={residual:.3e})"
                )
            }
        }
    }
}
