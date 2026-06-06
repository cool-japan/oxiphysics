//! # MsgType - Trait Implementations
//!
//! This module contains trait implementations for `MsgType`.
//!
//! ## Implemented Traits
//!
//! - `TryFrom`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MsgType;

impl TryFrom<u8> for MsgType {
    type Error = ();
    fn try_from(v: u8) -> Result<Self, ()> {
        match v {
            1 => Ok(MsgType::JointState),
            2 => Ok(MsgType::Wrench),
            3 => Ok(MsgType::GripperState),
            4 => Ok(MsgType::PointCloud),
            255 => Ok(MsgType::Generic),
            _ => Err(()),
        }
    }
}
