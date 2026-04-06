#[allow(unused_imports)]
use super::functions::*;
// # PhysicsWorld - Trait Implementations
//
// This module contains trait implementations for `PhysicsWorld`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::PhysicsConfig;
use std::collections::VecDeque;

use super::types::{BodyArena, PhysicsWorld};

impl Default for PhysicsWorld {
    fn default() -> Self {
        Self {
            time: 0.0,
            config: PhysicsConfig::default(),
            fixed_dt: 1.0 / 60.0,
            accumulator: 0.0,
            bodies: BodyArena::new(),
            constraints: Vec::new(),
            events: VecDeque::new(),
            islands: Vec::new(),
            max_events: 1024,
        }
    }
}
