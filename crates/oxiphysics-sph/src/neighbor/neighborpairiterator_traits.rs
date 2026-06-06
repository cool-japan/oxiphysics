//! # NeighborPairIterator - Trait Implementations
//!
//! This module contains trait implementations for `NeighborPairIterator`.
//!
//! ## Implemented Traits
//!
//! - `Iterator`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::NeighborPairIterator;

impl<'a> Iterator for NeighborPairIterator<'a> {
    type Item = (usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        while self.i < self.lists.len() {
            let nbrs = &self.lists[self.i];
            while self.j_idx < nbrs.len() {
                let j = nbrs[self.j_idx];
                self.j_idx += 1;
                if j > self.i {
                    return Some((self.i, j));
                }
            }
            self.i += 1;
            self.j_idx = 0;
        }
        None
    }
}
